
#####################################
#
# author: Huzaifa Hassan
# Date: 07/17/2023
# Script to run gRNA design, parse RNAplfold output results and split each FASTA sequence into
# sliding windows of the gRNA length (set with -l; default 22).
# The scores are first averaged by row and then again averaged for the scores of individual bases
# of each window. (Formula designed by Ariel Bazzini)
#
# run 'python gRNADesign.py -h' for details
#
# 02/23: Added parallelization, optimized the code and added more documentation
# 11/21: Optimized the code for mismatch analysis
# 08/26: Fixed off-by-one in sliding window, hardcoded line-index parsing in
#         parse_water_align/water_out_parse, global variable scoping for output dirs,
#         sb.call → sb.run with error handling in run_RNAplfold_top, RNAplfold temp-file isolation
#         via tempfile.TemporaryDirectory, input validation in read_delimited_file, type hints,
#         pathlib throughout, logging instead of print
# 10/26: Water parse failures now raise instead of defaulting; alignment blocks paired by
#         position (empty blocks no longer shift pairs); '-m 0' honored; empty alignment results
#         handled; identity/similarity_perc '%' stripped (matches pre-08/26 output); missing
#         RNAplfold output raises; FASTA path shell-quoted; unused imports removed;
#         transcript IDs containing '|' preserved (EMBOSS renamed them to EMBOSS_###);
#         gRNA_gaps now counts gaps in the gRNA only (was header total gaps), cross-checked
#         against header; water stderr (e.g. skipped zero-length records) logged once;
#         window scoring vectorized (pandas rolling mean, output identical after rounding);
#         added seq_<L>nt_revcomp column (target site and its reverse complement both
#         reported; target stays at column index 1 for the mismatch module);
#         added --fold_mode {local,global} and --winsize. global uses RNAfold -p over the
#         whole sequence and reproduces the ViennaRNA RNAfold web server readout exactly
#         (verified 1532/1532 windows identical on cth1). local is unchanged and remains
#         the default.
#####################################

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.Emboss.Applications import WaterCommandline
import subprocess as sb
import pandas as pd
import argparse
import re
import sys
import multiprocessing as mp
import shutil
import shlex
import tempfile
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Output columns of parse_water_align (kept in one place so empty results keep the schema)
WATER_COLUMNS = [
    'gRNA_name', 'seq2_name', 'aligned_length', 'identity', 'identity_perc',
    'similarity', 'similarity_perc', 'gRNA_align_start', 'gRNA_align_end',
    'seq_align_start', 'seq_align_end', 'gRNA_gaps', 'seq_gaps', 'total_mismatches',
]


def mkdir(subcommand: str):
    """Create output directories based on subcommand. Returns Path(s) for created dirs."""
    results = Path("results")
    results.mkdir(exist_ok=True)

    if subcommand == "design":
        grna_dir = results / "grna"
        grna_dir.mkdir(exist_ok=True)
        return grna_dir
    elif subcommand == "mismatch":
        all_align_dir = results / "all_alignments"
        all_align_dir.mkdir(exist_ok=True)
        mismatches_dir = results / "mismatches"
        mismatches_dir.mkdir(exist_ok=True)
        return all_align_dir, mismatches_dir
    else:
        grna_dir = results / "grna"
        grna_dir.mkdir(exist_ok=True)
        all_align_dir = results / "all_alignments"
        all_align_dir.mkdir(exist_ok=True)
        mismatches_dir = results / "mismatches"
        mismatches_dir.mkdir(exist_ok=True)
        return grna_dir, all_align_dir, mismatches_dir


# Dot-plot 'ubox' line: "i j sqrt(P(i,j)) ubox"
DP_UBOX_RE = re.compile(r'^(\d+)\s+(\d+)\s+([\d.eE+-]+)\s+ubox')

# Above this length a global fold gets slow (RNAfold is ~O(n^3) time, O(n^2) memory)
GLOBAL_FOLD_WARN_LEN = 5000


def global_unpaired_probs(seq_fasta_abs: str, tmpdir: str, seq_lengths: dict) -> dict:
    """
    Run RNAfold -p (global partition function) and return, per sequence, the per-base
    probability of being unpaired.

    This is the quantity the ViennaRNA RNAfold web server uses to colour the structure:
    "For unpaired regions the color denotes the probability of being unpaired." Using it
    makes the design score reproduce the web server's readout locally.

    unpaired[i] = 1 - sum_j P(i pairs with j), read from the <id>_dp.ps dot plot, where
    each 'ubox' entry stores sqrt(P).

    Parameters:
        seq_fasta_abs: Absolute path to the input FASTA.
        tmpdir: Directory to run in (keeps .ps files isolated).
        seq_lengths: {safe_id: sequence length}.

    Returns:
        {safe_id: [per-base unpaired probability, 1-based order]}
    """
    if shutil.which("RNAfold") is None:
        raise FileNotFoundError("RNAfold (ViennaRNA) is not installed or not in the PATH.")

    long_seqs = {k: v for k, v in seq_lengths.items() if v > GLOBAL_FOLD_WARN_LEN}
    if long_seqs:
        logger.warning(
            f"Global fold requested for {len(long_seqs)} sequence(s) longer than "
            f"{GLOBAL_FOLD_WARN_LEN} nt (longest {max(long_seqs.values())} nt). "
            "RNAfold scales ~O(n^3) in time and O(n^2) in memory; this may be slow."
        )

    cmd = f"RNAfold -p --noPS < {shlex.quote(seq_fasta_abs)}"
    try:
        sb.run(cmd, shell=True, check=True, stdout=sb.DEVNULL, stderr=sb.PIPE,
               cwd=tmpdir, text=True)
    except sb.CalledProcessError as e:
        raise RuntimeError(f"RNAfold failed:\n{e.stderr}") from e

    probs = {}
    for safe_id, n in seq_lengths.items():
        dp_path = Path(tmpdir) / (safe_id + "_dp.ps")
        if not dp_path.exists():
            produced = sorted(p.name for p in Path(tmpdir).glob("*_dp.ps"))
            raise FileNotFoundError(
                f"RNAfold dot plot '{dp_path.name}' not found. RNAfold produced: "
                f"{produced[:10]}{' ...' if len(produced) > 10 else ''}. "
                "Check FASTA header characters."
            )

        paired = [0.0] * (n + 1)
        with dp_path.open() as fh:
            for line in fh:
                m = DP_UBOX_RE.match(line.strip())
                if m:
                    i, j, p = int(m.group(1)), int(m.group(2)), float(m.group(3)) ** 2
                    if i <= n:
                        paired[i] += p
                    if j <= n:
                        paired[j] += p

        # A base cannot pair with total probability > 1. If it does, we are reading the
        # wrong file - RNAplfold writes a <id>_dp.ps with the same name and would silently
        # overwrite RNAfold's. Fail loudly rather than score on invalid numbers.
        worst = max(paired[1:], default=0.0)
        if worst > 1.0 + 1e-6:
            raise ValueError(
                f"Dot plot for '{safe_id}' gives total pairing probability {worst:.4f} for a "
                "single base (maximum is 1.0). The file is not a valid RNAfold partition-"
                "function dot plot - it may have been overwritten by another ViennaRNA tool."
            )
        probs[safe_id] = [min(1.0, max(0.0, 1.0 - x)) for x in paired[1:]]

    return probs


def run_RNAplfold_top(seq_fasta: str, top_N: int = None, len_grna: int = 23,
                      fold_mode: str = "local", winsize: int = 70) -> pd.DataFrame:
    """
    Design gRNAs by scoring sliding windows by average unpaired probability.

    Two folding modes:
      local  (default) - RNAplfold with a sliding window of `winsize` nt. Base pairs are
                         limited to that span. This is the original pipeline behaviour.
      global           - RNAfold partition function over the whole sequence, using the
                         per-base unpaired probability. Reproduces the ViennaRNA RNAfold
                         web server readout (verified: Spearman 1.000 against the web
                         server quantity on real zebrafish transcripts).

    The two modes rank sites differently on long mRNAs, because global folding allows
    long-range base pairs that a local window cannot see. Neither is "correct" a priori.

    Parameters:
        seq_fasta: Path to input FASTA file.
        top_N: Number of top-scoring gRNAs to retain per sequence (None = all).
        len_grna: gRNA length in bp.
        fold_mode: "local" or "global".
        winsize: RNAplfold window size (local mode only).

    Returns:
        DataFrame of scored gRNAs with base composition columns.
    """
    if sys.version_info[0] < 3:
        raise Exception("Please use Python 3.7 or higher.")
    if fold_mode not in ("local", "global"):
        raise ValueError(f"fold_mode must be 'local' or 'global', got {fold_mode!r}")

    # Resolve to absolute path before changing working dir inside tempdir
    seq_fasta_abs = str(Path(seq_fasta).resolve())

    # Materialized once: needed for sequence lengths in global mode, and the design input
    # (genes of interest) is normally small
    records = list(SeqIO.parse(seq_fasta, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in '{seq_fasta}'.")

    all_top_rows_df = []

    # Run the folding tool in a temp dir to keep _lunp / _dp.ps files isolated
    with tempfile.TemporaryDirectory() as tmpdir:
        if fold_mode == "global":
            seq_lengths = {r.id.replace('|', '_'): len(r.seq) for r in records}
            global_probs = global_unpaired_probs(seq_fasta_abs, tmpdir, seq_lengths)
            logger.info(
                f"Global fold (RNAfold -p): scoring by per-base unpaired probability "
                f"over {len(records)} sequence(s) - matches the RNAfold web server."
            )
        else:
            # Quote the path so spaces/metacharacters don't break the shell redirect
            cmd = f"RNAplfold -W {winsize} -u {len_grna} < {shlex.quote(seq_fasta_abs)}"
            try:
                sb.run(cmd, shell=True, check=True, stderr=sb.PIPE, cwd=tmpdir, text=True)
            except sb.CalledProcessError as e:
                raise RuntimeError(f"RNAplfold failed:\n{e.stderr}") from e

        for seq_record in records:
            seq_id = seq_record.id
            # Replace | with _ to match the filename RNAplfold/RNAfold write
            safe_id = seq_id.replace('|', '_')
            fasta = seq_record.seq

            if fold_mode == "global":
                averages_list = global_probs[safe_id]
            else:
                lunp_path = Path(tmpdir) / (safe_id + "_lunp")
                # Fail loudly: silently skipping would drop genes from the design output
                if not lunp_path.exists():
                    produced = sorted(p.name for p in Path(tmpdir).glob("*_lunp"))
                    raise FileNotFoundError(
                        f"RNAplfold output '{lunp_path.name}' not found for sequence "
                        f"'{seq_id}'. RNAplfold produced: {produced[:10]}"
                        f"{' ...' if len(produced) > 10 else ''}. "
                        "Check FASTA header characters."
                    )

                df_lunp = pd.read_csv(lunp_path, skiprows=1, sep="\t", header=0)
                # Row-average the unpaired probability columns (one per position)
                df_lunp['average'] = df_lunp.iloc[:, 1:len_grna + 1].mean(axis=1)
                averages_list = list(df_lunp['average'])

            # Sliding windows — fixed off-by-one (was: range(len - len_grna - 1)).
            # Bounded by both the sequence and the lunp table, as the old zip() did.
            n_win = max(0, min(len(fasta), len(averages_list)) - len_grna + 1)
            if n_win == 0:
                logger.warning(
                    f"Sequence '{seq_id}' is shorter than the gRNA length ({len(fasta)} < "
                    f"{len_grna} nt) — no gRNAs designed for it."
                )

            # Vectorized sliding-window mean. Identical to statistics.mean after round(,3):
            # verified over 52,536 windows of real transcripts, max raw difference 1.1e-16.
            # pandas' rolling mean uses compensated summation, so it does not drift on long
            # sequences the way a cumsum-based window would.
            scores = pd.Series(averages_list).rolling(len_grna).mean().dropna().to_numpy()

            # Both strands are reported: the target site as it appears in the input FASTA,
            # and its reverse complement. reverse_complement() follows the input alphabet
            # (DNA stays T, RNA stays U) and preserves case.
            # Column ORDER MATTERS: read_delimited_file() takes column 0 as the name and
            # column 1 as the sequence to align, so the target site must stay at index 1
            # and the reverse complement goes after it.
            rows_all = [
                [
                    safe_id,
                    str(fasta[i:i + len_grna]),
                    str(fasta[i:i + len_grna].reverse_complement()),
                    str(i + 1),                 # 1-based start
                    str(i + len_grna),          # 1-based end
                    str(round(float(scores[i]), 3)),
                ]
                for i in range(n_win)
            ]

            seq_col = f"seq_{len_grna}nt"
            rc_col = f"seq_{len_grna}nt_revcomp"
            df_rows = pd.DataFrame(
                rows_all,
                columns=["Gene|Transcript", seq_col, rc_col,
                         "Seq_start", "Seq_end", "Average_score"]
            )
            df_rows_sort = df_rows.sort_values(by='Average_score', ascending=False)
            df_rows_sort['gRNA_name'] = (
                df_rows_sort["Gene|Transcript"] + "_"
                + df_rows_sort['Seq_start'] + "_"
                + df_rows_sort['Seq_end']
            )

            all_top_rows_df.append(df_rows_sort.head(top_N) if top_N else df_rows_sort)

    result = pd.concat(all_top_rows_df)
    seq_col = f"seq_{len_grna}nt"
    rc_col = f"seq_{len_grna}nt_revcomp"
    # Base composition is reported for the target site (column seq_col), as before.
    # For the reverse complement these swap: G%<->C% and A%<->T%.
    result['G_Percentage'] = round(result[seq_col].str.count('G', re.I) * 100 / len_grna, 3)
    result['C_Percentage'] = round(result[seq_col].str.count('C', re.I) * 100 / len_grna, 3)
    result['A_Percentage'] = round(result[seq_col].str.count('A', re.I) * 100 / len_grna, 3)
    result['T_Percentage'] = round(result[seq_col].str.count('T', re.I) * 100 / len_grna, 3)

    return result[[
        'gRNA_name', seq_col, rc_col, "Gene|Transcript",
        'Seq_start', 'Seq_end', 'Average_score',
        'G_Percentage', 'C_Percentage', 'A_Percentage', 'T_Percentage'
    ]]


def water_cmdline(g_fasta: str, tran_fasta_fh: str) -> tuple:
    """
    Run the EMBOSS water commandline for sequence alignment.

    Parameters:
        g_fasta: gRNA sequence as a string.
        tran_fasta_fh: File path to the FASTA file for the cDNA sequence.

    Returns:
        (stdout, stderr): alignment text, and water's warnings (e.g. skipped zero-length
        records) so the caller can report them instead of silently discarding them.

    Raises:
        FileNotFoundError: If the water tool is not installed or not in PATH.
        subprocess.CalledProcessError: If the water command fails.
    """
    if shutil.which("water") is None:
        raise FileNotFoundError("The 'water' tool from EMBOSS is not installed or not in the PATH.")

    try:
        water_cmd = WaterCommandline(gapopen=10, gapextend=0.5)
        water_cmd.asequence = "asis:" + g_fasta
        water_cmd.bsequence = tran_fasta_fh
        water_cmd.outfile = "stdout"

        result = sb.run(
            str(water_cmd),
            stdout=sb.PIPE,
            stderr=sb.PIPE,
            text=True,
            shell=True,
            check=True
        )
        return result.stdout, result.stderr

    except sb.CalledProcessError as e:
        logger.error(f"Error running water command: {e.stderr}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in water_cmdline: {e}")
        raise


def parse_water_align(text: list, g_name: str, len_grna: int = 23) -> dict:
    """
    Parse a single EMBOSS water alignment block using regex.
    Replaces brittle hardcoded line-index access with pattern matching.

    Parameters:
        text: List of lines from one alignment block (header + alignment lines).
        g_name: gRNA name.
        len_grna: Length of gRNA.

    Returns:
        dict of alignment metrics.
    """
    full_text = '\n'.join(text)

    # Extract header metadata — robust to line-order changes across EMBOSS versions
    seq2_m   = re.search(r'# 2:\s+(\S+)', full_text)
    len_m    = re.search(r'# Length:\s+(\d+)', full_text)
    id_m     = re.search(r'# Identity:\s+(\d+)/\d+\s+\(([^)]+)\)', full_text)
    sim_m    = re.search(r'# Similarity:\s+(\d+)/\d+\s+\(([^)]+)\)', full_text)
    gap_m    = re.search(r'# Gaps:\s+(\d+)/\d+', full_text)

    # Required header fields: raise rather than default, so parse failures can't masquerade
    # as 0-mismatch hits in the filtered output
    missing = [name for name, m in [("2", seq2_m), ("Length", len_m), ("Identity", id_m),
                                    ("Similarity", sim_m), ("Gaps", gap_m)] if m is None]
    if missing:
        raise ValueError(
            f"Could not parse water header field(s) {missing} for gRNA '{g_name}'. "
            f"Block start:\n{full_text[:500]}"
        )

    seq2_name     = seq2_m.group(1)
    length        = int(len_m.group(1))
    identity      = int(id_m.group(1))
    # Strip padding and '%' to keep numeric output format (matches pre-08/26 behavior)
    identity_perc = id_m.group(2).strip().rstrip('%')
    similarity    = int(sim_m.group(1))
    sim_perc      = sim_m.group(2).strip().rstrip('%')
    # Header '# Gaps' is the TOTAL gap count (gRNA + transcript); used below as a cross-check
    total_gaps    = int(gap_m.group(1))

    # Alignment lines (may wrap over several lines at 50 columns): (start, aligned_seq, end)
    asis_lines = re.findall(r'^asis\s+(\d+)\s+(\S+)\s+(\d+)', full_text, re.MULTILINE)
    # Subject lines: not "asis", not "#", not match lines (which start with whitespace)
    subj_lines = re.findall(
        r'^(?!asis\b|#|\s)\S+\s+(\d+)\s+(\S+)\s+(\d+)', full_text, re.MULTILINE
    )

    if asis_lines and subj_lines:
        seq1_start = int(asis_lines[0][0])
        seq1_end   = int(asis_lines[-1][2])
        seq2_start = int(subj_lines[0][0])
        seq2_end   = int(subj_lines[-1][2])
    elif length == 0:
        # Genuine empty local alignment: treat the whole gRNA as unaligned
        # (formula below then yields total_mismatches = len_grna)
        seq1_start, seq1_end = 1, 0
        seq2_start, seq2_end = None, None
    else:
        raise ValueError(
            f"water reports Length={length} for gRNA '{g_name}' vs '{seq2_name}' "
            "but alignment coordinate lines could not be parsed."
        )
    # Gaps per sequence, counted from the aligned sequence strings themselves
    grna_gaps = sum(s.count('-') for _, s, _ in asis_lines)
    seq_gaps  = sum(s.count('-') for _, s, _ in subj_lines)
    # Cross-check against water's own total; a mismatch means the alignment lines were misparsed
    if grna_gaps + seq_gaps != total_gaps:
        raise ValueError(
            f"Gap count mismatch for gRNA '{g_name}' vs '{seq2_name}': header Gaps={total_gaps}, "
            f"parsed gRNA_gaps={grna_gaps} + seq_gaps={seq_gaps}."
        )

    # Mismatch calculation: account for unaligned gRNA overhangs at start/end
    add_start       = seq1_start - 1
    add_end         = len_grna - seq1_end
    added_aligned   = length + add_start + add_end
    total_mismatches = added_aligned - identity

    return {
        'gRNA_name':        g_name,
        'seq2_name':        seq2_name,
        'aligned_length':   length,
        'identity':         identity,
        'identity_perc':    identity_perc,
        'similarity':       similarity,
        'similarity_perc':  sim_perc,
        'gRNA_align_start': seq1_start,
        'gRNA_align_end':   seq1_end,
        'seq_align_start':  seq2_start,
        'seq_align_end':    seq2_end,
        'gRNA_gaps':        grna_gaps,
        'seq_gaps':         seq_gaps,
        'total_mismatches': total_mismatches,
    }


def water_out_parse(waterlines: str, g_name: str, len_grna: int = 23) -> pd.DataFrame:
    """
    Parse multi-alignment water stdout into a DataFrame.
    Replaces hardcoded line-index skipping with regex-based header detection.

    Parameters:
        waterlines: Full stdout string from the water command.
        g_name: gRNA name.
        len_grna: Length of gRNA.

    Returns:
        DataFrame with one row per transcript alignment.
    """
    # Strip global program header (the ####... block at the top)
    body = re.split(r'#{40,}.*?#{40,}', waterlines, maxsplit=1, flags=re.DOTALL)[-1]

    # Split remaining text into blocks by the #=== separators. Do NOT drop empty blocks:
    # pairing must stay positional, or one empty alignment shifts every later pair.
    raw_blocks = re.split(r'#={39,}', body)

    # Each metadata block (contains '# Aligned_sequences') is followed by its alignment block
    records = []
    for i, block in enumerate(raw_blocks):
        if '# Aligned_sequences' not in block:
            continue
        aln = raw_blocks[i + 1] if i + 1 < len(raw_blocks) else ""
        combined_lines = (block + "\n" + aln).split('\n')
        records.append(parse_water_align(combined_lines, g_name, len_grna))

    # Fixed columns so an empty result still has 'total_mismatches' for downstream filtering
    df = pd.DataFrame(records, columns=WATER_COLUMNS)
    # Nullable int keeps coordinates integer when an empty alignment contributes None
    df[['seq_align_start', 'seq_align_end']] = df[['seq_align_start', 'seq_align_end']].astype("Int64")
    return df


def run_water_parallel(args: tuple) -> tuple:
    """Wrapper for multiprocessing — aligns one gRNA against the transcript FASTA.
    Returns (gRNA name, alignment DataFrame, water stderr)."""
    g_name, g_seq, transSeq_fh = args
    stdout, stderr = water_cmdline(g_seq, transSeq_fh)
    return g_name, water_out_parse(stdout, g_name, len_grna=len(g_seq)), stderr


def prepare_transcripts(transcript_fasta: str, workdir: Path) -> tuple:
    """
    Scan the transcript FASTA once before alignment.

    - EMBOSS replaces IDs containing '|' with EMBOSS_001, EMBOSS_002, ... so hits can't be
      traced back. If any ID contains '|', write a copy with '|' -> '_' into workdir and
      return a map to restore the original IDs in the output.
    - Zero-length records are skipped by water; report them once here.

    Returns:
        (path to FASTA for water, {water-visible ID: original ID})
    """
    try:
        records = list(SeqIO.parse(transcript_fasta, "fasta"))
    except (OSError, ValueError) as e:
        raise RuntimeError(f"Could not read transcript FASTA '{transcript_fasta}': {e}") from e
    if not records:
        raise ValueError(f"No FASTA records found in transcript file '{transcript_fasta}'.")

    empty = [r.id for r in records if len(r.seq) == 0]
    if empty:
        logger.warning(
            f"{len(empty)} zero-length transcript record(s) will be skipped by water "
            f"(no rows in output): {empty[:10]}{' ...' if len(empty) > 10 else ''}"
        )

    piped = [r.id for r in records if '|' in r.id]
    if not piped:
        return transcript_fasta, {}

    # Sanitize IDs; refuse if sanitizing would merge two different transcripts
    original_ids = {r.id for r in records}
    name_map = {}
    for r in records:
        safe = r.id.replace('|', '_')
        # Collision with an existing ID, or two piped IDs sanitizing to the same name
        if (safe != r.id and safe in original_ids) or name_map.get(safe, r.id) != r.id:
            raise ValueError(
                f"Transcript ID '{r.id}' becomes '{safe}' after replacing '|', which collides "
                "with an existing ID. Rename one of them in the transcript FASTA."
            )
        name_map[safe] = r.id

    clean_path = workdir / "transcripts_sanitized.fa"
    try:
        with clean_path.open("w") as fh:
            for r in records:
                fh.write(f">{r.id.replace('|', '_')}\n{str(r.seq)}\n")
    except OSError as e:
        raise RuntimeError(f"Could not write sanitized transcript FASTA to '{clean_path}': {e}") from e

    logger.info(f"{len(piped)} transcript ID(s) contain '|'; aligned with '|'->'_' and restored in output.")
    return str(clean_path), name_map


def read_delimited_file(file_path: str, sep: str = "\t", header: int = 0) -> dict:
    """
    Read a delimited file and return a {name: sequence} dictionary.
    Expects at least two columns: name in column 0, sequence in column 1.
    """
    df = pd.read_csv(file_path, sep=sep, header=header)
    if df.shape[1] < 2:
        raise ValueError(
            f"Expected at least 2 columns in {file_path}, got {df.shape[1]}. "
            "Column 0 = gRNA name, column 1 = sequence."
        )
    return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))


def read_fasta_file(fasta_fh: str) -> dict:
    """Read a FASTA file and return a {id: sequence} dictionary."""
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_fh, "fasta")}


def run_mismatch(
    args: argparse.Namespace,
    grnafile: str = None,
    all_align_dir: Path = None,
    mismatches_dir: Path = None,
) -> None:
    """
    Run parallelized mismatch analysis for all gRNAs against a transcript FASTA.

    Parameters:
        args: Parsed argparse namespace.
        grnafile: Optional path override for gRNA input file (used when called from 'all').
        all_align_dir: Output directory for all alignment CSVs.
        mismatches_dir: Output directory for mismatch-filtered CSVs.
    """
    grna_input   = grnafile or args.grna_file
    # 'is not None' so an explicit '-m 0' (perfect matches only) is honored
    mismatches_n = int(args.mismatches) if args.mismatches is not None else 5
    num_processes = int(args.num_processes) if args.num_processes else 1

    # Parse gRNA input file
    if grna_input.endswith(".txt"):
        grna_dict = read_delimited_file(grna_input, sep="\t", header=0)
    elif grna_input.endswith(".csv"):
        grna_dict = read_delimited_file(grna_input, sep=",", header=0)
    elif grna_input.endswith((".fasta", ".fa")):
        grna_dict = read_fasta_file(grna_input)
    else:
        logger.error(f"Unsupported file format for gRNA file: {grna_input}")
        sys.exit(1)

    # Sanitized copy (if needed) lives under results/ rather than /tmp, since a
    # transcriptome FASTA can be large; removed automatically when alignment finishes
    with tempfile.TemporaryDirectory(dir=Path(mismatches_dir).parent, prefix=".tmp_transcripts_") as tmpdir:
        transSeq_fh, name_map = prepare_transcripts(args.transcript_fasta, Path(tmpdir))
        try:
            with mp.Pool(processes=num_processes) as pool:
                results = pool.map(
                    run_water_parallel,
                    [(name, seq.replace('U', 'T'), transSeq_fh) for name, seq in grna_dict.items()]
                )
        except Exception as e:
            logger.error(f"Error during parallel alignment: {e}")
            raise

    # Report each distinct water warning once (they repeat identically for every gRNA)
    warn_counts = {}
    for _, _, stderr in results:
        for line in filter(None, (l.strip() for l in stderr.splitlines())):
            warn_counts[line] = warn_counts.get(line, 0) + 1
    # EMBOSS prefixes real problems with Warning/Error/Died; other stderr (e.g. the program
    # banner) is informational and goes to DEBUG
    problems = {l: n for l, n in warn_counts.items() if l.startswith(("Warning", "Error", "Died"))}
    for line, n in warn_counts.items():
        if line not in problems:
            logger.debug(f"water stderr ({n}/{len(results)} gRNA runs): {line}")
    for i, (line, n) in enumerate(problems.items()):
        if i == 20:
            logger.warning(f"... {len(problems) - 20} more distinct water warning(s) not shown.")
            break
        logger.warning(f"water ({n}/{len(results)} gRNA runs): {line}")

    known_ids = set(name_map.values()) | set(name_map.keys())
    for g_name, result_df, _ in results:
        # Restore original transcript IDs that contained '|'
        if name_map:
            result_df['seq2_name'] = result_df['seq2_name'].map(lambda s: name_map.get(s, s))
        # Any remaining EMBOSS_### placeholder means water renamed an ID we didn't anticipate
        renamed = result_df['seq2_name'].astype(str).str.fullmatch(r'EMBOSS_\d+') & \
            ~result_df['seq2_name'].isin(known_ids)
        if renamed.any():
            logger.warning(
                f"gRNA '{g_name}': {int(renamed.sum())} hit(s) have water placeholder names "
                "(EMBOSS_###); check transcript FASTA IDs for special characters."
            )

        # Write all alignments
        out_all = all_align_dir / f"all_alignments_{g_name}.csv"
        result_df.to_csv(out_all, index=False)

        # Write mismatch-filtered alignments
        out_mm = mismatches_dir / f"mismatches_{mismatches_n}_alignments_{g_name}.csv"
        result_df[result_df['total_mismatches'] <= mismatches_n].to_csv(out_mm, index=False)


def parse_arguments() -> tuple:
    parser = argparse.ArgumentParser(
        usage="\n",
        description="Script to design gRNA's using RNAplfold and formula designed by Ariel Bazzini",
        epilog="For any questions Contact Huzaifa Hassan: hhassan@stowers.org"
    )
    parser._optionals.title = "Optional Arguments"
    subparsers = parser.add_subparsers(title="Subcommands", dest="subcommand")

    # 'all' subcommand
    parser_all = subparsers.add_parser("all", help="Arguments to run full pipeline")
    parser_all._optionals.title = "Arguments"
    parser_all.add_argument("-f", "--fasta", type=str, required=True,
                            help="Fasta file with sequences to design gRNA")
    parser_all.add_argument("-t", "--transcript_fasta", type=str, required=True,
                            help="cDNA fasta file")
    parser_all.add_argument("-o", "--output_file", type=str, required=True,
                            help="Output CSV file")
    parser_all.add_argument("-n", "--top_N", type=int,
                            help="Top N gRNAs by score per sequence (default: 10)")
    parser_all.add_argument("-l", "--length", type=int,
                            help="gRNA length in bp (default: 22)")
    parser_all.add_argument("-m", "--mismatches", type=int,
                            help="Mismatch threshold for filtering (default: 5)")
    parser_all.add_argument("-p", "--num_processes", type=int, default=1,
                            help="Number of parallel processes (default: 1)")
    parser_all.add_argument("--fold_mode", choices=["local", "global"], default="local",
                            help="local = RNAplfold sliding window (default); "
                                 "global = RNAfold whole-sequence fold, matches the "
                                 "RNAfold web server")
    parser_all.add_argument("--winsize", type=int, default=70,
                            help="RNAplfold window size, local mode only (default: 70)")

    # 'design' subcommand
    parser_design = subparsers.add_parser("design", help="'design' specific arguments")
    parser_design._optionals.title = "Arguments"
    parser_design.add_argument("-f", "--fasta", type=str, required=True,
                               help="Fasta file with sequences to design gRNA")
    parser_design.add_argument("-o", "--output_file", type=str, required=True,
                               help="Output CSV file")
    parser_design.add_argument("-n", "--top_N", type=int,
                               help="Top N gRNAs by score per sequence (default: all windows)")
    parser_design.add_argument("-l", "--length", type=int,
                               help="gRNA length in bp (default: 22)")
    parser_design.add_argument("--fold_mode", choices=["local", "global"], default="local",
                               help="local = RNAplfold sliding window (default); "
                                    "global = RNAfold whole-sequence fold, matches the "
                                    "RNAfold web server")
    parser_design.add_argument("--winsize", type=int, default=70,
                               help="RNAplfold window size, local mode only (default: 70)")

    # 'mismatch' subcommand
    parser_mismatch = subparsers.add_parser("mismatch", help="'mismatch' specific arguments")
    parser_mismatch._optionals.title = "Arguments"
    parser_mismatch.add_argument("-g", "--grna_file", type=str, required=True,
                                 help="gRNA file (.fa, .fasta, .txt tab-delimited, or .csv)")
    parser_mismatch.add_argument("-t", "--transcript_fasta", type=str, required=True,
                                 help="cDNA fasta file")
    parser_mismatch.add_argument("-p", "--num_processes", type=int, default=1,
                                 help="Number of parallel processes (default: 1)")
    parser_mismatch.add_argument("-m", "--mismatches", type=int,
                                 help="Mismatch threshold for filtering (default: 5)")

    return parser, parser.parse_args()


if __name__ == "__main__":
    parser, args = parse_arguments()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.subcommand == "design":
        top_N  = args.top_N  or None
        length = args.length or 22
        grna_dir = mkdir(args.subcommand)
        run_RNAplfold_top(
            args.fasta, top_N, length,
            fold_mode=args.fold_mode, winsize=args.winsize
        ).to_csv(grna_dir / args.output_file, sep=",", header=True, index=False)

    elif args.subcommand == "mismatch":
        all_align_dir, mismatches_dir = mkdir(args.subcommand)
        run_mismatch(args, all_align_dir=all_align_dir, mismatches_dir=mismatches_dir)

    elif args.subcommand == "all" or args.subcommand is None:
        grna_dir, all_align_dir, mismatches_dir = mkdir(args.subcommand)
        top_N  = args.top_N  or 10
        length = args.length or 22
        grna_output = grna_dir / args.output_file
        run_RNAplfold_top(
            args.fasta, top_N, length,
            fold_mode=args.fold_mode, winsize=args.winsize
        ).to_csv(grna_output, sep=",", header=True, index=False)
        run_mismatch(
            args,
            grnafile=str(grna_output),
            all_align_dir=all_align_dir,
            mismatches_dir=mismatches_dir
        )
