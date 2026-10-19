

# gRNADesign: A pipeline to design gRNA and calculate RNA fold scores and mismatches to cDNA's

**v1.1.0** · Python 3.12 · It is specifically designed for Bazzini Lab at Stowers Institute for Medical Research.

## Information

There are 3 submodules in this pipeline
* design
* mismatch
* all

### design module

The design module designs gRNAs from an input fasta file. This modules runs RNAplfold and computes a score using method designed by Ariel Bazzini

1. Run RNA fold for each gene (Complete sequence). RNA fold takes a fasta sequence as input along with the length of gRNAs.
2. Take sliding window of gRNAs for a specific length from gene fasta. (The length of gRNA will be equal length from step 1. e.g., 20, 22, 23)
3. Compute the mean score of all columns (The number of columns will be equal to length of the gRNA) in each row. This gives one value for each bp of a gene.
4. For each sequence from step 2, take the scores from step 3 (The number of scores should be equal to the length of the gRNA sequence). Compute a mean of these scores to give one single score for each gRNA

#### Folding mode: `--fold_mode`

The design module can fold the input two ways. This changes which sites score highest, so pick deliberately.

| `--fold_mode` | Tool | Base pairs allowed | Use when |
|---|---|---|---|
| `local` (default) | `RNAplfold -W <winsize>` | only within the window (default 70 nt) | you want the original pipeline behaviour |
| `global` | `RNAfold -p` over the whole sequence | any span, including long-range | you want to reproduce the ViennaRNA RNAfold **web server** |

**`global` reproduces the web server exactly.** The RNAfold web server states: *"For unpaired regions the color denotes the probability of being unpaired."* `global` mode scores each window by that same per-base unpaired probability. Verified on cth1 (1553 nt): all 1532 windows identical to the web-server quantity, max difference 0.000000, Spearman 1.000.

**The two modes disagree substantially on long mRNAs.** On cth1 they share **0 of their top 10** sites (Spearman 0.31), because a 70 nt window cannot see the long-range pairs a global fold forms. Neither is automatically correct: global folding assumes long-range pairs actually form, while local folding reflects the fact that long mRNAs fold co-transcriptionally and are covered by ribosomes and RBPs. Only knockdown data can settle which predicts activity better.

`--winsize` sets the RNAplfold window in `local` mode (default 70, which is RNAplfold's own default). It is ignored in `global` mode. **`--winsize` must be at least `--length`**; the run is refused otherwise. RNAplfold accepts `-u` greater than `-W` without complaining but fills the columns past the window with placeholders (`NA` and `-0`), which would be averaged into the score and silently corrupt it. This matters if you raise `--length` without raising `--winsize`.

Note: `RNAfold -p` scales roughly as O(n^3) in time and O(n^2) in memory, so `global` mode is slow on very long transcripts; a warning is printed above 5000 nt. (`RNAplfold` with a very large `-W` is not a substitute — it overflows and silently produces no output on longer sequences.)

**Input file**

Fasta File: The input file is a FASTA file containing sequences for which gRNAs should be designed.

**Sequence IDs must be unique after `|` is replaced by `_`.** RNAplfold and RNAfold name their output files after the sequence ID, so two records whose IDs collide (`dup`/`dup`, or `a|b`/`a_b`) overwrite each other's fold and every colliding record would be scored from a single structure. The run is refused with the offending pair named.

### mismatch module

The mismatch module is used to compute the number of mismatches for gRNAs against the transcriptome.
The module runs the water pairwise aligner (Smith-Waterman local alignment) from EMBOSS.

How it works:
1. Each gRNA (U is converted to T) is aligned with one `water` run against **every** transcript in the transcript FASTA (`-p` runs several gRNAs in parallel).
2. `water` reports only the **single best local alignment per transcript**, so a second site in the same transcript is not reported.
3. Only the **same strand** as the transcript FASTA is searched (no reverse complement).
4. The intended target transcript is not excluded; it appears in the output with 0 mismatches.
5. Transcript IDs containing `|` are aligned as `_` internally (EMBOSS would otherwise rename them `EMBOSS_001`, ...) and restored in the output. Zero-length records are skipped by `water` and reported as a warning.

**Input file**
The input files are:
* gRNA File: The gRNA file can be provided in one of the following formats:
  * FASTA format: The file should end with `.fa` or `.fasta`, with the FASTA headers representing the gRNA names.
  * Tab-delimited format: The file should end with `.txt`, with the first column containing the gRNA names and the second column containing the gRNA sequences. The file must include headers.
  * Comma-delimited format: The file should end with `.csv`, same layout as `.txt` (name, sequence; with headers). The `gRNA.csv` written by the design module can be used directly.
* Transcript/cDNA FASTA File: A FASTA file containing transcript or cDNA sequences.

### all module

The 'all' module runs both the 'design' and 'mismatch' module

**Input file**
The input files are:
* Fasta File: The input file is a FASTA file containing sequences for which gRNAs should be designed.
* Transcript/cDNA FASTA File: A FASTA file containing transcript or cDNA sequences.


## Getting started

`gRNADesign.yml` installs everything the pipeline needs, including the external tools:
Python 3.12, Biopython, pandas, ViennaRNA (`RNAplfold`, `RNAfold`) and EMBOSS (`water`).

**Create and activate the environment**

```bash
mamba env create -f gRNADesign.yml    # or: conda env create -f gRNADesign.yml
conda activate gRNADesign
```

Check it worked — all three must resolve:

```bash
python -c "import Bio, pandas; print(Bio.__version__, pandas.__version__)"
which RNAplfold RNAfold water
```

EMBOSS is now part of the environment file. If you would rather use the system build on
Stowers HPC, it is available as a module instead:

```bash
ml emboss/6.6.0
```


## Usage

Check the different modules in the pipeline

``` bash
python gRNADesign.py --help

usage:

Script to design gRNA's using RNAplfold and formula designed by Ariel Bazzini

Optional Arguments:
  -h, --help            show this help message and exit

Subcommands:
  {all,design,mismatch}
    all                 Arguments to run full pipeline
    design              'design' specific arguments
    mismatch            'mismatch' specific arguments

For any questions Contact Huzaifa Hassan: hhassan@stowers.org

```

To check the arguments for each module, run

```bash

python gRNADesign.py <submodule name> --help


python gRNADesign.py design --help

python gRNADesign.py mismatch --help

python gRNADesign.py all --help

usage:  all [-h] -f FASTA -t TRANSCRIPT_FASTA -o OUTPUT_FILE [-n TOP_N]
            [-l LENGTH] [-m MISMATCHES] [-p NUM_PROCESSES]

Arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Fasta file with sequences to design gRNA
  -t TRANSCRIPT_FASTA, --transcript_fasta TRANSCRIPT_FASTA
                        cDNA fasta file
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output CSV file; must end in '.csv' (it is re-read by
                        the mismatch step)
  -n TOP_N, --top_N TOP_N
                        Top N gRNAs by score per sequence (default: 10)
  -l LENGTH, --length LENGTH
                        gRNA length in bp (default: 22)
  -m MISMATCHES, --mismatches MISMATCHES
                        Mismatch threshold for filtering (default: 5)
  -p NUM_PROCESSES, --num_processes NUM_PROCESSES
                        Number of parallel processes (default: 1)
  --fold_mode {local,global}
                        local = RNAplfold sliding window (default); global =
                        RNAfold whole-sequence fold, matches the RNAfold web
                        server
  --winsize WINSIZE     RNAplfold window size, local mode only (default: 70)


python gRNADesign.py design --help

usage:  design [-h] -f FASTA -o OUTPUT_FILE [-n TOP_N] [-l LENGTH]

Arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Fasta file with sequences to design gRNA
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output CSV file
  -n TOP_N, --top_N TOP_N
                        Top N gRNAs by score per sequence (default: all windows)
  -l LENGTH, --length LENGTH
                        gRNA length in bp (default: 22)


python gRNADesign.py mismatch --help

usage:  mismatch [-h] -g GRNA_FILE -t TRANSCRIPT_FASTA [-p NUM_PROCESSES]
                 [-m MISMATCHES]

Arguments:
  -h, --help            show this help message and exit
  -g GRNA_FILE, --grna_file GRNA_FILE
                        gRNA file (.fa, .fasta, .txt tab-delimited, or .csv)
  -t TRANSCRIPT_FASTA, --transcript_fasta TRANSCRIPT_FASTA
                        cDNA fasta file
  -p NUM_PROCESSES, --num_processes NUM_PROCESSES
                        Number of parallel processes (default: 1)
  -m MISMATCHES, --mismatches MISMATCHES
                        Mismatch threshold for filtering (default: 5)

```

Note: `-m 0` keeps only perfect matches. The mismatch output files are named `mismatches_<m>_alignments_<gRNA>.csv`.

These are examples for each module

```bash

python gRNADesign.py all -f gene.fasta -t cDNA_transcript.fasta -o gRNA.csv -n 10 -l 22 -m 5

python gRNADesign.py design -f gene.fasta -o gRNA.csv -n 10 -l 22 

python gRNADesign.py mismatch -g  gRNA.csv -t cDNA_transcript.fasta -m 5 -p 4

# reproduce the RNAfold web server readout instead of the local sliding window
python gRNADesign.py design -f gene.fasta -o gRNA.csv -n 10 -l 22 --fold_mode global

```


## Results

The results are written to a folder 'results'

### Results Folders
results/grna : The gRNAs designed are in this folder <br/>
results/all_alignments : The alignments of gRNA to the transcripts are in this folder <br/>
results/mismatches : The alignments of gRNA to the transcripts having mismatches less than or equal to the specified mismatches. <br/>

### Files
The results/grna will have one file with the name used while running the pipeline (eg gRNA.csv)

The columns of gRNA.csv are described below:

* gRNA_name : The name here will be the sequence id of the input sequences followed by the start and the end coordinate where this gRNA is from on the sequence.
* seq_22nt : The target-site sequence, as it appears in the input FASTA (same strand). Default length is 22, unless user defined.
* seq_22nt_revcomp : The reverse complement of the above. Follows the input alphabet (DNA stays T, RNA stays U). Use whichever strand your system requires.
* Gene|Transcript : The name here will be the sequence id of the input sequences
* Seq_start : The start coordinate where this gRNA is from on the sequence
* Seq_end : The end coordinate where this gRNA is from on the sequence
* Average_score : Average is the score defined in the design module
* G_Percentage : Percentage of Base 'G' in the gRNA
* C_Percentage : Percentage of Base 'C' in the gRNA
* A_Percentage : Percentage of Base 'A' in the gRNA
* T_Percentage : Percentage of Base 'T' in the gRNA

Note on base composition: these are computed on `seq_22nt` (the target site). For the reverse complement they swap, i.e. G%<->C% and A%<->T%.

Note on column order: the mismatch module reads a gRNA file by column position (column 0 = name, column 1 = sequence). `seq_22nt` is therefore kept at column 1 so that `gRNA.csv` can be passed straight to the mismatch module, which aligns the target-site strand. Do not reorder the first two columns.


The results/all_alignments and results/mismatches will have one file for each gRNA. results/all_alignments will have all the alignments of gRNA mapping to each gene/transcript from the TRANSCRIPT_FASTA file and  results/mismatches will have the same alignments filtered for gene/transcript with total mismatch of specified mismatches. 

The columns of the alignment files (eg all_alignments_gRNA_1.csv) are described below. There is one row per transcript (zero-length transcripts are skipped):

* gRNA_name : Name of the gRNA
* seq2_name : The ID of the transcript from the input TRANSCRIPT_FASTA file
* aligned_length : Length of the local alignment in columns (including gap columns)
* identity : Number of identical positions in the alignment
* identity_perc : identity / aligned_length, as a percentage (number only, no `%` sign)
* similarity : Number of similar positions as reported by water
* similarity_perc : similarity / aligned_length, as a percentage
* gRNA_align_start : Alignment start on the gRNA (1-based)
* gRNA_align_end : Alignment end on the gRNA
* seq_align_start : Alignment start on the transcript (1-based)
* seq_align_end : Alignment end on the transcript
* gRNA_gaps : Number of gap positions in the gRNA row of the alignment
* seq_gaps : Number of gap positions in the transcript row of the alignment
* total_mismatches : aligned_length + (gRNA_align_start - 1) + (gRNA length - gRNA_align_end) - identity. This counts mismatches, gaps, and the gRNA bases left outside the local alignment at either end.

Note: Please feel free to compute your own mismatch score from files in results/all_alignments if needed.



**This page is still being update**



**Helper Scripts**
1. Filter gRNAs based on the distance between the gRNAs

This script is used to select gRNAs based on the distance between the gRNAs. The script takes the gRNA.csv file as input and the distance between the gRNAs. The script will select the gRNAs based on the distance and write the selected gRNAs to a new file.

Note: the script accepts either `Gene|Transcript` (what the design module writes) or the older `Gene|Longest_transcript`, so no column renaming is needed. Run the design module without `-n` so that all windows are available for the distance selection.

Usage:

```bash
python select_gRNA_dist.py -h
usage: select_gRNA_dist.py -f <df_fh> [-d <dist_grna>] [-n <top_n>]

Select the top gRNAs per gene by 'Average_score' while keeping them at least a minimum
distance apart.

options:
  -h, --help            show this help message and exit
  -f DF_FH, --df_fh DF_FH
                        File path of the input gRNA CSV file.
  -d DIST_GRNA, --dist_grna DIST_GRNA
                        Minimum distance between selected gRNAs. Default is 30.
  -n TOP_N, --top_n TOP_N
                        Number of gRNAs to select within each group. Default is 10.
  -v, --version         show program's version number and exit

```

