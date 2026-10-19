#!/usr/bin/env python3

"""
Author: hhassan
Date: 03/29/2024

Select gRNAs that are a minimum distance apart, based on the start/end coordinates
written by the design module, keeping the highest-scoring ones.

Input is the gRNA CSV from `gRNADesign.py design` (run it without -n so that all
windows are available to choose from).
"""
import argparse
import logging
import sys

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

__version__ = "1.2.0"

# The design module writes 'Gene|Transcript'. Older inputs used
# 'Gene|Longest_transcript'; accept either so the column does not have to be renamed
# by hand before running this script.
GROUP_COLUMNS = ("Gene|Transcript", "Gene|Longest_transcript")

REQUIRED_COLUMNS = ("gRNA_name", "Seq_start", "Seq_end", "Average_score")


def resolve_group_column(df: pd.DataFrame) -> str:
    """Return whichever accepted grouping column is present, else raise."""
    for col in GROUP_COLUMNS:
        if col in df.columns:
            return col
    raise ValueError(
        f"No grouping column found. Expected one of {list(GROUP_COLUMNS)}; "
        f"the file has {list(df.columns)}."
    )


def select_gRNA(df: pd.DataFrame, dist_grna: int = 10, top_n: int = 10) -> pd.DataFrame:
    """
    Select up to top_n gRNAs from one group, highest 'Average_score' first, keeping only
    those whose start and end coordinates are at least dist_grna from every coordinate
    already selected.
    """
    df_sorted = df.sort_values(by="Average_score", ascending=False)

    coords: list = []        # start and end positions of everything selected so far
    selected_ids: list = []

    for _, row in df_sorted.iterrows():
        if not selected_ids:
            selected_ids.append(row['gRNA_name'])
            coords.extend([row['Seq_start'], row['Seq_end']])
            continue

        # Keep the gRNA only if both of its ends clear every selected coordinate
        far_enough = all(
            abs(pos - c) >= dist_grna
            for pos in (row['Seq_start'], row['Seq_end'])
            for c in coords
        )
        if far_enough:
            selected_ids.append(row['gRNA_name'])
            coords.extend([row['Seq_start'], row['Seq_end']])
            if len(selected_ids) >= top_n:
                break

    return df_sorted[df_sorted['gRNA_name'].isin(selected_ids)]


def select_gRNA_per_gene(df_fh: str, dist_grna: int = 10, top_n: int = 10) -> pd.DataFrame:
    """
    Apply select_gRNA to every gene/transcript group in the input CSV.

    Parameters:
        df_fh: Path to the gRNA CSV written by the design module.
        dist_grna: Minimum distance between selected gRNAs.
        top_n: Maximum number of gRNAs to keep per group.

    Returns:
        DataFrame with the selected rows from every group.
    """
    try:
        df = pd.read_csv(df_fh, sep=',', header=0)
    except (OSError, ValueError) as e:
        raise RuntimeError(f"Could not read gRNA CSV '{df_fh}': {e}") from e

    group_col = resolve_group_column(df)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Input '{df_fh}' is missing required column(s): {missing}")

    # Groups are iterated explicitly rather than via groupby().apply(). In pandas 2.2+
    # that call warns because it operates on the grouping column, and pandas 3 drops that
    # column from the result - which would silently remove the gene name from the output.
    selected = [
        select_gRNA(group, dist_grna=dist_grna, top_n=top_n)
        for _, group in df.groupby(group_col, sort=False)
    ]
    result = pd.concat(selected, ignore_index=True) if selected else df.head(0)
    logger.info(
        f"Selected {len(result)} gRNA(s) from {df.groupby(group_col).ngroups} group(s) "
        f"in '{group_col}' (min distance {dist_grna} bp, up to {top_n} per group)."
    )
    return result


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Select the top gRNAs per gene by 'Average_score' while keeping them "
                    "at least a minimum distance apart.",
        usage=f'{sys.argv[0]} -f <df_fh> [-d <dist_grna>] [-n <top_n>]',
        epilog='Contact Huzaifa Hassan at hhassan@stowers.org for help or questions.'
    )
    parser.add_argument('-f', '--df_fh', type=str, required=True,
                        help='File path of the input gRNA CSV file.')
    parser.add_argument('-d', '--dist_grna', type=int, default=30,
                        help='Minimum distance between selected gRNAs. Default is 30.')
    parser.add_argument('-n', '--top_n', type=int, default=10,
                        help='Number of gRNAs to select within each group. Default is 10.')
    parser.add_argument('-v', '--version', action='version',
                        version=f'select_gRNA_dist {__version__}')
    return parser.parse_args()


def main() -> None:
    """Select gRNAs by distance and score, and write them next to the input file."""
    args = parse_arguments()

    result_df = select_gRNA_per_gene(args.df_fh, args.dist_grna, args.top_n)
    out_fh = args.df_fh.replace('.csv', f'.{args.dist_grna}bp.top_{args.top_n}.csv')
    result_df.to_csv(out_fh, index=False)
    logger.info(f"Wrote {out_fh}")


if __name__ == '__main__':
    main()
