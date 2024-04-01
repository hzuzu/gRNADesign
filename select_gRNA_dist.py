#!/usr/bin/env python

"""
Author: hhassan
Date: 03/29/2024
python script to select gRNA which are a certain distance apart based on the seq_start of the g_RNA and select the top 10 gRNAs based on the score

"""
import pandas as pd
import sys


import pandas as pd

def select_gRNA(df, dist_grna=10, top_n=10):
    """
    Selects the top n rows based on 'Average_score' within each group of
    'Gene|Longest_transcript' column, while ensuring that the gRNAs are at least dist_grna apart.
    """
    df_sorted = df.sort_values(by="Average_score", ascending=False)
    
    start = [] # store the start positions of selected IDs
    end = [] # store the end positions of selected IDs
    selected_ids = [] # store the selected IDs

    for index, row in df_sorted.iterrows():
        # store 1st id 
        if not selected_ids:
            selected_ids.append(row['gRNA_name'])
            start.append(row['Seq_start'])
            end.append(row['Seq_end'])
        else:
            # check if current id is at least 20 bp away from any previous start and any previous end in the list
            if all(abs(row['Seq_start'] - i) >= dist_grna for i in start) and all(abs(row['Seq_start'] - i) >= dist_grna for i in end) and all(abs(row['Seq_end'] - i) >= dist_grna for i in start) and all(abs(row['Seq_end'] - i) >= dist_grna for i in end):
                selected_ids.append(row['gRNA_name'])
                start.append(row['Seq_start'])
                end.append(row['Seq_end'])
                if len(selected_ids) >= top_n:  
                    break
    return df_sorted[df_sorted['gRNA_name'].isin(selected_ids)]

def select_gRNA_per_gene(df_fh, dist_grna=10, top_n=10):
    """
    Selects the top n rows based on 'Average_score' within each group of
    'Gene|Longest_transcript' column, while ensuring that the 'Seq_start'
    values are at least dist_grna apart.
    
    Parameters:
        df (DataFrame): Input DataFrame containing columns 'Gene|Longest_transcript', 'gRNA_name', 'Seq_start', 'Seq_end', and 'Average_score'.
        dist_grna (int): Minimum distance between gRNAs.
        top_n (int): Number of top rows to select.
    
    Returns:
        DataFrame: Merged DataFrame containing the selected rows for each group.
    """
    # read the input DataFrame
    df = pd.read_csv(df_fh, sep=',', header=0)
    # Group the DataFrame by 'Gene|Longest_transcript' and apply the select_gRNA function to each group
    grouped_df = df.groupby('Gene|Longest_transcript').apply(select_gRNA, dist_grna=dist_grna, top_n=top_n)

    # Reset the index of the resulting DataFrame
    grouped_df.reset_index(drop=True, inplace=True)
    
    return grouped_df


def args():
    """
    Parses the command line arguments.

    Returns:
        argparse.Namespace: An object containing the parsed command line arguments.
    """
    import argparse

    parser = argparse.ArgumentParser(description='Selects the top n rows based on \'Average_score\' within each group of \'Gene|Longest_transcript\' column, while ensuring that the \'Seq_start\' values are at least dist_grna apart.',
                                     usage=f'{sys.argv[0]} -f <df_fh> -d <dist_grna> -n <top_n>',
                                     epilog='Contact Huzaifa Hassan at hhassan@stowers.org for help or questions.')
    parser.add_argument('-f', '--df_fh', type=str, help='File path of the input CSV file.')
    parser.add_argument('-d', '--dist_grna', type=int, default=30, help='Minimum distance between \'Seq_start\' values of selected rows. Default is 30.')
    parser.add_argument('-n', '--top_n', type=int, default=10, help='Number of top rows to select within each group. Default is 10.')
    args = parser.parse_args()
    return args

def main():
    """
    Main function that selects gRNA based on distance and top n criteria.

    Usage: python select_gRNA_dist.py -f <df_fh> -d <dist_grna> -n <top_n>

    Args:
        None

    Returns:
        None
    """
    args_ = args()
    # check if  df_fh is not provided, if not print usage
    if args_.df_fh is None:
        print("Usage: python select_gRNA_dist.py -f <df_fh> -d <dist_grna> -n <top_n>")
        sys.exit(1)
    
    df_fh = args_.df_fh
    dist_grna = args_.dist_grna
    top_n = args_.top_n

    result_df = select_gRNA_per_gene(df_fh, dist_grna, top_n)
    result_df.to_csv(df_fh.replace('.csv', '.' +  str(dist_grna) + 'bp.top_' + str(top_n) + '.csv'), index=False)


if __name__ == '__main__':
    main()
