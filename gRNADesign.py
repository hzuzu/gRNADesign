
#####################################
#
#author: Huzaifa Hassan
#Date: 07/17/2023
# Script to run,gRNA design, parse RNAplfold output results and FASTA sequence into 23bp sequences.
# The scores are first averaged by row and then again averaged for the scores of individual bases of the 23bp sequence. (Designed by Ariel Bazzini)

# run 'python RNAplfold_top.py -h' for details 

# 02/23
# Added parallelization, optimized the code and added more documentation
#11/21
# optimized the code for mismatch analysis.
#####################################

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Emboss.Applications import WaterCommandline
import subprocess as sb
import pandas as pd
from statistics import mean
import argparse
import re
import sys
import os
import glob
import multiprocessing as mp
import shutil
from collections import defaultdict


# Create directories
def mkdir(subcommand):
    """
    Create a directory if it does not exist
    """
    if not os.path.exists("results"):
        os.makedirs("results", exist_ok=True)
    
    if subcommand == "design":
        grna_dir = os.path.join("results", "grna")
        os.makedirs(grna_dir, exist_ok = True)
        return grna_dir
    elif subcommand == "mismatch":
        all_align_dir = os.path.join("results", "all_alignments")
        os.makedirs(all_align_dir,exist_ok = True)
        mismatches_dir = os.path.join("results", "mismatches")
        os.makedirs(mismatches_dir, exist_ok = True)
        return all_align_dir, mismatches_dir
    else:
        grna_dir = os.path.join("results", "grna")
        os.makedirs(grna_dir, exist_ok = True)
        all_align_dir = os.path.join("results", "all_alignments")
        os.makedirs(all_align_dir,exist_ok = True)
        mismatches_dir = os.path.join("results", "mismatches")
        os.makedirs(mismatches_dir, exist_ok = True)
        return grna_dir, all_align_dir, mismatches_dir



def run_RNAplfold_top(seq_fasta, top_N=None, len_grna=23):
    if sys.version_info[0] < 3:
        raise Exception("Please use Python 3.7 or higher \n To change python version on Stowers servers Run 'pyenv shell 3.7.2' ")

    cmd_RNAplFold = f"RNAplfold -W 70 -u {len_grna} < {seq_fasta}"
    sb.call(cmd_RNAplFold, shell=True)

    all_top_rows_df = [] #list to hold all genes and 22bp seqeunces in the fasta file.
    for seq_record in SeqIO.parse(seq_fasta, "fasta"): # read fasta file using biopython
        id_seq = seq_record.id # variable for seq header which will be used as name column in output file
        if '|' in id_seq:
            id_seq = id_seq.replace('|', '_')
        fasta = seq_record.seq # variable for sequence
        b = [(fasta[i:i+ len_grna],z) for i,z in enumerate(range(len(fasta)- len_grna -1))] #get seq's of 23 bp and their index
        df = pd.read_csv('_'.join(id_seq.split("|")) + "_lunp", skiprows = 1, sep = "\t", header = 0) #Read RNA plfold output file for the sequence being processed
        df['average'] = df.iloc[:, 1:len_grna +1].mean(axis =1) #Compute the mean score of all columns in each row
        averages_list = list(df['average']) # take the averages column into a list variable to maintain the order 
        d = [(averages_list[i:i+len_grna],z) for i,z in enumerate(range(len(averages_list)-len_grna -1))] # Compute the score for base i using averages of i+22 averages.
        rows_all = [] #Create an empty list to hold all the sequences.
        for i,z in  zip(b,d):
            l_str = [id_seq , str(i[0]) , str(i[1] + 1) , str(i[1] + len_grna) , str(round(mean(z[0]), 3))] #Concat all columns and reound the scores to 3 decimal places.
            rows_all.append(l_str) # append all rows to list
        #print(rows_all)
        df_rows = pd.DataFrame(rows_all) #convert the aboves list to pandas df for sorting.
        df_rows.columns = ["Gene|Transcript" ,  "seq_" + str(len_grna) + "nt" , "Seq_start" , "Seq_end" , "Average_score"] #Select only necessary columns for sorting.
        df_rows_sort = df_rows.sort_values(by=['Average_score'], ascending= False) #Sort the column score by descending order
        df_rows_sort['gRNA_name'] = df_rows_sort["Gene|Transcript"] + "_" + df_rows_sort['Seq_start'] + "_" + df_rows_sort['Seq_end'] #Create a new column with gRNA name
        if top_N: #If only the top rows are requested this will select only the top rows or else will output all.
            all_top_rows_df.append(df_rows_sort.head(top_N))
        else:
            all_top_rows_df.append(df_rows_sort)

    all_top_rows_df = pd.concat(all_top_rows_df) #Convert all seq's to pandas df

    all_top_rows_df['G_Percentage'] = round(all_top_rows_df["seq_" + str(len_grna) + "nt"].str.count('G', re.I) *100 / len_grna, 3) # get percentage of each base, this is case-insensitive so using package 're'
    all_top_rows_df['C_Percentage'] = round(all_top_rows_df["seq_" + str(len_grna) + "nt"].str.count('C', re.I) *100 / len_grna, 3)
    all_top_rows_df['A_Percentage'] = round(all_top_rows_df["seq_" + str(len_grna) + "nt"].str.count('A', re.I) *100 / len_grna, 3)
    all_top_rows_df['T_Percentage'] = round(all_top_rows_df["seq_" + str(len_grna) + "nt"].str.count('T', re.I) *100 / len_grna, 3)
    all_top_rows_df = all_top_rows_df[['gRNA_name', 'seq_' + str(len_grna) + 'nt', "Gene|Transcript", 'Seq_start', 'Seq_end',  'Average_score', 'G_Percentage', 'C_Percentage', 'A_Percentage', 'T_Percentage']] #rearrange columns
    return all_top_rows_df



# Create a water commandline object and run the command
def water_cmdline(g_fasta, tran_fasta_fh):
    """
    Run the EMBOSS water commandline for sequence alignment.
    
    Parameters:
        g_fasta (str): gRNA sequence as a string.
        tran_fasta_fh (str): File path to the FASTA file for the cDNA sequence.
    
    Returns:
        str: The stdout of the water command containing the alignment result.
    
    Raises:
        FileNotFoundError: If the water tool is not installed or not in PATH.
        CalledProcessError: If the water command fails.
    """
    # Check if the water tool is installed
    if shutil.which("water") is None:
        raise FileNotFoundError("The 'water' tool from EMBOSS is not installed or not in the PATH.")

    try:
        # Initialize the water command
        water_cmd = WaterCommandline(gapopen=10, gapextend=0.5)
        water_cmd.asequence = "asis:" + g_fasta  # Query sequence directly
        water_cmd.bsequence = tran_fasta_fh      # Subject sequence as a file
        water_cmd.outfile = "stdout"            # Redirect output to stdout

        # Run the command using subprocess
        result = sb.run(
            str(water_cmd),
            stdout=sb.PIPE,
            stderr=sb.PIPE,
            text=True,
            shell=True,
            check=True  # Raises CalledProcessError for non-zero exit codes
        )

        return result.stdout  # Return the alignment result

    except sb.CalledProcessError as e:
        print(f"Error running water command: {e.stderr}")
        raise
    except Exception as e:
        print(f"Unexpected error: {e}")
        raise


#parse the output of water

def parse_water_align(text, g_name, len_grna = 23) -> dict:
    seq1_name = g_name
    if "# 2:" in text[3]:
        seq2_name = text[3].strip().rsplit(' ', 1)[-1]
    if "# Length:" in text[8]:
        length = int(text[8].strip().rsplit(' ', 1)[-1])
    if "# Identity:" in text[9]:
        identity_perc = text[9].strip().rsplit(' ', 1)[-1].replace('(', '').replace(')', '').replace('%', '')
        identity = int(text[9].strip().rsplit(' ', 2)[-2].split('/')[0])
    if "# Similarity:" in text[10]:
        similarity_perc = text[10].strip().rsplit(' ', 1)[-1].replace('(', '').replace(')', '').replace('%', '')
        similarity = int(text[10].strip().rsplit(' ', 2)[-2].split('/')[0])
    if "# Gaps:" in text[11]:
        gRNA_gaps = text[11].strip().split('/')[0].strip().split(' ')[-1]

    seq_gaps = text[19].count('-')
    seq1_start = int([item for item in text[17].strip().split(' ') if item][1])
    seq2_start = int([item for item in text[19].strip().split(' ') if item][1])
    seq1_end = int([item for item in text[17].strip().split(' ') if item][-1])
    seq2_end = int([item for item in text[19].strip().split(' ') if item][-1])


    # Calculate additional gaps at start and end of sequence
    add_start = seq1_start - 1 # gap in start of seq
    add_end = len_grna - seq1_end # gap in end of seq
    total_add = add_start + add_end # total gap at start and end
    added_aligned_len = length + total_add # Add the total length of aligned sequence + total gaps at start and end 

    # Calculate mismatch length
    mismatch_len = added_aligned_len - identity

    return {'gRNA_name': seq1_name,
                    'seq2_name': seq2_name,
                    'aligned_length': length,
                    'identity': identity,
                    'identity_perc': identity_perc,
                    'similarity': similarity,
                    'similarity_perc': similarity_perc,
                    'gRNA_align_start': seq1_start,
                    'gRNA_align_end': seq1_end,
                    'seq_align_start': seq2_start,
                    'seq_align_end': seq2_end,
                    'gRNA_gaps': gRNA_gaps,
                    'seq_gaps': seq_gaps,
                    'total_mismatches' : mismatch_len
                    }

def water_out_parse(waterlines,g_name, len_grna = 23):
    waterlines = waterlines.split('\n')
    waterlines = waterlines[13:]
    waterlines = waterlines[:-1]
    seq1_name = "gRNA_1"
    chunks = defaultdict(list)
    count = 1
    for i in waterlines:
        if i.strip() == "#=======================================":
            chunks[count].append(i)
            count += 1
        else:
            chunks[count].append(i)
    chunk_dict = {}
    for i in range(2, len(chunks)+1, 2):
        chunk_dict[str(i) + "|" + str(i+1)] = parse_water_align(chunks[i]+ chunks[i+1], g_name,len_grna)
    return pd.DataFrame.from_dict(chunk_dict).T


def run_water_parallel(args) -> dict:
    g_name, g_seq, transSeq_fh = args
    return {g_name: water_out_parse(water_cmdline(g_seq, transSeq_fh), g_name,len_grna = len(g_seq))}

def read_delimited_file(file_path, sep = "\t",header = 0):
    """
    Read delimited file and return a dictionary
    """
    df = pd.read_csv(file_path, sep = sep, header=header)
    grna_dict =  dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
    return grna_dict

def read_fasta_file(fasta_fh):
    """
    Read fasta file and return a dictionary
    """
    fasta_dict = {}
    for seq_record in SeqIO.parse(fasta_fh, "fasta"):
        fasta_dict[seq_record.id] = str(seq_record.seq)
    return fasta_dict

def delete_files_with_endings(extensions):
    for extension in extensions:
        files = glob.glob(f"*{extension}")
        for file in files:
            try:
                os.remove(file)
            except FileNotFoundError:
                continue

def run_mismatch(args, grnafile = None, num_processes = 1):
        if grnafile:
            grnafile = grnafile
        else:
            grnafile = args.grna_file

        # Number of mismatches to report
        mismatches_n = int(args.mismatches) if args.mismatches else int(5)

        # Read gRNA fasta file
        if grnafile.endswith(".txt"):
            grna_dict = read_delimited_file(grnafile, sep = "\t",header = 0)
        elif grnafile.endswith(".csv"):
            grna_dict = read_delimited_file(grnafile, sep = ",", header = 0)
        elif grnafile.endswith((".fasta", ".fa")):
            grna_dict = read_fasta_file(grnafile)
        else:
            print("ERROR: Unsupported file format for gRNA file: " + grnafile)
            exit(1)

        # Read transcript fasta file
        transSeq_fh = args.transcript_fasta

        # run parallely for each gRNA
        if 'num_processes' not in locals() or num_processes is None:        # Ensure the number of processes is defined
            num_processes = 4  # Default to 4 processes if not defined
        else :
            num_processes = int(args.num_processes)

        try:
            # Create a pool of processes
            with mp.Pool(processes=num_processes) as pool:
                # Map the function to the list of arguments and obtain results
                results = pool.map(
                    run_water_parallel,
                    [(g_name,g_seq.replace('U', 'T'),transSeq_fh) for g_name, g_seq in grna_dict.items()]
                )

        except Exception as e:
            print(f"An error occurred during parallel processing: {e}")
            raise

        for result in results:
            for g_name, result_df in result.items():
                outfile_all = os.path.join(all_align_dir,f"all_alignments_{g_name}.csv")
                result_df.to_csv(outfile_all, index=False)

                if mismatches_n:
                    outfile_mismatches = os.path.join(mismatches_dir,f"mismatches_{mismatches_n}_alignments_{g_name}.csv")
                    df_grna_align_mismatch = result_df[result_df['total_mismatches'] <=mismatches_n]
                    df_grna_align_mismatch.to_csv(outfile_mismatches, index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(usage="\n",
                                        description='Script to design gRNA\'s using RNAplfold and formula designed by Ariel Bazzini',
                                        epilog="For any questions Contact Huzaifa Hassan: hhassan@stowers.org")

    # Create subparsers
    subparsers = parser.add_subparsers(title="Subcommands", dest="subcommand")

    # Create a parser for the 'all' subcommand
    parser_all = subparsers.add_parser("all", help="Arguments to run full pipeline")
    parser_all.add_argument("-f", "--fasta", type = str,  help='Fasta file with sequences to design gRNA',required=True)
    parser_all.add_argument("-t", "--transcript_fasta", help='cDNA fasta file', type=str, required=True)
    parser_all.add_argument("-o", "--output_file",  type = str, help='Output CSV file', required=True)
    parser_all.add_argument("-n", "--top_N", type = int, help='The number of sequences to write to file (based on highest score) (Optional). Default is top 10')
    parser_all.add_argument("-l", "--length", type= int, help='Length of gRNAs to design (Optional). Default is 23')
    parser_all.add_argument("-m", "--mismatches", type=int, help='The number of mismatches to report in a separate file. Default is 5')
    parser_all.add_argument("-p", "--num_processes", type=int,default=1, help='The number of processes/threads to use. Default is 1')

    parser_all._optionals.title = "Arguments"

    # Create a parser for the 'design' subcommand
    parser_design = subparsers.add_parser("design", help="'design' specific arguments")
    parser_design.add_argument("-f", "--fasta", type = str, help='Fasta file with sequences to design gRNA', required=True)
    parser_design.add_argument("-o", "--output_file",  type = str, help='Output CSV file', required=True)
    parser_design.add_argument("-n", "--top_N", type = int, help='The number of sequences to write to file (based on highest score) (Optional). Default is top 10')
    parser_design.add_argument("-l", "--length", type= int, help='Length of gRNAs to design (Optional). Default is 22')
    parser_design._optionals.title = "Arguments"


    # Create a parser for the 'mismatch' subcommand
    parser_mismatch = subparsers.add_parser("mismatch", help="'mismatch' specific arguments")
    parser_mismatch.add_argument("-g", "--grna_file",type=str, help='Fasta file with sequences (recommended fasta file)', required=True)
    parser_mismatch.add_argument("-t", "--transcript_fasta",type=str, help='cDNA fasta file', required=True)
    parser_mismatch.add_argument("-p", "--num_processes", type=int,default=1, help='The number of processes/threads to use. Default is 1')
    parser_mismatch.add_argument("-m", "--mismatches", type=int, help='The number of mismatches to report in a separate file. Default is 5')
    parser_mismatch._optionals.title = "Arguments"
    parser._optionals.title = "Optional Arguments"

    # Parse the command-line arguments
    return parser, parser.parse_args()

if __name__ == "__main__":
    parser, args = parse_arguments()

    if len(sys.argv) == 1: 
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.subcommand == "design":
        top_N = int(args.top_N) if args.top_N else None
        length = int(args.length) if args.length else 22
        grna_dir = mkdir(args.subcommand)

        run_RNAplfold_top(args.fasta, top_N, length).to_csv(os.path.join(grna_dir, args.output_file), sep=",", header=True, index=False)
        delete_files_with_endings(["_lunp", "_dp.ps"])

    elif args.subcommand == "mismatch":
        all_align_dir, mismatches_dir = mkdir(args.subcommand)

        run_mismatch(args)

    elif (args.subcommand == "all") or (args.subcommand == None):
        grna_dir, all_align_dir, mismatches_dir = mkdir(args.subcommand)
        top_N = int(args.top_N) if args.top_N else 10
        length = int(args.length) if args.length else 22

        run_RNAplfold_top(args.fasta, top_N, length).to_csv(os.path.join(grna_dir, args.output_file), sep=",", header=True, index=False)
        delete_files_with_endings(["_lunp", "_dp.ps"])
        run_mismatch(args, os.path.join(grna_dir, args.output_file))