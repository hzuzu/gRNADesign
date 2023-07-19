
#####################################
#
#author: Huzaifa Hassan
#Date: 07/17/2023
# Script to run,gRNA design, parse RNAplfold output results and FASTA sequence into 22bp sequences.
# The scores are first averaged by row and then again averaged for the scores of individual bases of the 22bp sequence. (Designed by Ariel Bazzini)

# run 'python RNAplfold_top.py -h' for details 
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



def run_RNAplfold_top(seq_fasta, top_N=None, len_grna=22):
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
        df_rows.columns = ["Gene|Longest_transcript" ,  "seq_" + str(len_grna) + "nt" , "Seq_start" , "Seq_end" , "Average_score"] #Select only necessary columns for sorting.
        df_rows_sort = df_rows.sort_values(by=['Average_score'], ascending= False) #Sort the column score by descending order
        df_rows_sort['gRNA_name'] = df_rows_sort['Gene|Longest_transcript'] + "_" + df_rows_sort['Seq_start'] + "_" + df_rows_sort['Seq_end'] #Create a new column with gRNA name
        if top_N: #If only the top rows are requested this will select only the top rows or else will output all.
            all_top_rows_df.append(df_rows_sort.head(top_N))
        else:
            all_top_rows_df.append(df_rows_sort)

    all_top_rows_df = pd.concat(all_top_rows_df) #Convert all seq's to pandas df

    all_top_rows_df['G_Percentage'] = round(all_top_rows_df["seq_" + str(len_grna) + "nt"].str.count('G', re.I) *100 / len_grna, 3) # get percentage of each base, this is case-insensitive so using package 're'
    all_top_rows_df['C_Percentage'] = round(all_top_rows_df["seq_" + str(len_grna) + "nt"].str.count('C', re.I) *100 / len_grna, 3)
    all_top_rows_df['A_Percentage'] = round(all_top_rows_df["seq_" + str(len_grna) + "nt"].str.count('A', re.I) *100 / len_grna, 3)
    all_top_rows_df['T_Percentage'] = round(all_top_rows_df["seq_" + str(len_grna) + "nt"].str.count('T', re.I) *100 / len_grna, 3)
    all_top_rows_df = all_top_rows_df[['gRNA_name', 'seq_' + str(len_grna) + 'nt', 'Gene|Longest_transcript', 'Seq_start', 'Seq_end',  'Average_score', 'G_Percentage', 'C_Percentage', 'A_Percentage', 'T_Percentage']] #rearrange columns
    return all_top_rows_df



# Create a water commandline object
def water_cmdline(g_fasta, tran_fasta):
    """
    run water commandline
    """
    water_cmd = WaterCommandline(gapopen=10, gapextend=0.5)
    water_cmd.asequence = "asis:" + g_fasta
    water_cmd.bsequence = "asis:" + tran_fasta
    water_cmd.outfile = "stdout"  # Redirect the output to stdout

    run_water = sb.run(str(water_cmd), capture_output=True, text=True, shell=True)
    stdout_output = run_water.stdout
    return stdout_output


#parse the output of water
def parse_water_output(stdout_output, len_grna):
    """
    parse the output of water and return the alignment length, identity, mismatch, and gap
    """
    if not isinstance(len_grna, int):
        len_grna = int(len_grna)

    lines = stdout_output.split('\n')
    line_1 = lines[23].strip().split('/')
    identity = int(line_1[0].strip().split(' ')[-1]) # aligned identity
    total_aligned_len = int(line_1[1].strip().split(' ')[0]) #total length of aligned sequcne including gaps and mismatches

    line_2 = lines[31].strip().split(' ')
    line_2 = list(filter(None, line_2))
    align_start = int(line_2[1]) #sequence start of query
    align_end = int(line_2[-1]) #sequence end of query

    add_start = align_start - 1 #gap in start of seq
    add_end = len_grna - align_end #gap in end of seq
    total_add = add_start + add_end #total gap at start and end
    added_aligned_len = total_aligned_len + total_add #Add the total length of aligned sequcne + total gaps at start and end 

    mismatch_len = added_aligned_len - identity
    return [align_start, align_end,identity, total_aligned_len, add_start, add_end, added_aligned_len, mismatch_len]
    


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

def run_mismatch(args, grnafile = None):
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
        trans_dict = read_fasta_file(args.transcript_fasta)
        

        for g_name, g_seq in grna_dict.items():
            #outfile_all = f"all_alignments_{g_name}.csv"
            outfile_all = os.path.join(all_align_dir,f"all_alignments_{g_name}.csv")
            y = open(outfile_all, 'w')
            y.write("Gene/transcript" + "," + "align_start" + "," + "align_end" + "," + "identity" + "," + "query_aligned_len" + "," + "missing_start_bp" + "," + "missing_end_bp" + "," + "total_missing_bp" + "," + "total_aligned_length" + ',' + "mismatches" + '\n')

            if mismatches_n:
                outfile_mismatches = os.path.join(mismatches_dir,f"mismatches_{mismatches_n}_alignments_{g_name}.csv")
                z = open(outfile_mismatches, 'w')
                z.write(f"Gene/transcript" + "," + "align_start" + "," + "align_end" + "," + "identity" + "," + "query_aligned_len" + "," + "missing_start_bp" + "," + "missing_end_bp" + "," + "total_missing_bp" + "," + "total_aligned_length" + ',' + "mismatches" + '\n')

            for seq_id, trans_seq in trans_dict.items():
                #print(seq_id, trans_seq)
                water_out= water_cmdline(g_seq, trans_seq)
                #print(water_out)
                rnafold_out= [str(seq_id)] + parse_water_output(water_out, str(len(g_seq)))
                y.write(",".join(str(i) for i in rnafold_out) + '\n')
                if int(rnafold_out[-1]) <= mismatches_n:
                    z.write(",".join(str(i) for i in rnafold_out) + '\n')
            y.close()
            if mismatches_n:
                z.close()

def parse_arguments():
    parser = argparse.ArgumentParser(usage="\n",
                                        description='Script to design gRNA\'s using RNAplfold and formula designed by Ariel Bazzini',
                                        epilog="For any questions Contact Huzaifa Hassan: hhassan@stowers.org")

    # Create subparsers
    subparsers = parser.add_subparsers(title="Subcommands", dest="subcommand")

    # Create a parser for the 'all' subcommand
    parser_all = subparsers.add_parser("all", help="Arguments to run full piepline")
    parser_all.add_argument("-f", "--fasta", type = str,  help='Fasta file with sequences to design gRNA',required=True)
    parser_all.add_argument("-t", "--transcript_fasta", help='cDNA fasta file', type=str, required=True)
    parser_all.add_argument("-o", "--output_file",  type = str, help='Output CSV file', required=True)
    parser_all.add_argument("-n", "--top_N", type = int, help='The number of sequences to write to file (based on highest score) (Optional). Default is top 10')
    parser_all.add_argument("-l", "--length", type= int, help='Length of gRNAs to design (Optional). Default is 22')
    parser_all.add_argument("-m", "--mismatches", type=int, help='The number of mismatches to report in a seperate file. Default is 5')

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

    parser_mismatch.add_argument("-m", "--mismatches", type=int, help='The number of mismatches to report in a seperate file. Default is 5')
    parser_mismatch._optionals.title = "Arguments"
    parser._optionals.title = "Optional Arguments"

    # Parse the command-line arguments
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

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


    # else:
    #     print("Invalid subcommand. Use 'foo' or 'bar'.")
    # #mode = args.sub
    # if mode == "design":
    #     top_N = int(args.top_N) if args.top_N else None
    #     length = int(args.length) if args.length else 22
    #     grna_dir = mkdir()

    #     run_RNAplfold_top(args.fasta, top_N, length).to_csv(os.path.join(grna_dir, args.output_file), sep=",", header=True, index=False)

