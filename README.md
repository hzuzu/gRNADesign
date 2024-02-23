

# gRNADesign: A pipeline to design gRNA and calculate RNA fold scores, mismatches to cDNA's
It is specifically designed for Bazzini Lab.

## Information

There are 3 submodules in this pipeline
* design
* msimatch
* all

### design module

The design module designs gRNAs from an input fasta file. This modules runs RNAplfold and computes a score using method designed by Ariel Bazzini

1. Run RNA fold for each gene (Complete sequence). RNA fold takes a fasta sequence as input along with the length of gRNAs.
2. Take sliding window of gRNAs for a specific length from gene fasta. (The length of gRNA will be equal length from step 1. e.g., 20, 22, 23)
3. Compute the mean score of all columns (The number of columns will be equal to length of the gRNA) in each row. This gives one value for each bp of a gene.
4. For each sequence from step 2, take the scores from step 3 (The number of scores should be equal to the length of the gRNA sequence). Compute a mean of these scores to give one single score for each gRNA.

### mismatch module

The mismatch module is used to compute the number of mismatches for gRNAs against the transcriptome.
The module runs water pairwise aligner from Emboss.

### all module

The 'all' module runs both the 'design' and 'mismatch' module


## Getting started
This pipeline runs in a conda environment and all the required packages are installed in the conda environment. <br/>
However additional tools mentioned below are required for running the pipeline. Make sure you have them installed and in path.

**Additional tools**


The mismatch mdule uses water aligner from [EMBOSS](https://ssbio.readthedocs.io/en/latest/instructions/emboss.html).
The installation instructions are mentioned in the links above. 

To run it on Stowers HPC, Emboss module is already installed and can be loaded using

```bash
ml emboss/6.6.0
```


**Create a Conda environment using the yaml file** 

```bash

conda env create -f gRNADesign.yml

```

**Activate the gRNADesign environment** 

```bash

conda activate gRNADesign

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
    all                 Arguments to run full piepline
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
                        Output CSV file
  -n TOP_N, --top_N TOP_N
                        The number of sequences to write to file (based on
                        highest score) (Optional). Default is top 10
  -l LENGTH, --length LENGTH
                        Length of gRNAs to design (Optional). Default is 23
  -m MISMATCHES, --mismatches MISMATCHES
                        The number of mismatches to report in a seperate file.
                        Default is 5
  -p NUM_PROCESSES, --num_processes NUM_PROCESSES
                        The number of processes/threads to use. Default is 1

```

These are examples for each module

```bash

python gRNADesign.py all -f gene.fasta -t cDNA_transcript.fasta -o gRNA.csv -n 10 -l 22 -m 5

python gRNADesign.py design -f gene.fasta -o gRNA.csv -n 10 -l 22 

python gRNADesign.py mismatch -g  gRNA.csv -t cDNA_transcript.fasta -m 5 -p 4

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

* gRNA_name : The name here will be the sequnce id of the input sequences followed by the start and the end cordinate where this gRNA is from on the seqenuce.
* seq_23nt : The lenght of the gRNA. Default is 23, unless user defined.
* Gene|Transcript : The name here will be the sequnce id of the input sequences.
* Seq_start : The start coordinate where this gRNA is from on the seqenuce.
* Seq_end : The end coordinate where this gRNA is from on the seqenuce.
* Average_score : Average is the score defined in the design module
* G_Percentage : Percenate of Base 'G' in the gRNA
* C_Percentage : Percenate of Base 'C' in the gRNA
* A_Percentage : Percenate of Base 'A' in the gRNA
* T_Percentage : Percenate of Base 'T' in the gRNA


The results/all_alignments and results/mismatches will have one file for each gRNA. results/all_alignments will have all the alignments of gRNA mapping to each gene/transcript from the TRANSCRIPT_FASTA file and  results/mismatches will have the same alignments filtered for gene/trasncript with total mismatch of specified mismatches. 

* The columns of all_alignments.csv (eg all_alignments_gRNA_1.csv) are described below:
* Sequence ID : The ID of sequence from the input TRANSCRIPT_FASTA file
* gRNA_align_start : Alignment start of the gRNA on the gene/transcript
* gRNA_align_end : Alignment end of the gRNA on the gene/transcript
* seq_align_start : Alignment start of the gene/transcript
* seq_align_end : Alignment end of the gene/transcript
* identity : Number of aligned bases inlcudin gaps in gRNA
* total_aligned_len : Number of bases which aligned between the two sequences
* missing_start : The number of unaligned bases at the start of gRNA 
* missing_end : The number of unaligned bases at the end of gRNA 
* added_aligned_len : This is the total_aligned_len + missing_start + missing_end (Which represnts a hypothetical aligned + unaligned length)
* gRNA_gaps : The number of gaps in the gRNA for this alignment
* seq_gaps : The number of gaps in the gene/transcript for this alignment
* total_mismatches : This is added_aligned_len - identity ( Which is the difference of aligned bases and the total alignment including the unaligned length)

Note: Please feel free to compute your our own mismatch score from files in results/all_alignments if needed.

For more information contact Huzaifa hassan at hhassan@stowers.org