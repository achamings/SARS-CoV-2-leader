# SARS-CoV-2-leader

This repository is for the scripts to find reads containing partial leader
sequence in the SARS-CoV-2 genome. The scripts here are writen in bash and 
awk for use on a linux system.

The scripts take a mapped SAM/BAM file and look for reads with the sequence
motif GTAGATCTGTTCTCT, which occurs at the 3' end (nucleotides 52-67) 
of the SARS-CoV-2 leader.

You will need to install samtools on the system for the scripts to work 
correctly. On Debian this is done by: sudo apt install samtools

The script will output a BAM and SAM file containing only the reads with 
the above motif. The script will also produce a summary csv file on the 
total number of reads, the number of reads containing the leader, how many 
of those reads belong to the 5'UTR and then ratios of leader reads between 
the 5' UTR and subgenomic regions.

It will also produce a csv file using samtools depth to show the coverage
at each position at the specified Q value.

