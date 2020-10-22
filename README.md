# SARS-CoV-2-leader

This repository is for the scripts to find reads containing partial leader
sequence in the SARS-CoV-2 genome. The scripts here are writen in bash and 
awk for use on a linux system.

The scripts take a mapped SAM/BAM file and look for reads with the sequence
motif GTAGATCTGTTCTCT, which occurs at the 3' end of the SARS-CoV-2 leader.

You will need to install samtools on the system for the scripts to work 
correctly. 
