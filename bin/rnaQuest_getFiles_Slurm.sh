#!/bin/bash

#SBATCH -p general  ## Partition
#SBATCH -q public  ## QOS
#SBATCH -N 1   ## Number of Nodes
#SBATCH --cpus-per-task=16
#SBATCH --time=900
#SBATCH --job-name=rnaQuest_getFiles
#SBATCH -o /rnaQuest/bin/out/rnaQuest_getFiles.%j.out
#SBATCH -e /rnaQuest/bin/out/rnaQuest_getFiles.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=USERNAME@EMAIL.edu

# This script will walk you through downloading the fastq files for the tutorial to your linux computer locally for downstream processing. 
# For this script each line should be copied into the terminal, line by line, for execution. An example executable bash script is also 
# available in the github repository if you have access to an HPC environment with SLURM (rnaQuest_getFiles_Slurm.sh).

# Be sure to make script executable using: $chmod a+x runQC.sh. Long-term storage of scripts should be in the /home/rpevey/tdpAD/bin/ directory.
# If you wrote your script in a windows environment convert files.txt to unix with $dos2unix files.txt.
# For the purposes of this tutorial I have created a project directory with three subdirectories; bin/, data/, results/. Unless stated 
# otherwise, the following code is executed in the data directory.
# bin/ has all of the script files and reference files necessary for their execution. 
# data/ houses all of the raw, intermediate, and processed data files for the project. 
# results/ is the location for any output results files for downstream analysis.
# Check on job using:$ squeue --job ########
# Check on job memory usage and status:$ seff ########

# Load necessary modules
module load sratoolkit
module load parallel

# Prefetch SRR accessions
# Configure prefetch to save to the current working directory. Move to desired directory for download.
#vdb-config --prefetch-to-cwd
cat ../SRR_Acc_List.txt | cut -d "," -f 1 > ../SRR.numbers

# Run prefetch in parallel. The -v tag produces verbose execution messages, you can exclude that tag by preference. 
# The job waits until the lines ending in '&' have all completed before moving onto the next part of the script.
cat ../SRR.numbers | parallel prefetch -v {} &
wait
# Or you can download each file individually
# prefetch SRR10446596
# prefetch SRR10446599
# prefetch SRR10446602
# prefetch SRR10446605
# prefetch SRR10446608
# prefetch SRR10446611

# Convert the file into fastq for downstream processing, this step may take awhile.
cat ../SRR.numbers | parallel fasterq-dump {}
# Or you can convert each file individually
# fasterq-dump SRR10446596
# fasterq-dump SRR10446599
# fasterq-dump SRR10446602
# fasterq-dump SRR10446605
# fasterq-dump SRR10446608
# fasterq-dump SRR10446611
