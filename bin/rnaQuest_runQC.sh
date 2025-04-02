#!/bin/bash

# This script will run fastQC and Trimmomatic on every pair of read files listed in SRR_Acc_List.txt, which should be in the projects /bin
# directory. First, it will read in SRR_Acc_List.txt, then access the first line in the file which contains the root of the file names for
# the paired read files. The actual file names are parsed from the root and declared separately as variables to be called later. Second,
# fastQC is run on the raw fastq input files. Third, Trimmomatic is used to trim illumina adapters and low quality reads. The paired and
# unpaired trimmed files are output to the /results/fastQC directory. Fourth, fastQC is run on each of the trimmed read files and output a
# multiQC report for before and after trimming. After completing that on the first set of files, it will then cycle back onto the next set
# of files listed in SRR_Acc_List.txt until all files have been trimmed and QC'd.

# Be sure to make script executable using: '$chmod a+x runQC.sh'. If the script was written on a windows computer convert it to unix with
# '$dos2unix rnaQuest_runQC.sh'.
# Use nohup to run the script just incase you get disconnected.
	# nohup bash ../bin/rnaQuest_runQC.sh > ../bin/out/runQC.out 2>&1 &
# if you get disconnected check on progress with 
	# tail -f output.log
		###OR###
 	# ps aux | grep my_script.sh

# Download fastqc
# sudo apt install fastqc
# sudo apt install multiqc

# I like running fastQC on one preliminary file to make sure that I have the proper adapater settings for trimmomatic.
# fastqc SRR10446596_1.fastq.gz -o ../results/fastQC/raw/
# Observe the output file to check the adapaters. It should be Illumina family adapters according to metadata from NCBI.

# Copy the relevant adapter file to the working directory for Trimmomatic
# cp ../data/TruSeq3-PE.fa TruSeq3-PE.fa ./
# nohup java -jar ~/refs/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 \
		# SRR10446596_1.fastq.gz SRR10446596_2.fastq.gz \
		# ../data/trimmed/SRR10446596_1_trimmed_paired.fastq.gz ../data/trimmed/SRR10446596_1_trimmed_unpaired.fastq.gz \
		# ../data/trimmed/SRR10446596_2_trimmed_paired.fastq.gz ../data/trimmed/SRR10446596_2_trimmed_unpaired.fastq.gz \
		# ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
		# LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>&1 &
		
# Declare files variable with file name roots
files="rnaQuest/bin/SRR_Acc_List.txt"
while IFS= read -r line
do
	########
        # Step 1
        ########
	# Declare file names for the files that are to be trimmed and QC'd.
	ROOT=$line
	TAIL1="_1.fastq.gz"
	TAIL2="_2.fastq.gz"
	SEQUENCE1=$ROOT$TAIL1
	echo $SEQUENCE1
	SEQUENCE2=$ROOT$TAIL2
        echo $SEQUENCE2
	# Declare output directory for trimmed fastq files and output file names.
	OUTQCDIR=rnaQuest/results/fastQC/
	OUTPUTDIR=rnaQuest/data/trimmed/

	TAIL1pair="_1_trimmed_paired.fq"
	TAIL1unpair="_1_trimmed_unpaired.fq"
	TAIL2pair="_2_trimmed_paired.fq"
	TAIL2unpair="_2_trimmed_unpaired.fq"

	OUTPUTFILE1paired=$ROOT$TAIL1pair
	OUTPUTFILE1unpaired=$ROOT$TAIL1unpair
	OUTPUTFILE2paired=$ROOT$TAIL2pair
	OUTPUTFILE2unpaired=$ROOT$TAIL2unpair

	echo $OUTPUTDIR$OUTPUTFILE1paired
	echo $OUTPUTDIR$OUTPUTFILE2paired

	########
	# Step 2
	########
	# Run fastQC on each file.
	fastqc $SEQUENCE1 -o $OUTQCDIR/raw/
	fastqc $SEQUENCE2 -o $OUTQCDIR/raw/

	########
	# Step 3
	########
	# Run Trimmomatic on the raw sequence data.
	java -jar $HOME/refs/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 \
		$SEQUENCE1 $SEQUENCE2 \
		$OUTPUTDIR$OUTPUTFILE1paired $OUTPUTDIR$OUTPUTFILE1unpaired \
		$OUTPUTDIR$OUTPUTFILE2paired $OUTPUTDIR$OUTPUTFILE2unpaired \
		ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

	########
	# Step 4
	########
	# Run fastQC on the trimmed data.
	fastqc $OUTPUTDIR$OUTPUTFILE1paired -o $OUTQCDIR/trimmed/
	fastqc $OUTPUTDIR$OUTPUTFILE2paired -o $OUTQCDIR/trimmed/

	echo "$line"
done < "$files"

# Run multiQC just once on each output fastQC directories to consolidate the fastQC output to one file.
	multiqc -o $OUTQCDIR/raw/ $OUTQCDIR/raw/
	multiqc -o $OUTQCDIR/trimmed/ $OUTQCDIR/trimmed/