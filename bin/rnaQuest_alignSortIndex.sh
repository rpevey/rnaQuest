#!/bin/bash

# This is an executable script with commands following the demo (https://github.com/alexdobin/STAR). First, files must be downloaded and 
# aligned using STAR (for short-reads, GMAP for long reads) to produce .bam files for downstream analysis. If you are looking for more
# control over read counts, you can follow this script up with HTSeq or a similar tool. For the purposes of this tutorial I will stick to
# using the --quantMode GeneCounts function in STAR for efficiency and speed and should produce the same output as HTSeq on defaults.
# The index directory contains indexes for GRCh38 release 113. The genome directory contains one file which is the full human genome
# sequence (GRCh38 from Ensembl, downloaded from Ensembl on 25/03/03). Untar the file before running the script using:
# '$tar xvzf 'file'.tar.gz'. The genes directory contains one file containing human gene an notations for GRCm38 from the Ensembl database
# (Ensembl release 113).

# Be sure to make script executable using: '$chmod a+x runQC.sh'. If the script was written on a windows computer convert it to unix with
# '$dos2unix rnaQuest_alignSortIndex.sh'.
# Use nohup to run the script just incase you get disconnected.
	# nohup bash ../bin/rnaQuest_alignSortIndex.sh > ../bin/out/alignSortIndex.out 2>&1 &
# if you get disconnected check on progress with 
	# tail -f output.log
		###OR###
 	# ps aux | grep my_script.sh

# Index reference genome, execute the following code in ~/refs/index/ directory. This step takes a long time but only needs executed once per project.
# mkdir ~/refs/index
# PATH=$PATH:~/refs/STAR/STAR-2.7.11b/bin/Linux_x86_64_static/
# nohup STAR --runThreadN 8 --runMode genomeGenerate \
     # --genomeDir ~/refs/index/ \
     # --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     # --sjdbGTFfile Homo_sapiens.GRCh38.113.gtf \
     # --sjdbOverhang 49 2>&1 &

### Download Ensemble GRCh38 reference genome (.fasta) and annotation file (.gtf) into index directory
#wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
#gunzip Homo_sapiens.GRCh38.113.gtf.gz

# Test out the software install with the following code before running the script.
# PATH=$PATH:~/refs/STAR/STAR-2.7.11b/bin/Linux_x86_64_static/
# nohup STAR --runThreadN 12 \
	# --readFilesIn trimmed/SRR10446596_1_trimmed_paired.fq trimmed/SRR10446596_2_trimmed_paired.fq \
	# --outFileNamePrefix ./SRR10446596. \
	# --outSAMtype BAM SortedByCoordinate \
	# --quantMode GeneCounts \
	# --genomeDir $HOME/refs/index \
	# --sjdbGTFfile $HOME/refs/index/Homo_sapiens.GRCh38.113.gtf 2>&1 &
	
# samtools index SRR10446596Aligned.sortedByCoord.out.bam


#Add STAR to PATH
PATH=$PATH:~/refs/STAR/STAR-2.7.11b/bin/Linux_x86_64_static/
# Open files.txt with file name roots
files="$HOME/evolio/posts/rnaQuest/bin/SRR_Acc_List.txt"
while IFS= read -r line
do
	########
	# Step 1
	########
	# Stop on error
	#set -uex

	# Assign reference to the variable REF
	ACC=Homo_sapiens.GRCh38.dna.primary_assembly
	REF=$HOME/refs/index/$ACC.fa

	# Declare file names for the files that are to be trimmed and QC'd.
	ROOT=$line
	TAIL1="_1_trimmed_paired.fq"
	TAIL2="_2_trimmed_paired.fq"
	READ1=$ROOT$TAIL1
	echo "Input file 1 is: $READ1"
	READ2=$ROOT$TAIL2
	echo "Input file 2 is: $READ2"
	# Declare output directory for trimmed fastq files and output file names.
	INPUTDIR=trimmed/

	echo $INPUTDIR$READ1
	echo $INPUTDIR$READ2


	# Align the paired end reads to the reference genome
	echo "Starting alignment of $line with STAR using $ACC reference genome."
	STAR --runThreadN 12 \
	--readFilesIn $INPUTDIR$READ1 $INPUTDIR$READ2 \
	--outFileNamePrefix ./$ROOT. \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts \
	--genomeDir $HOME/refs/index \
	--sjdbGTFfile $HOME/refs/index/Homo_sapiens.GRCh38.113.gtf
	echo "Alignment complete, output as BAM file, sorted by coordinate."

	# Use samtools to convert alignment to bam file. Sort alignment by leftmost coordinate.
	#echo "Convert sam to bam."
	#samtools view -bS $ROOT.bwa.sam > $ROOT.bwa.bam
	#echo "Sort bam file by leftmost coordinate."
	#samtools sort $ROOT.bwa.bam > $ROOT.bwa.bam
	# index the BAM file
	echo "Index bam file for fast random access."
	samtools index $ROOT.Aligned.sortedByCoord.out.bam

	echo "Completed processing file: $line"
done < "$files"

#mv processed bam files and counts files to a new folder.
mkdir ../results/bam
mv *.out.bam ../results/bam/
mv *.out.bam.bai ../results/bam/
mkdir ../results/counts
mv *.ReadsPerGene.out.tab ../results/counts/
