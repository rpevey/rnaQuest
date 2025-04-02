#!/bin/bash

# The following block of code will provide an estimate of the sequencing depth and output it in a separate file named depth.txt.
# This script calculates sequencing depth statistics for all sorted and indexed BAM files in the ../results/bam/ directory.
# Outputs: Mean depth, and percent coverage for each sample, then both of the same statistics for all of the samples overall,
# then it outputs those results in 'depth.txt'. The samtools depth command takes a while to run, so be sure to 

# Be sure to make script executable using: '$chmod a+x runQC.sh'. If the script was written on a windows computer convert it to unix with
# '$dos2unix rnaQuest_alignSortIndex.sh'.
# Use nohup to run the script just incase you get disconnected.
	# nohup bash ../bin/rnaQuest_readDepth.sh > ../bin/out/readDepth.out 2>&1 &
# if you get disconnected check on progress with 
	# tail -f output.log
		###OR###
 	# ps aux | grep my_script.sh

# Output file
OUTPUT_FILE="../results/depth.txt"
echo -e "Sample\tMean Depth\t% Coverage" > "$OUTPUT_FILE"


# Temporary file for collecting all depth values across samples
ALL_DEPTHS="all_depths.tmp"
> "$ALL_DEPTHS"  # Empty the file in case it exists

# Process each BAM file
for file in ../results/bam/*.bam; do
    echo "Processing $file..."
    
    # Ensure the BAM file is indexed
    if [ ! -f "${file}.bai" ]; then
        echo "Index file missing for $file. Skipping..."
        continue
    fi
    
    # Extract coverage statistics
    coverage_stats=$(samtools coverage "$file" | awk 'NR>1 {sum+=$7; total+=$6} END {print sum/NR, (total/NR)*100}')
    
    mean_depth=$(echo "$coverage_stats" | awk '{print $1}')
    percent_coverage=$(echo "$coverage_stats" | awk '{print $2}')
    
    # Store mean depth for overall statistics
    echo "$mean_depth" >> "$ALL_DEPTHS"
    
    # Append per-sample results to the output file
    echo -e "$file\t$mean_depth\t$percent_coverage\t$max_depth" >> "$OUTPUT_FILE"
	echo "Completed processing file: $file"
done

# Compute overall statistics from all samples
echo "Calculating overall statistics..."
overall_mean=$(awk '{sum+=$1} END {print sum/NR}' "$ALL_DEPTHS")

# Append overall statistics to the output file
echo -e "Overall\t$overall_mean\t-" >> "$OUTPUT_FILE"

# Cleanup temporary file
rm -f "$ALL_DEPTHS"

echo "All read depth statistics saved to $OUTPUT_FILE"