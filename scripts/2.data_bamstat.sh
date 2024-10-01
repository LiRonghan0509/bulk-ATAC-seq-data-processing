#!/bin/bash

# output file name
output_file="bam_mapped_stats_summary.csv"

echo "Filename,Total_Reads,Total_Mapped_Reads,Properly_Paired,Average_Coverage" > "$output_file"

for bam_file in *.bam; do
  # samtools flatstat
  stats=$(samtools flagstat "$bam_file")

  total_reads=$(echo "$stats" | grep "in total" | cut -d ' ' -f 1)
  mapped_reads=$(echo "$stats" | grep "mapped (" | cut -d ' ' -f 1)
  properly_paired=$(echo "$stats" | grep "properly paired (" | cut -d ' ' -f 1)

  # read
  samtools depth "$bam_file" | awk '{sum+=$3} END {print sum/NR}' > tmp_coverage.txt
  avg_coverage=$(cat tmp_coverage.txt)

  # calculate properly mapping ratio
  properly_map_ratio=$(awk "BEGIN {printf \"%.4f\", $properly_paired/$total_reads}")

  # write
  echo "$bam_file,$total_reads,$mapped_reads,$properly_paired,$avg_coverage,$properly_map_ratio" >> "$output_file"
  rm tmp_coverage.txt
done

echo "Done stat!"
