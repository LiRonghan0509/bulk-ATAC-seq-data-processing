#!/bin/bash

# TN5 shifting of tagaligns for ATACseq

mkdir tagalign_bed_files
mkdir tn5_bed_files

for bam in filtered_bam_files/*.bam
do
    prefix=$(basename $bam .bam)
    bedtools bamtobed -i $bam > tagalign_bed_files/${prefix}.tagAlign.bed
    awk 'BEGIN {OFS = "\t"} {if ($6 == "+") $2 = $2 + 4; else if ($6 == "-") $3 = $3 - 5; print}' tagalign_bed_files/${prefix}.tagAlign.bed > tn5_bed_files/${prefix}.tn5.tagAlign.bed
    gzip tn5_bed_files/${prefix}.tn5.tagAlign.bed
done


# Split .bam files into .chr.bam files
samples=("$@")
output_dir="split_bam_files"
mkdir -p $output_dir
for sample in "${samples[@]}"; do
    echo "Indexing $chr... in $sample"
    chromosomes=$(samtools idxstats filtered_bam_files/$sample.bam | cut -f 1 | grep -v '\*')
    echo "Done indexing $chr... in $sample"
    for chr in $chromosomes; do
        echo "Processing $chr... in $sample"
        output_bam="$output_dir/${sample%.bam}_${chr}.bam"
        samtools view -b filtered_bam_files/$sample.bam $chr > $output_bam

        samtools index $output_bam
    done

done
