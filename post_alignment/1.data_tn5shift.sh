mkdir tagalign_bed_files
mkdir tn5_bed_files

for bam in filtered_bam_files/*.bam
do
    prefix=$(basename $bam .bam)
    bedtools bamtobed -i $bam > tagalign_bed_files/${prefix}.tagAlign.bed
    awk 'BEGIN {OFS = "\t"} {if ($6 == "+") $2 = $2 + 4; else if ($6 == "-") $3 = $3 - 5; print}' tagalign_bed_files/${prefix}.tagAlign.bed > tn5_bed_files/${prefix}.tn5.tagAlign.bed
    gzip tn5_bed_files/${prefix}.tn5.tagAlign.bed
done
