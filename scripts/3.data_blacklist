samples=("$@")
blacklist="/storage1/fs1/yeli/Active/l.ronghan/data/hg38/hg38-blacklist.v2.bed"
# get the blacklist gene files: https://github.com/Boyle-Lab/Blacklist/tree/master/lists
  
# mkdir filtered_bam_files
for sample in "${samples[@]}"; do
    echo "Start processing ${sample}"
    bedtools intersect -v -abam bam_files/${sample}.sorted.filt.noMT.bam -b ${blacklist} > filtered_bam_files/${sample}.bam 
done
