# Split .bam files into .chr.bam files
#!/bin/bash
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
