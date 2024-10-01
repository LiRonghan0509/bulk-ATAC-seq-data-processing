samples=("$@")

reference="/storage1/fs1/yeli/Active/yeli/genome/hg38/GRCh38.p12.genome.fa"

# Function to perform FastQC on a pair of files
qc_files() {
    fastqc "$1_R1_001.fastq.gz"
    fastqc "$1_R2_001.fastq.gz"
}

# Function to trim files
trim_files() {
    trim_galore --paired -o trimmed/ "$1_R1_001_val_1.fq.gz" "$1_R2_R2_001_val_2.fq.gz"
}

# QC raw reads
echo "Performing QC on raw reads..."
for sample in "${samples[@]}"; do
    qc_files $sample
done

# Trim reads
echo "Trimming reads..."
for sample in "${samples[@]}"; do
    trim_files $sample
done

# QC trimmed reads
echo "Performing QC on trimmed reads..."
for sample in "${samples[@]}"; do
    # qc_files "${sample}_trimmed"
    qc_files trimmed/"${sample}"
done

Combine QC data
python multiqc . # in the fastqc output directory; repeat for raw and trimmed
echo "Performing MultiQC..."
multiqc .

# Align to reference genome
# first have to prepare genome index; if not first, please remove this line
echo "Mapping to reference genome..."
bwa index $reference

# Align to reference genome
for sample in "${samples[@]}"; do
    # BWA alignment (assuming you are using BWA MEM, which is suitable for longer reads)
    bwa mem -t 8 $reference trimmed/${sample}_R1_001_val_1.fq.gz trimmed/${sample}_R2_001_val_2.fq.gz | \
    samtools view -bS - > ${sample}.bam

    # Sort and index BAM files
    samtools sort -o ${sample}.sorted.bam ${sample}.bam
    samtools index ${sample}.sorted.bam

    # Remove mitochondrial DNA reads
    samtools view -h ${sample}.sorted.bam | python removeChrom.py - - chrM | samtools view -bh - > ${sample}.noMT.bam

    # Coordinate sort and index the bam without mitochondrial DNA
    samtools sort -o ${sample}.sorted.noMT.bam ${sample}.noMT.bam
    samtools index ${sample}.sorted.noMT.bam

    # Restrict to properly-paired reads only
    samtools view -bh -f 2 ${sample}.sorted.noMT.bam > ${sample}.filt.noMT.bam

    # Sort and index the filtered bam
    samtools sort -o ${sample}.sorted.filt.noMT.bam ${sample}.filt.noMT.bam
    samtools index ${sample}.sorted.filt.noMT.bam

done
