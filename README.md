# Bulk ATAC-seq data preprocessing pipeline
A pipeline for bulk ATAC-seq data preprocessing. From .fasta files to peak calling. Using the docker image `ronghanli2002/bulkatac:&lt;tag`

## Pre-alignment and alignment
- Data trimming on the .fastq files
- Cut off adaptors: `trim_galore`
- Quality control: `fastqc`, `multiqc`
- Align to the reference genome:
  ```{bash}
  bwa mem -t 8 $reference trimmed/${sample}_R1_001_val_1.fq.gz trimmed/${sample}_R2_001_val_2.fq.gz | samtools view -bS - > ${sample}.bam
  ```
- Remove mitochondrial reads: `samtools`

## Statistic results of the alignment results
- Count the ratio of properly paired reads: `samtools flagstat`
## Post-alignment analysis
- ATACseqQC
```{bash}
alanayuchen/cut_tag_bulk:4.7
```
- Peak calling and annotation
