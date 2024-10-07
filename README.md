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

```{bash}
ronghanli2002/atacseqqc:hg38
```
  ATACseqQC requires a very large amount of CPU memory. (Every .bam file~50GB. Even when CPU mem=512GB, sessions were killed as well). Using `0.data_splitbam.sh` to split the whole .bam file into small .chr.bam files
- Peak calling and annotation
