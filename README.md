# Bulk ATAC-seq data preprocessing pipeline
A pipeline for bulk ATAC-seq data preprocessing. From .fasta files to peak calling. Using the docker image `ronghanli2002/bulkatac:&lt;tag`
## Workflow for bulk ATAC-seq data processing and analysis
![atac_analysis](https://github.com/user-attachments/assets/f7ef2dca-341a-47ac-b15a-5900c279d83a)
Created in  https://BioRender.com
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
ronghanli2002/atacseqqc:hg38
```

```{bash}
bsub -Is -q general-interactive -G compute-yeli -R 'rusage[mem=128GB]' -a 'docker(ronghanli2002/atacseqqc:hg38)' /bin/bash    
```
  ATACseqQC requires a very large amount of CPU memory. (Every .bam file~50GB. Even when CPU mem=512GB, sessions were killed as well). Using `0.data_splitbam.sh` to split the whole .bam file into small .chr.bam files.
  For the splited .chr.bam files, they just require `CPUmem=128GB`, and use the docker image `ronghanli2002/atacseqqc:hg38`, `R=4.2.1`, `Bioconductor=3.16`.
  Required packages:
```{R}
.libPaths("/usr/lib/R/site-library")
library(BiocManager)
library(ATACseqQC)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(phastCons100way.UCSC.hg38)
library(Rsamtools)
library(ChIPpeakAnno)
```
- Peak calling and annotation
