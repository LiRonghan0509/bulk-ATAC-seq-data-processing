setwd("/storage1/fs1/yeli/Active/l.ronghan/projects/1.bulkATAC_0901")
.libPaths("/usr/lib/R/site-library")
library(BiocManager)
library(ATACseqQC)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(phastCons100way.UCSC.hg38)

library(Rsamtools)
library(ChIPpeakAnno)

outPath <- "splited"
dir.create(outPath)

possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                 "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                               "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                               "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                               "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                               "U2"))

bamfile_list

for (bamfile in bamfile_list) {
  bamfile.labels <- gsub(".bam", "", basename(bamfile))
  source(system.file("extdata", "IGVSnapshot.R", package = "ATACseqQC"))
  # bamQC(bamfile, outPath=NULL)
  # dupFreq <- readsDupFreq(bamfile)
  pdf(paste0(bamfile.labels,"library_complexity.pdf"))
  # estimateLibComplexity(dupFreq)
  estimateLibComplexity(readsDupFreq(bamfile))
  dev.off()

  pdf(paste0(bamfile.labels,"fragsize.pdf"))
  fragSize <- fragSizeDist(bamfile, bamfile.labels)
  dev.off()
  bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                     param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
  tags <- names(bamTop100)[lengths(bamTop100)>0] 

  seqlev <- "chr1" ## subsample data for quick run
  seqinformation <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
  which <- as(seqinformation[seqlev], "GRanges")
  gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
  # paste0(bamfile.labels,"fragsize.pdf"))
  shiftedBamfile <- file.path(paste0("splited/", bamfile.labels,"_shifted.bam"))
  gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)

  txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
  pt <- PTscore(gal1, txs)

  pdf(paste0(bamfile.labels,"_meancoverage.pdf"))
  plot(pt$log2meanCoverage, pt$PT_score, 
      xlab="log2 mean coverage",
      ylab="Promoter vs Transcript") 
  dev.off()

  nfr <- NFRscore(gal1, txs)
  pdf(paste0(bamfile.labels,"_nfrscore.pdf"))
  plot(nfr$log2meanCoverage, nfr$NFR_score, 
      xlab="log2 mean coverage",
      ylab="Nucleosome Free Regions score",
      main="NFRscore for 200bp flanking TSSs",
      xlim=c(-10, 0), ylim=c(-5, 5))
  dev.off()

  tsse <- TSSEscore(gal1, txs)
  tsse$TSSEscore
  pdf(paste0(bamfile.labels,"_tssescore.pdf"))
  plot(100*(-9:10-.5), tsse$values, type="b", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")
  dev.off()
  
  txs <- txs[seqnames(txs) %in% "chr1"]
  genome <- Hsapiens

  objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath,
                              conservation=phastCons100way.UCSC.hg38)
  dir(outPath)

  objs <- splitBam(bamfile, tags=tags, outPath=outPath,
                 txs=txs, genome=genome,
                 conservation=phastCons100way.UCSC.hg38)
  objs <- splitGAlignmentsByCut(gal1, txs=txs, outPath = outPath)

  bamfiles <- file.path(outPath,
                     c("NucleosomeFree.bam",
                     "mononucleosome.bam",
                     "dinucleosome.bam",
                     "trinucleosome.bam"))
  cumulativePercentage(bamfiles[1:2], as(seqinformation["chr1"], "GRanges"))
  TSS <- promoters(txs, upstream=0, downstream=1)
  TSS <- unique(TSS)
  ## estimate the library size for normalization
  (librarySize <- estLibSize(bamfiles))

  
}

