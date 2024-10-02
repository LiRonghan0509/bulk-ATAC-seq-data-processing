setwd("/Users/ronghanli/Downloads/peaks")
# Install and load required R libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)

# Load BroadPeak files
peak_files <- c("LIB047011-DIL01_22CYYVLT4_S129_L005.peaks.bed", 
                "LIB047013-DIL01_22CYYHLT4_S747_L007.peaks.bed", 
                "LIB047012-DIL01_22C7G2LT4_S1_L002.peaks.bed", 
                "LIB047014-DIL01_22CYYHLT4_S712_L004.peaks.bed")

# Plot annotated results
for (file in peak_files) {
    # Read peak files
    peaks <- read.table(file, header = FALSE)
    colnames(peaks) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")

    # Convert peaks to GRanges objects
    peaks_gr <- GRanges(seqnames = peaks$chr,
                        ranges = IRanges(start = peaks$start, end = peaks$end))

    # Annotate peaks
    MacsCalls_Anno <- annotatePeak(peaks_gr, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)

    # Anotated results
    print(MacsCalls_Anno)

    # Plot pie chart
    pie_chart_filename <- paste0("peak_annotation_", sub(".bed", "", file), "_pie_chart.jpg")
    jpeg(pie_chart_filename)
    plotAnnoPie(MacsCalls_Anno)
    dev.off()

    # Plot bar chart
    bar_chart_filename <- paste0("peak_annotation_", sub(".bed", "", file), "_bar_chart.jpg")
    jpeg(bar_chart_filename)
    plotAnnoBar(MacsCalls_Anno)
    dev.off()
}
