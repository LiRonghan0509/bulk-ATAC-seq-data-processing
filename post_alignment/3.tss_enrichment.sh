# Define paths
bw_dir="/storage1/fs1/yeli/Active/l.ronghan/projects/0.postalign/bw_files/"
tss_dir="/storage1/fs1/yeli/Active/l.ronghan/projects/0.postalign/tss_files/"
fil_dir="/storage1/fs1/yeli/Active/l.ronghan/projects/1.bulkATAC_0901/filtered_bam_files/"
corr_dir="/storage1/fs1/yeli/Active/l.ronghan/projects/0.postalign/correlation_files/"


# Create directories
# mkdir -p $bw_dir
# mkdir -p $tss_dir
mkdir -p $corr_dir
# Convert bam files into bigwig files for sequence peak visualization 
for bam in ${fil_dir}*.bam
do
    sample=$(basename $bam .bam)
    # Convert .bam files into .bigwig
    echo "Start converting $sample"
    bamCoverage -b ${fil_dir}${sample}.bam \
                -o ${bw_dir}${sample}.bw \
                --binSize 10 \
                --normalizeUsing RPKM \
                -p=max  
                # PKM normalizes the coverage to RPKM, which is commonly used for ATAC-seq data. 
done

for bam in ${fil_dir}*.bam
do
    sample=$(basename $bam .bam)
    TSS enrichment analysis
    echo "Start TSS enrichment analysis for $sample"
    computeMatrix reference-point \
                --referencePoint TSS \
                -S ${bw_dir}${sample}.bw \
                -R /storage1/fs1/yeli/Active/l.ronghan/data/hg38/hg38.refGene.TSS.bed \
                -a 2000 \
                -b 2000 \
                --skipZeros \
                -p=max \
                -o ${tss_dir}${sample}.matrix_TSS.gz \
                --binSize 10 \
                
    plotHeatmap -m ${tss_dir}${sample}.matrix_TSS.gz \
                -out ${tss_dir}${sample}.TSS_enrichment.png \
                --colorMap RdYlBu \
                --zMin 0 \
                --zMax 15 \
                --heatmapHeight 7 \
                --heatmapHeight 10 \
                --plotTitle 'TSS Enrichment' \
                --xAxisLabel 'TSS' 
                
    echo "Analysis for $sample is done"
done
echo "TSS enrichment analysis is done"

# Correlation analysis

echo "Starting correlation analysis"
multiBigwigSummary bins -b ${bw_dir}*.bw \
                        -o ${corr_dir}results.npz \
                        --binSize 10 \
# Compute and plot correlation heatmap
plotCorrelation -in ${corr_dir}results.npz \
                --corMethod pearson \
                --skipZeros \
                --plotHeatmap \
                -o ${corr_dir}heatmap_correlation.png \
                --outFileCorMatrix ${corr_dir}correlation_matrix.txt
echo "Correlation analysis completed!"
