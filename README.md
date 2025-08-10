# TingChen_NATURE_2025
## Bulkseq analysis_mRNA and retroelement analysis 
### System Requirements
Both Python and R are required to replicate the analysis pipeline. We used:
Python 3.11.1
R  4.2.2
in our analyses.
### Installation guide
install.packages(c("ggplot2", "stringr", "clusterProfiler", "DESeq2", "org.Mm.eg.db", "openxlsx", "dplyr","readr"))
ggplot2 ：3.5.1
stringr ：1.5.1
clusterProfiler ：4.6.2
DESeq2 ：1.38.3
org.Mm.eg.db ：3.16.0
openxlsx ：4.2.5.2
dplyr ：1.1.4
readr :2.1.4

### Work pepline
1.Quality Control
trim_galore -q 25 --phred33 --length 30 --stringency 3 --fastqc --gzip --paired DT_2_HFSC.1.fq DT_2_HFSC.2.fq -o example/output/1.cleandata/DT_2_HFSC

2.Alinement
NDEX=/reference/mm10_index
STAR --runThreadN 20 \
  --genomeDir ${INDEX} \
  --readFilesIn DT_2_HFSC.1.fq DT_2_HFSC.2.fq \
  --outFileNamePrefix "$output/2.mapping/DT_2_HFSC/DT_2_HFSC_" \
  --readFilesCommand zcat \
  --outFilterMultimapNmax 20 \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \


3.Expression Gene Quantification
#genome
GTF=ref/mm10.refGene.gtf
#retroelement
GTF=HOMER/gEVE/Mouse/Mmus38.geve.v1.gtf
featureCounts -T 5 -t exon -g gene_id --countReadPairs -p -a $GTF -o example/output/counts/DT_2_HFSC/DT_2_HFSC.featurecounts.txt  example/output/2.mapping/DT_2_HFSC/Aligned.sortedByCoord.out.bam

4.Differential Gene Expression Analysis between Different Groups  ,work pipline in bulkRNA.R .

## Sc-RNAseq analysis 
### System Requirements
RStudio 2023.06.2+561
R 4.4.2
### Installation guide
#Installation Guide for Required R Packages
#To install the packages needed for our single - cell analysis workflow, follow these steps from an R terminal:
#Run the following command in an R terminal to install all required packages from CRAN:
install.packages(c("Seurat","dplyr","patchwork","tidyverse","DoubletFinder","readxl","ggplot2"))
### Instruction for use  
The analysis workflow is divided into three separate R scripts, each corresponding to a specific cell type analysis. All scripts contain detailed inline comments to explain key steps and parameters:
scRNAseq_total cell.R
scRNAseq_immune cell.R
scRNAseq_keratinocyte.R
