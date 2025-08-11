# TingChen_NATURE_2025
## Bulkseq analysis_mRNA and retroelement analysis 
### System Requirements
Both Python and R are required to replicate the analysis pipeline. We used:
```
Python 3.11.1
trim_galore v0.6.10
STAR v2.7.10a
subread (FeatureCounts v2.0.6)
R  4.2.2
```
in our analyses.
### Installation guide
```
install.packages(c("ggplot2", "stringr", "clusterProfiler", "DESeq2", "openxlsx", "dplyr","readr"))
```
```
ggplot2 ：3.5.1
stringr ：1.5.1
clusterProfiler ：4.6.2
DESeq2 ：1.38.3
openxlsx ：4.2.5.2
dplyr ：1.1.4
readr :2.1.4
```

On a normal desktop computer with standard internet speed, the installation of these packages typically takes 2 hours.
### Work pepeline

For the specific workflow and demo, please refer to the bulkRNAseq folder.

Use this command to run pipeline.sh：
```
nohup bash pipeline.sh > pipeline.log 2>&1 &
```

1.Quality Control

2.Alinement

For alignment, the mouse mm10 genome assembly was used as the reference. 

3.Expression Gene Quantification

For mRNA analysis, we utilized the mm10 genome (refdata-gex-mm10-2020-A). For retroelement analysis, we employed the Mmus38.geve.v1 annotation (http://geve.med.u-tokai.ac.jp). 

4.Differential Gene Expression Analysis between Different Groups, work pipline in bulkRNA.R.


## Sc-RNAseq analysis 
### System Requirements
```
RStudio 2023.06.2+561
R 4.4.2
```
### Installation guide
```
install.packages(c("Seurat","dplyr","patchwork","tidyverse","DoubletFinder","readxl","ggplot2"))
```
On a normal desktop computer with standard internet speed, the installation of these packages typically takes 2 hours.
### Instruction for use  
The analysis workflow is divided into three separate R scripts, each corresponding to a specific cell type analysis. All scripts contain detailed inline comments to explain key steps and parameters:
scRNAseq_total cell.R;
scRNAseq_immune cell.R;
scRNAseq_keratinocyte.R.
