#!/bin/bash

current_dir=$(pwd)
parent_dir=$(dirname "$current_dir")
output=$(find "$parent_dir" -mindepth 1 -maxdepth 1 -type d -name "output")
reference=$(find "$parent_dir" -mindepth 1 -maxdepth 1 -type d -name "reference")
INDEX=$reference/mm10_index/

cat fq.liist  | while read -r line; do
id=$(echo "$line" | cut -d '	' -f1)  # 获取样本
fq1=$(echo "$line" | cut -d '	' -f2) # 获取fastq文件路径
fq2=$(echo "$line" | cut -d '	' -f3) # 获取fastq文件路径
mkdir -p $output/1.cleandata/$id
cd $output/1.cleandata/$id
trim_galore -q 25 --gzip --phred33 --length 30 --stringency 3  --fastqc --gzip   --paired $fq1 $fq2 -o $output/1.cleandata/$id


fq1=$output/1.cleandata/$id/${id}.1_val_1.fq.gz
fq2=$output/1.cleandata/$id/${id}.2_val_2.fq.gz
mkdir -p $output/2.mapping/$id
cd $output/2.mapping/$id
STAR --runThreadN 20 \
  --genomeDir ${INDEX} \
  --readFilesIn $fq1 $fq2 \
  --outFileNamePrefix "$output/2.mapping/$id/$id_" \
  --readFilesCommand zcat \
  --outFilterMultimapNmax 20 \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \


mkdir -p $output/3.featurecounts
GTF_1=/public/home/changying/project/mouse_bulk_RNA_cuijun_PD1_antivirus_2024_12_20/ref/mm10.refGene.gtf
GTF_2=/public/home/changying/softwares/HOMER/gEVE/Mouse/Mmus38.geve.v1.gtf

featureCounts -T 5 -t exon -g gene_id --countReadPairs -p -a $GTF_1 -o $output/3.featurecounts/${id}.mRNA.featurecounts.txt  $output/2.mapping/$id/Aligned.sortedByCoord.out.bam
awk 'BEGIN{FS=OFS="\t"}!/#/{print $1,$6,$7}' $output/3.featurecounts/${id}.mRNA.featurecounts.txt > $output/3.featurecounts/${id}.mRNA.count.txt

featureCounts -T 5 -t CDS -g gene_id --countReadPairs -p -a $GTF_2 -o $output/3.featurecounts/${id}.erv.featurecounts.txt $output/2.mapping/$id/Aligned.sortedByCoord.out.bam
awk 'BEGIN{FS=OFS="\t"}!/#/{print $1,$6,$7}' $output/3.featurecounts/${id}.erv.featurecounts.txt > $output/3.featurecounts/${id}.erv.count.txt
done

for  id in 6OHDA_DT_1_HFSC  6OHDA_DT_3_HFSC  PBS_HFSC_a11 PBS_HFSC_a12 ;do
echo "$id	$output/3.featurecounts/${id}.erv.count.txt" >>featurecount.file.txt
done

for  id in 6OHDA_DT_1_HFSC  6OHDA_DT_3_HFSC  PBS_HFSC_a11 PBS_HFSC_a12 ;do
echo "$id	$output/3.featurecounts/${id}.mRNA.count.txt" >>featurecount.file.txt
done

# 运行 bulkRNA.R 脚本
Rscript bulkRNA.R > bulkRNA.R.log 2>&1 &&

