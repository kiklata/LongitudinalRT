#!/bin/bash
# STAR env
set -e # Exit immediately if a command exits with a non-zero status

datapath=/home/zhepan/Project/MultiOmics/data/snRNA/P1013S2/Rawdata

refer=/home/zhepan/Reference/STAR
whitelist=/home/zhepan/Reference/STAR/3M-february-2018.txt

bam=$datapath/cellranger/P1013S2/outs/possorted_genome_bam.bam
mkdir -p $datapath/star
STAR --runThreadN 8 \
     --genomeDir $refer \
     --soloType CB_UMI_Simple \
     --readFilesType SAM SE --readFilesIn $bam  \
     --readFilesCommand samtools view -F 0x100 \
     --soloInputSAMattrBarcodeSeq CR UR \
     --soloInputSAMattrBarcodeQual CY UY \
     --soloCBwhitelist $whitelist \
     --soloCellFilter EmptyDrops_CR \
     --soloFeatures Gene SJ Velocyto \
     --outFileNamePrefix $datapath/star \
     --outSAMtype None
