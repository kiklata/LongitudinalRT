## in velocyto env

set -e # Exit immediately if a command exits with a non-zero status

datapath=/home/zhepan/Project/MultiOmics/data/snRNA/P1013S2/Rawdata/

ref=/home/zhepan/Reference/Cellranger

cd $datapath/cellranger

fastq=$datapath/P1013S2

cellranger count --id=P1013S2 \
                 --sample='P1013S2-1','P1013S2-2' \
                 --localcores=16 \
                 --localmem=180 \
                 --transcriptome=$ref \
                 --nosecondary --fastqs=$fastq

