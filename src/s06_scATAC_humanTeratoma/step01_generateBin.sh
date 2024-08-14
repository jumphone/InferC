CHROM_SIZE=/home/database/reference/hg19/hg19.fa.size
OUT=./CHROM/

RES=10000
TAG=10k
bedtools makewindows -g $CHROM_SIZE -w $RES > $OUT\/hg19.$TAG\_bin.bed

