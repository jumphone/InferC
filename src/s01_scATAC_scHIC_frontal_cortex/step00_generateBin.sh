

CHROM_SIZE=/home/database/reference/hg19/hg19.fa.size

BINSIZE=100000
cooler makebins $CHROM_SIZE $BINSIZE > $CHROM_SIZE\.100k_coolBin.bed

BINSIZE=10000
cooler makebins $CHROM_SIZE $BINSIZE > $CHROM_SIZE\.10k_coolBin.bed

