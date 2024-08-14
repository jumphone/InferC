BED=./rds/scmul_pbmc_usedType_atacMat.rds.bedAll.bed
python src/addID.py $BED $BED\.withID.bed

PHYLOP=/home/database/annotation/hg38/hg38.phyloP100way.bw
bigWigAverageOverBed $PHYLOP $BED\.withID.bed $BED\.hg38.phyloP100way.txt

PHAST=/home/database/annotation/hg38/hg38.phastCons100way.bw
bigWigAverageOverBed $PHAST $BED\.withID.bed $BED\.hg38.phastCons100way.txt

