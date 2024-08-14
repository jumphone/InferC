
PHYLOP='/home/database/annotation/hg38/hg38.phyloP100way.bw'
BED='../data/GSE184386/ALL_PEAK.bed'
bwAve='/home/toolkit/tools/ucsc/bigWigAverageOverBed'
addID='/home/database/data/InferX/data/humanRetinaDev/scATAC/src/addID.py'

python $addID $BED $BED\.withID.bed
$bwAve $PHYLOP $BED\.withID.bed $BED\.hg38.phyloP100way.txt


