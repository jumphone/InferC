PHYLOP='/home/database/annotation/hg19/hg19.100way.phyloP100way.bw'
BED='../data/GSE159557_RAW/ALL_PEAK.bed'
bwAve='/home/toolkit/tools/ucsc/bigWigAverageOverBed'
addID='/home/database/data/InferX/data/humanRetinaDev/scATAC/src/addID.py'

python $addID $BED $BED\.withID.bed
$bwAve $PHYLOP $BED\.withID.bed $BED\.hg19.phyloP100way.txt

