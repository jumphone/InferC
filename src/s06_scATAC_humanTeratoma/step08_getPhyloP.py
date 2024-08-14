
import sys
import subprocess as sp


PHYLOP='/home/database/annotation/hg19/hg19.100way.phyloP100way.bw'
BED='../data/ALL_PEAK.bed'
bwAve='/home/toolkit/tools/ucsc/bigWigAverageOverBed'
addID='./addID.py'
sp.Popen('python '+addID+' '+BED+' '+BED+'.withID.bed',shell=True).wait()
sp.Popen(bwAve+' '+PHYLOP+' '+BED+'.withID.bed '+BED+'.hg19.phyloP100way.txt',shell=True).wait()

