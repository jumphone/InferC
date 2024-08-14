
import sys
import subprocess as sp

A_PATH={}
F_PATH={}
for line in open('LST_ADULT_ID.txt'):
    seq=line.rstrip().split('\t')
    A_PATH[seq[0]]=seq[1]
for line in open('LST_FETAL_ID.txt'):
    seq=line.rstrip().split('\t')
    F_PATH[seq[0]]=seq[1]


fi=open('COM_ID.txt')
for line in fi:
    seq=line.rstrip().split('\t')
    seq_a=seq[1].split(',')
    seq_f=seq[2].split(',')
    this_name=seq[0]
    PHYLOP='/home/database/annotation/hg38/hg38.phyloP100way.bw'
    BED='/home/database/data/COM_scATAC/data/'+this_name+'.bed'
    bwAve='/home/toolkit/tools/ucsc/bigWigAverageOverBed'
    addID='/home/database/data/COM_scATAC/src/addID.py'
    sp.Popen('python '+addID+' '+BED+' '+BED+'.withID.bed',shell=True).wait()    
    sp.Popen(bwAve+' '+PHYLOP+' '+BED+'.withID.bed '+BED+'.hg38.phyloP100way.txt',shell=True).wait()    




