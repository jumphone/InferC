
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
    all_bed_file=[]
    for one in seq_a:
        all_bed_file.append(A_PATH[one]+'/peaks_macs2.rds.bed')
    for one in seq_f:
        all_bed_file.append(F_PATH[one]+'/peaks_macs2.rds.bed')
    bedtools='/home/toolkit/local/bin/bedtools'
    print(all_bed_file)
    sp.Popen('cat '+' '.join(all_bed_file)+' | '+bedtools+' sort -i - | '+bedtools+' merge -i - > '+'/home/database/data/COM_scATAC/data/'+this_name+'.bed',shell=True).wait()
        




