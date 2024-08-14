

fi=open('LST.txt')


all_bed_file=[]
for line in fi:
    seq=line.rstrip().split('\t')
    this_name=seq[0]+'.format'
    all_bed_file.append(this_name+'/peaks_macs2.rds.bed')


print(all_bed_file)


bedtools='/home/toolkit/local/bin/bedtools'

import subprocess as sp

sp.Popen('cat '+' '.join(all_bed_file)+' | '+bedtools+' sort -i - | '+bedtools+' merge -i - > '+'../data/ALL_PEAK.bed',shell=True).wait()




