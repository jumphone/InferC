
import subprocess
import multiprocessing
import os


def Work(input_path,a):
    print input_path
    subprocess.Popen('Rscript ./RSCRIPT/scATAC_mergeCount.R '+input_path,shell=True).wait()


fa=open('./LST.txt')
PROC_LIMIT=6
jobs=[]
i=1

for line in fa:
    print i;i+=1
    seq=line.rstrip().split('\t')
    input_path=seq[0]
    print input_path
    if 1==1:
        p=multiprocessing.Process(target=Work, args=(input_path,1))
        p.start()
        jobs.append(p)
        if len(jobs)>=PROC_LIMIT:
            for p in jobs:
                p.join()
            jobs=[]
for p in jobs:
    p.join()



