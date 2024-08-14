import subprocess
import multiprocessing
import os


def Work(input_paths,a):
    print input_paths
    subprocess.Popen('/home/toolkit/tools/R-4.2.0/bin/Rscript  ./RSCRIPT/scATAC_transform.R '+input_paths[0]+'.format/'+' '+input_paths[1],shell=True).wait()

fa=open('./LST.txt')
PROC_LIMIT=5
jobs=[]
i=1

for line in fa:
    print i;i+=1
    seq=line.rstrip().split('\t')
    input_paths=[seq[0],seq[1]]
    print input_paths
    if 1==1:
        p=multiprocessing.Process(target=Work, args=(input_paths,1))
        p.start()
        jobs.append(p)
        if len(jobs)>=PROC_LIMIT:
            for p in jobs:
                p.join()
            jobs=[]
for p in jobs:
    p.join()

