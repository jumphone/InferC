MAX=500000
MIN=1000
fi=open('scHIC/Merged_loop.bedpe')
fo=open('scHIC/Merged_loop_within500k.bedpe','w')

for line in fi:
    seq=line.rstrip().split('\t')
    l1=int(seq[1])
    l2=int(seq[2])
    l3=int(seq[4])
    l4=int(seq[5])
    this_dist=max([l1,l2,l3,l4])-min([l1,l2,l3,l4])
    if this_dist > MIN and this_dist < MAX:
        fo.write(line)

fi.close()
fo.close()

