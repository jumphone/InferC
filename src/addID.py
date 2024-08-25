import sys
fi=open(sys.argv[1])
fo=open(sys.argv[2],'w')
i=1
for line in fi:
    seq = line.rstrip().split('\t')
    out = [seq[0],seq[1],seq[2]]+[str(i)]
    i=i+1
    fo.write('\t'.join(out)+'\n')

fo.close()
fi.close()
