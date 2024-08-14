fi=open('../data/GSE216323/GSM6668798_scATAC_metadata.txt')
fo=open('../data/GSE216323/GSM6668798_scATAC_metadata_new.txt','w')
for line in fi:
    fo.write(line.replace('#','_'))
fo.close()
fi.close()
