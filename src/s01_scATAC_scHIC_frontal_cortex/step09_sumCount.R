


COUNT.Astro=read.table('./rds/ciceroFrame_conns_uniq.bedpe.Astro.count',sep='\t',row.names=NULL,header=F)
COUNT.MG=read.table('./rds/ciceroFrame_conns_uniq.bedpe.MG.count',sep='\t',row.names=NULL,header=F)
COUNT.Neuron.Ex=read.table('./rds/ciceroFrame_conns_uniq.bedpe.Neuron.Ex.count',sep='\t',row.names=NULL,header=F)
COUNT.Neuron.In=read.table('./rds/ciceroFrame_conns_uniq.bedpe.Neuron.In.count',sep='\t',row.names=NULL,header=F)
COUNT.ODC=read.table('./rds/ciceroFrame_conns_uniq.bedpe.ODC.count',sep='\t',row.names=NULL,header=F)
COUNT.OPC=read.table('./rds/ciceroFrame_conns_uniq.bedpe.OPC.count',sep='\t',row.names=NULL,header=F)



COUNT_ALL=cbind(COUNT.Astro[,7],
                COUNT.MG[,7],
                COUNT.Neuron.Ex[,7],
                COUNT.Neuron.In[,7],
                COUNT.ODC[,7],
                COUNT.OPC[,7])

colnames(COUNT_ALL)=c('AS','MG','NE','NI','OD','OP')
rownames(COUNT_ALL)=paste0(COUNT.Astro[,1],'-',COUNT.Astro[,2],'-',COUNT.Astro[,3],'.And.',
                           COUNT.Astro[,4],'-',COUNT.Astro[,5],'-',COUNT.Astro[,6])


COUNT_ALL=cbind(COUNT_ALL, rowSums(COUNT_ALL))
colnames(COUNT_ALL)[7]='ALL'
saveRDS(COUNT_ALL,'./rds/COUNT_ALL.rds')


CONNS=c()

conns=readRDS('./rds/ciceroFrame_conns_D_uniq.rds')
CONNS=cbind(CONNS,conns[,3])
conns=readRDS('./rds/ciceroFrame_conns_cicero_uniq.rds')
CONNS=cbind(CONNS,conns[,3])
conns=readRDS('./rds/ciceroFrame_conns_cor_uniq.rds')
CONNS=cbind(CONNS,conns[,3])
conns=readRDS('./rds/ciceroFrame_conns_spearman_uniq.rds')
CONNS=cbind(CONNS,conns[,3])

rownames(CONNS)=stringr::str_replace_all(rownames(conns),'_','-')
colnames(CONNS)=c('D','cicero','cor','spearman')


saveRDS(CONNS, file='./rds/CONNS_ALL.rds')












