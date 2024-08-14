
source('InferC.R')


COUNT=readRDS('./rds/COUNT_ALL.rds')
CONNS=readRDS('./rds/CONNS_ALL.rds')


TMP=read.table('./scHIC/Merged_loop_within500k_overPeakPair_allinfo.bedpe',sep='\t',row.names=NULL, header=F)
LOOP=paste0(TMP[,7],'-',TMP[,8],'-',TMP[,9],'.And.',TMP[,10],'-',TMP[,11],'-',TMP[,12])


LOOP_INDEX=which(rownames(CONNS) %in% LOOP)
nonLOOP_INDEX=which(!rownames(CONNS) %in% LOOP)
set.seed(123)
nonLOOP_INDEX_RANDOM=sample(nonLOOP_INDEX,length(LOOP_INDEX))


t.test(CONNS[LOOP_INDEX,1],CONNS[nonLOOP_INDEX_RANDOM,1])
#2.2e-16
t.test(CONNS[LOOP_INDEX,2],CONNS[nonLOOP_INDEX_RANDOM,2])
#2.2e-16
t.test(CONNS[LOOP_INDEX,3],CONNS[nonLOOP_INDEX_RANDOM,3])
#2.2e-16
t.test(CONNS[LOOP_INDEX,4],CONNS[nonLOOP_INDEX_RANDOM,4])
#2.2e-16

pdf('./plot/p05_gls_loop_boxplot.pdf',width=4.5,height=5)
boxplot(
        CONNS[LOOP_INDEX,1],CONNS[nonLOOP_INDEX_RANDOM,1],
        CONNS[LOOP_INDEX,2],CONNS[nonLOOP_INDEX_RANDOM,2],
        CONNS[LOOP_INDEX,3],CONNS[nonLOOP_INDEX_RANDOM,3],
        CONNS[LOOP_INDEX,4],CONNS[nonLOOP_INDEX_RANDOM,4],
        col=c('indianred1','indianred1',rep('grey60',6)),
        outline=TRUE, cex=0.1, pch='+')
abline(v=2.5,lty=2)
abline(v=4.5,lty=2)
abline(v=6.5,lty=2)

dev.off()

AUC=c()

i=1
while(i<=ncol(CONNS)){
    this_resp=rep(0,nrow(CONNS))
    this_resp[LOOP_INDEX]=1

    this_pred=CONNS[,i]
    this_roc=pROC::roc(response=this_resp, predictor=this_pred, quiet =TRUE,direction='<')
    this_auc=pROC::auc(this_roc)
    AUC=c(AUC, this_auc)
    i=i+1
    }
names(AUC)=colnames(CONNS)

saveRDS(AUC, './scHIC/Merged_loop_within500k_overPeakPair_allinfo.bedpe.AUC.rds')

print(AUC)
#        D    cicero       cor  spearman
#0.6696789 0.6306404 0.6309917 0.6239903



this_pred=CONNS[,1]
this_roc=pROC::roc(response=this_resp, predictor=this_pred, quiet =TRUE,direction='<')
this_auc=pROC::auc(this_roc)

this_best=pROC::coords(this_roc, "best", ret = "threshold")

this_best
#0.282869


apply(CONNS[LOOP_INDEX,],2,median)

#      D    cicero       cor  spearman
#0.3127556 0.0000000 0.1212900 0.1224800

apply(CONNS[nonLOOP_INDEX,],2,median)

#         D     cicero        cor   spearman
#0.19707181 0.00000000 0.03910837 0.05016869

apply(CONNS[LOOP_INDEX,],2,median)-apply(CONNS[nonLOOP_INDEX,],2,median)

# D     cicero        cor   spearman
#0.11568379 0.00000000 0.08218159 0.07231133






