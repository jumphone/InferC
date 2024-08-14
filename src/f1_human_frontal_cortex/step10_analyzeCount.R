
source('InferC.R')


COUNT=readRDS('./rds/COUNT_ALL.rds')
CONNS=readRDS('./rds/CONNS_ALL.rds')





cor(COUNT, CONNS,method='spearman')


set.seed(123)
random_select=sample(c(1:nrow(COUNT)),10000)

CONNS_Q=apply(CONNS,2,rank,ties.method='random')/nrow(CONNS)

plot(CONNS_Q[random_select,1],COUNT[random_select,7])
points(CONNS_Q[random_select,2],COUNT[random_select,7],col='red')

cor(CONNS[random_select,1],COUNT[random_select,7])
cor(CONNS[random_select,2],COUNT[random_select,7])



plot(CONNS[random_select,1],CONNS[random_select,2])




COR=cor(COUNT, CONNS, method='spearman')

library('ComplexHeatmap')
library('circlize')

o.mat=COR
###############################
vvv=as.numeric(o.mat)
col_fun =colorRamp2(
             c(0.06,0.07, 0.08, 0.15),
             c('royalblue3','skyblue1','grey95','indianred1') 
             )

ht1=Heatmap(o.mat,row_title='',name="v",
        cluster_rows=T, cluster_columns=F,
        show_column_dend = T, show_row_dend = T,
        show_column_names= T, show_row_names= T,
        #row_split=R.SPLIT,
        #column_split=C.SPLIT,
        col=col_fun, border = TRUE,
        row_names_side='right',
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", o.mat[i, j]), x, y, gp = gpar(fontsize = 10))
            }
        )

ht2=Heatmap(o.mat,row_title='',name="v",
        cluster_rows=T, cluster_columns=F,
        show_column_dend = T, show_row_dend = T,
        show_column_names= T, show_row_names= T,
        #row_split=R.SPLIT,
        #column_split=C.SPLIT,
        col=col_fun, border = TRUE,
        row_names_side='right',
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        )


pdf('./plot/p01_compare.pdf',width=5,height=5)
print(ht1)
print(ht2)
dev.off()


pdf('./plot/p02_compare_boxplot.pdf',width=5,height=2.5)
boxplot(COR,col=c('indianred1',rep('grey60',3)),boxwex=0.5)
dev.off()


apply(COR,2,mean)
#   D     cicero        cor   spearman
#0.11993970 0.09284959 0.09376649 0.09394565


t.test(COR[,1],COR[,2],paired=T)
#0.005963
t.test(COR[,1],COR[,3],paired=T)
#0.0007307
t.test(COR[,1],COR[,4],paired=T)
#0.0007863






COUNT_ALL=COUNT[,7]
COUNT_ALL[which(COUNT_ALL>10)]=10

TAB_ALL=table(COUNT_ALL)
barplot(TAB_ALL)

print(TAB_ALL/nrow(COUNT))
#         0           1           2           3           4           5
#0.542664650 0.197982326 0.096960468 0.051645859 0.029955251 0.018780949
#          6           7           8           9          10
#0.012488899 0.008802709 0.006392849 0.004810387 0.029515653

pdf('./plot/p03_taball_barplot.pdf',width=6,height=5)
tmp=barplot(TAB_ALL,col=c('white',rep('grey60',10)),space=0.5)
segments(x0=(tmp[1]+tmp[2])/2,x1=(tmp[1]+tmp[2])/2,y0=-0.1,y1=max(TAB_ALL)*1.2,lty=2,lwd=2)
pie(c(TAB_ALL[1],sum(TAB_ALL[2:11])),col=c('white','grey60'),labels=NA)
dev.off()


sum(TAB_ALL[2:11])/nrow(COUNT)


########################################
AUC=c()
this_cut=1
max_cut=10
while(this_cut<=max_cut){
    print(paste0(this_cut,' / ',max_cut))
    TAG1=c()
    TAG2=c()
    this_auc_all=c()
    i=7
    #while(i<=ncol(COUNT)){
        this_resp=rep(0,nrow(COUNT))
        this_resp[which(COUNT[,i]>=this_cut)]=1
        j=1
        while(j<=ncol(CONNS)){
            print(colnames(COUNT)[i])
            print(colnames(CONNS)[j])
            TAG1=c(TAG1, colnames(COUNT)[i])
            TAG2=c(TAG2, colnames(CONNS)[j])
            this_conns=CONNS[,j]
            this_pred=this_conns
            this_roc=pROC::roc(response=this_resp, predictor=this_pred, quiet =TRUE,direction='<')
            this_auc=pROC::auc(this_roc)
            this_auc_all=c(this_auc_all, this_auc)
            j=j+1}
    #   i=i+1
    #   }
    AUC=cbind(AUC, this_auc_all)
    this_cut=this_cut+1
    }
####################################

saveRDS(AUC, file='./rds/AUC_AUC.rds')
saveRDS(TAG1, file='./rds/AUC_TAG1.rds')
saveRDS(TAG2, file='./rds/AUC_TAG2.rds')



AUC=readRDS('./rds/AUC_AUC.rds')
TAG1=readRDS('./rds/AUC_TAG1.rds')
TAG2=readRDS('./rds/AUC_TAG2.rds')

SHOW_INDEX=c(1:10)

pdf('./plot/p04_auc_plot.pdf',width=4.5,height=5)
################
CEX=1.5
plot(x=SHOW_INDEX, y=rep(0,length(SHOW_INDEX)),ylim=c(0.54,0.7),col='white',xlim=c(1,10))
points(x=SHOW_INDEX,AUC[which(TAG2=='D' & TAG1=='ALL'),SHOW_INDEX],col='indianred1',type='b',pch=16,cex=CEX)
points(x=SHOW_INDEX,AUC[which(TAG2=='cicero' & TAG1=='ALL'),SHOW_INDEX],col='black',type='b',pch=17,cex=CEX)
points(x=SHOW_INDEX,AUC[which(TAG2=='cor' & TAG1=='ALL'),SHOW_INDEX],col='black',type='b',pch=1,cex=CEX)
points(x=SHOW_INDEX,AUC[which(TAG2=='spearman' & TAG1=='ALL'),SHOW_INDEX],col='black',type='b',pch=2,cex=CEX)
#################
dev.off()


mean(AUC[which(TAG2=='D' & TAG1=='ALL'),SHOW_INDEX])
#[1] 0.650046

mean(AUC[which(TAG2=='cicero' & TAG1=='ALL'),SHOW_INDEX])
#[1] 0.6099617

mean(AUC[which(TAG2=='cor' & TAG1=='ALL'),SHOW_INDEX])
#[1] 0.6242422

mean(AUC[which(TAG2=='spearman' & TAG1=='ALL'),SHOW_INDEX])
#[1] 0.6235652


t.test(AUC[which(TAG2=='D' & TAG1=='ALL'),SHOW_INDEX],AUC[which(TAG2=='cicero' & TAG1=='ALL'),SHOW_INDEX],paired=T)
#p-value = 7.6e-09


















