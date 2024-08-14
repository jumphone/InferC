LIST_COM=read.table('COM_ID.txt')
LIST_ADULT=read.table('LST_ADULT_ID.txt',row.names=1)
LIST_FETAL=read.table('LST_FETAL_ID.txt',row.names=1)

options(scipen = 999999999)


.calauc <- function(this_pred,this_resp){
    this_roc=pROC::roc(response=this_resp, predictor=this_pred, quiet =TRUE,direction='<')
    this_auc=pROC::auc(this_roc)[1]
    return(this_auc)
    }

COR=matrix(0,ncol=7,nrow=nrow(LIST_COM))
AUC=matrix(0,ncol=7,nrow=nrow(LIST_COM))

rownames(COR)=LIST_COM[,1]
rownames(AUC)=LIST_COM[,1]
colnames(COR)=c('raw','inferc','esc','nfeature','ncount','epi','epi_init')
colnames(AUC)=c('raw','inferc','esc','nfeature','ncount','epi','epi_init')

i=1
while(i<=nrow(LIST_COM)){
    this_organ=LIST_COM[i,1]
    ###########################
    OUT=readRDS(paste0('../rds/devOut.',this_organ,'.rds'))
    OUT1=readRDS(paste0('../rds/devOut.GA.',this_organ,'.rds'))
    #############################    
    AUC[i,1]=.calauc(OUT$ss1,OUT$time)
    COR[i,1]=cor(OUT$ss1,OUT$time,method='spearman')
    AUC[i,2]=.calauc(OUT$ss2,OUT$time)
    COR[i,2]=cor(OUT$ss2,OUT$time,method='spearman')
    AUC[i,3]=.calauc(OUT1$es,OUT$time)
    COR[i,3]=cor(OUT1$es,OUT$time,method='spearman')
    AUC[i,4]=.calauc(OUT$ss3,OUT$time)
    COR[i,4]=cor(OUT$ss3,OUT$time,method='spearman')
    AUC[i,5]=.calauc(OUT$ss4,OUT$time)
    COR[i,5]=cor(OUT$ss4,OUT$time,method='spearman')
    #########################
    EPI=readRDS(paste0('../rds/epitrace.obj.',this_organ,'.rds'))
    #OUT2=readRDS(paste0('../rds/devOut.EPI.',this_organ,'.rds'))    
    S1=EPI$EpiTraceAge_iterative
    S2=EPI$EpiTraceAge_Clock_initial
    CELL=names(OUT$ss3)
    names(OUT$time)=CELL
    ###############################
    this_score=-S1
    this_used=which(!is.na(this_score))
    COR[i,6]=cor(this_score[this_used],OUT$time[names(this_score)][this_used],method='spearman')
    AUC[i,6]=.calauc(this_score[this_used],OUT$time[names(this_score)][this_used])
    #################################
    this_score=-S2
    this_used=which(!is.na(this_score))
    COR[i,7]=cor(this_score[this_used],OUT$time[names(this_score)][this_used],method='spearman')
    AUC[i,7]=.calauc(this_score[this_used],OUT$time[names(this_score)][this_used])
    print(i)
    i=i+1
   }

print(AUC)
print(COR)

#print(mean(AUC[,2]))
print(mean(COR[,2]))
print( t.test(AUC[,2],AUC[,1],paired=TRUE) )
print( t.test(COR[,2],COR[,1],paired=TRUE) )

saveRDS(AUC,'../rds/Eva_AUC_EPI.rds')
saveRDS(COR,'../rds/Eva_COR_EPI.rds')




###########
library('ComplexHeatmap')
library('circlize')


MAT=COR[,c(1,3,6,5)]
MAT=MAT[order(-MAT[,1]),]

o.mat=MAT
col_fun=colorRamp2(
            c(-0.2, 0, 0.2,0.4),
             c('skyblue1',
              'white','indianred1','red3'))
ht1=Heatmap(o.mat,row_title='',name="v",
        cluster_rows=F, cluster_columns=F,
        #show_column_dend = T, show_row_dend = T,
        show_column_names= T, show_row_names= T,
        col=col_fun, border = TRUE,
        row_names_side='right',
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.3f", o.mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })



ht2=Heatmap(o.mat,row_title='',name="v",
        cluster_rows=F, cluster_columns=F,
        #show_column_dend = T, show_row_dend = T,
        show_column_names= T, show_row_names= T,
        col=col_fun, border = TRUE,
        row_names_side='right',
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        )

pdf('../plot/s15p01_heatmap_epitrace.pdf',height=8,width=5)
ht1
ht2
dev.off()

apply(o.mat,2,mean)
#        raw         esc         epi      ncount
# 0.34695804  0.13992435 -0.19132428 -0.01011802










