LIST_COM=read.table('COM_ID.txt')
LIST_ADULT=read.table('LST_ADULT_ID.txt',row.names=1)
LIST_FETAL=read.table('LST_FETAL_ID.txt',row.names=1)

options(scipen = 999999999)


.calauc <- function(this_pred,this_resp){
    this_roc=pROC::roc(response=this_resp, predictor=this_pred, quiet =TRUE,direction='<')
    this_auc=pROC::auc(this_roc)[1]
    return(this_auc)
    }

COR=matrix(0,ncol=5,nrow=nrow(LIST_COM))
AUC=matrix(0,ncol=5,nrow=nrow(LIST_COM))

rownames(COR)=LIST_COM[,1]
rownames(AUC)=LIST_COM[,1]
colnames(COR)=c('raw','inferc','esc','nfeature','ncount')
colnames(AUC)=c('raw','inferc','esc','nfeature','ncount')

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
    print(i)
    i=i+1
   }

print(AUC)
print(COR)

#print(mean(AUC[,2]))
print(mean(COR[,2]))
print( t.test(AUC[,2],AUC[,1],paired=TRUE) )
print( t.test(COR[,2],COR[,1],paired=TRUE) )

saveRDS(AUC,'../rds/Eva_AUC.rds')
saveRDS(COR,'../rds/Eva_COR.rds')


############################


LIST_COM=read.table('COM_ID.txt')
LIST_ADULT=read.table('LST_ADULT_ID.txt',row.names=1)
LIST_FETAL=read.table('LST_FETAL_ID.txt',row.names=1)

options(scipen = 999999999)

i=1
while(i<=nrow(LIST_COM)){
    this_organ=LIST_COM[i,1]
    ###########################
    this_adult=stringr::str_split(LIST_COM[i,2],',')[[1]]
    this_fetal=stringr::str_split(LIST_COM[i,3],',')[[1]]
    this_adult_path_list=LIST_ADULT[this_adult,1]
    this_fetal_path_list=LIST_FETAL[this_fetal,1]
    this_all_path_list=c(this_adult_path_list,this_fetal_path_list)
    this_all_tag=c(rep('adult',length(this_adult_path_list)),rep('fetal',length(this_fetal_path_list)))
    #print(this_all_path_list)
    print(this_organ)
    ###############
    pbmc=readRDS(paste0('/home/database/data/COM_scATAC/rds/',this_organ,'.seuratWithUmap.rds'))
    TIME=rep(0,length(pbmc$time))
    TIME[which(pbmc$time=='fetal')]=1
    print(table(TIME))
    i=i+1
    } 


#"Muscle"
#TIME
#    0     1
#60391 52886
#[1] "Adrenal"
#TIME
#     0      1
# 10537 101435
#[1] "Lung"
#TIME
#     0      1
# 41089 145131
#[1] "Liver"
#TIME
#     0      1
# 10557 293375
#[1] "Stomach"
#TIME
#    0     1
#29255  5797
#[1] "Pancreas"
#TIME
#    0     1
#43493  7144













