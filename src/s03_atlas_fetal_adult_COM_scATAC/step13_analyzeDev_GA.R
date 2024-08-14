library(Seurat)
library(Signac)
source('/home/toolkit/src/fitdevo.R')
source('../../frontal_cortex/scATAC/InferC.R')


LIST_COM=read.table('COM_ID.txt')
LIST_ADULT=read.table('LST_ADULT_ID.txt',row.names=1)
LIST_FETAL=read.table('LST_FETAL_ID.txt',row.names=1)

options(scipen = 999999999)

i=1
#i=4
while(i<=nrow(LIST_COM)){
    this_organ=LIST_COM[i,1]
    ###########################
    this_adult=stringr::str_split(LIST_COM[i,2],',')[[1]]
    this_fetal=stringr::str_split(LIST_COM[i,3],',')[[1]]
    this_adult_path_list=LIST_ADULT[this_adult,1]
    this_fetal_path_list=LIST_FETAL[this_fetal,1]
    this_all_path_list=c(this_adult_path_list,this_fetal_path_list)
    this_all_tag=c(rep('adult',length(this_adult_path_list)),rep('fetal',length(this_fetal_path_list)))
    print(this_all_path_list)
    ###############
    pbmc=readRDS(paste0('/home/database/data/COM_scATAC/data/',this_organ,'.seurat.GeneActivity.rds'))
    TIME=rep(0,length(pbmc$time))
    TIME[which(pbmc$time=='fetal')]=1
    #####################    
    GA=pbmc[['RNA']]@data
    ES=read.table('./GS/BHATTACHARYA_EMBRYONIC_STEM_CELL.gmt',sep='\t')[,1]
    GA_ES=GA[which(rownames(GA) %in% ES),]
    SS=apply(GA_ES,2,mean)
    print('ES:')
    print(cor(SS,TIME,method='spearman'))    
    OUT=list()
    OUT$es=SS    
    saveRDS(OUT, file=paste0('../rds/devOut.GA.',this_organ,'.rds'))
    #########################
    print(i)
    i=i+1
   }

