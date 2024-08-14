library(Seurat)
library(Signac)

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
    print(this_all_path_list)
    PBMC=c()
    j=1
    while(j<=length(this_all_path_list)){
        this_path=this_all_path_list[j]
        this_seurat_file=paste0(this_path,'/','seurat_merge_peak.rds')
        this_pbmc=readRDS(this_seurat_file)
        this_pbmc <- RunTFIDF(this_pbmc)
        this_pbmc$tag=rep(this_all_tag[j],ncol(this_pbmc))
        colnames(this_pbmc)=paste0(this_pbmc$tag,'_',j,'_',colnames(this_pbmc))
        PBMC=c(PBMC, this_pbmc)
        print(this_path)
        j=j+1}
    pbmc=merge(PBMC[[1]],y=PBMC[2:length(PBMC)])
    TMP=t(matrix(unlist(stringr::str_split(colnames(pbmc),'_')),nrow=3))
    pbmc$time=TMP[,1]
    #tapply(colSums(pbmc@assays$macs2@data),pbmc$time,mean) 
    saveRDS(pbmc, paste0('/home/database/data/COM_scATAC/data/',this_organ,'.seurat.rds'))
    print(i)
    i=i+1
    }

