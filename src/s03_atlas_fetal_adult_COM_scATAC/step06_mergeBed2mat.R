library(Seurat)
library(Signac)

LIST_COM=read.table('COM_ID.txt')
LIST_ADULT=read.table('LST_ADULT_ID.txt',row.names=1)
LIST_FETAL=read.table('LST_FETAL_ID.txt',row.names=1)

options(scipen = 999999999)

i=1
while(i<=nrow(LIST_COM)){
    this_organ=LIST_COM[i,1]
    this_merge_bed_file=paste0('/home/database/data/COM_scATAC/data/',this_organ,'.bed')
    BED=read.table(this_merge_bed_file,sep='\t')
    this_merge_peak=GenomicRanges::GRanges(as.character(BED[,1]), ranges=IRanges::IRanges(as.numeric(BED[,2]),as.numeric(BED[,3])))
    ##########################################################
    this_adult=stringr::str_split(LIST_COM[i,2],',')[[1]]
    this_fetal=stringr::str_split(LIST_COM[i,3],',')[[1]]
    this_adult_path_list=LIST_ADULT[this_adult,1]
    this_fetal_path_list=LIST_FETAL[this_fetal,1]
    this_all_path_list=c(this_adult_path_list,this_fetal_path_list)
    print(this_all_path_list)
    j=1
    while(j<=length(this_all_path_list)){
        this_path=this_all_path_list[j]
        this_seurat_file=paste0(this_path,'/','seurat_object.rds')
        pbmc=readRDS(this_seurat_file)
        macs2_counts <- FeatureMatrix(fragments = Fragments(pbmc),features = this_merge_peak,cells = colnames(pbmc))
        saveRDS(macs2_counts, file=paste0(this_path,'/','macs2_counts_merge_peak.rds'))
        print(this_path)
        j=j+1}
        print(i)
    i=i+1
    }

