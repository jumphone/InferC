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
    print(this_all_path_list)
    j=1
    while(j<=length(this_all_path_list)){
        this_path=this_all_path_list[j]
        this_seurat_file=paste0(this_path,'/','seurat_object.rds')
        pbmc=readRDS(this_seurat_file)
        macs2_counts <- readRDS(paste0(this_path,'/','macs2_counts_merge_peak.rds'))
        chrom_assay <- CreateChromatinAssay(
            counts = macs2_counts,
            sep = c(":", "-"),
            genome = 'hg38',
            fragments = Fragments(pbmc)[[1]]@path,
            min.cells = 0,
            min.features = 200
            )
        this_tmp=stringr::str_split(this_path,'/')[[1]]
        pbmc_macs2 <- CreateSeuratObject(
            project=this_tmp[length(this_tmp)],
            counts = macs2_counts,
            assay = "macs2",
            )

        saveRDS(pbmc_macs2, file=paste0(this_path,'/','seurat_merge_peak.rds'))
        print(this_path)
        j=j+1}
        print(i)
    i=i+1
    }

