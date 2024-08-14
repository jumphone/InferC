library(Seurat)
library(Signac)
source('/home/toolkit/src/fitdevo.R')
source('../../frontal_cortex/scATAC/InferC.R')

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
    ###############
    pbmc=readRDS(paste0('/home/database/data/COM_scATAC/data/',this_organ,'.seurat.rds'))
    DefaultAssay(pbmc)='macs2'
    pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
    pbmc <- RunSVD(pbmc)
    pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
    pdf(paste0('/home/database/data/COM_scATAC/plot/step10_',this_organ,'.seurat.umap.pdf'),width=10,height=10)
    DimPlot(pbmc, group.by='time',label=T,raster=TRUE)+NoLegend()
    dev.off()
    saveRDS(pbmc, paste0('/home/database/data/COM_scATAC/rds/',this_organ,'.seuratWithUmap.rds'))
    ####################################################
    atac_mat=pbmc[['macs2']]@data
    UMAP=pbmc@reductions$umap@cell.embeddings
    set.seed(123)
    random_index=sample(1:ncol(atac_mat),5000)
    atac_mat_sub=atac_mat[,random_index]
    UMAP_sub=UMAP[random_index,]
    used_coords=UMAP_sub
    indata=atac_mat_sub
    genome.df=inferloop.getGenomeDF.hg38()
    window=500000
    sample_num=100
    set.seed(123)
    cicero_cds=inferc.ciceroFrame_step01_prepareInput(indata, used_coords=used_coords,genome.df=genome.df,window=window,sample_num=sample_num)
    saveRDS(cicero_cds,paste0('../rds/ciceroFrame_cicero_cds_',this_organ,'.rds'))
    used_function=inferc.rowCalD
    conns <- inferc.ciceroFrame_step02_runUsedFuntion(cicero_cds, genome.df, window, sample_num, used_function)
    saveRDS(conns,paste0('../rds/ciceroFrame_conns_D_',this_organ,'.rds'))
    #############################
    BATCH=as.character(pbmc$orig.ident)
    AGG=inferc.aggMat(atac_mat, UMAP, clstNum=1000,seed=123, BATCH=BATCH)
    saveRDS(AGG, file=paste0('../rds/scmul_pbmc_usedType_atacMat_agg_',this_organ,'.rds'))
    ####################################################
    print(i)
    i=i+1
    }

