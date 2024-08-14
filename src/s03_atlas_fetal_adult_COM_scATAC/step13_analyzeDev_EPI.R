
library(dplyr)
library(tidyr)
library(readr)
library(GenomicRanges)
library(reshape2)
#library(openxlsx)
library(ggplot2)
library(Matrix)
library(EpiTrace)
library(Seurat)
library(SeuratObject)
library(ggtree)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(ArchR)
library(parallel)


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
    pbmc=readRDS(paste0('/home/database/data/COM_scATAC/rds/',this_organ,'.seuratWithUmap.rds'))
    TIME=rep(0,length(pbmc$time))
    TIME[which(pbmc$time=='fetal')]=1
    ####################
    atac_hg38=pbmc@assays$macs2@counts
    human_clock_hg19 <- plyranges::reduce_ranges(c(clock_gr_list[[1]],clock_gr_list[[2]]))
    #easyLift::easyLiftOver(plyranges::reduce_ranges(c(clock_gr_list[[1]],clock_gr_list[[2]])),map = '/home/database/reference/hg19/hg19ToHg38.over.chain') -> human_clock_hg38
    this_dat_ranges <- data.frame(peaks=rownames(atac_hg38)) %>% separate(col=1,into=c('chr','start','end'),remove=F,convert=T) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    ################
    if(ncol(atac_hg38)<=100000){
        this_index=c(1:ncol(atac_hg38))
    }else{
        set.seed(123)
        this_index=sort(sample(c(1:ncol(atac_hg38)),100000))
    }
    ################
    EpiTrace::Init_Peakset(this_dat_ranges) -> init_gr
    EpiTrace::Init_Matrix(peakname = this_dat_ranges$peaks,cellname = colnames(atac_hg38)[this_index],matrix = atac_hg38[,this_index]) -> init_mm
    CNAME=colnames(pbmc)
    rm(pbmc)
    rm(atac_hg38)
    gc()
    EpiTrace::EpiTraceAge_Convergence(peakSet = init_gr,matrix=init_mm,ref_genome='hg38',clock_gr=human_clock_hg19, iterative_time = 5,min.cutoff = 0,non_standard_clock = T,qualnum = 10,ncore_lim = 3,mean_error_limit = 0.1, sep_str=c("-","-"),parallel=T) -> epitrace_obj_age_estimated    

    #####################    
    print('EPI:')
    SS=rep(NA, length(TIME))
    names(SS)=CNAME[this_index]
    SS[names(epitrace_obj_age_estimated$EpiTraceAge_iterative)]=epitrace_obj_age_estimated$EpiTraceAge_iterative
    USED=which(!is.na(SS))
    print(cor(SS[USED],TIME[this_index][USED],method='spearman'))    
    OUT=list()
    OUT$epi=SS
    OUT$this_index=this_index    
    saveRDS(OUT, file=paste0('../rds/devOut.EPI.',this_organ,'.rds'))
    saveRDS(epitrace_obj_age_estimated, file=paste0('../rds/epitrace.obj.',this_organ,'.rds'))
    #########################
    print(i)
    i=i+1
   }

