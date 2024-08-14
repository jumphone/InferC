library(dplyr)
library(tidyr)
library(readr)
library(GenomicRanges)
library(reshape2)
library(openxlsx)
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
#library(ChIPseeker)
#library(BSgenome.Hsapiens.UCSC.hg19)
#library(BSgenome.Hsapiens.UCSC.hg38)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(org.Hs.eg.db)



###############################
library(Seurat)
library(Signac)
source('/home/toolkit/src/fitdevo.R')
source('../../frontal_cortex/scATAC/InferC.R')

META=readRDS('./rds/scmul_pbmc_usedType_meta.rds')
UMAP=readRDS('./rds/scmul_pbmc_usedType_umap.rds')
PCA=readRDS('./rds/scmul_pbmc_usedType_pca.rds')
LSI=readRDS('./rds/scmul_pbmc_usedType_lsi.rds')
DP=readRDS('./rds/scmul_pbmc_usedType_dp.rds')
AGE=readRDS('./rds/scmul_pbmc_usedType_age.rds')


pbmc_sub=readRDS('./rds/scmul_pbmc_usedType.rds')
TYPE=pbmc_sub$celltype
atac_hg38=pbmc_sub@assays$ATAC@counts
#########################################

human_clock_hg19 <- plyranges::reduce_ranges(c(clock_gr_list[[1]],clock_gr_list[[2]]))

this_dat_ranges_hg38 <- data.frame(peaks=rownames(atac_hg38)) %>% separate(col=1,into=c('chr','start','end'),remove=F,convert=T) %>% makeGRangesFromDataFrame(keep.extra.columns = T)


##########
EpiTrace::Init_Peakset(this_dat_ranges) -> init_gr
EpiTrace::Init_Matrix(peakname = this_dat_ranges$peaks,cellname = colnames(atac_hg38),matrix = atac_hg38) -> init_mm
############

############
EpiTrace::EpiTraceAge_Convergence(peakSet = init_gr,matrix=init_mm,ref_genome='hg38',clock_gr=human_clock_hg19, iterative_time = 5,min.cutoff = 0,non_standard_clock = T,qualnum = 10,ncore_lim = 12,mean_error_limit = 0.1, sep_str=c("-","-"),parallel=T) -> epitrace_obj_age_estimated
#good quality cells 39235 and peaks 17189
#############################

EPI=epitrace_obj_age_estimated$EpiTraceAge_iterative

saveRDS(epitrace_obj_age_estimated,file='./rds/epitrace_obj_age_estimated.rds')
saveRDS(EPI,file='./rds/epi_all.rds')

names(AGE)=rownames(META)

cor(AGE[names(EPI)],EPI,method='spearman')


EPI=epitrace_obj_age_estimated$EpiTraceAge_Clock_initial

cor(AGE[names(EPI)],EPI,method='spearman')















