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
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(ArchR)
library(parallel)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)



library(Seurat)
library(Signac)

library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
annotations=readRDS('/home/database/annotation/hg38/hg38_signac_ucsc_annotations.rds')

source('../../../frontal_cortex/scATAC/InferC.R')

UMAP=readRDS('../rds/pbmc_umap.rds')
BATCH=readRDS('../rds/BATCH.rds')
TIME=readRDS('../rds/TIME.rds')

pbmc=readRDS('../rds/pbmc.rds')

atac_hg38=pbmc@assays$RNA@counts

human_clock_hg19 <- plyranges::reduce_ranges(c(clock_gr_list[[1]],clock_gr_list[[2]]))

#easyLift::easyLiftOver(plyranges::reduce_ranges(c(clock_gr_list[[1]],clock_gr_list[[2]])),map = '/home/database/reference/hg19/hg19ToHg38.over.chain') -> human_clock_hg38

this_dat_ranges <- data.frame(peaks=rownames(atac_hg38)) %>% separate(col=1,into=c('chr','start','end'),remove=F,convert=T) %>% makeGRangesFromDataFrame(keep.extra.columns = T)

##########
EpiTrace::Init_Peakset(this_dat_ranges) -> init_gr
EpiTrace::Init_Matrix(peakname = this_dat_ranges$peaks,cellname = colnames(atac_hg38),matrix = atac_hg38) -> init_mm
############

############
EpiTrace::EpiTraceAge_Convergence(peakSet = init_gr,matrix=init_mm,ref_genome='hg38',clock_gr=human_clock_hg19, iterative_time = 5,min.cutoff = 0,non_standard_clock = T,qualnum = 10,ncore_lim = 3,mean_error_limit = 0.1, sep_str=c("-","-"),parallel=T) -> epitrace_obj_age_estimated
#############################
#good quality cells 32849 and peaks 22506

EPI=rep(NA,length(TIME))
names(EPI)=colnames(pbmc)
EPI[names(epitrace_obj_age_estimated$EpiTraceAge_iterative)]=epitrace_obj_age_estimated$EpiTraceAge_iterative

saveRDS(epitrace_obj_age_estimated,file='../rds/epitrace_obj_age_estimated.rds')
saveRDS(EPI,file='../rds/epi_all.rds')

#EPI is predicted age


EPI=-EPI
#negative EPI is developmental potential

cor(EPI, TIME, method='spearman')
# 0.09992432
# Stage is max(TIME)-TIME



pdf('../plot/s9p0x_time_epitrace.pdf',width=4,height=4)
boxplot(EPI~TIME,pch='+',cex=0.5,col='grey60',main='nFeature')
dev.off()

tapply(EPI,TIME,mean)
#1          2          3          4
#-0.5976259 -0.4344369 -0.4949837 -0.4819484









