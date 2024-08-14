library(Seurat)
library(Signac)

library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
annotations=readRDS('/home/database/annotation/hg38/hg38_signac_ucsc_annotations.rds')


ReadPeak<-function(PATH){
    counts <- Matrix::readMM(paste0(PATH, "/matrix.mtx"))
    barcodes <- readLines(paste0(PATH, "/barcodes.tsv"))
    peaks <- read.table(paste0(PATH,"/peaks.bed"), sep="\t")
    peaknames <- paste(peaks$V1, peaks$V2, peaks$V3, sep="-")
    colnames(counts) <- barcodes
    rownames(counts) <- peaknames
    return(counts)
    }


counts_d53=ReadPeak('../data/GSE184386/d53_fetal_filtered_peak_bc_matrix')
counts_d59=ReadPeak('../data/GSE184386/d59_fetal_filtered_peak_bc_matrix')
counts_d78=ReadPeak('../data/GSE184386/d78_fetal_filtered_peak_bc_matrix')
counts_d78_2=ReadPeak('../data/GSE184386/d78_fetal_2_filtered_peak_bc_matrix')
counts_d89C=ReadPeak('../data/GSE184386/d89C_fetal_filtered_peak_bc_matrix')
counts_d89P=ReadPeak('../data/GSE184386/d89P_fetal_filtered_peak_bc_matrix')

count2seurat<-function(counts, frag){
    chrom_assay <- CreateChromatinAssay(
        counts = counts,
        sep = c(":", "-"),
        fragments = frag,
        min.cells = 10,
        min.features = 200
        )
    pbmc <- CreateSeuratObject(
        counts = chrom_assay,
        assay = "peaks"
       )
    Annotation(pbmc) <- annotations
    return(pbmc)
    }

seurat_d53=count2seurat(counts_d53, '../data/GSE184386/hg38/GSM5585602_d53_fetal_fragments.tsv.gz')
seurat_d59=count2seurat(counts_d59, '../data/GSE184386/hg38/GSM5585601_d59_fetal_fragments.tsv.gz')
seurat_d78=count2seurat(counts_d78, '../data/GSE184386/hg38/GSM5585603_d78_fetal_fragments.tsv.gz')
seurat_d78_2=count2seurat(counts_d78_2, '../data/GSE184386/hg38/GSM5585604_d78_fetal_2_fragments.tsv.gz')
seurat_d89C=count2seurat(counts_d89C, '../data/GSE184386/hg38/GSM5585611_d89C_fetal_fragments.tsv.gz')
seurat_d89P=count2seurat(counts_d89P, '../data/GSE184386/hg38/GSM5585612_d89P_fetal_fragments.tsv.gz')



saveRDS(seurat_d53, '../rds/f1_d53/pbmc_count.rds')
saveRDS(seurat_d59, '../rds/f2_d59/pbmc_count.rds')
saveRDS(seurat_d78, '../rds/f3_d78/pbmc_count.rds')
saveRDS(seurat_d78_2, '../rds/f4_d78_2/pbmc_count.rds')
saveRDS(seurat_d89C, '../rds/f5_d89C/pbmc_count.rds')
saveRDS(seurat_d89P, '../rds/f6_d89P/pbmc_count.rds')























