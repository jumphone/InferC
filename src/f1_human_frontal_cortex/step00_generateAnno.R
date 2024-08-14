

library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

annotations@seqnames@values=paste0('chr',annotations@seqnames@values)
annotations@seqinfo@seqnames=paste0('chr',annotations@seqinfo@seqnames)
annotations@seqinfo@genome=rep('hg19',length(annotations@seqinfo@genome))
annotations@seqnames@values[which(annotations@seqnames@values=='chrMT')]='chrM'
annotations@seqinfo@seqnames[which(annotations@seqinfo@seqnames=='chrMT')]='chrM'
annotations@seqnames@values=factor(annotations@seqnames@values,levels=annotations@seqinfo@seqnames)

saveRDS(annotations, '/home/database/annotation/hg19/hg19_signac_ucsc_annotations.rds')



