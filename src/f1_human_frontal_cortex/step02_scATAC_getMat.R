
library(Signac)
library(GenomicRanges)

opath='/home/database/data/InferX/data/frontal_cortex/scATAC/GSM5589415_UMB4540_snATAC_frontal_cortex_rep2_fragments.hg19.bed.gz.format'

fpath=paste0(opath, '/sorted.bed.gz')
bpath=paste0('/home/database/reference/hg19/hg19.fa.size.10k_coolBin.bed')


BIN=read.table(bpath,sep='\t',header=F,row.names=NULL)

TMP=BIN
PEAK=GRanges(seqnames = TMP$V1,
            ranges = IRanges(start = TMP$V2+1,
                             end = TMP$V3,
                             names = paste0('bin',c(1:nrow(TMP)))
                             )
             )


fragments <- CreateFragmentObject(fpath,validate.fragments=FALSE)


MAT=FeatureMatrix(
      fragments,
      features=PEAK,
      cells = NULL,
      process_n = 2000,
      sep = c(":", "-"),
      verbose = TRUE
      )


saveRDS(MAT, paste0(opath,'/mat_10k.rds'))

