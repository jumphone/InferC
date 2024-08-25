<img src="https://fzhang.bioinfo-lab.com/img/INFERC_LOGO.png" width="250">

**InferC: inferring the chromatin structure and developmental potential of cells using single-cell chromatin accessibility**

**Paper Link:** coming soon

This tool is designed for inferring the chromatin interaction strength and developmental potential of cells using scATAC-seq data

# Updates:

**2024.08.25, v1.0.0 - Paper version.** details about this version is described in our paper.

</br>


# Content:

* [Requirements](#requirements)
* [Usage](#usage)
* [Benchmark Datasets](#benchmark-datasets)

</br>

# Requirements:

    R: 4.2.0
    Seurat: 4.3.0
    Signac: 1.9.0
    Cicero: 1.16.1
    Monocle: 2.24.0
    bigWigAverageOverBed: downloaded from the UCSC Genome Browser
    
</br>
    
# Usage:

## Calculate GLS：

    # Input: peak-count matrix

    library(Seurat)
    library(Signac)
    library(Signac)
    source('InferC.R')
    
    atac_count=readRDS('peak_count_matrix.rds')
    chrom_assay <- CreateChromatinAssay(
        counts = atac_count,
        sep = c(":", "-"),
        min.cells = 10,
        min.features = 200
        )
    seurat_object <- CreateSeuratObject(
        counts = chrom_assay,
        assay = "peaks"
       )
    seurat_object <- RunTFIDF(seurat_object) 
    seurat_object <- FindTopFeatures(seurat_object, min.cutoff = 'q0')
    seurat_object <- RunSVD(seurat_object,n = 30)
    seurat_object <- RunUMAP(object = seurat_object, reduction = 'lsi', dims = 2:30 )

    DATA=seurat_object[['peaks']]@data
    UMAP=seurat_object@reductions$umap@cell.embeddings

    genome.df=inferloop.getGenomeDF.hg19()
    window=500000
    sample_num=100

    set.seed(123)
    indata=DATA
    used_coords=UMAP
    cicero_cds=inferc.ciceroFrame_step01_prepareInput(indata, used_coords=used_coords,genome.df=genome.df,window=window,sample_num=sample_num)
    used_function=inferc.rowCalD
    conns <- inferc.ciceroFrame_step02_runUsedFuntion(cicero_cds, genome.df, window, sample_num, used_function)
    
    head(conns)
    
    #"conns" stores the GLS of all peak-pairs.

## Calculate LLS：   

    AGG=inferc.aggMat(DATA, UMAP, clstNum=1000,seed=123)
    conns_uniq=inferloop.getUniqLoop(conns)
    
    NET=apply(conns_uniq[,c(1,2)],2,stringr::str_replace_all,'_','-')
    PEAK_USED=c(NET[,1],NET[,2])
    MAT=AGG$agg[which(rownames(AGG$agg) %in% PEAK_USED),]

    infercOut = inferc.calLoopScore(MAT, NET)
    LLS=infercOut$lls

    LLS[1:5,1:5]
    
## Draw contact map, loops, and virtual 4C:
    
    CHR='chr4'
    VIEW_POINT=55095264 # hg19, PDGFRA
    SLOP=2000000
    START=VIEW_POINT-SLOP
    END=VIEW_POINT+SLOP
    
    CELLTYPE=seurat_object$celltype # prepare the cell-type information before visualizing.
    
    s1_agg = inferc.aggMat(MAT, UMAP, clstNum=100, seed=123)

    s2_loop = inferc.calCorLoop(s1_agg$agg, CHR, START, END,
                            Dcut=0.5, SPLIT1='-', SPLIT2='-')

    LOOP = inferc.agg2cell(AGG=s2_loop$loopScore$lls, CLST=s1_agg$clst)
    LOOP_TYPE = inferc.generate_mean(LOOP, CELLTYPE)

    i=1 # draw the figure for CELLTYPE 1
    LOOP_SIGNAL=LOOP_TYPE[,i]
    
    s3_draw=inferc.preContact(LOOP_SIGNAL,exMin=0,MAIN=colnames(LOOP_TYPE)[i],Zmax=2.5)
    s3_draw=inferc.preV4C(s3_draw, VIEW_POINT)
    inferc.drawContact(s3_draw)
    inferc.addV4C(s3_draw)

## Calculate developmental score:

#### Step 1, Linux shell:
    
    PHYLOP='hg19.100way.phyloP100way.bw' # downloaded from the UCSC Genome Browser
    BED='ALL_PEAK.bed'
    bwAve='./src/bigWigAverageOverBed'
    addID='./src/addID.py'
   
    python $addID $BED $BED\.withID.bed
    $bwAve $PHYLOP $BED\.withID.bed $BED\.hg19.phyloP100way.txt

#### Step 2, R:

    library(Seurat)
    library(Signac)
    library(Signac)
    source('InferC.R')
    
    DATA=seurat_object[['peaks']]@data
    UMAP=seurat_object@reductions$umap@cell.embeddings
    conns=readRDS('conns.rds') #from "GLS"
    PHYLOP=read.table('ALL_PEAK.bed.hg19.phyloP100way.txt',sep='\t')
    AGG=inferc.aggMat(DATA, UMAP, clstNum=1000,seed=123)
    MAT=AGG$agg
    
    devScore=inferc.calDevScore(conns, MAT, PHYLOP)
    devScore_cell=devScore[AGG$clst]
    

## Draw developmental directions on UMAP:

   DS=devScore_cell
   COL=inferc.colMap(DS,c(min(DS),median(DS),max(DS)),c('blue','grey90','red'))
   FIELD=inferc.field(DP=DS, VEC=UMAP, SHOW=TRUE,N=11,P=0.9,AL=0.5,CUT=0,CEX=0.3,COL=COL)

</br>

# Benchmark datasets:

Original Data: **Supplementary Table 1**

Processed Data: [BaiduNetdisk](https://pan.baidu.com/s/11o2fzAkUmNuIKMg0sk8KYg?pwd=34q2)

</br>

