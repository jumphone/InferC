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
    Monocle: 
    bigWigAverageOverBed: downloaded from UCSC Genome Browser
    
</br>
    
# Usage:

## GLS

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
    
    
    
</br>

# Benchmark datasets:

Original Data: **Supplementary Table 1**

Processed Data: [BaiduNetdisk](https://pan.baidu.com/s/11o2fzAkUmNuIKMg0sk8KYg?pwd=34q2)

</br>

