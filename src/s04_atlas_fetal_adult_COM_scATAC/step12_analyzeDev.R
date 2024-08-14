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
    AGG=readRDS(paste0('../rds/scmul_pbmc_usedType_atacMat_agg_',this_organ,'.rds'))
    TIME=rep(0,length(pbmc$time))
    TIME[which(pbmc$time=='fetal')]=1
    atac_mat=pbmc[['macs2']]@data
    UMAP=pbmc@reductions$umap@cell.embeddings
    #####################
    PHYLOP=read.table(paste0('/home/database/data/COM_scATAC/data/',this_organ,'.bed.hg38.phyloP100way.txt'),sep='\t')
    rownames(PHYLOP)=PHYLOP[,1]
    PHYLOP_=PHYLOP
    rownames(PHYLOP_)=stringr::str_replace_all(PHYLOP[,1],'-','_')
    PP_=rank(PHYLOP_[,6],ties.method='average')
    names(PP_)=rownames(PHYLOP_)

    #####################
    conns=readRDS(paste0('../rds/ciceroFrame_conns_D_',this_organ,'.rds'))
    conns[,1]=as.character(conns[,1])
    conns[,2]=as.character(conns[,2])    

    #################
    R1=rank(conns[,3],ties.method='average')
    R2=rank(PP_[conns[,1]]+PP_[conns[,2]],ties.method='average')
    RR=R1+R2
    this_cut=sort(RR,decreasing=T)[100000*2]
    used_index=which(RR>=this_cut)

    #################################
    conns_used=conns[used_index,]
    print(dim(conns_used))
    conns_uniq=inferloop.getUniqLoop(conns_used)

    #################
    NET=apply(conns_uniq[,c(1,2)],2,stringr::str_replace_all,'_','-')
    PEAK_USED=c(NET[,1],NET[,2])
    MAT=AGG$agg[which(rownames(AGG$agg) %in% PEAK_USED),]

    #################
    infercOut = inferc.calLoopScore(MAT, NET) 
    LLS=infercOut$lls
    GLS=infercOut$gls
    #################
    LOOP=inferloop.splitLoop(rownames(LLS))
    P1=PHYLOP[LOOP[,1],6]
    P2=PHYLOP[LOOP[,2],6]
    #################
    
    ##########################
    absLLS=abs(LLS)
    A=absLLS*P1*P2
    B=absLLS*(abs(P1)+abs(P2))
    C=colSums(A)/colSums(B)
    S1=C
    SS1=S1[AGG$clst]
    cor(SS1,TIME, method='spearman')

    ############################
    CUMAP=t(.generate_mean(t(UMAP),AGG$clst))
    CUMAP=CUMAP[order(as.numeric(rownames(CUMAP))),]
    DIST=as.matrix(dist(CUMAP))
    RMAT=apply(DIST,2,rank)
    #####
    N=nrow(CUMAP)
    Q=0.99
    S2=( S1 %*% ( Q**RMAT ) ) * (1-Q)/Q
    SS2=S2[AGG$clst]
    cor(SS2,TIME, method='spearman')

    ############################
    SS3=pbmc$nFeature_macs2
    cor(SS3,TIME, method='spearman')
    SS4=pbmc$nCount_macs2
    cor(SS4,TIME, method='spearman')
    #########################
    print('raw:')
    print(cor(SS1,TIME, method='spearman'))
    print('smooth:')
    print(cor(SS2,TIME, method='spearman'))
    print('nFeature:')
    print(cor(SS3,TIME, method='spearman'))
    print('nCount:')
    print(cor(SS4,TIME, method='spearman'))
    #############################
    OUT=list()
    OUT$time=TIME
    OUT$ss1=SS1
    OUT$ss2=SS2
    OUT$ss3=SS3
    OUT$ss4=SS4
    saveRDS(OUT, file=paste0('../rds/devOut.',this_organ,'.rds'))
    #########################
    print(i)
    i=i+1
   }

