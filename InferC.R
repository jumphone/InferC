
CORES=5

#####################################################
# InferC: infer contact map from scATAC-seq data
# Feng Zhang
# 2024.01
######################################################

library(stringr)
library(hash)


.norm1=function(x){
    if(min(x)!=max(x)){
        y=(x-min(x))/(max(x)-min(x))
       }else{
        y=rep(0,length(x))
       }
    return(y)
    }


.simple_combine <- function(exp_sc_mat1, exp_sc_mat2, FILL=FALSE){
    FILL=FILL
    exp_sc_mat=exp_sc_mat1
    exp_ref_mat=exp_sc_mat2
    ##############################################
    if(FILL==TRUE){
        gene1=rownames(exp_sc_mat)
        gene2=rownames(exp_ref_mat)
        gene12=gene2[which(!gene2 %in% gene1)]
        gene21=gene1[which(!gene1 %in% gene2)]
        exp_sc_mat_add=matrix(0,ncol=ncol(exp_sc_mat),nrow=length(gene12))
        rownames(exp_sc_mat_add)=gene12
        colnames(exp_sc_mat_add)=colnames(exp_sc_mat)
        exp_ref_mat_add=matrix(0,ncol=ncol(exp_ref_mat),nrow=length(gene21))
        rownames(exp_ref_mat_add)=gene21
        colnames(exp_ref_mat_add)=colnames(exp_ref_mat)
        exp_sc_mat=rbind(exp_sc_mat, exp_sc_mat_add)
        exp_ref_mat=rbind(exp_ref_mat, exp_ref_mat_add)
    }
    ############################################
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    OUT=list()
    OUT$combine=cbind(exp_sc_mat,exp_ref_mat)
    return(OUT)
    }








################################################
# InferLoop functions, old version: https://github.com/jumphone/InferLoop 

inferloop.bed2granges<-function(bed){
    library(GenomicRanges)
    colnames(bed)=c('seqname','start','end')
    bed=as.data.frame(bed)
    # split and convert per region
    res <- makeGRangesFromDataFrame(bed)
    return(res)
    }


inferloop.rmOut<-function(X){
    X=X
    Q3=quantile(X,0.75)
    Q1=quantile(X,0.25)
    RANGE=Q3-Q1
    UP=Q3+1.5*RANGE
    LW=Q1-1.5*RANGE
    OUT=X[which(X<UP)]
    OUT=OUT[which(OUT>LW)]
    return(OUT)
    }


inferloop.calABCD<-function(X, Y, X_base, Y_base){
    X=X
    Y=Y
    X_base=X_base
    Y_base=Y_base
    X_delta=X-X_base
    Y_delta=Y-Y_base
    A = sum(X_delta * Y_delta)
    B = sum(X_delta ** 2)
    C = sum(Y_delta ** 2)
    ############################
    D = A / sqrt( B * C )
    OUT=list()
    OUT[['A']]=A
    OUT[['B']]=B
    OUT[['C']]=C
    OUT[['D']]=D
    return(OUT)
    }


inferloop.calILS<-function(X, Y, r=0, only_pos=FALSE,xD=FALSE){
    X=X
    Y=Y
    r=r
    only_pos=only_pos
    xD=xD
    ############################
    # Ensure positive value
    if( min(X)<0 ){
        X = X - min(inferloop.rmOut(X))
        X[which(X<0)]=0
        }
    if( min(Y)<0 ){
        Y = Y - min(inferloop.rmOut(Y))
        Y[which(Y<0)]=0
        }
    #############################
    r=r
    N=length(X)
    ############################
    X_mean = mean(X)
    Y_mean = mean(Y)
    ###########################
    X_base = X_mean * r
    Y_base = Y_mean * r
    X_delta = X - X_base
    Y_delta = Y - Y_base
    ###########################
    ABCD = inferloop.calABCD(X, Y, X_base, Y_base)
    ###########################
    D = ABCD[['D']]
    D[which(is.na(D))]=0
    ############################
    D_plus = ( ABCD[['A']] + X_delta * Y_delta ) / sqrt( (ABCD[['B']] + X_delta**2) * (ABCD[['C']] + Y_delta**2) )
    ###########################
    M = D_plus - D
    S = (1-D**2)/(N-1)
    ###########################
    ILS = M / S
    ILS[which(is.na(ILS))]=0
    if(xD==TRUE){ILS = ILS * D}
    if(only_pos==TRUE){ILS[which(ILS<0)]=0 }
    return( ILS )
    }




inferloop.getUniqLoop <-function(net){
    #######################
    library(hash)
    tmp1=t(as.matrix(net[,c(1,2)]))
    print('sorting ends of each loop...')
    tmp2=apply(tmp1,2,sort)
    tag=apply(tmp2,2,paste0,collapse='.And.')
    utag=unique(tag)
    ####################
    print('hashing...')
    h=hash(keys=utag, values=rep(0,length(utag)))
    ###################
    print('getting unique loop...')
    flag=rep(0,nrow(net))
    i=1
    while(i<=nrow(net)){
        this_tag=tag[i]
        this_v=as.numeric(hash::values(x=h, keys=this_tag))
        if(this_v>0){flag[i]=1}
        .set(h, keys=this_tag, values=this_v+1)
        if(i %%50000==1){print(paste0(i,' / ',nrow(net)))}
        i=i+1}
    out=net[which(flag==0),]
    print(paste0('Number of unique loops: ',nrow(out)))
    print('finished!')
    return(out)
    }



inferloop.inferLoopSignal<-function(mat, net, r=0,sep='.And.',xD=FALSE){
    library(hash)
    r=r # default r=0
    sep=sep
    xD=xD
    ####################
    mat=as.matrix(mat)
    net=as.matrix(net)
    tag=rownames(mat)
    net=net[which(net[,1] %in% tag & net[,2] %in% tag),]
    #######################################
    print('hashing...')
    h=hash(keys=tag,values=rep(0,nrow(mat)))
    i=1
    while(i<=nrow(mat)){
        this_tag=tag[i]
        this_v=as.numeric(mat[i,])
        .set(h, keys=this_tag, values=this_v)
        i=i+1}
    #######################################
    print('calculating ILS...')
    out=matrix(0,ncol=ncol(mat),nrow=nrow(net))
    colnames(out)=colnames(mat)
    rownames(out)=paste0(net[,1],sep,net[,2])
    i=1
    while(i<=nrow(net)){
        this_tag1=net[i,1]
        this_tag2=net[i,2]
        #######################
        x=as.vector(hash::values(x=h, keys=this_tag1)[,1])
        y=as.vector(hash::values(x=h, keys=this_tag2)[,1])
        #######################
        z=inferloop.calILS(x,y,r=r,xD=xD)
        #######################
        out[i,]=z
        ########################
        if(i%%10000==1){print(paste0(i,' / ',nrow(net)))}
        i=i+1}
    print('finished!')
    return(out)
    }



inferloop.loadSignal <-function(input_path){
   library(data.table)
   HEADER=as.character(fread(input_path,header=FALSE, nrows=1))
   NET=fread(input_path,header=FALSE,sep='\t',select=c(1))$V1
   SIGNAL=fread(input_path,header=FALSE,sep='\t',select=c(2:(length(HEADER)+1)))
   SIGNAL=as.matrix(SIGNAL)
   rownames(SIGNAL)=NET
   colnames(SIGNAL)=HEADER
   return(SIGNAL)
   }



inferloop.loadSignalNoGap <-function(input_path){
   library(data.table)
   HEADER=as.character(fread(input_path,header=FALSE, nrows=1))
   HEADER=HEADER[2:length(HEADER)]
   NET=fread(input_path,header=FALSE,sep='\t',select=c(1), skip=1)$V1
   SIGNAL=fread(input_path,header=FALSE,sep='\t',select=c(2:(length(HEADER)+1)), skip=1)
   SIGNAL=as.matrix(SIGNAL)
   rownames(SIGNAL)=NET
   colnames(SIGNAL)=HEADER
   return(SIGNAL)
   }



inferloop.splitLoop <-function(loop,tag='.And.',n=2){
   library(stringr)
   tag=tag
   n=n
   pair=t(matrix(unlist(str_split(loop,tag)),nrow=n))
   return(pair)
   }



inferloop.writeNet <-function(conns, output_path,  cut=400000){
    conns=conns
    CUT=cut
    LOOP=conns[which(rank(-conns[,3], ties.method='random')<= CUT ),]
    LOOP[,1]=stringr::str_replace_all(LOOP[,1],'_','-')
    LOOP[,2]=stringr::str_replace_all(LOOP[,2],'_','-')
    write.table(LOOP,  file= output_path, row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
    }



#########################################################################################################
.generate_mean <- function(exp_sc_mat, TAG, print_step=50){
    print_step=print_step
    exp_sc_mat=exp_sc_mat
    TAG=TAG

    NewRef=matrix(0,ncol=length(unique(TAG)),nrow=nrow(exp_sc_mat))

    TAG=as.character(TAG)
    refnames=unique(TAG)
    total_num=length(refnames)
    outnames=c()
    i=1
    while(i<=length(refnames)){
        one=refnames[i]
        this_col=which(TAG==one)
        outnames=c(outnames,one)
        if(length(this_col) >1){
            #this_new_ref=apply(exp_sc_mat[,this_col],1,mean)
            this_new_ref=rowSums(exp_sc_mat[,this_col])/length(this_col)
            }else{
            this_new_ref = exp_sc_mat[,this_col]
            }
        NewRef[,i]=this_new_ref
        if(i%%print_step==1){print(paste0(i,' / ' ,total_num ))}
        i=i+1
        }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    return(NewRef)
    }


.generate_func <- function(exp_sc_mat, TAG, func=mean, print_step=50){
    print_step=print_step
    exp_sc_mat=exp_sc_mat
    TAG=TAG

    NewRef=matrix(0,ncol=length(unique(TAG)),nrow=nrow(exp_sc_mat))

    TAG=as.character(TAG)
    refnames=unique(TAG)
    total_num=length(refnames)
    outnames=c()
    i=1
    while(i<=length(refnames)){
        one=refnames[i]
        this_col=which(TAG==one)
        outnames=c(outnames,one)
        if(length(this_col) >1){
            #this_new_ref=apply(exp_sc_mat[,this_col],1,mean)
            this_new_ref=apply(exp_sc_mat[,this_col],1,func)
            #this_new_ref=rowSums(exp_sc_mat[,this_col])/length(this_col)
            }else{
            this_new_ref = exp_sc_mat[,this_col]
            }
        NewRef[,i]=this_new_ref
        if(i%%print_step==1){print(paste0(i,' / ' ,total_num ))}
        i=i+1
        }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    return(NewRef)
    }


.fastRowScale <- function (mat, center = TRUE, scale = TRUE, scale_max = 10){
    if (center) {
        rm <- matrixStats::rowMeans2(x = mat, na.rm = TRUE)
    }
    if (scale) {
        if (center) {
            rsd <- matrixStats::rowSds(mat, center = rm)
        }
        else {
            rsd <- sqrt(x = matrixStats::rowSums2(x = mat^2)/(ncol(x = mat) -
                1))
        }
    }
    if (center) {
        mat <- mat - rm
    }
    if (scale) {
        mat <- mat/rsd
    }
    if (scale_max != Inf) {
        mat[mat > scale_max] <- scale_max
    }
    return(mat)
    }






##############################################################################


inferloop.generateBin <- function(indata, used_coords, n=100, seed=123, type=NULL){
    DATA=indata
    used_coords=used_coords
    if(is.null(type)==FALSE){
        type_coords=as.numeric(as.factor(type))*10
        used_coords=cbind(used_coords, type_coords)
        }
    n=n
    set.seed(seed)
    KM=kmeans(used_coords,centers=n)
    CLST=KM$cluster
    names(CLST)=colnames(indata)
    OUT=.generate_mean(DATA, CLST)
    ##########################
    RETURN=list()
    RETURN$mat=OUT
    RETURN$clst=CLST
    return(RETURN)
    }



inferloop.bin2cell <-function(signal_mat, clst){
    CLST=clst
    SSN=signal_mat
    isLoop=matrix(0,nrow=nrow(SSN),ncol=length(CLST))
    rownames(isLoop)=rownames(SSN)
    colnames(isLoop)=names(CLST)
    UC=unique(CLST)
    i=1
    while(i<=length(UC)){
        this_clst=UC[i]
        this_index1=which(CLST==this_clst)
        this_index2=which(colnames(SSN)==this_clst)
        isLoop[,this_index1]=SSN[,this_index2]
        if(i%%20==1){print(i)}
        i=i+1
        }
    return(isLoop)
    }




######################################
# Cicero

inferloop.getGenomeDF.mm10 <-function(){
    library(GenomicRanges)
    library(BSgenome.Mmusculus.UCSC.mm10)
    genome <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
    genome <- genome[1:21]
    genome.df <- data.frame("chr" = names(genome), "length" = genome)
    return(genome.df)
    }

inferloop.getGenomeDF.hg38 <-function(){
    library(GenomicRanges)
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
    genome <- genome[1:24]
    genome.df <- data.frame("chr" = names(genome), "length" = genome)
    return(genome.df)
    }

inferloop.getGenomeDF.hg19 <-function(){
    library(GenomicRanges)
    library(BSgenome.Hsapiens.UCSC.hg19)
    genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
    genome <- genome[1:24]
    genome.df <- data.frame("chr" = names(genome), "length" = genome)
    return(genome.df)
    }






assemble_connections<-function (cicero_model_list, silent = FALSE){
    types <- vapply(cicero_model_list, FUN = class, FUN.VALUE = "character")
    char_hbn <- cicero_model_list[types == "character"]
    gl_only <- cicero_model_list[types == "list"]
    if (!silent) {
        print(paste("Successful cicero models: ", length(gl_only)))
        print("Other models: ")
        print(table(unlist(char_hbn)))
        print(paste("Models with errors: ", sum(is.null(cicero_model_list))))
    }
    cors <- lapply(gl_only, function(gl) {
        cors <- stats::cov2cor(gl$w)
        data.table::melt(data.table::as.data.table(cors, keep.rownames = TRUE),
            measure = patterns("[0-9]"))
    })
    cors <- data.table::rbindlist(cors)
    names(cors) <- c("Var1", "Var2", "value")
    data.table::setkey(cors, "Var1", "Var2")
    cors_rec <- as.data.frame(cors[, list(mean_coaccess = mean(value)),
        by = "Var1,Var2"])
    names(cors_rec) <- c("Peak1", "Peak2", "coaccess")
    cors_rec <- cors_rec[cors_rec$Peak1 != cors_rec$Peak2, ]
    return(cors_rec)
}


generate_windows <- function(window, genomic_coords) {
  if(!is(genomic_coords, "data.frame")) {
    chr_maxes <- read.table(genomic_coords)
  } else {
    chr_maxes <- genomic_coords
  }
  names(chr_maxes) <- c("V1", "V2")
  win_ranges <- plyr::ddply(chr_maxes, plyr::.(V1), function(x) {
    r <- seq(from = 1, to = x$V2[1], by = window/2)
    l <- r + window - 1
    data.frame(start = r, end = l)
  })
  gr <- GenomicRanges::GRanges(win_ranges$V1,
                               ranges=IRanges::IRanges(win_ranges$start,
                                                       win_ranges$end))
  return(gr)
}



get_genomic_range <- function(grs, cds, win) {
  end1 <- as.numeric(as.character(GenomicRanges::end(grs[win])))
  end2 <- as.numeric(as.character(GenomicRanges::start(grs[win])))
  win_range <- cds[(fData(cds)$bp1 < end1 &
                                fData(cds)$bp1 > end2) |
                               (fData(cds)$bp2 < end1 &
                                  fData(cds)$bp2 > end2), ]
  win_range <-
    win_range[as.character(fData(win_range)$chr) ==
                gsub("chr", "",
                     as.character(GenomicRanges::seqnames(grs[win]))),]
  fData(win_range)$mean_bp <-
    (as.numeric(as.character(fData(win_range)$bp1)) +
       as.numeric(as.character(fData(win_range)$bp2)))/2

  return(win_range)
}


calc_dist_matrix <- function(gene_range) {
  dist_mat <- as.matrix(dist(fData(gene_range)$mean_bp))
  row.names(dist_mat) <- colnames(dist_mat) <- row.names(fData(gene_range))

  return(dist_mat)
}


get_rho_mat <- function(dist_matrix, distance_parameter, s) {
  xmin <- 1000
  out <- (1-(xmin/dist_matrix)^s) * distance_parameter
  out[!is.finite(out)] <- 0
  out[out < 0] <- 0
  return(out)
}



generate_cicero_models<-function (cds, distance_parameter, s = 0.75, window = 5e+05,
    max_elements = 200, genomic_coords = cicero::human.hg19.genome)
{
    assertthat::assert_that(assertthat::is.number(s), s < 1, s > 0)
    grs <- generate_windows(window, genomic_coords)
    fData(cds)$chr <- gsub("chr", "", fData(cds)$chr)
    fData(cds)$bp1 <- as.numeric(as.character(fData(cds)$bp1))
    fData(cds)$bp2 <- as.numeric(as.character(fData(cds)$bp2))
    print('go')
    outlist <- parallel::mclapply(seq_len(length(grs)), mc.cores = CORES,
        function(win) {
            #print(win)
            GL <- "Error"
            win_range <- get_genomic_range(grs, cds, win)
            if (nrow(exprs(win_range)) <= 1) {
                return("Zero or one element in range")
            }
            if (nrow(exprs(win_range)) > max_elements) {
                return("Too many elements in range")
            }
            dist_matrix <- calc_dist_matrix(win_range)
            rho_mat <- get_rho_mat(dist_matrix, distance_parameter,
                s)
            vals <- exprs(win_range)
            cov_mat <- cov(t(vals))
            diag(cov_mat) <- diag(cov_mat) + 1e-04

            #############################################################
            #GL=glasso::glasso(cov_mat, rho_mat)
            GL=glassoFast::glassoFast(cov_mat, rho_mat)

            colnames(GL$w) <- row.names(GL$w) <- row.names(vals)
            colnames(GL$wi) <- row.names(GL$wi) <- row.names(vals)
            return(GL)
        })
    names_df <- as.data.frame(grs)
    names(outlist) <- paste(names_df$seqnames, names_df$start,
        names_df$end, sep = "_")
    return(outlist)
}


.run_cicero_new<-function(cds,genomic_coords , window=5e+05,sample_num=100){
    cds=cds
    genomic_coords=genomic_coords
    window = window
    silent = FALSE
    sample_num = sample_num
    ########################
    distance_parameters <- estimate_distance_parameter(cds, window = window,
        maxit = 100, sample_num = sample_num, distance_constraint = 250000,
        distance_parameter_convergence = 1e-22, genomic_coords = genomic_coords)
    mean_distance_parameter <- mean(unlist(distance_parameters))
    #####################
    cicero_out <- generate_cicero_models(cds, distance_parameter = mean_distance_parameter,
        window = window, genomic_coords = genomic_coords)
    all_cons <- assemble_connections(cicero_out, silent = silent)
    return(all_cons)
          }



inferloop.cicero<-function(indata, used_coords, genome.df, k=50, window=5e+05,sample_num=100){
    indata=as.matrix(indata)
    used_coords=used_coords
    genome.df=genome.df
    window=window
    k=k
    sample_num=sample_num
    #########################
    library(monocle)
    library(cicero)
    ###########################
    rownames(used_coords)=stringr::str_replace_all(rownames(used_coords),'-','_')
    colnames(indata)=stringr::str_replace_all(colnames(indata),'-','_')
    rownames(indata)=stringr::str_replace_all(rownames(indata),'-','_')
    peakinfo=as.data.frame(t(matrix(unlist(stringr::str_split(rownames(indata),'_')),nrow=3)))
    names(peakinfo) <- c("chr", "bp1", "bp2")
    peakinfo$site_name=rownames(indata)
    rownames(peakinfo) <- peakinfo$site_name 
    fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
    input_cds <-  suppressWarnings(newCellDataSet(indata,
                            phenoData = NULL,
                            featureData = fd,
                            expressionFamily=VGAM::binomialff(),
                            lowerDetectionLimit=0))
    cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = used_coords, k =k, size_factor_normalize = FALSE)
    conns <- .run_cicero_new(cicero_cds, genome.df, window, sample_num)
    return(conns)
    }



########################################################
#ciceroFrame start
##########################################################


assemble_connections_ciceroFrame<-function (cicero_model_list, silent = FALSE){
    types <- vapply(cicero_model_list, FUN = class, FUN.VALUE = "character")
    char_hbn <- cicero_model_list[types == "character"]
    gl_only <- cicero_model_list[types == "list"]
    cors <- lapply(gl_only, function(gl) {
        cors <- gl$cor
        data.table::melt(data.table::as.data.table(cors, keep.rownames = TRUE),
            measure = patterns("[0-9]"))
        })
    cors <- data.table::rbindlist(cors)
    names(cors) <- c("Var1", "Var2", "value")
    data.table::setkey(cors, "Var1", "Var2")
    cors_rec <- as.data.frame(cors[, list(mean_coaccess = mean(value)),
        by = "Var1,Var2"])
    names(cors_rec) <- c("Peak1", "Peak2", "coaccess")
    cors_rec <- cors_rec[cors_rec$Peak1 != cors_rec$Peak2, ]
    return(cors_rec)
    }


generate_ciceroFrame_models<-function (cds, window = 5e+05, s=0.75, max_elements = 200, genomic_coords = cicero::human.hg19.genome, used_function=used_function){
    assertthat::assert_that(assertthat::is.number(s), s < 1, s > 0)
    grs <- generate_windows(window, genomic_coords)
    fData(cds)$chr <- gsub("chr", "", fData(cds)$chr)
    fData(cds)$bp1 <- as.numeric(as.character(fData(cds)$bp1))
    fData(cds)$bp2 <- as.numeric(as.character(fData(cds)$bp2))
    print('go')
    outlist <- parallel::mclapply(seq_len(length(grs)), mc.cores = CORES,
        function(win) {
            #print(win)
            GL <- "Error"
            win_range <- get_genomic_range(grs, cds, win)
            if (nrow(exprs(win_range)) <= 1) {
                return("Zero or one element in range")
            }
            if (nrow(exprs(win_range)) > max_elements) {
                return("Too many elements in range")
            }
            win_range <- get_genomic_range(grs, cds, win)
            vals <- exprs(win_range)
            cor <- used_function(vals)
            GL$cor=cor
            colnames(GL$cor) <- row.names(GL$cor) <- row.names(vals)
            return(GL)
        })
    names_df <- as.data.frame(grs)
    names(outlist) <- paste(names_df$seqnames, names_df$start,
        names_df$end, sep = "_")
    return(outlist)
    }

.run_ciceroFrame<-function(cds,genomic_coords , window=5e+05,sample_num=100,used_function){
    cds=cds
    genomic_coords=genomic_coords
    window = window
    silent = FALSE
    sample_num = sample_num
    ########################
    cicero_out <- generate_ciceroFrame_models(cds, window = window, genomic_coords = genomic_coords,used_function=used_function)
    all_cons <- assemble_connections_ciceroFrame(cicero_out, silent = silent)
    return(all_cons)
          }


inferc.ciceroFrame<-function(indata, used_function,used_coords, genome.df, k=50, window=5e+05,sample_num=100){
    indata=as.matrix(indata)
    used_coords=used_coords
    genome.df=genome.df
    window=window
    k=k
    sample_num=sample_num
    #########################
    library(monocle)
    library(cicero)
    ###########################
    rownames(used_coords)=stringr::str_replace_all(rownames(used_coords),'-','_')
    colnames(indata)=stringr::str_replace_all(colnames(indata),'-','_')
    rownames(indata)=stringr::str_replace_all(rownames(indata),'-','_')
    peakinfo=as.data.frame(t(matrix(unlist(stringr::str_split(rownames(indata),'_')),nrow=3)))
    names(peakinfo) <- c("chr", "bp1", "bp2")
    peakinfo$site_name=rownames(indata)
    rownames(peakinfo) <- peakinfo$site_name
    fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
    input_cds <-  suppressWarnings(newCellDataSet(indata,
                            phenoData = NULL,
                            featureData = fd,
                            expressionFamily=VGAM::binomialff(),
                            lowerDetectionLimit=0))
    print('step01: prepare input...')
    cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = used_coords, k =k, size_factor_normalize = FALSE)
    print('step02: run ciceroFrame...')
    conns <- .run_ciceroFrame(cicero_cds, genome.df, window, sample_num, used_function)
    return(conns)
    }

inferc.ciceroFrame_step01_prepareInput<-function(indata, used_coords, genome.df, k=50, window=5e+05,sample_num=100){
    indata=as.matrix(indata)
    used_coords=used_coords
    genome.df=genome.df
    window=window
    k=k
    sample_num=sample_num
    #########################
    library(monocle)
    library(cicero)
    ###########################
    rownames(used_coords)=stringr::str_replace_all(rownames(used_coords),'-','_')
    colnames(indata)=stringr::str_replace_all(colnames(indata),'-','_')
    rownames(indata)=stringr::str_replace_all(rownames(indata),'-','_')
    peakinfo=as.data.frame(t(matrix(unlist(stringr::str_split(rownames(indata),'_')),nrow=3)))
    names(peakinfo) <- c("chr", "bp1", "bp2")
    peakinfo$site_name=rownames(indata)
    rownames(peakinfo) <- peakinfo$site_name
    fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
    input_cds <-  suppressWarnings(newCellDataSet(indata,
                            phenoData = NULL,
                            featureData = fd,
                            expressionFamily=VGAM::binomialff(),
                            lowerDetectionLimit=0))
    print('step01: prepare input...')
    cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = used_coords, k =k, size_factor_normalize = FALSE)
    return(cicero_cds)
    }


inferc.ciceroFrame_step02_runUsedFuntion<-function(cicero_cds, genome.df, window, sample_num, used_function){
    cicero_cds=cicero_cds
    genome.df=genome.df
    window=window
    sample_num=sample_num
    used_function=used_function
    ########################
    library(monocle)
    library(cicero)
    #########################
    print('step02: run ciceroFrame...')
    conns <- .run_ciceroFrame(cicero_cds, genome.df, window, sample_num, used_function)
    return(conns)
    }


inferc.ciceroFrame_step02_runCicero<-function(cicero_cds, genome.df, window, sample_num){
    cicero_cds=cicero_cds
    genome.df=genome.df
    window=window
    sample_num=sample_num
    print('step02: run cicero...')
    conns <- .run_cicero_new(cicero_cds, genome.df, window, sample_num)
    return(conns)
    }
#####################################




assemble_sd_connections<-function (cicero_model_list, silent = FALSE){
    types <- vapply(cicero_model_list, FUN = class, FUN.VALUE = "character")
    char_hbn <- cicero_model_list[types == "character"]
    gl_only <- cicero_model_list[types == "list"]
    if (!silent) {
        print(paste("Successful cicero models: ", length(gl_only)))
        print("Other models: ")
        print(table(unlist(char_hbn)))
        print(paste("Models with errors: ", sum(is.null(cicero_model_list))))
    }
    cors <- lapply(gl_only, function(gl) {
        cors <- gl$d
        data.table::melt(data.table::as.data.table(cors, keep.rownames = TRUE),
            measure = patterns("[0-9]"))
    })
    cors <- data.table::rbindlist(cors)
    names(cors) <- c("Var1", "Var2", "value")
    data.table::setkey(cors, "Var1", "Var2")
    cors_rec <- as.data.frame(cors[, list(mean_coaccess = mean(value)),
        by = "Var1,Var2"])
    names(cors_rec) <- c("Peak1", "Peak2", "coaccess")
    cors_rec <- cors_rec[cors_rec$Peak1 != cors_rec$Peak2, ]
    return(cors_rec)
}


generate_sd_models<-function (cds, distance_parameter, s = 0.75, window = 5e+05,
    max_elements = 200, genomic_coords = cicero::human.hg19.genome)
{
    assertthat::assert_that(assertthat::is.number(s), s < 1, s > 0)
    grs <- generate_windows(window, genomic_coords)
    fData(cds)$chr <- gsub("chr", "", fData(cds)$chr)
    fData(cds)$bp1 <- as.numeric(as.character(fData(cds)$bp1))
    fData(cds)$bp2 <- as.numeric(as.character(fData(cds)$bp2))
    print('go')
    outlist <- parallel::mclapply(seq_len(length(grs)), mc.cores = CORES,
        function(win) {
            #print(win)
            GL <- "Error"
            win_range <- get_genomic_range(grs, cds, win)
            if (nrow(exprs(win_range)) <= 1) {
                return("Zero or one element in range")
            }
            if (nrow(exprs(win_range)) > max_elements) {
                return("Too many elements in range")
            }
            dist_matrix <- calc_dist_matrix(win_range)
            rho_mat <- get_rho_mat(dist_matrix, distance_parameter,
                s)
            vals <- exprs(win_range)
            #cov_mat <- cov(t(vals))
            vals=as.matrix(vals)
            cov_mat=vals %*% t(vals) / ncol(vals)
            ##################################### 
            diag(cov_mat) <- diag(cov_mat) + 1e-04

            #############################################################
            #GL=glasso::glasso(cov_mat, rho_mat)
            GL=glassoFast::glassoFast(cov_mat, rho_mat)
           
            colnames(GL$w) <- row.names(GL$w) <- row.names(vals)
            colnames(GL$wi) <- row.names(GL$wi) <- row.names(vals)
            #############################################################
            this_a=GL$w * ncol(vals)
            this_b=matrix(rep(rowSums(vals**2),nrow(vals)),ncol=nrow(vals))
            this_c= matrix(rep(rowSums(vals**2),each=nrow(vals)),ncol=nrow(vals))
            this_d= this_a/sqrt(this_b*this_c)
             ###########################################################
            #this_bc=cbind(as.numeric(this_b),as.numeric(this_c)) 
            #this_bc=rowMax(this_bc)
            #this_bc=matrix(this_bc,ncol=ncol(this_b))
            #this_d= this_a/this_bc
            #####################################################
            GL$d=this_d
            colnames(GL$d)<-row.names(GL$d)<-row.names(vals)

            return(GL)
        })
    names_df <- as.data.frame(grs)
    names(outlist) <- paste(names_df$seqnames, names_df$start,
        names_df$end, sep = "_")
    return(outlist)
}



.run_sd_new<-function(cds,genomic_coords , window=5e+05,sample_num=100){
    cds=cds
    genomic_coords=genomic_coords
    window = window
    silent = FALSE
    sample_num = sample_num
    ########################
    distance_parameters <- estimate_distance_parameter(cds, window = window,
        maxit = 100, sample_num = sample_num, distance_constraint = 250000,
        distance_parameter_convergence = 1e-22, genomic_coords = genomic_coords)
    mean_distance_parameter <- mean(unlist(distance_parameters))
    #####################
    cicero_out <- generate_sd_models(cds, distance_parameter = mean_distance_parameter,
        window = window, genomic_coords = genomic_coords)
    all_cons <- assemble_sd_connections(cicero_out, silent = silent)
    return(all_cons)
          }



inferc.ciceroFrame_step02_runSD<-function(cicero_cds, genome.df, window, sample_num){
    cicero_cds=cicero_cds
    genome.df=genome.df
    window=window
    sample_num=sample_num
    print('step02: run cicero...')
    conns <- .run_sd_new(cicero_cds, genome.df, window, sample_num)
    return(conns)
    }




####################################







used_function_cor<-function(vals){
    this_out=cor(t(vals),t(vals))
    this_out[which(is.na(this_out))]=0
    return(this_out)
    }

used_function_spearman<-function(vals){
    this_out=cor(t(vals),t(vals),method='spearman')
    this_out[which(is.na(this_out))]=0
    return(this_out)
    }



#####################
#ciceroFrame end
#####################









######################################
# Archr

archr.computeKNN <- function(
  data = NULL,
  query = NULL,
  k = 50,
  includeSelf = FALSE,
  ...
  ){
  if(is.null(query)){
    query <- data
    searchSelf <- TRUE
  }else{
    searchSelf <- FALSE
  }
  if(searchSelf & !includeSelf){
    knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
    knnIdx <- knnIdx[,-1,drop=FALSE]
  }else{
    knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
  }
  knnIdx
}


archr.getQuantiles <- function(v = NULL, len = length(v)){
  if(length(v) < len){
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  }else{
    v2 <- v
  }
  p <- trunc(rank(v2))/length(v2)
  if(length(v) < len){
    p <- p[seq_along(v)]
  }
  return(p)
}


determineOverlapCpp <- function(m, overlapCut) {
    .Call('_ArchR_determineOverlapCpp', PACKAGE = 'ArchR', m, overlapCut)
}

rowCorCpp <- function(idxX, idxY, X, Y) {
    .Call('_ArchR_rowCorCpp', PACKAGE = 'ArchR', idxX, idxY, X, Y)
}


inferloop.archrCoA<-function(MAT,VEC,k=50,maxDist=500000,SEED=123){
    MAT=MAT
    VEC=VEC
    SEED=SEED
    set.seed(SEED)
    k = k
    knnIteration = 500
    overlapCutoff = 0.8
    maxDist = maxDist
    scaleTo = 10^4
    log2Norm = TRUE
    threads = 10
    verbose = TRUE
    #################

    rD=VEC
    idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)
    knnObj <- archr.computeKNN(data = rD, query = rD[idx,], k = k)
    keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))
    knnObj <- knnObj[keepKnn==0,]
    knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
        rownames(rD)[knnObj[x, ]]
          }) %>% SimpleList
    peak=rownames(MAT)
    bed=inferloop.splitLoop(peak,'-',3)
    bed=as.data.frame(bed)
    colnames(bed)=c('chr','start','end')
    bed$start=as.integer(bed$start)
    bed$end=as.integer(bed$end)
    peakSet <- with(bed, GRanges(chr, IRanges(start, end)))
    peakSet$idx=c(1:nrow(bed))
    peakSummits <- resize(peakSet, 1, "center")
    peakWindows <- resize(peakSummits, 2*maxDist + 1, "center")

    o <- DataFrame(findOverlaps(peakSummits, peakWindows, ignore.strand = TRUE))
    o <- o[o[,1] != o[,2],]
    o$seqnames <- seqnames(peakSet)[o[,1]]
    o$idx1 <- peakSet$idx[o[,1]]
    o$idx2 <- peakSet$idx[o[,2]]
    o$correlation <- -999.999
    o$Variability1 <- 0.000
    o$Variability2 <- 0.000

    umat=MAT[,unlist(knnObj)]
    gmat=.generate_mean(umat, rep(1:length(knnObj),each=k))
    gS=colSums(gmat)
    groupMat=gmat
    groupMat <- t(t(groupMat) / gS) * scaleTo
        if(log2Norm){
            groupMat <- log2(groupMat + 1)
            }

    CHR=bed[,1]
    chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))

    for(x in seq_along(chri)){
        this_chr=chri[x]
        print(this_chr)
        this_row_index=which(CHR %in% this_chr)
        #Correlations
        idx <- BiocGenerics::which(o$seqnames==chri[x])
        corVals <- rowCorCpp(idxX = o[idx,]$idx1, idxY = o[idx,]$idx2, X = as.matrix(groupMat), Y = as.matrix(groupMat))
        rowVars <- as.numeric(matrixStats::rowVars(groupMat))
        o[idx,]$correlation <- as.numeric(corVals)
        o[idx,]$Variability1 <- rowVars[o[idx,]$idx1]
        o[idx,]$Variability2 <- rowVars[o[idx,]$idx2]
        }
    o$idx1 <- NULL
    o$idx2 <- NULL
    o <- o[!is.na(o$correlation),]
    o$TStat <- (o$correlation / sqrt((pmax(1-o$correlation^2, 0.00000000000000001, na.rm = TRUE))/(length(knnObj)-2))) #T-statistic P-value
    o$Pval <- 2*pt(-abs(o$TStat), length(knnObj) - 2)
    o$FDR <- p.adjust(o$Pval, method = "fdr")

    o$VarQuantile1 <- archr.getQuantiles(o$Variability1)
    o$VarQuantile2 <- archr.getQuantiles(o$Variability2)
    mcols(peakSet) <- NULL
    o@metadata$peakSet <- peakSet

    o$chr1=o$seqnames
    o$start1=bed$start[o$queryHits]
    o$end1=bed$end[o$queryHits]
    o$chr2=o$seqnames
    o$start2=bed$start[o$subjectHits]
    o$end2=bed$end[o$subjectHits]
    o$distance = abs((o$start1+o$end1)/2-(o$start2+o$end2)/2)

    return(o)
    }


###################################################

.cart2clock<-function(x,y,circle){
    phi <- (atan(x/y)/2/pi * circle + ifelse(y >= 0, circle,1.5 * circle))
    phi <- phi %% circle
    output=data.frame(rho = sqrt(x * x + y * y), phi = phi)
    output=as.matrix(output)
    return(output)
    }

.clock2cart<-function(rho, phi, circle){
    output=data.frame(x = rho * sin(phi/circle * 2 * pi), y = rho *
        cos(phi/circle * 2 * pi))
    output=as.matrix(output)
    return(output)
    }



inferloop.splitILS<-function(X, Y, r=0){
    X=X
    Y=Y
    r=r
    #######################
    ############################
    # Ensure positive value
    if( min(X)<0 ){
        X = X - min(inferloop.rmOut(X))
        X[which(X<0)]=0
        }
    if( min(Y)<0 ){
        Y = Y - min(inferloop.rmOut(Y))
        Y[which(Y<0)]=0
        }
    #######################
    X_mean = mean(X)
    Y_mean = mean(Y)
    X_base = X_mean * r
    Y_base = Y_mean * r
    X_delta = X - X_base
    Y_delta = Y - Y_base
    ######################
    ILS=inferloop.calILS(X,Y,r)
    #######################
    CVEC=as.matrix(.cart2clock(X_delta,Y_delta,360))
    CVEC[which(is.na(CVEC))]=0
    ANGLE=CVEC[,2]
    ##############################
    ILS_FLAG=rep(-1,length(ILS))
    ILS_FLAG[which(ILS>0)]=1
    #############################
    this_order=order(ANGLE)
    ANGLE_order=ANGLE[this_order]
    ILS_FLAG_order=ILS_FLAG[this_order]
    ##############################
    SM_ILS_FLAG_order=smooth.spline(ILS_FLAG_order)$y
    SM_ILS_FLAG_order_fat=c(SM_ILS_FLAG_order[length(SM_ILS_FLAG_order)],
                         SM_ILS_FLAG_order,
                         SM_ILS_FLAG_order[1])
    #################################################
    ILS_FLAG_order_fat=SM_ILS_FLAG_order_fat
    ILS_FLAG_order_fat[which(SM_ILS_FLAG_order_fat>0)]=1
    ILS_FLAG_order_fat[which(SM_ILS_FLAG_order_fat<=0)]=-1
    ######################################################
    FLAG_order_fat=rep(0,length(ILS_FLAG_order_fat))
    i=2
    while(i<length(ILS_FLAG_order_fat)){
        this_before=ILS_FLAG_order_fat[i-1]
        this_after=ILS_FLAG_order_fat[i+1]
        if(this_before!=this_after & FLAG_order_fat[i-1]!=1){FLAG_order_fat[i]=1}
        i=i+1}
    FLAG_order=FLAG_order_fat[2:(length(ILS_FLAG_order)+1)]
    ##############################
    CUT_ANGLE=ANGLE_order[which(FLAG_order>0)]
    ############################
    CLST=rep(1,length(ANGLE))
    if(max(CUT_ANGLE)-min(CUT_ANGLE)>180){
        CLST[which(ANGLE>max(CUT_ANGLE))]=1
        }else{
        CLST[which(ANGLE>max(CUT_ANGLE))]=length(CUT_ANGLE)+1
        }
    i=1
    while(i<length(CUT_ANGLE)){
        this_lw=CUT_ANGLE[i]
        this_up=CUT_ANGLE[i+1]
        CLST[which(ANGLE>this_lw & ANGLE <= this_up)]=i+1
        i=i+1}
    #########################
    OUT=list()
    OUT$ils=ILS
    OUT$clst=CLST
    return(OUT)
    }

#######################################
#20230519


inferloop.calD <-function(X, Y, r=0){
    X=X
    Y=Y
    r=r
    ############################
    # Ensure positive value
    if( min(X)<0 ){
        X = X - min(inferloop.rmOut(X))
        X[which(X<0)]=0
        }
    if( min(Y)<0 ){
        Y = Y - min(inferloop.rmOut(Y))
        Y[which(Y<0)]=0
        }
    #############################
    r=r
    N=length(X)
    ############################
    X_mean = mean(X)
    Y_mean = mean(Y)
    ###########################
    X_base = X_mean * r
    Y_base = Y_mean * r
    X_delta = X - X_base
    Y_delta = Y - Y_base
    ###########################
    ABCD = inferloop.calABCD(X, Y, X_base, Y_base)
    ###########################
    D = ABCD[['D']]
    return(D)
    }


inferloop.inferD<-function(mat, net, r=0,sep='.And.'){
    library(hash)
    r=r # default r=0
    sep=sep
    ###################
    mat=as.matrix(mat)
    net=as.matrix(net)
    tag=rownames(mat)
    net=net[which(net[,1] %in% tag & net[,2] %in% tag),]
    #######################################
    print('hashing...')
    h=hash(keys=tag,values=rep(0,nrow(mat)))
    i=1
    while(i<=nrow(mat)){
        this_tag=tag[i]
        this_v=as.numeric(mat[i,])
        .set(h, keys=this_tag, values=this_v)
        i=i+1}
    #######################################
    print('calculating ILS...')
    out=matrix(0,ncol=1,nrow=nrow(net))
    colnames(out)='D'
    rownames(out)=paste0(net[,1],sep,net[,2])
    i=1
    while(i<=nrow(net)){
        this_tag1=net[i,1]
        this_tag2=net[i,2]
        #######################
        x=as.vector(hash::values(x=h, keys=this_tag1)[,1])
        y=as.vector(hash::values(x=h, keys=this_tag2)[,1])
        z=inferloop.calD(x,y,r=r)
        out[i,]=z
        ########################
        if(i%%10000==1){print(paste0(i,' / ',nrow(net)))}
        i=i+1}
    print('finished!')
    return(out)
    }








##########################################################
# InferC functions

inferc.logit<-function(p){return(log(p/(1-p)))}

inferc.sigmoid<-function(x,a=0.5,b=10){return(1/(1+exp(-(x-a)*b)))}



inferc.getUniqConns<-function(net){
    #######################
    library(hash)
    #######################
    h=hash()
    #######################
    print('getting unique loop...')
    flag=rep(0,nrow(net))
    i=1
    while(i<=nrow(net)){
        this_tag1=tag1[i]
        this_tag2=tag2[i]
        if(has.key(this_tag1,h) | has.key(this_tag2,h)){
            flag[i]=1
            #del(this_tag1,h)
            #del(this_tag2,h)
            }else{
            .set(h, keys=this_tag1, values=1)
            .set(h, keys=this_tag2, values=1)
            }
        if(i %%50000==1){print(paste0(i,' / ',nrow(net)))}
        i=i+1}
    #########################
    out=net[which(flag==0),]
    print(paste0('Number of unique loops: ',nrow(out)))
    print('finished!')
    return(out)
    }




inferc.generate_mean <- function(exp_sc_mat, TAG, print_step=10){
    print_step=print_step
    exp_sc_mat=exp_sc_mat
    TAG=TAG
    ########
    NewRef=matrix(0,ncol=length(unique(TAG)),nrow=nrow(exp_sc_mat))
    ########
    TAG=as.character(TAG)
    refnames=unique(TAG)
    total_num=length(refnames)
    outnames=c()
    i=1
    while(i<=length(refnames)){
        one=refnames[i]
        this_col=which(TAG==one)
        outnames=c(outnames,one)
        if(length(this_col) >1){
            #this_new_ref=apply(exp_sc_mat[,this_col],1,mean)
            this_new_ref=rowMeans(exp_sc_mat[,this_col])
            }else{
            this_new_ref = exp_sc_mat[,this_col]
            }
        NewRef[,i]=this_new_ref
        if(i%%print_step==1){print(paste0(i,' / ' ,total_num ))}
        i=i+1
        }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    return(NewRef)
    }


inferc.getChrLoc<-function(TAG, SPLIT1='-', SPLIT2='-'){
    TAG=TAG
    SPLIT1=SPLIT1
    SPLIT2=SPLIT2
    OUT=list()
    TMP1=TAG
    if(SPLIT1==SPLIT2){
        TMP2=t(matrix(unlist(stringr::str_split(TMP1,SPLIT1)),nrow=3))
        OUT$chr=TMP2[,1]
        OUT$start=as.numeric(TMP2[,2])
        OUT$end=as.numeric(TMP2[,3])
        }else{
        TMP2=t(matrix(unlist(stringr::str_split(TMP1,SPLIT1)),nrow=2))
        OUT$chr=TMP2[,1]
        TMP3=t(matrix(unlist(stringr::str_split(TMP2[,2],SPLIT2)),nrow=2))
        OUT$start=as.numeric(TMP3[,1])
        OUT$end=as.numeric(TMP3[,2])
        }
    return(OUT)
    }



inferc.colMap<-function(x, value, color){
    value=value
    color=color
    #####################
    colMat=col2rgb(color)
    RED=colMat[1,]    
    GRE=colMat[2,]    
    BLU=colMat[3,]    
    ##########################     
    red.fit=approxfun(RED~value)
    gre.fit=approxfun(GRE~value)
    blu.fit=approxfun(BLU~value)
    ##########################
    OUT=x
    OUT[which(x<=min(value))]=color[1]
    OUT[which(x>=max(value))]=color[length(color)]
    used_index=which(x>min(value) & x <max(value))
    tmp=x[used_index]
    tmp_col=rgb(red.fit(tmp),gre.fit(tmp),blu.fit(tmp),maxColorValue=255)
    OUT[used_index]=tmp_col
    return(OUT)
    }


inferc.runTFIDF <- function(mat, method=2, scaleTo=10000){
    mat=mat
    method=method
    scaleTo=scaleTo
    #############
    #Adapted from ArchR
    #############
    mat=Matrix::Matrix(mat)
    colSm <- Matrix::colSums(mat)
    rowSm <- Matrix::rowSums(mat)
    idx <- which(rowSm > 0)
    mat <- mat[idx, ]
    rowSm <- rowSm[idx]
    #TF - Normalize
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))
    #########################
    if(method == 1 | tolower(method) == "tf-logidf"){
      print('tf-logidf | Adapted from Casanovich et al.')
      #Adapted from Casanovich et al.
      #LogIDF
      idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
      #TF-LogIDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
    }else if(method == 2 | tolower(method) == "log(tf-idf)"){
      print('log(tf-idf) | Adapted from Stuart et al.')
      #Adapted from Stuart et al.
      #IDF
      idf   <- as(ncol(mat) / rowSm, "sparseVector")
      #TF-IDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
      mat@x <- log(mat@x * scaleTo + 1)
    }else if(method == 3 | tolower(method) == "logtf-logidf"){
      print('logtf-logidf | Adapted from the method 3 of ArchR.')
      #LogTF
      mat@x <- log(mat@x + 1)
      #LogIDF
      idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
      #TF-IDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
      }
    return(mat)
    }



inferc.aggMat<-function(MAT, VEC, clstNum=100, seed=123, BATCH=NULL){
    MAT=MAT
    VEC=VEC
    BATCH=BATCH
    clstNum=clstNum
    seed=seed
    ##################
    set.seed(seed)
    N=clstNum
    #########################
    if(is.null(BATCH)){
        KM=kmeans(VEC, centers=N)
        CLST=KM$cluster
    }else{
        TAB=table(BATCH)
        BATCH.N=(N*TAB)%/%sum(TAB)
        BATCH.N[which(BATCH.N==min(BATCH.N))[1]]=BATCH.N[which(BATCH.N==min(BATCH.N))[1]]+N-sum(BATCH.N)
        print(BATCH.N)
        CLST=rep(0,length(BATCH))
        j=0
        i=1
        while(i<=length(BATCH.N)){
           this_n=BATCH.N[i]
           this_index=which(BATCH==names(BATCH.N)[i])
           this_vec=VEC[this_index,]
           this_km=kmeans(this_vec, centers=this_n)
           this_clst=this_km$cluster
           CLST[this_index]=j+this_km$cluster
           j=j+max(this_clst)
           i=i+1}
        }
    names(CLST)=colnames(MAT)
    AGG=inferc.generate_mean(MAT, CLST)
    AGG=AGG[,order(as.numeric(colnames(AGG)))]
    ##################
    OUT=list()
    #OUT$km=KM
    OUT$clst=CLST
    OUT$agg=AGG
    return(OUT)
    }

inferc.agg2cell<-function(AGG, CLST){
    CLST=CLST
    SSN=AGG
    ###################
    isLoop=matrix(0,nrow=nrow(SSN),ncol=length(CLST))
    rownames(isLoop)=rownames(SSN)
    colnames(isLoop)=names(CLST)
    UC=unique(CLST)
    i=1
    while(i<=length(UC)){
        this_clst=UC[i]
        this_index1=which(CLST==this_clst)
        this_index2=which(colnames(SSN)==this_clst)
        isLoop[,this_index1]=SSN[,this_index2]
        if(i%%20==1){print(i)}
        i=i+1
        }
    return(isLoop)
    }


inferc.posCor<-function(mat1,mat2,minN=3){
    mat1=mat1
    mat2=mat2
    minN=minN
    ################################
    out=matrix(0,nrow=ncol(mat1),ncol=ncol(mat2))
    rownames(out)=colnames(mat1)
    colnames(out)=colnames(mat2)
    i=1
    while(i<=ncol(mat1)){
        j=1
        while(j<=ncol(mat2)){
            x=mat1[,i]
            y=mat2[,j]
            used_index=which(x>0 & y>0)
            if(length(used_index)>=minN){
                this_cor=cor(x[used_index],y[used_index])
                out[i,j]=this_cor}
            j=j+1}
        i=i+1}
    #################################
    return(out)
    }


inferc.rowCalD<-function(mat){
    UMAT=mat
    DMAT=matrix(0,nrow=nrow(UMAT),ncol=nrow(UMAT))
    rownames(DMAT)=rownames(UMAT)
    colnames(DMAT)=rownames(UMAT)
    i=1
    while(i<=nrow(DMAT)){
        this_x=UMAT[i,]
        j=i
        while(j<=ncol(DMAT)){
            this_y=UMAT[j,]
            this_D=inferloop.calABCD(X=this_x, Y=this_y, X_base=0, Y_base=0)$D
            DMAT[i,j]=this_D
            j=j+1}
        i=i+1}
    #################################
    DMAT=(DMAT+t(DMAT))
    diag(DMAT)=diag(DMAT)/2
    ################################
    DMAT[which(is.na(DMAT))]=0
    return(DMAT)
    }





inferc.calScore<-function(X, Y){
    X=X
    Y=Y
    ############################
    # Ensure positive value
    if( min(X)<0 ){
        X = X - min(inferloop.rmOut(X))
        X[which(X<0)]=0
        }
    if( min(Y)<0 ){
        Y = Y - min(inferloop.rmOut(Y))
        Y[which(Y<0)]=0
        }
    #############################
    N=length(X)
    ############################
    X_mean = mean(X)
    Y_mean = mean(Y)
    ###########################
    X_base = 0
    Y_base = 0
    X_delta = X - X_base
    Y_delta = Y - Y_base
    ###########################
    ABCD = inferloop.calABCD(X, Y, X_base, Y_base)
    ###########################
    D = ABCD[['D']]
    D_plus = ( ABCD[['A']] + X_delta * Y_delta ) / sqrt( (ABCD[['B']] + X_delta**2) * (ABCD[['C']] + Y_delta**2) )
    ###########################
    M = D_plus - D
    S = (1-D**2)/(N-1)
    ###########################
    LS = M / S
    ###########################
    LS[which(is.na(LS))]=0
    D[which(is.na(D))]=0
    ###########################
    OUT=list()
    OUT$local=LS
    OUT$global=D
    return( OUT )
    }




inferc.calLoopScore<-function(mat, net,sep='.And.', a=0.5, b=10){
    library(hash)
    sep=sep
    a=a
    b=b
    ####################
    mat=as.matrix(mat)
    net=as.matrix(net)
    tag=rownames(mat)
    net=net[which(net[,1] %in% tag & net[,2] %in% tag),]
    #######################################
    print('hashing...')
    h=hash(keys=tag,values=rep(0,nrow(mat)))
    i=1
    while(i<=nrow(mat)){
        this_tag=tag[i]
        this_v=as.numeric(mat[i,])
        .set(h, keys=this_tag, values=this_v)
        i=i+1}
    #######################################
    print('calculating GLS & LLS...')
    gls=rep(0,nrow(net))
    names(gls)=paste0(net[,1],sep,net[,2])
    #########################################
    out=matrix(0,ncol=ncol(mat),nrow=nrow(net))
    colnames(out)=colnames(mat)
    rownames(out)=paste0(net[,1],sep,net[,2])
    i=1
    while(i<=nrow(net)){
        this_tag1=net[i,1]
        this_tag2=net[i,2]
        #######################
        x=as.vector(hash::values(x=h, keys=this_tag1)[,1])
        y=as.vector(hash::values(x=h, keys=this_tag2)[,1])
        #######################
        this_out=inferc.calScore(x,y)
        gls[i]=this_out$global
        out[i,]=this_out$local
        ########################
        if(i%%10000==1){print(paste0(i,' / ',nrow(net)))}
        i=i+1}
    #############################
    print('utilizing weight...')
    W=inferc.sigmoid(gls,a=a,b=b)
    wout=out * W
    #############################
    print('finished!')
    ############################
    OUT=list()
    OUT$gls=gls
    OUT$weight=W
    OUT$ori.lls=out
    OUT$lls=wout
    return(OUT)
    }







inferc.calCorLoop<-function(MAT, CHR, START, END, Dcut=0.3, SPLIT1='-',SPLIT2='-',SPLIT_NEW='_x_'){
    MAT=MAT
    SPLIT_NEW=SPLIT_NEW
    #####################
    CHR_INPUT=CHR
    START_INPUT=as.numeric(START)
    END_INPUT=as.numeric(END)
    ########################
    #####################
    Dcut=Dcut
    SPLIT1=SPLIT1
    SPLIT2=SPLIT2
    #####################
    AGG=MAT
    TMP1=inferc.getChrLoc(rownames(AGG),SPLIT1=SPLIT1,SPLIT2=SPLIT2)
    CHR=TMP1$chr
    START=TMP1$start
    END=TMP1$end
    LOC=(START+END)/2
    ######################
    rownames(AGG)=paste0(CHR,SPLIT_NEW,START,SPLIT_NEW,END)
    #####################
    USED_INDEX=which(CHR==CHR_INPUT & LOC >= START_INPUT & LOC <= END_INPUT)
    print(length(USED_INDEX))
    UMAT=AGG[USED_INDEX,]
    ############################
    if(min(as.numeric(UMAT))<0){print('Remove negative values in MAT !!!');UMAT[which(UMAT<0)]=0}
    ############################
    UCHR=CHR[USED_INDEX]
    ULOC=LOC[USED_INDEX]
    UNAME=rownames(UMAT)
    ############################
    DMAT=inferc.rowCalD(UMAT)
    #################################
    ##############################
    CUT=Dcut
    mat=AGG[USED_INDEX,]
    vvv=as.numeric(DMAT)
    ################################
    p1=rep(rownames(DMAT),times=ncol(DMAT),each=1)
    p2=rep(colnames(DMAT),times=1,each=nrow(DMAT))
    ##################################
    over_cut=which(vvv>CUT)
    ##################################
    net=cbind(p1[over_cut],p2[over_cut])
    net_uniq=inferloop.getUniqLoop(net)
    over_cut_uniq=which(paste0(net[,1],'.And.',net[,2]) %in% paste0(net_uniq[,1],'.And.',net_uniq[,2]))
    net_uniq_index=over_cut[over_cut_uniq]
    ##########################
    ################################
    mat_plus=matrix(0,nrow=2,ncol=ncol(mat))
    rownames(mat_plus)=c('1','2')
    rownames(mat_plus)[1]=paste0(CHR_INPUT[1],SPLIT_NEW,START_INPUT-1,SPLIT_NEW,START_INPUT)
    rownames(mat_plus)[2]=paste0(CHR_INPUT[1],SPLIT_NEW,END_INPUT-1,SPLIT_NEW,END_INPUT)
    net_uniq_plus=c(rownames(mat_plus)[1],rownames(mat_plus)[2])
    #############################
    inferc.out=inferc.calLoopScore(rbind(mat_plus,mat), rbind(net_uniq_plus,net_uniq))
    #############################
    OUT=list()
    OUT$dmat=DMAT
    OUT$mat=mat
    OUT$net=net_uniq
    OUT$net_index=net_uniq_index
    OUT$loopScore=inferc.out
    OUT$SPLIT_NEW=SPLIT_NEW
    OUT$mat_plus=mat_plus
    OUT$net_plus=net_uniq_plus
    ##############################
    return(OUT)
    }



inferc.preContact<-function(LOOP_SIGNAL,Hh=1,Zh=1,Ch=1,Zmax=2.5,maxCol.Z='.',lineCol='grey50',minCol='grey90',maxCol='red1', exMin=0, MAIN='',SPLIT_NEW='_x_'){
    LOOP_SIGNAL=LOOP_SIGNAL
    Hh=Hh
    Zh=Zh
    Ch=Ch
    Zmax=Zmax
    maxCol.Z=maxCol.Z
    if(maxCol.Z=='.'){maxCol.Z=Zmax}
    lineCol=lineCol
    minCol=minCol
    maxCol=maxCol
    exMin=exMin
    MAIN=MAIN
    SPLIT_NEW=SPLIT_NEW
    ##########################
    LOOP=inferloop.splitLoop(names(LOOP_SIGNAL))
    ##########################
    OUT1=inferc.getChrLoc(LOOP[,1],SPLIT_NEW,SPLIT_NEW)
    OUT2=inferc.getChrLoc(LOOP[,2],SPLIT_NEW,SPLIT_NEW)
    CHR1=OUT1$chr
    CHR2=OUT2$chr
    LOC1=(OUT1$start+OUT1$end )/2
    LOC2=(OUT2$start+OUT2$end )/2
    ############################
    ############################
    Z=LOOP_SIGNAL
    ###########################
    X=(LOC1+LOC2)/2
    Y=abs(LOC2-LOC1)/2
    Xall=c(LOC1,LOC2)
    ###########################
    Hbase=Ch+Zh
    Y=Y/max(Y)*Hh
    ###########################
    COL=inferc.colMap(Z,c(0, maxCol.Z),c(minCol,maxCol))
    #############
    Zpos=Z-0
    Zpos[which(Zpos<0)]=0
    Zpos[which(Zpos>Zmax)]=Zmax
    ZZZ= Hbase - Zpos/Zmax * Zh
    #############
    OUT=list()
    OUT$LOOP_SIGNAL=LOOP_SIGNAL
    OUT$Hh=Hh
    OUT$Zh=Zh
    OUT$Ch=Ch
    OUT$Hbase=Hbase
    OUT$Z=Z
    OUT$X=X
    OUT$Y=Y
    OUT$chr=CHR1[1]
    OUT$Zpos=Zpos
    OUT$ZZZ=ZZZ
    OUT$COL=COL 
    OUT$CHR1=CHR1
    OUT$CHR2=CHR2
    OUT$LOC1=LOC1
    OUT$LOC2=LOC2
    OUT$Xall=Xall
    OUT$Zmax=Zmax
    OUT$maxCol.Z=maxCol.Z
    OUT$lineCol=lineCol
    OUT$minCol=minCol
    OUT$maxCol=maxCol
    OUT$exMin=exMin
    OUT$MAIN=MAIN
    ###############
    return(OUT)
    }




inferc.drawContact<-function(OUT){
    OUT=OUT
    ################################
    plot(c(min(OUT$Xall),max(OUT$Xall)),c(0,0),col='white',
        ylim=c(OUT$exMin,OUT$Hh+OUT$Zh+OUT$Ch ),xlim=c(min(OUT$Xall),max(OUT$Xall)),
        yaxt='n',xlab=OUT$CHR1[1],ylab='',main=OUT$MAIN)
        #############
    points(OUT$X[order(OUT$Z)],OUT$Y[order(OUT$Z)]+OUT$Hbase,col=OUT$COL[order(OUT$Z)],pch=18)
    #############
    segments(x0=OUT$LOC1[order(OUT$Z)],y0=OUT$ZZZ[order(OUT$Z)],x1=OUT$LOC2[order(OUT$Z)],y1=OUT$ZZZ[order(OUT$Z)],col=OUT$COL[order(OUT$Z)])
    ##############
    segments(x0=min(OUT$Xall),y0=OUT$Hbase, x1=max(OUT$Xall), y1=OUT$Hbase, col=OUT$lineCol)
    segments(x0=min(OUT$Xall),y0=OUT$Hbase, x1=(max(OUT$Xall)+min(OUT$Xall))/2, y1= OUT$Hbase+OUT$Hh, col=OUT$lineCol)
    segments(x0=max(OUT$Xall),y0=OUT$Hbase, x1=(max(OUT$Xall)+min(OUT$Xall))/2, y1= OUT$Hbase+OUT$Hh, col=OUT$lineCol)
    abline(h=OUT$Hbase-OUT$Zh,lty=2,col=OUT$lineCol)
    ##############
    }



inferc.preV4C<-function(OUT, VIEW_POINT, VPCOL='red',COL1='black',COL2='black',PCH1=3,PCH2=3,CEX1=0.5,CEX2=0.5,Napp=1000,Wapp=2,appCol='grey90'){
    OUT=OUT
    VIEW_POINT=as.numeric(VIEW_POINT)
    PCH1=PCH1
    PCH2=PCH2
    CEX1=CEX1
    CEX2=CEX2
    COL1=COL1
    COL2=COL2
    Napp=Napp
    Wapp=Wapp
    appCol=appCol
    #########
    POS=VIEW_POINT
    ########################
    D=abs(OUT$Xall-POS)
    MIN_INDEX=which(D==min(D))
    MIN_LOC=OUT$Xall[MIN_INDEX][1]
    ##########################
    if(MIN_LOC==min(OUT$Xall)){
        XL=c(max(OUT$Xall))
        ZL=c(0)
       }else if(MIN_LOC==max(OUT$Xall)){
        XL=c(min(OUT$Xall))
        ZL=c(0)
       }else{
        XL=c(min(OUT$Xall),max(OUT$Xall))
        ZL=c(0,0)
       }
    ######################
    ######################
    Xp=c()
    Zp=c()
    Xh=c()
    Zh=c()
    Z_index=c()
    i=1
    while(i<=length(OUT$LOC1)){
        this_loc1=OUT$LOC1[i]
        this_loc2=OUT$LOC2[i]
        this_z=OUT$Zpos[i]
        if(this_loc1==MIN_LOC){
            this_x=this_loc2
            XL=c(XL, this_x)
            ZL=c(ZL, this_z)
            Xh=c(Xh, this_loc2)
            Zh=c(Zh, OUT$Hbase-this_z/OUT$Zmax * OUT$Zh)
            Xp=c(Xp,this_loc2)
            Zp=c(Zp,this_z/OUT$Zmax * OUT$Ch)
            Z_index=c(Z_index,i)
        }else if(this_loc2==MIN_LOC){
            this_x=this_loc1
            XL=c(XL, this_x)
            ZL=c(ZL, this_z)
            Xh=c(Xh, this_loc1)
            Zh=c(Zh, OUT$Hbase-this_z/OUT$Zmax * OUT$Zh)
            Xp=c(Xp,this_loc1)
            Zp=c(Zp,this_z/OUT$Zmax * OUT$Ch)
            Z_index=c(Z_index,i)
            }
        i=i+1
        }    
    #####################
    ZL=ZL[order(XL)]
    XL=XL[order(XL)]
    ZL[which(ZL<0)]=0
    ZL=ZL/OUT$Zmax*OUT$Ch
    ########################
    Napp=1000
    Lapp=c(1:Napp)/Napp * (max(XL)-min(XL))+min(XL)
    Lapp=c(Lapp,XL)
    Lapp=Lapp[order(Lapp)]
    Zapp=approxfun(x=XL,y=ZL,ties='mean')(Lapp)
    #########################
    this_chr=OUT$CHR1[1]
    #########################################
    OUT$v4c_Lapp=Lapp
    OUT$v4c_Zapp=Zapp
    OUT$v4c_XL=XL
    OUT$v4c_ZL=ZL
    OUT$v4c_Xp=Xp
    OUT$v4c_Zp=Zp
    OUT$v4c_Xh=Xh
    OUT$v4c_Zh=Zh
    OUT$v4c_Z_index=Z_index
    OUT$v4c_MIN_LOC=MIN_LOC
    OUT$v4c_VIEW_POINT=POS    
    OUT$v4c_POS=POS    
    OUT$v4c_VPCOL=VPCOL
    OUT$v4c_COL1=COL1
    OUT$v4c_COL2=COL2
    OUT$v4c_PCH1=PCH1
    OUT$v4c_PCH2=PCH2
    OUT$v4c_CEX1=CEX1
    OUT$v4c_CEX2=CEX2
    OUT$v4c_Wapp=Wapp
    OUT$v4c_appCol=appCol
    #########################################
    OUT$v4c_pos=POS
    OUT$v4c_mindist=min(D)
    OUT$v4c_minloc=MIN_LOC
    OUT$v4c_chr=this_chr
    OUT$v4c_x=XL
    OUT$v4c_z=ZL
    return(OUT)
    }


inferc.addV4C<-function(OUT){
    Xtmp=OUT$X[OUT$v4c_Z_index]
    Ytmp=OUT$Y[OUT$v4c_Z_index]
    Ztmp=OUT$Z[OUT$v4c_Z_index]
    COLtmp=OUT$COL[OUT$v4c_Z_index]
    #points(Xtmp[order(Ztmp)],Ytmp[order(Ztmp)]+OUT$Hbase,col=OUT$lineCol,lwd=1, pch=5)
    points(Xtmp[order(Ztmp)],Ytmp[order(Ztmp)]+OUT$Hbase,col='black',lwd=1.5, pch=5)
    points(Xtmp[order(Ztmp)],Ytmp[order(Ztmp)]+OUT$Hbase,col=COLtmp[order(Ztmp)],cex=1.5,pch=18)
    #points(Xtmp[order(Ztmp)],Ytmp[order(Ztmp)]+OUT$Hbase,col='black',cex=1, pch='*')
    #######################
    points(OUT$v4c_Xh,OUT$v4c_Zh,pch=OUT$v4c_PCH1,col=OUT$v4c_COL1,cex=OUT$v4c_CEX1)
    points(OUT$v4c_Lapp,OUT$v4c_Zapp,type='h',lwd=OUT$v4c_Wapp,col=OUT$v4c_appCol)
    points(OUT$v4c_XL,OUT$v4c_ZL,type='b',col='black',cex=0)
    points(OUT$v4c_Xp,OUT$v4c_Zp,pch=OUT$v4c_PCH2,col=OUT$v4c_COL2, cex=OUT$v4c_CEX2)
    segments(x0=OUT$v4c_MIN_LOC,x1=OUT$v4c_MIN_LOC,y0=0,y1=OUT$Ch+OUT$Zh,lty=2,col=OUT$lineCol)
    segments(x0=OUT$v4c_POS,x1=OUT$v4c_POS,y0=0,y1=OUT$Ch+OUT$Zh,lty=1,col=OUT$v4c_VPCOL)
    }



inferc.drawV4C<-function(OUT){
    plot(OUT$v4c_Lapp,OUT$v4c_Zapp,type='h',lwd=OUT$v4c_Wapp,col=OUT$v4c_appCol,ylim=c(0,1),
         yaxt='n',xlab=OUT$CHR1[1],ylab='',main=OUT$MAIN)
    points(OUT$v4c_XL,OUT$v4c_ZL,type='b',col='black',cex=0)
    points(OUT$v4c_Xp,OUT$v4c_Zp,pch=OUT$v4c_PCH2,col=OUT$v4c_COL2, cex=OUT$v4c_CEX2)
    segments(x0=OUT$v4c_MIN_LOC,x1=OUT$v4c_MIN_LOC,y0=0,y1=OUT$Ch,lty=2,col=OUT$lineCol)
    segments(x0=OUT$v4c_POS,x1=OUT$v4c_POS,y0=0,y1=OUT$Ch,lty=1,col=OUT$v4c_VPCOL)
    }















###############


inferc.calCSP<-function(X, Y, r=0, only_pos=FALSE){
    X=X
    Y=Y
    r=r
    only_pos=only_pos
    ############################
    # Ensure positive value
    if( min(X)<0 ){
        X = X - min(inferloop.rmOut(X))
        X[which(X<0)]=0
        }
    if( min(Y)<0 ){
        Y = Y - min(inferloop.rmOut(Y))
        Y[which(Y<0)]=0
        }
    #############################
    r=r
    N=length(X)
    ############################
    X_mean = mean(X)
    Y_mean = mean(Y)
    ###########################
    X_base = X_mean * r
    Y_base = Y_mean * r
    X_delta = X - X_base
    Y_delta = Y - Y_base
    ###########################
    ABCD = inferloop.calABCD(X, Y, X_base, Y_base)
    ###########################
    D = ABCD[['D']]
    D[which(is.na(D))]=0
    ############################
    D_plus = ( ABCD[['A']] + X_delta * Y_delta ) / sqrt( (ABCD[['B']] + X_delta**2) * (ABCD[['C']] + Y_delta**2) )
    ###########################
    M = D_plus - D
    ###########################
    CSP = M 
    CSP[which(is.na(CSP))]=0
    if(only_pos==TRUE){CSP[which(CSP<0)]=0 }
    return( CSP )
    }



inferc.inferLoopSignal<-function(mat, net, r=0,sep='.And.'){
    library(hash)
    r=r # default r=0
    sep=sep
    ####################
    mat=as.matrix(mat)
    net=as.matrix(net)
    tag=rownames(mat)
    net=net[which(net[,1] %in% tag & net[,2] %in% tag),]
    #######################################
    print('hashing...')
    h=hash(keys=tag,values=rep(0,nrow(mat)))
    i=1
    while(i<=nrow(mat)){
        this_tag=tag[i]
        this_v=as.numeric(mat[i,])
        .set(h, keys=this_tag, values=this_v)
        i=i+1}
    #######################################
    print('calculating ILS...')
    out=matrix(0,ncol=ncol(mat),nrow=nrow(net))
    colnames(out)=colnames(mat)
    rownames(out)=paste0(net[,1],sep,net[,2])
    i=1
    while(i<=nrow(net)){
        this_tag1=net[i,1]
        this_tag2=net[i,2]
        #######################
        x=as.vector(hash::values(x=h, keys=this_tag1)[,1])
        y=as.vector(hash::values(x=h, keys=this_tag2)[,1])
        #######################
        z=inferc.calCSP(x,y,r=r)
        #######################
        out[i,]=z
        ########################
        if(i%%10000==1){print(paste0(i,' / ',nrow(net)))}
        i=i+1}
    print('finished!')
    return(out)
    }









