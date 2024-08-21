
############
.generate_mean <- function(exp_sc_mat, TAG, print_step=10){
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


.getChrLoc<-function(TAG, SPLIT1='-', SPLIT2='-'){
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








.colMap<-function(x, value, color){
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


