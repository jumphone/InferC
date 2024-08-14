LST=read.table('LST.txt')



options(scipen = 999999999)

    this_all_path_list=LST[,1]
    print(this_all_path_list)
    j=1
    while(j<=length(this_all_path_list)){
        this_path=paste0(this_all_path_list[j],'.format')
        this_peak_file=paste0(this_path,'/','peaks_macs2.rds')
        this_peak=readRDS(this_peak_file)
        this_chr=as.character(this_peak@seqnames)
        this_start=this_peak@ranges@start
        this_end=this_peak@ranges@start+this_peak@ranges@width-1
        BED=cbind(this_chr,this_start,this_end)
        this_bed_file=paste0(this_path,'/','peaks_macs2.rds.bed')
        write.table(BED, this_bed_file,quote=F,sep='\t',row.names=F,col.names=F)
        j=j+1}




