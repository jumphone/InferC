LIST_COM=read.table('COM_ID.txt')
LIST_ADULT=read.table('LST_ADULT_ID.txt',row.names=1)
LIST_FETAL=read.table('LST_FETAL_ID.txt',row.names=1)

options(scipen = 999999999)

i=1
while(i<=nrow(LIST_COM)){
    this_organ=LIST_COM[i,1]
    this_adult=stringr::str_split(LIST_COM[i,2],',')[[1]]
    this_fetal=stringr::str_split(LIST_COM[i,3],',')[[1]]
    this_adult_path_list=LIST_ADULT[this_adult,1]
    this_fetal_path_list=LIST_FETAL[this_fetal,1]
    this_all_path_list=c(this_adult_path_list,this_fetal_path_list)
    print(this_all_path_list)
    j=1
    while(j<=length(this_all_path_list)){
        this_path=this_all_path_list[j]
        this_peak_file=paste0(this_path,'/','peaks_macs2.rds')
        this_peak=readRDS(this_peak_file)
        this_chr=as.character(this_peak@seqnames)
        this_start=this_peak@ranges@start
        this_end=this_peak@ranges@start+this_peak@ranges@width-1
        BED=cbind(this_chr,this_start,this_end)
        this_bed_file=paste0(this_path,'/','peaks_macs2.rds.bed')
        write.table(BED, this_bed_file,quote=F,sep='\t',row.names=F,col.names=F)
        j=j+1}
        print(i)
    i=i+1
    }

