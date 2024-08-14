
bedtools pairtopair -type both -a ./scHIC/Merged_loop_within500k.bedpe -b  ./rds/ciceroFrame_conns_uniq.bedpe >  ./scHIC/Merged_loop_within500k_overPeakPair_allinfo.bedpe


bedtools pairtopair -type both -a ./scHIC/Merged_loop_within500k.bedpe -b  ./rds/ciceroFrame_conns_uniq.bedpe | cut -f 1,2,3,4,5,6 - | uniq -> ./scHIC/Merged_loop_within500k_overPeakPair_uniq.bedpe

