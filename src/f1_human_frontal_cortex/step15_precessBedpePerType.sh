
scHIC_BEDPE=./scHIC/Astro_all_brain.txt_1kb_contacts.cool.fithic/FitHiC.spline_pass1.res1000.significances.txt.gz.0.05.filtered.bedpe
bedtools pairtopair -type both -a $scHIC_BEDPE -b  ./rds/ciceroFrame_conns_uniq.bedpe >  $scHIC_BEDPE\.overPeakPair_uniq.bedpe

scHIC_BEDPE=./scHIC/MG_all_brain.txt_1kb_contacts.cool.fithic/FitHiC.spline_pass1.res1000.significances.txt.gz.0.05.filtered.bedpe
bedtools pairtopair -type both -a $scHIC_BEDPE -b  ./rds/ciceroFrame_conns_uniq.bedpe >  $scHIC_BEDPE\.overPeakPair_uniq.bedpe

scHIC_BEDPE=./scHIC/Neuron.Ex_all_brain.txt_1kb_contacts.cool.fithic/FitHiC.spline_pass1.res1000.significances.txt.gz.0.05.filtered.bedpe
bedtools pairtopair -type both -a $scHIC_BEDPE -b  ./rds/ciceroFrame_conns_uniq.bedpe >  $scHIC_BEDPE\.overPeakPair_uniq.bedpe

scHIC_BEDPE=./scHIC/Neuron.In_all_brain.txt_1kb_contacts.cool.fithic/FitHiC.spline_pass1.res1000.significances.txt.gz.0.05.filtered.bedpe
bedtools pairtopair -type both -a $scHIC_BEDPE -b  ./rds/ciceroFrame_conns_uniq.bedpe >  $scHIC_BEDPE\.overPeakPair_uniq.bedpe

scHIC_BEDPE=./scHIC/ODC_all_brain.txt_1kb_contacts.cool.fithic/FitHiC.spline_pass1.res1000.significances.txt.gz.0.05.filtered.bedpe
bedtools pairtopair -type both -a $scHIC_BEDPE -b  ./rds/ciceroFrame_conns_uniq.bedpe >  $scHIC_BEDPE\.overPeakPair_uniq.bedpe

scHIC_BEDPE=./scHIC/OPC_all_brain.txt_1kb_contacts.cool.fithic/FitHiC.spline_pass1.res1000.significances.txt.gz.0.05.filtered.bedpe
bedtools pairtopair -type both -a $scHIC_BEDPE -b  ./rds/ciceroFrame_conns_uniq.bedpe >  $scHIC_BEDPE\.overPeakPair_uniq.bedpe
