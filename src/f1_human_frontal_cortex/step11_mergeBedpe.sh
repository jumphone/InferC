
B1=./scHIC/Astro_all_brain.txt_1kb_contacts.cool.fithic/FitHiC.spline_pass1.res1000.significances.txt.gz.0.05.filtered.bedpe
B2=./scHIC/MG_all_brain.txt_1kb_contacts.cool.fithic/FitHiC.spline_pass1.res1000.significances.txt.gz.0.05.filtered.bedpe
B3=./scHIC/Neuron.Ex_all_brain.txt_1kb_contacts.cool.fithic/FitHiC.spline_pass1.res1000.significances.txt.gz.0.05.filtered.bedpe
B4=./scHIC/Neuron.In_all_brain.txt_1kb_contacts.cool.fithic/FitHiC.spline_pass1.res1000.significances.txt.gz.0.05.filtered.bedpe
B5=./scHIC/ODC_all_brain.txt_1kb_contacts.cool.fithic/FitHiC.spline_pass1.res1000.significances.txt.gz.0.05.filtered.bedpe
B6=./scHIC/OPC_all_brain.txt_1kb_contacts.cool.fithic/FitHiC.spline_pass1.res1000.significances.txt.gz.0.05.filtered.bedpe

cat $B1 $B2 $B3 $B4 $B5 $B6 | cut -f 1,2,3,4,5,6 - | sort - | uniq - > ./scHIC/Merged_loop.bedpe 




