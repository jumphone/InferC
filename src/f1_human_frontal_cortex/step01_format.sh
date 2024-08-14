#FRAG_HG38=/home/disk/database/data/scHIC/GSE130711_analysis/scATAC/GSM5589415_UMB4540_snATAC_frontal_cortex_rep2_fragments.bed.gz
#FRAG_HG19=/home/disk/database/data/scHIC/GSE130711_analysis/scATAC/GSM5589415_UMB4540_snATAC_frontal_cortex_rep2_fragments.hg19.bed
FRAG_HG38=/home/database/data/InferX/data/frontal_cortex/scATAC/GSM5589415_UMB4540_snATAC_frontal_cortex_rep2_fragments.bed.gz
FRAG_HG19=/home/database/data/InferX/data/frontal_cortex/scATAC/GSM5589415_UMB4540_snATAC_frontal_cortex_rep2_fragments.hg19.bed

CHAIN=/home/database/reference/hg38/hg38ToHg19.over.chain
liftOver $FRAG_HG38 $CHAIN $FRAG_HG19 $FRAG_HG38\.NAtoHg19.bed
gzip $FRAG_HG19
INPUT=$FRAG_HG19\.gz
OUTPUT=$INPUT\.format
mkdir $OUTPUT
gunzip -c $INPUT | cut -f 1,2,3,4,5 -  | bedtools sort -i - > $OUTPUT/sorted.bed
bgzip -c $OUTPUT/sorted.bed > $OUTPUT/sorted.bed.gz
rm -rf $OUTPUT/sorted.bed
tabix -p bed $OUTPUT/sorted.bed.gz
