
cat ../data/GSE184386/*/peaks.bed | bedtools sort -i - | bedtools merge -i - > ../data/GSE184386/ALL_PEAK.bed


