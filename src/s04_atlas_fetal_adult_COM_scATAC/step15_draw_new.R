

AUC=readRDS('../rds/Eva_AUC.rds')
COR=readRDS('../rds/Eva_COR.rds')


library('ComplexHeatmap')
library('circlize')

MAT=COR[,c(1,3:5)]
MAT=MAT[order(-MAT[,1]),]

o.mat=MAT
col_fun=colorRamp2(
            c(-0.2, 0, 0.2,0.4),
             c('skyblue1',
              'white','indianred1','red3'))
ht1=Heatmap(o.mat,row_title='',name="v",
        cluster_rows=F, cluster_columns=F,
        #show_column_dend = T, show_row_dend = T,
        show_column_names= T, show_row_names= T,
        col=col_fun, border = TRUE,
        row_names_side='right',
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.3f", o.mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })


ht2=Heatmap(o.mat,row_title='',name="v",
        cluster_rows=F, cluster_columns=F,
        #show_column_dend = T, show_row_dend = T,
        show_column_names= T, show_row_names= T,
        col=col_fun, border = TRUE,
        row_names_side='right',
        row_names_gp = gpar(fontsize = 10, lineheight=NULL),
        column_names_gp = gpar(fontsize = 10, lineheight=NULL),
        row_names_max_width = max_text_width(rownames(o.mat),gp = gpar(fontsize = 15)),
        column_names_max_height = max_text_width(colnames(o.mat),gp = gpar(fontsize = 15)),
        )

pdf('../plot/s15p01_heatmap.pdf',height=8,width=5)
ht1
ht2
dev.off()

apply(o.mat,2,mean)
#     inferc         esc    nfeature      ncount
# 0.38650499  0.13992435  0.02598008 -0.01025985


ORGAN=c('Muscle','Adrenal','Lung','Liver','Stomach','Pancreas')
ADULT=c(60391,10537,41089,10557,29255,43493)
FETAL=c(52886,101435,145131,293375,5797,7144)
RATIO=cbind(ADULT, FETAL)
rownames(RATIO)=ORGAN

AUC=readRDS('../rds/Eva_AUC.rds')
COR=readRDS('../rds/Eva_COR.rds')

RATIO=RATIO[order(-COR[,2]),]

pdf('../plot/s15p02_barplot.pdf',height=5,width=5)
barplot(t(RATIO/rowSums(RATIO)),col=c('skyblue1','indianred1'))
dev.off()

#              ADULT     FETAL
#Muscle   0.53312676 0.4668732
#Pancreas 0.85891739 0.1410826
#Stomach  0.83461714 0.1653829
#Adrenal  0.09410388 0.9058961
#Lung     0.22064762 0.7793524
#Liver    0.03473474 0.9652653

#        ADULT  FETAL
#Muscle   60391  52886
#Pancreas 43493   7144
#Stomach  29255   5797
#Adrenal  10537 101435
#Lung     41089 145131
#Liver    10557 293375

barplot(RATIO[6,])











