
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

#        raw         esc    nfeature      ncount
# 0.34695804  0.13992435  0.02597306 -0.01011802



