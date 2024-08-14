
data <- readRDS("./rds/Pancreas_10x_downsampled.rds")
# extract expression data
expression_data <- data$expression_data
# running CytoTRACE 2 main function - cytotrace2 - with default parameters
cytotrace2_result <- CytoTRACE2::cytotrace2(expression_data)



