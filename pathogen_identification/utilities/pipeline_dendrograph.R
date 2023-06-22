#!/usr/bin/env Rscript

library(collapsibleTree) 
library(ellipsis)

# data frame providade by the user as argument
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

df_path= args[1]
output_path= args[2]

### read the data frame
df = read.table(df_path, header = TRUE, sep = "\t")

### create the dendrogram

p <- collapsibleTree( df, c(
    "QC", "Assembly", "Contig.classification", "Read.classification", "Remapping"), collapsed= FALSE
)

### save the dendrogram as html file

htmltools::save_html(p, file = output_path) 