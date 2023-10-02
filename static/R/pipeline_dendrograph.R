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
columns= args[3]
## columns come separated by comma, so we need to split them
columns= strsplit(columns, ",")[[1]]
## remove spaces
columns= gsub(" ", ".", columns)

### read the data frame
input = read.table(df_path, header = TRUE, sep = "\t")
### create the dendrogram

p <- collapsibleTree( input, columns, 
    collapsed= FALSE, fontSize = 15, fill = "lightsteelblue", linkLength= 150
)

### save the dendrogram as html file

htmltools::save_html(p, file = output_path)

