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

cat(columns, file = stderr())

##
pipeline_steps = c(
  "input",
  "X",
  "Extra.QC",
  "Viral.enrichment",
  "Host.depletion",
  "Assembly",
  "Contig.classification",
  "Read.classification",
  "Combined.analysis",
  "Remapping",
  "Request.Mapping",
  "Remap.filtering",
  "Map.filtering",
  "Screening",
  "Reporting",
  "leaves"
)

step_to_colour_map = c(
  "input" = "#FFFFFF",
  "X" = "#000000",
  "Extra.QC" = "#FF0000",
  "Viral.enrichment" = "#FFA500",
  "Host.depletion" = "#FFFF00",
  "Assembly" = "#00FF00",
  "Contig.classification" = "#0000FF",
  "Read.classification" = "#4B0082",
  "Remapping" = "#EE82EE",
  "Request.Mapping"= "#FFC0CB",
  "Remap.filtering" = "#FFC0CB",
  "Map.filtering" = "#FFC0CB",
  "Screening" = "#FFC0CB",
  "Reporting" = "#FFC0CB",
  "Combined.analysis" = "#FFC0CB",
  "leaves" = "#FFFFFF"
)

## create a vector with the colours for each column (pipeline step)
colours = rep("#FFFFFF", length(columns))

for (i in 1:length(columns)) {
  colours[i] = step_to_colour_map[columns[i]]
}

cat("######################################## colors\n", file = stderr())
cat(colours, file = stderr())
cat("\n", file = stderr())
### read the data frame

input = read.table(df_path, header = TRUE, sep = "\t")

write.table(input, file = stderr(), sep = "\t", row.names = FALSE)

cat("\n", file = stderr())
# print column names to stderr

#cat("column names: ", colnames(input), file = stderr())
#cat("\n", file = stderr())

p <- collapsibleTreeNetwork(
  input, 
  attribute= "module",
  fill= "colour",
  tooltipHtml = "module",
  collapsed= FALSE,
  fontSize = 15,
)
### save the dendrogram as html file

htmltools::save_html(p, file = output_path)

