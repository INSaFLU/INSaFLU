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


##
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

cat(colours, file = stderr())
cat("\n", file = stderr())
### read the data frame
input = read.table(df_path, header = TRUE, sep = "\t")
# replace empty values with NA
input[input == ""] <- NA
### prep the fill vector. for each column, repeat the colour of the pipeline step by the number of unique elements in that column
# print input shape to stderr
cat("input shape: ", dim(input), file = stderr())

# jump line
cat("\n", file = stderr())
# print input table to stderr
write.table(input, file = stderr(), sep = "\t", row.names = FALSE)

cat("\n", file = stderr())
cat(columns, file = stderr())
cat("\n", file = stderr())

fill = c(
  "seashell"
)
cat("\n", file = stderr())

# for each column (pipeline step)
annotations = c("input")

for (i in 1:length(columns)) {
  # get the unique values in that column
  unique_values = unique(input[, columns[i]])
  # remove NA  values
  unique_values = unique_values[!is.na(unique_values)]
  # print column
  cat(columns[i], file = stderr())
  cat("\n", file = stderr())
  # print unique values to stderr
  cat("unique values: ", unique_values, file = stderr())
  # jump line
  cat("\n", file = stderr())
  cat(length(unique_values))
  cat("\n", file = stderr())

  # for each unique value
  fill = c(fill, rep(colours[i], length(unique_values)))
  annotations = c(annotations, rep(columns[i], length(unique_values)))
}
cat(fill, file = stderr())
cat("\n", file = stderr())

### create the dendrogram

p <- collapsibleTree(input,  
    collapsed= FALSE, 
    fontSize = 15, 
    fill = fill,
    linkLength= 150,
    hierarchy= columns,
    tooltip = FALSE)

### save the dendrogram as html file

species <- read.csv("https://apps.fs.usda.gov/fia/datamart/CSV/REF_SPECIES_GROUP.csv")

fill <- c(
    # The root
    "seashell",
    # Unique regions
    rep("brown", length(unique(species$REGION))),
    # Unique classes per region
    rep("khaki", length(unique(paste(species$REGION, species$CLASS)))),
    # Unique names per region
    rep("forestgreen", length(unique(paste(species$NAME, species$REGION))))
  )

#cat(fill, file = stderr())
#cat("\n", file = stderr())
p <- collapsibleTree(
  species,
  hierarchy = c("REGION", "CLASS", "NAME"),
  fill = fill
)

htmltools::save_html(p, file = output_path)

