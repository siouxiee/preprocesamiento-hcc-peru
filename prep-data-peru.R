library(affy)
library(GEOquery)
library(tidyverse)
library(oligo)
library(R.utils)
require(hta20transcriptcluster.db)

options(timeout = 1600)
getGEOSuppFiles("GSE136247",baseDir="D:/d o c u/PUCP/TESIS")

untar("D:/d o c u/PUCP/TESIS/GSE136247/GSE136247_RAW.tar",exdir = 'D:/d o c u/PUCP/TESIS/GSE136247/data/')

# Lists .gz files on the directory
gz_files <- list.files("D:/d o c u/PUCP/TESIS/GSE136247/data/", pattern = "\\.CEL.gz$", full.names = TRUE)

# unzip files
sapply(gz_files, gunzip, remove = FALSE)  # 'remove = FALSE' para mantener los archivos originales comprimidos

celfiles <- list.files("D:/d o c u/PUCP/TESIS/GSE136247/data/", pattern = "\\.CEL$", full.names = TRUE)

raw.data <- read.celfiles(celfiles)

# performing RMA normalization
normalized.data <- rma(raw.data)
#normalized.data <- rma(raw.data, normalize = FALSE)

# get expression estimates
normalized.expr <- as.data.frame(exprs(normalized.data))

# map probe IDs to gene symbols
gse <- getGEO("GSE136247", GSEMatrix = TRUE)

# fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE136247_series_matrix.txt.gz@featureData@data

probe_ids <- rownames(normalized.expr)

mapping <- mapIds(hta20transcriptcluster.db,keys = probe_ids,column = 'SYMBOL',keytype = 'PROBEID')

head(mapping)

head(mapping[!is.na(mapping)])

normalized.expr$GeneSymbol <- mapping[rownames(normalized.expr)]

# Remove rows without gene symbol mappings
normalized.expr <- normalized.expr[!is.na(normalized.expr$GeneSymbol), ]

# Save the mapped and normalized expression data
write.csv(normalized.expr, "D:/d o c u/PUCP/TESIS/expression_mapped.csv", row.names = TRUE)

normalized.expr