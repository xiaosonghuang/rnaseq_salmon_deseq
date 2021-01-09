## Gene-level differential expression analysis using DESeq2
## Setup
### Bioconductor and CRAN libraries used
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
#library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
## List all directories containing data
samples <- list.files(path = "./data", full.names = T, pattern = "salmon$")
## Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")
## Since all quant files have the same name it is useful to have names for each element
names(files) <- str_replace(samples, "./data/", "") %>% str_replace(".salmon", "")
# Load the annotation table for GrCh38
tx2gene <- read.delim("tx2gene_grch38_ens94.txt")
# Take a look at it
tx2gene %>% View()
# Run tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("tx_id", "ensgene")], countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)
# Look at the counts
txi$counts %>% View()
# Write the counts to an object
data <- txi$counts %>% round() %>% data.frame()
## Create a sampletable/metadata
sampletype <- factor(c(rep("control",3), rep("E2F2_knockdown", 3)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))
# plot a histogram of the counts for a single sample
ggplot(data) +
  geom_histogram(aes(x = NC_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")
mean_counts <- apply(data[,1:3], 1, mean)
variance_counts <- apply(data[,1:3], 1, var)
df <- data.frame(mean_counts, variance_counts)
ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) +
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color = "red")

### Check that sample names match in both files
all(colnames(txi$counts) %in% rownames(meta))
all(colnames(txi$counts) == rownames(meta))
