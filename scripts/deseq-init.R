#log <- file(snakemake@log[[1]], open="wt")
#sink(log, type="message")

library("DESeq2")

parallel <- FALSE
#if (snakemake@threads > 1) {
    #library("BiocParallel")
    # setup parallelization
    #register(MulticoreParam(snakemake@threads))
    #parallel <- TRUE
#}

args = commandArgs(trailingOnly=TRUE)
output=paste(as.character(args[3]), "/", sep='')

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
#cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene", check.names=FALSE)
cts= read.table(args[1],header=TRUE, row.names="gene_id", check.names=FALSE)
#coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample", check.names=FALSE)
coldata <- read.table(args[2], header=TRUE, row.names="sample", check.names=FALSE)

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=~ condition)

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]

#normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)
saveRDS(dds, file=paste(output,"all.rds",sep=""))
