library(ACE)
library(QDNAseq)
library(Biobase)
library(future)
library(ggplot2)

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = c("output", "message"))

future::plan("multisession", workers = snakemake@threads)

counts_tsv <- snakemake@input[["counts"]]
binSize <- snakemake@params[["binsize"]]
sample_name <- snakemake@wildcards[["sample"]]
output_plot <- snakemake@output[["plot"]]
output_rds <- snakemake@output[["rds"]]

bins <- getBinAnnotations(binSize = binSize)
fd <- Biobase::pData(bins)

bin_chroms <- as.character(fd$chromosome)
bin_keys <- paste(bin_chroms, as.integer(fd$start), sep = ":")
counts_vec <- integer(nrow(fd))

if (file.exists(counts_tsv) && file.info(counts_tsv)$size > 0) {
    tsv <- read.table(
        counts_tsv,
        header = FALSE,
        sep = "\t",
        col.names = c("chrom", "start", "count"),
        colClasses = c("character", "integer", "integer"),
        stringsAsFactors = FALSE
    )

    tsv_chroms <- tsv$chrom
    bins_have_chr <- any(grepl("^chr", bin_chroms))
    tsv_have_chr <- any(grepl("^chr", tsv_chroms))
    if (!bins_have_chr && tsv_have_chr) {
        tsv_chroms <- sub("^chr", "", tsv_chroms)
        tsv_chroms[tsv_chroms == "M"] <- "MT"
    } else if (bins_have_chr && !tsv_have_chr) {
        tsv_chroms <- ifelse(tsv_chroms == "MT", "chrM", paste0("chr", tsv_chroms))
    }

    tsv_keys <- paste(tsv_chroms, tsv$start, sep = ":")
    idx <- match(tsv_keys, bin_keys)
    valid <- !is.na(idx)
    counts_vec[idx[valid]] <- tsv$count[valid]
}

counts_mat <- matrix(
    counts_vec,
    ncol = 1,
    dimnames = list(rownames(fd), sample_name)
)

total_reads <- sum(as.numeric(counts_vec))
used_bins <- if ("use" %in% colnames(fd)) as.logical(fd$use) else rep(TRUE, nrow(fd))
used_reads <- sum(as.numeric(counts_vec[used_bins]))

pheno <- Biobase::AnnotatedDataFrame(
    data.frame(
        name = sample_name,
        total.reads = total_reads,
        used.reads = used_reads,
        row.names = sample_name,
        stringsAsFactors = FALSE
    )
)

readCounts <- new("QDNAseqReadCounts", bins = bins, counts = counts_mat, phenodata = pheno)
if (exists("expectedVariance", envir = asNamespace("QDNAseq"), inherits = FALSE)) {
    readCounts$expected.variance <- QDNAseq:::expectedVariance(readCounts)
} else {
    mean_used <- mean(counts_vec[used_bins])
    readCounts$expected.variance <- if (mean_used > 0) 1 / mean_used else NA_real_
}

saveRDS(readCounts, file = output_rds)

readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE)
readCountsFiltered <- estimateCorrection(readCountsFiltered)

png(filename = output_plot, width = 3000, height = 1000, res = 300)
plot(readCountsFiltered, ylim = c(10, 20))
dev.off()

sink()
