# SalmonTElocus differential expression (called by `SalmonTElocus.py test`).
# Args: counts_table condition_file out_dir level conditions min_count
#   conditions = "ref,other" for a two-group comparison, or "" when the
#   condition column is numeric (continuous covariate, e.g. age).
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) stop("usage: SalmonTElocus_DE.R table conditions out level conditions min_count")
table_path <- args[1]; cond_path <- args[2]; out_dir <- args[3]; level <- args[4]
cond_levels <- if (nzchar(args[5])) strsplit(args[5], ",", fixed = TRUE)[[1]] else NULL
min_count <- as.numeric(args[6])

annot_cols <- if (level == "locus") {
  c("locus_id", "chrom", "start", "end", "strand", "transcript_id", "subfamily", "family", "class", "length",
    "gene_context", "overlapping_genes", "gene_types", "gene_strand")  # last four only with --gene-gtf
} else {
  c("subfamily", "family", "class")
}

message("Reading ", table_path)
tab <- read.csv(table_path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "")
annot_cols <- intersect(annot_cols, colnames(tab))
annot <- tab[, annot_cols, drop = FALSE]
ids <- if (level == "locus") annot$locus_id else paste(annot$subfamily, annot$family, annot$class, sep = ":")
counts <- as.matrix(tab[, setdiff(colnames(tab), annot_cols), drop = FALSE])
rownames(counts) <- ids; rownames(annot) <- ids

cond <- read.csv(cond_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample", "condition") %in% colnames(cond))) stop(cond_path, " needs columns sample,condition")
unknown <- setdiff(cond$sample, colnames(counts))
if (length(unknown)) stop("Samples in ", cond_path, " not in the count table: ", paste(unknown, collapse = ", "))
cond <- cond[!is.na(cond$condition) & cond$condition != "NA" & cond$condition != "", ]

if (!is.null(cond_levels)) {
  if (length(cond_levels) != 2) stop("--conditions needs exactly two labels, reference first")
  cond <- cond[cond$condition %in% cond_levels, ]
  n <- table(factor(cond$condition, levels = cond_levels))
  if (any(n < 2)) stop("Need >= 2 samples per condition; have: ", paste(names(n), n, sep = "=", collapse = ", "))
  cond$condition <- factor(cond$condition, levels = cond_levels)
  message(sprintf("Comparing %s (n=%d) vs reference %s (n=%d)", cond_levels[2], n[2], cond_levels[1], n[1]))
} else {
  num <- suppressWarnings(as.numeric(cond$condition))
  if (anyNA(num)) stop("Condition column is not numeric; pass --conditions ref,other for a group comparison")
  if (length(num) < 3) stop("Need >= 3 samples for a numeric covariate")
  cond$condition <- num
  message(sprintf("Testing association with a numeric covariate (n=%d samples)", length(num)))
}

counts <- round(counts[, cond$sample, drop = FALSE])  # Salmon counts are fractional
storage.mode(counts) <- "integer"
keep <- rowSums(counts) >= min_count
message(sprintf("%d of %d %s features have >= %g reads", sum(keep), length(keep), level, min_count))
if (sum(keep) < 2) stop("Too few features pass --min-count to test")
counts <- counts[keep, , drop = FALSE]

col_data <- data.frame(condition = cond$condition, row.names = cond$sample)
dds <- DESeqDataSetFromMatrix(counts, col_data, design = ~ condition)
dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) {
  # Too few features to fit a dispersion trend (typical for small or
  # low-coverage data): use gene-wise estimates, as DESeq2 itself advises.
  message("Dispersion trend fit failed; using gene-wise dispersion estimates (", conditionMessage(e), ")")
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  nbinomWaldTest(dds)
})
res <- results(dds)

out <- cbind(annot[rownames(res), , drop = FALSE], as.data.frame(res))
out <- out[order(out$padj, out$pvalue, na.last = TRUE), ]
write.csv(out, file.path(out_dir, sprintf("results_%s.csv", level)), row.names = FALSE)

norm <- cbind(annot[rownames(dds), , drop = FALSE], counts(dds, normalized = TRUE))
write.csv(norm, file.path(out_dir, sprintf("normalized_counts_%s.csv", level)), row.names = FALSE)

sig <- !is.na(out$padj) & out$padj < 0.05
by_class <- do.call(rbind, lapply(split(seq_len(nrow(out)), out$class), function(i) {
  lfc <- out$log2FoldChange[i]
  data.frame(class = out$class[i[1]], n_tested = length(i), n_padj_0.05 = sum(sig[i]),
             n_up = sum(sig[i] & lfc > 0, na.rm = TRUE), n_down = sum(sig[i] & lfc < 0, na.rm = TRUE),
             mean_log2FoldChange = mean(lfc, na.rm = TRUE))
}))
write.csv(by_class, file.path(out_dir, sprintf("summary_by_class_%s.csv", level)), row.names = FALSE)

if ("gene_context" %in% colnames(out)) {
  # same summary split by genomic context: intergenic hits are the best
  # candidates for TE-driven (rather than host-gene-driven) expression
  key <- paste(out$class, out$gene_context, sep = "|")
  by_ctx <- do.call(rbind, lapply(split(seq_len(nrow(out)), key), function(i) {
    lfc <- out$log2FoldChange[i]
    data.frame(class = out$class[i[1]], gene_context = out$gene_context[i[1]], n_tested = length(i),
               n_padj_0.05 = sum(sig[i]), n_up = sum(sig[i] & lfc > 0, na.rm = TRUE),
               n_down = sum(sig[i] & lfc < 0, na.rm = TRUE), mean_log2FoldChange = mean(lfc, na.rm = TRUE))
  }))
  write.csv(by_ctx, file.path(out_dir, "summary_by_class_and_context_locus.csv"), row.names = FALSE)
}

pdf(file.path(out_dir, sprintf("ma_plot_%s.pdf", level)))
plotMA(res, main = sprintf("SalmonTElocus: %s level", level))
invisible(dev.off())

message(sprintf("%d %s features with padj < 0.05. Results in %s", sum(sig), level, out_dir))
