#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GenomeInfoDb))

usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript peak_annotation.R --peak <macs3_peak> --sample <sample_id> --species <mouse|human> --annotation-dir <dir>",
    "",
    "Required arguments:",
    "  --peak             MACS3 narrowPeak, broadPeak, or BED-like peak file.",
    "  --sample           Sample ID used as output filename prefix.",
    "  --species          Species name: mouse/mm/mmu or human/hs/hsa.",
    "  --annotation-dir   Directory for annotated peak, summary, enrichment-gene-list tables, and plots.",
    "",
    "Optional arguments:",
    "  --tss-upstream     TSS upstream window. Default: 3000.",
    "  --tss-downstream   TSS downstream window. Default: 3000.",
    sep = "\n"
  ))
}

parse_args <- function(argv) {
  args <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected positional argument: %s", key))
    }
    name <- substring(key, 3)
    if (name %in% c("help", "h")) {
      usage()
      quit(save = "no", status = 0)
    }
    if (i == length(argv) || startsWith(argv[[i + 1]], "--")) {
      args[[name]] <- TRUE
      i <- i + 1
    } else {
      args[[name]] <- argv[[i + 1]]
      i <- i + 2
    }
  }
  args
}

get_arg <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(as.character(value))) {
    if (required) {
      stop(sprintf("Missing required argument --%s", name))
    }
    return(default)
  }
  as.character(value)
}

species_txdb <- function(species) {
  normalized <- tolower(species)
  if (normalized %in% c("mouse", "mus_musculus", "mm", "mmu")) {
    return(list(species = "mouse", txdb_package = "TxDb.Mmusculus.UCSC.mm39.knownGene"))
  }
  if (normalized %in% c("human", "homo_sapiens", "hs", "hsa")) {
    return(list(species = "human", txdb_package = "TxDb.Hsapiens.UCSC.hg38.knownGene"))
  }
  stop("--species must be mouse/mm/mmu or human/hs/hsa")
}

load_species_txdb <- function(txdb_package) {
  suppressWarnings(suppressPackageStartupMessages(require(txdb_package, character.only = TRUE)))
  get(txdb_package)
}

write_table <- function(data, path) {
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

normalize_peak_seqnames <- function(peak_path, txdb) {
  txdb_seqlevels <- seqlevels(txdb)
  txdb_has_chr <- any(startsWith(txdb_seqlevels, "chr"))
  peaks <- read.delim(peak_path, header = FALSE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(peaks) == 0 || ncol(peaks) < 3) {
    stop(sprintf("Peak file must contain at least 3 BED columns: %s", peak_path))
  }
  peak_has_chr <- any(startsWith(as.character(peaks[[1]]), "chr"))
  if (txdb_has_chr && !peak_has_chr) {
    peaks[[1]] <- ifelse(peaks[[1]] %in% c("M", "MT", "Mt", "mt"), "chrM", paste0("chr", peaks[[1]]))
  } else if (!txdb_has_chr && peak_has_chr) {
    peaks[[1]] <- sub("^chr", "", as.character(peaks[[1]]))
    peaks[[1]] <- ifelse(peaks[[1]] == "M", "MT", peaks[[1]])
  } else {
    return(peak_path)
  }
  normalized_peak_path <- tempfile(pattern = "peak_annotation.", fileext = ".bed")
  write.table(peaks, normalized_peak_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  normalized_peak_path
}

argv <- commandArgs(trailingOnly = TRUE)
args <- parse_args(argv)

peak_path <- get_arg(args, "peak", required = TRUE)
sample_id <- get_arg(args, "sample", required = TRUE)
species_config <- species_txdb(get_arg(args, "species", required = TRUE))
species <- species_config$species
txdb_package <- species_config$txdb_package
annotation_dir <- get_arg(args, "annotation-dir", required = TRUE)
plot_dir <- annotation_dir
tss_upstream <- as.integer(get_arg(args, "tss-upstream", "3000"))
tss_downstream <- as.integer(get_arg(args, "tss-downstream", "3000"))

if (!file.exists(peak_path) || file.info(peak_path)$size == 0) {
  stop(sprintf("Peak file does not exist or is empty: %s", peak_path))
}
dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

txdb <- load_species_txdb(txdb_package)
normalized_peak_path <- normalize_peak_seqnames(peak_path, txdb)
peak_anno <- annotatePeak(normalized_peak_path, TxDb = get(txdb_package), tssRegion = c(-tss_upstream, tss_downstream), verbose = FALSE)
peak_df <- as.data.frame(peak_anno)
write_table(peak_df, file.path(annotation_dir, paste0(sample_id, ".peaks.annotated.tsv")))

annotation_summary <- as.data.frame(table(peak_df$annotation), stringsAsFactors = FALSE)
colnames(annotation_summary) <- c("annotation", "count")
annotation_summary$ratio <- if (sum(annotation_summary$count) > 0) annotation_summary$count / sum(annotation_summary$count) else 0
write_table(annotation_summary, file.path(annotation_dir, paste0(sample_id, ".annotation_summary.tsv")))

gene_ids <- unique(as.character(na.omit(peak_df$geneId)))
gene_ids <- gene_ids[gene_ids != ""]
enrichment_genes <- data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)
write_table(enrichment_genes, file.path(annotation_dir, paste0(sample_id, ".peak_genes.txt")))
write_table(enrichment_genes, file.path(annotation_dir, paste0(sample_id, ".enrichment_genes.tsv")))

peak_stats <- data.frame(
  sample_id = sample_id,
  species = species,
  peak_file = peak_path,
  txdb_package = txdb_package,
  normalized_peak_file = normalized_peak_path,
  total_peaks = nrow(peak_df),
  enrichment_gene_count = nrow(enrichment_genes),
  stringsAsFactors = FALSE
)
write_table(peak_stats, file.path(annotation_dir, paste0(sample_id, ".peak_annotation_stats.tsv")))

pdf(file.path(plot_dir, paste0(sample_id, ".annotation_pie.pdf")), width = 8, height = 8)
print(plotAnnoPie(peak_anno))
dev.off()

pdf(file.path(plot_dir, paste0(sample_id, ".tss_distance.pdf")), width = 8, height = 6)
print(plotDistToTSS(peak_anno, title = paste0(sample_id, " peak distance to TSS")))
dev.off()

pdf(file.path(plot_dir, paste0(sample_id, ".genomic_feature_distribution.pdf")), width = 9, height = 6)
print(plotAnnoBar(peak_anno) + ggtitle(paste0(sample_id, " genomic feature distribution")))
dev.off()
