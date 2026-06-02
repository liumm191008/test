#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
SPECIES="${SPECIES:-mouse}"
TXDB_PATH="${TXDB_PATH:-${WORK_ROOT}/pipeline/database/mm39/TxDb.Mmusculus.GRCm39.115.sqlite}"
GENOME_FASTA="${GENOME_FASTA:-${WORK_ROOT}/pipeline/database/mm39/Mus_musculus.GRCm39.dna.toplevel.fa}"
ENRICHMENT_GENE_ID_TYPE="${ENRICHMENT_GENE_ID_TYPE:-ENTREZID}"
DIFF_PEAK_PVALUE="${DIFF_PEAK_PVALUE:-0.05}"
THREADS="${THREADS:-8}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
BIOCONDUCTOR_IMAGE="${BIOCONDUCTOR_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/bioconductor}"
HOMER_IMAGE="${HOMER_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/homer}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)
SAMPLE_GROUPS=(con con con exo exo exo LPS LPS LPS)

peak_paths=(); bam_paths=()
for sample in "${SAMPLES[@]}"; do
  peak_paths+=("${OUTPUT_DIR}/peaks/${sample}/${sample}.peaks.bed")
  bam_paths+=("${OUTPUT_DIR}/filtered/${sample}/${sample}.filtered.bam")
done
sample_csv=$(IFS=,; echo "${SAMPLES[*]}")
group_csv=$(IFS=,; echo "${SAMPLE_GROUPS[*]}")
peak_csv=$(IFS=,; echo "${peak_paths[*]}")
bam_csv=$(IFS=,; echo "${bam_paths[*]}")
${DOCKER_RUN} "${BIOCONDUCTOR_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/differential_peaks' '${OUTPUT_DIR}/differential_peaks/annotation' '${OUTPUT_DIR}/differential_peaks/enrichment' '${OUTPUT_DIR}/differential_peaks/motif'; Rscript - <<RSCRIPT
suppressPackageStartupMessages(library(DiffBind))
suppressPackageStartupMessages(library(ChIPseeker))
species <- tolower('${SPECIES}')
txdb <- loadDb('${TXDB_PATH}')
outdir <- file.path('${OUTPUT_DIR}', 'differential_peaks')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outdir, 'annotation'), recursive=TRUE, showWarnings=FALSE)
sample_ids <- strsplit('${sample_csv}', ',')[[1]]; groups <- strsplit('${group_csv}', ',')[[1]]; peaks <- strsplit('${peak_csv}', ',')[[1]]; bams <- strsplit('${bam_csv}', ',')[[1]]
expand_sample_groups <- function(ids, groups) {
  if (length(groups) == length(ids)) return(groups)
  unique_groups <- unique(groups[nzchar(groups)])
  if (length(unique_groups) > 0) {
    assigned <- vapply(ids, function(id) {
      hits <- unique_groups[startsWith(tolower(id), tolower(unique_groups))]
      if (length(hits) == 1) hits[1] else NA_character_
    }, character(1))
    if (all(!is.na(assigned))) return(assigned)
    if (length(ids) %% length(unique_groups) == 0) return(rep(unique_groups, each=length(ids) / length(unique_groups)))
  }
  stop(sprintf('groups length (%s) does not match sample_ids length (%s): %s', length(groups), length(ids), paste(groups, collapse=',')))
}
groups <- expand_sample_groups(sample_ids, groups)
if (length(peaks) != length(sample_ids) || length(bams) != length(sample_ids)) stop('sample_ids, peaks, and bams must have the same length')
sheet <- data.frame(SampleID=sample_ids, Tissue=groups, Factor='CutTag', Condition=groups, Replicate=ave(seq_along(sample_ids), groups, FUN=seq_along), bamReads=bams, Peaks=peaks, PeakCaller='bed', stringsAsFactors=FALSE)
write.csv(sheet, file.path(outdir, 'diffbind_sample_sheet.csv'), row.names=FALSE)
safe_plot <- function(path, plot_fun) { tryCatch({ pdf(path); plot_fun(); dev.off() }, error=function(e) { if (dev.cur() != 1) dev.off(); writeLines(conditionMessage(e), paste0(path, '.error.txt')) }) }
write_peak_table <- function(gr, path) { write.table(as.data.frame(gr), path, sep='	', quote=FALSE, row.names=FALSE) }
sanitize_name <- function(value) gsub('[^A-Za-z0-9_.-]+', '_', value)
get_fold_change <- function(report) {
  for (col in c('Fold', 'log2FoldChange', 'Log2FoldChange', 'log2FC', 'logFC')) {
    if (col %in% colnames(report)) return(suppressWarnings(as.numeric(report[[col]])))
  }
  rep(0, nrow(report))
}
write_peak_annotation_subset <- function(report_subset, prefix) {
  gene_path <- file.path(outdir, 'annotation', paste0(prefix, '.enrichment_genes.tsv'))
  if (nrow(report_subset) > 0) {
    bed_path <- file.path(outdir, 'annotation', paste0(prefix, '.bed'))
    write.table(report_subset[, c('seqnames', 'start', 'end')], bed_path, sep='	', quote=FALSE, row.names=FALSE, col.names=FALSE)
    anno <- annotatePeak(bed_path, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
    anno_df <- as.data.frame(anno)
    write.table(anno_df, file.path(outdir, 'annotation', paste0(prefix, '.annotated.tsv')), sep='	', quote=FALSE, row.names=FALSE)
    gene_ids <- unique(as.character(na.omit(anno_df[["geneId"]])))
    write.table(data.frame(gene_id=gene_ids), gene_path, sep='	', quote=FALSE, row.names=FALSE)
  } else {
    write.table(data.frame(gene_id=character()), gene_path, sep='	', quote=FALSE, row.names=FALSE)
  }
}
write_contrast_annotation <- function(report, prefix) {
  write_peak_annotation_subset(report, prefix)
  fold_change <- get_fold_change(report)
  gain_report <- report[!is.na(fold_change) & fold_change > 0, , drop=FALSE]
  loss_report <- report[!is.na(fold_change) & fold_change < 0, , drop=FALSE]
  write.table(gain_report, file.path(outdir, paste0(prefix, '.gain.tsv')), sep='	', quote=FALSE, row.names=FALSE)
  write.table(loss_report, file.path(outdir, paste0(prefix, '.loss.tsv')), sep='	', quote=FALSE, row.names=FALSE)
  write_peak_annotation_subset(gain_report, paste0(prefix, '.gain'))
  write_peak_annotation_subset(loss_report, paste0(prefix, '.loss'))
  write.table(data.frame(category=c('gain', 'loss'), peak_count=c(nrow(gain_report), nrow(loss_report))), file.path(outdir, 'annotation', paste0(prefix, '.gain_loss_summary.tsv')), sep='	', quote=FALSE, row.names=FALSE)
}
dba_obj <- dba(sampleSheet=sheet)
consensus <- dba.peakset(dba_obj, bRetrieve=TRUE)
if (length(consensus) > 0) {
  consensus_df <- data.frame(chrom=as.character(seqnames(consensus)), start=pmax(start(consensus)-1, 0), end=end(consensus), name=paste0('consensus_peak_', seq_along(consensus)), stringsAsFactors=FALSE)
} else {
  consensus_df <- data.frame(chrom=character(), start=integer(), end=integer(), name=character())
}
write.table(consensus_df, file.path(outdir, 'consensus_peaks.bed'), sep='	', quote=FALSE, row.names=FALSE, col.names=FALSE)
dba_obj <- dba.count(dba_obj)
saveRDS(dba_obj, file.path(outdir, 'diffbind.rds'))
write_peak_table(dba.peakset(dba_obj, bRetrieve=TRUE), file.path(outdir, 'counts_matrix.tsv'))
normalized_obj <- dba.count(dba_obj, score=DBA_SCORE_NORMALIZED)
write_peak_table(dba.peakset(normalized_obj, bRetrieve=TRUE), file.path(outdir, 'normalized_matrix.tsv'))
if (length(unique(groups)) >= 2) {
  dba_obj <- dba.contrast(dba_obj, categories=DBA_CONDITION, minMembers=2); dba_obj <- dba.analyze(dba_obj)
  saveRDS(dba_obj, file.path(outdir, 'diffbind.rds'))
  safe_plot(file.path(outdir, 'PCA.pdf'), function() dba.plotPCA(dba_obj, label=DBA_ID))
  safe_plot(file.path(outdir, 'correlation_heatmap.pdf'), function() dba.plotHeatmap(dba_obj))
  contrast_pairs <- combn(unique(groups), 2, simplify=FALSE)
  n_contrasts <- if (!is.null(dba_obj[["contrasts"]])) length(dba_obj[["contrasts"]]) else length(contrast_pairs)
  for (contrast_index in seq_len(n_contrasts)) {
    pair <- if (contrast_index <= length(contrast_pairs)) contrast_pairs[[contrast_index]] else c('groupA', paste0('groupB_', contrast_index))
    contrast_prefix <- paste0('contrast_', sanitize_name(pair[1]), '_vs_', sanitize_name(pair[2]))
    report <- as.data.frame(dba.report(dba_obj, contrast=contrast_index, th=${DIFF_PEAK_PVALUE}))
    write.table(report, file.path(outdir, paste0(contrast_prefix, '.tsv')), sep='	', quote=FALSE, row.names=FALSE)
    write_contrast_annotation(report, contrast_prefix)
    safe_plot(file.path(outdir, paste0(contrast_prefix, '.MA_plot.pdf')), function() dba.plotMA(dba_obj, contrast=contrast_index))
    safe_plot(file.path(outdir, paste0(contrast_prefix, '.volcano_plot.pdf')), function() dba.plotVolcano(dba_obj, contrast=contrast_index))
    safe_plot(file.path(outdir, paste0(contrast_prefix, '.boxplot.pdf')), function() dba.plotBox(dba_obj, contrast=contrast_index))
  }
  report <- as.data.frame(dba.report(dba_obj, contrast=1, th=${DIFF_PEAK_PVALUE}))
  write.table(report, file.path(outdir, 'diffbind_report.tsv'), sep='	', quote=FALSE, row.names=FALSE)
} else { write.table(data.frame(status='skipped', reason='less_than_two_groups'), file.path(outdir, 'diffbind_report.tsv'), sep='	', quote=FALSE, row.names=FALSE) }
RSCRIPT
 test -s '/kegg_data/gene_cluster_enrich.R'; shopt -s nullglob; for gene_file in '${OUTPUT_DIR}/differential_peaks/annotation'/*.enrichment_genes.tsv; do prefix=\$(basename "\${gene_file}" .enrichment_genes.tsv); if awk 'NR > 1 && NF > 0 {found=1} END {exit found ? 0 : 1}' "\${gene_file}"; then Rscript '/kegg_data/gene_cluster_enrich.R' "\${gene_file}" '${ENRICHMENT_GENE_ID_TYPE}' '${SPECIES}' '${OUTPUT_DIR}/differential_peaks/enrichment' "\${prefix}"; else printf 'status\treason\nskipped\tno_genes\n' > '${OUTPUT_DIR}/differential_peaks/enrichment/'"\${prefix}"'.enrichment.skipped.tsv'; fi; done"

${DOCKER_RUN} "${HOMER_IMAGE}" bash -c "set -euo pipefail; shopt -s nullglob; mkdir -p '${OUTPUT_DIR}/differential_peaks/motif'; for bed_file in '${OUTPUT_DIR}/differential_peaks/annotation'/*.gain.bed '${OUTPUT_DIR}/differential_peaks/annotation'/*.loss.bed; do prefix=\$(basename "\${bed_file}" .bed); motif_dir='${OUTPUT_DIR}/differential_peaks/motif/'"\${prefix}"; mkdir -p "\${motif_dir}"; findMotifsGenome.pl "\${bed_file}" '${GENOME_FASTA}' "\${motif_dir}" -size given -len 6,8,10,12 -p ${THREADS}; done"
