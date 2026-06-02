version 1.0

## Paired-end RNA-seq sample definition.
## All paths must be absolute host paths available through the docker_run mount.
struct RnaSeqSample {
  String sample_id
  String group
  String read1_path
  String read2_path
}

workflow RnaSeq {
  input {
    Array[RnaSeqSample] samples
    String genome_fasta_path = "/home/data/vip01/work/pipeline/database/mm39/Mus_musculus.GRCm39.dna.toplevel.fa"
    String annotation_gtf_path = "/home/data/vip01/work/pipeline/database/mm39/Mus_musculus.GRCm39.115.gtf"
    String star_index_path = ""
    String output_dir = "/home/data/vip01/work/rnaseq_results"
    String enrichment_annotation_path = ""
    Int threads = 8
    Int sjdb_overhang = 149
    Int min_count = 10
    Float padj_cutoff = 0.05
    Float log2fc_cutoff = 1.0

    String docker_run = "docker run --rm --security-opt seccomp=unconfined  -v /home/data/vip01/work:/home/data/vip01/work"
    String fastqc_image = "registry.cn-guangzhou.aliyuncs.com/origen/fastqc"
    String trim_galore_image = "registry.cn-guangzhou.aliyuncs.com/origen/trim-galore"
    String star_image = "registry.cn-guangzhou.aliyuncs.com/origen/star"
    String subread_image = "registry.cn-guangzhou.aliyuncs.com/origen/subread"
    String multiqc_image = "registry.cn-guangzhou.aliyuncs.com/origen/multiqc"
    String differential_analysis_image = "registry.cn-guangzhou.aliyuncs.com/origen/deseq2"
  }

  if (star_index_path == "") {
    call BuildStarIndex {
      input:
        genome_fasta_path = genome_fasta_path,
        annotation_gtf_path = annotation_gtf_path,
        output_dir = output_dir,
        threads = threads,
        sjdb_overhang = sjdb_overhang,
        docker_run = docker_run,
        image = star_image
    }
  }

  String resolved_star_index_path = if star_index_path != "" then star_index_path else select_first([BuildStarIndex.genome_dir])

  scatter (sample in samples) {
    String current_sample_id = sample.sample_id
    String current_sample_group = sample.group

    call FastQc as RawFastQc {
      input:
        sample_id = sample.sample_id,
        read1_path = sample.read1_path,
        read2_path = sample.read2_path,
        output_dir = output_dir,
        threads = threads,
        docker_run = docker_run,
        image = fastqc_image
    }

    call TrimReads {
      input:
        sample_id = sample.sample_id,
        read1_path = sample.read1_path,
        read2_path = sample.read2_path,
        output_dir = output_dir,
        threads = threads,
        docker_run = docker_run,
        image = trim_galore_image
    }

    call FastQc as TrimmedFastQc {
      input:
        sample_id = sample.sample_id + ".trimmed",
        read1_path = TrimReads.trimmed_read1_path,
        read2_path = TrimReads.trimmed_read2_path,
        output_dir = output_dir,
        threads = threads,
        docker_run = docker_run,
        image = fastqc_image
    }

    call StarAlign {
      input:
        sample_id = sample.sample_id,
        read1_path = TrimReads.trimmed_read1_path,
        read2_path = TrimReads.trimmed_read2_path,
        genome_dir = resolved_star_index_path,
        output_dir = output_dir,
        threads = threads,
        docker_run = docker_run,
        image = star_image
    }
  }

  call FeatureCounts {
    input:
      bam_paths = StarAlign.sorted_bam_path,
      annotation_gtf_path = annotation_gtf_path,
      output_dir = output_dir,
      threads = threads,
      docker_run = docker_run,
      image = subread_image
  }

  call DifferentialExpressionEnrichment {
    input:
      count_matrix_path = FeatureCounts.count_matrix_path,
      sample_ids = current_sample_id,
      sample_groups = current_sample_group,
      output_dir = output_dir,
      enrichment_annotation_path = enrichment_annotation_path,
      min_count = min_count,
      padj_cutoff = padj_cutoff,
      log2fc_cutoff = log2fc_cutoff,
      docker_run = docker_run,
      image = differential_analysis_image
  }

  call MultiQc {
    input:
      output_dir = output_dir,
      count_matrix_path = FeatureCounts.count_matrix_path,
      docker_run = docker_run,
      image = multiqc_image
  }

  output {
    String star_index_dir = resolved_star_index_path
    Array[String] raw_fastqc_dirs = RawFastQc.fastqc_dir_path
    Array[String] trimmed_fastqc_dirs = TrimmedFastQc.fastqc_dir_path
    Array[String] trimmed_read1 = TrimReads.trimmed_read1_path
    Array[String] trimmed_read2 = TrimReads.trimmed_read2_path
    Array[String] sorted_bams = StarAlign.sorted_bam_path
    Array[String] alignment_logs = StarAlign.final_log_path
    String count_matrix = FeatureCounts.count_matrix_path
    String count_summary = FeatureCounts.summary_path
    String sample_metadata = DifferentialExpressionEnrichment.sample_metadata_path
    String de_summary = DifferentialExpressionEnrichment.de_summary_path
    String de_results_dir = DifferentialExpressionEnrichment.de_results_dir
    String enrichment_results_dir = DifferentialExpressionEnrichment.enrichment_results_dir
    String visualization_dir = DifferentialExpressionEnrichment.visualization_dir
    String pca_plot = DifferentialExpressionEnrichment.pca_plot_path
    String sample_distance_heatmap = DifferentialExpressionEnrichment.sample_distance_heatmap_path
    String multiqc_report = MultiQc.report_path
  }
}

task BuildStarIndex {
  input {
    String genome_fasta_path
    String annotation_gtf_path
    String output_dir
    Int threads
    Int sjdb_overhang
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/star_index'"
    ~{docker_run} ~{image} bash -c "STAR --runThreadN ~{threads} --runMode genomeGenerate --genomeDir '~{output_dir}/star_index' --genomeFastaFiles '~{genome_fasta_path}' --sjdbGTFfile '~{annotation_gtf_path}' --sjdbOverhang ~{sjdb_overhang}"
  >>>

  output {
    String genome_dir = output_dir + "/star_index"
  }
}

task FastQc {
  input {
    String sample_id
    String read1_path
    String read2_path
    String output_dir
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/fastqc/~{sample_id}'"
    ~{docker_run} ~{image} bash -c "fastqc --threads ~{threads} --outdir '~{output_dir}/fastqc/~{sample_id}' '~{read1_path}' '~{read2_path}'"
  >>>

  output {
    String fastqc_dir_path = output_dir + "/fastqc/" + sample_id
  }
}

task TrimReads {
  input {
    String sample_id
    String read1_path
    String read2_path
    String output_dir
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/trimmed/~{sample_id}'"
    ~{docker_run} ~{image} bash -c "trim_galore --paired --cores ~{threads} --gzip --basename '~{sample_id}' --output_dir '~{output_dir}/trimmed/~{sample_id}' '~{read1_path}' '~{read2_path}'"
    ~{docker_run} ~{image} bash -c "test -s '~{output_dir}/trimmed/~{sample_id}/~{sample_id}_val_1.fq.gz'"
    ~{docker_run} ~{image} bash -c "test -s '~{output_dir}/trimmed/~{sample_id}/~{sample_id}_val_2.fq.gz'"
  >>>

  output {
    String trimmed_read1_path = output_dir + "/trimmed/" + sample_id + "/" + sample_id + "_val_1.fq.gz"
    String trimmed_read2_path = output_dir + "/trimmed/" + sample_id + "/" + sample_id + "_val_2.fq.gz"
    String trim_dir_path = output_dir + "/trimmed/" + sample_id
  }
}

task StarAlign {
  input {
    String sample_id
    String read1_path
    String read2_path
    String genome_dir
    String output_dir
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/star/~{sample_id}'"
    ~{docker_run} ~{image} bash -c "STAR --runThreadN ~{threads} --genomeDir '~{genome_dir}' --readFilesIn '~{read1_path}' '~{read2_path}' --readFilesCommand zcat --outFileNamePrefix '~{output_dir}/star/~{sample_id}/' --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts"
    ~{docker_run} ~{image} bash -c "mv '~{output_dir}/star/~{sample_id}/Aligned.sortedByCoord.out.bam' '~{output_dir}/star/~{sample_id}/~{sample_id}.sorted.bam'"
  >>>

  output {
    String sorted_bam_path = output_dir + "/star/" + sample_id + "/" + sample_id + ".sorted.bam"
    String final_log_path = output_dir + "/star/" + sample_id + "/Log.final.out"
    String gene_counts_path = output_dir + "/star/" + sample_id + "/ReadsPerGene.out.tab"
  }
}

task FeatureCounts {
  input {
    Array[String] bam_paths
    String annotation_gtf_path
    String output_dir
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/counts'"
    ~{docker_run} ~{image} bash -c "featureCounts -T ~{threads} -p -B -C -a '~{annotation_gtf_path}' -o '~{output_dir}/counts/gene_counts.tsv' ~{sep=" " bam_paths}"
  >>>

  output {
    String count_matrix_path = output_dir + "/counts/gene_counts.tsv"
    String summary_path = output_dir + "/counts/gene_counts.tsv.summary"
  }
}

task DifferentialExpressionEnrichment {
  input {
    String count_matrix_path
    Array[String] sample_ids
    Array[String] sample_groups
    String output_dir
    String enrichment_annotation_path
    Int min_count
    Float padj_cutoff
    Float log2fc_cutoff
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/differential_expression' '~{output_dir}/enrichment' '~{output_dir}/plots'"
    ~{docker_run} ~{image} Rscript - <<'RSCRIPT'
    suppressPackageStartupMessages(library(DESeq2))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(pheatmap))

    count_file <- "~{count_matrix_path}"
    output_dir <- "~{output_dir}"
    annotation_file <- "~{enrichment_annotation_path}"
    min_count <- ~{min_count}
    padj_cutoff <- ~{padj_cutoff}
    log2fc_cutoff <- ~{log2fc_cutoff}
    sample_ids <- c("~{sep='","' sample_ids}")
    sample_groups_raw <- c("~{sep='","' sample_groups}")
    sample_groups <- make.names(sample_groups_raw)

    de_dir <- file.path(output_dir, "differential_expression")
    enrich_dir <- file.path(output_dir, "enrichment")
    plot_dir <- file.path(output_dir, "plots")
    dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    metadata <- data.frame(
      sample_id = sample_ids,
      group = factor(sample_groups),
      original_group = sample_groups_raw,
      stringsAsFactors = FALSE
    )
    rownames(metadata) <- metadata$sample_id
    write.table(metadata, file.path(de_dir, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    feature_counts <- read.delim(count_file, comment.char = "#", check.names = FALSE)
    if (!"Geneid" %in% colnames(feature_counts)) {
      stop("featureCounts output must contain a Geneid column")
    }
    if (ncol(feature_counts) < 7) {
      stop("featureCounts output must contain annotation columns plus at least one count column")
    }

    gene_ids <- feature_counts$Geneid
    count_data <- feature_counts[, 7:ncol(feature_counts), drop = FALSE]
    if (ncol(count_data) != length(sample_ids)) {
      stop(sprintf("count column number (%s) does not match sample number (%s)", ncol(count_data), length(sample_ids)))
    }
    rownames(count_data) <- gene_ids
    colnames(count_data) <- sample_ids
    count_data <- round(as.matrix(count_data))
    keep <- rowSums(count_data) >= min_count
    count_data <- count_data[keep, , drop = FALSE]
    write.table(count_data, file.path(de_dir, "filtered_count_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)

    if (length(unique(metadata$group)) < 2) {
      writeLines("Only one group was provided; pairwise differential expression and enrichment were skipped.", file.path(de_dir, "pairwise_de_summary.tsv"))
      quit(save = "no", status = 0)
    }

    dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ group)
    dds <- DESeq(dds)
    vst_counts <- varianceStabilizingTransformation(dds, blind = FALSE)

    pdf(file.path(plot_dir, "pca.pdf"), width = 7, height = 6)
    print(plotPCA(vst_counts, intgroup = "group") + theme_bw())
    dev.off()

    sample_dist <- dist(t(assay(vst_counts)))
    pdf(file.path(plot_dir, "sample_distance_heatmap.pdf"), width = 8, height = 7)
    pheatmap(as.matrix(sample_dist), clustering_distance_rows = sample_dist, clustering_distance_cols = sample_dist)
    dev.off()

    annotation <- NULL
    if (nzchar(annotation_file) && file.exists(annotation_file) && file.info(annotation_file)$size > 0) {
      annotation <- read.delim(annotation_file, check.names = FALSE, stringsAsFactors = FALSE)
    }

    run_enrichment <- function(sig_genes, universe_genes, annotation, id_col, term_col, name_col, output_prefix) {
      if (is.null(annotation) || !(id_col %in% colnames(annotation)) || !(term_col %in% colnames(annotation))) {
        writeLines(sprintf("Missing annotation columns: %s and/or %s", id_col, term_col), paste0(output_prefix, ".skipped.txt"))
        return(invisible(NULL))
      }
      term_name <- if (name_col %in% colnames(annotation)) annotation[[name_col]] else annotation[[term_col]]
      anno <- data.frame(
        gene_id = annotation[[id_col]],
        term_id = annotation[[term_col]],
        term_name = term_name,
        stringsAsFactors = FALSE
      )
      anno <- anno[!is.na(anno$gene_id) & !is.na(anno$term_id) & anno$gene_id != "" & anno$term_id != "", ]
      anno <- unique(anno)
      anno <- anno[anno$gene_id %in% universe_genes, ]
      sig_genes <- intersect(sig_genes, universe_genes)
      if (length(sig_genes) == 0 || nrow(anno) == 0) {
        writeLines("No significant genes or no usable annotations for enrichment.", paste0(output_prefix, ".skipped.txt"))
        return(invisible(NULL))
      }
      terms <- sort(unique(anno$term_id))
      res <- lapply(terms, function(term) {
        term_genes <- unique(anno$gene_id[anno$term_id == term])
        a <- length(intersect(sig_genes, term_genes))
        b <- length(sig_genes) - a
        c <- length(setdiff(term_genes, sig_genes))
        d <- length(universe_genes) - a - b - c
        pvalue <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
        data.frame(
          term_id = term,
          term_name = anno$term_name[match(term, anno$term_id)],
          overlap = a,
          term_size = length(term_genes),
          query_size = length(sig_genes),
          pvalue = pvalue,
          genes = paste(intersect(sig_genes, term_genes), collapse = ","),
          stringsAsFactors = FALSE
        )
      })
      res <- do.call(rbind, res)
      res$padj <- p.adjust(res$pvalue, method = "BH")
      res <- res[order(res$padj, res$pvalue), ]
      write.table(res, paste0(output_prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
      top <- head(res[res$overlap > 0, ], 20)
      if (nrow(top) > 0) {
        top$term_label <- factor(top$term_name, levels = rev(top$term_name))
        pdf(paste0(output_prefix, ".barplot.pdf"), width = 9, height = 6)
        print(ggplot(top, aes(x = term_label, y = -log10(padj), fill = overlap)) + geom_col() + coord_flip() + theme_bw() + xlab(NULL) + ylab("-log10 adjusted P value"))
        dev.off()
      }
    }

    groups <- levels(metadata$group)
    comparisons <- combn(groups, 2, simplify = FALSE)
    summary_list <- list()
    universe_genes <- rownames(count_data)

    for (comparison in comparisons) {
      group_a <- comparison[1]
      group_b <- comparison[2]
      comparison_name <- paste0(group_b, "_vs_", group_a)
      res <- as.data.frame(results(dds, contrast = c("group", group_b, group_a)))
      res$gene_id <- rownames(res)
      res <- res[, c("gene_id", setdiff(colnames(res), "gene_id"))]
      res <- res[order(res$padj), ]
      all_path <- file.path(de_dir, paste0(comparison_name, ".all.tsv"))
      deg_path <- file.path(de_dir, paste0(comparison_name, ".deg.tsv"))
      write.table(res, all_path, sep = "\t", quote = FALSE, row.names = FALSE)

      deg <- res[!is.na(res$padj) & res$padj <= padj_cutoff & abs(res$log2FoldChange) >= log2fc_cutoff, ]
      write.table(deg, deg_path, sep = "\t", quote = FALSE, row.names = FALSE)
      summary_list[[comparison_name]] <- data.frame(
        comparison = comparison_name,
        group_a = group_a,
        group_b = group_b,
        total_tested = nrow(res),
        deg_total = nrow(deg),
        deg_up = sum(deg$log2FoldChange > 0),
        deg_down = sum(deg$log2FoldChange < 0),
        stringsAsFactors = FALSE
      )

      plot_df <- res
      plot_df$significant <- !is.na(plot_df$padj) & plot_df$padj <= padj_cutoff & abs(plot_df$log2FoldChange) >= log2fc_cutoff
      volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) + geom_point(alpha = 0.6, size = 1) + theme_bw() + ggtitle(comparison_name) + xlab("log2 fold change") + ylab("-log10 adjusted P value")
      ggsave(file.path(plot_dir, paste0(comparison_name, ".volcano.pdf")), volcano, width = 7, height = 6)

      ma <- ggplot(plot_df, aes(x = baseMean, y = log2FoldChange, color = significant)) + geom_point(alpha = 0.6, size = 1) + scale_x_log10() + theme_bw() + ggtitle(comparison_name) + xlab("mean normalized count") + ylab("log2 fold change")
      ggsave(file.path(plot_dir, paste0(comparison_name, ".ma.pdf")), ma, width = 7, height = 6)

      top_genes <- head(deg$gene_id, 50)
      if (length(top_genes) > 1) {
        pdf(file.path(plot_dir, paste0(comparison_name, ".top_deg_heatmap.pdf")), width = 9, height = 9)
        pheatmap(assay(vst_counts)[top_genes, , drop = FALSE], scale = "row", annotation_col = metadata[, "group", drop = FALSE])
        dev.off()
      }

      run_enrichment(deg$gene_id, universe_genes, annotation, "gene_id", "go_id", "go_term", file.path(enrich_dir, paste0(comparison_name, ".GO")))
      run_enrichment(deg$gene_id, universe_genes, annotation, "gene_id", "pathway_id", "pathway_name", file.path(enrich_dir, paste0(comparison_name, ".pathway")))
    }

    summary_table <- do.call(rbind, summary_list)
    write.table(summary_table, file.path(de_dir, "pairwise_de_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    RSCRIPT
  >>>

  output {
    String sample_metadata_path = output_dir + "/differential_expression/sample_metadata.tsv"
    String de_summary_path = output_dir + "/differential_expression/pairwise_de_summary.tsv"
    String de_results_dir = output_dir + "/differential_expression"
    String enrichment_results_dir = output_dir + "/enrichment"
    String visualization_dir = output_dir + "/plots"
    String pca_plot_path = output_dir + "/plots/pca.pdf"
    String sample_distance_heatmap_path = output_dir + "/plots/sample_distance_heatmap.pdf"
  }
}


task MultiQc {
  input {
    String output_dir
    String count_matrix_path
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c "test -s '~{count_matrix_path}'"
    ~{docker_run} ~{image} bash -c "mkdir -p '~{output_dir}/multiqc'"
    ~{docker_run} ~{image} bash -c "multiqc --outdir '~{output_dir}/multiqc' --filename multiqc_report.html '~{output_dir}'"
  >>>

  output {
    String report_path = output_dir + "/multiqc/multiqc_report.html"
  }
}
