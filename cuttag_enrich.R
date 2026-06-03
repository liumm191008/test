#!/usr/bin/env Rscript

############################################################
# CUT&Tag enrichment pipeline
# Author: ChatGPT
# Version: stable offline KEGG
############################################################

suppressPackageStartupMessages({
    library(ChIPseeker)
    library(clusterProfiler)
    library(enrichplot)
    library(GenomeInfoDb)
    library(gson)
    library(org.Mm.eg.db)
    library(TxDb.Mmusculus.UCSC.mm39.knownGene)
    library(ggplot2)
})

usage <- function() {
    cat("Usage:\n\n",
        "Rscript cuttag_enrich.R \\\n",
        "    peak.narrowPeak \\\n",
        "    kegg_mmu.gson \\\n",
        "    kegg_class.tsv \\\n",
        "    output_dir \\\n",
        "    prefix\n",
        sep = "")
}

log_step <- function(...) {
    message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

validate_file <- function(path, label) {
    if (!file.exists(path)) {
        stop(label, " not found: ", path, call. = FALSE)
    }
    if (file.info(path)$size == 0) {
        stop(label, " is empty: ", path, call. = FALSE)
    }
    invisible(path)
}

write_result <- function(result_df, path) {
    write.csv(result_df, path, row.names = FALSE)
    log_step("Wrote ", nrow(result_df), " rows: ", path)
}

empty_kegg_class_df <- function() {
    data.frame(
        ID = character(0),
        category = character(0),
        subcategory = character(0),
        Description = character(0),
        stringsAsFactors = FALSE
    )
}

normalize_gson_term2name <- function(kk_gson) {
    term2name <- kk_gson@gsid2name

    if (is.data.frame(term2name) && ncol(term2name) >= 2) {
        return(data.frame(
            ID = as.character(term2name[[1]]),
            Description = as.character(term2name[[2]]),
            stringsAsFactors = FALSE
        ))
    }

    if (is.list(term2name) && !is.null(names(term2name))) {
        return(data.frame(
            ID = names(term2name),
            Description = as.character(unlist(term2name, use.names = FALSE)),
            stringsAsFactors = FALSE
        ))
    }

    if (is.atomic(term2name) && !is.null(names(term2name))) {
        return(data.frame(
            ID = names(term2name),
            Description = as.character(term2name),
            stringsAsFactors = FALSE
        ))
    }

    empty_kegg_class_df()[, c("ID", "Description")]
}

normalize_gson_term_ids <- function(kk_gson) {
    term2gene <- kk_gson@gsid2gene

    if (is.data.frame(term2gene) && ncol(term2gene) >= 1) {
        return(unique(as.character(term2gene[[1]])))
    }

    if (is.list(term2gene) && !is.null(names(term2gene))) {
        return(names(term2gene))
    }

    character(0)
}

make_kegg_class_from_gson <- function(kk_gson) {
    term_ids <- normalize_gson_term_ids(kk_gson)
    term2name <- normalize_gson_term2name(kk_gson)

    if (length(term_ids) == 0 && nrow(term2name) > 0) {
        term_ids <- term2name$ID
    }

    if (length(term_ids) == 0) {
        stop("Unable to infer KEGG pathway IDs from gson; please provide a non-empty kegg_class.tsv", call. = FALSE)
    }

    class_df <- data.frame(
        ID = unique(term_ids),
        category = "Unclassified",
        subcategory = "Unclassified",
        stringsAsFactors = FALSE
    )
    class_df <- merge(class_df, term2name, by = "ID", all.x = TRUE, sort = FALSE)
    missing_description <- is.na(class_df$Description) | !nzchar(class_df$Description)
    class_df$Description[missing_description] <- class_df$ID[missing_description]
    class_df[, c("ID", "category", "subcategory", "Description")]
}

prepare_kegg_class_df <- function(class_file, kk_gson) {
    validate_file(class_file, "KEGG class file")
    class_df <- read.delim(
        class_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if (!"ID" %in% colnames(class_df)) {
        stop("KEGG class file must contain an ID column", call. = FALSE)
    }

    if (nrow(class_df) == 0) {
        warning(
            "KEGG class file contains 0 rows; using gson-derived descriptions with Unclassified categories.",
            call. = FALSE
        )
        class_df <- make_kegg_class_from_gson(kk_gson)
    }

    class_df
}

html_escape <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub('"', "&quot;", x, fixed = TRUE)
    x
}

extract_kegg_pathway_id <- function(pathway_id, organism = "mmu") {
    pathway_id <- sub("^path:", "", pathway_id)
    pathway_id <- sub("^map", organism, pathway_id)

    if (grepl("^[0-9]{5}$", pathway_id)) {
        pathway_id <- paste0(organism, pathway_id)
    }

    pathway_id
}

split_gene_ids <- function(gene_ids) {
    gene_ids <- as.character(gene_ids)
    if (is.na(gene_ids) || !nzchar(gene_ids)) {
        return(character(0))
    }
    unique(strsplit(gene_ids, "/", fixed = TRUE)[[1]])
}

make_kegg_pathway_url <- function(pathway_id, gene_ids, organism = "mmu",
                                  bg_color = "#FF9999", fg_color = "#000000") {
    pathway_id <- extract_kegg_pathway_id(pathway_id, organism)
    genes <- split_gene_ids(gene_ids)
    genes <- genes[!is.na(genes) & nzchar(genes)]

    if (length(genes) == 0) {
        return(paste0("https://www.kegg.jp/kegg-bin/show_pathway?", pathway_id))
    }

    genes <- ifelse(grepl(":", genes, fixed = TRUE), genes, paste0(organism, ":", genes))
    color_lines <- paste0(genes, " ", bg_color, ",", fg_color)
    color_query <- utils::URLencode(paste(color_lines, collapse = "\r\n"), reserved = TRUE)
    color_query <- gsub("%20", "+", color_query, fixed = TRUE)

    paste0(
        "https://www.kegg.jp/kegg-bin/show_pathway?map=",
        pathway_id,
        "&multi_query=",
        color_query
    )
}

add_kegg_pathway_links <- function(kegg_df, organism = "mmu") {
    if (nrow(kegg_df) == 0) {
        kegg_df$Pathway_Link <- character(0)
        return(kegg_df)
    }

    if (!all(c("ID", "geneID") %in% colnames(kegg_df))) {
        stop("KEGG result must contain ID and geneID columns to create pathway links", call. = FALSE)
    }

    kegg_df$Pathway_Link <- mapply(
        make_kegg_pathway_url,
        pathway_id = kegg_df$ID,
        gene_ids = kegg_df$geneID,
        MoreArgs = list(organism = organism),
        USE.NAMES = FALSE
    )

    kegg_df
}

write_kegg_html_result <- function(kegg_df, path) {
    if (nrow(kegg_df) == 0) {
        writeLines(
            c(
                "<!doctype html>",
                "<html><head><meta charset=\"UTF-8\"><title>KEGG enrichment</title></head>",
                "<body><p>No KEGG enrichment terms found.</p></body></html>"
            ),
            path
        )
        log_step("Wrote KEGG HTML result: ", path)
        return(invisible(path))
    }

    display_df <- kegg_df
    link_col <- display_df$Pathway_Link
    display_df$Pathway_Link <- paste0(
        "<a href=\"", html_escape(link_col),
        "\" target=\"_blank\" rel=\"noopener noreferrer\">View highlighted pathway</a>"
    )

    header <- paste0("<th>", html_escape(colnames(display_df)), "</th>", collapse = "")
    rows <- apply(display_df, 1, function(row) {
        cells <- mapply(function(value, name) {
            if (name == "Pathway_Link") {
                paste0("<td>", value, "</td>")
            } else {
                paste0("<td>", html_escape(value), "</td>")
            }
        }, row, names(row), USE.NAMES = FALSE)
        paste0("<tr>", paste0(cells, collapse = ""), "</tr>")
    })

    html <- c(
        "<!doctype html>",
        "<html>",
        "<head>",
        "<meta charset=\"UTF-8\">",
        "<title>KEGG enrichment</title>",
        "<style>",
        "body{font-family:Arial,sans-serif;margin:24px;}",
        "table{border-collapse:collapse;font-size:13px;}",
        "th,td{border:1px solid #ddd;padding:6px 8px;vertical-align:top;}",
        "th{background:#f2f2f2;position:sticky;top:0;}",
        "tr:nth-child(even){background:#fafafa;}",
        "</style>",
        "</head>",
        "<body>",
        "<h1>KEGG enrichment result</h1>",
        "<p>Click <strong>Pathway_Link</strong> to open the KEGG pathway map. Genes from this sample are highlighted on the KEGG interactive map; hover nodes for details and click nodes to open KEGG database entries.</p>",
        "<table>",
        paste0("<thead><tr>", header, "</tr></thead>"),
        "<tbody>",
        rows,
        "</tbody>",
        "</table>",
        "</body></html>"
    )

    writeLines(html, path)
    log_step("Wrote KEGG HTML result: ", path)
    invisible(path)
}

plot_enrichment <- function(enrich_result, result_df, plot_dir, prefix, label,
                            show_category = 15, width = 10, height = 8,
                            draw_barplot = TRUE) {
    if (nrow(result_df) == 0) {
        log_step("No ", label, " enrichment terms found; skipping plots")
        return(invisible(FALSE))
    }

    dotplot_file <- file.path(plot_dir, paste0(prefix, "_", label, "_dotplot.pdf"))
    barplot_file <- file.path(plot_dir, paste0(prefix, "_", label, "_barplot.pdf"))

    log_step("Plotting ", label, " dotplot: ", dotplot_file)
    pdf(dotplot_file, width = width, height = height)
    print(enrichplot::dotplot(enrich_result, showCategory = show_category))
    dev.off()

    if (draw_barplot) {
        log_step("Plotting ", label, " barplot: ", barplot_file)
        pdf(barplot_file, width = width, height = height)
        print(barplot(enrich_result, showCategory = show_category))
        dev.off()
    }

    invisible(TRUE)
}

go_ontology_labels <- c(
    BP = "biological_process",
    CC = "cellular_component",
    MF = "molecular_function"
)

go_ontology_colors <- c(
    biological_process = "#F8766D",
    cellular_component = "#00BA38",
    molecular_function = "#619CFF"
)

run_go_classification <- function(genes, level = 2) {
    go_class_list <- lapply(names(go_ontology_labels), function(ontology) {
        go_class <- groupGO(
            gene = genes,
            OrgDb = org.Mm.eg.db,
            keyType = "ENTREZID",
            ont = ontology,
            level = level,
            readable = TRUE
        )
        go_class_df <- as.data.frame(go_class)
        go_class_df$ONTOLOGY <- go_ontology_labels[[ontology]]
        go_class_df
    })

    go_class_df <- do.call(rbind, go_class_list)
    go_class_df$Count <- as.numeric(go_class_df$Count)
    go_class_df <- go_class_df[!is.na(go_class_df$Count) & go_class_df$Count > 0, ]
    row.names(go_class_df) <- NULL

    go_class_df
}

plot_go_classification_bar <- function(go_class_df, plot_dir, prefix,
                                       width = 10, height = 9) {
    if (nrow(go_class_df) == 0) {
        log_step("No GO classification terms found; skipping GO classification barplot")
        return(invisible(FALSE))
    }

    go_class_df$ONTOLOGY <- factor(
        go_class_df$ONTOLOGY,
        levels = unname(go_ontology_labels)
    )
    go_class_df$Description <- factor(
        go_class_df$Description,
        levels = rev(unique(go_class_df$Description))
    )

    barplot_file <- file.path(plot_dir, paste0(prefix, "_GO_barplot.pdf"))
    log_step("Plotting GO classification barplot: ", barplot_file)

    pdf(barplot_file, width = width, height = height)
    print(
        ggplot(go_class_df, aes(x = Count, y = Description, fill = ONTOLOGY)) +
            geom_col(width = 0.55) +
            facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
            scale_fill_manual(values = go_ontology_colors, guide = "none") +
            labs(x = "Number", y = "Description", caption = "GO 分类注释结果") +
            theme_gray(base_size = 11) +
            theme(
                plot.caption = element_text(hjust = 0.5, size = 14, margin = margin(t = 18)),
                axis.text.y = element_text(size = 8),
                strip.text.y = element_text(size = 9),
                panel.spacing.y = grid::unit(0.15, "lines")
            )
    )
    dev.off()

    invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
    usage()
    quit(save = "no", status = 1)
}

peak_file <- args[1]
gson_file <- args[2]
class_file <- args[3]
outdir <- args[4]
prefix <- args[5]

validate_file(peak_file, "Peak file")
validate_file(gson_file, "KEGG gson file")

if (!nzchar(prefix)) {
    stop("prefix must not be empty", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(outdir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

log_step("Reading peaks: ", peak_file)
peaks <- readPeakFile(peak_file)

log_step("Normalizing chromosome names to UCSC style")
seqlevelsStyle(peaks) <- "UCSC"
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

if (length(peaks) == 0) {
    stop("No standard-chromosome peaks remained after filtering", call. = FALSE)
}

log_step("Annotating ", length(peaks), " peaks")
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
peak_anno <- annotatePeak(
    peaks,
    TxDb = txdb,
    tssRegion = c(-3000, 3000),
    annoDb = "org.Mm.eg.db"
)
anno_df <- as.data.frame(peak_anno)

write_result(
    anno_df,
    file.path(outdir, paste0(prefix, "_peak_annotation.csv"))
)

genes <- unique(anno_df$geneId)
genes <- genes[!is.na(genes) & nzchar(genes)]
log_step("Unique annotated genes: ", length(genes))

if (length(genes) == 0) {
    stop("No genes found in peak annotation", call. = FALSE)
}

log_step("Running GO BP enrichment")
ego <- enrichGO(
    gene = genes,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)
go_df <- as.data.frame(ego)
write_result(go_df, file.path(outdir, paste0(prefix, "_GO.csv")))

log_step("Running GO classification")
go_class_df <- run_go_classification(genes, level = 2)
write_result(
    go_class_df,
    file.path(outdir, paste0(prefix, "_GO_classification.csv"))
)

log_step("Running offline KEGG enrichment")
kk_gson <- read.gson(gson_file)
kegg <- enricher(
    gene = genes,
    TERM2GENE = kk_gson@gsid2gene,
    TERM2NAME = kk_gson@gsid2name,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
)
kegg_df <- as.data.frame(kegg)

log_step("Adding KEGG categories")
class_df <- prepare_kegg_class_df(
    class_file = class_file,
    kk_gson = kk_gson
)

kegg_df$.enrich_order <- seq_len(nrow(kegg_df))
kegg_df <- merge(class_df, kegg_df, by = "ID", all.y = TRUE, sort = FALSE)
kegg_df <- kegg_df[order(kegg_df$.enrich_order), , drop = FALSE]
kegg_df$.enrich_order <- NULL
kegg_df <- add_kegg_pathway_links(kegg_df, organism = "mmu")

write_result(kegg_df, file.path(outdir, paste0(prefix, "_KEGG.csv")))
write_kegg_html_result(
    kegg_df,
    file.path(outdir, paste0(prefix, "_KEGG.html"))
)

plot_enrichment(ego, go_df, plot_dir, prefix, "GO", draw_barplot = FALSE)
plot_go_classification_bar(go_class_df, plot_dir, prefix)
plot_enrichment(kegg, kegg_df, plot_dir, prefix, "KEGG")

log_step("================================")
log_step("Analysis finished")
log_step("================================")
