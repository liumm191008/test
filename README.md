# RNA-seq WDL workflow

This repository contains a WDL 1.0 paired-end RNA-seq workflow in `workflows/rna_seq.wdl`.

The workflow performs:

1. STAR genome-index generation.
2. Raw-read FastQC.
3. Adapter/quality trimming with Trim Galore.
4. Trimmed-read FastQC.
5. STAR alignment to coordinate-sorted BAM and per-gene STAR counts.
6. Gene-level summarization with featureCounts.
7. Pairwise differential-expression analysis across all sample groups with DESeq2.
8. GO and pathway over-representation enrichment for each pairwise DEG set.
9. PCA, sample-distance heatmap, volcano, MA, DEG heatmap, and enrichment bar-plot visualization.
10. Aggregated QC reporting with MultiQC.

## Docker execution model

Each tool image parameter includes the Docker runner and the image name, for example:

```wdl
String fastqc_image = "docker run --rm -v /data:/data registry.cn-guangzhou.aliyuncs.com/origen/fastqc"
```

Tasks use that combined image command directly, so no separate Docker runner parameter needs to be passed to every task:

```bash
~{image} bash -c "fastqc --threads ~{threads} --outdir '~{output_dir}/fastqc/~{sample_id}' '~{read1_path}' '~{read2_path}'"
```

Because only `/data` is mounted into each container by the default image commands, all input and output paths in the WDL inputs must be absolute paths under `/data` unless you override the image command strings with additional mounts.

## Sample groups and downstream analysis

Each sample object must include `sample_id`, `group`, `read1_path`, and `read2_path`. If two or more groups are present, the workflow runs all pairwise group comparisons and writes:

- `differential_expression/*.all.tsv`: full DESeq2 result table for each comparison.
- `differential_expression/*.deg.tsv`: filtered DEG table for each comparison using `padj_cutoff` and `log2fc_cutoff`.
- `differential_expression/pairwise_de_summary.tsv`: DEG counts for all pairwise comparisons.
- `enrichment/*.GO.tsv` and `enrichment/*.pathway.tsv`: GO/pathway enrichment results when an enrichment annotation table is supplied.
- `plots/*.pdf`: PCA, sample distance heatmap, volcano plots, MA plots, top-DEG heatmaps, and enrichment bar plots.

The optional `enrichment_annotation_path` should point to a tab-delimited file with these columns when GO/pathway enrichment is desired:

```text
gene_id	go_id	go_term	pathway_id	pathway_name
```

The `gene_id` values must match the `Geneid` column emitted by featureCounts. If `enrichment_annotation_path` is empty or missing required columns, DESeq2 and visualization still run, and enrichment is skipped with a `.skipped.txt` note. The `differential_analysis_image` must provide R with DESeq2, ggplot2, and pheatmap installed.

## Example inputs

See `docs/rna_seq.inputs.json` for a minimal paired-end example. Update the sample FASTQ paths, sample groups, reference FASTA, annotation GTF, enrichment annotation table, thread count, output directory, thresholds, and image commands before running. The root-level `input.json` is a larger project-specific example with multiple groups.

## Example Cromwell run

```bash
java -jar cromwell.jar run workflows/rna_seq.wdl --inputs docs/rna_seq.inputs.json
```

The default output directory is `/data/rnaseq_results`.
