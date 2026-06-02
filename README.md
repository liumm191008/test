# RNA-seq WDL workflow

This repository contains a WDL 1.0 paired-end RNA-seq workflow in `workflows/rna_seq.wdl`.

The workflow performs:

1. STAR genome-index generation, or reuse of a supplied STAR index directory.
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

The Docker runner command is configured once with `docker_run`:

```wdl
String docker_run = "docker run --rm --security-opt seccomp=unconfined  -v /home/data/vip01/work:/home/data/vip01/work"
```

Each tool image parameter contains only the image name, for example:

```wdl
String fastqc_image = "registry.cn-guangzhou.aliyuncs.com/origen/fastqc"
```

Tasks combine these two values directly:

```bash
~{docker_run} ~{image} bash -c "fastqc --threads ~{threads} --outdir '~{output_dir}/fastqc/~{sample_id}' '~{read1_path}' '~{read2_path}'"
```

Because the default Docker runner mounts `/home/data/vip01/work`, all default input and output paths should be absolute paths under `/home/data/vip01/work` unless you override `docker_run` with a different mount. The workflow defaults to the mm39 FASTA `/home/data/vip01/work/pipeline/database/mm39/Mus_musculus.GRCm39.dna.toplevel.fa` and GTF `/home/data/vip01/work/pipeline/database/mm39/Mus_musculus.GRCm39.115.gtf`.

## STAR index reuse

Set `RnaSeq.star_index_path` to an existing STAR index directory to skip `BuildStarIndex` and use that index directly. For the provided mm39 setup, use:

```json
"RnaSeq.star_index_path": "/home/data/vip01/work/pipeline/database/mm39/star_index_seqlen150"
```

If `RnaSeq.star_index_path` is empty, the workflow builds a new STAR index under `output_dir/star_index` using `RnaSeq.genome_fasta_path` and `RnaSeq.annotation_gtf_path`.

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

See `docs/rna_seq.inputs.json` for a minimal paired-end example. Update the sample FASTQ paths, sample groups, reference FASTA, annotation GTF, STAR index path, enrichment annotation table, thread count, output directory, thresholds, `docker_run`, and image names before running. The root-level `input.json` is a larger project-specific example with multiple groups.

## Example Cromwell run

```bash
java -jar cromwell.jar run workflows/rna_seq.wdl --inputs docs/rna_seq.inputs.json
```

The default output directory is `/home/data/vip01/work/rnaseq_results`.
