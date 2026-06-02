# CUT&Tag WDL workflow

This repository contains a WDL 1.0 workflow for paired-end CUT&Tag data. The workflow follows the same host-path and `docker_run` style as the reference pipeline: all input and output paths are strings, all files are expected to be absolute paths visible inside the Docker mount, and each task invokes a configured container image through `docker run`.

## Workflow

`workflows/cuttag.wdl` implements these steps:

1. **Raw read quality control** with `FastQC`.
2. **Adapter and quality trimming** with `Trim Galore`.
3. **Trimmed read quality control** with `FastQC`.
4. **Reference alignment** with `Bowtie2 --end-to-end --very-sensitive -I 10 -X 700`, followed by coordinate sorting and indexing with `samtools`.
5. **Duplicate handling** with `Picard AddOrReplaceReadGroups` inside the MarkDuplicates task, followed by `Picard MarkDuplicates`; set `deduplicate=true` to remove duplicates after collecting duplication metrics.
6. **BAM filtering** with `samtools` for properly paired reads, mapping quality, SAM flags, optional mitochondrial contig removal, optional primary-contig filtering for no-`chr` references, and optional blacklist removal.
7. **Alignment and insert-size metrics** with Picard, including an insert-size histogram PDF.
8. **Track generation** with `samtools` fragment BED output and `deepTools bamCoverage` bigWig output.
9. **Peak calling** with `MACS3` in paired-end mode. Narrow peak calling is the default, and broad peak calling can be enabled with `call_broad_peaks`.
10. **Peak annotation** with ChIPseeker, including annotation pie chart, TSS distance plot, genomic feature distribution, and GO/pathway enrichment through `/kegg_data/gene_cluster_enrich.R`.
11. **Statistics visualizations** including peak functional-element distribution, genome-wide signal distribution, sample fingerprint curves, peak/gene heatmaps and profiles, and peak-length histograms.
12. **E. coli spike-in alignment** using the default Bowtie2 index when `spikein_index_path` is not empty.
13. **Sample sheet and QC aggregation** with a metadata table and `MultiQC` report.

## Default reference paths

The workflow defaults are set to the requested mm39 resources:

| Input | Default |
| --- | --- |
| `bowtie2_index_path` | `/home/data/vip01/work/pipeline/database/mm39/mm39` |
| `genome_fasta_path` | `/home/data/vip01/work/pipeline/database/mm39/Mus_musculus.GRCm39.dna.toplevel.fa` |
| `annotation_gtf_path` | `/home/data/vip01/work/pipeline/database/mm39/Mus_musculus.GRCm39.115.gtf` |
| `txdb_path` | `/home/data/vip01/work/pipeline/database/mm39/TxDb.Mmusculus.GRCm39.115.sqlite` |
| `spikein_index_path` | `/home/data/vip01/work/pipeline/database/E.coli/MG1655` |

The mm39 FASTA is expected to use contig names without a `chr` prefix and may include alternate/unplaced contigs such as `JH584299.1`. By default, `keep_primary_contigs_only=true` retains contigs matching `primary_contig_regex` so these alternate contigs are removed before downstream tracks and peaks.

If `bowtie2_index_path` is empty, the workflow builds a Bowtie2 index from `genome_fasta_path` under `output_dir/bowtie2_index`.

## Inputs

See `examples/cuttag.inputs.json` for a complete template.

Required inputs:

| Input | Description |
| --- | --- |
| `samples` | Array of paired-end samples with `sample_id`, `group`, `read1_path`, and `read2_path`. |
| `bowtie2_index_path` or `genome_fasta_path` | Existing Bowtie2 index prefix, or a FASTA used to build an index when `bowtie2_index_path` is empty. |
| `output_dir` | Absolute output directory visible through the `docker_run` mount. |
| `docker_run` | Docker command and bind mounts used by every task; defaults include `docker run -i --rm`. |

Common optional inputs:

| Input | Default | Description |
| --- | --- | --- |
| `blacklist_bed_path` | empty | BED intervals removed from the filtered BAM when provided. |
| `spikein_index_path` | `/home/data/vip01/work/pipeline/database/E.coli/MG1655` | Optional Bowtie2 index prefix for spike-in QC. Set to empty to skip spike-in alignment. |
| `threads` | `8` | CPU threads per task. |
| `mapq` | `30` | Minimum mapping quality for retained alignments. |
| `exclude_mito` | `true` | Remove mitochondrial alignments before final outputs. |
| `keep_primary_contigs_only` | `true` | Remove alternate/unplaced contigs before downstream outputs. |
| `primary_contig_regex` | `^([1-9]|1[0-9]|X|Y|MT)$` | Regex for primary no-`chr` contig names retained when `keep_primary_contigs_only` is true. |
| `deduplicate` | `false` | Keep duplicate reads by default; set true to remove duplicates after Picard metrics are generated. |
| `read_group_platform` | `ILLUMINA` | `RGPL` value passed to Picard `AddOrReplaceReadGroups` inside the duplicate-handling task. |
| `call_broad_peaks` | `false` | Add MACS3 broad peak mode. |
| `macs3_qvalue` | `0.01` | MACS3 q-value cutoff. |
| `bigwig_normalization` | `CPM` | deepTools normalization method. Use `RPGC` with a positive `effective_genome_size` for genome-size scaling. |
| `effective_genome_size` | `2652783500` | Effective genome size used by deepTools when `bigwig_normalization` is `RPGC`. |
| `bigwig_bin_size` | `10` | bigWig bin size in base pairs. |
| `species` | `mouse` | Species used to derive MACS3 genome size (`mm`/`hs`), OrgDb (`org.Mm.eg.db`/`org.Hs.eg.db`), KEGG organism (`mmu`/`hsa`), and peak annotation behavior. |
| `txdb_path` | mm39 TxDb SQLite path | TxDb object used by downstream differential-peak annotation; per-sample `peak_annotation.R` chooses its TxDb package from `species`. |
| `enrichment_gene_id_type` | `ENTREZID` | Gene ID type passed to `/kegg_data/gene_cluster_enrich.R` for sample-level and differential-peak enrichment. |
| `motif_kmer_size` | `6` | HOMER motif length; the workflow also requests 8, 10, and 12 bp motifs. |
| `diff_peak_pvalue` | `0.05` | DiffBind report threshold for differential peaks. |
| `diff_peak_min_fraction_delta` | `0.5` | Legacy input retained for compatibility; DiffBind uses `diff_peak_pvalue`. |

## Container images

The main alignment, QC, and peak-calling images are configured as requested:

| Input | Default image |
| --- | --- |
| `bowtie2_image` | `registry.cn-guangzhou.aliyuncs.com/origen/bowtie2` |
| `samtools_image` | `registry.cn-guangzhou.aliyuncs.com/origen/samtools` |
| `picard_image` | `registry.cn-guangzhou.aliyuncs.com/origen/picard` |
| `deeptools_image` | `registry.cn-guangzhou.aliyuncs.com/origen/deeptools` |
| `macs3_image` | `registry.cn-guangzhou.aliyuncs.com/origen/macs3` |
| `bioconductor_image` | `registry.cn-guangzhou.aliyuncs.com/origen/bioconductor` |
| `homer_image` | `registry.cn-guangzhou.aliyuncs.com/origen/homer` |

The workflow also keeps configurable images for read-level QC and report aggregation:

| Input | Default image |
| --- | --- |
| `fastqc_image` | `registry.cn-guangzhou.aliyuncs.com/origen/fastqc` |
| `trim_galore_image` | `registry.cn-guangzhou.aliyuncs.com/origen/trim-galore` |
| `multiqc_image` | `registry.cn-guangzhou.aliyuncs.com/origen/multiqc` |


## Example test dataset

`examples/cuttag.inputs.json` is configured for the provided MJ20260515140 test dataset with three groups and three biological replicates per group:

| Group | Samples | FASTQ directory |
| --- | --- | --- |
| `con` | `con1`, `con2`, `con3` | `/home/data/vip01/work/bioproject/MJ20260515140/seq_data/rawdata` |
| `exo` | `exo1`, `exo2`, `exo3` | `/home/data/vip01/work/bioproject/MJ20260515140/seq_data/rawdata` |
| `LPS` | `LPS1`, `LPS2`, `LPS3` | `/home/data/vip01/work/bioproject/MJ20260515140/seq_data/rawdata` |

The example writes results to `/home/data/vip01/work/bioproject/MJ20260515140/cuttag_results`.

## Peak annotation helper script

Differential peak contrasts are split into gain (`Fold`/log fold-change > 0) and loss (`Fold`/log fold-change < 0) peak sets; each subset is annotated, enriched with `/kegg_data/gene_cluster_enrich.R`, and analyzed with HOMER motif discovery.

`scripts/peak_annotation.R` accepts MACS3 `narrowPeak`, `broadPeak`, or BED-like peak files, supports `--species mouse` and `--species human`, chooses the TxDb package from `--species` (`TxDb.Mmusculus.UCSC.mm39.knownGene` for mouse; `TxDb.Hsapiens.UCSC.hg38.knownGene` for human) instead of accepting `--txdb`, and writes annotated peaks, annotation summaries/statistics, annotation plots directly under `--annotation-dir`, and the gene list (`*.enrichment_genes.tsv`). The WDL and direct shell task always call this script from the hard-coded path `/kegg_data/peak_annotation.R` inside `registry.cn-guangzhou.aliyuncs.com/origen/bioconductor`, then run `/kegg_data/gene_cluster_enrich.R` on each sample gene list and on each differential-peak contrast gene list to produce GO, KEGG/pathway, and Reactome enrichment outputs.

## Direct shell scripts for the test dataset

In addition to the WDL, each workflow task has a directly runnable shell script under `scripts/tasks/`. Each script is self-contained and only defines the sample list plus the reference paths, container images, runtime parameters, and output paths needed by that specific task, so there is no separate shared config file to source. Run one task at a time, for example:

```bash
scripts/tasks/03_trim_reads.sh
scripts/tasks/04_align_reads.sh
scripts/tasks/10_call_peaks.sh
```

To execute the whole task sequence without a WDL engine, run:

```bash
scripts/run_all_tasks.sh
```

You can override defaults by exporting variables before running a script, for example `THREADS=16 DEDUPLICATE=true scripts/tasks/05_mark_duplicates.sh`.

## Running with miniwdl

```bash
miniwdl run workflows/cuttag.wdl \
  --input examples/cuttag.inputs.json
```

## Running with Cromwell

```bash
java -jar cromwell.jar run workflows/cuttag.wdl \
  --inputs examples/cuttag.inputs.json
```

## Outputs

Important workflow outputs include:

| Output | Description |
| --- | --- |
| `reference_fasta` / `annotation_gtf` / `bowtie2_index` | Reference paths used by the run. |
| `raw_fastqc_dirs` / `trimmed_fastqc_dirs` | FastQC output directories before and after trimming. |
| `trimmed_read1` / `trimmed_read2` | Trimmed paired FASTQ files. |
| `sorted_bams` / `alignment_logs` | Coordinate-sorted alignment BAMs and Bowtie2 logs. |
| `duplicate_marked_bams` / `duplication_metrics` | Picard duplicate-marked or deduplicated BAMs and duplication metrics. |
| `filtered_bams` / `filtered_bam_indexes` | Final filtered BAM files and indexes. |
| `filter_stats` / `picard_alignment_metrics` | samtools and Picard alignment QC outputs. |
| `insert_size_metrics` / `insert_size_plots` | Picard insert-size statistics and histogram PDFs. |
| `fragments` | Fragment BED files inferred from paired-end alignments. |
| `bigwigs` | Normalized bigWig signal tracks. |
| `peak_files` / `peak_summaries` | MACS3 peak BED files and per-sample peak summaries. |
| `annotated_peaks` / `peak_annotation_summaries` / `peak_annotation_stats` | ChIPseeker peak annotation tables, annotation summaries, and per-sample annotation statistics from `peak_annotation.R`. |
| `peak_enrichment_gene_lists` / `peak_enrichment_dirs` | Gene lists exported by `peak_annotation.R` and per-sample GO/KEGG/Reactome enrichment outputs from `gene_cluster_enrich.R`. |
| `annotation_pie_charts` / `tss_distance_plots` / `genomic_feature_distribution_plots` | ChIPseeker annotation pie, TSS distance, and genomic feature plots. |
| `peak_feature_distribution_plots` / `all_peak_feature_distribution_plot` | Per-sample and all-sample peak functional-element distribution plots. |
| `peak_length_distribution_plots` | Per-sample peak length histogram plots. |
| `genome_signal_distribution_plots` | Per-sample genome-wide signal distribution plots by chromosome. |
| `fingerprint_plot` / `fingerprint_counts` | deepTools fingerprint plot and raw count table across samples from the renamed `PlotPeak` task. |
| `peak_reference_point_heatmap` / `peak_reference_point_profile` / `gene_body_heatmap` / `gene_body_profile` | deepTools heatmap/profile plots over consensus peaks and gene bodies from the `PlotPeak` task. |
| `motif_tables` / `motif_plots` | HOMER known motif table and HOMER motif HTML report generated from all called peaks. |
| `differential_peak_*` / `diffbind_report` / `differential_consensus_peak_set` | DiffBind differential peak results, consensus peak BED set, saved DiffBind object, count matrices, PCA/correlation plots (PCA labels samples) plus per-contrast MA/volcano/box plots written directly in the differential peak output directory, gain/loss-classified differential peak annotations, separate gain/loss HOMER motif results, GO/KEGG/Reactome enrichment directory, and visualization directory. |
| `spikein_bams` / `spikein_logs` | Optional spike-in alignment BAMs and Bowtie2 logs. |
| `sample_sheet` / `reference_paths` | Sample metadata table and reference-path metadata table. |
| `multiqc_report` | Aggregated HTML QC report. |
