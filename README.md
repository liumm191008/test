# RNA-seq WDL workflow

This repository contains a WDL 1.0 paired-end RNA-seq workflow in `workflows/rna_seq.wdl`.

The workflow performs:

1. STAR genome-index generation.
2. Raw-read FastQC.
3. Adapter/quality trimming with Trim Galore.
4. Trimmed-read FastQC.
5. STAR alignment to coordinate-sorted BAM and per-gene STAR counts.
6. Gene-level summarization with featureCounts.
7. Aggregated QC reporting with MultiQC.

## Docker execution model

Each tool image parameter includes the Docker runner and the image name, for example:

```wdl
String fastqc_image = "docker run --rm -v /data:/data registry.cn-guangzhou.aliyuncs.com/origen/fastqc"
```

Tasks use that combined image command directly, so no separate Docker runner parameter needs to be passed to every task:

```bash
~{image} bash -c "fastqc --threads ~{threads} --outdir '~{output_dir}/fastqc/~{sample_id}' '~{read1_path}' '~{read2_path}'"
```

Because only `/data` is mounted into each container by the default image commands, all input and output paths in the WDL inputs must be absolute paths under `/data`.

## Example inputs

See `docs/rna_seq.inputs.json` for a minimal paired-end example. Update the sample FASTQ paths, reference FASTA, annotation GTF, thread count, output directory, and image commands before running.

## Example Cromwell run

```bash
java -jar cromwell.jar run workflows/rna_seq.wdl --inputs docs/rna_seq.inputs.json
```

The default output directory is `/data/rnaseq_results`.
