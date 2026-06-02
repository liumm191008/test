#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
SPECIES="${SPECIES:-mouse}"
ENRICHMENT_GENE_ID_TYPE="${ENRICHMENT_GENE_ID_TYPE:-ENTREZID}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
BIOCONDUCTOR_IMAGE="${BIOCONDUCTOR_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/bioconductor}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)

for sample in "${SAMPLES[@]}"; do
  ${DOCKER_RUN} "${BIOCONDUCTOR_IMAGE}" bash -c "set -euo pipefail; test -s '/kegg_data/peak_annotation.R'; test -s '/kegg_data/gene_cluster_enrich.R'; mkdir -p '${OUTPUT_DIR}/annotation/${sample}' '${OUTPUT_DIR}/annotation/${sample}/enrichment'; Rscript '/kegg_data/peak_annotation.R' --peak '${OUTPUT_DIR}/peaks/${sample}/${sample}.peaks.bed' --sample '${sample}' --species '${SPECIES}' --annotation-dir '${OUTPUT_DIR}/annotation/${sample}'; Rscript '/kegg_data/gene_cluster_enrich.R' '${OUTPUT_DIR}/annotation/${sample}/${sample}.enrichment_genes.tsv' '${ENRICHMENT_GENE_ID_TYPE}' '${SPECIES}' '${OUTPUT_DIR}/annotation/${sample}/enrichment' '${sample}'"
done
