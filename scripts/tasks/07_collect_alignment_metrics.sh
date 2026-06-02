#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
GENOME_FASTA="${GENOME_FASTA:-${WORK_ROOT}/pipeline/database/mm39/Mus_musculus.GRCm39.dna.toplevel.fa}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
PICARD_IMAGE="${PICARD_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/picard}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)

for sample in "${SAMPLES[@]}"; do
  ${DOCKER_RUN} "${PICARD_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/picard/${sample}'; picard CollectAlignmentSummaryMetrics I='${OUTPUT_DIR}/filtered/${sample}/${sample}.filtered.bam' R='${GENOME_FASTA}' O='${OUTPUT_DIR}/picard/${sample}/${sample}.alignment_metrics.txt'"
done
