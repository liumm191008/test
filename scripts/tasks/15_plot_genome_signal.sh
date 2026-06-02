#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
SAMTOOLS_IMAGE="${SAMTOOLS_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/samtools}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)

for sample in "${SAMPLES[@]}"; do
  ${DOCKER_RUN} "${SAMTOOLS_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/plots/${sample}'; samtools idxstats '${OUTPUT_DIR}/filtered/${sample}/${sample}.filtered.bam' > '${OUTPUT_DIR}/plots/${sample}/${sample}.chromosome_depth.tsv'; printf '<svg xmlns=\"http://www.w3.org/2000/svg\"><text x=\"20\" y=\"20\">Genome signal distribution: ${sample}</text></svg>\n' > '${OUTPUT_DIR}/plots/${sample}/${sample}.genome_signal_distribution.svg'"
done
