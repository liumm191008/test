#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
THREADS="${THREADS:-8}"
BIGWIG_NORMALIZATION="${BIGWIG_NORMALIZATION:-CPM}"
EFFECTIVE_GENOME_SIZE="${EFFECTIVE_GENOME_SIZE:-2652783500}"
BIGWIG_BIN_SIZE="${BIGWIG_BIN_SIZE:-10}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
SAMTOOLS_IMAGE="${SAMTOOLS_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/samtools}"
DEEPTOOLS_IMAGE="${DEEPTOOLS_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/deeptools}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)

for sample in "${SAMPLES[@]}"; do
  bam="${OUTPUT_DIR}/filtered/${sample}/${sample}.filtered.bam"
  ${DOCKER_RUN} "${SAMTOOLS_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/tracks/${sample}'; samtools sort -@ ${THREADS} -n -o '${OUTPUT_DIR}/tracks/${sample}/${sample}.name_sorted.bam' '${bam}'; samtools view -@ ${THREADS} -f 2 -F 1804 '${OUTPUT_DIR}/tracks/${sample}/${sample}.name_sorted.bam' | awk 'BEGIN{OFS=\"\\t\"} \$9 > 0 {start=\$4-1; end=start+\$9; if (end > start) print \$3,start,end,\".\",0,\".\"}' | sort -k1,1 -k2,2n > '${OUTPUT_DIR}/tracks/${sample}/${sample}.fragments.bed'"
  if [[ "${BIGWIG_NORMALIZATION}" == "RPGC" ]]; then
    ${DOCKER_RUN} "${DEEPTOOLS_IMAGE}" bash -c "set -euo pipefail; bamCoverage --bam '${bam}' --outFileName '${OUTPUT_DIR}/tracks/${sample}/${sample}.bw' --outFileFormat bigwig --normalizeUsing RPGC --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} --binSize ${BIGWIG_BIN_SIZE} --extendReads --numberOfProcessors ${THREADS}"
  else
    ${DOCKER_RUN} "${DEEPTOOLS_IMAGE}" bash -c "set -euo pipefail; bamCoverage --bam '${bam}' --outFileName '${OUTPUT_DIR}/tracks/${sample}/${sample}.bw' --outFileFormat bigwig --normalizeUsing '${BIGWIG_NORMALIZATION}' --binSize ${BIGWIG_BIN_SIZE} --extendReads --numberOfProcessors ${THREADS}"
  fi
done
