#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
CALL_BROAD_PEAKS="${CALL_BROAD_PEAKS:-false}"
MACS3_QVALUE="${MACS3_QVALUE:-0.01}"
SPECIES="${SPECIES:-mouse}"
case "${SPECIES}" in
  human|hsa|hs) MACS3_GENOME_SIZE="hs" ;;
  *) MACS3_GENOME_SIZE="mm" ;;
esac
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
MACS3_IMAGE="${MACS3_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/macs3}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)

for sample in "${SAMPLES[@]}"; do
  broad_flag=""
  if [[ "${CALL_BROAD_PEAKS}" == "true" ]]; then broad_flag="--broad --broad-cutoff ${MACS3_QVALUE}"; fi
  ${DOCKER_RUN} "${MACS3_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/peaks/${sample}'; macs3 callpeak -t '${OUTPUT_DIR}/filtered/${sample}/${sample}.filtered.bam' -f BAMPE -g '${MACS3_GENOME_SIZE}' -n '${sample}' --outdir '${OUTPUT_DIR}/peaks/${sample}' -q ${MACS3_QVALUE} --keep-dup all ${broad_flag}; if [[ -f '${OUTPUT_DIR}/peaks/${sample}/${sample}_peaks.broadPeak' ]]; then peak_file='${OUTPUT_DIR}/peaks/${sample}/${sample}_peaks.broadPeak'; else peak_file='${OUTPUT_DIR}/peaks/${sample}/${sample}_peaks.narrowPeak'; fi; awk 'BEGIN{OFS=\"\\t\"} {count += 1; bp += (\$3 - \$2)} END{print \"sample\",\"peak_count\",\"peak_bp\"; print \"${sample}\",count+0,bp+0}' \"\${peak_file}\" > '${OUTPUT_DIR}/peaks/${sample}/${sample}.peak_summary.tsv'; cp \"\${peak_file}\" '${OUTPUT_DIR}/peaks/${sample}/${sample}.peaks.bed'"
done
