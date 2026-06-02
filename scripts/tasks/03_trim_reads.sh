#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
RAW_DIR="${RAW_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/seq_data/rawdata}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
THREADS="${THREADS:-8}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
TRIM_GALORE_IMAGE="${TRIM_GALORE_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/trim-galore}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)

for sample in "${SAMPLES[@]}"; do
  ${DOCKER_RUN} "${TRIM_GALORE_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/trimmed/${sample}'; trim_galore --paired --cores ${THREADS} --gzip --basename '${sample}' --output_dir '${OUTPUT_DIR}/trimmed/${sample}' '${RAW_DIR}/${sample}.R1.raw.fastq.gz' '${RAW_DIR}/${sample}.R2.raw.fastq.gz'; test -s '${OUTPUT_DIR}/trimmed/${sample}/${sample}_val_1.fq.gz'; test -s '${OUTPUT_DIR}/trimmed/${sample}/${sample}_val_2.fq.gz'"
done
