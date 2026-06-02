#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
RAW_DIR="${RAW_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/seq_data/rawdata}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
THREADS="${THREADS:-8}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
FASTQC_IMAGE="${FASTQC_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/fastqc}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)
MODE="${1:-raw}"

for sample in "${SAMPLES[@]}"; do
  if [[ "${MODE}" == "trimmed" ]]; then
    r1="${OUTPUT_DIR}/trimmed/${sample}/${sample}_val_1.fq.gz"
    r2="${OUTPUT_DIR}/trimmed/${sample}/${sample}_val_2.fq.gz"
    sample_label="${sample}.trimmed"
  else
    r1="${RAW_DIR}/${sample}.R1.raw.fastq.gz"
    r2="${RAW_DIR}/${sample}.R2.raw.fastq.gz"
    sample_label="${sample}"
  fi
  ${DOCKER_RUN} "${FASTQC_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/fastqc/${sample_label}'; fastqc --threads ${THREADS} --outdir '${OUTPUT_DIR}/fastqc/${sample_label}' '${r1}' '${r2}'; for zip_file in '${OUTPUT_DIR}/fastqc/${sample_label}'/*.zip; do [[ -f \"\${zip_file}\" ]] && unzip -o \"\${zip_file}\" -d '${OUTPUT_DIR}/fastqc/${sample_label}'; done"
done
