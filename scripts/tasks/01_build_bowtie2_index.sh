#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
GENOME_FASTA="${GENOME_FASTA:-${WORK_ROOT}/pipeline/database/mm39/Mus_musculus.GRCm39.dna.toplevel.fa}"
BOWTIE2_INDEX="${BOWTIE2_INDEX:-${WORK_ROOT}/pipeline/database/mm39/mm39}"
THREADS="${THREADS:-8}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
BOWTIE2_IMAGE="${BOWTIE2_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/bowtie2}"

if [[ -n "${BOWTIE2_INDEX}" ]]; then
  echo "Using existing Bowtie2 index prefix: ${BOWTIE2_INDEX}"
  exit 0
fi

${DOCKER_RUN} "${BOWTIE2_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/bowtie2_index'; bowtie2-build --threads ${THREADS} '${GENOME_FASTA}' '${OUTPUT_DIR}/bowtie2_index/genome'"
