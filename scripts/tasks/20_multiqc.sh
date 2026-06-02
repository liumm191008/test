#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
MULTIQC_IMAGE="${MULTIQC_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/multiqc}"

${DOCKER_RUN} "${MULTIQC_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/multiqc'; multiqc --outdir '${OUTPUT_DIR}/multiqc' --filename multiqc_report.html '${OUTPUT_DIR}'"
