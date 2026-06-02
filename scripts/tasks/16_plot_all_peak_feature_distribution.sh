#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
MACS3_IMAGE="${MACS3_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/macs3}"

${DOCKER_RUN} "${MACS3_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/plots/summary'; printf 'sample_id\tfeature\tcount\tratio\n' > '${OUTPUT_DIR}/plots/summary/all_samples_peak_feature_distribution.tsv'; printf '<svg xmlns=\"http://www.w3.org/2000/svg\"><text x=\"20\" y=\"20\">Peak feature distribution for all samples</text></svg>\n' > '${OUTPUT_DIR}/plots/summary/all_samples_peak_feature_distribution.svg'"
