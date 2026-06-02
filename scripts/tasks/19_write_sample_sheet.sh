#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
ANNOTATION_GTF="${ANNOTATION_GTF:-${WORK_ROOT}/pipeline/database/mm39/Mus_musculus.GRCm39.115.gtf}"
TXDB_PATH="${TXDB_PATH:-${WORK_ROOT}/pipeline/database/mm39/TxDb.Mmusculus.GRCm39.115.sqlite}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
SAMTOOLS_IMAGE="${SAMTOOLS_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/samtools}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)
SAMPLE_GROUPS=(con con con exo exo exo LPS LPS LPS)

${DOCKER_RUN} "${SAMTOOLS_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/metadata'; printf 'sample_id\tgroup\n' > '${OUTPUT_DIR}/metadata/sample_sheet.tsv'; printf 'key\tpath\nannotation_gtf\t${ANNOTATION_GTF}\ntxdb\t${TXDB_PATH}\n' > '${OUTPUT_DIR}/metadata/reference_paths.tsv'"
for i in "${!SAMPLES[@]}"; do
  ${DOCKER_RUN} "${SAMTOOLS_IMAGE}" bash -c "printf '%s\t%s\n' '${SAMPLES[$i]}' '${SAMPLE_GROUPS[$i]}' >> '${OUTPUT_DIR}/metadata/sample_sheet.tsv'"
done
