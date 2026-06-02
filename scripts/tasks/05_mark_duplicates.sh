#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
DEDUPLICATE="${DEDUPLICATE:-false}"
READ_GROUP_PLATFORM="${READ_GROUP_PLATFORM:-ILLUMINA}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
PICARD_IMAGE="${PICARD_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/picard}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)

for sample in "${SAMPLES[@]}"; do
  if [[ "${DEDUPLICATE}" == "true" ]]; then
    out_bam="${OUTPUT_DIR}/duplicates/${sample}/${sample}.dedup.bam"
    remove_duplicates=true
  else
    out_bam="${OUTPUT_DIR}/duplicates/${sample}/${sample}.markdup.bam"
    remove_duplicates=false
  fi
  ${DOCKER_RUN} "${PICARD_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/duplicates/${sample}'; read_group_bam='${OUTPUT_DIR}/duplicates/${sample}/${sample}.rg.bam'; picard AddOrReplaceReadGroups I='${OUTPUT_DIR}/alignment/${sample}/${sample}.sorted.bam' O=\"\${read_group_bam}\" RGID='${sample}' RGLB='cuttag_${sample}' RGPL='${READ_GROUP_PLATFORM}' RGPU='${sample}' RGSM='${sample}' SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT; picard MarkDuplicates I=\"\${read_group_bam}\" O='${out_bam}' M='${OUTPUT_DIR}/duplicates/${sample}/${sample}.duplication_metrics.txt' REMOVE_DUPLICATES='${remove_duplicates}' ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT; picard BuildBamIndex I='${out_bam}'; rm -f \"\${read_group_bam}\" \"\${read_group_bam}.bai\""
done
