#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
THREADS="${THREADS:-8}"
MAPQ="${MAPQ:-30}"
EXCLUDE_MITO="${EXCLUDE_MITO:-true}"
KEEP_PRIMARY_CONTIGS_ONLY="${KEEP_PRIMARY_CONTIGS_ONLY:-true}"
PRIMARY_CONTIG_REGEX="${PRIMARY_CONTIG_REGEX:-^([1-9]|1[0-9]|X|Y|MT)$}"
DEDUPLICATE="${DEDUPLICATE:-false}"
BLACKLIST_BED="${BLACKLIST_BED:-}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
SAMTOOLS_IMAGE="${SAMTOOLS_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/samtools}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)

for sample in "${SAMPLES[@]}"; do
  input_markdup="${OUTPUT_DIR}/duplicates/${sample}/${sample}.markdup.bam"
  [[ "${DEDUPLICATE}" == "true" ]] && input_markdup="${OUTPUT_DIR}/duplicates/${sample}/${sample}.dedup.bam"
  ${DOCKER_RUN} "${SAMTOOLS_IMAGE}" bash -c "set -euo pipefail
mkdir -p '${OUTPUT_DIR}/filtered/${sample}'
samtools view -@ ${THREADS} -b -q ${MAPQ} -f 2 -F 1804 '${input_markdup}' > '${OUTPUT_DIR}/filtered/${sample}/${sample}.mapq${MAPQ}.proper.bam'
samtools index -@ ${THREADS} '${OUTPUT_DIR}/filtered/${sample}/${sample}.mapq${MAPQ}.proper.bam'
input_bam='${OUTPUT_DIR}/filtered/${sample}/${sample}.mapq${MAPQ}.proper.bam'
samtools idxstats \"\${input_bam}\" | cut -f 1 | awk 'NF > 0 && \$1 != \"*\" {print \$1}' > '${OUTPUT_DIR}/filtered/${sample}/keep_contigs.txt'
if [[ '${KEEP_PRIMARY_CONTIGS_ONLY}' == 'true' ]]; then awk -v regex='${PRIMARY_CONTIG_REGEX}' '\$1 ~ regex {print \$1}' '${OUTPUT_DIR}/filtered/${sample}/keep_contigs.txt' > '${OUTPUT_DIR}/filtered/${sample}/primary_contigs.txt'; mv '${OUTPUT_DIR}/filtered/${sample}/primary_contigs.txt' '${OUTPUT_DIR}/filtered/${sample}/keep_contigs.txt'; fi
if [[ '${EXCLUDE_MITO}' == 'true' ]]; then grep -v -E '^(MT|M)$' '${OUTPUT_DIR}/filtered/${sample}/keep_contigs.txt' > '${OUTPUT_DIR}/filtered/${sample}/non_mito_contigs.txt'; mv '${OUTPUT_DIR}/filtered/${sample}/non_mito_contigs.txt' '${OUTPUT_DIR}/filtered/${sample}/keep_contigs.txt'; fi
if [[ -s '${OUTPUT_DIR}/filtered/${sample}/keep_contigs.txt' ]]; then samtools view -@ ${THREADS} -b \"\${input_bam}\" \$(cat '${OUTPUT_DIR}/filtered/${sample}/keep_contigs.txt') > '${OUTPUT_DIR}/filtered/${sample}/${sample}.contig_filtered.bam'; input_bam='${OUTPUT_DIR}/filtered/${sample}/${sample}.contig_filtered.bam'; fi
if [[ -n '${BLACKLIST_BED}' ]]; then samtools view -@ ${THREADS} -b -U '${OUTPUT_DIR}/filtered/${sample}/${sample}.filtered.bam' -L '${BLACKLIST_BED}' \"\${input_bam}\" > /dev/null; else cp \"\${input_bam}\" '${OUTPUT_DIR}/filtered/${sample}/${sample}.filtered.bam'; fi
samtools index -@ ${THREADS} '${OUTPUT_DIR}/filtered/${sample}/${sample}.filtered.bam'
samtools flagstat -@ ${THREADS} '${OUTPUT_DIR}/filtered/${sample}/${sample}.filtered.bam' > '${OUTPUT_DIR}/filtered/${sample}/${sample}.filtered.flagstat.txt'
samtools idxstats '${OUTPUT_DIR}/filtered/${sample}/${sample}.filtered.bam' > '${OUTPUT_DIR}/filtered/${sample}/${sample}.filtered.idxstats.txt'"
done
