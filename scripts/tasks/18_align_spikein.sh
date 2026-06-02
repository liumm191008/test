#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
SPIKEIN_INDEX="${SPIKEIN_INDEX:-${WORK_ROOT}/pipeline/database/E.coli/MG1655}"
THREADS="${THREADS:-8}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
BOWTIE2_IMAGE="${BOWTIE2_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/bowtie2}"
SAMTOOLS_IMAGE="${SAMTOOLS_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/samtools}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)

[[ -z "${SPIKEIN_INDEX}" ]] && { echo "SPIKEIN_INDEX is empty; skip."; exit 0; }
for sample in "${SAMPLES[@]}"; do
  ${DOCKER_RUN} "${BOWTIE2_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/spikein/${sample}'; bowtie2 --end-to-end --very-sensitive -I 10 -X 700 --no-mixed --no-discordant --threads ${THREADS} -x '${SPIKEIN_INDEX}' -1 '${OUTPUT_DIR}/trimmed/${sample}/${sample}_val_1.fq.gz' -2 '${OUTPUT_DIR}/trimmed/${sample}/${sample}_val_2.fq.gz' -S '${OUTPUT_DIR}/spikein/${sample}/${sample}.spikein.sam' 2> '${OUTPUT_DIR}/spikein/${sample}/${sample}.spikein.bowtie2.log'"
  ${DOCKER_RUN} "${SAMTOOLS_IMAGE}" bash -c "set -euo pipefail; samtools sort -@ ${THREADS} -o '${OUTPUT_DIR}/spikein/${sample}/${sample}.spikein.sorted.bam' '${OUTPUT_DIR}/spikein/${sample}/${sample}.spikein.sam'; samtools index -@ ${THREADS} '${OUTPUT_DIR}/spikein/${sample}/${sample}.spikein.sorted.bam'; samtools flagstat -@ ${THREADS} '${OUTPUT_DIR}/spikein/${sample}/${sample}.spikein.sorted.bam' > '${OUTPUT_DIR}/spikein/${sample}/${sample}.spikein.flagstat.txt'; rm -f '${OUTPUT_DIR}/spikein/${sample}/${sample}.spikein.sam'"
done
