#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
GENOME_FASTA="${GENOME_FASTA:-${WORK_ROOT}/pipeline/database/mm39/Mus_musculus.GRCm39.dna.toplevel.fa}"
THREADS="${THREADS:-8}"
MOTIF_KMER_SIZE="${MOTIF_KMER_SIZE:-6}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
HOMER_IMAGE="${HOMER_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/homer}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)

for sample in "${SAMPLES[@]}"; do
  ${DOCKER_RUN} "${HOMER_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/motif/${sample}'; awk 'BEGIN{OFS=\"\\t\"} NF >= 3 {print \$1,\$2,\$3,\"peak_\"NR}' '${OUTPUT_DIR}/peaks/${sample}/${sample}.peaks.bed' > '${OUTPUT_DIR}/motif/${sample}/${sample}.homer_peaks.bed'; findMotifsGenome.pl '${OUTPUT_DIR}/motif/${sample}/${sample}.homer_peaks.bed' '${GENOME_FASTA}' '${OUTPUT_DIR}/motif/${sample}/homer' -size given -len ${MOTIF_KMER_SIZE},8,10,12 -p ${THREADS}; if [[ -f '${OUTPUT_DIR}/motif/${sample}/homer/knownResults.txt' ]]; then cp '${OUTPUT_DIR}/motif/${sample}/homer/knownResults.txt' '${OUTPUT_DIR}/motif/${sample}/${sample}.known_motifs.tsv'; fi; if [[ -f '${OUTPUT_DIR}/motif/${sample}/homer/homerResults.html' ]]; then cp '${OUTPUT_DIR}/motif/${sample}/homer/homerResults.html' '${OUTPUT_DIR}/motif/${sample}/${sample}.homer_motifs.html'; fi"
done
