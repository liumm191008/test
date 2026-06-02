#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-/home/data/vip01/work}"
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_ROOT}/bioproject/MJ20260515140/cuttag_results}"
ANNOTATION_GTF="${ANNOTATION_GTF:-${WORK_ROOT}/pipeline/database/mm39/Mus_musculus.GRCm39.115.gtf}"
THREADS="${THREADS:-8}"
DOCKER_RUN="${DOCKER_RUN:-docker run -i --rm --security-opt seccomp=unconfined -v ${WORK_ROOT}:${WORK_ROOT}}"
DEEPTOOLS_IMAGE="${DEEPTOOLS_IMAGE:-registry.cn-guangzhou.aliyuncs.com/origen/deeptools}"
SAMPLES=(con1 con2 con3 exo1 exo2 exo3 LPS1 LPS2 LPS3)

bams=(); bigwigs=()
for sample in "${SAMPLES[@]}"; do
  bams+=("${OUTPUT_DIR}/filtered/${sample}/${sample}.filtered.bam")
  bigwigs+=("${OUTPUT_DIR}/tracks/${sample}/${sample}.bw")
done
sample_labels="${SAMPLES[*]}"
bam_list="${bams[*]}"
bigwig_list="${bigwigs[*]}"
consensus_peaks="${OUTPUT_DIR}/differential_peaks/consensus_peaks.bed"

${DOCKER_RUN} "${DEEPTOOLS_IMAGE}" bash -c "set -euo pipefail; mkdir -p '${OUTPUT_DIR}/plots/summary'; plotFingerprint --bamfiles ${bam_list} --labels ${sample_labels} --plotFile '${OUTPUT_DIR}/plots/summary/fingerprint.svg' --outRawCounts '${OUTPUT_DIR}/plots/summary/fingerprint_counts.tsv' --numberOfProcessors ${THREADS} --plotTitle 'Fingerprints of different samples'; computeMatrix reference-point -S ${bigwig_list} -R '${consensus_peaks}' -a 3000 -b 3000 -bs 10 -p ${THREADS} -o '${OUTPUT_DIR}/plots/summary/peak_matrix.gz' --skipZeros; plotHeatmap -m '${OUTPUT_DIR}/plots/summary/peak_matrix.gz' -o '${OUTPUT_DIR}/plots/summary/peak_heatmap.png' --colorMap Reds Blues --whatToShow 'plot, heatmap and colorbar' --heatmapHeight 15 --heatmapWidth 4 -x 'peak distance(bp)' --refPointLabel 'center' --samplesLabel ${sample_labels}; plotProfile -m '${OUTPUT_DIR}/plots/summary/peak_matrix.gz' -out '${OUTPUT_DIR}/plots/summary/peak_profile.png' --samplesLabel ${sample_labels}; computeMatrix scale-regions -S ${bigwig_list} -R '${ANNOTATION_GTF}' -a 3000 -b 3000 --regionBodyLength 5000 -bs 10 -p ${THREADS} -o '${OUTPUT_DIR}/plots/summary/gene_matrix.gz' --skipZeros; plotHeatmap -m '${OUTPUT_DIR}/plots/summary/gene_matrix.gz' -o '${OUTPUT_DIR}/plots/summary/gene_heatmap.png' --colorMap Reds Blues --whatToShow 'plot, heatmap and colorbar' --heatmapHeight 15 --heatmapWidth 4 --samplesLabel ${sample_labels}; plotProfile -m '${OUTPUT_DIR}/plots/summary/gene_matrix.gz' -out '${OUTPUT_DIR}/plots/summary/gene_profile.png' --samplesLabel ${sample_labels}"
