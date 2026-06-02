#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TASK_DIR="${SCRIPT_DIR}/tasks"

"${TASK_DIR}/01_build_bowtie2_index.sh"
"${TASK_DIR}/02_fastqc.sh" raw
"${TASK_DIR}/03_trim_reads.sh"
"${TASK_DIR}/02_fastqc.sh" trimmed
"${TASK_DIR}/04_align_reads.sh"
"${TASK_DIR}/05_mark_duplicates.sh"
"${TASK_DIR}/06_filter_bam.sh"
"${TASK_DIR}/07_collect_alignment_metrics.sh"
"${TASK_DIR}/08_collect_insert_size_metrics.sh"
"${TASK_DIR}/09_make_tracks.sh"
"${TASK_DIR}/10_call_peaks.sh"
"${TASK_DIR}/11_annotate_peaks.sh"
"${TASK_DIR}/12_motif_analysis.sh"
"${TASK_DIR}/13_differential_peaks.sh"
"${TASK_DIR}/14_plot_peak_statistics.sh"
"${TASK_DIR}/15_plot_genome_signal.sh"
"${TASK_DIR}/16_plot_all_peak_feature_distribution.sh"
"${TASK_DIR}/17_plot_peak.sh"
"${TASK_DIR}/18_align_spikein.sh"
"${TASK_DIR}/19_write_sample_sheet.sh"
"${TASK_DIR}/20_multiqc.sh"
