version 1.0

## Paired-end CUT&Tag sample definition.
## All paths must be absolute host paths available through the docker_run mount.
struct CutTagSample {
  String sample_id
  String group
  String read1_path
  String read2_path
}

workflow CutTag {
  input {
    Array[CutTagSample] samples
    String genome_fasta_path = "/home/data/vip01/work/pipeline/database/mm39/Mus_musculus.GRCm39.dna.toplevel.fa"
    String annotation_gtf_path = "/home/data/vip01/work/pipeline/database/mm39/Mus_musculus.GRCm39.115.gtf"
    String txdb_path = "/home/data/vip01/work/pipeline/database/mm39/TxDb.Mmusculus.GRCm39.115.sqlite"
    String species = "mouse"
    String bowtie2_index_path = "/home/data/vip01/work/pipeline/database/mm39/mm39"
    String blacklist_bed_path = ""
    String spikein_index_path = "/home/data/vip01/work/pipeline/database/E.coli/MG1655"
    String output_dir = "/home/data/vip01/work/cuttag_results"
    Int threads = 8
    Int mapq = 30
    Boolean exclude_mito = true
    Boolean keep_primary_contigs_only = true
    String primary_contig_regex = "^([1-9]|1[0-9]|X|Y|MT)$"
    Boolean deduplicate = false
    String read_group_platform = "ILLUMINA"
    Boolean call_broad_peaks = false
    Float macs3_qvalue = 0.01
    String bigwig_normalization = "CPM"
    Int effective_genome_size = 2652783500
    Int bigwig_bin_size = 10
    String enrichment_annotation_path = ""
    String enrichment_gene_id_type = "ENTREZID"
    Int motif_kmer_size = 6
    Float diff_peak_pvalue = 0.05
    Float diff_peak_min_fraction_delta = 0.5

    String docker_run = "docker run -i --rm --security-opt seccomp=unconfined  -v /home/data/vip01/work:/home/data/vip01/work"
    String fastqc_image = "registry.cn-guangzhou.aliyuncs.com/origen/fastqc"
    String trim_galore_image = "registry.cn-guangzhou.aliyuncs.com/origen/trim-galore"
    String bowtie2_image = "registry.cn-guangzhou.aliyuncs.com/origen/bowtie2"
    String samtools_image = "registry.cn-guangzhou.aliyuncs.com/origen/samtools"
    String picard_image = "registry.cn-guangzhou.aliyuncs.com/origen/picard"
    String deeptools_image = "registry.cn-guangzhou.aliyuncs.com/origen/deeptools"
    String macs3_image = "registry.cn-guangzhou.aliyuncs.com/origen/macs3"
    String bioconductor_image = "registry.cn-guangzhou.aliyuncs.com/origen/bioconductor"
    String homer_image = "registry.cn-guangzhou.aliyuncs.com/origen/homer"
    String multiqc_image = "registry.cn-guangzhou.aliyuncs.com/origen/multiqc"
  }

  if (bowtie2_index_path == "") {
    call BuildBowtie2Index {
      input:
        genome_fasta_path = genome_fasta_path,
        output_dir = output_dir,
        threads = threads,
        docker_run = docker_run,
        image = bowtie2_image
    }
  }

  String resolved_bowtie2_index_path = if bowtie2_index_path != "" then bowtie2_index_path else select_first([BuildBowtie2Index.index_prefix_path])

  scatter (sample in samples) {
    String current_sample_id = sample.sample_id
    String current_sample_group = sample.group

    call FastQc as RawFastQc {
      input:
        sample_id = sample.sample_id,
        read1_path = sample.read1_path,
        read2_path = sample.read2_path,
        output_dir = output_dir,
        threads = threads,
        docker_run = docker_run,
        image = fastqc_image
    }

    call TrimReads {
      input:
        sample_id = sample.sample_id,
        read1_path = sample.read1_path,
        read2_path = sample.read2_path,
        output_dir = output_dir,
        threads = threads,
        docker_run = docker_run,
        image = trim_galore_image
    }

    call FastQc as TrimmedFastQc {
      input:
        sample_id = sample.sample_id + ".trimmed",
        read1_path = TrimReads.trimmed_read1_path,
        read2_path = TrimReads.trimmed_read2_path,
        output_dir = output_dir,
        threads = threads,
        docker_run = docker_run,
        image = fastqc_image
    }

    call AlignReads {
      input:
        sample_id = sample.sample_id,
        read1_path = TrimReads.trimmed_read1_path,
        read2_path = TrimReads.trimmed_read2_path,
        bowtie2_index_path = resolved_bowtie2_index_path,
        output_dir = output_dir,
        threads = threads,
        docker_run = docker_run,
        bowtie2_image = bowtie2_image,
        samtools_image = samtools_image
    }

    call MarkDuplicates {
      input:
        sample_id = sample.sample_id,
        sorted_bam_path = AlignReads.sorted_bam_path,
        output_dir = output_dir,
        read_group_platform = read_group_platform,
        deduplicate = deduplicate,
        docker_run = docker_run,
        image = picard_image
    }

    call FilterBam {
      input:
        sample_id = sample.sample_id,
        bam_path = MarkDuplicates.output_bam_path,
        output_dir = output_dir,
        threads = threads,
        mapq = mapq,
        exclude_mito = exclude_mito,
        keep_primary_contigs_only = keep_primary_contigs_only,
        primary_contig_regex = primary_contig_regex,
        blacklist_bed_path = blacklist_bed_path,
        docker_run = docker_run,
        image = samtools_image
    }

    call CollectAlignmentMetrics {
      input:
        sample_id = sample.sample_id,
        filtered_bam_path = FilterBam.filtered_bam_path,
        genome_fasta_path = genome_fasta_path,
        output_dir = output_dir,
        docker_run = docker_run,
        image = picard_image
    }

    call CollectInsertSizeMetrics {
      input:
        sample_id = sample.sample_id,
        filtered_bam_path = FilterBam.filtered_bam_path,
        output_dir = output_dir,
        docker_run = docker_run,
        image = picard_image
    }

    call MakeTracks {
      input:
        sample_id = sample.sample_id,
        filtered_bam_path = FilterBam.filtered_bam_path,
        output_dir = output_dir,
        threads = threads,
        bigwig_normalization = bigwig_normalization,
        effective_genome_size = effective_genome_size,
        bigwig_bin_size = bigwig_bin_size,
        docker_run = docker_run,
        samtools_image = samtools_image,
        deeptools_image = deeptools_image
    }

    call CallPeaks {
      input:
        sample_id = sample.sample_id,
        filtered_bam_path = FilterBam.filtered_bam_path,
        output_dir = output_dir,
        call_broad_peaks = call_broad_peaks,
        macs3_qvalue = macs3_qvalue,
        species = species,
        docker_run = docker_run,
        image = macs3_image
    }

    call AnnotatePeaks {
      input:
        sample_id = sample.sample_id,
        peaks_path = CallPeaks.peaks_path,
        species = species,
        enrichment_gene_id_type = enrichment_gene_id_type,
        output_dir = output_dir,
        docker_run = docker_run,
        image = bioconductor_image
    }

    call MotifAnalysis {
      input:
        sample_id = sample.sample_id,
        peaks_path = CallPeaks.peaks_path,
        genome_fasta_path = genome_fasta_path,
        output_dir = output_dir,
        motif_kmer_size = motif_kmer_size,
        threads = threads,
        docker_run = docker_run,
        image = homer_image
    }

    call PlotPeakStatistics {
      input:
        sample_id = sample.sample_id,
        peaks_path = CallPeaks.peaks_path,
        annotation_summary_path = AnnotatePeaks.annotation_summary_path,
        output_dir = output_dir,
        docker_run = docker_run,
        image = macs3_image
    }

    call PlotGenomeSignal {
      input:
        sample_id = sample.sample_id,
        filtered_bam_path = FilterBam.filtered_bam_path,
        output_dir = output_dir,
        threads = threads,
        keep_primary_contigs_only = keep_primary_contigs_only,
        primary_contig_regex = primary_contig_regex,
        docker_run = docker_run,
        image = samtools_image
    }

    if (spikein_index_path != "") {
      call AlignSpikeIn {
        input:
          sample_id = sample.sample_id,
          read1_path = TrimReads.trimmed_read1_path,
          read2_path = TrimReads.trimmed_read2_path,
          spikein_index_path = spikein_index_path,
          output_dir = output_dir,
          threads = threads,
          docker_run = docker_run,
          bowtie2_image = bowtie2_image,
          samtools_image = samtools_image
      }
    }
  }

  call DifferentialPeaks {
    input:
      peak_paths = CallPeaks.peaks_path,
      bam_paths = FilterBam.filtered_bam_path,
      sample_ids = current_sample_id,
      sample_groups = current_sample_group,
      species = species,
      txdb_path = txdb_path,
      enrichment_gene_id_type = enrichment_gene_id_type,
      genome_fasta_path = genome_fasta_path,
      output_dir = output_dir,
      diff_peak_pvalue = diff_peak_pvalue,
      threads = threads,
      docker_run = docker_run,
      image = bioconductor_image,
      homer_image = homer_image
  }

  call PlotAllPeakFeatureDistribution {
    input:
      sample_ids = current_sample_id,
      annotation_summary_paths = AnnotatePeaks.annotation_summary_path,
      output_dir = output_dir,
      docker_run = docker_run,
      image = macs3_image
  }

  call PlotPeak {
    input:
      sample_ids = current_sample_id,
      filtered_bam_paths = FilterBam.filtered_bam_path,
      bigwig_paths = MakeTracks.bigwig_path,
      consensus_peaks_path = DifferentialPeaks.consensus_peak_set_path,
      annotation_gtf_path = annotation_gtf_path,
      output_dir = output_dir,
      threads = threads,
      docker_run = docker_run,
      image = deeptools_image
  }

  call WriteSampleSheet {
    input:
      sample_ids = current_sample_id,
      sample_groups = current_sample_group,
      annotation_gtf_path = annotation_gtf_path,
      output_dir = output_dir,
      docker_run = docker_run,
      image = samtools_image
  }

  call MultiQc {
    input:
      output_dir = output_dir,
      docker_run = docker_run,
      image = multiqc_image
  }

  output {
    String bowtie2_index = resolved_bowtie2_index_path
    String reference_fasta = genome_fasta_path
    String annotation_gtf = annotation_gtf_path
    String txdb = txdb_path
    Array[String] raw_fastqc_dirs = RawFastQc.fastqc_dir_path
    Array[String] trimmed_fastqc_dirs = TrimmedFastQc.fastqc_dir_path
    Array[String] trimmed_read1 = TrimReads.trimmed_read1_path
    Array[String] trimmed_read2 = TrimReads.trimmed_read2_path
    Array[String] sorted_bams = AlignReads.sorted_bam_path
    Array[String] alignment_logs = AlignReads.alignment_log_path
    Array[String] duplicate_marked_bams = MarkDuplicates.output_bam_path
    Array[String] duplication_metrics = MarkDuplicates.duplication_metrics_path
    Array[String] filtered_bams = FilterBam.filtered_bam_path
    Array[String] filtered_bam_indexes = FilterBam.filtered_bam_index_path
    Array[String] filter_stats = FilterBam.flagstat_path
    Array[String] picard_alignment_metrics = CollectAlignmentMetrics.alignment_metrics_path
    Array[String] insert_size_metrics = CollectInsertSizeMetrics.insert_size_metrics_path
    Array[String] insert_size_plots = CollectInsertSizeMetrics.insert_size_histogram_path
    Array[String] fragments = MakeTracks.fragments_bed_path
    Array[String] bigwigs = MakeTracks.bigwig_path
    Array[String] peak_files = CallPeaks.peaks_path
    Array[String] peak_summaries = CallPeaks.peak_summary_path
    Array[String] annotated_peaks = AnnotatePeaks.annotated_peaks_path
    Array[String] peak_annotation_summaries = AnnotatePeaks.annotation_summary_path
    Array[String] peak_annotation_stats = AnnotatePeaks.peak_annotation_stats_path
    Array[String] peak_enrichment_gene_lists = AnnotatePeaks.enrichment_genes_path
    Array[String] peak_enrichment_dirs = AnnotatePeaks.peak_enrichment_dir_path
    Array[String] annotation_pie_charts = AnnotatePeaks.annotation_pie_plot_path
    Array[String] tss_distance_plots = AnnotatePeaks.tss_distance_plot_path
    Array[String] genomic_feature_distribution_plots = AnnotatePeaks.genomic_feature_distribution_plot_path
    Array[String] peak_feature_distribution_plots = PlotPeakStatistics.feature_distribution_plot_path
    Array[String] peak_length_distribution_plots = PlotPeakStatistics.length_distribution_plot_path
    Array[String] genome_signal_distribution_plots = PlotGenomeSignal.genome_signal_plot_path
    String all_peak_feature_distribution_plot = PlotAllPeakFeatureDistribution.feature_distribution_plot_path
    String all_peak_feature_distribution_table = PlotAllPeakFeatureDistribution.feature_distribution_table_path
    String fingerprint_plot = PlotPeak.fingerprint_plot_path
    String fingerprint_counts = PlotPeak.fingerprint_counts_path
    String peak_reference_point_matrix = PlotPeak.peak_reference_point_matrix_path
    String peak_reference_point_heatmap = PlotPeak.peak_reference_point_heatmap_path
    String peak_reference_point_profile = PlotPeak.peak_reference_point_profile_path
    String gene_body_matrix = PlotPeak.gene_body_matrix_path
    String gene_body_heatmap = PlotPeak.gene_body_heatmap_path
    String gene_body_profile = PlotPeak.gene_body_profile_path
    Array[String] motif_tables = MotifAnalysis.motif_table_path
    Array[String] motif_plots = MotifAnalysis.motif_plot_path
    String differential_peak_summary = DifferentialPeaks.diff_summary_path
    String differential_peak_results_dir = DifferentialPeaks.diff_results_dir
    String differential_peak_annotation_dir = DifferentialPeaks.diff_annotation_dir
    String differential_consensus_peak_set = DifferentialPeaks.consensus_peak_set_path
    String differential_counts_matrix = DifferentialPeaks.counts_matrix_path
    String differential_normalized_matrix = DifferentialPeaks.normalized_matrix_path
    String differential_diffbind_rds = DifferentialPeaks.diffbind_rds_path
    String differential_pca_plot = DifferentialPeaks.pca_plot_path
    String differential_correlation_heatmap = DifferentialPeaks.correlation_heatmap_path
    String differential_peak_enrichment_dir = DifferentialPeaks.diff_enrichment_dir
    String differential_peak_motif_dir = DifferentialPeaks.diff_motif_dir
    String differential_peak_visualization_dir = DifferentialPeaks.diff_visualization_dir
    String diffbind_report = DifferentialPeaks.diffbind_report_path
    Array[String] spikein_bams = select_all(AlignSpikeIn.spikein_bam_path)
    Array[String] spikein_logs = select_all(AlignSpikeIn.alignment_log_path)
    String sample_sheet = WriteSampleSheet.sample_sheet_path
    String reference_paths = WriteSampleSheet.reference_paths_path
    String multiqc_report = MultiQc.report_path
  }
}

task BuildBowtie2Index {
  input {
    String genome_fasta_path
    String output_dir
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/bowtie2_index"
      bowtie2-build --threads ~{threads} "~{genome_fasta_path}" "~{output_dir}/bowtie2_index/genome"
    '
  >>>

  output {
    String index_prefix_path = output_dir + "/bowtie2_index/genome"
  }
}

task FastQc {
  input {
    String sample_id
    String read1_path
    String read2_path
    String output_dir
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/fastqc/~{sample_id}"
      fastqc --threads ~{threads} --outdir "~{output_dir}/fastqc/~{sample_id}" "~{read1_path}" "~{read2_path}"
      for zip_file in "~{output_dir}/fastqc/~{sample_id}"/*.zip; do
        [[ -f "${zip_file}" ]] && unzip -o "${zip_file}" -d "~{output_dir}/fastqc/~{sample_id}"
      done
    '
  >>>

  output {
    String fastqc_dir_path = output_dir + "/fastqc/" + sample_id
  }
}

task TrimReads {
  input {
    String sample_id
    String read1_path
    String read2_path
    String output_dir
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/trimmed/~{sample_id}"
      trim_galore --paired --cores ~{threads} --gzip --basename "~{sample_id}" --output_dir "~{output_dir}/trimmed/~{sample_id}" "~{read1_path}" "~{read2_path}"
      test -s "~{output_dir}/trimmed/~{sample_id}/~{sample_id}_val_1.fq.gz"
      test -s "~{output_dir}/trimmed/~{sample_id}/~{sample_id}_val_2.fq.gz"
    '
  >>>

  output {
    String trimmed_read1_path = output_dir + "/trimmed/" + sample_id + "/" + sample_id + "_val_1.fq.gz"
    String trimmed_read2_path = output_dir + "/trimmed/" + sample_id + "/" + sample_id + "_val_2.fq.gz"
    String trim_dir_path = output_dir + "/trimmed/" + sample_id
  }
}

task AlignReads {
  input {
    String sample_id
    String read1_path
    String read2_path
    String bowtie2_index_path
    String output_dir
    Int threads
    String docker_run
    String bowtie2_image
    String samtools_image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{bowtie2_image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/alignment/~{sample_id}"
      bowtie2 \
        --end-to-end \
        --very-sensitive \
        -I 10 \
        -X 700 \
        --no-mixed \
        --no-discordant \
        --threads ~{threads} \
        -x "~{bowtie2_index_path}" \
        -1 "~{read1_path}" \
        -2 "~{read2_path}" \
        -S "~{output_dir}/alignment/~{sample_id}/~{sample_id}.sam" \
        2> "~{output_dir}/alignment/~{sample_id}/~{sample_id}.bowtie2.log"
    '
    ~{docker_run} ~{samtools_image} bash -c '
      set -euo pipefail
      samtools sort -@ ~{threads} -o "~{output_dir}/alignment/~{sample_id}/~{sample_id}.sorted.bam" "~{output_dir}/alignment/~{sample_id}/~{sample_id}.sam"
      samtools index -@ ~{threads} "~{output_dir}/alignment/~{sample_id}/~{sample_id}.sorted.bam"
      rm -f "~{output_dir}/alignment/~{sample_id}/~{sample_id}.sam"
    '
  >>>

  output {
    String sorted_bam_path = output_dir + "/alignment/" + sample_id + "/" + sample_id + ".sorted.bam"
    String sorted_bam_index_path = output_dir + "/alignment/" + sample_id + "/" + sample_id + ".sorted.bam.bai"
    String alignment_log_path = output_dir + "/alignment/" + sample_id + "/" + sample_id + ".bowtie2.log"
  }
}


task MarkDuplicates {
  input {
    String sample_id
    String sorted_bam_path
    String output_dir
    String read_group_platform
    Boolean deduplicate
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/duplicates/~{sample_id}"
      read_group_bam="~{output_dir}/duplicates/~{sample_id}/~{sample_id}.rg.bam"
      picard AddOrReplaceReadGroups \
        I="~{sorted_bam_path}" \
        O="${read_group_bam}" \
        RGID="~{sample_id}" \
        RGLB="cuttag_~{sample_id}" \
        RGPL="~{read_group_platform}" \
        RGPU="~{sample_id}" \
        RGSM="~{sample_id}" \
        SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=SILENT
      if [[ "~{deduplicate}" == "true" ]]; then
        output_bam="~{output_dir}/duplicates/~{sample_id}/~{sample_id}.dedup.bam"
        remove_duplicates=true
      else
        output_bam="~{output_dir}/duplicates/~{sample_id}/~{sample_id}.markdup.bam"
        remove_duplicates=false
      fi
      picard MarkDuplicates \
        I="${read_group_bam}" \
        O="${output_bam}" \
        M="~{output_dir}/duplicates/~{sample_id}/~{sample_id}.duplication_metrics.txt" \
        REMOVE_DUPLICATES="${remove_duplicates}" \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT
      picard BuildBamIndex I="${output_bam}"
      rm -f "${read_group_bam}" "${read_group_bam}.bai"
    '
  >>>

  output {
    String output_bam_path = if deduplicate then output_dir + "/duplicates/" + sample_id + "/" + sample_id + ".dedup.bam" else output_dir + "/duplicates/" + sample_id + "/" + sample_id + ".markdup.bam"
    String output_bam_index_path = output_bam_path + ".bai"
    String duplication_metrics_path = output_dir + "/duplicates/" + sample_id + "/" + sample_id + ".duplication_metrics.txt"
  }
}

task FilterBam {
  input {
    String sample_id
    String bam_path
    String output_dir
    Int threads
    Int mapq
    Boolean exclude_mito
    Boolean keep_primary_contigs_only
    String primary_contig_regex
    String blacklist_bed_path
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/filtered/~{sample_id}"
      samtools view \
        -@ ~{threads} \
        -b \
        -q ~{mapq} \
        -f 2 \
        -F 1804 \
        "~{bam_path}" \
      > "~{output_dir}/filtered/~{sample_id}/~{sample_id}.mapq~{mapq}.proper.bam"
      samtools index -@ ~{threads} "~{output_dir}/filtered/~{sample_id}/~{sample_id}.mapq~{mapq}.proper.bam"
      input_bam="~{output_dir}/filtered/~{sample_id}/~{sample_id}.mapq~{mapq}.proper.bam"

      samtools idxstats "${input_bam}" | cut -f 1 | awk 'NF > 0 && $1 != "*" {print $1}' > "~{output_dir}/filtered/~{sample_id}/all_contigs.txt"
      cp "~{output_dir}/filtered/~{sample_id}/all_contigs.txt" "~{output_dir}/filtered/~{sample_id}/keep_contigs.txt"
      if [[ "~{keep_primary_contigs_only}" == "true" ]]; then
        awk -v regex="~{primary_contig_regex}" '$1 ~ regex {print $1}' "~{output_dir}/filtered/~{sample_id}/keep_contigs.txt" > "~{output_dir}/filtered/~{sample_id}/primary_contigs.txt"
        mv "~{output_dir}/filtered/~{sample_id}/primary_contigs.txt" "~{output_dir}/filtered/~{sample_id}/keep_contigs.txt"
      fi
      if [[ "~{exclude_mito}" == "true" ]]; then
        grep -v -E "^(MT|M)$" "~{output_dir}/filtered/~{sample_id}/keep_contigs.txt" > "~{output_dir}/filtered/~{sample_id}/non_mito_contigs.txt"
        mv "~{output_dir}/filtered/~{sample_id}/non_mito_contigs.txt" "~{output_dir}/filtered/~{sample_id}/keep_contigs.txt"
      fi
      if [[ -s "~{output_dir}/filtered/~{sample_id}/keep_contigs.txt" ]]; then
        samtools view -@ ~{threads} -b "${input_bam}" $(cat "~{output_dir}/filtered/~{sample_id}/keep_contigs.txt") > "~{output_dir}/filtered/~{sample_id}/~{sample_id}.contig_filtered.bam"
        samtools index -@ ~{threads} "~{output_dir}/filtered/~{sample_id}/~{sample_id}.contig_filtered.bam"
        input_bam="~{output_dir}/filtered/~{sample_id}/~{sample_id}.contig_filtered.bam"
      fi

      if [[ -n "~{blacklist_bed_path}" ]]; then
        samtools view -@ ~{threads} -b -U "~{output_dir}/filtered/~{sample_id}/~{sample_id}.filtered.bam" -L "~{blacklist_bed_path}" "${input_bam}" > /dev/null
      else
        cp "${input_bam}" "~{output_dir}/filtered/~{sample_id}/~{sample_id}.filtered.bam"
      fi

      samtools index -@ ~{threads} "~{output_dir}/filtered/~{sample_id}/~{sample_id}.filtered.bam"
      samtools flagstat -@ ~{threads} "~{output_dir}/filtered/~{sample_id}/~{sample_id}.filtered.bam" > "~{output_dir}/filtered/~{sample_id}/~{sample_id}.filtered.flagstat.txt"
      samtools idxstats "~{output_dir}/filtered/~{sample_id}/~{sample_id}.filtered.bam" > "~{output_dir}/filtered/~{sample_id}/~{sample_id}.filtered.idxstats.txt"
    '
  >>>

  output {
    String filtered_bam_path = output_dir + "/filtered/" + sample_id + "/" + sample_id + ".filtered.bam"
    String filtered_bam_index_path = output_dir + "/filtered/" + sample_id + "/" + sample_id + ".filtered.bam.bai"
    String flagstat_path = output_dir + "/filtered/" + sample_id + "/" + sample_id + ".filtered.flagstat.txt"
    String idxstats_path = output_dir + "/filtered/" + sample_id + "/" + sample_id + ".filtered.idxstats.txt"
  }
}

task CollectAlignmentMetrics {
  input {
    String sample_id
    String filtered_bam_path
    String genome_fasta_path
    String output_dir
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/picard/~{sample_id}"
      picard CollectAlignmentSummaryMetrics \
        I="~{filtered_bam_path}" \
        R="~{genome_fasta_path}" \
        O="~{output_dir}/picard/~{sample_id}/~{sample_id}.alignment_metrics.txt"
    '
  >>>

  output {
    String alignment_metrics_path = output_dir + "/picard/" + sample_id + "/" + sample_id + ".alignment_metrics.txt"
  }
}


task CollectInsertSizeMetrics {
  input {
    String sample_id
    String filtered_bam_path
    String output_dir
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/insert_size/~{sample_id}"
      picard CollectInsertSizeMetrics \
        I="~{filtered_bam_path}" \
        O="~{output_dir}/insert_size/~{sample_id}/~{sample_id}.insert_size_metrics.txt" \
        H="~{output_dir}/insert_size/~{sample_id}/~{sample_id}.insert_size_histogram.pdf" \
        M=0.5 \
        VALIDATION_STRINGENCY=SILENT
    '
  >>>

  output {
    String insert_size_metrics_path = output_dir + "/insert_size/" + sample_id + "/" + sample_id + ".insert_size_metrics.txt"
    String insert_size_histogram_path = output_dir + "/insert_size/" + sample_id + "/" + sample_id + ".insert_size_histogram.pdf"
  }
}

task MakeTracks {
  input {
    String sample_id
    String filtered_bam_path
    String output_dir
    Int threads
    String bigwig_normalization
    Int effective_genome_size
    Int bigwig_bin_size
    String docker_run
    String samtools_image
    String deeptools_image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{samtools_image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/tracks/~{sample_id}"
      samtools sort -@ ~{threads} -n -o "~{output_dir}/tracks/~{sample_id}/~{sample_id}.name_sorted.bam" "~{filtered_bam_path}"
      samtools view -@ ~{threads} -f 2 -F 1804 "~{output_dir}/tracks/~{sample_id}/~{sample_id}.name_sorted.bam" \
      | awk '\''BEGIN{OFS="\t"} $9 > 0 {start=$4-1; end=start+$9; if (end > start) print $3,start,end,".",0,"."}'\'' \
      | sort -k1,1 -k2,2n > "~{output_dir}/tracks/~{sample_id}/~{sample_id}.fragments.bed"
    '

    if [[ "~{bigwig_normalization}" == "RPGC" ]]; then
      if [[ ~{effective_genome_size} -le 0 ]]; then
        echo "effective_genome_size must be greater than 0 when bigwig_normalization is RPGC" >&2
        exit 1
      fi
      ~{docker_run} ~{deeptools_image} bash -c '
        set -euo pipefail
        bamCoverage \
          --bam "~{filtered_bam_path}" \
          --outFileName "~{output_dir}/tracks/~{sample_id}/~{sample_id}.bw" \
          --outFileFormat bigwig \
          --normalizeUsing RPGC \
          --effectiveGenomeSize ~{effective_genome_size} \
          --binSize ~{bigwig_bin_size} \
          --extendReads \
          --numberOfProcessors ~{threads}
      '
    else
      ~{docker_run} ~{deeptools_image} bash -c '
        set -euo pipefail
        bamCoverage \
          --bam "~{filtered_bam_path}" \
          --outFileName "~{output_dir}/tracks/~{sample_id}/~{sample_id}.bw" \
          --outFileFormat bigwig \
          --normalizeUsing "~{bigwig_normalization}" \
          --binSize ~{bigwig_bin_size} \
          --extendReads \
          --numberOfProcessors ~{threads}
      '
    fi
  >>>

  output {
    String fragments_bed_path = output_dir + "/tracks/" + sample_id + "/" + sample_id + ".fragments.bed"
    String bigwig_path = output_dir + "/tracks/" + sample_id + "/" + sample_id + ".bw"
  }
}

task CallPeaks {
  input {
    String sample_id
    String filtered_bam_path
    String output_dir
    Boolean call_broad_peaks
    Float macs3_qvalue
    String species
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/peaks/~{sample_id}"
      case "~{species}" in
        human|hsa|hs) macs3_genome_size="hs" ;;
        *) macs3_genome_size="mm" ;;
      esac
      broad_flag=""
      if [[ "~{call_broad_peaks}" == "true" ]]; then
        broad_flag="--broad --broad-cutoff ~{macs3_qvalue}"
      fi

      macs3 callpeak \
        -t "~{filtered_bam_path}" \
        -f BAMPE \
        -g "${macs3_genome_size}" \
        -n "~{sample_id}" \
        --outdir "~{output_dir}/peaks/~{sample_id}" \
        -q ~{macs3_qvalue} \
        --keep-dup all \
        ${broad_flag}

      if [[ -f "~{output_dir}/peaks/~{sample_id}/~{sample_id}_peaks.broadPeak" ]]; then
        peak_file="~{output_dir}/peaks/~{sample_id}/~{sample_id}_peaks.broadPeak"
      else
        peak_file="~{output_dir}/peaks/~{sample_id}/~{sample_id}_peaks.narrowPeak"
      fi
      awk '\''BEGIN{OFS="\t"} {count += 1; bp += ($3 - $2)} END{print "sample","peak_count","peak_bp"; print "~{sample_id}",count+0,bp+0}'\'' "${peak_file}" > "~{output_dir}/peaks/~{sample_id}/~{sample_id}.peak_summary.tsv"
      cp "${peak_file}" "~{output_dir}/peaks/~{sample_id}/~{sample_id}.peaks.bed"
    '
  >>>

  output {
    String peaks_path = output_dir + "/peaks/" + sample_id + "/" + sample_id + ".peaks.bed"
    String peak_summary_path = output_dir + "/peaks/" + sample_id + "/" + sample_id + ".peak_summary.tsv"
    String peaks_dir = output_dir + "/peaks/" + sample_id
  }
}

task AnnotatePeaks {
  input {
    String sample_id
    String peaks_path
    String species
    String enrichment_gene_id_type
    String output_dir
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      test -s "/kegg_data/peak_annotation.R"
      mkdir -p "~{output_dir}/annotation/~{sample_id}"
      Rscript "/kegg_data/peak_annotation.R" \
        --peak "~{peaks_path}" \
        --sample "~{sample_id}" \
        --species "~{species}" \
        --annotation-dir "~{output_dir}/annotation/~{sample_id}"
      test -s "/kegg_data/gene_cluster_enrich.R"
      mkdir -p "~{output_dir}/annotation/~{sample_id}/enrichment"
      Rscript "/kegg_data/gene_cluster_enrich.R" \
        "~{output_dir}/annotation/~{sample_id}/~{sample_id}.enrichment_genes.tsv" \
        "~{enrichment_gene_id_type}" \
        "~{species}" \
        "~{output_dir}/annotation/~{sample_id}/enrichment" \
        "~{sample_id}"
    '
  >>>

  output {
    String annotated_peaks_path = output_dir + "/annotation/" + sample_id + "/" + sample_id + ".peaks.annotated.tsv"
    String peak_genes_path = output_dir + "/annotation/" + sample_id + "/" + sample_id + ".peak_genes.txt"
    String enrichment_genes_path = output_dir + "/annotation/" + sample_id + "/" + sample_id + ".enrichment_genes.tsv"
    String peak_enrichment_dir_path = output_dir + "/annotation/" + sample_id + "/enrichment"
    String annotation_summary_path = output_dir + "/annotation/" + sample_id + "/" + sample_id + ".annotation_summary.tsv"
    String peak_annotation_stats_path = output_dir + "/annotation/" + sample_id + "/" + sample_id + ".peak_annotation_stats.tsv"
    String annotation_pie_plot_path = output_dir + "/annotation/" + sample_id + "/" + sample_id + ".annotation_pie.pdf"
    String tss_distance_plot_path = output_dir + "/annotation/" + sample_id + "/" + sample_id + ".tss_distance.pdf"
    String genomic_feature_distribution_plot_path = output_dir + "/annotation/" + sample_id + "/" + sample_id + ".genomic_feature_distribution.pdf"
  }
}
task MotifAnalysis {
  input {
    String sample_id
    String peaks_path
    String genome_fasta_path
    String output_dir
    Int motif_kmer_size
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/motif/~{sample_id}"
      awk '\''BEGIN{OFS="\t"} NF >= 3 {print $1,$2,$3,"peak_"NR}'\'' "~{peaks_path}" > "~{output_dir}/motif/~{sample_id}/~{sample_id}.homer_peaks.bed"
      findMotifsGenome.pl \
        "~{output_dir}/motif/~{sample_id}/~{sample_id}.homer_peaks.bed" \
        "~{genome_fasta_path}" \
        "~{output_dir}/motif/~{sample_id}/homer" \
        -size given \
        -len ~{motif_kmer_size},8,10,12 \
        -p ~{threads}
      if [[ -f "~{output_dir}/motif/~{sample_id}/homer/knownResults.txt" ]]; then
        cp "~{output_dir}/motif/~{sample_id}/homer/knownResults.txt" "~{output_dir}/motif/~{sample_id}/~{sample_id}.known_motifs.tsv"
      else
        printf "status\tmessage\nmissing\tknownResults.txt not generated\n" > "~{output_dir}/motif/~{sample_id}/~{sample_id}.known_motifs.tsv"
      fi
      if [[ -f "~{output_dir}/motif/~{sample_id}/homer/homerResults.html" ]]; then
        cp "~{output_dir}/motif/~{sample_id}/homer/homerResults.html" "~{output_dir}/motif/~{sample_id}/~{sample_id}.homer_motifs.html"
      else
        printf "<html><body>HOMER motif report was not generated.</body></html>\n" > "~{output_dir}/motif/~{sample_id}/~{sample_id}.homer_motifs.html"
      fi
    '
  >>>

  output {
    String homer_peak_bed_path = output_dir + "/motif/" + sample_id + "/" + sample_id + ".homer_peaks.bed"
    String homer_results_dir = output_dir + "/motif/" + sample_id + "/homer"
    String motif_table_path = output_dir + "/motif/" + sample_id + "/" + sample_id + ".known_motifs.tsv"
    String motif_plot_path = output_dir + "/motif/" + sample_id + "/" + sample_id + ".homer_motifs.html"
  }
}

task DifferentialPeaks {
  input {
    Array[String] peak_paths
    Array[String] bam_paths
    Array[String] sample_ids
    Array[String] sample_groups
    String species
    String txdb_path
    String enrichment_gene_id_type
    String genome_fasta_path
    String output_dir
    Float diff_peak_pvalue
    Int threads
    String docker_run
    String image
    String homer_image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/differential_peaks" "~{output_dir}/differential_peaks/annotation" "~{output_dir}/differential_peaks/enrichment" "~{output_dir}/differential_peaks/motif"
      Rscript - <<'\''RSCRIPT'\''
      suppressPackageStartupMessages(library(DiffBind))
      suppressPackageStartupMessages(library(ChIPseeker))
      peak_paths <- c("~{sep='","' peak_paths}")
      bam_paths <- c("~{sep='","' bam_paths}")
      sample_ids <- c("~{sep='","' sample_ids}")
      sample_groups <- c("~{sep='","' sample_groups}")
      expand_sample_groups <- function(ids, groups) {
        if (length(groups) == length(ids)) {
          return(groups)
        }
        unique_groups <- unique(groups[nzchar(groups)])
        if (length(unique_groups) > 0) {
          assigned <- vapply(ids, function(id) {
            hits <- unique_groups[startsWith(tolower(id), tolower(unique_groups))]
            if (length(hits) == 1) hits[1] else NA_character_
          }, character(1))
          if (all(!is.na(assigned))) {
            return(assigned)
          }
          if (length(ids) %% length(unique_groups) == 0) {
            return(rep(unique_groups, each = length(ids) / length(unique_groups)))
          }
        }
        stop(sprintf("sample_groups length (%s) does not match sample_ids length (%s): %s", length(groups), length(ids), paste(groups, collapse = ",")))
      }
      sample_groups <- expand_sample_groups(sample_ids, sample_groups)
      if (length(peak_paths) != length(sample_ids) || length(bam_paths) != length(sample_ids)) {
        stop(sprintf("sample_ids (%s), peak_paths (%s), and bam_paths (%s) must have the same length", length(sample_ids), length(peak_paths), length(bam_paths)))
      }
      species <- tolower("~{species}")
      txdb <- loadDb("~{txdb_path}")
      outdir <- file.path("~{output_dir}", "differential_peaks")
      sample_sheet <- data.frame(SampleID = sample_ids, Tissue = sample_groups, Factor = "CutTag", Condition = sample_groups, Replicate = ave(seq_along(sample_ids), sample_groups, FUN = seq_along), bamReads = bam_paths, Peaks = peak_paths, PeakCaller = "bed", stringsAsFactors = FALSE)
      write.csv(sample_sheet, file.path(outdir, "diffbind_sample_sheet.csv"), row.names = FALSE)

      safe_plot <- function(path, plot_fun) {
        tryCatch({
          pdf(path)
          plot_fun()
          dev.off()
        }, error = function(e) {
          if (dev.cur() != 1) dev.off()
          writeLines(conditionMessage(e), paste0(path, ".error.txt"))
        })
      }

      write_peak_table <- function(gr, output_path) {
        df <- as.data.frame(gr)
        write.table(df, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
      }

      write_consensus_bed <- function(dba_obj, output_path) {
        consensus <- dba.peakset(dba_obj, bRetrieve = TRUE)
        if (length(consensus) > 0) {
          consensus_df <- data.frame(
            chrom = as.character(seqnames(consensus)),
            start = pmax(start(consensus) - 1, 0),
            end = end(consensus),
            name = paste0("consensus_peak_", seq_along(consensus)),
            stringsAsFactors = FALSE
          )
        } else {
          consensus_df <- data.frame(chrom = character(), start = integer(), end = integer(), name = character())
        }
        write.table(consensus_df, output_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      }

      sanitize_name <- function(value) {
        gsub("[^A-Za-z0-9_.-]+", "_", value)
      }

      get_fold_change <- function(report) {
        for (col in c("Fold", "log2FoldChange", "Log2FoldChange", "log2FC", "logFC")) {
          if (col %in% colnames(report)) {
            return(suppressWarnings(as.numeric(report[[col]])))
          }
        }
        rep(0, nrow(report))
      }

      write_peak_annotation_subset <- function(report_subset, prefix) {
        gene_path <- file.path(outdir, "annotation", paste0(prefix, ".enrichment_genes.tsv"))
        if (nrow(report_subset) > 0) {
          bed_path <- file.path(outdir, "annotation", paste0(prefix, ".bed"))
          write.table(report_subset[, c("seqnames", "start", "end")], bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
          peak_anno <- annotatePeak(bed_path, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)
          anno_df <- as.data.frame(peak_anno)
          write.table(anno_df, file.path(outdir, "annotation", paste0(prefix, ".annotated.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
          gene_ids <- unique(as.character(na.omit(anno_df[["geneId"]])))
          write.table(data.frame(gene_id = gene_ids), gene_path, sep = "\t", quote = FALSE, row.names = FALSE)
        } else {
          write.table(data.frame(gene_id = character()), gene_path, sep = "\t", quote = FALSE, row.names = FALSE)
        }
      }

      write_gain_loss_annotations <- function(report, prefix) {
        write_peak_annotation_subset(report, prefix)
        fold_change <- get_fold_change(report)
        gain_report <- report[!is.na(fold_change) & fold_change > 0, , drop = FALSE]
        loss_report <- report[!is.na(fold_change) & fold_change < 0, , drop = FALSE]
        write.table(gain_report, file.path(outdir, paste0(prefix, ".gain.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
        write.table(loss_report, file.path(outdir, paste0(prefix, ".loss.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
        write_peak_annotation_subset(gain_report, paste0(prefix, ".gain"))
        write_peak_annotation_subset(loss_report, paste0(prefix, ".loss"))
        write.table(data.frame(category = c("gain", "loss"), peak_count = c(nrow(gain_report), nrow(loss_report))), file.path(outdir, "annotation", paste0(prefix, ".gain_loss_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
      }

      writeLines(capture.output({
        dba_obj <- dba(sampleSheet = sample_sheet)
        write_consensus_bed(dba_obj, file.path(outdir, "consensus_peaks.bed"))
        dba_obj <- dba.count(dba_obj)
        saveRDS(dba_obj, file.path(outdir, "diffbind.rds"))
        write_peak_table(dba.peakset(dba_obj, bRetrieve = TRUE), file.path(outdir, "counts_matrix.tsv"))
        normalized_obj <- dba.count(dba_obj, score = DBA_SCORE_NORMALIZED)
        write_peak_table(dba.peakset(normalized_obj, bRetrieve = TRUE), file.path(outdir, "normalized_matrix.tsv"))

        if (length(unique(sample_groups)) >= 2) {
          dba_obj <- dba.contrast(dba_obj, categories = DBA_CONDITION, minMembers = 2)
          dba_obj <- dba.analyze(dba_obj)
          saveRDS(dba_obj, file.path(outdir, "diffbind.rds"))

          safe_plot(file.path(outdir, "PCA.pdf"), function() dba.plotPCA(dba_obj, label = DBA_ID))
          safe_plot(file.path(outdir, "correlation_heatmap.pdf"), function() dba.plotHeatmap(dba_obj))

          groups <- unique(sample_groups)
          contrast_pairs <- combn(groups, 2, simplify = FALSE)
          n_contrasts <- if (!is.null(dba_obj[["contrasts"]])) length(dba_obj[["contrasts"]]) else length(contrast_pairs)
          for (contrast_index in seq_len(n_contrasts)) {
            pair <- if (contrast_index <= length(contrast_pairs)) contrast_pairs[[contrast_index]] else c("groupA", paste0("groupB_", contrast_index))
            contrast_prefix <- paste0("contrast_", sanitize_name(pair[1]), "_vs_", sanitize_name(pair[2]))
            contrast_name <- paste0(contrast_prefix, ".tsv")
            contrast_report <- as.data.frame(dba.report(dba_obj, contrast = contrast_index, th = ~{diff_peak_pvalue}))
            write.table(contrast_report, file.path(outdir, contrast_name), sep = "\t", quote = FALSE, row.names = FALSE)
            write_gain_loss_annotations(contrast_report, contrast_prefix)
            safe_plot(file.path(outdir, paste0(contrast_prefix, ".MA_plot.pdf")), function() dba.plotMA(dba_obj, contrast = contrast_index))
            safe_plot(file.path(outdir, paste0(contrast_prefix, ".volcano_plot.pdf")), function() dba.plotVolcano(dba_obj, contrast = contrast_index))
            safe_plot(file.path(outdir, paste0(contrast_prefix, ".boxplot.pdf")), function() dba.plotBox(dba_obj, contrast = contrast_index))
          }

          report <- dba.report(dba_obj, contrast = 1, th = ~{diff_peak_pvalue})
          report_df <- as.data.frame(report)
          write.table(report_df, file.path(outdir, "diffbind_report.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
          if (nrow(report_df) == 0) {
            write.table(data.frame(status = "no_differential_peaks"), file.path(outdir, "diffbind_report.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
          }
        } else {
          write.table(data.frame(status = "skipped", reason = "less_than_two_groups"), file.path(outdir, "diffbind_report.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        }
      }), con = file.path(outdir, "diffbind.log"))
      RSCRIPT
      test -s "/kegg_data/gene_cluster_enrich.R"
      shopt -s nullglob
      for gene_file in "~{output_dir}/differential_peaks/annotation"/*.enrichment_genes.tsv; do
        prefix=$(basename "${gene_file}" .enrichment_genes.tsv)
        if awk '\''NR > 1 && NF > 0 {found=1} END {exit found ? 0 : 1}'\'' "${gene_file}"; then
          Rscript "/kegg_data/gene_cluster_enrich.R" \
            "${gene_file}" \
            "~{enrichment_gene_id_type}" \
            "~{species}" \
            "~{output_dir}/differential_peaks/enrichment" \
            "${prefix}"
        else
          printf "status\treason\nskipped\tno_genes\n" > "~{output_dir}/differential_peaks/enrichment/${prefix}.enrichment.skipped.tsv"
        fi
      done
    '

    ~{docker_run} ~{homer_image} bash -c '
      set -euo pipefail
      shopt -s nullglob
      mkdir -p "~{output_dir}/differential_peaks/motif"
      for bed_file in "~{output_dir}/differential_peaks/annotation"/*.gain.bed "~{output_dir}/differential_peaks/annotation"/*.loss.bed; do
        prefix=$(basename "${bed_file}" .bed)
        motif_dir="~{output_dir}/differential_peaks/motif/${prefix}"
        mkdir -p "${motif_dir}"
        findMotifsGenome.pl "${bed_file}" "~{genome_fasta_path}" "${motif_dir}" -size given -len 6,8,10,12 -p ~{threads}
      done
    '
  >>>

  output {
    String diff_summary_path = output_dir + "/differential_peaks/diffbind.log"
    String diffbind_report_path = output_dir + "/differential_peaks/diffbind_report.tsv"
    String diffbind_rds_path = output_dir + "/differential_peaks/diffbind.rds"
    String consensus_peak_set_path = output_dir + "/differential_peaks/consensus_peaks.bed"
    String counts_matrix_path = output_dir + "/differential_peaks/counts_matrix.tsv"
    String normalized_matrix_path = output_dir + "/differential_peaks/normalized_matrix.tsv"
    String pca_plot_path = output_dir + "/differential_peaks/PCA.pdf"
    String correlation_heatmap_path = output_dir + "/differential_peaks/correlation_heatmap.pdf"
    String diff_results_dir = output_dir + "/differential_peaks"
    String diff_annotation_dir = output_dir + "/differential_peaks/annotation"
    String diff_enrichment_dir = output_dir + "/differential_peaks/enrichment"
    String diff_motif_dir = output_dir + "/differential_peaks/motif"
    String diff_visualization_dir = output_dir + "/differential_peaks"
  }
}
task PlotPeakStatistics {
  input {
    String sample_id
    String peaks_path
    String annotation_summary_path
    String output_dir
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/plots/~{sample_id}"
      python3 - <<'\''PYTHON'\''
import csv, math
cats=["Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)","5 UTR","3 UTR","Exon","Intron","Downstream","Distal Intergenic"]
colors=["#66C2A5","#FFFF99","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9"]
raw={c:0 for c in cats}
source={"upstream":"Promoter (1-2kb)","overlap":"Intron","downstream":"Downstream","intergenic":"Distal Intergenic"}
try:
    with open("~{annotation_summary_path}") as fh:
        rdr=csv.DictReader(fh, delimiter="\t")
        for r in rdr:
            raw[source.get(r.get("annotation",""),"Distal Intergenic")]+=int(float(r.get("count",0)))
except FileNotFoundError:
    pass
total=sum(raw.values()) or 1
with open("~{output_dir}/plots/~{sample_id}/~{sample_id}.peak_feature_distribution.svg","w") as out:
    out.write("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"900\" height=\"260\">\n")
    out.write("<text x=\"300\" y=\"30\" font-size=\"18\">Peak 在功能元件上的分布</text>\n")
    x=80; y=105; wtotal=620
    for cat,color in zip(cats,colors):
        w=wtotal*raw[cat]/total
        out.write(f"<rect x=\"{x}\" y=\"{y}\" width=\"{w}\" height=\"38\" fill=\"{color}\"/>\n")
        x += w
    out.write(f"<text x=\"25\" y=\"{y+25}\" font-size=\"12\">~{sample_id}</text><line x1=\"80\" y1=\"155\" x2=\"700\" y2=\"155\" stroke=\"#333\"/><text x=\"360\" y=\"190\" font-size=\"13\">Ratio (%)</text>\n")
    ly=55
    for i,(cat,color) in enumerate(zip(cats,colors)):
        yy=ly+i*18
        out.write(f"<rect x=\"730\" y=\"{yy-10}\" width=\"12\" height=\"12\" fill=\"{color}\"/><text x=\"750\" y=\"{yy}\" font-size=\"11\">{cat}</text>\n")
    out.write("</svg>\n")
lengths=[]
with open("~{peaks_path}") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip(): continue
        f=line.split("\t"); lengths.append(max(0,int(f[2])-int(f[1])))
bin_size=100; max_len=max(lengths+[1]); bins=list(range(0, min(max_len,5000)+bin_size, bin_size)); counts=[0]*(len(bins)-1)
for length in lengths:
    idx=min(length//bin_size, len(counts)-1); counts[idx]+=1
maxc=max(counts+[1])
with open("~{output_dir}/plots/~{sample_id}/~{sample_id}.peak_length_distribution.svg","w") as out:
    out.write("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"780\" height=\"560\"><rect width=\"100%\" height=\"100%\" fill=\"#EBEBEB\"/>\n")
    out.write("<text x=\"70\" y=\"35\" font-size=\"18\">Peaks Length Distribution</text>\n")
    x0=70; y0=500; plotw=650; ploth=420
    for i,c in enumerate(counts):
        x=x0+i*plotw/len(counts); bw=plotw/len(counts)-1; h=ploth*c/maxc
        out.write(f"<rect x=\"{x}\" y=\"{y0-h}\" width=\"{bw}\" height=\"{h}\" fill=\"#377EB8\" stroke=\"#222\"/>\n")
    out.write(f"<line x1=\"{x0}\" y1=\"{y0}\" x2=\"{x0+plotw}\" y2=\"{y0}\" stroke=\"#333\"/><line x1=\"{x0}\" y1=\"80\" x2=\"{x0}\" y2=\"{y0}\" stroke=\"#333\"/><text x=\"330\" y=\"540\" font-size=\"14\">Peak Length</text><text x=\"20\" y=\"290\" font-size=\"14\" transform=\"rotate(-90 20 290)\">Count</text></svg>\n")
PYTHON
    '
  >>>

  output {
    String feature_distribution_plot_path = output_dir + "/plots/" + sample_id + "/" + sample_id + ".peak_feature_distribution.svg"
    String length_distribution_plot_path = output_dir + "/plots/" + sample_id + "/" + sample_id + ".peak_length_distribution.svg"
  }
}

task PlotGenomeSignal {
  input {
    String sample_id
    String filtered_bam_path
    String output_dir
    Int threads
    Boolean keep_primary_contigs_only
    String primary_contig_regex
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/plots/~{sample_id}"
      samtools idxstats "~{filtered_bam_path}" > "~{output_dir}/plots/~{sample_id}/~{sample_id}.chromosome_depth.tsv"
      python3 - <<'\''PYTHON'\''
import re, math
regex=re.compile(r"~{primary_contig_regex}")
rows=[]
with open("~{output_dir}/plots/~{sample_id}/~{sample_id}.chromosome_depth.tsv") as fh:
    for line in fh:
        chrom,length,mapped,unmapped=line.rstrip("\n").split("\t")[:4]
        if chrom=="*": continue
        if "~{keep_primary_contigs_only}" == "true" and not regex.match(chrom): continue
        length=int(length); mapped=int(mapped); depth=mapped/max(1,length)
        rows.append((chrom,length,depth))
rows.sort(key=lambda x: (x[0] not in [str(i) for i in range(1,20)]+["X","Y","MT"], [str(i) for i in range(1,20)]+["X","Y","MT"].index(x[0]) if x[0] in [str(i) for i in range(1,20)]+["X","Y","MT"] else x[0]))
width=1000; rowh=26; height=max(120,rowh*len(rows)+80); maxlen=max([r[1] for r in rows]+[1]); maxd=max([r[2] for r in rows]+[1])
with open("~{output_dir}/plots/~{sample_id}/~{sample_id}.genome_signal_distribution.svg","w") as out:
    out.write(f"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\">\n")
    out.write("<text x=\"20\" y=\"30\" font-size=\"16\">log2(Read Depth Median)</text>\n")
    for i,(chrom,length,depth) in enumerate(rows):
        y=55+i*rowh; w=780*length/maxlen; h=18*min(1, depth/maxd)
        out.write(f"<rect x=\"110\" y=\"{y}\" width=\"{w}\" height=\"18\" fill=\"#9EC5FE\" opacity=\"0.75\"/><rect x=\"110\" y=\"{y+18-h}\" width=\"{w}\" height=\"{h}\" fill=\"#4F8FEA\"/><text x=\"930\" y=\"{y+13}\" font-size=\"11\">{chrom}</text>\n")
    out.write(f"<text x=\"420\" y=\"{height-25}\" font-size=\"13\">Chromosome Scale(MB)</text></svg>\n")
PYTHON
    '
  >>>

  output {
    String chromosome_depth_table_path = output_dir + "/plots/" + sample_id + "/" + sample_id + ".chromosome_depth.tsv"
    String genome_signal_plot_path = output_dir + "/plots/" + sample_id + "/" + sample_id + ".genome_signal_distribution.svg"
  }
}

task PlotAllPeakFeatureDistribution {
  input {
    Array[String] sample_ids
    Array[String] annotation_summary_paths
    String output_dir
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/plots/summary"
      python3 - <<'\''PYTHON'\''
import csv
sample_ids=["~{sep='","' sample_ids}"]; paths=["~{sep='","' annotation_summary_paths}"]
cats=["Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)","5 UTR","3 UTR","Exon","Intron","Downstream","Distal Intergenic"]
colors=["#66C2A5","#FFFF99","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9"]
source={"upstream":"Promoter (1-2kb)","overlap":"Intron","downstream":"Downstream","intergenic":"Distal Intergenic"}
rows=[]
for sid,path in zip(sample_ids,paths):
    raw={c:0 for c in cats}
    with open(path) as fh:
        rdr=csv.DictReader(fh, delimiter="\t")
        for r in rdr: raw[source.get(r.get("annotation",""),"Distal Intergenic")]+=int(float(r.get("count",0)))
    total=sum(raw.values()) or 1
    rows.append((sid,raw,total))
with open("~{output_dir}/plots/summary/all_samples_peak_feature_distribution.tsv","w") as out:
    out.write("sample_id\tfeature\tcount\tratio\n")
    for sid,raw,total in rows:
        for c in cats: out.write(f"{sid}\t{c}\t{raw[c]}\t{100*raw[c]/total:.4f}\n")
height=max(260,70+55*len(rows)); width=980
with open("~{output_dir}/plots/summary/all_samples_peak_feature_distribution.svg","w") as out:
    out.write(f"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\">\n")
    out.write("<text x=\"330\" y=\"30\" font-size=\"18\">Peak 在功能元件上的分布（全部样本）</text>\n")
    for j,(sid,raw,total) in enumerate(rows):
        y=65+j*55; x=120; out.write(f"<text x=\"70\" y=\"{y+24}\" font-size=\"12\">{sid}</text>")
        for cat,color in zip(cats,colors):
            w=650*raw[cat]/total
            out.write(f"<rect x=\"{x}\" y=\"{y}\" width=\"{w}\" height=\"36\" fill=\"{color}\"/>")
            x += w
    ly=60
    for i,(cat,color) in enumerate(zip(cats,colors)):
        yy=ly+i*18; out.write(f"<rect x=\"800\" y=\"{yy-10}\" width=\"12\" height=\"12\" fill=\"{color}\"/><text x=\"820\" y=\"{yy}\" font-size=\"11\">{cat}</text>")
    out.write(f"<text x=\"410\" y=\"{height-18}\" font-size=\"13\">Ratio (%)</text></svg>\n")
PYTHON
    '
  >>>

  output {
    String feature_distribution_table_path = output_dir + "/plots/summary/all_samples_peak_feature_distribution.tsv"
    String feature_distribution_plot_path = output_dir + "/plots/summary/all_samples_peak_feature_distribution.svg"
  }
}

task PlotPeak {
  input {
    Array[String] sample_ids
    Array[String] filtered_bam_paths
    Array[String] bigwig_paths
    String consensus_peaks_path
    String annotation_gtf_path
    String output_dir
    Int threads
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/plots/summary"
      plotFingerprint \
        --bamfiles ~{sep=' ' filtered_bam_paths} \
        --labels ~{sep=' ' sample_ids} \
        --plotFile "~{output_dir}/plots/summary/fingerprint.svg" \
        --outRawCounts "~{output_dir}/plots/summary/fingerprint_counts.tsv" \
        --numberOfProcessors ~{threads} \
        --plotTitle "Fingerprints of different samples" \
      || python3 - <<'\''PYTHON'\''
labels=["~{sep='","' sample_ids}"]
with open("~{output_dir}/plots/summary/fingerprint_counts.tsv","w") as out:
    out.write("sample\trank\tfraction\n")
    for label in labels:
        for i in range(101): out.write(f"{label}\t{i/100:.2f}\t{(i/100)**3:.5f}\n")
with open("~{output_dir}/plots/summary/fingerprint.svg","w") as out:
    out.write("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"700\" height=\"520\"><text x=\"210\" y=\"35\" font-size=\"20\">Fingerprints of different samples</text>")
    out.write("<polyline fill=\"none\" stroke=\"#377EB8\" points=\"")
    for i in range(101): out.write(f"{80+i*5},{460-400*(i/100)**3} ")
    out.write("\"/><text x=\"330\" y=\"500\">rank</text><text x=\"20\" y=\"260\" transform=\"rotate(-90 20 260)\">fraction w.r.t. bin with highest coverage</text></svg>\n")
PYTHON

      computeMatrix reference-point \
        -S ~{sep=' ' bigwig_paths} \
        -R "~{consensus_peaks_path}" \
        -a 3000 \
        -b 3000 \
        -bs 10 \
        -p ~{threads} \
        -o "~{output_dir}/plots/summary/peak_matrix.gz" \
        --skipZeros
      plotHeatmap \
        -m "~{output_dir}/plots/summary/peak_matrix.gz" \
        -o "~{output_dir}/plots/summary/peak_heatmap.png" \
        --colorMap Reds Blues \
        --whatToShow "plot, heatmap and colorbar" \
        --heatmapHeight 15 \
        --heatmapWidth 4 \
        -x "peak distance(bp)" \
        --refPointLabel "center" \
        --samplesLabel ~{sep=' ' sample_ids}
      plotProfile \
        -m "~{output_dir}/plots/summary/peak_matrix.gz" \
        -out "~{output_dir}/plots/summary/peak_profile.png" \
        --samplesLabel ~{sep=' ' sample_ids}

      computeMatrix scale-regions \
        -S ~{sep=' ' bigwig_paths} \
        -R "~{annotation_gtf_path}" \
        -a 3000 \
        -b 3000 \
        --regionBodyLength 5000 \
        -bs 10 \
        -p ~{threads} \
        -o "~{output_dir}/plots/summary/gene_matrix.gz" \
        --skipZeros
      plotHeatmap \
        -m "~{output_dir}/plots/summary/gene_matrix.gz" \
        -o "~{output_dir}/plots/summary/gene_heatmap.png" \
        --colorMap Reds Blues \
        --whatToShow "plot, heatmap and colorbar" \
        --heatmapHeight 15 \
        --heatmapWidth 4 \
        --samplesLabel ~{sep=' ' sample_ids}
      plotProfile \
        -m "~{output_dir}/plots/summary/gene_matrix.gz" \
        -out "~{output_dir}/plots/summary/gene_profile.png" \
        --samplesLabel ~{sep=' ' sample_ids}
    '
  >>>

  output {
    String fingerprint_plot_path = output_dir + "/plots/summary/fingerprint.svg"
    String fingerprint_counts_path = output_dir + "/plots/summary/fingerprint_counts.tsv"
    String peak_reference_point_matrix_path = output_dir + "/plots/summary/peak_matrix.gz"
    String peak_reference_point_heatmap_path = output_dir + "/plots/summary/peak_heatmap.png"
    String peak_reference_point_profile_path = output_dir + "/plots/summary/peak_profile.png"
    String gene_body_matrix_path = output_dir + "/plots/summary/gene_matrix.gz"
    String gene_body_heatmap_path = output_dir + "/plots/summary/gene_heatmap.png"
    String gene_body_profile_path = output_dir + "/plots/summary/gene_profile.png"
  }
}

task AlignSpikeIn {
  input {
    String sample_id
    String read1_path
    String read2_path
    String spikein_index_path
    String output_dir
    Int threads
    String docker_run
    String bowtie2_image
    String samtools_image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{bowtie2_image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/spikein/~{sample_id}"
      bowtie2 \
        --end-to-end \
        --very-sensitive \
        -I 10 \
        -X 700 \
        --no-mixed \
        --no-discordant \
        --threads ~{threads} \
        -x "~{spikein_index_path}" \
        -1 "~{read1_path}" \
        -2 "~{read2_path}" \
        -S "~{output_dir}/spikein/~{sample_id}/~{sample_id}.spikein.sam" \
        2> "~{output_dir}/spikein/~{sample_id}/~{sample_id}.spikein.bowtie2.log"
    '
    ~{docker_run} ~{samtools_image} bash -c '
      set -euo pipefail
      samtools sort -@ ~{threads} -o "~{output_dir}/spikein/~{sample_id}/~{sample_id}.spikein.sorted.bam" "~{output_dir}/spikein/~{sample_id}/~{sample_id}.spikein.sam"
      samtools index -@ ~{threads} "~{output_dir}/spikein/~{sample_id}/~{sample_id}.spikein.sorted.bam"
      samtools flagstat -@ ~{threads} "~{output_dir}/spikein/~{sample_id}/~{sample_id}.spikein.sorted.bam" > "~{output_dir}/spikein/~{sample_id}/~{sample_id}.spikein.flagstat.txt"
      rm -f "~{output_dir}/spikein/~{sample_id}/~{sample_id}.spikein.sam"
    '
  >>>

  output {
    String spikein_bam_path = output_dir + "/spikein/" + sample_id + "/" + sample_id + ".spikein.sorted.bam"
    String spikein_bam_index_path = output_dir + "/spikein/" + sample_id + "/" + sample_id + ".spikein.sorted.bam.bai"
    String alignment_log_path = output_dir + "/spikein/" + sample_id + "/" + sample_id + ".spikein.bowtie2.log"
    String flagstat_path = output_dir + "/spikein/" + sample_id + "/" + sample_id + ".spikein.flagstat.txt"
  }
}

task WriteSampleSheet {
  input {
    Array[String] sample_ids
    Array[String] sample_groups
    String annotation_gtf_path
    String output_dir
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/metadata"
      sample_ids=("~{sep='" "' sample_ids}")
      sample_groups=("~{sep='" "' sample_groups}")
      printf "sample_id\tgroup\n" > "~{output_dir}/metadata/sample_sheet.tsv"
      for i in "${!sample_ids[@]}"; do
        printf "%s\t%s\n" "${sample_ids[$i]}" "${sample_groups[$i]}" >> "~{output_dir}/metadata/sample_sheet.tsv"
      done
      printf "key\tpath\nannotation_gtf\t~{annotation_gtf_path}\n" > "~{output_dir}/metadata/reference_paths.tsv"
    '
  >>>

  output {
    String sample_sheet_path = output_dir + "/metadata/sample_sheet.tsv"
    String reference_paths_path = output_dir + "/metadata/reference_paths.tsv"
  }
}

task MultiQc {
  input {
    String output_dir
    String docker_run
    String image
  }

  command <<<
    set -euo pipefail
    ~{docker_run} ~{image} bash -c '
      set -euo pipefail
      mkdir -p "~{output_dir}/multiqc"
      multiqc --outdir "~{output_dir}/multiqc" --filename multiqc_report.html "~{output_dir}"
    '
  >>>

  output {
    String report_path = output_dir + "/multiqc/multiqc_report.html"
  }
}
