version 1.0

## Paired-end RNA-seq sample definition.
## All paths must be absolute host paths under /data because each image command
## includes a Docker runner that mounts only /data inside containers.
struct RnaSeqSample {
  String sample_id
  String read1_path
  String read2_path
}

workflow RnaSeq {
  input {
    Array[RnaSeqSample] samples
    String genome_fasta_path
    String annotation_gtf_path
    String output_dir = "/data/rnaseq_results"
    Int threads = 8
    Int sjdb_overhang = 149

    String fastqc_image = "docker run --rm -v /data2:/data registry.cn-guangzhou.aliyuncs.com/origen/fastqc"
    String trim_galore_image = "docker run --rm -v /data:/data registry.cn-guangzhou.aliyuncs.com/origen/trim-galore"
    String star_image = "docker run --rm -v /data:/data registry.cn-guangzhou.aliyuncs.com/origen/star"
    String subread_image = "docker run --rm -v /data:/data registry.cn-guangzhou.aliyuncs.com/origen/subread"
    String multiqc_image = "docker run --rm -v /data:/data registry.cn-guangzhou.aliyuncs.com/origen/multiqc"
  }

  call BuildStarIndex {
    input:
      genome_fasta_path = genome_fasta_path,
      annotation_gtf_path = annotation_gtf_path,
      output_dir = output_dir,
      threads = threads,
      sjdb_overhang = sjdb_overhang,
      image = star_image
  }

  scatter (sample in samples) {
    call FastQc as RawFastQc {
      input:
        sample_id = sample.sample_id,
        read1_path = sample.read1_path,
        read2_path = sample.read2_path,
        output_dir = output_dir,
        threads = threads,
        image = fastqc_image
    }

    call TrimReads {
      input:
        sample_id = sample.sample_id,
        read1_path = sample.read1_path,
        read2_path = sample.read2_path,
        output_dir = output_dir,
        threads = threads,
        image = trim_galore_image
    }

    call FastQc as TrimmedFastQc {
      input:
        sample_id = sample.sample_id + ".trimmed",
        read1_path = TrimReads.trimmed_read1_path,
        read2_path = TrimReads.trimmed_read2_path,
        output_dir = output_dir,
        threads = threads,
        image = fastqc_image
    }

    call StarAlign {
      input:
        sample_id = sample.sample_id,
        read1_path = TrimReads.trimmed_read1_path,
        read2_path = TrimReads.trimmed_read2_path,
        genome_dir = BuildStarIndex.genome_dir,
        output_dir = output_dir,
        threads = threads,
        image = star_image
    }
  }

  call FeatureCounts {
    input:
      bam_paths = StarAlign.sorted_bam_path,
      annotation_gtf_path = annotation_gtf_path,
      output_dir = output_dir,
      threads = threads,
      image = subread_image
  }

  call MultiQc {
    input:
      output_dir = output_dir,
      count_matrix_path = FeatureCounts.count_matrix_path,
      image = multiqc_image
  }

  output {
    String star_index_dir = BuildStarIndex.genome_dir
    Array[String] raw_fastqc_dirs = RawFastQc.fastqc_dir_path
    Array[String] trimmed_fastqc_dirs = TrimmedFastQc.fastqc_dir_path
    Array[String] trimmed_read1 = TrimReads.trimmed_read1_path
    Array[String] trimmed_read2 = TrimReads.trimmed_read2_path
    Array[String] sorted_bams = StarAlign.sorted_bam_path
    Array[String] alignment_logs = StarAlign.final_log_path
    String count_matrix = FeatureCounts.count_matrix_path
    String count_summary = FeatureCounts.summary_path
    String multiqc_report = MultiQc.report_path
  }
}

task BuildStarIndex {
  input {
    String genome_fasta_path
    String annotation_gtf_path
    String output_dir
    Int threads
    Int sjdb_overhang
    String image
  }

  command <<<
    set -euo pipefail
    ~{image} bash -c "mkdir -p '~{output_dir}/star_index'"
    ~{image} bash -c "STAR --runThreadN ~{threads} --runMode genomeGenerate --genomeDir '~{output_dir}/star_index' --genomeFastaFiles '~{genome_fasta_path}' --sjdbGTFfile '~{annotation_gtf_path}' --sjdbOverhang ~{sjdb_overhang}"
  >>>

  output {
    String genome_dir = output_dir + "/star_index"
  }
}

task FastQc {
  input {
    String sample_id
    String read1_path
    String read2_path
    String output_dir
    Int threads
    String image
  }

  command <<<
    set -euo pipefail
    ~{image} bash -c "mkdir -p '~{output_dir}/fastqc/~{sample_id}'"
    ~{image} bash -c "fastqc --threads ~{threads} --outdir '~{output_dir}/fastqc/~{sample_id}' '~{read1_path}' '~{read2_path}'"
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
    String image
  }

  command <<<
    set -euo pipefail
    ~{image} bash -c "mkdir -p '~{output_dir}/trimmed/~{sample_id}'"
    ~{image} bash -c "trim_galore --paired --cores ~{threads} --gzip --basename '~{sample_id}' --output_dir '~{output_dir}/trimmed/~{sample_id}' '~{read1_path}' '~{read2_path}'"
    ~{image} bash -c "test -s '~{output_dir}/trimmed/~{sample_id}/~{sample_id}_val_1.fq.gz'"
    ~{image} bash -c "test -s '~{output_dir}/trimmed/~{sample_id}/~{sample_id}_val_2.fq.gz'"
  >>>

  output {
    String trimmed_read1_path = output_dir + "/trimmed/" + sample_id + "/" + sample_id + "_val_1.fq.gz"
    String trimmed_read2_path = output_dir + "/trimmed/" + sample_id + "/" + sample_id + "_val_2.fq.gz"
    String trim_dir_path = output_dir + "/trimmed/" + sample_id
  }
}

task StarAlign {
  input {
    String sample_id
    String read1_path
    String read2_path
    String genome_dir
    String output_dir
    Int threads
    String image
  }

  command <<<
    set -euo pipefail
    ~{image} bash -c "mkdir -p '~{output_dir}/star/~{sample_id}'"
    ~{image} bash -c "STAR --runThreadN ~{threads} --genomeDir '~{genome_dir}' --readFilesIn '~{read1_path}' '~{read2_path}' --readFilesCommand zcat --outFileNamePrefix '~{output_dir}/star/~{sample_id}/' --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts"
    ~{image} bash -c "mv '~{output_dir}/star/~{sample_id}/Aligned.sortedByCoord.out.bam' '~{output_dir}/star/~{sample_id}/~{sample_id}.sorted.bam'"
  >>>

  output {
    String sorted_bam_path = output_dir + "/star/" + sample_id + "/" + sample_id + ".sorted.bam"
    String final_log_path = output_dir + "/star/" + sample_id + "/Log.final.out"
    String gene_counts_path = output_dir + "/star/" + sample_id + "/ReadsPerGene.out.tab"
  }
}

task FeatureCounts {
  input {
    Array[String] bam_paths
    String annotation_gtf_path
    String output_dir
    Int threads
    String image
  }

  command <<<
    set -euo pipefail
    ~{image} bash -c "mkdir -p '~{output_dir}/counts'"
    ~{image} bash -c "featureCounts -T ~{threads} -p -B -C -a '~{annotation_gtf_path}' -o '~{output_dir}/counts/gene_counts.tsv' ~{sep=" " bam_paths}"
  >>>

  output {
    String count_matrix_path = output_dir + "/counts/gene_counts.tsv"
    String summary_path = output_dir + "/counts/gene_counts.tsv.summary"
  }
}

task MultiQc {
  input {
    String output_dir
    String count_matrix_path
    String image
  }

  command <<<
    set -euo pipefail
    ~{image} bash -c "test -s '~{count_matrix_path}'"
    ~{image} bash -c "mkdir -p '~{output_dir}/multiqc'"
    ~{image} bash -c "multiqc --outdir '~{output_dir}/multiqc' --filename multiqc_report.html '~{output_dir}'"
  >>>

  output {
    String report_path = output_dir + "/multiqc/multiqc_report.html"
  }
}
