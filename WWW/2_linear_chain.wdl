version 1.0
## Testing if alignement for 1 sample will work
## Samples:
## NCIH2172_LUNG : Cell-line with EGFR L858R mutation
##
## Input requirements:
## - fastq files for reference and alternative
##
## Output Files:
## - An aligned, marked duplicated bam for 1 sample
## 
## Workflow developed by Sitapriya Moorthi @ Fred Hutch LMD: 11/22/23 for use by DaSL @ Fred Hutch.

workflow minidata_test_alignment {
  input {
    # Sample info
    File sampleFastq
    # Reference Genome information
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
  }

  #  Map reads to reference
  call BwaMem {
    input:
      input_fastq = sampleFastq,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa
  }
   
  # Mark duplicates
  call MarkDuplicatesSpark {
    input:
      input_bam = BwaMem.analysisReadySorted
  }

  # Outputs that will be retained when execution is complete
  output {
    File alignedBamSorted = BwaMem.analysisReadySorted
    File markDuplicates = MarkDuplicatesSpark.output_bam
  }

  parameter_meta {
    sampleFastq: "Filepath to sample .fastq"
    ref_fasta: "Filepath to reference genome"
    ref_fasta_index: "Filepath to reference genome index"
    ref_dict: "Filepath to reference genome dictionary"
    ref_amb: "Filepath to reference genome info"
    ref_ann: "Filepath to reference genome info"
    ref_bwt: "Filepath to reference genome info"
    ref_pac: "Filepath to reference genome info"
    ref_sa: "Filepath to reference genome info"
  }
# End workflow
}

#### TASK DEFINITIONS
# align to genome
## Currently uses -M but GATK uses -Y and no -M
task BwaMem {
  input {
    File input_fastq
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
  }
  
  String base_file_name = basename(input_fastq, ".fastq")
  String ref_fasta_local = basename(ref_fasta)

  # CL added some fake read groups for MarkDuplicates to run. 
  command <<<
    set -eo pipefail

    mv ~{ref_fasta} .
    mv ~{ref_fasta_index} .
    mv ~{ref_dict} .
    mv ~{ref_amb} .
    mv ~{ref_ann} .
    mv ~{ref_bwt} .
    mv ~{ref_pac} .
    mv ~{ref_sa} .

    bwa mem \
      -p -v 3 -t 16 -M -R '@RG\tID:foo\tSM:foo2' \
      ~{ref_fasta_local} ~{input_fastq} > ~{base_file_name}.sam 
    samtools view -1bS -@ 15 -o ~{base_file_name}.aligned.bam ~{base_file_name}.sam
    samtools sort -n -@ 15 -o ~{base_file_name}.sorted_query_aligned.bam ~{base_file_name}.aligned.bam

  >>>
  output {
    File analysisReadyBam = "~{base_file_name}.aligned.bam"
    File analysisReadySorted = "~{base_file_name}.sorted_query_aligned.bam"
  }
  runtime {
    memory: "48 GB"
    cpu: 16
    docker: "fredhutch/bwa:0.7.17"
    disks: "local-disk 100 SSD"
  }
}

task MarkDuplicatesSpark {
  input {
    File input_bam
  }
  
  String base_file_name = basename(input_bam, ".sorted_query_aligned.bam")
  String output_bam_string = "~{base_file_name}.duplicates_marked.bam"
  String metrics_file = "~{base_file_name}.duplicate_metrics"

  # Later use: --verbosity WARNING
  # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
  # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
  # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Dsamjdk.compression_level=5 -Xms32g" \
      MarkDuplicatesSpark \
      --input ~{input_bam} \
      --output ~{output_bam_string} \
      --metrics-file ~{metrics_file} \
      --optical-duplicate-pixel-distance 2500 
  }
  runtime {
    docker: "broadinstitute/gatk:4.1.4.0"
    memory: "48 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
  }
  output {
    File output_bam = "~{output_bam_string}"
    File output_bai = "~{output_bam_string}.bai"
    File duplicate_metrics = "~{metrics_file}"
  }
}



