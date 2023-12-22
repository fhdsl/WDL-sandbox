version 1.0
## Testing if alignement for 1 sample will work
## Samples:
## NCIH2172_LUNG : Cell-line with EGFR L858R mutation
##
## Input requirements:
## - fastq files for reference and alternative]
##
## Output Files:
## - An aligned bam for 1 sample
## 
## Workflow developed by Sitapriya Moorthi @ Fred Hutch LMD: 11/22/23 for use by DaSL @ Fred Hutch.

workflow minidata_test_alignment {
  input {
    # Batch and cohort information
    File sampleFastq
    String base_file_name = "SRR8618962_combined"
    # Reference Genome information
    String ref_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative". Leave blank in JSON for legacy
    # references such as b37 and hg19.
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    }

    # Docker containers this workflow has been designed for
    String bwadocker = "fredhutch/bwa:0.7.17"
    String GATKdocker = "broadinstitute/gatk:4.1.4.0"
    
    #Array[Object] batchInfo = read_objects(batchFile)

  #  Map reads to reference
  call BwaMem as sampleBwaMem {
    input:
      input_fastq = sampleFastq,
      base_file_name = base_file_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_alt = ref_alt,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      taskDocker = bwadocker
  }
  # Outputs that will be retained when execution is complete
  output {
    File analysisReadyBam = "~{base_file_name}" +".aligned.bam"
    File analysisReadyBamSorted = "~{base_file_name}" +"_sorted_query_aligned.bam"
  }
   
   # Mark duplicates
  call MarkDuplicatesSpark as MarkDuplicates {
    input:
      input_bam = base_file_name + "_sorted_query_aligned.bam",
      output_bam_basename = base_file_name + ".aligned.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      taskDocker = GATKdocker
  }
# End workflow
}
#### TASK DEFINITIONS

# align to genome
## Currently uses -M but GATK uses -Y and no -M
task BwaMem {
  input {
    File input_fastq
    String base_file_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    String taskDocker
  }
  command <<<
    set -eo pipefail

    bwa mem \
      -p -v 3 -t 16 -M \
      ~{ref_fasta} ~{input_fastq} > ~{base_file_name}.sam 
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
    docker: taskDocker
    walltime: "2:00:00"
  }
}

task MarkDuplicatesSpark {
  input {
    File input_bam
    String output_bam_basename
    String metrics_filename
    String taskDocker
  }
  # Later use: --verbosity WARNING
 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Dsamjdk.compression_level=5 -Xms32g" \
      MarkDuplicatesSpark \
      --input ~{input_bam} \
      --output ~{output_bam_basename}.bam \
      --metrics-file ~{metrics_filename} \
      --optical-duplicate-pixel-distance 2500 
  }
  runtime {
    docker: taskDocker
    memory: "48 GB"
    cpu: 4
    walltime: "6:00:00"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bai = "~{output_bam_basename}.bam.bai"
    File duplicate_metrics = "~{metrics_filename}"
  }
}



