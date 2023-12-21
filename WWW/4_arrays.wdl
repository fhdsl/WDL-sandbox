version 1.0
## Testing if alignement for a few samples will work
## Samples:
## NCIH2172_LUNG : Cell-line with EGFR L858R mutation
##
## Input requirements:
## - fastq files for reference and alternative]
##
## Output Files:
## - An aligned bam for many samples
## 
## Workflow developed by Sitapriya Moorthi & Chris Lo @ Fred Hutch LMD.

# Sample info
struct sampleInputs {
  File fastq
  String id
}

# Reference Genome information
struct referenceGenome {
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


workflow minidata_test_alignment {
  input {
    Array[sampleInputs] allSamples
    referenceGenome refGenome
  }

  # Docker containers this workflow has been designed for
  String bwadocker = "fredhutch/bwa:0.7.17"
  String GATKdocker = "broadinstitute/gatk:4.1.4.0"
  
  scatter (sample in allSamples){
    #  Map reads to reference
    call BwaMem {
      input:
        input_fastq = sample.fastq,
        base_file_name = sample.id,
        ref_fasta = refGenome.ref_fasta,
        ref_fasta_index = refGenome.ref_fasta_index,
        ref_dict = refGenome.ref_dict,
        ref_alt = refGenome.ref_alt,
        ref_amb = refGenome.ref_amb,
        ref_ann = refGenome.ref_ann,
        ref_bwt = refGenome.ref_bwt,
        ref_pac = refGenome.ref_pac,
        ref_sa = refGenome.ref_sa,
        taskDocker = bwadocker
    }
    
    # Mark duplicates
    call MarkDuplicatesSpark {
      input:
        input_bam = BwaMem.sorted_bam,
        output_bam_basename = "~{sample.id}.duplicates_marked",
        taskDocker = GATKdocker
    }
  #End scatter
  } 

  # Outputs that will be retained when execution is complete
  output {
    Array[File] alignedBamSorted = BwaMem.sorted_bam
    Array[File] markDuplicates = MarkDuplicatesSpark.markDuplicates_bam
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
  # CL added some fake read groups for MarkDuplicates to run. 
  command <<<
    set -eo pipefail

    bwa mem \
      -p -v 3 -t 16 -M -R '@RG\tID:foo\tSM:foo2' \
      ~{ref_fasta} ~{input_fastq} > ~{base_file_name}.sam 
    samtools view -1bS -@ 15 -o ~{base_file_name}.aligned.bam ~{base_file_name}.sam
    samtools sort -n -@ 15 -o ~{base_file_name}.sorted_query_aligned.bam ~{base_file_name}.aligned.bam

  >>>
  output {
    File bam = "~{base_file_name}.aligned.bam"
    File sorted_bam = "~{base_file_name}.sorted_query_aligned.bam"
    
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
      --metrics-file ~{output_bam_basename}.duplicate_metrics \
      --optical-duplicate-pixel-distance 2500 
  }
  runtime {
    docker: taskDocker
    memory: "48 GB"
    cpu: 4
    walltime: "6:00:00"
  }
  output {
    File markDuplicates_bam = "~{output_bam_basename}.bam"
    File output_bai = "~{output_bam_basename}.bam.bai"
    File duplicate_metrics = "~{output_bam_basename}.duplicate_metrics"
  }
}



