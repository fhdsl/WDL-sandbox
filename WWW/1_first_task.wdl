version 1.0
## Testing if alignement for 1 sample will work
## Samples:
## NCIH2172_LUNG : Cell-line with EGFR L858R mutation
##
## Input requirements:
## - fastq files for reference and alternative]
##
## Output Files:
## - An aligned, markduplicates bam for 1 sample
## 
## Workflow developed by Sitapriya Moorthi & Chris Lo @ Fred Hutch LMD.

workflow minidata_test_alignment {
  input {
    # Sample info
    File input_fastq
    String base_file_name
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
  call BwaMem {
    input:
      input_fastq = input_fastq,
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
    File alignedBamSorted = BwaMem.sorted_bam
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
