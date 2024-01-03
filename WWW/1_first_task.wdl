version 1.0
## Testing if alignement for 1 sample will work
## Samples:
## NCIH2172_LUNG : Cell-line with EGFR L858R mutation
##
## Input requirements:
## - fastq files for reference and alternative
##
## Output Files:
## - An aligned bam for 1 sample
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
  }

  #  Map reads to reference
  call BwaMem {
    input:
      input_fastq = sampleFastq,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict
  }
   
  # Outputs that will be retained when execution is complete
  output {
    File alignedBamSorted = BwaMem.analysisReadySorted
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
  }
  
  String base_file_name = basename(input_fastq, ".fastq")

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
    File analysisReadyBam = "~{base_file_name}.aligned.bam"
    File analysisReadySorted = "~{base_file_name}.sorted_query_aligned.bam"
  }
  runtime {
    memory: "48 GB"
    cpu: 16
    docker: "fredhutch/bwa:0.7.17"
  }
}
