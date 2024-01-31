version 1.0

task BwaMem {
  input {
    File input_fastq
    File ref_fasta
    File ref_fasta_index
    #File ref_dict
    #File ref_amb
    #File ref_ann
    #File ref_bwt
    #File ref_pac
    #File ref_sa
  }
  
  String base_file_name = basename(input_fastq, ".fastq")
  String ref_fasta_local = basename(ref_fasta)
  String read_group_id = "ID:" + base_file_name
  String sample_name = "SM:" + base_file_name
  String platform = "illumina"
  String platform_info = "PL:" + platform   # Create the platform information


  command <<<
    set -eo pipefail
    
    mv ~{ref_fasta} .
    
    bwa index -p goober ~{ref_fasta_local}
    ls
    
    bwa mem \
      -p -P goober -v 3 -t 3 -M -R '@RG\t~{read_group_id}\t~{sample_name}\t~{platform_info}' \
      ~{ref_fasta_local} ~{input_fastq} > ~{base_file_name}.sam 
    samtools view -1bS -o ~{base_file_name}.aligned.bam ~{base_file_name}.sam
    samtools sort -n -o ~{base_file_name}.sorted_query_aligned.bam ~{base_file_name}.aligned.bam
  >>>

  output {
    File analysisReadySorted = "~{base_file_name}.sorted_query_aligned.bam"
  }
  
  runtime {
    memory: "48 GB"
    cpu: 16
    docker: "fredhutch/bwa:0.7.17"
  }
}

workflow test_bwa_mem {
  input {
    # Sample info
    File sampleFastq
    # Reference Genome information
    File ref_fasta
    File ref_fasta_index
    #File ref_dict
    #File ref_amb
    #File ref_ann
    #File ref_bwt
    #File ref_pac
    #File ref_sa
  }

  #  Map reads to reference
  call BwaMem {
    input:
      input_fastq = sampleFastq,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      #ref_dict = ref_dict,
      #ref_amb = ref_amb,
      #ref_ann = ref_ann,
      #ref_bwt = ref_bwt,
      #ref_pac = ref_pac,
      #ref_sa = ref_sa
  }
}