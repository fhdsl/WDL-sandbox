version 1.0
## Notes: Testing if alignement for 1 sample will work
## Samples:NCIH2172_LUNG : Cell-line with EGFR L858R mutation
##
## Input requirements:
## - fastq files for reference and alternative
##
## Output Files:
## - An aligned bam for 1 sample
## 
## Workflow developed by Sitapriya Moorthi @ Fred Hutch LMD: 01/10/24 for use by DaSL @ Fred Hutch.

workflow minidata_mutation_calling_v1 {
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
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    File af_only_gnomad
    File af_only_gnomad_index
    File annovarTAR
    String annovar_protocols
    String annovar_operation
    String ref_name
  }
    
  # Map reads to reference
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
  
  call MarkDuplicates {
    input:
      input_bam = BwaMem.analysisReadySorted
  }

  call ApplyBaseRecalibrator {
    input:
      input_bam = MarkDuplicates.markDuplicates_bam,
      input_bam_index = MarkDuplicates.markDuplicates_bai,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
    }

  call Mutect2TumorOnly {
      input:
        input_bam = ApplyBaseRecalibrator.recalibrated_bam,
        input_bam_index = ApplyBaseRecalibrator.recalibrated_bai,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        genomeReference = af_only_gnomad,
        genomeReferenceIndex = af_only_gnomad_index
    }
  
  call annovar {
      input:
        input_vcf = Mutect2TumorOnly.output_vcf,
        ref_name = ref_name,
        annovarTAR = annovarTAR,
        annovar_operation = annovar_operation,
        annovar_protocols = annovar_protocols
    }
  
  # Outputs that will be retained when execution is complete
  output {
    File alignedBamSorted = BwaMem.analysisReadySorted
    File markDuplicates_bam = MarkDuplicates.markDuplicates_bam
    File markDuplicates_bai = MarkDuplicates.markDuplicates_bai
    File analysisReadyBam = ApplyBaseRecalibrator.recalibrated_bam 
    File analysisReadyIndex = ApplyBaseRecalibrator.recalibrated_bai
    File Mutect_Vcf = Mutect2TumorOnly.output_vcf
    File Mutect_VcfIndex = Mutect2TumorOnly.output_vcf_index
    File Mutect_AnnotatedVcf = annovar.output_annotated_vcf
    File Mutect_AnnotatedTable = annovar.output_annotated_table
  }
}






# TASK DEFINITIONS

# Align fastq file to the reference genome
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
  String read_group_id = "ID:" + base_file_name
  String sample_name = "SM:" + base_file_name
  String platform = "illumina"
  String platform_info = "PL:" + platform   # Create the platform information


  command <<<
    set -eo pipefail

    bwa mem \
      -p -v 3 -t 16 -M -R '@RG\t~{read_group_id}\t~{sample_name}\t~{platform_info}' \
      ~{ref_fasta} ~{input_fastq} > ~{base_file_name}.sam 
    samtools view -1bS -@ 15 -o ~{base_file_name}.aligned.bam ~{base_file_name}.sam
    samtools sort -@ 15 -o ~{base_file_name}.sorted_query_aligned.bam ~{base_file_name}.aligned.bam
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

# Mark duplicates (not SPARK, for some reason that does something weird)
task MarkDuplicates {
  input {
    File input_bam
  }

  String base_file_name = basename(input_bam, ".sorted_query_aligned.bam")
  String output_bam = "~{base_file_name}.duplicates_marked.bam"
  String output_bai = "~{base_file_name}.duplicates_marked.bai"
  String metrics_file = "~{base_file_name}.duplicate_metrics"

  command <<<
    gatk MarkDuplicates \
      --INPUT ~{input_bam} \
      --OUTPUT ~{output_bam} \
      --METRICS_FILE ~{metrics_file} \
      --CREATE_INDEX true \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 \
      --VALIDATION_STRINGENCY SILENT
  >>>

  runtime {
    docker: "broadinstitute/gatk:4.1.4.0"
    memory: "48 GB"
    cpu: 4
  }

  output {
    File markDuplicates_bam = "~{output_bam}"
    File markDuplicates_bai = "~{output_bai}"
    File duplicate_metrics = "~{metrics_file}"
  }
}

# Base quality recalibration
task ApplyBaseRecalibrator {
  input {
    File input_bam
    File input_bam_index
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index
  }
  
  String base_file_name = basename(input_bam, ".duplicates_marked.bam")

  command <<<
  set -eo pipefail

  samtools index ~{input_bam}

  gatk --java-options "-Xms8g" \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -O ~{base_file_name}.recal_data.csv \
      --known-sites ~{dbSNP_vcf} \
      --known-sites ~{sep=" --known-sites " known_indels_sites_VCFs} \
      

  gatk --java-options "-Xms8g" \
      ApplyBQSR \
      -bqsr ~{base_file_name}.recal_data.csv \
      -I ~{input_bam} \
      -O ~{base_file_name}.recal.bam \
      -R ~{ref_fasta} \
      

  #finds the current sort order of this bam file
  samtools view -H ~{base_file_name}.recal.bam | grep @SQ | sed 's/@SQ\tSN:\|LN://g' > ~{base_file_name}.sortOrder.txt
>>>

  output {
    File recalibrated_bam = "~{base_file_name}.recal.bam"
    File recalibrated_bai = "~{base_file_name}.recal.bai"
    File sortOrder = "~{base_file_name}.sortOrder.txt"
  }
  runtime {
    memory: "36 GB"
    cpu: 2
    docker: "broadinstitute/gatk:4.1.4.0"
  }
}






# Mutect 2 calling

task Mutect2TumorOnly {
  input {
    File input_bam
    File input_bam_index
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File genomeReference
    File genomeReferenceIndex
  }

    String base_file_name = basename(input_bam, ".recal.bam")

command <<<
    set -eo pipefail

    gatk --java-options "-Xms16g" Mutect2 \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -O preliminary.vcf.gz \
      --germline-resource ~{genomeReference} \
     
    gatk --java-options "-Xms16g" FilterMutectCalls \
      -V preliminary.vcf.gz \
      -O ~{base_file_name}.mutect2.vcf.gz \
      -R ~{ref_fasta} \
      --stats preliminary.vcf.gz.stats \
     
>>>

runtime {
    docker: "broadinstitute/gatk:4.1.4.0"
    memory: "24 GB"
    cpu: 1
  }

output {
    File output_vcf = "${base_file_name}.mutect2.vcf.gz"
    File output_vcf_index = "${base_file_name}.mutect2.vcf.gz.tbi"
  }

}




# annotate with annovar
task annovar {
  input {
  File input_vcf
  String ref_name
  File annovarTAR
  String annovar_protocols
  String annovar_operation
}
  String base_vcf_name = basename(input_vcf, ".vcf.gz")
  
  command <<<
  set -eo pipefail
  
  tar -xzvf ~{annovarTAR}
  
  perl annovar/table_annovar.pl ~{input_vcf} annovar/humandb/ \
    -buildver ~{ref_name} \
    -outfile ~{base_vcf_name} \
    -remove \
    -protocol ~{annovar_protocols} \
    -operation ~{annovar_operation} \
    -nastring . -vcfinput
>>>
  runtime {
    docker : "perl:5.28.0"
    cpu: 1
    memory: "2GB"
  }
  output {
    File output_annotated_vcf = "${base_vcf_name}.${ref_name}_multianno.vcf"
    File output_annotated_table = "${base_vcf_name}.${ref_name}_multianno.txt"
  }
}


