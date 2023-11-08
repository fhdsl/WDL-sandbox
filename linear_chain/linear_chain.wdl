version 1.0
workflow LinearChainExample {
  input {
    File originalBAM
  }
  call MarkDuplicates { 
    input: 
      bamFile = originalBAM 
  }
  call HaplotypeCaller { 
    input: 
      bamFile = MarkDuplicates.dedupBam 
  }
  call SelectVariants { 
    input: 
      VCF = HaplotypeCaller.rawVCF, 
      type="INDEL" 
  }
}

task MarkDuplicates {
  input {
    File bamFile
  }
  # command <<<
  #   java -jar picard.jar MarkDuplicates \
  #   I=~{bamFile} O=dedupped.bam M= dedupped.metrics
  # >>>
  command <<<
    cat ~{bamFile} > dedupped.bam
  >>>
  output {
    File dedupBam = "dedupped.bam"
    #File metrics = "dedupped.metrics"
  }
}

task HaplotypeCaller {
  input {
    File bamFile
  }
  # command <<<
  #   java -jar GenomeAnalysisTK.jar \
  #     -T HaplotypeCaller \
  #     -R reference.fasta \
  #     -I ~{bamFile} \
  #     -o rawVariants.vcf
  # >>>
  command <<<
    cat ~{bamFile} > rawVariants.vcf
  >>>
  output {
    File rawVCF = "rawVariants.vcf"
  }
}

task SelectVariants {
  input {
    File VCF
    String type
  }
  # command <<<
  #   java -jar GenomeAnalysisTK.jar \
  #     -T SelectVariants \
  #     -R reference.fasta \
  #     -V ~{VCF} \
  #     --variantType ~{type} \
  #     -o rawIndels.vcf
  # >>>
  command <<<
    cat ~{VCF} > rawIndels.vcf
  >>>
  output {
    File subsetVCF = "rawIndels.vcf"
  }
}