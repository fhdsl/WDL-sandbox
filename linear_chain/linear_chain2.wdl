version 1.0

#do not hard-code output file names. 

workflow LinearChainExample {
  input {
    File originalBAM
    String id
  }
  call MarkDuplicates { 
    input: 
      bamFile = originalBAM, 
      outputName = id + ".markDuplicates.bam"
  }
  call HaplotypeCaller { 
    input: 
      bamFile = MarkDuplicates.dedupBam,
      outputName = id + ".vcf" 
  }
  call SelectVariants { 
    input: 
      VCF = HaplotypeCaller.rawVCF, 
      outputName = id + ".indels.vcf", 
      type="INDEL" 
  }
}

task MarkDuplicates {
  input {
    File bamFile
    String outputName
  }
  # command <<<
  #   java -jar picard.jar MarkDuplicates \
  #   I=~{bamFile} O=dedupped.bam M= dedupped.metrics
  # >>>
  command <<<
    cat ~{bamFile} > ~{outputName}
  >>>
  output {
    File dedupBam = outputName
    #File metrics = "dedupped.metrics"
  }
}

task HaplotypeCaller {
  input {
    File bamFile
    String outputName
  }
  # command <<<
  #   java -jar GenomeAnalysisTK.jar \
  #     -T HaplotypeCaller \
  #     -R reference.fasta \
  #     -I ~{bamFile} \
  #     -o rawVariants.vcf
  # >>>
  command <<<
    cat ~{bamFile} >  ~{outputName}
  >>>
  output {
    File rawVCF = outputName
  }
}

task SelectVariants {
  input {
    File VCF
    String type
    String outputName

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
    cat ~{VCF} > ~{outputName}
  >>>
  output {
    File subsetVCF = outputName
  }
}