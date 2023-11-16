version 1.0

#do not hard-code output file names. 
#now allows an array of inputs. 

workflow LinearChainExample {
  input {
    Array[File] originalBAM
    Array[String] id
  }
  scatter (idx in range(length(originalBAM))) {
    call MarkDuplicates { 
    input: 
        bamFile = originalBAM[idx], 
        outputName = id[idx] + ".markDuplicates.bam"
    }
  }
  scatter (idx in range(length(MarkDuplicates.dedupBam))) {
    call HaplotypeCaller { 
        input: 
        bamFile = MarkDuplicates.dedupBam[idx],
        outputName = id[idx] + ".vcf"
    }
  }
  scatter (idx in range(length(HaplotypeCaller.rawVCF))) {
    call SelectVariants { 
        input: 
        VCF = HaplotypeCaller.rawVCF[idx], 
        outputName = id[idx] + "indels.vcf",
        type="INDEL" 
    }
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
    String outputName
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
    cat ~{VCF} > ~{outputName}
  >>>
  output {
    File subsetVCF = outputName
  }
}