version 1.0
workflow somaticMutationCalilng {
	input {
		File tumorBam
	}
	call SamToFastq {
		input: 
			input_bam = tumorBam,
      base_file_name = basename(tumorBam, ".unmapped.bam")
	}
	output {
		File fastq = SamToFastq.output_fastq
	}
}

task SamToFastq {
  input {
    File input_bam
    String base_file_name
  }
  command {
    set -eo pipefail

    gatk \
      SamToFastq \
      --INPUT=~{input_bam} \
      --FASTQ=~{base_file_name}.fastq.gz \
      --INTERLEAVE=true \
      --INCLUDE_NON_PF_READS=true 
  }
  output {
    File output_fastq = "~{base_file_name}.fastq.gz"
  }
  runtime {
    memory: "16 GB"
    docker: "broadinstitute/gatk:4.1.4.0"
    cpu: 1
    walltime: "2:00:00"
  }
}