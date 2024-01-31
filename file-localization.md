

## Dude, where's my files?

```
task diffdiff_backmask {
    input {
        Array[File] diffs
    }
    command <<<
    wget https://raw.githubusercontent.com/aofarrel/diffdiff/main/diffdiff.py
    
    echo "~{sep='\n' diffs}" >> diff_paths.txt
    
    python3 diffdiff.py diff_paths.txt -b
    
    >>>
    runtime {
		cpu: 4
		disks: "local-disk " + 10 + " HDD"
		docker: "ashedpotatoes/sranwrp:1.1.15"
		memory: "8 GB"
		preemptible: 2
	}
    output {
        Array[File] backmasked_diffs = glob("*.backmask.diff")
    }
}
```

End result: backmasked_diffs is an empty array (oddly enough on miniwdl this does not error?)