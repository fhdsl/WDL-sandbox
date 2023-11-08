# Example workflow
# Declare WDL version 1.0 if working in Terra
version 1.0
workflow myWorkflow {
    call myTask
}

task myTask {
    command <<<
        echo "hello world"
    >>>
    output {
        String out = read_string(stdout())
    }
}