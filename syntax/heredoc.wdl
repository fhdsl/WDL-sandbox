version 1.0

task heredoc_syntax {
    input {
        String some_input = "foo"
    }
    
    String declared_elsewhere = "bar"
    
    command <<<
    BASH_VARIABLE="bizz"
    
    echo "This is heredoc syntax."
    echo "Here is the input variable: ~{some_input}"
    echo "Here is a variable declared outside the input section: ~{declared_elsewhere}"
    echo "Here is a variable defined in the command section: $BASH_VARIABLE"
    >>>

}