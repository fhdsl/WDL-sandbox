version 1.0

workflow Test_File_Localization {
  input {
    Array[File] many_files  # two or more files
    File move_this_file     # any file
    File bar                # a file literally called "bar.txt" that contains the text "bar!"
    File? copy_this_file    # any file
    File? undefined         # leave this one undefined
  }
  
  call test_localization {
    input:
      many_files = many_files,
      move_this_file = move_this_file,
      bar = bar,
      copy_this_file = copy_this_file,
      undefined = undefined
  }
}

task test_localization {
  input {
    Array[File] many_files  # two or more files
    File move_this_file     # any file
    File bar                # a file literally called "bar.txt" that contains the text "bar!"
    File? copy_this_file    # any file
    File? undefined         # leave this one undefined
  }

  command <<<
  seperator="----------------------------------------------------------"
  echo "workdir is $(pwd)" &>> outfile
  
  echo $seperator &>> outfile
  ls -lha &>> outfile
  echo $seperator &>> outfile
  
  echo "move_this_file: ~{move_this_file}" &>> outfile
  echo "bar: ~{bar}" &>> outfile
  echo "copy_this_file: ~{copy_this_file}" &>> outfile
  echo "undefined: ~{undefined}" &>> outfile
  echo "array: ~{sep=' ' many_files}" &>> outfile
  echo $seperator &>> outfile
  
  echo "cat bar as if it's in workdir" &>> outfile
  cat bar.txt &>> outfile
  echo $seperator &>> outfile
  
  echo "copy copy_this_file into workdir" &>> outfile
  cp "~{copy_this_file}" . &>> outfile  # FYI: stderr redirect changes "Read-only file system" into "Device or resource busy"
  ls -lha &>> outfile
  echo $seperator &>> outfile
  
  echo "move move_this_file into workdir" &>> outfile
  mv "~{move_this_file}" . &>> outfile
  ls -lha &>> outfile
  echo $seperator &>> outfile
  
  echo "directly edit the original copy_this_file" &>> outfile
  echo "HELLO WORLD" > "~{copy_this_file}"
  cat "~{copy_this_file}" &>> outfile
  >>>

  runtime { docker: "ubuntu:latest" } # bad practice, but okay for a test
  
  output { File stderr_and_stdout = "outfile" }
}