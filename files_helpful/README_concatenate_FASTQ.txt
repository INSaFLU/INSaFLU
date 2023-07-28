###### Concatenate FASTQ files on WINDOWS or MAC #######

- WINDOWS

Download one of these files = copy/paste the content to a NOTEPAD and save with the respective name:

  WINDOWS_concat_fastq.bat (if your files are fastq)

  WINDOWS_concat_fastqgz.bat (if your files are fastq.gz)

To concatenate several fastq/fastq.gz files into a single fastq/fastq.gz, do the following:

  - create a folder where you only put the fastq/fastq.gz you wish to concatenate

  - copy&paste the adequate “.bat” file (whether you have fastq or fastq.gz) to the same folder and double-click on the “.bat” file
  
This will automatically create a single file named “concat.fastq.gz” (or “concat.fastq”) inside the same folder. You can then rename this file as needed. (Note that this will not eliminate or change the original fastq inside the folder.)


- MAC

Download one of these files = copy/paste the content to a NOTEPAD and save with the respective name:

  MAC_concat_fastq.command (if your files are fastq)

  MAC_concat_fastqgz.command (if your files are fastq.gz)

To concatenate several fastq/fastq.gz files into a single fastq/fastq.gz, do the following:
  - create a folder where you only put the fastq/fastq.gz you wish to concatenate
  - copy&paste the adequate “.command” file (whether you have fastq or fastq.gz) to the same folder and double-click on the “.command” file
  
This will automatically create a single file named “concat.fastq.gz” (or “concat.fastq”) inside the same folder. You can then rename this file as needed. (Note that this will not eliminate or change the original fastq inside the folder.)


###### Concatenate FASTQ files on UNIX #######

You can use "findONTtime:

https://github.com/INSaFLU/findONTime (see example 1)

$ findontime -i input_directory -o output_directory --tag suffix --max_size 100000000 --merge --upload none
