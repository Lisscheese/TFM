# ****************************
# Subsampling Module
# ****************************
library(ShortRead)

# Subsamples the first N reads from large FASTQ files using function of the library 
# ShortRead: FastqStreamer()
run_subsampling <- function( readfile1, readfile2 = NULL, sample_name = "sample", n = 10000, out_dir = file.path("..","Outputs", "SubSampling"))
{
  message("Subsampling ", n, " reads from FASTQ")
  
  # Creating the output directory
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Creating the name of the output file
  output1 <- file.path(out_dir, paste0(sample_name, "_1_subsample.fastq"))
  output2 <- NULL
  
  # Read and write R1 using streaming, how the n parameter works is that it is the block-size that will take 
  # in each yield call.
  fq_stream1 <- FastqStreamer(readfile1, n = n)
  fq1_sub <- yield(fq_stream1)
  writeFastq(fq1_sub, output1)
  close(fq_stream1)
  
  # If the second file is not null then it means that is paired-end
  # so, we can process readfile2
  if (!is.null(readfile2)) {
    fq_stream2 <- FastqStreamer(readfile2, n = n)
    fq2_sub <- yield(fq_stream2)
    output2 <- file.path(out_dir, paste0(sample_name, "_2_subsample.fastq"))
    writeFastq(fq2_sub, output2)
    close(fq_stream2)
  }
  
  message("Subsampling completed. Output files are in ", out_dir)
  
  # Returning to the main flow, the file paths of generated fastq files
  return(list(fastq1 = output1, fastq2 = output2))
}
