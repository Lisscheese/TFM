# ***********************
# Quality Control Module
# ***********************

library(ShortRead)

# This function generates a html report with many metrics provided by the qa function of ShortRear
# most of the steps done in this function are to there to handle the input and output files 
run_qc <- function(readfile1, readfile2 = NULL, sample_name, read_type = "single", file_type) 
{
  # Creating path of output folder.
  qc_dir <- file.path("..", "Outputs", "Quality_Control", sample_name)
  
  # Creating folder if this does not exist (for new instalations)
  if (!dir.exists(qc_dir)) {
    dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")  
  # Calling function qa
  qa_summary1 <- qa(readfile1, type = file_type)
  report(qa_summary1,dest = file.path(qc_dir, paste0(sample_name, "_1_QC_", timestamp)))

  # In case type of read is paired, the second file is processed.
  if (tolower(read_type) == "paired" && !is.null(readfile2)) {
    qa_summary2 <- qa(readfile2, type = file_type)
    report( qa_summary2, dest = file.path( qc_dir, paste0(sample_name, "_2_QC_", timestamp)))
  }
}
