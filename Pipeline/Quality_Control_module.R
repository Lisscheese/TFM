# --------------------------
# Quality Control Module
# --------------------------

library(ShortRead)

# Performs quality assessment on FASTQ files and generates HTML/CSV reports
run_qc <- function(
    readfile1,
    readfile2 = NULL,
    sample_name,
    read_type = "single",
    file_type
) 
{
  # Create timestamp for filenames
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  
  # Define QC output directory
  qc_dir <- file.path("..", "Outputs", "Quality_Control", sample_name)
  
  # Create output directory recursively if it does not exist
  if (!dir.exists(qc_dir)) {
    dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Perform QA on readfile1 and generate report
  qa_summary1 <- qa(readfile1, type = file_type)
  report(
    qa_summary1,
    dest = file.path(
      qc_dir,
      paste0(sample_name, "_1_QC_", timestamp)
    )
  )
  
  # If paired-end, perform QA on readfile2 and generate report
  if (tolower(read_type) == "paired" && !is.null(readfile2)) {
    qa_summary2 <- qa(readfile2, type = file_type)
    report(
      qa_summary2,
      dest = file.path(
        qc_dir,
        paste0(sample_name, "_2_QC_", timestamp)
      )
    )
  }
}
