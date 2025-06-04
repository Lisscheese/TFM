# ******************************
# Alignment Evaluation Module 
# ******************************

library(Rsamtools)
library(GenomicAlignments)

# Evaluates an aligned BAM file and executes idxstats to generate statistics
run_alignment_evaluation <- function( bam_file, sample_name = tools::file_path_sans_ext(basename(bam_file))) 
{
  # Init of crea
  eval_dir <- file.path("..", "Outputs", "Alignment_Evaluation", sample_name)
  if (!dir.exists(eval_dir)) {
    dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)
  }
  # Calling sortBam and indexBam
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  sorted_bam_base <- file.path(eval_dir, paste0(sample_name, "_sorted"))
  sortBam(bam_file, sorted_bam_base)
  sorted_bam <- paste0(sorted_bam_base, ".bam")
  indexBam(sorted_bam)
  
  # Callning idxstats to generate statistics
  # idxstatsBam returns dfs with seqname, seqlength, mapped, unmapped counts
  idxstats_df <- idxstatsBam(sorted_bam)
  # Saving to csv statistics of idxstats
  idxstats_file <- file.path(eval_dir, paste0(sample_name, "_idxstats_", timestamp, ".csv"))
  write.csv(as.data.frame(idxstats_df), idxstats_file, row.names = FALSE)

  return(list(sorted_bam = sorted_bam, idxstats_file = idxstats_file))
}
# run_alignment_evaluation (
#   "C:/Users/newli/OneDrive/Documentos/TFM/Outputs/Alignment/sample26/sample26_aligned_20250601_0008.bam",
#     sample_name = "sample27"
# )
