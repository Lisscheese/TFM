# *****************************************
# Main Script: Bioinformatics Pipeline
# *****************************************

# Load modules
source("Reduce_Sampling_module.R")
source("Quality_Control_module.R")
source("Filtering_module.R")
source("Alignment_module.R")
source("Alignment_Evaluation_module.R")
source("Variant_Calling_module.R")
source("Annotation_Module.R")

# Preparing arguments
args <- commandArgs(trailingOnly = TRUE)
get_param <- function(i, text, default = NULL) {
  if (interactive()) {
    resp <- readline(paste0(text, if (!is.null(default)) paste0(" [", default, "]"), ": "))
    return(ifelse(resp == "", default, resp))
  }
  if (length(args) >= i) return(args[[i]])
  return(default)
}

# Reading parameters
readfile1    <- get_param(1, "FASTQ R1 path")
readfile2    <- get_param(2, "FASTQ R2 path", default = NULL)
sample_name  <- get_param(3, "Sample name", default = "sample")
read_type    <- tolower(get_param(4, "Read type (single/paired)", default = "single"))
file_type    <- tolower(get_param(5, "File type (fastq/fq)",    default = "fastq"))
ref_genome   <- get_param(6, "Reference genome FASTA path")
do_subsamp   <- tolower(get_param(7, "Do you want to subsample? (yes/no)", default="no"))

# Subsampling
if (do_subsamp == "yes") {
  subs <- run_subsampling(
    readfile1, readfile2,
    sample_name = sample_name,
    n           = 100000,
    out_dir     = file.path("Outputs","SubSampling")
  )
  readfile1 <- subs$fastq1
  readfile2 <- subs$fastq2
}

# QC
run_qc(
  readfile1, readfile2,
  sample_name = sample_name,
  read_type   = read_type,
  file_type   = file_type
)

# Filtering FASTQ
do_filter_qc <- tolower(get_param(8, "Filter reads by quality? (yes/no)", default="no"))
if (do_filter_qc == "yes") {
  q_thr <- as.numeric(get_param(9, "Min average quality", default="30"))
  fastq_res <- run_filtering(
    fastq1           = readfile1,
    fastq2           = readfile2,
    sample_name      = sample_name,
    read_type        = read_type,
    quality_threshold= q_thr,
    do_bam_filter    = FALSE
  )
  readfile1 <- fastq_res$fastq$fastq1
  readfile2 <- fastq_res$fastq$fastq2
}

# Alignment
aligned_bam <- run_alignment(
  readfile1, readfile2,
  ref_genome  = ref_genome,
  sample_name = sample_name,
  read_type   = read_type
)

# Alignemnt Evaluation
eval_res   <- run_alignment_evaluation(aligned_bam, sample_name)
sorted_bam <- eval_res$sorted_bam

# Filtering BAM
do_filter_bam <- tolower(get_param(11, "Filter BAM by coverage? (yes/no)", default="no"))
if (do_filter_bam == "yes") {
  min_cov <- as.numeric(get_param(12, "Min coverage", default="10"))
  sorted_bam <- filter_bam_by_coverage(
    bam_file    = sorted_bam,
    sample_name = sample_name,
    min_coverage= min_cov
  )
}

# Variant Calling
vcf <- run_variant_calling(
  bam_file       = sorted_bam,
  reference_fasta= ref_genome,
  sample_name    = sample_name
)

# Annotation
annot_res <- run_annotation(
  vcf_file    = vcf,
  genome_build= "hg38",
  sample_name = sample_name
)
