# ***************************************
# Main Script: Bioinformatics Pipeline
# ***************************************

# Load modules
source("Reduce_Sampling_module.R")
source("Quality_Control_module.R")
source("Filtering_module.R")
source("Alignment_module.R")
source("Alignment_Evaluation_module.R")
source("Variant_Calling_module.R")
source("Annotation_Module.R")


# Activating args[]
args <- commandArgs(trailingOnly = TRUE)

# Getting args[]
get_param <- function(i, line_text, default = NULL) {
  if (interactive()) {
    input <- readline(paste0(line_text, if (!is.null(default)) paste0(" [", default, "]"), ": "))
    return(ifelse(input == "", default, input))
  } else {
    return(ifelse(length(args) >= i, args[[i]], default))
  }
}

# Saving args[] in variables
readfile1   <- get_param(1, "File path 1")
readfile2   <- get_param(2, "File path 2", default = NULL)
sample_name <- get_param(3, "Sample name", default = "sample")
read_type   <- get_param(4, "Read type (single/paired)", default = "single")
file_type   <- get_param(5, "Type of file (fastq, fq)", default = "fastq")
ref_genome  <- get_param(6, "Genome of reference")
do_subsample <- get_param(7, "Do subsampling? (yes/no)", default = "no")

# Derived names
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
subsampled_fastq1 <- file.path("../Data/SubSampling", paste0(sample_name, "_1_subsample.fastq"))
subsampled_fastq2 <- file.path("../Data/SubSampling", paste0(sample_name, "_2_subsample.fastq"))


# Step 0: Subsampling if requested
if (tolower(do_subsample) == "yes") {
  run_subsampling(readfile1, readfile2, sample_name, 100000, out_dir = "../Data/SubSampling")
  readfile1 <- subsampled_fastq1
  if (!is.null(readfile2)) {
    readfile2 <- subsampled_fastq2
  }
}

if (tolower(file_type) == "fq") {
  file_type <- "fastq"
}

# Run QC
run_qc(readfile1, readfile2, sample_name, read_type, file_type)

# Ask user if they want to filter after quality control
do_filtering_qc <- get_param(8, "Do you want to filter the reads? (yes/no)", default = "no")
if (tolower(do_filtering_qc) == "yes") {
  quality_threshold <- as.numeric(get_param(9, "Minimum average quality (Phred)", default = "30"))
  min_length <- as.numeric(get_param(10, "Minimum read length", default = "30"))
  filtered <- run_filtering(readfile1, readfile2, sample_name, read_type, quality_threshold, min_length)
  readfile1 <- filtered$fastq1
  if (read_type == "paired") readfile2 <- filtered$fastq2
}

# Alignment
aligned_bam <- run_alignment(readfile1, readfile2, ref_genome, sample_name, read_type)
run_alignment_evaluation(aligned_bam)



# Ask user if they want to filter after alignment
do_filtering_bam <- get_param(11, "Do you want to filter the reads? (yes/no)", default = "no")
if (tolower(do_filtering_bam) == "yes") {
  min_cov <- as.numeric(get_param(12, "Minimum coverage", default = "10"))
  aligned_bam <- filter_bam_by_coverage(aligned_bam, sample_name, min_coverage = min_cov)
}

# Variant Calling and Annotation
vcf <- run_variant_calling(aligned_bam, ref_genome)
## annot <- run_annotation(vcf)

