# *************************************************
# Filtering Module for Read Quality and Coverage
# *************************************************

library(ShortRead)
library(Rsamtools)
library(GenomicRanges)

# FILTER FASTQ BY QUALITY
filter_fastq_by_quality <- function(
    fastq_file1,
    fastq_file2 = NULL,
    sample_name = "sample",
    read_type = "single",
    min_quality = 30
) {
  # Define output directory for filtered FASTQ files
  out_dir <- file.path("..", "Outputs", "Filtered_Fastq")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create timestamp for output filenames
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  
  # Initialize output paths
  output1 <- file.path(out_dir, paste0(sample_name, "_1_filtered_", timestamp, ".fastq"))
  output2 <- NULL
  
  # Read FASTQ1 and filter by average quality
  fq1 <- readFastq(fastq_file1)
  qual1 <- alphabetScore(fq1) / width(fq1)
  fq1_filtered <- fq1[qual1 >= min_quality]
  writeFastq(fq1_filtered, output1)
  
  # If paired-end, read FASTQ2 and apply same filtering
  if (tolower(read_type) == "paired" && !is.null(fastq_file2)) {
    fq2 <- readFastq(fastq_file2)
    qual2 <- alphabetScore(fq2) / width(fq2)
    fq2_filtered <- fq2[qual2 >= min_quality]
    output2 <- file.path(out_dir, paste0(sample_name, "_2_filtered_", timestamp, ".fastq"))
    writeFastq(fq2_filtered, output2)
  }
  
  message("FASTQ quality filtering complete. Outputs saved to ", out_dir)
  return(list(fastq1 = output1, fastq2 = output2))
}

# FILTER BAM BY COVERAGE
filter_bam_by_coverage <- function(bam_file,sample_name   = "sample26",min_coverage  = 10)
{
  # Set ups of outputs folders and output file names
  out_dir <- file.path("..", "Outputs", "Filtered_Bam")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  bam <- BamFile(bam_file)
  idxstats <- idxstatsBam(bam)
  
  # Erasing row with seqname == "*" (not corresponding to any valid chromosome)
  idxstats <- idxstats[idxstats$seqnames != "*", ]
  
  # Seleccionar únicamente aquellos cromosomas cuya cantidad de lecturas mapeadas >= min_coverage
  keep_regions <- idxstats[idxstats$mapped >= min_coverage, "seqnames"]
  
  if (length(keep_regions) == 0) {
    warning("No regions meet the minimum coverage threshold.")
    return(NULL)
  }
  gr_keep <- GRanges(seqnames = keep_regions, ranges   = IRanges(start = 1, end = .Machine$integer.max))
  param <- ScanBamParam(which = gr_keep, what  = "qname")
  aln <- scanBam(bam, param = param)[[1]]
  
  if (length(aln$qname) == 0) {
    warning("No reads passed the coverage threshold.")
    return(NULL)
  }
  
  # Buildin filtered BAM
  filtered_bam <- file.path(out_dir,paste0(sample_name, "_filtered_", timestamp, ".bam"))
  
  # Using Samtools to keep only Qnames that passed filter
  qnames_pattern <- paste0(aln$qname, collapse = "|")
  cmd <- paste(
    "samtools view -h", shQuote(bam_file), "|",
    "grep -E '^@|", qnames_pattern, "' |",
    "samtools view -b -o", shQuote(filtered_bam)
  )
  system(cmd)
  message("BAM coverage filtering complete. Output saved to ", filtered_bam)
  return(filtered_bam)
}

# WRAPPER FUNCTIONS
# Wrapper for FASTQ filtering only
run_fastq_filtering <- function(fastq1,fastq2 = NULL,sample_name = "sample",read_type = "single",quality_threshold = 30
) {
  return(filter_fastq_by_quality(fastq_file1 = fastq1,fastq_file2 = fastq2,sample_name = sample_name,read_type   = read_type,
    min_quality = quality_threshold
  ))
}

# Wrapper for BAM filtering only
run_bam_filtering <- function(bam_file,sample_name = "sample",min_coverage = 10) {
  return(filter_bam_by_coverage(bam_file= bam_file,sample_name=sample_name,min_coverage=min_coverage))
}
# run_bam_filtering (
#   "C:/Users/newli/OneDrive/Documentos/TFM/Outputs/Alignment_Evaluation/sample26/sample26_sorted.bam",
#    "sample12",
#   10
# )
