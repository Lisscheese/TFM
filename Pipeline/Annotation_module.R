# *****************************************
# Annotation_module.R
# ******************************************

# Loading specific libraries 
library(S4Vectors)
library(IRanges)
library(GenomicRanges)
library(VariantAnnotation)
library(biomaRt)

# Path converter function same as variant calling module (/mnt/c/…).
windows_to_wsl <- function(win_path, must_exist = TRUE) {
  abs_win <- normalizePath(win_path, winslash = "/", mustWork = must_exist)
  drive_letter <- tolower(substr(abs_win, 1, 1))      
  rest_path <- substr(abs_win, 3, nchar(abs_win))     
  rest_path <- gsub("\\\\", "/", rest_path)           
  rest_path <- sub("^/", "", rest_path)               
  paste0("/mnt/", drive_letter, "/", rest_path)       
}

run_snpEff <- function( vcf_file, sample_name, genome_db   = "GRCh38.86", snpeff_dir  = "/mnt/c/Users/newli/Downloads/snpEff_latest_core/snpEff")
{
  # Checking if vcf exists
  if (!file.exists(vcf_file)) {
    stop("error: VCF file not found:\n  ", vcf_file)
  }
  
  # Creating output folder if not created
  out_dir <- file.path("..", "Outputs", "Variant_Annotation_SnpEff", sample_name)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Outputs names of files, building
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  annotated_vcf <- file.path(out_dir, paste0(sample_name, "_snpeff_annotated_", timestamp, ".vcf"))
  stats_file <- file.path(out_dir,paste0(sample_name, "_snpeff_report_", timestamp, ".html"))
  
  # Converting paths of windows to wsl paths
  vcf_wsl       <- windows_to_wsl(vcf_file, must_exist = TRUE)
  annotated_wsl <- windows_to_wsl(annotated_vcf, must_exist = FALSE)
  stats_wsl     <- windows_to_wsl(stats_file, must_exist = FALSE)
  
  # Building path to JAR of SnpEff in WSL
  snpeff_jar <- file.path(snpeff_dir, "snpEff.jar")
  
  # Checking that snpEff.jar is installed in wsl
  check_cmd <- paste("wsl", shQuote(paste0("[ -f ", snpeff_jar, " ] && echo 'OK' || echo 'MISSING'")))
  chk <- system(check_cmd, intern = TRUE)
  if (!any(grepl("OK", chk))) {
    stop("ERROR: snpEff.jar not found in path:\n  ", snpeff_jar)
  }
  
  # Building and calling the command to make the anotation
  # with command java -Xmx4g -jar snpEff.jar <genome_db> <VCF> > <annotated.vcf> -> example in docs of snpEff
  snpEff_cmd <- paste("java -Xmx4g -jar", snpeff_jar, genome_db, vcf_wsl, ">", annotated_wsl)
  full_cmd <- paste("wsl bash -lc", shQuote(snpEff_cmd))
  message("Calling SnpEff for annotation...")
  exit_code <- system(full_cmd)
  if (exit_code != 0) {
    stop("error: SnpEff exit code", exit_code)
  }
  
  # Generating the html report
  # with command java -jar snpEff.jar stats <annotated.vcf> > <stats.html> -> example in docs of snpEff
  stats_cmd <- paste("java -Xmx4g -jar", snpeff_jar, "stats", annotated_wsl, ">", stats_wsl)
  full_stats_cmd <- paste("wsl bash -lc", shQuote(stats_cmd))
  system(full_stats_cmd)
  
  message(
    "SnpEff completed\n",
    "Annotated VCF saved in", annotated_wsl, "\n",
    "HTML report saved in: ", stats_wsl
  )

  return(list(annotated_vcf = annotated_wsl, stats_file = stats_wsl))
}


# 
# resultado <- run_snpEff(
#   vcf_file     = "C:/Users/newli/OneDrive/Documentos/TFM/Outputs/Variant_Calling/sample29/sample29_variants_20250602_1254.vcf",
#   sample_name = "sample29",
#   genome_db   = "GRCh38.86",    # Base de datos descargada en snpEff/data/
#   snpeff_dir  = "/mnt/c/Users/newli/Downloads/snpEff_latest_core/snpEff"
# )

