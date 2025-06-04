# ****************************************
# Variant_Calling_module.R
# ****************************************

run_variant_calling <- function(
    bam_file,
    reference_fasta,
    sample_name
) {
  # Output folders set ups
  out_dir <- file.path("..", "Outputs", "Variant_Calling", sample_name)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Timestamp win path to generate vcf as result
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  output_vcf <- file.path( out_dir, paste0(sample_name, "_variants_", timestamp, ".vcf")
  )
  
  # Error handeling when relevant files does not exist
  if (!file.exists(bam_file)) {
    stop("error: BAM file does not exist:\n  ", bam_file)
  }
  if (!file.exists(reference_fasta)) {
    stop("error: reference file does not exist:\n  ", reference_fasta)
  }
  
  # Error handeling if wsl not installed in windows
  wsl_exe <- Sys.which("wsl")
  if (wsl_exe == "") {
    stop("ERROR: wsl not available in system.\n",
         "Please install wsl and test from R console.")
  }
  
  # Helper function to convert window path \ to wsl path /
  windows_to_wsl <- function(win_path, must_exist = TRUE) {
    # "C:\\a\\b" in "C:/a/b"
    abs_win <- normalizePath(win_path, winslash = "/", mustWork = must_exist)
    drive_letter <- tolower(substr(abs_win, 1, 1))
    rest_path <- substr(abs_win, 3, nchar(abs_win))   
    rest_path <- gsub("\\\\", "/", rest_path)         
    rest_path <- sub("^/", "", rest_path)            
    paste0("/mnt/", drive_letter, "/", rest_path)     
  }
  
  # Applying the changes to bam and fastq paths
  ref_wsl <- windows_to_wsl(reference_fasta, must_exist = TRUE)
  bam_wsl <- windows_to_wsl(bam_file, must_exist = TRUE)
  
  # Turning path of VCF output to wsl
  # using mustWork = FALSE to avoid error when vcf is not created yet, important to keep 
  abs_win_vcf <- normalizePath(output_vcf, winslash = "/", mustWork = FALSE)
  drive_letter <- tolower(substr(abs_win_vcf, 1, 1))
  rest_vcf <- substr(abs_win_vcf, 3, nchar(abs_win_vcf))
  rest_vcf <- gsub("\\\\", "/", rest_vcf)
  rest_vcf <- sub("^/", "", rest_vcf)
  vcf_wsl <- paste0("/mnt/", drive_letter, "/", rest_vcf)
  
  # Checking that freebayes is installed 
  check_cmd <- "which freebayes"
  cmd_check <- paste("wsl bash -lc", shQuote(check_cmd))
  path_freebayes <- tryCatch(
    system(cmd_check, intern = TRUE),
    error   = function(e) character(0),
    warning = function(w) character(0)
  )
  if (length(path_freebayes) == 0 || grepl("^which: no freebayes", path_freebayes[1])) {
    stop(
      "error: FreeBayes not installed in wsl.\n",
    )
  }
  
  # Building FreeBayes command to pass to the wsl
  # Command to execute :wsl bash -lc "freebayes -f '/mnt/.../genome.fa' '/mnt/.../sample.bam' > '/mnt/.../sample.vcf'"
  freebayes_cmd <- sprintf( "freebayes -f %s %s > %s", shQuote(ref_wsl), shQuote(bam_wsl), shQuote(vcf_wsl))
  full_cmd <- paste("wsl bash -lc", shQuote(freebayes_cmd))
  
  # Calling FreeBayes inside wsl and getting the exit code, 
  # this is different when working in Linux environemt
  message("Ejecutando FreeBayes dentro de WSL...")
  exit_code <- system(full_cmd)
  if (exit_code != 0) {
    stop("ERROR: FreeBayes exit with code: ", exit_code, ".\n")
  }
  
  message("Variant Calling completed. VCF generated: ", output_vcf)
  return(output_vcf)
}

# 
# run_variant_calling(
# "C:/Users/newli/OneDrive/Documentos/TFM/Outputs/Alignment_Evaluation/sample26/sample26_sorted.bam",
#  "C:/Users/newli/OneDrive/Documentos/TFM/Data/genome.fa",
#  "sample26"
# )
#  