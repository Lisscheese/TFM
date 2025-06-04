# **************************
# Alignment Module
# **************************

library(Rsubread)

# Align FASTQ reads to a reference genome locally stored and return the path to the resulting BAM file.
# Uses an index in Data/Prebuilt_Indexes/ref_index, building it only if missing.
run_alignment <- function(readfile1, readfile2 = NULL, ref_genome,sample_name,read_type = "single")
{
  # Init of folder where we store shared index, I added this step because it is usefull when
  # we are doing different analysis against the same genome_reference.
  index_dir <- file.path("..", "Data", "Prebuilt_Indexes")
  dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)
  index_prefix <- file.path(index_dir, "ref_index")
  index_prefix <- normalizePath(index_prefix, winslash = "/", mustWork = FALSE)
  
  # Checking if the index is already built, to build it only when necessary
  ejemplo_array <- paste0(index_prefix, ".00.b.array")
  if (!file.exists(ejemplo_array)) {
    message("Índice no encontrado en '", ejemplo_array, "'. Construyendo índice compartido...")
    buildindex(basename = index_prefix, reference = ref_genome)
    message("Índice construido correctamente en: ", index_prefix)
  } else {
    message("Índice ya existe en: ", index_prefix, ". No será reconstruido.")
  }
  
  # Initializing outputs folder for Alignment
  output_dir <- file.path("..", "Outputs", "Alignment", sample_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Setting the timestamp and name of output file, timestamp are added
  # to differ between different analysi of same sample.
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  output_file <- file.path(
    output_dir,
    paste0(sample_name, "_aligned_", timestamp, ".bam")
  )
  
  # Calling Alig function of Rsubread
  if (tolower(read_type) == "paired" && !is.null(readfile2)) {
    message("Running alignment in paired-end mode...")
    align( index = index_prefix,readfile1 = readfile1, readfile2   = readfile2,output_file = output_file)
  } else {
    message("Running alignment in single-end mode...")
    align(index = index_prefix, readfile1 = readfile1,output_file = output_file)
  }
  
  message("Alignment completed. Output saved to ", output_file)
  return(output_file)
}


#run_alignment(
#  "C:/Users/newli/OneDrive/Documentos/TFM/Outputs/SubSampling/sample_7_1_subsample.fastq",
#  "C:/Users/newli/OneDrive/Documentos/TFM/Outputs/SubSampling/sample_7_2_subsample.fastq",
#  "C:/Users/newli/OneDrive/Documentos/TFM/Data/genome.fa",
#  "sample24",
#  "paired"
# )

