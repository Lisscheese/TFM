# -------------------------------------
# Shiny Interface for Genomic Pipeline
# -------------------------------------

options(shiny.maxRequestSize = 88 * 1024 * 1024 * 1024 * 2)

library(shiny)
library(fs)
library(shinythemes)

# Load modules
source("Reduce_Sampling_module.R")
source("Quality_Control_module.R")
source("Filtering_module.R")
source("Alignment_module.R")
source("Alignment_Evaluation_module.R")
source("Variant_Calling_module.R")
source("Annotation_Module.R")

ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Genomic Pipeline - Interactive Interface"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("readfile1", "Read File 1 (R1)", accept = c(".fq", ".fastq")),
      fileInput("readfile2", "Read File 2 (R2, optional)", accept = c(".fq", ".fastq")),
      textInput("sample_name", "Sample Name", value = "sample"),
      selectInput("read_type", "Read Type", choices = c("single", "paired")),
      fileInput("ref_genome", "Reference Genome (FASTA)", accept = ".fa"),
      checkboxInput("do_subsample", "Perform Subsampling?", value = FALSE),
      numericInput("subsample_n", "Reads for Subsampling", value = 100000),
      checkboxInput("do_filter_qc", "Filter after QC", value = FALSE),
      numericInput("min_quality", "Min Average Quality", value = 30),
      numericInput("min_length", "Min Read Length", value = 30),
      checkboxInput("do_filter_bam", "Filter BAM by Coverage", value = FALSE),
      numericInput("min_coverage", "Min Coverage for BAM Filter", value = 10),
      actionButton("run_qc", "Step 1: Run Quality Control"),
      actionButton("run_filter_qc", "Step 2: Filter after QC"),
      actionButton("run_align", "Step 3: Alignment"),
      actionButton("run_filter_bam", "Step 4: Filter BAM"),
      actionButton("run_variant_call", "Step 5: Variant Calling"),
      verbatimTextOutput("file_status")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Console Log", verbatimTextOutput("log")),
        tabPanel("Alignment Summary", verbatimTextOutput("align_summary")),
        tabPanel("Variant Calls", verbatimTextOutput("variant_output"))
      )
    )
  )
)

server <- function(input, output, session) {
  log_text <- reactiveVal("")
  append_log <- function(msg) {
    log_text(paste(log_text(), msg, sep = "\n"))
  }
  
  output$file_status <- renderText({
    if (is.null(input$readfile1)) {
      "⏳ Waiting for Read File 1 upload..."
    } else {
      paste("✔ File uploaded:", input$readfile1$name)
    }
  })
  
  observeEvent(input$run_qc, {
    req(input$readfile1)
    withProgress(message = "Running Quality Control...", value = 0.1, {
      readfile1 <- input$readfile1$datapath
      readfile2 <- if (!is.null(input$readfile2)) input$readfile2$datapath else NULL
      run_qc(readfile1, readfile2, input$sample_name, input$read_type)
      incProgress(1)
    })
    append_log("✔ Quality Control completed.")
  })
  
  observeEvent(input$run_filter_qc, {
    req(input$readfile1)
    withProgress(message = "Running QC Filtering...", value = 0.1, {
      readfile1 <- input$readfile1$datapath
      readfile2 <- if (!is.null(input$readfile2)) input$readfile2$datapath else NULL
      run_filtering(readfile1, readfile2, input$sample_name, input$read_type, input$min_quality, input$min_length)
      incProgress(1)
    })
    append_log("✔ QC filtering completed.")
  }) 
  
  observeEvent(input$run_align, {
    req(input$ref_genome)
    withProgress(message = "Running Alignment...", value = 0.1, {
      readfile1 <- input$readfile1$datapath
      readfile2 <- if (!is.null(input$readfile2)) input$readfile2$datapath else NULL
      aligned <- run_alignment(readfile1, readfile2, input$ref_genome$datapath, input$sample_name, input$read_type)
      incProgress(1)
      append_log(paste("✔ Alignment completed. Output:", aligned))
    })
  })
  
  observeEvent(input$run_filter_bam, {
    withProgress(message = "Filtering BAM file...", value = 0.1, {
      bam_file <- paste0("../Data/Alignment_output/", input$sample_name, "_aligned.bam")
      filtered <- filter_bam_by_coverage(bam_file, input$sample_name, min_coverage = input$min_coverage)
      incProgress(1)
      append_log(paste("✔ BAM filtering completed. Output:", filtered))
    })
  })
  
  observeEvent(input$run_variant_call, {
    withProgress(message = "Calling Variants...", value = 0.1, {
      bam_file <- paste0("../Data/Alignment_output/", input$sample_name, "_aligned.bam")
      vcf_file <- run_variant_calling(bam_file, input$ref_genome$datapath)
      incProgress(1)
      append_log(paste("✔ Variant calling completed. Output:", vcf_file))
    })
  })
  
  output$log <- renderText({
    log_text()
  })
  
  output$align_summary <- renderText({
    file <- "alignment_summary.csv"
    if (file.exists(file)) paste(readLines(file, warn = FALSE), collapse = "\n") else "Alignment summary not available."
  })
  
  output$variant_output <- renderText({
    file <- "variants_called.vcf"
    if (file.exists(file)) paste(head(readLines(file, warn = FALSE), 20), collapse = "\n") else "VCF output not available."
  })
}

shinyApp(ui = ui, server = server)