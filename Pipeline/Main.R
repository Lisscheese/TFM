# ***************************************
# Shiny Interface for Genomic Pipeline
# ***************************************

options(shiny.maxRequestSize = 8 * 1024^10)

library(shiny)
library(shinythemes)

# Load modules 
source("Reduce_Sampling_module.R")
source("Quality_Control_module.R")
source("Filtering_module.R")
source("Alignment_module.R")
source("Alignment_Evaluation_module.R")
source("Variant_Calling_module.R")
source("Annotation_Module.R")

# UI definition
ui <- fluidPage(
  # Custom CSS
  tags$head(tags$style(HTML(
    "
    body { background-color: white; }
    .btn-primary { background-color: #2b8cbe; color: white; }
    .btn { background-color: #2b8cbe; color: white; }
    .btn-danger { background-color: #d73027; color: white; }
    .form-control, textarea { background-color: #e6f5e6; }
    .shiny-input-container { margin-bottom: 15px; }
    h2, h3, h4 { color: #2b8cbe; }
    #console { 
              background: #000; 
              color: #0f0; 
              padding: 10px; 
              font-family: monospace; 
              height: 1200px; 
              overflow-y: auto; 
              white-space: pre-wrap; 
    }
    "
  ))),
  theme = shinytheme("flatly"),
  titlePanel("Genomic Pipeline - Console Interface"),
  sidebarLayout(
    sidebarPanel(
      fileInput("readfile1", "Read File 1 (R1)", accept = c(".fq", ".fastq")),
      fileInput("readfile2", "Read File 2 (R2, optional)", accept = c(".fq", ".fastq")),
      fileInput("sorted_bam_upload", "Sorted BAM (optional)", accept = ".bam"),
      fileInput("ref_genome", "Reference Genome (FASTA)", accept = ".fa"),
      textInput("sample_name", "Sample Name", value = "sample"),
      selectInput("read_type", "Read Type", choices = c("single", "paired")),
      selectInput("file_type", "File Format", choices = c("fastq", "fq")),
      hr(),
      checkboxInput("do_subsample", "Enable Subsampling", value = FALSE),
      numericInput("subsample_n", "Subsample reads N", value = 100000),
      actionButton("run_subsample", "Run Subsampling", class = "btn btn-primary"),
      hr(),
      actionButton("run_qc", "Run Quality Control", class = "btn btn-primary"),
      hr(),
      checkboxInput("do_filter_qc", "Enable FASTQ Filtering", value = FALSE),
      numericInput("min_quality", "Min Average Quality", value = 30),
      actionButton("run_filter_qc", "Run FASTQ Filtering", class = "btn btn-primary"),
      hr(),
      actionButton("run_align", "Run Alignment", class = "btn btn-primary"),
      actionButton("run_eval", "Evaluate Alignment", class = "btn btn-primary"),
      hr(),
      checkboxInput("do_filter_bam", "Enable BAM Filtering", value = FALSE),
      numericInput("min_coverage", "Min Coverage", value = 10),
      actionButton("run_filter_bam", "Run BAM Filtering", class = "btn btn-primary"),
      hr(),
      actionButton("run_variant_call", "Run Variant Calling", class = "btn btn-primary"),
      actionButton("run_annotation", "Run Annotation", class = "btn btn-primary"),
      hr(),
      actionButton("exit", "Exit App", class = "btn btn-danger")
    ),
    mainPanel(
      tags$div(id = "console", verbatimTextOutput("console", placeholder = TRUE))
    )
  )
)

# Server logic
server <- function(input, output, session) {
  rv <- reactiveValues(
    read1      = NULL,
    read2      = NULL,
    sorted_bam = NULL,
    vcf        = NULL,
    annot      = NULL,
    console    = "Console log started"
  )
  
  # helper to append console
  appendConsole <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    rv$console <- paste(rv$console, paste0("[", timestamp, "] ", msg), sep = "\n")
  }
  
  # helper to run an expression once, capture stdout + message()
  runAndLog <- function(expr) {
    expr_code <- substitute(expr)

    txt <- capture.output(
      withCallingHandlers(
        {
          eval(expr_code, envir = parent.frame())
        },
        message = function(m) {

          cat(conditionMessage(m), "\n")
          invokeRestart("muffleMessage")
        }
      ),
      type = "output"
    )
    lines <- unlist(strsplit(paste(txt, collapse = "\n"), "\n"))
    lapply(lines, appendConsole)
  }
  
  # monitor uploads
  observeEvent(input$readfile1, {
    rv$read1 <- input$readfile1$datapath
    appendConsole(paste("Loaded Read 1:", basename(rv$read1)))
  })
  observeEvent(input$readfile2, {
    rv$read2 <- input$readfile2$datapath
    appendConsole(paste("Loaded Read 2:", basename(rv$read2)))
  })
  observeEvent(input$sorted_bam_upload, {
    rv$sorted_bam <- input$sorted_bam_upload$datapath
    appendConsole(paste("Loaded sorted BAM:", basename(rv$sorted_bam)))
  })
  observeEvent(input$ref_genome, {
    appendConsole(paste("Loaded reference FASTA:", basename(input$ref_genome$datapath)))
  })
  
  output$console <- renderText({ rv$console })
  
  # Subsampling
  observeEvent(input$run_subsample, {
    req(rv$read1)
    appendConsole("** Subsampling step started **")
    tryCatch({
      runAndLog({ subs <- run_subsampling(
        rv$read1, rv$read2,
        sample_name = input$sample_name,
        n           = input$subsample_n
      )
      rv$read1 <- subs$fastq1; rv$read2 <- subs$fastq2
      })
      appendConsole(paste("Subsampling completed:", basename(rv$read1), basename(rv$read2)))
    }, error = function(e) appendConsole(paste("ERROR in subsampling:", e$message)))
  })
  
  # Quality Control
  observeEvent(input$run_qc, {
    req(rv$read1, input$ref_genome)
    appendConsole("** QC step started **")
    tryCatch({
      fmt <- ifelse(input$file_type == "fq", "fastq", input$file_type)
      runAndLog(run_qc(
        readfile1   = rv$read1,
        readfile2   = rv$read2,
        sample_name = input$sample_name,
        read_type   = input$read_type,
        file_type   = fmt
      ))
      appendConsole("QC completed successfully.")
    }, error = function(e) appendConsole(paste("ERROR in QC:", e$message)))
  })
  
  # FASTQ Filtering
  observeEvent(input$run_filter_qc, {
    req(rv$read1)
    appendConsole("** FASTQ filtering step started **")
    tryCatch({
      if (input$do_filter_qc) {
        runAndLog({
          ff <- run_fastq_filtering(
            fastq1            = rv$read1,
            fastq2            = rv$read2,
            sample_name       = input$sample_name,
            read_type         = input$read_type,
            quality_threshold = input$min_quality
          )
          rv$read1 <- ff$fastq1
          rv$read2 <- ff$fastq2
        })
        appendConsole("FASTQ filtering completed.")
      } else {
        appendConsole("FASTQ filtering skipped (checkbox disabled)")
      }
    }, error = function(e) appendConsole(paste("ERROR in FASTQ filtering:", e$message)))
  })
  
  
  # Alignment
  observeEvent(input$run_align, {
    req(rv$read1, input$ref_genome)
    appendConsole("** Alignment step started **")
    tryCatch({
      runAndLog({
        aligned <- run_alignment(
          readfile1   = rv$read1,
          readfile2   = rv$read2,
          ref_genome  = input$ref_genome$datapath,
          sample_name = input$sample_name,
          read_type   = input$read_type
        )
        rv$sorted_bam <- aligned
      })
      appendConsole("Alignment completed.")
    }, error = function(e) appendConsole(paste("ERROR in alignment:", e$message)))
  })
  
  # Evaluate Alignment
  observeEvent(input$run_eval, {
    req(rv$sorted_bam)
    appendConsole("** Evaluation step started **")
    tryCatch({
      runAndLog({ rs <- run_alignment_evaluation(rv$sorted_bam, sample_name = input$sample_name)
      rv$sorted_bam <- rs$sorted_bam
      })
      appendConsole("Evaluation completed.")
    }, error = function(e) appendConsole(paste("ERROR in evaluation:", e$message)))
  })
  
  # BAM Filtering
  observeEvent(input$run_filter_bam, {
    req(rv$sorted_bam)
    appendConsole("** BAM filtering step started **")
    tryCatch({
      if (input$do_filter_bam) {
        runAndLog({
          rv$sorted_bam <- run_bam_filtering(
            bam_file     = rv$sorted_bam,
            sample_name  = input$sample_name,
            min_coverage = input$min_coverage
          )
        })
        appendConsole("BAM filtering completed.")
      } else {
        appendConsole("BAM filtering skipped (checkbox disabled)")
      }
    }, error = function(e) appendConsole(paste("ERROR in BAM filtering:", e$message)))
  })
  
  
  # Variant Calling
  observeEvent(input$run_variant_call, {
    req(rv$sorted_bam, input$ref_genome)
    appendConsole("** Variant calling step started **")
    tryCatch({
      runAndLog({ rv$vcf <- run_variant_calling(
        bam_file        = rv$sorted_bam,
        reference_fasta = input$ref_genome$datapath,
        sample_name     = input$sample_name
      )
      })
      appendConsole("Variant calling completed.")
    }, error = function(e) appendConsole(paste("ERROR in variant calling:", e$message)))
  })
  
  # Annotation
  observeEvent(input$run_annotation, {
    req(rv$vcf)
    appendConsole("** Annotation step started **")
    tryCatch({
      runAndLog({ ann <- run_annotation(
        vcf_file     = rv$vcf,
        genome_build = "hg38",
        sample_name  = input$sample_name
      )
      rv$annot <- ann
      })
      appendConsole("Annotation completed.")
    }, error = function(e) appendConsole(paste("ERROR in annotation:", e$message)))
  })
  
  # Exit
  observeEvent(input$exit, { stopApp() })
}

# Launch the app
shinyApp(ui = ui, server = server)
