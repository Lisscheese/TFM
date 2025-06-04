# ***************************************
# Shiny Interface for Genomic Pipeline
# ***************************************

options(shiny.maxRequestSize = 8 * 1024^10)

library(shiny)
library(shinythemes)

# Helper to capture startup messages and warnings
runAndLogStartup <- function(expr) {
  expr_code <- substitute(expr)
  output <- capture.output(
    withCallingHandlers(
      eval(expr_code, envir = parent.frame()),
      message = function(m) {
        cat("stating message", conditionMessage(m), "\n")
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        cat("warn", conditionMessage(w), "\n")
        invokeRestart("muffleWarning")
      }
    ),
    type = "output"
  )
  cat(paste(output, collapse = "\n"), "\n")
}

# Load libraries and modules with logging
runAndLogStartup(library(ShortRead))
runAndLogStartup(library(Rsamtools))
runAndLogStartup(library(GenomicRanges))
runAndLogStartup(source("Reduce_Sampling_module.R"))
runAndLogStartup(source("Quality_Control_module.R"))
runAndLogStartup(source("Filtering_module.R"))
runAndLogStartup(source("Alignment_module.R"))
runAndLogStartup(source("Alignment_Evaluation_module.R"))
runAndLogStartup(source("Variant_Calling_module.R"))
runAndLogStartup(source("Annotation_Module.R"))

# UI definition
ui <- fluidPage(
  tags$head(
    tags$style(HTML(
      "
      body { background-color: white; }
      .btn-primary { background-color: #2b8cbe; color: white; }
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
    ))
  ),
  theme = shinytheme("flatly"),
  titlePanel("Genomic Pipeline"),
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
      numericInput("subsample_n", "Block size N", value = 100000),
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
      div(id = "console", verbatimTextOutput("console", placeholder = TRUE))
    )
  )
)

# Server logic
server <- function(input, output, session) {
  rv <- reactiveValues(
    read1 = NULL,
    read2 = NULL,
    sorted_bam = NULL,
    vcf = NULL,
    annot = NULL,
    console = "Console log started",
    ref = NULL
  )
  
  # Helper to append console
  appendConsole <- function(msg) {
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    rv$console <- paste(rv$console, paste0("[", ts, "] ", msg), sep = "\n")
  }
  
  # Run and log expressions
  runAndLog <- function(expr) {
    expr_code <- substitute(expr)
    txt <- capture.output(
      withCallingHandlers(
        eval(expr_code, envir = parent.frame()),
        message = function(m) { cat("[MSG]", conditionMessage(m), "\n"); invokeRestart("muffleMessage") },
        warning = function(w) { cat("[WARN]", conditionMessage(w), "\n"); invokeRestart("muffleWarning") }
      ), type = "output"
    )
    sapply(strsplit(paste(txt, collapse = "\n"), "\n")[[1]], appendConsole)
  }
  
  # Observe file uploads
  observeEvent(input$readfile1, {
    req(input$readfile1)
    rv$read1 <- input$readfile1$datapath
    appendConsole(paste("Loaded Read 1:", input$readfile1$name))
  })
  observeEvent(input$readfile2, {
    rv$read2 <- input$readfile2$datapath
    appendConsole(paste("Loaded Read 2:", input$readfile2$name))
  })
  observeEvent(input$sorted_bam_upload, {
    rv$sorted_bam <- input$sorted_bam_upload$datapath
    appendConsole(paste("Loaded sorted BAM:", input$sorted_bam_upload$name))
  })
  observeEvent(input$ref_genome, {
    rv$ref <- input$ref_genome$datapath
    appendConsole(paste("Loaded reference FASTA:", input$ref_genome$name))
  })
  
  # Render console output
  output$console <- renderText({ rv$console })
  
  # Subsampling
  observeEvent(input$run_subsample, {
    req(rv$read1)
    appendConsole("** Subsampling started **")
    tryCatch({
      runAndLog({
        subs <- run_subsampling(readfile1 = rv$read1,
                                readfile2 = rv$read2,
                                sample_name = input$sample_name,
                                n = input$subsample_n)
        rv$read1 <- subs$fastq1
        rv$read2 <- subs$fastq2
      })
      appendConsole("Subsampling completed.")
    }, error = function(e) {
      appendConsole(paste("ERROR in Subsampling:", conditionMessage(e)))
    })
  })
  
  # Quality Control
  observeEvent(input$run_qc, {
    req(rv$read1, rv$ref)
    appendConsole("** QC started **")
    tryCatch({
      fmt <- ifelse(input$file_type == "fq", "fastq", input$file_type)
      runAndLog({
        run_qc(readfile1 = rv$read1,
               readfile2 = rv$read2,
               sample_name = input$sample_name,
               read_type = input$read_type,
               file_type = fmt)
      })
      appendConsole("QC completed.")
    }, error = function(e) {
      appendConsole(paste("ERROR in QC:", conditionMessage(e)))
    })
  })
  
  # FASTQ Filtering
  observeEvent(input$run_filter_qc, {
    req(rv$read1)
    appendConsole("** FASTQ filtering started **")
    tryCatch({
      if (input$do_filter_qc) {
        runAndLog({
          ff <- run_fastq_filtering(
            fastq1 = rv$read1,
            fastq2 = rv$read2,
            sample_name = input$sample_name,
            read_type = input$read_type,
            quality_threshold = input$min_quality
          )
          rv$read1 <- ff$fastq1
          rv$read2 <- ff$fastq2
        })
        appendConsole("FASTQ filtering completed.")
      } else {
        appendConsole("FASTQ filtering skipped.")
      }
    }, error = function(e) {
      appendConsole(paste("ERROR in FASTQ filtering:", conditionMessage(e)))
    })
  })
  
  # Alignment
  observeEvent(input$run_align, {
    req(rv$read1, rv$ref)
    appendConsole("** Alignment step started **")
    aligned_bam <- run_alignment(
      readfile1   = rv$read1,
      readfile2   = rv$read2,
      ref_genome  = rv$ref,
      sample_name = input$sample_name,
      read_type   = input$read_type
    )
    rv$sorted_bam <- aligned_bam
    appendConsole("Alignment completed.")
  })
  
  
  # Alignment Evaluation
  observeEvent(input$run_eval, {
    req(rv$sorted_bam)
    appendConsole("** Evaluation started **")
    tryCatch({
      runAndLog({
        res <- run_alignment_evaluation(bam_file = rv$sorted_bam,
                                        sample_name = input$sample_name)
        rv$sorted_bam <- res$sorted_bam
      })
      appendConsole("Evaluation completed.")
    }, error = function(e) {
      appendConsole(paste("ERROR in evaluation:", conditionMessage(e)))
    })
  })
  
  
  # BAM Filtering
  observeEvent(input$run_filter_bam, {
    req(rv$sorted_bam)
    appendConsole("** BAM filtering started **")
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
        appendConsole("BAM filtering skipped.")
      }
    }, error = function(e) {
      appendConsole(paste("ERROR in BAM filtering:", conditionMessage(e)))
    })
  })
  
  # Variant Calling
  observeEvent(input$run_variant_call, {
    req(rv$sorted_bam, rv$ref)
    appendConsole("** Variant calling started **")
    tryCatch({
      runAndLog({
        rv$vcf <- run_variant_calling(
          bam_file        = rv$sorted_bam,
          reference_fasta = rv$ref,
          sample_name     = input$sample_name
        )
      })
      appendConsole("Variant calling completed.")
    }, error = function(e) {
      appendConsole(paste("ERROR in variant calling:", conditionMessage(e)))
    })
  })
  
  # Annotation
  observeEvent(input$run_annotation, {
    req(rv$vcf)
    appendConsole("** Annotation started **")
    tryCatch({
      runAndLog({
        rv$annot <- run_snpEff(vcf_file = rv$vcf, sample_name = input$sample_name)
      })
      appendConsole("Annotation completed.")
    }, error = function(e) {
      appendConsole(paste("ERROR in annotation:", conditionMessage(e)))
    })
  })
  
  # Exit App
  observeEvent(input$exit, { stopApp() })
}

# Launch the application
shinyApp(ui, server)