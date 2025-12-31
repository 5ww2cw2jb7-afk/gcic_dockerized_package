# app.R (GCIC single-gene pipeline; single-page view; Log top; show GCIC island FASTA; overwrite summary each run)
library(shiny)
library(readr)
library(DT)
library(fs)

ui <- fluidPage(
  titlePanel("GCIC single-gene pipeline (Rice / Arabidopsis)"),
  sidebarLayout(
    sidebarPanel(
      selectInput("species", "Species",
                  choices = c("Rice (Oryza sativa)" = "rice",
                              "Arabidopsis thaliana" = "arabidopsis"),
                  selected = "rice"),
      textInput("gene_id", "Gene ID", value = "Os03g0215200"),

      selectInput("strand_mode", "strand-mode",
                  choices = c("fwd", "rev", "total"), selected = "fwd"),

      numericInput("window_bp", "Window (bp)", value = 200, min = 10, step = 1),
      numericInput("step_bp",   "Step (bp)",   value = 50,  min = 1,  step = 1),

      textInput("python_bin", "python 実行コマンド", value = "python3"),
      actionButton("run_btn", "RUN", class = "btn-primary"),

      tags$hr(),
      verbatimTextOutput("input_paths"),
      tags$small("※ アプリフォルダに gene_regions.5prime.fa / gene_regions.bed（Arabidopsis は ATgene_regions.*）がある前提")
    ),
    mainPanel(
      h3("Log"),
      tags$pre(textOutput("log_out")),
      tags$hr(),

      h3("Summary"),
      DTOutput("summary_tbl"),
      tags$hr(),

      h3("GCIC islands (FASTA)"),
      tags$pre(textOutput("fa_out"))
    )
  )
)

server <- function(input, output, session) {

  rv <- reactiveValues(
    log = "",
    summary_df = NULL,
    fa_text = ""
  )

  append_log <- function(x) rv$log <- paste0(rv$log, x, "\n")

  # 起動時の作業ディレクトリをアプリルートとして固定（cdしてrunAppする前提に対応）
  app_root <- normalizePath(getwd())

  master_fasta <- reactive({
    if (input$species == "arabidopsis") {
      path(app_root, "ATgene_regions.5prime.fa")
    } else {
      path(app_root, "gene_regions.5prime.fa")
    }
  })

  bed_path <- reactive({
    if (input$species == "arabidopsis") {
      path(app_root, "ATgene_regions.bed")
    } else {
      path(app_root, "gene_regions.bed")
    }
  })

  output$input_paths <- renderText({
    paste0(
      "app_root: ", app_root, "\n",
      "master_fasta: ", master_fasta(), "\n",
      "bed: ", bed_path(), "\n"
    )
  })

  output$log_out <- renderText(rv$log)
  output$fa_out  <- renderText(rv$fa_text)

  observeEvent(input$species, {
    if (input$species == "arabidopsis") {
      updateTextInput(session, "gene_id", value = "AT4G18960")
      updateNumericInput(session, "window_bp", value = 100)
      updateNumericInput(session, "step_bp", value = 25)
    } else {
      updateTextInput(session, "gene_id", value = "Os03g0215200")
      updateNumericInput(session, "window_bp", value = 200)
      updateNumericInput(session, "step_bp", value = 50)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$run_btn, {

    # 実行のたびに UI 表示をクリア（「溜めない」）
    rv$log <- ""
    rv$summary_df <- NULL
    rv$fa_text <- ""

    if (!file_exists(master_fasta())) {
      append_log(paste0("[ERROR] master-fasta not found: ", master_fasta()))
      return()
    }
    if (!file_exists(bed_path())) {
      append_log(paste0("[ERROR] bed not found: ", bed_path()))
      return()
    }

    out_root <- path(app_root, "results")
    dir_create(out_root)

    append_log(paste0("[INFO] app_root: ", app_root))

    args <- c(
      "run_single_gene.py",
      "--species", input$species,
      "--gene-id", input$gene_id,
      "--master-fasta", master_fasta(),
      "--bed", bed_path(),
      "--strand-mode", input$strand_mode,
      "--window", as.character(input$window_bp),
      "--step", as.character(input$step_bp)
    )
    append_log(paste0("[INFO] CMD: ", input$python_bin, " ", paste(shQuote(args), collapse = " ")))

    out <- tryCatch({
      system2(input$python_bin, args = args, stdout = TRUE, stderr = TRUE)
    }, error = function(e) {
      paste0("[R ERROR] ", conditionMessage(e))
    })

    if (length(out) > 0) append_log(paste(out, collapse = "\n"))

    # summary は「毎回上書き」される前提で読み込み
    summary_path <- path(app_root, "results", "all_genes_gcic_summary.tsv")
    if (file_exists(summary_path)) {
      df <- tryCatch({
        readr::read_tsv(summary_path, show_col_types = FALSE)
      }, error = function(e) NULL)
      rv$summary_df <- df
    } else {
      append_log(paste0("[ERROR] summary missing: ", summary_path))
      rv$summary_df <- NULL
    }

    # GCIC islands FASTA 表示（入力FASTAではなく island 配列）
    island_fa <- path(app_root, "results", "per_gene", input$gene_id, "GCIC.multi.island.fa")
    if (file_exists(island_fa)) {
      rv$fa_text <- paste(readLines(island_fa, warn = FALSE), collapse = "\n")
      append_log(paste0("[INFO] island FASTA: ", island_fa))
    } else {
      rv$fa_text <- ""
      append_log(paste0("[WARN] GCIC island FASTA not found: ", island_fa, " (no islands or generation failed)."))
    }
  })

  output$summary_tbl <- renderDT({
    if (is.null(rv$summary_df)) return(NULL)
    datatable(rv$summary_df, options = list(pageLength = 10, scrollX = TRUE))
  })

}

shinyApp(ui, server)
