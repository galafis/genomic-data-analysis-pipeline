# dashboard_genomic.R - Template R Shiny para análise de variantes VCF
# Autor: Gabriel
# Descrição: Dashboard básico com upload de VCF, tabela de variantes,
# gráficos de frequência alélica e heatmap de qualidade (Phred),
# usando módulos simples do Shiny. Este é um template inicial e pode ser
# expandido conforme necessário.

# =====================
# Instruções de uso
# =====================
# 1) Instale dependências em R:
#    install.packages(c("shiny","shinydashboard","DT","ggplot2","reshape2"))
#    # Para leitura de VCF via vroom/readr, este template aceita VCF sem BGZF/Tabix.
#    # Para VCFs grandes/compactados, considere usar VariantAnnotation (Bioconductor).
#    # Bioc: if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#    # BiocManager::install("VariantAnnotation")
# 2) Salve este arquivo em src/visualization/interactive/shiny_app/dashboard_genomic.R
# 3) Rode localmente: shiny::runApp("src/visualization/interactive/shiny_app")
# 4) Faça upload de um .vcf (ou .vcf.gz já descompactado) para visualizar.

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(DT)
  library(ggplot2)
  library(reshape2)
})

# =====================
# Funções auxiliares
# =====================
# Leitura simples de VCF: parse de linhas não iniciadas por '#',
# separando colunas padrão: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT samples...
# Observação: Para VCF completos com genótipos, recomendas-se VariantAnnotation.

read_vcf_simple <- function(path) {
  # Lê linhas e filtra cabeçalho
  lines <- readLines(path, warn = FALSE)
  data_lines <- lines[!startsWith(lines, "#")]
  if (length(data_lines) == 0) return(data.frame())
  # Divide por tabulação
  split_rows <- strsplit(data_lines, "\t", fixed = TRUE)
  # Garante no mínimo 8 colunas (VCF mínimo)
  max_cols <- max(vapply(split_rows, length, integer(1)))
  mat <- do.call(rbind, lapply(split_rows, function(x) {
    length(x) <- max_cols
    x
  }))
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  # Nomeia colunas básicas
  base_cols <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  if (ncol(df) >= length(base_cols)) {
    colnames(df)[1:length(base_cols)] <- base_cols
  }
  # Tipos numéricos
  if ("POS" %in% names(df)) df$POS <- suppressWarnings(as.integer(df$POS))
  if ("QUAL" %in% names(df)) df$QUAL <- suppressWarnings(as.numeric(df$QUAL))
  df
}

# Extrai frequência alélica (AF) do campo INFO quando disponível (ex.: AF=0.12)
extract_af <- function(info) {
  if (is.na(info) || info == "") return(NA_real_)
  m <- regmatches(info, regexpr("AF=([0-9.eE-]+)", info))
  if (length(m) == 0 || m == "") return(NA_real_)
  val <- sub("AF=", "", m)
  suppressWarnings(as.numeric(val))
}

# =====================
# Módulos Shiny
# =====================
# Módulo de upload e parsing de VCF
vcfUploadUI <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("vcf_file"), "Upload VCF", accept = c(".vcf","text/vcf")),
    checkboxInput(ns("has_header"), "Arquivo contém cabeçalho (#)", TRUE),
    helpText("Suporta VCF texto simples. Para VCFs grandes/compactados, considere VariantAnnotation.")
  )
}

vcfUploadServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    vcf_data <- reactive({
      req(input$vcf_file)
      df <- read_vcf_simple(input$vcf_file$datapath)
      if (nrow(df) == 0) return(NULL)
      # Extrai AF quando possível
      if ("INFO" %in% names(df)) {
        df$AF <- vapply(df$INFO, extract_af, numeric(1))
      }
      df
    })
    return(vcf_data)
  })
}

# Módulo de tabela de variantes
variantsTableUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(width = 12, title = "Tabela de Variantes", status = "primary", solidHeader = TRUE,
          DTOutput(ns("tbl")))
    )
  )
}

variantsTableServer <- function(id, data_r) {
  moduleServer(id, function(input, output, session) {
    output$tbl <- renderDT({
      df <- data_r()
      validate(need(!is.null(df) && nrow(df) > 0, "Carregue um VCF válido."))
      datatable(df[, intersect(names(df), c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","AF"))],
                rownames = FALSE, filter = "top", options = list(pageLength = 10))
    })
  })
}

# Módulo de gráficos de frequência alélica
alleleFreqUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(width = 6, title = "Distribuição de AF", status = "success", solidHeader = TRUE,
          plotOutput(ns("af_hist"), height = 300)),
      box(width = 6, title = "AF vs QUAL", status = "success", solidHeader = TRUE,
          plotOutput(ns("af_scatter"), height = 300))
    )
  )
}

alleleFreqServer <- function(id, data_r) {
  moduleServer(id, function(input, output, session) {
    output$af_hist <- renderPlot({
      df <- data_r()
      validate(need(!is.null(df) && "AF" %in% names(df), "Campo AF não encontrado no INFO."))
      ggplot(df, aes(x = AF)) +
        geom_histogram(color = "white", fill = "#3182bd", bins = 30, na.rm = TRUE) +
        labs(x = "Frequência Alélica (AF)", y = "Contagem") +
        theme_minimal()
    })
    output$af_scatter <- renderPlot({
      df <- data_r()
      validate(need(!is.null(df) && all(c("AF","QUAL") %in% names(df)), "AF e/ou QUAL indisponíveis."))
      ggplot(df, aes(x = AF, y = QUAL)) +
        geom_point(alpha = 0.6, color = "#2ca02c", na.rm = TRUE) +
        labs(x = "AF", y = "QUAL (Phred)") +
        theme_minimal()
    })
  })
}

# Módulo de heatmap de qualidade (SOLID)
# Interpretação: criamos um heatmap simples de QUAL por cromossomo vs janelas de posição.
# "solid" aqui refere-se a estilo sólido (geom_tile com cores sólidas).
qualityHeatmapUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(width = 12, title = "Heatmap de Qualidade (Sólido)", status = "warning", solidHeader = TRUE,
          sliderInput(ns("bins"), "Número de janelas por cromossomo:", min = 10, max = 200, value = 50, step = 10),
          plotOutput(ns("heatmap"), height = 380))
    )
  )
}

qualityHeatmapServer <- function(id, data_r) {
  moduleServer(id, function(input, output, session) {
    output$heatmap <- renderPlot({
      df <- data_r()
      validate(need(!is.null(df) && all(c("CHROM","POS","QUAL") %in% names(df)), "Campos CHROM, POS e QUAL são necessários."))
      # Remove NAs
      df <- df[!is.na(df$QUAL) & !is.na(df$POS) & !is.na(df$CHROM), ]
      validate(need(nrow(df) > 0, "Sem dados suficientes para heatmap."))
      # Binning por cromossomo
      bin_df <- do.call(rbind, lapply(split(df, df$CHROM), function(sub) {
        if (nrow(sub) < 1) return(NULL)
        rng <- range(sub$POS, na.rm = TRUE)
        if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) return(NULL)
        cuts <- cut(sub$POS, breaks = input$bins, include.lowest = TRUE)
        agg <- aggregate(QUAL ~ cuts, data = transform(sub, cuts = cuts), FUN = median, na.rm = TRUE)
        data.frame(CHROM = unique(sub$CHROM)[1], BIN = agg$cuts, QUAL = agg$QUAL)
      }))
      validate(need(!is.null(bin_df) && nrow(bin_df) > 0, "Não foi possível agregar por janelas."))
      # Ordena cromossomos de forma natural se possível
      chrom_levels <- unique(bin_df$CHROM)
      bin_df$CHROM <- factor(bin_df$CHROM, levels = chrom_levels)
      ggplot(bin_df, aes(x = BIN, y = CHROM, fill = QUAL)) +
        geom_tile(color = NA) +
        scale_fill_gradient(low = "#fee5d9", high = "#a50f15", na.value = "grey85") +
        labs(x = "Janelas de posição", y = "Cromossomo", fill = "QUAL (mediana)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
  })
}

# =====================
# UI e Server principal
# =====================
ui <- dashboardPage(
  dashboardHeader(title = "Genomic Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload", tabName = "upload", icon = icon("file-arrow-up")),
      menuItem("Variantes", tabName = "variants", icon = icon("table")),
      menuItem("Frequência Alélica", tabName = "af", icon = icon("chart-area")),
      menuItem("Heatmap Qualidade", tabName = "heatmap", icon = icon("th"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "upload",
              fluidRow(box(width = 12, vcfUploadUI("u1")))) ,
      tabItem(tabName = "variants", variantsTableUI("t1")),
      tabItem(tabName = "af", alleleFreqUI("a1")),
      tabItem(tabName = "heatmap", qualityHeatmapUI("h1"))
    )
  )
)

server <- function(input, output, session) {
  vcf_data <- vcfUploadServer("u1")
  variantsTableServer("t1", vcf_data)
  alleleFreqServer("a1", vcf_data)
  qualityHeatmapServer("h1", vcf_data)
}

# Executa app
shinyApp(ui, server)
