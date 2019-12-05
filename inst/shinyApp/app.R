# The following script links to the runAnnotator function in
# runAnnotationApp.R. It handles presentation of Gviz plots through
# a Shiny app, in the same way calling the individual functions of
# annotator3D would (see attached vignette).
#
# To build this app, I followed RStudio's tutorial
# (found at https://shiny.rstudio.com/tutorial/).
#
# Author: Nada Elnour

library(annotator3D)
library(shiny)

gui <- fluidPage(
  titlePanel(h1(
    "Classify Motifs into Chromatin Loops", align = "center"
  )),

  sidebarPanel(
    h4("Main Control"),
    fileInput("motifs", h5("Upload matched motifs:")),
    fileInput("chromLoops", h5("Upload predicted loops:")),
    sliderInput(
      "threshold",
      "Score Threshold:",
      min = 0,
      max = 1,
      value = 0
    ),
    textInput("genome", h5("Genome:"), value = "Enter genome"),
    textInput("chr", h5("Choromosome no.:"), value = "Enter chromosome number"),
    #tags$hr(),
    h4("Gviz paramaters"),
    sliderInput(
      "from",
      "From:",
      min = 1,
      max = 248999999,
      value = 0
    ),
    sliderInput(
      "to",
      "To:",
      min = 2,
      max = 249000000,
      value = 0
    ),
    actionButton("submit", "Submit")
  ),

  mainPanel(plotOutput("GvizLayers", height="700px"))
)

server <- function(input, output) {
  plot_gviz <- eventReactive(input$submit, {
    plot(runif(1000))
    req(input$motifs)
    req(input$chromLoops)
    req(input$genome)
    req(input$chr)

    matched_motifs <- read.delim(
      input$motifs$datapath,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )

    matched_motifs <-
      matched_motifs[complete.cases(matched_motifs), ]

    check1 <- assertthat::are_equal(
      colnames(matched_motifs),
      c(
        "motif_id",
        "motif_alt_id",
        "sequence_name",
        "start",
        "stop",
        "strand",
        "score",
        "p.value",
        "q.value",
        "matched_sequence"
      )
    )

    if (!check1) {
      stop(
        "Motif file must have the following columns in order: motif_id,
        motif_alt_id, sequence_name, start, stop, strand, score, p.value,
        q.value, matched_sequence "
      )
    }

    loops <- read.delim(
      input$chromLoops$datapath,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )

    loops <- loops[complete.cases(loops), ]

    check2 <-
      assertthat::are_equal(colnames(loops), c("chr", "start", "stop"))

    if (!check2) {
      stop("Chromatin loops TSV must have the following columns in order: chr, start, stop")
    }

    motifs.cl <- classify(matched_motifs, loops)

    colnames(motifs.cl) <- c("chromosome",
                             "start",
                             "end",
                             "symbol",
                             "score")

    motifs.cl <- filter_motifs(motifs.cl, thresh = input$threshold)

    chrom <- unique(as.character(motifs.cl$chromosome))
    gen <- c("mm10")
    names(gen) <- chrom

    antrack <- Gviz::GeneRegionTrack(
      motifs.cl,
      genome = gen,
      chromosome = chrom,
      name = "Mapped Features",
      transcript = "symbol",
      transcriptAnnotation = "symbol"
    )

    plot_layered_gviz(antrack, gen, from = input$from, to = input$to)

  })

  output$GvizLayers <- renderPlot({
    plot_gviz()
  })
}

shinyApp (gui, server)

# [END]
