#' Run annotator3D as a Shiny app
#'
#' runAnnotationApp runs annotator3D as a Shiny app.
#' The user needs to upload FIMO matched motif TSV file and a TSV file
#' containing predicted chromatin loops. The user can then customize parameters
#' of visualization such as lower bound for loop membership score, and Gviz
#' plotting options.
#'
#' Example TSV files are included in \code{extdata} to run annotator3D.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#'   runAnnotator()
#' }
#'
runAnnotator <- function() {
  annDir <- system.file("shinyApp",
              package= "annotator3D")

  shiny::runApp(annDir, display.mode = "normal")
}

# [END]
