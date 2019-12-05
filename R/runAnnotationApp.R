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
#' @export
#'
#' @examples
#' \dontrun{
#'   runAnnotator()
#' }
#'
#' @importFrom shiny runApp
#'
runAnnotator <- function() {
  path_tp_app <- system.file("shinyApp",
                        package = "annotator3D")

  shiny::runApp(appDir=path_tp_app, display.mode = "normal")
}

# [END]
