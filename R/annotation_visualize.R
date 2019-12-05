#' Title
#'
#' @param motif_bed A string of the full path to the BED file
#' containing the classified motifs
#' @param gen A named character containing the chromosome name
#' @param thresh An optional numeric value denoting the threshold of
#' significance. Motifs with membership scores greater than or equal to thresh
#' will be returned.
#'
#' @return A Formal Class AnnotationTrack for enriched motif with GRanges
#' corresponding to motif's genomic ranges,
#' @export
#'
#' @examples
#' \dontrun{
#' data("A549ChromLoops")
#' data("matchedMotifs")
#' classify(matchedMotifs, A549ChromLoops)
#' gen <- c("mm10")
#' names(gen) <- "chr18"
#' antrack <- classifiedAnnotationTrack("chromLoops.cl"./BEDs/H3K9me3.bed", gen)
#' }
#'
#'@import Gviz
#'
classifiedAnnotationTrack <- function(motif_bed, gen, thresh = 0) {
  if (missing(motif_bed) | missing(gen)) {
    stop("Missing mapped motifs file path")
  }

  cols <- c("chromosome",
            "start",
            "end",
            "symbol",
            "score")

  motifs <- utils::read.delim(
    motif_bed,
    stringsAsFactors = FALSE,
    col.names = cols,
    header = FALSE
  )

  motifs <- filter_motifs(motifs, thresh = thresh)

  chrom <- unique(as.character(motifs$chromosome))

  antrack <- Gviz::GeneRegionTrack(
    motifs,
    genome = gen,
    chromosome = chrom,
    name = "Mapped Features",
    transcript = "symbol",
    transcriptAnnotation = "symbol"
  )

  return(antrack)

}


#' Filter BED Dataframe by Threshold
#'
#' Filter classified motifs based on a lower bound threshold on their membership
#' score.
#'
#' @param motifs A dataframe with column names corresponding to BED file format.
#' @param thresh An optional numeric value that serves as lower bound for scores
#' in the dataframe
#'
#' @return A dataframe with column names corresponding to BED format and scores
#' greater than or equal to lower bound threshold.

filter_motifs <- function(motifs, thresh = 0) {
  tmp <- motifs[motifs$score >= thresh, ]
  return(tmp)

}

#' Plot a Layered Gviz Figure
#'
#' Plot the classified motifs, along with the ideogram of the chromosome of
#' interest and and genome axis.
#'
#' @param antrack A classified annotation track obtained with
#' \code{classifiedAnnotationTrack}.
#' @param gen A named character containing genome of interest
#' @param from An optional integer value of start of viewing window.
#' @param to An optional integer value of end of viewing window.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("A549ChromLoops")
#' data("matchedMotifs")
#' classify(matchedMotifs, A549ChromLoops)
#' gen <- c("mm10")
#' names(gen) <- "chr18"
#' antrack <- classifiedAnnotationTrack("chromLoops.cl"./BEDs/H3K9me3.bed", gen)
#' plot_layered_gviz(antrack, gen)
#' }
#'
#'@import Gviz
#'
plot_layered_gviz <- function(antrack, gen, from=1, to=1000000) {
  ideogram <- Gviz::IdeogramTrack(genome = gen, chromosome = names(gen))
  genometrack <- Gviz::GenomeAxisTrack()
  Gviz::plotTracks(list(ideogram, genometrack, antrack),
                   from=from,
                   to=to,
                   col=NULL)
}

# [END]
