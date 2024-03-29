% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\name{classify}
\alias{classify}
\title{Motif Soft Classification into Topologically-associating Domains}
\usage{
classify(motifs, chromLoops)
}
\arguments{
\item{motifs}{A dataframe of matched motifs with at least the following
attributes:
\itemize{
\item \emph{motif_id}: a string array containing the motifs' names
\item \emph{start}: an integer array cotaining the starting coordinates of
matched motif
\item \emph{stop}: an integer array containing the end coordinate of matched
motif
}}

\item{chromLoops}{A dataframe of hierarchically clustered chromatin loops
with at least the
following attributes:
#' \itemize{
\item \emph{start}: an integer array containing the starting coordinates
of the predicted chromatin loops
\item \emph{stop}: an integer array containing the end coordinates of
the predicted chromatin loops
}}
}
\value{
A dataframe of classified of motifs, their membership scores, and
locations
}
\description{
Classifies motifs by location into hierarchically clustered chromatin loops.
It assumes that each base position has an equal chance of occurring in a
chromatin loop.
}
\examples{
data("A549ChromLoops")
data("matchedMotifs")
classify(matchedMotifs, A549ChromLoops)

}
