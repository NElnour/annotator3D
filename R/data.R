#' Matched Motifs that Regulate H5-K9 Methylation
#'
#' The data set is obtained by scanning murine PPMs provided by Ngo \emph{et al.} with FIMO against mm10 chromsome 18. The top 50 matches by \emph{q}=value were selected and organized into a dataframe.
#'
#' @format A dataframe with 50 rows and 10 variables:
#' \describe{
#'   \item{motif_id}{chromosome 18 fragments start coordinates}
#'   \item{motif_alt_id}{chromosome 18 fragments end coordinates}
#'   \item{sequence_name}{chromsome name at which the matches were found}
#'   \item{start}{start coordinates of the matched motif on chromsome 18}
#'   \item{stop}{end coordinate of the matched motif on chromosome 18}
#'   \item{strand}{Either "+" or "-" denoting the strand orientation of the match}
#'   \item{score}{Maximum likelihood score}
#'   \item{p.value}{The probability that the motif matchs to the chromsome region bounded by start and stop coordinates.}
#'   \item{q.value}{The \emph{p}-value corrected by Benjamini and Hochsberg false discvery rate for multiple testing}
#'   \item{matched_sequence}{Consensus sequence of matched motif}
#' }
#'
#' @references {
#' Ngo V, Chen Z, Zhang K, Whitaker J, Wang M, Wang W. (2019) Epigenomic analysis reveals DNA motifs regulating histone modifications in human and mouse. \emph{PNAS} \strong{116}(9): 3668--3677.
#' }
#' @source \url{https://www.pnas.org/content/116/9/3668}
"matchedMotifs"

#' Chromsome interactions for 40kb segments of mES cells' chromsome 18
#'
#' A subset 1000-by-1000 symmetric matrix of interchromsome interactions murine chromsome 18 embryonic stem cells. Murine genome was sequenced on Illumina HiSeq 2000 after HindIII and NCoI digestion, achiving resolution of 40 kb.
#' @format A matrix with 1000 rows and 1000 columns:
#' \describe{
#'   \item{row}{chromosome 18 fragments start coordinates}
#'   \item{column}{chromosome 18 fragments end coordinates}
#' }
#' @references {
#' Dixon JR, Selvaraj S, Yue F, Kim A \emph{et al.} (2012) Topological domains in mammalian genomes identified by analysis of chromatin interactions. \emph{Nature} \strong{485}(7398):376--80. PMID: \href{https://www.ncbi.nlm.nih.gov/pubmed/22495300}{22495300}
#' }
#' @source \url{http://sysbio.rnet.missouri.edu/3dgenome/GSDB/details.php?id=TE1402WS}
"mESChromInteractions"
