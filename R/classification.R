#' Motif Soft Classification into Topologically-associating Domains
#'
#' Classifies motifs by location into hierarchically clustered Hi-C domains.
#' It assumes that each base position has an equal chance of occurring in a
#' topologically-associating domain (TAD)
#'
#' @param motifs A dataframe of matched motifs with at least the following
#' attributes:
#' \itemize{
#' \item \emph{motif_id}: a string array containing the motifs' names
#' \item \emph{start}: an integer array cotaining the starting coordinates of
#' matched motif
#' \item \emph{stop}: an integer array containing the end coordinate of matched
#' motif
#' }
#' @param tads A dataframe of hierarchically clustered TADs with at least the
#' following attributes:
#' #' \itemize{
#' \item \emph{segment1}: an integer array containing the starting coordinates
#' of the predicted TADs
#' \item \emph{segment2}: an integer array containing the end coordinates of
#' the predicted TADs
#' \item \emph{hcluster}: a "-" separated string array containing the
#' hierarchical classification of predicted TADs.
#' }
#' @param offset A numeric denoting the resolution of Hi-C data bins.
#'
#' @return A list of motifs and their corresponding non-zero membership
#' probabilities to TADs.
#' @export
#'
#' @examples
#' data("mESChromInteractions")
#' data("matchedMotifs")
#' tads <- cluster(mESChromInteractions, offset=40000)
#' tads.cl <- classify(matchedMotifs, tads, offset=40000)
#'
classify <- function(motifs, tads, offset = 40000) {
    # check parameters are provided

    if (missing(motifs) | missing(tads)) {
        stop("classify needs both motifs and tads")
    } else if (!missing(tads) &&
               !(all(c("segment1", "segment2", "hcluster") %in%
                     colnames(tads)))) {
        stop("tads dataframe has improper column names")
    } else if (!missing(motifs) &&
               !(all(c("motif_id", "start", "stop") %in% colnames(motifs)))) {
        stop("motifs dataframe has improper column names")
    }

    # create a nmotifs x ntads dynamic matrix of 0s
    m <- dim(tads)[1]
    n <- dim(motifs)[1]
    M <- matrix(0, n, m)
    colnames(M) <- tads$hcluster
    rows <- motifs$motif_id
    rownames(M) <- rows

    # bin_motifs
    for (j in 1:(m - 1)) {
        M <-  bin_motifs(motifs, tads[j, ], j, M, offset = offset)
    }

    tmp <- rowsum(M, rownames(M))

    # return non-zero memberships for each motif
    membership <- vector("list", dim(tmp)[1])
    names(membership) <- rownames(tmp)

    # format output
    for (idx in 1:length(membership)) {
        curr <- tmp[idx,]
        nonzeros <-  curr[curr != 0]
        if (!is.null(nonzeros))
            membership[idx][[1]] <-
            (paste0(names(nonzeros), ": ", nonzeros))
    }

    return(membership)
}

#' Calculate Motif Membership Probabilty per TAD
#'
#' Assuming equal chance of assignment given the matched sequence, this function
#' calculates membership
#' proportion of motifs per TAD.
#'
#' @param motifs A dataframe of matched motifs with at least the following
#' attributes:
#' \itemize{
#' \item \emph{start}: an integer denoting the starting coordinate of matched
#' motif
#' \item \emph{stop}: an integer denoting the end coordinate of matched motif
#' }
#' @param j An integer denoting the index of the TAD for classification
#' @param tad A row of a dataframe containing hierarchically clustered TADs.
#' Compatible with output of the annotator3D::cluster function. The row should
#' have at least the following attributes:
#' \itemize{
#' \item \emph{segment1}: an integer denoting the starting coordinate of the
#' predicted TAD
#' \item \emph{segment2}: an integer denoting the end coordinate of the
#' predicted TAD
#' }
#' @param M A n-by-m numeric matrix to be used as a dynamic array of one TAD.
#' @param offset A numeric denoting the resolution of Hi-C data bins.
#'
#' @return M A n-by-m numeric matrix storing calculated membership scores.
#'
bin_motifs <- function(motifs, tad, j, M, offset) {
    n <- dim(motifs)[1]
    next_tad = c((tad$segment1 + offset), (tad$segment2 + offset)) # boundaries
    # of next TAD
    # handle case TAD has coordinates (0, offset)
    if (tad$segment1 != 0) {
        prev_tad = c((tad$segment1 -  offset), (tad$segment2 - offset))
    } else {
        prev_tad = c(0, tad$segment2)
    }

    # calculate membership score relative to TAD/sub-TAD boundaries.
    # Assuming that each TAD has an equal chance containing a motif in full,
    # on average a TAD has probability 1/n of matching to the motif. Since
    # motifs can cross TAD boundaries, however, what the TADs match to are
    # fractions of the motif - i.e. the motif belongs to 2+ TADs but not
    # necessarily to the same extent. Therefore, the probability either TAD
    # matches to the motif is (fraction of motif in TAD) * (1/n).
    #
    for (i in 1:n) {

        motif_len <- motifs$stop[i] - motifs$start[i]

        if (motifs$start[i] >= tad$segment1 &&
            motifs$stop[i] <= tad$segment2) {
            # motif is entirely contained in current TAD
            M[i, j] <- M[i, j] + (1 / n) # divide by n
        } else if (motifs$stop[i] >= tad$segment1 &&
                   motifs$stop[i] <= tad$segment2 &&
                   motifs$start[i] < (tad$segment1)) {
            # motif contained in current and previous TADs
            if (motifs$start[i] >= prev_tad[1]) {
                # motif contained entirely across current and previous TAD
                M[i, j] <-
                    M[i, j] + calculate_overflow(motifs$stop[i],
                                                 tad$segment1, motif_len) /
                    (n)
                M[i, j - 1] <-
                    M[i, j - 1] + calculate_overflow(motifs$start[i],
                                                     tad$segment1, motif_len) /
                    (n)
            }
        } else if (motifs$stop[i] > (tad$segment2) &&
                   motifs$start[i] >= tad$segment1 &&
                   motifs$start[i] <= tad$segment2) {
            # motif contained in current and next TAD
            if (motifs$stop[i] < next_tad[2])
            {
                # motif contained entirely across current and next TAD
                M[i, j + 1] <-
                    M[i, j + 1] + calculate_overflow(motifs$stop[i],
                                                     tad$segment2, motif_len) /
                    (n)
                M[i, j - 1] <-
                    M[i, j] + calculate_overflow(motifs$start[i],
                                                 tad$segment2, motif_len) /
                    (n)
            }
        }
    }

    return(M)
}

#' Calculate Amount of Motif Length to Hard TAD Boundary
#'
#' Calculates 'overflow' of a motif's length across TAD boundaries. It returns
#' the proportion of motif's length that overlaps with the TAD specified by
#' tad_thresh.
#'
#' @param motif_thresh An integer denoting the motif's coordinate.
#' @param tad_thresh An integer denoting the TAD boundary (hard)
#' @param motif_len An integer denoting the length of the motif
#'
#' @return A double measuring the proportion of the motif's length to the TAD
#' boundary.
#'
calculate_overflow <-
    function(motif_thresh, tad_thresh, motif_len) {
        # calculate proportion of motif in a TAD
        spill <- abs(motif_thresh - tad_thresh) / motif_len
        return(spill)
    }
