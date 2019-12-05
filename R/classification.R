#' Motif Soft Classification into Topologically-associating Domains
#'
#' Classifies motifs by location into hierarchically clustered chromatin loops.
#' It assumes that each base position has an equal chance of occurring in a
#' chromatin loop.
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
#' @param chromLoops A dataframe of hierarchically clustered chromatin loops
#' with at least the
#' following attributes:
#' #' \itemize{
#' \item \emph{start}: an integer array containing the starting coordinates
#' of the predicted chromatin loops
#' \item \emph{stop}: an integer array containing the end coordinates of
#' the predicted chromatin loops
#' }
#'
#' @export
#'
#' @return A dataframe of classified of motifs, their membership scores, and
#' locations
#'
#' @examples
#' data("A549ChromLoops")
#' data("matchedMotifs")
#' classify(matchedMotifs, A549ChromLoops)
#'
#'@importFrom stats complete.cases
#'
classify <- function(motifs, chromLoops) {
    # check parameters are provided
    if (missing(motifs) | missing(chromLoops)) {
        stop("classify needs both motifs and chromLoops")
    } else if (!missing(chromLoops) &&
               !(all(c("start", "stop") %in%
                     colnames(chromLoops)))) {
        stop("chromLoops dataframe has improper column names")
    } else if (!missing(motifs) &&
               !(all(c("motif_id", "start", "stop") %in% colnames(motifs)))) {
        stop("motifs dataframe has improper column names")
    }

    chrom <- unique(chromLoops$chr)
    chrom <- chrom[stats::complete.cases(chrom)]

    # create a nmotifs x nchromLoops dynamic matrix of 0s
    m <- dim(chromLoops)[1]
    n <- dim(motifs)[1]
    M <- matrix(0, n, m)
    colnames(M) <- paste0(chromLoops$start, "-", chromLoops$stop)
    if (all(is.na(motifs$motif_alt_id)))
    {
        rows <- motifs$motif_id
    } else {
        rows <- motifs$motif_alt_id
    }

    rownames(M) <- rows

    # bin_motifs
    for (j in 1:(m - 1)) {
        M <-  bin_motifs(motifs, chromLoops[j, ], j, M, m)
    }

    tmp <- rowsum(M, rownames(M))

    # return non-zero memberships for each motif
    all_motifs <- data.frame(
        chrom = character(),
        start = integer(),
        end = integer(),
        symbol = character(),
        score = numeric()
    )

    membership <- vector("list", dim(tmp)[1])
    names(membership) <- rownames(tmp)

    for (idx in 1:length(membership)) {
        name <- names(membership[idx])
        curr <- tmp[idx, ]
        nonzeros <-  curr[curr != 0]
        if (!is.null(nonzeros))
        {
            motif <- write_to_BED(nonzeros, name, chrom)
            all_motifs <- rbind(all_motifs, motif)
        }
    }

    return(all_motifs)
}

#' Calculate Motif Membership Probabilty per chromatin loop
#'
#' Assuming equal chance of assignment given the matched sequence, this function
#' calculates membership
#' proportion of motifs per chromatin loop.
#'
#' @param motifs A dataframe of matched motifs with at least the following
#' attributes:
#' \itemize{
#' \item \emph{start}: an integer denoting the starting coordinate of matched
#' motif
#' \item \emph{stop}: an integer denoting the end coordinate of matched motif
#' }
#' @param j An integer denoting the index of the chromatin loop for classification
#' @param chromLoops A row of a dataframe containing hierarchically clustered chromatin loops.
#' Compatible with output of the annotator3D::cluster function. The row should
#' have at least the following attributes:
#' \itemize{
#' \item \emph{start}: an integer denoting the starting coordinate of the
#' predicted chromatin loop
#' \item \emph{stop}: an integer denoting the end coordinate of the
#' predicted chromatin loop
#' }
#' @param M An n-by-m numeric matrix to be used as a dynamic array of one chromatin loop.
#'
#' @param m An integer denoting the number of predicted loops.
#'
#' @return M An n-by-m numeric matrix storing calculated membership scores.
#'
bin_motifs <- function(motifs, chromLoops, j, M, m) {
    n <- dim(motifs)[1]

    # calculate membership score relative to chromatin loop/sub-chromatin loop boundaries.
    # Assuming that each chromatin loop has an equal chance containing a motif in full,
    # on average a chromatin loop has probability 1/n of matching to the motif. Since
    # motifs can cross chromatin loop boundaries, however, what the chromatin loops match to are
    # fractions of the motif - i.e. the motif belongs to 2+ chromatin loops but not
    # necessarily to the same extent. Therefore, the probability either chromatin loop
    # matches to the motif is (fraction of motif in chromatin loop) * (1/n).
    #
    for (i in 1:n) {
        motif_len <- motifs$stop[i] - motifs$start[i]

        if (motifs$start[i] >= chromLoops$start &&
            motifs$stop[i] <= chromLoops$stop) {
            # motif is entirely contained in current chromatin loop
            M[i, j] <- M[i, j] + (1 / m)

        } else if (motifs$stop[i] >= chromLoops$start &&
                   motifs$stop[i] <= chromLoops$stop &&
                   motifs$start[i] < (chromLoops$start)) {
            # motif contained in current and previous chromatin loops
            M[i, j] <-
                M[i, j] + calculate_overflow(motifs$stop[i],
                                             chromLoops$start, motif_len) /
                (m)
            M[i, j - 1] <-
                M[i, j - 1] + calculate_overflow(motifs$start[i],
                                                 chromLoops$start, motif_len) /
                (m)
        } else if (motifs$stop[i] > (chromLoops$stop) &&
                   motifs$start[i] >= chromLoops$start &&
                   motifs$start[i] <= chromLoops$stop) {
            # motif contained in current and next chromatin loop

            M[i, j + 1] <-
                M[i, j + 1] + calculate_overflow(motifs$stop[i],
                                                 chromLoops$stop, motif_len) /
                (m)
            M[i, j - 1] <-
                M[i, j] + calculate_overflow(motifs$start[i],
                                             chromLoops$stop, motif_len) /
                (m)

        }
    }

    return(M)
}

#' Calculate Amount of Motif Length to Hard chromatin loop Boundary
#'
#' Calculates 'overflow' of a motif's length across chromatin loop boundaries. It returns
#' the proportion of motif's length that overlaps with the chromatin loop specified by
#' chromLoop_thresh.
#'
#' @param motif_thresh An integer denoting the motif's coordinate.
#' @param chromLoop_thresh An integer denoting the chromatin loop boundary (hard)
#' @param motif_len An integer denoting the length of the motif
#'
#' @return A double measuring the proportion of the motif's length to the chromatin loop
#' boundary.
#'
calculate_overflow <-
    function(motif_thresh,
             chromLoop_thresh,
             motif_len) {
        # calculate proportion of motif in a chromatin loop
        spill <- abs(motif_thresh - chromLoop_thresh) / motif_len
        return(spill)
    }


#' Write Motif Vector to BED File
#'
#' Write a motif and its location to a BED file stored locally.
#'
#' @param motif_ranges A numeric vector of motif's membership scores. The names
#' of the entries correspond to the motif's location on the chromsome chrom.
#' @param motif_name A string containing the name of the motif to be stored.
#' @param chrom A string denoting the chromosome containing the motif and
#' chromatin loop.
#'
#' @importFrom stringr str_split_fixed str_extract
#' @importFrom utils write.table
#'
write_to_BED <- function(motif_ranges, motif_name, chrom) {
    # check if storage directory is available. Create it, if not.

    wd <- getwd()

    if (!dir.exists("BEDs")) {
        dir.create("BEDs")
    }

    # parse classified motifs list in BED format
    range <- stringr::str_split_fixed(names(motif_ranges), "[-]", 2)

    storage <- paste0(wd, "/", "BEDs")

    bed <- data.frame(
        chrom = chrom,
        start = as.integer(range[, 1]),
        end = as.integer(range[, 2]),
        symbol = motif_name,
        score = as.numeric(motif_ranges)
    )

    utils::write.table(
        bed,
        file = paste0(storage, "/", motif_name, ".bed"),
        quote = F,
        sep = "\t",
        row.names = F,
        col.names = F,
        append = T
    )

    return(bed)
}

# [END]
