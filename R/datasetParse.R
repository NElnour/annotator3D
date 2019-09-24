#' Title
#'
#' @return
#' @export
#'
#' @examples
readMotifs <- function(){

}

#' Title
#'
#' @return
#' @export
#'
#' @examples
motifParse <- function(){

}

#' Title
#'
#' @return
#' @export
#'
#' @examples
readHiC <- function(HSA_filename){
    data <- read.delim(HSA_filename, header=FALSE)
    end <- dim(data)[2] - 1

    contact_mat <- as.matrix(data[, 3:end])
    rownames(contact_mat) <- data$V1
    colnames(contact_mat) <- data$V2

    return(contact_mat)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
hicParse <- function(){

}
