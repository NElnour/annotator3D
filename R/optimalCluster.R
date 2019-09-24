#' soften
#'
#' @return
#' @export
#' @import mclust
#'
#' @examples
soften <-
    function(data, eigen_mat) {
        clustered <- Mclust(eigen_mat)

        loci <- paste(rownames(data), colnames(data), sep = "-")
        cluster <- clustered$classification
        pval <- clustered$z
        names <- colnames(pval)

        for(i in 1:(dim(pval)[2])){
            col <- paste0("p", i)
            names[i] <- col
        }

        colnames(pval) <- names

        output <- data.frame(loci,
                             cluster,
                             pval)

        return(output)
    }
