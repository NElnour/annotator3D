#' soften
#'
#' @return
#' @export
#' @import cluster
#' @import igraph
#'
#' @examples
fuzzify <-
    function(data, eigen_mat, k, distance, fit) {
        g <- igraph::graph.adjacency(data, mode="undirected", weighted = TRUE)

        clustered <- cluter::fanny(eigen_mat, k=k, metric = distance, memb.exp = fit)

        loci <- paste(rownames(data), colnames(data), sep = "-")
        cluster <- clustered$classification
        pval <- clustered$membership
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
