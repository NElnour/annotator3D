#' Fuzzy Cluster Membership with C-means
#' 
#' Adds probability of cluster membership to gene sequence
#'
#' @return A dataframe of chrosome regions in contact, their predicted TAD
#' clusters, and the probability of their membership.
#' @export
#' @import cluster
#'
fuzzify <-
    function(data, eigen_mat, k, distance, fit) {
        
        clustered <- cluster::fanny(eigen_mat, k=k, metric = distance, memb.exp = fit)

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
