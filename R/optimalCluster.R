#' Fuzzy Cluster Membership with C-means
#'
#' Adds probability of contact cluster membership to TADs
#'
#' @return A dataframe of chrosome regions in contact, their predicted TAD
#' clusters, and the probability of their membership.
#' @export
#'
fuzzify <-
    function(hard_clusters) {
        loci <- paste(rownames(data), colnames(data), sep = "-")
        cluster <- clustered$classification
        pval <- clustered$membership
        names <- colnames(pval)
        
        for (i in 1:(dim(pval)[2])) {
            col <- paste0("p", i)
            names[i] <- col
        }
        
        colnames(pval) <- names
        
        output <- data.frame(loci,
                             cluster,
                             pval)
        
        return(output)
    }

cluster <- function(data, cluster_method="ward.D2") {
    k <- floor(dim(data)[1] / 2) - 1
    data_scaled <- scale(data)
    data.d <- stats::dist(data_scaled, method = "euclidean")
    data.cl <- stats::hclust(data.d, method = cluster_method)
    
    print(data.cl)
    
    dendrogram <- as.dendrogram(data.cl)
    
    tads <- dendextend:::cutree.dendrogram(dendrogram, k = k)
    dendrogram <- dendextend::color_branches(dendrogram, k = k)
    dendrogram <- dendextend::color_labels(dendrogram, k = k)
    
    plot(dendrogram)
    subtads <- tad.group(dendrogram, tads=tads, k)
    
}

tad.group <- function(dendrogram, tads, k, plot_subtads = FALSE) {
    subtads <- dendextend::get_subdendrograms(dendrogram, k)
    
    if (plot_subtads)
        sapply(subtads, plot)
    
    labelled_tads <- names(LinearizeNestedList(subtads))
    labelled_tads <- stringr::str_replace_all(labelled_tads, "/", "-")
    
    collate <- data.frame(segment = names(tads),
                          flat_cluster = tads,
                          hcluster = labelled_tads)
    return(collate)
}
