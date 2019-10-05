#' Hierarchical Hard-Clustering of Hi-C HSA Matrix
#' 
#' This is a wrapper function streamlining hierarchical clustering of Hi-C contact matrices into topologically-associating domains(TADs). It uses hclust and dendrogram cutting to predict TADs and sub-TADs.
#'
#' @param hsa a Hi-C numeric matrix of chromatin interactions
#' @param offset An integer denoting the binning resolution of the Hi-C dataset
#' @param cluster_method a string indicating the agglomeration method to be used for hierarchical clustering. It is "ward.D2" by default, but can also be set to "ward.D", "single", "complete", "average" or "centroid".
#' @param plot_subtads A boolean denoting whether predicted sub-TADs should be plotted individually
#' @param dend_labels a string one of "none", "perpindicular" or "textlike" for no, vertical, or horizontal labels, respectively.
#'
#' @return dataframe
#' @export
#'
#' @examples
#' data(hsa)
#' data.cl <- cluster(hsa, offset=40000)
#' 
#' @importFrom stats hclust dist as.dendrogram
#' @importFrom dendextend cutree_1k.dendrogram color_branches color_labels
#' @importFrom graphics plot
#' @import utils
#' 
cluster <- function(hsa, offset, cluster_method="ward.D2", plot_subtads=FALSE, dend_labels="none") {
    
    if(missing(hsa) | missing(offset)){
        stop("cluster needs both HSA and ffset arguments")
    }
    
    k <- floor(dim(hsa)[1] / 2) - 1
    hsa_scaled <- scale(hsa)
    hsa.d <- stats::dist(hsa_scaled, method = "euclidean")
    hsa.cl <- stats::hclust(hsa.d, method = cluster_method)
    
    print(hsa.cl)
    
    dendrogram <- stats::as.dendrogram(hsa.cl)
    
    tads <- dendextend::cutree_1k.dendrogram(dendrogram, k = k)
    dendrogram <- dendextend::color_branches(dendrogram, k = k)
    
    if (dend_labels != "none")
        dendrogram <- dendextend::color_labels(dendrogram, k = k)
    
    graphics::plot(dendrogram, leaflab=dend_labels)
    subtads <- tad.group(dendrogram, tads=tads, k, plot_subtads = plot_subtads, offset = offset)
    return(subtads)
}

#' Label TADs into Flat and Hierarchical Clusters
#' 
#' Finds and groups sub-TADs and TADs into hierarchical clusters, with the option of the plotting them separatelyy.
#'
#' @param dendrogram A dendrogram (list) produced by hclust.
#' @param tads A cut dendrogram object into $$k$$ objects
#' @param k The number of groups into which the dendrogram is cut
#' @param plot_subtads A boolean specifying whether sub-TADs should be plotted
#' @param offset An integer denoting the binning resolution of the Hi-C dataset
#'
#' @return A dataframe with four attributes:
#' \itemize{
#' \item \emph{segment1}: the start coordinate of a Hi-C bin
#' \item \emph{segment2}: the end coordinate of a Hi-C bin
#' \item \emph{flat_cluster}: predicted cluster group of TADs output by hcluster
#' \item \emph{hcluster}: full hierarchical subclusters of the TAD given the flat_cluster
#' }
#' 
#' @importFrom dendextend get_subdendrograms
#' @importFrom stringr str_replace_all
#' @importFrom graphics plot
#'
tad.group <- function(dendrogram, tads, k, plot_subtads, offset) {
    subtads <- dendextend::get_subdendrograms(dendrogram, k)
    
    if (plot_subtads)
        sapply(subtads, graphics::plot)
    
    labelled_tads <- names(LinearizeNestedList(subtads))
    labelled_tads <- stringr::str_replace(labelled_tads, "/", "-")
    labelled_tads <- stringr::str_replace(labelled_tads, "^[0-9]*", as.character(tads))
    segment1 <- as.numeric(names(tads)) - 40000
    segment2 <- segment1 + offset
    
    collate <- data.frame(segment1 = segment1,
                          segment2 = segment2,
                          flat_cluster = tads,
                          hcluster = labelled_tads)
    return(collate)
}
