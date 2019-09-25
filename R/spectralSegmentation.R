#' Spectral segmentation of chromatin interaction/contact matrices
#'
#' Creates eigenvectors of clustered components in the chromatin interaction/contact matrix. This function is authored by Nura Kawa (2018).
#' @author Nura Kawa
#' @references{\cite{Kawa, N. 2018. Spectral Clustering. RPubs. URL: \url{https://rpubs.com/nurakawa/spectral-clustering}}}
#'
#' @param x the chromatin interaction/contact matrix
#' @param nn number of k-nearest neighbours to consider
#' @param n_eig number of eigenvectors to keep
#'
#' @return eigenvectors of clustered components in X for
#'  kmeans analysis
#'
#' @examples
#' clustered <- spectral_clustering(X)
#'
#' @export

spectral_clustering <- function(x, nn = 10, n_eig = 2)
{

    W = mutual_knn_graph(x) # 1. matrix of similarities
    L = graph_laplacian(W) # 2. compute graph laplacian
    ei = eigen(L, symmetric = TRUE) # 3. Compute the eigenvectors and values of L
    n = nrow(L)
    return(ei$vectors[,(n - n_eig):(n - 1)]) # return the eigenvectors of the n_eig smallest eigenvalues

}

#' Nearest neighbourdetection based on Euclidean distances
#'
#' Finds nn nearest Euclidean neighbours. This function is authored by Nura Kawa (2018).
#' @author Nura Kawa
#' @references{\cite{Kawa, N. 2018. Spectral Clustering. RPubs. URL: \url{https://rpubs.com/nurakawa/spectral-clustering}}}
#'
#'
#' @param x the chromatin interaction/contact matrix
#' @param nn number of k-nearest neighbours to consider
#'
#' @return matrix
#'

mutual_knn_graph <- function(x, nn = 10)
{
    D <- as.matrix( dist(x) ) # matrix of euclidean distances between data points in X

    # intialize the knn matrix
    knn_mat <- matrix(0,
                      nrow = nrow(x),
                      ncol = nrow(x))

    # find the 10 nearest neighbors for each point
    for (i in 1: nrow(x)) {
        neighbor_index <- order(D[i,])[2:(nn + 1)]
        knn_mat[i,][neighbor_index] <- 1
    }

    # Now we note that i,j are neighbors iff K[i,j] = 1 or K[j,i] = 1
    knn_mat <- knn_mat + t(knn_mat) # find mutual knn

    knn_mat[ knn_mat == 2 ] = 1

    return(knn_mat)
}

#' Laplacian of contact matrix
#'
#' Calculates the Laplacian transformation of the contact matrix. This function is authored by Nura Kawa (2018).
#' @author Nura Kawa
#' @author Nada Elnour
#' @references{\cite{Kawa, N. 2018. Spectral Clustering. RPubs. URL: \url{https://rpubs.com/nurakawa/spectral-clustering}}}
#'
#' @param W similarity matrix
#'
#' @return matrix Laplacian transformation of the contact matrix
#'

graph_laplacian <- function(W)
{
    stopifnot(nrow(W) == ncol(W))

    g = colSums(W) # degrees of vertices
    n = nrow(W)

    return( diag(g) - W )

}
