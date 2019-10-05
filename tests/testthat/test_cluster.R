context("Hierarchical Clustering")
library(annotator3D)

data(hsa)
res <- cluster(hsa, 40000)

test_that("cluster has four columns", {
    expect_equal(dim(res)[2], 4)
    expect_type(res$segment1, "double")
    expect_type(res$segment2, "double")
    expect_type(res$flat_cluster, "integer")
    expect_type(res$hcluster, "integer")
})

test_that("cluster produces 1000 rows on hsa", {
    expect_equal(dim(res)[1], 1000)
})

test_that("cluster needs both hsa and offset", {
    expect_error(cluster())
    expect_error(cluster(hsa))
    expect_error(offset = 40000)
})

test_that("cluster hcluster has hierarchical info", {
    expect_true(any(sapply(res$hcluster, function(x) grepl("-", res$hcluster[x]))))
})

test_that("cluster produces a dendrogram image", {
    expect_true("Rplots.pdf" %in% list.files(".") )
})
