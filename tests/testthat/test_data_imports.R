context("Data Imports Correctly")
library(annotator3D)

data(hsa)
data("matchedMotifs")

test_that("chr18 HSA dataset is loaded correctly", {
    expect_silent(data(hsa))
    expect_equal(dim(hsa)[1], 1000)
    expect_equal(dim(hsa)[2], 1000)
    expect_type(hsa, "double")
    expect_true(class(hsa)== "matrix")
})

test_that("FIMO matched motifs' dataset is loaded correctly", {
    expect_silent(data("matchedMotifs"))
    expect_equal(dim(matchedMotifs)[1], 50)
    expect_equal(dim(matchedMotifs)[2], 10)
    expect_true(class(matchedMotifs)== "data.frame")
    expect_type(matchedMotifs$motif_id, "integer")
    expect_type(matchedMotifs$start, "integer")
    expect_type(matchedMotifs$stop, "integer")
    
})