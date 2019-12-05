context("Data Imports Correctly")
library(annotator3D)

data("A549ChromLoops")
data("matchedMotifs")

test_that("chr18 chromLoops dataset is loaded correctly", {
    expect_silent(data(A549ChromLoops))
    expect_equal(dim(A549ChromLoops)[1], 216)
    expect_equal(dim(A549ChromLoops)[2], 3)
    expect_true(class(A549ChromLoops)== "data.frame")
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
