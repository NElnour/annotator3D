context("Motif Classification")
library(annotator3D)

data("mESChromInteractions")
data("matchedMotifs")

tads <- cluster(mESChromInteractions, offset = 40000)
res <- classify(matchedMotifs, tads)

motif <- matchedMotifs
colnames(motif) <-
    c(
        "x",
        "motif_alt_id",
        "sequence_name",
        "start",
        "stop",
        "strand",
        "score",
        "p.value",
        "q.value",
        "matched_sequence"
    )

test_that("classify needs motif_id column names outlined in docs", {
    expect_error(classify(motif, tads))
})

motif <- matchedMotifs
colnames(motif) <-
    c(
        "motif_id",
        "motif_alt_id",
        "sequence_name",
        "x",
        "stop",
        "strand",
        "score",
        "p.value",
        "q.value",
        "matched_sequence"
    )

test_that("classify needs start column name outlined in docs", {
    expect_error(classify(motif, tads))
})

motif <- matchedMotifs
colnames(motif) <-
    c(
        "motif_id",
        "motif_alt_id",
        "sequence_name",
        "start",
        "x",
        "strand",
        "score",
        "p.value",
        "q.value",
        "matched_sequence"
    )

test_that("classify needs stop column name outlined in docs", {
    expect_error(classify(motif, tads))
})

tad <- tads
colnames(tad) <- c("x", "segment2", "flat_cluster", "hcluster")

test_that("classify needs segment1 column name outlined in docs", {
    expect_error(classify(matchedMotifs, tad))
})

tad <- tads
colnames(tad) <- c("segmetn1", "x", "flat_cluster", "hcluster")

test_that("classify needs segment2 column name outlined in docs", {
    expect_error(classify(matchedMotifs, tad))
})

tad <- tads
colnames(tad) <- c("segmetn1", "segment2", "flat_cluster", "x")

test_that("classify needs hcluster column name outlined in docs", {
    expect_error(classify(matchedMotifs, tad))
})

test_that("classify needs both motifs and tads in order", {
    expect_error(classify())
    expect_error(classify(matchedMotifs))
    expect_error(classif(tads))
    expect_error(classify(tads, matchedMotifs))
})


test_that("classify outputs a list", {
    expect_type(res, "list")
})
