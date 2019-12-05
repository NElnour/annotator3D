context("Motif Classification")
library(annotator3D)

data("A549ChromLoops")
data("matchedMotifs")

results <- classify(matchedMotifs, A549ChromLoops)

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
    expect_error(classify(motif, clusteredAnnotationTrack))
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
    expect_error(classify(motif, clusteredAnnotationTrack))
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
    expect_error(classify(motif, A549ChromLoops))
})

tad <- A549ChromLoops
colnames(tad) <- c("x", "stop", "flat_cluster")

test_that("classify needs segment1 column name outlined in docs", {
    expect_error(classify(matchedMotifs, tad))
})

tad <- A549ChromLoops
colnames(tad) <- c("start", "x", "flat_cluster")

test_that("classify needs segment2 column name outlined in docs", {
    expect_error(classify(matchedMotifs, tad))
})


test_that("classify needs both motifs and A549ChromLoops in order", {
    expect_error(classify())
    expect_error(classify(matchedMotifs))
    expect_error(classif(A549ChromLoops))
    expect_error(classify(A549ChromLoops, matchedMotifs))
})
