test_that("bccenter returns numeric output", {
    peak <- c(100, 200, 300, 400, 500)
    qc <- c(150, 250, 350)

    result_log <- bccenter(peak, qc, log = TRUE)
    expect_type(result_log, "double")
    expect_length(result_log, length(peak))

    result_nolog <- bccenter(peak, qc, log = FALSE)
    expect_type(result_nolog, "double")
    expect_length(result_nolog, length(peak))

    # Log vs non-log should differ
    expect_false(identical(result_log, result_nolog))
})

test_that("bccenter without QC works", {
    peak <- c(100, 200, 300, 400, 500)
    result <- bccenter(peak, qc = NULL, log = TRUE)
    expect_type(result, "double")
    expect_length(result, length(peak))
})

test_that("bcscaling returns numeric output", {
    peak <- c(100, 200, 300, 400, 500)
    qc <- c(150, 250, 350)

    result <- bcscaling(peak, qc, log = TRUE)
    expect_type(result, "double")
    expect_length(result, length(peak))
})

test_that("bcpareto returns numeric output", {
    peak <- c(100, 200, 300, 400, 500)
    qc <- c(150, 250, 350)

    result <- bcpareto(peak, qc, log = TRUE)
    expect_type(result, "double")
    expect_length(result, length(peak))
})

test_that("bcrange returns numeric output", {
    peak <- c(100, 200, 300, 400, 500)
    qc <- c(150, 250, 350)

    result <- bcrange(peak, qc, log = TRUE)
    expect_type(result, "double")
    expect_length(result, length(peak))
})

test_that("bcvast returns numeric output", {
    peak <- c(100, 200, 300, 400, 500)
    qc <- c(150, 250, 350)

    result <- bcvast(peak, qc, log = TRUE)
    expect_type(result, "double")
    expect_length(result, length(peak))
})

test_that("bclevel returns numeric output", {
    peak <- c(100, 200, 300, 400, 500)
    qc <- c(150, 250, 350)

    result <- bclevel(peak, qc, log = TRUE)
    expect_type(result, "double")
    expect_length(result, length(peak))
})
