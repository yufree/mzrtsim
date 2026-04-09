test_that("simmzml generates mzML and CSV files", {
    skip_on_cran()
    data(monahrms1)
    tmpdir <- tempdir()
    outname <- file.path(tmpdir, "test_sim")

    simmzml(db = monahrms1, name = outname, n = 3,
            rtrange = c(0, 10), scanrate = 1, seed = 42)

    mzml_file <- paste0(outname, ".mzML")
    csv_file <- paste0(outname, ".csv")

    expect_true(file.exists(mzml_file))
    expect_true(file.exists(csv_file))

    # Check CSV contents
    df <- utils::read.csv(csv_file)
    expect_true("mz" %in% colnames(df))
    expect_true("rt" %in% colnames(df))
    expect_true("ins" %in% colnames(df))
    expect_true("sim_ins" %in% colnames(df))
    expect_true("name" %in% colnames(df))
    expect_true(nrow(df) > 0)
    expect_true(all(df$sim_ins >= 0))

    # Check mzML content
    msdata <- RaMS::grabMSdata(mzml_file, grab_what = "MS1", verbosity = 0L)$MS1
    expect_true(nrow(msdata) > 0)
    expect_true(all(msdata$mz >= 100 & msdata$mz <= 1000))

    # Clean up
    unlink(c(mzml_file, csv_file))
})

test_that("simmzml errors on mismatched rtime length", {
    data(monahrms1)
    tmpdir <- tempdir()
    outname <- file.path(tmpdir, "test_sim_err")

    expect_error(
        simmzml(db = monahrms1, name = outname, n = 3,
                rtime = c(10, 20), seed = 42),
        "same number of compounds"
    )
})

test_that("simmzml errors on mismatched baseline vector length", {
    data(monahrms1)
    tmpdir <- tempdir()
    outname <- file.path(tmpdir, "test_sim_bl")

    expect_error(
        simmzml(db = monahrms1, name = outname, n = 3,
                baseline = c(100, 200, 300), seed = 42),
        "baseline"
    )
})

test_that("simmzml_blank generates files", {
    skip_on_cran()
    tmpdir <- tempdir()
    outname <- file.path(tmpdir, "test_blank")

    simmzml_blank(name = outname, rtrange = c(0, 5), scanrate = 1)

    mzml_file <- paste0(outname, ".mzML")
    csv_file <- paste0(outname, ".csv")

    expect_true(file.exists(mzml_file))
    expect_true(file.exists(csv_file))

    # Clean up
    unlink(c(mzml_file, csv_file))
})
