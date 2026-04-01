test_that("mzrtsim returns expected structure", {
    data(hmdbcms)
    sim <- mzrtsim(ncomp = 10, ncond = 2, ncpeaks = 0.1,
                   nbatch = 2, nbpeaks = 0.1, npercond = 5,
                   nperbatch = c(5, 5), seed = 42, batchtype = "mb",
                   db = hmdbcms)

    expect_type(sim, "list")
    expect_true("data" %in% names(sim))
    expect_true("mz" %in% names(sim))
    expect_true("rt" %in% names(sim))
    expect_true("group" %in% names(sim))
    expect_true("name" %in% names(sim))
    expect_true("conp" %in% names(sim))
    expect_true("batchp" %in% names(sim))

    # data should have 10 columns (5 per condition * 2 conditions)
    expect_equal(ncol(sim$data), 10)
    # group should have 10 rows
    expect_equal(nrow(sim$group), 10)
})

test_that("mzrtsim errors without database", {
    expect_error(mzrtsim(ncomp = 10, db = NULL),
                 "You need database to simulate")
})

test_that("mzrtsim errors with mismatched batch/condition columns", {
    data(hmdbcms)
    expect_error(
        mzrtsim(ncomp = 10, ncond = 2, npercond = 5,
                nbatch = 3, nperbatch = c(3, 3, 3), db = hmdbcms),
        "same numbers"
    )
})
