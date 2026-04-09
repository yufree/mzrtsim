test_that("mzrtsim_se returns a SummarizedExperiment", {
    data(hmdbcms)
    se <- mzrtsim_se(ncomp = 10, ncond = 2, ncpeaks = 0.1,
                     nbatch = 2, nbpeaks = 0.1, npercond = 5,
                     nperbatch = c(5, 5), seed = 42,
                     batchtype = "mb", db = hmdbcms)

    expect_s4_class(se, "SummarizedExperiment")

    # Check dimensions
    expect_equal(ncol(se), 10)  # 5 per condition * 2 conditions
    expect_true(nrow(se) > 0)

    # Check assay
    expect_true("counts" %in% SummarizedExperiment::assayNames(se))

    # Check rowData
    rd <- SummarizedExperiment::rowData(se)
    expect_true("mz" %in% colnames(rd))
    expect_true("rt" %in% colnames(rd))
    expect_true("condition_peak" %in% colnames(rd))
    expect_true("batch_peak" %in% colnames(rd))

    # Check colData
    cd <- SummarizedExperiment::colData(se)
    expect_true("condition" %in% colnames(cd))
    expect_true("batch" %in% colnames(cd))
    expect_true("sample_name" %in% colnames(cd))
})
