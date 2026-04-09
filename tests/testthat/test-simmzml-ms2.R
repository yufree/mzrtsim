test_that("simmzml supports MS2 parameters", {
    # Setup mock database
    mock_db <- list(
        list(
            name = "FragmentComp",
            spectra = list(
                mz = c(50.1, 100.2, 150.3),
                ins = c(500, 1000, 500)
            )
        )
    )

    tmpdir <- tempdir()
    outname <- file.path(tmpdir, "sim_ms2")
    mzml_file <- paste0(outname, ".mzML")

    # Run simmzml with MS2 parameters
    simmzml(
        db = mock_db,
        name = outname,
        n = 1,
        compound = 1,
        ms_level = 2L,
        precursor_mz = 250.5,
        precursor_charge = 1L,
        collision_energy = 30,
        rtrange = c(0, 5),
        scanrate = 5,
        seed = 42
    )

    expect_true(file.exists(mzml_file))

    # Verify content with RaMS
    ms2_data <- RaMS::grabMSdata(mzml_file, grab_what = "MS2", verbosity = 0L)$MS2
    expect_true(nrow(ms2_data) > 0)
    expect_equal(unique(ms2_data$premz), 250.5, tolerance = 1e-4)

    # Clean up
    unlink(c(mzml_file, paste0(outname, ".csv")))
})
