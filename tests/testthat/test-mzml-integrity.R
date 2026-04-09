test_that("mzML integrity is preserved", {
    # 1. Setup mock database
    mock_db <- list(
        list(
            name = "Comp1",
            spectra = list(
                mz = c(100.12345, 200.54321),
                ins = c(1000, 500)
            )
        )
    )

    tmpdir <- tempdir()
    outname <- file.path(tmpdir, "integrity_test")
    mzml_file <- paste0(outname, ".mzML")

    # 2. Run simulation with fixed parameters
    # Using very simple parameters: 1 scan total, fixed RT
    simmzml(
        db = mock_db,
        name = outname,
        n = 1,
        rtrange = c(5, 5), # Results in one scan at 5s
        scanrate = 1,
        baseline = 0,
        baselinesd = 0,
        ppm = 0,
        sampleppm = 0,
        seed = 42,
        compound = 1,
        rtime = 5
    )

    expect_true(file.exists(mzml_file))

    # 3. Read back using RaMS
    # Note: RaMS is a dependency of the package
    msdata <- RaMS::grabMSdata(mzml_file, grab_what = "MS1", verbosity = 0L)$MS1

    # 4. Assertions
    # mzdigit defaults to 5 in simmzml
    expected_mz <- round(c(100.12345, 200.54321), 5)

    # Check length
    expect_equal(nrow(msdata), 2)

    # Check RT (RaMS returns RT in minutes by default, but write_mzml writes in seconds)
    # Actually RaMS grabMSdata documentation: rt: retention time in minutes
    # Let's check simmzml code: rtime0 <- seq(rtrange[1], rtrange[2], scanrate)
    # rtime0 is in seconds. write_mzml writes it as "second" unit.
    expect_equal(unique(msdata$rt) * 60, 5, tolerance = 1e-5)

    # Check m/z
    expect_equal(sort(msdata$mz), sort(expected_mz), tolerance = 1e-5)

    # Check Intensity
    # In simmzml: nret <- matrix(intensity[[i]])%*%re[i,]
    # re[i,] is dnorm(...) scaled to 100*rf[i]*peakheight[i]
    # For a single scan at the mean of dnorm, re[i,] should be 100*rf*pheight
    # Default rf=100, pheight=10. So re = 100*100*10 = 1e5.
    # Intensity should be ins * 1e5
    # Comp1 ins = 1000 and 500. So 1e8 and 5e7.
    expect_equal(sort(msdata$int), sort(c(1000, 500) * 1e5), tolerance = 1e-2)

    # Clean up
    unlink(c(mzml_file, paste0(outname, ".csv")))
})

test_that("mzML header contains correct metadata", {
    mock_db <- list(list(name="test", spectra=list(mz=100, ins=100)))
    tmpdir <- tempdir()
    outname <- file.path(tmpdir, "header_test")
    mzml_file <- paste0(outname, ".mzML")

    simmzml(db = mock_db, name = outname, n = 1, seed = 42, compound = 1, rtrange = c(0,0))

    lines <- readLines(mzml_file, n = 40)
    # Check for mzrtsim software tag
    expect_true(any(grepl('<software id="mzrtsim" version="0.99.0">', lines)))
    expect_true(any(grepl('custom unreleased software tool', lines)))

    unlink(c(mzml_file, paste0(outname, ".csv")))
})
