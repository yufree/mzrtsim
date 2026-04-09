test_that("write_mzml supports MS2 data", {
    # 1. Setup mixed MS1/MS2 data
    mz_list <- list(
        c(100.1, 200.2, 300.3),          # Scan 1 (MS1)
        c(50.05, 100.10, 150.15)         # Scan 2 (MS2 of 200.2)
    )
    int_list <- list(
        c(1000, 2000, 3000),
        c(500, 1000, 500)
    )
    rt <- c(1, 1.2)
    ms_lvls <- c(1L, 2L)
    prec_mz <- c(NA, 200.2)
    prec_charge <- c(NA, 1L)
    prec_int <- c(NA, 2000)
    coll_eng <- c(NA, 35.0)

    tmpdir <- tempdir()
    mzml_file <- file.path(tmpdir, "ms2_test.mzML")

    # 2. Extract internal write_mzml (it's not exported, so we use ::: or just load_all)
    # Since we are in testthat and devtools::test() loads the package, 
    # we can use mzrtsim:::write_mzml if needed, 
    # but I'll use the package name to be safe if it's not exported.
    
    mzrtsim:::write_mzml(
        mz_list = mz_list,
        intensity_list = int_list,
        rtime = rt,
        file = mzml_file,
        ms_level = ms_lvls,
        precursor_mz = prec_mz,
        precursor_charge = prec_charge,
        precursor_intensity = prec_int,
        collision_energy = coll_eng
    )

    expect_true(file.exists(mzml_file))

    # 3. Verify using RaMS
    # RaMS grabMSdata for MS1
    ms1_data <- RaMS::grabMSdata(mzml_file, grab_what = "MS1", verbosity = 0L)$MS1
    expect_equal(nrow(ms1_data), 3)
    expect_equal(unique(ms1_data$rt) * 60, 1)

    # RaMS grabMSdata for MS2
    ms2_data <- RaMS::grabMSdata(mzml_file, grab_what = "MS2", verbosity = 0L)$MS2
    expect_equal(nrow(ms2_data), 3)
    expect_equal(unique(ms2_data$rt) * 60, 1.2)
    
    # RaMS returns precursor m/z as 'premz' column
    expect_equal(unique(ms2_data$premz), 200.2, tolerance = 1e-4)

    # Clean up
    unlink(mzml_file)
})

test_that("write_mzml generates valid isolation window metadata for MS2", {
    mz_list <- list(c(50, 100))
    int_list <- list(c(100, 100))
    rt <- 1.0
    ms_lvls <- 2L
    prec_mz <- 200.2
    
    tmpdir <- tempdir()
    mzml_file <- file.path(tmpdir, "ms2_meta_test.mzML")

    mzrtsim:::write_mzml(
        mz_list = mz_list,
        intensity_list = int_list,
        rtime = rt,
        file = mzml_file,
        ms_level = ms_lvls,
        precursor_mz = prec_mz,
        isolation_window = 1.0 # +/- 1.0 Da
    )

    lines <- readLines(mzml_file)
    # Check for isolation window target, lower and upper offsets
    expect_true(any(grepl('name="isolation window target m/z" value="200.200000"', lines)))
    expect_true(any(grepl('name="isolation window lower offset" value="1.000000"', lines)))
    expect_true(any(grepl('name="isolation window upper offset" value="1.000000"', lines)))
    # Check for collision energy
    # We didn't provide collision_energy in this specific call, so it shouldn't be there or should be NA handled.
    # In write_mzml: if (!is.null(collision_energy) && !is.na(collision_energy[i])) ...
    expect_false(any(grepl('name="collision energy"', lines)))

    unlink(mzml_file)
})
