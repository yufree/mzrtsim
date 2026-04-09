# mzrtsim 0.99.1

## New Features

* Implemented a **native mzML writer** using `base64enc`, removing the heavyweight dependency on `mzR` for data generation.
* Extended `simmzml()` with support for **MS2/MSn data simulation**, including parameters for precursor ion details, collision energy, and isolation windows.

# mzrtsim 0.99.0

## New Features

* Added `mzrtsim_se()` to return peak list simulation results as a
  `SummarizedExperiment` for seamless integration with Bioconductor workflows.
* Added `simmzml_blank()` for generating blank/matrix-only mzML files.
* Added `sim_ins` column to CSV output of `simmzml()`, providing the absolute
  ground-truth maximum simulated intensity for each peak.
* `baseline` parameter in `simmzml()` now accepts a numeric vector to simulate
  gradient elution baseline shifts over time.

## Bug Fixes

* Updated GitHub Actions workflows to use current action versions.
