# mzrtsim 0.99.0

## New Features

* Added `simmzml_blank()` for generating blank/matrix-only mzML files.
* Added `sim_ins` column to CSV output of `simmzml()`, providing the absolute
  ground-truth maximum simulated intensity for each peak.
* `baseline` parameter in `simmzml()` now accepts a numeric vector to simulate
  gradient elution baseline shifts over time.

## Bug Fixes

* Updated GitHub Actions workflows to use current action versions.
