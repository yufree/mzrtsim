#' Shiny application for batch simulation and correction based on real data
#' @export
runmzrtsim <- function() {
        file <- system.file("sim.Rmd",
                            package = "mzrtsim")
        if (file == "") {
                stop("Could not find directory. Try re-installing `mzrtsim`.",
                     call. = FALSE)
        }
        rmarkdown::run(file)
}

#' Shiny application for raw data simulation
#' @export
runmzMLsim <- function() {
        file <- system.file("mzMLsim.Rmd",
                            package = "mzrtsim")
        if (file == "") {
                stop("Could not find directory. Try re-installing `mzrtsim`.",
                     call. = FALSE)
        }
        rmarkdown::run(file)
}
