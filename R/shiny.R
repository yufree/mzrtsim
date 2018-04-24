#' Shiny application for batch simulation and correction based on real data
#' @export
runmzrtsim <- function() {
        file <- system.file("shinyapp", "sim.Rmd",
                            package = "mzrtsim")
        if (file == "") {
                stop("Could not find directory. Try re-installing `mzrtsim`.",
                     call. = FALSE)
        }
        rmarkdown::run(file)
}
