% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qcbc.R
\name{bccenter}
\alias{bccenter}
\title{Centering batch correction}
\usage{
bccenter(peak, qc = NULL, log = TRUE)
}
\arguments{
\item{peak}{peaks intensity across samples}

\item{qc}{peaks intensity pooled QC samples}

\item{log}{log transformation}
}
\value{
corrected peaks intensity
}
\description{
Centering batch correction
}
