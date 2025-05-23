% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mzmlsim.R
\name{simmzml}
\alias{simmzml}
\title{Generate simulated mzml data and compounds list}
\usage{
simmzml(
  db,
  name,
  n = 100,
  inscutoff = 0.05,
  mzrange = c(100, 1000),
  rtrange = c(0, 600),
  ppm = 5,
  sampleppm = 5,
  mzdigit = 5,
  noisesd = 0.5,
  scanrate = 0.2,
  pwidth = 10,
  pheight = 10,
  baseline = 100,
  baselinesd = 30,
  rf = 100,
  tailingfactor = 1.2,
  compound = NULL,
  rtime = NULL,
  tailingindex = NULL,
  seed = 42,
  unique = FALSE,
  matrix = FALSE,
  matrixmz = NULL
)
}
\arguments{
\item{db}{compound database with MS1 data. e.g. hmdbcms or monams1}

\item{name}{file name of mzml}

\item{n}{compound numbers from database, if compound is not NULL, n will be compound number, default 100}

\item{inscutoff}{intensity cutoff for MS1 spectra, default 0.05}

\item{mzrange}{m/z range for simulation, peaks out of the range will be removed, default c(100,1000)}

\item{rtrange}{retention time range for simulation, default c(0,600)}

\item{ppm}{m/z shift in ppm for input compounds, default 5}

\item{sampleppm}{m/z shift in ppm within one sample, default 5}

\item{mzdigit}{m/z digits, default 5}

\item{noisesd}{standard deviation for normal distribution of m/z shift, default 0.5}

\item{scanrate}{time for each full scan, default 0.2 second or 5 spectra per second}

\item{pwidth}{peak width for the compound. If it's a single value, simulated peaks' width will use this number as peak width for all compounds. If it's a numeric vector, it will be used as the peak widths for every compounds.}

\item{pheight}{peak height for the compound. If it's a single value, simulated peaks' height will use this number as peak height for all compounds. If it's a numeric vector, it will be used as the peak heights for every compounds.}

\item{baseline}{noise baseline, default 100}

\item{baselinesd}{standard deviation for noise, default 30}

\item{rf}{response factor of each compound, default 100 for all compounds}

\item{tailingfactor}{tailing factor for peaks, larger means larger tailing, default 1.2}

\item{compound}{numeric compound index in the database for targeted analysis, default NULL}

\item{rtime}{retention time for the compounds, if NULL, retention time will be simulated. Default NULL}

\item{tailingindex}{numeric index for tailing compounds, if NULL, all peaks will tailing. Default NULL}

\item{seed}{Random seed for reproducibility}

\item{unique}{if TRUE, one compound will have one spectra. Default FALSE}

\item{matrix}{if TRUE, m/z from experimental data will be used for background m/z simulation.Default FALSE}

\item{matrixmz}{custom matrix m/z vector, default NULL and predefined list from serum blank will be used.}
}
\value{
one mzML file for simulated data and one csv file the simulated compounds with retention time, m/z and name
}
\description{
Generate simulated mzml data and compounds list
}
\examples{
data(monams1)
simmzml(db=monams1, name = 'test')
}
