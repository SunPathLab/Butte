% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnmut.R
\name{scnaInput}
\alias{scnaInput}
\title{Reading and sorting scnaFile}
\usage{
scnaInput(scnaFile, skipchunk = 19)
}
\arguments{
\item{scnaFile}{the SCNA segmentation file}

\item{skipchunk}{segments with number of data points (probes) no more than this number will be skipped, set it to zero if not skipped}
}
\value{
sorted scna segmentation data frame
}
\description{
Reading and sorting scnaFile
}
