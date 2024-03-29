% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnmutHistoryMatrix.R
\name{cnmutHistory}
\alias{cnmutHistory}
\title{generate history matrix in relating SCNA to SSNVs}
\usage{
cnmutHistory(nt, nb)
}
\arguments{
\item{nt}{total copy number}

\item{nb}{copy number of the minor allele}
}
\value{
A list of possible matrices
}
\description{
Given a SCNA configuration Nt (total copy) and Nb (minor copy)
this function produces the possible history matrices in relating 
CN timing (time period for each stage) to the burden of SSNVs
at distinct allele states.
}
