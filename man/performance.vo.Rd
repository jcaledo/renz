% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance.R
\name{performance.vo}
\alias{performance.vo}
\title{Performance of Methods Based on Initial Velocities}
\usage{
performance.vo(Km, Vm, unit_S = 'mM', unit_v = 'au', error, sd)
}
\arguments{
\item{Km}{Michaelis contant.}

\item{Vm}{maximal velocity.}

\item{unit_S}{concentration unit.}

\item{unit_v}{time unit.}

\item{error}{it should be one among c('absolute', 'relative').}

\item{sd}{standard deviation of the error.}
}
\value{
Returns a dataframe containing the fitted parameters under different conditions.
}
\description{
Evaluates the performance of the methods to estimate Km and Vm based on initial velocities.
}
\seealso{
sE.progress(), int.MM()
}
\author{
Juan Carlos Aledo
}
