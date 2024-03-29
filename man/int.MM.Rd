% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LintMM.R
\name{int.MM}
\alias{int.MM}
\title{Linearization of The Integrated Michaelis-Mente Equation}
\usage{
int.MM(data, unit_S = 'mM', unit_t = 'min')
}
\arguments{
\item{data}{a dataframe with two columns. The first column contains the values of the independent variable time, t, and the second column contains the substrate concentrations.}

\item{unit_S}{concentration unit.}

\item{unit_v}{time unit.}
}
\value{
A list of two elements. The first element is named vector containing the Km and Vm. The second element is a dataframe where the first two columns are the original data and the last two columns are the transformed variables. Also a linear plot the transformed variables together with the parameters values are provided.
}
\description{
Estimates the kinetic parameters using an linearized form of the integrated Michaelins-Menten equation.
}
\author{
Juan Carlos Aledo
}
