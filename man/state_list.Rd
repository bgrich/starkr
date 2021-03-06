% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/state_list.R
\name{state_list}
\alias{state_list}
\title{List of States for Stark Matrix.}
\usage{
state_list(nmin, nmax, mj, n_add_min = 7, n_add_max = 7)
}
\arguments{
\item{nmin}{A numeric. The minimum principle quantum number for the whole
manifold}

\item{nmax}{A numeric. THe maximum principle quantum number for the whole
manifold}

\item{mj}{A numeric. The magnetic momentum quantum number.}

\item{n_add_min}{A numeric. The minimum principle quantum number for low
angular momentum states.}

\item{n_add_max}{A numeric. The maximum principle quantum number for low
angular momentum states.}
}
\description{
\code{state_list} creates a matrix containing the states to be calculated in
the Stark matrix calculation.
}
\details{
This function outputs a matrix with the n, l, and j states for a given m_j
and range of n's. It gives teh full manifold between nmin and max. It also
returns the low angular momentum states (l < 6) for a range of nmin to
n_add_min and nmax to n_add_max. If n_add_min or n_add_max are <= 0, then
those sections will be ignored and no extra states will be added.

The output takes the form of a matrix with columns n, l, and j.
}
