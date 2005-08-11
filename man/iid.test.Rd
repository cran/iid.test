\name{iid.test}
\alias{iid.test}
\title{iid test}
\description{
Test for whether a variable is independent and identically distributed (iid). Used in \code{\link{daily.station.records}}.

Reference:

Benestad, R.E., 2003: How often can we expect a record-event? Climate Research. 23, 3-13 (pdf)

Benestad, R.E., 2004: Record values, nonstationarity tests and extreme value distributions, Global and Planetary Change, vol 44, p. 11-26

The papers are available in the pdf format from \url{http://regclim.met.no/results_iii_artref.html}.

Note, gaps of missing data (NA) can bias the results and produce an under-count. The sign of non-iid behaviour is when the 'forward' analysis indicated higher number of record-events than the confidence region and the backward analysis gives lower than the confidence region.
}
\usage{
iid.test(Y,plot=TRUE,Monte.Carlo=TRUE,N.test=200)
}
\arguments{
    \item{Y}{A data matrix or a vector.}
    \item{plot}{Flag: plot the diagnostics.}
    \item{Monte.Carlo}{Flag: for estimating confidence limits.}
    \item{N.test}{Number of Monre-Carlo runs.}
}
\value{list: 'record.density' and 'record.density.rev' for the reverse analysis.}
\author{R.E. Benestad}}
\examples{
dat <- rnorm(100*30)
dim(dat) <- c(100,30)
iid.test(dat)
}
\keyword{manip}

