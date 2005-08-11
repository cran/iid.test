\name{daily.station.records}
\alias{daily.station.records}
\title{daily.station.records}
\description{
Test for whether daily station data are independent and identically distributed (iid). 

Reference:

Benestad, R.E., 2003: How often can we expect a record-event? Climate Research. 23, 3-13 (pdf)

Benestad, R.E., 2004: Record values, nonstationarity tests and extreme value distributions, Global and Planetary Change, vol 44, p. 11-26

The papers are available in the pdf format from \url{http://regclim.met.no/results_iii_artref.html}.

}
\usage{
daily.station.records(obs,element="precip",subsample=5,tolerance=2)
}
\arguments{
    \item{obs}{A daily.station.record (The clim.pact-package).}
    \item{element}{Decides which element is analysed ("t2m" or "precip").}
    \item{subsample}{Subsampling to reduce the effects of dependencies.}
    \item{tolerance}{How many days with missing value is tolerable before the results are significantly biased.}
}
\value{list: 'record.density' and 'record.density.rev' for the reverse analysis.}
\author{R.E. Benestad}}
\examples{
data(obs)
daily.station.records(obs,subsample=30)
}
\keyword{manip}

