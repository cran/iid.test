\name{n.records}
\alias{n.records}
\title{n.records}
\description{
Number of record-events. 

Reference:

Benestad, R.E., 2003: How often can we expect a record-event? Climate Research. 23, 3-13 (pdf)

Benestad, R.E., 2004: Record values, nonstationarity tests and extreme value distributions, Global and Planetary Change, vol 44, p. 11-26

The papers are available in the pdf format from \url{http://regclim.met.no/results_iii_artref.html}.

}
\usage{
n.records(y)
}
\arguments{
    \item{y}{A vector.}
}
\value{list: N= number of records, t=vector of time index corresponding to when a new record is set, events= a vector with TRUE where there is a record-event. The same fields with a suffix '.rev' represent the same analysis but for a reversed chronological order.} 
\author{R.E. Benestad}}
\examples{
y <- rnorm(100)
n.records(y)
}
\keyword{manip}

