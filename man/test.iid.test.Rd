\name{test.iid.test}
\alias{test.iid.test}
\title{test.iid.test}
\description{
Testing routine for the iid-test.. 

Reference:

Benestad, R.E., 2003: How often can we expect a record-event? Climate Research. 23, 3-13 (pdf)

Benestad, R.E., 2004: Record values, nonstationarity tests and extreme value distributions, Global and Planetary Change, vol 44, p. 11-26

The papers are available in the pdf format from \url{http://regclim.met.no/results_iii_artref.html}.

}
\usage{
test.iid.test(distr="rnorm",d=c(70, 50),plot=TRUE,Monte.Carlo=TRUE)
}
\arguments{
    \item{distr}{What type of random numbers to be used, .eg normally distributed..}
    \item{d}{Dimension of the test matrix.}
    \item{plot}{For plotting the diagnostics.}
    \item{Monte.Carlo}{Flag for running Monte-Carlo simulations to estimate 95\% C.I.}
}
\value{list: 'record.density' and 'record.density.rev' for the reverse analysis.}
\author{R.E. Benestad}
\examples{
test.iid.test()
}
\keyword{manip}

