\name{obs}
\alias{obs}
\alias{giss}
\alias{aberdeen.t2m}
\alias{bismarck.t2m}
\alias{buenos.aires.t2m}
\alias{concord.t2m}
\alias{funchal.t2m}
\alias{honululu.t2m}
\alias{ishigakijima.t2m}
\alias{lisboa.t2m}
\alias{moskva.t2m}
\alias{nassau.t2m}
\alias{portland.t2m}
\alias{saentis.t2m}
\alias{sao.paulo.t2m}
\alias{seychelles.t2m}
\alias{thessaloniki.t2m}
\alias{thiruvanantha.t2m}
\alias{turuhansk.t2m}
\alias{cet}
\alias{GPC2004}
\title{Data sets.}
\description{
\code{obs} contains daily temperature and precipitation from Oslo. Stored as a 'daily.station.record' (clim.pact).

The other data sets are used for producing the plots in Benestad (2004)
Record-values, non-stationarity tests and extreme value distributions
Global and Planetary Change vol 44, issue 1-4, p.11-26 and Benestad,
R.E. (2008) 'A Simple Test for Changes in Statistical Distributions', Eos, 89 (41), 7 October 2008, p. 389-390: \code{demo(GPC2004a)}, \code{demo(GPC2004b)}, \& \code{demo(Eos2008)}:  
\code{giss} contains meta-data ofmonthly mean temperatures from
NASA/GISS (\url{http://data.giss.nasa.gov/gistemp/})  
\code{aberdeen.t2m}... \code{turuhansk.t2m}.

\code{cet} contains the central england temperature.

\code{GPC2004} is a list containing results from Monte-Carlo simulations
carried out by \code{demo(GPC2004a)} is \code{do.monte.carlo} is set to
\code{TRUE}. These data are used in \code{demo(GPC2004b)} for
plotting. The data are from 7 different Monte-Carlo simulations
(MC1...MC7), where MC1-3 are for random gamma-distributed data (\code{rgamma}) with
scale of 1, 10, and 0.1 respectively. MC4-6 are for GEV-distributed data
(\code{rgev}) with shape of 1, 10, and -1 respectively, while MC7 is for random
binomially distributed data (\code{rbinom}). The \code{GPC2004} contains:
\tabular{rlll}{
         [,1] \tab b.fit  \tab numeric \tab coefficients from LS-fit of
	 exp(actual number of records) against the length of series. \cr
         [,2-8] \tab b.fitX \tab numeric \tab coefficients from LS-fit of
	 exp(number of records from MC1-7, where X signifies 1-7) against the length of series. \cr
         [,9-10] \tab b.l.fit/b.u.fit  \tab numeric \tab coefficients from LS-fit of
	 exp(5/95-percentile of MC-simulated number of records) against the length of series. \cr
         [,11-17] \tab E.mc.X    \tab numeric \tab MC-simulated number
	 of record densities for the MC-simulations 1-7. \cr
         [,18] \tab NS   \tab numeric \tab Number of parrallel series \cr
       }
       

}
\usage{
obs
giss
aberdeen.t2m
bismarck.t2m
buenos.aires.t2m
concord.t2m
funchal.t2m
honululu.t2m
ishigakijima.t2m
lisboa.t2m
moskva.t2m
nassau.t2m
portland.t2m
saentis.t2m
saopaulo.t2m
seychelles.t2m
thessaloniki.t2m
thiruvanantha.t2m
turuhansk.t2m
cet
GPC2004
}

\format{The data is a 'daily.station.record' object (see \code{\link[clim.pact]{station.obj}}).}
\examples{
data(obs)
}
\source{
  http://data.giss.nasa.gov/gistemp/
  http://hadobs.metoffice.com/hadcet/
}
\references{
Benestad (2004) Record-values, non-stationarity tests and extreme value
distributions Global and Planetary Change vol 44, issue 1-4, p.11-26

Benestad,R.E. (2008) 'A Simple Test for Changes in Statistical Distributions', Eos, 89 (41), 7 October 2008, p. 389-390
}
\keyword{datasets}
