\name{hmp.stat}
\alias{hmp.stat}
\title{
Compute the harmonic mean \emph{p}-value
}
\description{
The harmonic mean \emph{p}-value (HMP) is defined as the inverse of the (possibly weighted) arithmetic mean of the inverse \emph{p}-values. When the HMP is small (e.g. less than 0.05), it is approximately well-calibrated, meaning that it can be directly interpreted. However, the function p.hmp calculates an exact \emph{p}-value from the HMP and is preferred.
}
\usage{
hmp.stat(p, w = NULL)
}
\arguments{
  \item{p}{
A numeric vector of one or more \emph{p}-values. Missing values (NAs) will cause a missing value to be returned.
}
  \item{w}{
An optional numeric vector of weights that can be interpreted as prior model probabilities for each of the alternative hypotheses represented by the individual \emph{p}-values. If specified, it should sum to one.
}
}
\value{
The harmonic mean \emph{p}-value is returned.
}
\references{
Daniel J. Wilson (2017) The harmonic mean \emph{p}-value and model averaging by mean maximum likelihood. In preparation.
}
\author{
Daniel J. Wilson
}

\seealso{
p.hmp
}
\examples{
p = rbeta(1000,1/1.5,1)
hmp.stat(p)
p.hmp(p)
}
\keyword{ ~hmp }
