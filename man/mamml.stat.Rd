\name{mamml.stat}
\alias{mamml.stat}
\title{
Compute the model-averaged mean maximized likelihood
}
\description{
The model-averaged mean maximized likelihood (MAMML) is defined as the (possibly weighted) arithmetic mean of the maximized likelihood ratios from a series of likelihood ratio tests comparing mutually exclusive alternative hypotheses with the same nested null hypothesis based on the exact same data.
}
\usage{
mamml.stat(R, w = NULL)
}
\arguments{
  \item{R}{
A numeric vector of one or more maximized likelihood ratios. Missing values (NAs) will cause a missing value to be returned.
}
  \item{w}{
An optional numeric vector of weights that can be interpreted as prior model probabilities for each of the alternative hypotheses represented by the individual \emph{p}-values. If specified, it should sum to one.
}
}
\value{
The model-averaged mean maximized likelihood ratio is returned.
}
\references{
Daniel J. Wilson (2017) The harmonic mean \emph{p}-value and model averaging by mean maximum likelihood. In preparation.
}
\author{
Daniel J. Wilson
}

\seealso{
p.mamml
}
\examples{
nu = 3
R = exp(0.5*rchisq(1000,nu))
mamml.stat(R)
p.mamml(R,nu)
}
\keyword{ ~MAMML }
