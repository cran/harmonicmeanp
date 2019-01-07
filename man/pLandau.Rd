\name{pLandau}
\alias{pLandau}
\title{
The cumulative distribution function of a Landau distribution.
}
\description{
If X has a Landau(mu,sigma) distribution, the function produces the probability Pr(X<x) when lower.tail=TRUE or Pr(X>x) when lower.tail=FALSE. Read the cited paper for the parameterization of the distribution.
}
\usage{
pLandau(x,mu=log(pi/2),sigma=pi/2,log=FALSE,lower.tail=TRUE)
}
\arguments{
  \item{x}{
The value of a Landau-distributed random variable.
}
  \item{mu}{
The location parameter of the Landau distribution. Defaults to log(pi/2) to give Landau's original distribution.
}
  \item{sigma}{
The scale parameter of the Landau distribution. Defaults to pi/2 to give Landau's original distribution.
}
  \item{log}{
If true the log probability is output.
}
  \item{lower.tail}{
If true (the default) the lower tail probability is returned. Otherwise the upper tail probability.
}
}
\value{
A tail probability is returned.
}
\references{
Daniel J. Wilson (2019) The harmonic mean \emph{p}-value for combining dependent tests. \emph{Proceedings of the National Academy of Sciences USA}, in press.
}
\author{
Daniel J. Wilson
}

\seealso{
p.hmp
}
\examples{
# For detailed examples type vignette("harmonicmeanp")
# Example: simulate from a non-uniform distribution mildly enriched for small \emph{p}-values. 
# Compare the significance of the combined p-value for Bonferroni, Benjamini-Hochberg (i.e. Simes), 
# HMP and (equivalently) MAMML with 2 degrees of freedom.
p = rbeta(1000,1/1.5,1)
min(p.adjust(p,"bonferroni"))
min(p.adjust(p,"BH"))
x = hmp.stat(p)
pLandau(1/x,log(length(p))+(1 + digamma(1) - log(2/pi)),pi/2,lower.tail=FALSE)
p.hmp(p)
p.mamml(1/p,2)
}
\keyword{ ~hmp }
