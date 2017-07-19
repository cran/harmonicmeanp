\name{p.hmp}
\alias{p.hmp}
\title{
Compute a combined \emph{p}-value via the harmonic mean \emph{p}-value
}
\description{
The harmonic mean \emph{p}-value (HMP) test combines \emph{p}-values and corrects for multiple testing while controlling the family-wise error rate in a way that is more powerful than common alternatives including Bonferroni and Simes procedures, more stringent than controlling the false discovery rate, and is robust to positive correlations between tests and unequal weights. It is a multi-level test in the sense that a superset of one or more significant tests is almost certain to be significant and conversely when the superset is non-significant, the constituent tests are almost certain to be non-significant. It is based on MAMML (model averaging by mean maximum likelihood), a frequentist analogue to Bayesian model averaging, and is theoretically grounded in generalized central limit theorem.
}
\usage{
p.hmp(p, w = NULL)
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
A combined \emph{p}-value is returned.
}
\references{
Daniel J. Wilson (2017) The harmonic mean \emph{p}-value and model averaging by mean maximum likelihood. In preparation.
}
\author{
Daniel J. Wilson
}

\seealso{
hmp.stat
}
\examples{
# Example: simulate from a non-uniform distribution mildly enriched for small \emph{p}-values. 
# Compare the significance of the combined p-value for Bonferroni, Benjamini-Hochberg (i.e. Simes), 
# HMP and (equivalently) MAMML with 2 degrees of freedom.
p = rbeta(1000,1/1.5,1)
min(p.adjust(p,"bonferroni"))
min(p.adjust(p,"BH"))
p.hmp(p)
p.mamml(1/p,2)
}
\keyword{ ~hmp }
