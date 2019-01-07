\name{harmonicmeanp}
\alias{harmonicmeanp}
\alias{hmp}
\docType{package}
\title{
Harmonic mean p-values and model averaging by mean maximum likelihood
}
\description{
The harmonic mean \emph{p}-value (HMP) test combines \emph{p}-values and corrects for multiple testing while controlling the family-wise error rate in a way that is more powerful than common alternatives including Bonferroni and Simes procedures, more stringent than controlling the false discovery rate, and is robust to positive correlations between tests and unequal weights. It is a multi-level test in the sense that a superset of one or more significant tests is almost certain to be significant and conversely when the superset is non-significant, the constituent tests are almost certain to be non-significant. It is based on MAMML (model averaging by mean maximum likelihood), a frequentist analogue to Bayesian model averaging, and is theoretically grounded in generalized central limit theorem.
}
\details{
\tabular{ll}{
Package: \tab harmonicmeanp\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2019-01-05\cr
License: \tab GPL-3\cr
}
The key function is p.hmp for combining p-values using the HMP. Type vignette("harmonicmeanp") for detailed examples.
}
\author{
Daniel J. Wilson

Maintainer: Daniel Wilson <hmp.R.package@gmail.com>
}
\references{
Daniel J. Wilson (2019) The harmonic mean \emph{p}-value for combining dependent tests. \emph{Proceedings of the National Academy of Sciences USA}, in press.
}
\seealso{
Package FMStable
}
\examples{
# For detailed examples type vignette("harmonicmeanp")
# Example: simulate from a non-uniform distribution mildly enriched for small \emph{p}-values. 
# Compare the significance of the combined p-value for Bonferroni, Benjamini-Hochberg (i.e. Simes), 
# HMP and (equivalently) MAMML with 2 degrees of freedom.
p = rbeta(1000,1/1.5,1)
min(p.adjust(p,"bonferroni"))
min(p.adjust(p,"BH"))
p.hmp(p)
p.mamml(1/p,2)
}
