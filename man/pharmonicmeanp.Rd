\name{pharmonicmeanp}
\alias{pharmonicmeanp}
\title{
The (asymptotic) probability that a harmonic mean \emph{p}-value (the test statistic of the HMP procedure) has value less than or equal to x.
}
\description{
The recommended usage of the HMP procedure is via the function p.hmp which internally computes the harmonic mean \emph{p}-value test statistic and outputs an asymptotically exact \emph{p}-value for the test statistic. If the test statistic (x) has already been computed, this function calculates the asymptotically exact cumulative distribution function (Pr(X<x)), i.e. the \emph{p}-value).
}
\usage{
pharmonicmeanp(x, L, log=FALSE, lower.tail=TRUE)
}
\arguments{
  \item{x}{
The previously computed harmonic mean \emph{p}-value test statistic, for example using function hmp.stat.
}
  \item{L}{
The number of constituent \emph{p}-values used in calculating x.
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
pharmonicmeanp(x,length(p))
p.hmp(p)
p.mamml(1/p,2)
}
\keyword{ ~hmp }
