library(FMStable)
hmp.stat = function(p, w = NULL) {
  p = as.numeric(p)
  if(is.null(w)) return(c("hmp.stat"=1/mean(1/p)))
  if(abs(sum(w)-1)>1e-6) warning("Weights do not sum to one, renormalizing")
  w = w/sum(w)
  return(c("hmp.stat"=1/sum(w/p)))
}
p.hmp = function(p, w = NULL) {
  HMP = hmp.stat(p, w)
  O.874 = 1+digamma(1)-log(2/pi)
  return(c("p.hmp"=pEstable(1/HMP,setParam(alpha=1,location=(log(length(p))+O.874),logscale=log(pi/2),pm=0),lower.tail=FALSE)))
}
mamml.stat = function(R, w = NULL) {
  R = as.numeric(R)
  if(any(R<1)) stop("Maximized likelihood ratios cannot be less than one")
  if(is.null(w)) return(c("mamml.stat"=mean(R)))
  if(abs(sum(w)-1)>1e-6) warning("Weights do not sum to one, renormalizing")
  w = w/sum(w)
  return(c("mamml.stat"=sum(w*R)))
}
p.mamml = function(R, nu, w = NULL) {
  if(length(nu)!=1 & length(nu)!=length(R)) stop("Degrees of freedom (nu) must have length one or length of R")
  Rbar = mamml.stat(R, w)
  nu = as.numeric(nu)
  if(any(nu<=0)) stop("Degrees of freedom (nu) must be positive")
  nu.max = max(nu)
  if(nu.max<2) {
    c = pgamma(log(Rbar),nu.max/2,1,lower.tail=FALSE)*Rbar
  } else {
    c = pgamma(log(length(R)*Rbar),nu.max/2,1,lower.tail=FALSE)*length(R)*Rbar
  }
  O.874 = 1+digamma(1)-log(2/pi)
  return(c("p.mamml"=pEstable(Rbar,setParam(alpha=1,location=c*(log(length(R))+O.874),logscale=log(pi/2*c),pm=0),lower.tail=FALSE)))
}
