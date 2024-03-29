## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(strip.white=FALSE,class.output="Routput",class.source="Rsource")

## ----installation, exercise=TRUE, eval=FALSE----------------------------------
#  install.packages("harmonicmeanp")

## ----require, exercise=TRUE---------------------------------------------------
library(harmonicmeanp)

## ----checkversion, exercise=TRUE, eval=FALSE----------------------------------
#  stopifnot(packageVersion("harmonicmeanp")>=3.0)

## ----download, exercise=TRUE--------------------------------------------------
system.time((gwas = read.delim("https://www.danielwilson.me.uk/files/Neuroticism_ch12.txt",
  header=TRUE,as.is=TRUE)))
head(gwas)

## ----HMP, exercise=TRUE-------------------------------------------------------
L = 6524432
gwas$w = 1/L
R = 1:nrow(gwas)
(HMP.R = sum(gwas$w[R])/sum(gwas$w[R]/gwas$p[R]))

## ----hmpthreshold, exercise=TRUE----------------------------------------------
# Specify the false positive rate
alpha = 0.05
# Compute the HMP significance threshold
(alpha.L = qharmonicmeanp(alpha, L))

## ----hmpthresholdadjust, exercise=TRUE----------------------------------------
# Test whether the HMP for subset R is significance
w.R = sum(gwas$w[R])
alpha.L * w.R

## ----phmp, exercise=TRUE------------------------------------------------------
# Use p.hmp instead to compute the HMP test statistic and
# calculate its asymptotically exact p-value in one step
# Note this line has changed because of a previous error.
w.R*pharmonicmeanp(HMP.R/w.R, L=L, lower.tail=TRUE)
# Compare it to the multiple testing threshold
w.R*alpha

## ----p.hmp, exercise=TRUE-----------------------------------------------------
# Note that the p.hmp function has been redefined to take argument L. Omitting L will issue a warning.
R = 1:nrow(gwas)
p.hmp(gwas$p[R],gwas$w[R],L)

## ----oddsevens, exercise=TRUE-------------------------------------------------
R = which(gwas$pos%%2==0)
p.hmp(gwas$p[R],gwas$w[R],L)
w.R = sum(gwas$w[R])
alpha*w.R
R = which(gwas$pos%%2==1)
p.hmp(gwas$p[R],gwas$w[R],L)
w.R = sum(gwas$w[R])
alpha*w.R

## ----oddsevens.adjust, exercise=TRUE------------------------------------------
R = which(gwas$pos%%2==0)
p.R = p.hmp(gwas$p[R],gwas$w[R],L)
w.R = sum(gwas$w[R])
(p.R.adjust = p.R/w.R)
R = which(gwas$pos%%2==1)
p.R = p.hmp(gwas$p[R],gwas$w[R],L)
w.R = sum(gwas$w[R])
(p.R.adjust = p.R/w.R)

## ----twohalves, exercise=TRUE-------------------------------------------------
R = 1:156229
p.R = p.hmp(gwas$p[R],gwas$w[R],L)
w.R = sum(gwas$w[R])
(p.R.adjust = p.R/w.R)
R = 156230:312457
p.R = p.hmp(gwas$p[R],gwas$w[R],L)
w.R = sum(gwas$w[R])
(p.R.adjust = p.R/w.R)

## ----win, exercise=TRUE-------------------------------------------------------
# Define overlapping sliding windows of 50 megabase at 10 megabase intervals
win.50M.beg = outer(0:floor(max(gwas$pos/50e6-1)),(0:4)/5,"+")*50e6
win.50M.beg = win.50M.beg[win.50M.beg+50e6<=max(gwas$pos)]
# Calculate the combined p-values for each window
system.time({
  p.50M = sapply(win.50M.beg,function(beg) {
    R = which(gwas$pos>=beg & gwas$pos<(beg+50e6))
    p.hmp(gwas$p[R],gwas$w[R],L)
  })
})
# Calculate sums of weights for each combined test
system.time({
  w.50M = sapply(win.50M.beg,function(beg) {
    R = which(gwas$pos>=beg & gwas$pos<(beg+50e6))
    sum(gwas$w[R])
  })
})
# Calculate adjusted p-value for each window
p.50M.adj = p.50M/w.50M 

## ----winplot, exercise=TRUE, fig.width=7, fig.height=5------------------------
# Took a few seconds, plotting over 312k points
gwas$p.adj = gwas$p/gwas$w
plot(gwas$pos/1e6,-log10(gwas$p.adj),pch=".",xlab="Position on chromosome 12 (megabases)",
  ylab="Adjusted significance (-log10 adjusted p-value)",
  ylim=sort(-log10(range(gwas$p.adj,p.50M.adj,na.rm=TRUE)))) 
arrows(win.50M.beg/1e6,-log10(p.50M.adj),(win.50M.beg+50e6)/1e6,len=0,col="#D3C991",lwd=2)
# Superimpose the significance threshold, alpha, e.g. alpha=0.05
abline(h=-log10(0.05),col="black",lty=2)
# When using the HMP to evaluate individual p-values, the HMP threshold must be used,
# which is slightly more stringent than Bonferroni for individual tests
abline(h=-log10(qharmonicmeanp(0.05,L)),col="grey",lty=2)
# For comparison, plot the conventional GWAS threshold of 5e-8. Need to convert
# this into the adjusted p-value scale. Instead of comparing each raw p-value
# against a Bonferonni threshold of alpha/L=0.05/6524432, we would be comparing each
# against 5e-8. So the adjusted p-values p/w=p*L would be compared against
# 5e-8*L = 5e-8 * 6524432 = 0.3262216
abline(h=-log10(0.3262216),col="grey",lty=3) 

## ----listpos, exercise=TRUE---------------------------------------------------
win.50M.beg[which(p.50M.adj<=0.05)]
# Also list the position of the most significant individual (adjusted) p-value
(peakpos = gwas$pos[gwas$p.adj==min(gwas$p.adj)])

## ----winlengths, exercise=TRUE------------------------------------------------
# Window of 100 base pairs
wlen = 100
R = which(abs(gwas$pos-peakpos)<wlen)
(p.R.adjust = p.hmp(gwas$p[R],gwas$w[R],L)/sum(gwas$w[R]))
# Window of 1 kilobase
wlen = 1e3
R = which(abs(gwas$pos-peakpos)<wlen)
(p.R.adjust = p.hmp(gwas$p[R],gwas$w[R],L)/sum(gwas$w[R]))
# Window of 10 kilobases
wlen = 1e4
R = which(abs(gwas$pos-peakpos)<wlen)
(p.R.adjust = p.hmp(gwas$p[R],gwas$w[R],L)/sum(gwas$w[R]))
# Window of 100 kilobases
wlen = 1e5
R = which(abs(gwas$pos-peakpos)<wlen)
(p.R.adjust = p.hmp(gwas$p[R],gwas$w[R],L)/sum(gwas$w[R]))
# Window of 1 megabase
wlen = 1e6
R = which(abs(gwas$pos-peakpos)<wlen)
(p.R.adjust = p.hmp(gwas$p[R],gwas$w[R],L)/sum(gwas$w[R]))
# Window of 10 megabases
wlen = 1e7
R = which(abs(gwas$pos-peakpos)<wlen)
(p.R.adjust = p.hmp(gwas$p[R],gwas$w[R],L)/sum(gwas$w[R]))
# Window of 20 megabases
wlen = 2e7
R = which(abs(gwas$pos-peakpos)<wlen)
(p.R.adjust = p.hmp(gwas$p[R],gwas$w[R],L)/sum(gwas$w[R]))
# Window of 30 megabases
wlen = 3e7
R = which(abs(gwas$pos-peakpos)<wlen)
(p.R.adjust = p.hmp(gwas$p[R],gwas$w[R],L)/sum(gwas$w[R]))
# Window of 40 megabases
wlen = 4e7
R = which(abs(gwas$pos-peakpos)<wlen)
(p.R.adjust = p.hmp(gwas$p[R],gwas$w[R],L)/sum(gwas$w[R]))
# Window of 50 megabases
wlen = 5e7
R = which(abs(gwas$pos-peakpos)<wlen)
(p.R.adjust = p.hmp(gwas$p[R],gwas$w[R],L)/sum(gwas$w[R]))

## ----optwin, exercise=TRUE----------------------------------------------------
# Find the smallest window centred on position 118876918 significant at alpha=0.05
f = function(wlen) {
  R = which(abs(gwas$pos-peakpos)<wlen)
  p.R.adjust = p.hmp(gwas$p[R],gwas$w[R],L)/sum(gwas$w[R])
  return(p.R.adjust - 0.05)
}
(wlen.opt = uniroot(f,c(3e7,4e7))$root)
# Show that the group of SNPs in this window is indeed significant
wlen = wlen.opt
R = which(abs(gwas$pos-peakpos)<wlen)
(p.R.adjust = p.hmp(gwas$p[R],gwas$w[R],L)/sum(gwas$w[R]))
# The number of individual SNPs included in this group
length(R)

## ----fiddler, exercise=TRUE, fig.width=4, fig.height=4------------------------
# Load the ape package for reading and plotting the tree
library(ape)
tree = read.tree(text=
"(((chlorophthalmus:1,crassipes:1)A:1,(inversa:1,sindensis:1)B:1):1,argillicola:3);")
plot(tree, show.node.label=TRUE)

log.carapace.breadth = c("chlorophthalmus"=1.02,"crassipes"=1.06,"inversa"=0.96,
"sindensis"=0.92,"argillicola"=0.89)
log.propodus.length = c("chlorophthalmus"=1.38,"crassipes"=1.41,"inversa"=1.36,
"sindensis"=1.22,"argillicola"=1.13)

plot(log.propodus.length ~ log.carapace.breadth)

## ----dataframe, exercise=TRUE-------------------------------------------------
# Convert branches in the tree into informative 'partitions'
informative.partitions = function(tree) {
  n = length(tree$tip.label)
  m = sapply(n+1:tree$Nnode,function(node) {
    1*is.na(match(tree$tip.label,extract.clade(tree,node)$tip.label))
  })
  rownames(m) = tree$tip.label
  colnames(m) = paste0("node.",tree$node.label)
  cs = colSums(m)
  is.informative = pmin(cs,n-cs)>1
  m[,is.informative]
}
# Extract phylogenetically informative partitions from the tree
partition = informative.partitions(tree)

# Create a data frame combining all the information
Uca = data.frame(log.propodus.length,log.carapace.breadth,partition)

## ----models, exercise=TRUE----------------------------------------------------
# Claw size does not vary by species
m0 = formula(log.propodus.length ~ 1) # grand null
# Claw size is associated with body size and there is no phylogenetic correlation
m1 = formula(log.propodus.length ~ log.carapace.breadth)
# Claw size isn't associated with body size but it is different in the descendents of ancestor A
m2 = formula(log.propodus.length ~ node.A)
# Claw size isn't associated with body size but it is different in the descendents of ancestor B
m3 = formula(log.propodus.length ~ node.B)
# Claw size is associated with body size and it is different in the descendants of ancestor A
m4 = formula(log.propodus.length ~ log.carapace.breadth + node.A)
# Claw size is associated with body size and it is different in the descendants of ancestor B
m5 = formula(log.propodus.length ~ log.carapace.breadth + node.B)
# Claw size isn't associated with body size but is different in descendents of ancestors A & B
m6 = formula(log.propodus.length ~ node.A + node.B)
# Claw size is associated with body size and is different in descendants of ancestors A & B
m7 = formula(log.propodus.length ~ log.carapace.breadth + node.A + node.B) # grand alternative
# List the alternatives together
mA = list(m1,m2,m3,m4,m5,m6,m7)

## ----pairwisetests, exercise=TRUE---------------------------------------------
# Output p-values from all tests for the inclusion of the primary regressor
pairwise.p = function(response,primary,data) {
  # Define a model space including the grand null
  rid = which(colnames(data)==response)
  if(length(rid)!=1) stop("Could not find response variable")
  # Define the 'primary' regressor
  pid = which(colnames(data)==primary)
  if(length(pid)!=1) stop("Could not find primary regressor")
  # Define the 'secondary' regressors, excluding the response and 'primary' regressor
  xid = (1:ncol(data))[-c(rid,pid)]
  if(length(xid)<1) stop("Could find only the primary regressor")
  # Create a table of every unique combination of models involving the secondary regressors
  delta = expand.grid(lapply(xid,function(j) 0:1))
  colnames(delta) = colnames(data)[xid]
  # Sort them by the number of regressors included, from fewest to most
  delta = delta[order(rowSums(delta)),]
  # Enumerate the models, adding the primary regressor to every one
  mpairs = apply(delta,1,function(x) {
    if(all(x==0)) {
      formula(paste0(colnames(data)[rid],"~",colnames(data)[pid]))
    } else {
      formula(paste0(colnames(data)[rid],"~",colnames(data)[pid],"+",
        paste(colnames(data)[xid][x==1],collapse="+")))
    }
  })
  names(mpairs) = gsub(colnames(data)[pid],paste0("[",colnames(data)[pid],"]"),
    as.character(mpairs),perl=TRUE)
  # Calculate a p-value for the inclusion of the primary regressor in each model
  lapply(mpairs,function(m) {
    fit = lm(m, data=data)
    drop1(fit,colnames(data)[pid],test="Chisq")[colnames(data)[pid],"Pr(>Chi)"]
  })
}
# Calculate the p-values from all tests for the inclusion of log.carapace.breadth
(p = pairwise.p(response="log.propodus.length",primary="log.carapace.breadth",data=Uca))

## ----hmpbodysize, exercise=TRUE-----------------------------------------------
# Specify the weight of each test, assuming equal weights
L = 12
(w = rep(1/L,length(p)))
# Calculate the model-averaged (asymptotically exact) HMP
(p.comb = p.hmp(p,w,L))
# Sum the weights of the constituent tests
(w.comb = sum(w))
# Calculate an adjusted model-averaged p-value for comparison to the ssFWER alpha
(p.comb.adj = p.comb/w.comb)

## ----bonferronibodysize, exercise=TRUE----------------------------------------
(p.HMP.adj = Vectorize(p.hmp)(p,w,L)/w)
(p.Bonf.adj = unlist(p)/w)
(p.Bonf = min(p.Bonf.adj))

## ----hmpnodeAB, exercise=TRUE-------------------------------------------------
# Is there a significant difference in claw size between the descendants of ancestor A
# and other species?
p = pairwise.p(response="log.propodus.length",primary="node.A",data=Uca)
w = rep(1/L,length(p))
p.hmp(p,w,L)/sum(w)
# Individual tests: HMP and Bonferroni
(p.hmp.adj = Vectorize(p.hmp)(p,w,L)/w)
(p.adj = unlist(p)/w)
(p.Bonf = min(p.adj))
# Is there a significant difference in claw size between the descendants of ancestor B
# and other species?
p = pairwise.p(response="log.propodus.length",primary="node.B",data=Uca)
w = rep(1/L,length(p))
p.hmp(p,w,L)/sum(w)
# Individual tests: HMP and Bonferroni
(p.hmp.adj = Vectorize(p.hmp)(p,w,L)/w)
(p.adj = unlist(p)/w)
(p.Bonf = min(p.adj))

## ----ppairwise, exercise=TRUE-------------------------------------------------
p = c(
  unlist(pairwise.p(response="log.propodus.length",primary="log.carapace.breadth",data=Uca)),
  unlist(pairwise.p(response="log.propodus.length",primary="node.A",data=Uca)),
  unlist(pairwise.p(response="log.propodus.length",primary="node.B",data=Uca)))

## ----terms, exercise=TRUE-----------------------------------------------------
terms = lapply(names(p),function(s) labels(terms(as.formula(gsub("\\[|\\]","",s,perl=TRUE)))))
(nterms = unlist(lapply(terms,length)))
(s = ncol(Uca)-1)
(m = 1/s)
(mu = m^nterms * (1-m)^(s-nterms))

## ----test.term, exercise=TRUE-------------------------------------------------
test.term = sapply(names(p),function(s) 
  gsub("\\[|\\]","",regmatches(s,regexpr("\\[.*?\\]",s,perl=TRUE)),perl=TRUE)
)
Uca.var = apply(Uca,2,var)
ssqb.over.ssqe = 2/Uca.var[test.term]
names(ssqb.over.ssqe) = names(p)
Var.beta.over.ssqe = sapply(names(p),function(s) {
  test.term = gsub("\\[|\\]","",
    regmatches(s,regexpr("\\[.*?\\]",s,perl=TRUE)),perl=TRUE)
  X = model.matrix(as.formula(gsub("\\[|\\]","",s,perl=TRUE)), data=Uca)
  solve(crossprod(X))[test.term,test.term]
})

# When the beta approximation performs poorly, best to evaluate
# near the likely value of the final threshold
smallp = 0.05/length(p)
# These are the test powers assuming threshold smallp
(smallp.pow = pchisq(qchisq(smallp,1,lower.tail=FALSE)/
  (1+ssqb.over.ssqe/Var.beta.over.ssqe),1,lower.tail=FALSE))
# Sanity checks
if(any(smallp.pow==1))
  warning("Perfect power test detected, check this is plausible")
if(any(smallp.pow<smallp))
  stop("Tests with worse power than smallp violate assumptions")
if(any(!is.finite(smallp.pow)) | any(smallp.pow<0) | any(smallp.pow>1))
  stop("Power cannot be outside range 0-1")
# Convert them into the parameter of the Beta(xi,1) distribution
xi = log(smallp.pow)/log(smallp)
if(any(!is.finite(xi)) | any(xi<0) | any(xi>1))
  stop("Beta(xi,1): xi cannot be outside range 0-1")

# Optimize the weights
wfunc = function(mu,xi,lambda,alpha) (mu*xi/lambda)^(1/(1-xi))/alpha
lambdafunc = function(lambda,alpha) sum(wfunc(mu,xi,lambda,alpha))-1
lambda = uniroot(lambdafunc,c(1e-6,1e6),0.05)$root
lambda = uniroot(lambdafunc,lambda*c(0.1,1.1),0.05)$root
(w = wfunc(mu,xi,lambda,0.05))
# Check the weights sum to one
sum(w)
if(abs(1-sum(w))>1e-4) stop("weights do not sum to one, check")

## ----wvsunw, exercise=TRUE----------------------------------------------------
# Compare the weighted and unweighted results
wequal = rep(1/length(w),length(w))
# For log.carapace.breadth
incl = which(test.term == "log.carapace.breadth")
c("weighted"=p.hmp(p[incl],w[incl],L)/sum(w[incl]),
  "unweighted"=p.hmp(p[incl],wequal[incl],L)/sum(wequal[incl]))
# For node.A
incl = which(test.term == "node.A")
c("weighted"=p.hmp(p[incl],w[incl],L)/sum(w[incl]),
  "unweighted"=p.hmp(p[incl],wequal[incl],L)/sum(wequal[incl]))
# For node.B
incl = which(test.term == "node.B")
c("weighted"=p.hmp(p[incl],w[incl],L)/sum(w[incl]),
  "unweighted"=p.hmp(p[incl],wequal[incl],L)/sum(wequal[incl]))
# Headline
incl = 1:length(p)
c("weighted"=p.hmp(p[incl],w[incl],L)/sum(w[incl]),
  "unweighted"=p.hmp(p[incl],wequal[incl],L)/sum(wequal[incl]))

## ----enumerateallmodels, exercise=TRUE----------------------------------------
enumerate.models = function(response,data) {
  # Define the response variable
  rid = which(colnames(data)==response)
  if(length(rid)!=1) stop("Could not find the response variable")
  # Define the regressors
  xid = (1:ncol(data))[-rid]
  # Create a table defining every unique combination of alternative hypotheses
  delta = expand.grid(lapply(xid,function(j) 0:1))
  colnames(delta) = colnames(data)[xid]
  # Sort them from fewest to most terms
  delta = delta[order(rowSums(delta)),]
  # Remove the grand null model
  delta = delta[rowSums(delta)>0,]
  # Define the grand null model separately
  m0 = formula(paste0(colnames(data)[rid],"~1"))
  # Define the alternative models
  mA = apply(delta,1,function(x) formula(paste0(colnames(data)[rid],"~",
    paste(colnames(data)[xid][x==1],collapse="+"))))
  names(mA) = as.character(mA)
  return(list("m0"=m0,"mA"=mA))
}
# E.g. on the Uca data
enumerate.models("log.propodus.length",Uca)

