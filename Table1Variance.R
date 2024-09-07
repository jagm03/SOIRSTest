#---
#title: "third paper simulations" Variance Correctio, FIxed LISA table 1
#author: "Jonatan A. Gonzalez "
#---
################################################################################
#Packages
library(spatstat)
library(GET)
library(doParallel)
################################################################################
#Preliminary functions
Rshift.J <- function(Z, radius = 0.5){
  jump <- runifdisc(1, radius = radius)
  X <- shift(Z, jump)
  Wf <- intersect.owin(Z$window, X$window)
  Xok <- inside.owin(X$x, X$y, Wf)
  return(ppp(x = X$x[Xok], y = X$y[Xok], marks = X$marks[Xok],window = Wf))
}

#Generating random field observations
GG <- function() attr(rLGCP("exp", mu = 0, var = 1, scale = 0.05), "Lambda")

#We set n=100 replicates as an easy example
RF <- replicate(n = 100, simplify = FALSE, expr = GG()) 

tpp <- function(S){
  rThomas(scale = S, kappa = 25, mu = 10)
}

################################################################################
onesimu <- function(GG, nsim = 99, S1 = 0.01, rmaxx = 0.15)
{
  #simulate a log-Gaussian with scale  =0.1
  pp <- tpp(S1)
  # Calculate LISA functions based on K
  TemplateLisa <- localL(pp, verbose = F, rmax = rmaxx, correction = "translate")
  r0 <- floor(seq(1,513, length.out = 51))
  rr <- TemplateLisa$r[r0]
  KLisasObs <- as.matrix(TemplateLisa)
  KLisasObs <- KLisasObs[r0, 1:(dim(KLisasObs)[2] - 2)]
  #Covariate values
  Covariate <- GG
  CovariateObs <- safelookup(Covariate, pp)
  
  Pearson <- function(L, Z){
    Pear <- apply(L, 1, cor, y = Z, method = "pearson")
    Pear[is.na(Pear)] <- 0
    return(Pear)
  }
  
  RhoObs <- Pearson(L = KLisasObs, Z = CovariateObs)
  
  pp <- pp %mark% 1:pp$n
  #Random Shiftings with variance
  randomloc <- function() {
    repeat {
      pp.shift <- Rshift.J(pp, radius = 0.5)
      nn <- npoints(pp.shift)
      if (nn > 5) break
    }
    Lisa.shift <- KLisasObs[, pp.shift$marks]
    CovariateSim <- safelookup(Covariate, pp.shift)
    corr <- Pearson(L = Lisa.shift, Z = CovariateSim)
    return(list(Corr = corr, N = nn))
  }
  simu <- replicate(nsim, randomloc())
  Rhosimu <- sapply(simu[1, ], "[")
  Nsimu <- sapply(simu[2, ], "[")
  
  #No correction but weighted variance sqrt(n) (Ti -T)
  Rhosimu <- cbind(Rhosimu, RhoObs)
  Rhomean <- apply(Rhosimu, 1, mean)
  Ti <- sweep(Rhosimu, 1, Rhomean, FUN = "-")
  Si <- c(Nsimu, npoints(pp))
  TT <- sweep(Ti, 2, sqrt(Si), FUN = "*")
  ####
  
  CS <- create_curve_set(list(r = rr, obs = TT[, nsim + 1], sim_m = TT[, 1:nsim]))
  attr(rank_envelope(CS, type = "erl"), "p")
}

system.time(A <- onesimu(RFn[[1]], nsim = 999, S1 = 0.05, rmaxx = 0.15))

nP <- function(s) mclapply(RFn[sample(1:4000, 1000)], mc.cores = 60,
                           FUN = function(x) onesimu(GG = x, nsim = 999, S1 = s, rmaxx = 0.15))

P <- sapply(c(0.0125, 0.025, 0.05, 0.1), nP, simplify = "array")

nominal.rejection <- function(P) mean(unlist(P) <= 0.05)
apply(P, 2, nominal.rejection)

