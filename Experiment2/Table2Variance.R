#---
#title: Table 2: Variance case (second line execution)
#Fixed LISA functions
#author: "Jonatan A. Gonzalez"
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
GG <- function() attr(rLGCP("exp", mu = 0, var = 0.05, scale = 0.2), "Lambda")

#We set n=100 replicates as an easy example
RF <- replicate(n = 100, simplify = FALSE, expr = GG()) 

#Generating Thomas pp
tpp <- function(S){
  rThomas(scale = S, kappa = 20, mu = 10)
}

################################################################################
onesimu <- function(GG, nsim = 99, S1 = 0.05, rmaxx = 0.15)
{
  pp0 <- tpp(S1)
  tG <- GG / max(GG) #Building a probability surface
  pp1 <- safelookup(tG, pp0)
  pp <- rthin(X = pp0, P = pp1) #A ppp still SOIRS
  #
  #Intens <- 200 * tG (TRUE) #True Intensity
  
  Intens <- predict(ppm(pp ~ Z, covariates = list(Z = tG)), type = "intensity")
  #Parametric estimate
  
  #Intens <- density.ppp(pp, at = "points") #Non-parametric estimate
  
  # Calculate LISA functions based on K
  TemplateLisa <- localLinhom(pp, lambda = Intens,
                              verbose = F, rmax = rmaxx, correction = "translate")
  r0 <- floor(seq(1,513, length.out = 51))
  rr <- TemplateLisa$r[r0]
  KLisasObs <- as.matrix(TemplateLisa)
  KLisasObs <- KLisasObs[r0, 1:(dim(KLisasObs)[2] - 2)]
  #Covariate values
  Covariate <- tG
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

#Checking the time for one simulation
system.time(onesimu(RF[[1]], nsim = 99, S1 = 0.05, rmaxx = 0.15))

#Using parallel computing for acceletating things
#Note that we only use 99 simulations as illustration
nP <- function(s) mclapply(RF, mc.cores = 14,
                           FUN = function(x) onesimu(GG = x, nsim = 99, S1 = s, rmaxx = 0.15))

#Executing parallel procedure with the diferent scales of Thomas pp
P <- sapply(c(0.0125, 0.025, 0.05, 0.1), nP, simplify = "array")

#Estimating nominal significance (this case is the fifth line of Table 2)
nominal.rejection <- function(P) mean(unlist(P) <= 0.05)
apply(P, 2, nominal.rejection)

