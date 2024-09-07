#---
#title: Table 1: 
#Fixed LISA functions
#author: "Jonatan A. Gonzalez"
#---
################################################################################
#Packages
library(spatstat)
library(GET)
library(doParallel)
################################################################################
#Generating random field observations
GG <- function() attr(rLGCP("exp", mu = 0, var = 1, scale = 0.05), "Lambda")

#We set n=100 replicates as an easy example
RF <- replicate(n = 100, simplify = FALSE, expr = GG()) 

#Generating Thomas pp
tpp <- function(S){
  rThomas(scale = S, kappa = 25, mu = 10)
}

################################################################################
#Generating one simulation function
onesimu <- function(GG, nsim = 99, S1 = 0.05, rmaxx = 0.15)
{
  pp <- tpp(S1)
  # Calculate LISA functions based on K
  TemplateLisa <- localLinhom(pp, verbose = F, rmax = rmaxx, correction = "translate")
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
  
  #Random Shiftings with torus
  randomloc <- function() {
    pp.shift <- rshift(pp, radius = 0.5, edge = "torus")
    CovariateSim <- safelookup(Covariate, pp.shift)
    return(Pearson(L = KLisasObs, Z = CovariateSim))
  }
  simu <- replicate(nsim, randomloc())
  CS <- create_curve_set(list(r = rr, obs = RhoObs, sim_m = simu))
  attr(rank_envelope(CS, type = "erl"), "p")
}

#Checking the time for one simulation
system.time(onesimu(RF[[1]], nsim = 999, S1 = 0.05, rmaxx = 0.15))

#Using parallel computing for acceletating things
#Note that we only use 99 simulations as illustration
nP <- function(s) mclapply(RF, mc.cores = 14,
                           FUN = function(x) onesimu(GG = x, nsim = 99, S1 = s, rmaxx = 0.15))

#Executing parallel procedure with the diferent scales of Thomas pp
P <- sapply(c(0.0125, 0.025, 0.05, 0.1), nP, simplify = "array")

#Estimating nominal significance (this case is the second line of Table 1)
nominal.rejection <- function(P) mean(unlist(P) <= 0.05)
apply(P, 2, nominal.rejection)
