#---
#title: Table 5: Torus case (second line execution)
#Fixed LISA functions
#author: "Jonatan A. Gonzalez"
#---
# Please load the file "LocalLfunctions.R"
################################################################################
#Packages
library(spatstat)
library(spatstat.utils)
library(GET)
library(doParallel)
################################################################################

#Preliminary functions

#Simulation of a inhomogeneus thomas pattern with inhomogenoeus variance
rThomasInhomParents <- function(scale = 0.2, k0, mu, sd, beta = 0.2){
  W <- owin()
  frame <- boundingbox(W = owin())
  dilated <- grow.rectangle(frame, 4 * scale)
  GG <- attr(rLGCP(win = dilated, scale = beta), "Lambda")
  k1 <- GG * exp(k0)
  A <- rThomas(kappa = k1, scale = sd, mu = mu, win = W, expand = 4 * 0.2)
  return(list(pp = A, rf = k1[W]))
}

#Generating Thomas pp
tpp <- function(S = 0.2){
  rThomasInhomParents(sd = S, k0 = 2.70, mu = 10, beta = 0.2) #1.8 = 100, 2.7 = 250, 3.45 = 500
}

################################################################################
onesimu <- function(nsim = 99, S1 = 0.025, rmaxx = 0.15)
{
  PPP <- tpp(S = S1)
  pp <- PPP$pp 
  # Calculate LISA functions based on K
  lambda <- density.ppp(pp, at = "points")
  TemplateLisa <- localL1inhom(pp, lambda = lambda, verbose = F, 
                               correction = "translate", rmax = rmaxx)
  r0 <- floor(seq(1,512, length.out = 51))
  rr <- TemplateLisa$r[r0]
  KLisasObs <- as.matrix(TemplateLisa)
  KLisasObs <- KLisasObs[r0, 1:(dim(KLisasObs)[2] - 2)]
  #Covariate values
  Covariate <- PPP$rf
  CovariateObs <- Covariate[pp, drop = F]
  
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
system.time(onesimu(nsim = 99, S1 = 0.05, rmaxx = 0.15))

#Using parallel computing for acceletating things
#Note that we only use 99 simulations as illustration
nP <- function(s) mclapply(1:100, mc.cores = 14,
                           FUN = function(x) onesimu(nsim = 99, S1 = s, rmaxx = 0.15))

#Executing parallel procedure with the diferent scales of Thomas pp
P <- sapply(c(0.0125, 0.025, 0.05, 0.1), nP, simplify = "array")

#Estimating nominal significance (this case is the second line of Table 5)
nominal.rejection <- function(P) mean(unlist(P) <= 0.05)
apply(P, 2, nominal.rejection)
