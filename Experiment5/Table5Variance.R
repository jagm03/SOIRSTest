#---
#title: Table 5: Variance case (fifth line execution)
#Fixed LISA functions
#author: "Jonatan A. Gonzalez"
#---
# Please load the file "LocalLfunctions.R"
################################################################################
#Packages
library(spatstat)
library(GET)
library(RandomFields)
library(doParallel)
################################################################################
#Preliminary functions
Rshift.J <- function(Z, radius = 0.5){
  jump <- runifdisc(1, radius = radius)
  X <- shift(Z, jump)
  Wf <- intersect.owin(Z$window, X$window)
  Xok <- inside.owin(X$x, X$y, Wf)
  return(list(PP = ppp(x = X$x[Xok], y = X$y[Xok], window = Wf), OK = Xok))
}

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

tpp <- function(S = 0.2){
  rThomasInhomParents(sd = S, k0 = 2.7, mu = 10) #1.8 = 100, 2.7 = 250, 3.45 = 500
}

################################################################################
onesimu <- function(nsim = 99, S1 = 0.01, rmaxx = 0.15)
{
  repeat {
    PPP <- tpp(S = S1)
    pp <- PPP$pp
    if (npoints(pp) > 20) break
  }
  PPP <- tpp(S = S1)
  pp <- PPP$pp 
  # Calculate LISA functions based on K
  lambda <- density.ppp(pp, at = "points")
  #TemplateLisa <- localL(pp, verbose = F, rmax = rmaxx, correction = "translate")
  TemplateLisa <- localL1inhom(pp, lambda = lambda, verbose = F, 
                               correction = "translate", rmax = rmaxx)
  r0 <- floor(seq(1,512, length.out = 51))
  rr <- TemplateLisa$r[r0]
  KLisasObs <- as.matrix(TemplateLisa)
  KLisasObs <- KLisasObs[r0, 1:(dim(KLisasObs)[2] - 2)]
  #Covariate values
  Covariate <- PPP$rf
  CovariateObs <- safelookup(Covariate, pp, warn = F)
  
  Pearson <- function(L, Z){
    Pear <- apply(L, 1, cor, y = Z, method = "pearson")
    Pear[is.na(Pear)] <- 0
    return(Pear)
  }
  
  RhoObs <- Pearson(L = KLisasObs, Z = CovariateObs)
  
  #Random Shiftings with torus
  randomloc <- function() {
    repeat {
      PP.shift <- Rshift.J(pp)
      pp.shift <- PP.shift$PP
      nn <- npoints(pp.shift)
      if (nn > 4) break
    }
    CovariateSim <- safelookup(Covariate, pp.shift, warn = F)
    corr <- Pearson(L = KLisasObs[,PP.shift$OK], Z = CovariateSim)
    return(list(Corr = corr, N = nn))
  }
  simu <- replicate(nsim, randomloc())
  Rhosimu <- sapply(simu[1, ], "[")
  Nsimu <- sapply(simu[2, ], "[")
  Rhosimu <- cbind(Rhosimu, RhoObs)
  Rhomean <- apply(Rhosimu, 1, mean)
  Ti <- sweep(Rhosimu, 1, Rhomean, FUN = "-")
  Si <- c(Nsimu, npoints(pp))
  TT <- sweep(Ti, 2, sqrt(Si), FUN = "*")
  
  CS <- create_curve_set(list(r = rr, obs = TT[, nsim + 1], sim_m = TT[, 1:nsim]))
  attr(rank_envelope(CS, type = "erl"), "p")
}

#Checking the time for one simulation
system.time(onesimu(nsim = 99, S1 = 0.05, rmaxx = 0.15))

#Using parallel computing for acceletating things
#Note that we only use 99 simulations as illustration
nP <- function(s) mclapply(1:100, mc.cores = 14,
                           FUN = function(x) onesimu(nsim = 99, S1 = s, rmaxx = 0.15))

#Executing parallel procedure with the diferent scales of Thomas pp
P <- sapply(c(0.0125, 0.025,0.05, 0.1), nP, simplify = "array")

#Estimating nominal significance (this case is the fifth line of Table 5)
nominal.rejection <- function(P) mean(unlist(P) <= 0.05)
apply(P, 2, nominal.rejection)
