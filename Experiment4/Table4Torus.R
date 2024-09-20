#---
#title: Table 3: Torus case (second line execution)
#Fixed LISA functions
#author: "Jonatan A. Gonzalez and Jiri Dvorak"
#---
################################################################################
#Packages
library(spatstat)
library(GET)
library(doParallel)
################################################################################
#Preliminary functions
#Generating random field observations
GRF <- function(beta = 0.2, W = owin()){
  attr(rLGCP(win = W, scale = beta), "Lambda")
}

#Simulation of a inhomogeneus thomas pattern with inhomogenoeus variance
rThomasInhomSigma <- function(scale, kappa, mu, s0, s1){
  x.out <- y.out <- NULL
  W <- owin()
  frame <- boundingbox(W)
  dilated <- grow.rectangle(frame, 4 * scale)
  parents <- rpoispp(kappa, win = dilated)
  Z0 <- abs(GRF(beta = scale, W = dilated))
  sigma <- eval.im(s0 + s1 * Z0)
  
  if (parents$n > 0){
    n.offsprings <- rpois(n = parents$n, lambda = mu)
    for (i in 1:parents$n){
      if (n.offsprings[i] == 0) {next}
      aux.sigma <- sigma[parents[i]]
      aux.dis <- rnorm(n = 2 * n.offsprings[i], mean = 0, sd = aux.sigma)
      x.aux <- parents$x[i] + aux.dis[1:n.offsprings[i]]
      y.aux <- parents$y[i] + aux.dis[(n.offsprings[i] + 1):(2 * n.offsprings[i])]
      x.out <- c(x.out, x.aux)
      y.out <- c(y.out, y.aux)
    }
  }
  else return(rpoispp(0))
  ok <- inside.owin(x.out, y.out, W)
  out <- ppp(x=x.out[ok], y=y.out[ok], window = W)
  return(list(pp = out, rf = Z0[W]))
}

#Generating Thomas pp
tpp <- function(S = 0.01){
  rThomasInhomSigma(scale = 0.2, kappa = 20, mu = 5, 
                    s0 = 0.01, s1 = S)
}

################################################################################
onesimu <- function(nsim = 99, S1 = 0.01, rmaxx = 0.15)
{
  PP <- tpp(S1)
  pp <- PP$pp
  # Calculate LISA functions based on K
  TemplateLisa <- localL(pp, verbose = F, rmax = rmaxx, correction = "translate")
  r0 <- floor(seq(1,513, length.out = 51))
  rr <- TemplateLisa$r[r0]
  KLisasObs <- as.matrix(TemplateLisa)
  KLisasObs <- KLisasObs[r0, 1:(dim(KLisasObs)[2] - 2)]
  #Covariate values
  Covariate <- PP$rf
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
system.time(onesimu(nsim = 99, S1 = 0.05, rmaxx = 0.15))

#Using parallel computing for acceletating things
#Note that we only use 99 simulations as illustration
nP <- function(s) mclapply(1:100, FUN = function(x) onesimu(nsim = 99, S1 = s, 
                                                            rmaxx = 0.15), mc.cores = 14)

#Executing parallel procedure with the diferent scales of Thomas pp
P <- sapply(c(0.03, 0.05), nP, simplify = "array")


#Estimating nominal significance (this case is the fifth line of Table 2)
nominal.rejection <- function(P) mean(unlist(P) <= 0.05)
apply(P, 2, nominal.rejection)
