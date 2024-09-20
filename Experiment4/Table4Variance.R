#---
#title: Table 3: Variance case (second line execution)
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
Rshift.J <- function(Z, radius = 0.5){
  jump <- runifdisc(1, radius = radius)
  X <- shift(Z, jump)
  Wf <- intersect.owin(Z$window, X$window)
  Xok <- inside.owin(X$x, X$y, Wf)
  return(ppp(x = X$x[Xok], y = X$y[Xok], marks = X$marks[Xok],window = Wf))
}

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

tpp <- function(S = 0.01){
  rThomasInhomSigma(scale = 0.2, kappa = 20, mu = 5, 
                    s0 = 0.01, s1 = S)
}
################################################################################
onesimu <- function(nsim = 99, S1 = 0.01, rmaxx = 0.15)
{
  #simulate a log-Gaussian with scale  =0.1
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
  pp <- pp %mark% 1:pp$n
  #Random Shiftings with variance
  randomloc <- function() {
    repeat {
      pp.shift <- Rshift.J(pp)
      nn <- npoints(pp.shift)
      if (nn > 4) break
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
system.time(onesimu(nsim = 99, S1 = 0.05, rmaxx = 0.15))

#Using parallel computing for acceletating things
#Note that we only use 99 simulations as illustration
nP <- function(s) mclapply(1:1000, FUN = function(x) onesimu(nsim = 999, S1 = s, 
                                                             rmaxx = 0.15), mc.cores = 14)

#Executing parallel procedure with the diferent scales of Thomas pp
P <- sapply(c(0.03, 0.05, 0.08), nP, simplify = "array")


#Estimating nominal significance (this case is the fifth line of Table 2)
nominal.rejection <- function(P) mean(unlist(P) <= 0.05)
apply(P, 2, nominal.rejection)
