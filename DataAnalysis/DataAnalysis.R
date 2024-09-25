#---
#title: Data analysis: torus case 
#Fixed LISA functions
#author: "Jonatan A. Gonzalez"
#---
# Please load the data "Data.RData"
#Please load the functions 2_LfunctionsData.R
################################################################################

# Data analysis

library(spatstat)
library(spatstat.utils)
library(GET)
library(doParallel)
library(sparr)
library(spatstat.utils)
library(sostatpp)

# load data ("Data.RData")
#Load functions 2_LfunctionsData.R

lambda <- density.ppp(PP, sigma = 50, diggle = T, at = "points")

################################################################################

pp <- PP
rmaxx <- rmax.rule("K", W = pp$window, lambda = intensity(pp))
pp <- pp %mark% lambda

# Calculate LISA functions based on K
TemplateLisa <- localL1inhom(pp, lambda = lambda, verbose = F, 
                             correction = "translate", rmax = rmaxx)
#r0 <- floor(seq(1,513, length.out = 51))
rr <- TemplateLisa$r
KLisasObs <- as.matrix(TemplateLisa)
KLisasObs <- KLisasObs[, 1:(dim(KLisasObs)[2] - 2)]

#Covariate values
CovariateObs <- lapply(Cov, function(Co) safelookup(Co, pp))

Pearson <- function(L, Z){
  Pear <- apply(L, 1, cor, y = Z, method = "pearson")
  Pear[is.na(Pear)] <- 0    
  return(Pear)
}

RhoObs <- lapply(CovariateObs, function(C) Pearson(L = KLisasObs, Z = C))

############### Random Shifting with torus #####################################
################################################################################
randomloc <- function(C) {
  pp.shift <- rshift(pp, width = 450, height = 5, edge = "torus")
  CovariateSim <- safelookup(C, pp.shift)
  return(Pearson(L = KLisasObs, Z = CovariateSim))
}

nsim <- 9999
Process <- function(R,C) {
  simu <- replicate(nsim, randomloc(C))
  create_curve_set(list(r = rr, obs = R, sim_m = simu))
}

simu <- mcmapply(Process, RhoObs, Cov, mc.cores = 15, SIMPLIFY = F)

#simuPearsonTorus <- simu
#simuKendallTorus <- simu
save(simuKendallTorus, simuPearsonTorus, file = "DataPerm.RData")


plot(rank_envelope(simu[[1]], type = "erl"))
plot(rank_envelope(simu[[2]], type = "erl"))
plot(rank_envelope(simu[[3]], type = "erl"))
plot(rank_envelope(simu[[4]], type = "erl"))


#simuPearsonTorus <- simu


############### Random Shifting with variance ##################################
################################################################################
#Preliminary functions
Rshift.J <- function(Z, width = 450, height = 5){
  #jump <- runifdisc(1, radius = radius)
  jump <- list(x = runif(1, min = 0, max = width), 
               y = runif(1, min = 0, max = height))
  X <- shift(Z, jump)
  Wf <- intersect.owin(Z$window, X$window)
  Xok <- inside.owin(X$x, X$y, Wf)
  return(list(PP = ppp(x = X$x[Xok], y = X$y[Xok], window = Wf), OK = Xok))
}

randomlocVariance <- function(C) {
  repeat {
    PP.shift <- Rshift.J(pp, width = 450, height = 5)
    pp.shift <- PP.shift$PP
    nn <- npoints(pp.shift)
    if (nn > 4) break
  }
  CovariateSim <- safelookup(C, pp.shift, warn = F)
  corr <- Pearson(L = KLisasObs[,PP.shift$OK], Z = CovariateSim)
  return(list(Corr = corr, N = nn))
}

nsim <- 9999
Process <- function(R,C) {
  simu <- replicate(nsim, randomlocVariance(C))
  Rhosimu <- sapply(simu[1, ], "[")
  Nsimu <- sapply(simu[2, ], "[")
  Rhosimu <- cbind(Rhosimu, R)
  Rhomean <- apply(Rhosimu, 1, mean)
  Ti <- sweep(Rhosimu, 1, Rhomean, FUN = "-")
  Si <- c(Nsimu, npoints(pp))
  TT <- sweep(Ti, 2, sqrt(Si), FUN = "*")
  create_curve_set(list(r = rr, obs = TT[, nsim + 1], 
                        sim_m = TT[, 1:nsim]))
}

simu <- mcmapply(Process, RhoObs, Cov, mc.cores = 15, SIMPLIFY = F)

simuPearsonVariance <- simu
simuKendallVariance <- simu

plot(rank_envelope(simuKendallTorus[[1]], type = "erl"))
plot(rank_envelope(simu[[2]], type = "erl"))
plot(rank_envelope(simu[[3]], type = "erl"))
plot(rank_envelope(simu[[4]], type = "erl"))

################################################################################
save(simuPearsonVariance,
     simuPearsonTorus,
     simuKendallVariance,
     simuKendallTorus, 
     file = "DataPerm.RData")


################################################################################
#Test UTE 
beiquads <- twoquadsets(PP, nx = 12, ny = 1, minpoints = 5)
beistyle <- list(hi = simplist(col="red",  alpha=.4, col.win="red", alpha.win=.4), 
                 lo = simplist(col="blue", alpha=.4, col.win="blue"))
quadratsplot(PP, beiquads, beistyle, pch = 16, cex = .4, 
             main = "Beilschmiedia quadrats for testing")
nperm <- 9999
sos.test(reweighted(PP, intensity = lambda),
         beiquads, rmax = 2.5, rlen = 512,
         use.tbar=TRUE, nperm = nperm)
