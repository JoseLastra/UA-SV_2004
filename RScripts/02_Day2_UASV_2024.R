# 1.- UASV 2025: Day 2 'MONTECARLO' ----------
## ------------------------------------------------------------------------#
## @date 2024-12-10
## @project C:/github/UA-SV_2004
## @R version R version 4.3.2 (2023-10-31 ucrt)
## @OS system Windows 10 x64
## @author Jose A. Lastra
## @email jose.lastramunoz@wur.nl | jose.lastra@pucv.cl
## ------------------------------------------------------------------------#
# 2.- Libraries --------
pacman::p_load(tidyverse, sf, terra)
## Clean environment ----
rm(list = ls(all = T))

## ------------------------------------------------------------------------#
# 3.- Exercises -----
## 3.1.- Monte Carlo approximation for solar energy flux --------
# Plot solar energy as a function of interception angle
phi <- seq(0, 50, 0.01)
hist(phi, col = "LightBlue", main = "phi")
z <- 185 * phi / (35 + phi)
hist(z, col = "Orange", main = "z")
plot(phi, z, type = "l", lwd = 2, col = "darkgreen")

# Monte Carlo uncertainty propagation
n <- 10000
phi_sample <- rnorm(n, 28, 5.5)
z_sample <- 185 * phi_sample / (35 + phi_sample)
mean(z_sample)
sd(z_sample)
hist(phi_sample, col = "LightBlue", main = "phi")
hist(z_sample, col = "Orange", main = "z")

## ----------------------------------------------------------------------#
## 3.2.- USLE erosion model --------
# Monte Carlo uncertainty propagation USLE model
n <- 5000
set.seed(12345)
R <- rnorm(n, 297, 72)
K <- rnorm(n, 0.10, 0.05)
L <- rnorm(n, 2.13, 0.05)
S <- rnorm(n, 1.17, 0.12)
C <- rnorm(n, 0.63, 0.15)
P <- rnorm(n, 0.50, 0.10)
E <- R * K * L * S * C * P # running the model
mean(E)
sd(E)
hist(E, col = "LightBlue")
# Setting each uncertainty source to its mean, one by one
E_R <- 297 * K * L * S * C * P
E_K <- R * 0.10 * L * S * C * P
E_L <- R * K * 2.13 * S * C * P
E_S <- R * K * L * 1.17 * C * P
E_C <- R * K * L * S * 0.63 * P
E_P <- R * K * L * S * C * 0.50
# Uncertainty source contributions
C_R <- 1 - var(E_R) / var(E)
C_K <- 1 - var(E_K) / var(E)
C_L <- 1 - var(E_L) / var(E)
C_S <- 1 - var(E_S) / var(E)
C_C <- 1 - var(E_C) / var(E)
C_P <- 1 - var(E_P) / var(E)
C_R
C_K
C_L
C_S
C_C
C_P
# Do they sum to one?
C_R + C_K + C_L + C_S + C_C + C_P

## -----------------------------------------------------------------------#
# 4.- Afternoon practical session: Veg Indices Monte Carlo ------
# Load the data
fc_lorraine <- rast("1 Monday/Practical/FCLorraine.tif")


s_RED <- 0.025 # sigma RED
s_NIR <- 0.03 # sigma NIR
rho <- 0.8 # correlation errors in NIR and RED
# covariance matrix and Choleski factorization
covmat <- matrix(c(s_RED^2, rep(rho * s_RED * s_NIR, 2), s_NIR^2), 2, 2)
M <- t(chol(covmat)) # Choleski factorization

all(M %*% t(M) == covmat) # just a check; %*% is matrix multiplication
set.seed(1234) # initialize pseudo random number generator
n <- 1000 # sample size
Z <- matrix(rnorm(2 * n), 2, n) # 2 rows, n columns; independent draws from

# the standard normal distribution
devs <- t(M %*% Z) # independent draws from multivariate normal

# distribution with covmat
plot(devs, pch = 19, cex = 0.7, xlab = "red error", ylab = "NIR error")

## --------------------------------------------------------------------#
## Single point
# compute mean and sd SR from a series of realised deviations of
# red and nir (in x) and central values for red and nir
MC_SR <- function(red, nir, devs) {
  RED <- red + devs[, 1]
  NIR <- nir + devs[, 2]
  # ignore red reflectances < 0.025% (division by values close to zero)
  RED[RED < 0.00025] <- NA
  # return SR realizations without missing values
  na.omit(NIR / RED)
}
set.seed(1234567) # initialize pseudo random number generator
n <- 200 # sample size
Z <- matrix(rnorm(2 * n), 2, n) # 2 rows, n columns; independent draws from
# standard normal distribution
devs <- t(M %*% Z)
SRsamp <- MC_SR(0.1, 0.6, devs)
sd(SRsamp)
mean(SRsamp)

## 4.1.- With the computed data: (i) Histogram, (ii) number of replicates --------
hist(SRsamp, col = "LightBlue", main = "SR") # skew to the left (positive skew)

n <- 3000 # sample size
Z <- matrix(rnorm(2 * n), 2, n) # 2 rows, n columns; independent draws from
# standard normal distribution
devs <- t(M %*% Z)
SRsamp <- MC_SR(0.1, 0.6, devs)
sd(SRsamp)
mean(SRsamp)
hist(SRsamp, col = "LightBlue", main = "SR") # skew to the left (positive skew)


## 4.2.- Automate the process with a function --------
N <- 1:100 * 1000 # a series of values for n
set.seed(123456) # set a seed for the PRN generator
sdSR <- function(n, red, nir) {
  Z <- matrix(rnorm(2 * n), 2, n)
  devs <- t(M %*% Z)
  sd(MC_SR(red, nir, devs))
}
funSR <- function(n, red, nir) replicate(50, sdSR(n, red, nir))
sdsSR <- sapply(N, function(n) funSR(n, 0.1, 0.6))
sd_var <- apply(sdsSR, 2, function(x) sd(x^2))
plot(sapply(N, rep, times = 50), sdsSR,
  xlab = "sample size",
  ylab = "sd fom sample",
  main = "red = 0.1, NIR = 0.6", pch = 16, cex = 0.4
)
plot(N, sd_var,
  xlab = "sample size", ylab = "sd(var(SR))",
  main = "red = 0.1, NIR = 0.6"
)
lines(N, sd_var[100] / sqrt(N / N[100]), col = "red") # inverse sqrt relation

## 4.3.- Repeat for NDVI --------
MC_NDVI <- function(red, nir, devs) {
  RED <- red + devs[, 1]
  NIR <- nir + devs[, 2]
  Den <- RED + NIR
  # ignore sampled denominators near zero
  Den[abs(Den) < 0.0025] <- NA
  samp <- (NIR - RED) / Den
  na.omit(samp)
}

sd(MC_NDVI(0.1, 0.6, devs))
mean(MC_NDVI(0.1, 0.6, devs))
hist(MC_NDVI(0.1, 0.6, devs), col = "LightBlue", main = "NDVI")

N <- 1:100 * 1000 # a series of values for n
set.seed(123456) # set a seed for the PRN generator
sd_NDVI <- function(n, red, nir) {
  Z <- matrix(rnorm(2 * n), 2, n)
  devs <- t(M %*% Z)
  sd(MC_NDVI(red, nir, devs))
}
fun_NDVI <- function(n, red, nir) replicate(50, sd_NDVI(n, red, nir))
sds_NDVI <- sapply(N, function(n) fun_NDVI(n, 0.1, 0.6))
sd_var <- apply(sds_NDVI, 2, function(x) sd(x^2))
plot(sapply(N, rep, times = 50), sds_NDVI,
  xlab = "sample size",
  ylab = "sd fom sample",
  main = "red = 0.1, NIR = 0.6", pch = 16, cex = 0.4
)
plot(N, sd_var,
  xlab = "sample size", ylab = "sd(var(SR))",
  main = "red = 0.1, NIR = 0.6"
)
lines(N, sd_var[100] / sqrt(N / N[100]), col = "red") # inverse sqrt relation

##----------------------------------------------------------------#
## Contour plot
red <- 1:500/500 # generate a sequence 0.002, 0.004, ., 1.00
nir <- 1:500/500
rednir <- expand.grid(red, nir)
mu_sdSR <- function(x) {
  samp <- MC_SR(x[1], x[2], devs=devs)
  c(mean(samp), sd(samp))
}

## SR
n = 20000
set.seed(1234567)
Z <- matrix(rnorm(2*n),2,n)
devs <- t(M %*% Z)
mu_sd <- apply(rednir, 1, mu_sdSR)
MCmu_SR <- matrix(mu_sd[1,], 500, 500)
MCsd_SR <- matrix(mu_sd[2,], 500, 500)
CV_SR <- (MCsd_SR / MCmu_SR)

## NDVI
mu_sdNDVI <- function(x) {
  samp <- MC_NDVI(x[1], x[2], devs=devs)
  c(mean(samp), sd(samp))
}


set.seed(1234567)
Z <- matrix(rnorm(2*n),2,n)
devs <- t(M %*% Z)
mu_sd <- apply(rednir, 1, mu_sdNDVI)
MCmu_NDVI <- matrix(mu_sd[1,], 500, 500)
MCsd_NDVI <- matrix(mu_sd[2,], 500, 500)
CV_NDVI <- (MCsd_SR / MCmu_SR)

## 4.5.- Contour plots for SR and NDVI --------
x11(width=15, height=5)

par(mfcol=c(1,3),mar=c(4,4,2,1), mai=rep(0.6,4), xaxs="i", yaxs="i")

contour(red, nir, MCmu_SR, levels=c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 50), 
        xlab="red", ylab="NIR", main="SR")
contour(red, nir, MCsd_SR, levels=c(0.025, 0.05, 0.1, 0.5, 1, 5, 50),
        xlab="red", ylab="nir", main="sd(SR) by Monte Carlo method")
contour(red, nir, CV_SR, levels=c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2),
        xlab="red", ylab="NIR", main="CV(SR) by Monte Carlo method")

# NDVI
x11(width=15, height=5)

par(mfcol=c(1,3),mar=c(4,4,2,1), mai=rep(0.6,4), xaxs="i", yaxs="i")

contour(red, nir, MCmu_NDVI, levels=c(-0.95, -0.5, -0.2, -0.1, -0.05, 
                                      0, 0.05, 0.1, 0.2, 0.5, 0.95), 
        xlab="red", ylab="NIR", main="NDVI")
contour(red, nir, MCsd_NDVI, levels=c(0.025, 0.05, 0.1, 0.5, 1, 5, 50),
        xlab="red", ylab="nir", main="sd(NDVI) by Monte Carlo method")
contour(red, nir, CV_NDVI, levels=c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2),
        xlab="red", ylab="NIR", main="CV(NDVI) by Monte Carlo method")

## --------------------------------------------------------------------#
## Image
## 4.6.- NDVI and SR images --------
false_color <- rast("1 Monday/Practical/FCLorraine.tif")

SRim <- rast("1 Monday/Practical/SR.tif")
sd_SRmc <- app(false_color, fun = function(x) sd(MC_SR(x[2], x[3], devs)),
               filename="2 Tuesday/Practical/sd_SR.tif", overwrite = T)
CV_SRmc <- app(c(sd_SRmc, SRim), fun=function(x) {x[1]/x[2]}, 
               filename="2 Tuesday/Practical/CV_SRmc.tif", overwrite=T)


x11(width=15, height=5)

par(mfcol=c(1,3),xaxs="i", yaxs="i",cex=1) # split window
plot(SRim, col=gray(1:100/100), main="SR", mar=c(2,2,2,4))

plot(sd_SRmc, col=gray(1:100/100), main= "sd(SR) Monte Carlo", range=c(0,50),
     mar=c(2,2,2,4))

plot(CV_SRmc, col=gray(1:100/100), main= "CV(SR) Monte Carlo", range=c(0, 2.5),
     mar=c(2,2,2,4))


## NDVI
NDVIim <- rast("1 Monday/Practical/NDVI.tif")

sd_NDVImc <- app(false_color, fun = function(x) sd(MC_NDVI(x[2], x[3], devs)),
                 filename="2 Tuesday/Practical/sd_NDVI_2.tif", overwrite = T)

CV_NDVImc <- app(c(sd_NDVImc, NDVIim), fun=function(x) {x[1]/x[2]}, 
                 filename="2 Tuesday/Practical/CV_NDVImc.tif", overwrite=T)

x11(width=15, height=5)

par(mfcol=c(1,3),xaxs="i", yaxs="i",cex=1) # split window
plot(NDVIim, col=gray(1:100/100), main="NDVI", mar=c(2,2,2,4))

plot(sd_NDVImc, col=gray(1:100/100), main= "sd(NDVI) Monte Carlo", range=c(0,1),
     mar=c(2,2,2,4))

plot(CV_SRmc, col=gray(1:100/100), main= "CV(NDVI) Monte Carlo", range=c(0, 2.5),
     mar=c(2,2,2,4))

## 4.7.- Assess the probability of NDVI < 0.3 (a typical value for differentiating between vegetated and non-vegetated areas) using the precompute data and interpret the result -----

MC_NDVI_p <- function(red, nir, devs, thr) {
  RED <- red + devs[, 1]
  NIR <- nir + devs[, 2]
  Den <- RED + NIR
  # ignore sampled denominators near zero
  Den[abs(Den) < 0.0025] <- NA
  samp <- (NIR - RED) / Den
  samp_count <- samp[samp < thr]
  pp <- length(samp_count) / length(samp)
  na.omit(pp)
}


pp_NDVImc <- app(false_color, fun = function(x) MC_NDVI_p(x[2], x[3], devs, thr = 0.3),
                 filename="2 Tuesday/Practical/pp_NDVI_2.tif", overwrite = T)

pp_NDVImc_sytze <- app(false_color, fun = function(x) mean(ifelse(MC_NDVI(x[2], x[3], devs) < 0.3,1,2)),
                 filename="2 Tuesday/Practical/pp_NDVI_3.tif", overwrite = T)

all.equal(pp_NDVImc, pp_NDVImc_sytze)
plot(pp_NDVImc)
plot(pp_NDVImc_sytze)
plotRGB(false_color, r=3, g=2, b=1, stretch="lin", main="False Color")







