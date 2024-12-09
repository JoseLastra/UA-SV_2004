# 1.- UASV 2025: Day 1 'Prob modelling and Taylor series' ----------
## ------------------------------------------------------------------------#
## @date 2024-12-09
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
## 3.1.- First order Taylor approximation  for vectors ------
### I.- Show that using this method the variance of SR is approximated by: ------
# SR = NIR/Red
# Var(SR) = Var(NIR/Red) = Var(NIR)/Red^2 + NIR^2*Var(Red)/Red^4 - 2*NIR*Cov(NIR, Red)/Red^3

### II.- Show that using this method the variance of NDVI is approximated by: ------
# NDVI = (NIR - Red)/(NIR + Red)
# Var(NDVI) = Var(NIR - Red)/(NIR + Red)^2 + Var(NIR + Red)*(NIR - Red)^2/(NIR + Red)^4 - 2*Cov(NIR, Red)*(NIR - Red)/(NIR + Red)^3


### III .- Taylor series function for SR ------

## function
# Uncertainty SR by Taylor series method, returns tau.
# nir and red are the measured reflectances in the two
# bands; s_red and s_nir are the standard deviations
# of the measurement errors; rho is the correlation
# between the two measurement errors.
TaylorSR <- function(red, nir, s_red, s_nir, rho) {
  tau_sq <- s_nir^2 / red^2 +
    s_red^2 * nir^2 / red^4 - 2 * rho * s_nir * s_red * nir / red^3
  sqrt(tau_sq) # return tau
}
# call the function with some data
TaylorSR(0.1, 0.6, 0.025, 0.03, 0.8)


### IV.- Taylor series function for NDVI ------
# Uncertainty NDVI by Taylor series method, returns tau.
# Inputs as TaylorSR.
TaylorNDVI <- function(red, nir, s_red, s_nir, rho) {
  tau_sq <- 4 * (s_nir^2 * red^2 + s_red^2 * nir^2 -
    2 * rho * nir * red * s_red * s_nir) / (nir + red)^4
  tau_sq[tau_sq < 0] <- 0
  sqrt(tau_sq)
}

# call the function with some data
TaylorNDVI(0.1, 0.6, 0.025, 0.03, 0.8)

### V.- Add the lines in the below box to your script file and study the code as well as the results of running it. Repeat the calculations for ρ = 0.8. Do the plots conform to your expectations? Save the plots (File | Save as | Png) ----------

red <- 1:500 / 500 # generate a "red" sequence 0.002, 0.004, ..., 1.00
nir <- 1:500 / 500 # ditto for the nir band
# compute SR for every combination of red and nir using the function outer
SR <- outer(red, nir, FUN = function(r, n) n / r)
# compute uncertainty propagation in SR with fixed s_nir and s_red
tau_SR <- outer(red, nir, TaylorSR, s_red = 0.02, s_nir = 0.03, rho = 0)
CV_SR <- (tau_SR / SR) # coefficient of variation, see Wikipedia
# Open graphics window
x11(width = 15, height = 5)
# Split window to display multiple plots
par(mfcol = c(1, 3), mar = c(4, 4, 2, 1), mai = rep(0.6, 4), xaxs = "i", yaxs = "i", cex = 1)
# contour plots with specified levels
contour(
  z = SR, levels = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 50),
  xlab = "RED", ylab = "NIR", main = "SR", labcex = 1
)
contour(
  z = tau_SR, levels = c(0.025, 0.05, 0.1, 0.5, 1, 5, 50),
  xlab = "RED", ylab = "NIR", main = "Tau(SR) by Taylor method", labcex = 1
)
contour(
  z = CV_SR, levels = c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2), labcex = 1,
  xlab = "RED", ylab = "NIR", main = "CV(SR) by Taylor method"
)


### VI.- Make (and save) similar plots for the uncertainty on NDVI given the same ranges of reflectance. Try to understand differences in the coefficients of variation of NDVI and SR. ----------

NDVI <- outer(red, nir, FUN = function(r, n) (n - r) / (n + r))
tau_NDVI <- outer(red, nir, TaylorNDVI, s_red = 0.02, s_nir = 0.03, rho = 0)
CV_NDVI <- (tau_NDVI / NDVI)

# Open graphics window
x11(width = 15, height = 5)
# Split window to display multiple plots
par(mfcol = c(1, 3), mar = c(4, 4, 2, 1), mai = rep(0.6, 4), xaxs = "i", yaxs = "i", cex = 1)
# contour plots with specified levels
contour(
  z = NDVI, levels = c(-0.9, -0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 0.9),
  xlab = "RED", ylab = "NIR", main = "NDVI", labcex = 1
)
contour(
  z = tau_NDVI, levels = c(0.025, 0.05, 0.1, 0.5, 1, 5, 50),
  xlab = "RED", ylab = "NIR", main = "Tau(NDVI) by Taylor method", labcex = 1
)
contour(
  z = CV_NDVI, levels = c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2), labcex = 1,
  xlab = "RED", ylab = "NIR", main = "CV(NDVI) by Taylor method"
)
## --------------------------------------------------------------------------------#
## 3.2.- Taylor approcimation for rasters -------
### Read tiff data ------
fc_lorraine <- rast("1 Monday/Practical/FCLorraine.tif")

### VII.- Plot info from lorraine ------
# close any separate R graphical devices
while (!is.null(dev.list())) dev.off()
op <- par(no.readonly = TRUE) # save default graphical settings
plotRGB(fc_lorraine, 3, 2, 1, stretch = "lin")
par(op) # restore default graphical settings
summary(fc_lorraine, maxsamp = 2e5)
hist(fc_lorraine, 1, maxcell = 2e5)
hist(fc_lorraine, 2, maxcell = 2e5)
hist(fc_lorraine, 3, maxcell = 2e5)


### VIII.- Calculate Cv_sr, tau_sr, and CV_sr for the lorraine data ------
SRim <- app(fc_lorraine,
  fun = function(x) x[3] / x[2], filename = "1 Monday/Practical/SR.tif",
  overwrite = T
) ## calculate SR
plot(SRim, col = gray(1:100 / 100), main = "SR")

tau_SR <- app(fc_lorraine, ## Calculate Tau sr
  fun = function(x) {
    TaylorSR(x[2], x[3], 0.025, 0.03, 0.8)
  }, filename = "1 Monday/Practical/tau_SR.tif",
  overwrite = T
)
plot(tau_SR, col = gray(1:100 / 100), main = "Tau(SR)")

CV_SR <- app(c(tau_SR, SRim),
  fun = function(x) {
    x[1] / x[2]
  },
  filename = "1 Monday/Practical/CV_SR.tif", overwrite = T
)

plot(CV_SR, col = gray(1:100 / 100), main = "CV(SR)", cex = 1)



### IX.- Repeat the above calculations for NDVI and save the results as NDVI.tif, tau_NDVI.tif, and CV_NDVI.tif. ----------
NDVIim <- app(fc_lorraine,
  fun = function(x) (x[3] - x[2]) / (x[3] + x[2]), filename = "1 Monday/Practical/NDVI.tif",
  overwrite = T
) ## calculate NDVI

plot(NDVIim, col = gray(1:100 / 100), main = "NDVI")

tau_NDVI <- app(fc_lorraine, ## Calculate Tau NDVI
  fun = function(x) {
    TaylorNDVI(x[2], x[3], 0.025, 0.03, 0.8)
  }, filename = "1 Monday/Practical/tau_NDVI.tif",
  overwrite = T
)
plot(tau_NDVI, col = gray(1:100 / 100), main = "Tau(NDVI)")

CV_NDVI <- app(c(tau_NDVI, NDVIim), ## Calculate CV NDVI
  fun = function(x) {
    x[1] / x[2]
  },
  filename = "1 Monday/Practical/CV_NDVI.tif", overwrite = T
)
plot(CV_NDVI, col = gray(1:100 / 100), main = "CV(NDVI)")
