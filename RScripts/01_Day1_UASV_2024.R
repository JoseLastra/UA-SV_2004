# 1.- UASV 2025: Day 1 'Prob modelling and Taylor series' ----------
##------------------------------------------------------------------------#
## @date 2024-12-09
## @project C:/github/UA-SV_2004
## @R version R version 4.3.2 (2023-10-31 ucrt)
## @OS system Windows 10 x64
## @author Jose A. Lastra
## @email jose.lastramunoz@wur.nl | jose.lastra@pucv.cl
##------------------------------------------------------------------------#
# 2.- Libraries --------
pacman::p_load(tidyverse, sf, terra)
## Clean environment ----
rm(list = ls(all = T))

##------------------------------------------------------------------------#
# 3.- Exercises -----
## 3.1.- First order Taylor approximation ------
### I.- Show that using this method the variance of SR is approximated by: ------
# SR = NIR/Red
# Var(SR) = Var(NIR/Red) = Var(NIR)/Red^2 + NIR^2*Var(Red)/Red^4 - 2*NIR*Cov(NIR, Red)/Red^3

### II.- Show that using this method the variance of NDVI is approximated by: ------
# NDVI = (NIR - Red)/(NIR + Red)
# Var(NDVI) = Var(NIR - Red)/(NIR + Red)^2 + Var(NIR + Red)*(NIR - Red)^2/(NIR + Red)^4 - 2*Cov(NIR, Red)*(NIR - Red)/(NIR + Red)^3


### III .- Taylor series function for SR ------
## Read tiff data ------
fc_lorraine <- rast('1 Monday/Practical/FCLorraine.tif')
plotRGB(fc_lorraine, stretch = 'lin')

## function
# Uncertainty SR by Taylor series method, returns tau.
# nir and red are the measured reflectances in the two
# bands; s_red and s_nir are the standard deviations
# of the measurement errors; rho is the correlation
# between the two measurement errors.
TaylorSR <- function(red, nir, s_red, s_nir, rho) {
  tau_sq <- s_nir^2/red^2 +
    s_red^2*nir^2/red^4-2*rho*s_nir*s_red*nir/red^3
  sqrt(tau_sq) # return tau
}
# call the function with some data
TaylorSR(0.1, 0.6, 0.025, 0.03, 0.8)


### IV.- Taylor series function for NDVI ------
# Uncertainty NDVI by Taylor series method, returns tau.
# Inputs as TaylorSR.
TaylorNDVI <- function(red, nir, s_red, s_nir, rho) {
  tau_sq <- 4*(s_nir^2*red^2+s_red^2*nir^2 - 
                 2*rho*nir*red*s_red*s_nir)/(nir+red)^4
  tau_sq[tau_sq < 0] <- 0
  sqrt(tau_sq)
}

# call the function with some data
TaylorNDVI(0.1, 0.6, 0.025, 0.03, 0.8)
