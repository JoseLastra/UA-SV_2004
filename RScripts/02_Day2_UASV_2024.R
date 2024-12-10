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

