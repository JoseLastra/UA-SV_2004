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

## 3.2.- USLE erosion model --------
# Monte Carlo uncertainty propagation USLE model
n <- 1000
R <- rnorm(n, 297, 72)
K <- rnorm(n, 0.10, 0.05)
L <- rnorm(n, 2.13, 0.05)
S <- rnorm(n, 1.17, 0.12)
C <- rnorm(n, 0.63, 0.15)
P <- rnorm(n, 0.50, 0.10)
E <- R * K * L * S * C * P
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
