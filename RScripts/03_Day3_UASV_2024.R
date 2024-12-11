# 1.- UASV 2025: Day 3 'Model uncertainty' ----------
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

# Try Bayes’ law yourself (tip: make also use of the ‘Law of
#                          Total Probability’):
# ▪ Thunderstorms happen on average 2 out of 100 days (2%)
# ▪ The probability of rain in case of a thunderstorm is 80%
# ▪ The probability of rain in case of no thunderstorm is 5%
# ▪ What is the probability of thunderstorm in case of rain?
#Calculate the probability of thunderstorm in case of rain
#P(Thunderstorm | Rain) = P(Rain | Thunderstorm) * P(Thunderstorm) / P(Rain)
#P(Thunderstorm | Rain) = 0.8 * 0.02 / P(Rain)

#Calculate the probability of rain
#P(Rain) = P(Rain | Thunderstorm) * P(Thunderstorm) + P(Rain | No Thunderstorm) * P(No Thunderstorm)
#P(Rain) = 0.8 * 0.02 + 0.05 * 0.98







