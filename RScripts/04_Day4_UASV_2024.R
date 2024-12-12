# 1.- UASV 2025: Day 4 'Validation' ----------
## ------------------------------------------------------------------------#
## @date 2024-12-12
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
## 3.1.- Map purity index Drenthe -----
drenthe_data <- read_csv('4 Thursday/Lecture/soil_obs_drenthe.csv')

### contingency table for columns Observation, Soil_old, Soil_new ------

# P=thick
# peat; mP=thick peat soil with mineral surface horizon; PY=thin peat soil; mPY=thin
# peat soil with mineral surface horizon; BF=brown forest soil; PZ=podzol; E=earth
# soil; PS=plaggen soil; T=till soil; S=sandy vague soil
soil_classes_count <- drenthe_data %>%
  group_by(Observation) %>%
  count()


conf_matrix_old <- drenthe_data %>%
  select(Observation, Soil_old) %>%
  table()
print(conf_matrix_old)

conf_matrix_new <- drenthe_data %>%
  select(Observation, Soil_new) %>%
  table()
print(conf_matrix_new)

### sum and sum diagonal ----
## Purity is the proportion of the total number of observations that were correctly classified 
purity_old <- sum(diag(conf_matrix_old)/sum(conf_matrix_old))
purity_new <- sum(diag(conf_matrix_new)/sum(conf_matrix_new))









