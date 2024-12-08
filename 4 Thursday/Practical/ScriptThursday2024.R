
# ===================================================================
#       UA&SV validation practical THURSDAY (Sytze de Bruin)
# ===================================================================

#### load libraries ####
library(terra)
library(ranger)
library(spatstat)
library(sf)

#### get and plot data ####
setwd("D:/PERC/Thursday") # path to data folder

OCSstack <- rast('OCSstack.tif')
# the first layer contains reference soil organic carbon stock (OCS) values
# the others are covariates, see https://doi.org/10.1016/j.ecoinf.2022.101665

plot(OCSstack[[1]], range=c(18, 155)) # holes are rock outcrops, water & urban


# *******************************************************************
####                 predict OCS with RF model                   ####
# *******************************************************************

# simple random sampling (training), sample size = 2500
set.seed(12345)
idx <- as.integer(spatSample(OCSstack[[1]], 2500, method="random", 
                             na.rm=T, cells=T)[,1])

# extract training data from the raster stack
OCStrain <- extract(OCSstack, idx)

# random forest training, mostly using defaults for hyper parameters
RFmodel <- ranger(ocs~., data = OCStrain, num.trees = 100)

# random forest prediction; takes several minutes
# the wrapper "predfun" is used to retrieve only the RF predictions and to 
# set NA outside the study area 
predfun <- function(model, data) {
  predict(model, data, num_threads=0,)$predictions * 
    ifelse(is.na(data[,2]), NA, 1)
}
RFpredict <- predict(OCSstack, RFmodel, fun=predfun, 
                     filename = "OCSpreds.tif", overwrite=TRUE)

# quick and dirty plot
plot(RFpredict, range=c(18, 155))


# *******************************************************************
####                     validation                             ####
# *******************************************************************
# avoid leakage training data for validation; set points to NA
OCSstack[[1]][idx] <- NA


#### validation using the entire map; we have exhaustive reference data ####
resids <- OCSstack[[1]] - RFpredict  # residual = reference - prediction
MEr  <- global(resids, "mean", na.rm=TRUE)[[1]]   # mean error
RMSE <- global(resids, function(x) {sqrt(mean(x^2, na.rm = TRUE))})[[1]] # RMSE
SSR  <- global(resids, function(x) {sum(x^2, na.rm = TRUE)})[[1]] # squared res
mref <- global(OCSstack[[1]], "mean", na.rm=TRUE)[[1]]  # mean reference OCS
SST  <- global(OCSstack[[1]], function(x) {sum((x-mref)^2, na.rm=TRUE)})[[1]]
MEC  <- 1 - SSR/SST


#### validation using SIMPLE RANDOM SAMPLE ####

#### first just once with sample size 350 ####
set.seed(6789)
idx <- as.integer(spatSample(OCSstack[[1]], 350, method="random", 
                             na.rm=T, cells=T)[,1])
sres <- extract(resids, idx)[,1] # sampled residuals

# compute the error metrics
srsME  <- mean(sres)           # mean error
srsSdME <- sd(sres)/350^0.5    # standard deviation of mean error
srsSdSE <- sd(sres^2)/350^0.5  # standard deviation of mean squared error
srsRMSE <- sqrt(mean(sres^2))  # RMSE simple random sampling
srsRMSElw <- sqrt(max(0, mean(sres^2) + qt(0.05, 349) * srsSdSE)) # low 90% CI RMSE
srsRMSEup <- sqrt(max(0, mean(sres^2) + qt(0.95, 349) * srsSdSE)) # upp 90% CI RMSE

#### now repeat many (100) times (only for ME) ####
getMetrics <- function(){
  idx <- as.integer(spatSample(OCSstack[[1]], 350, method="random", 
                               na.rm=T, cells=T)[,1])
  sres <-  extract(resids, idx)[,1] # sampled residuals
  c(mean(sres), sd(sres)/350^0.5)
}
sMEs <- replicate(100, getMetrics()) # takes several minutes

hist(sMEs[1,])

# how  often is true ME within confidence interval (uses t-distribution)?
inside <- ifelse(sMEs[1,] + qt(0.05, 349) * sMEs[2,] < MEr & 
                   sMEs[1,] + qt(0.95, 349) * sMEs[2,] > MEr, 1, 0)
table(inside)


#### if time allows you can also implement for RMSE ####


#### validation using STRATIFIED RANDOM SAMPLE; stratification based on ####
# layer 9 in OCSstack: GLC2017 (generalized 2017 Copernicus land cover map 
# categories) there are 7 strata; first EQUAL SAMPLE SIZE (50) per stratum

# to save time, we do it just once
set.seed(12345)
idx <- as.integer(spatSample(OCSstack[[9]], 50, method="stratified", 
                             na.rm=T, cells=T)[,1])

# strata weights
w <- as.numeric(global(OCSstack[[9]], table))
w <- w/sum(w)

# set (temporary) variables to zero
stratME_1 <- 0
tmpVar    <- 0
VarSqr    <- 0
MSE       <- 0

for (strat in 1:7){
  sidx <- 1:50 + (strat-1) * 50
  sres <-  extract(resids, idx[sidx])[,1]        # sampled residuals
  stratME_1 <- stratME_1 + w[strat] * mean(sres) # mean error
  tmpVar <- tmpVar + w[strat]^2 * var(sres)/50   # variance of mean residuals
  VarSqr <- VarSqr + w[strat]^2 * var(sres^2)/50 # variance of squared mean res.
  MSE <- MSE + w[strat] * mean(sres^2)           # mean squared error
}
stratSD_1   <- sqrt(tmpVar)
stratRMSE_1 <- sqrt(MSE)
stratRMSElw_1 <- sqrt(max(0, MSE + qt(0.05, 349) * sqrt(VarSqr))) # low 90% CI RMSE
stratRMSEup_1 <- sqrt(max(0, MSE + qt(0.95, 349) * sqrt(VarSqr))) # upp 90% CI RMSE
        
#### same but now with PROPORTINAL SAMPLING; at least 2 points per stratum ####
sizes <- pmax(2, floor(350 * w))
sum(sizes) # great! same sample size

# set (temporary) variables to zero
stratME_2 <- 0
tmpVar    <- 0
VarSqr    <- 0
MSE       <- 0

set.seed(12345)
for (strat in 1:7){
  idx <- as.integer(spatSample(OCSstack[[1]], sizes[strat], na.rm=T, 
                               cells=T)[,1])
  sres <-  extract(resids, idx)[,1]                 # sampled residuals
  stratME_2 <- stratME_2 + w[strat] * mean(sres)    # mean error
  tmpVar <- tmpVar + w[strat]^2 * var(sres)/sizes[strat]   # var. mean res.
  VarSqr <- VarSqr + w[strat]^2 * var(sres^2)/sizes[strat] # var. sq. mean res.
  MSE <- MSE + w[strat] * mean(sres^2)              # mean squared error
}

stratSD_2   <- sqrt(tmpVar)
stratRMSE_2 <- sqrt(MSE)
stratRMSElw_2 <- sqrt(max(0, MSE + qt(0.05, 349) * sqrt(VarSqr))) # low 90% CI RMSE
stratRMSEup_2 <- sqrt(max(0, MSE + qt(0.95, 349) * sqrt(VarSqr))) # upp 90% CI RMSE


# ****************************************************************************
####  PICP: prediction interval coverage probability, for random sample,  ####
#     n = 1000; need to setup the RF model for quantile prediction
# ****************************************************************************
QRFmodel <- ranger(ocs~., data = OCStrain, num.trees = 100, quantreg = TRUE)
set.seed(6789)

# random sample of size 1000
idx <- as.integer(spatSample(OCSstack[[1]], 1000, method="random", 
                             na.rm=T, cells=T)[,1])

# we only need to predict the quantiles at the sample points
dat <- extract(OCSstack, idx)

# predict the quantiles
quants <- predict(QRFmodel, dat, type = "quantiles", 
                  quantiles = 0:20 * 0.05)$predictions

# frequency of points falling in CIs of different widths
CIW <- 0:10/10
frq <- numeric(11)

for(i in 0:10){
  frq[i+1] <- sum(dat$ocs > quants[,11-i] & dat$ocs < quants[,11+i])
}

plot(CIW, frq/1000, asp = 1, ylab="proportion")
abline(0,1, col="red")


# *******************************************************
#### Cross validation using SRS and clustered sample ####
# *******************************************************

# we use one strongly clustered sample (n = 5000) from 
# https://doi.org/10.1016/j.ecoinf.2022.101665
load("OCSdata050.Rdata")  # strongly clustered OCS sample
names(OCSdata)

plot(OCSdata$xcoord, OCSdata$ycoord, asp=1, pch=19, cex=0.5)

# random 10-fold cross validation
# make the random partition
set.seed(12345)
idx <- sample(5000)

resids <- c()

# compute the residuals, one fold at a time
for (i in 1:10){
  istart <-  1 + (i-1) * 500
  istop  <- istart + 499
  isub <- idx[istart:istop]
  RFmodelCV <- ranger(ocs~., data = OCSdata[-isub,], num.trees = 100)
  preds <- predict(RFmodelCV,  OCSdata[isub,])$predictions
  resids <- c(resids, OCSdata$ocs[isub] - preds)             
}

#### simple random CV ####
ME_rCV   <- mean(resids)
RMSE_rCV <- sqrt(mean(resids^2))
MEC_rCV <- 1 - sum(resids^2)/sum((OCSdata$ocs - mean(OCSdata$ocs))^2)

#### give less weight to areas having greater data density ####
#    which are likely to have greater accuracy

# define area over which density has to be computed
studarea <- rast("mask.tif")
pols <- as.polygons(studarea)
polssf <- st_as_sf(pols)
polbuf <- st_buffer(polssf, 2500) # expand a bit
rm(studarea, pols, polssf)

#### sample inytensity calculation using spatstat; don't mind the details ####
getDensity <- function(x, y, sf_pol, rsl){
  CRSlaea <- paste0("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 ",
                      "+ellps=GRS80 +units=m +no_defs") # coordinate ref. syst.
  win <- as.owin(sf_pol)  # analysis window
  spp <- ppp(x, y, win)   # create point pattern
  s <- bw.CvL(spp)        # bandwidth selection for kernel density
  den <- density.ppp(spp, eps=rsl, sigma=s, positive=T) # kernel dens. est.
  # below den is converted into a SpatRaster
  denmat <- as.matrix.im(den)
  denrot <- transmat(denmat, from="spatstat", to="Europe")
  denext <- ext(round(c(den$xrange, den$yrange),0))
  denrast<- rast(denrot, crs=CRSlaea)
  ext(denrast) <- denext
  return(denrast)
}

# calling the function
dens <- getDensity(OCSdata$xcoord*1000, OCSdata$ycoord*1000,
                   polbuf, 2500)
plot(dens)

#### weigh residuals by the inverse of the sample density ####
w <- 1/(extract(dens, OCSdata[, c("xcoord","ycoord")]*1000)[,2])

ME_wCV   <- weighted.mean(resids, w)
RMSE_wCV <- sqrt(sum(resids^2*w)/sum(w))
muref <- weighted.mean(OCSdata$ocs, w)
MEC_wCV <- 1 - sum(w * resids^2)/sum(w * (OCSdata$ocs - muref)^2)
