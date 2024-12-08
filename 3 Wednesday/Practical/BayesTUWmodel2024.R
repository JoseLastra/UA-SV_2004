# ********************************************************************
# Uncertainty 2024: Bayesian calibration selected parameters TUWmodel,
# see https://cran.r-project.org/web/packages/TUWmodel/index.html
# ************ For education/demonstration purposes only *************
# ****** 2016, 2020, 2022 Sytze de Bruin, Wageningen University ******
# ** 2018, Alexandre Wadoux & Sytze de Bruin, Wageningen University **
# ********************************************************************

#### Load and attach "TUWmodel" library ####
# Check whether it has been installed; if not install it
options(warn=-1)
if (!require(TUWmodel)){
  install.packages("TUWmodel")
  library(TUWmodel)
}
options(warn=0)

#### Load Thur data, similar to those used in DOI 10.7717/peerj.9558 ####
#    you may have to set the working directory with setwd() 
setwd("C:/PERC/Wednesday")
prec <- as.numeric(readLines("agERA5prec.txt"))
airt <- as.numeric(readLines("temp.txt"))
epot <- as.numeric(readLines("evap.txt"))
runoff <- as.numeric(readLines("runoff.txt"))
dates <- as.Date(readLines("dates.txt"), "%Y%m%d")

# **********************************************************
#### ****    Default run with initial values       **** ####
# **********************************************************

defaultPars <- c(1.0, 1.5, 2, 0, -0.5, 1.0, 60, 2.52, 0.4, 
                 5.0, 120, 60, 6.0, 9, 15)
initState <- c(60, 0, 37, 300)

simDefault <- TUWmodel(prec, airt, epot, 1, defaultPars, initState)
plot(dates, runoff, type="l", xlab="", ylab="Discharges [mm/day]", lwd=1.5)
lines(dates, simDefault$q, col=2, lwd=1.5)
legend("topright", legend=c("Observations","Simulations"), 
       col=c(1,2), lty=1, bty="n", lwd=1.5)


#### Prepare time series for calibration &validation ####
# calIds: 01/01/2004 - 31/12/2008 --> the calibration data
calIds <- 1:1827
runoff.cal <- runoff[calIds]
prec.cal <- prec[calIds]
airt.cal <- airt[calIds]
epot.cal <- epot[calIds]

# valIds: 01/01/2009 - 31/12/2011 --> the validation data
valIds <- 1828:2922
runoff.val <- data.matrix(runoff[valIds])
prec.val <- data.matrix(prec[valIds])
airt.val <- data.matrix(airt[valIds])
epot.val <- data.matrix(epot[valIds])

# RMSE of the (relative) residuals: validation period, default run
relResidsDflt <- 100 *
  (as.numeric(simDefault$q)[valIds] - runoff.val)/runoff.val
summary(relResidsDflt)
relRMSE_dflt <- sqrt(mean(relResidsDflt^2))  # RMSE(relResids)
cat("RMSE(relative residuals) default run: ", round(relRMSE_dflt,0), "%", sep= "")

# **********************************************************
#### ******       DEFINE AND PLOT PRIORS        ******* #### 
# Selected parameters: k1, lsuz, cperc, croute + error pars.
# **********************************************************
ipar <- c(10, 12, 13, 15) # indices of the selected model pars

# labels for axis plotting
labs <- c("k1 day", "lsuz mm", "cperc mm/day",  
          expression(croute ~ day^2 ~ mm^-1), 
          expression(beta[1]), 
          expression(sigma[delta]^2), 
          expression(sigma[eta]^2))

# Function for priors as used in Parajka et al., 2007; see fig 2.
# Let's call these "scaled beta distributions"
dBetaPrior <- function(x, shape1, shape2, min, max){
  xx <- (x-min)/(max-min)
  dbeta(xx, shape1, shape2)/(max-min)
}

# Parameters prior distributions of selected parameters
pModpars <- rbind(c(2.0, 4.0, 2, 30.0), c(3.0, 3.0, 1.0, 100.0), 
                  c(2.0, 4.0, 0, 8.0), c(1.05, 1.05, 0, 75.0))

# Uniform priors of the additional error model parameters
beta1 <- c(0,1)
sigma.delta <- c(0,1)
sigma.eta <- c(0,1)

# Plot prior distribution k1
plot(4:60/2, dBetaPrior(4:60/2, 2.0, 4.0, 2, 30), ty="l", 
     xlab=labs[1], ylab="Density")

# Plot prior distribution lsuz;
plot(0:101, dBetaPrior(0:101, 3.0, 3.0, 1, 100), ty="l", 
     xlab=labs[2], ylab="Density")

# Plot prior distribution cperc;
plot(0:80/10, dBetaPrior(0:80/10, 2.0, 4.0, 0, 8), ty="l", 
     xlab=labs[3], ylab="Density")

# prior distribution croute
plot(0:75, dBetaPrior(0:75, 1.05, 1.05, 0, 75), ty="l", 
     xlab=labs[4], ylab="Density")

# Plot prior distribution beta_1
plot(c(0.0,0,1, 1), c(0, 1, 1, 0), ty="l", ylim=c(0, 1.5),
     xlab=labs[5], ylab="Density")

# Plot prior distribution sigma_delta
plot(c(0.0,0,1, 1), c(0, 1, 1, 0), ty="l", ylim=c(0, 1.5),
     xlab=labs[6], ylab="Density")

# Plot prior distribution sigma_eta
plot(c(0.0,0,1, 1), c(0, 1, 1, 0), ty="l", ylim=c(0, 1.5),
     xlab=labs[7], ylab="Density")


# ******************************************************
# ******                                          ******
#### ***     Bayesian calibration functions     *** ####
# ******                                          ******
# ******************************************************

#### Log likelihood function
loglik = function(params, init, beta1, sigma.delta, sigma.eta, 
                  prec, airt, epot, runoff){
  # Computes the log likelihood of the data given sampled parameters 

  # the expectation of the structural uncertainty and discharge measurement 
  # error terms is forced to 1 to avoid bias, hence
  beta0 <- -0.5 * ((1-beta1)/(1-beta1^2))*sigma.delta
  
  pred.log <- as.numeric(log(TUWmodel(prec, airt, epot, 1, params, init)$q))
  
  Y <- log(runoff) - pred.log

  second.term <- c()
  negloglik <- c()
  first.term <- c()
  # Kalman approach Eqs. 9-18 Wadoux et al. 2020
  hat.vaeps.min <- c()
  hat.vaeps.plus <- c()
  sigma.plus <- c()
  sigma.min <-c()
  k <- c()
  d <- c()
  for (t in 1:length(runoff)){
    if (t==1){
      hat.vaeps.plus[1] <- 10
      sigma.plus[1] <- 10
      
    }else{
      hat.vaeps.min[t] <- beta0 + beta1*hat.vaeps.plus[t-1] 
      sigma.min[t] <- (beta1^2)*sigma.plus[t-1] + sigma.delta
      k[t] <- sigma.min[t]/(sigma.min[t]+ sigma.eta)
      d[t] <- Y[t]- hat.vaeps.min[t]
      hat.vaeps.plus[t] <- hat.vaeps.min[t] + k[t]*d[t]
      sigma.plus[t] <- (1-k[t])*sigma.min[t]
      
      tt1 <- 0.5*(log(sigma.min[t]+sigma.eta))
      
      tt2 <- 0.5*(((Y[t]-hat.vaeps.min[t])^2)/(sigma.min[t]+sigma.eta))
      
      # likelihood function
      negloglik <- append(negloglik, -((1/2)*log(2*pi) + tt1 + tt2))
    }
  }
  return(sum(negloglik))
}


metropolis <- function(N, params, init, ipar, pModpars, beta1, sigma.delta, 
                       sigma.eta, prop0, step, prec, airt, epot, runoff){
  # Metropolis algorithm; N specifies the number of runs, including burning in;
  # params, prec, airt, and epot are TUWmodel inputs while runoff is reference 
  # runoff; pModpars specifies the priors of model parameters identified by ipar; 
  # beta1 and sigma.delta are the two free parameters of the multiplicative 
  # first order auto regressive structural error model and sigma.eta is the prior
  # of measurement error; prop0 is the initial proposal and step defines the jump
  # parameters for exploring parameter space. 
  # The function returns a sample of size N of the model parameters identified by 
  # ipar and prints the acceptance rate (after burning in).
  
  if(N < 1)
    stop("N should be larger than 0")
  
  np <- length(ipar) # number of parameters
  ns1 <- np + 1      # beta1 structural error slope
  ns2 <- ns1 + 1     # sigma.delta structural error variance
  ns3 <- ns2+1       # sigma.eta measurement error variance

  
  # log(~probability) initial proposal
  logdenom <- sum(log(dBetaPrior(prop0[1:np], pModpars[,1], pModpars[,2], 
                                 pModpars[,3], pModpars[,4]))) + 
    dunif(prop0[ns1], beta1[1], beta1[2], log=T) + 
    dunif(prop0[ns2], sigma.delta[1], sigma.delta[2], log=T) +
    dunif(prop0[ns3], sigma.eta[1], sigma.eta[2], log=T) +
    loglik(replace(params, ipar, prop0[1:np]), init, prop0[ns1], prop0[ns2],
           prop0[ns3], prec, airt, epot, runoff)
  
  if (logdenom == -Inf)
    stop("Initial proposal has zero probability")
  
  subMetropolis <- function(jump){
    # Actual Metropolis sampling using predefined jumps
    
    propt <- prop0 + jump

    # everything on log scale to avoid numerical problems
    lognum <- sum(log(dBetaPrior(propt[1:np], pModpars[,1], pModpars[,2], 
                                 pModpars[,3], pModpars[,4]))) + 
      dunif(propt[ns1], beta1[1], beta1[2], log=T) + 
      dunif(propt[ns2], sigma.delta[1], sigma.delta[2], log=T)+
      dunif(propt[ns3], sigma.eta[1], sigma.eta[2], log=T)+
      loglik(replace(params, ipar, propt[1:np]), init, propt[ns1], propt[ns2],
             propt[ns3], prec, airt, epot, runoff)
    
    a <- lognum - logdenom
    accept <- 0
    if (!is.nan(a) & ((a>0) | (a>log(runif(1,0,1))))) {
      prop0 <<- propt
      accept <- 1
      logdenom <<- lognum
    }
    c(prop0, accept)
  }
  
  apply_pb <- function(X, MARGIN, FUN, ...)
  {
    env <- environment()
    pb_Total <- sum(dim(X)[MARGIN])
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total,
                         style = 3)
    
    wrapper <- function(...)
    {
      curVal <- get("counter", envir = env)
      assign("counter", curVal +1 ,envir= env)
      setTxtProgressBar(get("pb", envir= env),
                        curVal +1)
      FUN(...)
    }
    res <- apply(X, MARGIN, wrapper, ...)
    close(pb)
    res
  }
  
  # apply subMetropolis
  jumps  <- replicate(N, rnorm(ns3, 0, step))
  result <- apply_pb(jumps, 2, subMetropolis)
  
  # transpose result
  result <- t(result)
  
  # report acceptance rate after burning in; aim for approx 0.25 
  # (Roberts et al. 1997)
  cat("acceptance rate:", round(mean(result[-(1:500),ns3+1]),3))
  
  return(result[,-(ns3+1)])
  # end Metropolis algorithm
}


# **********************************************************
#### ****   Preparation & calling the functions    **** ####
# **********************************************************

# initial proposal selected parameters & error model pars
prop0 <- c(defaultPars[ipar], 0.03, 0.005,0.005)

# step parameters 
step <- c(0.3, 1.2, 0.2, 2.0, 0.02,0.002,0.0015)

#### calling the Metropolis sampling function ####
set.seed(1234567)
N <- 10000
#### BELOW: replace 2000 by N once the acceptance rate ~ 0.25 !!!! ####
s <- metropolis(2000, defaultPars, initState, ipar, pModpars, beta1, sigma.delta, 
                sigma.eta, prop0, step, prec.cal, airt.cal, epot.cal, runoff.cal) 

# trace plots
oldpar <- par(no.readonly = TRUE)
par(mai=rep(1.5,4))
plot(s[,1], type = "l", ylab=labs[1])
plot(s[,2], type = "l", ylab=labs[2])
plot(s[,3], type = "l", ylab=labs[3])
plot(s[,4], type = "l", ylab=labs[4])
plot(s[,5], type = "l", ylab=labs[5])
plot(s[,6], type = "l", ylab=labs[6])
plot(s[,7], type = "l", ylab=labs[7])
par(oldpar)


# subsampling (thinning) after burning in (first set N!!!)
nsub <- (N - 1500)/50
ssub <- s[1500+(1:nsub)*50,]

hist(ssub[,1], freq=F, xlim=c(2,30), xlab=labs[1], main="")
lines(4:60/2, dBetaPrior(4:60/2, 2.0, 4.0, 2, 30), col="red")
legend("topright", c("prior", "posterior"), col=c("red", "black"), 
       lty=1, bty="n")

hist(ssub[,2], freq=F, xlim=c(1,100), xlab=labs[2], main="")
lines(0:100, dBetaPrior(0:100, 3.0, 3.0, 1, 100), col="red")
legend("topright", c("prior", "posterior"), col=c("red", "black"), 
       lty=1, bty="n")

hist(ssub[,3], freq=F, xlim=c(0,8), xlab=labs[3], main="")
lines(0:80/10, dBetaPrior(0:80/10, 2.0, 4.0, 0, 8), col="red")
legend("topright", c("prior", "posterior"), col=c("red", "black"), 
       lty=1, bty="n")

hist(ssub[,4], freq=F, xlim=c(0,80), xlab=labs[4], main="")
lines(0:75, dBetaPrior(0:75, 1.05, 1.05, 0, 75), col="red")
legend("topleft", c("prior", "posterior"), col=c("red", "black"), 
       lty=1, bty="n")

hist(ssub[,5], freq=F, xlim=c(0,1), xlab=labs[5], main="")
lines(c(0.0,0,1, 1), c(0, 1, 1, 0), col="red")
legend("topleft", c("prior", "posterior"), col=c("red", "black"), 
       lty=1, bty="n")

hist(ssub[,6], freq=F, xlim=c(0,1), xlab=labs[6], main="")
lines(c(0.0,0,1, 1), c(0, 1, 1, 0), col="red")
legend("toprigh", c("prior", "posterior"), col=c("red", "black"), 
       lty=1, bty="n")

hist(ssub[,7], freq=F, xlim=c(0,1), xlab=labs[7], main="")
lines(c(0.0,0,1, 1), c(0, 1, 1, 0), col="red")
legend("toprigh", c("prior", "posterior"), col=c("red", "black"), 
       lty=1, bty="n")


colnames(ssub) <- c("k1", "lsuz", "cperc", "croute", 'beta1', 'sigma.delta', 'sigma.eta')
summary(ssub)

# scatterplots and correlations parameters
plot(as.data.frame(ssub[,1:4]))
cor(ssub[,1:4])


# *********************************************************
#### ****            Uncertainty analysis       ****** ####
# Series of runs computed by sampling from the posteriors *
# ******               for the validation period     ******
# *********************************************************

# TUWmodel runs with posteriors
outTUW <- numeric()
for (i in 1:nsub){ # nsub = sample size from the parameter pdf
  tmpPars <- replace(defaultPars, ipar, ssub[i, 1:4])
  outTUW <- rbind(outTUW,TUWmodel(prec, airt, epot, 
                                  param=tmpPars, initState)$q)
}

# multiplicative error # all considered components
multAR <- numeric()    # auto regressive model structural error
multMS <- numeric()    # measurement error Q
for (i in 1:nsub){
  AR.error <- c() 
  for (t in 1:length(runoff.val)){
    if (t==1){AR.error[t] <- 0.3 }else{
      beta0 <- -0.5 * ((1-ssub[i,5])/(1-ssub[i,5]^2))*ssub[i,6]
      AR.error[t] <- beta0 +  ssub[i,5]*(AR.error[t-1]) + 
        rnorm(1, 0, sqrt(ssub[i,6]))
    }}
  multAR <- rbind(multAR, exp(AR.error))
  multMS <- rbind(multMS, exp(rnorm(length(prec.val), 0, sqrt(ssub[i, 7]))))
}


# all considered uncertainty components included
qs <- outTUW[, valIds] * multAR * multMS

# mean over sample
mu  <- apply(qs, 2, mean, na.rm=T)

# visual comparison of parameterizations
x11(width=10, height=7)
plot(as.Date(dates[valIds]), runoff[valIds], type="l", lwd=1.5, 
     col="grey50", xlab="", ylab="Discharges [mm/day]", ylim=c(0,30))
lines(dates[valIds], simDefault$q[valIds], col="blue", lwd=1.5)
lines(dates[valIds], mu, col="red", lwd=1.5)
legend("topleft", legend=c("Observations","Original parameterization", 
                            "Updated parameterization"), 
       col=c("grey50", "blue", "red"), lty=1, lwd=1.5)


# statistics residuals over validation period
# RMSE of the relative residuals: validation period, default run
relResidsClt <- 100 * (mu - runoff.val)/runoff.val
summary(relResidsClt)
relRMSE_clt <- sqrt(mean(relResidsClt^2))  # RMSE(relResids)

cat("RMSE(relative residuals) calibrated parameters: ", round(relRMSE_clt,0), 
    "%", sep = "")

cat("RMSE(relative residuals) default run: ", round(relRMSE_dflt,0), "%", sep= "")

# *************************************
####   UNCERTAINTY CONTRIBUTIONS   ####
# *************************************

# 1 - Model parameters, structure & measurement error
# Compute and plot 5 and 95 percentiles and mean for predictions
q5 <- apply(qs, 2, quantile, 0.05, na.rm=T)
q95 <- apply(qs, 2, quantile, 0.95, na.rm=T)
qPolx <- dates[valIds]

x11(width=10, height=7)
plot(c(min(dates[valIds]), max(dates[valIds])), y=c(NA,NA), ty="n",  
     ylim=c(0,50), xlab="Date", ylab="Discharge [mm/day]", 
     main = "Model parameters, structure & measurement error")

polygon(c(qPolx, rev(qPolx)),c(q5,rev(q95)),
        col = "grey60", border = FALSE)

lines(dates[valIds], runoff.val, col="red")
lines(dates[valIds], mu, col="black")

legend("topright", c("predicted","90% prediction interval", "measured"), 
       col=c("black", "grey60", "red"), lwd=c(1,5,1))


# 2 - Only model parameter uncertainty
mu_mult2 <- mean((multAR * multMS))
mu2 <- apply(outTUW[, valIds] * mu_mult2, 2, mean, na.rm=T)
q5_2 <- apply(outTUW[, valIds] * mu_mult2, 2, quantile, 0.05, na.rm=T)
q95_2 <- apply(outTUW[, valIds] * mu_mult2, 2, quantile, 0.95, na.rm=T)

# the plotting
x11(width=10, height=7)
plot(c(min(dates[valIds]), max(dates[valIds])), y=c(NA,NA), ty="n",  
     ylim=c(0,50), xlab="Date", ylab="Discharge [mm/day]", 
     main = "Only model parameters")

polygon(c(qPolx, rev(qPolx)),c(q5_2,rev(q95_2)),
        col = "grey60", border = FALSE)

lines(dates[valIds], runoff.val, col="red")
lines(dates[valIds], mu2, col="black")

legend("topright", c("predicted","90% prediction interval", "measured"), 
       col=c("black", "grey60", "red"), lwd=c(1,5,1))


# 3 - Model parameter uncertainty and measurement error
mu_mult3 <- mean(multAR)
mu3 <- apply(outTUW[, valIds] * multMS * mu_mult3, 2, mean, na.rm=T)
q5_3 <- apply(outTUW[, valIds] * multMS * mu_mult3, 2, quantile, 0.05, na.rm=T)
q95_3 <- apply(outTUW[, valIds] * multMS * mu_mult3, 2, quantile, 0.95, na.rm=T)

# the plotting
x11(width=10, height=7)
plot(c(min(dates[valIds]), max(dates[valIds])), y=c(NA,NA), ty="n",  
     ylim=c(0,50), xlab="Date", ylab="Discharge [mm/day]", 
     main = "Model parameters & measurement error")

polygon(c(qPolx, rev(qPolx)),c(q5_3,rev(q95_3)),
        col = "grey60", border = FALSE)

lines(dates[valIds], runoff.val, col="red")
lines(dates[valIds], mu3, col="black")

legend("topright", c("predicted","90% prediction interval", "measured"), 
       col=c("black", "grey60", "red"), lwd=c(1,5,1))



