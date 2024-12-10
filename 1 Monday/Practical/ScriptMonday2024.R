
# =================================================================
#                       practical UA&SV MONDAY
# =================================================================


TaylorSR <- function(red, nir, s_red, s_nir, rho) {
  # Uncertainty SR by Taylor series method, returns tau.
  # nir and red are the measured reflectances in the two
  # bands; s_red and s_nir are the standard deviations
  # of the measurement errors; rho is the correlation
  # between the two measurement errors.
  tau_sq <- s_nir^2/red^2 +
    s_red^2*nir^2/red^4-2*rho*s_nir*s_red*nir/red^3
  sqrt(tau_sq) # return tau
}


#### call the function TaylorSR with some data ####
TaylorSR(0.1, 0.6, 0.025, 0.03, 0.8)


#### ===================== Exercise 4 ======================== ####

TaylorNDVI <- function(red, nir, s_red, s_nir, rho) {
  # Uncertainty NDVI by Taylor series method, returns tau.
  # Inputs as those for TaylorSR.
  tau_sq <- 4*(s_nir^2*red^2+s_red^2*nir^2 - 
                 2*rho*nir*red*s_red*s_nir)/(nir+red)^4
  tau_sq[tau_sq < 0] <- 0
  sqrt(tau_sq)
}


# call the function with some data
TaylorNDVI(0.1, 0.6, 0.025, 0.03, 0.8)


#### =================== Exercise 5 ========================== ####

red <- 1:500/500 # generate a "red" sequence 0.002, 0.004, ..., 1.00
nir <- 1:500/500 # ditto for the "nir" band

# compute SR for every combination of red and nir using the function outer
SR <- outer(red, nir, function(r, n){n/r})

# compute uncertainty propagation in SR with fixed s_nir and s_red
tau_SR <- outer(red, nir, TaylorSR, s_red=0.02, s_nir=0.03, rho=0.8)
CV_SR <- (tau_SR / SR) # coefficient of variation, see Wikipedia

# Open graphics window
x11(width=15, height=5)

# Split window to display multiple plots
par(mfcol=c(1,3),mai=rep(0.85,4), xaxs="i", yaxs="i", cex=1)

# contour plots with specified levels 
contour(z=SR, levels=c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 50),
        xlab="RED", ylab="NIR", main="SR", labcex=1)
contour(z=tau_SR, levels=c(0.025, 0.05, 0.1, 0.5, 1, 5, 50),
        xlab="RED", ylab="NIR", main="Tau(SR) by Taylor method", labcex=1)
contour(z=CV_SR, levels=c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2), labcex=1,
        xlab="RED", ylab="NIR", main="CV(SR) by Taylor method")



##### ================== Exercise 6 =========================== ####
# compute NDVI for every combination of red and nir using the function outer
NDVI <- outer(red, nir, function(r, n) {(n-r)/(n+r)})

# compute uncertainty propagation in NDVI with fixed s_nir and s_red
tau_NDVI <- outer(red, nir, TaylorNDVI, s_red=0.025, s_nir=0.03, rho=0.8)

CV_NDVI <- tau_NDVI / NDVI  # coefficient of variation, see Wikipedia

# Open graphics window
x11(width=15, height=5)

# Split window to display multiple plots
par(mfcol=c(1,3),mai=rep(0.85,4), xaxs="i", yaxs="i", cex=1)

# contour plots with specified levels (negative = red)
contour(z=NDVI, levels=c(-0.95, -0.5, -0.2, -0.1, -0.05, 
                         0, 0.05, 0.1, 0.2, 0.5, 0.95), 
        col=c(rep("red", 5), "blue", rep("black",5)), xlab="RED", 
        ylab="NIR", main="NDVI", labcex=1)
contour(z=tau_NDVI, levels=c(0.01, 0.015, 0.02, 0.05, 0.1, 0.5),
        xlab="RED", ylab="NIR", main="Tau(NDVI) by Taylor method",
        labcex=1)
contour(z=CV_NDVI, levels=c(-0.5, -0.2, -0.1, -0.05, 0.05, 0.1, 0.2, 0.5, 2),
        xlab="RED", ylab="NIR", main="CV(NDVI) by Taylor method",
        col=c(rep("red", 4), rep("black",5)), labcex=1)

##### =================== Exercise 7 ========================== ####
library(terra)
setwd("C:/PERC/Monday") # set to the folder where downloaded data are
false_color <- rast("FCLorraine.tif")  # reads in a SpatRaster

# close any separate R graphical devices
while (!is.null(dev.list())) dev.off()
    
op <- par(no.readonly = TRUE) # save default graphical settings
plotRGB(false_color, 3, 2, 1, stretch='lin', mar=2)
par(op) # restore default graphical settings

summary(false_color + 0, size=2e5) 
hist(false_color, 1)
hist(false_color, 2)
hist(false_color, 3)

SRim <- app(false_color, fun=function(x) x[3]/x[2], filename="SR.tif", 
             overwrite = T)
plot(SRim, col=gray(1:100/100), main="SR")

tau_SR <- app(false_color, fun=function(x) 
  TaylorSR(x[2], x[3], 0.025, 0.03, 0.8), filename="tau_SR.tif", 
  overwrite = T)

plot(tau_SR, col=gray(1:100/100), main= "Tau(SR)")


##### ===================== Exercise 8 ============================= ####
CV_SR <- app(c(tau_SR, SRim), fun=function(x) {x[1]/x[2]}, 
              filename="CV_SR.tif", overwrite=T)
plot(CV_SR, col=gray(1:100/100), main= "CV(SR)", cex=1)


##### ===================== Exercise 9 ============================= ####
NDVIim <- app(false_color, fun=function(x) (x[3]-x[2])/(x[3]+x[2]), 
               filename="NDVI.tif", overwrite = T)
plot(NDVIim, col=gray(1:100/100), main="NDVI")

tau_NDVI <- app(false_color, fun=function(x) 
  TaylorNDVI(x[2], x[3], 0.025, 0.03, 0.8), 
  filename="tau_NDVI.tif", overwrite = T)
plot(tau_NDVI, col=gray(c(1:20/50, 41:65/100,262:400/400)), 
     main= "Tau(NDVI)")

CV_NDVI <- app(c(tau_NDVI, NDVIim), fun=function(x) {x[1]/x[2]}, 
                   filename="CV_NDVI.tif", overwrite=T)
plot(CV_NDVI, col=gray(c(1:20/50, 41:65/100,262:400/400)), 
     main= "CV(NDVI)")


