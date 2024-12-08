# Exercise PERC Uncertainty Friday: spatial aggregation
library(MASS)

# define parameters
dim <- 25
mean <- rep(70,dim)
var <- 30^2
rho <- 0.0  # change this later
threshold <- 100
varcov <- matrix(rho*var, nrow=dim, ncol=dim)
diag(varcov) <- var
varcov

# simulate nsim realisations and compute average of each
nsim <- 2000
data <- mvrnorm(nsim, mean, varcov)  

# plot the first four
plot(data[1,], type="l", col="darkgreen", lwd=2)
lines(data[2,], col="blue", lwd=2)
lines(data[3,], col="red", lwd=2)
lines(data[4,], col="purple", lwd=2)

# compute averages for all and estimate probability > threshold  
average <- rowMeans(data)
hist(average, col="Lightblue")
prob <- mean(average>threshold); prob

# correct answer rho=1 (same as one location)
1-pnorm(100,70,30)