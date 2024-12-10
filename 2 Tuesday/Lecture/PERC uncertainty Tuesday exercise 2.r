# PERC uncertainty Tuesday exercise 2

# Monte Carlo uncertainty propagation USLE model
n <- 5000
set.seed(12345)
R <- rnorm(n,297,72); K <-rnorm(n,0.10,0.05)
L <- rnorm(n,2.13,0.05); S <- rnorm(n,1.17,0.12)
C <- rnorm(n,0.63,0.15); P <- rnorm(n,0.50,0.10)
E <- R*K*L*S*C*P

mean(E); sd(E)
hist(E, col="LightBlue")

# Setting each uncertainty source to its mean, one by one
E_R <- 297*K*L*S*C*P
E_K <- R*0.10*L*S*C*P
E_L <- R*K*2.13*S*C*P
E_S <- R*K*L*1.17*C*P
E_C <- R*K*L*S*0.63*P
E_P <- R*K*L*S*C*0.50

var(E); var(E_R); var(E_K); var(E_L)
var(E_S); var(E_C); var(E_P)

# Uncertainty source contributions
C_R <- 1-var(E_R)/var(E)
C_K <- 1-var(E_K)/var(E)
C_L <- 1-var(E_L)/var(E)
C_S <- 1-var(E_S)/var(E)
C_C <- 1-var(E_C)/var(E)
C_P <- 1-var(E_P)/var(E)

# Can scale to sum to one:
sum <- C_R+C_K+C_L+C_S+C_C+C_P; sum
C_R <- C_R/sum; C_K <- C_K/sum
C_L <- C_L/sum; C_S <- C_S/sum
C_C <- C_C/sum; C_P <- C_P/sum

C_R; C_K; C_L; C_S; C_C; C_P
