library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(tidyverse)
library(grid)
library(gridExtra)
library(truncnorm)
rm(list=ls(all=TRUE))

sourceCpp('/Users/yinziwei/Documents/成功大學(碩)/碩士論文的東西（二）/程式碼/PathSamplingInvariant_v3.cpp')

IntegrationSpline = function(ExpU, start, end)
{
  fit = smooth.spline(ExpU[, 1], ExpU[, 2])
  smooth = function(x) predict(fit, x)$y
  value = integrate(smooth, start, end)$value
  return(value)
}

#### Set up ####
num.of.rows = 20
num.of.cols = 20
n = num.of.obs = num.of.rows*num.of.cols
T = 9
p = num.of.covariates = 2
x.easting <- 1:num.of.rows
x.northing <- 1:num.of.cols
Grid <- expand.grid(x.easting, x.northing)
distance <- as.matrix(dist(Grid, method="euclidean"))
K <- nrow(Grid)

W = array(0, c(K,K))
W[distance==1] = 1

ExpU = ExpU_Evaluation(50, 3, 2, 40, W)

#### Potts Image ####
PottsImage = matrix(0, n, p) #1d
#2d PottsImage
PottsImage2d = array(0, c(num.of.rows, num.of.cols, 2))
PottsImage2d[c(4:8, 13:17), 13:17, 1] = 1
PottsImage2d[7:14, 4:8, 1] = 2
PottsImage2d[c(4:8, 13:17), 13:17, 2] = 2
PottsImage2d[7:14, 4:8, 2] = 1

PottsImage = cbind(c(PottsImage2d[,,1]), c(PottsImage2d[,,2]))
image(matrix(PottsImage[,1], num.of.rows, num.of.cols))
image(matrix(PottsImage[,2], num.of.rows, num.of.cols))

#### Generate Data ####
## phi
rho = 0.5
lambda = 2
Q = rho*( diag(c( W%*%matrix(1, num.of.obs, 1) ))-W) + (1-rho)*diag(1, num.of.obs)
phi = mvrnorm(1, rep(0, num.of.obs), lambda*solve(Q))

## psi
zeta = 0.2
psi = arima.sim(model=list(ar=zeta),n=T, sd=sqrt(0.1))

## variables
beta.given = PottsImage
beta.given[PottsImage[, 1]==0, 1] = 0
beta.given[PottsImage[, 1]==1, 1] = rtruncnorm(sum(PottsImage[, 1]==1), a=0, b=Inf, mean = 2, sd = 0.1)
beta.given[PottsImage[, 1]==2, 1] = rtruncnorm(sum(PottsImage[, 1]==2), a=-Inf, b=0, mean = -1, sd = 0.1)

beta.given[PottsImage[, 2]==0, 2] = 0
beta.given[PottsImage[, 2]==1, 2] = rtruncnorm(sum(PottsImage[, 2]==1), a=0, b=Inf, mean = 1, sd = 0.1)
beta.given[PottsImage[, 2]==2, 2] = rtruncnorm(sum(PottsImage[, 2]==2), a=-Inf, b=0, mean = -2, sd = 0.1)
beta.given = t(beta.given)
x = array(rnorm(n*p*T), c(T, p, n))
eta = matrix(0, n, T)
for(i in 1:n){
  eta[i, ] = x[, , i]%*%beta.given[, i, drop=F] + psi
}
eta = eta+phi	
## Generate y
y = apply(exp(eta), c(1, 2), rpois, n=1)

#### Estimation ####
### 3 cluster
IntImage = array(sample(0:2, size = (n*p), replace = TRUE), c(n,p))
IntBeta = dplyr::case_when(
  as.vector(IntImage)==0 ~ 0,
  as.vector(IntImage)==1 ~ 5,
  as.vector(IntImage)==2 ~ -5,
  TRUE ~ NA
) |> array(c(p,n))
InititalValues = list(ActImg = IntImage, beta = beta.given, 
                      theta=c(0., 0.), phi = phi, psi = psi, rho=rho, lambda=1)
# Eit = matrix(mean(y), nrow = T, ncol = n)
Eit = matrix(1, nrow = T, ncol = n)
G = 2
alpha_gj_range = array(NA, c(G, num.of.covariates))
alpha_gj_range[1,] <- c(-Inf,0); alpha_gj_range[2,] <- c(0, Inf)
HyperParam = list(alpha_gj_range = alpha_gj_range)
output1 = Update_Gamma_Beta_Theta(5000, 3, W, y, x, ExpU, InititalValues, HyperParam, Eit)
output1$acceptance_rate_beta_gamma %>% dim
output1$acceptance_rate_theta
### 2 cluster
sourceCpp('/Users/yinziwei/Documents/成功大學(碩)/碩士論文的東西（二）/IsingVersion/Codes/rcppCode/SpVCisingModel.cpp')
InititalValues = list(ActImg = matrix(0,p,n), beta = matrix(0,p,n), 
                      theta=c(0., 0.), phi = phi, psi = psi, rho=rho, lambda=1)

output2 = Update_Gamma_Beta_Theta_ising(5000, 2, W, y, x, ExpU, InititalValues, Eit)
#### DIC ####
dim(output$beta_samples)
dim(output$phi_samples)
dim(output$psi_samples)
dim(y); dim(x); dim(Eit)
source("Simulation/DIC_calculation.R")
calculate_DIC(output1$beta_samples, output1$phi_samples, output1$psi_samples, Eit, y, x)
calculate_DIC(output2$beta_samples, output2$phi_samples, output2$psi_samples, Eit, y, x)