if (!requireNamespace("BayeSTMwPotts", quietly = TRUE)) {
  devtools::install('/Users/yinziwei/Documents/成功大學(碩)/BayeSTMwPotts')
} else {
  print("Package is ready to use!")
}
library(BayeSTMwPotts)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(tidyverse)
library(grid)
library(gridExtra)
library(truncnorm)
library(MASS)
library(latex2exp)
rm(list=ls(all=TRUE))
path = '/Users/yinziwei/Documents/成功大學(碩)/碩士論文總整理'
# sourceCpp(paste(path, '/thesisCode/cppCode/PathSamplingInvariant_v2copy.cpp', sep = ""))

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
Grid <- expand.grid(x.easting, x.northing)## 左下到右上
distance <- as.matrix(dist(Grid, method="euclidean"))
K <- nrow(Grid)

W = array(0, c(K,K))
W[distance==1] = 1

ExpU = ExpU_Evaluation(50, 3, 2, 40, W)


#### Potts Image ####
#2d PottsImage
PottsImage2d = array(0, c(num.of.rows, num.of.cols, 2))
PottsImage2d[c(4:8, 13:17), 13:17, 1] = 1
PottsImage2d[7:14, 4:8, 1] = 2
PottsImage2d[c(4:8, 13:17), 13:17, 2] = 2
PottsImage2d[7:14, 4:8, 2] = 1
##1d PottsImage
PottsImage = cbind(c(PottsImage2d[,,1]), c(PottsImage2d[,,2]))
image(matrix(PottsImage[,1], num.of.rows, num.of.cols))
image(matrix(PottsImage[,2], num.of.rows, num.of.cols))

#### Generate Data ####
iter.num <- 50
output1.list <- list()
output2.list <- list()
data.list <- list()
for(iter in 1:iter.num){
  ## Data Generation
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
  
  #1. theta neq 0
  ## Estimation
  IntImage = array(sample(0:2, size = (n*p), replace = TRUE), c(n,p))
  IntBeta = dplyr::case_when(
    as.vector(IntImage)==0 ~ 0,
    as.vector(IntImage)==1 ~ 5,
    as.vector(IntImage)==2 ~ -5,
    TRUE ~ NA
  ) |> array(c(p,n))
  InititalValues = list(ActImg = IntImage, beta = beta.given,
                        theta=c(0., 0.), phi = phi, psi = rep(0, T), rho=rho, lambda=1)
  # Eit = matrix(mean(y), nrow = T, ncol = n)
  Eit = matrix(1, nrow = T, ncol = n)
  output = Update_Gamma_Beta_Theta(5000, 3, W, y, x, ExpU, InititalValues, Eit, UpdateTheta = TRUE, UpdatePsi = TRUE)
  
  #2. theta = 0
  
  output2 = Update_Gamma_Beta_Theta(5000, 3, W, y, x, ExpU, InititalValues, Eit, UpdateTheta = TRUE, UpdatePsi = FALSE)
  
  output1.list[[iter]] <- output
  output2.list[[iter]] <- output2
  data.list[[iter]] <- list(y = y, x = x)
}

save(output1.list, output2.list, data.list, file = paste(path, "/thesisCode/Simulation/temporalComp_50rep.RData", sep = ""))
## DIC
source(paste(path,"/thesisCode/Simulation/DIC_calculation.R", sep = ""))
load(paste(path, "/thesisCode/Simulation/temporalComp_50rep.RData", sep = ""))

dic.table <- data.frame(NA,ncol = 2)
for(iter in 1:iter.num){
  print(iter)
  output <- output1.list[[iter]]
  output2 <- output2.list[[iter]]
  y <- data.list[[iter]]$y
  x <- data.list[[iter]]$x
  dic.table[iter, 1] <- calculate_DIC(output$beta_samples, output$phi_samples, output$psi_samples, Eit, y, x)
  dic.table[iter, 2] <- calculate_DIC(output2$beta_samples, output2$phi_samples, output2$psi_samples, Eit, y, x)
}

colMeans(dic.table)
dic.table[,1] > dic.table[,2]

## consequences
output <- output1.list[[50]]
output2 <- output2.list[[50]]
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
activation.plot <- function(jj,Grid,z.est.mat,title.string){
  z.Data <- Grid %>% dplyr::mutate(value = as.factor(z.est.mat[,jj]))
  ## z.Data has to include the Grid map and the value of each grid
  ggplot(z.Data, aes(x = Var1,y = Var2))+ 
    geom_tile(aes(fill = value),colour = "grey50", size = 0.1, alpha = 0.5)+
    scale_shape_manual(name = "category", values = c(0, 1, 2)) + 
    scale_fill_manual(name = "category", values = c('0' = "grey95", '1' = "#B40404", '2' = "#0B6121")) +
    scale_x_discrete(breaks=c(1, 5, 10, 15, 20))+
    scale_y_discrete(breaks=c(1, 5, 10, 15, 20))+
    theme(legend.text=element_text(size=7)) + 
    ggtitle(title.string) + xlab("x") + ylab("y")
}
EstPottsImage = apply(output$SimActivationImg, c(1, 2), Mode)
p1 <- activation.plot(1,Grid,EstPottsImage, expression("cov1 with" ~theta))
p2 <- activation.plot(2,Grid,EstPottsImage, expression("cov2 with" ~theta))

EstPottsImage2 = apply(output2$SimActivationImg, c(1, 2), Mode)
p3 <- activation.plot(1,Grid,EstPottsImage2, expression("cov1 without" ~theta))
p4 <- activation.plot(2,Grid,EstPottsImage2, expression("cov2 without" ~theta))

p.lis <- list(p1,p2,p3,p4)
#### ACA ####
acc.mat1 <- matrix(NA,ncol = 2,nrow = iter.num)
acc.mat2 <- matrix(NA,ncol = 2,nrow = iter.num)
accTemporal1.count <- 0
accTemporal2.count <- 0
for(i in 1:iter.num){
  print(i)
  output1 <- output1.list[[i]]
  EstPottsImage = apply(output1$SimActivationImg, c(1, 2), Mode)
  # apply(PottsImage==EstPottsImage, 2, mean) %>% print()
  acc.mat1[i,] <- apply(PottsImage==EstPottsImage, 2, mean)
  
  output2 <- output2.list[[i]]
  EstPottsImage2 = apply(output2$SimActivationImg, c(1,2), Mode)
  # apply(PottsImage==EstPottsImage2, 2, mean) %>% print()
  acc.mat2[i, ] <- apply(PottsImage==EstPottsImage2, 2, mean)
  # print(acc.mat1[i,1] > acc.mat2[i, 1])
  # print(acc.mat1[i,2] > acc.mat2[i, 2])
  if(acc.mat1[i,1] > acc.mat2[i, 1]) accTemporal1.count = accTemporal1.count + 1
  if(acc.mat1[i,2] > acc.mat2[i, 2]) accTemporal2.count = accTemporal2.count + 1
}

apply(acc.mat1, c(2), mean)
apply(acc.mat2, c(2), mean)
accTemporal1.count;accTemporal2.count
# --- Function to extract legend ---
get_legend <- function(myplot) {
  tmp <- ggplotGrob(myplot)
  gtable_filter(tmp, "guide-box")
}

# --- Extract one legend ---
legend <- get_legend(p1)

# Remove legends from all plots
p.lis <- lapply(p.lis, function(p) p + theme(legend.position="none"))

# Build combined layout: plots + legend
combined <- arrangeGrob(
  do.call("arrangeGrob", c(p.lis, ncol = 2)),
  legend,  # second row
  ncol = 2,
  widths = c(10, 1)
)

svg(filename = paste(path, "/thesisCode/Simulation/plots/thetaSim.svg", sep = ""), width = 12, height = 10)
grid.draw(combined)
dev.off()