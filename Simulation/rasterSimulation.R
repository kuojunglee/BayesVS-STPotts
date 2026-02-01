if (!requireNamespace("BayeSTMwPotts", quietly = TRUE)) {
  devtools::install('/Users/yinziwei/Documents/成功大學(碩)/BayeSTMwPotts')
} else {
  print("Package is ready to use!")
}
library(BayeSTMwPotts)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(DescTools)
library(fields)
library(truncnorm)
library(grid)
library(gridExtra)
library(stringr)
rm(list=ls(all=TRUE))
path = '/Users/yinziwei/Documents/成功大學(碩)/碩士論文總整理'
# sourceCpp('/Users/yinziwei/Documents/成功大學(碩)/碩士論文的東西（二）/STBayesVCM/teachCode/ver2/PathSamplingInvariant.cpp')
# sourceCpp('/Users/yinziwei/Documents/成功大學(碩)/碩士論文的東西（二）/程式碼/PathSamplingInvariant_v2.cpp')

IntegrationSpline = function(ExpU, start, end)
{
  fit = smooth.spline(ExpU[, 1], ExpU[, 2])
  smooth = function(x) predict(fit, x)$y
  value = integrate(smooth, start, end)$value
  return(value)
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
#### Settings for Data Generation####
## phi
rho = 0.5
lambda = 2
Q = rho*( diag(c( W%*%matrix(1, num.of.obs, 1) ))-W) + (1-rho)*diag(1, num.of.obs)
phi = mvrnorm(1, rep(0, num.of.obs), lambda*solve(Q))

## psi
zeta = 0.2
psi = arima.sim(model=list(ar=zeta),n=T, sd=sqrt(0.1))

## covariates
SimImage = sample(0:2, num.of.rows*num.of.rows, replace = TRUE)

theta = runif(p, 0.1, 1.0)

PottsImage = matrix(0, num.of.rows*num.of.rows, p)


for(j in 1:p){
  PottsImage[, j] = SimulatePottsModelCompWise(1000, theta[j], 3, SimImage, W)
}

## Setting up the gamma value for each covariate
PottsImage1 = array(0, c(num.of.rows, num.of.rows, 2))
PottsImage1[c(4:8, 13:17), 13:17, 1] = 1
PottsImage1[7:14, 4:8, 1] = 2
PottsImage1[c(4:8, 13:17), 13:17, 2] = 2
PottsImage1[7:14, 4:8, 2] = 1

PottsImage = cbind(c(PottsImage1[, ,1]), c(PottsImage1[, ,2]))
svg("simulation/Raster/Raster_realSetting.svg", width = 20, height = 8)
grid.arrange(
  activation.plot(1, Grid, PottsImage,"Covariate1"),
  activation.plot(2, Grid, PottsImage,"Covariate2"), ncol=2
)
dev.off()
#### Simulate 50 times ####
iter.num = 50
output.list <- list()
true_data.list <- list()
for(simulation in 1:iter.num){
  phi = mvrnorm(1, rep(0, num.of.obs), lambda*solve(Q))
  psi = arima.sim(model=list(ar=zeta),n=T, sd=sqrt(0.1))
  ## covariates
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
  true_data <- list(phi = phi, psi = psi, beta = beta.given, y = y)
  ## Estimation
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
  output = Update_Gamma_Beta_Theta(5000, 3, W, y, x, ExpU, InititalValues, Eit, UpdateTheta = TRUE, UpdatePsi = TRUE)
  output.list[[simulation]] <- output
  true_data.list[[simulation]] <- true_data
}
save(output.list, true_data.list,file = paste(path, "/thesisCode/Simulation/rasterSimulation_50Times.RData", sep = ""))

## est ###
load( paste(path, "/thesisCode/Simulation/rasterSimulation_50Times.RData", sep = ""))
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
EstPottsImage = apply(output.list[[50]]$SimActivationImg, c(1, 2), Mode)

# --- Plots ---
p1 <- activation.plot(1, Grid, PottsImage, "True Cov1")
p2 <- activation.plot(2, Grid, PottsImage, "True Cov2")
p3 <- activation.plot(1, Grid, EstPottsImage, "Est Cov1")
p4 <- activation.plot(2, Grid, EstPottsImage, "Est Cov2")
p.lis <- list(p1, p2, p3, p4)
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
svg(filename = paste(path, "/thesisCode/Simulation/plots/rasterSim.svg", sep = ""), width = 12, height = 10)
grid.draw(combined)
dev.off()

#### Print Acc ####
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
acc.mat <- matrix(NA,ncol = 2,nrow = 50)
for(i in 1:iter.num){
  print(i)
  output1 <- output.list[[i]]
  EstPottsImage = apply(output1$SimActivationImg, c(1, 2), Mode)
  apply(PottsImage==EstPottsImage, 2, mean) %>% print()
  acc.mat[i,] <- apply(PottsImage==EstPottsImage, 2, mean)
}

apply(acc.mat, c(2), mean)
