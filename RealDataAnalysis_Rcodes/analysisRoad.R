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
library(DescTools) #for Mode
library(grid)
library(gridExtra)
library(MASS)
library(truncnorm)
library(sf)
rm(list=ls(all=TRUE))

path = '/Users/yinziwei/Documents/成功大學(碩)/碩士論文總整理'
# sourceCpp(paste(path, '/thesisCode/cppCode/PathSamplingInvariant_v2.cpp', sep = ""))

IntegrationSpline = function(ExpU, start, end)
{
  fit = smooth.spline(ExpU[, 1], ExpU[, 2])
  smooth = function(x) predict(fit, x)$y
  value = integrate(smooth, start, end)$value
  return(value)
}

#### Read Data ####
load(paste(path, "/thesisCode/Analysis/交通意外分析/data/data.RData", sep = ""))

#### Read map data ####

map = st_read(paste(path, '/thesisCode/Analysis/交通意外分析/data/鄉(鎮、市、區)界線檔(TWD97經緯度)1130719/TOWN_MOI_1130718.shp', sep = ""))
map.main = map[((!(map$COUNTYNAME %in% c('金門縣', '連江縣', '澎湖縣')))&(!(map$TOWNNAME %in% c('琉球鄉', '綠島鄉', '蘭嶼鄉')))),]
which(!(map.main$TOWNID == towns$TOWNID))


cov.of.interest.fixed <- c("rain","temp","pm25","M_F_RAT","P_DEN","A65_A0A14_RAT")
p = length(cov.of.interest.fixed)
t = 12
num.of.obs = n = dim(map.main)[1]
#### accumulate month to season ####
## 如果需要更換時間維度可使用，若使用月份為單位無需使用
# aggregate_to_season_time_dim1 <- function(arr) {
#   stopifnot(length(dim(arr)) == 3)
#   time_len <- dim(arr)[1]
#   stopifnot(time_len %% 3 == 0)  # 必須能被 3 整除
#   
#   # 新的時間維度（季數）
#   num_seasons <- time_len / 3
#   out_array <- array(NA, dim = c(num_seasons, dim(arr)[2], dim(arr)[3]))
#   
#   for (i in 1:num_seasons) {
#     idx <- (3 * (i - 1) + 1):(3 * i)  # 該季的月份索引
#     out_array[i,,] <- apply(arr[idx,,], c(2, 3), mean)  # 對空間位置聚合
#   }
#   
#   return(out_array)
# }
# 
# aggregate_to_season_time_dim1_2d <- function(arr) {
#   stopifnot(length(dim(arr)) == 2)
#   time_len <- dim(arr)[2]
#   stopifnot(time_len %% 3 == 0)  # 必須能被 3 整除
#   
#   # 新的時間維度（季數）
#   num_seasons <- time_len / 3
#   out_array <- array(NA, dim = c(dim(arr)[1],num_seasons))
#   
#   for (i in 1:num_seasons) {
#     idx <- (3 * (i - 1) + 1):(3 * i)  # 該季的月份索引
#     out_array[,i] <- apply(arr[,idx], c(1), sum)  # 對空間位置聚合
#   }
#   
#   return(out_array)
# }

#### Preprocess Data ####
x = array(NA, c(t,p,n))

for(pp in 1:p){
  x[,pp,] = towns[, grepl(cov.of.interest.fixed[pp], colnames(towns))] %>% st_drop_geometry() %>% t()
}
x = as.numeric(x) %>% array(c(t,p,n))
## aggregation of x ##
# x <- aggregate_to_season_time_dim1(x)
##########
x.norm = x
x.mean = apply(x, 2, mean)
x.sd = apply(x, 2, sd)
for(j in 1:length(cov.of.interest.fixed)){
  x.norm[, j, ] = (x[,j,]-x.mean[j])/x.sd[j]
}
x = x.norm

num.in.risk = ref.map[, grepl("^freq", colnames(ref.map))] %>% st_drop_geometry() %>% as.matrix()
num.of.death = ref.map[, grepl("^death_freq", colnames(ref.map))] %>% st_drop_geometry() %>% as.matrix()
num.of.population = ref.map[, grepl("^population_count", colnames(ref.map))] %>% st_drop_geometry() %>% as.matrix()

y = num.of.death + num.in.risk
# y.out1 <- aggregate_to_season_time_dim1_2d(num.in.risk)
# y.out2 <- aggregate_to_season_time_dim1_2d(num.of.death)
# y <- y.out1/(y.out1+y.out2)*100
# y <- y.out2

p_star <- sum(y) / sum(num.of.population/100)
# Eit <- (num.of.death+num.in.risk) * p_star
Eit <- num.of.population/100 * p_star
Eit <- Eit %>% t()
Eit <- Eit/mean(Eit)

#### estimation ####
# t=4
## prelim
W = neighbor.matrix.origin <- matrix(as.numeric(st_touches(map.main,sparse = FALSE)), ncol = nrow(map.main))
ExpU = ExpU_Evaluation(50, 3, 2, 40, W)
# Eit <- matrix(1, nrow = t, ncol = n)
iter.num = 5
est.list <- list()
int.list <- list()
for (run in 1:iter.num) {
  ## initial
  rho =0.5
  lambda = 2
  Q = rho*( diag(c( W%*%matrix(1, num.of.obs, 1) ))-W) + (1-rho)*diag(1, num.of.obs)
  ## phi
  phi = mvrnorm(1, rep(0, num.of.obs), lambda*solve(Q))
  ## psi
  psi = arima.sim(model=list(ar=0.2),n=t, sd=sqrt(0.1))
  ## Covariate
  # PottsImageRandom = matrix(sample(0:2, num.of.obs*p, replace = TRUE), num.of.obs, p)
  PottsImageRandom = matrix(0, num.of.obs, p)
  beta.ini = PottsImageRandom
  beta.ini[PottsImageRandom==0] = 0
  beta.ini[PottsImageRandom==1] = rtruncnorm(sum(PottsImageRandom==1), a=0, b=Inf, mean = 2, sd = 0.1)
  beta.ini[PottsImageRandom==2] = rtruncnorm(sum(PottsImageRandom==2), a=-Inf, b=0, mean = -2, sd = 0.1)
  beta.ini = t(beta.ini)
  
  InititalValues = list(ActImg = PottsImageRandom, beta = beta.ini, 
                        theta=runif(p, 0, 1), phi = phi, psi = psi, rho=rho, lambda=1)
  
  output = Update_Gamma_Beta_Theta(10000, 3, W, y, x, ExpU, InititalValues, Eit, UpdateTheta = TRUE, UpdatePsi = TRUE)
  
  est.list[[run]] <- output
  int.list[[run]] <- InititalValues
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
EstPottsImage.list <- list()
for(i in 1:iter.num){
  EstPottsImage = apply(est.list[[i]]$SimActivationImg, c(1, 2), Mode)
  EstPottsImage.list[[i]] <- EstPottsImage
}

save(est.list, int.list, file = paste(path, "/thesisCode/Analysis/交通意外分析/output_5timesMCMC_smallEit_fixInit.RData", sep = ""))