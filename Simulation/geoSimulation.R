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
sourceCpp(paste(path, '/thesisCode/cppCode/PathSamplingInvariant_v2.cpp', sep = ""))

IntegrationSpline = function(ExpU, start, end)
{
  fit = smooth.spline(ExpU[, 1], ExpU[, 2])
  smooth = function(x) predict(fit, x)$y
  value = integrate(smooth, start, end)$value
  return(value)
}

#### Read Data ####
load(paste(path, "/thesisCode/Analysis/交通意外分析/data/data.RData", sep = ""))
urbanize.data <- read.csv(paste(path,"/thesisCode/Simulation/urbanize.csv",sep = ""))
#### Read map data ####

map = st_read(paste(path, '/thesisCode/Analysis/交通意外分析/data/鄉(鎮、市、區)界線檔(TWD97經緯度)1130719/TOWN_MOI_1130718.shp', sep = ""))
map.main = map[((!(map$COUNTYNAME %in% c('金門縣', '連江縣', '澎湖縣')))&(!(map$TOWNNAME %in% c('琉球鄉', '綠島鄉', '蘭嶼鄉')))),]
which(!(map.main$TOWNID == towns$TOWNID))
map.townName <- paste(map.main$COUNTYNAME, map.main$TOWNNAME, sep = "")
idx <- sapply(map.townName, function(x){
  which(x == urbanize.data$TOWNNAME)
})
urbanize.data.reorder = urbanize.data[unlist(idx), ]

cov.of.interest.fixed <- c("temp","pm25")
p = length(cov.of.interest.fixed)
t = T = 12
num.of.obs = n = dim(map.main)[1]

#### Preprocess Data ####
x = array(NA, c(t,p,n))

for(pp in 1:p){
  x[,pp,] = towns[, grepl(cov.of.interest.fixed[pp], colnames(towns))] %>% st_drop_geometry() %>% t()
}
x = as.numeric(x) %>% array(c(t,p,n))
##########
x.norm = x
x.mean = apply(x, 2, mean)
x.sd = apply(x, 2, sd)
for(j in 1:length(cov.of.interest.fixed)){
  x.norm[, j, ] = (x[,j,]-x.mean[j])/x.sd[j]
}
x = x.norm

#### estimation ####
# t=4
## prelim
W = neighbor.matrix.origin <- matrix(as.numeric(st_touches(map.main,sparse = FALSE)), ncol = nrow(map.main))
ExpU = ExpU_Evaluation(50, 3, 2, 40, W)

city.index = urbanize.data.reorder$Urbanization %in% c(1, 2)
suburb.index = urbanize.data.reorder$Urbanization %in% c(7)
PottsImage = matrix(0, num.of.obs, p)
PottsImage[city.index, 1] = 1
PottsImage[suburb.index, 1] = 2

PottsImage[city.index, 2] = 2
PottsImage[suburb.index, 2] = 1

activation.plot = function(index.num,map.main,z.est.mat,title.string){
  map.main.use <- st_simplify(map.main, dTolerance = 1000)
  pl <- map.main.use %>% 
    ggplot(aes(fill = factor(z.est.mat[, index.num]))) + 
    # geom_sf() +
    geom_sf(size = 0.5, alpha = 0.5) +
    scale_shape_manual(name = "category", values = c(0, 1, 2)) + 
    scale_fill_manual(name = "category", values = c('0' = "grey95", '1' = "#B40404", '2' = "#0B6121")) +
    # lims(x = c(120, 122), y = c(21,26)) +
    scale_x_continuous(breaks = seq(120,122), limits = c(120,122))+
    scale_y_continuous(breaks = seq(21,26), limits = c(21.5,25.5))+
    coord_sf() +
    ggtitle(cov.of.interest.fixed[index.num])+
    ggtitle(title.string)
  return(pl)
}
grid.arrange(
activation.plot(1, map.main, PottsImage, "1"),
activation.plot(2, map.main, PottsImage, "2")
)
#### Simulate 50 times ####
iter.num = 50
output.list <- list()
true_data.list <- list()
for(sim in 1:iter.num){
  beta.given = PottsImage
  beta.given[PottsImage[, 1]==0, 1] = 0
  beta.given[PottsImage[, 1]==1, 1] = rtruncnorm(sum(PottsImage[, 1]==1), a=0, b=Inf, mean = 2, sd = 0.1)
  beta.given[PottsImage[, 1]==2, 1] = rtruncnorm(sum(PottsImage[, 1]==2), a=-Inf, b=0, mean = -1, sd = 0.1)
  
  beta.given[PottsImage[, 2]==0, 2] = 0
  beta.given[PottsImage[, 2]==1, 2] = rtruncnorm(sum(PottsImage[, 2]==1), a=0, b=Inf, mean = 1, sd = 0.1)
  beta.given[PottsImage[, 2]==2, 2] = rtruncnorm(sum(PottsImage[, 2]==2), a=-Inf, b=0, mean = -2, sd = 0.1)
  
  
  
  beta.given = t(beta.given)
  rho =0.5
  lambda = 2
  Q = rho*( diag(c( W%*%matrix(1, num.of.obs, 1) ))-W) + (1-rho)*diag(1, num.of.obs)
  phi = mvrnorm(1, rep(0, num.of.obs), lambda*solve(Q))
  
  psi = arima.sim(model=list(ar=0.2),n=T, sd=sqrt(0.1))
  eta = matrix(0, n, T)
  for(i in 1:n){
    eta[i, ] = x[, , i]%*%beta.given[, i, drop=F] + psi
  }
  
  
  eta = eta+phi	
  
  y = apply(exp(eta), c(1, 2), rpois, n=1)  #nxT
  Eit <- matrix(1,  nrow = T, ncol = n)
  
  
  PottsImageRandom = matrix(sample(0:2, num.of.obs*p, replace = TRUE), num.of.obs, p)
  beta.ini = PottsImageRandom
  beta.ini[PottsImageRandom==0] = 0
  beta.ini[PottsImageRandom==1] = rtruncnorm(sum(PottsImageRandom==1), a=0, b=Inf, mean = 2, sd = 0.1)
  beta.ini[PottsImageRandom==2] = rtruncnorm(sum(PottsImageRandom==2), a=-Inf, b=0, mean = -2, sd = 0.1)
  beta.ini = t(beta.ini)
  
  # InititalValues = list(ActImg = PottsImageRandom, beta = beta.ini, 
  #                       theta=runif(p, 0, 1), phi = phi, psi = psi, rho=rho, lambda=1)
  InititalValues = list(ActImg = PottsImageRandom, beta = beta.ini, 
                        theta=rep(0,p), phi = phi, psi = psi, rho=rho, lambda=1)
  
  output = Update_Gamma_Beta_Theta(5000, 3, W, y, x, ExpU, InititalValues, Eit, UpdateTheta = TRUE, UpdatePsi = TRUE)
  true_data <- list(phi = phi, psi = psi, beta = beta.given, y = y)
  output.list[[sim]] = output
  true_data.list[[sim]] = true_data
}
save(output.list, true_data.list,file = paste(path, "/thesisCode/Simulation/geoSimulation_50Times.RData", sep = ""))

 ## est ###
load(paste(path, "/thesisCode/Simulation/geoSimulation_50Times.RData", sep = ""))
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
EstPottsImage = apply(output.list[[50]]$SimActivationImg, c(1, 2), Mode)
# --- Plots ---
p1 <- activation.plot(1, map.main, PottsImage, "True Cov1")
p2 <- activation.plot(2, map.main, PottsImage, "True Cov2")
p3 <- activation.plot(1, map.main, EstPottsImage, "Est Cov1")
p4 <- activation.plot(2, map.main, EstPottsImage, "Est Cov2")
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

svg(filename = paste(path, "/thesisCode/Simulation/plots/geoSim.svg", sep = ""), width = 10, height = 15)
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
