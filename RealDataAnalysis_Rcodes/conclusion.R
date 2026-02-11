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
library(gtable)
rm(list=ls(all=TRUE))
path = '/Users/yinziwei/Documents/成功大學(碩)/碩士論文總整理'
#### Read Data ####
load(paste(path, "/thesisCode/Analysis/交通意外分析/data/data.RData", sep = ""))
load(paste(path, "/thesisCode/Analysis/交通意外分析/output_5timesMCMC_smallEit_fixInit.RData", sep = ""))
#### Read map data ####

map = st_read('/Users/yinziwei/Documents/成功大學(碩)/碩士論文的東西（二）/STBayesVCM/data/鄉(鎮、市、區)界線檔(TWD97經緯度)1130719/TOWN_MOI_1130718.shp')
map.main = map[((!(map$COUNTYNAME %in% c('金門縣', '連江縣', '澎湖縣')))&(!(map$TOWNNAME %in% c('琉球鄉', '綠島鄉', '蘭嶼鄉')))),]
# map.main <- st_simplify(map.main, dTolerance = 500)
which(!(map.main$TOWNID == towns$TOWNID))


cov.of.interest.fixed <- c("rain","temp","pm25","M_F_RAT","P_DEN","A65_A0A14_RAT")
p = length(cov.of.interest.fixed)
t = 12
num.of.obs = n = dim(map.main)[1]


#### Preprocess Data ####
x = array(NA, c(t,p,n))

for(pp in 1:p){
  x[,pp,] = towns[, grepl(cov.of.interest.fixed[pp], colnames(towns))] %>% st_drop_geometry() %>% t()
}
x = as.numeric(x) %>% array(c(t,p,n))
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
risk = y/(num.of.population/100)

p_star <- sum(y) / sum(num.of.population/100)
# Eit <- (num.of.death+num.in.risk) * p_star
Eit <- num.of.population/100 * p_star
Eit <- Eit %>% t()
Eit.avg <- mean(Eit)
Eit <- Eit/Eit.avg

#### Estimation ####
output <- est.list[[3]]
## theta ##
theta.est <- output$theta_samples %>% apply(c(1), mean)
theta.table <- rbind(cov.of.interest.fixed, round(theta.est, 3))

# credible interval #
theta.UL <- output$theta_samples %>% apply(c(1), quantile, probs = 0.975)
theta.LL <- output$theta_samples %>% apply(c(1), quantile, probs = 0.025)

## alpha ##
alpha.est <- output$alpha_samples %>% apply(c(1,2), mean)
alpha.table <- rbind(cov.of.interest.fixed, round(t(alpha.est), 3))

# credible interval #
alpha.UL <- output$alpha_samples %>% apply(c(1,2), quantile, probs = 0.975)
alpha.LL <- output$alpha_samples %>% apply(c(1,2), quantile, probs = 0.025)

#### Activation Plot ####
cat.gtr.r <- function(x){
  ux <- unique(x)
  ifelse(length(which(tabulate(match(x, ux))/length(x) > 0.5)) != 0, ux[which.max(tabulate(match(x, ux)))], 0)
}
EstPottsImage = apply(output$SimActivationImg, c(1,2), cat.gtr.r)
activation.plot <- function(index.num ,map.main, EstPottsImage){
  map.main.use <- st_simplify(map.main, dTolerance = 1000)
  pl <- map.main.use %>% 
    ggplot(aes(fill = factor(EstPottsImage[, index.num]))) + 
    # geom_sf() +
    geom_sf(size = 0.5, alpha = 0.5) +
    scale_shape_manual(name = "category", values = c(0, 1, 2)) + 
    scale_fill_manual(name = "category", values = c('0' = "grey95", '1' = "#B40404", '2' = "#0B6121")) +
    # lims(x = c(120, 122), y = c(21,26)) +
    scale_x_continuous(breaks = seq(120,122), limits = c(120,122))+
    scale_y_continuous(breaks = seq(21,26), limits = c(21.5,25.5))+
    coord_sf() +
    ggtitle(cov.of.interest.fixed[index.num]) +
    theme_minimal()
  return(pl)
}
activation.plot.list <- lapply(1:length(cov.of.interest.fixed), activation.plot, map.main = map.main, EstPottsImage = EstPottsImage)
# --- Function to extract legend ---
get_legend <- function(myplot) {
  tmp <- ggplotGrob(myplot)
  gtable_filter(tmp, "guide-box")
}
# --- Extract one legend ---
legend <- get_legend(activation.plot.list[[1]])

# Remove legends from all plots
activation.plot.list <- lapply(activation.plot.list, function(p) p + theme(legend.position="none"))

# Build combined layout: plots + legend
combined <- arrangeGrob(
  do.call("arrangeGrob", c(activation.plot.list, ncol = 3)),  # your 12-column layout
  legend,  # second row
  ncol = 2,
  widths = c(10, 1)
)
svg(filename = paste(path,"/thesisCode/Analysis/交通意外分析/plots/activation_plot.svg", sep = ""), height = 10, width = 10)
grid.arrange(combined)
dev.off()
