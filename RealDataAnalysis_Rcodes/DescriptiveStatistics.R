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
library(spdep)
library(gtable)
rm(list=ls(all=TRUE))

path = '/Users/yinziwei/Documents/成功大學(碩)/碩士論文總整理'
load(paste(path, "/thesisCode/Analysis/交通意外分析/data/data.RData", sep = ""))

map = st_read(paste(path, '/thesisCode/Analysis/交通意外分析/data/鄉(鎮、市、區)界線檔(TWD97經緯度)1130719/TOWN_MOI_1130718.shp', sep = ""))
map.main = map[((!(map$COUNTYNAME %in% c('金門縣', '連江縣', '澎湖縣')))&(!(map$TOWNNAME %in% c('琉球鄉', '綠島鄉', '蘭嶼鄉')))),]
which(!(map.main$TOWNID == towns$TOWNID))

W = neighbor.matrix.origin <- matrix(as.numeric(st_touches(map.main,sparse = FALSE)), ncol = nrow(map.main))
neighbor.list <- map.main %>% poly2nb(queen = TRUE) %>% nb2listw(style = "W", zero.policy = TRUE)

cov.of.interest.fixed <- c("rain","temp","pm25_","M_F_RAT","P_DEN","A65_A0A14_RAT")
p = length(cov.of.interest.fixed)
t = 12
num.of.obs = n = dim(map.main)[1]


num.in.risk = ref.map[, grepl("^freq", colnames(ref.map))] %>% st_drop_geometry() %>% as.matrix()
num.of.death = ref.map[, grepl("^death_freq", colnames(ref.map))] %>% st_drop_geometry() %>% as.matrix()
num.of.population = ref.map[, grepl("^population_count", colnames(ref.map))] %>% st_drop_geometry() %>% as.matrix()

y = num.of.death + num.in.risk
p_star = sum(y)/sum(num.of.population/100)
Eit = num.of.population/100*p_star
SMR = y/Eit

#### LISA ####
MCi <- localmoran(SMR[, 1], neighbor.list)
map.main$Ii1 <- hotspot(MCi, Prname="Pr(z != E(Ii))", cutoff = 0.05, p.adjust = "none")
map.main$Ii1 <- factor(map.main$Ii1, levels=c("High-High","Low-Low", "Low-High", "High-Low", "Not Significant"))
map.main$Ii1[is.na(map.main$Ii1)] <- "Not Significant"


lmoran <- localmoran(SMR[, 1], neighbor.list)
map.main$var <- SMR[, 1]
map.main$Ii <- lmoran[, "Ii"]
map.main$Z_Ii <- lmoran[, "Z.Ii"]
map.main$p_value <- lmoran[, "Pr(z != E(Ii))"]

mean_val <- mean(SMR[, 1])
map.main <- map.main %>%
  mutate(
    cluster = case_when(
      var >= mean_val & Ii >= 0 & p_value < 0.05 ~ "High-High",
      var <  mean_val & Ii >= 0 & p_value < 0.05 ~ "Low-Low",
      var >= mean_val & Ii <  0 & p_value < 0.05 ~ "High-Low",
      var <  mean_val & Ii <  0 & p_value < 0.05 ~ "Low-High",
      TRUE ~ "Not Significant"
    )
  )

LISA.function <- function(idx, map.main){
  lmoran <- localmoran(SMR[, idx], neighbor.list)
  map.main$var <- SMR[, idx]
  map.main$Ii <- lmoran[, "Ii"]
  map.main$Z_Ii <- lmoran[, "Z.Ii"]
  map.main$p_value <- lmoran[, "Pr(z != E(Ii))"]
  
  mean_val <- mean(SMR[, idx])
  map.main <- map.main %>%
    mutate(
      cluster = case_when(
        var >= mean_val & Ii >= 0 & p_value < 0.05 ~ "High-High",
        var <  mean_val & Ii >= 0 & p_value < 0.05 ~ "Low-Low",
        var >= mean_val & Ii <  0 & p_value < 0.05 ~ "High-Low",
        var <  mean_val & Ii <  0 & p_value < 0.05 ~ "Low-High",
        TRUE ~ "Not Significant"
      )
    )
  shp <- st_simplify(map.main, dTolerance = 1000)
  pl <- ggplot(shp) +
    geom_sf(aes(fill = cluster), color = "white", size = 0.2) +
    scale_fill_manual(
      values = c(
        "High-High" = "#E41A1C",
        "Low-Low" = "#377EB8",
        "High-Low" = "#984EA3",
        "Low-High" = "#4DAF4A",
        "Not Significant" = "grey80"
      ),
      drop = FALSE
    ) +
    labs(title = paste("Month",idx), fill = "Cluster Type") +
    scale_x_continuous(breaks = seq(120,122), limits = c(120,122))+
    scale_y_continuous(breaks = seq(21,26), limits = c(21.5,25.5))+
    coord_sf() +
    theme_minimal()
  return(pl)
}
p.lis <- lapply(1:12, LISA.function, map.main = map.main)
# --- Function to extract legend ---
get_legend <- function(myplot) {
  tmp <- ggplotGrob(myplot)
  gtable_filter(tmp, "guide-box")
}

# --- Extract one legend ---
legend <- get_legend(p.lis[[1]])


# Remove legends from all plots
p.lis <- lapply(p.lis, function(p) p + theme(legend.position="none"))

# Build combined layout: plots + legend
combined <- arrangeGrob(
  do.call("arrangeGrob", c(p.lis, ncol = 4)),  # your 12-column layout
  legend,  # second row
  ncol = 2,
  top = textGrob(
    "Lisa Cluster Map for 12 Month",
    gp = gpar(fontsize = 16, fontface = "bold")  # adjust font size/weight
  ),
  widths = c(10, 1)
)

# --- Save to SVG ---
svg(filename = file.path(path, "thesisCode/Analysis/交通意外分析/plots/LISA_SMR.svg"),
    width = 20, height = 21)  # adjust size as needed
grid.draw(combined)
dev.off()


#### SMR ####
SMR.plot <- function(month.idx, map.main){
  map.main <- map.main %>% mutate(smr = SMR[,month.idx])
  shp <- st_simplify(map.main, dTolerance = 1000)
  m1 <- ggplot(shp, color = "white", size = 0.2) +
    geom_sf(aes(fill = smr)) +
    scale_fill_distiller(palette = "Spectral") +
    # lims(x = c(120, 122), y = c(21,26)) +
    scale_x_continuous(breaks = seq(120,122), limits = c(120,122))+
    scale_y_continuous(breaks = seq(21,26), limits = c(21.5,25.5))+
    labs(title = paste("month", month.idx),fill = "Accdient SMR  ") +
    theme_minimal()
  return(m1)
}
SMR.plot.list <- lapply(1:12, SMR.plot, map.main = map.main)
svg(filename = paste(path,"/thesisCode/Analysis/交通意外分析/plots/SMR_plot.svg", sep = ""), width = 20, height = 16)
do.call("grid.arrange", c(SMR.plot.list, ncol = 4))
dev.off()