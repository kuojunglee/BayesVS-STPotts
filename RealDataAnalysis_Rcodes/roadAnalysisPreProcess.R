library(sf)
library(sp)
library(terra)
library(stringr)
library(tidyverse)
library(gstat)
library(stringr)
rm(list = ls(all=TRUE))
path = '/Users/yinziwei/Documents/成功大學(碩)/碩士論文總整理/thesisCode/Analysis/交通意外分析/data'
setwd(path)
#### Map Data ####
taiwan_shap <- st_read('鄉(鎮、市、區)界線檔(TWD97經緯度)1130719/TOWN_MOI_1130718.shp')
towns <- taiwan_shap[((!(taiwan_shap$COUNTYNAME %in% c('金門縣', '連江縣', '澎湖縣')))&(!(taiwan_shap$TOWNNAME %in% c('琉球鄉', '綠島鄉', '蘭嶼鄉')))),]
# Transform CRS to WGS84
towns <- st_transform(towns, crs = 4326)
ref.map <- towns
townNames <- paste(substr(towns$COUNTYNAME,1,3),substr(towns$TOWNNAME,1,2),sep = "")

#### Covariate Data ####
## rainfall, temperature
rain.file <- list.files(path = '雨量_2023')
temp.file <- list.files(path = "氣溫_2023")

allRainData <- lapply(paste("雨量_2023/", rain.file,sep = ""), rast)
allTempData <- lapply(paste("氣溫_2023/", temp.file, sep = ""), rast)
plot(allRainData[[1]])
plot(allTempData[[1]])
temperature <- allTempData[[1]]
temperature[temperature == -999] <- NA
plot(temperature)

t = length(rain.file)

towns <- st_transform(towns, crs(temperature)) ## align the map data and the raster data
for(tt in 1:t){
  temperature = allTempData[[tt]]
  rainfall = allRainData[[tt]]
  temperature[temperature == -999] <- NA
  rainfall[rainfall == -999] <- NA
  towns[,paste("temp",tt,sep = "")] = terra::extract(temperature, towns, fun = mean, na.rm = TRUE)[,2]
  towns[,paste("rain",tt,sep = "")] <- terra::extract(rainfall, towns, fun = mean, na.rm = TRUE)[,2]
}


## Air pollution
### Prepare data
pollution.data.list <- list.files(path = "空氣污染_2023",pattern = "月平均值查詢") %>% paste("空氣污染_2023/",.,sep = "")
pollution.data <- lapply(pollution.data.list,read.csv,skip = 2,header = TRUE,fileEncoding = 'big5')
station.name <- lapply(pollution.data.list, readLines, n = 1,encoding = 'big5') %>% unlist() %>% 
  iconv(from = "BIG5", to = "UTF-8") %>% str_sub(start = 8)
pollution_base.df <- read.csv('空氣污染_2023/即時資料(測站位置).csv') %>%
  filter(sitename%in%station.name) %>%
  dplyr::select(sitename, longitude, latitude) %>%
  `rownames<-`(.$sitename) %>%
  .[station.name,]
### Calculate each town's pm2.5 by id2
get_pollution_data <- function(df, pollute_type){
  out <- df %>% filter(測項 == pollute_type) %>% dplyr::select("日期","平均值") %>% arrange(日期)
  return(out$平均值)
}
data.out <- lapply(pollution.data,FUN = get_pollution_data,pollute_type = "PM2.5") %>% 
  unlist() %>%
  array(dim = c(12,78)) %>% t() %>% `colnames<-`(c("Jan","Feb","Mar","Apr","May","Jun",
                                                   "Jul","Aug","Sep","Oct","Dec","Feb"))


pollution.PM25.data <- cbind(pollution_base.df,data.out)
# Convert pollution data to spatial points
pollution_pm25.sf <- st_as_sf(pollution.PM25.data, coords = c("longitude", "latitude"), crs = crs(temperature))

grd <- as.data.frame(temperature, xy=T, na.rm=FALSE)
names(grd) = c("x","y","out")
coordinates(grd) <- c("x","y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- crs(temperature)

out.idw <- gstat::idw(formula = `Jan` ~ 1, locations =  as(pollution_pm25.sf, "Spatial"), newdata = grd, idp = 2.0)
out.idw$var1.pred
plot(out.idw)
pm25.raster <- rast(out.idw)
plot(pm25.raster)


coln <- c("Jan","Feb","Mar","Apr","May","Jun",
          "Jul","Aug","Sep","Oct","Dec","Feb")
for(tt in 1:t){
  formu <- formula(paste(coln[tt],"~ 1"))
  out.idw <- gstat::idw(formula = formu, locations =  as(pollution_pm25.sf, "Spatial"), newdata = grd, idp = 2.0)
  pm25.raster <- rast(out.idw)
  towns[,paste("pm25_",tt,sep = "")] = terra::extract(pm25.raster, towns, fun = mean, na.rm = TRUE)[,2]
}

## 人口
population.filelist <- list.files("人口統計_年齡/人口指標", pattern = "112年") %>% paste("人口統計_年齡/人口指標/",.,sep = "")
population.data <- lapply(population.filelist, read.csv, head=T)
countyName <- substr(ref.map$COUNTYNAME,1,2)
townName <- substr(ref.map$TOWNNAME,1,2)
maping.area_name <- paste(countyName,townName,sep = "")

reorder.row <- function(dat){
  dat <- dat[-1,]
  out <- dat[((!(dat$COUNTY %in% c('金門縣', '連江縣', '澎湖縣')))&(!(dat$TOWN %in% c('琉球鄉', '綠島鄉', '蘭嶼鄉')))),]
  datCounty <- substr(out$COUNTY,1,2)
  datTown <- substr(out$TOWN,1,2)
  datMapping.area_name <- paste(datCounty,datTown,sep = "")
  idx <- sapply(maping.area_name, function(x) {
    which(datMapping.area_name == x)
  })
  out <- out[unlist(idx), ]
  return(out)
}

## P_DEN:人口密度, M_F_RAT: 性別比, A65_A0A14_RAT: 老化指數
for(i in 1:4){
  out1 <- population.data[[i]] %>% reorder.row()
  add_col <- out1 %>% dplyr::select(M_F_RAT, P_DEN, A65_A0A14_RAT)
  data_bind <- cbind(add_col,add_col,add_col)
  colnames(data_bind) <- c(paste("M_F_RAT",i*3-2, sep = ""), paste("P_DEN",i*3-2, sep = ""), paste("A65_A0A14_RAT",i*3-2, sep = ""),
                           paste("M_F_RAT",i*3-1, sep = ""), paste("P_DEN",i*3-1, sep = ""), paste("A65_A0A14_RAT",i*3-1, sep = ""),
                           paste("M_F_RAT",i*3-0, sep = ""), paste("P_DEN",i*3-0, sep = ""), paste("A65_A0A14_RAT",i*3-0, sep = ""))
  towns <- cbind(towns, data_bind)
}

#### 床數資料 ####
bed <- read.csv("/Users/yinziwei/Documents/成功大學(碩)/碩士論文的東西（二）/分析_道路傷亡/病床數資料/STAT (2)/112年12月行政區醫療院所統計_鄉鎮市區.csv")
bed.out <- reorder.row(bed)

#### Outcome Data ####
a2.data.list <- list.files(path = "112年傷亡道路交通事故資料", pattern = "112年度A2交通事故資料")
a2.data.list <- paste("112年傷亡道路交通事故資料/",a2.data.list,sep = "")
a2Accidient.data <- lapply(a2.data.list,read.csv)

for(tt in 1:t){
  df1 <- a2Accidient.data[[tt]]
  t1 <- substr(df1$發生地點,1,5) %>% table() %>% as.data.frame()
  colnames(t1) <- c("townName", "freq")
  acc.freq <- towns %>% dplyr::select(COUNTYNAME,TOWNNAME) %>%
    mutate(townName = townNames) %>% mutate(match_value = left_join(.,t1, by="townName")$freq) %>% 
    mutate(match_value = ifelse(is.na(match_value), 0, match_value))
  ref.map[,paste("freq",tt,sep = "")] = acc.freq$match_value
}


a1Accidient.data <- read.csv('112年傷亡道路交通事故資料/112年度A1交通事故資料.csv')
date.select.v <- c("^202301","^202302","^202303","^202304","^202305","^202306",
                   "^202307","^202308","^202309","^202310","^202311","^202312")
for(tt in 1:t){
  df1 <- a1Accidient.data %>% filter(str_detect(`發生日期`,date.select.v[[tt]]))
  t1 <- substr(df1$發生地點,1,5) %>% table() %>% as.data.frame()
  colnames(t1) <- c("townName", "freq")
  acc.freq <- towns %>% dplyr::select(COUNTYNAME,TOWNNAME) %>%
    mutate(townName = townNames) %>% mutate(match_value = left_join(.,t1, by="townName")$freq) %>% 
    mutate(match_value = ifelse(is.na(match_value), 0, match_value))
  ref.map[,paste("death_freq",tt,sep = "")] = acc.freq$match_value
}

#### Outcome Data2 ####
num_of_population.file <- list.files("人口統計_人口數",pattern = "112年") %>% paste("人口統計_人口數/",.,sep = "")
num_of_population.data <- lapply(num_of_population.file, read.csv)

for(i in 1:4){
  out1 <- num_of_population.data[[i]] %>% reorder.row()
  population_count <- out1$P_CNT %>% as.numeric()
  data_bind <- cbind(population_count, population_count, population_count)
  colnames(data_bind) <- c(paste("population_count",i*3-2, sep = "" ), paste("population_count",i*3-1, sep = "" ), paste("population_count",i*3, sep = "" ))
  ref.map <- cbind(ref.map, data_bind)
}

save(ref.map, towns, bed.out, file = "data.RData")
