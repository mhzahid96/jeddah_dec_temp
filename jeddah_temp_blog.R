###### Mustafa Zahid 
###### 01/02/2024
###### An R script to produce the figures in the Jeddah December post 
###### (https://mustafahzahid.netlify.app/post/jeddah_temp/)

## 0) prep libraries and color palettes 
rm(list = ls())
gc()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
require(pacman)
pacman::p_load(sf, feather, raster, rgdal, exactextractr, dplyr, ggplot2,
               tidyverse, tools, terra, testit, rlang, readxl, tabularaster, 
               rworldmap, sp, rworldxtra, data.table, tictoc,  
               doParallel, foreach, ncdf4, fasterize, devtools, 
               remotes, future.apply, doFuture, data.table, doFuture, 
               bigmemory, devtools, MetBrewer)
cols <- MetBrewer::met.brewer("OKeeffe2", type = "continuous")


## 1) read and prep the data 

# jeddah station temperature data 
jeddah_daily <- read_csv(paste0(getwd(), "/data/SA000041024.csv"))
jeddah_daily$year <- substr(jeddah_daily$DATE, 1, 4)
jeddah_daily$month <- substr(jeddah_daily$DATE, 6, 7)
jeddah_daily <- subset(jeddah_daily, month == 12)
jeddah_daily$TAVG[is.na(jeddah_daily$TAVG)] <- (jeddah_daily$TMAX[is.na(jeddah_daily$TAVG)]+jeddah_daily$TMIN[is.na(jeddah_daily$TAVG)])/2
jeddah_daily$TAVG <- jeddah_daily$TAVG/10

jeddah_daily$DATE <- as.Date(jeddah_daily$DATE)
jeddah_daily$year <- as.numeric(jeddah_daily$year)
years <- unique(jeddah_daily$year)
jeddah_daily$day <- as.numeric(substr(jeddah_daily$DATE,9,10))
jeddah_daily$decade[jeddah_daily$year < 1990] <- "1980s" 
jeddah_daily$decade[jeddah_daily$year >= 1990 & jeddah_daily$year < 2000] <- "1990s" 
jeddah_daily$decade[jeddah_daily$year >= 2000 & jeddah_daily$year < 2010] <- "2000s" 
jeddah_daily$decade[jeddah_daily$year >= 2010 & jeddah_daily$year < 2020] <- "2010s" 
jeddah_daily$decade[jeddah_daily$year >= 2021 & jeddah_daily$year < 2023] <- "2021-22"
jeddah_daily$decade[jeddah_daily$year == 2023] <- "2023"

jeddah_daily_decadal <- jeddah_daily %>% 
  dplyr::group_by(decade, day) %>% 
  dplyr::summarise(average = mean(TAVG, na.rm = T))

# Jeddah ERA5 reanalysis data
region_dec <- raster::stack(paste0(getwd(), "/data/data-3.nc"))

# Saudi Shapefile 2nd division
#sa <- read_sf("~/Downloads/data-2/polbnda_sau.shp")
sau_p <- read_sf(paste0(getwd(), "/data/gadm41_SAU_shp/gadm41_SAU_2.shp"))

# get jeddah
jeddah_p <- read_sf(paste0(getwd(), "/data/gadm41_SAU_shp/gadm41_SAU_2.shp"))
jeddah_p <- subset(jeddah_p, NAME_2 == "Jiddah")

# extract and reshape data
region_dec <- readAll(region_dec)
jeddah_long <- exactextractr::exact_extract(region_dec, 
                                            jeddah_p, 
                                            fun = "mean")
jeddah_long <- melt(jeddah_long)
jeddah_long$year <- substr(jeddah_long$variable, 7, 10)
jeddah_long$month <- substr(jeddah_long$variable, 12, 13)
jeddah_long_12 <- subset(jeddah_long, month == 12)
jeddah_long_12$value <- jeddah_long_12$value - 273.15
jeddah_long_12 <- subset(jeddah_long_12, year < 2023)
jeddah_long_12$year <- as.numeric(jeddah_long_12$year)
jeddah_long_12$value <- as.numeric(jeddah_long_12$value)

# Region monthly average 1950-2023
region_yrs <- raster::stack(paste0(getwd(), "/data/data-4.nc"))
region_yrs <- readAll(region_yrs)
sau_p_temp <- exactextractr::exact_extract(region_yrs, 
                                           sau_p, 
                                           fun = "mean", 
                                           append_cols = c("GID_1", "NAME_1", 
                                                           "NAME_2", "geometry"))
sau_p_temp_long <- melt(sau_p_temp,
                        id = c("GID_1","NAME_1","NAME_2"))
sau_p_temp_long$value <- sau_p_temp_long$value - 273.15
sau_p_temp_long$year <- substr(sau_p_temp_long$variable, 7, 10)
sau_p_temp_long$decade[sau_p_temp_long$year <1970] <- "1950-60"
sau_p_temp_long$decade[sau_p_temp_long$year >2009] <- "2010-23"
sau_p_temp_long$month <- substr(sau_p_temp_long$variable, 12, 13)
sau_p_temp_long_12 <- subset(sau_p_temp_long, month == 12)
sau_p_temp_long_12 <- sau_p_temp_long_12 %>% 
  dplyr::group_by(GID_1, NAME_1, NAME_2, decade) %>% 
  dplyr::summarise(temp = mean(value, na.rm = T))
sau_p_temp_long_12_1 <- subset(sau_p_temp_long_12, decade == "1950-60")
colnames(sau_p_temp_long_12_1)[5] <- "temp_1950"
sau_p_temp_long_12_2 <- subset(sau_p_temp_long_12, decade == "2010-23")
colnames(sau_p_temp_long_12_2)[5] <- "temp_2023"
sau_p_temp_long_12_2<- left_join(sau_p_temp_long_12_2, 
                                 sau_p_temp_long_12_1,
                                 by = c("GID_1","NAME_1","NAME_2"))
sau_p_temp_long_12_2$change <- sau_p_temp_long_12_2$temp_2023 - sau_p_temp_long_12_2$temp_1950
sau_p <- left_join(sau_p, 
                   sau_p_temp_long_12_2,
                   by = c("GID_1","NAME_1","NAME_2"))

## 2) plot

############################################################################### fig1: 2023 December vs previous decades 
  plot("",
       "",
       col = alpha("black", 0.25), type = "l", ylim = range(20,30), 
       xlim = range(0:30), frame.plot = F, ylab = "Avg. temperature (°C)", xlab = "Day")
  
  
  # Assuming jeddah_daily_decadal$decade is a factor variable
  for (i in unique(jeddah_daily_decadal$decade)) {
    lines(jeddah_daily_decadal$day[jeddah_daily_decadal$decade == i], 
          jeddah_daily_decadal$average[jeddah_daily_decadal$decade == i], 
          col = alpha("black", 0.25))
  }
  
  n_decades <- length(unique(jeddah_daily_decadal$decade))
  shades <- seq(0.2, 1, length.out = n_decades)
  i <- 1
  
  for (i in 1:length(unique(jeddah_daily_decadal$decade))) {
    lines(jeddah_daily_decadal$day[jeddah_daily_decadal$decade == unique(jeddah_daily_decadal$decade)[i]], 
          jeddah_daily_decadal$average[jeddah_daily_decadal$decade == unique(jeddah_daily_decadal$decade)[i]], 
          col = alpha("black", shades[i]))
    
  }
  mtext("Data Source: King Abdulaziz Weather Station", side = 1, line = 4, cex = 0.8, adj = 1)
  
  text(5, 21, "1980s", col = alpha("black", 0.2000000))
  text(5, 21.5, "1990s", col = alpha("black", 0.3333333))
  text(5, 22, "2000s", col = alpha("black", 0.4666667))
  text(5, 22.5, "2010s", col = alpha("black", 0.6000000))
  text(5, 23, "2021-22", col = alpha("black", 0.7333333))
  text(5, 23.5, "2023", col = "red")
  
  lines(jeddah_daily_decadal$day[jeddah_daily_decadal$decade == "2023"], 
        jeddah_daily_decadal$average[jeddah_daily_decadal$decade == "2023"], 
        col = "red", lwd = 1.25)
  
  mtext("Jeddah December in the past decades vs. today", side = 3, cex = 1.05, adj = 0)


############################################################################### fig2: long term december average 
  plot(jeddah_long_12$year,
       jeddah_long_12$value, 
       type = "l", frame.plot = F, 
       ylab = "Avg December temperature (°C)", 
       xlab = "Year", 
       xlim = range(1950,2028))
  
  mtext("Jeddah long-term December average temp", side = 3, cex = 1.05, adj = 0)
  mtext("Data Source: ERA5-Land reanalysis dataset", side = 1, line = 4, cex = 0.8, adj = 1)
  
  # Add a decade running average
  #decade_avg <- rollmean(jeddah_long_12$value, k = 10, align = "right", fill = NA)
  lines(jeddah_long_12$year, decade_avg, col = "red", lwd = 2)
  #segments(x0 = 2000, x1 = 2000, y0 = 21, y1 = 26, lty = 4)
  #segments(x0 = 2022, x1 = 2022, y0 = 21, y1 = 26, lty = 4)
  segments(x0 = 2000, x1 = 2022, y0 = decade_avg[51], y1 = decade_avg[51], lty = 4)
  segments(x0 = 2022, x1 = 2022, y0 = decade_avg[51], y1 = decade_avg[73], lty = 4)
  text(2024, 24.25, "+0.98°C", adj = 0)
  segments(x0 = 2000, x1 =2003, y0 = 22, y1 = 22, col = "red", lwd = 1.75)
  text(2004.25, 22, "Rolling average", adj = 0, cex = 0.8)
  
############################################################################### fig3: long term december average 
  # Create the plot
  p <- ggplot() +
    geom_sf(data = sau_p, aes(fill = change)) +
    scale_fill_gradientn(name = "Degree Change (°C)", colors = cols) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0, size = 14),
      legend.title = element_text(size = 10)
    ) +
    ggtitle("Change in average December temperature\nbetween 1950-60s and 2010-20s across governorates") +
    labs(
      caption = "Data Source: ERA5 Reanalysis Dataset"
    )
  p
  
  # Print the plot
  print(p)

############################################################################### fig4: seasonal average over time
  jed_riy_long <- subset(sau_p_temp_long, NAME_2 == "Jiddah" | NAME_2 == "Riyadh")
  jed_riy_long_avg <- jed_riy_long %>% 
    dplyr::group_by(NAME_2, month, decade) %>% 
    dplyr::summarise(temp = mean(value, na.rm = T))
  jed_riy_long_avg_recent <- subset(jed_riy_long_avg, decade == "2010-23")
  jed_riy_long_avg_recent$month <- as.numeric(jed_riy_long_avg_recent$month)
  plot(jed_riy_long_avg_recent$month[jed_riy_long_avg_recent$NAME_2 == "Riyadh"], 
       jed_riy_long_avg_recent$temp[jed_riy_long_avg_recent$NAME_2 == "Riyadh"], 
       type = "l", frame.plot = F, xlab = "Month", ylab = "Monthly average temp (°C)", lwd = 1.25, 
       ylim = range(12,37), xaxt = "n")
  axis(1, at = c(1,2,3,4,5,6,7,8,9,10,11,12),
       label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
  lines(jed_riy_long_avg_recent$month[jed_riy_long_avg_recent$NAME_2 == "Jiddah"], 
        jed_riy_long_avg_recent$temp[jed_riy_long_avg_recent$NAME_2 == "Jiddah"], 
        col = "red", lwd = 1.25)
  lines(jed_riy_long_avg_old$month[jed_riy_long_avg_old$NAME_2 == "Riyadh"], 
        jed_riy_long_avg_old$temp[jed_riy_long_avg_old$NAME_2 == "Riyadh"], type = "l", col = "grey", lwd = 1.25)
  lines(jed_riy_long_avg_old$month[jed_riy_long_avg_old$NAME_2 == "Jiddah"], 
        jed_riy_long_avg_old$temp[jed_riy_long_avg_old$NAME_2 == "Jiddah"], col = "pink", lwd = 1.25)
  segments(x0 = 5, x1 = 5.5, y0 = 25, y1 = 25, col = "red")
  segments(x0 = 5, x1 = 5.5, y0 = 24, y1 = 24, col = "pink")
  text(6, 25, "Jeddah: 2010-20s", adj = 0, cex = 0.75)
  text(6, 24, "Jeddah: 1950-60s", adj = 0, cex = 0.75)
  segments(x0 = 5, x1 = 5.5, y0 = 22, y1 = 22, col = "black")
  segments(x0 = 5, x1 = 5.5, y0 = 21, y1 = 21, col = "grey")
  text(6, 22, "Riyadh: 2010-20s", adj = 0, cex = 0.75)
  text(6, 21, "Riyadh: 1950-60s", adj = 0, cex = 0.75)

## 3) check some stats
quantile(sau_p_temp_long_12_2$change, 0.4325)



## end of script