################################################################################
# This code is written by Yizhan Gu (and later modified by Dazhi Yang)
# Harbin Institute of Technology
# emails: guyizhan33@gmail.com
################################################################################

####this code is aimed to get the three climatology value of the 126 stations####

####load packages and clear workspace####
rm(list = ls(all = TRUE))
libs <- c("RANN", "dplyr")
invisible(lapply(libs, library, character.only = TRUE))

dir0 <- "/Users/DYang/Dropbox/Working papers/Separation" # paper directory 

# read station data 
setwd(file.path(dir0, "Data"))
loc <- read.csv("location.csv", header=T, na.strings=c("NA"))
stn <- loc[, c("stn", "lon", "lat")]

# AOD
aod <- read.csv("AOD.csv", header = T)
# albedo
albedo <- read.csv("Albedo.csv", header = T)
# cloud
cloud <- read.csv("Cloud.csv", header = T)

points <- tibble(aod, "albedo" = albedo$ALB, "cloud" = cloud$CC)
collocate.idx <- nn2(data = points[,1:2], query = stn[,2:3], k = 1)$nn.idx
result <- points[collocate.idx,]
result <- result %>%
  mutate(stn = loc$stn) %>%
  mutate(CZ = loc$CZ)

# save data
write.csv(result, file="climatology variables of 126 sites.csv", row.names = F)
write.csv(points, file = "normalized climatology global datapoints.csv", row.names = F)

