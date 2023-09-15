################################################################################
# This code is written by D. Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# email: yangdazhi.nus@gmail.com
################################################################################
#Clear all workspace
rm(list = ls(all = TRUE))
libs <- c("raster", "dplyr", "ggplot2", "ncdf4", "RANN", "kgc", "rnaturalearth", "maptools")
invisible(lapply(libs, library, character.only = TRUE))

#################################################################################
# global variables
################################################################################
dir0 <- "/Users/DYang/Dropbox/Working papers/Separation" # paper directory 
dir.input <- "/Volumes/Macintosh Research/Data/High-res climatology/AOD"
# plot parameters
plot.size = 7; line.size = 0.1; point.size = 0.5; legend.size = 0.4; text.size = plot.size*5/14
#################################################################################


################################################################################
# get KGC grid
# this grid covers some islands, as well as all land areas
# these grid points are to be used as the universal grid, to which all other grids are to be mapped
################################################################################
data("climatezones")
grid.final <- tibble(lon = climatezones$Lon, lat = climatezones$Lat)

# among the 126 ground stations, 6 are far away from any grid cell in "grid.final", and their coordinates are to be manually add to "grid.final"

# read the locations of the 126 stations
setwd(file.path(dir0, "Data"))
loc.126 <- tibble(read.csv("location.csv")) %>%
  dplyr::select(one_of("lon", "lat"))
# find collocated cell within a certain radius
collocate <- nn2(data = grid.final[, 1:2], query = loc.126, k = 1, searchtype = "radius", radius = 0.5*sqrt(2))
# here are the stations without a matching cell
loc.126[which(collocate$nn.idx == 0),]
# round the coordinates to nearest 0.25 or 0.75 
grid.add <- tibble(round((loc.126[which(collocate$nn.idx == 0),1] + 0.25) * 2) / 2 - 0.25,
                   round((loc.126[which(collocate$nn.idx == 0),2] + 0.25) * 2) / 2 - 0.25)
# add to "grid.final"
grid.final <- bind_rows(grid.final, grid.add)
# order the grids according first to latitude second to longitude
grid.final <- grid.final[order(grid.final[,2],grid.final[,1],decreasing=FALSE),]


################################################################################
# read MERRA-2 data
################################################################################
setwd(file.path(dir.input, "MERRA-2"))
files.merra2 <- dir(pattern = "*.nc")

# loop to extract data
MERRA2 <- list()
pb <- utils::txtProgressBar(min = 0, max = length(files.merra2), style = 3)
for(i in 1:length(files.merra2))
{
  ncin <- nc_open(files.merra2[i])
  # get location 
  lon <- ncvar_get(ncin,"lon")
  lat <- ncvar_get(ncin,"lat")
  # MERRA-2 grid
  grid.merra2 <- expand.grid(lon, lat)
  
  #get Total Aerosol Extinction AOT [550 nm], AOT or AOD
  AOD550 <- ncvar_get(ncin,"TOTEXTTAU") 
  fillvalue <- ncatt_get(ncin,"TOTEXTTAU","_FillValue")
  AOD550[AOD550==fillvalue$value] <- NA # replace netCDF fill values with NA's
  
  # close nc file
  nc_close(ncin)
  
  # save each monthly AOD map to a list
  MERRA2[[i]] <- AOD550
  utils::setTxtProgressBar(pb, i)
}  
close(pb)

# aggregate monthly AOD to AOD climatology
AOD.merra2 <- apply(simplify2array(MERRA2), 1:2, mean)

# construct a data frame for MERRA-2 AOD climatology
data.merra2 <- tibble(lon = grid.merra2$Var1, lat = grid.merra2$Var2, AOD = as.vector(AOD.merra2))

# quick plot to check whether the matrix to vector conversion is correct
ggplot(data.merra2) +
  geom_raster(aes(x = lon, y = lat, fill = AOD))


################################################################################
# read Yang-Gueymard data
################################################################################
setwd(file.path(dir.input, "Yang-Gueymard"))
files.YG <- dir(pattern = "*.nc")

# loop to extract data
YG <- list()
pb <- utils::txtProgressBar(min = 0, max = length(files.YG), style = 3)
for(i in 1:length(files.YG))
{
  ncin <- nc_open(files.YG[i])
  # get location 
  lon <- ncvar_get(ncin,"longitude")
  lat <- ncvar_get(ncin,"latitude")
  # Yang-Gueymard grid
  grid.YG <- expand.grid(lon, lat)
  
  #get median AOD of Yang-Gueymard product
  AOD550 <- ncvar_get(ncin,"AOD550 (p = 0.500)") 
 
  # close nc file
  nc_close(ncin)
  
  # save each monthly AOD map to a list
  YG[[i]] <- AOD550
  utils::setTxtProgressBar(pb, i)
}  
close(pb)

# aggregate monthly AOD to AOD climatology
AOD.YG <- apply(simplify2array(YG), 1:2, mean)

# construct a data frame for Yang-Gueymard AOD climatology
data.YG <- tibble(lon = grid.YG$Var1, lat = grid.YG$Var2, AOD = as.vector(AOD.YG))

# quick plot to check whether the matrix to vector conversion is correct
ggplot(data.YG) +
  geom_raster(aes(x = lon, y = lat, fill = AOD))

################################################################################
# combine the two data frames
################################################################################
# if latitude is within \pm 60, use YG, otherwise, use MERRA-2
data.merra2 <- data.merra2 %>%
  filter(lat < -60 | lat > 60)

# combine the two data frames
data.aod.raw <- bind_rows(data.merra2, data.YG)

# use kd-tree to fast-map the raw and final grid, only within "radius" to ensure the found neighbor is in fact a collocated pixel
mapping <- nn2(data = data.aod.raw[, 1:2], query = grid.final, k = 1, searchtype = "radius", radius = 0.5*sqrt(2))
wrong.ind <- which(mapping$nn.idx == 0)
if(length(wrong.ind)>0)
{
  mapping$nn.idx[wrong.ind] <- NA
}

# make the final data frame
data.aod.final <- tibble(grid.final[,1], grid.final[,2], data.aod.raw[mapping$nn.idx, 3]) %>%
  mutate(AOD = round(AOD, 5))

################################################################################
# normalize AOD value to 0-1 range
################################################################################
# get the AOD value at the 126 sites
AOD <- pull(data.aod.final, "AOD")
AOD.min <- min(AOD)
AOD.max <- max(AOD)

# normalize
data.aod.final <- data.aod.final %>%
  mutate(AOD = (AOD-AOD.min)/(AOD.max-AOD.min))
# check the normalized data
range(data.aod.final$AOD)


################################################################################
# plot to compare the raw and processed data
################################################################################
# coastline 
coastlines_sim2_df <- ne_coastline(scale = "small", returnclass = "sp")

# raw AOD map
no_classes <- 12
quantiles <- c(0.0, seq(0.08,0.35, length.out = no_classes), max(data.aod.raw$AOD))
quantiles <- scales::rescale(quantiles, to = c(0,1))

p1 <- ggplot() +
  geom_raster(data=data.merra2, aes(lon,lat, fill = AOD)) +
  geom_raster(data=data.YG, aes(lon,lat, fill = AOD)) +
  geom_path(data = coastlines_sim2_df, aes(x = long, y = lat, group = group), size = line.size) +
  scale_fill_distiller(name = "AOD550 [AOD unit]", palette = "Spectral", values = quantiles, limits = c(0, 0.92)) +
  coord_fixed(xlim = c(-180,180), ylim = c(-90,90), expand = FALSE) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0,0), "lines"), panel.spacing = unit(0.02, "lines"), panel.grid = element_blank(),text = element_text(family = "Times", size = plot.size), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), strip.text.x = element_text(margin = unit(c(0.05,0,0.05,0), "lines"), size = plot.size), legend.direction = "horizontal", legend.text = element_text(family = "Times", size = plot.size), legend.title = element_text(family = "Times", size = plot.size), legend.key.height = unit(0.2, "lines"), legend.key.width = unit(1.8, "lines"), legend.box.margin = unit(c(-1.2,0,-0.6,0), "lines"), legend.background = element_rect(fill = "transparent", colour = "transparent"), legend.position = "bottom", plot.background = element_rect(fill = "transparent", colour = "transparent")) +
  guides(fill = guide_colorbar(title.position = "top"))

p1


# Discrete classes with quantile scale
no_classes <- 10
quantiles <- quantile(data.aod.final$AOD, probs = seq(0, 1, length.out = no_classes + 1))
quantiles <- scales::rescale(quantiles, to = c(0,1))

p2 <- ggplot() +
  geom_raster(data=data.aod.final, aes(lon,lat, fill = AOD)) +
  geom_path(data = coastlines_sim2_df, aes(x = long, y = lat, group = group), size = line.size) +
  viridis::scale_fill_viridis(name = "AOD climatology [0-1]", option = "D", direction = 1, values = quantiles) +
  #scale_fill_distiller(name = "Climatological mean AOD@550nm", palette = "Spectral", values = quantiles, limits = c(0, 1)) +
  coord_fixed(xlim = c(-180,180), ylim = c(-90,90), expand = FALSE) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0,0), "lines"), panel.spacing = unit(0.02, "lines"), panel.grid = element_blank(),text = element_text(family = "Times", size = plot.size), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), strip.text.x = element_text(margin = unit(c(0.05,0,0.05,0), "lines"), size = plot.size), legend.direction = "horizontal", legend.text = element_text(family = "Times", size = plot.size), legend.title = element_text(family = "Times", size = plot.size), legend.key.height = unit(0.2, "lines"), legend.key.width = unit(1.8, "lines"), legend.box.margin = unit(c(-1.2,0,-0.6,0), "lines"), legend.background = element_rect(fill = "transparent", colour = "transparent"), legend.position = "bottom", plot.background = element_rect(fill = "transparent", colour = "transparent")) +
  guides(fill = guide_colorbar(title.position = "top"))

p2

p <- ggpubr::ggarrange(p1, p2, ncol = 1, align = "v", labels = c("(a)", "(b)"), heights = c(1, 1), font.label = list(size = plot.size, color = "black", face = "plain", family = "Times"))

# output plot and csv files
setwd(file.path(dir0, "Revision 1"))
ggsave("AOD.pdf", p, width = 85, height = 110, unit = "mm")

setwd(file.path(dir0, "Data"))
write.csv(data.aod.final, file = "AOD.csv", quote = FALSE, row.names = FALSE)










  
