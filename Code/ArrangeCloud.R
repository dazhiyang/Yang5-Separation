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
dir.input <- "/Volumes/Macintosh Research/Data/High-res climatology/Cloud"
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
grid.final <- grid.final[order(pull(grid.final, "lat"), pull(grid.final, "lon"), decreasing=FALSE),]

################################################################################
# read MODCF cloud cover data
################################################################################
# read the tiff file
setwd(dir.input)
MODCF <- raster("MODCF_meanannual.tif")
# reduce the resolution by a factor of 50, so fewer data points for plotting
MODCF <- aggregate(MODCF, fact = 50)
# get the inherent grid
grid.modcf <- as_tibble(coordinates(MODCF)) %>%
  dplyr::rename(lon = x, lat = y)
# construct a data frame for MODCF cloud climatology
data.modcf <- rasterToPoints(MODCF)
colnames(data.modcf) <- c("lon", "lat", "CC")
data.modcf <- as_tibble(data.modcf)

################################################################################
# normalize CC value to 0-1 range
################################################################################
# use kd-tree to fast-map the raw and final grid, only within "radius" to ensure the found neighbor is in fact a collocated pixel
mapping <- nn2(data = data.modcf[, 1:2], query = grid.final, k = 1, searchtype = "radius", radius = 0.5*sqrt(2))
wrong.ind <- which(mapping$nn.idx == 0)
if(length(wrong.ind)>0)
{
  mapping$nn.idx[wrong.ind] <- NA
}

# make the final data frame
data.modcf.final <- tibble(grid.final[,1], grid.final[,2], data.modcf[mapping$nn.idx, 3])

# get the CC value at the 126 sites
CC <- pull(data.modcf, "CC")
CC.min <- min(CC)
CC.max <- max(CC)

# normalize
data.modcf.final <- data.modcf.final %>%
  mutate(CC = (CC-CC.min)/(CC.max-CC.min))
# check the normalized data
range(data.modcf.final$CC)

################################################################################
# plot to compare the raw and processed data
################################################################################
# coastline 
coastlines_sim2_df <- ne_coastline(scale = "small", returnclass = "sp")

# raw CC map
p1 <- ggplot() +
  geom_raster(data=data.modcf, aes(lon,lat, fill = CC/100)) +
  geom_path(data = coastlines_sim2_df, aes(x = long, y = lat, group = group), size = line.size) +
  viridis::scale_fill_viridis(name = "Cloudy days [%]", direction = 1, option = "B") +
  coord_fixed(xlim = c(-180,180), ylim = c(-90,90), expand = FALSE) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0,0), "lines"), panel.spacing = unit(0.02, "lines"), panel.grid = element_blank(),text = element_text(family = "Times", size = plot.size), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), strip.text.x = element_text(margin = unit(c(0.05,0,0.05,0), "lines"), size = plot.size), legend.direction = "horizontal", legend.text = element_text(family = "Times", size = plot.size), legend.title = element_text(family = "Times", size = plot.size), legend.key.height = unit(0.2, "lines"), legend.key.width = unit(1.8, "lines"), legend.box.margin = unit(c(-1.2,0,-0.6,0), "lines"), legend.background = element_rect(fill = "transparent", colour = "transparent"), legend.position = "bottom", plot.background = element_rect(fill = "transparent", colour = "transparent")) +
  guides(fill = guide_colorbar(title.position = "top"))

p1


# Discrete classes with quantile scale
no_classes <- 10
quantiles <- quantile(data.modcf.final$CC, probs = seq(0, 1, length.out = no_classes + 1))
quantiles <- scales::rescale(quantiles, to = c(0,1))

p2 <- ggplot() +
  geom_raster(data=data.modcf.final, aes(lon,lat, fill = CC)) +
  geom_path(data = coastlines_sim2_df, aes(x = long, y = lat, group = group), size = line.size) +
  viridis::scale_fill_viridis(name = "Cloud climatology [0-1]", option = "D", direction = 1, values = quantiles) +
  coord_fixed(xlim = c(-180,180), ylim = c(-90,90), expand = FALSE) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0,0), "lines"), panel.spacing = unit(0.02, "lines"), panel.grid = element_blank(),text = element_text(family = "Times", size = plot.size), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), strip.text.x = element_text(margin = unit(c(0.05,0,0.05,0), "lines"), size = plot.size), legend.direction = "horizontal", legend.text = element_text(family = "Times", size = plot.size), legend.title = element_text(family = "Times", size = plot.size), legend.key.height = unit(0.2, "lines"), legend.key.width = unit(1.8, "lines"), legend.box.margin = unit(c(-1.2,0,-0.6,0), "lines"), legend.background = element_rect(fill = "transparent", colour = "transparent"), legend.position = "bottom", plot.background = element_rect(fill = "transparent", colour = "transparent")) +
  guides(fill = guide_colorbar(title.position = "top"))

p2

p <- ggpubr::ggarrange(p1, p2, ncol = 1, align = "v", labels = c("(a)", "(b)"), heights = c(1, 1), font.label = list(size = plot.size, color = "black", face = "plain", family = "Times"))

# output plot and csv files
setwd(file.path(dir0, "Revision 1"))
ggsave("Cloud.pdf", p, width = 85, height = 110, unit = "mm")

setwd(file.path(dir0, "Data"))
write.csv(data.modcf.final, file = "Cloud.csv", quote = FALSE, row.names = FALSE)










  
