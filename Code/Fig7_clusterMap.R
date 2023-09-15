#################################################################################
# This code is written by Yizhan Gu 
# Harbin Institute of Technology
# emails: guyizhan33@gmail.com
#################################################################################

#Clear all workspace
rm(list = ls(all = TRUE))
libs <- c("dplyr", "maptools", "ggplot2", "rnaturalearth", "ggthemes", "ggrepel", "fdm2id", "DMwR2", "ncdf4", "raster", "cluster", "RANN", "factoextra", "openxlsx", "kgc")
invisible(lapply(libs, library, character.only = TRUE))

#################################################################################
# Inputs
#################################################################################
dir0 <- "/Users/DYang/Dropbox/Working papers/Separation" # paper directory 
plot.size = 7; line.size = 0.05; point.size = 0.4; legend.size = 0.4; text.size = 1.5

#################################################################################
# Get stn
#################################################################################
setwd(file.path(dir0, "Data"))
loc <- read.csv("climatology variables of 126 sites.csv",header = T)
stn.plot <- as_tibble(loc[,c("lon", "lat", "cluster")])
stn.plot$cluster <- as.factor(stn.plot$cluster)

# coastline 
coastlines_sim2_df <- ne_coastline(scale = "small", returnclass = "sp")
#################################################################################
# Get corresponding cluster type for each point on map, which doesn't need running
#################################################################################

# create lat/lon grid for plotting
#lon <- seq(-180, 180, by = 0.1)
#lat <- seq(-90, 90, by = 0.1)
#point <- expand.grid(lon, lat)
points <- read.csv("RCC.csv", header = T)

# reproduce the cluster model based on 126 stations
xx <- loc[, 3:5]
set.seed(111)
yy <- kmeans(xx, 5, nstart=24)

# predict cluster type of every point on the map
tp <- predict(yy, newdata = points[,3:5])

# prepare plotting data
data.plot <- as_tibble(points)%>%
  #rename(Lon = Var1, Lat = Var2) %>%
  mutate(tp = tp)
data.plot$tp <- factor(data.plot$tp, levels = c(1:5))

# save climatology map data
readr::write_csv(data.plot, "RCC.csv")

#################################################################################
# plot SCC map
#################################################################################
p1 <- ggplot() + 
  geom_raster(data = data.plot, aes(x = lon, y = lat, fill = tp), show.legend = T)+
  #geom_point(data = stn.plot, aes(x = lon, y = lat, color = cluster), shape = 4, show.legend = T, size = point.size) +
  labs(fill = "Solar climate classification")+
  scale_color_manual(values = colorblind_pal()(8)[2:7])+
  scale_fill_manual(values = colorblind_pal()(8)[2:7]) +
  geom_path(data = coastlines_sim2_df, aes(x = long, y = lat, group = group), size = line.size*2)+
  coord_fixed(xlim = c(-180,180), ylim = c(-90,90), expand = FALSE) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0,0), "lines"), 
        panel.spacing = unit(0.02, "lines"), 
        panel.grid = element_blank(),
        text = element_text(family = "Times", size = plot.size), 
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        legend.direction = "horizontal",
        legend.text = element_text(family = "Times", size = plot.size), 
        #legend.key = element_rect(fill = "transparent", size = plot.size),
        legend.key.height = unit(0.2, "lines"), 
        legend.key.width = unit(1.2, "lines"),
        legend.box.background = element_rect(fill = "transparent", color = "transparent"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"), 
        legend.position = "bottom", 
        legend.title = element_text(family = "Times", size = plot.size), 
        legend.box.margin = unit(c(-1.2,0,-0.4,0), "lines"), 
        plot.background = element_rect(fill = "transparent", colour = "transparent"))


p1
setwd(file.path(dir0, "Revision 1"))
ggsave(filename = "cluster_map.pdf", plot = p1, scale = 1, width = 90, height = 50, unit = "mm", dpi = 300)


  