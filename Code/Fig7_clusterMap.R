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
points <- read.csv("normalized climatology global datapoints.csv", header = T)

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
readr::write_csv(data.plot, "normalized climatology global datapoints.csv")

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


p2 <- ggplot() + 
  #geom_raster(data = data.plot, aes(x = lon, y = lat, fill = tp), show.legend = T)+
  geom_point(data = stn.plot, aes(x = lon, y = lat, color = cluster), shape = 4, show.legend = T, size = point.size) +
  labs(fill = "Cluster of 126 Stations")+
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
        legend.key.width = unit(1.8, "lines"),
        legend.box.background = element_rect(fill = "transparent", color = "transparent"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"), 
        legend.position = "bottom", 
        legend.title = element_text(family = "Times", size = plot.size), 
        legend.box.margin = unit(c(-1.2,0,-0.4,0), "lines"), 
        plot.background = element_rect(fill = "transparent", colour = "transparent"))
  
  
p2
ggsave(filename = "station_map.pdf", plot = p2, path = dir, scale = 1, width = 90, height = 50, unit = "mm", dpi = 300)


#################################################################################
# plot KGC map
#################################################################################
data(climatezones)
kp.color31 <- c(
  "#1100FF", "#0F72FA", "#3397E4", "#69B0E5",
  "#FC0200", "#F4998F", "#F0A000", "#EBC35B",
  "#FDFF03", "#D0C708", "#8E8F0A",
  "#8BFF96", "#5BC861", "#349431",
  "#C3FB48", "#62FD33", "#37C800",
  "#FC00FF", "#CB00C4", "#98329F", "#8F5A92",
  "#9FB5FF", "#4377DC", "#4851AE", "#2D028C",
  "#00FAFF", "#43C2FF", "#057B7B", "#004665",
  "#9F9F9F", "#5F615E"
)

climatezones$Cls = factor(climatezones$Cls, levels = c("Af", "Am", "As", "Aw", "BWh", "BWk", "BSh", "BSk", "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Cfa", "Cfb", "Cfc", "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd", "Dfa", "Dfb", "Dfc", "Dfd", "ET", "EF"))

p <- ggplot() + 
  geom_raster(data = climatezones, aes(Lon, Lat, fill = Cls))+
  geom_point(data = stn.plot, aes(x = lon, y = lat), size = point.size*2, shape = 4) +
  #geom_polygon(data=wmap_df, aes(long,lat, group=group))+
  #geom_polygon(data=world_gSimplify,aes(x=long, y=lat, group = group), size = line.size, color = "grey20", fill = "grey20")+
  coord_quickmap(xlim=range(-180,180), ylim=range(-90,90), expand = 0) +
  #geom_point(aes(y=lat, x=lon, color = error), size = point.size, stroke = 0) +
  #scale_color_distiller(name = "RMSE of Angstrom exponent", limits=c(0, 1), palette='Spectral') +
  scale_fill_manual(name = "Koppen-Geiger climate classification", values = kp.color31)+
  theme_minimal() +
  theme(plot.margin = unit(c(0,0,0,-0.1), "lines"), panel.spacing = unit(0, "lines"), panel.grid = element_blank(), text = element_text(family = "Times", size = plot.size), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), strip.text.x = element_text(margin = margin(0,0,0,0, "lines"), size = plot.size), legend.position = c(0.15, 0.3), legend.direction = "vertical", legend.text = element_text(family = "Times", size = plot.size, color = "black"), legend.title = element_text(family = "Times", size = plot.size, color = "black"), legend.key.height = unit(0.2, "lines"), legend.box.margin = unit(c(-0.7,0,0,0), "lines"), legend.background = element_rect(fill = "transparent", colour = "transparent"), legend.title.align = "top") +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5))

p
ggsave(filename = "KoppenGeiger.pdf", plot = p, path = dir, scale = 1, width = 180, height = 90, unit = "mm", dpi = 300)


#################################################################################
# KGC as row and SCC as col
#################################################################################
#sim <- tibble(tp = data.plot$tp, Cls = climatezones$Cls)
sim <- matrix(nrow = length(c("Af", "Am", "As", "Aw", "BWh", "BWk", "BSh", "BSk", "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Cfa", "Cfb", "Cfc", "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd", "Dfa", "Dfb", "Dfc", "Dfd", "ET", "EF")), ncol = nrow(yy$centers))
rownames(sim) <- c("Af", "Am", "As", "Aw", "BWh", "BWk", "BSh", "BSk", "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Cfa", "Cfb", "Cfc", "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd", "Dfa", "Dfb", "Dfc", "Dfd", "ET", "EF")
colnames(sim) <- seq(1:nrow(yy$centers))
for(i in 1:nrow(sim))
  {
  for (j in 1:ncol(sim)) 
    {
      KGC_ind <- which(as.character(climatezones$Cls)==rownames(sim)[i])
      SCC_ind <- which(as.character(data.plot$tp)==colnames(sim)[j])
      sim[i,j] <- length(intersect(KGC_ind, SCC_ind))
    }
  }# end i
if(sum(sim)==length(points)){print("Correct")}
sim <- as.data.frame(sim)
write.csv(sim, "commom datapoints of KGC and SCC.csv", row.names = T)
  
  
  