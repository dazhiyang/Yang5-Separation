#################################################################################
# This code is written by Dazhi Yang
# Singapore Institute of Manufacturing Technology (SIMTech)
# Agency for Science, Technology and Research (A*STAR)
# emails: yangdazhi.nus@gmail.com
#################################################################################

#Clear all workspace
rm(list = ls(all = TRUE))
#devtools::install_github("GuangchuangYu/chinamap")
libs <- c("dplyr", "ggplot2", "chinamap", "ggthemes", "raster", "maptools", "RColorBrewer", "tiff")
invisible(lapply(libs, library, character.only = TRUE))

#################################################################################
# Inputs
#################################################################################
dir0 <- "/Users/DYang/Dropbox/Working papers/Separation" # paper directory 
plot.size = 7; line.size = 0.15; point.size = 0.7; legend.size = 0.4; text.size = plot.size*5/14
#PROJ <- "+proj=eck4 +ellps=WGS84 +datum=WGS84"
#################################################################################

# read station locations
setwd(file.path(dir0, "Data"))
loc <- read.csv("location.csv")

#################################################################################
# Koeppen-Geiger Climatic Zones
#################################################################################
#gpclibPermit() # turn on the permit
countriesSP <- rworldmap::getMap()
countriesSP <- cleangeo::clgeo_Clean(countriesSP)
countriesSP <- fortify(countriesSP, region = "REGION")

kgc <- raster(file.path(dir0, "Data/Beck_KG_V1_present_0p083.tif")) 
kgc <- as_tibble(rasterToPoints(kgc))
names(kgc) <- c("lon", "lat", "class")
kgc <- kgc %>%
  filter(class > 0)
  
kgc$class <- factor(kgc$class, levels = c(1:30))
levels(kgc$class) <- c("Af", "Am", "As", 
                       "BWh", "BWk", "BSh", "BSk", 
                       "Csa", "Csb", "Csc", 
                       "Cwa", "Cwb", "Cwc", 
                       "Cfa", "Cfb", "Cfc", 
                       "Dsa", "Dsb", "Dsc", "Dsd", 
                       "Dwa", "Dwb", "Dwc", "Dwd", 
                       "Dfa", "Dfb", "Dfc", "Dfd", 
                       "ET", "EF")

kp.color30 <- c(
  "#1100FF", "#0F72FA", "#3397E4", 
  "#FC0200", "#F4998F", "#F0A000", "#EBC35B",
  "#FDFF03", "#D0C708", "#8E8F0A",
  "#8BFF96", "#5BC861", "#349431",
  "#C3FB48", "#62FD33", "#37C800",
  "#FC00FF", "#CB00C4", "#98329F", "#8F5A92",
  "#9FB5FF", "#4377DC", "#4851AE", "#2D028C",
  "#00FAFF", "#43C2FF", "#057B7B", "#004665",
  "#9F9F9F", "#5F615E"
) #"#69B0E5",

p <- ggplot() + 
  geom_raster(data = kgc, aes(lon, lat, fill = class))+
  geom_point(data = loc, aes(x = lon, y = lat), size = point.size*2, shape = 4) +
  #geom_polygon(data = cn, aes(long, lat, group=group), fill = NA, color = "gray50", size = line.size, alpha = 0.7) +
  coord_quickmap(expand = 0) +
  scale_fill_manual(name = "Koppen-Geiger climate classification system", values = kp.color30)+
  theme_minimal() +
  theme(plot.margin = unit(c(0,0,0,-0.1), "lines"), panel.spacing = unit(0, "lines"), panel.grid = element_blank(), text = element_text(family = "Times", size = plot.size), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), strip.text.x = element_text(margin = margin(0,0,0,0, "lines"), size = plot.size), legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(family = "Times", size = plot.size, color = "black"), legend.title = element_text(family = "Times", size = plot.size, color = "black"), legend.key.height = unit(0.2, "lines"), legend.box.margin = unit(c(-0.7,0,0,0), "lines"), legend.background = element_rect(fill = "transparent", colour = "transparent")) +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5))

p

setwd(file.path(dir0, "tex"))
ggsave(filename = "mapStn.pdf", plot = p, path = getwd(), scale = 1, width = 180, height = 110, unit = "mm", dpi = 300)




