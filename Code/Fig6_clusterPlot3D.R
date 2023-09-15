################################################################################
# This code is written by Yizhan Gu 
# Harbin Institute of Technology
# emails: guyizhan33@gmail.com
################################################################################

#Clear all workspace
rm(list = ls(all = TRUE))
libs <- c("ggplot2", "dplyr", "insol", "lubridate", "scatterplot3d", "scales", "ggthemes")
invisible(lapply(libs, library, character.only = TRUE))

dir0 <- "/Users/DYang/Dropbox/Working papers/Separation" # paper directory 

plot.size = 7; line.size = 0.15; point.size = 0.2; legend.size = 0.4; text.size = 1.5;
#################################################################################
# plot 3D
#################################################################################
setwd(file.path(dir0, "Data"))
data.plot <- read.csv("climatology variables of 126 sites.csv", header = TRUE)

data.plot$cluster <- colorblind_pal()(8)[data.plot$cluster+1]
setwd(file.path(dir0, "Revision 1"))
pdf("cluster3D.pdf", width = 7, height = 7)
par(family = "Times", oma = c(0,0,0,0), cex = 1)
p <- scatterplot3d(data.plot$albedo, data.plot$AOD, data.plot$cloud, pch = 16, type = "p", xlab = "Albedo [0-1]", ylab = "Aerosol [0-1]", zlab = "Cloud [0-1]", scale.y = 0.8, color = data.plot$cluster, cex.axis = 1.4, cex.symbols = 1, cex.lab = 1.4, angle = 60, grid = T, lab.z = 4, lab = c(2, 2, 0.2), font.axis = 1, font.lab = 1, lty.axis = 1, zlim = c(0,1))


legend(p$xyz.convert(0.25, -0.2, -0.15), legend = seq(1:5), bty = "n", col = colorblind_pal()(8)[2:7], pch = 16, horiz = T, xpd = T, cex = 1.4)

#text(p$xyz.convert(data.plot[,c("albedo","AOD","cloud")] + 0.02))

dev.off()




