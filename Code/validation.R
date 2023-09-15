#################################################################################
# This code is written by Yizhan Gu (later modified by Dazhi Yang) 
# Harbin Institute of Technology
# emails: guyizhan33@gmail.com
#################################################################################

#Clear all workspace
rm(list = ls(all = TRUE))
libs <- c("dplyr", "lubridate", "xtable", "ggplot2", "ggthemes", "forecast", "openxlsx")
invisible(lapply(libs, library, character.only = TRUE))

#################################################################################
# Inputs
#################################################################################
dir0 <- "/Users/DYang/Dropbox/Working papers/Separation" # paper directory 
source(file.path(dir0, "Code/functions.R"))
# models of interest
models <- c("Engerer2", "Starke2", "Starke3", "Yang4", "Yang5")
plot.size = 7; line.size = 0.1; point.size = 0.5; legend.size = 0.4; text.size = plot.size*5/14
#################################################################################

# read metadata file
setwd(file.path(dir0, "Data"))
loc <- read.csv("climatology variables of 126 sites.csv")

setwd(file.path(dir0, "Data/testing data"))
files <- dir()

# some empty matrix to hold the error
Dh.mean = Bn.mean = array(NA, length(files))
MBE.Dh = MAE.Dh = RMSE.Dh = MBE.Bn = MAE.Bn = RMSE.Bn <- matrix(NA, nrow = length(files), ncol = length(models))

# an empty matrix to hold the result of DM tests
n.model <- length(models)
DM.Dh = DM.Bn <- matrix(0, n.model, n.model)

pb = txtProgressBar(min = 0, max = length(files), initial = 0, style = 3)
for(i in 1:length(files))
{
  # load(files[i])
  load(paste(substr(files[i], 1,3), ".RData", sep=""))
  data <- data.test
  clim1 <- loc$cluster[match(substr(files[i], 1, 3), loc$stn)]
  clim2 <- loc$CZ[match(substr(files[i], 1, 3), loc$stn)]
  
  # mean of Gh and Bn
  Dh.mean[i] <- mean(data$Dh)
  Bn.mean[i] <- mean(data$Bn)
  
  # various decomposition models
  kd.list <- list()
  for(m in 1:length(models))
  {
    #cat("Separation modeling using", models[m], "\n")
    if(models[m] %in% c("Yang5"))
    {
      # cluster-dependent models
      clim <- clim1
      assign(paste0("kd.", models[m]), value = get(models[m])(data, clim))
    } else if(models[m] %in% c("Starke3"))
    {
      # climate-dependent models
      clim <- clim2
      assign(paste0("kd.", models[m]), value = get(models[m])(data, clim))
    } else{
      # other models
      assign(paste0("kd.", models[m]), value = get(models[m])(data))
    }
    
    kd.list[[m]] <- get(paste0("kd.", models[m]))
  }
  kd.mat <- matrix(unlist(kd.list), ncol = length(models), byrow = FALSE)
  colnames(kd.mat) <- models
  
  # compute errors
  Dh.pred <- apply(kd.mat, 2, function(x) x*data$Gh)
  Bn.pred <- apply(kd.mat, 2, function(x) (data$Gh-x*data$Gh)/cos(insol::radians(data$Z)))
  
  # errors
  MBE.Dh[i,] <- apply(Dh.pred, 2, function(x) mbe(x, data$Dh))
  MAE.Dh[i,] <- apply(Dh.pred, 2, function(x) mae(x, data$Dh))
  RMSE.Dh[i,] <- apply(Dh.pred, 2, function(x) rmse(x, data$Dh))
  MBE.Bn[i,] <- apply(Bn.pred, 2, function(x) mbe(x, data$Bn))
  MAE.Bn[i,] <- apply(Bn.pred, 2, function(x) mae(x, data$Bn))
  RMSE.Bn[i,] <- apply(Bn.pred, 2, function(x) rmse(x, data$Bn))

  # DM tests
  for(j in 1:(n.model-1))
  {
    for(k in (j+1):n.model)
    {
      # for Dh
      e1.Dh <- data$Dh - Dh.pred[,j]
      e2.Dh <- data$Dh - Dh.pred[,k]
      DM_stat <- dm.test(e1.Dh, e2.Dh, alternative = "two.sided", h = 1, power = 2)$statistic
      #m1 is better or worse than m2
      if(DM_stat < qnorm(0.025))
        DM.Dh[k,j] <- DM.Dh[k,j] + 1
      if(DM_stat > qnorm(0.975))
        DM.Dh[j,k] <- DM.Dh[j,k] + 1
      
      # for Bn
      e1.Bn <- data$Bn - Bn.pred[,j]
      e2.Bn <- data$Bn - Bn.pred[,k]
      DM_stat <- dm.test(e1.Bn, e2.Bn, alternative = "two.sided", h = 1, power = 2)$statistic
      #m1 is better or worse than m2
      if(DM_stat < qnorm(0.025))
        DM.Bn[k,j] <- DM.Bn[k,j] + 1
      if(DM_stat > qnorm(0.975))
        DM.Bn[j,k] <- DM.Bn[j,k] + 1
    }
  }
  
  # for kt-kd plotting
  if(i==64)
  {
    kt.all <- data$kt
    kd.meas <- data$kd
    kd.all <- as.vector(kd.mat)
  }
  
  setTxtProgressBar(pb,i)
}
close(pb)

setwd(file.path(dir0, "Data"))
save(DM.Dh, DM.Bn, Dh.mean, Bn.mean, MBE.Dh, MAE.Dh, RMSE.Dh, MBE.Bn, MAE.Bn, RMSE.Bn, file = "error.RData")
colnames(MBE.Bn) = colnames(MBE.Dh) = colnames(MAE.Bn) = colnames(MAE.Dh) = colnames(RMSE.Bn) = colnames(RMSE.Dh) = models
sheets1 = data.frame("RMSE.Bn"=round(RMSE.Bn, 2), "MBE.Bn"=round(MBE.Bn, 2))
write.csv(sheets1,"Bn_error.csv", )
sheets2 = list("RMSE.Dh"=round(RMSE.Dh, 2), "MBE.Dh" = round(MBE.Dh, 2))
write.csv(sheets2,"Dh_error.csv")

#################################################################################
# Print error tables
#################################################################################
# table 3 
table3.RMSE <- cbind(as.matrix(Bn.mean), RMSE.Bn)
table3.MBE <- cbind(as.matrix(Bn.mean), MBE.Bn)

table3.print <- matrix(NA, nrow = 10, ncol = 6)
for(i in 1:5)
{
 select <- which(loc$cluster == i)
 table3.print[i,] <- round(colMeans(table3.RMSE[select,]), 1)
 table3.print[i+5,] <- round(colMeans(table3.MBE[select,]), 1)
}
table3.print
# compute overall
round(colMeans(table3.RMSE), 1)
round(colMeans(table3.MBE), 1)

# table 4 
table4.RMSE <- cbind(as.matrix(Dh.mean), RMSE.Dh)
table4.MBE <- cbind(as.matrix(Dh.mean), MBE.Dh)

table4.print <- matrix(NA, nrow = 10, ncol = 6)
for(i in 1:5)
{
  select <- which(loc$cluster == i)
  table4.print[i,] <- round(colMeans(table4.RMSE[select,]), 1)
  table4.print[i+5,] <- round(colMeans(table4.MBE[select,]), 1)
}
table4.print
# compute overall
round(colMeans(table4.RMSE), 1)
round(colMeans(table4.MBE), 1)

#################################################################################
# Print ranking tables
#################################################################################
score <- c(1:n.model)
ranking.Bn <- apply(RMSE.Bn, 1, function(x) match(score, order(x))) 
ranking.Dh <- apply(RMSE.Dh, 1, function(x) match(score, order(x)))
mean.rank.Bn <- rowMeans(ranking.Bn)
mean.rank.Dh <- rowMeans(ranking.Dh)

# write.csv(ranking.Bn, file="rank of Bn.csv", row.names = T)
# write.csv(ranking.Dh, file="rank of Dh.csv", row.names = T)

# table 5
comp <- "Dh" # Dh or Bn
error.table <- tibble(model = paste0("\\textsc{", models, "}", sep = ""), as_tibble(get(paste("ranking", comp, sep = "."))[,c(1, 2, 3)]), "$\\cdots$", get(paste("ranking", comp, sep = "."))[,126], get(paste("mean.rank", comp, sep = ".")))

print(xtable(error.table), include.rownames = FALSE, sanitize.text.function=function(x){sanitize.numbers(x, "latex", TRUE, TRUE)})

# table 6
comp <- "Bn" # Dh or Bn
error.table <- tibble(model = paste0("\\textsc{", models, "}", sep = ""), as_tibble(get(paste("ranking", comp, sep = "."))[,c(1, 2, 3)]), "$\\cdots$", get(paste("ranking", comp, sep = "."))[,126], get(paste("mean.rank", comp, sep = ".")))

print(xtable(error.table), include.rownames = FALSE, sanitize.text.function=function(x){sanitize.numbers(x, "latex", TRUE, TRUE)})


#################################################################################
# DM test result
#################################################################################
rownames(DM.Dh) = colnames(DM.Dh) = rownames(DM.Bn) = colnames(DM.Bn) <- models

# plot for Dh
data.plot <- tibble(x = rep(rownames(DM.Dh), each = n.model), y = rep(rownames(DM.Dh), n.model), win = as.vector(DM.Dh))
data.plot$x <- factor(data.plot$x, levels = rownames(DM.Dh))
data.plot$y <- factor(data.plot$y, levels = rownames(DM.Dh))

p1 <- ggplot(data.plot) +
  geom_tile(aes(x = x, y = y, fill = win)) +
  viridis::scale_fill_viridis(name = "# wins",limits = c(0, 126), direction = -1, begin = 0.2) +
  #scale_fill_gradientn(colours = c("lightblue","yellow","red"), limits = c(0,100)) +
  geom_text(aes(x = x, y = y, label = win), size = plot.size*5/14, family = "Times") +
  scale_x_discrete(name = "Model A", expand = c(0,0)) +
  scale_y_discrete(name = "Model B", expand = c(0,0)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0,0), "lines"), legend.position = "bottom", text = element_text(family = "Times"), legend.text = element_text(size = plot.size), legend.title = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key.height = unit(0.3, "lines"), legend.key.width = unit(1.4, "lines"), legend.box.margin = unit(c(-0.7,0,0,0), "lines"), axis.title = element_text(size = plot.size), axis.text.x = element_text(size = plot.size, angle = 90, vjust = 0.5), axis.text.y = element_text(size = plot.size), strip.text.x = element_text(size = plot.size, margin = margin(0.02,0,0.02,0, "lines")), strip.text.y = element_text(size = plot.size, margin = margin(0,0.02,0,0.02, "lines")), panel.spacing = unit(0 , "lines"), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent"))

# plot for Dh
data.plot <- tibble(x = rep(rownames(DM.Bn), each = n.model), y = rep(rownames(DM.Bn), n.model), win = as.vector(DM.Bn))
data.plot$x <- factor(data.plot$x, levels = rownames(DM.Bn))
data.plot$y <- factor(data.plot$y, levels = rownames(DM.Bn))

p2 <- ggplot(data.plot) +
  geom_tile(aes(x = x, y = y, fill = win)) +
  viridis::scale_fill_viridis(name = "# wins",limits = c(0, 126), direction = -1, begin = 0.2) +
  #scale_fill_gradientn(colours = c("lightblue","yellow","red"), limits = c(0,100)) +
  geom_text(aes(x = x, y = y, label = win), size = plot.size*5/14, family = "Times") +
  scale_x_discrete(name = "Model A", expand = c(0,0)) +
  scale_y_discrete(name = "Model B", expand = c(0,0)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0,0), "lines"), legend.position = "bottom", text = element_text(family = "Times"), legend.text = element_text(size = plot.size), legend.title = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key.height = unit(0.3, "lines"), legend.key.width = unit(1.4, "lines"), legend.box.margin = unit(c(-0.7,0,0,0), "lines"), axis.title = element_text(size = plot.size), axis.text.x = element_text(size = plot.size, angle = 90, vjust = 0.5), axis.text.y = element_text(size = plot.size), strip.text.x = element_text(size = plot.size, margin = margin(0.02,0,0.02,0, "lines")), strip.text.y = element_text(size = plot.size, margin = margin(0,0.02,0,0.02, "lines")), panel.spacing = unit(0 , "lines"), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent"))

setwd(file.path(dir0, "Revision 1"))
ggsave(filename = "DMDh.pdf", plot = p1, path = getwd(), scale = 1, width = 60, height = 67, unit = "mm", dpi = 300)
ggsave(filename = "DMBn.pdf", plot = p2, path = getwd(), scale = 1, width = 60, height = 67, unit = "mm", dpi = 300)

#################################################################################
# Box plot
#################################################################################
metrics <- c( "RMSE", "MBE")
comps <- c("Bn", "Dh")
# select only three models for comparison
models.select <- c("Engerer2", "Starke2", "Starke3", "Yang4", "Yang5")
plot.col <- which(models %in% models.select)
# remove stations with high RMSE
# which(RMSE.Bn)

data.plot <- NULL
for(i in 1:length(metrics))
{
  for(j in 1:length(comps))
  {
    error <- get(paste(metrics[i], comps[j], sep = "."))
    for(k in 1:length(models.select))
    {
      data.plot <- bind_rows(data.plot, tibble(error = error[, plot.col][,k], model = models.select[k], comp = comps[j], metric = metrics[i], clim = substr(loc$cluster, 1, 1)))
    }
  }
}
data.plot$metric <- factor(data.plot$metric, levels = metrics)
levels(data.plot$metric) <- c("nRMSE", "nMBE")
data.plot$model <- factor(data.plot$model, levels = models)
data.plot$clim <- factor(data.plot$clim, levels = c(1,2,3,4,5))
levels(data.plot$clim) <- c("1", "2", "3", "4", "5")
data.plot$comp <- factor(data.plot$comp, levels = comps)
levels(data.plot$comp) <- c("italic(B)[italic(n)]", "italic(D)[italic(h)]")

p <- ggplot(data.plot) + 
  geom_hline(yintercept = 0, color="red", linetype="dotted") +
  geom_boxplot(aes(x = clim, y = error, fill = model), size = line.size, outlier.size = point.size, varwidth = TRUE, alpha = 0.8) +
  scale_fill_manual(values = colorblind_pal()(8)[2:7]) +
  #ggrepel::geom_text_repel(aes(x = continent, y = error, label = outlier), na.rm = TRUE, size = 1.8, family = "Times", point.padding = unit(0.15, "lines"), force = 4, segment.alpha = 0.3, segment.size = line.size, color = "blue", position = "dodge") +
  facet_wrap(~ metric + comp, ncol = 2, scales = "free_x", labeller = label_parsed) +
  coord_flip()+
  theme_minimal()+
  theme(plot.margin = unit(c(0.5,0.2,0,0.2), "lines"), panel.spacing = unit(0.5, "lines"), text = element_text(family = "Times", size = plot.size), axis.title = element_blank(), axis.text = element_text(size = plot.size), strip.text.x = element_text(margin = margin(0,0,0,0, "lines"), size = plot.size), plot.title = element_text(family = "Times", size=plot.size, face="bold"), legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = plot.size), strip.text.y = element_text(margin = margin(0,0,0,0, "lines"), size = plot.size))


p

setwd(file.path(dir0, "Revision 1"))
ggsave(filename = "bar.pdf", plot = p, path = getwd(), scale = 1, width = 180, height = 180, unit = "mm", dpi = 300)



#################################################################################
# kt-kd plot
#################################################################################

data.plot.bg <- data.frame(x = kt.all, y = kd.meas)
data.plot <- data.frame(x = rep(kt.all,5), y = kd.all, group = rep(c("Engerer2", "Starke2", "Starke3", "Yang4","Yang5"), each = length(kt.all)))
x = seq(0, 1.1, by = 0.01)
y = sapply(x, Erbs)
data.plot$group <- factor(data.plot$group, levels = c("Engerer2", "Starke2", "Starke3", "Yang4", "Yang5"))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
plot.size = 7
line.size = 0.5
point.size = 0.6
p1 <- ggplot(data.plot) +  
  stat_binhex(data = data.plot.bg, bins = 65, aes(x = x, y = y), fill = "grey60") +  
  stat_binhex(bins = 100, aes(x = x, y = y)) +
  facet_wrap(~group, ncol = 5) +
  viridis::scale_fill_viridis(name = "count",trans = "log", direction = 1, option = "D", breaks = c(1, 10, 100, 1000, 10000)) +
  scale_x_continuous(name = expression(paste("Clearness index, ", italic(k)[italic(t)], " [dimensionless]")), breaks = seq(0, 1, by = 0.2), limits = c(0,1.15), expand = c(0,0)) +
  scale_y_continuous(name = expression(paste("Diffuse fraction, ", italic(k), " [dimensionless]")), limits = c(0,1.1), expand = c(0,0)) +
  #xlab(expression(paste("Clearness index, ", k[t], " [dimensionless]", sep = ""))) +
  #ylab(expression(paste("Diffuse fraction, ", k[d], " [dimensionless]", sep = ""))) +
  #xlim(c(0,1.1))+ylim(c(0,1.1))+
  theme_bw() +
  theme(plot.margin = unit(c(0.1,0.1,0,0.5), "lines"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent", color = NA), legend.position = "right", text = element_text(family = "Times"), legend.text = element_text(size = plot.size), legend.title = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key.width = unit(0.3, "lines"), legend.key.height = unit(1.2, "lines"), legend.box.margin = unit(c(0,0,0,-0.4), "lines"), axis.title = element_text(size = plot.size), axis.text.x = element_text(size = plot.size), axis.text.y = element_text(size = plot.size), strip.text.x = element_text(size = plot.size, margin = margin(0.02,0,0.02,0, "lines")), strip.text.y = element_text(size = plot.size, margin = margin(0,0.02,0,0.02, "lines")), panel.spacing = unit(0 , "lines"))

p1
setwd(file.path(dir0, "Revision 1"))
ggsave(filename = "scatter output.pdf", plot = p1, path = getwd(), scale = 1, width = 180, height = 50, unit = "mm", dpi = 300)











