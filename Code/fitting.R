#################################################################################
# This code is written by Yizhan Gu 
# Harbin Institute of Technology
# emails: guyizhan33@gmail.com
#################################################################################

####this code is aimed to fit the new models called Yang5, using the divided locations generated in "dividing.R"####

####load packages and clear workspace####
rm(list = ls(all = TRUE))
libs <- c("ggplot2", "dplyr", "insol", "xtable","reshape2","lubridate", "NbClust","cluster","zoo","vegan","factoextra","ggrepel","tidyverse","ggthemes","scatterD3","scatterplot3d")
invisible(lapply(libs, library, character.only = TRUE))
plot.size = 7; line.size = 0.1; point.size = 0.5; legend.size = 0.4; text.size = plot.size*5/14

dir0 <- "/Users/DYang/Dropbox/Working papers/Separation" # paper directory 
#################################################################################
# clustering
#################################################################################
setwd(file.path(dir0, "Data"))
loc <- read.csv("climatology variables of 126 sites.csv")
#loc <- read.csv("normalized climatology global datapoints.csv")
xx <- loc[, 3:5]

####################################################################
# perform clustering based on three variables from 126 stations
####################################################################
# NbClust computes using 26 indexes the "best" number of clusters
nc <- NbClust(xx,min.nc = 2,max.nc = 10,method = "kmeans",distance = "euclidean",index ="all")
# convert indexes to a data.frame
tmp <- as.data.frame(nc$All.index)
# since not all indexes are of the same scale, for better plotting, normalize all indexes
tmp <- apply(tmp, 2, function(x) (x - min(x))/(max(x)-min(x)))
# select five indexes
index.choices <- c("Friedman","Ratkowsky","Rubin","Hartigan","Scott","Cindex","McClain")
data.plot <- tibble(x = rep(c(2:10), length(index.choices)), y = as.numeric(tmp[,index.choices]), index = rep(index.choices, each = 9))
# draw picture of index variation for choosing the best cluster number
p1 <- ggplot(data = data.plot,aes(x=x,y=y,group = index,color=index))+
  geom_point()+
  geom_line()+
  labs(#title="Index Variation of NbCLust",
  x="Cluster Number [dimensionless]",
  y="Index Value [dimensionless]")+
  ylim(c(0,1))+
  theme_bw()+
  scale_color_manual(values = colorblind_pal()(8)[2:11])+
  theme(panel.grid.major=element_line(colour=NA),
        plot.margin = unit(c(0.1,0.1,0,0), "lines"),
        plot.title = element_text(hjust = 0.5),
        #legend.box.background = element_rect(color="black"),
        panel.spacing = unit(0.02, "lines"),
        text = element_text(family = "Times", size = plot.size), 
        strip.text.x = element_text(margin = unit(c(0.05,0,0.05,0), "lines"), size = plot.size), 
        legend.direction = "horizontal",
        legend.text = element_text(family = "Times", size = plot.size), 
        legend.title = element_blank(),
        #legend.title = element_text(family = "Times", size = plot.size), 
        legend.key.height = unit(0.2, "lines"),
        legend.key.width = unit(1.0, "lines"), 
        legend.box.margin = unit(c(-0.5,0,-0.4,0), "lines"), 
        legend.background = element_rect(fill = "transparent", colour = "transparent"), 
        legend.position = "bottom", plot.background = element_rect(fill = "transparent", colour = "transparent"))
p1
setwd(file.path(dir0, "Revision 1"))
ggsave(filename = "cluster index.pdf", plot = p1, scale = 1, width = 85, height = 50, unit = "mm", dpi = 300)

# Based on the plot above, the data can be best separated into 6 clusters, see paper, 5 should be used.
# PCA
#pca_xx = prcomp(xx, center = TRUE, scale = F)
#summary(pca_xx)

# set seed, for reproducibility
set.seed(111)
yy <- kmeans(xx,5,nstart=24)
loc$cluster <- yy$cluster
write.csv(loc, file="climatology variables of 126 sites.csv", row.names = F)
# plot cluster PCA result 
cluster_map <-
             fviz_cluster(yy, data = xx,
             #fviz_cluster(kmeans_xx12, data = xx_transform12,           
             palette = colorblind_pal()(8)[2:7],
             ellipse.type = "euclid",
             star.plot = FALSE, 
             geom = "point",
             repel = TRUE,
             show.clust.cent = TRUE,
             main=NULL,
             pointsize = 1,
             xlab = "Principal component 1 [dimensionless]",
             ylab = "Principal component 2 [dimensionless]")+
             scale_x_continuous(breaks = c(-4,-2,0,2,4))+
               theme_bw()+
               theme(panel.grid.major=element_line(colour=NA),
                              #panel.grid=element_blank(),
                              plot.margin = unit(c(0.1,0.1,0,0), "lines"),
                              plot.title = element_text(hjust = 0.5),
                              #legend.box.background = element_rect(color="black"),
                              panel.spacing = unit(0.02, "lines"),
                              #axis.ticks = element_blank(), 
                              #axis.text = element_blank(),
                              #axis.title = element_blank(),
                              axis.title.x = element_text(size = plot.size, family = "Times"),
                              axis.title.y = element_text(size = plot.size, family = "Times"),
                              axis.text.x = element_text(size = plot.size, family = "Times"),
                              axis.text.y = element_text(size = plot.size, family = "Times"),
                              text = element_text(family = "Times", size = plot.size), 
                              strip.text.x = element_text(margin = unit(c(0.05,0,0.05,0), "lines"), size = plot.size), 
                              legend.box = "horizontal",
                              legend.box.background = element_blank(),
                              legend.text = element_text(family = "Times", size = plot.size), 
                              legend.title = element_blank(),
                              #legend.title = element_text(family = "Times", size = plot.size), 
                              legend.key.height = unit(0.1, "lines"),
                              legend.key.width = unit(1.3, "lines"), 
                              legend.box.margin = unit(c(-0.6,0,-0.4,0), "lines"), 
                              legend.background = element_rect(fill = "transparent", colour = "transparent"), 
                              legend.position = "bottom",
                              #plot.background = element_blank()))
                              plot.background = element_rect(fill = "transparent", colour = "transparent"))
cluster_map
setwd(file.path(dir0, "Revision 1"))
ggsave(filename = "cluster_output.pdf", plot = cluster_map, scale = 1, width = 62.5, height = 62.5, unit = "mm", dpi = 300)

#################################################################################
# Fitting
#################################################################################

path <- "/Volumes/Macintosh Research/Completed works/RSER_2022d_Separation/Data/txt"
setwd(path)
fileNames <- dir(path)

# make empty matrices to hold the parameters of the Yang5 models
pars.nls.Yang5 <- matrix(nrow=nrow(yy$centers),ncol=8,byrow=TRUE)
colnames(pars.nls.Yang5) <- c("C","beta0","beta1","beta2","beta3","beta4","beta5","beta6")

for(k in 1:nrow(yy$centers))
{
  data.fit <- NULL
  stn.fit <- which(loc["cluster"]==k)
  cat("Performing train/test split for cluster", k, "\n")
  pb <- txtProgressBar(min = 0, max = length(stn.fit), style = 3) 
  for (i in 1:length(stn.fit))
  {
    data.tmp <- read.table(fileNames[stn.fit[i]], header = TRUE)
    data.tmp <- as_tibble(data.tmp) %>% # convert to tibble
      dplyr::mutate(Time = lubridate::ymd_hms(Time)) %>% # change time format
      dplyr::mutate(kd = Dh/Gh) %>% # compute diffuse fraction
      dplyr::filter(Z < 80) 

    # fitting/testing data split
    # set.seed(100)
    set.seed(k)
    fit.ind <- sample(1:nrow(data.tmp), nrow(data.tmp)*0.5, replace = FALSE)
    test.ind <- sample(c(1:nrow(data.tmp))[-fit.ind], length(c(1:nrow(data.tmp))[-fit.ind]), replace = FALSE)
    # append data from all stations
    data.fit <- bind_rows(data.fit,add_column(data.tmp[fit.ind,],"stn"=substr(fileNames[stn.fit[i]], 1,3)))
    setTxtProgressBar(pb, i)
    # save testing data for DM test in validating.R
    setwd(file.path(dir0, "Data/testing data"))
    data.test <- data.tmp[test.ind,]
    save(data.test, file =  paste(substr(fileNames[stn.fit[i]], 1,3), ".RData", sep = ""))
    #write.csv(data.tmp[test.ind,], file = paste(substr(fileNames[stn.fit[i]], 1,3), ".csv", sep=''))
    setwd(path)
  }# end i
  close(pb) # Close the connection
  # pars of new model
  #Yang5
  pars_Yang5 <- nls(kd~ (C + (1-C)/(1+exp(beta0+beta1*kt + beta2*AST + beta3*Z + beta4*delta_ktc+ beta6*kds60)) + beta5*kde),data=data.fit,start=list(C=0.036101762,beta0=-0.5744039,beta1=4.3184193,beta2=-0.001119887,beta3=0.0003659426,beta4=-4.795195,beta5=1.441433,beta6=-2.839606))
  pars_Yang5 <- coef(pars_Yang5)
  pars.nls.Yang5[k,] <- as.matrix(pars_Yang5)
}#end k

setwd(file.path(dir0, "Data"))
write.csv(pars.nls.Yang5,file="pars of Yang5.csv",row.names = F)


