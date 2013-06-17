# Copyright 2013 Finlay Scott. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, JRC
# drfinlayscott@gmail.com

# Figures for paper

rm(list=ls())
library(plyr)
library(reshape2)
library(ggplot2)
setwd("~/Work/MEC_analysis/OA_workshop")


#-------------------------------------------------------------------------
# Preparing species data 

# Loading original data
load("../MEC_data/Datv16b.rdata")

# Pulling in the species sensitivity scores from spreadsheet
species_dat <- read.csv("species_data_12_1130.csv")
head(species_dat)

# Add species sens score to data set
dat <- join(dat,species_dat[,c("Species","Species.Score")],by="Species")
dat$Total.Score <- dat$Ai * dat$Species.Score

# Add mean abundance, and mean Bpp columns
dat <- ddply(dat, .(Species), transform, mean_Ai = mean(Ai), mean_Bpp = mean(Ai) * Bpi)

# Just pull out species data
species <- unique(dat[,c("Species","Bpi","Bi","Species.Score","mean_Ai", "mean_Bpp","R","M")])
# Compare to plots on page 139 of book chapter - matches
#par(mfrow=c(1,3))
#plot(x=species$mean_Ai, y = species$Bpi, log="xy")
#plot(x=species$mean_Ai, y = species$mean_Bpp, log="xy")
#plot(x=species$Bi, y = species$mean_Bpp, log="xy")

# Remove unscored species
species <- species[species$Species.Score>0,]

#-------------------------------------------------------------------------
# Barplot - just for our own iterest
# Order it
species <- species[order(species$Bpi, decreasing=TRUE),]

# Plot horizontal
par(mar=c(15,2,2,2))
barplot(height=species$Species.Score, names.arg=species$Species, las=3)

# Plot vertical
species <- species[order(species$Species.Score, decreasing=FALSE),]
par(mar=c(2,15,0,2))
barplot(height=species$Species.Score, names.arg=species$Species, las=1, horiz=TRUE, lab.cex = 0.2)

#-------------------------------------------------------------------------
# Correlation figures
# Plot:
# OA sens by Size
#            Bpi (Bioturbation of an individual)
#            Mean abundance across sites
#            Mean Bpp (mean bioturbation of the population across sites)
#            Reworking
#            Movement

plot_function <- function(x,y="Species.Score", xlab, ylab = "", text_x = min(log(species[,x])), text_y = min(log(species[,y])), main, pch = 16){
    plot(x=log(species[,x]), y = log(species[,y]), ylab = ylab, xlab=xlab, main = main, pch = pch)
    fit <- lm(log(species[,y]) ~ log(species[,x]))
    abline(fit)
    r2_text <- bquote(R^{2} ~ "=" ~ .(signif(summary(fit)$r.squared,2)))
    text(text_x, text_y, labels = r2_text, pos = 4)
}

# Will need to move text position about
setEPS() # or use ps.options to fine tune it
postscript("correlation_plots.eps")
par(mfrow=c(2,3))
plot_function(x="Bi",ylab ="Log sensitivity", xlab="Log size (g)", main="a")
plot_function(x="Bpi",xlab="Log BPI", main="b")
plot_function(x="mean_Ai",xlab="Log mean abundance", main="c")
plot_function(x="mean_Bpp",xlab="Log mean BPP", main="d", ylab ="Log sensitivity")
plot_function(x="R",xlab="Log reworking score", main="e")
plot_function(x="M",xlab="Log movement score", main="f")
dev.off()


#-----------------------------------------------------------------------
# Map total sens
library(mapdata)
library(maps)

# Sum over the stations to get 'station sensitivity score'
fid_sens <- ddply(dat, .(FID, Long, Lat), summarise, fid_sens = sum(Total.Score))

range(dat$Long)
range(dat$Lat)
latrange <- c(51, 59) 
lonrange <- c(-4, 9)
plotborder <- 0
latplot <- c(latrange[1]-plotborder,latrange[2]+plotborder)
lonplot <- c(lonrange[1]-plotborder,lonrange[2]+plotborder)

plot(1,1,type="n",xlim=lonplot, ylim=latplot, xlab="Longitude", ylab="Latitude")
map("worldHires", fill=T, col="black",add=T)

fids <- unique(fid_sens$FID)
for (f in fids){
    points(
           x = fid_sens[fid_sens$FID == f,"Long"],
           y = fid_sens[fid_sens$FID == f,"Lat"],
           cex = fid_sens[fid_sens$FID == f,"fid_sens"] * 1e-4,
           pch = 16,
           col="red")
}

#-----------------------------------------------------------
library(mapdata)
library(maps)
# Image test
range(dat$Long)
range(dat$Lat)
latrange <- c(51.5, 58.5) 
lonrange <- c(-3, 8.5)
plotborder <- 0
latplot <- c(latrange[1]-plotborder,latrange[2]+plotborder)
lonplot <- c(lonrange[1]-plotborder,lonrange[2]+plotborder)

plot(1,1,type="n",xlim=lonplot, ylim=latplot, xlab="Longitude", ylab="Latitude")
map("worldHires", fill=T, col="black",add=T)

#image_x <- seq(from= 2, to = 4, length=20)
#image_y <- seq(from= 54, to = 56, length=20)
#image_z <- matrix(rnorm(400), nrow = length(image_y), ncol = length(image_x))
#image(x=image_x, y=image_y, z = image_z, add=TRUE)
# Ok...
comp <- FALSE
measure <- "OrgC"
# width and height in terms of lat and long
box_width <-  0.75
box_height <- 0.75
# Have same dims for each box - so max nsp at any station
max_nsp <- max(ddply(dat, .(FID), summarise, nsp = length(Species))$nsp)
bins <- 70+1
nsp_breaks <- seq(from=1,to=max_nsp)
# What is max of all measures?
min_measure <-2 
max_measure <-20 
breaks <- seq(from=min_measure,to=max_measure,length=bins)
fid <- 1
for (fid in 1:109){
    # centre
    fid_lat_long <- unlist(unique(dat[dat$FID==fid, c("Lat","Long")]))
    nsp <- nrow(dat[dat$FID==fid,])
    image_x <- seq(from= fid_lat_long["Long"],
                   to = fid_lat_long["Long"] + box_width,
                   length=max_nsp)
    image_y <- seq(from= fid_lat_long["Lat"], to = fid_lat_long["Lat"] + box_height, length=bins)

    # Generate data here
    load(paste("~/Work/MEC_analysis/OA_workshop/Outputs/fast_op/output",fid,".RData",sep=""))
    # For the histogram, need to add some rows of 0 when 
    tempop <- output[output$comp == comp,c("species_left","sim",measure)]
    if (nsp < max_nsp){
        extra_bit <- data.frame(species_left = (nsp+1):max_nsp, sim=1,measure=NA)
        names(extra_bit)[names(extra_bit)=="measure"] <- measure
        tempop <- rbind(tempop,extra_bit)
    }

    image_z <- daply(tempop, .(species_left), function(x){
          return(hist(x[,measure],plot=FALSE, breaks=breaks)$counts) })

    image_z[image_z == 0] <- NA

    #image_z <- t(matrix(rnorm(bins*nsp), nrow = length(image_y), ncol = length(image_x)))

    image(x=image_x, y=image_y, z = image_z, add=TRUE)
    rect(xleft = image_x[1],
         ybottom = image_y[1],
         xright = image_x[nsp],
         ytop = image_y[max(which(max(tempop[,measure], na.rm=TRUE) > breaks))],
         border="grey")
    # tick marks?
    # scale marks at least
}

