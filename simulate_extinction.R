# Copyright 2013 Finlay Scott. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, JRC
# drfinlayscott@gmail.com

# Code to run extinction simulations based on sensitivity (e.g. larger guys go first)
# Includes biomass compensation

library(plyr)
library(reshape)
library(ggplot2)
rm(list=ls())
setwd("~/Work/MEC_analysis")

# Load the data
load("MEC_data/Datv16b.rdata")

# cluster is single value
# bpc can be greater than 1
get_orgc_chla <- function(bpc, cluster){
    if (length(cluster) != 1)
        stop("cluster can only be a single value")
    if (cluster %in% c(1,2,3,5)){
        orgc <- sin(0.108291 - 0.017583 * log10(bpc)) * 100
    }
    if (cluster == 4){
        orgc <- sin(0.121421 - 0.019874 * log10(bpc)) * 100
    }
    if (cluster %in% c(1,2,3,4)){
        chla <- 10^(-1.14372 + 0.50808*log10(bpc))
    }
    if (cluster == 5){
        chla <- 10^(-1.1282 + 0.5046*log10(bpc))
    }
    return(data.frame(OrgC = orgc, Chla = chla))
}


nsim <- 2000
# Pick a FID
#comp <- FALSE 
#fid <- 84
for (fid in 1:109){
    cat("FID: ", fid, "\n")
    fid_data <- dat[dat$FID==fid,]
    nsp <- nrow(fid_data)
    fid_data$BPp <- fid_data$Ai * fid_data$Bpi
    fid_data[,c("Species","Ai","Bi","BPp")]
    total_biomass = sum(fid_data$Bi * fid_data$Ai)
    # Biomass extinction
    extinct_prob <- fid_data$Bi / sum(fid_data$Bi)
     
    output <- data.frame()
    for (comp in c(TRUE,FALSE)){
        extinction_species_array <- array(NA, dim=c(nsp,nsim))
        for (i in 1:nsim){
            extinction_species_array[,i] <- sample(1:nsp, replace=FALSE, prob=extinct_prob)
        }

        lower_tri <- rep(lower.tri(array(NA, dim=c(nsp,nsp)), diag=TRUE),nsim)
        order_array <- array(0, dim=c(nsp,nsp,nsim))
        order_array[lower_tri] <- 1
        order_array <- sweep(order_array, c(1,3), extinction_species_array, "*")

        ai_array <- array(0, dim=c(nsp,nsp,nsim))
        ai_array[lower_tri] <- sweep(order_array, c(1,3), fid_data$Ai, function(x,y) y[x])
        # Quick check of sim 3
        #extinction_species_array[,3]
        #fid_data$Ai[extinction_species_array[,3]]
        #ai_array[,,3]

        bi_array <- array(0, dim=c(nsp,nsp,nsim))
        bi_array[lower_tri] <- sweep(order_array, c(1,3), fid_data$Bi, function(x,y) y[x])

        bpi_array <- array(0, dim=c(nsp,nsp,nsim))
        bpi_array[lower_tri] <- sweep(order_array, c(1,3), fid_data$Bpi, function(x,y) y[x])

        # Correct abundances if compensation is on
        if (comp == TRUE){
            aibi_array <- array(0, dim=c(nsp,nsp,nsim))
            aibi_array[lower_tri] <- sweep(order_array, c(1,3), fid_data$Ai * fid_data$Bi, function(x,y) y[x])
            # So abundances must be corrected so that total biomass is same as species go extinct
            biomass_correction <- total_biomass / apply(aibi_array,c(2,3),sum) 
            # Correct ai
            ai_array <- sweep(ai_array, c(2,3), biomass_correction, "*")
            # Just checking we have corrected for total biomass
            aibi_corrected <- ai_array * bi_array
            if (!(any(apply(aibi_corrected,c(2,3),sum) == total_biomass))){
                stop("We've lost some biomass")
            }
        }
        bpp_array <- ai_array * bpi_array
        bpc <- apply(bpp_array,c(2,3),sum)
        dimnames(bpc) <- list(species_left = nsp:1, sim = 1:nsim)

        bpcm <- melt(bpc)
        output <- rbind(output,cbind(comp=comp,bpcm))
    }
    names(output)[names(output)=="value"] <- "Bpc"
    output <- cbind(output, get_orgc_chla(output$Bpc, cluster=fid_data$Cluster[1]))
    output <- cbind(output, FID = fid)
    save(output,file=paste("OA_workshop/Outputs/fast_op/output",fid,".RData",sep=""))
}


#library(ggplot2)
#ggplot(output) + geom_point(aes(x=species_left, y=Bpc), alpha = 0.1) + facet_wrap(~comp, scales="free")
#
#--------------------------------------------------------------------------
# Compare to original method
load("~/Work/MEC_analysis/OA_workshop/Outputs/fast_op/output84.RData")
newop <- output
load("~/Work/MEC_analysis/OA_workshop/Outputs/output84.RData")

op <- rbind(cbind(newop[,c("comp","species_left","sim","Bpc")], set="new"),
            data.frame(set="old",species_left = output$Nsp, Bpc = output$BPC, sim=output$Simulation, comp=output$Compensation))

ggplot(op) + geom_point(aes(x=species_left, y=Bpc),alpha=0.1) + facet_grid(comp~set, scales="free")




