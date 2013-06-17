# Copyright 2013 Finlay Scott. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, JRC
# drfinlayscott@gmail.com

# Putting together test data set for PML
# See email from M Solan 16/06/2013
# Used simulate_extinction.R to generate the results

library(plyr)
library(reshape)
library(ggplot2)
rm(list=ls())
setwd("~/Work/MEC_analysis")
load("MEC_data/Datv16b.rdata")

# Load up datasets - rbind into a huge dataset
# join with dat to add Lat and Long of each FID
# Only include without compensation for the moment
fid <- 1
load(paste("OA_workshop/Outputs/fast_op/output",fid,".RData",sep=""))
output <- output[output$comp==FALSE,c("species_left","sim","Bpc","OrgC","Chla","FID")]
output <- join(output,unique(dat[dat$FID==fid,c("FID","Long","Lat")]), by="FID")
write.table(output, file = "OA_workshop/Outputs/combined_op/output.csv", row.names=FALSE, sep=",")
for (fid in 2:109){
    cat("FID: ", fid, "\n")
    load(paste("OA_workshop/Outputs/fast_op/output",fid,".RData",sep=""))
    output <- output[output$comp==FALSE,c("species_left","sim","Bpc","OrgC","Chla","FID")]
    output <- join(output,unique(dat[dat$FID==fid,c("FID","Long","Lat")]), by="FID")
    write.table(output, file = "OA_workshop/Outputs/combined_op/output.csv",append=TRUE, row.names=FALSE, col.names=FALSE, sep=",")
}

