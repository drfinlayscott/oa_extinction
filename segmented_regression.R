# Copyright 2013 Finlay Scott. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, JRC
# drfinlayscott@gmail.com

# Using segmented regression to identify a breakpoint

library(plyr)
library(reshape)
library(ggplot2)
rm(list=ls())
setwd("~/Work/MEC_analysis")

fid <-84 
load(paste("OA_workshop/Outputs/fast_op/output",fid,".RData",sep=""))

# what does it look like?
ggplot(output) + geom_point(aes(x=species_left, y = Bpc), alpha=0.1) + facet_wrap(~comp, scales="free")

# Old segreg code from MEC analysis
# The segmented regression function
fn <- function(th, x, y) {
	X <- cbind(x, pmax(0, x - th))
	sum(lsfit(X, y)$resid^2)
}

segreg <- function(x,y){
  th.m  <- optimize(fn, range(x), x = x, y = y)$minimum
  fm.m  <- lm(y ~ x + pmax(0, x - th.m))
  m.est <- fm.m$coefficients[2]*th.m + fm.m$coefficients[1]
    output <- list(bp = th.m, y.bp = m.est, y.int = fm.m$coefficients[[1]],
                  y.right = fm.m$fitted.values[1], grad.left = fm.m$coefficients[[2]],
                  #grad.right = (fm.m$fitted.values[1] - m.est) / (max(out$SpeciesLeft) - th.m))
                  grad.right = (fm.m$fitted.values[1] - m.est) / (max(x) - th.m))
  return(output)
}

# Or use segmented package?
# Needs starting estimate of breakpoint parameter
# We only have one predictor so might not be necessary


#segmented(obj, seg.Z, psi, control = seg.control(), model = TRUE, ...)
# Demo with segmented
library(segmented)
set.seed(12)
xx<-1:100
zz<-runif(100)
yy<-2+1.5*pmax(xx-35,0)-1.5*pmax(xx-70,0)+15*pmax(zz-.5,0)+rnorm(100,0,2)
dati<-data.frame(x=xx,y=yy,z=zz)
out.lm<-lm(y~x,data=dati)
plot(dati$x,dati$y)
# two breakpoints
o<-segmented(out.lm,seg.Z=~x,psi=list(x=c(30,60)), control=seg.control(display=FALSE))
slope(o)
plot(o)

# Just one break point
o<-segmented(out.lm,seg.Z=~x,psi=list(x=c(60)), control=seg.control(display=FALSE))
plot(o)

# We only want to fit one
dat <- output[output$comp==FALSE,]
nsp <- max(dat$species_left)
plot(Bpc ~ species_left, data = dat)
# Fit old school seg reg
old_seg_reg <- segreg(x = dat$species_left, y = dat$Bpc)
# Using segmented
lm_fit <- lm(Bpc ~ species_left, data = dat)
abline(lm_fit)
new_seg_reg <- segmented(lm_fit,seg.Z=~species_left,psi=list(species_left=c(nsp/2)), control=seg.control(display=TRUE))
plot(new_seg_reg)
summary(new_seg_reg)
# Compare breakpoints
new_seg_reg$psi[1,"Est."]
old_seg_reg$bp


# Does start point really matter?
start_points <- seq(from = 2, to = nsp-1, length = 20)
new_bp <- rep(NA, length(start_points))
for (i in 1:length(start_points)){
    cat("start point: ", i, "\n")
    new_bp[i] <- segmented(lm_fit,seg.Z=~species_left,psi=list(species_left=start_points[i]), control=seg.control(display=FALSE))$psi[1,"Est."]
}
# It does make a difference - but a small one
# How do we pick which one?
# R2?

fit1 <- segmented(lm_fit,seg.Z=~species_left,psi=list(species_left=start_points[1]), control=seg.control(display=FALSE))
fit2 <- segmented(lm_fit,seg.Z=~species_left,psi=list(species_left=start_points[20]), control=seg.control(display=FALSE))
# Which is 'better'?
summary(fit1)$r.squared
summary(fit2)$r.squared
(fit1)$psi 
(fit2)$psi
# Difference in R2 is tiny, but small change in BP
old_seg_reg
# Original has same BP as highest R2

# Should test if segreg is better fit than the linear
# Then there is a 'collapse' in ecosystem function
# List which stations have a collapse?
# And what is their of declines, before and after collapse
# What about compensation and bifurcations


