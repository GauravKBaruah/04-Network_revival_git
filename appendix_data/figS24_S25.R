rm(list=ls())
source("01_functions_ode.R")
library(statmod)
require(deSolve) ## for integrating ordinary differential equations
library(tidyr) ## for efficient data manipulation & plotting
library(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(viridis)
library(beepr)
library(viridis)
library(ggplot2)
library(network)
library(sna)
library(GGally)

theme_set(theme_classic()) 


mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:154]

load("gaussian_width.RData")
net_dat<-sp_dat<-NULL
for(i in 1:4416){
  net_dat <- rbind(net_dat,outt[[i]]$output)
  sp_dat <- rbind(sp_dat, outt[[i]]$ddf)
  
}
require(deSolve) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
require(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(viridis)
library(ggplot2)
require(reshape2)
require(flux)
require(akima)
library(GGally)
library(network)
library(sna)
library(ggdist)

theme_set(theme_classic()) 

str(net_dat)
net_dat$mut_strength<-as.factor(net_dat$mut_strength)
net_dat$individual_variation <- plyr::revalue(net_dat$individual_variation, c( "high" = "high variation", "low" = "low variation"))

n1<-(net_dat %>% 
       ggplot(aes(x=Nestedness, y = recovery_richness, color= factor(width)))+
       geom_point(position=position_jitter(height=0.0,width=0.0),
                  alpha = 0.05, size = 3)+
       scale_color_viridis(discrete = TRUE)+
       ylab("Recovery richness")+
       #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
       labs(color=expression(paste("Gaussian width", ",", omega )))+
       # theme(legend.position="none")+
       xlab("Nestedness")+
       # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
       stat_smooth(method = "glm", size=2,alpha=0.1,
                   method.args = list(family = "quasibinomial"),
                   se =FALSE,
                   aes(color=factor(width))) +
       facet_grid(individual_variation~mut_strength))

n2<-(net_dat %>% 
       ggplot(aes(x=connectance, y = recovery_richness, color= factor(width)))+
       geom_point(position=position_jitter(height=0.0,width=0.0),
                  alpha = 0.05, size = 3)+
       scale_color_viridis(discrete = TRUE)+
       ylab("Recovery richness")+
       #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
       labs(color=expression(paste("Gaussian width", ",", omega )))+
       # theme(legend.position="none")+
       xlab("Connectance")+
       # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
       stat_smooth(method = "glm", size=2,alpha=0.1,
                   method.args = list(family = "quasibinomial"),
                   se =FALSE,
                   aes(color=factor(width))) +
       facet_grid(individual_variation~mut_strength))



ggpubr::ggarrange(n1, n2,labels=c("A","B"),
                  nrow=1)

