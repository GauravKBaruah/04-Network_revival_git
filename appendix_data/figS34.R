rm(list=ls())



load("nestedness_matrices_net.RData")
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

str(net_dat)
#net_dat$mut_strength<-as.factor(net_dat$mut_strength)
net_dat$individual_variation <- plyr::revalue(net_dat$individual_variation, c( "high" = "high variation",
                                                                               "low" = "low variation"))


n1<-(net_dat %>% filter(mut_strength > 0.8) %>% 
       ggplot(aes(x=Nestedness, y = recovery_richness, color=individual_variation))+
       geom_point(position=position_jitter(height=0.0,width=0.0),
                  alpha = 0.05, size = 2)+
       ylab("Recovery richness")+
       scale_color_viridis_d()+
       #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
       labs(color=expression("Individual variation"))+
       # theme(legend.position="none")+
       
       xlab("Nestedness")+
       # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
       stat_smooth(method = "glm", size=1.5,alpha=0.1,
                   method.args = list(family = "binomial"),
                   se =FALSE) +
       facet_wrap(.~mut_strength,nrow = 1))

n1

