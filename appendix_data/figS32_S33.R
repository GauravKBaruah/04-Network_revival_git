rm(list=ls())


load("init_trait_dist_net.RData")
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
net_dat$individual_variation <- plyr::revalue(net_dat$individual_variation, c( "high" = "high variation",
                                                                               "low" = "low variation"))

net_dat$init_trait<-as.factor(net_dat$init_trait)
net_dat$init_trait <- plyr::revalue(net_dat$init_trait, c( "c(-0.25,0.25)" = "U(-0.25,0.25)", "c(-0.75,0.75)"="U(-0.75,0.75)"))

n1<-(net_dat %>% 
       ggplot(aes(x=Nestedness, y = recovery_richness, color= factor(init_trait)))+
       geom_point(position=position_jitter(height=0.0,width=0.0),
                  alpha = 0.05, size = 2)+
       ylab("Recovery richness")+
       scale_color_viridis_d()+
       #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
       labs(color=expression("Initial trait distribution"))+
       # theme(legend.position="none")+
       
       xlab("Nestedness")+
       # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
       stat_smooth(method = "glm", size=2,alpha=0.1,
                   method.args = list(family = "quasibinomial"),
                   se =FALSE,
                   aes(color=factor(init_trait))) +
       facet_grid(individual_variation~mut_strength))

n1
n2<-(net_dat %>% 
       ggplot(aes(x=connectance, y = recovery_richness, color= factor(init_trait)))+
       geom_point(position=position_jitter(height=0.0,width=0.0),
                  alpha = 0.05, size = 2)+
       scale_color_viridis(discrete = TRUE)+
       ylab("Recovery richness")+
       #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
       labs(color=expression("Initial trait distribution"))+
       # theme(legend.position="none")+
       xlab("Connectance")+
       # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
       stat_smooth(method = "glm", size=2,alpha=0.1,
                   method.args = list(family = "quasibinomial"),
                   se =FALSE,
                   aes(color=factor(init_trait))) +
       facet_grid(individual_variation~mut_strength))



