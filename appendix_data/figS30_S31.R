rm(list=ls())

load("growth_rates_net.RData")
require(deSolve) ## for integrating ordinary differential equations
require(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(viridis)
library(ggplot2)
require(reshape2)
library(GGally)
library(network)
library(sna)
library(ggdist)


str(net_dat)
net_dat$mut_strength<-as.factor(net_dat$mut_strength)
net_dat$individual_variation <- plyr::revalue(net_dat$individual_variation, c( "high" = "high variation", "low" = "low variation"))
net_dat$network_size <- as.numeric(as.character(net_dat$network_size))
n1<-(net_dat %>% 
       ggplot(aes(x=Nestedness, y = recovery_richness, 
                  color= factor(growth_rate)))+
       geom_point(position=position_jitter(height=0.0,width=0.0),
                  alpha = 0.05, size = 3)+
       scale_color_viridis(discrete = TRUE)+
       ylab("Recovery richness")+
       #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
       labs(color=expression("growth rate"))+
       # theme(legend.position="none")+
       xlab("Nestedness")+
       # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
       stat_smooth(method = "glm", size=2,alpha=0.1,
                   method.args = list(family = "quasibinomial"),
                   se =FALSE,
                   aes(color=factor(growth_rate))) +
       facet_grid(individual_variation~mut_strength))
n1
n2<-(net_dat %>%
       ggplot(aes(x=connectance, y = recovery_richness, color= factor(growth_rate)))+
       geom_point(position=position_jitter(height=0.0,width=0.0),
                  alpha = 0.05, size = 3)+
       scale_color_viridis(discrete = TRUE)+
       ylab("Recovery richness")+
       labs(color=expression("growth rate"))+
       # theme(legend.position="none")+
       xlab("Connectance")+
       # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
       stat_smooth(method = "glm", size=2,alpha=0.1,
                   method.args = list(family = "quasibinomial"),
                   se =FALSE,
                   aes(color=factor(growth_rate))) +
       facet_grid(individual_variation~mut_strength))




