rm(list=ls())
source("01_ODE_Function.R", echo=F)

load("S_targeted_degree_random_net.RData")

#save(net_dat,file="S_targeted_random_degree.RData")
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


str(net_dat)
net_dat$mut_strength<-as.factor(net_dat$mut_strength)
net_dat$individual_variation <- plyr::revalue(net_dat$individual_variation, 
                                              c( "high" = "high variation", "low" = "low variation"))
net_dat$network_size <- as.numeric(as.character(net_dat$network_size))
n1<-(net_dat %>% 
       ggplot(aes(x=Nestedness, y = recovery_richness, 
                  color= factor(perturbation)))+
       geom_point(position=position_jitter(height=0.0,width=0.0),
                  alpha = 0.05, size = 3)+
       scale_color_viridis(discrete = TRUE)+
       ylab("Recovery richness")+
       #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
       labs(color=expression("Recovery perturbation"))+
       # theme(legend.position="none")+
       xlab("Nestedness")+
       # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
       stat_smooth(method = "glm", size=2,alpha=0.1,
                   method.args = list(family = "quasibinomial"),
                   se =FALSE,
                   aes(color=factor(perturbation))) +
       facet_grid(individual_variation~mut_strength))
n1
n2<-(net_dat %>%
       ggplot(aes(x=connectance, y = recovery_richness, color= factor(perturbation)))+
       geom_point(position=position_jitter(height=0.0,width=0.0),
                  alpha = 0.05, size = 3)+
       scale_color_viridis(discrete = TRUE)+
       ylab("Recovery richness")+
       labs(color=expression("Recovery perturbation"))+
       # theme(legend.position="none")+
       xlab("Connectance")+
       # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
       stat_smooth(method = "glm", size=2,alpha=0.1,
                   method.args = list(family = "quasibinomial"),
                   se =FALSE,
                   aes(color=factor(perturbation))) +
       facet_grid(individual_variation~mut_strength))

n3<-(net_dat %>% filter(network_size < 160) %>% 
       ggplot(aes(x=network_size, y = recovery_richness, color= factor(perturbation)))+
       geom_point(position=position_jitter(height=0.0,width=0.0),
                  alpha = 0.05, size = 3)+
       scale_color_viridis(discrete = TRUE)+
       ylab("Recovery richness")+
       labs(color=expression("Recovery perturbation"))+
       # theme(legend.position="none")+
       xlab("Network size")+
       # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
       stat_smooth(method = "glm", size=2,alpha=0.1,
                   method.args = list(family = "quasibinomial"),
                   se =FALSE,
                   aes(color=factor(perturbation))) +
       facet_grid(individual_variation~mut_strength))

ggpubr::ggarrange(n1, n2,labels=c("A","B"),
                  nrow=1)



