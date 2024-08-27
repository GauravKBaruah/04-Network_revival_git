
rm(list=ls())
source("01_hyst_functions.R")

require(tidyverse) ## for efficient data manipulation & plotting
require(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(viridis)
library(beepr)
library(GGally)
library(network)
library(sna)
library(nlme)
library(lmerTest)
library(ggplot2)
library(ggnet)
library(gganimate)
library(gifski)
#\############################################# proportion exhibited hysteresis : supplementary figure
load("FigS1.RData")

# reading all the datasets
# calculating nestedness and connectance
mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles


for(r in 1 :nrow(hyst)){
  g<-adj.mat(as.character(newfiles[which(newfiles == hyst$web.name[r])])) #network web names
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant spec
  hyst$network.size[r]<- Aspecies+Plantspecies
}


#at high mutualistic strength, proportion showing histeresis
prop_net<-hyst %>% filter(hysteresis== "yes", mut_strength == 4.35)


for(i in 1:nrow(prop_net)){
  
  
  prop<-prop_net$richness[i]/prop_net$network.size[i] 
  if(prop > 0.7){
    hysteresis_presence<- 0
  }else{
    hysteresis_presence<- 1
  }
  
  prop_net$hysteresis_presence[i]<- hysteresis_presence
  prop_net$proportion_recovery[i]<-prop
  
}

sum(prop_net$hysteresis_presence)/nrow(prop_net)
h0<-ggplot(prop_net, aes(y=proportion_recovery, x= hysteresis, color=network.size))+
  geom_boxplot(fill = "white", outlier.position=NA)+
  geom_jitter(size=3,alpha=0.5,  position=position_jitter(0.2))+
  scale_colour_gradientn(colours = terrain.colors(10))+
  theme_bw()+
  ylim(c(0,1))+
  xlab("Presence of hysteresis")+
  ylab("Recovery richness at highest mutualistic strength")

#

(h1<-hyst %>% filter(hysteresis== "yes") %>% 
    ggplot(aes(x=mut_strength , y = biomass,color=factor(hysteresis)))+
    geom_point(size=5,alpha=0.7)+
    xlim(c(0,4))+
    theme_classic()+
    theme(legend.position = "")+
    scale_color_manual(values= c("#CC79A7"))+
    annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = 1500,
             alpha = .2)+
    labs(color="")+
    geom_segment(aes(x = 1.05, y =75, xend = 2.5, yend = 75),
                 arrow = arrow(length = unit(0.025, "npc")), 
                 size=1.15,  col = "#CC79A7")+
    annotate(geom = "text", x = 1.5, y = 115, label = "No recovery", color="#CC79A7",
             size=3.5,  angle=0)+
    
    #ggtitle("D")+
    xlab( expression(paste("Avg. mutualistic strength,", gamma[0])))+
    ylab("Network biomass"))


(h2<-hyst %>% filter(hysteresis== "yes") %>% 
    ggplot(aes(x=mut_strength , y = richness,color=factor(hysteresis)))+
    geom_point(size=5,alpha=0.7)+
    xlim(c(0,4))+
    theme_classic()+
    #ggtitle("E")+
    theme(legend.position = "")+
    scale_color_manual(values= c("#CC79A7"))+
    annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = 200,
             alpha = .2)+
    labs(color="")+
    geom_segment(aes(x = 1.05, y = 10, xend = 2.25, yend = 10),
                 arrow = arrow(length = unit(0.025, "npc")), 
                 size=1.15,  col = "#CC79A7")+
    annotate(geom = "text", x = 1.5, y = 15, label = "No recovery", color="#CC79A7",
             size=3.5,  angle=0)+
    
    xlab( expression(paste("Avg. mutualistic strength,", gamma[0])))+
    ylab("Richness"))

ggpubr::ggarrange(h1,h2, h0,ncol=3,nrow=1,
                  labels = c("A","B", "C"))
