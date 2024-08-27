rm(list=ls())
source("01_ODE_Function.R")
library(statmod)
library(dplyr)
library(readr)
library(beepr)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(tidyr)
theme_set(theme_classic()) 


# reading all the datasets
# calculating nestedness and connectance
mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:154]


#load("04_figure_data_1.RData")


load("fig_S2.RData")

sp_dat<-net_dat<-NULL
for(i in 1:1472){
  
  sp_dat<-rbind(sp_dat,outt[[i]]$ddf)
  net_dat<-rbind(net_dat,outt[[i]]$output)
  
}


(sfig_1<-net_dat %>% filter(forcing_strength == 0.5, mut_strength<1.4) %>% 
    ggplot(aes(x=Nestedness, y = recovery_richness, color= factor(mut_strength)))+
    geom_point(position=position_jitter(height=0.0,width=0.0),
               alpha = 0.05, size = 3)+
    scale_color_viridis(discrete = TRUE)+
    ylab("Recovery richness")+
    #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
    labs(color="Mutualism strength")+
    ylim(c(0,0.25))+
    # theme(legend.position="none")+
    xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
    theme_classic()+
    # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
    stat_smooth(formula = 'y~x',method = "lm", size=2,alpha=0.1,
                se =FALSE,
                aes(color=factor(mut_strength))) +
    facet_wrap(~individual_variation))


(sfig_2<-net_dat %>% filter(forcing_strength == 0.5, mut_strength<1.4) %>% 
    ggplot(aes(x=connectance, y = recovery_richness, color= factor(mut_strength)))+
    geom_point(position=position_jitter(height=0.0,width=0.0),
               alpha = 0.05, size = 3)+
    scale_color_viridis(discrete = TRUE)+
    ylab("Recovery richness")+
    #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
    labs(color="Mutualism strength")+
    ylim(c(0,0.25))+
    # theme(legend.position="none")+
    xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
    theme_classic()+
    # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
    stat_smooth(formula = 'y~x',method = "lm", size=2,alpha=0.1,
                se =FALSE,
                aes(color=factor(mut_strength)))  +
    facet_wrap(~individual_variation))
