rm(list=ls())
source("01_ODE_Function.R", echo=F)
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


load(file="figS18_net.RData")
load("figS19_sp.RData")

str(net_dat)
net_dat$mut_strength<-as.factor(net_dat$mut_strength)

n1<-(net_dat %>% filter(mut_strength == 1.5 | mut_strength == 1.425 | mut_strength==1.350 | mut_strength == 1.275  
                    | mut_strength == 1.2) %>% 
    ggplot(aes(x=Nestedness, y = recovery_richness, color= factor(mut_strength)))+
    geom_point(position=position_jitter(height=0.0,width=0.0),
               alpha = 0.05, size = 3)+
    scale_color_viridis(discrete = TRUE)+
    ylab("Recovery richness")+
    #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
      labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
      # theme(legend.position="none")+
    xlab("Nestedness")+
    # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
    stat_smooth(method = "glm", size=2,alpha=0.1,
                method.args = list(family = "quasibinomial"),
                se =FALSE,
                aes(color=factor(mut_strength))) +
    facet_wrap(~individual_variation))

n2<-(net_dat %>% filter(mut_strength == 1.5 | mut_strength == 1.425 | mut_strength==1.350 | mut_strength == 1.275  
                        | mut_strength == 1.2) %>% 
    ggplot(aes(x=connectance, y = recovery_richness, color= factor(mut_strength)))+
    geom_point(position=position_jitter(height=0.0,width=0.0),
               alpha = 0.05, size = 3)+
    scale_color_viridis(discrete = TRUE)+
    ylab("Recovery richness")+
   labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
    # theme(legend.position="none")+
    xlab("Connectance")+
    # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
    stat_smooth(method = "glm", size=2,alpha=0.1,
                method.args = list(family = "quasibinomial"),
                se =FALSE,
                aes(color=factor(mut_strength))) +
    facet_wrap(~individual_variation))



ggpubr::ggarrange(n1, n2,labels=c("A","B"),
                  nrow=1)




str(sp_dat)

#sp_dat$w<-as.factor(sp_dat$web)
sp_dat$Degree<-as.numeric(as.character(sp_dat$Degree))
sp_dat$Species<-as.numeric(as.character(sp_dat$Species))
sp_dat$response_thr<-as.numeric(as.character(sp_dat$response_thr))
sp_dat$Nestedness<-as.numeric(as.character(sp_dat$Nestedness))
sp_dat$Connectance<-as.numeric(as.character(sp_dat$Connectance))
sp_dat$response_time<-as.numeric(as.character(sp_dat$response_time))
sp_dat$Individual_variation<-as.factor(as.character(sp_dat$Individual_variation))
sp_dat$Network_size<-as.numeric(as.character(sp_dat$Network_size))
sp_dat$mutualism_strength<-as.numeric(as.character(sp_dat$mutualism_strength))
sp_dat$abundance<-as.numeric(as.character(sp_dat$abundance))
sp_dat$Network_respose_time<-as.numeric(as.character(sp_dat$Network_respose_time))
str(sp_dat)


sp_dat$Individual_variation <- plyr::revalue(sp_dat$Individual_variation, c( "high" = "high variation"))



sp_dat$fraction_abundance<- 0

sp_dat$fraction_abundance[which(sp_dat$abundance>= 0.5) ]<-1
sp_dat$fraction_abundance[which(sp_dat$abundance< 0.5) ]<-0


appender <- function(string) 
  latex2exp::TeX(paste("$\\gamma_{\\0} = $", string))  

sp_dat$Individual_variation <- plyr::revalue(spdat_3$Individual_variation, c( "high" = "high variation", "low" = "low variation"))


(a2<-sp_dat %>% group_by(Nestedness,Connectance, forcing_strength,
                          Individual_variation, Network_size , mutualism_strength) %>%
    summarise(count=n(),
              count_abundance =  sum(abundance >= 0.5),
              proportion_abundance = count_abundance/count ) %>%
    filter(forcing_strength == 0.5,mutualism_strength > 1.1) %>% 
    ggplot( aes(y = (proportion_abundance),
                x = Nestedness,
                colour = Individual_variation))+
    geom_point(size=4)+
    ylim(c(0,1))+
    geom_smooth(method = "glm", 
                method.args = list(family = "quasibinomial"), 
                se = T, size=1.5) +
    #geom_smooth(method = lm, formula = y ~x, se = F)+
    theme_classic()+
    theme(legend.position = "right")+
    xlab("Nestedness (NODF)")+
    scale_color_manual(values=c( "#E69F00", "#56B4E9"))+
    labs(color= "individual variation")+
    labs(y = expression(paste("Proportion of species with ", N[i] > 0.5 )))+ 
    facet_wrap(.~mutualism_strength, nrow = 2,ncol = 4, labeller = as_labeller(appender, default = label_parsed)))



(r2<-sp_dat %>% filter(mutualism_strength>1) %>% 
    select(response_time, Individual_variation, 
           response_thr ,Network_size, abundance, mutualism_strength) %>% 
    group_by(Individual_variation, Network_size) %>% 
    summarise(Fraction = mean(abundance, na.rm=T)) %>% 
    ggplot(aes(x = Individual_variation, y = Fraction, fill= Individual_variation ))+
    ylab("Mean species desnity")+
    xlab("Individual variation")+
    xlab("")+
    guides(fill="none")+
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1,
      # move to the right
      justification = -0.15,
      # remove the slub interval
      .width = 0,
      point_colour = NA
    ) +
    theme_classic()+
    geom_boxplot(
      width = 0.25,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5
    )+ scale_fill_manual(values=c( "#E69F00", "#56B4E9")))


ggpubr::ggarrange(a2,r2, labels=c("A","B"), nrow = 2,ncol=1)
