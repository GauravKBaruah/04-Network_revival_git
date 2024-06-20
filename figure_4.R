rm(list=ls())
source("01_ODE_Function.R")
library(foreach)
library(doMC)
library(statmod)
library(ggplot2)
#require(deSolve) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
library(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(viridis)
library(ggdist)
theme_set(theme_classic()) 


# reading all the datasets
# calculating nestedness and connectance
mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:154]


#load("04_figure_data_1.RData")


load("figure_4_network_data.RData")
webdat_2$individual_variation<- factor(webdat_2$individual_variation, 
                                       labels = c("high variation", 
                                                  "low variation"))
webdat_2$network_group <- cut(webdat_2$network_size, 
                            breaks =c(20, 40, 60, 80 ,100, 120, 140, 160 ,180),
                            labels = c("<40", "40-60","60-80" ,"80-100", "100-120", "120-140",
                                       "140-160", ">160"))


webdat_2$Nestedness_group<- cut(webdat_2$Nestedness,
                              breaks = c(0.15, 0.30, 0.45,  0.6, 0.75, 0.9),
                              labels=c("<0.3", "0.3-0.45", "0.45-0.6", "0.6-0.75", ">0.75"))


webdat_2$Connectance_group<- cut(webdat_2$connectance,
                              breaks = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5), 
                              labels=c("<0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4",">0.4"))

tempdat<-webdat_2 %>% filter(mut_strength>1, individual_variation == "high")
# glm_fit_var<-(glm(recovery_richness~ network_size*(mut_strength) , family = "binomial",data=tempdat))
# plot(predictorEffects(glm_fit_var, ~ network_size*(mut_strength)),
#      axes=list(grid=F,
#                x=list(rug=FALSE),
#                y=list(type="response"),
#                main = "", xlab ="Nestedness", ylab="Probability of abrupt collapse"))
# 
# plot(allEffects(glm_fit_var))

(n1<-webdat_2 %>% filter(forcing_strength == 0.5, mut_strength>1) %>% 
  ggplot(aes(x=mut_strength, y = recovery_richness, color= factor(network_group)))+
  geom_point(position=position_jitter(height=0.0,width=0.0),
             alpha = 0.05, size = 3)+
  scale_color_viridis(discrete = TRUE)+
  ylab("Recovery richness")+
  #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
    labs(color="Network size bins")+
   # theme(legend.position="none")+
    xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
 # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
  stat_smooth(method = "glm", size=2,alpha=0.1,
              method.args = list(family = "quasibinomial"),
                se =FALSE,
              aes(color=factor(network_group))) +
  facet_wrap(~individual_variation))


(n2<-webdat_2 %>% filter(forcing_strength == 0.5,mut_strength>1) %>% 
  ggplot(aes(x=mut_strength, y = recovery_richness, 
             color= factor(Nestedness_group)))+
  geom_point(position=position_jitter(height=0.0,width=0.),
             alpha = 0.05, size = 3)+
  scale_color_viridis(discrete = TRUE)+
  ylab("Recovery richness")+
    labs(color="Nestedness bins")+
   # labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
    #xlab("Nestedness (NODF)")+
    xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
    stat_smooth(method = "glm", size=2,alpha=0.1,
              method.args = list(family = "quasibinomial"),
                se =FALSE,
              aes(color=factor(network_group)))+
  facet_wrap(~individual_variation))


 (n3<-webdat_2 %>% filter(forcing_strength == 0.5, mut_strength>1) %>% 
  ggplot(aes(x=mut_strength, y = recovery_richness, 
             color= factor(Connectance_group)))+
  geom_point(alpha = 0.5, size = 3)+
  scale_color_viridis(discrete = TRUE)+
  ylab("Recovery richness")+
     labs(color="Connectance bins")+
     xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
     stat_smooth(method = "glm", size=2,alpha=0.1,
                 method.args = list(family = "quasibinomial"),
                 se =FALSE,
                 aes(color=factor(Connectance_group)))+ 
     facet_wrap(~individual_variation))

#figure 4
ggpubr::ggarrange(n1,n2,n3,
                  nrow=3,ncol=1, labels =c("A", "B", "C"))


###### species level data #################################################

load("figure_4_species_level_data.RData")


str(spdat_3)

#spdat_3$w<-as.factor(spdat_3$web)
spdat_3$Degree<-as.numeric(as.character(spdat_3$Degree))
spdat_3$Species<-as.numeric(as.character(spdat_3$Species))
spdat_3$response_thr<-as.numeric(as.character(spdat_3$response_thr))
spdat_3$Nestedness<-as.numeric(as.character(spdat_3$Nestedness))
spdat_3$Connectance<-as.numeric(as.character(spdat_3$Connectance))
spdat_3$response_time<-as.numeric(as.character(spdat_3$response_time))
spdat_3$Individual_variation<-as.factor(as.character(spdat_3$Individual_variation))
spdat_3$Network_size<-as.numeric(as.character(spdat_3$Network_size))
spdat_3$mutualism_strength<-as.numeric(as.character(spdat_3$mutualism_strength))
spdat_3$abundance<-as.numeric(as.character(spdat_3$abundance))
spdat_3$Network_respose_time<-as.numeric(as.character(spdat_3$Network_respose_time))
str(spdat_3)




spdat_3$fraction_abundance<- 0

spdat_3$fraction_abundance[which(spdat_3$abundance>= 0.5) ]<-1
spdat_3$fraction_abundance[which(spdat_3$abundance< 0.5) ]<-0


(r2<-spdat_3 %>% filter(mutualism_strength>1) %>% 
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



appender <- function(string) 
  latex2exp::TeX(paste("$\\gamma_{\\0} = $", string))  



 (a2<-spdat_3 %>% group_by(Nestedness,Connectance, forcing_strength,
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
     #labs(color= expression(paste(gamma[0])))+
    labs(y = expression(paste("Proportion of species with ", N[i] > 0.5 )))+ 
     facet_wrap(.~mutualism_strength, nrow = 2,ncol = 4, labeller = as_labeller(appender, default = label_parsed)))


#figure 5
ggpubr::ggarrange(a2,
                  ggpubr::ggarrange(r2, ncol=2, labels=c("B")),
                  nrow=2, labels ="A")
