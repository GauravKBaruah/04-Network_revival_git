

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



#theme_set(theme_classic()) 


# reading all the datasets
# calculating nestedness and connectance
mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles

trait_data_list<-list()
trait_data<-expand.grid(web=c(newfiles)) %>% 
  as_tibble %>% mutate(networksize=0) 



for(i in 1:nrow(trait_data)){
  
  g<-adj.mat(myfiles[which(myfiles == trait_data$web[i])]) #network web names
  # g<-g[-1,-1] 
  
  
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species
  
  
  nestedness<-nestedness_NODF(g)
  C<-Connectance(g)
  
  trait_data$networksize[i] <- Aspecies+Plantspecies
  muinit <-runif((Aspecies+Plantspecies), -0.5,0.5)
  
  trait_data_list[i]<-list(muinit)
  
  names(trait_data_list)[i] <- as.character(trait_data$web[i])
  
}


fact2<- trait_data %>% filter(networksize < 185)
webfiles<-fact2$web



load("hysteresis.RData")

hyst<-NULL
for(i in 1:2800){
  
  hyst <- rbind(hyst,outt[[i]])
  
  
}


c1<-hyst %>% filter(hysteresis== "no") %>% 
ggplot(aes(x=mut_strength , y = biomass,color=factor(hysteresis)))+
  geom_point(size=5,alpha=0.7)+
  xlim(c(0,4))+
  theme_classic()+
  theme(legend.position = "")+
  scale_color_manual(values= c("#CC79A7"))+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = 1500,
          alpha = .2)+
  labs(color="")+
  #ggtitle("D")+
  xlab( expression(paste("Avg. mutualistic strength,", gamma[0])))+
  ylab("Network biomass")


c2<-hyst %>% filter(hysteresis== "no") %>% 
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
  xlab( expression(paste("Avg. mutualistic strength,", gamma[0])))+
  ylab("Richness")


hyst %>% filter(hysteresis== "yes") %>% 
  ggplot(aes(x=mut_strength , y = richness,color=factor(hysteresis)))+
  geom_point(size=5,alpha=0.7)+
  xlim(c(0,4.5))+
  theme_classic()+
  theme(legend.position = "")+
  #ggtitle("B")+
  scale_color_manual(values= c("#CC79A7"))+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = 180,
           alpha = .2)+
  labs(color="")+
  xlab("Average mutualistic strength")+
  ylab("Richness")


(w1<-hyst %>% filter(web.name == "datasets_1/M_PL_018.csv" )%>% 
  ggplot(aes(x=mut_strength , y = richness,color=factor(hysteresis)))+
  geom_point(size=5,alpha=0.7)+
  xlim(c(0,4))+
  theme_classic()+
 # ggtitle("B")+
  theme(legend.position = "")+
  scale_color_manual(values= c("#009E73", "#CC79A7"))+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = 150,
           alpha = .2)+
  labs(color="")+
  geom_curve(aes(x = 2.5, y = 145, xend = 1.45, yend = 75),
             arrow = arrow(length = unit(0.025, "npc")), 
             size=1.15, alpha = 0.25, curvature = 0.3, col="#CC79A7" )+
  annotate(geom = "text", x = 2.5, y = 150, label = "Collapse regime", color="#CC79A7", 
           size=3.5,  hjust = "left")+
  geom_curve(aes(x = 1.5, y = 5, xend = 2.5, yend = 40),
             arrow = arrow(length = unit(0.025, "npc")), 
             size=1.15, curvature = 0.3, col = "#009E73")+
    annotate(geom = "text", x = 1.25, y = 25 ,label ="Recovery regime",
             color="#009E73", size=3.5,  hjust = "left")+
    
    #scale_color_manual(values= c("#009E73", "#CC79A7"))+
  xlab( expression(paste("Avg. mutualistic strength,", gamma[0])))+
  ylab("Richness"))

(w2<-hyst %>% filter(web.name == "datasets_1/M_PL_018.csv" )%>% 
  ggplot(aes(x=mut_strength , y = biomass,color=factor(hysteresis)))+
  geom_point(size=5,alpha=0.5)+
  xlim(c(0,4))+
  theme_classic()+
 # ggtitle("C")+
  theme(legend.position = "")+
  scale_color_manual(values= c("#009E73", "#CC79A7"))+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = 1150,
           alpha = .2)+
  labs(color="")+
  geom_segment(aes(x = 2.5, y = 550, xend = 1.45, yend = 55),
             arrow = arrow(length = unit(0.025, "npc")), 
             size=1.15, alpha = 0.25, col="#CC79A7" )+
  annotate(geom = "text", x = 1.7, y = 300, label = "Collapse regime", color="#CC79A7",
           size=3.5,  angle=55)+
  geom_segment(aes(x = 1.45, y = 40, xend = 2.5, yend = 40),
             arrow = arrow(length = unit(0.025, "npc")), 
             size=1.15,  col = "#009E73")+
    annotate(geom = "text", x = 2.45, y = 95, label = "Recovery regime", color="#009E73",
             size=3.5,  angle=0)+
    
  xlab( expression(paste("Avg. mutualistic strength,", gamma[0])))+
  ylab("Biomass"))

g<-  adj.mat(as.character(webfiles[which(webfiles == "datasets_1/M_PL_018.csv" )])) 

Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

#degree of plants and anichmals
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}

net1 = network(g, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net1)
net1 %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net1 %v% "color" = ifelse(net1 %v% "groups" == "plants", "#0072B2", "#E69F00" )
#ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)

deg<-c(degree.plants,degree.animals)
w6<-ggnet2(net1, mode="circle", size=deg, edge.size = 1,max_size =8, color ="color",edge.alpha = 1, legend.position = "",
           label ="A")
w6




#########################################################################################################################

load("Hysteresis_species_data.RData")


appender <- function(string) 
  latex2exp::TeX(paste("$\\gamma_{\\0} = $", string))  


(a1<-sp_dat %>% filter(Individual_variation=="high variation", mutualism_strength == 1.8 |
                         mutualism_strength == 2,
                       Nestedness < 0.8) %>%
    group_by(Nestedness,Connectance, forcing_strength,Individual_variation, Network_size , mutualism_strength) %>%
    summarise(count=n(),
              count_abundance = mean(abundance),
              proportion_abundance = count_abundance/count ) %>%
    ggplot( aes(y = (count_abundance),
                x = Nestedness,
                colour = factor(mutualism_strength)))+
    geom_point(size=5,alpha=0.6)+
    geom_smooth(method = lm, formula = y ~x, se = F)+
    theme_classic()+
    xlab("Nestedness")+scale_color_manual(values=c( "#0072B2", "#D55E00"))+
    labs(color= expression(paste( gamma[0])))+
    ylab("Mean pollinator density"))
# labs(y = expression(paste("Proportion of species with ", N[i] > 0.5 )))+)# nrow = 2,ncol = 4, labeller = as_labeller(appender, default = label_parsed)))

(a2<-sp_dat %>% filter( mutualism_strength > 1.4, Individual_variation == "high variation",mutualism_strength == 1.7 | mutualism_strength == 1.9, Nestedness < 0.8) %>%
    group_by(Nestedness,Connectance, forcing_strength,
             Individual_variation, Network_size , mutualism_strength) %>%
    summarise(count=n(),
              count_abundance = mean(abundance),
              proportion_abundance = count_abundance/count ) %>%
    filter(Network_size < 200) %>% 
    ggplot( aes(y = (proportion_abundance),
                x = Nestedness,
                colour = factor(mutualism_strength)))+
    geom_point(size=4)+
    ggtitle("B")+
    geom_smooth(method = lm, formula = y ~x, se = F,lwd=1.5)+
    theme_classic()+
    xlab("Nestedness (NODF)")+
    labs(color= expression(paste(gamma[0])))+
    labs(y = expression(paste("Proportion of pollinator species with ", N[i] > 0.5 ))))# nrow = 2,ncol = 4, labeller = as_labeller(appender, default = label_parsed)))


load("Empirical_data.RData")
#all_dat<-rbind(fact2,fact)


summary(model1<-lmer(mean_visits~nestedness + Treatment + (1|site) + (1|month), data = all_dat))
summary(model2<-lmer(mean_visitation_rate~nestedness + Treatment + (1|site) + (1|month), data = all_dat))



#save(all_dat,file="Empirical_data.RData")
(g1<-all_dat %>%  ggplot(aes(x=nestedness,y=mean_visits,color = Treatment))+
    geom_point(size =5,alpha=0.5 )+
    theme_classic()+
    xlab("Nestedness")+
    ggtitle("Seychelles Data")+
    ylab(expression("Mean pollinator visits per network"))+
    scale_color_manual(values=c(  "#0072B2", "#D55E00"))+
    geom_smooth(method = "lm", se=F))



(g2<-all_dat %>%  ggplot(aes(x=nestedness,y=mean_visitation_rate,color = Treatment))+
    geom_point(size =5,alpha=0.5 )+
    theme_classic()+
    ylim(c(0,6))+
    xlab("Nestedness")+
    ggtitle("Seychelles Data")+
    ylab(expression("Mean visitation rate per network"))+
    scale_color_manual(values=c(  "#0072B2", "#D55E00"))+
    geom_smooth(method = "lm", se=F))


#Figure 1 
ggpubr::ggarrange(w6,w1,w2,
                  a1,g2,g1, labels=c("A","B", "C", "D", "E", "F"), nrow = 2,ncol=3)









#\############################################# proportion exhibited hysteresis : supplementary figure

for(r in 1 :nrow(hyst)){
g<-adj.mat(as.character(webfiles[which(webfiles == hyst$web.name[r])])) #network web names
# g<-g[-1,-1] 


Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant spec

hyst$network.size[r]<- Aspecies+Plantspecies
}




prop_net<-hyst %>% filter(hysteresis== "yes", mut_strength == 4.35)


for(i in 1:nrow(prop_net)){
  
  
  prop<-prop_net$richness[i]/prop_net$network.size[i] 
 if(prop > 0.8){
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
