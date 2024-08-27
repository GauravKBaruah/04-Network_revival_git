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


load("figS5.RData")
fact$h2 <- factor(fact$h2, labels = c("heritability=0.0",  
                                      "heritability=0.4"))

(d1<-fact %>% filter(web=="datasets_1/M_PL_045.csv") %>% 
    ggplot(aes(x=Strength_mutualism,y=forcing_strength,z=recovery_richness))+
    geom_raster(aes(fill=recovery_richness),show.legend =TRUE)+ 
    scale_fill_gradient(limits=range(fact$recovery_richness), high = 'yellow', low = 'red')+
    xlab(expression(gamma[0]~",Mutualism strength"))+ylab("Forcing ") + 
    theme_cowplot()+
    facet_wrap(.~h2)+
    labs(fill="Recovery richness"))


d2<-(fact %>% filter(web=="datasets_1/M_PL_060_19.csv") %>% 
       ggplot(aes(x=Strength_mutualism,y=forcing_strength,z=recovery_richness))+
       geom_raster(aes(fill=recovery_richness),show.legend =TRUE)+ 
       scale_fill_gradient(limits=range(fact$recovery_richness), high = 'yellow', low = 'red')+
       xlab(expression(gamma[0]~",Mutualism strength"))+
       ylab("Forcing ") + 
       theme_cowplot()+
       facet_wrap(.~h2)+
       labs(fill="Recovery richness"))


d3<-(fact %>% filter(web=="datasets_1/M_PL_003.csv") %>% 
       ggplot(aes(x=Strength_mutualism,y=forcing_strength,z=recovery_richness))+
       geom_raster(aes(fill=recovery_richness),show.legend =TRUE)+ 
       scale_fill_gradient(limits=range(fact$recovery_richness), high = 'yellow', low = 'red')+
       xlab(expression(gamma[0]~",Mutualism strength"))+
       ylab("Forcing ") + 
       theme_cowplot()+
       facet_wrap(.~h2)+
       labs(fill="Recovery richness"))



g<-adj.mat("datasets_1/M_PL_045.csv") #network web names
Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

net = network(g, directed = FALSE)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net)
net %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
#ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}

deg<-c(degree.plants,degree.animals)
webg<-ggnet2(net, mode="circle", size=deg, 
             edge.size = 1.1,max_size =12, 
             color ="color",edge.alpha = 1, legend.position = "")






g<-adj.mat("datasets_1/M_PL_060_19.csv") #network web names
Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

net = network(g, directed = FALSE)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net)
net %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
#ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}

deg<-c(degree.plants,degree.animals)
webg1<-ggnet2(net, mode="circle", size=deg, 
              edge.size = 1.1,max_size =12, 
              color ="color",edge.alpha = 1, legend.position = "")







g<-adj.mat("datasets_1/M_PL_003.csv") #network web names
Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

net = network(g, directed = FALSE)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net)
net %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
#ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}

deg<-c(degree.plants,degree.animals)
webg3<-ggnet2(net, mode="circle", size=deg, 
              edge.size = 1.1,max_size =12, 
              color ="color",edge.alpha = 1, legend.position = "")





ggpubr::ggarrange(
  webg,d1,
  webg1,d2,
  webg3,d3,
  ncol=2,nrow=3,
  labels = c("A", "B", "C"))


