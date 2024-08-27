source("01_functions_ode.R")
library(statmod)
#require(deSolve) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
library(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(viridis)
library(ggdist)
library(ggplot2)
library(akima)
library(network)
library(GGally)
library(sna)

load("fig_S8.RData")

############################## NETWORK NO 1 ####################################
web1<- net_dat %>% filter(web.name  == "datasets_1/M_PL_003.csv", individual_variation == "high")


g<-adj.mat("datasets_1/M_PL_003.csv") #network web names
Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

#degree of plants and anichmals
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}
net = network(g, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net)
net %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
ggnet2(net,  mode="circle",  color ="color", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)


deg<-c(degree.plants,degree.animals)
net_1<-ggnet2(net, mode="circle", size=deg, max_size =8, color ="color",edge.alpha = 1.5, legend.position = "")



temp<-with(web1,interp(web1$mut_strength,web1$forcing_duration,web1$recovery_richness,xo=seq(min(web1$mut_strength),max(web1$mut_strength),length=15), 
                       yo=seq(min(web1$forcing_duration),max(web1$forcing_duration), length=8), duplicate ='mean'))


temp<-interp2xyz(temp, data.frame=TRUE)
colnames(temp)<-c("x","y","Richness")
hvar_richness_w1<-ggplot(temp , aes(x=x,y=y,z=Richness))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp$Richness), high = 'yellow', low = 'red')+
  xlab("Mutualistic strength")+ylab("Forcing duration") + ggtitle("High variation")+
  theme_cowplot()+
  xlim(c(min(web1$mut_strength),max(web1$mut_strength)))+ 
  ylim(c(min(web1$forcing_duration),max(web1$forcing_duration)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
hvar_richness_w1



web2<- net_dat %>% filter(web.name  == "datasets_1/M_PL_003.csv", individual_variation == "low")

temp2<-with(web2,interp(web2$mut_strength,web2$forcing_duration,web2$recovery_richness,xo=seq(min(web2$mut_strength),max(web2$mut_strength),length=15), 
                       yo=seq(min(web2$forcing_duration),max(web2$forcing_duration), length=8), duplicate ='mean'))


temp2<-interp2xyz(temp2, data.frame=TRUE)
colnames(temp2)<-c("x","y","Richness")
lvar_richness_w2<-ggplot(temp2 , aes(x=x,y=y,z=Richness))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp2$Richness), high = 'yellow', low = 'red')+
  xlab("Mutualistic strength")+ylab("Forcing duration") + ggtitle("Low variation")+
  theme_cowplot()+
  xlim(c(min(web2$mut_strength),max(web2$mut_strength)))+ 
  ylim(c(min(web2$forcing_duration),max(web2$forcing_duration)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
lvar_richness_w2




########################################### Network 2 ##########################################

web3<- net_dat %>% filter(web.name  == "datasets_1/M_PL_061_01.csv", individual_variation == "high")


g<-adj.mat("datasets_1/M_PL_061_01.csv") #network web names
Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

#degree of plants and anichmals
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}
net2 = network(g, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net)
net2 %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net2 %v% "color" = ifelse(net2 %v% "groups" == "plants", "#0072B2", "#E69F00" )
ggnet2(net2,  mode="circle",  color ="color", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)


deg<-c(degree.plants,degree.animals)
net_2<-ggnet2(net2, mode="circle", size=deg, max_size =8, color ="color",edge.alpha = 1.5, legend.position = "")




temp3<-with(web3,interp(web3$mut_strength,web3$forcing_duration,web3$recovery_richness,xo=seq(min(web3$mut_strength),max(web3$mut_strength),length=15), 
                        yo=seq(min(web3$forcing_duration),max(web3$forcing_duration), length=8), duplicate ='mean'))


temp3<-interp2xyz(temp3, data.frame=TRUE)
colnames(temp3)<-c("x","y","Richness")
hvar_richness_w3<-ggplot(temp3 , aes(x=x,y=y,z=Richness))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp3$Richness), high = 'yellow', low = 'red')+
  xlab("Mutualistic strength")+ylab("Forcing duration") + ggtitle("High variation")+
  theme_cowplot()+
  xlim(c(min(web3$mut_strength),max(web3$mut_strength)))+ 
  ylim(c(min(web3$forcing_duration),max(web3$forcing_duration)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
hvar_richness_w3


web4<- net_dat %>% filter(web.name  == "datasets_1/M_PL_061_01.csv", individual_variation == "low")

temp4<-with(web4,interp(web4$mut_strength,web4$forcing_duration,web4$recovery_richness,xo=seq(min(web4$mut_strength),max(web4$mut_strength),length=15), 
                        yo=seq(min(web4$forcing_duration),max(web4$forcing_duration), length=8), duplicate ='mean'))


temp4<-interp2xyz(temp4, data.frame=TRUE)
colnames(temp4)<-c("x","y","Richness")
lvar_richness_w4<-ggplot(temp4 , aes(x=x,y=y,z=Richness))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp4$Richness), high = 'yellow', low = 'red')+
  xlab("Mutualistic strength")+ylab("Forcing duration") + ggtitle("Low variation")+
  theme_cowplot()+
  xlim(c(min(web4$mut_strength),max(web4$mut_strength)))+ 
  ylim(c(min(web4$forcing_duration),max(web4$forcing_duration)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
lvar_richness_w4









#########################################    network 3  ############################################################
web5<- net_dat %>% filter(web.name  == "datasets_1/M_PL_043.csv", individual_variation == "high")


g<-adj.mat("datasets_1/M_PL_043.csv") #network web names
Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

#degree of plants and anichmals
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}
net3 = network(g, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net)
net3 %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net3 %v% "color" = ifelse(net3 %v% "groups" == "plants", "#0072B2", "#E69F00" )
ggnet2(net3,  mode="circle",  color ="color", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)


deg<-c(degree.plants,degree.animals)
net_3<-ggnet2(net3, mode="circle", size=deg, max_size =8, color ="color",edge.alpha = 1.5, legend.position = "")




temp5<-with(web5,interp(web5$mut_strength,web5$forcing_duration,web5$recovery_richness,xo=seq(min(web5$mut_strength),max(web5$mut_strength),length=15), 
                        yo=seq(min(web5$forcing_duration),max(web5$forcing_duration), length=8), duplicate ='mean'))


temp5<-interp2xyz(temp5, data.frame=TRUE)
colnames(temp5)<-c("x","y","Richness")
hvar_richness_w5<-ggplot(temp5 , aes(x=x,y=y,z=Richness))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp5$Richness), high = 'yellow', low = 'red')+
  xlab("Mutualistic strength")+ylab("Forcing duration") + ggtitle("High variation")+
  theme_cowplot()+
  xlim(c(min(web5$mut_strength),max(web5$mut_strength)))+ 
  ylim(c(min(web5$forcing_duration),max(web5$forcing_duration)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
hvar_richness_w5




web6<- net_dat %>% filter(web.name  == "datasets_1/M_PL_043.csv", individual_variation == "low")

temp6<-with(web6,interp(web6$mut_strength,web6$forcing_duration,web6$recovery_richness,xo=seq(min(web6$mut_strength),max(web6$mut_strength),length=15), 
                        yo=seq(min(web6$forcing_duration),max(web6$forcing_duration), length=8), duplicate ='mean'))


temp6<-interp2xyz(temp6, data.frame=TRUE)
colnames(temp6)<-c("x","y","Richness")
lvar_richness_w6<-ggplot(temp6 , aes(x=x,y=y,z=Richness))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp6$Richness), high = 'yellow', low = 'red')+
  xlab("Mutualistic strength")+ylab("Forcing duration") + ggtitle("Low variation")+
  theme_cowplot()+
  xlim(c(min(web6$mut_strength),max(web6$mut_strength)))+ 
  ylim(c(min(web6$forcing_duration),max(web6$forcing_duration)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
lvar_richness_w6




library(ggpubr)

ggarrange(net_1,hvar_richness_w1,lvar_richness_w2,
          net_3,hvar_richness_w5,lvar_richness_w6,
          net_2,hvar_richness_w3,lvar_richness_w4, nrow=3,ncol=3)










