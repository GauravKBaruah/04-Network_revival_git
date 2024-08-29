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
theme_set(theme_classic()) 

load("S_targeted_random_data.RData")
############################## NETWORK NO 1  nad perturbation regime of random ####################################
web1<- fact %>% filter(web == "datasets_1/M_PL_060_11.csv",individual.variation=="high",perturbation == "random")


g<-adj.mat( "datasets_1/M_PL_060_11.csv") #network web names
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



temp<-with(web1,interp(web1$Strength_mutualism,web1$forcing_strength,
                       web1$recovery_richness,xo=seq(min(web1$Strength_mutualism),max(web1$Strength_mutualism),length=9), 
                       yo=seq(min(web1$forcing_strength),max(web1$forcing_strength), length=10), duplicate ='mean'))


temp<-interp2xyz(temp, data.frame=TRUE)
colnames(temp)<-c("x","y","Richness")
hvar_richness_w1_random<-ggplot(temp , aes(x=x,y=y,z=round(Richness,3)))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp$Richness), high = 'yellow', low = 'red')+
  xlab(expression( paste(gamma[0],",", "Mutualism strength")))+ylab("Forcing strength") + 
  ggtitle("Random perturbation")+
  theme_cowplot()+
  xlim(c(min(web1$Strength_mutualism),max(web1$Strength_mutualism)))+ 
  ylim(c(min(web1$forcing_strength),max(web1$forcing_strength)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
hvar_richness_w1_random


########################### WITH low trait variation ######################################
web2<- fact %>% filter(web == "datasets_1/M_PL_060_11.csv",individual.variation=="low", perturbation == "random")

temp2<-with(web2,interp(web2$Strength_mutualism,web2$forcing_strength,
                        web2$recovery_richness,xo=seq(min(web2$Strength_mutualism),max(web2$Strength_mutualism),length=9), 
                        yo=seq(min(web2$forcing_strength),max(web2$forcing_strength), length=10), duplicate ='mean'))


temp2<-interp2xyz(temp2, data.frame=TRUE)
colnames(temp2)<-c("x","y","Richness")
lvar_richness_w1_random<-ggplot(temp2 , aes(x=x,y=y,z=round(Richness,3)))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp2$Richness), high = 'yellow', low = 'red')+
  xlab(expression( paste(gamma[0],",", "Mutualism strength")))+ylab("Forcing strength") + 
  ggtitle("Random perturbation")+
  theme_cowplot()+
  xlim(c(min(web2$Strength_mutualism),max(web2$Strength_mutualism)))+ 
  ylim(c(min(web2$forcing_strength),max(web2$forcing_strength)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
lvar_richness_w1_random




#################################################################
############################ 2nd network ########################

web2<- fact %>% filter(web == "datasets_1/M_PL_003.csv",individual.variation=="high", perturbation == "random" )

g<-adj.mat( "datasets_1/M_PL_003.csv") #network web names
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
net_2<-ggnet2(net, mode="circle", size=deg, max_size =8, color ="color",edge.alpha = 1.5, legend.position = "")


temp3<-with(web2,interp(web2$Strength_mutualism,web2$forcing_strength,
                        web2$recovery_richness,xo=seq(min(web2$Strength_mutualism),max(web2$Strength_mutualism),length=9), 
                        yo=seq(min(web2$forcing_strength),max(web2$forcing_strength), length=10), duplicate ='mean'))


temp3<-interp2xyz(temp3, data.frame=TRUE)
colnames(temp3)<-c("x","y","Richness")
hvar_richness_w2_random<-ggplot(temp3 , aes(x=x,y=y,z=round(Richness,3)))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp3$Richness), high = 'yellow', low = 'red')+
  xlab(expression( paste(gamma[0],",", "Mutualism strength")))+ylab("Forcing strength") + 
  ggtitle("Random perturbation")+
  theme_cowplot()+
  xlim(c(min(web2$Strength_mutualism),max(web2$Strength_mutualism)))+ ylim(c(min(web2$forcing_strength),max(web2$forcing_strength)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
hvar_richness_w2_random




web3<- fact %>% filter(web == "datasets_1/M_PL_003.csv",individual.variation=="low", perturbation == "random" )


temp4<-with(web3,interp(web3$Strength_mutualism,web3$forcing_strength,
                        web3$recovery_richness,xo=seq(min(web3$Strength_mutualism),max(web3$Strength_mutualism),length=9), 
                        yo=seq(min(web3$forcing_strength),max(web3$forcing_strength), length=10), duplicate ='mean'))


temp4<-interp2xyz(temp4, data.frame=TRUE)
colnames(temp4)<-c("x","y","Richness")
lvar_richness_w3_random<-ggplot(temp4 , aes(x=x,y=y,z=round(Richness,3)))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp4$Richness), high = 'yellow', low = 'red')+
  xlab(expression( paste(gamma[0],",", "Mutualism strength")))+ylab("Forcing strength") + 
  ggtitle("Random perturbation")+
  theme_cowplot()+
  xlim(c(min(web3$Strength_mutualism),max(web3$Strength_mutualism)))+ ylim(c(min(web3$forcing_strength),max(web3$forcing_strength)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
lvar_richness_w3_random





############################ Targeted forcing #####################################

web4<- fact %>% filter(web ==  "datasets_1/M_PL_060_11.csv", individual.variation == "high", perturbation == "targeted")


temp5<-with(web4,interp(web4$Strength_mutualism,web4$forcing_strength,
                        web4$recovery_richness,xo=seq(min(web4$Strength_mutualism),max(Strength_mutualism),length=10), 
                        yo=seq(min(web4$forcing_strength),max(web4$forcing_strength), length=10), duplicate ='mean'))


temp5<-interp2xyz(temp5, data.frame=TRUE)
colnames(temp5)<-c("x","y","Richness")
hvar_richness_w4_targeted<-ggplot(temp5 , aes(x=x,y=y,z=round(Richness,3)))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp5$Richness), high = 'yellow', low = 'red')+
  xlab(expression( paste(gamma[0],",", "Mutualism strength")))+ylab("Forcing strength") + ggtitle("Targeted perturbation")+
  theme_cowplot()+
  xlim(c(min(web4$Strength_mutualism),max(web4$Strength_mutualism)))+ ylim(c(min(web4$forcing_strength),max(web4$forcing_strength)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
hvar_richness_w4_targeted



## low trait variation



web5<- fact %>% filter(web ==  "datasets_1/M_PL_060_11.csv", individual.variation == "low", perturbation == "targeted")


temp6<-with(web5,interp(web5$Strength_mutualism,web5$forcing_strength,
                        web5$recovery_richness,xo=seq(min(web5$Strength_mutualism),max(Strength_mutualism),length=10), 
                        yo=seq(min(web5$forcing_strength),max(web5$forcing_strength), length=10), duplicate ='mean'))


temp6<-interp2xyz(temp6, data.frame=TRUE)
colnames(temp6)<-c("x","y","Richness")
lvar_richness_w5_targeted<-ggplot(temp6 , aes(x=x,y=y,z=round(Richness,3)))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp6$Richness), high = 'yellow', low = 'red')+
  xlab(expression( paste(gamma[0],",", "Mutualism strength")))+ylab("Forcing strength") + ggtitle("Targeted perturbation")+
  theme_cowplot()+
  xlim(c(min(web4$Strength_mutualism),max(web4$Strength_mutualism)))+ ylim(c(min(web4$forcing_strength),max(web4$forcing_strength)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
lvar_richness_w5_targeted


################################# Network 2  targeted perturbation ###############################



web6<- fact %>% filter(web == "datasets_1/M_PL_003.csv",individual.variation=="high", perturbation == "targeted" )


temp7<-with(web6,interp(web6$Strength_mutualism,web6$forcing_strength,
                        web6$recovery_richness,xo=seq(min(web6$Strength_mutualism),max(Strength_mutualism),length=10), 
                        yo=seq(min(web6$forcing_strength),max(web6$forcing_strength), length=10), duplicate ='mean'))


temp7<-interp2xyz(temp7, data.frame=TRUE)
colnames(temp7)<-c("x","y","Richness")
hvar_richness_w6_targeted<-ggplot(temp7 , aes(x=x,y=y,z=round(Richness,3)))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp7$Richness), high = 'yellow', low = 'red')+
  xlab(expression( paste(gamma[0],",", "Mutualism strength")))+ylab("Forcing strength") + ggtitle("Targeted perturbation")+
  theme_cowplot()+
  xlim(c(min(web4$Strength_mutualism),max(web4$Strength_mutualism)))+ ylim(c(min(web4$forcing_strength),max(web4$forcing_strength)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
hvar_richness_w6_targeted

#low trait variation


web7<- fact %>% filter(web == "datasets_1/M_PL_003.csv",individual.variation=="low", perturbation == "targeted" )


temp8<-with(web7,interp(web7$Strength_mutualism,web7$forcing_strength,
                        web7$recovery_richness,xo=seq(min(web7$Strength_mutualism),max(Strength_mutualism),length=10), 
                        yo=seq(min(web7$forcing_strength),max(web7$forcing_strength), length=10), duplicate ='mean'))


temp8<-interp2xyz(temp8, data.frame=TRUE)
colnames(temp8)<-c("x","y","Richness")
lvar_richness_w6_targeted<-ggplot(temp8 , aes(x=x,y=y,z=round(Richness,3)))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp8$Richness), high = 'yellow', low = 'red')+
  xlab(expression( paste(gamma[0],",", "Mutualism strength")))+ylab("Forcing strength") + ggtitle("Targeted perturbation")+
  theme_cowplot()+
  xlim(c(min(web4$Strength_mutualism),max(web4$Strength_mutualism)))+ ylim(c(min(web4$forcing_strength),max(web4$forcing_strength)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
lvar_richness_w6_targeted








pdf(file = "s_targeted_degree_example_high_var.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 10) # The height of the plot in inches

ggpubr::ggarrange(net_1, hvar_richness_w1_random,
                  net_1, hvar_richness_w4_targeted,
                  net_2, hvar_richness_w2_random,
                  net_2,hvar_richness_w6_targeted,
                  
                  labels = c("A","B","C","D","E","F","G","H"),
                  nrow=4,ncol=2)

dev.off()


