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


# reading all the datasets
# calculating nestedness and connectance
mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:154]



load("figS9.RData")

############################## NETWORK NO 1 ####################################
web1<- fact %>% filter(web == "datasets_1/M_PL_060_11.csv")


g<-adj.mat( newfiles[70]) #network web names
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



temp<-with(web1,interp(web1$h2,web1$h2_sd,
                       web1$recovery_richness,xo=seq(min(web1$h2),max(web1$h2),length=9), 
                          yo=seq(min(web1$h2_sd),max(web1$h2_sd), length=10), duplicate ='mean'))


temp<-interp2xyz(temp, data.frame=TRUE)
colnames(temp)<-c("x","y","Richness")
hvar_richness_w1<-ggplot(temp , aes(x=x,y=y,z=round(Richness,3)))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp$Richness), high = 'yellow', low = 'red')+
  xlab("Mean heritability")+ylab("sd in heritability") + ggtitle("High variation")+
  theme_cowplot()+
  xlim(c(min(web1$h2),max(web1$h2)))+ ylim(c(min(web1$h2_sd),max(web1$h2_sd)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
hvar_richness_w1





############################ 2nd network ########################

web2<- fact %>% filter(web == "datasets_1/M_PL_003.csv")

g<-adj.mat( newfiles[3]) #network web names
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


temp2<-with(web2,interp(web2$h2,web2$h2_sd,
                       web2$recovery_richness,xo=seq(min(web2$h2),max(web2$h2),length=9), 
                       yo=seq(min(web2$h2_sd),max(web2$h2_sd), length=10), duplicate ='mean'))


temp2<-interp2xyz(temp2, data.frame=TRUE)
colnames(temp2)<-c("x","y","Richness")
hvar_richness_w2<-ggplot(temp2 , aes(x=x,y=y,z=round(Richness,3)))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp2$Richness), high = 'yellow', low = 'red')+
  xlab("Mean heritability")+ylab("sd in heritability") + ggtitle("High variation")+
  theme_cowplot()+
  xlim(c(min(web2$h2),max(web2$h2)))+ ylim(c(min(web2$h2_sd),max(web2$h2_sd)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
hvar_richness_w2




############################ 3rd network #####################################

web3<- fact %>% filter(web == "datasets_1/M_PL_030.csv")

g<-adj.mat( newfiles[30]) #network web names
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
net_3<-ggnet2(net, mode="circle", size=deg, max_size =8, color ="color",edge.alpha = 1.5, legend.position = "")


temp3<-with(web3,interp(web3$h2,web3$h2_sd,
                        web3$recovery_richness,xo=seq(min(web3$h2),max(web3$h2),length=10), 
                        yo=seq(min(web3$h2_sd),max(web3$h2_sd), length=10), duplicate ='mean'))


temp3<-interp2xyz(temp3, data.frame=TRUE)
colnames(temp3)<-c("x","y","Richness")
hvar_richness_w3<-ggplot(temp3 , aes(x=x,y=y,z=round(Richness,3)))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp3$Richness), high = 'yellow', low = 'red')+
  xlab("Mean heritability")+ylab("sd in heritability") + ggtitle("High variation")+
  theme_cowplot()+
  xlim(c(min(web3$h2),max(web3$h2)))+ ylim(c(min(web3$h2_sd),max(web3$h2_sd)))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
hvar_richness_w3

ggpubr::ggarrange(net_1, hvar_richness_w1,
                  net_2, hvar_richness_w2,
                  net_3, hvar_richness_w3,
                  labels = c("A","B","C","D","E","F"),
                  nrow=3,ncol=2)


