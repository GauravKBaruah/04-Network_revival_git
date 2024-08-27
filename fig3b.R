rm(list=ls())
source("01_ODE_Function.R")
library(statmod)
require(deSolve) ## for integrating ordinary differential equations
library(cowplot) ## for arranging plots in a grid
library(dplyr)
library(beepr)
library(viridis)
library(GGally)
library(network)
library(sna)
library(ggplot2)
require(reshape2)
require(flux)
require(akima)


theme_set(theme_classic()) 


mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:154]




# loading the dataset

load("03_figure_data.RData")

web1<- fact %>% filter(web == "datasets_1/M_PL_043.csv")

w1.hvar<-web1 %>% filter(individual.variation =="high")

temp<-with(w1.hvar,interp(w1.hvar$Strength_mutualism,w1.hvar$forcing_strength,
                          w1.hvar$recovery_richness,xo=seq(min(w1.hvar$Strength_mutualism),max(w1.hvar$Strength_mutualism),length=10), 
                          yo=seq(min(w1.hvar$forcing_strength),max(w1.hvar$forcing_strength), length=20), duplicate ='mean'))


temp<-interp2xyz(temp, data.frame=TRUE)
colnames(temp)<-c("x","y","Richness")
hvar_richness_w1<-ggplot(temp , aes(x=x,y=y,z=round(Richness,3)))+
  geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp$Richness), high = 'yellow', low = 'red')+
  xlab("Mutualistic strength")+ylab("Forcing strength") + ggtitle("High variation")+
  theme_cowplot()+
  xlim(c(0.52,1.45))+ ylim(c(0,0.99))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
hvar_richness_w1


lvar<-web1 %>% filter(individual.variation =="low")

akima.lvar<-with(lvar,interp(lvar$Strength_mutualism,lvar$forcing_strength,
                             lvar$recovery_richness,xo=seq(min(lvar$Strength_mutualism),max(lvar$Strength_mutualism),length=50), 
                             yo=seq(min(lvar$forcing_strength),max(lvar$forcing_strength), length=50), duplicate ='mean'))

gdat_lvar<-interp2xyz(akima.lvar, data.frame=TRUE)
colnames(gdat_lvar)<-c("x","y","Richness")
lvar_richness_w1<-ggplot(gdat_lvar , aes(x=x,y=y,z=round(Richness,3)))+geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(gdat_lvar$Richness),high = 'yellow', low = 'red', na.value="NA")+
  xlab("Mutualistic strength")+ylab("Forcing strength") + ggtitle("Low variation")+
  theme_cowplot()+
  xlim(c(0.52,1.45))+ ylim(c(0,0.99))#+scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
lvar_richness_w1


w1<-fact %>% filter(web == 'datasets_1/M_PL_043.csv', forcing_strength == 0) %>% 
  ggplot( aes(x=Strength_mutualism, y = recovery_richness*network.size, color=web ))+
  geom_point(size =4, alpha =0.75)+
  scale_color_manual(values ="#56B4E9")+
  theme_cowplot()+
  ylab("Species richness")+
  xlab("Mutualistic strength")+
  ggtitle("")+
  theme(legend.position="")+ 
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = 0.5, xmax = 1.5,
           ymin = 0, ymax = 105) +
  labs(color=NULL)




newfiles<-myfiles[1:153]

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









##################################WEB NO 2


web2<- fact %>% filter(web == "datasets_1/M_PL_045.csv")

w2.hvar<-web2 %>% filter(individual.variation =="high")

temp<-with(w2.hvar,interp(w2.hvar$Strength_mutualism,w2.hvar$forcing_strength,
                          w2.hvar$recovery_richness,xo=seq(min(w2.hvar$Strength_mutualism),max(w2.hvar$Strength_mutualism),length=50), 
                          yo=seq(min(w2.hvar$forcing_strength),max(w2.hvar$forcing_strength), length=50), duplicate ='mean'))


temp<-interp2xyz(temp, data.frame=TRUE)
colnames(temp)<-c("x","y","Richness")
hvar_richness_w2<-ggplot(temp , aes(x=x,y=y,z=round(Richness,3)))+geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(temp$Richness), high = 'yellow', low = 'red',na.value = "NA")+
  xlab("Mutualistic strength")+ylab("Forcing strength") + ggtitle("")+
  theme_cowplot()+
  xlim(c(0.52,1.45))+ ylim(c(0,0.99))
#scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
hvar_richness_w2


lvar<-web2 %>% filter(individual.variation =="low")

akima.lvar<-with(lvar,interp(lvar$Strength_mutualism,lvar$forcing_strength,
                             lvar$recovery_richness,xo=seq(min(lvar$Strength_mutualism),max(lvar$Strength_mutualism),length=50), 
                             yo=seq(min(lvar$forcing_strength),max(lvar$forcing_strength), length=50), duplicate ='mean'))

gdat_lvar<-interp2xyz(akima.lvar, data.frame=TRUE)
colnames(gdat_lvar)<-c("x","y","Richness")
lvar_richness_w2<-ggplot(gdat_lvar , aes(x=x,y=y,z=round(Richness,3)))+geom_raster(aes(fill=Richness),show.legend =TRUE)+ 
  scale_fill_gradient(limits=range(gdat_lvar$Richness), high = 'yellow', low = 'red', na.value="NA")+
  xlab("Mutualistic strength")+ylab("Forcing strength") + ggtitle("")+
  theme_cowplot()+
  xlim(c(0.52,1.45))+ ylim(c(0,0.99))#+scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
lvar_richness_w2


w2<-fact %>% filter(web == 'datasets_1/M_PL_045.csv', forcing_strength==0) %>% 
  ggplot( aes(x=Strength_mutualism, 
              y = recovery_richness*network.size, color=web ))+
  geom_point(size =4, alpha =1)+
  scale_color_manual(values ="#56B4E9")+
  ylab("Species richness")+
  xlab("Mutualistic strength")+
  
  theme_cowplot()+
  theme(legend.position="")+ 
  annotate("rect", fill = "grey", alpha = 0.5, 
           xmin = 0.5, xmax = 1.5,
           ymin = 0, ymax = 45) +
  labs(color=NULL)




newfiles<-myfiles[1:153]

g<-adj.mat(myfiles[45]) #network web names
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
net_2<-ggnet2(net,  mode="circle",  color ="color", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)


deg<-c(degree.plants,degree.animals)
net_2<-ggnet2(net, mode="circle", size=deg, max_size =8, color ="color",edge.alpha = 1.5, legend.position = "")

####

ggpubr::ggarrange(net_1, hvar_richness_w1, lvar_richness_w1, 
          net_2, hvar_richness_w2, lvar_richness_w2, 
          labels = c("I", "J", "K", "L", "M", "N"),
          ncol = 3, nrow = 2)

