rm(list=ls())
source("01_functions_ode.R")
library(tidyr) ## for efficient data manipulation & plotting
library(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(viridis)
library(matlib)
library(ggplot2)


theme_set(theme_classic()) 

mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:154]

#make a empty dataframe of all the mutualism networks
fact_2<- expand.grid(`web` =newfiles[1:153]) %>%
  as_tibble %>%
  mutate(`Nestedness`= 0,
         `Connectance` = 0,
         `network_size`=0)
# a loop to go over the rows of the dataframe to calculate nestedness and network size for each network
for (r in 1:nrow(fact_2)){
  g<-adj.mat(myfiles[which(myfiles == fact_2$web[r])]) #network web names
  # g<-g[-1,-1] 
  
  
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1]
  
  fact_2$network_size[r]<-Aspecies+Plantspecies
  fact_2$Nestedness[r]<-  round(nestedness_NODF(g),2)
  fact_2$Connectance[r] <- Connectance(g)
}

fact_2 <- fact_2 %>% filter(network_size < 185 )

webfiles<-as.character(fact_2$web)


fact<- expand.grid(`web` =webfiles,
                   `random_seed`=4327+(1:1)*100) %>%
  as_tibble  


out<-NULL
load("Upper_bound_abundance_response_time.RData")
for(r in 1:nrow(fact)){
  
  g<-adj.mat(myfiles[which(myfiles == fact$web[r])]) #network web names
  # g<-g[-1,-1] 
  
  
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species
  degree.animals<-degree.plants<-numeric()
  
  #degree of plants and anichmals
  for(i in 1:Plantspecies){
    degree.plants[i]<-sum(g[i,])} # degree of plants
  for(j in 1:Aspecies){
    degree.animals[j]<-sum(g[,j]) # degree of animals
  }
  

  degrees<-c(degree.animals,degree.plants)
  
  
  ddf <-as.data.frame(cbind( rep(as.character(fact$web),each=(Aspecies+Plantspecies)), 
                      degrees))
  
  colnames(ddf)<-c("Web","Degree")
  
  out<-rbind(out,ddf)
  
  
  
  print(r)
}

str(out)
out$Degree<- as.numeric(as.character(out$Degree))

webfiles[1:20]
(r1<-out %>% filter(Web ==webfiles[1:20]) %>% 
  ggplot(aes(x=Degree))+
  geom_histogram(bins=40)+
  theme_classic()+
  facet_wrap(.~Web))

(r2<-out %>% filter(Web ==webfiles[21:40]) %>% 
    ggplot(aes(x=Degree))+
    geom_histogram(bins=40)+
    theme_classic()+
    facet_wrap(.~Web))

(r3<-out %>% filter(Web ==webfiles[41:60]) %>% 
    ggplot(aes(x=Degree))+
    geom_histogram(bins=40)+
    theme_classic()+
    facet_wrap(.~Web))

(r4<-out %>% filter(Web ==webfiles[61:101]) %>% 
    ggplot(aes(x=Degree))+
    geom_histogram(bins=40)+
    theme_classic()+
    facet_wrap(.~Web))

