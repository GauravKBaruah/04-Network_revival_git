rm(list=ls())
source("01_ODE_Function.R", echo=F)
library(foreach)
library(doMC)
library(statmod)
require(deSolve) ## for integrating ordinary differential equations
library(tidyr) ## for efficient data manipulation & plotting
#library(cowplot) ## for arranging plots in a grid
library(dplyr)
library(bipartite)
library(beepr)
#library(viridis)
numcores<- 20
registerDoMC(numcores)



# reading all the datasets
# calculating nestedness and connectance
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
#141

# creating another dataframe where trait values of species will be kept for final simulations
trait_data_list<-list()
trait_data<-expand.grid(web=webfiles) %>% 
  as_tibble  

for (i in 1:nrow(trait_data)){
  
  g<-adj.mat(myfiles[which(myfiles == trait_data$web[i])]) #network web names
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species
  muinit <-runif((Aspecies+Plantspecies), -0.25,0.25)
  trait_data_list[i]<-list(muinit)
  trait_data$networksize[i]<-Aspecies+Plantspecies
  names(trait_data_list)[i] <- as.character(trait_data$web[i])
}

fact<- expand.grid(`Strength_mutualism`=c(0.8, 0.9, 1, 1.1, 1.2),
                   `web` = webfiles[1:92],
                   forcing_strength=0.5,
                   `forcing_duration`=500,
                    h2=0.4,
                    perturbation=c("targeted","degree", "random"),
                   `model`="abundance",
                   `interaction_type`= "trade_off", #no_trade_off
                   `individual.variation` = c("high","low"),
                   `random_seed`=4327+(1:1)*100) %>%
  as_tibble
model.t<-list()

set.seed(1234)

new_ddf<-NULL
load("Upper_bound_abundance_response_time.RData")
outt<-foreach(r = 1:nrow(fact))%dopar%{
  #print(r)
  
  
  g<-adj.mat(myfiles[which(myfiles == fact$web[r])]) #network web names
 Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species
  degree.animals<-degree.plants<-numeric()
  
  #degree of plants and anichmals
  for(i in 1:Plantspecies){
    degree.plants[i]<-sum(g[i,])} # degree of plants
  for(j in 1:Aspecies){
    degree.animals[j]<-sum(g[,j]) # degree of animals
  }
  
  ##control loop for selecting whether variation is high or low
  if(fact$individual.variation[r] == "low"){
    sig <-runif((Aspecies+Plantspecies),0.001,0.001)}else if(fact$individual.variation[r] == "high"){
      sig <-runif((Aspecies+Plantspecies),0.02,0.02)}
  
  
  h2 <- runif((Aspecies+Plantspecies),0.4,0.4)
  
  
  ## vector of species trait standard deviations
  
  index<-which(names(trait_data_list)==fact$web[r])
  
  muinit<- trait_data_list[[index]]        #runif((Aspecies+Plantspecies), -0.25,0.25)
  mainit<-muinit[1:Aspecies]
  mpinit<-muinit[(Aspecies+1): (Aspecies+Plantspecies)]
  
  Amatrix<-mat.comp(g,degree.animals,degree.plants)$Amatrix
  Pmatrix<-mat.comp(g,degree.animals,degree.plants)$Pmatrix
  gamma=0.35#fact_lessvar$Strength_mutualism[r]
  nestedness<-nestedness_NODF(g)
  C<-Connectance(g)
  web.name<-fact$web[r]
  ba<-runif(Aspecies, 0.0,0.0)
  bp<-runif(Plantspecies,0.0,0.00)
  dganimals<-degree.animals
  dgplants<-degree.plants
  
  #species_index_animals<-seq(1,fact$no_species_forced[r],1)
  N <- runif( (Aspecies+Plantspecies) , 0,0.005)
  nainit<-N[1:Aspecies]
  npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
  
  degree_1<-c(degree.animals,degree.plants)
 # index_maxdegree_animal<-which(degree.animals == max(degree.animals))
  #index_maxdegree_plants<-which(degree.plants == max(degree.plants))
  if(fact$perturbation[r] == "degree"){
    
    indexes<- which(degree_1 == max(degree_1))
  }else if(fact$perturbation[r] == "targeted"){
    
    btwn_plantsp<- (bipartite::BC(g, rescale = T, weighted = FALSE))$lower
    btwn_Asp<- (bipartite::BC(g, rescale = T, weighted = FALSE))$higher
    between_centrality <- c(btwn_Asp,btwn_plantsp)
    species_1_index<-which(between_centrality == max(between_centrality))
    second_largest<-sort(between_centrality,decreasing = T)[2]
    species_2_index<-which(between_centrality == second_largest)
    indexes<-c(species_1_index,species_2_index)
  }else{
    indexes<-sample(seq(1:(Aspecies+Plantspecies)),1,replace = F)
}
  
  species_index<-indexes
  
  ba<-runif(Aspecies,0.00,0.00)
  
  # recovery mutualistic network simulation------------------
  mut.strength<-runif((Aspecies+Plantspecies), fact$Strength_mutualism[r],fact$Strength_mutualism[r])
  mut.strength[species_index]<- fact$Strength_mutualism[r] 
  tmax<-900
  time_range<-c(0,tmax)
  deltat<- (time_range[2]-time_range[1])/1 + 1
  duration <- fact$forcing_duration[1]
  d<- c(rep(1,duration),rep(0,(deltat-duration)))
  duration_mat_A<-(replicate(Aspecies,d))
  duration_mat_P<-(replicate(Plantspecies,d))
  times<-seq(0, tmax, 1)
  t1_A<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_P<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_A$import<-duration_mat_A[,1]
  t1_P$import<-duration_mat_P[,1]
  forcing_strength <- rep(fact$forcing_strength[r], (Aspecies+Plantspecies))
  ic_f<-c(nainit, npinit, mainit,mpinit)
  
  params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                         Pmatrix=Pmatrix,w=gamma,
                         ic=ic_f,
                         dat=dat,
                         perturbation=fact$perturbation[r],
                         individual_variation=fact$individual.variation[r],
                         Strength_mutualism=fact$Strength_mutualism[r],
                         mut.strength=mut.strength,m=mainit,C=C,nestedness=nestedness,
                         web.name=web.name,h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                         dgplants=dgplants, forcing_strength=forcing_strength,
                         ua=mainit,up=mpinit,
                         species_index=species_index,
                         interaction_type=fact$interaction_type[r],
                         na=nainit,np=npinit, duration=duration,
                         t1_A=t1_A,t1_P=t1_P)
  
  dat<-cluster_run_func(params_forcing = params_forcing,ic_f = ic_f,tmax = tmax)
  
  
  
}

save(fact, file="S_targeted_degree_data.RData")


load("S_targeted_degree_data.RData")



net_dat<-sp_dat<-NULL
for(i in 1:3312){
 net_dat<-rbind(net_dat, outt[[i]]$output)
  sp_dat<- rbind(sp_dat,outt[[i]]$ddf)
  print(i)
}

#save(net_dat,file="S_targeted_random_degree.RData")
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


str(net_dat)
net_dat$mut_strength<-as.factor(net_dat$mut_strength)
net_dat$individual_variation <- plyr::revalue(net_dat$individual_variation, 
                                              c( "high" = "high variation", "low" = "low variation"))
net_dat$network_size <- as.numeric(as.character(net_dat$network_size))
n1<-(net_dat %>% 
       ggplot(aes(x=Nestedness, y = recovery_richness, 
                  color= factor(perturbation)))+
       geom_point(position=position_jitter(height=0.0,width=0.0),
                  alpha = 0.05, size = 3)+
       scale_color_viridis(discrete = TRUE)+
       ylab("Recovery richness")+
       #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
       labs(color=expression("Recovery perturbation"))+
       # theme(legend.position="none")+
       xlab("Nestedness")+
       # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
       stat_smooth(method = "glm", size=2,alpha=0.1,
                   method.args = list(family = "quasibinomial"),
                   se =FALSE,
                   aes(color=factor(perturbation))) +
       facet_grid(individual_variation~mut_strength))
n1
n2<-(net_dat %>%
       ggplot(aes(x=connectance, y = recovery_richness, color= factor(perturbation)))+
       geom_point(position=position_jitter(height=0.0,width=0.0),
                  alpha = 0.05, size = 3)+
       scale_color_viridis(discrete = TRUE)+
       ylab("Recovery richness")+
       labs(color=expression("Recovery perturbation"))+
       # theme(legend.position="none")+
       xlab("Connectance")+
       # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
       stat_smooth(method = "glm", size=2,alpha=0.1,
                   method.args = list(family = "quasibinomial"),
                   se =FALSE,
                   aes(color=factor(perturbation))) +
       facet_grid(individual_variation~mut_strength))

n3<-(net_dat %>% filter(network_size < 160) %>% 
       ggplot(aes(x=network_size, y = recovery_richness, color= factor(perturbation)))+
       geom_point(position=position_jitter(height=0.0,width=0.0),
                  alpha = 0.05, size = 3)+
       scale_color_viridis(discrete = TRUE)+
       ylab("Recovery richness")+
       labs(color=expression("Recovery perturbation"))+
       # theme(legend.position="none")+
       xlab("Network size")+
       # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
       stat_smooth(method = "glm", size=2,alpha=0.1,
                   method.args = list(family = "quasibinomial"),
                   se =FALSE,
                   aes(color=factor(perturbation))) +
       facet_grid(individual_variation~mut_strength))

ggpubr::ggarrange(n1, n2,labels=c("A","B"),
                  nrow=1)



