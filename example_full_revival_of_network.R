rm(list=ls())

#load libraries and functions related to this R script.

source("01_ODE_Function.R")
library(statmod)
require(deSolve) 
library(dplyr)
library(beepr)
library(GGally)
library(network)
library(sna)

set.seed(124)

#set directory of all the empirical plant-pollinator network incidence data matrices
mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:153]

#create an empty data frame where different parameters could be varied and state variables could be stored.
fact<- expand.grid(`Strength_mutualism`=4, 
                   `web` = "datasets_1/M_PL_003.csv",
                   `forcing_strength`= 0,
                   `forcing_duration`=0,
                   `no_species_forced`=0,
                   `model`="abundance",
                   `interaction_type`= "trade_off", 
                   `individual.variation` ="high",
                   `random_seed`=4327+(1:1)*100) %>%
  as_tibble()

#loop over the empty dataframe

  
  #adjacency matrix
  g<-adj.mat(myfiles[which(myfiles == fact$web[1])]) #network web names
  
  
  
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
  sig <-runif((Aspecies+Plantspecies),0.005,0.005)
  N <- runif( (Aspecies+Plantspecies) , 1,1) 
  nainit<-N[1:Aspecies]
  npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
  
  #randomly sample mean trait values bounded between -0.5,0.5
  ma<-  runif((Aspecies+Plantspecies), -0.5, 0.5) 
  mainit<-ma[1:Aspecies]
  mpinit<-ma[(Aspecies+1): (Aspecies+Plantspecies)]
  
  index_maxdegree_animal<-which(degree.animals == max(degree.animals))
  index_maxdegree_plants<-which(degree.plants == max(degree.plants))
  
  degree.vector<-c(degree.animals,degree.plants)
  index_max_degree<-which(degree.vector == max(degree.vector))
  
  species_index<-index_max_degree
  Amatrix<-mat.comp(g,degree.animals,degree.plants)$Amatrix
  Pmatrix<-mat.comp(g,degree.animals,degree.plants)$Pmatrix
  gamma=0.35
  tmax<-1000 #total time point for simulation
  nestedness<-nestedness_NODF(g) #nestedness calculation
  C<-Connectance(g) #connectance calculation
  web.name<-fact$web[1] #name of the web
  ba<-runif(Aspecies, 0.,0.) #all zero growth rate, b
  bp<-runif(Plantspecies,0.,0.) # zero growth rate for plants, b
  dganimals<-degree.animals
  dgplants<-degree.plants
  fact$Strength_mutualism[1]<- 4 #gamma_0 or average mutualistic strength
  mut.strength<-runif((Aspecies+Plantspecies), fact$Strength_mutualism[1],fact$Strength_mutualism[1])
  mut.strength[index_max_degree]<- fact$Strength_mutualism[1] 
  
  #this below lines are if you want to add species specific perturbation, which we dont do in this Rscript  
  time_range<-c(0,tmax)
  deltat<- (time_range[2]-time_range[1])/1 + 1
  duration <- fact$forcing_duration[1]<- 0 #duration of the perturbation which is zero here.
  d<- c(rep(1,duration),rep(0,(deltat-duration)))
  duration_mat_A<-(replicate(Aspecies,d))
  duration_mat_P<-(replicate(Plantspecies,d))
  times<-seq(0, tmax, 1)
  t1_A<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_P<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_A$import<-duration_mat_A[,1]
  t1_P$import<-duration_mat_P[,1]
  fact$forcing_strength[1]<-0 #forcing strenght is zero here in this example hystereis simulation
  forcing_strength <- rep(0, (Aspecies+Plantspecies))# zero species specific forcing strength here.
  ic_f<-c(nainit, npinit, mainit,mpinit) #initial N, and mean trait values
  
  #list of variables
  params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                         Pmatrix=Pmatrix,w=gamma,
                         ic=ic_f,
                         individual_variation=fact$individual.variation[1],
                         mut.strength=mut.strength,m=mainit,C=C,nestedness=nestedness,
                         web.name=web.name,h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                         dgplants=dgplants, forcing_strength=forcing_strength,
                         ua=mainit,up=mpinit,
                         species_index=species_index,
                         interaction_type=fact$interaction_type[1],
                         na=nainit,np=npinit, duration=duration,
                         t1_A=t1_A,t1_P=t1_P)
  
  #solving eco-evo dynamics
  sol1<-ode(func=eqs, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
    organize_results(pars = params_forcing)
  
 final_mean_trait_values<- (sol1 %>% filter(time == tmax, type %in%c("ma", "mp")) )$v
  #plotting N for all speices
 (density<-sol1 %>% filter(type %in%c("N", "P")) %>% 
     ggplot +
     geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
     scale_y_continuous(name="population density", limits=c(0, NA)) +
     scale_color_viridis_d()+
     theme_classic()+
     theme(legend.position="none"))
 


 ######################## collapse of plant-pollinator networks at low average mutualistic strength - low variation ##############



 
 N <- runif( (Aspecies+Plantspecies) , 1,1) # #mimicking a system with low population density near the collapse state
 
 nainit<-N[1:Aspecies]
 npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
 
 sig <-runif((Aspecies+Plantspecies),0.02,0.02)
 ma<- final_mean_trait_values
 mainit<-ma[1:Aspecies]
 mpinit<-ma[(Aspecies+1): (Aspecies+Plantspecies)]
 
 
 fact$Strength_mutualism[1]<- 1.225 #gamma_0 or average mutualistic strength
 mut.strength<-runif((Aspecies+Plantspecies), fact$Strength_mutualism[1],fact$Strength_mutualism[1])
 mut.strength[index_max_degree]<- fact$Strength_mutualism[1] 
 
 tmax<-1000
 #this below lines are if you want to add species specific perturbation, which we dont do in this Rscript  
 time_range<-c(0,tmax)
 deltat<- (time_range[2]-time_range[1])/1 + 1
 duration <- fact$forcing_duration[1]<- 0 #duration of the perturbation which is zero here.
 d<- c(rep(1,duration),rep(0,(deltat-duration)))
 duration_mat_A<-(replicate(Aspecies,d))
 duration_mat_P<-(replicate(Plantspecies,d))
 times<-seq(0, tmax, 1)
 t1_A<-as.data.frame(list(times=times, import = rep(0,length(times))))
 t1_P<-as.data.frame(list(times=times, import = rep(0,length(times))))
 t1_A$import<-duration_mat_A[,1]
 t1_P$import<-duration_mat_P[,1]
 fact$forcing_strength[1]<-0 #forcing strenght is zero here in this example hystereis simulation
 forcing_strength <- rep(0, (Aspecies+Plantspecies))# zero species specific forcing strength here.
 ic_f<-c(nainit, npinit, mainit,mpinit) #initial N, and mean trait values
 
 #list of variables
 params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                        Pmatrix=Pmatrix,w=gamma,
                        ic=ic_f,
                        individual_variation=fact$individual.variation[1],
                        mut.strength=mut.strength,m=mainit,C=C,nestedness=nestedness,
                        web.name=web.name,h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                        dgplants=dgplants, forcing_strength=forcing_strength,
                        ua=mainit,up=mpinit,
                        species_index=species_index,
                        interaction_type=fact$interaction_type[1],
                        na=nainit,np=npinit, duration=duration,
                        t1_A=t1_A,t1_P=t1_P)
 

 #solving eco-evo dynamics - no perturbation- zero forcing
 sol2<-ode(func=eqs, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
   organize_results(pars = params_forcing)
 
  #plotting N for all speices
 (netwr_no_forcing_lvar<-sol2 %>% filter(type != "ma", type != "mp") %>% 
     ggplot +
     geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
     scale_y_continuous(name="population density", limits=c(0, NA)) +
     scale_color_viridis_d()+
     theme_classic()+
     theme(legend.position="none")+
     annotate("text", x = 250, y = 1.5, label = "gamma[0] == 1.2 + no_forcing",
              parse = TRUE)+
     annotate("text", x = 1100, y = 1.9, label = "sigma[i] == 0.02",
              parse = TRUE,size=4 , color="darkgreen")+
     annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2,
              alpha = .2))
 
 
 
######################## Revival of plant-pollinator networks at low average mutualistic strength - High variation + forcing ##############
#################################### forcing strength = 0.5, for duration = 500 #####################################
 
 N <- runif( (Aspecies+Plantspecies) , 0,0.005) # #mimicking a system with low population density near the collapse state
 nainit<-N[1:Aspecies]
 npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
 
 sig <-runif((Aspecies+Plantspecies),0.02,0.02)
 ma<-  (sol1 %>% filter(time == 1000, type %in%c("ma", "mp")) )$v
 mainit<-ma[1:Aspecies]
 mpinit<-ma[(Aspecies+1): (Aspecies+Plantspecies)]
 
 
 fact$Strength_mutualism[1]<- 1.225 #gamma_0 or average mutualistic strength
 mut.strength<-runif((Aspecies+Plantspecies), fact$Strength_mutualism[1],fact$Strength_mutualism[1])
 mut.strength[index_max_degree]<- fact$Strength_mutualism[1] 
 
 tmax<-1200
 #this below lines are if you want to add species specific perturbation, which we dont do in this Rscript  
 time_range<-c(0,tmax)
 deltat<- (time_range[2]-time_range[1])/1 + 1
 duration <- fact$forcing_duration[1]<- 500 #duration of the perturbation which is zero here.
 d<- c(rep(1,duration),rep(0,(deltat-duration)))
 duration_mat_A<-(replicate(Aspecies,d))
 duration_mat_P<-(replicate(Plantspecies,d))
 times<-seq(0, tmax, 1)
 t1_A<-as.data.frame(list(times=times, import = rep(0,length(times))))
 t1_P<-as.data.frame(list(times=times, import = rep(0,length(times))))
 t1_A$import<-duration_mat_A[,1]
 t1_P$import<-duration_mat_P[,1]
 fact$forcing_strength[1]<-0.6 #forcing strenght is zero here in this example hystereis simulation
 forcing_strength <- rep(0.6, (Aspecies+Plantspecies))# zero species specific forcing strength here.
 ic_f<-c(nainit, npinit, mainit,mpinit) #initial N, and mean trait values
 
 #list of variables
 params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                        Pmatrix=Pmatrix,w=gamma,
                        ic=ic_f,
                        individual_variation=fact$individual.variation[1],
                        mut.strength=mut.strength,m=mainit,C=C,nestedness=nestedness,
                        web.name=web.name,h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                        dgplants=dgplants, forcing_strength=forcing_strength,
                        ua=mainit,up=mpinit,
                        species_index=species_index,
                        interaction_type=fact$interaction_type[1],
                        na=nainit,np=npinit, duration=duration,
                        t1_A=t1_A,t1_P=t1_P)
 
 
 #solving eco-evo dynamics - no perturbation- zero forcing
 sol3<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
   organize_results(pars = params_forcing)
 
 #plotting N for all speices
 (netwr_forcing_lvar<-sol3 %>% filter(type != "ma", type != "mp") %>% 
     ggplot +
     geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
     scale_y_continuous(name="population density", limits=c(0, NA)) +
     scale_color_viridis_d()+
     theme_classic()+
     theme(legend.position="none")+
     annotate("text", x = 250, y = 1.5, label = "gamma[0] == 1.2+ forcing",
              parse = TRUE)+
     annotate("text", x = 1300, y = 1.9, label = "sigma[i] == 0.02",
              parse = TRUE,size=4 , color="darkgreen")+
     annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2,
              alpha = .2))

 