rm(list=ls())
source("01_ODE_Function.R")
require(deSolve) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
require(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(viridis)

library(network)
library(sna)
library(ggplot2)
library(ggnet)

theme_set(theme_classic()) 


# reading all the datasets
# calculating nestedness and connectance
mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:154]





#newfiles<-ex_webs
fact<- expand.grid(`Strength_mutualism`=1.15, 
                 `web` = c(newfiles[70],newfiles),
                   `forcing_strength`=0.5,
                   `forcing_duration`=500,
                   `no_species_forced`=1, #no_trade_off
                 interaction_type = "trade_off",
                   `individual.variation` = c("high"),
                   `random_seed`=4327+(1:1)*100) %>%
  as_tibble

trait_data_list<-list()
trait_data<-expand.grid(web=c(newfiles[1:152])) %>%
  as_tibble

for(i in 1:nrow(trait_data)){

  g<-adj.mat(myfiles[which(myfiles == trait_data$web[i])]) #network web names
  # g<-g[-1,-1]


  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species



  muinit <-runif((Aspecies+Plantspecies), -0.25,0.25)

  trait_data_list[i]<-list(muinit)
  trait_data$network_size[i]<- Aspecies+Plantspecies
  names(trait_data_list)[i] <- as.character(trait_data$web[i])

}

load("Mean_trait_data.RData")
set.seed(1234)
new_ddf<-NULL
load("Upper_bound_abundance_response_time.RData")


  
  g<-adj.mat(myfiles[which(myfiles == fact$web[1])]) #network web names
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
  
  ##control loop for selecting whether variation is high or low
  if(fact$individual.variation[1] == "low"){
    sig <-runif((Aspecies+Plantspecies),0.005,0.005)}else if(fact$individual.variation[1] == "high"){
      sig <-runif((Aspecies+Plantspecies),0.02,0.02)}
  

      h2 <- runif((Aspecies+Plantspecies),0.4,0.4)
  
  
  ## vector of species trait standard deviations
  
  N <- runif( (Aspecies+Plantspecies) , 1,1)
  nainit<- N[1:Aspecies]
  npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
  
  index<-which(names(trait_data_list)==fact$web[1])
  
  muinit<-trait_data_list[[index]]
  mainit<-muinit[1:Aspecies]
  mpinit<-muinit[(Aspecies+1): (Aspecies+Plantspecies)]
  
  Amatrix<-mat.comp(g,degree.animals,degree.plants)$Amatrix
  Pmatrix<-mat.comp(g,degree.animals,degree.plants)$Pmatrix
  gamma=0.35#fact_lessvar$Strength_mutualism[r]
  nestedness<-nestedness_NODF(g)
  C<-Connectance(g)
  web.name<-fact$web[1]
  ba<-runif(Aspecies, 0.00,0.00)
  bp<-runif(Plantspecies,0.00,0.00)
  dganimals<-degree.animals
  dgplants<-degree.plants
  mut.strength<-runif((Aspecies+Plantspecies), 1.15,1.15)
  ic <-c(nainit, npinit, mainit,mpinit)
 tmax=1000
  
  params_noforcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                 Pmatrix=Pmatrix,w=gamma,
                 ic=ic,
                 dat=dat,
                 individual_variation=fact$individual.variation[1],
                 mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                 web.name=web.name,h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                 dgplants=dgplants, forcing_strength=0,
                 ua=mainit,up=mpinit,
                 interaction_type=fact$interaction_type[1],
                 na=nainit,np=npinit, duration=0)
  
  
  sol1<-ode(func=eqs, y=ic, parms=params_noforcing, times=seq(0, tmax, by=1)) %>% 
    organize_results(pars = params_noforcing) 
    
  
    (netwr_no_forcing_hvar_n1<-sol1 %>% filter(type != "ma", type != "mp") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 0.5, label = "gamma[0] == 1.15",
             parse = TRUE)+
        annotate("text", x = 465, y = 0.505, label = ", with no forcing",
                 parse = F) )
    
  ##### the simulation with with forcing  
  N <- runif( (Aspecies+Plantspecies) , 0,0.005) #mimicking a system with low population density near the collapse state
  nainit<-N[1:Aspecies]
  npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
  
 
  mainit<- (sol1 %>% filter(type=="ma",time ==tmax))$v
  mpinit<-(sol1 %>% filter(type=="mp", time ==tmax))$v
  
  index_maxdegree_animal<-which(degree.animals == max(degree.animals))
  index_maxdegree_plants<-which(degree.plants == max(degree.plants))
  
  degree.vector<-c(degree.animals,degree.plants)
  index_max_degree<-which(degree.vector == max(degree.vector))
  
  species_index<-index_max_degree
  
  mut.strength<-runif((Aspecies+Plantspecies), fact$Strength_mutualism[1],fact$Strength_mutualism[1])
  mut.strength[index_max_degree]<- fact$Strength_mutualism[1] 
  time_range<-c(0,tmax)
  deltat<- (time_range[2]-time_range[1])/1 + 1
  duration <- 500#fact$forcing_duration[1]
  d<- c(rep(1,duration),rep(0,(deltat-duration)))
  duration_mat_A<-(replicate(Aspecies,d))
  duration_mat_P<-(replicate(Plantspecies,d))
  times<-seq(0, tmax, 1)
  t1_A<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_P<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_A$import<-duration_mat_A[,1]
  t1_P$import<-duration_mat_P[,1]
  #time_func_A<-approxfun(t1_A, method="linear", rule =2)
  #time_func_P<-approxfun(t1_P, method="linear", rule =2)
  forcing_strength <- rep(fact$forcing_strength[1], (Aspecies+Plantspecies))
  ic_f<-c(nainit, npinit, mainit,mpinit)
  
  params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                           Pmatrix=Pmatrix,w=gamma,
                           ic=ic_f,
                           dat=dat,
                           individual_variation=fact$individual.variation[1],
                           mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                           web.name=web.name,h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                           dgplants=dgplants, forcing_strength=forcing_strength,
                           ua=mainit,up=mpinit,
                           species_index=species_index,
                          interaction_type=fact$interaction_type[1],
                           na=nainit,np=npinit, duration=duration,
                         t1_A=t1_A,t1_P=t1_P)
  
  
  sol<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
    organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs
  
    (netwr_forcing_hvar_n1<-sol %>% filter(type != "ma", type != "mp") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme(legend.position="none")+
        annotate("text", x = 250, y = 1.5, label = "gamma[0] == 1.15",
                 parse = TRUE)+
        annotate("text", x = 465, y = 1.505, label = ", with forcing",
                 parse = F)+
    annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2,
             alpha = .2))

  
################################for second network plotting for figure 2################################
#########################################################################################################
    
    
    (g<-adj.mat(myfiles[which(myfiles == fact$web[51])])) #network web names
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
    
    ##control loop for selecting whether variation is high or low
    if(fact$individual.variation[1] == "low"){
      sig <-runif((Aspecies+Plantspecies),0.005,0.005)}else if(fact$individual.variation[1] == "high"){
        sig <-runif((Aspecies+Plantspecies),0.02,0.02)}
    
    
    h2 <- runif((Aspecies+Plantspecies),0.4,0.4)
    
    
    ## vector of species trait standard deviations
    
    N <- runif( (Aspecies+Plantspecies) , 1,1)
    nainit<- N[1:Aspecies]
    npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
    
    index<-which(names(trait_data_list)==fact$web[51])
    
    muinit<-trait_data_list[[index]]
    mainit<-muinit[1:Aspecies]
    mpinit<-muinit[(Aspecies+1): (Aspecies+Plantspecies)]
    
    Amatrix<-mat.comp(g,degree.animals,degree.plants)$Amatrix
    Pmatrix<-mat.comp(g,degree.animals,degree.plants)$Pmatrix
    gamma=0.35#fact_lessvar$Strength_mutualism[r]
    nestedness<-nestedness_NODF(g)
    C<-Connectance(g)
    web.name<-fact$web[2]
    ba<-runif(Aspecies, 0.00,0.00)
    bp<-runif(Plantspecies,0.00,0.00)
    dganimals<-degree.animals
    dgplants<-degree.plants
    mut.strength<-runif((Aspecies+Plantspecies), 1.15,1.15)
    ic <-c(nainit, npinit, mainit,mpinit)
    tmax =1000
    
    params_noforcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                             Pmatrix=Pmatrix,w=gamma,
                             ic=ic,
                             dat=dat,
                             individual_variation=fact$individual.variation[1],
                             mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                             web.name=web.name,h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                             dgplants=dgplants, forcing_strength=0,
                             ua=mainit,up=mpinit,
                             interaction_type=fact$interaction_type[1],
                             na=nainit,np=npinit, duration=0)
    
    
    sol2<-ode(func=eqs, y=ic, parms=params_noforcing, times=seq(0, tmax, by=1)) %>% 
      organize_results(pars = params_noforcing) 
    
    
    (netwr_no_forcing_hvar_n2<-sol2 %>% filter(type != "ma", type != "mp") %>% 
      ggplot +
      geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
      scale_y_continuous(name="population density", limits=c(0, NA)) +
      scale_color_viridis_d()+
      theme(legend.position="none")+
        annotate("text", x = 250, y = 0.5, label = "gamma[0] == 1.15",
                 parse = TRUE)+
        annotate("text", x = 465, y = 0.505, label = ", with no forcing",
                 parse = F)) #annotate("rect", xmin = 0, xmax = 300, ymin = 0, ymax = 0.1,
    #        alpha = .2)
    #%>% filter(time == 500)  #%>% plot_all() ## solve ODEs
    
    N <- runif( (Aspecies+Plantspecies) , 0,0.005)
    nainit<-N[1:Aspecies]
    npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
    
    
    mainit<- (sol2 %>% filter(type=="ma",time ==tmax))$v
    mpinit<-(sol2 %>% filter(type=="mp", time ==tmax))$v
    
    index_maxdegree_animal<-which(degree.animals == max(degree.animals))
    index_maxdegree_plants<-which(degree.plants == max(degree.plants))
    
    degree.vector<-c(degree.animals,degree.plants)
    index_max_degree<-which(degree.vector == max(degree.vector))
    
    species_index<-index_max_degree
    
    mut.strength<-runif((Aspecies+Plantspecies), fact$Strength_mutualism[1],fact$Strength_mutualism[1])
    mut.strength[index_max_degree]<- fact$Strength_mutualism[1] 
    time_range<-c(0,tmax)
    deltat<- (time_range[2]-time_range[1])/1 + 1
    duration <- 500#fact$forcing_duration[1]
    d<- c(rep(1,duration),rep(0,(deltat-duration)))
    duration_mat_A<-(replicate(Aspecies,d))
    duration_mat_P<-(replicate(Plantspecies,d))
    times<-seq(0, tmax, 1)
    t1_A<-as.data.frame(list(times=times, import = rep(0,length(times))))
    t1_P<-as.data.frame(list(times=times, import = rep(0,length(times))))
    t1_A$import<-duration_mat_A[,1]
    t1_P$import<-duration_mat_P[,1]
    time_func_A<-approxfun(t1_A, method="linear", rule =2)
    time_func_P<-approxfun(t1_P, method="linear", rule =2)
    forcing_strength <- rep(fact$forcing_strength[1], (Aspecies+Plantspecies))
    ic_f<-c(nainit, npinit, mainit,mpinit)
    
    params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                           Pmatrix=Pmatrix,w=gamma,
                           ic=ic_f,
                           dat=dat,
                           individual_variation=fact$individual.variation[1],
                           mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                           web.name=web.name,h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                           dgplants=dgplants, forcing_strength=forcing_strength,
                           ua=mainit,up=mpinit,
                           species_index=species_index,
                           interaction_type=fact$interaction_type[1],
                           na=nainit,np=npinit, duration=duration,
                           t1_A=t1_A,t1_P=t1_P)
    
    
    sol3<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
      organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs
    
    (netwr_forcing_hvar_n3<-sol3 %>% filter(type != "ma", type != "mp") %>% 
      ggplot +
      geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
      scale_y_continuous(name="population density", limits=c(0, NA)) +
      scale_color_viridis_d()+
      theme(legend.position="none")+
        annotate("text", x = 250, y = 1.5, label = "gamma[0] == 1.15",
                 parse = TRUE)+
        annotate("text", x = 465, y = 1.505, label = ", with forcing",
                 parse = F)+
      annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2,
               alpha = .2))
    
    
    
###### plotting the two networks and every population dynamics together############################################
################################################################################################################    
    
    
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
    
    #vertex names 
    net = network(g, bipartite = T, directed = FALSE)
    
    names<-network.vertex.names(net)
    net %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
    net %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
    #ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)
    
    deg<-c(degree.plants,degree.animals)
    (w1<-ggnet2(net, mode="circle", size=deg, 
               edge.size = 1.1,max_size =12, 
               color ="color",edge.alpha = 1, legend.position = "",label = "A"))
    w1
    
    
################ network no. 2 #####################
    
    g1<-adj.mat(myfiles[which(myfiles == fact$web[51])]) #network web names
    Aspecies<- dim(g1)[2] # no of animal species
    Plantspecies<- dim(g1)[1] # no of plant species
    degree.animals<-degree.plants<-numeric()
    
    #degree of plants and anichmals
    for(i in 1:Plantspecies){
      degree.plants[i]<-sum(g1[i,])} # degree of plants
    for(j in 1:Aspecies){
      degree.animals[j]<-sum(g1[,j]) # degree of animals
    }
    
    #vertex names 
    net = network(g1, bipartite = T, directed = FALSE)
    
    names<-network.vertex.names(net)
    net %v% "groups" = ifelse( names[1:sum(dim(g1))] %in% c( as.character(seq(1:dim(g1)[1])) ), "plants", "animals")
    net %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
    #ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)
    
    deg<-c(degree.plants,degree.animals)
    (w2<-ggnet2(net, mode="circle", size=deg, 
                edge.size = 1.1,max_size =12, 
                color ="color",edge.alpha = 1, 
                legend.position = "",label = "A"))
    w2
    
    
    
    ggpubr::ggarrange(w1,
                        netwr_no_forcing_hvar_n1,
                        netwr_forcing_hvar_n1,
                        w2,
                        netwr_no_forcing_hvar_n2,
                        netwr_forcing_hvar_n3, 
                        nrow=2, ncol=3,
                      widths = c(0.7,1,1),
                       
                        labels = c("A", "B", "C",
                                   "D","E", "F"))
    
