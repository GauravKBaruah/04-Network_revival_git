rm(list=ls())
source("01_functions_ode.R")
library(statmod)
require(deSolve) ## for integrating ordinary differential equations
library(tidyr) ## for efficient data manipulation & plotting
library(cowplot) ## for arranging plots in a grid
library(dplyr)
#library(readr)
library(beepr)
library(viridis)
library(ggplot2)

library(network)
library(sna)
library(GGally)

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

webfiles<-fact_2$web





  g<-adj.mat(myfiles[which(myfiles == "datasets_1/M_PL_061_01.csv")]) #network web names
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species
  degree.animals<-degree.plants<-numeric()
  
  #degree of plants and anichmals
  for(i in 1:Plantspecies){
    degree.plants[i]<-sum(g[i,])} # degree of plants
  for(j in 1:Aspecies){
    degree.animals[j]<-sum(g[,j]) # degree of animals
  }
  sig <-runif((Aspecies+Plantspecies),0.02,0.02)
  h2 <- runif((Aspecies+Plantspecies),0.4,0.4)
  
  
  ## vector of species trait standard deviations
  
  N <- runif( (Aspecies+Plantspecies) , 1,1)
  nainit<- N[1:Aspecies]
  npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
  muinit<-runif((Aspecies+Plantspecies), -0.5,0.5)
  mainit<-muinit[1:Aspecies]
  mpinit<-muinit[(Aspecies+1): (Aspecies+Plantspecies)]
  
  Amatrix<-mat.comp(g,degree.animals,degree.plants)$Amatrix
  Pmatrix<-mat.comp(g,degree.animals,degree.plants)$Pmatrix
  width<-0.35
  nestedness<-nestedness_NODF(g)
  C<-Connectance(g)
  web.name<-"datasets_1/M_PL_061_01.csv"
  ba<-runif(Aspecies, 0.00,0.00)
  bp<-runif(Plantspecies,0.00,0.00)
  dganimals<-degree.animals
  dgplants<-degree.plants
  
  tmax=1000
  
  ##### the simulation with with forcing  
  index_maxdegree_animal<-which(degree.animals == max(degree.animals))
  index_maxdegree_plants<-which(degree.plants == max(degree.plants))
  degree.vector<-c(degree.animals,degree.plants)
  index_max_degree<-which(degree.vector == max(degree.vector))[1]
  species_index<-index_max_degree
  # sig[index_max_degree]<- sig[index_max_degree] + fact$variaiton_perturbation[r] 
  strength_mutualism<- 4
  mut.strength<-runif((Aspecies+Plantspecies), strength_mutualism,strength_mutualism)
  mut.strength[index_max_degree]<- strength_mutualism
    ic<-c(nainit, npinit, mainit,mpinit)
  
  params_nforcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                         Pmatrix=Pmatrix,w=width,
                         ic=ic,
                         w1=0.1,
                         individual_variation="high",
                         mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                         web.name="datasets_1/M_PL_061_01.csv",h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                         dgplants=dgplants, forcing_strength=0,
                         ua=mainit,up=mpinit,
                         species_index=species_index,
                         interaction_type="trade_off",
                         interaction="assym",
                         na=nainit,np=npinit, duration=0)
  
  #no forcing simulations
  sol1<-ode(func=eqs, y=ic, parms=params_nforcing, times=seq(0, tmax, by=1)) %>% 
    organize_results(pars = params_nforcing) 
  
  
  (netwr_no_forcing_hvar_n1<-sol1 %>% filter(type != "ma", type != "mp") %>% 
      ggplot +
      geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
      scale_y_continuous(name="population density", limits=c(0, NA)) +
      scale_color_viridis_d()+
      theme(legend.position="none")+
      annotate("text", x = 250, y = 0.5, label = "gamma[0] == 5",
               parse = TRUE))

  
  ################# forcing of 0.5 for a duration of 500 time points with assymetric interaction kernel #################
  N <- runif( (Aspecies+Plantspecies) ,1,1 )
  nainit<- N[1:Aspecies]
  npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
  mainit<-(sol1 %>% filter(type == "ma", time==tmax))$v
  mpinit<-(sol1 %>% filter(type == "mp", time==tmax))$v
  time_range<-c(0,tmax)
  deltat<- (time_range[2]-time_range[1])/1 + 1
  duration <- 500
  d<- c(rep(1,duration),rep(0,(deltat-duration)))
  duration_mat_A<-(replicate(Aspecies,d))
  duration_mat_P<-(replicate(Plantspecies,d))
  times<-seq(0, tmax, 1)
  t1_A<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_P<-as.data.frame(list(times=times, import = rep(0,length(times))))
  t1_A$import<-duration_mat_A[,1]
  t1_P$import<-duration_mat_P[,1]
  # time_func_A<-approxfun(t1_A, method="linear", rule =2)
  #time_func_P<-approxfun(t1_P, method="linear", rule =2)
  forcing_strength <- rep(0, (Aspecies+Plantspecies))
  strength_mutualism<- 1.2
  mut.strength<-runif((Aspecies+Plantspecies), strength_mutualism,strength_mutualism)
  mut.strength[index_max_degree]<- strength_mutualism
  
  ic_f<- c(nainit, npinit, mainit,mpinit)
  params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                          Pmatrix=Pmatrix,w=width,
                          ic=ic_f,
                          w1=0.1,
                          individual_variation="high",
                          mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                          web.name="datasets_1/M_PL_061_01.csv",h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                          dgplants=dgplants, forcing_strength=forcing_strength,
                          ua=mainit,up=mpinit,
                          species_index=species_index,
                          interaction_type="trade_off",
                          interaction="assym",
                          na=nainit,np=npinit, duration=duration,t1_A=t1_A,
                          t1_P=t1_P)

  sol<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
    organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs
  
  
  (f1_no_forcing<-sol %>% filter(type != "ma", type != "mp") %>% 
      ggplot +
      geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
      scale_y_continuous(name="population density", limits=c(0, NA)) +
      scale_color_viridis_d()+
      theme(legend.position="none")+
    annotate("text", x = 250, y = 0.5, label = "gamma[0] == 1.2",
               parse = TRUE)+
    annotate("text", x = 250, y = 0.45, label = "No forcing")+
    annotate("rect", xmin = 0, xmax =500, ymin = 0.44, ymax = 0.52,
               alpha = .2)+  
    annotate("text", x = 250, y = 0.47, label = "Asymmetric interaction"))
  
  
  
  ############# forcing of 0.5 for 500 time points ###################
   N <- runif( (Aspecies+Plantspecies) ,0,0.005 )
  nainit<- N[1:Aspecies]
  npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
  
  ic_f<- c(nainit, npinit, mainit,mpinit)
 
  forcing_strength <- rep(0.5, (Aspecies+Plantspecies))
  strength_mutualism<- 1.2
  mut.strength<-runif((Aspecies+Plantspecies), strength_mutualism,strength_mutualism)
  mut.strength[index_max_degree]<- strength_mutualism
  
  ic_f<- c(nainit, npinit, mainit,mpinit)
  
  params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                         Pmatrix=Pmatrix,w=width,
                         ic=ic_f,
                         w1=0.1,
                         individual_variation="high",
                         mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                         web.name="datasets_1/M_PL_061_01.csv",h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                         dgplants=dgplants, forcing_strength=forcing_strength,
                         ua=mainit,up=mpinit,
                         species_index=species_index,
                         interaction_type="trade_off",
                         interaction="assym",
                         na=nainit,np=npinit, duration=duration,t1_A=t1_A,
                         t1_P=t1_P)
  
  sol1<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
    organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs
  
  
  (f2_forcing_1<-sol1 %>% filter(type != "ma", type != "mp") %>% 
      ggplot +
      geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
      scale_y_continuous(name="population density", limits=c(0, NA)) +
      scale_color_viridis_d()+
      theme(legend.position="none")+
      annotate("rect", xmin = 0, xmax =500, ymin = 2, ymax = 0,
               alpha = .2)+
      annotate("text", x = 250, y = 1.65, label = "gamma[0] == 1.2",
               parse = TRUE)+
     annotate("text", x = 250, y = 1.58, label = "Asymmetric interaction")+
     annotate("text", x = 250, y = 1.5, label = "forcing"))
  
  
  ############# forcing of 0.5 for 500 time points with 1.4 ###################
  ic_f<- c(nainit, npinit, mainit,mpinit)
  
  forcing_strength <- rep(0.5, (Aspecies+Plantspecies))
  strength_mutualism<- 1.3
  mut.strength<-runif((Aspecies+Plantspecies), strength_mutualism,strength_mutualism)
  mut.strength[index_max_degree]<- strength_mutualism
  
  ic_f<- c(nainit, npinit, mainit,mpinit)
  
  params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                         Pmatrix=Pmatrix,w=width,
                         ic=ic_f,
                         w1=0.1,
                         individual_variation="high",
                         mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                         web.name="datasets_1/M_PL_061_01.csv",h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                         dgplants=dgplants, forcing_strength=forcing_strength,
                         ua=mainit,up=mpinit,
                         species_index=species_index,
                         interaction_type="trade_off",
                         interaction="assym",
                         na=nainit,np=npinit, duration=duration,t1_A=t1_A,
                         t1_P=t1_P)
  
  sol2<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
    organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs
  
  
  (f2_forcing_2<-sol2 %>% filter(type != "ma", type != "mp") %>% 
      ggplot +
      geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
      scale_y_continuous(name="population density", limits=c(0, NA)) +
      scale_color_viridis_d()+
      theme(legend.position="none")+
      annotate("rect", xmin = 0, xmax =500, ymin = 2, ymax = 0,
               alpha = .2)+
      annotate("text", x = 250, y = 1.65, label = "gamma[0] == 1.3",
               parse = TRUE)+
      annotate("text", x = 250, y = 1.58, label = "Asymmetric interaction")+
      annotate("text", x = 250, y = 1.5, label = "forcing"))
  
  
  
  
  ############ symmetric interaction #####################  
  ic_f<- c(nainit, npinit, mainit,mpinit)
  
  forcing_strength <- rep(0, (Aspecies+Plantspecies))
  strength_mutualism<- 1.2
  mut.strength<-runif((Aspecies+Plantspecies), strength_mutualism,strength_mutualism)
  mut.strength[index_max_degree]<- strength_mutualism
  
  ic_f<- c(nainit, npinit, mainit,mpinit)
  
  params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                         Pmatrix=Pmatrix,w=width,
                         ic=ic_f,
                         w1=0.1,
                         individual_variation="high",
                         mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                         web.name="datasets_1/M_PL_061_01.csv",h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                         dgplants=dgplants, forcing_strength=forcing_strength,
                         ua=mainit,up=mpinit,
                         species_index=species_index,
                         interaction_type="trade_off",
                         interaction="ssymetric",
                         na=nainit,np=npinit, duration=duration,t1_A=t1_A,
                         t1_P=t1_P)
  
  sol3<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
    organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs
  
  
  (f3_forcing_3<-sol3 %>% filter(type != "ma", type != "mp") %>% 
      ggplot +
      geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
      scale_y_continuous(name="population density", limits=c(0, NA)) +
      scale_color_viridis_d()+
      theme(legend.position="none")+
      annotate("text", x = 250, y = 1.65, label = "gamma[0] == 1.2",
               parse = TRUE)+
      annotate("text", x = 250, y = 1.58, label = "Gaussian interaction")+
      annotate("rect", xmin = 0, xmax =500, ymin = 1.4, ymax = 1.7,
               alpha = .2)+ 
      annotate("text", x = 250, y = 1.5, label = "No forcing"))
  
  
  
  ############ symmetric interaction with forcing involved #####################  
  ic_f<- c(nainit, npinit, mainit,mpinit)
  
  forcing_strength <- rep(0.5, (Aspecies+Plantspecies))
  strength_mutualism<- 1.2
  mut.strength<-runif((Aspecies+Plantspecies), strength_mutualism,strength_mutualism)
  mut.strength[index_max_degree]<- strength_mutualism
  
  ic_f<- c(nainit, npinit, mainit,mpinit)
  
  params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                         Pmatrix=Pmatrix,w=width,
                         ic=ic_f,
                         w1=0.1,
                         individual_variation="high",
                         mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                         web.name="datasets_1/M_PL_061_01.csv",h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                         dgplants=dgplants, forcing_strength=forcing_strength,
                         ua=mainit,up=mpinit,
                         species_index=species_index,
                         interaction_type="trade_off",
                         interaction="ssymetric",
                         na=nainit,np=npinit, duration=duration,t1_A=t1_A,
                         t1_P=t1_P)
  
  sol4<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
    organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs
  
  
  (f3_forcing_4<-sol4 %>% filter(type != "ma", type != "mp") %>% 
      ggplot +
      geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
      scale_y_continuous(name="population density", limits=c(0, NA)) +
      scale_color_viridis_d()+
      theme(legend.position="none")+
      annotate("text", x = 250, y = 1.75, label = "gamma[0] == 1.2",
               parse = TRUE)+
      annotate("rect", xmin = 0, xmax =500, ymin = 2, ymax = 0,
               alpha = .2)+
      annotate("text", x = 250, y = 1.65, label = "Gaussian interaction")+
      annotate("text", x = 250, y = 1.55, label = "forcing"))
  
  
  
  
  ############ symmetric interaction with forcing involved with strength 0.3 #####################  
  ic_f<- c(nainit, npinit, mainit,mpinit)
  
  forcing_strength <- rep(0.5, (Aspecies+Plantspecies))
  strength_mutualism<- 1.3
  mut.strength<-runif((Aspecies+Plantspecies), strength_mutualism,strength_mutualism)
  mut.strength[index_max_degree]<- strength_mutualism
  
  ic_f<- c(nainit, npinit, mainit,mpinit)
  
  params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                         Pmatrix=Pmatrix,w=width,
                         ic=ic_f,
                         w1=0.1,
                         individual_variation="high",
                         mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                         web.name="datasets_1/M_PL_061_01.csv",h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                         dgplants=dgplants, forcing_strength=forcing_strength,
                         ua=mainit,up=mpinit,
                         species_index=species_index,
                         interaction_type="trade_off",
                         interaction="ssymetric",
                         na=nainit,np=npinit, duration=duration,t1_A=t1_A,
                         t1_P=t1_P)
  
  sol5<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
    organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs
  
  
  (f3_forcing_5<-sol5 %>% filter(type != "ma", type != "mp") %>% 
      ggplot +
      geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
      scale_y_continuous(name="population density", limits=c(0, NA)) +
      scale_color_viridis_d()+
      theme(legend.position="none")+
      annotate("text", x = 250, y = 2.7, label = "gamma[0] == 1.3",
               parse = TRUE)+
      annotate("rect", xmin = 0, xmax =500, ymin = 3, ymax = 0,
               alpha = .2)+
      annotate("text", x = 250, y = 2.6, label = "Gaussian interaction")+
      annotate("text", x = 250, y = 2.5, label = "forcing"))
  
    
  
  
  ###### plotting the two networks and every population dynamics together############################################
  ################################################################################################################    
  
  
 
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
  
  
  ggpubr::ggarrange(w1, f1_no_forcing,f2_forcing_1,f2_forcing_2,
                    f3_forcing_3,f3_forcing_4,f3_forcing_5,
                    labels = c("A", "B" , "C", "D",
                               "E","F","G"),
                    nrow = 2,ncol=4)
  
  
  
  
