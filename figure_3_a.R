rm(list=ls())
source("01_ODE_Function.R")
library(statmod)
require(deSolve) ## for integrating ordinary differential equations
#require(tidyverse) ## for efficient data manipulation & plotting
#library(cowplot) ## for arranging plots in a grid
library(dplyr)
#library(readr)
library(beepr)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(ggnet)
library(gganimate)
library(gifski)

mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:153]
#make a empty dataframe of all the mutualism networks
fact_2<- expand.grid(`web` =newfiles[1:153],
                     #`initial.trait`=ru
                     `model`="abundance",
                     `interaction_type`= "trade_off", #no_trade_off
                     `random_seed`=4327+1*100) %>%
  as_tibble %>%
  mutate(`Biomass.animal`=0,
         `Biomass.plant` =0,
         `Nestedness`= 0,
         `Connectance` = 0,
         `network_size`=0)
model.t<-list()


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


fact_2 <- fact_2 %>% filter(network_size < 270 )

webfiles<-as.character(fact_2$web)





#	datasets_1/M_PL_006.csv, nestedness =0.5, network size = 78
# datasets_1/M_PL_061_17.cs

fact<- expand.grid(`Strength_mutualism`=1.15, 
                   `web` = "datasets_1/M_PL_037.csv", #"datasets_1/M_PL_045.csv",
                   mut_strength_perturbation=0, 
                   `forcing_strength`= 0.5,
                   `forcing_duration`=500,
                   `no_species_forced`=1,
                   `model`="abundance",
                   `interaction_type`= "trade_off", 
                   `individual.variation` ="high",
                   `random_seed`=4327+(1:1)*100) %>%
  as_tibble()

new_ddf<-NULL
load("Mean_trait_data.RData")
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
sig <-runif((Aspecies+Plantspecies),0.02,0.02)
h2 <- runif((Aspecies+Plantspecies),0.4,0.4)


## vector of species trait standard deviations
N <- runif( (Aspecies+Plantspecies) , 0,0.005) #mimicking a system with low population density near the collapse state
nainit<-N[1:Aspecies]
npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]

index<-which(names(trait_data_list)==fact$web[1])

muinit<-trait_data_list[[index]]
mainit<-muinit[1:Aspecies]
mpinit<-muinit[(Aspecies+1): (Aspecies+Plantspecies)]
#mainit<- runif(Aspecies,-0.25,0.25)
#mpinit<-runif(Plantspecies, -0.25,0.25)
index_maxdegree_animal<-which(degree.animals == max(degree.animals))
index_maxdegree_plants<-which(degree.plants == max(degree.plants))

degree.vector<-c(degree.animals,degree.plants)
index_max_degree<-which(degree.vector == max(degree.vector))

species_index<-index_max_degree
Amatrix<-mat.comp(g,degree.animals,degree.plants)$Amatrix
Pmatrix<-mat.comp(g,degree.animals,degree.plants)$Pmatrix
gamma=0.35#fact_lessvar$Strength_mutualism[r]
tmax<-1200
nestedness<-nestedness_NODF(g)
C<-Connectance(g)
web.name<-fact$web[1]
ba<-runif(Aspecies, 0,0)
bp<-runif(Plantspecies,0,0)
dganimals<-degree.animals
dgplants<-degree.plants

mut.strength<-runif((Aspecies+Plantspecies), fact$Strength_mutualism[1],fact$Strength_mutualism[1])
mut.strength[index_max_degree]<- fact$Strength_mutualism[1] 
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
forcing_strength <- rep(fact$forcing_strength[1], (Aspecies+Plantspecies))
ic_f<-c(nainit, npinit, mainit,mpinit)

params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                       Pmatrix=Pmatrix,w=gamma,
                       ic=ic_f,
                       dat=dat,
                       individual_variation=fact$individual.variation[1],
                       mut.strength=mut.strength,m=mainit,C=C,nestedness=nestedness,
                       web.name=web.name,h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                       dgplants=dgplants, forcing_strength=forcing_strength,
                       ua=mainit,up=mpinit,
                       species_index=species_index,
                       interaction_type=fact$interaction_type[1],
                       na=nainit,np=npinit, duration=duration,
                       t1_A=t1_A,t1_P=t1_P)


sol1<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
  organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs



(netwr_forcing_hvar_n1<-sol1 %>% filter(type != "ma", type != "mp") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme(legend.position="none")+
    annotate("text", x = 350, y = 2.35, label = "gamma[0] == 1.15",
             parse = TRUE)+
    annotate("text", x = 750, y = 2.359, label = ", with forcing",
             parse = FALSE)+
    annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2.5,
             alpha = .2))


#animate(netwr_forcing_lvar_n1, height = 500, width = 800, fps = 30, duration = 15,
#        end_pause = 60, res = 100, renderer = gifski_renderer())
#anim_save("revival_net_low.gif")



########## time =50 ###############
newg_recovery1<-network_structure(Na=(sol1 %>% filter(time==50,type =="N"))$v ,
                                  Np=(sol1 %>% filter(time==50,type =="P"))$v , g=g )
net1 = network(newg_recovery1, directed = FALSE)


#net = network(model.t$``, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net1)
net1 %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net1 %v% "color" = ifelse(net1 %v% "groups" == "plants", "#0072B2", "#E69F00" )
#ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)

deg<-c(degree.plants,degree.animals)
w_h_t_50<-ggnet2(net1, mode="circle", size=deg, 
                   edge.size = 1.1,max_size =12, 
                   color ="color",edge.alpha = 1, legend.position = "")




############ time = 300 ##############

newg_recovery_t300<-network_structure(Na=(sol1 %>% filter(time==300,type =="N"))$v ,
                                      Np=(sol1 %>% filter(time==300,type =="P"))$v , g=g )
net2 = network(newg_recovery_t300, directed = FALSE)


#net = network(model.t$``, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net2)
net2 %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net2 %v% "color" = ifelse(net2 %v% "groups" == "plants", "#0072B2", "#E69F00" )
#ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)

deg<-c(degree.plants,degree.animals)
w_h_t_300<-ggnet2(net2, mode="circle", size=deg, 
                    edge.size = 1.1,max_size =15, 
                    color ="color",edge.alpha = 1, legend.position = "")



########### time = 500 ##########
newg_recovery_t1000<-network_structure(Na=(sol1 %>% filter(time==1000,type =="N"))$v ,
                                       Np=(sol1 %>% filter(time==1000,type =="P"))$v , g=g )
net3 = network(newg_recovery_t1000, directed = FALSE)


#net = network(model.t$``, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net3)
net3 %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net3 %v% "color" = ifelse(net3 %v% "groups" == "plants", "#0072B2", "#E69F00" )
#ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)

deg<-c(degree.plants,degree.animals)
(w_h_t_1000<-ggnet2(net3, mode="circle", size=deg, 
                      edge.size = 1.1,max_size =12, 
                      color ="color",edge.alpha = 1, legend.position = ""))













###################### LOW INTRASpecific Variation and resurrection of a network##############################
############################################################################################################



sig <-runif((Aspecies+Plantspecies),0.005,0.005)
N <- runif( (Aspecies+Plantspecies) , 0,0.005) #mimicking a system with low population density near the collapse state
nainit<-N[1:Aspecies]
npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]


#mainit<- runif(Aspecies, -0.25, 0.25)
#mpinit<-runif(Plantspecies, -0.25, 0.25)
index_maxdegree_animal<-which(degree.animals == max(degree.animals))
index_maxdegree_plants<-which(degree.plants == max(degree.plants))

degree.vector<-c(degree.animals,degree.plants)
index_max_degree<-which(degree.vector == max(degree.vector))

species_index<-index_max_degree

mut.strength<-runif((Aspecies+Plantspecies), fact$Strength_mutualism[1],fact$Strength_mutualism[1])
mut.strength[index_max_degree]<- fact$Strength_mutualism[1] 
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
time_func_A<-approxfun(t1_A, method="linear", rule =2)
time_func_P<-approxfun(t1_P, method="linear", rule =2)
forcing_strength <- rep(fact$forcing_strength[1], (Aspecies+Plantspecies))
ic_f<-c(nainit, npinit, mainit,mpinit)

params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                       Pmatrix=Pmatrix,w=gamma,
                       ic=ic_f,
                       dat=dat,
                       individual_variation="low",
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

(netwr_forcing_lvar_n2<-sol %>% filter(type != "ma", type != "mp") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme(legend.position="none")+
    annotate("text", x = 350, y = 1.5, label = "gamma[0] == 1.15",
             parse = TRUE)+
    annotate("text", x = 750, y = 1.509, label = ", with forcing",
             parse = FALSE)+
    annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2.5,
             alpha = .2))



########## time =50 ###############
newg_recovery<-network_structure(Na=(sol %>% filter(time==50,type =="N"))$v ,Np=(sol %>% filter(time==50,type =="P"))$v , g=g )
net = network(newg_recovery, directed = FALSE)


#net = network(model.t$``, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net)
net %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
#ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)

deg<-c(degree.plants,degree.animals)
w_low_t_50<-ggnet2(net, mode="circle", size=deg, 
           edge.size = 1.1,max_size =12, 
           color ="color",edge.alpha = 1, legend.position = "")




############ time = 300 ##############

newg_recovery_t300<-network_structure(Na=(sol %>% filter(time==300,type =="N"))$v ,
                                 Np=(sol %>% filter(time==300,type =="P"))$v , g=g )
net2 = network(newg_recovery_t300, directed = FALSE)


#net = network(model.t$``, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net2)
net2 %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net2 %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
#ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)

deg<-c(degree.plants,degree.animals)
w_low_t_300<-ggnet2(net2, mode="circle", size=deg, 
                   edge.size = 1.1,max_size =12, 
                   color ="color",edge.alpha = 1, legend.position = "")



########### time = 500 ##########
newg_recovery_t1000<-network_structure(Na=(sol %>% filter(time==750,type =="N"))$v ,
                                      Np=(sol %>% filter(time==750,type =="P"))$v , g=g )
net3 = network(newg_recovery_t1000, directed = FALSE)


#net = network(model.t$``, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net3)
net3 %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net3 %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
#ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)

deg<-c(degree.plants,degree.animals)
(w_low_t_1000<-ggnet2(net3, mode="circle", size=deg, 
                    edge.size = 1.1,max_size =12, 
                    color ="color",edge.alpha = 1, legend.position = ""))




###################### PLOTTING EVERYTHING #####################################


ggpubr::ggarrange(w_h_t_50,w_h_t_300,w_h_t_1000,netwr_forcing_hvar_n1,
                  w_low_t_50,w_low_t_300,w_low_t_1000,netwr_forcing_lvar_n2,
                  ncol=4,nrow=2,
                  labels = c("A. t = 50", "B. t = 300", "C. t = 1000",
                             "D. High variation",
                             "E.", "F.", "G.", "H. Low variation"))


