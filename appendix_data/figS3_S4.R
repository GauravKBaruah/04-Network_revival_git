rm(list=ls())
source("01_ODE_Function.R")
library(statmod)
require(deSolve) ## for integrating ordinary differential equations
library(dplyr)
library(beepr)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(gganimate)

set.seed(1234)
mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:153]
load("trait_data_allnets.RData")

fact<- expand.grid(`Strength_mutualism`=1.15, 
                   `web` = "datasets_1/M_PL_045.csv",
                   mut_strength_perturbation=0, 
                   `forcing_strength`= 0,
                   `forcing_duration`=500,
                   `no_species_forced`=1,
                   `model`="abundance",
                   `interaction_type`= "trade_off", 
                   `individual.variation` ="high",
                   `random_seed`=4327+(1:1)*100) %>%
  as_tibble()



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
sig <-runif((Aspecies+Plantspecies),0.005,0.005)
h2 <- runif((Aspecies+Plantspecies),0,0)


## vector of species trait standard deviations
N <- runif( (Aspecies+Plantspecies) , 0,0.001) #mimicking a system with low population density near the collapse state
nainit<-N[1:Aspecies]
npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]

ma<-(outt %>% filter(webname == fact$web[1]))$m
mainit<-ma[1:Aspecies]
mpinit<-ma[(Aspecies+1): (Aspecies+Plantspecies)]
index_maxdegree_animal<-which(degree.animals == max(degree.animals))
index_maxdegree_plants<-which(degree.plants == max(degree.plants))

degree.vector<-c(degree.animals,degree.plants)
index_max_degree<-which(degree.vector == max(degree.vector))

species_index<-index_max_degree
Amatrix<-mat.comp(g,degree.animals,degree.plants)$Amatrix
Pmatrix<-mat.comp(g,degree.animals,degree.plants)$Pmatrix
gamma=0.35#fact_lessvar$Strength_mutualism[r]
tmax<-1500
nestedness<-nestedness_NODF(g)
C<-Connectance(g)
web.name<-fact$web[1]
ba<-runif(Aspecies, 0.,0.)
bp<-runif(Plantspecies,0.,0.)
dganimals<-degree.animals
dgplants<-degree.plants
fact$Strength_mutualism[1]<-1.15
mut.strength<-runif((Aspecies+Plantspecies), fact$Strength_mutualism[1],fact$Strength_mutualism[1])
mut.strength[index_max_degree]<- fact$Strength_mutualism[1] 
time_range<-c(0,tmax)
deltat<- (time_range[2]-time_range[1])/1 + 1
duration <- fact$forcing_duration[1]<- 500
d<- c(rep(1,duration),rep(0,(deltat-duration)))
duration_mat_A<-(replicate(Aspecies,d))
duration_mat_P<-(replicate(Plantspecies,d))
times<-seq(0, tmax, 1)
t1_A<-as.data.frame(list(times=times, import = rep(0,length(times))))
t1_P<-as.data.frame(list(times=times, import = rep(0,length(times))))
t1_A$import<-duration_mat_A[,1]
t1_P$import<-duration_mat_P[,1]
fact$forcing_strength[1]<-0.5
forcing_strength <- rep(fact$forcing_strength[1], (Aspecies+Plantspecies))
ic_f<-c(nainit, npinit, mainit,mpinit)

params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                       Pmatrix=Pmatrix,w=gamma,
                       ic=ic_f,
                       individual_variation=fact$individual.variation[1],
                       mut.strength=mut.strength,m=mainit,C=C,nestedness=nestedness,
                       web.name=web.name,h2=0.0, ba=ba,bp=bp,dganimals=dganimals,
                       dgplants=dgplants, forcing_strength=forcing_strength,
                       ua=mainit,up=mpinit,
                       species_index=species_index,
                       interaction_type=fact$interaction_type[1],
                       na=nainit,np=npinit, duration=duration,
                       t1_A=t1_A,t1_P=t1_P)


sol1<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
  organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs



(netwr_forcing_lvar<-sol1 %>% filter(type != "ma", type != "mp") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    theme(legend.position="none")+
     annotate("text", x = 250, y = 1.5, label = "gamma[0] == 1.15 + forcing",
      parse = TRUE)+
    annotate("text", x = 1300, y = 2, label = "evolution == No",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 1.9, label = "sigma[i] == 0.005",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2,
     alpha = .2))


(netwr_forcing_lvar_m<-sol1 %>% filter(type != "N", type != "P") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    ylab("mean trait values")+
    ylim(c(-0.5,0.5))+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 0.45, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 0.5, label = "evolution == No",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 0.45, label = "sigma[i] == 0.005",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = -0.5, ymax = 0.5,
             alpha = .2))



##########high variation forcing- no evolutiona ###############


sig <-runif((Aspecies+Plantspecies),0.02,0.02)
params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                       Pmatrix=Pmatrix,w=gamma,
                       ic=ic_f,
                       individual_variation=fact$individual.variation[1],
                       mut.strength=mut.strength,m=mainit,C=C,nestedness=nestedness,
                       web.name=web.name,h2=0, ba=ba,bp=bp,dganimals=dganimals,
                       dgplants=dgplants, forcing_strength=forcing_strength,
                       ua=mainit,up=mpinit,
                       species_index=species_index,
                       interaction_type=fact$interaction_type[1],
                       na=nainit,np=npinit, duration=duration,
                       t1_A=t1_A,t1_P=t1_P)


solhvar<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
  organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs

(netwr_forcing_hvar<-solhvar %>% filter(type != "ma", type != "mp") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 1.5, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 2, label = "evolution == No",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 1.9, label = "sigma[i] == 0.02",
             parse = TRUE,size=4, color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2,
             alpha = .2))

(netwr_forcing_hvar_m<-solhvar %>% filter(type != "N", type != "P") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    ylim(c(-0.5,0.5))+
    ylab("mean trait values")+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 0.45, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 0.5, label = "evolution == No",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 0.45, label = "sigma[i] == 0.02",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = -0.5, ymax = 0.5,
             alpha = .2))
##########high variation forcing- yes evolution ###############


sig <-runif((Aspecies+Plantspecies),0.02,0.02)
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


solhvarevo<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
  organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs



(netwr_forcing_hvar_ev0<-solhvarevo %>% filter(type != "ma", type != "mp") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 1.5, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 2, label = "evolution == Yes",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 1.9, label = "sigma[i] == 0.02",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2,
             alpha = .2))

(netwr_forcing_hvar_ev0_m<-solhvarevo %>% filter(type != "N", type != "P") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    ylim(c(-0.5,0.5))+
    ylab("mean trait values")+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 0.45, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 0.5, label = "evolution == Yes",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 0.45, label = "sigma[i] == 0.02",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = -0.5, ymax = 0.5,
             alpha = .2))
########## evolution + low variation ########


sig <-runif((Aspecies+Plantspecies),0.005,0.005)
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


sollvarevo<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
  organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs



(netwr_forcing_lvar_ev0<-sollvarevo %>% filter(type != "ma", type != "mp") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 1.5, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 2, label = "evolution == Yes",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 1.9, label = "sigma[i] == 0.005",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2,
             alpha = .2))

(netwr_forcing_lvar_ev0_m<-sollvarevo %>% filter(type != "N", type != "P") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    ylim(c(-0.5,0.5))+
    ylab("mean trait values")+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 0.45, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 0.5, label = "evolution == Yes",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 0.45, label = "sigma[i] == 0.005",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = -0.5, ymax = 0.5,
             alpha = .2))



########## time =50 ###############
#newg_recovery<-network_structure(Na=(sol %>% filter(time==50,type =="N"))$v ,Np=(sol %>% filter(time==50,type =="P"))$v , g=g )
net = network(g, directed = FALSE)


#net = network(model.t$``, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net)
net %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
#ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)

deg<-c(degree.plants,degree.animals)
webg<-ggnet2(net, mode="circle", size=deg, 
                   edge.size = 1.1,max_size =12, 
                   color ="color",edge.alpha = 1, legend.position = "")



###################### PLOTTING EVERYTHING #####################################


ggpubr::ggarrange(webg,netwr_forcing_lvar,netwr_forcing_lvar_m,
                  netwr_forcing_hvar,
                  netwr_forcing_hvar_m,netwr_forcing_lvar_ev0,netwr_forcing_lvar_ev0_m,
                  netwr_forcing_hvar_ev0,netwr_forcing_hvar_ev0_m,
                  ncol=3,nrow=3,
                  labels = c("A", "B", "C",
                             "D",
                             "E", "F","G","H"))



## NETWORK 2###




fact<- expand.grid(`Strength_mutualism`=1.15, 
                     `web` = "datasets_1/M_PL_060_19.csv",
                   mut_strength_perturbation=0, 
                   `forcing_strength`= 0,
                   `forcing_duration`=400,
                   `no_species_forced`=1,
                   `model`="abundance",
                   `interaction_type`= "trade_off", 
                   `individual.variation` ="high",
                   `random_seed`=4327+(1:1)*100) %>%
  as_tibble()

new_ddf<-NULL
load("Upper_bound_abundance_response_time.RData")



g<-adj.mat(myfiles[which(myfiles == fact$web[1])]) #network web names
dim(g)
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
sig <-runif((Aspecies+Plantspecies),0.005,0.005)
h2 <- runif((Aspecies+Plantspecies),0,0)


## vector of species trait standard deviations
N <- runif( (Aspecies+Plantspecies) , 0,0.001) #mimicking a system with low population density near the collapse state
nainit<-N[1:Aspecies]
npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]


mainit<- runif(Aspecies,-0.25,0.25)
mpinit<-runif(Plantspecies, -0.25,-0.25)
index_maxdegree_animal<-which(degree.animals == max(degree.animals))
index_maxdegree_plants<-which(degree.plants == max(degree.plants))

degree.vector<-c(degree.animals,degree.plants)
index_max_degree<-which(degree.vector == max(degree.vector))

species_index<-index_max_degree
Amatrix<-mat.comp(g,degree.animals,degree.plants)$Amatrix
Pmatrix<-mat.comp(g,degree.animals,degree.plants)$Pmatrix
gamma=0.35#fact_lessvar$Strength_mutualism[r]
tmax<-1500
nestedness<-nestedness_NODF(g)
C<-Connectance(g)
web.name<-fact$web[1]
ba<-runif(Aspecies, 0.,0.)
bp<-runif(Plantspecies,0.,0.)
dganimals<-degree.animals
dgplants<-degree.plants
fact$Strength_mutualism[1]<-1.15
mut.strength<-runif((Aspecies+Plantspecies), fact$Strength_mutualism[1],fact$Strength_mutualism[1])
mut.strength[index_max_degree]<- fact$Strength_mutualism[1] 
time_range<-c(0,tmax)
deltat<- (time_range[2]-time_range[1])/1 + 1
duration <- fact$forcing_duration[1]<- 500
d<- c(rep(1,duration),rep(0,(deltat-duration)))
duration_mat_A<-(replicate(Aspecies,d))
duration_mat_P<-(replicate(Plantspecies,d))
times<-seq(0, tmax, 1)
t1_A<-as.data.frame(list(times=times, import = rep(0,length(times))))
t1_P<-as.data.frame(list(times=times, import = rep(0,length(times))))
t1_A$import<-duration_mat_A[,1]
t1_P$import<-duration_mat_P[,1]
fact$forcing_strength[1]<-0.5
forcing_strength <- rep(fact$forcing_strength[1], (Aspecies+Plantspecies))
ic_f<-c(nainit, npinit, mainit,mpinit)

params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                       Pmatrix=Pmatrix,w=gamma,
                       ic=ic_f,
                       individual_variation=fact$individual.variation[1],
                       mut.strength=mut.strength,m=mainit,C=C,nestedness=nestedness,
                       web.name=web.name,h2=0.0, ba=ba,bp=bp,dganimals=dganimals,
                       dgplants=dgplants, forcing_strength=forcing_strength,
                       ua=mainit,up=mpinit,
                       species_index=species_index,
                       interaction_type=fact$interaction_type[1],
                       na=nainit,np=npinit, duration=duration,
                       t1_A=t1_A,t1_P=t1_P)


sol1<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
  organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs



(netwr_forcing_lvar<-sol1 %>% filter(type != "ma", type != "mp") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 1.5, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 2, label = "evolution == No",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 1.9, label = "sigma[i] == 0.005",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2,
             alpha = .2))


(netwr_forcing_lvar_m<-sol1 %>% filter(type != "N", type != "P") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    ylab("mean trait values")+
    ylim(c(-0.5,0.5))+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 0.45, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 0.5, label = "evolution == No",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 0.45, label = "sigma[i] == 0.005",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = -0.5, ymax = 0.5,
             alpha = .2))



##########high variation forcing- no evolutiona ###############


sig <-runif((Aspecies+Plantspecies),0.02,0.02)
params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                       Pmatrix=Pmatrix,w=gamma,
                       ic=ic_f,
                       individual_variation=fact$individual.variation[1],
                       mut.strength=mut.strength,m=mainit,C=C,nestedness=nestedness,
                       web.name=web.name,h2=0, ba=ba,bp=bp,dganimals=dganimals,
                       dgplants=dgplants, forcing_strength=forcing_strength,
                       ua=mainit,up=mpinit,
                       species_index=species_index,
                       interaction_type=fact$interaction_type[1],
                       na=nainit,np=npinit, duration=duration,
                       t1_A=t1_A,t1_P=t1_P)


solhvar<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
  organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs

(netwr_forcing_hvar<-solhvar %>% filter(type != "ma", type != "mp") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 1.5, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 2, label = "evolution == No",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 1.9, label = "sigma[i] == 0.02",
             parse = TRUE,size=4, color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2,
             alpha = .2))

(netwr_forcing_hvar_m<-solhvar %>% filter(type != "N", type != "P") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    ylim(c(-0.5,0.5))+
    ylab("mean trait values")+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 0.45, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 0.5, label = "evolution == No",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 0.45, label = "sigma[i] == 0.02",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = -0.5, ymax = 0.5,
             alpha = .2))
##########high variation forcing- yes evolution ###############


sig <-runif((Aspecies+Plantspecies),0.02,0.02)
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


solhvarevo<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
  organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs



(netwr_forcing_hvar_ev0<-solhvarevo %>% filter(type != "ma", type != "mp") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 1.5, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 2, label = "evolution == Yes",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 1.9, label = "sigma[i] == 0.02",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2,
             alpha = .2))

(netwr_forcing_hvar_ev0_m<-solhvarevo %>% filter(type != "N", type != "P") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    ylim(c(-0.5,0.5))+
    ylab("mean trait values")+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 0.45, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 0.5, label = "evolution == Yes",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 0.45, label = "sigma[i] == 0.02",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = -0.5, ymax = 0.5,
             alpha = .2))
########## evolution + low variation ########


sig <-runif((Aspecies+Plantspecies),0.005,0.005)
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


sollvarevo<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
  organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs



(netwr_forcing_lvar_ev0<-sollvarevo %>% filter(type != "ma", type != "mp") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 1.5, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 2, label = "evolution == Yes",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 1.9, label = "sigma[i] == 0.005",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = 0, ymax = 2,
             alpha = .2))

(netwr_forcing_lvar_ev0_m<-sollvarevo %>% filter(type != "N", type != "P") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=2) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme_classic()+
    ylim(c(-0.5,0.5))+
    ylab("mean trait values")+
    theme(legend.position="none")+
    annotate("text", x = 250, y = 0.45, label = "gamma[0] == 1.15 + forcing",
             parse = TRUE)+
    annotate("text", x = 1300, y = 0.5, label = "evolution == Yes",
             parse = TRUE,size=4,color="darkgreen")+
    annotate("text", x = 1300, y = 0.45, label = "sigma[i] == 0.005",
             parse = TRUE,size=4 , color="darkgreen")+
    annotate("rect", xmin = 0, xmax = 500, ymin = -0.5, ymax = 0.5,
             alpha = .2))



########## time =50 ###############
#newg_recovery<-network_structure(Na=(sol %>% filter(time==50,type =="N"))$v ,Np=(sol %>% filter(time==50,type =="P"))$v , g=g )
net = network(g, directed = FALSE)


#net = network(model.t$``, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net)
net %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net %v% "color" = ifelse(net %v% "groups" == "plants", "#0072B2", "#E69F00" )
#ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)

deg<-c(degree.plants,degree.animals)
webg<-ggnet2(net, mode="circle", size=deg, 
             edge.size = 1.1,max_size =12, 
             color ="color",edge.alpha = 1, legend.position = "")



ggpubr::ggarrange(webg,netwr_forcing_lvar,netwr_forcing_lvar_m,
                  netwr_forcing_hvar,
                  netwr_forcing_hvar_m,netwr_forcing_lvar_ev0,netwr_forcing_lvar_ev0_m,
                  netwr_forcing_hvar_ev0,netwr_forcing_hvar_ev0_m,
                  ncol=3,nrow=3,
                  labels = c("A", "B", "C",
                             "D",
                             "E", "F","G","H"))



