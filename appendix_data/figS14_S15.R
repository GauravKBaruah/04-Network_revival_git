rm(list=ls())
source("01_functions_ode.R")
library(statmod)
require(deSolve) ## for integrating ordinary differential equations
library(tidyr) ## for efficient data manipulation & plotting
library(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(viridis)
library(matlib)
library(ggplot2)

set.seed(1234)
theme_set(theme_classic()) 


load("trait_data_allnets.RData")

load("figure_4_network_data.RData")


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

# creating another dataframe where trait values of species will be kept for final simulations
trait_data_list<-list()
trait_data<-expand.grid(web=webfiles) %>% 
  as_tibble  

for (i in 1:nrow(trait_data)){
  
  g<-adj.mat(myfiles[which(myfiles == trait_data$web[i])]) #network web names
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species
  muinit <-runif((Aspecies+Plantspecies), -0.2,0.2)
  trait_data_list[i]<-list(muinit)
  trait_data$networksize[i]<-Aspecies+Plantspecies
  names(trait_data_list)[i] <- as.character(trait_data$web[i])
}



fact<- expand.grid(`Strength_mutualism`=unique(webdat_2$mut_strength), 
                   `web` =unique(webdat_2$web.name),
                   `forcing_strength`= 0.5,
                   `forcing_duration`=500,
                   `no_species_forced`=1,
                   `model`="abundance",
                   `interaction_type`= "trade_off", #no_trade_off
                   `individual.variation` = c("high","low"),
                   `random_seed`=4327+(1:1)*100) %>%
  as_tibble  


new_ddf<-NULL
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
  
  ##control loop for selecting whether variation is high or low
  if(fact$individual.variation[r] == "low"){
    sig <-runif((Aspecies+Plantspecies),0.002,0.002)}else if(fact$individual.variation[r] == "high"){
      sig <-runif((Aspecies+Plantspecies),0.02,0.02)}
  
  
  h2 <- runif((Aspecies+Plantspecies),0.4,0.4)
  
  
  ## vector of species trait standard deviations
  
  N <- runif( (Aspecies+Plantspecies) , 0.000,0.005)
  nainit<- N[1:Aspecies]
  npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
  index<-which(names(trait_data_list)==fact$web[r])
  muinit<-(outt %>% filter(webname == fact$web[r]))$m
  mainit<-muinit[1:Aspecies]
  mpinit<-muinit[(Aspecies+1): (Aspecies+Plantspecies)]
  
  Amatrix<-mat.comp(g,degree.animals,degree.plants)$Amatrix
  Pmatrix<-mat.comp(g,degree.animals,degree.plants)$Pmatrix
  width<-0.35
  nestedness<-nestedness_NODF(g)
  C<-Connectance(g)
  web.name<-fact$web[r]
  ba<-runif(Aspecies, 0.00,0.00)
  bp<-runif(Plantspecies,0.00,0.00)
  dganimals<-degree.animals
  dgplants<-degree.plants
  
  tmax=900
  
  ##### the simulation with with forcing  
  index_maxdegree_animal<-which(degree.animals == max(degree.animals))
  index_maxdegree_plants<-which(degree.plants == max(degree.plants))
  degree.vector<-c(degree.animals,degree.plants)
  index_max_degree<-which(degree.vector == max(degree.vector))[1]
  species_index<-index_max_degree
  # sig[index_max_degree]<- sig[index_max_degree] + fact$variaiton_perturbation[r] 
  
  mut.strength<-runif((Aspecies+Plantspecies), fact$Strength_mutualism[r],fact$Strength_mutualism[r])
  mut.strength[index_max_degree]<- fact$Strength_mutualism[r] 
  time_range<-c(0,tmax)
  deltat<- (time_range[2]-time_range[1])/1 + 1
  duration <- fact$forcing_duration[r]
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
  forcing_strength <- rep(fact$forcing_strength[r], (Aspecies+Plantspecies))
  ic_f<-c(nainit, npinit, mainit,mpinit)
  
  params_forcing <- list(time=tmax,matrix=g,sig=sig,Amatrix=Amatrix,
                         Pmatrix=Pmatrix,w=width,
                         ic=ic_f,
                         dat=dat,
                         individual_variation=fact$individual.variation[r],
                         mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                         web.name=web.name,h2=0.4, ba=ba,bp=bp,dganimals=dganimals,
                         dgplants=dgplants, forcing_strength=forcing_strength,
                         ua=mainit,up=mpinit,
                         species_index=species_index,
                         interaction_type=fact$interaction_type[r],
                         na=nainit,np=npinit, duration=duration,t1_A=t1_A,
                         t1_P=t1_P)
  
  
  params_indirect<-list(mut_strength=mut.strength[1],w=0.35, comp_matrixA = Amatrix,comp_matrixP =Pmatrix,trait_vals=c(mainit,mpinit),Na=nainit,Np=npinit, g=g,
        A=Aspecies,P=Plantspecies, s=sig)
  
  
  fact$sum_indirect_effects[r]<- indirect_effects(parameters = params_indirect)$sum_indirect_effects
  fact$network_size[r] <- Aspecies+Plantspecies
  fact$nestedness[r] <- nestedness_NODF(g)
  fact$connectance[r] <- Connectance(g)
  fact$mean_indirect_effects[r] <- indirect_effects(parameters = params_indirect)$mean_indirect_effects
  fact$sum_indirect_effects[r]<-indirect_effects(parameters = params_indirect)$sum_indirect_effects
  fact$spectral_radius[r]<-indirect_effects(parameters = params_indirect)$spectral_radius
  fact$sum_indirect_effects3[r]<-indirect_effects(parameters = params_indirect)$indirect_effects_order3
  fact$sum_indirect_effects4[r]<-indirect_effects(parameters = params_indirect)$indirect_effects_order4
  fact$mean_indirect_effects3[r]<-indirect_effects(parameters = params_indirect)$mean_indirect_effects3
  fact$mean_indirect_effects4[r]<-indirect_effects(parameters = params_indirect)$mean_indirect_effects4
  fact$mean_indirect_effects5[r]<-indirect_effects(parameters = params_indirect)$mean_indirect_effects5
  fact$sum_indirect_effects5[r]<-indirect_effects(parameters = params_indirect)$indirect_effects_order5
  

  #dat<-cluster_run_func(params_forcing = params_forcing,ic_f = ic_f,tmax = tmax)
  
  print(r)
  
}

(a1<-fact %>% filter(spectral_radius <=1) %>%  ggplot(aes(x=Strength_mutualism, y=mean_indirect_effects3,
                    color=factor(individual.variation )))+
  geom_point(size=3,alpha=0.5)+
  scale_color_manual(values=c("#000000", "#E69F00"))+
  theme_classic()+
    ylab(expression(paste(beta^3, " ", ",", "Mean indirect effects of order 3")))+
  xlab(expression(paste(gamma[0], " ", ", ", "Mutualism strength")))+
  labs(color="Individual variation")+
  geom_hline(yintercept = 0,linetype = "dashed"))
  
(a2<-fact %>% filter(spectral_radius <=1) %>% ggplot(aes(x=Strength_mutualism, y=mean_indirect_effects,
                    color=factor(individual.variation )))+
  geom_point(size=3,alpha=0.5)+
  ylim(c(-100,200))+
  scale_color_manual(values=c("#000000", "#E69F00"))+
  theme_classic()+
  ylab(expression(paste(beta, " ", ",", "Net mean indirect effects")))+
  xlab(expression(paste(gamma[0], " ", ", ", "Mutualism strength")))+
  labs(color="Individual variation")+
  geom_hline(yintercept = 0,linetype = "dashed"))


ggpubr::ggarrange(a1,a2, nrow=1,ncol=2,labels = c("A","B"))


 
nrow(webdat_2)
new_data<-webdat_2 
for(i in 1:nrow(new_data)){
print(i)  
  new_dat <- fact %>% filter(  web == webdat_2$web.name[i], Strength_mutualism == webdat_2$mut_strength[i],
                  individual.variation == webdat_2$individual_variation[i])
  new_data$sum_indirect_effects[i] <-  new_dat$sum_indirect_effects
  new_data$mean_indirect_effects[i] <-  new_dat$mean_indirect_effects
  new_data$spectral_radius[i] <-  new_dat$spectral_radius
  new_data$mean_indirect_effects_1[i] <-  new_dat$mean_indirect_effects_1
  new_data$mean_indirect_effects_1[i] <-  new_dat$mean_indirect_effects_1
  new_data$mean_indirect_effects_1[i] <-  new_dat$mean_indirect_effects_1
  new_data$mean_indirect_effects4[i]<- new_dat$mean_indirect_effects4
  new_data$mean_indirect_effects3[i]<-new_dat$sum_indirect_effects3
  new_data$sum_indirect_effects3[i]<-new_dat$sum_indirect_effects3
  new_data$sum_indirect_effects4[i]<-new_dat$sum_indirect_effects4
  
  
}


appender <- function(string) 
  latex2exp::TeX(paste("$\\gamma_{\\0} = $", string))  

new_data$mut_strength <-as.factor(new_data$mut_strength)
n1<-(new_data %>%filter(spectral_radius <=1) %>% 
       ggplot(aes(y=((mean_indirect_effects)), 
                  x = log(recovery_biomass) , 
                  color=mut_strength))+
       geom_hline(yintercept = 0,linetype ="dashed")+
       geom_point(alpha = 0.5, size = 2)+
       ylim(c(-50,10))+
       scale_color_viridis_d()+
       ylab(expression(paste(beta,","," ", "Mean net indirect effects")))+
       theme_classic()+
         labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
       #labs(color=("Individual variation"))+
       # theme(legend.position="none")+
       xlab("log(Recovery biomass)"))+
  stat_smooth( method = "lm", size=1,formula = y ~ x, se = FALSE)+
  facet_wrap(.~individual_variation)


n1

n2<-(new_data %>%filter(spectral_radius <=1) %>% 
       ggplot(aes(y=((mean_indirect_effects3)), 
                  x = log(recovery_biomass) , 
                  color=mut_strength))+
       geom_point(alpha = 0.5, size = 2)+
       ylab(expression(paste(beta^3,","," ", "Mean indirect effects order 3")))+
       theme_classic()+
      scale_color_viridis_d()+
         labs(color=expression(paste(gamma[0], "," ,"mutualism strength")))+
       #labs(color=("Individual variation"))+
      # scale_color_manual(values=c("#000000", "#E69F00"))+
       geom_hline(yintercept = 0,linetype ="dashed")+
       # theme(legend.position="none")+
       xlab("log(Recovery biomass)"))+
  stat_smooth( method = "lm", formula = y ~ x, se = FALSE)+
  facet_wrap(.~individual_variation)

ggpubr::ggarrange(n1,n2, nrow=2,ncol=1, labels=c("A","B"))
n2
n3<-(new_data %>%filter(mut_strength %in% c(1.0, 1.1, 1.2, 1.3, 1.5)) %>% 
       ggplot(aes(y=((mean_indirect_effects4)), 
                  x = log(recovery_biomass) , 
                  color=individual_variation))+
       geom_point(alpha = 0.5, size = 2)+
       ylab(expression(paste(beta^4,","," ", "Mean indirect effects order 4")))+
       theme_classic()+
       #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
       labs(color=("Individual variation"))+
       scale_color_manual(values=c("#000000", "#E69F00"))+
       # theme(legend.position="none")+
       xlab("log(Recovery biomass)"))+
  stat_smooth( method = "lm", formula = y ~ x, se = FALSE)+
  facet_wrap(.~mut_strength)
n3


ggpubr::ggarrange(n1,n2,n3,nrow = 2,ncol=2)
