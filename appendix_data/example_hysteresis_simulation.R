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
load("Mean_trait_data.RData")

fact<- expand.grid(`Strength_mutualism`=seq(0.8,2.5,0.1), 
                   `web` ="datasets_1/M_PL_045.csv",
                   hysteresis=c("yes","no"), 
                   `forcing_strength`= 0,
                   `forcing_duration`=0,
                   `no_species_forced`=0,
                   `model`="abundance",
                   `interaction_type`= "trade_off", 
                   `individual.variation` ="high",
                   `random_seed`=4327+(1:1)*100) %>%
  as_tibble()

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
sig <-runif((Aspecies+Plantspecies),0.005,0.005)


## vector of species trait standard deviations
if(fact$hysteresis[r] == "yes"){
N <- runif( (Aspecies+Plantspecies) , 0,0.005)
}else{
  N <- runif( (Aspecies+Plantspecies) , 1,1)
} #mimicking a system with low population density near the collapse state
nainit<-N[1:Aspecies]
npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]

index_m<-which(names(trait_data_list)==fact$web[1])

ma<-  runif((Aspecies+Plantspecies), -0.5, 0.5) #(outt %>% filter(webname == fact$web[1]))$m
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
tmax<-1000
nestedness<-nestedness_NODF(g)
C<-Connectance(g)
web.name<-fact$web[1]
ba<-runif(Aspecies, 0.,0.)
bp<-runif(Plantspecies,0.,0.)
dganimals<-degree.animals
dgplants<-degree.plants
fact$Strength_mutualism[1]<-fact$Strength_mutualism[r]
mut.strength<-runif((Aspecies+Plantspecies), fact$Strength_mutualism[1],fact$Strength_mutualism[1])
mut.strength[index_max_degree]<- fact$Strength_mutualism[1] 
time_range<-c(0,tmax)
deltat<- (time_range[2]-time_range[1])/1 + 1
duration <- fact$forcing_duration[1]<- 0
d<- c(rep(1,duration),rep(0,(deltat-duration)))
duration_mat_A<-(replicate(Aspecies,d))
duration_mat_P<-(replicate(Plantspecies,d))
times<-seq(0, tmax, 1)
t1_A<-as.data.frame(list(times=times, import = rep(0,length(times))))
t1_P<-as.data.frame(list(times=times, import = rep(0,length(times))))
t1_A$import<-duration_mat_A[,1]
t1_P$import<-duration_mat_P[,1]
fact$forcing_strength[1]<-0
forcing_strength <- rep(0, (Aspecies+Plantspecies))# rep(fact$forcing_strength[1], (Aspecies+Plantspecies))
ic_f<-c(nainit, npinit, mainit,mpinit)

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

print(netwr_forcing_lvar)
proportion_richness<-length(which((sol1 %>% filter(time == tmax, type %in%c("N","P")))$v>0.5))/(Aspecies+Plantspecies)
richness<-length(which((sol1 %>% filter(time == tmax, type %in%c("N","P")))$v>0.5))

fact$prop_richness[r]<- proportion_richness
fact$richness[r]<- richness
print(r)
}



#plotting the hysteresis figure for a smaller network of 44 species
fact %>% ggplot(aes(x=Strength_mutualism, y = richness,color = hysteresis))+
  geom_point(size=5, alpha=0.5)+
  theme_classic()+
  annotate("rect", xmin = 1, xmax = 1.5, ymin = 0, ymax = 50,
           alpha = .2)+
  labs(color="hysteresis")+
  geom_curve(aes(x = 2.5, y = 45, xend = 1.45, yend = 0),
             arrow = arrow(length = unit(0.025, "npc")), 
             size=1.15, alpha = 0.25, curvature = 0.3, col="#CC79A7" )+
  annotate(geom = "text", x = 1.75, y = 44, label = "Collapse regime", color="#CC79A7", 
           size=3.5,  hjust = "left")+
  geom_curve(aes(x = 1.5, y = 1, xend = 2.5, yend = 30),
             arrow = arrow(length = unit(0.025, "npc")), 
             size=1.15, curvature = 0.3, col = "#009E73")+
  annotate(geom = "text", x = 2, y = 15 ,label ="Recovery regime",
           color="#009E73", size=3.5,  hjust = "left")+
  xlab( expression(paste("Avg. mutualistic strength,", gamma[0])))+
  ylab("Richness")






