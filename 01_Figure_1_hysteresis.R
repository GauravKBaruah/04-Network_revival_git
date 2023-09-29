

rm(list=ls())
source("01_hyst_functions")

#require(deSolve) ## for integrating ordinary differential equations
#require(tidyverse) ## for efficient data manipulation & plotting
#require(cowplot) ## for arranging plots in a grid
#library(dplyr)
#library(readr)
#library(beepr)
#library(viridis)

library(foreach)
library(doMC)
library(statmod)
#require(deSolve) ## for integrating ordinary differential equations
#require(tidyverse) ## for efficient data manipulation & plotting
#library(cowplot) ## for arranging plots in a grid
library(dplyr)
#library(readr)
library(beepr)
#library(viridis)
numcores<- 20
registerDoMC(numcores)


#theme_set(theme_classic()) 


# reading all the datasets
# calculating nestedness and connectance
mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles

trait_data_list<-list()
trait_data<-expand.grid(web=c(newfiles)) %>% 
  as_tibble %>% mutate(networksize=0) 



for(i in 1:nrow(trait_data)){
  
  g<-adj.mat(myfiles[which(myfiles == trait_data$web[i])]) #network web names
  # g<-g[-1,-1] 
  
  
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species
  
  
  nestedness<-nestedness_NODF(g)
  C<-Connectance(g)
  
#  trait_data$nestedness[i]<- nestedness
 # trait_data$connectance[i]<- C
  trait_data$networksize[i] <- Aspecies+Plantspecies
  muinit <-runif((Aspecies+Plantspecies), -0.5,0.5)
  
  trait_data_list[i]<-list(muinit)
  
  names(trait_data_list)[i] <- as.character(trait_data$web[i])
  
}


fact2<- trait_data %>% filter(networksize < 185)
webfiles<-fact2$web
#creating the final dataframe over which all the data will be simulated and collected
fact<- expand.grid(`Strength_mutualism`=seq(0.1, 5, 0.25),
                   `web` = webfiles[1:70],
                   `h2`= 0.25,
                   `model`="abundance",
                    hysteresis_check=c("yes", "no"),
                   `interaction_type`= "trade_off", 
                   `random_seed`=4327+(1:1)*100) %>%
  as_tibble %>%
  mutate(biomass=0,
         richness =0)
model.t<-list()


set.seed(1234)

new_ddf<-NULL



outt<-foreach(r = 1:nrow(fact))%dopar%{
  #print(r)
  
  
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
  
  h2<-fact$h2[r]
  
  ## vector of species trait standard deviations
  
  if(fact$hysteresis_check[r] == "yes"){
  N <- runif( (Aspecies+Plantspecies) , 0,0.005)}else if(fact$hysteresis_check[r] == "no"){
    N <- runif( (Aspecies+Plantspecies) , 1,1)
  }  ## initial species densities
  index<-which(names(trait_data_list)==fact$web[r])
  
  muinit<-trait_data_list[[index]]
  mainit<-muinit[1:Aspecies]
  mpinit<-muinit[(Aspecies+1): (Aspecies+Plantspecies)]
  
  nainit<- N[1:Aspecies]
  npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
  
  Amatrix<-mat.comp(g,degree.animals,degree.plants)$Amatrix
  Pmatrix<-mat.comp(g,degree.animals,degree.plants)$Pmatrix
  gamma=0.35#fact_lessvar$Strength_mutualism[r]
  mut.strength<-runif( (Aspecies+Plantspecies), fact$Strength_mutualism[r],fact$Strength_mutualism[r])
  nestedness<-nestedness_NODF(g)
  C<-Connectance(g)
  web.name<-fact$web[r]
  ba<-runif(Aspecies, 0.05,0.05)
  bp<-runif(Plantspecies,0.05,0.05)
  dganimals<-degree.animals
  dgplants<-degree.plants
  
  
  ic <-c(nainit, npinit, mainit,mpinit)
  
  #fact$noise
  params <- list(time=time,matrix=g,sig=sig,Amatrix=Amatrix,
                 Pmatrix=Pmatrix,w=gamma,model=fact$model[r],
                 interaction_type=fact$interaction_type[r],
                 hysteresis_check=fact$hysteresis_check[r],
                 mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                 web.name=web.name,h2=h2, ba=ba,bp=bp,dganimals=dganimals,
                 dgplants=dgplants)
  
  
  
  
  start.time =3000
  
   model.t<-eqs(time = start.time,state = ic,pars = params)
  #system.time( model.t<-lapply(1, Mcommunity,time=start.time,state=ic,
   #                            pars=params))
  
  
  
  #par(mfrow=c(2,1))
 # ts.plot(model.t$Animals,gpars= list(col=viridis(40)),lwd=2.5,   ylab = "Plant Abundance")
  #ts.plot(model.t$Plant.trait,gpars= list(col=viridis(40)),lwd=2.5,   ylab = "Plant Abundance")
  
  
  
  #ts.plot(model.t[[1]]$Animals,gpars= list(col=viridis(40)),lwd=2.5,
  #       xlab = "Time", ylab = "Animal Abundance")
  # 
  #plot_snapshot(Na = model.t[[1]]$Animals[1000,],
  #             Np = model.t[[1]]$Plants[1000,],
  #             m = c(model.t[[1]]$Animal.trait[1000,], model.t[[1]]$Plant.trait[1000,]),
  #           sigma =sig, moment=0, limits=c(-1, 1), res=1001)
  # #
  # pbiomass<-sum( colMeans(model.t$Plants[200:800,]))
  # abiomass<-sum( colMeans(model.t$Animals[200:800,]))
  # 
  # 
  # fact$Biomass[r] = abiomass+pbiomass
  # fact$richness[r] =  (length(which(model.t$Plants[799,] > 0.5))+length(which(model.t$Animals[799,] > 0.5)))
  # fact$Nestedness[r] = nestedness_NODF(g)
  # fact$Connectance[r] = Connectance(g)
  # fact$mean.trait.matching[r] = trait.matching(mA=model.t$Animal.trait[799,],
  #                                              mP = model.t$Plant.trait[799,], 
  #                                              adj.mat = g, gamma=gamma
  # )
  # fact$community.weighted.mean[r] = community.w.mean(mA=model.t[[1]]$Animal.trait[800,],
  #                                                    mP = model.t[[1]]$Plant.trait[800,],
  #                                                    Na =  model.t[[1]]$Animals[800,],
  #                                                    Np = model.t[[1]]$Plants[800,])$cwm_community
  # 
  # fact$community.weighted.mean.plants[r] = community.w.mean(mA=model.t[[1]]$Animal.trait[800,],
  #                                                    mP = model.t[[1]]$Plant.trait[800,],
  #                                                    Na =  model.t[[1]]$Animals[800,],
  #                                                    Np = model.t[[1]]$Plants[800,])$cwm_plants
  # 
  # 
  # fact$community.weighted.mean.animals[r] = community.w.mean(mA=model.t[[1]]$Animal.trait[800,],
  #                                                           mP = model.t[[1]]$Plant.trait[800,],
  #                                                           Na =  model.t[[1]]$Animals[800,],
  #                                                           Np = model.t[[1]]$Plants[800,])$cwm_animals
  # 
  
 # print(r)
  
  
}




save(outt,file="hysteresis.RData")
#






load("hysteresis.RData")

hyst<-NULL
for(i in 1:2800){
  
  hyst <- rbind(hyst,outt[[i]])
  
  
}


c1<-hyst %>% filter(hysteresis== "no") %>% 
ggplot(aes(x=mut_strength , y = biomass,color=factor(hysteresis)))+
  geom_point(size=5,alpha=0.7)+
  xlim(c(0,4))+
  theme_classic()+
  theme(legend.position = "")+
  scale_color_manual(values= c("#CC79A7"))+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = 1500,
          alpha = .2)+
  labs(color="")+
  #ggtitle("D")+
  xlab( expression(paste("Avg. mutualistic strength,", gamma[0])))+
  ylab("Network biomass")


c2<-hyst %>% filter(hysteresis== "no") %>% 
  ggplot(aes(x=mut_strength , y = richness,color=factor(hysteresis)))+
  geom_point(size=5,alpha=0.7)+
  xlim(c(0,4))+
  theme_classic()+
  #ggtitle("E")+
  theme(legend.position = "")+
  scale_color_manual(values= c("#CC79A7"))+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = 200,
           alpha = .2)+
  labs(color="")+
  xlab( expression(paste("Avg. mutualistic strength,", gamma[0])))+
  ylab("Richness")

gridExtra::grid.arrange(c1,c2, nrow=1,ncol=2)



hyst %>% filter(hysteresis== "yes") %>% 
  ggplot(aes(x=mut_strength , y = richness,color=factor(hysteresis)))+
  geom_point(size=5,alpha=0.7)+
  xlim(c(0,4.5))+
  theme_classic()+
  theme(legend.position = "")+
  #ggtitle("B")+
  scale_color_manual(values= c("#CC79A7"))+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = 180,
           alpha = .2)+
  labs(color="")+
  xlab("Average mutualistic strength")+
  ylab("Richness")


(w1<-hyst %>% filter(web.name == "datasets_1/M_PL_018.csv" )%>% 
  ggplot(aes(x=mut_strength , y = richness,color=factor(hysteresis)))+
  geom_point(size=5,alpha=0.7)+
  xlim(c(0,4))+
  theme_classic()+
 # ggtitle("B")+
  theme(legend.position = "")+
  scale_color_manual(values= c("#009E73", "#CC79A7"))+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = 150,
           alpha = .2)+
  labs(color="")+
  geom_curve(aes(x = 2.5, y = 145, xend = 1.45, yend = 75),
             arrow = arrow(length = unit(0.025, "npc")), 
             size=1.15, alpha = 0.25, curvature = 0.3, col="#CC79A7" )+
  annotate(geom = "text", x = 2.5, y = 150, label = "Collapse regime", color="#CC79A7", 
           size=3.5,  hjust = "left")+
  geom_curve(aes(x = 1.5, y = 5, xend = 2.5, yend = 40),
             arrow = arrow(length = unit(0.025, "npc")), 
             size=1.15, curvature = 0.3, col = "#009E73")+
    annotate(geom = "text", x = 1.25, y = 25 ,label ="Recovery regime",
             color="#009E73", size=3.5,  hjust = "left")+
    
    #scale_color_manual(values= c("#009E73", "#CC79A7"))+
  xlab( expression(paste("Avg. mutualistic strength,", gamma[0])))+
  ylab("Richness"))

(w2<-hyst %>% filter(web.name == "datasets_1/M_PL_018.csv" )%>% 
  ggplot(aes(x=mut_strength , y = biomass,color=factor(hysteresis)))+
  geom_point(size=5,alpha=0.5)+
  xlim(c(0,4))+
  theme_classic()+
 # ggtitle("C")+
  theme(legend.position = "")+
  scale_color_manual(values= c("#009E73", "#CC79A7"))+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = 1150,
           alpha = .2)+
  labs(color="")+
  geom_segment(aes(x = 2.5, y = 550, xend = 1.45, yend = 55),
             arrow = arrow(length = unit(0.025, "npc")), 
             size=1.15, alpha = 0.25, col="#CC79A7" )+
  annotate(geom = "text", x = 1.7, y = 300, label = "Collapse regime", color="#CC79A7",
           size=3.5,  angle=55)+
  geom_segment(aes(x = 1.45, y = 40, xend = 2.5, yend = 40),
             arrow = arrow(length = unit(0.025, "npc")), 
             size=1.15,  col = "#009E73")+
    annotate(geom = "text", x = 2.45, y = 95, label = "Recovery regime", color="#009E73",
             size=3.5,  angle=0)+
    
  xlab( expression(paste("Avg. mutualistic strength,", gamma[0])))+
  ylab("Biomass"))

g<-  adj.mat(as.character(webfiles[which(webfiles == "datasets_1/M_PL_018.csv" )])) 

Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant species
degree.animals<-degree.plants<-numeric()

#degree of plants and anichmals
for(i in 1:Plantspecies){
  degree.plants[i]<-sum(g[i,])} # degree of plants
for(j in 1:Aspecies){
  degree.animals[j]<-sum(g[,j]) # degree of animals
}

net1 = network(g, bipartite = T, directed = FALSE)

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#vertex names 
names<-network.vertex.names(net1)
net1 %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net1 %v% "color" = ifelse(net1 %v% "groups" == "plants", "#0072B2", "#E69F00" )
#ggnet2(net,  mode="circle",  color ="groups", edge.size = 1,edge.alpha = 1, edge.color = "black", edge.lty = 1)

deg<-c(degree.plants,degree.animals)
w6<-ggnet2(net1, mode="circle", size=deg, edge.size = 1,max_size =8, color ="color",edge.alpha = 1, legend.position = "",
           label ="A")
w6


#grid::ggarrnge(w6, w1,w2, c1,c2, nrow=2,ncol=3)
ggpubr::ggarrange(w6,w1,w2,c1,c2,ncol=3,nrow=2,
                  labels = c("A","B","C","D","E"))
#18. 28, 



# proportion exhibited hysteresis

for(r in 1 :nrow(hyst)){
g<-adj.mat(as.character(webfiles[which(webfiles == hyst$web.name[r])])) #network web names
# g<-g[-1,-1] 


Aspecies<- dim(g)[2] # no of animal species
Plantspecies<- dim(g)[1] # no of plant spec

hyst$network.size[r]<- Aspecies+Plantspecies
}




prop_net<-hyst %>% filter(hysteresis== "yes", mut_strength == 4.35)


for(i in 1:nrow(prop_net)){
  
  
  prop<-prop_net$richness[i]/prop_net$network.size[i] 
 if(prop > 0.8){
   hysteresis_presence<- 0
 }else{
   hysteresis_presence<- 1
 }
  
  prop_net$hysteresis_presence[i]<- hysteresis_presence
  prop_net$proportion_recovery[i]<-prop
  
}

sum(prop_net$hysteresis_presence)/nrow(prop_net)
h0<-ggplot(prop_net, aes(y=proportion_recovery, x= hysteresis, color=network.size))+
  geom_boxplot(fill = "white", outlier.position=NA)+
  geom_jitter(size=3,alpha=0.5,  position=position_jitter(0.2))+
  scale_colour_gradientn(colours = terrain.colors(10))+
  theme_bw()+
  ylim(c(0,1))+
  xlab("Presence of hysteresis")+
ylab("Recovery richness at highest mutualistic strength")

#

(h1<-hyst %>% filter(hysteresis== "yes") %>% 
  ggplot(aes(x=mut_strength , y = biomass,color=factor(hysteresis)))+
  geom_point(size=5,alpha=0.7)+
  xlim(c(0,4))+
  theme_classic()+
  theme(legend.position = "")+
  scale_color_manual(values= c("#CC79A7"))+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = 1500,
           alpha = .2)+
  labs(color="")+
    geom_segment(aes(x = 1.05, y =75, xend = 2.5, yend = 75),
                 arrow = arrow(length = unit(0.025, "npc")), 
                 size=1.15,  col = "#CC79A7")+
    annotate(geom = "text", x = 1.5, y = 115, label = "No recovery", color="#CC79A7",
             size=3.5,  angle=0)+
    
  #ggtitle("D")+
  xlab( expression(paste("Avg. mutualistic strength,", gamma[0])))+
  ylab("Network biomass"))


(h2<-hyst %>% filter(hysteresis== "yes") %>% 
  ggplot(aes(x=mut_strength , y = richness,color=factor(hysteresis)))+
  geom_point(size=5,alpha=0.7)+
  xlim(c(0,4))+
  theme_classic()+
  #ggtitle("E")+
  theme(legend.position = "")+
  scale_color_manual(values= c("#CC79A7"))+
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = 200,
           alpha = .2)+
  labs(color="")+
    geom_segment(aes(x = 1.05, y = 10, xend = 2.25, yend = 10),
                 arrow = arrow(length = unit(0.025, "npc")), 
                 size=1.15,  col = "#CC79A7")+
    annotate(geom = "text", x = 1.5, y = 15, label = "No recovery", color="#CC79A7",
             size=3.5,  angle=0)+
    
  xlab( expression(paste("Avg. mutualistic strength,", gamma[0])))+
  ylab("Richness"))

ggpubr::ggarrange(h1,h2, h0,ncol=3,nrow=1,
                  labels = c("A","B", "C"))
