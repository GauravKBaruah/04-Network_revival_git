
library(foreach)
library(doMC)
library(statmod)
library(sna)
library(bipartite)
require(tidyverse) ## for efficient data manipulation & plotting
library(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(beepr)
library(viridis)
library(ggdist)

# reading all the datasets
# calculating nestedness and connectance
mydir = 'datasets_1'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:154]


load("figs10.RData")

webnames<- unique(webdat_2$web.name)

for(r in 1:nrow(webdat_2)){
  
  g<-adj.mat(newfiles[which(newfiles == as.character(webdat_2$web.name[r]))])
  w_nestedness<-networklevel(g,c("weighted nestedness"))
  w_nodf<-networklevel(g,c("weighted NODF"))
  betn_temp<-betweenness(g)
  normalised_betweeness<-(betn_temp - min(betn_temp))/(max(betn_temp)- min(betn_temp))
  
  webdat_2$modularity[r] <- NOS(web=g)$mod
  webdat_2$median_centrality[r] <- median(normalised_betweeness,na.rm=T)
  webdat_2$sd_centrality[r] <- sd(normalised_betweeness,na.rm=T)
  webdat_2$weighted_nestedness[r] <- w_nestedness
  webdat_2$weighted_NODF[r] <-w_nodf
  
  spdat_3$modularity[r] <- NOS(web=g)$mod
  spdat_3$median_centrality[r] <-median(normalised_betweeness,na.rm=T)
  spdat_3$sd_centrality[r] <- sd(normalised_betweeness,na.rm=T)
  spdat_3$weighted_nestedness[r] <- w_nestedness
  spdat_3$weighted_NODF[r] <- w_nodf
  print(r)
}

#save(webdat_2,file = "figs10.RData")
##################### Supplementary analysis- of modularity, centrality, and weighted nestedness ##########
webdat_2$modularity_group<- cut(webdat_2$modularity,
                                breaks = c(0.48, 0.5, 0.6, 0.7, 0.8, 0.9), 
                                labels=c("<0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8",">0.8"))


webdat_2$centrality_group<- cut(webdat_2$median_centrality,
                                breaks = c(0.005, 0.01, 0.03, 0.05, 0.07, 0.08), 
                                labels=c("<0.01", "0.01-0.03", "0.03-0.05", "0.05-0.07", ">0.07"))


webdat_2$weighted_nestedness_group<- cut(webdat_2$weighted_nestedness,
                                         breaks = c(0.15, 0.30, 0.45,  0.6, 0.75, 0.9), 
                                         labels=c("<0.3", "0.3-0.45", "0.45-0.6", "0.6-0.75", ">0.75"))


webdat_2$sd_centrality_group<- cut(webdat_2$sd_centrality,
                                   breaks = c(50, 100, 300, 500, 700, 900, 1000), 
                                   labels=c("<100", "100-300", "300-500", "500-700","700-900", ">900"))




###### Supplementary figures: Modularity, Weighted Nestedness, Avg. Centrality #########


####### supplementary data figures on weightedness, modularity and centrality ################

(sfig_1<-webdat_2 %>% filter(forcing_strength == 0.5, mut_strength>1) %>% 
   ggplot(aes(x=mut_strength, y = recovery_richness, color= factor(modularity_group)))+
   geom_point(position=position_jitter(height=0.0,width=0.0),
              alpha = 0.05, size = 3)+
   scale_color_viridis(discrete = TRUE)+
   ylab("Recovery richness")+
   #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
   labs(color="Modularity bins")+
   theme_classic()+
   # theme(legend.position="none")+
   xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
   # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
   stat_smooth(method = "glm", size=2,alpha=0.1,
               method.args = list(family = "quasibinomial"),
               se =FALSE,
               aes(color=factor(modularity_group))) +
   facet_wrap(~individual_variation))


(sfig_2<-webdat_2 %>% filter(forcing_strength == 0.5, mut_strength>1) %>% 
    ggplot(aes(x=mut_strength, y = recovery_richness, color= factor(centrality_group)))+
    geom_point(position=position_jitter(height=0.0,width=0.0),
               alpha = 0.05, size = 3)+
    scale_color_viridis(discrete = TRUE)+
    ylab("Recovery richness")+
    #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
    labs(color="median betweeness centrality bins")+
    theme_classic()+
    # theme(legend.position="none")+
    xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
    # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
    stat_smooth(method = "glm", size=2,alpha=0.1,
                method.args = list(family = "quasibinomial"),
                se =FALSE,
                aes(color=factor(centrality_group))) +
    facet_wrap(~individual_variation))



(sfig_3<-webdat_2 %>% filter(forcing_strength == 0.5, mut_strength>1) %>% 
    ggplot(aes(x=mut_strength, y = recovery_richness, color= factor(weighted_nestedness_group)))+
    geom_point(position=position_jitter(height=0.0,width=0.0),
               alpha = 0.05, size = 3)+
    scale_color_viridis(discrete = TRUE)+
    ylab("Recovery richness")+
    #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
    labs(color="Wighted Nestedness bins")+
    # theme(legend.position="none")+
    xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
    theme_classic()+
    # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
    stat_smooth(method = "glm", size=2,alpha=0.1,
                method.args = list(family = "quasibinomial"),
                se =FALSE,
                aes(color=factor(weighted_nestedness_group))) +
    facet_wrap(~individual_variation))

(sfig_4<-webdat_2 %>% filter(forcing_strength == 0.5, mut_strength>1) %>% 
    ggplot(aes(x=mut_strength, y = recovery_richness, color= factor(sd_centrality_group)))+
    geom_point(position=position_jitter(height=0.0,width=0.0),
               alpha = 0.05, size = 3)+
    scale_color_viridis(discrete = TRUE)+
    ylab("Recovery richness")+
    #  labs(color=expression(paste(gamma[0],"," ,"mutualism strength")))+
    labs(color="SD centrality bins")+
    # theme(legend.position="none")+
    xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
    theme_classic()+
    # xlab(expression(paste(gamma[0],"," ,"mutualism strength")))+
    stat_smooth(method = "glm", size=2,alpha=0.1,
                method.args = list(family = "quasibinomial"),
                se =FALSE,
                aes(color=factor(sd_centrality_group))) +
    facet_wrap(~individual_variation))




ggpubr::ggarrange(sfig_1,sfig_2,sfig_3,
                  nrow=3,ncol=1, labels =c("A", "B", "C"))




###### species level data #################################################
load("figS11_S13.RData")



str(spdat_3)

#spdat_3$w<-as.factor(spdat_3$web)
spdat_3$Degree<-as.numeric(as.character(spdat_3$Degree))
spdat_3$Species<-as.numeric(as.character(spdat_3$Species))
spdat_3$response_thr<-as.numeric(as.character(spdat_3$response_thr))
spdat_3$Nestedness<-as.numeric(as.character(spdat_3$Nestedness))
spdat_3$Connectance<-as.numeric(as.character(spdat_3$Connectance))
spdat_3$response_time<-as.numeric(as.character(spdat_3$response_time))
spdat_3$Individual_variation<-as.factor(as.character(spdat_3$Individual_variation))
spdat_3$Network_size<-as.numeric(as.character(spdat_3$Network_size))
spdat_3$mutualism_strength<-as.numeric(as.character(spdat_3$mutualism_strength))
spdat_3$abundance<-as.numeric(as.character(spdat_3$abundance))
spdat_3$Network_respose_time<-as.numeric(as.character(spdat_3$Network_respose_time))

str(spdat_3)



spdat_3$Individual_variation <- plyr::revalue(spdat_3$Individual_variation, c( "high" = "high variation", "low" = "low variation"))



spdat_3$fraction_abundance<- 0

spdat_3$fraction_abundance[which(spdat_3$abundance>= 0.5) ]<-1
spdat_3$fraction_abundance[which(spdat_3$abundance< 0.5) ]<-0



appender <- function(string) 
  latex2exp::TeX(paste("$\\gamma_{\\0} = $", string))  



(a2<-spdat_3 %>% group_by(Nestedness,Connectance, forcing_strength,median_centrality,weighted_nestedness,modularity,
                          Individual_variation, Network_size , mutualism_strength) %>%
    summarise(count=n(),
              count_abundance =  sum(abundance >= 0.5),
              proportion_abundance = count_abundance/count ) %>%
    filter(forcing_strength == 0.5,mutualism_strength > 1.1) %>% 
    ggplot( aes(y = (proportion_abundance),
                x = weighted_nestedness,
                colour = Individual_variation))+
    geom_point(size=4)+
    ylim(c(0,1))+
    geom_smooth(method = "glm", 
                method.args = list(family = "quasibinomial"), 
                se = T, size=1.5) +
    #geom_smooth(method = lm, formula = y ~x, se = F)+
    theme_classic()+
    theme(legend.position = "right")+
    xlab("Weighted nestedness")+
    scale_color_manual(values=c( "#E69F00", "#56B4E9"))+
    #labs(color= expression(paste(gamma[0])))+
    labs(y = expression(paste("Proportion of species with ", N[i] > 0.5 )))+ 
    facet_wrap(.~mutualism_strength, nrow = 2,ncol = 4, labeller = as_labeller(appender, default = label_parsed)))


(a3<-spdat_3 %>% group_by(Nestedness,Connectance, forcing_strength,median_centrality,weighted_nestedness,modularity,
                          Individual_variation, Network_size , mutualism_strength) %>%
    summarise(count=n(),
              count_abundance =  sum(abundance >= 0.5),
              proportion_abundance = count_abundance/count ) %>%
    filter(forcing_strength == 0.5,mutualism_strength > 1.1) %>% 
    ggplot( aes(y = (proportion_abundance),
                x = median_centrality,
                colour = Individual_variation))+
    geom_point(size=4)+
    ylim(c(0,1))+
    geom_smooth(method = "glm", 
                method.args = list(family = "quasibinomial"), 
                se = T, size=1.5) +
    #geom_smooth(method = lm, formula = y ~x, se = F)+
    theme_classic()+
    theme(legend.position = "right")+
    xlab("median betweeness centrality")+
    scale_color_manual(values=c( "#E69F00", "#56B4E9"))+
    #labs(color= expression(paste(gamma[0])))+
    labs(y = expression(paste("Proportion of species with ", N[i] > 0.5 )))+ 
    facet_wrap(.~mutualism_strength, nrow = 2,ncol = 4, labeller = as_labeller(appender, default = label_parsed)))




(a4<-spdat_3 %>% group_by(Nestedness,Connectance, forcing_strength,median_centrality,weighted_nestedness,modularity,
                          Individual_variation, Network_size , mutualism_strength) %>%
    summarise(count=n(),
              count_abundance =  sum(abundance >= 0.5),
              proportion_abundance = count_abundance/count ) %>%
    filter(forcing_strength == 0.5,mutualism_strength > 1.1) %>% 
    ggplot( aes(y = (proportion_abundance),
                x = modularity,
                colour = Individual_variation))+
    geom_point(size=4)+
    ylim(c(0,1))+
    geom_smooth(method = "glm", 
                method.args = list(family = "quasibinomial"), 
                se = T, size=1.5) +
    #geom_smooth(method = lm, formula = y ~x, se = F)+
    theme_classic()+
    theme(legend.position = "right")+
    xlab("Modularity")+
    scale_color_manual(values=c( "#E69F00", "#56B4E9"))+
    #labs(color= expression(paste(gamma[0])))+
    labs(y = expression(paste("Proportion of species with ", N[i] > 0.5 )))+ 
    facet_wrap(.~mutualism_strength, nrow = 2,ncol = 4, labeller = as_labeller(appender, default = label_parsed)))


