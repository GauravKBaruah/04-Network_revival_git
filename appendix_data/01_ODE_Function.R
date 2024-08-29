#Functions used in the paper " Reviving collapsed plant-polliantor networks from a single species" 2024. Gaurav Baruah & Meike Witmann. Some of these functions are borrowed from various open access codes and Rscripts
#including Baruah 2022 Eco. Letts paper.
#email: gbaruahecoevo@gmail.com

require(statmod)

cutoff <- function(x) ifelse(x<1, (1*(x>0))*(x*x*x*(10+x*(-15+6*x))), 1)

network_structure<-function(Na,Np, g ){
  new_g<- matrix(0, nrow= length(Np), ncol =length(Na) )
  
  Na[which(Na < 0.05)]<-0
  Np[which(Np < 0.05)]<-0
  
  for(i in 1:length(Np)){
    for(j in 1:length(Na)){
      new_g[i,j]<-g[i,j]*Na[j]*Np[i]    
      
    }
  }
  new_g[which(new_g > 0)]<-1
  
  return(new_g) 
}


#conversion of .csv incidence adjacency matrices to a matrix format of plant-pollinator networks
adj.mat<-function(data){
  #dat <- paste('network.csv',sep='')
  d <- read.csv(file=data,header=FALSE )
  dat<-as.matrix(d)
  dat[dat > 0] = 1
  dat[dat < 0] = 
    dat<-apply(dat,2,as.numeric)
  return(dat)}



#function for estimating trait matching for a network
trait.matching<-function(mA,mP,adj.mat_1,gamma){
  tm<-numeric()
  for(i in 1:nrow(adj.mat_1)){
    tm[i] <- mean(adj.mat_1[i,]*exp(-(mA-mP[i])^2)/gamma)
    
  }
  return(tm=mean(tm))
}



#gaussian quadrature approximation used to numerically determine the integrals used in the paper 
gausquad.animals<-function(m,sigma,w,h,np,na,mut.strength,points,mat,degree.animal,interaction_type){
  
  
  temp2<-dat2<-x2<-x3<-array(dim=c(points))
  if(mat == 0){
    return(list(G= 0, B = 0))
  }
  else if(mat == 1){
    #nodes oir points in the abscissa where the integral will be evaluated numerically
    if(interaction_type == "trade_off"){
      z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z'
      z2<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w1<-gauss.quad.prob(points, dist = "normal", 
                          mu=m$ma,sigma =sigma$sa)$weights #pi(z')
      w2<-gauss.quad.prob(points, dist = "normal", 
                          mu=m$mp,sigma =sigma$sp)$weights #pj(z'')
      
      
      #for the pairwise model however there are only two species interacting and hence i and j
      #or in other words the integral goes over z and z'
      for (i in 1: points){
        
        
        f <-  exp(-(z1[i]- z2)^2/w^2) # + 2*alpha*(sign(z1[i] - z2))*(1- exp(-(z1[i]-z2)^2/w^2)) + sign(alpha))
        
        temp2[i]<- sum(np*(mut.strength/degree.animal)*f/(1+h*np*(mut.strength/degree.animal)*f)*w2*w1[i])
        
        #temp2[i]<- sum(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])
        #x2[i]<- sum(h*w2*f)
        #x3[i]<-sum(h*exp(-(z1[i]-z2)^2/w^2)*w2)
        #dat2[i]<- sum((z1[i]-m$mp)*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]))
        
        dat2[i]<- sum( ((z1[i]-m$ma)*f*(mut.strength/degree.animal)*np/(1+h*np*(mut.strength/degree.animal)*f) )*w2*w1[i])
        
        # x3[i] <-sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
        
      }
      G = sum(temp2)
      B = sum(dat2) 
      
    } else {
      
      z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z'
      z2<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w1<-gauss.quad.prob(points, dist = "normal", 
                          mu=m$ma,sigma =sigma$sa)$weights #pi(z')
      w2<-gauss.quad.prob(points, dist = "normal", 
                          mu=m$mp,sigma =sigma$sp)$weights #pj(z'')
      
      
      #for the pairwise model however there are only two species interacting and hence i and j
      #or in other words the integral goes over z and z'
      for (i in 1: points){
        
        f <-  exp(-(z1[i]- z2)^2/w^2) # + 2*alpha*(sign(z1[i] - z2))*(1- exp(-(z1[i]-z2)^2/w^2)) + sign(alpha))
        
        temp2[i]<- sum(np*(mut.strength)*f/(1+h*np*(mut.strength)*f)*w2*w1[i])
        
        #temp2[i]<- sum(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])
        #x2[i]<- sum(h*w2*f)
        #x3[i]<-sum(h*exp(-(z1[i]-z2)^2/w^2)*w2)
        #dat2[i]<- sum((z1[i]-m$mp)*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]))
        
        dat2[i]<- sum( ((z1[i]-m$ma)*f*(mut.strength)*np/(1+h*np*(mut.strength)*f) )*w2*w1[i])
        
        
      }
      G = sum(temp2)
      B = sum(dat2) 
    }
    
    
    return(list(G= G, B = B))
  }
}

#gaussian quadrature approximation used to numerically determine the integrals.
gausquad.plants<-function(m,sigma,w,h,np,na,mut.strength,points,mat,degree.plant,interaction_type){
  
  temp2<-dat2<-x3<-x4<-array(dim=c(points))
  
  if(mat == 0){
    
    return(list(G= 0, 
                B = 0))
    
  }
  else if (mat==1){
    if(interaction_type == "trade_off"){
      #nodes oir points in the abscissa where the integral will be evaluated numerically
      z1<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z'
      z2<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w1<-gauss.quad.prob(points, dist = "normal", 
                          mu=m$mp,sigma =sigma$sp)$weights #pi(z')
      w2<-gauss.quad.prob(points, dist = "normal", 
                          mu=m$ma,sigma =sigma$sa)$weights #pj(z'')
      
      
      #for the pairwise model however there are only two species interacting and hence i and j
      #or in other words the integral goes over z and z'
      for (i in 1: points){
        
        f <-  exp(-(z1[i]- z2)^2/w^2) # + 2*alpha*(sign(z1[i] - z2))*(1- exp(-(z1[i]-z2)^2/w^2)) + sign(alpha))
        
        temp2[i]<- sum(na*(mut.strength/degree.plant)*f/(1+h*na*(mut.strength/degree.plant)*f)*w2*w1[i])
        
        dat2[i]<- sum( ((z1[i]-m$mp)*f*(mut.strength/degree.plant)*na/(1+h*na*(mut.strength/degree.plant)*f) )*w2*w1[i])
        
        # x4[i]<-sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
      }
      
      G= sum(temp2)
      B = sum(dat2)
      
    } else {
      
      z1<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z'
      z2<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z''
      
      #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
      w1<-gauss.quad.prob(points, dist = "normal", 
                          mu=m$mp,sigma =sigma$sp)$weights #pi(z')
      w2<-gauss.quad.prob(points, dist = "normal", 
                          mu=m$ma,sigma =sigma$sa)$weights #pj(z'')
      
      
      #for the pairwise model however there are only two species interacting and hence i and j
      #or in other words the integral goes over z and z'
      for (i in 1: points){
        f <-  exp(-(z1[i]- z2)^2/w^2) # + 2*alpha*(sign(z1[i] - z2))*(1- exp(-(z1[i]-z2)^2/w^2)) + sign(alpha))
        
        temp2[i]<- sum(na*(mut.strength)*f/(1+h*na*(mut.strength)*f)*w2*w1[i])

        
        dat2[i]<- sum( ((z1[i]-m$mp)*f*(mut.strength)*na/(1+h*na*(mut.strength)*f) )*w2*w1[i])
        
        # x4[i]<-sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
      }
      
      G= sum(temp2)
      B = sum(dat2)
      
    }
    
    
    return(list(G= G, 
                B = B))
  }
  
}


# this is the function that numerically simulates the eco-evo dynamics of plant-pollinator networks without any species specific perturbation
#time: time
#state: initial state values
#pars: list of parameters

eqs <- function(time, state, pars) {
  A <- dim(pars$matrix)[2]  ## number of animal species
  P <-dim(pars$matrix)[1]
  Np<-state[(A+1):(A+P)]
  Na<-state[1:A]
  s <- pars$sig ## species' trait standard deviations
  ma<-state[(A+P+1):(A+P+A)]
  mp<-state[(A+P+A+1):(A+P+A+P)]
  ## define g, where g[i] is the selection pressure on species i from growth
  alpha.a<-pars$Amatrix ## alpha matrix
  alpha.p<-pars$Pmatrix
  
  aij<-bij<-matrix(0, nrow=A,ncol=P) 
  aji<-bji<-matrix(0, nrow=P,ncol=A) 
  muA<-ma
  muP<-mp
  aj<-bj<-ai<-bi<-numeric()
  
  
  for(r in 1:A){
    for(l in 1:P){
      #
      m.temp<-list(ma=muA[r],mp=muP[l])
      sigma1<-list(sa=s[r],sp=s[(A)+l])
      temp1<-gausquad.animals(m=m.temp,sigma=sigma1,w=pars$w,h=0.25,np=Np[l],na=Na[r],
                              mut.strength=pars$mut.strength[r], points=5
                              ,mat=pars$matrix[l,r],
                              degree.animal = pars$dganimals[r],
                              interaction_type=pars$interaction_type)
      aij[r,l] <-temp1$G
      bij[r,l] <-temp1$B
      
    }
    ai[r]<-sum(aij[r,])
    bi[r]<-sum(bij[r,])
  }
  for(k in 1:P){
    for(m in 1:A){
      m2.temp<-list(ma=muA[m],mp=muP[k])
      sigma2<-list(sa=s[m],sp=s[(A+k)])
      temp2<-gausquad.plants(m=m2.temp,sigma=sigma2,w=pars$w,h=0.25,np=Np[k],na=Na[m],
                             mut.strength=pars$mut.strength[A+k],
                             points=5,mat=pars$matrix[k,m], 
                             degree.plant =pars$dgplants[k],
                             interaction_type=pars$interaction_type)
      aji[k,m] <-temp2$G
      bji[k,m]<-temp2$B
    }
    aj[k]<-sum(aji[k,])
    bj[k]<-sum(bji[k,])
  }
  #print(t)
  
  
  
  dndt_a<- Na*(pars$ba-alpha.a%*%Na+ai)*cutoff(Na/(1e-8)) #  na*(ba-alpha.a%*%na+ai)*cutoff(na/(1e-8)) #population dynamics
  dndt_p<- Np*(pars$bp-alpha.p%*%Np+aj)*cutoff(Np/(1e-8))  #population dynamics
  dudt_A<- pars$h2[1]*(bi) #mean trait dynamics
  dudt_P<- pars$h2[1]*(bj) #mean trait dynamics
  
  ## return equations by first flattening them back into a single vector
  return(list(c(dndt_a, dndt_p,dudt_A,dudt_P)))
}


# this is the function that numerically simulates the eco-evo dynamics of plant-pollinator networks WITH  species specific perturbation
#time: time
#state: initial state values
#pars: list of parameters
eqs_perturbation <- function(t, state, pars) {
  A <- dim(pars$matrix)[2]  ## number of animal species
  P <-dim(pars$matrix)[1]
  Np<-state[(A+1):(A+P)]
  Na<-state[1:A]
  s <- pars$sig ## species' trait standard deviations
  ma<-state[(A+P+1):(A+P+A)]
  mp<-state[(A+P+A+1):(A+P+A+P)]
  alpha.a<-pars$Amatrix ## alpha matrix
  alpha.p<-pars$Pmatrix
  aij<-bij<-matrix(0, nrow=A,ncol=P) 
  aji<-bji<-matrix(0, nrow=P,ncol=A) 
  muA<-ma
  muP<-mp
  aj<-bj<-ai<-bi<-numeric()
  time_func_A<-approxfun(pars$t1_A, method="linear", rule =2 )
  time_func_P<-approxfun(pars$t1_P, method="linear", rule= 2)
  e_A<-replicate(A,time_func_A(t))
  e_P<-replicate(P,time_func_P(t))
  
  forcing_strength_A<-pars$forcing_strength[1:A]
  forcing_strength_P<-pars$forcing_strength[ (A+1):(A+P)]
  
  forcing_index<-rep(0,(A+P))
  forcing_index[pars$species_index]<-1
  forcing_index_sp_A<-forcing_index[1:A]
  forcing_index_sp_P<-forcing_index[(A+1):(A+P)]
  for(r in 1:A){
    for(l in 1:P){
      #
      m.temp<-list(ma=muA[r],mp=muP[l])
      sigma1<-list(sa=s[r],sp=s[(A)+l])
      temp1<-gausquad.animals(m=m.temp,sigma=sigma1,w=pars$w,h=0.25,np=Np[l],na=Na[r],
                              mut.strength=pars$mut.strength[r], points=6
                              ,mat=pars$matrix[l,r],
                              degree.animal = pars$dganimals[r],
                              interaction_type=pars$interaction_type)
      aij[r,l] <-temp1$G
      bij[r,l] <-temp1$B
      
    }
    ai[r]<-sum(aij[r,])
    bi[r]<-sum(bij[r,])
  }
  for(k in 1:P){
    for(m in 1:A){
      m2.temp<-list(ma=muA[m],mp=muP[k])
      sigma2<-list(sa=s[m],sp=s[(A+k)])
      temp2<-gausquad.plants(m=m2.temp,sigma=sigma2,w=pars$w,h=0.25,np=Np[k],na=Na[m],
                             mut.strength=pars$mut.strength[A+k],
                             points=6,mat=pars$matrix[k,m], 
                             degree.plant =pars$dgplants[k],
                             interaction_type=pars$interaction_type)
      aji[k,m] <-temp2$G
      bji[k,m]<-temp2$B
    }
    aj[k]<-sum(aji[k,])
    bj[k]<-sum(bji[k,])
  }  
  
  
  dndt_a<- Na*(pars$ba-alpha.a%*%Na+ai+forcing_index_sp_A*forcing_strength_A*e_A)*cutoff(Na/(1e-8))  #  na*(ba-alpha.a%*%na+ai)*cutoff(na/(1e-8)) #population dynamics
    dndt_p<- Np*(pars$bp-alpha.p%*%Np+aj+forcing_index_sp_P*forcing_strength_P*e_P)*cutoff(Np/(1e-8))   #population dynamics
  dudt_A<- pars$h2[1]*(bi) #mean trait dynamics
  dudt_P<- pars$h2[1]*(bj) #mean trait dynamics
  
  ## return equations by first flattening them back into a single vector
  return(list(c(dndt_a, dndt_p,dudt_A,dudt_P)))
}



# this is the function that numerically simulates the eco-evo dynamics of plant-pollinator networks WITH species specific perturbation but a specific density and not rate
#time: time
#state: initial state values
#pars: list of parameters
eqs_perturbation_c <- function(t, state, pars) {
  A <- dim(pars$matrix)[2]  ## number of animal species
  P <-dim(pars$matrix)[1]
  Np<-state[(A+1):(A+P)]
  Na<-state[1:A]
  s <- pars$sig ## species' trait standard deviations
  ma<-state[(A+P+1):(A+P+A)]
  mp<-state[(A+P+A+1):(A+P+A+P)]
  alpha.a<-pars$Amatrix ## alpha matrix
  alpha.p<-pars$Pmatrix
  aij<-bij<-matrix(0, nrow=A,ncol=P) 
  aji<-bji<-matrix(0, nrow=P,ncol=A) 
  muA<-ma
  muP<-mp
  aj<-bj<-ai<-bi<-numeric()
  time_func_A<-approxfun(pars$t1_A, method="linear", rule =2 )
  time_func_P<-approxfun(pars$t1_P, method="linear", rule= 2)
  e_A<-replicate(A,time_func_A(t))
  e_P<-replicate(P,time_func_P(t))
  
  forcing_strength_A<-pars$forcing_strength[1:A]
  forcing_strength_P<-pars$forcing_strength[ (A+1):(A+P)]
  
  forcing_index<-rep(0,(A+P))
  forcing_index[pars$species_index]<-1
  forcing_index_sp_A<-forcing_index[1:A]
  forcing_index_sp_P<-forcing_index[(A+1):(A+P)]

  #mutualistic strength of animals on plants
  for(r in 1:A){
    for(l in 1:P){
      #
      m.temp<-list(ma=muA[r],mp=muP[l])
      sigma1<-list(sa=s[r],sp=s[(A)+l])
      temp1<-gausquad.animals(m=m.temp,sigma=sigma1,w=pars$w,h=0.25,np=Np[l],na=Na[r],
                              mut.strength=pars$mut.strength[r], points=6
                              ,mat=pars$matrix[l,r],
                              degree.animal = pars$dganimals[r],
                              interaction_type=pars$interaction_type)
      aij[r,l] <-temp1$G
      bij[r,l] <-temp1$B
      
    }
    ai[r]<-sum(aij[r,])
    bi[r]<-sum(bij[r,])
  }

  #mutualistic strength of plants on animals
  for(k in 1:P){
    for(m in 1:A){
      m2.temp<-list(ma=muA[m],mp=muP[k])
      sigma2<-list(sa=s[m],sp=s[(A+k)])
      temp2<-gausquad.plants(m=m2.temp,sigma=sigma2,w=pars$w,h=0.25,np=Np[k],na=Na[m],
                             mut.strength=pars$mut.strength[A+k],
                             points=6,mat=pars$matrix[k,m], 
                             degree.plant =pars$dgplants[k],
                             interaction_type=pars$interaction_type)
      aji[k,m] <-temp2$G
      bji[k,m]<-temp2$B
    }
    aj[k]<-sum(aji[k,])
    bj[k]<-sum(bji[k,])
  }

  
  
  
  dndt_a<- Na*(pars$ba-alpha.a%*%Na+ai)*cutoff(Na/(1e-8)) + forcing_index_sp_A*forcing_strength_A*e_A*cutoff(Na/(1e-8)) #  na*(ba-alpha.a%*%na+ai)*cutoff(na/(1e-8)) #population dynamics
  dndt_p<- Np*(pars$bp-alpha.p%*%Np+aj)*cutoff(Np/(1e-8)) + forcing_index_sp_P*forcing_strength_P*e_P*cutoff(Np/(1e-8))  #population dynamics
  dudt_A<- pars$h2[1]*(bi) #mean trait dynamics
  dudt_P<- pars$h2[1]*(bj) #mean trait dynamics
  
  ## return equations by first flattening them back into a single vector
  return(list(c(dndt_a, dndt_p,dudt_A,dudt_P)))
}



# this function is designed for the THEOBIOTA cluster at University of Bielefeld.
#this function includes 
#params_Forcing: list of parameters
# ic_f : initlat state variables
#tmax - max time
# also includes ode solver and uses the above eqs_perturbation() function to simulate eco-evo dynamics with species specific forcing and
# lets out recovery species richness at max time point, biomass, degree, nestedness, connectance etc. 
 cluster_run_func<-function(params_forcing, ic_f, tmax ){
  
  
  sol<-ode(func=eqs_perturbation, y=ic_f, parms=params_forcing, times=seq(0, tmax, by=1)) %>% 
    organize_results(pars = params_forcing)#%>% plot_all() ## solve ODEs
  
  
  
  newg_recovery<-network_structure(Na=(sol %>% filter(time ==tmax, type !="ma",type !="mp", type =="N"))$v,
                                   Np=(sol %>% filter(time ==tmax, type !="ma",type !="mp", type =="P"))$v, 
                                   g=params_forcing$matrix )
  
  A<-dim(params_forcing$matrix)[2]
  P<-dim(params_forcing$matrix)[1]
  
  rr<-response_time(sol=sol,
                    Aspecies = A, Plantspecies = P, g = params_forcing$matrix,t = tmax, 
                    dat = params_forcing$dat, webname =  as.character(params_forcing$web.name))
  
  
  richness_perturb<- length(which((sol %>%  filter(time == tmax, type !="ma",type !="mp"))$v >= 5e-1))
  
  na_1<-  sol %>%  filter(time > params_forcing$duration, type !="ma",type !="mp", type =="N")
  np_1<-   sol  %>%  filter(time > params_forcing$duration, type !="ma",type !="mp", type =="P") 
  na <- spread(na_1, key = species, v = v)
  np<- spread(np_1, key = species, v = v)
  
  # total equilibrium biomass after perturbation during the recovery simulation
  pbiomass_perturb<-sum( colMeans(np[,3:P]))
  abiomass_perturb<-sum( colMeans(na[,3:A]))
  
  
  
  
  
  
  
  recovery_biomass <- (pbiomass_perturb+abiomass_perturb) 
  recovery_richness <- richness_perturb/(A+P)
  network_response_time<- mean(rr$Response_time,na.rm=T)
  
  perturb_biomass.animal = abiomass_perturb
  perturb_biomass.plant = pbiomass_perturb
  Nestedness =nestedness_NODF(params_forcing$matrix)
  connectance = Connectance(params_forcing$matrix)
  recovery_nestedness = abs(nestedness_NODF(newg_recovery)-0)
  recovery_connectance = abs(Connectance(newg_recovery) - 0)
  #  
  
  
  #  
  out<- data.frame(recovery_biomass=recovery_biomass,
                   recovery_richness=recovery_richness,
                   network_response_time=network_response_time,
                   perturb_biomass.animal=perturb_biomass.animal,
                   perturb_biomass.plant=perturb_biomass.plant,
                   Nestedness=Nestedness,
                   connectance=connectance,
                   web.name=as.character(params_forcing$web.name),
                   recovery_nestedness=recovery_nestedness,
                   recovery_connectance=recovery_connectance,
                   individual_variation=params_forcing$individual_variation,
                   forcing_strength= params_forcing$forcing_strength,
                   mut_strength=params_forcing$mut.strength[1],
                   forcing_duration = params_forcing$duration,
                   interaction_type=params_forcing$interaction_type)
  
  ## return equations by first flattening them back into a single vector
  output=out 
  #list(Plants = Np[1:time,],Animals=Na[1:time,], Plant.trait = muP[1:time,], 
  #          Animal.trait=muA[1:time,],s=pars$s,aij=aij,ba=pars$ba)
  
  
  
  
  ddf <-as.data.frame(cbind(
    rep(seq(1,(nrow(params_forcing$matrix)+ncol(params_forcing$matrix)),1)),
    params_forcing$forcing_strength,
    as.numeric(c( (na %>% filter(time == tmax))[3:A],(np %>% filter(time == tmax))[3:P])),
    params_forcing$mut.strength[1],
    rr$Response_time,
    c(rr$Degree),
    rr$response,
    rep(network_response_time,each=(A+P)),
    rep(nestedness_NODF(params_forcing$matrix), each=((A+P)) ),
    rep(Connectance(params_forcing$matrix), each=((A+P)) ),
    rep( (A+P),each=((A+P))),
    rep(as.character(params_forcing$individual_variation),each=(A+P)),
    interaction_type = params_forcing$interaction_type))
  
  
  colnames(ddf)<-c("Species",'forcing_strength',"abundance","mutualism_strength", 
                   "response_time","Degree","response_thr", "Network_respose_time",
                   "Nestedness", "Connectance","Network_size",
                   "Individual_variation", "interaction_type")
  
  
  final_output<- list(output=output, ddf=ddf, sol=sol)
  return(final_output)
  
}


## Organize simulation results into tidy table (code adapted from Barabas and D'Andrea 2016 Eco.Letts paper)
## Input:
## - sol: output produced by the function ode()
## - pars: list of parameters, with the following elements:
## Output:
## - a tibble with columns: time, species, n (density), m (trait mean),
##   sigma
organize_results <- function(sol, pars) {
  S <- length(pars$sigma) ## number of species
  A<-dim(pars$matrix)[2] # no. of animals
  P<-dim(pars$matrix)[1] # no. of plants
  temp<- sol %>% as.data.frame %>% as_tibble ## convert to tibble
  ## name the first column "time"
  # temp<- temp %>% filter(time >= pars$cutoff.time)
  names(temp)[2:(A+1)] <- paste0("N_", 1:(A)) ## name abundance columns (n_k)
  names(temp)[1] <- "time"
  names(temp)[(A+1+1):(A+P+1)] <- paste0("P_",  A+1:P) ## name trait mean columns
  names(temp)[(A+P+2):(A+P+A+1)] <- paste0("ma_",1:A)
  names(temp)[(A+P+A+2):(A+P+A+P+1)]<-paste0("mp_",  A+1:P)
  temp <- temp %>%
    tidyr:: gather("variable", "v", 2:ncol(temp)) %>% ## normalize the data
    tidyr::separate(variable, c("type", "species"), sep="_") %>%
    # spread(type, v) %>% ## separate columns for animal densities n and plant densities m
    dplyr::select(time, type, species, v) %>% ## rearrange columns
    mutate(species=as.integer(species), Temp=pars$Temp) ## add params
  return(as_tibble(temp))
}


## Plot time series of densities, time series of trait values, and
## snapshot of the trait distributions at time = moment (code adapted from Barabas and D'Andrea 2016 Eco.Letts paper)
## Input:
## - dat: data generated by organize_results()
## - moment: time at which trait distribution should be plotted
## - limits: a vector of two entries (x_low, x_high) for the x-axis limits
## - res: number of evenly spaced sampling points along the trait axis
##               for the trait distribution plot
## Output:
## - a ggplot2 plot with three panels in one column: abundance time series,
##   trait value time seties, and snapshot of trait distribution
plot_all <- function(dat, moment=0, limits=c(-0.5, 0.5), res=1001) {
  plot_grid(plot_density(dat), ncol=1, align="hv") %>%
    return
}



## Plot species densities through time
## Input:
## - dat: data generated by organize_results()
## Output:
## - a ggplot2 plot
## used to produce figure 1.

## code adapted from Barabas and D'Andrea 2016 Eco.Letts paper
plot_density<- function(dat) {
  dat %>% filter(type != "ma", type != "mp") %>% 
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species)),size=3) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    scale_color_viridis_d()+
    theme(legend.position="none") %>%  #+ facet_wrap(.~type, scales = "free",) %>%
    return
}



#multiplot of ggplot2 figures with a common shared legend. Code taken from :https://rpubs.com/sjackman/grid_arrange_shared_legend
grid_arrange_shared_legend <- function(..., ncol, nrow, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol =ncol , nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}


# plots the density distribution  at a particular timepoint.
# Na: Abundance of animals at equilibrium
# Np: Abundance of plants at equilibrium
# m: mean traits at equilibrium
# sigma: variance of traits
# moment: mean
# limits: limits of the mean trait axis which in the study are -1,1

plot_snapshot <- function(Na, Np, m, sigma, moment=0, limits=c(-1, 1), res=1001) {
  Sa <- length(Na) ## number of species
  Sp <- length(Np)
  ma<- m[1:(Sa)]
  mp<- m[Sa+1:Sp]
  sigma_a <-sigma[1:(Sa)]
  sigma_p <- sigma[Sa+1:Sp]
  traitaxis <- seq(limits[1], limits[2], l=res) ## sampling the trait axis
  #snap <- dat %>% filter(time==moment) %>% select(-time) ## time = moment
  traits_a <- expand.grid(species=1:Sa, trait=traitaxis) %>% as_tibble ## trait table
  traits_p <- expand.grid(species=Sa+1:Sp, trait=traitaxis) %>% as_tibble ## trait table
  
  traits_a["density"] <- 0 ## add column for population densities
  traits_p["density"] <- 0
  
  for (i in 1:Sa) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_a$density[(traits_a$species==i)] <- Na[i]*
      dnorm(traits_a$trait[(traits_a$species==i)], ma[i], sigma_a[i]) ## times density
  }
  traits_a$density[traits_a$density<max(traits_a$density)/1e3] <- NA
  
  for (i in 1:Sp) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_p$density[(traits_p$species==(Sa+i))] <- Np[i]*dnorm(traits_p$trait[(traits_p$species==(Sa+i))], 
                                                                mp[i], sigma_p[i]) ## times density
  }
  traits_p$density[traits_p$density<max(traits_p$density)/1e3] <- NA
  
  
  traits<-data.frame(rbind(traits_a,traits_p), 
                     species_group=c(rep("Animals", nrow(traits_a)),
                                     rep("Plants", nrow(traits_p))))
  
  ggplot(traits) + ## generate plot
    geom_line(aes(x=trait, y=density, colour=factor(species)), na.rm=TRUE) +
    geom_ribbon(aes(x=trait, ymin=0, ymax=density, fill=factor(species)),
                alpha=0.15, colour=NA)+scale_fill_viridis_d()+
    facet_wrap(.~species_group, nrow = 2)+
    theme(legend.title = element_text(size = 14, face = "bold"), 
          legend.position = "right", panel.background = element_blank(), 
          axis.text = element_text(colour = "black", size = 14, face = "bold"), 
          axis.title = element_text(size = 14, face = "bold"), 
          legend.text = element_text(size = 14), legend.key = element_blank(),
          strip.text.x = element_text(size= 14, face ="bold"))+
    #geom_line(data=landscape, aes(x=trait, y=r), linetype="dashed",
    #         colour="darkred", alpha=0.5, na.rm=TRUE) +
    scale_x_continuous(name="trait value", limits=limits) +
    scale_y_continuous(name="density", limits=c(0, NA)) +
    theme(legend.position="none") %>%
    return }




#computes the raw NODF taken from Song et al 2017 J. Animal Ecology
#input: web = mutualistic network
#output: raw NODF of the given network
nestedness_NODF <- function(web){
  web[web > 0] = 1
  SA <- nrow(web)
  SP <- ncol(web)
  N <- t(web) %*% web
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SP,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SP)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele == 0] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)] = 1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n1 <- sum(nes)
  
  N <- web %*% t(web)
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SA,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SA)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele ==0 ] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)]=1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n2 <- sum(nes)
  out <- 2*(n1 + n2) / (SA*(SA-1)+SP*(SP-1))
  return(out)
}



# measures connectance of a web network
# web: interaction network
Connectance<-function(web)
{
  return(sum(web)/(ncol(web)*nrow(web)))}

# function for sampling competitive coefficients from random uniform distribution 
# matrix: network of interactions which are 0 or 1. 
# strength: average competition strength
mat.comp<-function(matrix,degree.animals,degree.plants){
  Aspecies<- dim(matrix)[2]
  Plantspecies<- dim(matrix)[1]
  
  Amatrix<-matrix(runif(Aspecies^2, 0.0001, 0.0005), nrow=Aspecies, ncol = Aspecies)
  diag(Amatrix)<-1
  #diag(Amatrix)<-  diag(Amatrix) #/degree.animals
  
  Pmatrix<-matrix(runif(Plantspecies^2, 0.0001, 0.0005), nrow=Plantspecies, ncol = Plantspecies)
  diag(Pmatrix)<-1
  #diag(Pmatrix)<-2-  diag(Pmatrix)/degree.plants
  out<-return(list(Amatrix=Amatrix,Pmatrix=Pmatrix))
  
}




# function for sampling competitive coefficients from exponential distribution 
# competitive interactions  were scaled by the total number of species within a guild as Dakos & Bascompte 2014 PNAS.
# matrix: network of interactions which are 0 or 1. 
# strength: average competition strength

mat.comp_exp<-function(matrix){
  Aspecies<- dim(matrix)[2]
  Plantspecies<- dim(matrix)[1]
  
  Amatrix<- matrix(rexp(Aspecies^2, 70),nrow=Aspecies,ncol=Aspecies)/Aspecies #scaling by no. of competitors within a guild
  
  diag(Amatrix)<-1
  Pmatrix<-matrix(rexp(Plantspecies^2, 70),nrow=Plantspecies,ncol=Plantspecies)/Plantspecies  #scaling by no. of competitors within a guild
  
  diag(Pmatrix)<-1
  
  out<-return(list(Amatrix=Amatrix,Pmatrix=Pmatrix))
  
}

# function for competition coefficients within a guild for feasibility analysis.
# competitive interactions  were scaled by the total number of species within a guild as Dakos & Bascompte 2014 PNAS.
# matrix: network of interactions which are 0 or 1. 
# strength: average competition strength
mat.comp_feasibility<-function(matrix,strength){
  Aspecies<- dim(matrix)[2]
  Plantspecies<- dim(matrix)[1]
  
  Amatrix<-strength*matrix(runif(Aspecies^2, 0.0001, 0.001), nrow=Aspecies, ncol = Aspecies)/Aspecies #scaled by number of competitors within a guild
  diag(Amatrix)<-1 #intraspecific competition for animals
  Pmatrix<-strength*matrix(runif(Plantspecies^2, 0.0001, 0.001), nrow=Plantspecies, ncol = Plantspecies)/Plantspecies ##scaled by number of competitors within a guild
  diag(Pmatrix)<-1 #intraspecific competion for plants
  
  out<-return(list(Amatrix=Amatrix,Pmatrix=Pmatrix))
  
}


#plot for feasibility analysis that takes data frame as input
feasibility_plot<-function(dat){
  akima.R<-with(dat,interp(dat$Strength_mutualism,log(dat$range_competition),dat$mean_feasibility,
                           xo=seq(min(dat$Strength_mutualism),max(dat$Strength_mutualism),length=20), 
                           yo=seq(min(log(dat$range_competition)),max(log(dat$range_competition)), length=20)))
  
  gdat1<-interp2xyz(akima.R, data.frame=TRUE)
  colnames(gdat1)<-c("mutualism_strength","range_competition","Feasibility")
  a<-ggplot(gdat1 , aes(x=mutualism_strength,y=range_competition,
                        z=Feasibility))+
    geom_raster(aes(fill=Feasibility),show.legend =T)+ 
    #geom_contour(aes(colour = ..level..),bins=1)+
    theme_cowplot()+
    theme(legend.title = element_text(size = 9, face = "bold"), 
          legend.text=element_text(size=rel(0.5)),
          legend.position = "bottom", panel.background = element_blank(), 
          axis.text = element_text(colour = "black", size = 9, face = "bold"), 
          axis.title = element_text(size = 9, face = "bold"), 
          legend.key = element_blank())+
    scale_fill_continuous(low = "#BFE1B0", high = "#137177") +
    labs(x = expression(paste("Mutualistic strength, ",gamma[0])),
         y = expression(paste(log(rho)))) +
    scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))
  
  a
  
  
  return(a)
}

####### MISCELLANEOUS FUNCTIONS not used in the manuscript ####################


#response time for each species in a network

#Na: timeseries of animals
#Np: time series of plants
#Aspecies: no of animals
#Plantspecies: no of plants
#g:adjacency matrix
#t: time

response_time<-function(sol, Aspecies, Plantspecies, g, t , dat, webname){
  
  thr_time_Aspecies<-thr_time_Plants<-degree.plants<-degree.animals<-thr_response_Aspecies<-thr_response_Pspecies<-numeric()
  
  Na = (sol %>% filter(type !="ma",type !="mp", type =="N") %>% select(time,v,species))
  Np = (sol %>% filter(type !="ma",type !="mp", type =="P") %>% select(time,v,species))
  
  
  for(i in 1:Plantspecies){
    degree.plants[i]<-sum(g[i,])} # degree of plants
  for(j in 1:Aspecies){
    degree.animals[j]<-sum(g[,j]) # degree of animals
  }
  
  
  for(i in 1:Aspecies) { 
    
    indexA<-which((Na %>% filter(species == i))$v >= 0.5)
    #0.25*upper_state)
    
    if(length(indexA) == 0){
      indexA <- NA
      thr_response_A<- 0
      
    }else { indexA  <- min(indexA)
    thr_response_A<- 1
    }
    
    thr_time_Aspecies[i]<- indexA
    thr_response_Aspecies[i]<- thr_response_A
  }
  
  for (j in 1:Plantspecies){
    indexP<- which((Np %>% filter(species == i))$v >= 0.5)
    
    if(length(indexP) == 0){
      indexP <- NA
      thr_response_P<-0
    }else { indexP  <- min(indexP)
    thr_response_P<-1
    
    }
    
    thr_time_Plants[j]<-indexP
    thr_response_Pspecies[j]<-thr_response_P
    
  }
  
  
  
  return(list(Response_time=c(thr_time_Aspecies,thr_time_Plants),
              response=c(thr_response_Aspecies,thr_response_Pspecies),
              Degree=c(degree.animals,degree.plants)))
  
}



# final state after perturbation
final_state<-function(forcing_time ){
  
  
  for(i in 1:Plantspecies){
    
    du<-exp(-(0.5)^2/(2*si[i]^2 + 2*sj^2 + w^2))
    tu[i]<- bi[i] - 2*(degree.plants[i]+0.1)+  sum(g[i,]*0.1*0.2*(1/sqrt(2*si[i]^2 + 2*sj^2 + w^2))*du)
    tau[i]<- - 1/tu[i]
    
  }
  
  
  for(j in 1:Aspecies){
    
    dua<-exp(-(0.5)^2/(2*si^2 + 2*sj[j]^2 + w^2))
    tua[j]<- bj[j] - 2*(degree.animals[j]+0.1)+  sum(g[,j]*0.1*0.2*(1/sqrt(2*si^2 + 2*sj[j]^2 + w^2))*dua)
    tau_a[j]<- - 1/tua[j]
    
  }
  
  
  
  
}

