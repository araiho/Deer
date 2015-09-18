var N[4,T,Park], N.3yr[4,T,Park], N.1yr[4,T,Park], N.ster[4,T,Park], N.cull[4,T,Park],  mu[4,T,Park], mu.3yr[4,T,Park], mu.1yr[4,T,Park], mu.ster[4,T,Park], mu.cull[4,T,Park]

model{
  
  #####
  ##### Priors ##### 
  #####
  
  a ~ dgamma(mu.a^2/sd.a^2, mu.a/sd.a^2)
  b ~ dgamma(mu.b^2/sd.b^2, mu.b/sd.b^2)
  
  mean.adult.female<- a/(a+b)
  
  a1 ~ dgamma(mu.a1^2/sd.a1^2, mu.a1/sd.a1^2)
  b1 ~ dgamma(mu.b1^2/sd.b1^2, mu.b1/sd.b1^2)
  
  mean.fawn <- a1/(a1+b1)
  
  a2 ~ dgamma(mu.a2^2/sd.a2^2, mu.a2/sd.a2^2)
  b2 ~ dgamma(mu.b2^2/sd.b2^2, mu.b2/sd.b2^2)
  
  mean.adult.male <- a2/(a2+b2)
  W<-65
  
  r ~ dnorm(2*3.09*W^-.33,1/((.1304)^2))
  
  beta ~ dunif(0,100)
  
  ratio ~ dbeta(312,312) # sex ratio # females:males
  
  sigma.p ~ dunif(0,2)
  tau.p <- 1/(sigma.p)^2
  
  #####
  ##### Initial Conditions #####
  #####
  
  for(p in 1:Park){
    
    s.adult.female[p] ~ dbeta(a,b)T(0,.99)
    s.fawn[p] ~ dbeta(a1,b1)T(0,.99)
    s.adult.male[p] ~ dbeta(a2,b2)T(0,.99)
    
    prop[1:3,p] ~ ddirch(N.init[1:3,1,p]+1)
    number.park[p] ~ dnorm(N.obs.init[p,5],1/(N.obs.init[p,8])^2)
    N[1:3,1,p] <- prop[1:3,p]*number.park[p]*area[p]
    sumNpark[1,p]<-sum(N[1:3,1,p])
    sumNpark.3yr[1,p]<-sum(N[1:3,1,p])
    sumNpark.1yr[1,p]<-sum(N[1:3,1,p])
    sumNpark.ster[1,p]<-sum(N[1:3,1,p])
    sumNpark.cull[1,p]<-sum(N[1:3,1,p])
    
    denNpark[1,p]<-sumNpark[1,p]/area[p]
    denNpark.3yr[1,p]<-sumNpark[1,p]/area[p]
    denNpark.1yr[1,p]<-sumNpark[1,p]/area[p]
    denNpark.ster[1,p]<-sumNpark[1,p]/area[p]
    denNpark.cull[1,p]<-sumNpark[1,p]/area[p]
    
    N.cull[1:3,1,p] <- N[1:3,1,p]
    N.ster[1:3,1,p] <- N[1:3,1,p]
    N.1yr[1:3,1,p] <- N[1:3,1,p]
    N.3yr[1:3,1,p] <- N[1:3,1,p]
    
  }
  
  sumN[1,1] <-sum(sumNpark[1,1:Park])
  sumN.3yr[1,1] <-sum(sumNpark[1,1:Park])
  sumN.1yr[1,1] <-sum(sumNpark[1,1:Park])
  sumN.ster[1,1] <-sum(sumNpark[1,1:Park])
  sumN.cull[1,1] <-sum(sumNpark[1,1:Park])
  
  sumNstage[1,1] <- sum(N[1,1,1:Park])
  sumNstage[2,1] <- sum(N[2,1,1:Park])
  sumNstage[3,1] <- sum(N[3,1,1:Park])
  
  for(p in 1:Park){
    
    f.adult[1,p] <- min(5,exp(r-(r/beta)*denNpark[1,p]))
    
    M[1,2,p,1]<-s.adult.female[p]*f.adult[1,p]
    M[2,1,p,1]<-s.fawn[p]*ratio
    M[2,2,p,1]<-s.adult.female[p]
    M[3,1,p,1]<-s.fawn[p]*(1-ratio)
    M[3,3,p,1]<-s.adult.male[p]
    
  }
  
  #####
  ##### Process Model #####
  #####
  
  for (i in 1:length(N.obsno1s[,1])-40){
    f.adult[N.obsno1s[i,4],N.obsno1s[i,3]] <- min(5,exp(r-(r/beta)*denNpark[N.obsno1s[i,4]-1,N.obsno1s[i,3]]))
    
    M[1,2,N.obsno1s[i,3],N.obsno1s[i,4]]<-s.adult.female[N.obsno1s[i,3]]*f.adult[N.obsno1s[i,4],N.obsno1s[i,3]]
    M[2,1,N.obsno1s[i,3],N.obsno1s[i,4]]<-s.fawn[N.obsno1s[i,3]]*ratio
    M[2,2,N.obsno1s[i,3],N.obsno1s[i,4]]<-s.adult.female[N.obsno1s[i,3]]
    M[3,1,N.obsno1s[i,3],N.obsno1s[i,4]]<-s.fawn[N.obsno1s[i,3]]*(1-ratio)
    M[3,3,N.obsno1s[i,3],N.obsno1s[i,4]]<-s.adult.male[N.obsno1s[i,3]]
    
    
    mu[,N.obsno1s[i,4],N.obsno1s[i,3]] <- M[,,N.obsno1s[i,3],N.obsno1s[i,4]]%*%N[,N.obsno1s[i,4]-1,N.obsno1s[i,3]]
    
    N[1,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu[1,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    N[2,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu[2,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    N[3,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu[3,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)	
    
    sumNpark[N.obsno1s[i,4],N.obsno1s[i,3]] <- sum(N[1:4,N.obsno1s[i,4],N.obsno1s[i,3]])
    sumNpark.3yr[N.obsno1s[i,4],N.obsno1s[i,3]] <- sum(N[1:4,N.obsno1s[i,4],N.obsno1s[i,3]])
    sumNpark.1yr[N.obsno1s[i,4],N.obsno1s[i,3]] <- sum(N[1:4,N.obsno1s[i,4],N.obsno1s[i,3]])
    sumNpark.ster[N.obsno1s[i,4],N.obsno1s[i,3]] <- sum(N[1:4,N.obsno1s[i,4],N.obsno1s[i,3]])
    sumNpark.cull[N.obsno1s[i,4],N.obsno1s[i,3]] <- sum(N[1:4,N.obsno1s[i,4],N.obsno1s[i,3]])
    
    denNpark[N.obsno1s[i,4],N.obsno1s[i,3]] <- sumNpark[N.obsno1s[i,4],N.obsno1s[i,3]]/area[N.obsno1s[i,3]]
    denNpark.3yr[N.obsno1s[i,4],N.obsno1s[i,3]] <- sumNpark[N.obsno1s[i,4],N.obsno1s[i,3]]/area[N.obsno1s[i,3]]
    denNpark.1yr[N.obsno1s[i,4],N.obsno1s[i,3]] <- sumNpark[N.obsno1s[i,4],N.obsno1s[i,3]]/area[N.obsno1s[i,3]]
    denNpark.ster[N.obsno1s[i,4],N.obsno1s[i,3]] <- sumNpark[N.obsno1s[i,4],N.obsno1s[i,3]]/area[N.obsno1s[i,3]]
    denNpark.cull[N.obsno1s[i,4],N.obsno1s[i,3]] <- sumNpark[N.obsno1s[i,4],N.obsno1s[i,3]]/area[N.obsno1s[i,3]]
    
    N.3yr[1:3,N.obsno1s[i,4],N.obsno1s[i,3]] <- N[1:3,N.obsno1s[i,4],N.obsno1s[i,3]]
    N.1yr[1:3,N.obsno1s[i,4],N.obsno1s[i,3]] <- N[1:3,N.obsno1s[i,4],N.obsno1s[i,3]]
    N.ster[1:3,N.obsno1s[i,4],N.obsno1s[i,3]] <- N[1:3,N.obsno1s[i,4],N.obsno1s[i,3]]
    N.cull[1:3,N.obsno1s[i,4],N.obsno1s[i,3]] <- N[1:3,N.obsno1s[i,4],N.obsno1s[i,3]] 
    
  }
  
  #####		
  ##### Forecasting and Model Experiments #####		
  #####	
  
  treated.fawns <- 0
  
  effect <- 1-(1/3)*exp((-1/3)*1) 
  
  for (i in length(N.obsno1s[,1])-39:length(N.obsno1s[,1])){
    
    #####
    ##### Average of 3 Year Effectiveness #####
    #####
    
    f.adult.3yr[N.obsno1s[i,4],N.obsno1s[i,3]] <- min(5,exp(r-(r/beta)*denNpark.3yr[N.obsno1s[i,4]-1,N.obsno1s[i,3]]))
    
    M.3yr[1,2,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*f.adult.3yr[N.obsno1s[i,4],N.obsno1s[i,3]]
    M.3yr[2,1,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.fawn[N.obsno1s[i,3]]*ratio*(1-treated.fawns)
    M.3yr[2,2,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*(1-treated.does)
    M.3yr[3,1,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.fawn[N.obsno1s[i,3]]*(1-ratio)
    M.3yr[3,3,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.male[N.obsno1s[i,3]]
    
    M.3yr[4,1,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.fawn[N.obsno1s[i,3]]*ratio*treated.fawns
    M.3yr[4,2,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*treated.does
    M.3yr[2,4,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*(1-effect)*(1-treated.does)
    M.3yr[4,4,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*(treated.does)+s.adult.female[N.obsno1s[i,3]]*(1-treated.does)*effect
    
    
    mu.3yr[,N.obsno1s[i,4],N.obsno1s[i,3]] <- M.3yr[,,N.obsno1s[i,3],N.obsno1s[i,4]]%*%N.3yr[,N.obsno1s[i,4]-1,N.obsno1s[i,3]]
    
    N.3yr[1,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.3yr[1,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    N.3yr[2,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.3yr[2,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    N.3yr[3,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.3yr[3,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    N.3yr[4,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.3yr[4,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    
    sumNpark.3yr[N.obsno1s[i,4],N.obsno1s[i,3]] <- sum(N.3yr[1:4,N.obsno1s[i,4],N.obsno1s[i,3]])
    denNpark.3yr[N.obsno1s[i,4],N.obsno1s[i,3]] <- sumNpark.3yr[N.obsno1s[i,4],N.obsno1s[i,3]]/area[N.obsno1s[i,3]]
    
    #####
    ##### Average of 1 Year Effectiveness #####
    #####
    
    f.adult.1yr[N.obsno1s[i,4],N.obsno1s[i,3]] <- min(5,exp(r-(r/beta)*denNpark.1yr[N.obsno1s[i,4]-1,N.obsno1s[i,3]]))
    
    M.1yr[1,2,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*f.adult.1yr[N.obsno1s[i,4],N.obsno1s[i,3]]
    M.1yr[2,1,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.fawn[N.obsno1s[i,3]]*ratio*(1-treated.fawns)
    M.1yr[2,2,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*(1-treated.does)
    M.1yr[3,1,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.fawn[N.obsno1s[i,3]]*(1-ratio)
    M.1yr[3,3,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.male[N.obsno1s[i,3]]
    
    M.1yr[2,4,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*(1-treated.does)
    M.1yr[4,1,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.fawn[N.obsno1s[i,3]]*ratio*treated.fawns
    M.1yr[4,2,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*treated.does
    M.1yr[4,4,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*treated.does
    
    
    mu.1yr[,N.obsno1s[i,4],N.obsno1s[i,3]] <- M.1yr[,,N.obsno1s[i,3],N.obsno1s[i,4]]%*%N.1yr[,N.obsno1s[i,4]-1,N.obsno1s[i,3]]
    
    N.1yr[1,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.1yr[1,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    N.1yr[2,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.1yr[2,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    N.1yr[3,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.1yr[3,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    N.1yr[4,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.1yr[4,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    
    sumNpark.1yr[N.obsno1s[i,4],N.obsno1s[i,3]] <- sum(N.1yr[1:4,N.obsno1s[i,4],N.obsno1s[i,3]])
    denNpark.1yr[N.obsno1s[i,4],N.obsno1s[i,3]] <- sumNpark.1yr[N.obsno1s[i,4],N.obsno1s[i,3]]/area[N.obsno1s[i,3]]
    
    #####
    ##### Sterilization #####
    #####
    
    f.adult.ster[N.obsno1s[i,4],N.obsno1s[i,3]] <- min(5,exp(r-(r/beta)*denNpark.ster[N.obsno1s[i,4]-1,N.obsno1s[i,3]]))
    
    M.ster[1,2,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*f.adult.ster[N.obsno1s[i,4],N.obsno1s[i,3]]*(1-treated.does)
    M.ster[2,1,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.fawn[N.obsno1s[i,3]]*ratio*(1-treated.fawns)
    M.ster[2,2,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*(1-treated.does)
    M.ster[3,1,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.fawn[N.obsno1s[i,3]]*(1-ratio)
    M.ster[3,3,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.male[N.obsno1s[i,3]]
    
    M.ster[4,1,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.fawn[N.obsno1s[i,3]]*ratio*treated.fawns 
    M.ster[4,2,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*treated.does
    M.ster[4,4,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]   
    
    mu.ster[,N.obsno1s[i,4],N.obsno1s[i,3]] <- M.ster[,,N.obsno1s[i,3],N.obsno1s[i,4]]%*%N.ster[,N.obsno1s[i,4]-1,N.obsno1s[i,3]]
    
    N.ster[1,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.ster[1,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    N.ster[2,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.ster[2,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    N.ster[3,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.ster[3,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    N.ster[4,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.ster[4,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    
    sumNpark.ster[N.obsno1s[i,4],N.obsno1s[i,3]] <- sum(N.ster[1:4,N.obsno1s[i,4],N.obsno1s[i,3]])
    denNpark.ster[N.obsno1s[i,4],N.obsno1s[i,3]] <- sumNpark.ster[N.obsno1s[i,4],N.obsno1s[i,3]]/area[N.obsno1s[i,3]]
    
    #####
    ##### Culling #####
    #####
    
    f.adult.cull[N.obsno1s[i,4],N.obsno1s[i,3]] <- min(5,exp(r-(r/beta)*denNpark.cull[N.obsno1s[i,4]-1,N.obsno1s[i,3]]))
    
    M.cull[1,2,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*f.adult.cull[N.obsno1s[i,4],N.obsno1s[i,3]]*(1-treated.does)
    M.cull[2,1,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.fawn[N.obsno1s[i,3]]*ratio*(1-treated.fawns)
    M.cull[2,2,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female[N.obsno1s[i,3]]*(1-treated.does)
    M.cull[3,1,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.fawn[N.obsno1s[i,3]]*(1-ratio)
    M.cull[3,3,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.male[N.obsno1s[i,3]]
    
    mu.cull[,N.obsno1s[i,4],N.obsno1s[i,3]] <- M.cull[,,N.obsno1s[i,3],N.obsno1s[i,4]]%*%N.cull[,N.obsno1s[i,4]-1,N.obsno1s[i,3]]
    
    N.cull[1,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.cull[1,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    N.cull[2,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.cull[2,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    N.cull[3,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu.cull[3,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
    
    sumNpark.cull[N.obsno1s[i,4],N.obsno1s[i,3]] <- sum(N.cull[1:3,N.obsno1s[i,4],N.obsno1s[i,3]])
    denNpark.cull[N.obsno1s[i,4],N.obsno1s[i,3]] <- sumNpark.cull[N.obsno1s[i,4],N.obsno1s[i,3]]/area[N.obsno1s[i,3]]
    
  }
  
  for(i in 2:T){
    
    sumN.3yr[1,i] <-mean(denNpark.3yr[i,])
    sumN.1yr[1,i] <-mean(denNpark.1yr[i,])
    sumN.ster[1,i] <-mean(denNpark.ster[i,])
    sumN.cull[1,i] <-mean(denNpark.cull[i,])
    
  }
  
  
  #####
  ##### Data Model #####
  #####
  
  for (i in 1:length(N.obs[,1])){
    
    density[N.obs[i,4],N.obs[i,3]] <- sumNpark[N.obs[i,4],N.obs[i,3]]/area[N.obs[i,3]]
    N.obs[i,5] ~ dnorm(density[N.obs[i,4],N.obs[i,3]],1/(N.obs[i,8])^2)T(0,500)
    N.obs.new[i] ~ dnorm(density[N.obs[i,4],N.obs[i,3]],1/(N.obs[i,8])^2)T(0,500)
    c[1,N.obs[i,4],N.obs[i,3]] <- N[1,N.obs[i,4],N.obs[i,3]]/sum(N[,N.obs[i,4],N.obs[i,3]])
    c[2,N.obs[i,4],N.obs[i,3]] <- N[2,N.obs[i,4],N.obs[i,3]]/sum(N[,N.obs[i,4],N.obs[i,3]])
    c[3,N.obs[i,4],N.obs[i,3]] <- N[3,N.obs[i,4],N.obs[i,3]]/sum(N[,N.obs[i,4],N.obs[i,3]])
    
  }
  
  for(i in 1:length(y.alpha[,1])){
    
    y.alpha[i,4:6] ~ dmulti(c[1:3,y.alpha[i,2],y.alpha[i,1]],y.alpha[i,8])
    
  }
  
  #####
  ###### Posterior Predictive Checks #####
  ######
  
  for(i in 1:length(N.obs[,1])){
    sq[i] <- (N.obs[i,5] - sumNpark[N.obs[i,4],N.obs[i,3]] / area[N.obs[i,3]])^2
    sq.new[i] <- (N.obs.new[i] - sumNpark[N.obs[i,4],N.obs[i,3]] / area[N.obs[i,3]])^2
  }
  
  fit <- sum(sq[])
  fit.new <- sum(sq.new[])
  pvalue.fit <- step(fit.new-fit)
  
  f.adult.lambda <- min(5,exp(r-(r/beta)*0))/2
  
  M.lambda.cull[1,2] <- mean.adult.female*f.adult.lambda*(1-treated.does)
  M.lambda.cull[2,1] <- mean.fawn*ratio
  M.lambda.cull[2,2] <- mean.adult.female*(1-treated.does)
  
  M.lambda.ster[1,2] <- mean.adult.female*f.adult.lambda*(1-treated.does)
  M.lambda.ster[2,1] <- mean.fawn*ratio
  M.lambda.ster[2,2] <- mean.adult.female*(1-treated.does)
  M.lambda.ster[3,2] <- mean.adult.female*treated.does
  M.lambda.ster[3,3] <- mean.adult.female
  
  M.lambda.1yr[1,2] <- mean.adult.female*f.adult.lambda
  M.lambda.1yr[2,1] <- mean.fawn*ratio
  M.lambda.1yr[2,2] <- mean.adult.female*(1-treated.does)
  M.lambda.1yr[2,3] <- mean.adult.female*(1-treated.does)
  M.lambda.1yr[3,2] <- mean.adult.female*treated.does
  M.lambda.1yr[3,3] <- mean.adult.female*treated.does
  
  M.lambda.3yr[1,2] <- mean.adult.female*f.adult.lambda
  M.lambda.3yr[2,1] <- mean.fawn*ratio
  M.lambda.3yr[2,2] <- mean.adult.female*(1-treated.does)
  M.lambda.3yr[2,3] <- mean.adult.female*(1-treated.does)*(1-effect)
  M.lambda.3yr[3,2] <- mean.adult.female*treated.does
  M.lambda.3yr[3,3] <- mean.adult.female*treated.does+mean.adult.female*(1-treated.does)*effect  
  
}
