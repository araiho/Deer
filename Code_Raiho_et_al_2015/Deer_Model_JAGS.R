var N[4,T,Park],  mu[4,T,Park]

model{
	
#####
##### Priors ##### 
#####
  	 
  s.fawn ~ dbeta(1,1)T(0,.99)
  s.adult.female ~ dbeta(1,1)T(0,.99)
  s.adult.male ~ dbeta(1,1)T(0,.99)
    
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
    
    prop[1:3,p] ~ ddirch(N.init[1:3,1,p]+1)
    number.park[p] ~ dnorm(N.obs.init[p,5],1/(N.obs.init[p,8])^2)
    N[1:3,1,p] <- prop[1:3,p]*number.park[p]*area[p]
    sumNpark[1,p]<-sum(N[1:3,1,p])
    denNpark[1,p]<-sumNpark[1,p]/area[p]
    }
    
    sumN[1,1] <-mean(denNpark[1,1:Park])
    
    sumNstage[1,1] <- sum(N[1,1,1:Park])
    sumNstage[2,1] <- sum(N[2,1,1:Park])
    sumNstage[3,1] <- sum(N[3,1,1:Park])

  for(p in 1:Park){

    f.adult[1,p] <- min(5,exp(r-(r/beta)*denNpark[1,p]))
 
    M[1,2,p,1]<-s.adult.female*f.adult[1,p]
    M[2,1,p,1]<-s.fawn*ratio
    M[2,2,p,1]<-s.adult.female
    M[3,1,p,1]<-s.fawn*(1-ratio)
    M[3,3,p,1]<-s.adult.male
  
}

#####
##### Process Model #####
#####

  for (i in 1:length(N.obsno1s[,1])-40){
    f.adult[N.obsno1s[i,4],N.obsno1s[i,3]] <- min(5,exp(r-(r/beta)*denNpark[N.obsno1s[i,4]-1,N.obsno1s[i,3]]))

    M[1,2,N.obsno1s[i,3],N.obsno1s[i,4]]<-s.adult.female*f.adult[N.obsno1s[i,4],N.obsno1s[i,3]]
    M[2,1,N.obsno1s[i,3],N.obsno1s[i,4]]<-s.fawn*ratio
    M[2,2,N.obsno1s[i,3],N.obsno1s[i,4]]<-s.adult.female
    M[3,1,N.obsno1s[i,3],N.obsno1s[i,4]]<-s.fawn*(1-ratio)
    M[3,3,N.obsno1s[i,3],N.obsno1s[i,4]]<-s.adult.male

			
	mu[,N.obsno1s[i,4],N.obsno1s[i,3]] <- M[,,N.obsno1s[i,3],N.obsno1s[i,4]]%*%N[,N.obsno1s[i,4]-1,N.obsno1s[i,3]]
	
	N[1,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu[1,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
	N[2,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu[2,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
	N[3,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu[3,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)

	sumNpark[N.obsno1s[i,4],N.obsno1s[i,3]] <- sum(N[1:4,N.obsno1s[i,4],N.obsno1s[i,3]])
	denNpark[N.obsno1s[i,4],N.obsno1s[i,3]] <- sumNpark[N.obsno1s[i,4],N.obsno1s[i,3]]/area[N.obsno1s[i,3]]

		}
		
#####		
##### Forecasting and Model Experiments #####		
#####
	
  for (i in length(N.obsno1s[,1])-39:length(N.obsno1s[,1])){
  	
  	#####
  	##### No Action #####
  	#####
  	
    f.adult[N.obsno1s[i,4],N.obsno1s[i,3]] <- min(5,exp(r-(r/beta)*denNpark[N.obsno1s[i,4]-1,N.obsno1s[i,3]]))
  
    M[1,2,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female*f.adult[N.obsno1s[i,4],N.obsno1s[i,3]]
    M[2,1,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.fawn*ratio
    M[2,2,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.female
    M[3,1,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.fawn*(1-ratio)
    M[3,3,N.obsno1s[i,3],N.obsno1s[i,4]] <- s.adult.male
    
    mu[,N.obsno1s[i,4],N.obsno1s[i,3]] <- M[,,N.obsno1s[i,3],N.obsno1s[i,4]]%*%N[,N.obsno1s[i,4]-1,N.obsno1s[i,3]]
	
	N[1,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu[1,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
	N[2,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu[2,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)
	N[3,N.obsno1s[i,4],N.obsno1s[i,3]] ~ dlnorm(log(max(1,mu[3,N.obsno1s[i,4],N.obsno1s[i,3]])),tau.p)

	sumNpark[N.obsno1s[i,4],N.obsno1s[i,3]] <- sum(N[1:4,N.obsno1s[i,4],N.obsno1s[i,3]])
	 denNpark[N.obsno1s[i,4],N.obsno1s[i,3]] <- sumNpark[N.obsno1s[i,4],N.obsno1s[i,3]]/area[N.obsno1s[i,3]]
    
			}

for(i in 2:T){
	
  sumN[1,i] <- mean(denNpark[i,])
  
  sumNstage[1,i] <- sum(N[1,i,])
  sumNstage[2,i] <- sum(N[2,i,])
  sumNstage[3,i] <- sum(N[3,i,])

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
  
    M.lambda[1,2] <- s.adult.female*f.adult.lambda
    M.lambda[2,1] <- s.fawn*ratio
    M.lambda[2,2] <- s.adult.female
   
  
 
  		
}
			