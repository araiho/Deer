
#####
#####
##### MODEL RUN
##### "Forecasting the Effects of Fertility Control on Overabundant Ungulates" 
##### Created by: Ann Raiho 
##### Last Edited on: 26 November 2014
#####
#####

library(popbio)

#####
##### Add empty rows for forecast #####
#####

add=matrix(0,40,17)
N.obs.no1s.new=rbind(N.obsno1s,add)
forecast.start = 97
N.obs.no1s.new[forecast.start:nrow(N.obs.no1s.new),3]=rep(1:8,5) # add factors for parks
N.obs.no1s.new[forecast.start:104,4]=14 # add factors for years
N.obs.no1s.new[105:112,4]=15
N.obs.no1s.new[113:120,4]=16
N.obs.no1s.new[121:128,4]=17
N.obs.no1s.new[129:136,4]=18
N.obsno1s=N.obs.no1s.new


#####
##### JAGS Implementation #####
#####

library(rjags)
T=13+5

#####
##### Elements in the matrix without zeros get NAs #####
#####

N=array(0,dim=c(4,T,Park)); mu=N;
M=array(0,dim=c(4,4,Park,T));
sumNpark=matrix(NA,T,Park); sumNpark.3yr=sumNpark; sumNpark.1yr=sumNpark ; sumNpark.ster=sumNpark; sumNpark.cull=sumNpark
sumN=matrix(NA,1,T); sumN.3yr=sumN; sumN.1yr=sumN ; sumN.ster=sumN; sumN.cull=sumN
sumNstage=matrix(NA,3,T)

denNpark=matrix(NA,T,Park)

##### Creating empty spots for calculating growth rate of different experiments
M[1,2,,14:18]<-NA ## Matrix rows will be labeled female.fawns, does, male.fawns, bucks
M[2,1,,14:18]<-NA
M[2,2,,14:18]<-NA
M[3,1,,14:18]<-NA
M[3,3,,14:18]<-NA
M.3yr=M; M.1yr=M; M.ster=M; M.cull=M; 
M[1,2,,]<-NA
M[2,1,,]<-NA
M[2,2,,]<-NA
M[3,1,,]<-NA
M[3,3,,]<-NA

M.3yr[4,1,,14:18]<-NA
M.3yr[4,2,,14:18]<-NA
M.3yr[2,4,,14:18]<-NA
M.3yr[4,4,,14:18]<-NA

M.1yr[4,1,,14:18]<-NA
M.1yr[4,2,,14:18]<-NA
M.1yr[2,4,,14:18]<-NA
M.1yr[4,4,,14:18]<-NA

M.ster[4,1,,14:18]<-NA
M.ster[4,2,,14:18]<-NA
M.ster[4,4,,14:18]<-NA

N[1:3,1:13,]=NA
N[1:3,14:T,]=NA;
N.cull=N;N.3yr=N
N.3yr[,14:T,]=NA;
N.1yr=N.3yr; N.ster=N.3yr;

N.init=array(0,dim=c(4,1,Park))
y.alpha1=y.alpha2[y.alpha2[,2]==1,]
y.alpha1=rbind(y.alpha1, y.alpha2[y.alpha2[,2]==2&y.alpha2[,1]==8,])

for(i in 1:Park){
  N.init[1,1,i]=y.alpha1[i,4]
  N.init[2,1,i]=y.alpha1[i,5]
  N.init[3,1,i]=y.alpha1[i,6]
}

N.obs.init=as.matrix(N.obs[N.obs[,4]==1,])

inits= list(list(s.fawn=.7,s.adult.female=.7,s.adult.male=.7,sigma.p=.5,r=.5,beta=10,ratio=.3),
             list(s.fawn=.5,s.adult.female=.5,s.adult.male=.5,sigma.p=1.9,r=.3,beta=25,ratio=.5),
             list(s.fawn=.3,s.adult.female=.9,s.adult.male=.9,sigma.p=.2,r=.5,beta=50,ratio=.7))
            
M.lambda=matrix(0,3,3)
M.lambda[1,2] <- NA
M.lambda[2,1] <- NA
M.lambda[2,2] <- NA

data = list(N.init=N.init, N.obs.init=N.obs.init, T=T, Park=Park, M=M, N.obs=N.obs, N.obsno1s=N.obsno1s, y.alpha=y.alpha2, area=area,N=N,sumNpark=sumNpark,sumN=sumN,M.lambda=M.lambda,denNpark=denNpark)

#####
##### Base Model #####
#####

jM=jags.model(paste(model.dir,"Deer_Model_JAGS.R",sep=""), data=data,inits=inits, n.chain=length(inits), n.adapt=n.adapt)

update(jM,n.iter=n.update)

zm=coda.samples(jM,variable.names=c("s.fawn","s.adult.female","s.adult.male","r","beta","sigma.p","ratio"),n.iter=n.iter,thin=1)

gelman.diag(zm)
#plot(zm)

summary.params=summary(zm)
if(SAVE==TRUE) save(summary.params,file=paste(dump.dir,"summary.params.Rdata",sep=""))

#sum.order=cbind(summary.params$statistics[,1],summary.params$quantiles[,3],summary.params$statistics[,2],summary.params$quantiles[,1],summary.params$quantiles[,5])
#rownames(sum.order)<-c("Carrying Capacity $K_{f}$","Maximum Fecundity $r_{f}$","Sex Ratio $\\phi$","Female Survival $(s_{2})$","Male Survival $(s_{3})$","Juvenile Survival $(s_{1})$","Process Variance $(\\sigma^2)$")
#colnames(sum.order)<- c("Mean","Median","SD","2.5$\\%$ BCI","97.5$\\%$ BCI")
#print(xtable(sum.order,digits=c(0,rep(2,5)),align=c("l","c","c","c","c","c"),display=c("s","f","f","f","f","f")),floating=FALSE, sanitize.text.function=function(x){x})

zmj3=jags.samples(jM,variable.names=c("pvalue.fit","fit","fit.new"),n.iter=n.iter,thin=10)
PB=mean(zmj3$pvalue.fit)
if(SAVE==TRUE) save(PB,file=paste(dump.dir,"denpvalue.Rdata",sep=""))

if(DRAW==TRUE) pdf(paste(dump.dir,"deer.ppc.pdf",sep=""))
plot(zmj3$fit,zmj3$fit.new,xlim=c(2000,10000),ylim=c(2000,10000),xlab=expression(paste(T^obs)),ylab=expression(paste(T^rep)),cex=.7)
abline(a=0,b=1)
text(x=9000,3000,labels=bquote(P[B]==.(format(PB,digits=3))))
if(DRAW==TRUE) dev.off()

zmj=jags.samples(jM,variable.names=c("s.fawn","s.adult.female","s.adult.male","r","beta","ratio","sigma.p"),n.iter=n.iter,thin=1)

if(DRAW==TRUE) pdf(paste(dump.dir,"deer.posteriors.pdf",sep=""))

par(mfrow=c(3,3))
hist(zmj$s.adult.female,xlab="Adult Female Survival",col=8,breaks=20,xlim=c(0,1),freq=FALSE,ylab="Probability Density",main=NA)
lines(seq(0,1,.005),dunif(seq(0,1,.005),0,1),lty=2)
hist(zmj$s.adult.male,xlab="Adult Male Survival",col=8,breaks=20,xlim=c(0,1),freq=FALSE,ylab="Probability Density",main=NA)
lines(seq(0,1,.005),dunif(seq(0,1,.005),0,1),lty=2)
hist(zmj$s.fawn,xlab="Juvenile Survival",col=8,breaks=20,xlim=c(0,1),freq=FALSE,ylab="Probability Density",main=NA)
lines(seq(0,1,.005),dunif(seq(0,1,.005),0,1),lty=2)
hist(zmj$r,xlab="Maximum Birth Rate (fawns per doe)",col=8,breaks=20,freq=FALSE,xlim=c(0,2),ylab="Probability Density",main=NA)
lines(seq(-10,10,.005),dnorm(seq(-10,10,.005),2*3.09*65^-.33,.1304),lty=2)
hist(zmj$ratio,xlab="Juvenile Sex Ratio (females:males)",col=8,breaks=20,xlim=c(.4,.6),ylim=c(0,25),freq=FALSE,ylab="Probability Density",main=NA)
lines(seq(0,1,.005),dbeta(seq(0,1,.005),312,312),lty=2)
hist(zmj$beta,xlab=expression(paste("Carrying Capacity (deer per ", km^2,")")),col=8,breaks=20,freq=FALSE,ylab="Probability Density",main=NA)
lines(seq(0,100,.005),dunif(seq(0,100,.005),0,100),lty=2)
hist(zmj$sigma.p,xlab="Process Variance",col=8,breaks=20,freq=FALSE,ylab="Probability Density",main=NA)
lines(seq(0,2,.005),dunif(seq(0,2,.005),0,2),lty=2)
if(DRAW==TRUE) dev.off()
tvec=as.character(unique(raw.noFOWA[,1]))

#if(DRAW==TRUE) pdf("deer.posteriors.pvar.pdf")
#par(mfrow=c(3,3))
#for(i in 1:Park){
#	hist(zmj$sigma.p[i,,],main=substitute(paste(a), #list(a=tvec[i])),col=8,breaks=20,xlab=NA,freq=FALSE)
#    lines(seq(0,2,.005),dunif(seq(0,2,.005),0,2),lty=2)
#}
#if(DRAW==TRUE) dev.off()


plot(seq(1,150,1),exp(median(zmj$r)-(median(zmj$r)/median(zmj$beta))*seq(1,150,1)),type="l")

if(DRAW==TRUE) pdf(paste(dump.dir,"rK.pdf",sep=""))
hist(zmj$r/zmj$beta,col=8,breaks=20,xlab=expression(paste("r/K (", km^2,")")),main=NA,freq=FALSE,ylab="Probability Density")
if(DRAW==TRUE) dev.off()

#plot(seq(1,median(zmj$beta),1)*-.03+median(zmj$r))

#for(i in 1:length(y.alpha2[,1])){
#	plot(y.alpha2[,8],y.alpha2[,4]/y.alpha2[,5])
#}


#####
##### Forecast Plots -- No Experiments #####
#####

zmj1.base=jags.samples(jM,variable.names=c("denNpark"),n.iter=n.iter,thin=1)
zmj1.N=jags.samples(jM,variable.names=c("N"),n.iter=n.iter,thin=1)

sumN.quant=matrix(0,18,3) ;sumN.3yr.quant=sumN.quant ;sumN.1yr.quant=sumN.quant ;sumN.ster.quant=sumN.quant ;sumN.cull.quant=sumN.quant; sumN.out.quant=sumN.quant

median.create=array(0,dim=c(18,n.iter,3))

for(t in 1:18){
	for(i in 1:n.iter){
		for(c in 1:3){
			median.create[t,i,c] = median(zmj1.base$denNpark[t,,i,c])
		}
	}	
	sumN.quant[t,1]=quantile(median.create[t,,],probs=c(.025))
	sumN.quant[t,2]=quantile(median.create[t,,],probs=c(.50))
	sumN.quant[t,3]=quantile(median.create[t,,],probs=c(.975))
}

save(median.create,file=paste(dump.dir,"median.create.Rdata",sep=""))

if(DRAW==TRUE) pdf(paste(dump.dir,"deer.forecast.pdf",sep=""))
Sum.area=mean(as.numeric(area))
plot(sumN.quant[,2],ylim=c(0,80),ylab=expression(paste("Deer Density ",(km^2))),xlab="Year",type="l",lwd=2, xaxt='n')
lines(sumN.quant[,1],lty=2,lwd=2)
lines(sumN.quant[,3],lty=2,lwd=2)

# for(t in 14:18){
	# points(t,quantile(median.create1[4,1,t,,],probs=c(.025)),col="green",pch=8)
	# points(t,quantile(median.create1[4,1,t,,],probs=c(.975)),col="green",pch=8)
	# points(t,quantile(median.create1[4,1,t,,],probs=c(.5)),col="green",pch=8)
	
	# points(t,quantile(median.create1[4,2,t,,],probs=c(.025)),col="purple",pch=8)
	# points(t,quantile(median.create1[4,2,t,,],probs=c(.975)),col="purple",pch=8)
	# points(t,quantile(median.create1[4,2,t,,],probs=c(.5)),col="purple",pch=8)
	
	# points(t,quantile(median.create1[4,3,t,,],probs=c(.025)),col="blue",pch=8)
	# points(t,quantile(median.create1[4,3,t,,],probs=c(.975)),col="blue",pch=8)
	
	# points(t,quantile(median.create1[4,3,t,,],probs=c(.025)),col="red",pch=8)
	# points(t,quantile(median.create1[4,3,t,,],probs=c(.975)),col="red",pch=8)
# }


#abline(v=13,lty=3,lwd=2,col="red")
axis(1,at=seq(2,18,2),labels=seq(2002,2018,2))

points=numeric(13)
for(i in 1:13){	
	points[i]=median(N.obs[N.obs[,4]==i,5])
	}
	points(points,pch=19)
	for(i in 1:13){
	arrows(i,sd(N.obs[N.obs[,4]==i,5])+points[i],i,points[i]-sd(N.obs[N.obs[,4]==i,5]),code=3, angle=90, length=.01)
}
if(DRAW==TRUE) dev.off()

sumNpark.quant=array(0,dim=c(18,3,10))
sum.den.park=sumNpark.quant
area.vec=as.numeric(area)
tvec=as.character(unique(raw.noFOWA[,1]))

if(DRAW==TRUE) pdf(paste(dump.dir,"deer.bypark.forecast.pdf",sep=""))
par(mfrow=c(3,3))

for(p in 1:Park){
	for(t in 1:18){
		sumNpark.quant[t,1,p]=quantile(zmj1.base$denNpark[t,p,,],probs=c(.025))
		sumNpark.quant[t,2,p]=quantile(zmj1.base$denNpark[t,p,,],probs=c(.5))
		sumNpark.quant[t,3,p]=quantile(zmj1.base$denNpark[t,p,,],probs=c(.975))
	
	}
		sum.den.park[,,p]=sumNpark.quant[,,p]
		
		plot(sum.den.park[,2,p],ylim=c(0,150),ylab=expression(paste("Deer Density ",(km^2))),xlab="Year",main=substitute(paste(a), list(a=tvec[p])),type="l",lwd=1, xaxt='n')
		lines(sum.den.park[,1,p],lty=2,lwd=1)
		lines(sum.den.park[,3,p],lty=2,lwd=1)
		points(N.obs[N.obs[,3]==p,5],pch=20,cex=.5)
		arrows(seq(1,13,1),N.obs[N.obs[,3]==p,5]+N.obs[N.obs[,3]==p,8],seq(1,13,1),N.obs[N.obs[,3]==p,5]-N.obs[N.obs[,3]==p,8],code=3, angle=90, length=.01,lwd=1.5)
		axis(1,at=seq(2,18,2),labels=seq(2002,2018,2))
}
if(DRAW==TRUE) dev.off()

zmjM.base=jags.samples(jM,variable.names=c("M.lambda"),n.iter=n.iter,thin=1)

library(popbio)
lambda=array(0,dim=c(n.iter,3,1))

for(i in 1:n.iter){
	for(c in 1:3){

				h<-eigen.analysis(zmjM.base$M.lambda[1:2,1:2,i,c])
				
				
				lambda[i,c,1] <- h$lambda1
				
	}
	if(i%%5000==0) cat(i," ");flush.console()
}

lambda.base.quant=matrix(0,1,3)

lambda.base.quant[1,1]=quantile(lambda,probs=c(.025))
lambda.base.quant[1,2]=quantile(lambda,probs=c(.5))
lambda.base.quant[1,3]=quantile(lambda,probs=c(.975))

lambda.noaction.mean=mean(lambda)
lambda.noaction.sd=sd(lambda)

save(lambda,file=paste(dump.dir,"lambda.noaction.Rdata",sep=""))

if(DRAW==TRUE) pdf(paste(dump.dir,"lambda.noaction.pdf",sep=""))

plot(density(lambda),xlab=expression(paste(lambda)),main=NA,xlim=c(mean(lambda)-.3,mean(lambda)+.3),lwd=2, xaxt='n')
abline(v=mean(lambda),col="black",lty=2,lwd=2)
axis(1,at=seq(0,1.5,.2))

if(DRAW==TRUE) dev.off()

base.sens=array(0,dim=c(2,2,n.iter,3))
for(i in 1:n.iter){
  for(c in 1:3){
    
    h<-eigen.analysis(zmjM.base$M.lambda[1:2,1:2,i,c])
    
    base.sens[,,i,c] <- h$sensitivities
  }
  if(i%%5000==0) cat(i," ");flush.console()
}


if(DRAW==TRUE) pdf(paste(dump.dir,"base_sens.pdf",sep=""))
par(mfrow=c(2,2))

hist(base.sens[1,2,,],xlim=c(0,1),col=8,breaks=20,main=NA,xlab=expression(paste(s[2],f[it])),freq=FALSE,ylab="Probability Density")
text(x=.8,5,labels=bquote(Mean==.(format(mean(base.sens[1,2,,]),digits=3))))
text(x=.8,4.6,labels=bquote(Median==.(format(median(base.sens[1,2,,]),digits=3))))
text(x=.8,4.2,labels=bquote(paste("2.5% BCI")==.(format(quantile(base.sens[1,2,,],.025),digits=3))))
text(x=.8,3.8,labels=bquote(paste("97.5% BCI")==.(format(quantile(base.sens[1,2,,],.975),digits=3))))

hist(base.sens[2,1,,],xlim=c(0,1),col=8,breaks=20,main=NA,xlab=expression(paste(s[1],m)),freq=FALSE,ylab="Probability Density")
text(x=.2,4.7,labels=bquote(Mean==.(format(mean(base.sens[2,1,,]),digits=3))))
text(x=.2,4.3,labels=bquote(Median==.(format(median(base.sens[2,1,,]),digits=3))))
text(x=.2,3.9,labels=bquote(paste("2.5% BCI")==.(format(quantile(base.sens[2,1,,],.025),digits=3))))
text(x=.2,3.5,labels=bquote(paste("97.5% BCI")==.(format(quantile(base.sens[2,1,,],.975),digits=3))))

hist(base.sens[2,2,,],xlim=c(0,1),col=8,breaks=20,main=NA,xlab=expression(paste(s[2])),freq=FALSE,ylab="Probability Density")
text(x=.2,8.7,labels=bquote(Mean==.(format(mean(base.sens[2,2,,]),digits=3))))
text(x=.2,8,labels=bquote(Median==.(format(median(base.sens[2,2,,]),digits=3))))
text(x=.2,7.2,labels=bquote(paste("2.5% BCI")==.(format(quantile(base.sens[2,2,,],.025),digits=3))))
text(x=.2,6.5,labels=bquote(paste("97.5% BCI")==.(format(quantile(base.sens[2,2,,],.975),digits=3))))
if(DRAW==TRUE) dev.off()


#####
##### Model Experiments #####
#####


#####
##### 20% Treated #####
#####

treated.does = .2

sumNpark.2=matrix(NA,T-5,Park)
sumN.2=matrix(NA,1,T-5)

M.lambda.cull=matrix(0,2,2)
M.lambda.cull[1,2] <- NA
  M.lambda.cull[2,1] <- NA
  M.lambda.cull[2,2] <- NA
  
  M.lambda.ster=matrix(0,3,3)
  M.lambda.ster[1,2] <- NA
  M.lambda.ster[2,1] <-NA
  M.lambda.ster[2,2] <- NA
  M.lambda.ster[3,2] <- NA
  M.lambda.ster[3,3] <- NA
  
   M.lambda.1yr=matrix(0,3,3)
  M.lambda.1yr[1,2] <- NA
  M.lambda.1yr[2,1] <- NA
  M.lambda.1yr[2,2] <- NA
  M.lambda.1yr[2,3] <- NA
  M.lambda.1yr[3,2] <- NA
  M.lambda.1yr[3,3] <- NA
  
   M.lambda.3yr=matrix(0,3,3)
  M.lambda.3yr[1,2] <- NA
  M.lambda.3yr[2,1] <- NA
  M.lambda.3yr[2,2] <- NA
  M.lambda.3yr[2,3] <- NA
  M.lambda.3yr[3,2] <- NA
  M.lambda.3yr[3,3] <- NA

data = list(treated.does=treated.does,N.init=N.init, N.obs.init=N.obs.init, T=T, Park=Park, M=M, M.3yr=M.3yr, M.1yr=M.1yr, M.ster=M.ster, M.cull=M.cull, N.obs=N.obs, N.obsno1s=N.obsno1s, y.alpha=y.alpha2, area=area,N=N,sumNpark=sumNpark.2,sumNpark.3yr=sumNpark.3yr,sumNpark.1yr=sumNpark.1yr,sumNpark.ster=sumNpark.ster,sumN=sumN.2,sumN.3yr=sumN.3yr,sumN.1yr=sumN.1yr,sumN.ster=sumN.ster,sumN.cull=sumN.cull,N.3yr=N.3yr,N.1yr=N.1yr,N.ster=N.ster,N.cull=N.cull,M.lambda.cull=M.lambda.cull,M.lambda.ster=M.lambda.ster,M.lambda.1yr=M.lambda.1yr,M.lambda.3yr=M.lambda.3yr)

jM.2=jags.model(paste(model.dir,"Deer_Model_Experiments_JAGS.R",sep=""), data=data,inits=inits, n.chain=length(inits), n.adapt=n.adapt)

update(jM.2,n.iter=n.update)

zmjM.2=jags.samples(jM.2,variable.names=c("M.lambda.3yr","M.lambda.1yr","M.lambda.ster","M.lambda.cull"),n.iter=n.iter,thin=1)

library(popbio)
lambda=array(0,dim=c(n.iter,3,1)); lambda.3yr.2=lambda; lambda.1yr.2=lambda; lambda.ster.2=lambda; lambda.cull.2=lambda

for(i in 1:n.iter){
	for(c in 1:3){

				h<-eigen.analysis(zmjM.base$M.lambda[1:2,1:2,i,c])
				h.3yr<-eigen.analysis(zmjM.2$M.lambda.3yr[,,i,c])
				h.1yr<-eigen.analysis(zmjM.2$M.lambda.1yr[,,i,c])
				h.ster<-eigen.analysis(zmjM.2$M.lambda.ster[,,i,c])
				h.cull<-eigen.analysis(zmjM.2$M.lambda.cull[,,i,c])
				
				lambda[i,c,1] <- h$lambda1
				lambda.3yr.2[i,c,1] <- h.3yr$lambda1
				lambda.1yr.2[i,c,1] <- h.1yr$lambda1
				lambda.ster.2[i,c,1] <- h.ster$lambda1
				lambda.cull.2[i,c,1] <- h.cull$lambda1
	}
	if(i%%5000==0) cat(i," ");flush.console()
}


lambda.matrix.2=matrix(0,5,3)

lambda.matrix.2[1,1]=quantile(lambda,probs=c(.025))
lambda.matrix.2[1,2]=quantile(lambda,probs=c(.5))
lambda.matrix.2[1,3]=quantile(lambda,probs=c(.975))

lambda.matrix.2[2,1]=quantile(lambda.3yr.2,probs=c(.025))
lambda.matrix.2[2,2]=quantile(lambda.3yr.2,probs=c(.5))
lambda.matrix.2[2,3]=quantile(lambda.3yr.2,probs=c(.975))

lambda.matrix.2[3,1]=quantile(lambda.1yr.2,probs=c(.025))
lambda.matrix.2[3,2]=quantile(lambda.1yr.2,probs=c(.5))
lambda.matrix.2[3,3]=quantile(lambda.1yr.2,probs=c(.975))

lambda.matrix.2[4,1]=quantile(lambda.ster.2,probs=c(.025))
lambda.matrix.2[4,2]=quantile(lambda.ster.2,probs=c(.5))
lambda.matrix.2[4,3]=quantile(lambda.ster.2,probs=c(.975))

lambda.matrix.2[5,1]=quantile(lambda.cull.2,probs=c(.025))
lambda.matrix.2[5,2]=quantile(lambda.cull.2,probs=c(.5))
lambda.matrix.2[5,3]=quantile(lambda.cull.2,probs=c(.975))

#####
##### 40% Treated #####
#####

treated.does = .4

sumNpark.4=matrix(NA,T-5,Park)
sumN.4=matrix(NA,1,T-5)

data = list(treated.does=treated.does,N.init=N.init, N.obs.init=N.obs.init, T=T, Park=Park, M=M, M.3yr=M.3yr, M.1yr=M.1yr, M.ster=M.ster, M.cull=M.cull, N.obs=N.obs, N.obsno1s=N.obsno1s, y.alpha=y.alpha2, area=area,N=N,sumNpark=sumNpark.4,sumNpark.3yr=sumNpark.3yr,sumNpark.1yr=sumNpark.1yr,sumNpark.ster=sumNpark.ster,sumN=sumN.4,sumN.3yr=sumN.3yr,sumN.1yr=sumN.1yr,sumN.ster=sumN.ster,sumN.cull=sumN.cull,N.3yr=N.3yr,N.1yr=N.1yr,N.ster=N.ster,N.cull=N.cull,M.lambda.cull=M.lambda.cull,M.lambda.ster=M.lambda.ster,M.lambda.1yr=M.lambda.1yr,M.lambda.3yr=M.lambda.3yr)

jM.4=jags.model(past(model.dir,"Deer_Model_Experiments_JAGS.R",sep=""), data=data,inits=inits, n.chain=length(inits), n.adapt=n.adapt)

update(jM.4,n.iter=n.update)

zmjM.4=jags.samples(jM.4,variable.names=c("M.lambda.3yr","M.lambda.1yr","M.lambda.ster","M.lambda.cull"),n.iter=n.iter,thin=1)

library(popbio)
lambda.3yr.4=array(0,dim=c(n.iter,3,1)); lambda.1yr.4=lambda.3yr.4; lambda.ster.4=lambda.3yr.4; lambda.cull.4=lambda.3yr.4

for(i in 1:n.iter){
	for(c in 1:3){

				h.3yr<-eigen.analysis(zmjM.4$M.lambda.3yr[,,i,c])
				h.1yr<-eigen.analysis(zmjM.4$M.lambda.1yr[,,i,c])
				h.ster<-eigen.analysis(zmjM.4$M.lambda.ster[,,i,c])
				h.cull<-eigen.analysis(zmjM.4$M.lambda.cull[,,i,c])
			
				lambda.3yr.4[i,c,1] <- h.3yr$lambda1
				lambda.1yr.4[i,c,1] <- h.1yr$lambda1
				lambda.ster.4[i,c,1] <- h.ster$lambda1
				lambda.cull.4[i,c,1] <- h.cull$lambda1
	}
	if(i%%5000==0) cat(i," ");flush.console()
}

lambda.matrix.4=matrix(0,5,3)

lambda.matrix.4[1,1]=quantile(lambda,probs=c(.025))
lambda.matrix.4[1,2]=quantile(lambda,probs=c(.5))
lambda.matrix.4[1,3]=quantile(lambda,probs=c(.975))

lambda.matrix.4[2,1]=quantile(lambda.3yr.4,probs=c(.025))
lambda.matrix.4[2,2]=quantile(lambda.3yr.4,probs=c(.5))
lambda.matrix.4[2,3]=quantile(lambda.3yr.4,probs=c(.975))

lambda.matrix.4[3,1]=quantile(lambda.1yr.4,probs=c(.025))
lambda.matrix.4[3,2]=quantile(lambda.1yr.4,probs=c(.5))
lambda.matrix.4[3,3]=quantile(lambda.1yr.4,probs=c(.975))

lambda.matrix.4[4,1]=quantile(lambda.ster.4,probs=c(.025))
lambda.matrix.4[4,2]=quantile(lambda.ster.4,probs=c(.5))
lambda.matrix.4[4,3]=quantile(lambda.ster.4,probs=c(.975))

lambda.matrix.4[5,1]=quantile(lambda.cull.4,probs=c(.025))
lambda.matrix.4[5,2]=quantile(lambda.cull.4,probs=c(.5))
lambda.matrix.4[5,3]=quantile(lambda.cull.4,probs=c(.975))

#####
##### 60% Treated #####
#####

treated.does = .6

sumNpark.6=matrix(NA,T-5,Park)
sumN.6=matrix(NA,1,T-5)

data = list(treated.does=treated.does,N.init=N.init, N.obs.init=N.obs.init, T=T, Park=Park, M=M, M.3yr=M.3yr, M.1yr=M.1yr, M.ster=M.ster, M.cull=M.cull, N.obs=N.obs, N.obsno1s=N.obsno1s, y.alpha=y.alpha2, area=area,N=N,sumNpark=sumNpark.6,sumNpark.3yr=sumNpark.3yr,sumNpark.1yr=sumNpark.1yr,sumNpark.ster=sumNpark.ster,sumN=sumN.6,sumN.3yr=sumN.3yr,sumN.1yr=sumN.1yr,sumN.ster=sumN.ster,sumN.cull=sumN.cull,N.3yr=N.3yr,N.1yr=N.1yr,N.ster=N.ster,N.cull=N.cull,M.lambda.cull=M.lambda.cull,M.lambda.ster=M.lambda.ster,M.lambda.1yr=M.lambda.1yr,M.lambda.3yr=M.lambda.3yr)

jM.6=jags.model(paste(model.dir,"Deer_Model_Experiments_JAGS.R",sep=""), data=data,inits=inits, n.chain=length(inits), n.adapt=n.adapt)

update(jM.6,n.iter=n.update)

zmjM.6=jags.samples(jM.6,variable.names=c("M.lambda.3yr","M.lambda.1yr","M.lambda.ster","M.lambda.cull"),n.iter=n.iter,thin=1)

library(popbio)
lambda.3yr.6=array(0,dim=c(n.iter,3,1)); lambda.1yr.6=lambda.3yr.6; lambda.ster.6=lambda.3yr.6; lambda.cull.6=lambda.3yr.6

for(i in 1:n.iter){
	for(c in 1:3){

				h<-eigen.analysis(zmjM.base$M.lambda[,,i,c])
				h.3yr<-eigen.analysis(zmjM.6$M.lambda.3yr[,,i,c])
				h.1yr<-eigen.analysis(zmjM.6$M.lambda.1yr[,,i,c])
				h.ster<-eigen.analysis(zmjM.6$M.lambda.ster[,,i,c])
				h.cull<-eigen.analysis(zmjM.6$M.lambda.cull[,,i,c])
				
				lambda[i,c,1] <- h$lambda1
				lambda.3yr.6[i,c,1] <- h.3yr$lambda1
				lambda.1yr.6[i,c,1] <- h.1yr$lambda1
				lambda.ster.6[i,c,1] <- h.ster$lambda1
				lambda.cull.6[i,c,1] <- h.cull$lambda1
	}
	if(i%%5000==0) cat(i," ");flush.console()
}

lambda.matrix.6=matrix(0,5,3)

lambda.matrix.6[1,1]=quantile(lambda,probs=c(.025))
lambda.matrix.6[1,2]=quantile(lambda,probs=c(.5))
lambda.matrix.6[1,3]=quantile(lambda,probs=c(.975))

lambda.matrix.6[2,1]=quantile(lambda.3yr.6,probs=c(.025))
lambda.matrix.6[2,2]=quantile(lambda.3yr.6,probs=c(.5))
lambda.matrix.6[2,3]=quantile(lambda.3yr.6,probs=c(.975))

lambda.matrix.6[3,1]=quantile(lambda.1yr.6,probs=c(.025))
lambda.matrix.6[3,2]=quantile(lambda.1yr.6,probs=c(.5))
lambda.matrix.6[3,3]=quantile(lambda.1yr.6,probs=c(.975))

lambda.matrix.6[4,1]=quantile(lambda.ster.6,probs=c(.025))
lambda.matrix.6[4,2]=quantile(lambda.ster.6,probs=c(.5))
lambda.matrix.6[4,3]=quantile(lambda.ster.6,probs=c(.975))

lambda.matrix.6[5,1]=quantile(lambda.cull.6,probs=c(.025))
lambda.matrix.6[5,2]=quantile(lambda.cull.6,probs=c(.5))
lambda.matrix.6[5,3]=quantile(lambda.cull.6,probs=c(.975))


#####
##### 90% Treated #####
#####

treated.does = .9

sumNpark.9=matrix(NA,T-5,Park)
sumN.9=matrix(NA,1,T-5)

data = list(treated.does=treated.does,N.init=N.init, N.obs.init=N.obs.init, T=T, Park=Park, M=M, M.3yr=M.3yr, M.1yr=M.1yr, M.ster=M.ster, M.cull=M.cull, N.obs=N.obs, N.obsno1s=N.obsno1s, y.alpha=y.alpha2, area=area,N=N,sumNpark=sumNpark.9,sumNpark.3yr=sumNpark.3yr,sumNpark.1yr=sumNpark.1yr,sumNpark.ster=sumNpark.ster,sumN=sumN.9,sumN.3yr=sumN.3yr,sumN.1yr=sumN.1yr,sumN.ster=sumN.ster,sumN.cull=sumN.cull,N.3yr=N.3yr,N.1yr=N.1yr,N.ster=N.ster,N.cull=N.cull,M.lambda.cull=M.lambda.cull,M.lambda.ster=M.lambda.ster,M.lambda.1yr=M.lambda.1yr,M.lambda.3yr=M.lambda.3yr)

jM.9=jags.model(paste(model.dir,"Deer_Model_Experiments_JAGS.R"), data=data,inits=inits, n.chain=length(inits), n.adapt=n.adapt)

zm.9=coda.samples(jM.9,variable.names=c("s.fawn","s.adult.female","s.adult.male","r","beta","ratio","sigma.p"),n.iter=n.iter,thin=1)

update(jM.9,n.iter=n.update)

zmjM.9=jags.samples(jM.9,variable.names=c("M.lambda.3yr","M.lambda.1yr","M.lambda.ster","M.lambda.cull"),n.iter=n.iter,thin=1)

library(popbio)
lambda.3yr.9=array(0,dim=c(n.iter,3,1)); lambda.1yr.9=lambda.3yr.9; lambda.ster.9=lambda.3yr.9; lambda.cull.9=lambda.3yr.9

for(i in 1:n.iter){
	for(c in 1:3){

				h<-eigen.analysis(zmjM.base$M.lambda[,,i,c])
				h.3yr<-eigen.analysis(zmjM.9$M.lambda.3yr[,,i,c])
				h.1yr<-eigen.analysis(zmjM.9$M.lambda.1yr[,,i,c])
				h.ster<-eigen.analysis(zmjM.9$M.lambda.ster[,,i,c])
				h.cull<-eigen.analysis(zmjM.9$M.lambda.cull[,,i,c])
				
				lambda[i,c,1] <- h$lambda1
				lambda.3yr.9[i,c,1] <- h.3yr$lambda1
				lambda.1yr.9[i,c,1] <- h.1yr$lambda1
				lambda.ster.9[i,c,1] <- h.ster$lambda1
				lambda.cull.9[i,c,1] <- h.cull$lambda1
	}
	if(i%%5000==0) cat(i," ");flush.console()
}

lambda.matrix.9=matrix(0,5,3)

lambda.matrix.9[1,1]=quantile(lambda,probs=c(.025))
lambda.matrix.9[1,2]=quantile(lambda,probs=c(.5))
lambda.matrix.9[1,3]=quantile(lambda,probs=c(.975))

lambda.matrix.9[2,1]=quantile(lambda.3yr.9,probs=c(.025))
lambda.matrix.9[2,2]=quantile(lambda.3yr.9,probs=c(.5))
lambda.matrix.9[2,3]=quantile(lambda.3yr.9,probs=c(.975))

lambda.matrix.9[3,1]=quantile(lambda.1yr.9,probs=c(.025))
lambda.matrix.9[3,2]=quantile(lambda.1yr.9,probs=c(.5))
lambda.matrix.9[3,3]=quantile(lambda.1yr.9,probs=c(.975))

lambda.matrix.9[4,1]=quantile(lambda.ster.9,probs=c(.025))
lambda.matrix.9[4,2]=quantile(lambda.ster.9,probs=c(.5))
lambda.matrix.9[4,3]=quantile(lambda.ster.9,probs=c(.975))

lambda.matrix.9[5,1]=quantile(lambda.cull.9,probs=c(.025))
lambda.matrix.9[5,2]=quantile(lambda.cull.9,probs=c(.5))
lambda.matrix.9[5,3]=quantile(lambda.cull.9,probs=c(.975))


zmjN.2=jags.samples(jM.2,variable.names=c("N.3yr","N.1yr","N.ster","N.cull"),n.iter=n.iter,thin=1)

zmjN.4=jags.samples(jM.4,variable.names=c("N.3yr","N.1yr","N.ster","N.cull"),n.iter=n.iter,thin=1)

zmjN.6=jags.samples(jM.6,variable.names=c("N.3yr","N.1yr","N.ster","N.cull"),n.iter=n.iter,thin=1)

zmjN.9=jags.samples(jM.9,variable.names=c("N.3yr","N.1yr","N.ster","N.cull"),n.iter=n.iter,thin=1)

number.treated=array(0,dim=c(17,6,n.iter,3))

save.image(file=paste(dump.dir,"deer.all.mcmc.inc.Rdata",sep=""))
#load("deer.all.mcmc.inc.Rdata")

for(t in 13:17){
	for(i in 1:n.iter){
		for(c in 1:3){
	number.treated[14,t-12,i,c]=sum(zmjN.2$N.cull[2,t,,i,c])*.2
	number.treated[10,t-12,i,c]=sum(zmjN.2$N.ster[2,t,,i,c])*.2
	number.treated[6,t-12,i,c]=(sum(zmjN.2$N.1yr[2,t,,i,c])+sum(zmjN.2$N.1yr[4,t,,i,c]))*.2
	number.treated[2,t-12,i,c]=(sum(zmjN.2$N.3yr[2,t,,i,c])+sum(zmjN.2$N.3yr[4,t,,i,c]))*.2
	
	number.treated[15,t-12,i,c]=sum(zmjN.4$N.cull[2,t,,i,c])*.4
	number.treated[11,t-12,i,c]=sum(zmjN.4$N.ster[2,t,,i,c])*.4
	number.treated[7,t-12,i,c]=(sum(zmjN.4$N.1yr[2,t,,i,c])+sum(zmjN.4$N.1yr[4,t,,i,c]))*.4
	number.treated[3,t-12,i,c]=(sum(zmjN.4$N.3yr[2,t,,i,c])+sum(zmjN.4$N.3yr[4,t,,i,c]))*.4
	
	number.treated[16,t-12,i,c]=sum(zmjN.6$N.cull[2,t,,i,c])*.6
	number.treated[12,t-12,i,c]=sum(zmjN.6$N.ster[2,t,,i,c])*.6
	number.treated[8,t-12,i,c]=(sum(zmjN.6$N.1yr[2,t,,i,c])+sum(zmjN.6$N.1yr[4,t,,i,c]))*.6
	number.treated[4,t-12,i,c]=(sum(zmjN.6$N.3yr[2,t,,i,c])+sum(zmjN.6$N.3yr[4,t,,i,c]))*.6
	
	number.treated[17,t-12,i,c]=sum(zmjN.9$N.cull[2,t,,i,c])*.9
	number.treated[13,t-12,i,c]=sum(zmjN.9$N.ster[2,t,,i,c])*.9
	number.treated[9,t-12,i,c]=(sum(zmjN.9$N.1yr[2,t,,i,c])+sum(zmjN.9$N.1yr[4,t,,i,c]))*.9
	number.treated[5,t-12,i,c]=(sum(zmjN.9$N.3yr[2,t,,i,c])+sum(zmjN.9$N.3yr[4,t,,i,c]))*.9
	
				
		}
	}
	if(t%%1==0) cat(t," ");flush.console()
}

save(number.treated,file=paste(dump.dir,"number.treated.Rdata",sep=""))
#load("number.treated.Rdata")

number.treated1=array(0,dim=c(17,5,2))
for(r in 1:17){
	for(c in 1:5){
		number.treated1[r,c,1]=median(number.treated[r,c,,])
		number.treated1[r,c,2]=sd(number.treated[r,c,,])
	}
}

number.treated2.numbers=number.treated1[,,1]

new.number.treated1=number.treated2.numbers[1,]
for(i in 2:17){
	new.number.treated1=c(new.number.treated1,number.treated2.numbers[i,])
}

number.treated2.numbers.sd=number.treated1[,,2]
new.number.treated1.sd=number.treated2.numbers.sd[1,]
for(i in 2:17){
	new.number.treated1.sd=c(new.number.treated1.sd,number.treated2.numbers.sd[i,])
}

growth.rate.prob=matrix(0,17,4)
growth.rate.prob[,1]=NA
growth.rate.prob[1,1]=0
growth.rate.prob[2,1]=20
growth.rate.prob[6,1]=40
growth.rate.prob[10,1]=60
growth.rate.prob[14,1]=90
growth.rate.prob[1,2:4]=lambda.base.quant

growth.rate.prob[2,2:4]=lambda.matrix.2[5,]
growth.rate.prob[3,2:4]=lambda.matrix.2[4,]
growth.rate.prob[4,2:4]=lambda.matrix.2[3,]
growth.rate.prob[5,2:4]=lambda.matrix.2[2,]

growth.rate.prob[6,2:4]=lambda.matrix.4[5,]
growth.rate.prob[7,2:4]=lambda.matrix.4[4,]
growth.rate.prob[8,2:4]=lambda.matrix.4[3,]
growth.rate.prob[9,2:4]=lambda.matrix.4[2,]

growth.rate.prob[10,2:4]=lambda.matrix.6[5,]
growth.rate.prob[11,2:4]=lambda.matrix.6[4,]
growth.rate.prob[12,2:4]=lambda.matrix.6[3,]
growth.rate.prob[13,2:4]=lambda.matrix.6[2,]

growth.rate.prob[14,2:4]=lambda.matrix.9[5,]
growth.rate.prob[15,2:4]=lambda.matrix.9[4,]
growth.rate.prob[16,2:4]=lambda.matrix.9[3,]
growth.rate.prob[17,2:4]=lambda.matrix.9[2,]


addon<-matrix(c("No Action",rep(c("Cull","Sterilize","1 Year","3 Year"),4)),17,1)
growth.rate.prob=data.frame(addon,growth.rate.prob)

colnames(growth.rate.prob)<-c("Treatment","Percent Treated","2.5$\\%$ BCI","Median","97.5$\\%$ BCI")

rel.diff.matrix = matrix(NA,17,3)
rel.diff.matrix[2,]=quantile(lambda-lambda.cull.2,probs=c(.05,.5,.975))
rel.diff.matrix[3,]=quantile(lambda-lambda.ster.2,probs=c(.05,.5,.975))
rel.diff.matrix[4,]=quantile(lambda-lambda.1yr.2,probs=c(.05,.5,.975))
rel.diff.matrix[5,]=quantile(lambda-lambda.3yr.2,probs=c(.05,.5,.975))

rel.diff.matrix[6,]=quantile(lambda-lambda.cull.4,probs=c(.05,.5,.975))
rel.diff.matrix[7,]=quantile(lambda-lambda.ster.4,probs=c(.05,.5,.975))
rel.diff.matrix[8,]=quantile(lambda-lambda.1yr.4,probs=c(.05,.5,.975))
rel.diff.matrix[9,]=quantile(lambda-lambda.3yr.4,probs=c(.05,.5,.975))

rel.diff.matrix[10,]=quantile(lambda-lambda.cull.6,probs=c(.05,.5,.975))
rel.diff.matrix[11,]=quantile(lambda-lambda.ster.6,probs=c(.05,.5,.975))
rel.diff.matrix[12,]=quantile(lambda-lambda.1yr.6,probs=c(.05,.5,.975))
rel.diff.matrix[13,]=quantile(lambda-lambda.3yr.6,probs=c(.05,.5,.975))

rel.diff.matrix[14,]=quantile(lambda-lambda.cull.9,probs=c(.05,.5,.975))
rel.diff.matrix[15,]=quantile(lambda-lambda.ster.9,probs=c(.05,.5,.975))
rel.diff.matrix[16,]=quantile(lambda-lambda.1yr.9,probs=c(.05,.5,.975))
rel.diff.matrix[17,]=quantile(lambda-lambda.3yr.9,probs=c(.05,.5,.975))

growth.rate.prob=cbind(growth.rate.prob,rel.diff.matrix)
colnames(growth.rate.prob)<-c("Treatment","Percent Treated",
"2.5$\\%$ BCI","Median","97.5$\\%$ BCI","2.5$\\%$ BCI RNA","Median Difference RNA",
"97.5$\\%$ BCI RNA")

save(growth.rate.prob,file=paste(dump.dir,"growth.rate.prob.Rdata",sep=""))
#load("growth.rate.prob.Rdata")

zmj1.2=jags.samples(jM.2,variable.names=c("denNpark.3yr","denNpark.1yr","denNpark.ster","denNpark.cull"),n.iter=n.iter,thin=1)

zmj1.4=jags.samples(jM.4,variable.names=c("denNpark.3yr","denNpark.1yr","denNpark.ster","denNpark.cull"),n.iter=n.iter,thin=1)

zmj1.6=jags.samples(jM.6,variable.names=c("denNpark.3yr","denNpark.1yr","denNpark.ster","denNpark.cull"),n.iter=n.iter,thin=1)

zmj1.9=jags.samples(jM.9,variable.names=c("denNpark.3yr","denNpark.1yr","denNpark.ster","denNpark.cull"),n.iter=n.iter,thin=1)

sumN.quant=matrix(0,18,3) ;sumN.3yr.quant=sumN.quant ;sumN.1yr.quant=sumN.quant ;sumN.ster.quant=sumN.quant ;sumN.cull.quant=sumN.quant; sumN.out.quant=sumN.quant

median.create1=array(0,dim=c(4,4,18,n.iter,3))

for(t in 13:18){
	for(i in 1:n.iter){
		for(c in 1:3){
			median.create1[1,1,t,i,c] = median(zmj1.2$denNpark.3yr[t,,i,c])
			median.create1[1,2,t,i,c] = median(zmj1.2$denNpark.1yr[t,,i,c])
			median.create1[1,3,t,i,c] = median(zmj1.2$denNpark.ster[t,,i,c])
			median.create1[1,4,t,i,c] = median(zmj1.2$denNpark.cull[t,,i,c])
			
			median.create1[2,1,t,i,c] = median(zmj1.4$denNpark.3yr[t,,i,c])
			median.create1[2,2,t,i,c] = median(zmj1.4$denNpark.1yr[t,,i,c])
			median.create1[2,3,t,i,c] = median(zmj1.4$denNpark.ster[t,,i,c])
			median.create1[2,4,t,i,c] = median(zmj1.4$denNpark.cull[t,,i,c])
			
			median.create1[3,1,t,i,c] = median(zmj1.6$denNpark.3yr[t,,i,c])
			median.create1[3,2,t,i,c] = median(zmj1.6$denNpark.1yr[t,,i,c])
			median.create1[3,3,t,i,c] = median(zmj1.6$denNpark.ster[t,,i,c])
			median.create1[3,4,t,i,c] = median(zmj1.6$denNpark.cull[t,,i,c])
			
			median.create1[4,1,t,i,c] = median(zmj1.9$denNpark.3yr[t,,i,c])
			median.create1[4,2,t,i,c] = median(zmj1.9$denNpark.1yr[t,,i,c])
			median.create1[4,3,t,i,c] = median(zmj1.9$denNpark.ster[t,,i,c])
			median.create1[4,4,t,i,c] = median(zmj1.9$denNpark.cull[t,,i,c])
		}
	}	
}

save(median.create1,file=paste(dump.dir,"median.create1.Rdata",sep=""))

alt.MGMT=matrix(0,80,6)
alt.MGMT[,1]=c(rep(1,20),rep(2,20),rep(3,20),rep(4,20))
alt.MGMT[,2]=rep(c(rep(1,5),rep(2,5),rep(3,5),rep(4,5)),4)
alt.MGMT[,3]=rep(c(14,15,16,17,18),16)

objective1=5
objective=15
library(xtable)
for(i in 1:length(alt.MGMT[,1])){
			alt.MGMT[i,4:6]=c(ecdf(median.create1[alt.MGMT[i,2],alt.MGMT[i,1],alt.MGMT[i,3],,])(objective1),ecdf(median.create1[alt.MGMT[i,2],alt.MGMT[i,1],alt.MGMT[i,3],,])(objective)-ecdf(median.create1[alt.MGMT[i,2],alt.MGMT[i,1],alt.MGMT[i,3],,])(objective1),1-ecdf(median.create1[alt.MGMT[i,2],alt.MGMT[i,1],alt.MGMT[i,3],,])(objective))

}

donothing=matrix(0,5,6)
for(t in 14:18){
donothing[t-13,4:6]=c(ecdf(median.create[t,,])(objective1),ecdf(median.create[t,,])(objective)-ecdf(median.create[t,,])(objective1),1-ecdf(median.create[t,,])(objective))
}
alt.MGMT=rbind(donothing,alt.MGMT)

alt.MGMTsave=alt.MGMT[,3:6]
newnames=matrix(NA,85,1)
newnames[1,1]="No Action"
newnames[6,1]="3 Year"
newnames[26,1]="1 Year"
newnames[46,1]="Sterilize"
newnames[66,1]="Cull"
alt.MGMT[,2]=NA
alt.MGMT[seq(6,85,5),2]=rep(c(20,40,60,90),4)
alt.MGMT[,3]=rep(c(2014,2015,2016,2017,2018),17)

alt.MGMT1=data.frame(newnames,alt.MGMT[,2:6])

alt.MGMT1=cbind(alt.MGMT1[,1:6],new.number.treated1)

colnames(alt.MGMT1)<-c("Treatment","Percent Treated","Year","P <","P in","P >","Number Treated")

print(xtable(alt.MGMT1,digits=c(0,0,2,4,rep(4,3),0),align=c("l","l","c","c","c","c","c","c"),display=c("s","s","f","f","f","f","f","f")),floating=FALSE,
sanitize.text.function=function(x){x},
include.rownames=FALSE,tabular.environment = 'longtable')

save(alt.MGMT1,file=paste(dump.dir,"alt.MGMT.Rdata",sep=""))


#####
##### Plotting #####
#####

yvec=as.character(seq(2001,2018,1))
labs = c("20%","40%","60%","90%")

if(DRAW==TRUE) pdf(paste(dump.dir,"experiments14-18.pdf",sep=""))
par(mfrow=c(4,2),mar=c(2,4.5,2,3))

for(r in 1:4){
	for(t in c(14,18)){
		if(t==14){
plot(density(median.create[t,,]),col="black",lwd=3,xlim=c(0,60),ylim=c(0,1.1),main=substitute(paste(a), list(a=yvec[t])),xlab=expression(paste("Deer Density ",(km^2))),ylab = "Probability Density")
} else{
plot(density(median.create[t,,]),col="black",lwd=3,xlim=c(0,60),ylim=c(0,1.1),main=substitute(paste(a), list(a=yvec[t])),xlab=expression(paste("Deer Density ",(km^2))),ylab=NA)
mtext(labs[r], las=1,
       side=2,line=3)
}
lines(density(median.create1[r,1,t,,]),col="green",lwd=2)
lines(density(median.create1[r,2,t,,]),col="purple",lwd=2)
lines(density(median.create1[r,3,t,,]),col="blue",lwd=2)
lines(density(median.create1[r,4,t,,]),col="red",lwd=2)
abline(v=5,lty=2)
abline(v=15,lty=2)

	}
}
if(DRAW==TRUE) dev.off()

#for(r in 1:4){
#plot(c(median(median.create[14,,]),median(median.create[15,,]),median(median.create[16,,]),median(median.create[17,,]),median(median.create[18,,])),col="black",lwd=3,main=substitute(paste(a), list(a=yvec[t])),xlab=expression(paste("Deer Density ",(km^2))),ylim=c(0,60),typ="l")
#lines(c(quantile(median.create[14,,],probs=c(.975)),quantile(median.create[15,,],probs=c(.975)),quantile(median.create[16,,],probs=c(.975)),quantile(median.create[17,,],probs=c(.975)),quantile(median.create[18,,],probs=c(.975))),col="black",lwd=3,lty=2)
#lines(c(quantile(median.create[14,,],probs=c(.05)),quantile(median.create[15,,],probs=c(.05)),quantile(median.create[16,,],probs=c(.05)),quantile(median.create[17,,],probs=c(.05)),quantile(median.create[18,,],probs=c(.05))),col="black",lwd=3,lty=2)
#lines(c(median(median.create1[r,1,14,,]),median(median.create1[r,1,15,,]),median(median.create1[r,1,16,,]),median(median.create1[r,1,17,,]),median(median.create1[r,1,18,,])),col="green",lwd=2)
#lines(c(quantile(median.create1[r,1,14,,],probs=c(.975)),quantile(median.create1[r,1,15,,],probs=c(.975)),quantile(median.create1[r,1,16,,],probs=c(.975)),quantile(median.create1[r,1,17,,],probs=c(.975)),quantile(median.create1[r,1,18,,],probs=c(.975))),col="green",lwd=2,lty=2)
#lines(c(quantile(median.create1[r,1,14,,],probs=c(.05)),quantile(median.create1[r,1,15,,],probs=c(.05)),quantile(median.create1[r,1,16,,],probs=c(.05)),quantile(median.create1[r,1,17,,],probs=c(.05)),quantile(median.create1[r,1,18,,],probs=c(.05))),col="green",lwd=2,lty=2)
#lines(c(median(median.create1[r,2,14,,]),median(median.create1[r,2,15,,]),median(median.create1[r,2,16,,]),median(median.create1[r,2,17,,]),median(median.create1[r,2,18,,])),col="purple",lwd=2)
#lines(c(median(median.create1[r,3,14,,]),median(median.create1[r,3,15,,]),median(median.create1[r,3,16,,]),median(median.create1[r,3,17,,]),median(median.create1[r,3,18,,])),col="blue",lwd=2)
#lines(c(median(median.create1[r,4,14,,]),median(median.create1[r,4,15,,]),median(median.create1[r,4,16,,]),median(median.create1[r,4,17,,]),median(median.create1[r,4,18,,])),col="red",lwd=2)

#}



#####
##### 20% Treated AFter Culling 90% #####
#####

N=array(0,dim=c(4,T,Park));
N[1:3,1:13,]=NA
N[1:3,14:T,]=0;

N.cull=N;
N.cull[1:3,14:T,]=NA;

N.3yr=array(0,dim=c(4,T,Park));
N.3yr[1:3,1:14,]=NA;
N.3yr[,15:T,]=NA;

N.1yr=N.3yr; N.ster=N.3yr;

treated.does = .2

sumNpark.2=matrix(NA,T-5,Park)
sumN.2=matrix(NA,1,T-5)

data = list(treated.does=treated.does,N.init=N.init, N.obs.init=N.obs.init, T=T, Park=Park, M=M, M.3yr=M.3yr, M.1yr=M.1yr, M.ster=M.ster, M.cull=M.cull, N.obs=N.obs, N.obsno1s=N.obsno1s, y.alpha=y.alpha2, area=area,N=N,sumNpark=sumNpark.2,sumNpark.3yr=sumNpark.3yr,sumNpark.1yr=sumNpark.1yr,sumNpark.ster=sumNpark.ster,sumN=sumN.2,sumN.3yr=sumN.3yr,sumN.1yr=sumN.1yr,sumN.ster=sumN.ster,sumN.cull=sumN.cull,N.3yr=N.3yr,N.1yr=N.1yr,N.ster=N.ster,N.cull=N.cull)

jM.2.cull1st=jags.model(paste(model.dir,"Deer_Model_Experiments_Cull1st_JAGS.R",sep=""), data=data,inits=inits, n.chain=length(inits), n.adapt=n.adapt)

update(jM.2.cull1st,n.iter=n.update)

treated.does = .4

sumNpark.4=matrix(NA,T-5,Park)
sumN.4=matrix(NA,1,T-5)

data = list(treated.does=treated.does,N.init=N.init, N.obs.init=N.obs.init, T=T, Park=Park, M=M, M.3yr=M.3yr, M.1yr=M.1yr, M.ster=M.ster, M.cull=M.cull, N.obs=N.obs, N.obsno1s=N.obsno1s, y.alpha=y.alpha2, area=area,N=N,sumNpark=sumNpark.4,sumNpark.3yr=sumNpark.3yr,sumNpark.1yr=sumNpark.1yr,sumNpark.ster=sumNpark.ster,sumN=sumN.4,sumN.3yr=sumN.3yr,sumN.1yr=sumN.1yr,sumN.ster=sumN.ster,sumN.cull=sumN.cull,N.3yr=N.3yr,N.1yr=N.1yr,N.ster=N.ster,N.cull=N.cull)

jM.4.cull1st=jags.model(paste(model.dir,"Deer_Model_Experiments_Cull1st_JAGS.R",sep=""), data=data,inits=inits, n.chain=length(inits), n.adapt=n.adapt)

update(jM.4.cull1st,n.iter=n.update)

treated.does = .6

sumNpark.6=matrix(NA,T-5,Park)
sumN.6=matrix(NA,1,T-5)

data = list(treated.does=treated.does,N.init=N.init, N.obs.init=N.obs.init, T=T, Park=Park, M=M, M.3yr=M.3yr, M.1yr=M.1yr, M.ster=M.ster, M.cull=M.cull, N.obs=N.obs, N.obsno1s=N.obsno1s, y.alpha=y.alpha2, area=area,N=N,sumNpark=sumNpark.6,sumNpark.3yr=sumNpark.3yr,sumNpark.1yr=sumNpark.1yr,sumNpark.ster=sumNpark.ster,sumN=sumN.6,sumN.3yr=sumN.3yr,sumN.1yr=sumN.1yr,sumN.ster=sumN.ster,sumN.cull=sumN.cull,N.3yr=N.3yr,N.1yr=N.1yr,N.ster=N.ster,N.cull=N.cull)

jM.6.cull1st=jags.model(paste(model.dir,"Deer_Model_Experiments_Cull1st_JAGS.R",sep=""), data=data,inits=inits, n.chain=length(inits), n.adapt=n.adapt)

update(jM.6.cull1st,n.iter=n.update)

treated.does = .9

sumNpark.9=matrix(NA,T-5,Park)
sumN.9=matrix(NA,1,T-5)

data = list(treated.does=treated.does,N.init=N.init, N.obs.init=N.obs.init, T=T, Park=Park, M=M, M.3yr=M.3yr, M.1yr=M.1yr, M.ster=M.ster, M.cull=M.cull, N.obs=N.obs, N.obsno1s=N.obsno1s, y.alpha=y.alpha2, area=area,N=N,sumNpark=sumNpark.9,sumNpark.3yr=sumNpark.3yr,sumNpark.1yr=sumNpark.1yr,sumNpark.ster=sumNpark.ster,sumN=sumN.9,sumN.3yr=sumN.3yr,sumN.1yr=sumN.1yr,sumN.ster=sumN.ster,sumN.cull=sumN.cull,N.3yr=N.3yr,N.1yr=N.1yr,N.ster=N.ster,N.cull=N.cull)

jM.9.cull1st=jags.model(paste(model.dir,"Deer_Model_Experiments_Cull1st_JAGS.R",sep=""), data=data,inits=inits, n.chain=length(inits), n.adapt=n.adapt)

update(jM.9.cull1st,n.iter=n.update)

zmjN.2.cull1st=jags.samples(jM.2.cull1st,variable.names=c("N.3yr","N.1yr","N.ster","N.cull"),n.iter=n.iter,thin=1)

zmjN.4.cull1st=jags.samples(jM.4.cull1st,variable.names=c("N.3yr","N.1yr","N.ster","N.cull"),n.iter=n.iter,thin=1)

zmjN.6.cull1st=jags.samples(jM.6.cull1st,variable.names=c("N.3yr","N.1yr","N.ster","N.cull"),n.iter=n.iter,thin=1)

zmjN.9.cull1st=jags.samples(jM.9.cull1st,variable.names=c("N.3yr","N.1yr","N.ster","N.cull"),n.iter=n.iter,thin=1)

number.treated.cull1st=array(0,dim=c(16,4,n.iter,3))

for(t in 14:17){
	for(i in 1:n.iter){
		for(c in 1:3){

	number.treated.cull1st[13,t-13,i,c]=sum(zmjN.2.cull1st$N.cull[2,t,,i,c])*.2
	number.treated.cull1st[9,t-13,i,c]=sum(zmjN.2.cull1st$N.ster[2,t,,i,c])*.2
	number.treated.cull1st[5,t-13,i,c]=(sum(zmjN.2.cull1st$N.1yr[2,t,,i,c])+sum(zmjN.2.cull1st$N.1yr[4,t,,i,c]))*.2
	number.treated.cull1st[1,t-13,i,c]=(sum(zmjN.2.cull1st$N.3yr[2,t,,i,c])+sum(zmjN.2.cull1st$N.3yr[4,t,,i,c]))*.2
	
	number.treated.cull1st[14,t-13,i,c]=sum(zmjN.4.cull1st$N.cull[2,t,,i,c])*.4
	number.treated.cull1st[10,t-13,i,c]=sum(zmjN.4.cull1st$N.ster[2,t,,i,c])*.4
	number.treated.cull1st[6,t-13,i,c]=(sum(zmjN.4.cull1st$N.1yr[2,t,,i,c])+sum(zmjN.4.cull1st$N.1yr[4,t,,i,c]))*.4
	number.treated.cull1st[2,t-13,i,c]=(sum(zmjN.4.cull1st$N.3yr[2,t,,i,c])+sum(zmjN.4.cull1st$N.3yr[4,t,,i,c]))*.4
	
	number.treated.cull1st[15,t-13,i,c]=sum(zmjN.6.cull1st$N.cull[2,t,,i,c])*.6
	number.treated.cull1st[11,t-13,i,c]=sum(zmjN.6.cull1st$N.ster[2,t,,i,c])*.6
	number.treated.cull1st[7,t-13,i,c]=(sum(zmjN.6.cull1st$N.1yr[2,t,,i,c])+sum(zmjN.6.cull1st$N.1yr[4,t,,i,c]))*.6
	number.treated.cull1st[3,t-13,i,c]=(sum(zmjN.6.cull1st$N.3yr[2,t,,i,c])+sum(zmjN.6.cull1st$N.3yr[4,t,,i,c]))*.6
	
	number.treated.cull1st[16,t-13,i,c]=sum(zmjN.9.cull1st$N.cull[2,t,,i,c])*.9
	number.treated.cull1st[12,t-13,i,c]=sum(zmjN.9.cull1st$N.ster[2,t,,i,c])*.9
	number.treated.cull1st[8,t-13,i,c]=(sum(zmjN.9.cull1st$N.1yr[2,t,,i,c])+sum(zmjN.9.cull1st$N.1yr[4,t,,i,c]))*.9
	number.treated.cull1st[4,t-13,i,c]=(sum(zmjN.9.cull1st$N.3yr[2,t,,i,c])+sum(zmjN.9.cull1st$N.3yr[4,t,,i,c]))*.9
				
		}
	}
	if(t%%1==0) cat(t," ");flush.console()
}

number.treated3=array(0,dim=c(16,4,2))
for(r in 1:16){
	for(c in 1:4){
		number.treated3[r,c,1]=median(number.treated.cull1st[r,c,,])
		number.treated3[r,c,2]=sd(number.treated.cull1st[r,c,,])
	}
}

number.treated3.numbers=number.treated3[,,1]
new.number.treated3=number.treated3.numbers[2,]
for(i in 2:16){
	new.number.treated3=c(new.number.treated3,number.treated3.numbers[i,])
}

number.treated3.numbers.sd=number.treated3[,,2]
new.number.treated3.sd=number.treated3.numbers.sd[1,]
for(i in 2:16){	new.number.treated3.sd=c(new.number.treated3.sd,number.treated3.numbers.sd[i,])
}

zmj1.2.cull1st=jags.samples(jM.2.cull1st,variable.names=c("denNpark.3yr","denNpark.1yr","denNpark.ster","denNpark.cull"),n.iter=n.iter,thin=1)

zmj1.4.cull1st=jags.samples(jM.4.cull1st,variable.names=c("denNpark.3yr","denNpark.1yr","denNpark.ster","denNpark.cull"),n.iter=n.iter,thin=1)

zmj1.6.cull1st=jags.samples(jM.6.cull1st,variable.names=c("denNpark.3yr","denNpark.1yr","denNpark.ster","denNpark.cull"),n.iter=n.iter,thin=1)

zmj1.9.cull1st=jags.samples(jM.9.cull1st,variable.names=c("denNpark.3yr","denNpark.1yr","denNpark.ster","denNpark.cull"),n.iter=n.iter,thin=1)


median.create2=array(0,dim=c(4,4,18,n.iter,3))

for(t in 14:18){
	for(i in 1:n.iter){
		for(c in 1:3){
			median.create2[1,1,t,i,c] = median(zmj1.2.cull1st$denNpark.3yr[t,,i,c])
			median.create2[1,2,t,i,c] = median(zmj1.2.cull1st$denNpark.1yr[t,,i,c])
			median.create2[1,3,t,i,c] = median(zmj1.2.cull1st$denNpark.ster[t,,i,c])
			median.create2[1,4,t,i,c] = median(zmj1.2.cull1st$denNpark.cull[t,,i,c])
			
			median.create2[2,1,t,i,c] = median(zmj1.4.cull1st$denNpark.3yr[t,,i,c])
			median.create2[2,2,t,i,c] = median(zmj1.4.cull1st$denNpark.1yr[t,,i,c])
			median.create2[2,3,t,i,c] = median(zmj1.4.cull1st$denNpark.ster[t,,i,c])
			median.create2[2,4,t,i,c] = median(zmj1.4.cull1st$denNpark.cull[t,,i,c])
			
			median.create2[3,1,t,i,c] = median(zmj1.6.cull1st$denNpark.3yr[t,,i,c])
			median.create2[3,2,t,i,c] = median(zmj1.6.cull1st$denNpark.1yr[t,,i,c])
			median.create2[3,3,t,i,c] = median(zmj1.6.cull1st$denNpark.ster[t,,i,c])
			median.create2[3,4,t,i,c] = median(zmj1.6.cull1st$denNpark.cull[t,,i,c])
			
			median.create2[4,1,t,i,c] = median(zmj1.9.cull1st$denNpark.3yr[t,,i,c])
			median.create2[4,2,t,i,c] = median(zmj1.9.cull1st$denNpark.1yr[t,,i,c])
			median.create2[4,3,t,i,c] = median(zmj1.9.cull1st$denNpark.ster[t,,i,c])
			median.create2[4,4,t,i,c] = median(zmj1.9.cull1st$denNpark.cull[t,,i,c])
		}
	}	
}

save(median.create2,file=paste(dump.dir,"median.create2.Rdata",sep=""))

alt.MGMT.cull1st=matrix(0,64,6)
alt.MGMT.cull1st[,1]=c(rep(1,16),rep(2,16),rep(3,16),rep(4,16))
alt.MGMT.cull1st[,2]=rep(c(rep(1,4),rep(2,4),rep(3,4),rep(4,4)),4)
alt.MGMT.cull1st[,3]=rep(c(15,16,17,18),16)

objective1=5
objective=15
library(xtable)
for(i in 1:length(alt.MGMT.cull1st[,1])){
			alt.MGMT.cull1st[i,4:6]=c(ecdf(median.create2[alt.MGMT.cull1st[i,2],alt.MGMT.cull1st[i,1],alt.MGMT.cull1st[i,3],,])(objective1),ecdf(median.create2[alt.MGMT.cull1st[i,2],alt.MGMT.cull1st[i,1],alt.MGMT.cull1st[i,3],,])(objective)-ecdf(median.create2[alt.MGMT.cull1st[i,2],alt.MGMT.cull1st[i,1],alt.MGMT.cull1st[i,3],,])(objective1),1-ecdf(median.create2[alt.MGMT.cull1st[i,2],alt.MGMT.cull1st[i,1],alt.MGMT.cull1st[i,3],,])(objective))

}

alt.MGMTsave.cull1st=alt.MGMT.cull1st[,3:6]
newnames.cull1st=matrix(NA,64,1)
newnames.cull1st[1,1]="3 Year"
newnames.cull1st[17,1]="1 Year"
newnames.cull1st[33,1]="Sterilize"
newnames.cull1st[49,1]="Cull"
alt.MGMT.cull1st[,2]=NA
alt.MGMT.cull1st[seq(1,64,4),2]=rep(c(20,40,60,90),4)
alt.MGMT.cull1st[,3]=rep(c(2015,2016,2017,2018),16)

alt.MGMT1.cull1st=data.frame(newnames.cull1st,alt.MGMT.cull1st[,2:6])

alt.MGMT1.cull1st=cbind(alt.MGMT1.cull1st[,1:6],new.number.treated3)

colnames(alt.MGMT1.cull1st)<-c("Treatment","Percent Treated","Year","P <","P in","P >","Number Treated")

save(alt.MGMT1.cull1st,file=paste(dump.dir,"alt.MGMT1.cull1st.Rdata",sep=""))

yvec=as.character(seq(2001,2018,1))

if(DRAW==TRUE) pdf(paste("experiments15-18.pdf",sep=""))
par(mfrow=c(4,2),mar=c(2,4.5,2,3))

for(r in 1:4){
	for(t in c(15,18)){
		if(t==15){
plot(density(median.create[t,,]),col="white",lwd=3,xlim=c(0,20),ylim=c(0,1.1),main=substitute(paste(a), list(a=yvec[t])),xlab=expression(paste("Deer Density ",(km^2))),ylab="Probability Density")
		} else {			plot(density(median.create[t,,]),col="white",lwd=3,xlim=c(0,20),ylim=c(0,1.1),main=substitute(paste(a), list(a=yvec[t])),xlab=expression(paste("Deer Density ",(km^2))),ylab=NA)
mtext(labs[r], las=1, side=2, line=3)
		}
lines(density(median.create2[r,1,t,,]),col="green",lwd=2)
lines(density(median.create2[r,2,t,,]),col="purple",lwd=2)
lines(density(median.create2[r,3,t,,]),col="blue",lwd=2)
lines(density(median.create2[r,4,t,,]),col="red",lwd=2)
abline(v=5,lty=2)
abline(v=15,lty=2)

	}
}
if(DRAW==TRUE) dev.off()



##################################################################################################################################################################################################################################################################################################################################################################################################
