library(rjags)
S=(read.csv("~/Documents/Deer/Data/Survival.csv"))

#####
##### Adult Female #####
#####

y.m=na.omit(S$SP.Adult.F)
y.m=as.list(y.m)
y.sd=na.omit(S$SD.Adult.F)
y.sd=as.list(y.sd)

data = list(y.m=y.m, y.sd=y.sd)

inits= list(list(a=5,b=5),list(a=10,b=1),list(a=1,b=10))

jM=jags.model(paste(model.dir,"Survival JAGS.R",sep=""), data=data,inits=inits, n.chain=length(inits), n.adapt=n.adapt)

update(jM,n.iter=n.update)

zm=coda.samples(jM,variable.names=c("a","b","fe"),n.iter=n.iter,n.thin=1)

gelman.diag(zm)
mu.a = summary(zm)$statistics[1,1]
sd.a = summary(zm)$statistics[1,2]

mu.b = summary(zm)$statistics[2,1]
sd.b = summary(zm)$statistics[2,2]

    shape.a <- mu.a^2/sd.a^2
    rate.a <- mu.a/sd.a^2

    shape.b <- mu.b^2/sd.b^2
    rate.b <- mu.b/sd.b^2

#####
##### Juveniles #####
#####

y.m=na.omit(c(S$SP.Yearling.F,S$SP.Yearling.M))
y.m=as.list(y.m)
y.sd=na.omit(c(S$SD.Yearling.F,S$SD.Yearling.M))
y.sd=as.list(y.sd)

data = list(y.m=y.m, y.sd=y.sd)

inits= list(list(a=355,b=125),list(a=355+.2*355,b=125+.2*125),list(a=355-.2*355,b=125-.2*125))

jM=jags.model(paste(model.dir,"Survival JAGS.R",sep=""), data=data,inits=inits, n.chain=length(inits), n.adapt=n.adapt)

update(jM,n.iter=n.update)

zm=coda.samples(jM,variable.names=c("a","b","fe"),n.iter=n.iter,n.thin=1)

gelman.diag(zm)
mu.a1 = summary(zm)$statistics[1,1]
sd.a1 = summary(zm)$statistics[1,2]

mu.b1 = summary(zm)$statistics[2,1]
sd.b1 = summary(zm)$statistics[2,2]

    shape.a1 <- mu.a1^2/sd.a1^2
    rate.a1 <- mu.a1/sd.a1^2

    shape.b1 <- mu.b1^2/sd.b1^2
    rate.b1 <- mu.b1/sd.b1^2

#####
##### Adult Male #####
#####


y.m=na.omit(S$SP.Adult.M)
y.m=as.list(y.m)
y.sd=na.omit(S$SD.Adult.M)
y.sd=as.list(y.sd)

data = list(y.m=y.m, y.sd=y.sd)

inits= NA

jM=jags.model(paste(model.dir,"survival JAGS 2.R",sep=""), data=data,n.chain = 3, n.adapt=n.adapt)

update(jM,n.iter=n.update)

zm=coda.samples(jM,variable.names=c("mu.true"),n.iter=n.iter,n.thin=10)

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

mu.a2 = as.numeric(estBetaParams(summary(zm)$statistics[1],summary(zm)$statistics[2]^2)$alpha)
mu.b2 = as.numeric(estBetaParams(summary(zm)$statistics[1],summary(zm)$statistics[2]^2)$beta)

sd.a2 = 5
sd.b2 = 5

    shape.a2 <- mu.a2^2/sd.a2^2
    rate.a2 <- mu.a2/sd.a2^2
    
    shape.b2 <- mu.b2^2/sd.b2^2
    rate.b2 <- mu.b2/sd.b2^2