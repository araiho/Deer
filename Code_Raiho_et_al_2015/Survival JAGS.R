
model{
	
	#priors
	a~dgamma(.001,.001)
	b~dgamma(.001,.001)
	
	for(i in 1:length(y.m[])){
		s[i] ~ dbeta(a,b)
		shape1[i] <- max(.001,((s[i])^2-(s[i]^3)-s[i]*y.sd[i]^2)/(y.sd[i]^2))
		shape2[i] <- max(.001,(s[i]-(2*s[i]^2)+s[i]^3-y.sd[i]^2+s[i]*y.sd[i]^2)/(y.sd[i]^2))
		y.m[i] ~ dbeta(shape1[i],shape2[i])
		
		
	}
	fe <- a/(a+b)
}

