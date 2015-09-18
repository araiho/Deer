model{

	for(i in 1:length(y.m)){
		shape1[i] <- max(.001,((y.m[i])^2-(y.m[i]^3)-y.m[i]*y.sd[i]^2)/(y.sd[i]^2))
		shape2[i] <- max(.001,(y.m[i]-(2*y.m[i]^2)+y.m[i]^3-y.sd[i]^2+y.m[i]*y.sd[i]^2)/(y.sd[i]^2))
		mu[i] ~ dbeta(shape1[i],shape2[i])	
	}
	
	mu.true <- sum(mu)/5
	
}