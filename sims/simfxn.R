balding_nichols <- function(n, f, p) {
  return(rbeta(n, (1-f)/f * p, (1-f)/f) * (1-p))
}

simfxn <- function(N, M, Fst, gamma, h2, tcor, ngen=1, Mcause=M/10) { 

  ps	<- runif( M, .05, .5 )
  freqs1 =  balding_nichols(M, Fst, ps) #allelel freqs in pop 1
  freqs2 =  balding_nichols(M, Fst, ps) #allele freqs in pop 2

  pop = c(rep(1,N),rep(2,N)) #population labels
  
  betas = rep(0,M)
  betas[sample(M,Mcause)]	<- sqrt(h2/(Mcause)) * rnorm(Mcause)
  
  # F-statistics is always about ratios of variances
  # FST - short for f-statistic
  # Ratio of in population variance to between population variance
  # FST = 0, doesn't mean same pop but same pop distribution
  
  # FST is some measure of distance arbitrarily determined
  # Only matters here when we're god but not in real pops

  #make genotypes in males and females of each population
  genosM1 = matrix(rbinom(N*M/2,2,freqs1),N/2,M,byrow=T) 
  genosF1 = matrix(rbinom(N*M/2,2,freqs1),N/2,M,byrow=T)
  genosM2 = matrix(rbinom(N*M/2,2,freqs2),N/2,M,byrow=T)
  genosF2 = matrix(rbinom(N*M/2,2,freqs2),N/2,M,byrow=T)

  #generate phenotypes with mean shift in pop 1
  phenosM1 = genosM1%*%betas + sqrt(1-h2) * rnorm(N/2) + gamma
  phenosF1 = genosF1%*%betas + sqrt(1-h2) * rnorm(N/2) + gamma
  phenosM2 = genosM2%*%betas + sqrt(1-h2) * rnorm(N/2)
  phenosF2 = genosF2%*%betas + sqrt(1-h2) * rnorm(N/2)
  
  for (cur in 1:ngen) {

    #make couples
    couple1	<- couple_fxn( phenosM1, phenosF1, tcor, N )
    couple2	<- couple_fxn( phenosM2, phenosF2, tcor, N )
    
    #make children
    genosC1	<- child_fxn( genosM1, genosF1, couple1$sM, couple1$sF, N, M )
    genosC2	<- child_fxn( genosM2, genosF2, couple2$sM, couple2$sF, N, M )
    
    #make phenotypes in children
    phenosC1 = genosC1%*%betas + rnorm(N,0,sqrt(0.5)) + gamma
    phenosC2 = genosC2%*%betas + rnorm(N,0,sqrt(0.5))
    
    m1 = sample(1:N,N/2)
    f1 = setdiff(1:N,m1) 
    m2 = sample(1:N,N/2)
    f2 = setdiff(1:N,m2) 
    
    genosM1 <- genosC1[m1 ,]
    genosF1 <- genosC1[f1 ,]
    genosM2 <- genosC2[m2 ,]
    genosF2 <- genosC2[f2 ,]
    
    phenosM1 <- phenosC1[m1]
    phenosF1 <- phenosC1[f1]
    phenosM2 <- phenosC2[m2]
    phenosF2 <- phenosC2[f2]
    
  }
  
  train = sample(N,N/2)
  test = setdiff(1:N,train)

  return(
    list( betas=betas,
        poptest=rep( 1:2, each=N/2 ),
        poptrain=rep( 1:2, each=N/2 ),
        ytest=c( phenosC1[test], phenosC2[test] ),
        ytrain=c( phenosC1[train], phenosC2[train] ),
        Gtest=rbind( genosC1[test,], genosC2[test,] ),
        Gtrain=rbind( genosC1[train,], genosC2[train,] ),
        ngen=ngen
    )
  )
}

couple_fxn	<- function( phenosM, phenosF, tcor, N ) {
  
  pcor = cor(phenosM,phenosF)
	iter = 1
	while(pcor < tcor){
	  # tcor is the target mate pair correlation
		sM = sort(phenosM+10/iter * rnorm(N/2),index.return=T)
		sF = sort(phenosF+10/iter * rnorm(N/2),index.return=T)
		
		pcor = cor(phenosM[sM$ix],phenosF[sF$ix])
		iter = iter+1

	}
	return(list( sM=sM, sF=sF ))
}

child_fxn	<- function( genosM, genosF, sM, sF, N, M ) {
	genosC = matrix(0,N,M)
	for(i in 1:(N/2)) {
		matpat = rbinom(M,1,0.5)
		genosC[i,matpat==1] = genosM[sM$ix[i],matpat==1]
		genosC[i,matpat==0] = genosF[sF$ix[i],matpat==0]
		matpat = rbinom(M,1,0.5)
		genosC[i+(N/2),matpat==1] = genosM[sM$ix[i],matpat==1]
		genosC[i+(N/2),matpat==0] = genosF[sF$ix[i],matpat==0]
	}
	return(genosC)
}
