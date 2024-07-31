rm( list=ls() )

setwd("/gpfs/data/ukb-share/dahl/ophelia/hairpin/sims")

source( 'simfxn.R' )
source( 'gwas_fxns.R' )
source( 'setup.R' )

for (it in sample(maxit)) {
  for (gen in gens) {
    for ( F_i in 1:nF ) {
    	savefile	<- paste0( 'Rdata/hairpins_', N, '_', M, '_', Fsts[F_i], '_', gamma, '_', gen, '_', it, '.Rdata' )
    	sinkfile	<- paste0(  'Rout/hairpins_', N, '_', M, '_', Fsts[F_i], '_', gamma, '_', gen, '_', it, '.Rout' )
    	
    	if ( file.exists( savefile ) | file.exists( sinkfile ) ) { next }
    	
    	sink( sinkfile )
    
    	ev	<- 1:(M/2)
    	od	<- 1:(M/2)+M/2
    
    	corves	<- array( NA, dim=c(2,5,nS), dimnames=list( c('cor','ve'), c( 'true', 'estraw', 'estpc', 'estpc10', 'estorac' ), maxS ) )
    	set.seed( it )
    	dat		<- simfxn( N, M, Fsts[F_i], gamma, h2, tcor, gen )
    	print('returned')
      
    	
    	pc_train <- svd(dat$Gtrain)$u
    	pctrain		<- pc_train[,1]
    	pctrain10	<- pc_train[,1:10]
      rm(pc_train)
    	
      pc_test <- svd(dat$Gtest)$u
    	pctest		<- pc_test[,1]
    	pctest10	<- pc_test[,1:10]
      rm(pc_test)
    	
    	print( summary(lm( dat$poptrain ~ pctrain ))$adj.r.squared )
    	print( summary(lm( dat$poptrain ~ pctrain10 ))$adj.r.squared )
    
    	gwas	<- run_gwas(dat$ytrain,dat$Gtrain,dat$poptrain,pctrain,pctrain10)
    
    	for( s_i in 1:nS ) {
    		print( s_i )
    
    		prsE	<- beta2prs( dat$betas[ev]				, -abs(dat$betas[ev])	, maxS[s_i], dat$Gtest[,ev] )
    		prsO	<- beta2prs( dat$betas[od]				, -abs(dat$betas[od])	, maxS[s_i], dat$Gtest[,od] )
    		corves[,'true'		,s_i]	<- corve_fxn( prsE, prsO, dat$ytest )
    
    		prsE	<- beta2prs( gwas$betas[1,ev]			, gwas$betas[4,ev]		, maxS[s_i], dat$Gtest[,ev] )
    		prsO	<- beta2prs( gwas$betas[1,od]			, gwas$betas[4,od]		, maxS[s_i], dat$Gtest[,od] )
    		corves[,'estraw'		,s_i]	<- corve_fxn( prsE, prsO, dat$ytest )
    
    		prsE	<- beta2prs( gwas$betas_pc [1,ev]	, gwas$betas_pc [4,ev], maxS[s_i], dat$Gtest[,ev] )
    		prsO	<- beta2prs( gwas$betas_pc [1,od]	, gwas$betas_pc [4,od], maxS[s_i], dat$Gtest[,od] )
    		corves[,'estpc'		,s_i]	<- corve_fxn( prsE, prsO, dat$ytest, z=pctest )
    
    		prsE	<- beta2prs( gwas$betas_pc10[1,ev], gwas$betas_pc10[4,ev],maxS[s_i], dat$Gtest[,ev] )
    		prsO	<- beta2prs( gwas$betas_pc10[1,od], gwas$betas_pc10[4,od],maxS[s_i], dat$Gtest[,od] )
    		corves[,'estpc10'	,s_i]	<- corve_fxn( prsE, prsO, dat$ytest, z=pctest10 )
    
    		prsE	<- beta2prs( gwas$betas_adj[1,ev]	, gwas$betas_adj[4,ev], maxS[s_i], dat$Gtest[,ev] )
    		prsO	<- beta2prs( gwas$betas_adj[1,od]	, gwas$betas_adj[4,od], maxS[s_i], dat$Gtest[,od] )
    		corves[,'estorac'	,s_i]	<- corve_fxn( prsE, prsO, dat$ytest, z=dat$poptest )
    	}
  	save( corves, file=savefile )
  	sink()
    }
  }
}
