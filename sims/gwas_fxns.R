run_gwas	<- function(y,G,pop,pc,pc10){

	betas			<- apply( G, 2, function(g) safe_lm( summary( lm( y ~ 1+g			) )$coef ) )
	betas_adj	<- apply( G, 2, function(g) safe_lm( summary( lm( y ~ 1+g+pop	) )$coef ) )
	betas_pc	<- apply( G, 2, function(g) safe_lm( summary( lm( y ~ 1+g+pc	) )$coef ) )
	betas_pc10<- apply( G, 2, function(g) safe_lm( summary( lm( y ~ 1+g+pc10) )$coef ) )
	#betas_pc10<- NA

	list( betas=betas, betas_adj=betas_adj, betas_pc=betas_pc, betas_pc10=betas_pc10 )
}

corve_fxn	<- function(x1,x2,y,z){
	if( !missing( z ) ){
		x1	<- resid( lm( x1 ~ z ) )
		x2	<- resid( lm( x2 ~ z ) )
		y		<- resid( lm( y  ~ z ) )
	}
	c( cor(x1,x2), summary(lm(y~I(x1 + x2)))$adj.r.squared )
}

beta2prs	<- function( betas, pvals, maxS, G ){
	if( maxS < length( betas ) ){
		nonzeros	<- sort.list( pvals, dec=F )[1:maxS]
		betas[-nonzeros]	<- 0
	}
	G %*% betas
}

safe_lm	<- function( lm.coef ){
	out	<- c( 0, NA, NA, 1 )
	try( out	<- lm.coef['g',] )
	out
}

# add bootstrap so we can see the GLS plots
  # slope of the line is the mate pair correlation
    # show that the true mate pare correlation falls within CI 
    # Showing the estimator is unbiased

# Want to be able to say that the p-value rejection is "correct"

# Fraction explained by population structure
  # gotten from the sim parameters
