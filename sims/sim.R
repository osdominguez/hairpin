rm( list=ls() )

setwd("/gpfs/data/ukb-share/dahl/ophelia/hairpin/sims")

source( 'simfxn.R' )
source( 'gwas_fxns.R' )
source( 'setup.R' )

for (it in sample(maxit)) {
    #boot_df <- data.frame(matrix(ncol=6, nrow=0, dimnames=list(NULL, c("Fst", "SNPS", "pc_num", "r2", "theta_eo", "replicate"))))
	#boot_save <- paste0( 'boots/hairpins_', N, '_', M, '_', gamma, '_', gen, '_', it, '_bootstrap.table' ) 
	for ( F_i in 1:nF ) {
    	savefile	<- paste0( 'Rdata/hairpins_', N, '_', M, '_', Fsts[F_i], '_', gamma, '_', it, '.Rdata' )
    	sinkfile	<- paste0(  'Rout/hairpins_', N, '_', M, '_', Fsts[F_i], '_', gamma, '_', it, '.Rout' )
    	
    	if ( file.exists( savefile ) | file.exists( sinkfile ) ) { next }
    	
    	sink( sinkfile )
    
    	ev	<- 1:(M/2)
    	od	<- 1:(M/2)+M/2
    
    	corves	<- array( NA, dim=c(2,6,nS), dimnames=list( c('cor','ve'), c( 'true', 'estraw', 'estpc', 'estpc4', 'estpc10', 'estorac' ), maxS ) )
    	set.seed( it )
    	dat		<- simfxn( N, M, Fsts[F_i], gamma, h2, tcor, 1)
    	
    	pc_train <- svd(dat$Gtrain)$u
    	pctrain	<- pc_train[,1]
		pctrain4 <- pc_train[,4]
    	pctrain10 <- pc_train[,1:10]
    	rm(pc_train)
    	
    	pc_test <- svd(dat$Gtest)$u
		pctest <- pc_test[,1]
		pctest4 <- pc_test[,4]
		pctest10 <- pc_test[,1:10]
    	rm(pc_test)
    	
    	print( summary(lm( dat$poptrain ~ pctrain ))$adj.r.squared )
    	print( summary(lm( dat$poptrain ~ pctrain10 ))$adj.r.squared )
    
    	gwas <- run_gwas(dat$ytrain,dat$Gtrain,dat$poptrain,pctrain,pctrain4,pctrain10)
    	for( s_i in 1:nS ) {
    		print( s_i )
    
    		prsE	<- beta2prs( dat$betas[ev], -abs(dat$betas[ev]), maxS[s_i], dat$Gtest[,ev] )
    		prsO	<- beta2prs( dat$betas[od], -abs(dat$betas[od]), maxS[s_i], dat$Gtest[,od] )
    		corves[,'true',s_i]	<- corve_fxn( prsE, prsO, dat$ytest )
    
    		prsE	<- beta2prs( gwas$betas[1,ev], gwas$betas[4,ev], maxS[s_i], dat$Gtest[,ev] )
    		prsO	<- beta2prs( gwas$betas[1,od], gwas$betas[4,od], maxS[s_i], dat$Gtest[,od] )
    		corves[,'estraw',s_i]	<- corve_fxn( prsE, prsO, dat$ytest )
    
    		prsE	<- beta2prs( gwas$betas_pc [1,ev], gwas$betas_pc [4,ev], maxS[s_i], dat$Gtest[,ev] )
    		prsO	<- beta2prs( gwas$betas_pc [1,od], gwas$betas_pc [4,od], maxS[s_i], dat$Gtest[,od] )
    		corves[,'estpc',s_i]	<- corve_fxn( prsE, prsO, dat$ytest, z=pctest )
    
			prsE	<- beta2prs( gwas$betas_pc4[1,ev], gwas$betas_pc4[4,ev],maxS[s_i], dat$Gtest[,ev] )
    		prsO	<- beta2prs( gwas$betas_pc4[1,od], gwas$betas_pc4[4,od],maxS[s_i], dat$Gtest[,od] )
    		corves[,'estpc4',s_i]	<- corve_fxn( prsE, prsO, dat$ytest, z=pctest4 )

    		prsE	<- beta2prs( gwas$betas_pc10[1,ev], gwas$betas_pc10[4,ev],maxS[s_i], dat$Gtest[,ev] )
    		prsO	<- beta2prs( gwas$betas_pc10[1,od], gwas$betas_pc10[4,od],maxS[s_i], dat$Gtest[,od] )
    		corves[,'estpc10',s_i]	<- corve_fxn( prsE, prsO, dat$ytest, z=pctest10 )
    
    		prsE	<- beta2prs( gwas$betas_adj[1,ev], gwas$betas_adj[4,ev], maxS[s_i], dat$Gtest[,ev] )
    		prsO	<- beta2prs( gwas$betas_adj[1,od], gwas$betas_adj[4,od], maxS[s_i], dat$Gtest[,od] )
    		corves[,'estorac',s_i]	<- corve_fxn( prsE, prsO, dat$ytest, z=dat$poptest )

			#for (b_i in 1:boot_n) {
		#		t_df <- data.frame(matrix(ncol=6, nrow=1, dimnames=list(NULL, c("Fst", "SNPS", "pc_num", "r2", "theta_eo", "replicate")))) 
		#		t_df$replicate <- b_i
		#		t_df$SNPS <- maxS[s_i]
		#		t_df$Fst <- Fsts[F_i]
		#		
		#		ev_t <- sample(ev, length(ev), replace=TRUE)
		#		od_t <- sample(od, length(od), replace=TRUE)
		#		
		#		t_df$pc_num <- "unadj"
		#		prsE	<- beta2prs( dat$betas[ev_t], -abs(dat$betas[ev_t])	, maxS[s_i], dat$Gtest[,ev_t] )
    	#		prsO	<- beta2prs( dat$betas[od_t], -abs(dat$betas[od_t])	, maxS[s_i], dat$Gtest[,od_t] )
    	#		cvs	<- corve_fxn( prsE, prsO, dat$ytest )
		#		
		#		teo <- cvs[1]
		#		r2 <- cvs[2]
		#		t_df$theta_eo <- teo
		#		t_df$r2 <- r2
#
#				boot_df <- rbind(boot_df, t_df)
#
#				t_df$pc_num <- "0"
#				prsE	<- beta2prs( gwas$betas[1,ev_t], gwas$betas[4,ev_t], maxS[s_i], dat$Gtest[,ev_t] )
 #   			prsO	<- beta2prs( gwas$betas[1,od_t], gwas$betas[4,od_t], maxS[s_i], dat$Gtest[,od_t] )
  #  			cvs	<- corve_fxn( prsE, prsO, dat$ytest )
	#			
	#			teo <- cvs[1]
	#			r2 <- cvs[2]
	#			t_df$theta_eo <- teo
	#			t_df$r2 <- r2

	#			boot_df <- rbind(boot_df, t_df)
#
#				t_df$pc_num <- "1"
#				prsE	<- beta2prs( gwas$betas_pc [1,ev_t]	, gwas$betas_pc [4,ev_t], maxS[s_i], dat$Gtest[,ev_t] )
 #   			prsO	<- beta2prs( gwas$betas_pc [1,od_t]	, gwas$betas_pc [4,od_t], maxS[s_i], dat$Gtest[,od_t] )
  #  			cvs	<- corve_fxn( prsE, prsO, dat$ytest )
	#			
	#			teo <- cvs[1]
	#			r2 <- cvs[2]
	#			t_df$theta_eo <- teo
	#			t_df$r2 <- r2

	#			boot_df <- rbind(boot_df, t_df)

	#			t_df$pc_num <- "4"
	#			prsE	<- beta2prs( gwas$betas_pc4[1,ev_t], gwas$betas_pc4[4,ev_t],maxS[s_i], dat$Gtest[,ev_t] )
    #			prsO	<- beta2prs( gwas$betas_pc4[1,od_t], gwas$betas_pc4[4,od_t],maxS[s_i], dat$Gtest[,od_t] )
    #			cvs	<- corve_fxn( prsE, prsO, dat$ytest )
				
	#			teo <- cvs[1]
	#			r2 <- cvs[2]
	#			t_df$theta_eo <- teo
	#			t_df$r2 <- r2

	#			boot_df <- rbind(boot_df, t_df)

	#			t_df$pc_num <- "10"
    #			prsE	<- beta2prs( gwas$betas_pc10[1,ev_t], gwas$betas_pc10[4,ev_t],maxS[s_i], dat$Gtest[,ev_t] )
    #			prsO	<- beta2prs( gwas$betas_pc10[1,od_t], gwas$betas_pc10[4,od_t],maxS[s_i], dat$Gtest[,od_t] )
    #			cvs	<- corve_fxn( prsE, prsO, dat$ytest )
				
	#			teo <- cvs[1]
	#			r2 <- cvs[2]
	#			t_df$theta_eo <- teo
	#			t_df$r2 <- r2

	#			boot_df <- rbind(boot_df, t_df)

	#			t_df$pc_num <- "true"
	#			prsE	<- beta2prs( gwas$betas_adj[1,ev_t]	, gwas$betas_adj[4,ev_t], maxS[s_i], dat$Gtest[,ev_t] )
    #			prsO	<- beta2prs( gwas$betas_adj[1,od_t]	, gwas$betas_adj[4,od_t], maxS[s_i], dat$Gtest[,od_t] )
    #			cvs	<- corve_fxn( prsE, prsO, dat$ytest )
				
	#			teo <- cvs[1]
	#			r2 <- cvs[2]
	#			t_df$theta_eo <- teo
	#			t_df$r2 <- r2

	#			boot_df <- rbind(boot_df, t_df)

	#		}
	#		rm(t_df)

    	}
  	save( corves, file=savefile )
  	sink()
    }
	#write.table(boot_df, file = boot_save, row.names = F, quote = F)
	#rm(boot_df)
  }
