#N = 2e3 #individuals in each population (note testing is done in half the sample size)
#M = 3e3 #SNPs (note half are causal)

N = 5e3 #individuals in each population (note testing is done in half the sample size)
M = 1e4 #SNPs (note half are causal)

gamma	<- 0.1	#pop mean shift in phenotype
h2		<- 0.5	#heritability
tcor	<- 0.6	#true correlation of parents at phenotype level

Fsts	<- c(1e-4,.001,.01,.1)#,.5
nF		<- length(Fsts)

maxS	<- c(1, 5, 10, 50, 100, 500, 1000, 2500, 5000, 10000)
nS		<- length(maxS)

gens <- c(1, 2, 3) #num generations to look at

maxit	<- 200
boot_n <- 100