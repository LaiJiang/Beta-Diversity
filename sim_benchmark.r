# Box 2: Simulate a Unifrac Distance matrix based on benchmark distance matrix from
# American Gut data with a scaling parameter sigma that adjust effects sizes
 
N_subjects <- 500 ; N_groups <- 10 ; Group_inds <- N_subjects/N_groups
Sigma.list<- seq(from = 1, to =10, by = 0.1) #candidate values for scaling parameter Sigma
load(paste0(PATH.load, "uUniFrac.RData"),verbose=TRUE) #load benchmark matrix from AG data
#the loaded Robject of benchmark distance matrix:  dist.matrix
sim.group   <- rep(c(1:N_groups), N_subjects/N_groups)  #simulated group ID
 
Func_simDM <-function(sigma){  #sigma: scaling parameter that boosts between-group distances
  #obtain nonzero density values
  density.vals <- as.vector(dist.matrix);
  density.vals.nonzero <- density.vals[density.vals!=0]
  adj.nonzero <- (density.vals.nonzero)*sigma	# scaled distance values
  #first simulate between group distance
  sim.between.group.distance <- sample(adj.nonzero,
                                   	size=N_subjects*N_subjects, replace=TRUE)
  Sim.DM <- matrix(sim.between.group.distance,N_subjects,N_subjects)
  #simulate within group distance using the original (unscaled) distance values
  for(Gp.id in c(1:N_groups)){
	temp.dist <- sample(density.vals.nonzero, size=Group_inds*Group_inds, replace=TRUE)
	temp.dist1 <- matrix(temp.dist, Group_inds,Group_inds );
	temp.dist1[lower.tri(temp.dist1)] <- 0
	temp.matrix <- (temp.dist1 + t(temp.dist1))/2; diag(temp.matrix) = 0
	Sim.DM[sim.group==Gp.id, sim.group==Gp.id] <- temp.matrix#within-grup distances}
	Sim.DM   #output: simulated Density matrix}
 
Sim.DM.list <-  lapply(Func_simDM, Sigma.list) #the simulated distance matrices
Save(Sigma.list, N_subjects, N_groups, sim.group, Sigma.list,  Sim.DM.list, file=”SimDM_Bechmark.RData”) #save the list of simulated DMs, sigmas, and parameters.
