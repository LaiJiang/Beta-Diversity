# Box 4: Estimate PERMANOVA power of a distance matrix using bootstrap sampling
 
N.groups <- 2	#number of simulated groups
N.perm.sub.group <- 10 #number of simulated subjects per group
 
#number of bootstrap samples for each Distance matrix (DM)
N.boots <- 100
 
typeI.error <- 0.05 #criteria for controlled type-I error
#################################################################
library(micropower)
library(phyloseq)
library(vegan)
#load one single DM (previously simulated) either with Benchmark method (Box2):
# Load ("SimDM_Bechmark.RData", verbose=TRUE); sim.DM = Sim.DM.list[[1]]
# Or with Kelly’s method (Box3):
# Load ("DMs.RData", verbose=TRUE); sim.DM = dist.mat
##########################################################################################
#function to calculate effect size
effect.size.cal <- function(dist.mat, group.id){
  N <- nrow(dist.mat)
  a <- length(unique(group.id))
  list.groupID <- unique(group.id)
  SSw <- 0
  for(temp.id in list.groupID){
	dist.group <- dist.mat[group.id==temp.id,group.id==temp.id] 
	ss.group <- sum(dist.group^2)/2/nrow(dist.group)
	SSw <- SSw + ss.group
	#print(ss.group)
  }
  SST <- sum(dist.mat^2)/2/N
  SSa <- SST – SSw             
  effect.size <- (SSa - (a-1)*(SSw/(N-a)) )/ (SST + SSw/(N-a))
  effect.size
}
###########################################################################################
  rec.pvalues <- NULL   #start loop for bootstrap samplings
  for(i.perm in 1:N.boots){	boots.sample.id <- as.vector(sapply(c(1:N.groups),
 	function(x){sample(which(sim.group.id==x), size=N.perm.sub.group, replace=TRUE)}))
 
	#the bootstrapped DM and subject group ID
	boots.sample.DM <- as.matrix(sim.DM)[boots.sample.id,boots.sample.id]
	boots.group.id  <- sim.group.id[boots.sample.id]
 
	#apply permanova to this bootstrapped DM
	permanova = adonis(formula=boots.sample.DM ~ boots.group.id, permutations=999, method='bray')
	#pvalue from PERMANOVA
	boots.pvalues <- permanova$aov.tab$`Pr(>F)`[1]
	rec.pvalues <- c(rec.pvalues,boots.pvalues)  #record
  }
  #################################################################
  #power: the proportion of bootstrap distance matrices for which
  #PERMANOVA P-values are less than the pre-specified threshold for type I error.
  boots.power <- mean(rec.pvalues<typeI.error)
  effect.size <- effect.size.cal(as.matrix(sim.DM),sim.group.id) #”true” effect size
  #record the power and the effect size
  save(boots.power, effect.size, file=”power_example.RData”)
