# Box 3: Simulate a Unifrac Distance matrix with pre-specified within-group distance level
# and effect sizes using techniques from Kelly (2015)
 
library(micropower)  #Kelly paper
library(phyloseq)  
 
#specify population parameters
N.total <- 100    	#total number of subjects
N.groups <- 2     	#number of groups
otu.number <- 20  	#number of OTUs
sequence_depth <- 10  #number of sequence counts per OTU bin
 
rare_depth <- 0.2 	#proportion of sequence counts to retain after subsampling step.
                  	#this parameter corresponds to average within-group distances
                  	#higher values corresponds to close/low distance
 
effect  <- 0.2    	#the parameter value for the proportions of unique community
                  	#membership in membership segregation step (see Kelly 2015).
                  	#this parameter corresponds to the ``true” effect sizes
 
Distance.metric <- “Jaccard” 	#user can choose to create Jaccard distance matrix
Distance.metric <- “Unifrac” 	#or Unifrac distance
#number of subjects in each group
group_size_vector <- rep(N.total/N.groups,N.groups)
#simulated group ID for each individual in N.total
sim.group.id   <- rep(c(1:N.groups),each=N.total/N.groups) #simulated group ID
 
#simulate OTU tables using the method from Kelly paper
#rarefaction and segragation is embedded in their function simStudy
  Sim.OTUs <- simStudy( group_size_vector, otu.number, sequence_depth,
                    	rare_depth, effect)
 
#no need to simulate tree file for Jaccard distance
  if(Distance.metric==”Jaccard”)dist.mat <- phyloseq::distance(otu_table(Sim.OTUs, taxa_are_rows = TRUE), "jaccard")
 
  if(Distance.metric==”Unifrac”){
  #we need to simulate a phylogenetic tree for unifrac distance
  Sim.tree <- simTreeTable(Sim.OTUs)
 
  #calculate the Unifrac Distance
  combine_file_samples =  phyloseq::phyloseq( otu_table(Sim.OTUs, taxa_are_rows = TRUE), Sim.tree)
  dist.mat <- phyloseq::distance(combine_file_samples, "jaccard")
}
 
#save distance matrix, group IDs 
save(dist.mat, sim.group.id, file= "DMs.RData")
