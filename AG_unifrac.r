# Box 1 : sample code for reading in American Gut data and estimating a Unifrac distance
# Reading in and extracting OTU tables for American Gut data
# https://github.com/biocore/American-Gut
library(phyloseq)
 
PATH.main <- "yourpathtodata/"
PATH.dat.trim <- paste0(PATH.main, "dat/03-otus/100nt/gg-13_8-97-percent/")
PATH.save <- paste0(PATH.main, "save/")
 
#Load tree file
tree_file = read_tree(paste0(PATH.dat.trim, "97_otus.tree"))
 
#load the selected biom data
biom_file1 = import_biom(paste0(PATH.main,”/AG.biom"))
 
#Read meta data for selection (requires read.delim – tab separated)
#define missing values to exclude
dat1 = read.delim(paste0((PATH.main,"AG.txt"), sep="\t", na.strings = c("no_data","unknown"))
 
#select samples from two small sample types with 165 or 150 samples for illustration
sel.sample.type = which((dat1$BODY_SITE %in% c("UBERON:hand", "UBERON:skin")))
 
#Extracting OTU tables from the phylum rank
phyloseq_phylum = tax_glom(biom_file1, "Rank2")
otu_tab_phylum = otu_table(phyloseq_phylum)
 
#select desired samples under phylum
sel_otu_tab_phylum = otu_tab_phylum[,sel.sample.type]
 
#combine file
combine_file_plylum =  phyloseq(sel_otu_tab_phylum, phy_tree(tree_file))
 
#calculate the unifac distance matrix, and save the result
dist.mat <- distance(combine_file_plylum, "uUniFrac")
save(sel.sample.type, dist.mat, file=paste0(PATH.save, "uUniFrac.RData"))
