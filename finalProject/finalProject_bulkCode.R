# BINF 8205 Computational Molecular Evolution Final Code
# Authored by Daisy Fry Brumit
# Due date December 14 2022
# This code is meant to accompany finalProject.Rmd and is intended to hold the bulk
# of code necessary to generate results without cluttering the markdown intended for publishing.

###### Part 0: imports and setup ######
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("msa", force = TRUE)
BiocManager::install("phangorn")
BiocManager::install("phytools")
BiocManager::install("bios2mds", force = TRUE)
BiocManager::install("pegas")
BiocManager::install("geiger")

library(msa)
library(phangorn)
library(phytools)
library(bios2mds)
library(knitr)
library(tidyverse)
library(pegas)
library(castor)
library(evobiR)
library(MCMCtreeR)
library(geiger)

###### Part 1: Sequence Alignments ######

# list loci
loci.list = c('atp6', 'atp8', 'cox1','cox2', 'cox3', 'cytb', 'nd1', 'nd2', 'nd3', 'nd4',
              'nd4L', 'nd5', 'nd6')

# get a list of all fasta files in my fish_fasta folder
my.files = list.files(path="./fish_fasta", pattern="*.fasta", full.names=TRUE)
my.alns = vector("list", length(my.files))

# loop through my files and align using muscle, output as .fa file for VerAlign
for (i in 1:length(my.files)) {
  x = readDNAStringSet(my.files[i])
  my.alns[[i]] = msa(x, method="Muscle", order="input") # perform alignment
  
  # did export fasta for all msa methods. No longer needed. 
  #tmp_aln_for_fasta = msaConvert(my.alns[[i]], "bios2mds::align") # convert to format 
  #export.fasta(tmp_aln_for_fasta, outfile=paste0(loci.list[i], "_muscle_aln.fa")) # convert to fasta file
}

###### Part 2: Test Clockwise and Substitution Models ######
# determine which loci are violating clockwise assumptions
num.nonclock = rep(0, length(my.alns))

for (i in 1:length(my.alns)) {
  temp = my.alns[[i]]
  dna = as.DNAbin(temp)
  outindex = grep("japonica", labels(dna))
  outgroup = dna[outindex,]
  FUN <- function(x)
    rr.test(dna[x[1], ], dna[x[2], ], outgroup)$Pval
  n = nrow(dna)
  cc = combn(2:n, 2)
  OUT <- apply(cc, 2, FUN)
  num.nonclock[i] = length(which(OUT<(0.05/(length(cc)))))
}
# which loci have 0 clock-wise violations
clock.like = which(num.nonclock==0) # none! (CO, CW, or muscle)
clock.like.aln = my.alns[clock.like] # none!

# test which substitution models fit which loci
my.sub.models = rep(NA, length(my.alns))

for (i in 1:length(my.alns)) {
  temp = my.alns[[i]]
  phang = msaConvert(temp, type="phangorn::phyDat")
  mt = modelTest(phang)
  my.sub.models[i] = mt[which(mt$BIC == min(mt$BIC)),1]
}

table(my.sub.models, loci.list)


###### Part 3: Concatenate for ML and Bayes Treebuilding ######
# concatenate for mrbayes
x = my.alns[[1]]
my.lengths = rep(0, length(my.alns))
my.lengths[1] = ncol(x)
for (i in 2:length(my.alns)) {
  y = my.alns[[i]]
  my.lengths[i] = ncol(y)
  temp = paste0(x,y)
  names(temp) = rownames(y)
  x = DNAStringSet(temp)
}
concat.aln = DNAMultipleAlignment(x)

write.nexus.data(as.DNAbin(concat.aln), file="FishConcat.nex") # mr bayes
write.phylip(concat.aln, "FishConcat.phy")

# get partition file (for RAxML)
my.starts = c(1)
for (i in 1:(length(my.lengths)-1)) {
  my.starts[i+1] = my.lengths[i] + my.starts[i]
}
my.ends = my.starts + (my.lengths-1)

info1 = paste(my.starts, my.ends, sep="-")
info2 = paste(loci.list, info1, sep="=")
info3 = paste("DNA,", info2, sep=" ")
write.table(info3, file="Fish-Raxml-part.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

### Create the partition file for MrBayes
info3 = gsub("DNA,", "charset", info3)
info3 = paste(info3, ";", sep="")
loc.list = paste(loci.list, collapse=",")
list.w.length = paste(c(length(loci.list), loc.list), collapse=":")
info4 = paste(c("partition", "byLocus", "=", list.w.length), collapse=" ")
info4 = paste(info4, ";", sep="")
write.table(c(info3,info4), file="Fish-MrB-part.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

###### Part 4: Consensus trees ######
#start with ML output and check phylogeny
my.tree = read.tree('treeBuilds_raxml/RAxML_bestTree.fishCon')
my.boot = read.tree('treeBuilds_raxml/RAxML_bootstrap.fishCon')

plotBS(tree=my.tree, BStrees=my.boot, cex=0.75, type="phylogram", p=0, bs.col="blue") ## doesn't look quite the way I want it to

## Bayesian
bayesTree = read.nexus("treeBuilds_mrBayes/FishConcat.nex.con.tre")
bayesTree1 = bayesTree[[1]] # contains branch length info
bayesTree2 = bayesTree[[2]] # only contains topology

bayesTree1$node.label = round(((as.numeric(bayesTree1$node.label))*100), digits=2)
bayesTree1 = root(bayesTree1, "japonica", resolve.root=T)

plotTree(bayesTree1, edge.width=2, font=1)
nodelabels(bayesTree2$node.label, cex=0.8)

###### Part 5: individual gene trees & ASTRAL ######
# Align each sequence one by one, output as .nexus file
for (i in 1:length(my.files)) {
  seqs = readDNAStringSet(my.files[i])
  alns = msa(seqs, method="Muscle", order="input")
  out = gsub("fasta", "nexus", my.files[i])
  write.nexus.data(as.DNAbin(alns), out)
}

# Download all of the run1.t files to get mcc trees
# Read all of this files into a loop,
# get the rooted MCC tree, and write to a newick file
# Get all of the max clade credibility trees
tree.files = list.files(path="./single-gene-trees", pattern="*run1.t", full.names = TRUE)
for (i in 1:length(tree.files)) {
  out = gsub("nexus\\.run1\\.t", "mcc.tre", tree.files[i])
  trees = read.nexus(tree.files[i])
  mcctree = root(mcc(trees), outgroup=c("japonica"))
  write.tree(mcctree, out)
}

### Get all of the max clade credibility trees
tree.files = list.files(path="./single-gene-treesets", pattern="*run1.t", full.names = TRUE)
for (i in 1:length(tree.files)) {
  out = gsub("nex\\.run1\\.t", "mcc.tre", gsub("single-gene-treesets\\/", "mcc_trees\\/", tree.files[i]))
  trees = read.nexus(tree.files[i])
  mcctree = root(mcc(trees), outgroup=c("japonica"))
  write.tree(mcctree, out)
}

# run in astral on cluster
# ...

# output final tree
astral = read.tree("fish_astral.tre")
astral$node.label = round(((as.numeric(astral$node.label))*100), digits=2)
astral$edge.length[which(is.na(astral$edge.length))] = 1
astral2 = root(astral, "japonica", resolve.root=TRUE)
jpeg("Fish-Astral-Tree.jpg", width=8, height=8, units="in", res=300)
plotTree(astral2, edge.width=2, font=1)
nodelabels(astral2$node.label, cex=0.8)
dev.off()

#get density tree 

all.trees = read.tree("mcc_trees/fish_allLoci_mcc.tre")
y = root(all.trees, outgroup=c("japonica"), resolve.root=TRUE)
densityTree(y, type = "cladogram", alpha = 0.05, compute.consensus=FALSE, use.edge.length=FALSE)

#compare trees, please!
treedist(my.tree,bayesTree1) # differ at 1 site
treedist(my.tree, astral2) # differ at 1 site
treedist(bayesTree1,astral2) # differ at 2 sites

###### Part 6: ASR ######
# list out species
TreeID = c('japonica','americanus','setigerus','caulinaris','ocellatum','pictus','brevicaudata',
           'jordani','pelagica','spinifer','murrayi','thompsoni','vanhoeffeni','macrothrix','mollis')
maxDepth = c(628,800,800,311,51,978,1024,1510,2500,1200,6370,2014,5300,2000,2250)
names(maxDepth)=TreeID
m = match(my.tree$tip.label, names(maxDepth))
maxDepth = maxDepth[m]
contMap(my.tree, maxDepth)
