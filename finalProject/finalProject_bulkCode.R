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
<<<<<<< HEAD
  tmp_aln_for_fasta = msaConvert(my.alns[[i]], "bios2mds::align") # convert to format 
  export.fasta(tmp_aln_for_fasta, outfile=paste0(loci.list[i], "_muscle_aln.fa")) # convert to fasta file
}

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
clock.like = which(num.nonclock==0) # none!
clock.like.aln = my.alns[clock.like]

=======
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
write.nexus.data(as.DNAbin(concat.aln), file="FishConcat.nex")

# get partition file (for RAxML)
my.starts = c(1)
for (i in 1:(length(my.lengths)-1)) {
  my.starts[i+1] = my.lengths[i] + my.starts[i]
}
my.ends = my.starts + (my.lengths-1)

info1 = paste(my.starts, my.ends, sep="-")
info2 = paste(loci.list, info1, sep="=")
info3 = paste("DNA,", info2, sep=" ")

### Create the partition file for MrBayes
info3 = gsub("DNA,", "charset", info3)
info3 = paste(info3, ";", sep="")
loc.list = paste(loci.list, collapse=",")
list.w.length = paste(c(length(loci.list), loc.list), collapse=":")
info4 = paste(c("partition", "byLocus", "=", list.w.length), collapse=" ")
info4 = paste(info4, ";", sep="")
write.table(c(info3,info4), file="Fish-MrB-part.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
>>>>>>> ecc829586540db04b459c425a23a10858f74d235
