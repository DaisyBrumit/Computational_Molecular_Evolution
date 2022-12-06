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

