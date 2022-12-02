# BINF 8205 Computational Molecular Evolution Final Code
# Authored by Daisy Fry Brumit
# Due date December 14 2022
# This code is meant to accompany finalProject.Rmd and is intended to hold the bulk
# of code necessary to generate results without cluttering the markdown intended for publishing.

###### Part 0: imports and setup ######
library(msa)
library(phangorn)
library(phytools)
library(bios2mds)
library(knitr)
library(tidyverse)



###### Part 1: Sequence Alignments ######

# list loci
loci.list = c('atp6', 'atp8', 'cox1','cox2', 'cox3', 'cytb', 'nd1', 'nd2', 'nd3', 'nd4',
              'nd4L', 'nd5', 'nd6')
# get a list of all fasta files in my fish_fasta folder
my.files = list.files(path="./fish_fasta", pattern="*.fasta", full.names=TRUE)

### get concatenated alignments for multiple aln methods to compare later
x = aln1
my.lengths = rep(0, length(my.files))
my.lengths[1] = ncol(x)

for (i in 2:length(my.files)) {
  seqs = readDNAStringSet(my.files[i])
  y = msa(seqs, method="ClustalOmega", order="input")
  my.lengths[i] = ncol(y)
  temp = paste0(x,y)
  names(temp) = rownames(y)
  x = DNAStringSet(temp)
}
concat.aln = DNAMultipleAlignment(x)
write.nexus.data(as.DNAbin(concat.aln), file="RosesConcat.nex")
