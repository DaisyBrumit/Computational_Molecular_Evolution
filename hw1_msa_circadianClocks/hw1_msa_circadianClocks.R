# BINF 8205 ASSIGNMENT 1
# 1-MSA and Ortholog Identification in Circadian Clock Genes 
# DUE SEPT 5 2022
# DAISY FRY BRUMIT

# data input: 'cryptochrome.fa'
# data output: 'clustalW_align.pdf' --> results of initial alignment
#              ''

# install and access relevant packages
install.packages('knitr')
install.packages("tinytex")

tinytex::install_tinytex()
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")
BiocManager::install("phangorn")
BiocManager::install("phytools")
BiocManager::install("bios2mds")

library(msa)
library(phangorn)
library(phytools)
library(bios2mds)
library(knitr)
library(tinytex)

# DON'T FORGET TO SETWD!

# Part 1: Align the Sequences with Different Methods
aaSeqs = readAAStringSet("cryptochrome.fa") # read in
aligns1 = msa(aaSeqs) # alignment, default = Clustal W
aligns2 = msa(aaSeqs, method = "ClustalOmega")

#print(aligns1, show="complete") # check work
# create more visually appealing output
msaPrettyPrint(aligns1, 
               output="pdf", 
               file="clustalW_align.pdf", 
               showNames="right",
               showNumbering="none", 
               showLogo="none", 
               askForOverwrite=TRUE, 
               verbose=TRUE)

msaPrettyPrint(aligns2, 
               output="pdf", 
               file="clustalOmega_align.pdf", 
               showNames="right",
               showNumbering="none", 
               showLogo="none", 
               askForOverwrite=TRUE, 
               verbose=TRUE)


