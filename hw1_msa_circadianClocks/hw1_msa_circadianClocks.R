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
BiocManager::install("phangorn", force = TRUE)
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


# Part 2: Compare Alignments with a Reference Alignment

aligns1_as_align = msaConvert(aligns1, "bios2mds::align") #generate clustalW alignment
export.fasta(aligns1_as_align, outfile="clustalW_aln.fa") #output evaluated alignment
aligns2_as_align = msaConvert(aligns2, "bios2mds::align") #generate clustal omega alignment
export.fasta(aligns2_as_align, outfile="clustalO_aln.fa") 

# Part 3: ID ortho and paralogs

aaPhydat = as.phyDat(aligns1)
d = dist.ml(aaPhydat, model="JTT")
my.tree = nj(d)
plot(my.tree)

my.tree2 = root(my.tree, outgroup="InsectJ")
plot(my.tree2)

# Part 4: examine diffs

cry1.seqs = aaSeqs[grep("(FishC|FishD|MouseA|BirdF|FrogH)", names(aaSeqs))]
cry2.seqs = aaSeqs[grep("(BirdG|MouseB|FrogI|FishE)", names(aaSeqs))]

cry1_align = msa(cry1.seqs) #generate clustalW alignment
cry2_align = msa(cry2.seqs) #generate clustalW alignment


msaPrettyPrint(cry1_align,
               output="pdf",
               file="cry1_align.pdf",
               showNames="right",
               showNumbering="none",
               shadingMode="functional",
               shadingModeArg="structure",
               verbose=TRUE,
               askForOverwrite=TRUE)

msaPrettyPrint(cry2_align,
               output="pdf",
               file="cry2_align.pdf",
               showNames="right",
               showNumbering="none",
               shadingMode="functional",
               shadingModeArg="structure",
               verbose=TRUE,
               askForOverwrite=TRUE)

x = msaConsensusSequence(cry1_align)
y = msaConsensusSequence(cry2_align)
cry.con = AAStringSet(c(x,y))
cry_align = msa(cry.con)

msaPrettyPrint(cry_align, 
               output="pdf", 
               file="cry_align.pdf", 
               showNames="right",
               showNumbering="none", 
               showLogo="none", 
               askForOverwrite=TRUE, 
               verbose=TRUE)
