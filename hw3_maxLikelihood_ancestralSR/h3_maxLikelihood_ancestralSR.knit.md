---
title: 'Homework 3: Max Likelihood & Ancestral State Reconstruction'
author: "Daisy Fry-Brumit"
date: "09/26/2022"
output:
  word_document: default
  pdf_document: default
---



## Background

For this exercise, you will use maximum likelihood and consensus tree methods to look at the phylogenetic relationships among several species of snakes, lizards, and these legless lizard intermediates to determine whether or not the evolution of limblessness has happened more than once.

The dataset for this exercise consists of 2 sets of sequences corresponding to 2 different genes: SLC8A1^1^ and R35^2^. Use the links associated with each gene name to download the fasta files. Each file contains DNA sequences from 12 different reptile species.

## Use Neighbor-Joining Trees to Determine Which Gene To Use



### Question 1 (2pts)

#####  Provide the plot for your R35 neighbor-joining tree. Make sure that the tip labels are clearly readable.


![](h3_maxLikelihood_ancestralSR_files/figure-docx/plot_R35-1.png)

### Question 2 (2pts)

##### Provide the plot for you SLC8A1 tree; make sure the labels are readable. 


![](h3_maxLikelihood_ancestralSR_files/figure-docx/plot_SLC8A1-1.png)

### Question 3 (1pt)

#####  How many clades (or bipartitions) are unique to EACH tree? Which score tells you this?


Table: Distance Metrics Output

|                          |         x|
|:-------------------------|---------:|
|symmetric.difference      | 16.000000|
|branch.score.difference   |  1.927149|
|path.difference           | 15.716234|
|quadratic.path.difference |  6.094534|

The Symmetric (Robinson Foulds Distance) tells us the **sum** of the number of clades unique to each tree. So there are 8 unique clades in each tree.

### Question 4 (1pt)

##### Which gene (R35 or SLC8A1) shows the correct relationships?

The R35 gene shows the correct relationships.

## Run RAxML with 1000 Bootstrap Replicates



## Question 5 (3 pts)

##### Copy and paste your 3 different RAxML commands (with all parameters), to show how you got 3 independent runs with the same data.


![script 1](slurmPic.jpg)


![script 2](slurmPic1.jpg)

## Question 6 (1 pt)

##### Did your results appear to converge? Why or why not?

run_1 likelihood: -4082.569435

run_2 likelihood: -4082.569434

run_3 likelihood: -4082.569434

The highly similar optimized likelihoods suggest convergence.

## Plot your RAxML Tree with Bootstrap Values

### Question 7 (5 pts) 

##### Provide the plot of your best RAxML tree with bootstrap support values. Make sure it is clearly labeled, and all labels are legible.


![](h3_maxLikelihood_ancestralSR_files/figure-docx/bootstrapping-1.png)

## Perform Ancestral State Reconstruction


### Question 8 (5 pts)

##### Provide the plot of your RAxML tree with reconstructed states based on the "all equal" transition cost model. Make sure tip labels are clearly visible, and provide a legend to indicate which state is represented by each color in your plot.  According to these results, how many times do you think limblessness evolved? Does this suggest multiple evolutions of "4-leggedness"?

![](h3_maxLikelihood_ancestralSR_files/figure-docx/ASR_plot-1.png)

It would appear that limblessness has evolved twice since the common ancestor (once by the Tuatara, and once between the worm lizard and the clade above it). The tree plot is showing no uncertainty (no splits in the pie charts)

### Question 9 (5 pts)

This time, set up a custom transition matrix, with the transition from 0 legs to 4 legs being 10 times more costly than the transition from 4 legs to 0. Then re-run the ancestral state reconstruction using maximum parsimony with the custom matrix, and plot your new results with pie charts for the reconstructed values.



##### Provide the plot of your final tree with ancestral state reconstructions using the custom matrix. Again, make sure you have clear labels and a legend.  How does this tree compare to the previous tree in terms of the number of transitions to limblessness or 4-leggedness? 

![](h3_maxLikelihood_ancestralSR_files/figure-docx/ASR_custom_plot-1.png)
This tree actually makes less sense to me as it would seem there are multiple transitions from leglessness, to having legs, and back again even within 2 completely legless clades. I suppose it's possible but it violates our notion of parsimonious tree building, I think.

## References

^1^ Reptiles_SLC8A1.fasta

^2^ Reptiles_R35.fasta
