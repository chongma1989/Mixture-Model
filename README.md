# Reproducible Screening For Large Scale Testing Problem
Screen differentially expressed genes and gene-gene linkage associations 

## Introduction 
In large scale testing problems, p-values are usually obtained by using the whole data, which are in turn used to conduct the significance screening for all of the hypotheses parametrically or nonparametrically. False discovery rate has become widely popular in massive hypothesis tests, because it controls the expected fraction of false positives among all detected significant cases. We combine the subsampling and false discovery rate to reproduce the discovery of significance hypotheses for gaining more power. For each subsampled data, the obtained p-values are fitted to a mixture of baseline and signal distributions by using a boosted EM algorithm, which is used to calculate the false discovery rate (FDR) for each p-value. Repeat this subsampling process many times, a sequence of FDR's for each hypothesis can produce reproducibility for truly significant hypotheses. We conduct a simulation study and real data analysis to show the power of our approach. 

## Steps For The Reproducible Screening 

1. Split the whole samples into $y==x^2$

## Boosted EM algorithm to fit the mixture of Uniform and Beta Distributions

