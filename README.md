# Reproducible Screening For Large Scale Testing Problem
Screen differentially expressed genes and gene-gene linkage associations 

## Introduction 
In large scale testing problems, p-values are usually obtained by using the whole data, which are in turn used to conduct the significance screening for all of the hypotheses parametrically or nonparametrically. False discovery rate has become widely popular in massive hypothesis tests, because it controls the expected fraction of false positives among all detected significant cases. We combine the subsampling and false discovery rate to reproduce the discovery of significance hypotheses for gaining more power. For each subsampled data, the obtained p-values are fitted to a mixture of baseline and signal distributions by using a boosted EM algorithm, which is used to calculate the false discovery rate (FDR) for each p-value. Repeat this subsampling process many times, a sequence of FDR's for each hypothesis can produce reproducibility for truly significant hypotheses. We conduct a simulation study and real data analysis to show the power of our approach. 

## Steps For The Reproducible Screening 

1. Randomly split the whole samples into two data sets with almost equal sample size. 
2. Calculate the marginal p-values of each hypothesis for each of the two data sets, respectively. Combine the two bunch of p-values together, and fit them to a mixture of uniform and beta distributions. 
3. Using the fitted mixture model, for each hypothesis, calculate the false discover rate (FDR) under the two subsamples, respectively. Take the minimum one as its adjusted false discovery rate. 
4. Repeat steps 1-3 100 times, calculate the integral of the empirical cdf of FDR for each hypothesis, which is essentially its expected true discovery rate. This statistics can become the stability selection. 

## Boosted EM algorithm to fit the mixture of Uniform and Beta Distributions

Algorithm:
1. Update the weights 
2. Update the scale parameter alpha in the beta distribution
3. Update the scale parameter beta in the beta distribution

Repeat steps 1-3 until converges 

<!DOCTYPE html>
<html>
<body>
<img src="https://latex.codecogs.com/gif.latex?O_t=\text { Onset event at time bin } t" width="100", height="20">
</body>
</html>

<dl>
  <dt>Definition list</dt>
  <dd>Is something people use sometimes.</dd>

  <dt>Markdown in HTML</dt>
  <dd>Does *not* work **very** well. Use HTML <em>tags</em>.</dd>
</dl>
