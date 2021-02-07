# UKin: an Unbiased kinship estimator

---
We develop algorithms of an improved kinship estimation method, UKin, which can reduce both bias 
and root mean square error (RMSE) in the estimation of genomic relationship matrix.
Please refer to our paper "Correcting statistical bias in correlation-based kinship estimators" 
for more details about the UKin method. (https://biorxiv.org/cgi/content/short/2021.01.13.426515v1)

## Content

---
The ukin.py file includes the core algorithms for the UKin method, see the Function section 
for more details. In the test.py file, we provide some algorithms for generating 
genotype data for correlated individuals in a population, which serves as 
a test dataset for the UKin method. The test.R file needs to be sourced in test.py.

Note that to run the test.py file, you also need to install rpy2. 



## Function

---
We develop algorithms of two scGRM methods and UKin method. For each method we provide an original version together 
with an alternative one which filtered the unrelated individual pairs by setting estimated correlation under specified threshold to be 0.

### Original scGRM method:
Usage:
```
pairCorr_raw(X)
```
Estimating pairwise correlations among individuals with classical GRM method (original scGRM method, unfiltered)

Input:
```
X -- a matrix (n*m dimension) of numbers of reference alleles, either 2d list or 2d nparray 
     n is the number of individuals, m is the number of SNPs
```
Output:
```
R -- Estimated relationship matrix
sigma2 -- Sample variance
N -- The number of available SNPs for all individual pairs
```
    
### Filtered original scGRM method:
Usage:
```
pairCorr_raw1(X, nLevel='auto', prefix=None)
```
Estimating pairwise correlations among individuals classical GRM method (original scGRM method, filtered)

Input:
```
X -- a matrix (n*m dimension) of numbers of reference alleles, either 2d list or 2d nparray 
     n is the number of individuals, m is the number of SNPs
nLevel -- The threshold value, can be specified or automatically determined
```
Output:
```
R -- Estimated relationship matrix
sigma2 -- Sample variance
N -- The number of available SNPs for all individual pairs
nLevel -- The threshold value we actually use
```

### Improved scGRM method:
Usage:
```
pairCorr_raw(X)
```
Estimating pairwise correlations among individuals with improved scGRM method where the sample variance is replaced by another unbiased estimation of variance 
(improved scGRM method, unfiltered)

Input:
```
X -- a matrix (n*m dimension) of numbers of reference alleles, either 2d list or 2d nparray 
     n is the number of individuals, m is the number of SNPs
```
Output:
```
R -- Estimated relationship matrix
sigma2Hat -- Improved sample variance
N -- The number of available SNPs for all individual pairs
```
    
### Filtered improved scGRM method:
Usage:
```
pairCorr_raw1(X,nLevel='auto', prefix=None)
```
Estimating pairwise correlations among individuals with improved scGRM method where the sample variance is replaced by another unbiased estimation of variance 
(improved scGRM method, filtered)

Input:
```
X -- a matrix (n*m dimension) of numbers of reference alleles, either 2d list or 2d nparray 
     n is the number of individuals, m is the number of SNPs
nLevel -- The threshold value, can be specified or automatically determined
```
Output:
```
R -- Estimated relationship matrix
sigma2Hat -- Improved sample variance
N -- The number of available SNPs for all individual pairs
nLevel -- The threshold value we actually use
```


### UKin method:
Usage:
```
ukinEst(X, alg=0)
```
Estimating pairwise correlations among individuals with the UKin method (UKin method, unfiltered)

Input:
```
X -- a matrix (n*m dimension), either 2d list or 2d nparray 
     n is the number of individuals, m is the number of SNPs
alg -- an indicator of which scGRM method we use in UKin, if alg=0, we use the improved scGRM method, otherwise the original one
```
Output:
```
Radjusted -- Estimated relationship matrix
sigma2Hat -- Improved sample variance (if alg=0) or sample variance (else)
N -- The number of available SNPs for all individual pairs
```
    
### Filtered UKin method:
Usage:
```
ukinEst1(X, nLevel='auto', prefix=None, alg=0)
```
Estimating pairwise correlations among individuals with the UKin method (UKin Method, filtered)

Input:
```
X -- a matrix (n*m dimension) of numbers of reference alleles, either 2d list or 2d nparray 
     n is the number of individuals, m is the number of SNPs
nLevel -- The threshold value, can be specified or automatically determined
alg -- an indicator of which scGRM method we use in UKin, if alg=0, we use the improved scGRM method, otherwise the original one
```
Output:
```
Radjusted_new -- Estimated Relationship matrix
sigma2Hat -- Improved sample variance (if alg=0) or sample variance (else)
N -- The number of available SNPs for all individual pairs
nLevel -- The threshold value we actually use
```



## Contact Information

---
For help or communication related to UKin, please contact Wei Jiang (w.jiang@yale.edu) or Xiangyu Zhang (zxy22320@mail.ustc.edu.cn).
