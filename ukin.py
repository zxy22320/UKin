import pandas as pd
import numpy as np
import rpy2
import rpy2.robjects as robjects
import pdb
from pandas import DataFrame
import gc
import pdb
import os
import math





def pairCorr_raw(X):
    '''
    Estimating pairwise correlations among individuals (scGRM Method, unfiltered)
    Input:
        X -- a matrix (n*m dimension) of numbers of reference alleles, either 2d list or 2d nparray 
             n is the number of individuals, m is the number of SNPs
    Output:
        R -- Estimated Relationship matrix
        sigma2 -- Sample Variance
    '''
    try:
        assert isinstance(X, (list, np.ndarray))
    except AssertionError:
        print('accUtils.accEst: X should be a list/Numpy.ndarray!')
        return
    try:
        assert all(isinstance(item, (int,float,np.int32,np.int64,np.float)) for item in X.flatten())
    except AssertionError:
        print('accUtils.accEst: The elements in X are not numeric!')
        return
    X = np.asarray(X)
    try: 
        assert X.ndim == 2
    except AssertionError:
        print('accUtils.accEst: X is not a matrix (2d array)')
        return
    
    # Reference allele frequencies
    fHat = np.nanmean(X, axis=0)/2.

    # Variances of genotypic values
    sigma2 = 2*fHat*(1-fHat)
    # Sample Variances
    s2 = np.nanvar(X, axis=0, ddof=1)  #Compute the variance along the specified axis, while ignoring NaNs.
    n = np.sum(~np.isnan(X), axis=0)

    sigma2Hat = 2*(sigma2-(1./2-1./(2*n))*s2)

    fHatMat = np.outer(np.ones(X.shape[0]), fHat)
    
    sigma2HatMat = np.outer(np.ones(X.shape[0]), sigma2)
    Xnorm = (X-2*fHatMat)/np.sqrt(sigma2HatMat)
    mask = (~np.isnan(Xnorm)).astype(int) 
    Xnorm[mask==0] = 0
    N = np.dot(mask, mask.T)
    R = np.dot(Xnorm, Xnorm.T)/N

    return R, sigma2, N





def pairCorr_raw1(X,nLevel='auto', prefix=None):
    '''
    Estimating pairwise correlations among individuals (scGRM Method, filtered)
    Input:
        X -- a matrix (n*m dimension) of numbers of reference alleles, either 2d list or 2d nparray 
             n is the number of individuals, m is the number of SNPs
        nLevel -- The threshold value, cane be specified or automatically determined
    Output:
        R -- Estimated Relationship matrix
        sigma2 -- Sample Variance
        nLevel -- The threshold value we actually use
    '''
    try:
        assert isinstance(X, (list, np.ndarray))
    except AssertionError:
        print('accUtils.accEst: X should be a list/Numpy.ndarray!')
        return
    try:
        assert all(isinstance(item, (int,float,np.int32,np.int64,np.float)) for item in X.flatten())
    except AssertionError:
        print('accUtils.accEst: The elements in X are not numeric!')
        return
    X = np.asarray(X)
    try: 
        assert X.ndim == 2
    except AssertionError:
        print('accUtils.accEst: X is not a matrix (2d array)')
        return
    
    l = X.shape[0]
    # Reference allele frequencies
    fHat = np.nanmean(X, axis=0)/2.

    # Variances of genotypic values
    sigma2 = 2*fHat*(1-fHat)
    # Sample Variances
    s2 = np.nanvar(X, axis=0, ddof=1)  #Compute the variance along the specified axis, while ignoring NaNs.
    n = np.sum(~np.isnan(X), axis=0)

    sigma2Hat = 2*(sigma2-(1./2-1./(2*n))*s2)

    fHatMat = np.outer(np.ones(X.shape[0]), fHat)
    
    sigma2HatMat = np.outer(np.ones(X.shape[0]), sigma2)
    Xnorm = (X-2*fHatMat)/np.sqrt(sigma2HatMat)
    mask = (~np.isnan(Xnorm)).astype(int) 
    Xnorm[mask==0] = 0
    N = np.dot(mask, mask.T)
    R = np.dot(Xnorm, Xnorm.T)/N
    
    if nLevel == 'auto':
        Rvec = R[np.triu_indices(l, 1)]
        hist, binEdges = np.histogram(Rvec, bins=100, range=(0, Rvec.max()))
        hist1 = np.append(hist[1:], hist[-1])
        if len(np.where((hist1-hist)>=0)[0])!=0:                                                                                                                                                                                                                             
            nLevel = binEdges[np.where((hist1-hist)>=0)[0][0]]
        else:
            nLevel = None
#设置过滤，低于某个水平设为0
#    nLevel = pow(2,-3.5)
        
    def isNum(value):
        try:
            value + 1
        except TypeError:
            return False
        else:
            return True
    
    if isNum(nLevel):
        diag = np.diag(R)
        R[R<=nLevel] = 0.
        np.fill_diagonal(R, diag)
        if prefix!= None:
            saveGRM2MTX(R, prefix)
 
    return R, sigma2, N, nLevel




def pairCorr(X):
    '''
    Estimating pairwise correlations among individuals (Small bias one, unfiltered)
    Input:
        X -- a matrix (n*m dimension) of numbers of reference alleles, either 2d list or 2d nparray 
             n is the number of individuals, m is the number of SNPs
    Output:
        R -- Estimated Relationship matrix
        sigma2Hat -- Improved sample variance
    '''
    try:
        assert isinstance(X, (list, np.ndarray))
    except AssertionError:
        print('accUtils.accEst: X should be a list/Numpy.ndarray!')
        return
    try:
        assert all(isinstance(item, (int,float,np.int32,np.int64,np.float)) for item in X.flatten())
    except AssertionError:
        print('accUtils.accEst: The elements in X are not numeric!')
        return
    X = np.asarray(X)
    try: 
        assert X.ndim == 2
    except AssertionError:
        print('accUtils.accEst: X is not a matrix (2d array)')
        return
    
    # Reference allele frequencies
    fHat = np.nanmean(X, axis=0)/2.

    # Variances of genotypic values
    sigma2 = 2*fHat*(1-fHat)
    # Sample Variances
    s2 = np.nanvar(X, axis=0, ddof=1)
    n = np.sum(~np.isnan(X), axis=0)

    sigma2Hat = 2*(sigma2-(1./2-1./(2*n))*s2)

    fHatMat = np.outer(np.ones(X.shape[0]), fHat)
    
    sigmaHatMat = np.outer(np.ones(X.shape[0]), np.sqrt(sigma2Hat))
    # Xnorm = (X-2*fHatMat)/np.sqrt(sigma2HatMat)
    # mask = (~np.isnan(Xnorm)).astype(int) 
    # Xnorm[mask==0] = 0
    # R = np.dot(Xnorm, Xnorm.T)/np.dot(mask, mask.T)

    Xcenter = (X-2*fHatMat)
    mask = (~np.isnan(Xcenter)).astype(int)
    Xcenter[mask==0] = 0
    sigmaHatMat[mask==0] = 0
    R = np.dot(Xcenter, Xcenter.T)/np.dot(sigmaHatMat, sigmaHatMat.T)
    N = np.dot(mask, mask.T)
#    R = np.dot(Xcenter, Xcenter.T)/(N*np.nanmean(sigma2Hat))

#    avgCorrRaw = np.mean(R[np.triu_indices(R.shape[0], 1)])
    return R, sigma2Hat, N





def pairCorr1(X, nLevel='auto', prefix=None):
    '''
    Estimating pairwise correlations among individuals (Small bias one, filtered)
    Input:
        X -- a matrix (n*m dimension) of numbers of reference alleles, either 2d list or 2d nparray 
             n is the number of individuals, m is the number of SNPs
        nLevel -- The threshold value, cane be specified or automatically determined
    Output:
        R -- Estimated relationship matrix
        sigma2Hat -- Improved sample variance
        nLevel -- The threshold value we actually use
    '''
    try:
        assert isinstance(X, (list, np.ndarray))
    except AssertionError:
        print('accUtils.accEst: X should be a list/Numpy.ndarray!')
        return
    try:
        assert all(isinstance(item, (int,float,np.int32,np.int64,np.float)) for item in X.flatten())
    except AssertionError:
        print('accUtils.accEst: The elements in X are not numeric!')
        return
    X = np.asarray(X)
    try: 
        assert X.ndim == 2
    except AssertionError:
        print('accUtils.accEst: X is not a matrix (2d array)')
        return
    
    l = X.shape[0]
    
    # Reference allele frequencies
    fHat = np.nanmean(X, axis=0)/2.

    # Variances of genotypic values
    sigma2 = 2*fHat*(1-fHat)
    # Sample Variances
    s2 = np.nanvar(X, axis=0, ddof=1)
    n = np.sum(~np.isnan(X), axis=0)

    sigma2Hat = 2*(sigma2-(1./2-1./(2*n))*s2)

    fHatMat = np.outer(np.ones(X.shape[0]), fHat)
    
    sigmaHatMat = np.outer(np.ones(X.shape[0]), np.sqrt(sigma2Hat))
    # Xnorm = (X-2*fHatMat)/np.sqrt(sigma2HatMat)
    # mask = (~np.isnan(Xnorm)).astype(int) 
    # Xnorm[mask==0] = 0
    # R = np.dot(Xnorm, Xnorm.T)/np.dot(mask, mask.T)

    Xcenter = (X-2*fHatMat)
    mask = (~np.isnan(Xcenter)).astype(int)
    Xcenter[mask==0] = 0
    sigmaHatMat[mask==0] = 0
    R = np.dot(Xcenter, Xcenter.T)/np.dot(sigmaHatMat, sigmaHatMat.T)
    N = np.dot(mask, mask.T)
#    R = np.dot(Xcenter, Xcenter.T)/(N*np.nanmean(sigma2Hat))

#    avgCorrRaw = np.mean(R[np.triu_indices(R.shape[0], 1)])
    
    if nLevel == 'auto':
        Rvec = R[np.triu_indices(l, 1)]
        hist, binEdges = np.histogram(Rvec, bins=100, range=(0, Rvec.max()))
        hist1 = np.append(hist[1:], hist[-1])
        if len(np.where((hist1-hist)>=0)[0])!=0:                                                                                                                                                                                                                             
            nLevel = binEdges[np.where((hist1-hist)>=0)[0][0]]
        else:
            nLevel = None
#    nLevel = pow(2,-3.5)
        
    def isNum(value):
        try:
            value + 1
        except TypeError:
            return False
        else:
            return True
    
    if isNum(nLevel):
        diag = np.diag(R)
        R[R<=nLevel] = 0.
        np.fill_diagonal(R, diag)
        if prefix!= None:
            saveGRM2MTX(R, prefix)
            
    return R, sigma2Hat, N, nLevel







def ukinEst1(X, nLevel='auto', prefix=None, alg=0):
    '''
    Estimating pairwise correlations among individuals (UKin Method, filtered)
    Input:
        X -- a matrix (n*m dimension) of numbers of reference alleles, either 2d list or 2d nparray 
             n is the number of individuals, m is the number of SNPs
        nLevel -- The threshold value, cane be specified or automatically determined
    Output:
        R -- Estimated Relationship matrix
        sigma2Hat -- Sample Variance
        nLevel -- The threshold value we actually use
    '''
    n = X.shape[0]
    if alg==0:
        R, sigma2Hat, N = pairCorr(X)
    else:
        R, sigma2Hat, N = pairCorr_raw(X)

    sumRel = np.sum(R, axis=1)-np.diagonal(R)
    sumRelMat1 = np.outer(np.ones(n), sumRel)
    sumRelMat2 = sumRelMat1.T
    Radjusted = R+sumRelMat1/2.+sumRelMat2/2.+1
    # saveGRM2MTX(Radjusted, 'ukin_nofilter')
    Radjusted_new = np.copy(Radjusted)
                                                                                                                                                                                                                                                                                                                                                          
    if nLevel == 'auto':
        Rvec = Radjusted[np.triu_indices(n, 1)]
        hist, binEdges = np.histogram(Rvec, bins=100, range=(0, Rvec.max()))
        hist1 = np.append(hist[1:], hist[-1])
        if len(np.where((hist1-hist)>=0)[0])!=0:                                                                                                                                                                                                                             
            nLevel = binEdges[np.where((hist1-hist)>=0)[0][0]]
        else:
            nLevel = None
    #nLevel = pow(2,-3.5)
        
    def isNum(value):
        try:
            value + 1
        except TypeError:
            return False
        else:
            return True
    
    if isNum(nLevel):
        diag = np.diag(Radjusted_new)
        Radjusted_new[Radjusted_new<=nLevel] = 0.
        np.fill_diagonal(Radjusted_new, diag)
        if prefix!= None:
            saveGRM2MTX(Radjusted_new, prefix)
    
    return Radjusted_new, sigma2Hat, N, nLevel



def ukinEst(X, alg=0):
    '''
    Estimating pairwise correlations among individuals (UKin method, unfiltered)
    Input:
        X -- a matrix (n*m dimension), either 2d list or 2d nparray 
             n is the number of individuals, m is the number of SNPs
    Output:
        R -- Estimated Relationship matrix
        sigma2Hat -- Sample Variance
        N -- 
    '''
    n = X.shape[0]
    if alg==0:
        R, sigma2Hat, N = pairCorr(X)
    else:
        R, sigma2Hat, N = pairCorr_raw(X)

    sumRel = np.sum(R, axis=1)-np.diagonal(R)
    sumRelMat1 = np.outer(np.ones(n), sumRel)
    sumRelMat2 = sumRelMat1.T
    Radjusted = R+sumRelMat1/2.+sumRelMat2/2.+1
    # saveGRM2MTX(Radjusted, 'ukin_nofilter')
 #   Radjusted_new = np.copy(Radjusted)
   
    return Radjusted, sigma2Hat, N



