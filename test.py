import pandas as pd
import numpy as np
#import rpy2.robjects as robjects
import pdb
from pandas import DataFrame
import rpy2
from rpy2 import robjects
import gc
import pdb
import os
import math

#print(os.getcwd())

#os.getcwd()
#'C:\\Users\\lenovo'
#os.chdir("C:\\Users\\lenovo\\Desktop\\UKin\\ukin")

import ukin



def simData1(n, m, rho, f, maxNRel=None):
    '''
    Generate Genotype Data for correlated individuals in a population
    Input:
        n -- Number of individuals
        m -- Number of SNPs
        rho -- Average pairwise correlation between individuals
        f -- minor allele frequency for each SNP (If f is given, m will not be effective; If f is not given, f will be generated based on U(0.05, 0.5))
    Output:
        X -- Genotype matrix
        f -- minor allele frequency for each SNP
    # Dependency:
        # binData in R package
    '''
   
    try:
        assert isinstance(f, (list, np.ndarray))
    except AssertionError:
        print('genoSim.simData: f should be a list/Numpy.ndarray!')           
    try:
        assert all(isinstance(item, (int,float,np.int32,np.int64,np.float)) for item in f)
    except AssertionError:
        print('genoSim.simData: The elements in f are not numeric!')
            
    f = np.asarray(f)
    m = f.size  
    if maxNRel==None: maxNRel=n
    
    if n<=maxNRel:
        try:
            assert (0<= rho) and (rho<=1)
        except AssertionError:
            print('genoSim.simData: rho out of range!')
            return
        rho_rel = rho
        n_rel = n
    else:
        try:
            assert (0<=rho) and (rho<=(maxNRel*(maxNRel-1))/(n*(n-1)))
        except AssertionError:
            print('genoSim.simData: rho out of range!')
            return
        rho_rel = rho*n*(n-1)/(maxNRel*(maxNRel-1))
        n_rel = maxNRel

    f_rvector = robjects.FloatVector(f)
    tmp = robjects.r.source('test.R')
    simGenoData = robjects.globalenv['simGenoData']
    X_rel_rvector = simGenoData(n_rel, f_rvector, rho_rel)
    X_rel = np.array(X_rel_rvector, dtype=np.int32)
    
    if n<=maxNRel:
        X=X_rel
    else:
        X_unrel_rvector = simGenoData(n-n_rel, f_rvector, 0)
        X = np.vstack((X_rel, X_unrel_rvector))
        
    #Clear memory
    robjects.r('rm(list=ls(all=TRUE))')
    robjects.r('gc()')
    # robjects.r('print(memory.size())')
    del X_rel_rvector
    if n>maxNRel: del X_unrel_rvector
    gc.collect()

    return X






def simulation1(n1,n2,n3,n4,n5,m,rho1,rho2,rho3,rho4,f):    
    
    X=simData1(2,m,rho=rho1,f=f)

    for i in range((n1-1)):
        a=simData1(2,m,rho=rho1,f=f)
        X=np.vstack((X,a))
    for i in range(n2):
        a=simData1(2,m,rho=rho2,f=f)
        X=np.vstack((X,a))
    for i in range(n3):
        a=simData1(2,m,rho=rho3,f=f)
        X=np.vstack((X,a))
    for i in range(n4):
        a=simData1(2,m,rho=rho4,f=f)
        X=np.vstack((X,a))    
    a=simData1(n5,m,0,f)
    X=np.vstack((X,a))
    return X





def SimData(m,n1,n2,n3,n4,n5,rho1,rho2,rho3,rho4,pi1=0.05,f0=None,f1=None):
    n=2*(n1+n2+n3+n4)+n5
    if f0!=None:
        try:
            assert isinstance(f0, (list, np.ndarray))
        except AssertionError:
            print('gwasSim.simData: f0 should be a list/Numpy.ndarray!')
            return
        try:
            assert all(isinstance(item, (int,float,np.int32,np.int64,np.float)) for item in f0)
        except AssertionError:
            print('gwasSim.simData: The elements in f0 are not numeric!')
            return
        f0 = np.asarray(f0)
        m = f0.size
    else:
        f0 = np.random.uniform(0.1, 0.5, m)
    
    if f1!=None:
        try:
            assert isinstance(f1, (list, np.ndarray))
        except AssertionError:
            print('gwasSim.simData: f1 should be a list/Numpy.ndarray!')
            return
        try:
            assert all(isinstance(item, (int,float,np.int32,np.int64,np.float)) for item in f1)
        except AssertionError:
            print('gwasSim.simData: The elements in f1 are not numeric!')
            return
        try:
            assert len(f0)==len(f1)
        except AssertionError:
            print('gwasSim.simData: The lengths of f0 and f1 are not equal!')
            return
        f1 = np.asarray(f1)
        idx = np.where(f1!=f0)[0]
        m1 = length(idx)
    else:
        m1 = round(m*pi1)
        idx = list(range(0, m1))
        f1 = f0 + np.concatenate((np.random.normal(0.0, 0.1, size=m1), np.zeros(m-m1))) #np.concatenate((np.choice([-1,1],m1)*0.1, np.zeros(m-m1)))
        f1[f1>0.95]=0.95
        f1[f1<0.05]=0.05
           
    X0=simulation1(n1,n2,n3,n4,n5,m,rho1,rho2,rho3,rho4,f0)
    X1=simulation1(n1,n2,n3,n4,n5,m,rho1,rho2,rho3,rho4,f1)
    X=np.vstack((X0,X1))
    Y=np.repeat([0,1], [n,n], axis=0)

    return X,Y,f0,f1,idx



X=SimData(1000,25,25,25,0,50,0.125,0.25,0.5,1,pi1=0.05,f0=None,f1=None)  

A1=ukin.ukinEst1(X[0])
A2=ukin.ukinEst(X[0])
A3=ukin.pairCorr_raw1(X[0])
A4=ukin.pairCorr_raw(X[0])