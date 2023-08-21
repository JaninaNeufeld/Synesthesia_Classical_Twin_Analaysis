#==========
# Program: univariate twin analysis
#
# Author: Mark Taylor
# Date: 2013-11-01
# modified by Janina Neufeld 2023
#==========

# set working directory:
#C:\Users\janneu\OneDrive - Karolinska Institutet\Skrivbordet\Analyses\CATSS_new
setwd('C:/Users/janneu/OneDrive - Karolinska Institutet/Skrivbordet/Analyses/CATSS_new')
list.files()

# load packages and functions:

require(OpenMx); require(psych)

## put miFunctions into working directory
#source('miFunctions.R')

# OR load from link
source('https://openmx.ssri.psu.edu/sites/default/files/myFunctions.R')

#=====
# Prepare data
#=====

### import and check the data:
## v6 = randomized twin order also for OS twins, but nothing changes - I don't get why
data <- read.csv(file='masked_data_for_sharing.csv', header=T, sep=';')
fix(data)


Vars <- 'synesthesia_score' 

# 1=MZ, 2=DZ
MZ<-subset(data, zygosity==1)
DZ<-subset(data, zygosity==2)


nv <- 1 # number of phenotypes
ntv <- nv*2 # total number of variables (nv times 2 twins)
(selVars <- paste(Vars, c(rep(1, nv), rep(2, nv)), sep=''))

### create subsets of data for mz and dz twins:
mz <- subset(data, zygosity==1, selVars) #only the variables of interest, otherwise problem to pcik certain colums
dz <- subset(data, zygosity==2, selVars)


dim(mz)
dim(dz)



#=====
# Assumptions testing 
#=====

#  DEFINING THE MATRIX STRUCTURES

# means: value = starting value which should be close but not exact actual mean
expMeanMZ <- mxMatrix(type='Full', nrow=1, ncol=ntv, free=T, values=0, labels=labFull('m_mz', 1, ntv), 
                      name='ExpMeanMZ')
expMeanDZ <- mxMatrix(type='Full', nrow=1, ncol=ntv, free=T, values=0, labels=labFull('m_dz', 1, ntv),
                      name='ExpMeanDZ')


# variance and covariance:

expCovMZ <- mxMatrix(type='Symm', nrow=ntv, ncol=ntv, free=T, values=c(1, .5, 1), labels=labSymm('mz_cov', ntv),
                     name='ExpCovMZ')
expCovDZ <- mxMatrix(type='Symm', nrow=ntv, ncol=ntv, free=T, values=c(1, .5, 1), labels=labSymm('dz_cov', ntv),
                     name='ExpCovDZ')



# data:
# WHAT FOR THIS STEP: where to find the data for the model, raw rather than matrix

dataMZ <- mxData(mz, type='raw')
dataDZ <- mxData(dz, type='raw')



# objectives:
# with capitals it is the matrix within the object which is called similar just small letters
objMZ <- mxExpectationNormal(covariance='ExpCovMZ', means='ExpMeanMZ', dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ', means='ExpMeanDZ', dimnames=selVars)
funcML <- mxFitFunctionML() #maximun likelyhood which is most commonly used (eg handles missing values well)

# data groups:

modelMZ <- mxModel('MZ', expMeanMZ, expCovMZ, dataMZ, objMZ, funcML)
modelDZ <- mxModel('DZ', expMeanDZ, expCovDZ, dataDZ, objDZ, funcML)

# combine the data groups:

multi <- mxFitFunctionMultigroup(c('MZ', 'DZ'))

# confidence intervals (how to build them)

ci <- mxCI(c('MZ.ExpMeanMZ', 'DZ.ExpMeanDZ', 'MZ.ExpCovMZ', 'DZ.ExpCovDZ'))

# combine all model objects:

satModel <- mxModel('Sat', modelMZ, modelDZ, multi, ci)

# fit the model:
# and add the CIs if wanted (much faster without if many variables)
SatModel <- mxRun(satModel, intervals=T)
(sumSat <- summary(SatModel))

# !!! use describe to check the means and variances for MZ and DZ subcohort to see if it matches with the model
describe(mz$synesthesia_score1)
describe(mz$synesthesia_score2)
describe(dz$synesthesia_score1)
describe(dz$synesthesia_score2)




### assumptions tests

# equal means over twin order:
# we rename all MZ twins the same and all DZ the same so we collapse both twins for means and variances
# then we check the model fit of this model against the saturated


satModel2 <- mxModel(SatModel, name='Sat2')
satModel2 <- omxSetParameters(satModel2, labels=labFull('m_mz', 1, ntv), free=T, values=0, newlabels='m_mz')
satModel2 <- omxSetParameters(satModel2, labels=labFull('m_dz', 1, ntv), free=T, values=0, newlabels='m_dz')
SatModel2 <- mxRun(satModel2, intervals=T)
(sumSat2 <- summary(SatModel2))

# equal means over zygosity:

satModel3 <- mxModel(satModel2, name='Sat3')
satModel3 <- omxSetParameters(satModel3, labels=c('m_mz', 'm_dz'), free=T, values=0, newlabel='m')
SatModel3 <- mxRun(satModel3, intervals=T)
(sumSat3 <- summary(SatModel3))

# equal variances over twin order:

satModel4 <- mxModel(satModel3, name='Sat4')
satModel4 <- omxSetParameters(satModel4, labels=labDiag('mz_cov', ntv), free=T, values=1, newlabel='var_mz')
satModel4 <- omxSetParameters(satModel4, labels=labDiag('dz_cov', ntv), free=T, values=1, newlabel='var_dz')
SatModel4 <- mxRun(satModel4, intervals=T)
# SatModel4 <- mxTryHard(satModel4, intervals=T) # if the starting value fails might just need more attempts
(sumSat4 <- summary(SatModel4))



# equal variances over zygosity:

satModel5 <- mxModel(SatModel4, name='Sat5')
satModel5 <- omxSetParameters(satModel5, labels=c('var_mz', 'var_dz'), free=T, values=1, newlabel='var')
SatModel5 <- mxRun(satModel5, intervals=T)
# SatModel5 <- mxTryHard(satModel5, intervals=T)
(sumSat5 <- summary(SatModel5))



# compare the models:
# checking if the models with assumptions are sign worse
# if twins unqual across twin order needs to be discussed in limiations
# if MZ DZ differ, sibling contrast effect needs to be modelled

mxCompare(SatModel, c(SatModel2, SatModel3, SatModel4, SatModel5))
# unequal means and variances in MZ and DZ can be due to contrast effect




#=====
#=====
# Estimate twin correlations
#=====

### fits a constrained saturated model # similar to model five but more constraints such as SD, intra-class correl

# starting values:

svMe <- mean(mz[, 1], na.rm=T)
svCorMZ <- vech(cov(mz[, 1], mz[, 2], use='complete')) #automatically removes data with missing values
svCorDZ <- vech(cov(dz[, 1], dz[, 2], use='complete'))

# mean:

expMean <- mxMatrix(type='Full', nrow=1, ncol=ntv, free=T, values=svMe, label='m', name='ExpMean')

# standard deviation:

sd <- mxMatrix(type='Diag', nrow=nv, ncol=nv, free=T, values=1.5, labels='sd', name='SD') #diagonal matrix
pad <- mxMatrix(type='Zero', nrow=nv, ncol=nv, name='Pad') # off-diagonal set to zero to remove correlations
expSD <- mxAlgebra(rbind(cbind(SD, Pad),
                         cbind(Pad, SD)), name='ExpSD')

# correlations:
# symmetrical matrix
# this time off-diagonals left in and diagonal replaced with 1
corMZ <- mxMatrix(type='Symm', nrow=nv, ncol=nv, free=T, values=svCorMZ, label='mz_cor', name='CorMZ')
corDZ <- mxMatrix(type='Symm', nrow=nv, ncol=nv, free=T, values=svCorDZ, label='dz_cor', name='CorDZ')

diag <- mxMatrix(type='Unit', nrow=nv, ncol=nv, name='Diag')

expCorMZ <- mxAlgebra(rbind(cbind(Diag, CorMZ),
                            cbind(CorMZ, Diag)), name='ExpCorMZ')
expCorDZ <- mxAlgebra(rbind(cbind(Diag, CorDZ),
                            cbind(CorDZ, Diag)), name='ExpCorDZ')

# calculate covariance:
# to check whether correlations are reasonable

expCovMZ <- mxAlgebra(ExpSD%*%ExpCorMZ%*%t(ExpSD), name='ExpCovMZ')
expCovDZ <- mxAlgebra(ExpSD%*%ExpCorDZ%*%t(ExpSD), name='ExpCovDZ')

# observed data:

dataMZ <- mxData(mz, type='raw')
dataDZ <- mxData(dz, type='raw')


# model objectives:

objMZ <- mxExpectationNormal(covariance='ExpCovMZ', means='ExpMean', dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ', means='ExpMean', dimnames=selVars)


# data groups:
# parameters that assumed to be the same across zygosity
pars <-  list(expMean, sd, pad, expSD)

modelMZ <- mxModel('MZ', pars, corMZ, diag, expCorMZ, expCovMZ, dataMZ, objMZ, funcML)
modelDZ <- mxModel('DZ', pars, corDZ, diag, expCorDZ, expCovDZ, dataDZ, objDZ, funcML)

# combine the data groups_

multi <- mxFitFunctionMultigroup(c('MZ', 'DZ'))


# confidence intervals:
ci <- mxCI(c('MZ.ExpMean', 'MZ.ExpSD', 'MZ.ExpCorMZ', 'DZ.ExpCorDZ'))

# combine all model objects:

satModelCor <- mxModel('SatCor', modelMZ, modelDZ, multi, ci)


# fit the model:

SatModelCor <- mxTryHard(satModelCor, intervals=T)
(sumCor <- summary(SatModelCor))



#=====
# ACE model WITHOUT OS TWINS
#=====

# path coefficients:

pathA <- mxMatrix(type='Full', nrow=nv, ncol=nv, free=T, values=1, label='a_1_1', name='a')
pathC <- mxMatrix(type='Full', nrow=nv, ncol=nv, free=T, values=1, label='c_1_1', name='c')
pathE <- mxMatrix(type='Full', nrow=nv, ncol=nv, free=T, values=1, label='e_1_1', name='e')

# variance components:

varA <- mxAlgebra(a%*%t(a), name='A')
varC <- mxAlgebra(c%*%t(c), name='C')
varE <- mxAlgebra(e%*%t(e), name='E')

# total variance:

varP <- mxAlgebra(A+C+E, name='V')

# means:

expMean <- mxMatrix(type='Full', nrow=1, ncol=ntv, free=T, values=0, label='m', name='ExpMean')

# variance/covariance:
#expected covariance in A+C for MZ and .5A+C for DZ

expCovMZ <- mxAlgebra(rbind(cbind(V, A+C),
                            cbind(A+C, V)), name='ExpCovMZ')
expCovDZ <- mxAlgebra(rbind(cbind(V, 0.5%x%A+C),
                            cbind(0.5%x%A+C, V)), name='ExpCovDZ')

# variance components as proportions of variance:

estVC <- mxAlgebra(cbind(A/V, C/V, E/V), name='EstVC')

# data:

dataMZ <- mxData(mz, type='raw')
dataDZ <- mxData(dz, type='raw')

# objectives:

objMZ <- mxExpectationNormal(covariance='ExpCovMZ', means='ExpMean', dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ', means='ExpMean', dimnames=selVars)
funcML <- mxFitFunctionML()

# data groups:

pars <- list(pathA, pathC, pathE, varA, varC, varE, varP, expMean)

modelMZ <- mxModel('MZ', pars, expCovMZ, estVC, dataMZ, objMZ, funcML)
modelDZ <- mxModel('DZ', pars, expCovDZ, estVC, dataDZ, objDZ, funcML)

# combine groups:

multi <- mxFitFunctionMultigroup(c('MZ', 'DZ'))

# confidence intervals:

ci <- mxCI(c('MZ.EstVC'))

# combine all model objects:

modelACE <- mxModel('ACE', modelMZ, modelDZ, multi, ci)

# fit the model:

ModelACE <- mxRun(modelACE, intervals=T)
(sumACE <- summary(ModelACE))

# compare fit of the model to the saturated model:
mxCompare(SatModel, ModelACE)

## it is significantly worse
#=====
# Nested models
#=====

# AE model:

modelAE <- mxModel(ModelACE, name='AE')
modelAE <- omxSetParameters(modelAE, labels='c_1_1', free=F, values=0)
ModelAE <- mxRun(modelAE, intervals=T)
(sumAE <- summary(ModelAE))
mxCompare(SatModel, ModelAE)

# CE model:

modelCE <- mxModel(ModelACE, name='CE')
modelCE <- omxSetParameters(modelCE, labels='a_1_1', free=F, values=0)
ModelCE <- mxRun(modelCE, intervals=T)
(sumCE <- summary(ModelCE))

# E model:

modelE <- mxModel(ModelAE, name='E')
modelE <- omxSetParameters(modelE, labels='a_1_1', free=F, values=0)
ModelE <- mxRun(modelE)
(sumE <- summary(ModelE))

# compare the fits of the models:
mxCompare(ModelACE, c(ModelAE, ModelCE, ModelE))

