#==========
# Program: bivariate twin model, including OS twins as DZ
# Author: Mark Taylor
# modified by Janina Neufeld, 2023

# Note: for saturated used cor instead of cov for twin correlations 
# reason: since synethesia measure so skewed, the covariance differs too much from the correlation and model can't handle that
  

setwd('C:/Users/janneu/OneDrive - Karolinska Institutet/Skrivbordet/Analyses/CATSS_new')
list.files()
  
#==========


require(OpenMx); require(xlsx)


## put miFunctions into working directory
#source('miFunctions.R')

# OR load from link
source('https://openmx.ssri.psu.edu/sites/default/files/myFunctions.R')


#===== 
# Prepare data
#=====

### import and check data:
# chose one of the options below, depending on the bivariate model to run
data <- read.csv(file='Supplementary_data_masked_for_sharing_overall_traits.csv', header=T, sep=';')
data <- read.csv(file='Supplementary_data_masked_for_sharing_RRBID.csv', header=T, sep=';')
data <- read.csv(file='Supplementary_data_masked_for_sharing_SIC.csv', header=T, sep=';')
fix(data)






### select variables for analysis:


Vars <- c('synesthesia_score', 'autistic_traits')
Vars <- c('synesthesia_score', 'RRBID_traits')
Vars <- c('synesthesia_score', 'SIC_traits')


nv <- 2
ntv <- nv*2
(selVars <- paste(Vars, c(rep(1, nv), rep(2, nv)), sep=''))



### subsets of data for analysis:
## here OS integrated in DZ
mz <- subset(data, zygosity==1, selVars)
dz <- subset(data, zygosity==2, selVars)




#=====
# Estimate twin correlations
#=====

### fit a constrained saturated model to estimate the twin correlations:

# starting values:
# should be close to but not too close to correct value
# adding a little bit of random noise might solve it

svMe <- colMeans(mz[, 1:2], na.rm=T)
#svMe <- svMe + .5

svPh <- vechs(cor(mz[, 1:2], use='complete'))

# use cor instead of cov because otherwise too different -- due to the skew of the synesthesia score the model fails with cov
svCorMZ <- vech(cor(mz[, 1:2], mz[, 3:4], use='complete'))
svCorDZ <- vech(cor(dz[, 1:2], dz[, 3:4], use='complete'))


# means:

expMean <- mxMatrix(type='Full', nrow=1, ncol=ntv, free=T, values=svMe, labels=labFull('m', 1, nv), 
                     name='ExpMean')

# standard deviations:

sd <- mxMatrix(type='Diag', nrow=nv, ncol=nv, free=T, values=.9999, labels=labDiag('sd', nv), name='SD')

pad <- mxMatrix(type='Zero', nrow=nv, ncol=nv, name='Pad')

expSD <- mxAlgebra(rbind(cbind(SD, Pad),
                         cbind(Pad, SD)), name='ExpSD')

# phenotypic correlation:

rph <- mxMatrix(type='Stand', nrow=nv, ncol=nv, free=T, values=svPh, label='rph', name='Rph')

# twin correlations:

corMZ <- mxMatrix(type='Symm', nrow=nv, ncol=nv, free=T, values=svCorMZ, labels=labSymm('mz_cor', nv),
                   name='CorMZ')
corDZ <- mxMatrix(type='Symm', nrow=nv, ncol=nv, free=T, values=svCorDZ, labels=labSymm('dz_cor', nv),
                   name='CorDZ')

expCorMZ <- mxAlgebra(rbind(cbind(Rph, CorMZ),
                            cbind(CorMZ, Rph)), name='ExpCorMZ')
expCorDZ <- mxAlgebra(rbind(cbind(Rph, CorDZ),
                            cbind(CorDZ, Rph)), name='ExpCorDZ')

# calculate covariances:

expCovMZ <- mxAlgebra(ExpSD%*%ExpCorMZ%*%t(ExpSD), name='ExpCovMZ')
expCovDZ <- mxAlgebra(ExpSD%*%ExpCorDZ%*%t(ExpSD), name='ExpCovDZ')

# data:

dataMZ <- mxData(mz, type='raw')
dataDZ <- mxData(dz, type='raw')

# objectives:

objMZ <- mxExpectationNormal(covariance='ExpCovMZ', means='ExpMean', dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ', means='ExpMean', dimnames=selVars)
funcML <- mxFitFunctionML()

# data groups:

pars <- list(expMean, sd, expSD, rph)

modelMZ <- mxModel('MZ', pars, pad, corMZ, expCorMZ, expCovMZ, dataMZ, objMZ, funcML)
modelDZ <- mxModel('DZ', pars, pad, corDZ, expCorDZ, expCovDZ, dataDZ, objDZ, funcML)

# combine the groups:

multi <- mxFitFunctionMultigroup(c('MZ', 'DZ'))

# confidence intervals:

ci <- mxCI(c('MZ.ExpCorMZ', 'DZ.ExpCorDZ'))

# combine all model objectives:

satModelCor <- mxModel('SatCor', modelMZ, modelDZ, multi, ci)

# fit the model:
SatModelCor <- mxTryHard(satModelCor, intervals=T) 
#In confidence intervals: look if lower bound avove 0 and note estimate
(sumCor <- summary(SatModelCor))

  


#=====
# Fully saturated model # this is mostly just to double check that nothing odd is going on
#=====

# starting values:

svCov <- c(1, rep(.5, 3), 1, rep(.5, 2), 1, .5, 1)

# means:

expMeanMZ <- mxMatrix(type='Full', nrow=1, ncol=ntv, free=T, values=.001, labels=labFull('m_mz', 1, ntv),
                       name='ExpMeanMZ')
expMeanDZ <- mxMatrix(type='Full', nrow=1, ncol=ntv, free=T, values=.001, labels=labFull('m_dz', 1, ntv),
                       name='ExpMeanDZ')

# covariances:

expCovMZ <- mxMatrix(type='Symm', nrow=ntv, ncol=ntv, free=T, values=svCov, labels=labSymm('mz_cov', ntv),
                      name='ExpCovMZ')
expCovDZ <- mxMatrix(type='Symm', nrow=ntv, ncol=ntv, free=T, values=svCov, labels=labSymm('dz_cov', ntv),
                      name='ExpCovDZ')

# data:

dataMZ <- mxData(mz, type='raw')
dataDZ <- mxData(dz, type='raw')

# objectives:

objMZ <- mxExpectationNormal(covariance='ExpCovMZ', means='ExpMeanMZ', dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ', means='ExpMeanDZ', dimnames=selVars)
funcML <- mxFitFunctionML()

# data groups:

modelMZ <- mxModel('MZ', expMeanMZ, expCovMZ, dataMZ, objMZ, funcML)
modelDZ <- mxModel('DZ', expMeanDZ, expCovDZ, dataDZ, objDZ, funcML)

# combine groups:

multi <- mxFitFunctionMultigroup(c('MZ', 'DZ'))

# combine all model objects:

satModel <- mxModel('Sat', modelMZ, modelDZ, multi)

# fit the model:

SatModel <- mxTryHard(satModel)
(sumSat <- summary(SatModel))





#=====
# ACE model
#=====

# path coefficients:

pathA <- mxMatrix(type='Lower', nrow=nv, ncol=nv, free=T, values=c(1, .5, 1), labels=labLower('a', nv), name='a')
pathC <- mxMatrix(type='Lower', nrow=nv, ncol=nv, free=T, values=c(1, .5, 1), labels=labLower('c', nv), name='c')
pathE <- mxMatrix(type='Lower', nrow=nv, ncol=nv, free=T, values=c(1, .5, 1), labels=labLower('e', nv), name='e')

# variance components:

varA <- mxAlgebra(a%*%t(a), name='A')
varC <- mxAlgebra(c%*%t(c), name='C')
varE <- mxAlgebra(e%*%t(e), name='E')

# total variance and standard deviations:

varP <- mxAlgebra(A+C+E, name='V')

matI <- mxMatrix(type='Iden', nrow=nv, ncol=nv, name='I')
isd <- mxAlgebra(solve(sqrt(I*V)), name='iSD')

# etiological correlations:

rph <- mxAlgebra(solve(sqrt(I*V))%&%V, name='rPh')
ra <- mxAlgebra(solve(sqrt(I*A))%&%A, name='rA')
rc <- mxAlgebra(solve(sqrt(I*C))%&%C, name='rC')
re <- mxAlgebra(solve(sqrt(I*E))%&%E, name='rE')

# variance components as proportions:

estVC <- mxAlgebra(cbind(A/V, C/V, E/V), name='EstVC')

# means:

expMean <- mxMatrix(type='Full', nrow=1, ncol=ntv, free=T, values=.001, labels=labFull('m', 1, nv), 
                    name='ExpMean')

# variance/covariance:

expCovMZ <- mxAlgebra(rbind(cbind(V, A+C),
                            cbind(A+C, V)), name='ExpCovMZ')
expCovDZ <- mxAlgebra(rbind(cbind(V, 0.5%x%A+C),
                            cbind(0.5%x%A+C, V)), name='ExpCovDZ')

# data:

dataMZ <- mxData(mz, type='raw')
dataDZ <- mxData(dz, type='raw')

# objectives:

objMZ <- mxExpectationNormal(covariance='ExpCovMZ', means='ExpMean', dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ', means='ExpMean', dimnames=selVars)
funcML <- mxFitFunctionML()

# data groups:

pars <- list(pathA, pathC, pathE, varA, varC, varE, varP, matI, isd, rph, ra, rc, re, estVC, expMean)

modelMZ <- mxModel('MZ', pars, expCovMZ, dataMZ, objMZ, funcML)
modelDZ <- mxModel('DZ', pars, expCovDZ, dataDZ, objDZ, funcML)

# combine groups:

multi <- mxFitFunctionMultigroup(c('MZ', 'DZ'))

# confidence intervals:

ci <- mxCI(c('MZ.EstVC', 'MZ.rPh', 'MZ.rA', 'MZ.rC', 'MZ.rE'))

# combine all model objects:

modelACE <- mxModel('ACE', modelMZ, modelDZ, multi, ci)

# fit the model:

ModelACE <- mxTryHard(modelACE, intervals=T)
# output confidence intervals: EstVC = estimated variance explained by A: 1,1= syn; 2,2 = OCD, 1,2/2,1 = syn*OCD
# C = 1,3; 2,4; 2,3/ 1,4; E = 1,5; 2,6; 1,6/2,5
(sumACE <- summary(ModelACE, verbose=F))

# compare model fits:

mxCompare(SatModel, ModelACE)


#=====
# Nested models
#=====

### drop individual components:

# AE model:

modelAE <- mxModel(ModelACE, name='AE')
modelAE <- omxSetParameters(modelAE, labels=labLower('c', nv), free=F, values=0)
ModelAE <- mxTryHard(modelAE, intervals=T)
(sumAE <- summary(ModelAE))


# CE model:

modelCE <- mxModel(ModelACE, name='CE')
modelCE <- omxSetParameters(modelCE, labels=labLower('a', nv), free=F, values=0)
ModelCE <- mxTryHard(modelCE, intervals=T)
(sumCE <- summary(ModelCE))

# E model:

modelE <- mxModel(ModelAE, name='E')
modelE <- omxSetParameters(modelE, labels=labLower('a', nv), free=F, values=0)
ModelE <- mxTryHard(modelE, intervals=T)
(sumE <- summary(ModelE))

# model fits:

mxCompare(ModelACE, c(ModelAE, ModelCE, ModelE))

# No E
modelNoE <- mxModel(ModelAE, name='noE')
modelNoE <- omxSetParameters(modelNoE, labels='e_2_1', free=F, values=0)
ModelNoE <- mxRun(modelNoE, intervals=T)
summary(ModelNoE)
coefficients(ModelNoE)

ModelNoE$MZ.EstVC

mxCompare(ModelAE, ModelNoE)



### variance components AE:
# colums 1&2 = A, 3&4 = C, 5&6 = E
varcomp <- rbind(
  mxEval(MZ.EstVC,ModelAE),
  mxSE(MZ.EstVC,ModelAE)
)

varcomp <- round(varcomp, 2)

#CIs AE
CIs<- rbind(mxEval(MZ.EstVC,ModelAE)-1.96*mxSE(MZ.EstVC,ModelAE),
            mxEval(MZ.EstVC,ModelAE)+1.96*mxSE(MZ.EstVC,ModelAE))
CIs <- round(CIs, 2)

varcomp <- rbind(varcomp, CIs)
rownames(CIs) <- c('lower', 'lower', 'upper', 'upper')
colnames(CIs) <- c('A', 'A', 'C', 'C', 'E', 'E')

