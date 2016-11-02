#### Copyright 2016 Andrew D Fox
#### 
#### Licensed under the Apache License, Version 2.0 (the "License");
#### you may not use this file except in compliance with the License.
#### You may obtain a copy of the License at
#### 
#### http://www.apache.org/licenses/LICENSE-2.0
#### 
#### Unless required by applicable law or agreed to in writing, software
#### distributed under the License is distributed on an "AS IS" BASIS,
#### WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#### See the License for the specific language governing permissions and
#### limitations under the License.
#### 

# Code adapted from:
# SAMPLE CODE FOR HOUSEMAN, ACCOMANDO, ET AL. (November 6, 2011)
#                  PubMed ID: 22568884
# Code adapted by: Andrew D Fox

library(nlme)
source("wbcInference.R")

load("metaDataMSmethyl.RData")    ## metadata
cellPropsFACS = read.csv("CellPropsFACS.txt", sep="\t")
BetaVal = read.table( "methylation_beta_values_001.txt", header = T, sep = '\t', quote='' )

NUM_COLUMNS = 9
NUM_CTRLS = 8
NUM_CASES = 7	
NUM_CELLTYPES = 6
TARG_RANGE_CTRL = 1:6
TARG_RANGE_CASE = 7:19
TARG_RANGE_CTRL_P1 = 1:8
TARG_RANGE_CASE_P1 = 9:15
NUM_CTRLS_TARGDATA = length( TARG_RANGE_CTRL_P1 )
NUM_CASES_TARGDATA = length( TARG_RANGE_CASE_P1 )
NUM_VALS_STORAGE = 5

m2b <- function(m){ return( 2^m/(1+2^m)  ) }
b2m <- function(b){ return( log2(b/(1-b)) ) }

BetaVal = as.matrix( BetaVal )
colnames(BetaVal) <- substr(colnames(BetaVal), start = 2 , stop = length(colnames(BetaVal)) )
Mv = b2m(BetaVal)

# Cell Type indexes:
i_n   =  c(1,10,19,28,37,45,54,63)
i_4   =  c(2,11,20,29,38,46,55,64)
i_8   =  c(3,12,21,30,39,47,56,65)
i_k   =  c(4,13,22,31,40,48,57,66)
i_b   =  c(5,14,23,32,41,49,58,67)
i_m   =  c(6,15,24,33,42,50,59,68)
i_wbc =  c(7,16,25,34,43,51,60,69)
i_wb  =  c(8,17,26,35,44,52,61,70, 72,73,74,75,76,77)
i_ms  =  c(9,18,27,36,   53,62,71, 78,79,80,81,82,83)
ind = list( i_n, i_4, i_8, i_k, i_b, i_m, i_wbc, i_wb, i_ms)

dataIndex <- function( ctrlSampleNum, ct_num){ return(  (ctrlSampleNum-1)*9 + ct_num ); }
dataIndexVal <- function( ctrlSampleNum, ct_num, vals){ return(  vals[(ctrlSampleNum-1)*9 + ct_num] ); }
cpg_data = as.matrix( read.csv("houseman_refSites_n1826.txt", sep="\t", header=F) )
cpgs = cpg_data[,1]
cpg_cts = cpg_data[,2]
cpg_dirs = cpg_data[,3]

##############################################
##### ##### Step 1: Fit Validation Model (S0)
##############################################

trainData = Mv[ cpgs , -c( i_wb, i_ms ) ]   ## i_wb are wholeblood(ctrl), 9==failedSample, i_ms are wholeblood(case)
targData = Mv[ cpgs , c(i_wb[1:NUM_CTRLS], i_ms[1:NUM_CASES]) ]	# Original controls (8), then original cases (7)
M = length(cpgs)
NTOP_CPGS = length(cpgs)
# Define the validation model:
theModel = y ~ Neut + CD4T + CD8T + NK + Bcell + Mono
sizeModel = 7
validationData_Assay = m2b(trainData)
validationData_Pheno = trainData_pheno
targetData_Assay = m2b(targData)
targetData_Covariates = targData_covariates
# Linear transformation of coefficient vector
# representing contrast to test F statistic
L.forFstat = diag(sizeModel)[-1,]  #All non-intercept coefficients
# Initialize various containers
sigmaResid = sigmaIcept = nObserved = nClusters = Fstat = rep(NA, M)
coefEsts = matrix(NA, M, sizeModel)
coefVcovs =list()

for(j in 1:M){ # For each CpG
  #Remove missing methylation values
  ii = !is.na(validationData_Assay[j,])
  nObserved[j] = sum(ii)
  validationData_Pheno$y = validationData_Assay[j,]
  if(j%%100==0) cat(j,"\n") # Report progress
  try({ # Try to fit a mixed model to adjust for plate
    fit = try(lme(theModel, random=~1|PLATE, data=validationData_Pheno[ii,]))
    if(inherits(fit,"try-error")){ # If LME can't be fit, just use OLS
      fit = lm(theModel, data=validationData_Pheno[ii,])
      fitCoef = fit$coef
      sigmaResid[j] = summary(fit)$sigma
      sigmaIcept[j] = 0
      nClusters[j] = 0
    }
    else{ 
      fitCoef = fit$coef$fixed
      sigmaResid[j] = fit$sigma
      sigmaIcept[j] = sqrt(getVarCov(fit)[1])
      nClusters[j] = length(fit$coef$random[[1]])
    }
    coefEsts[j,] = fitCoef
    coefVcovs[[j]] = vcov(fit)
    useCoef = L.forFstat %*% fitCoef
    useV = L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
    Fstat[j] = (t(useCoef) %*% solve(useV, useCoef))/sizeModel
  })
}
# Name the rows so that they can be easily matched to the target data set
rownames(coefEsts) = rownames(validationData_Assay)
colnames(coefEsts) = names(fitCoef)
# Get P values corresponding to F statistics
Pval = pf(Fstat, sizeModel, nObserved - nClusters - sizeModel + 1, lower.tail=FALSE)
Fstat_mtx = as.matrix(Fstat)
Pval_mtx = as.matrix(Pval)
rownames(Fstat_mtx) <- rownames(validationData_Assay)
rownames(Pval_mtx) <- rownames(validationData_Assay)

###############################################
######## Step 2: Fit Target Model (S1) ########
###############################################
# Contrast matrix:
Lwbc = diag(7)[-1,]
Lwbc[,1]=1
#colnames( coefEsts )
rownames(Lwbc) = colnames(coefEsts)[-1]
colnames(Lwbc) = colnames(coefEsts)

# Denominator degrees-of-freedom for parametric bootstrap
degFree = nObserved - nClusters - (sizeModel-1)
CpGSelection = rownames(coefEsts)[1:NTOP_CPGS]
#Note:  if the CpGs were scattered throughout the array,
#    you would want to select them by name as is indicated here.
#    For this sample version, it would be easier just to use "[1:NTOP_CPGS]"

targetEst = inferWBCbyLme(
  targetData_Assay[1:NTOP_CPGS,],     # Target methylation (cpGs x subjects)
  targetData_Covariates,        # Target phenotype frame (subjects x covariates)
  y~case,                		######+gender+ageCtr,  # Target model (fixed effects)
  ~1|BeadChip,           		# Target adjustment (random effects) [*Footnote 3*]
  coefEsts[CpGSelection,],      # Raw coefficient estimates for WBC 
  Lwbc                   		# Contrast matrix [*Footnote 2*]
)

##############################################
##### ##### Step 3: View projections
##############################################
# Contrast matrix

Lwbc = diag(NUM_CELLTYPES+1)[-1,]
Lwbc[,1]=1
rownames(Lwbc) = colnames(coefEsts)[-1]
colnames(Lwbc) = colnames(coefEsts)
#Lwbc # View contrast matrix
CpGSelection = rownames(coefEsts)[1:NTOP_CPGS]

####### Projections for target data
constrainedCoefs = projectWBC(
  targetData_Assay[1:NTOP_CPGS,],
  coefEsts[CpGSelection,],    
  Lwbc)


cellPropTarget = cellPropsFACS
cellPropsFACS = as.matrix(cellPropsFACS)
rownames(cellPropsFACS) = rownames(constrainedCoefs)[1:NUM_CTRLS]
colnames(cellPropsFACS)<- c("Neut_Gold", "CD4T_Gold", "CD8T_Gold", "NK_Gold", "Bcell_Gold", "Mono_Gold")
ctrlProps = constrainedCoefs[TARG_RANGE_CTRL_P1,]
caseProps = constrainedCoefs[TARG_RANGE_CASE_P1,]

#View the 8 control sample cell proportion predictions:
ctrlProps

cor( c(as.matrix(props_expt)), c(as.matrix(ctrlProps)) )^2
sqrt(mean((props_expt-ctrlProps)^2))

