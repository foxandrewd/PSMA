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

m2b <- function(m){ return( 2^m/(1+2^m)  ) }	# Turn methylation M values to beta values
b2m <- function(b){ return( log2(b/(1-b)) ) }	# Turn methylation beta values to M values

props_expt = read.csv("CellPropsFACS.txt", sep="\t", head=F)
BetaVal = read.table( "methylation_beta_values_001.txt", header = T, sep = '\t', quote='' )
BetaVal = as.matrix( BetaVal )
colnames(BetaVal) <- substr(colnames(BetaVal), start = 2 , stop = length(colnames(BetaVal)) )
Mv = b2m(BetaVal)

N = dim(BetaVal)[1]
m = dim(BetaVal)[2]
linear_loocv_n = 8

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

# The Cell-type - specific Reference markers that identify the different cell types
# These are strongly demethylated in the cell type for which they are specific (M < ~-2);
# and they are strongly methylated (M > ~2) in the other cell types (exception is CD4 & CD8 T cells which
# could not be separated from each other given the data we had. However, these T cell subsets can be
# well-separated from each of the 4 other non-T-cell types, and the markers shown do just that:

CTS_NEUTRO = c('cg13618969', 'cg01699630', 'cg06270401', 'cg25600606', 'cg11070172', 'cg03146219' )
CTS_CD4T = c('cg07545925', 'cg16452866', 'cg05160234', 'cg07015803', 'cg15880738', 'cg10111816' )
CTS_CD8T = c('cg24841244', 'cg07728874', 'cg05160234', 'cg13750061', 'cg24612198', 'cg10111816')
CTS_CD3T = c('cg24841244', 'cg07728874', 'cg05160234', 'cg13750061', 'cg24612198', 'cg10111816')
CTS_NK = c('cg23015664', 'cg19915997', 'cg27274718', 'cg25386954', 'cg26275360')
CTS_BCELL = c('cg11661493', 'cg19260718', 'cg02087075', 'cg21743182', 'cg03402235', 'cg22907103' )
CTS_MONO = c('cg23244761','cg10480329','cg02244028','cg18066690','cg04468741','cg04045544')
CTS_IN_USE = NA
CTSes_ALL = list( CTS_NEUTRO, CTS_CD4T, CTS_CD8T, CTS_NK, CTS_BCELL, CTS_MONO)

res_pred = c(); res_pred2 = c();
colnames(props_expt) = c("NT", "T4", "T8", "NK", "BC", "MO")
rownames(props_expt) = colnames(Mv)[i_wb][1:dim(props_expt)[1]]

THE_CTS = 0
for( CTSes_i in  CTSes_ALL  ) {
	THE_CTS = THE_CTS + 1								# Cell type index (1,2,3,4,5,6) <=> (N,T4,T8,NK,B,M)
	CTS_IN_USE = CTSes_i								# Cell type currently under examination	
	L_WhBlood = Mv[,i_wb]								# WB data (Control)
	S_WhBlood = Mv[,i_ms]								# WB data (Case)
	L_CTSavg =  colMeans( Mv[ CTS_IN_USE , i_wb  ])		# CTS reference signal (mean of CTS ref. markers) for controls
	S_CTSavg =  colMeans( Mv[ CTS_IN_USE , i_ms  ])		# CTS reference signal (mean of CTS ref. markers) for cases
	data = Mv[CTSes_i , i_wb[c(1:linear_loocv_n)] ]
	data = colMeans( data )
	LOO_cts = c() ;	LOO_prop = c() ; LOO_pred = c()
	all_intercepts = c() ; all_slopes = c()
	all_full_int = c() ; all_full_slope = c() ;	all_full_pred = c() ;
	for( LOO in c(1:8) ){
		the_rest_data = data[ -LOO ]
		the_rest_prop = props_expt[ -LOO, THE_CTS]
		the_all_data = data[]
		the_all_prop = props_expt[ , THE_CTS]
		the_LOO_cts_data = data[ LOO ]					# CTS 'REF' methyl for LOO sample
		the_LOO_prop_expt = props_expt[LOO, THE_CTS]	# CTS true proportion of LOO sample
		the_LOO_prop_actual_true_value = props_expt[LOO , THE_CTS]
		#### Build the LOO linear model
		curr_LM = lm( the_rest_prop ~ the_rest_data)
		curr_LM_intercept = coef(curr_LM)[1]
		curr_LM_slope = coef(curr_LM)[2]
		the_LOO_prop_prediction = curr_LM_intercept + curr_LM_slope*the_LOO_cts_data
		LOO_cts = c(LOO_cts, the_LOO_cts_data)
		LOO_prop = c(LOO_prop, the_LOO_prop_expt)
		LOO_pred = c( LOO_pred, the_LOO_prop_prediction)
		all_intercepts = c( all_intercepts, curr_LM_intercept )
		all_slopes = c( all_slopes , curr_LM_slope )
		LM_fulldata = lm( the_all_prop ~ the_all_data)
		curr_fulldata_LM_intercept = coef(LM_fulldata)[1]
		curr_fulldata_LM_slope = coef(LM_fulldata)[2]
		curr_fulldata_pred = curr_fulldata_LM_intercept + curr_fulldata_LM_slope*the_LOO_cts_data
		all_full_int = c(all_full_int, curr_fulldata_LM_intercept)
		all_full_slope = c(all_full_slope, curr_fulldata_LM_slope)
		all_full_pred = c(all_full_pred, curr_fulldata_pred)
	}
	
	results = data.frame( LOO_cts=LOO_cts, LOO_prop=LOO_prop, LOO_pred=LOO_pred, LOO_slope=all_slopes, 
						  LOO_intercept=all_intercepts, FULL_pred=all_full_pred, FULL_slope=all_full_slope,
						  FULL_intercept=all_full_int)
	
	expected_intercept = mean(results$LOO_cts + results$LOO_pred)
	results$LOO_pred2 = expected_intercept - results$LOO_cts
	res_pred = cbind( res_pred,  results$LOO_pred )
	res_pred2 = cbind( res_pred2,  results$LOO_pred2 )

}
colnames(res_pred) = c("NT", "T4", "T8", "NK", "BC", "MO")
rownames(res_pred) = colnames(Mv)[i_wb][1:dim(props_expt)[1]]

res_pred
cor( c(as.matrix(props_expt)), c(as.matrix(res_pred)) )^2
sqrt(mean((props_expt-res_pred)^2))
