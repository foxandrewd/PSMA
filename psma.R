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
error_sites = read.table("error_sites.txt")
BetaVal = as.matrix( BetaVal )
colnames(BetaVal) <- substr(colnames(BetaVal), start = 2 , stop = length(colnames(BetaVal)) )
Mv = b2m(BetaVal)

N = dim(BetaVal)[1]
m = dim(BetaVal)[2]

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

colnames(props_expt) = c("NT", "T4", "T8", "NK", "BC", "MO")
rownames(props_expt) = colnames(Mv)[i_wb][1:dim(props_expt)[1]]
PLIM = 0.05

sink("psma_output.txt", append = F)

THE_CTS = 0
for( CTSes_i in  list( CTS_NEUTRO, CTS_CD8T, CTS_NK, CTS_BCELL, CTS_MONO)  ) {

	THE_CTS = THE_CTS + 1								# Cell type index (1,2,3,4,5,6) <=> (N,T4,T8,NK,B,M)
	CTS_IN_USE = CTSes_i								# Cell type currently under examination	
	L_WhBlood = Mv[,i_wb]								# WB data (Control)
	S_WhBlood = Mv[,i_ms]								# WB data (Case)
	L_CTSavg =  colMeans( Mv[ CTS_IN_USE , i_wb  ])		# CTS reference signal (mean of CTS ref. markers) for controls
	S_CTSavg =  colMeans( Mv[ CTS_IN_USE , i_ms  ])		# CTS reference signal (mean of CTS ref. markers) for cases
	
	for( i in 1:N  ){
		current_cpg_site_name = rownames(Mv)[i]
		if( is.element(current_cpg_site_name, error_sites$V1)){next}				# Ignore if this is an invalid (error) site on the chip
		ctrl_i = L_WhBlood[i,]													# WB data for this CpG site (Controls)
		case_i = S_WhBlood[i,]													# WB data for this CpG site (Cases)
		r2_ctrl_i = cor( ctrl_i , L_CTSavg)^2									# Pearson R^2 for control vs. Celltype Reference signal
		r2_case_i = cor( case_i , S_CTSavg)^2									# Pearson R^2 for case vs. Celltype Reference signal
		L_lm = lm(ctrl_i ~ L_CTSavg)											# Linear model (LM) of this CpG vs. Celltype Reference signal (Controls)
		S_lm = lm(case_i ~ S_CTSavg)											# Linear model (LM) of this CpG vs. Celltype Reference signal (Cases)
		L_intercept = L_lm$coefficients[1]										# Control LM intercept
		L_slope = L_lm$coefficients[2]											# Control LM slope (CTS methylation level for this CpG in the current cell type)
		S_intercept = S_lm$coefficients[1]										# Case LM intercept
		S_slope = S_lm$coefficients[2]											# Case LM slope (CTS methylation level for this CpG in the current cell type)
		L_stderr = summary(L_lm)$coefficients["L_CTSavg","Std. Error"]			# Std error of slope (methylation level) estimate (Controls)
		S_stderr = summary(S_lm)$coefficients["S_CTSavg","Std. Error"]			# Std error of slope (methylation level) estimate (Cases)
		L_slope_range = c(L_slope-1.96*L_stderr, L_slope+1.96*L_stderr)
		S_slope_range = c(S_slope-1.96*S_stderr, S_slope+1.96*S_stderr)
		slopes_diff_i = abs( L_slope / S_slope)									# Ratio of case to control methylation
		slopes_diff_dist_i = abs(S_slope - L_slope)
		L_Fstat = summary(L_lm)$fstatistic
		S_Fstat = summary(S_lm)$fstatistic
		L_Fpval_i = pf( L_Fstat[1], L_Fstat[2], L_Fstat[3], lower.tail=FALSE)	# P-value of goodness-of-fit of the linear model (Controls)
		S_Fpval_i = pf( S_Fstat[1], S_Fstat[2], S_Fstat[3], lower.tail=FALSE)	# P-value of goodness-of-fit of the linear model (Cases)

		dt <- data.frame(L_CTSavg, ctrl_i, c("ctrl") )							# Contruct table of CTS ref signal, CpG site signal, and whether the sample is a case or control
		colnames(dt) = c("CTS_marker_mean", "GOI_value", "Aff")
		dt_case <- data.frame(S_CTSavg, case_i, c("case") )
		colnames(dt_case) = c("CTS_marker_mean", "GOI_value", "Aff")
		dt <- rbind( dt, dt_case)

		lmfit <- summary(lm(GOI_value ~ CTS_marker_mean * Aff, data=dt))		# Fit combined case/control linear model
		aff_by_cts_pVal = lmfit$coefficients[4,4]
		ctsMarkers_GOI_pVal = lmfit$coefficients[2,4]

    	# WHOLE-BLOOD TESTING:
     	WB_pVal = t.test( ctrl_i, case_i)
		WB_Ctrl_Methyl_Mv_mean = mean(ctrl_i)
		WB_Case_Methyl_Mv_mean = mean(case_i)
		WB_Ctrl_Methyl_Mv_SEM = sd(ctrl_i)/sqrt(length(ctrl_i))
		WB_Case_Methyl_Mv_SEM = sd(case_i)/sqrt(length(case_i))
		
		cat( rownames(Mv)[i], i, THE_CTS, L_slope, S_slope, L_stderr, S_stderr, L_Fpval_i, S_Fpval_i,
			slopes_diff_i, slopes_diff_dist_i, aff_by_cts_pVal, r2_ctrl_i, r2_case_i, WB_pVal$p.value, 
			WB_Ctrl_Methyl_Mv_mean, WB_Case_Methyl_Mv_mean, WB_Ctrl_Methyl_Mv_SEM, WB_Case_Methyl_Mv_SEM,
			"\n")

	}
}

sink()
