# Script Name: Phenotypic correlations explorations as prelude to MRDoC
# Adapted by Oginni OA (2022) from Oginni OA and Rijsdijk FV (2021)

# To determine Causal relationships from Sasme-sex attraction (SSA) to Psychological Distress (PD) and vice versa.
# Latent SSA; 	1 Indicator: Sexual orientation - SO (Liability threshold variable)
# Latent PD;	2 Indicators:  dep=Depressive seymptoms ; anx=Anxiety symptoms     
# Mediator:  Victimisation
# Moderators: Childhood Gender Nonconformity (CGN), retrospective Early-life Adversity (ELA)  
# Heterogeneity variable: 	Sex
# Covariate on Thresholds:	Age (mean-centered) and sex
# Zygosity variable: 		zyg (1=MZM, 2=DZM, 3=MZF, 4=DZF, 5=OSDZ)
# ------------------------------------------------------------------------------------

#a. Predictor: SSA
# i. SO

#b. Outcomes: PD
# i. Anx and Dep

#c. Mediators
# i. Victimisation

#d. Moderators
# i. CGN
# ii. ELA

# TABLE OF CONTENTS
#__________________

# Data Preparation

# PART 1: PHENOTYPIC MODELS By Zygosity groups: MZ/DZ 
#____________________________________________________

# Model Ia - Phenotypic correlations between variables including mediator and moderators
# Model Ib - Phenotypic correlations between all PRSs and variables of interest and mediator
# Model IIa - Phenotypic Covariance Model across Latent Constructs  with Vict
# Model III - Phenotypic Mediation Model

# PART 2: BIOMETRIC GENETIC MODELS By Zygosity groups: MZ/DZ
#___________________________________________________________

# Model IV - ACE MODEL for latent factors (SSA, PD, Vict)
# Model Va - Mendelian Randomisation Direction of Causation (MRDoC) MODEL for SO>MH
# Model Vb - Mendelian Randomisation Direction of Causation (MRDoC) MODEL for PD>SO
# Model VIa - Mendelian Randomisation Direction of Causation (MRDoC) MODEL for SO>PD incorporating moderation by ELA
# Model VIb - Mendelian Randomisation Direction of Causation (MRDoC) MODEL for SO>PD incorporating moderation by CGN
# Model VIc - Mendelian Randomisation Direction of Causation (MRDoC) MODEL for PD>SO incorporating moderation by ELA
# Model VId - Mendelian Randomisation Direction of Causation (MRDoC) MODEL for PD>SO incorporating moderation by CGN
# Model VIIa - Mendelian Randomisation Direction of Causation (MRDoC) MODEL for SO>PD with mediation by Vict
# Model VIIb - Mendelian Randomisation Direction of Causation (MRDoC) MODEL for PD>SO with mediation by Vict
# Model VIIIa - Mendelian Randomisation Direction of Causation (MRDoC) MODEL for SO>VICT>PD with mediation by Vict and confounding by CGN
# Model VIIIb - Mendelian Randomisation Direction of Causation (MRDoC) MODEL for PD>VICT>SO with mediation by Vict and confounding by CGN

# PART 3: SEX DIFFERENCES MODELS (Phenotypic and Biometric Genetic) by-Zygosity groups: MZM/DZM MZF/DZF DZO 
#_______________________________________________________________________

# Model IIa - Sex Diff - Phenotypic Covariance Model across Latent Constructs  including Vict with sex differences
# Model IIb - Sex Diff - Phenotypic Covariance Model across Latent Constructs  including ELA, CGN and Vict with sex differences
# Model IV - Sex Diff - ACE Factor MODEL (including victimisation)
# Model Vai - sex diff - Mendelian Randomisation Direction of Causation (MRDoC) MODEL for SO on PD with sex differences
# Model Vaii - Sex diff - Mendelian Randomisation Direction of Causation (MRDoC) MODEL for PD on SO with sex differences

# Note: In all models Age and Sex covariate effects are regressed out of continuous variables before modelling

rm(list=ls())
ls()
Sys.setenv(OMP_NUM_THREADS=parallel::detectCores()) #this has to be before loading the library OpenMx

#check you have the packages:
search()

library(psych)
library(OpenMx)
library(Hmisc)			
library(dplyr)

mxOption(NULL, "Default optimizer", "SLSQP") # SLSQP is a better optimizer for ordinal data

#to add in another file for directory path
#setwd("C:/Users/Kunle/Desktop/Project 6/Analysis/Data and scripts use/")

# ******************************************************************************************************************************************
# (1) Read in data file and check descriptive statistics
# ******************************************************************************************************************************************

# Reads data 

library(foreign)
Data <- read.dta ("3. Long - 524 Olakunle Oginni sexual orientation Aug2021.dta", convert.dates = TRUE, convert.underscore=TRUE, convert.factors=FALSE)
psych::describe(Data)

#exclude those with medical conditions
Data	<- subset(Data, exclude1==0) 
psych::describe(Data)

#fill in missing co-twin's age as proxy
Data <- mutate(Data, age21 = ifelse(is.na(u1cage1),u1cage2 , u1cage1))
Data <- mutate(Data, age21 = ifelse(is.na(u1cage2),u1cage1 , u1cage2))

Data <- mutate(Data, age4 = ifelse(is.na(drepage1),drepage2 , drepage1))
Data <- mutate(Data, age4 = ifelse(is.na(drepage2),drepage1 , drepage2))

#Set male dzo==twin 1 and female dzo==twin 2
Data$twin[Data$sex1==1 & Data$sexzyg==5] <- "1"		#Male
Data$twin[Data$sex1==0 & Data$sexzyg==5] <- "2"		#Female

#To view what I have done
tab	<-table(Data$sexzyg, Data$twin, Data$sex1, dnn=c('Zygosity','Birth order','Sex'))
colnames(tab) <- c('Twin 1','Twin 2')
rownames(tab) <- c('mzm','dzm','mzf','dzf','dzo')
tab

attr(Data, "var.labels")
summary(Data)
psych::describe(Data)
str(Data)
dim(Data)
names(Data)

# ******************************************************************************************************************************************
# (2) Explore variable of interest: recode, regress-out age and sex, check distribution and transform if necessary.
# ******************************************************************************************************************************************

#---------

Data$Rec.SexOri<-mxFactor(Data$Rec.SexOri, levels=c(1:5) )

#---------

Data$Dep_resid <- residuals(lm(Data$Sum.MFQ ~ Data$sex1 + Data$age21, na.action="na.exclude"))
psych::describe(Data$Sum.MFQ)
psych::describe(Data$Dep_resid)
multi.hist(Data$Sum.MFQ)
multi.hist(Data$Dep_resid)

Data$Dep_residlogT <- ((log(Data$Dep_resid + 8))*1.8)+1
psych::describe(Data$Dep_residlogT)
multi.hist(Data$Dep_residlogT)

#---------

Data$Depct_resid <- residuals(lm(Data$Sum.MFQ ~ Data$sex1 + Data$age21 + Data$ACE.Tot, na.action="na.exclude"))
psych::describe(Data$Sum.MFQ)
psych::describe(Data$Depct_resid)
multi.hist(Data$Sum.MFQ)
multi.hist(Data$Depct_resid)

Data$Depct_residlogT <- (((Data$Depct_resid + 10))/3)
psych::describe(Data$Depct_residlogT)
multi.hist(Data$Depct_residlogT)

#---------

Data$Depcg4_resid <- residuals(lm(Data$Sum.MFQ ~ Data$sex1 + Data$age21 + Data$CGN4.stand, na.action="na.exclude"))
psych::describe(Data$Sum.MFQ)
psych::describe(Data$Depcg4_resid)
multi.hist(Data$Sum.MFQ)
multi.hist(Data$Depcg4_resid)

Data$Depcg4_residlogT <- ((log(Data$Depcg4_resid + 8))*3)
psych::describe(Data$Depcg4_residlogT)
multi.hist(Data$Depcg4_residlogT)

#---------

Data$Anx_resid <- residuals(lm(Data$Sum.GenAnx ~ Data$sex1 + Data$u2cage1, na.action="na.exclude"))
psych::describe(Data$Sum.GenAnx)
psych::describe(Data$Anx_resid)
multi.hist(Data$Sum.GenAnx)
multi.hist(Data$Anx_resid)

Data$Anx_residlogT <- (log(Data$Anx_resid + 10)) *2
psych::describe(Data$Anx_residlogT)
multi.hist(Data$Anx_residlogT)

#---------

Data$Anxct_resid <- residuals(lm(Data$Sum.GenAnx ~ Data$sex1 + Data$u2cage1 + Data$ACE.Tot, na.action="na.exclude"))
psych::describe(Data$Sum.GenAnx)
psych::describe(Data$Anxct_resid)
multi.hist(Data$Sum.GenAnx)
multi.hist(Data$Anxct_resid)

#Data$Anxct_residlogT <- ((sqrt(Data$Anxct_resid + 30)) *2)-6
Data$Anxct_residlogT <- ((log(Data$Anxct_resid + 30)) *6)-12
psych::describe(Data$Anxct_residlogT)
multi.hist(Data$Anxct_residlogT)

#---------

Data$Anxcg4_resid <- residuals(lm(Data$Sum.GenAnx ~ Data$sex1 + Data$u2cage1 + Data$CGN4.stand, na.action="na.exclude"))
psych::describe(Data$Sum.GenAnx)
psych::describe(Data$Anxcg4_resid)
multi.hist(Data$Sum.GenAnx)
multi.hist(Data$Anxcg4_resid)

Data$Anxcg4_residlogT <- (log(Data$Anxcg4_resid + 11)) *2
psych::describe(Data$Anxcg4_residlogT)
multi.hist(Data$Anxcg4_residlogT)

#---------

Data$ELA_resid <- residuals(lm(Data$ACE.Tot ~ Data$sex1 + Data$u2cage1, na.action="na.exclude"))
psych::describe(Data$ACE.Tot)
psych::describe(Data$ELA_resid)
multi.hist(Data$ACE.Tot)
multi.hist(Data$ELA_resid)

Data$ELA_residlogT <- log(Data$ELA_resid + 7)*2
multi.hist(Data$ELA_residlogT)
psych::describe(Data$ELA_residlogT)

#---------

Data$VictTot_resid <- residuals(lm(Data$Tot.vict ~ Data$sex1 + Data$age21, na.action="na.exclude"))
psych::describe(Data$Tot.vict)
psych::describe(Data$VictTot_resid)
multi.hist(Data$Tot.vict)
multi.hist(Data$VictTot_resid)

Data$VictTot_residlogT <- ((log(Data$VictTot_resid + 5)*1))
Data$VictTot_residlogT <- ((Data$VictTot_residlogT) /0.8)+1
multi.hist(Data$VictTot_residlogT)
psych::describe(Data$VictTot_residlogT)

#---------

Data$CGN4_resid <- residuals(lm(Data$CGN4.stand ~ Data$sex1 + Data$age4, na.action="na.exclude"))
psych::describe(Data$CGN4.stand)
psych::describe(Data$CGN4_resid)
multi.hist(Data$CGN4.stand)
multi.hist(Data$CGN4_resid)

Data$CGN4_residlogT <- (Data$CGN4_resid + 4)*1
multi.hist(Data$CGN4_residlogT)
psych::describe(Data$CGN4_residlogT)

#---------

Data$PRSSO_resid <- residuals(lm(Data$Same.sex.behaviour.Ganna2019.FR6 ~ Data$sex1 + Data$age4, na.action="na.exclude"))
psych::describe(Data$Same.sex.behaviour.Ganna2019.FR6)
psych::describe(Data$PRSSO_resid)
multi.hist(Data$Same.sex.behaviour.Ganna2019.FR6)
multi.hist(Data$PRSSO_resid)

Data$PRSSO_residlogT <- (Data$PRSSO_resid * -25)+ 5 #I multiple by minus to reverse the effect - Andrea says LDPred may have made an error
multi.hist(Data$PRSSO_residlogT)
psych::describe(Data$PRSSO_residlogT)

#---------

Data$PRSDep_resid <- residuals(lm(Data$Depression.Howard2019.UKBonly.FR ~ Data$sex1 + Data$age4, na.action="na.exclude"))
psych::describe(Data$Depression.Howard2019.UKBonly.FR)
psych::describe(Data$PRSDep_resid)
multi.hist(Data$Depression.Howard2019.UKBonly.FR)
multi.hist(Data$PRSDep_resid)

Data$PRSDep_residlogT <- (Data$PRSDep_resid *8)+5
multi.hist(Data$PRSDep_residlogT)
psych::describe(Data$PRSDep_residlogT)

#---------

Data$PRSAnx_resid <- residuals(lm(Data$Anxiety.Purves2019.FRCT0.3 ~ Data$sex1 + Data$age4, na.action="na.exclude"))
psych::describe(Data$Anxiety.Purves2019.FRCT0.3)
psych::describe(Data$PRSAnx_resid)
multi.hist(Data$Anxiety.Purves2019.FRCT0.3)
multi.hist(Data$PRSAnx_resid)

Data$PRSAnx_residlogT <- (Data$PRSAnx_resid *-10)+5 #I multiple by -1 to reverse the effect - Andrea says LDPred may have made an error
multi.hist(Data$PRSAnx_residlogT)
psych::describe(Data$PRSAnx_residlogT)

#---------

names(Data)
psych::describe(Data)

# ******************************************************************************************************************************************
# (3) Generate Pair-Wise File
# ******************************************************************************************************************************************
# Subset only variables of interest

SubData <-Data[,c("randomfamid","sexzyg","sex1","age21","twin","Rec.SexOri","Dep_residlogT","Depct_residlogT","Depcg4_residlogT","Anx_residlogT","Anxct_residlogT","Anxcg4_residlogT",
			"ELA_residlogT","CGN4_residlogT","VictTot_residlogT","PRSSO_residlogT","PRSDep_residlogT","PRSAnx_residlogT")]

psych::describe(SubData)
summary(SubData)
str(SubData)

#Rename variables 
colnames(SubData) <- c('family','zyg','sex','age','TwNo','SO','Dep','Depct','Depcg4','Anx','Anxct','Anxcg4',
				'ELA','CGN4','Vict','PRSSO','PRSDep','PRSAnx') 
names(SubData)

TWINdata <- reshape(SubData, idvar = c("family"), timevar = "TwNo", direction = "wide")
#Some variables were NA, these are the siblings who are beyond 2nd female and male siblings (n=49)
names(TWINdata) 
psych::describe(TWINdata)

#Pick only twins for the time being
TWINdata2 <-TWINdata[,c('family','zyg.2','sex.2','age.2','SO.2','Dep.2','Depct.2','Depcg4.2','Anx.2','Anxct.2','Anxcg4.2','ELA.2','CGN4.2','Vict.2','PRSSO.2','PRSDep.2','PRSAnx.2',
					   'zyg.1','sex.1','age.1','SO.1','Dep.1','Depct.1','Depcg4.1','Anx.1','Anxct.1','Anxcg4.1','ELA.1','CGN4.1','Vict.1','PRSSO.1','PRSDep.1','PRSAnx.1')]

psych::describe(TWINdata2)

# ******************************************************************************************************************************************
# (4) Rename Variables
# Rename variables using shorter names and to get rid of the . in the names (OpenMx cannot handle . in names)
# ******************************************************************************************************************************************
colnames(TWINdata2) <- c('family','zyg2','sex2','age2','SO2','Dep2','Depct2','Depcg42','Anx2','Anxct2','Anxcg42','ELA2','CGN42','Vict2','PRSSO2','PRSDep2','PRSAnx2',
					    'zyg1','sex1','age1','SO1','Dep1','Depct1','Depcg41','Anx1','Anxct1','Anxcg41','ELA1','CGN41','Vict1','PRSSO1','PRSDep1','PRSAnx1') 

names(TWINdata2)
psych::describe(TWINdata2)

# To obtain complete number of twin pairs per variable
#SO
 selVars	<- c('SO1', 'SO2')
 matrix(
	c(
	table(TWINdata2$zyg2),
	table(TWINdata2$zyg2[complete.cases(TWINdata2[,selVars])]),
	table(TWINdata2$zyg2[!is.na(TWINdata2$SO1) & is.na(TWINdata2$SO2)|!is.na(TWINdata2$SO2) & is.na(TWINdata2$SO1)])
	),
	ncol=5,
	byrow=TRUE,
	dimnames=list(c("total", "complete","incomplete"), c("mzm","dzm","mzf","dzf","dzo"))
	)

#Dep
 selVars	<- c('Dep1', 'Dep2')
 matrix(
	c(
	table(TWINdata2$zyg2),
	table(TWINdata2$zyg2[complete.cases(TWINdata2[,selVars])]),
	table(TWINdata2$zyg2[!is.na(TWINdata2$Dep1) & is.na(TWINdata2$Dep2)|!is.na(TWINdata2$Dep2) & is.na(TWINdata2$Dep1)])
	),
	ncol=5,
	byrow=TRUE,
	dimnames=list(c("total", "complete","incomplete"), c("mzm","dzm","mzf","dzf","dzo"))
	)

#Anx
 selVars	<- c('Anx1', 'Anx2')
 matrix(
	c(
	table(TWINdata2$zyg2),
	table(TWINdata2$zyg2[complete.cases(TWINdata2[,selVars])]),
	table(TWINdata2$zyg2[!is.na(TWINdata2$Anx1) & is.na(TWINdata2$Anx2)|!is.na(TWINdata2$Anx2) & is.na(TWINdata2$Anx1)])
	),
	ncol=5,
	byrow=TRUE,
	dimnames=list(c("total", "complete","incomplete"), c("mzm","dzm","mzf","dzf","dzo"))
	)

#Vict
 selVars	<- c('Vict1', 'Vict2')
 matrix(
	c(
	table(TWINdata2$zyg2),
	table(TWINdata2$zyg2[complete.cases(TWINdata2[,selVars])]),
	table(TWINdata2$zyg2[!is.na(TWINdata2$Vict1) & is.na(TWINdata2$Vict2)|!is.na(TWINdata2$Vict2) & is.na(TWINdata2$Vict1)])
	),
	ncol=5,
	byrow=TRUE,
	dimnames=list(c("total", "complete","incomplete"), c("mzm","dzm","mzf","dzf","dzo"))
	)

#PRS - Note that only one of a pair of mz twins was genotyped
 selVars	<- c('PRSSO1', 'PRSSO2')
	table(TWINdata2$zyg2)
	table(TWINdata2$zyg2[complete.cases(TWINdata2[,selVars])])
	table(TWINdata2$zyg2[!is.na(TWINdata2$PRSSO1) & is.na(TWINdata2$PRSSO2)|!is.na(TWINdata2$PRSSO2) & is.na(TWINdata2$PRSSO1)])
	
#CGN4
 selVars	<- c('CGN41', 'CGN42')
 matrix(
	c(
	table(TWINdata2$zyg2),
	table(TWINdata2$zyg2[complete.cases(TWINdata2[,selVars])]),
	table(TWINdata2$zyg2[!is.na(TWINdata2$CGN41) & is.na(TWINdata2$CGN42)|!is.na(TWINdata2$CGN42) & is.na(TWINdata2$CGN41)])
	),
	ncol=5,
	byrow=TRUE,
	dimnames=list(c("total", "complete","incomplete"), c("mzm","dzm","mzf","dzf","dzo"))
	)

#ELA
 selVars	<- c('ELA1', 'ELA2')
 matrix(
	c(
	table(TWINdata2$zyg2),
	table(TWINdata2$zyg2[complete.cases(TWINdata2[,selVars])]),
	table(TWINdata2$zyg2[!is.na(TWINdata2$ELA1) & is.na(TWINdata2$ELA2)|!is.na(TWINdata2$ELA2) & is.na(TWINdata2$ELA1)])
	),
	ncol=5,
	byrow=TRUE,
	dimnames=list(c("total", "complete","incomplete"), c("mzm","dzm","mzf","dzf","dzo"))
	)

# ******************************************************************************************************************************************
# (5) Prepare data for modelling
# ******************************************************************************************************************************************

#*****************************************************************
# Model 1a - Phenotypic correlations between variables including mediator and moderators
#*****************************************************************

nv			<- 6				# number of variables for a twin = 1 in Univariate
nvo 			<- 1     			#number of ordinal variables per twin
nvc 			<- nv-nvo  			#number of continuous variables per twin
poso 			<- nvo 			#position where ordinal variables start
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nth			<- 4				# number of max thresholds
nlower		<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor			<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
ninc 			<- nth-1 			#number of max increments
ncovariates 	<- 2 				#number of covariates

Groups		<-c("mz", "dz")
Vars			<- c('SO','Dep','Anx','Vict','ELA','CGN4')
selVars		<- c('SO1','Dep1','Anx1','Vict1','ELA1','CGN41',
			     'SO2','Dep2','Anx2','Vict2','ELA2','CGN42')
useVars		<- c('SO1','Dep1','Anx1','Vict1','ELA1','CGN41',
		     	     'SO2','Dep2','Anx2','Vict2','ELA2','CGN42','sex1','age1','sex2','age2')

mzData		<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , useVars)
dzData		<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , useVars)

psych::describe(mzData)
psych::describe(dzData)

mzData$SO1[is.na(mzData$age1)] <- NA
mzData$SO2[is.na(mzData$age2)] <- NA
dzData$SO1[is.na(dzData$age1)] <- NA
dzData$SO2[is.na(dzData$age2)] <- NA

mzData$SO1[is.na(mzData$sex1)] <- NA
mzData$SO2[is.na(mzData$sex2)] <- NA
dzData$SO1[is.na(dzData$sex1)] <- NA
dzData$SO2[is.na(dzData$sex2)] <- NA

mzData$age1[is.na(mzData$age1)] <- 999
mzData$age2[is.na(mzData$age2)] <- 999
dzData$age1[is.na(dzData$age1)] <- 999
dzData$age2[is.na(dzData$age2)] <- 999

mzData$sex1[is.na(mzData$sex1)] <- 999
mzData$sex2[is.na(mzData$sex2)] <- 999
dzData$sex1[is.na(dzData$sex1)] <- 999
dzData$sex2[is.na(dzData$sex2)] <- 999

psych::describe(mzData)
psych::describe(dzData)

mzData_st		<- mzData%>% select(all_of(selVars))%>% mutate_if(is.factor, as.numeric)#I convert factors to numeric data so I can create starting values for correlations
dzData_st		<- dzData%>% select(all_of(selVars))%>% mutate_if(is.factor, as.numeric)

psych::describe(mzData_st)
psych::describe(dzData_st)

cor(mzData_st, use=  "complete")
cor(dzData_st, use = "complete")

# ******************************************************************************************************************************************
# Multivariate phenotypic Constrained MODEL 
# Fits a constrained Polychoric correlation model 
# Ordinal combined with Continuous
# Means and variance constrained to be equal across birth order and zyg groups
# ******************************************************************************************************************************************

# CREATE LABELS & START VALUES as objects(to ease specification)

# Labels For correlations
(rphLabs	<- paste("r",1:ncor,sep="")) #Within person
(MZbLabs <- paste("rmz", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep=""))
(DZbLabs <- paste("rdz", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep=""))

# Labels For means  #NAs because the first variable is ordinal
(LabM	<- paste("m",2:nv,sep=""))
MLabs		<-c(NA,LabM,NA,LabM) 

# Labels for Standard deviations (SDs), #NAs because the first variable is ordinal
(LabSD	<- paste("sd",2:nv,sep=""))
SDLabs	<-c(NA,LabSD,NA,LabSD) #NA because the first variable is categorical

# Starting values for correlations
jiggle	<-rnorm(nlower, mean = 0, sd = .1)

(StWithinperson  	<-vechs(cor(mzData_st[,1:nv],use="complete")))
(StBetweenMZ  	<-vech(cor(mzData_st[,1:nv],mzData_st[,(nv+1):ntv],use="complete")))
(StBetweenDZ  	<-vech(cor(dzData_st[,1:nv],dzData_st[,(nv+1):ntv],use="complete")))

# Starting values for M & SDs
(Stmean	<-colMeans(mzData_st[,2:nv],na.rm=TRUE))
(Stsd 	<-sapply(mzData_st[,2:nv],sd, na.rm=TRUE))

StM 		<-c(0, Stmean, 0, Stmean)
StSD  	<-c(1, Stsd, 1, Stsd)

# Free parameters
(Pat  	<- c( rep(F,nvo), rep(TRUE, nvc)))

# -------------------------------------------------------------------------------------------------------------------
#
# MODEL 1a:	Constrained Correlation Model
#
# -------------------------------------------------------------------------------------------------------------------
# 1) Fits a constrained Polychoric correlation model
# TH same across twins and across zyg groups
# Age effect is different acros variables, but same across thresholds within variables (if c>2)
# There is one overall set of rPH and the x-trait x-twin correlations are symmetric
# ------------------------------------------------------------------------------------------------------------------------------

# CREATE LABELS & START VALUES as objects(to ease specification)
# I allow the threshold and increments not to vary by birth order or zygosity and specify only one set
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 for mz twin 1 (mz)

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')

ThPat		<-c(TRUE,TRUE,TRUE,TRUE)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=T, values=.2, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=T, values=.2, labels=LabCovS, name="BsexTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=T, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=F, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2), name="expThres")

# Matrix for expected Means
# First & 10th means are false, i.e. means for SO  variable are set to 0
Mean		<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(Pat, Pat), values=c(StM), labels=c(MLabs), name="expm" ) #first variable is not free to be estimated

# Matrices to store SD and cor matrices
SD		<-mxMatrix( type="Diag", nrow=ntv, ncol=ntv, free=c(Pat,Pat), values=c(StSD), labels=c(SDLabs), lbound=0, name="sd")

Rph		<-mxMatrix("Stand", nv, nv, free = TRUE, values = StWithinperson, labels=rphLabs, name="within") 
MZb		<-mxMatrix("Symm", nv, nv, free = TRUE, values = StBetweenMZ, labels=MZbLabs, name="BetweenMZ") 
DZb		<-mxMatrix("Symm", nv, nv, free = TRUE, values = StBetweenDZ, labels=DZbLabs, name="BetweenDZ") 
CorMZ		<-mxAlgebra(rbind(cbind(within,BetweenMZ), cbind(BetweenMZ, within)), name="expCorMZ")
CorDZ		<-mxAlgebra(rbind(cbind(within,BetweenDZ), cbind(BetweenDZ, within)), name="expCorDZ")

# Expected covariance matrices
# SD & Correlation matrices from above
covMZ		<-mxAlgebra( expression=sd %&% expCorMZ, name="expCovMZ")
covDZ		<-mxAlgebra( expression=sd %&% expCorDZ, name="expCovDZ")

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objmz  	<- mxExpectationNormal( covariance="expCovMZ", means="expm", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2") )
objdz		<- mxExpectationNormal( covariance="expCovDZ", means="expm", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2") )

fitFunction 	<- mxFitFunctionML()
#fitFunction 	<- mxFitFunctionWLS()

# Combine Groups
modelMZ	<- mxModel(obsage1, obsage2, obssex1, obssex2, betaA, betaS, Thr, inc, Thres, Mean,  SD, Rph, MZb, CorMZ, covMZ, dataMZ, objmz, fitFunction, name="MZ")
modelDZ	<- mxModel(obsage1, obsage2, obssex1, obssex2, betaA, betaS, Thr, inc, Thres, Mean,  SD, Rph, DZb, CorDZ, covDZ, dataDZ, objdz, fitFunction, name="DZ")
ciW		<- mxCI('MZ.within')
ciBmz		<- mxCI('MZ.BetweenMZ')
ciBdz		<- mxCI('DZ.BetweenDZ')
minus2ll	<- mxAlgebra(expression=MZ.objective + DZ.objective, name="m2LL")
obj		<- mxFitFunctionAlgebra("m2LL")
corModel	<- mxModel('cor', modelMZ, modelDZ, minus2ll, obj, ciW, ciBmz, ciBdz) 

# ------------------------------------------------------------------------------
# RUN Multivariate cor Model
corFit    <- mxTryHardOrdinal(corModel, intervals=F)
#corFit    <- mxRun(corModel, intervals=F)
(corSumm  <- summary(corFit))

# Generate Output
mxEval(MZ.within, corFit)
mxEval(MZ.BetweenMZ, corFit)
mxEval(DZ.BetweenDZ, corFit)
mxEval(MZ.expCovMZ, corFit)
mxEval(DZ.expCovDZ, corFit)


#**************************************************************************
# Model 1b - Phenotypic correlations between all PRSs and variables of interest and mediator
#**************************************************************************

nv			<- 7				# number of variables for a twin = 1 in Univariate
nvo 			<- 1     			# number of ordinal variables per twin
nvc 			<- nv-nvo  			# number of continuous variables per twin
poso 			<- nvo 			# position where ordinal variables start
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nth			<- 4				# number of max thresholds
nlower		<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor			<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
ninc 			<- nth-1 			# number of max increments
ncovariates 	<- 2 				# number of covariates

Groups		<- c("mz", "dz")
#Vars			<- c('SO','Dep','Anx','Vict','PRSSO','PRSDep','PRSAnx'')
selVars		<- c('SO1','Dep1','Anx1','Vict1','PRSSO1','PRSDep1','PRSAnx1',
			     'SO2','Dep2','Anx2','Vict2','PRSSO2','PRSDep2','PRSAnx2')
useVars		<- c('SO1','Dep1','Anx1','Vict1','PRSSO1','PRSDep1','PRSAnx1',
			     'SO2','Dep2','Anx2','Vict2','PRSSO2','PRSDep2','PRSAnx2','sex1','age1','sex2','age2')

mzData		<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , useVars)
dzData		<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , useVars)

psych::describe(mzData)
psych::describe(dzData)

mzData$SO1[is.na(mzData$age1)] <- NA
mzData$SO2[is.na(mzData$age2)] <- NA
dzData$SO1[is.na(dzData$age1)] <- NA
dzData$SO2[is.na(dzData$age2)] <- NA

mzData$SO1[is.na(mzData$sex1)] <- NA
mzData$SO2[is.na(mzData$sex2)] <- NA
dzData$SO1[is.na(dzData$sex1)] <- NA
dzData$SO2[is.na(dzData$sex2)] <- NA

mzData$age1[is.na(mzData$age1)] <- 999
mzData$age2[is.na(mzData$age2)] <- 999
dzData$age1[is.na(dzData$age1)] <- 999
dzData$age2[is.na(dzData$age2)] <- 999

mzData$sex1[is.na(mzData$sex1)] <- 999
mzData$sex2[is.na(mzData$sex2)] <- 999
dzData$sex1[is.na(dzData$sex1)] <- 999
dzData$sex2[is.na(dzData$sex2)] <- 999

psych::describe(mzData)
psych::describe(dzData)

mzData_st		<- mzData%>% select(all_of(selVars))%>% mutate_if(is.factor, as.numeric)#I convert factors to numeric data so I can create starting values for correlations
dzData_st		<- dzData%>% select(all_of(selVars))%>% mutate_if(is.factor, as.numeric)

psych::describe(mzData_st)
psych::describe(dzData_st)

cor(mzData_st[,1:nv], use = "complete")
cor(dzData_st, use = "complete")

# ******************************************************************************************************************************************
# Multivariate phenotypic Constrained MODEL 
# Fits a constrained Polychoric correlation model 
# Ordinal combined with Continuous
# Means and variance constrained to be equal across birth order and zyg groups
# ******************************************************************************************************************************************

# CREATE LABELS & START VALUES as objects(to ease specification)

# Labels For correlations
(rphvLabs	<- paste("rv",1:6,sep="")) #Within person - for the 4 measured variables - different in mz
(rphPLabs	<- paste("rP",1:3,sep="")) #Within person - for the last 3 polygenic scores - same in mz constraint because only one mz twin was genotyped per pair
(rphvPLabs	<- paste("rvP",1:12,sep="")) #Within person - for the cross PRS-variable correlations - same in mz as above
(MZbLabs 	<- paste("rmz", do.call(c, sapply(seq(1, 4), function(x){ paste(x:4, x,sep="") })), sep="")) # first 3 vars in mz
(DZbLabs 	<- paste("rdz", do.call(c, sapply(seq(1, nv), function(x){ paste(x:nv, x,sep="") })), sep="")) # all vars in dz

# Labels For means  #NAs because the first variable is ordinal
(LabM	<- paste("m",2:nv,sep=""))
MLabs		<-c(NA,LabM,NA,LabM) 

# Labels for Standard deviations (SDs), #NAs because the first variable is ordinal
(LabSD	<- paste("sd",2:nv,sep=""))
SDLabs	<-c(NA,LabSD,NA,LabSD) #NA because the first variable is categorical

# Starting values for correlations
jiggle	<-rnorm(nlower, mean = 0, sd = .1)

(StWithinpersonv 	<-vechs(cor(mzData_st[,1:4],use="complete")))
(StWithinpersonP 	<-vechs(cor(mzData_st[,5:nv],use="complete")))
(StWithinpersonvP	<-as.vector(cor(mzData_st[,5:nv],mzData_st[,1:4],use="complete")))
(StBetweenMZ  	<-vech(cor(mzData_st[,1:4],mzData_st[,(nv+1):(nv+4)],use="complete")))
(StBetweenDZ  	<-vech(cor(dzData_st[,1:nv],dzData_st[,(nv+1):ntv],use="complete")))

# Starting values for M & SDs
(Stmean	<-colMeans(mzData_st[2:7],na.rm=TRUE))
(Stsd 	<-sapply(mzData_st[2:7],sd, na.rm=TRUE))

StM 		<-c(0, Stmean, 0, Stmean)
StSD  	<-c(1, Stsd, 1, Stsd)

# Free parameters
(Pat  	<- c( rep(F,nvo), rep(TRUE, nvc)))

# -------------------------------------------------------------------------------------------------------------------
#
# MODEL 1:	Constrained Correlation Model
#
# -------------------------------------------------------------------------------------------------------------------
# 1) Fits a constrained Polychoric correlation model
# TH same across twins and across zyg groups
# Age effect is different acros variables, but same across thresholds within variables (if c>2)
# There is one overall set of rPH and the x-trait x-twin correlations are symmetric
# ------------------------------------------------------------------------------------------------------------------------------

# CREATE LABELS & START VALUES as objects(to ease specification)
# I allow the threshold and increments not to vary by birth order or zygosity and specify only one set
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 for mz twin 1 (mz)

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')

ThPat		<-c(TRUE,TRUE,TRUE,TRUE)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=T, values=.2, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=T, values=.2, labels=LabCovS, name="BsexTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=T, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=F, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2), name="expThres")

# Matrix for expected Means
# First & 10th means are false, i.e. means for SO  variable are set to 0
Mean	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(Pat, Pat), values=c(StM), labels=c(MLabs), name="expm" ) #first variable is not free to be estimated

# Matrices to store SD and cor matrices
SD	<-mxMatrix( type="Diag", nrow=ntv, ncol=ntv, free=c(Pat,Pat), values=c(StSD), labels=c(SDLabs), lbound=0, name="sd")
Rphv	<-mxMatrix("Stand", 4, 4, free = TRUE, values = StWithinpersonv, labels=rphvLabs, name="withinv") 
RphP	<-mxMatrix("Stand", 3, 3, free = TRUE, values = StWithinpersonP, labels=rphPLabs, name="withinP") 
RphvP	<-mxMatrix("Full", 3, 4, free = TRUE, values = StWithinpersonvP, labels=rphvPLabs, name="withinvP") 
RphPv	<-mxMatrix("Full", 4, 3, free = TRUE, values = StWithinpersonvP, labels=rphvPLabs, byrow=T, name="withinPv") 
Rph	<-mxAlgebra(rbind(cbind(withinv,withinPv), cbind(withinvP, withinP)), name="within")

MZbv	<-mxMatrix("Symm", 4, 4, free = TRUE, values = StBetweenMZ, labels=MZbLabs, name="BetweenMZv") 
MZb	<-mxAlgebra(rbind(cbind(BetweenMZv,withinPv), cbind(withinvP, withinP)), name="BetweenMZ")

DZb	<-mxMatrix("Symm", nv, nv, free = T, values = StBetweenDZ, labels=DZbLabs, name="BetweenDZ") 

CorMZ	<-mxAlgebra(rbind(cbind(within,BetweenMZ), cbind(BetweenMZ, within)), name="expCorMZ")
CorDZ	<-mxAlgebra(rbind(cbind(within,BetweenDZ), cbind(BetweenDZ, within)), name="expCorDZ")

# Expected covariance matrices
# SD & Correlation matrices from above
covMZ	<-mxAlgebra( expression=sd %&% expCorMZ, name="expCovMZ")
covDZ	<-mxAlgebra( expression=sd %&% expCorDZ, name="expCovDZ")

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objmz  	<- mxExpectationNormal( covariance="expCovMZ", means="expm", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2") )
objdz		<- mxExpectationNormal( covariance="expCovDZ", means="expm", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2") )

fitFunction 	<- mxFitFunctionML()
#fitFunction 	<- mxFitFunctionWLS()

# Combine Groups
modelMZ	<- mxModel(obsage1, obsage2, obssex1, obssex2, betaA, betaS, Thr, inc, Thres, Mean,  SD, Rph, Rphv, RphP, RphvP, RphPv, MZb, MZbv, CorMZ, covMZ, dataMZ, objmz, fitFunction, name="MZ")
modelDZ	<- mxModel(obsage1, obsage2, obssex1, obssex2, betaA, betaS, Thr, inc, Thres, Mean,  SD, Rph, Rphv, RphP, RphvP, RphPv, DZb, CorDZ, covDZ, dataDZ, objdz, fitFunction, name="DZ")
ciW		<- mxCI('MZ.within')
ciBmz		<- mxCI('MZ.BetweenMZ')
ciBdz		<- mxCI('DZ.BetweenDZ')
minus2ll	<- mxAlgebra(expression=MZ.objective + DZ.objective, name="m2LL")
obj		<- mxFitFunctionAlgebra("m2LL")
corgModel	<- mxModel('corg', modelMZ, modelDZ, minus2ll, obj, ciW, ciBmz, ciBdz) 

# ------------------------------------------------------------------------------
# RUN Multivariate cor Model
corgFit    <- mxTryHardOrdinal(corgModel, extraTries=1, intervals=T)
(corgSumm  <- summary(corgFit))

# Generate Output
mxEval(MZ.within, corgFit)
mxEval(MZ.BetweenMZ, corgFit)
mxEval(DZ.BetweenDZ, corgFit)
mxEval(MZ.expCovMZ, corgFit)
mxEval(DZ.expCovDZ, corgFit)

#*******************************************************************************************************
# __(IIa)_________________________________________________________________________________________________
# Phenotypic Covariance Model across Latent Constructs  with Vict
# Restrictions: means and variances equated across birth-order & zygosity groups;
# One set of factor loadings; one set of correltional paths between the factors; one set of error terms
# We estimate the factor variances, giving them a scale by fixing the loading on the 1st variable to 1
# This model specifies a full var/cov structure between the latent factors for MZ and DZ twins 
#______________________________________________________________________________________________________

nv			<- 7				# number of variables for a twin = 1 in Univariate
nvo 			<- 1     			# number of ordinal variables per twin
nvc 			<- nv-nvo  			# number of continuous variables per twin
poso 			<- nvo 			# position where ordinal variables start
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nth			<- 4				# number of max thresholds
nlower		<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor			<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
ninc 			<- nth-1 			# number of max increments
ncovariates 	<- 2 				# number of covariates

nfact			<- 5				# number of Latent Factors for Mediation Model per twin
nfact2		<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nfcor			<-(nfact*(nfact+1)/2)-nfact	# number of free elements in a correlation matrix nfact*nfcat

Groups		<- c("mz", "dz")
#Vars			<- c('SO','Dep','Anx','Vict','PRSSO','PRSDep','PRSAnx')
Vars			<- c('SO','Dep','Anx','Vict','PRSSO','PRSDep','PRSAnx')
selVars		<- c('SO1','Dep1','Anx1','Vict1','PRSSO1','PRSDep1','PRSAnx1',
			     'SO2','Dep2','Anx2','Vict2','PRSSO2','PRSDep2','PRSAnx2')
useVars		<- c('SO1','Dep1','Anx1','Vict1','PRSSO1','PRSDep1','PRSAnx1',
			     'SO2','Dep2','Anx2','Vict2','PRSSO2','PRSDep2','PRSAnx2','age1','sex1','age2','sex2')

mzData		<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , useVars)
dzData		<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , useVars)

psych::describe(mzData)
psych::describe(dzData)

mzData$SO1[is.na(mzData$age1)] <- NA
mzData$SO2[is.na(mzData$age2)] <- NA
dzData$SO1[is.na(dzData$age1)] <- NA
dzData$SO2[is.na(dzData$age2)] <- NA

mzData$SO1[is.na(mzData$sex1)] <- NA
mzData$SO2[is.na(mzData$sex2)] <- NA
dzData$SO1[is.na(dzData$sex1)] <- NA
dzData$SO2[is.na(dzData$sex2)] <- NA

mzData$age1[is.na(mzData$age1)] <- 999
mzData$age2[is.na(mzData$age2)] <- 999
dzData$age1[is.na(dzData$age1)] <- 999
dzData$age2[is.na(dzData$age2)] <- 999

mzData$sex1[is.na(mzData$sex1)] <- 999
mzData$sex2[is.na(mzData$sex2)] <- 999
dzData$sex1[is.na(dzData$sex1)] <- 999
dzData$sex2[is.na(dzData$sex2)] <- 999

psych::describe(mzData)
psych::describe(dzData)


# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)

(Stmean	<-colMeans(mzData[,2:nv],na.rm=TRUE))
StM 		<-c(0, Stmean, 0, Stmean)

(LabM	<- paste("m",2:nv,sep=""))
MLabs		<-c(NA,LabM,NA,LabM) 

#(LabEr	<- paste("e",1:nv,sep=""))
(LabEr	<-c("e1","e2","e2","e4","e5","e6","e6"))

# Create Labels for the Factor parameters
(sdLabs	<- paste("sd",1:nfact,sep=""))	# SD
(rphvLabs	<- paste("rv",1:3,sep="")) #Within person - for the first 3 factors - different in mz
(rphPLabs	<- paste("rP",1,sep="")) #Within person - for the last 2 polygenic factors - same in mz
(rphvPLabs	<- paste("rvP",1:6,sep="")) #Within person - for the cross PRS-variable correlations - same in mz
(MZbLabs 	<- paste("rmz", do.call(c, sapply(seq(1, 3), function(x){ paste(x:3, x,sep="") })), sep="")) # first 2 vars in mz
(DZbLabs 	<- paste("rdz", do.call(c, sapply(seq(1, nfact), function(x){ paste(x:nfact, x,sep="") })), sep="")) # all vars in dz

# Create Labels for the Factor Loadings (1st loadings fixed to 1)

PatFl	<- c(F,F,F,F,F,F,F,			
	     F,T,F,F,F,F,F,
	     F,F,F,F,F,F,F,
	     F,F,F,F,F,F,F,
	     F,F,F,F,F,F,T)

StFl	<- c(1,0,0,0,0,0,0,
	     0,.5,1,0,0,0,0,
	     0,0,0,1,0,0,0,
	     0,0,0,0,1,0,0,
	     0,0,0,0,0,1,.5)

LabFl	<- c('l1',NA,NA,NA,NA,NA,NA,
	      NA,'l2','l3',NA,NA,NA,NA,
	      NA,NA,NA,'l4',NA,NA,NA,
	      NA,NA,NA,NA,'l5',NA,NA,
	      NA,NA,NA,NA,NA,'l6','l7')

# Free parameters
(Pat  <- c( rep(FALSE,nvo), rep(TRUE, nvc)))

(StWithinpersonv 	<-c(.2))
(StWithinpersonP 	<-c(.2))
(StWithinpersonvP	<-c(.1))
(StBetweenMZ  	<-c(.5,.2,.2,.5,.2,.5))
(StBetweenDZ  	<-c(.3,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.5,.1,.5))
(PatBetweenDZ  	<-c(T,T,T,T,T,T,T,T,T,T,T,T,F,T,F))

# ______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
# ______________________________________________________________________________________________________

Mean	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(Pat, Pat), values=c(StM), labels=c(MLabs), name="expm" ) #first variable is not free to be estimated

# I constrain the threshold and increments to be equal across birth order or zygosity and specify only one set
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 for mz twin 1 (mz)

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')

ThPat		<-c(TRUE,TRUE,TRUE,TRUE)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.05, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.05, labels=LabCovS, name="BsexTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=F, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2), name="expThres")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Ze75		<-mxMatrix("Zero", nv, nfact, free=F, name="Z75")
LoadTw	<-mxAlgebra(rbind(cbind(FactL,Z75), cbind(Z75, FactL)), name="FactLTw")

ErPath	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=c(F,T,T,F,F,T,T), values=c(0,.5,.5,0,0,.5,.5), labels=LabEr, name="Erp" )
Er		<-mxAlgebra(Erp %*% t(Erp), name="Error")
Ze7		<-mxMatrix("Zero", nv, nv, free=F, name="Z7")
ErTw		<-mxAlgebra(rbind(cbind(Error,Z7), cbind(Z7, Error)), name="ErrorTw")
 
# elements for the SD of Factors
Id2	<-mxMatrix("Iden", 2, 2, free=F, name="I2")
sdF	<-mxMatrix("Diag", nfact, nfact, free=c(F,T,T,T,T), values=1, labels=sdLabs, name="SDf") 
sdFTw	<-mxAlgebra(I2 %x% SDf, name="SDftwin")

# elements for the correlations of Factors
Rphv	<-mxMatrix("Stand", 3, 3, free = TRUE, values = StWithinpersonv, labels=rphvLabs, lbound=-.999, ubound=.999, name="withinv") 
RphP	<-mxMatrix("Stand", 2, 2, free = TRUE, values = StWithinpersonP, labels=rphPLabs, lbound=-.999, ubound=.999, name="withinP") 
RphvP	<-mxMatrix("Full", 2, 3, free = TRUE, values = StWithinpersonvP, labels=rphvPLabs, lbound=-.999, ubound=.999, name="withinvP") 
RphPv	<-mxMatrix("Full", 3, 2, free = TRUE, values = StWithinpersonvP, labels=rphvPLabs, lbound=-.999, ubound=.999, byrow=T, name="withinPv") 
Rph	<-mxAlgebra(rbind(cbind(withinv,withinPv), cbind(withinvP, withinP)), name="Rwithin")

MZbv	<-mxMatrix("Symm", 3, 3, free = TRUE, values = StBetweenMZ, labels=MZbLabs, lbound=-.999, ubound=.999, name="BetweenMZv") 
MZb	<-mxAlgebra(rbind(cbind(BetweenMZv,withinPv), cbind(withinvP, withinP)), name="RbetweenMZ")

DZb		<-mxMatrix("Symm", nfact, nfact, free=PatBetweenDZ, values=StBetweenDZ, labels=DZbLabs, lbound=-.999, ubound=.999, name="RbetweenDZ") 
FactCorMZ	<-mxAlgebra(rbind(cbind(Rwithin,RbetweenMZ), cbind(RbetweenMZ, Rwithin)), name="RMZ")
FactCorDZ	<-mxAlgebra(rbind(cbind(Rwithin,RbetweenDZ), cbind(RbetweenDZ, Rwithin)), name="RDZ")

# Generate expected Covariance matrices of Factors
FactCovMZ	<-mxAlgebra(SDftwin %&% RMZ , name="expFactCovMZ")
FactCovDZ	<-mxAlgebra(SDftwin %&% RDZ , name="expFactCovDZ")

## This second step then derives the var/cov matrix of the observed/measured variables in terms of the variance/covariances of the latent factors and the Factor Loadings
covMZ		<-mxAlgebra( expression= FactLTw  %&% expFactCovMZ , name="ExpCovMZ" )
covDZ		<-mxAlgebra( expression= FactLTw  %&% expFactCovDZ , name="ExpCovDZ" )

## Finally, we derive the total expected variance/covariances for the measured variables which go in the models
TOTcovMZ	<-mxAlgebra( expression= ExpCovMZ + ErrorTw , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= ExpCovDZ + ErrorTw , name="TOTexpCovDZ" )

# Standardizing parameters **********************

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% expFactCovMZ[1:5,1:5] / TOTexpCovMZ[1:7,1:7])), name="StandFact")

# Standardise error terms of the measured variables
StEr		<-mxAlgebra( expression= sqrt(diag2vec( Error/TOTexpCovMZ[1:7,1:7])), name="StandEr")

# ************************************

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expm", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expm", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))

fitFunction <- mxFitFunctionML()
#fitFunction <- mxFitFunctionWLS()

# Combine Groups
pars1		<-list(Mean, Load, Ze75, Ze7, LoadTw, ErPath, Er, ErTw, Id2, sdF, sdFTw, Rph, Rphv, RphP, RphvP, RphPv)
pars2		<-list(obsage1, obsage2, obssex1, obssex2, betaA, betaS, Thr, inc, Thres)
modelMZ	<-mxModel(pars1, pars2, MZb, MZbv, FactCorMZ, FactCovMZ, covMZ, TOTcovMZ, dataMZ, objMZ, fitFunction, StFL, StEr, name="MZ" )
modelDZ	<-mxModel(pars1, pars2, DZb, FactCorDZ, FactCovDZ, covDZ, TOTcovDZ, dataDZ, objDZ, fitFunction, name="DZ" )
minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cist1		<-mxCI (c ('MZ.Rwithin[2,1]','MZ.Rwithin[3,1]','MZ.Rwithin[4,1]','MZ.Rwithin[5,1]','MZ.Rwithin[3,2]','MZ.Rwithin[4,2]','MZ.Rwithin[5,2]','MZ.Rwithin[4,3]','MZ.Rwithin[5,3]','MZ.Rwithin[5,4]'))
cist2		<-mxCI (c ('MZ.StandFact[1,1]','MZ.StandFact[2,1]','MZ.StandFact[3,1]','MZ.StandFact[4,1]','MZ.StandFact[5,1]','MZ.StandFact[6,1]','MZ.StandFact[7,1]'))
cist3		<-mxCI (c ('MZ.StandEr[1,1]','MZ.StandEr[2,1]','MZ.StandEr[3,1]','MZ.StandEr[4,1]','MZ.StandEr[5,1]','MZ.StandEr[6,1]','MZ.StandEr[7,1]'))
cist4		<-mxCI (c ('MZ.RbetweenMZ[1,1]','MZ.RbetweenMZ[2,1]','MZ.RbetweenMZ[3,1]','MZ.RbetweenMZ[4,1]','MZ.RbetweenMZ[5,1]','MZ.RbetweenMZ[2,2]','MZ.RbetweenMZ[3,2]','MZ.RbetweenMZ[4,2]','MZ.RbetweenMZ[5,2]','MZ.RbetweenMZ[3,3]','MZ.RbetweenMZ[4,3]','MZ.RbetweenMZ[5,3]','MZ.RbetweenMZ[4,4]','MZ.RbetweenMZ[5,4]','MZ.RbetweenMZ[5,5]'))
cist5		<-mxCI (c ('DZ.RbetweenDZ[1,1]','DZ.RbetweenDZ[2,1]','DZ.RbetweenDZ[3,1]','DZ.RbetweenDZ[4,1]','DZ.RbetweenDZ[5,1]','DZ.RbetweenDZ[2,2]','DZ.RbetweenDZ[3,2]','DZ.RbetweenDZ[4,2]','DZ.RbetweenDZ[5,2]','DZ.RbetweenDZ[3,3]','DZ.RbetweenDZ[4,3]','DZ.RbetweenDZ[5,3]','DZ.RbetweenDZ[4,4]','DZ.RbetweenDZ[5,4]','DZ.RbetweenDZ[5,5]'))
PhCModel	<-mxModel("PhC", modelMZ, modelDZ, minus2ll, obj, cist1, cist2, cist3, cist4, cist5) 

# --------------------------------------------------------------------------------------------------------------------------------
# 2 RUN Phenotypic Fact Covariance Model by Zygosity

PhCFit	<-mxTryHardOrdinal(PhCModel, intervals=F)
(PhCSumm	<-summary(PhCFit))


# Generate confidence intervals
PhCCIModel	<-mxModel(PhCModel)
PhCCIFit	<-mxRun(PhCCIModel, intervals=TRUE)
(PhCCISumm	<-summary(PhCCIFit, verbose=T))

capture.output(print(PhCCISumm,row.names=F), file = "summary.PhCModel", append = FALSE)

mxEval(MZ.Rwithin, PhCFit)
mxEval(MZ.RbetweenMZ, PhCFit)
mxEval(DZ.RbetweenDZ, PhCFit)

mxEval(MZ.expFactCovMZ, PhCFit)
mxEval(DZ.expFactCovDZ, PhCFit)

mxEval(MZ.ExpCovMZ, PhCFit)
mxEval(DZ.ExpCovDZ, PhCFit)

mxEval(MZ.TOTexpCovMZ, PhCFit)
mxEval(DZ.TOTexpCovDZ, PhCFit)

mxEval(MZ.FactL, PhCFit)
mxEval(MZ.StandFact, PhCFit)

mxEval(MZ.Error, PhCFit)
mxEval(MZ.StandEr, PhCFit)

####

#*******************************************************************************************************
# __(III)________________________________________________________________________________________________
# Phenotypic Mediation Model, partly constrained, but allowing Twin correlations (MZ/DZ) across factors  
# Restrictions: means and variances equated across birth-order & zygosity groups;
# One set of factor loadings; one set of causal paths between the factors; one set of error terms
# We estimate the factor variances by giving them a scale by fixing the loading on the 1st and 10th variable to 1
# and corresponding residuals to 0
# F1=SO, F2=VICT, F3=PD; F1>F2>F3, F1>F3
#________________________________________________________________________________________________________

nv			<- 4				# number of variables for a twin = 1 in Univariate
nvo 			<- 1     			#number of ordinal variables per twin
nvc 			<- nv-nvo  			#number of continuous variables per twin
poso 			<- nvo 			#position where ordinal variables start
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact			<- 3				# number of Latent Factors for Mediation Model per twin
nfact2		<- 2*nfact			# number of Latent Factors for Mediation Model per pair
nth			<- 4				# number of max thresholds
nlower		<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor			<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
ninc 			<- nth-1 			#number of max increments
ncovariates 	<- 2 				#number of covariates

Groups		<-c("mz", "dz")
Vars			<- c('SO','Vict','Dep','Anx')
selVars		<- c('SO1','Vict1','Dep1','Anx1',
			     'SO2','Vict2','Dep2','Anx2')
useVars		<- c('SO1','Vict1','Dep1','Anx1',
		     	     'SO2','Vict2','Dep2','Anx2','sex1','age1','sex2','age2')

mzData		<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , useVars)
dzData		<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , useVars)

psych::describe(mzData)
psych::describe(dzData)

mzData$SO1[is.na(mzData$age1)] <- NA
mzData$SO2[is.na(mzData$age2)] <- NA
dzData$SO1[is.na(dzData$age1)] <- NA
dzData$SO2[is.na(dzData$age2)] <- NA

mzData$SO1[is.na(mzData$sex1)] <- NA
mzData$SO2[is.na(mzData$sex2)] <- NA
dzData$SO1[is.na(dzData$sex1)] <- NA
dzData$SO2[is.na(dzData$sex2)] <- NA

mzData$age1[is.na(mzData$age1)] <- 999
mzData$age2[is.na(mzData$age2)] <- 999
dzData$age1[is.na(dzData$age1)] <- 999
dzData$age2[is.na(dzData$age2)] <- 999

mzData$sex1[is.na(mzData$sex1)] <- 999
mzData$sex2[is.na(mzData$sex2)] <- 999
dzData$sex1[is.na(dzData$sex1)] <- 999
dzData$sex2[is.na(dzData$sex2)] <- 999

psych::describe(mzData)
psych::describe(dzData)

mzData_st		<- mzData%>% select(all_of(selVars))%>% mutate_if(is.factor, as.numeric)#I convert factors to numeric data so I can create starting values for correlations
dzData_st		<- dzData%>% select(all_of(selVars))%>% mutate_if(is.factor, as.numeric)

psych::describe(mzData_st)
psych::describe(dzData_st)

cor(mzData_st, use=  "complete")
cor(dzData_st, use = "complete")

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)

# Create Labels for Column and Diagonal Matrices
(Stmean	<-colMeans(mzData[,2:nv],na.rm=TRUE))
(LabM		<- paste("m",2:nv,sep=""))
mLabs		<-c(NA,LabM) 
(Pat  	<- c( rep(FALSE,nvo), rep(TRUE, nvc)))

#(LabEr	<- paste("e",1:nv,sep=""))
(LabEr	<- c("e1","e2","e3","e3"))
PatEr		<- c(F,F,TRUE,TRUE)
StEr		<- c(0,0,.5,.5)

# Create Labels for the Factor Loadings (1st loadings fixed to 1)
PatFl	<- c(F,F,F,F,			
	     F,F,F,F,
	     F,F,F,T)

StFl	<- c(1,0,0,0,
	     0,1,0,0,
	     0,0,1,.8)

LabFl	<- c('l1',NA,NA,NA,
	      NA,'l2',NA,NA,
	      NA,NA,'l3','l4')

PatPhC	<- c(F,TRUE,TRUE,
		     F,F,TRUE,
		     F,F,F)

StPhC		<- c(0,.3,.3,
		     0,0,.3,
		     0,0,0)

LabPhC	<- c(NA,'c1on2','c1on3',
		     NA,NA,'c2on3',	
		     NA,NA,NA)	 

LabFactrMZ	<- c('f1rmz','f2rmz','f3rmz')
LabFactrDZ	<- c('f1rdz','f2rdz','f3rdz')

# ______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
# ______________________________________________________________________________________________________

Means		<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(Pat,Pat), values=c(0,Stmean,0,Stmean), labels=c(mLabs,mLabs), name="expMean" ) #first variable is not free to be estimated

# I constrain the threshold and increments to be equal across birth order or zygosity and specify only one set
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 for mz twin 1 (mz)

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')

ThPat		<-c(TRUE,TRUE,TRUE,TRUE)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.05, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.05, labels=LabCovS, name="BsexTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=F, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2), name="expThres")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Ze43		<-mxMatrix("Zero", nv, nfact, free=F, name="Z43")
LoadTw	<-mxAlgebra(rbind(cbind(FactL,Z43), cbind(Z43, FactL)), name="FactLTw")

ErPath	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatEr, values=c(StEr), labels=LabEr, name="Erp" )
Er		<-mxAlgebra(Erp %*% t(Erp), name="Error")
Ze4		<-mxMatrix("Zero", nv, nv, free=F, name="Z4")
ErTw		<-mxAlgebra(rbind(cbind(Error,Z4), cbind(Z4, Error)), name="ErrorTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 4 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )
Ze3		<-mxMatrix("Zero", nfact, nfact, free=F, name="Z3")
PhCausTw	<-mxAlgebra( expression= rbind( (cbind(PhC,Z3)), (cbind(Z3,PhC))  ), name="PhCTw" )

# Define the matrix to hold the Double headed Arrows (Variances) of the 3 latent variables to themselves 
FactSD	<-mxMatrix(type="Diag",	nrow=nfact, ncol=nfact, free=c(F,TRUE,TRUE), values=c(1,1,1), labels=c('F1sd','F2res','F3res'),name="Fsd" )
FactVar	<-mxAlgebra(expression=Fsd %*% t(Fsd),name="Fcov" )
FactSDTw	<-mxAlgebra( expression= rbind( (cbind(Fsd,Z3)), (cbind(Z3,Fsd))  ), name="FsdTw" )

# Define the matrix to hold the Double headed Arrows (CoVariances) between the 3 latent variables (across twins) 
FactrMZ	<-mxMatrix(type="Diag",	nrow=nfact, ncol=nfact, free=TRUE, values=.4, labels=LabFactrMZ, lbound=-.999, ubound=.999, name="FcorMZ" )
FactrDZ	<-mxMatrix(type="Diag",	nrow=nfact, ncol=nfact, free=TRUE, values=.2, labels=LabFactrDZ, lbound=-.999, ubound=.999, name="FcorDZ" )

# Generate Covariance of Latent factor model: MZ
Id3		<-mxMatrix(type="Iden",	nrow=nfact, ncol=nfact, name="I3" )
FactCorMZ	<-mxAlgebra( expression= rbind( (cbind(I3,FcorMZ)), (cbind(FcorMZ,I3))  ), name="TwFcorMZ" )
FactCovMZ	<-mxAlgebra( expression= FsdTw %&% TwFcorMZ, name="FCovMZ" )

# Generate Covariance of Latent factor model: DZ
FactCorDZ	<-mxAlgebra( expression= rbind( (cbind(I3,FcorDZ)), (cbind(FcorDZ,I3))  ), name="TwFcorDZ" )
FactCovDZ	<-mxAlgebra( expression= FsdTw %&% TwFcorDZ, name="FCovDZ" )

Id6		<-mxMatrix(type="Iden",	nrow=nfact2, ncol=nfact2, name="I6" )
PhFactCovMZ	<-mxAlgebra( expression= solve(I6-PhCTw) %&% FCovMZ , name="ExpFactCovMZ" )
PhFactCovDZ	<-mxAlgebra( expression= solve(I6-PhCTw) %&% FCovDZ , name="ExpFactCovDZ" )

covMZ		<-mxAlgebra( expression= FactLTw  %&% ExpFactCovMZ , name="ExpCovMZ" )
covDZ		<-mxAlgebra( expression= FactLTw  %&% ExpFactCovDZ , name="ExpCovDZ" )

## Finally, we derive the total expected variance/covariances for the measured variables which go in the models
TOTcovMZ	<-mxAlgebra( expression= ExpCovMZ + ErrorTw , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= ExpCovDZ + ErrorTw , name="TOTexpCovDZ" )

# Standardizing parameters **********************

# Standardize the causal effects
Stcp1on2	<-mxAlgebra( expression= (PhC[2,1]* sqrt(ExpFactCovMZ[1,1]))/sqrt(ExpFactCovMZ[2,2]) , name="Stand_1on2" )
Stcp1on3	<-mxAlgebra( expression= (PhC[3,1]* sqrt(ExpFactCovMZ[1,1]))/sqrt(ExpFactCovMZ[3,3]) , name="Stand_1on3" )
Stcp2on3	<-mxAlgebra( expression= (PhC[3,2]* sqrt(ExpFactCovMZ[2,2]))/sqrt(ExpFactCovMZ[3,3]) , name="Stand_2on3" )

# Standardize the covariances between latent factors
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I6*ExpFactCovMZ)) %&% ExpFactCovMZ		, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I6*ExpFactCovDZ)) %&% ExpFactCovDZ		, name="FactcorDZ" )

# Standardized Factor Loadings
StFL	<-mxAlgebra( expression= diag2vec( FactL %&% ExpFactCovMZ[1:3,1:3] / TOTexpCovMZ[1:4,1:4]) , name="StandFact" )

# Standardise error terms of the measured variables
StEr		<-mxAlgebra( expression= diag2vec( Error/TOTexpCovMZ[1:4,1:4])		, name="StandEr" )

#Mediated effects
indir_eff1	<-mxAlgebra(expression = Stand_1on2%*%Stand_2on3, name="Indi_eff1")
prop_med1	<-mxAlgebra(expression = Indi_eff1/(Indi_eff1 + Stand_1on3), name="Prop_med1")
tot_eff	<-mxAlgebra(expression = Stand_1on3 + Indi_eff1, name="Tot_eff")

# ***************************************

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2") )
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2") )

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1		<-list(Means, Load, Ze43, Ze4, LoadTw, ErPath, Er, ErTw, PhCaus, Ze3, PhCausTw, FactSD, FactVar, FactSDTw, Id3, Id6)
pars2		<-list(obsage1, obsage2, obssex1, obssex2, betaA, betaS, Thr, inc, Thres)
parsmed	<-list(Stcp1on2, Stcp1on3, Stcp2on3, indir_eff1, prop_med1, tot_eff)
modelMZ	<-mxModel(pars1, pars2, FactrMZ, FactCorMZ, FactCovMZ, PhFactCovMZ, covMZ, TOTcovMZ, dataMZ, objMZ, Rfactmz, fitFunction, StFL, StEr, parsmed, name="MZ" )
modelDZ	<-mxModel(pars1, pars2, FactrDZ, FactCorDZ, FactCovDZ, PhFactCovDZ, covDZ, TOTcovDZ, dataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )
minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cistC		<-mxCI(c('MZ.Stand_1on2','MZ.Stand_1on3','MZ.Stand_2on3','MZ.Indi_eff1','MZ.Prop_med1','MZ.Tot_eff'))
cistC2	<-mxCI(c('MZ.StandFact[3,1]','MZ.StandFact[4,1]',
			   'MZ.StandEr[3,1]','MZ.StandEr[4,1]'))
cistCor	<-mxCI(c('MZ.FcorMZ[1,1]','MZ.FcorMZ[2,2]','MZ.FcorMZ[3,3]','DZ.FcorDZ[1,1]','DZ.FcorDZ[2,2]','DZ.FcorDZ[3,3]'))
PhMModel	<-mxModel("PhM", modelMZ, modelDZ, minus2ll, obj, cistC, cistC2) 
#PhMModel	<-mxModel("PhM", modelMZ, modelDZ, minus2ll, obj, cistC) 
#PhMModel	<-mxModel("PhM", modelMZ, modelDZ, obj, cistC) 

# --------------------------------------------------------------------------------------------------------------------------------
# 3 RUN Phenotypic Mediation Model by Zygosity

PhMFit	<-mxTryHardOrdinal(PhMModel, intervals=F)
#PhMFit	<-mxTryHard(PhMModel, intervals=F,bestInitsOutput=T,showInits=F)
(PhMSumm	<-summary(PhMFit))
#(PhMSumm	<-summary(PhMFit1, verbose=T))


# Generate confidence intervals
PhMCIModel	<-mxModel(PhMModel)
PhMCIFit	<-mxRun(PhMCIModel, intervals=TRUE)
(PhMCISumm	<-summary(PhMCIFit, verbose=F))

capture.output(print(PhMCISumm,row.names=F), file = "summary.PhMModel", append = FALSE)


# Get some output

mxEval(MZ.FactcorMZ, PhMFit)
mxEval(DZ.FactcorDZ, PhMFit)

mxEval(MZ.Fsd, PhMFit)
mxEval(MZ.Fcov, PhMFit)
mxEval(MZ.ExpFactCovMZ, PhMFit)
mxEval(DZ.ExpFactCovDZ, PhMFit)

mxEval(MZ.FcorMZ, PhMFit)
mxEval(DZ.FcorDZ, PhMFit)
mxEval(MZ.FCovMZ, PhMFit)
mxEval(DZ.FCovDZ, PhMFit)

mxEval(MZ.ExpCovMZ, PhMFit)
mxEval(DZ.ExpCovDZ, PhMFit)
mxEval(MZ.TOTexpCovMZ, PhMFit)
mxEval(DZ.TOTexpCovDZ, PhMFit)

mxEval(MZ.FactL, PhMFit)
mxEval(MZ.StandFact, PhMFit)

mxEval(MZ.Erp, PhMFit)
mxEval(MZ.StandEr, PhMFit)

mxEval(MZ.Stand_1on2, PhMFit)
mxEval(MZ.Stand_1on3, PhMFit)
mxEval(MZ.Stand_2on3, PhMFit)

mxEval(MZ.Tot_eff, PhMFit)
mxEval(MZ.Indi_eff1, PhMFit)
mxEval(MZ.Prop_med1, PhMFit)


#****************************************************************************************************************************
# __(IV)______________________________________________________________________________________________
# ACE MODEL for latent factors (SSA, PD, Vict) by zygosity
# NO causal paths between Phenotypic Factors; A, C and E latent factors have Cholesky Structure
# + Asp, Csp and Esp in the bottom with constraints to identify the model on top
# Correlation between Phenotypic Factors only due to shared A, C and E influences
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
# this because applying a constraint on the factor variances of 1 is problematic especially when we model the causal paths.
# To identify the model we constrain Asp, Csp and Esp variance components loading on variables 1 and 10 to zero
#_____________________________________________________________________________________________________________________________


nv			<- 4				# number of variables for a twin = 1 in Univariate
nvo 			<- 1     			# number of ordinal variables per twin
nvc 			<- nv-nvo  			# number of continuous variables per twin
poso 			<- nvo 			# position where ordinal variables start
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nth			<- 4				# number of max thresholds
nlower		<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor			<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
ninc 			<- nth-1 			# number of max increments
ncovariates 	<- 2 				# number of covariates

nfact			<- 3				# number of Latent Factors for Mediation Model per twin
nfact2		<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nfcor			<-(nfact*(nfact+1)/2)-nfact	# number of free elements in a correlation matrix nfact*nfcat

Groups		<- c("mz", "dz")

Vars			<- c('SO','Dep','Anx','Vict')
selVars		<- c('SO1','Dep1','Anx1','Vict1',
			     'SO2','Dep2','Anx2','Vict2')
useVars		<- c('SO1','Dep1','Anx1','Vict1',
			     'SO2','Dep2','Anx2','Vict2','age1','sex1','age2','sex2')

mzData		<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , useVars)
dzData		<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , useVars)

psych::describe(mzData)
psych::describe(dzData)

mzData$SO1[is.na(mzData$age1)] <- NA
mzData$SO2[is.na(mzData$age2)] <- NA
dzData$SO1[is.na(dzData$age1)] <- NA
dzData$SO2[is.na(dzData$age2)] <- NA

mzData$SO1[is.na(mzData$sex1)] <- NA
mzData$SO2[is.na(mzData$sex2)] <- NA
dzData$SO1[is.na(dzData$sex1)] <- NA
dzData$SO2[is.na(dzData$sex2)] <- NA

mzData$age1[is.na(mzData$age1)] <- 999
mzData$age2[is.na(mzData$age2)] <- 999
dzData$age1[is.na(dzData$age1)] <- 999
dzData$age2[is.na(dzData$age2)] <- 999

mzData$sex1[is.na(mzData$sex1)] <- 999
mzData$sex2[is.na(mzData$sex2)] <- 999
dzData$sex1[is.na(dzData$sex1)] <- 999
dzData$sex2[is.na(dzData$sex2)] <- 999

psych::describe(mzData)
psych::describe(dzData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- colMeans(mzData[,2:nv],na.rm=TRUE))
(Stsd 	<- sapply(mzData[,2:nv],sd, na.rm=TRUE))
(PatM		<- c(F,TRUE,TRUE,TRUE))

# Create Labels for Diagonal Matrices
# To identify this model we need to equate the sp effects of var and 2 and fix the Sp of last variable to 0)
(LabEs	<- c(NA,'es2','es2',NA))
(LabAs	<- c(NA,'as2','as2',NA))
(LabCs	<- c(NA,'cs2','cs2',NA))

PatSp		<- c(F,TRUE,TRUE,F)
StSpa		<- c(0,.5,.5,0)
StSpc		<- c(0,.5,.5,0)
StSpe		<- c(0,.5,.5,0)

# all 1st loadings fixed to 1
PatFl		<- c(F,F,F,F,			
		     F,F,T,F,
		     F,F,F,F)

StFl		<- c(1,0,0,0,
		     0,1,.5,0,
		     0,0,0,1)


LabFl		<- c('l1',NA,NA,NA,
		     NA,'l2','l3',NA,
		     NA,NA,NA,'l4')

PatPhC	<- c(F,F,F,
		     F,F,F,
		     F,F,F)

StPhC		<- c(0,0,0,
		     0,0,0,
		     0,0,0)

LabPhC	<- c(NA,NA,NA,
		     NA,NA,NA,	
		     NA,NA,NA)	 

# ______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
# ______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(0,Stmean,0,Stmean), labels=c(mLabs,mLabs), name="expMean") 

# Threshold and covariates
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')

ThPat		<-c(T,T,T,T)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.3, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.4, labels=LabCovS, name="BsexTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2), name="expThres")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )

# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSp, values=StSpa, labels=LabAs, name="as" )
PathsCs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSp, values=StSpc, labels=LabCs, name="cs" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSp, values=StSpe, labels=LabEs, name="es" )
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
covCs		<-mxAlgebra( expression= cs %*% t(cs), name="Cs" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Cs+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAc	<-mxMatrix(type="Lower", nrow=nfact, ncol=nfact, free=TRUE, values=.6, labels=c("a11","a21","a31","a22","a32","a33"), name="a_c" )
PathsCc	<-mxMatrix(type="Lower", nrow=nfact, ncol=nfact, free=TRUE, values=.3, labels=c("c11","c21","c31","c22","c32","c33"), name="c_c" )
PathsEc	<-mxMatrix(type="Lower", nrow=nfact, ncol=nfact, free=TRUE, values=.6, labels=c("e11","e21","e31","e22","e32","e33"), name="e_c" )
covAc		<-mxAlgebra( expression= a_c %*% t(a_c), name="Ac" )
covCc		<-mxAlgebra( expression= c_c %*% t(c_c), name="Cc" )
covEc		<-mxAlgebra( expression= e_c %*% t(e_c), name="Ec" )
covPc		<-mxAlgebra( expression= Ac+Cc+Ec, name="Vc" )

# Generate Covariance of Latent factor model Including Causal Paths between factors
Id3		<-mxMatrix(type="Iden",	nrow=nfact, ncol=nfact, name="I3" )
covFAc	<-mxAlgebra( expression= solve(I3-PhC) %&% Ac, name ="FAc") #(I3-PhC) gives the expression for the removal of the loop effect of causal relationships between the factors (1-4).
covFCc	<-mxAlgebra( expression= solve(I3-PhC) %&% Cc, name ="FCc")
covFEc	<-mxAlgebra( expression= solve(I3-PhC) %&% Ec, name ="FEc")
covFc		<-mxAlgebra( expression= FAc+FCc+FEc, name="FVc" )

# Constraint on total variance of Ordinal variable (A+C+E=1)
varL		<- mxConstraint( expression=FVc[1,1]==1, name="L" )

# Var-Cov of measured vars in terms of latent factors and AC, Cc, and Ec
FcovMZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, FAc+FCc), cbind(FAc+FCc, FVc) )) , name="expFCovMZ" )#This traces the path from vars to factors and back to vars
FcovDZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, .5%x%FAc+FCc), cbind(.5%x%FAc+FCc, FVc) )) , name="expFCovDZ" )

SpcovMZ	<-mxAlgebra( expression= rbind (cbind(Vs, As+Cs), cbind(As+Cs, Vs) ) , name="expSpCovMZ" )
SpcovDZ	<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As+Cs), cbind(.5%x%As+Cs, Vs) ) , name="expSpCovDZ" )

TOTcovMZ	<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )

# *******************************************************************************************************
# Calculator

# Standardize the Total var/covariances matrices of the observed variables
Id8		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I8" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )

# Standardize the Common Effects
stcovAc	<-mxAlgebra( expression= FAc/FVc, name="stAc" )
stcovCc	<-mxAlgebra( expression= FCc/FVc, name="stCc" )
stcovEc	<-mxAlgebra( expression= FEc/FVc, name="stEc" )

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= As/( (FactL %&% FVc) +Vs), name="stAs" )
stcovCs	<-mxAlgebra( expression= Cs/( (FactL %&% FVc) +Vs), name="stCs" )
stcovEs	<-mxAlgebra( expression= Es/( (FactL %&% FVc) +Vs), name="stEs" )

# Standardized Effects of Individual variables from the factors (Variance components) above
stAvar	<-mxAlgebra( expression= (FactL %&% FAc)/( (FactL %&% FVc) +Vs), name="stAvariables" )
stCvar	<-mxAlgebra( expression= (FactL %&% FCc)/( (FactL %&% FVc) +Vs), name="stCvariables" )
stEvar	<-mxAlgebra( expression= (FactL %&% FEc)/( (FactL %&% FVc) +Vs), name="stEvariables" )

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% FVc / TOTexpCovMZ[1:4,1:4])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2") )
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2") )

fitFunction <- mxFitFunctionML()
#fitFunction <- mxFitFunctionWLS()
 
# Combine Groups
pars1		<-list(Means,Load,LoadTw,PhCaus,PathsAs,PathsAs,PathsCs,PathsEs,covAs,covCs,covEs,covPs,Id3,Id2,Id8, stAvar, stCvar, stEvar)
pars2		<-list(PathsAc,PathsCc,PathsEc,covAc,covCc,covEc,covPc,covFAc,covFCc,covFEc,covFc,stcovAc,stcovCc,stcovEc, stcovAs, stcovCs, stcovEs)
pars3		<-list(obsage1, obsage2, obssex1, obssex2, betaA, betaS, Thr, inc, Thres)
modelMZ	<-mxModel(pars1, pars2, pars3, FcovMZ, SpcovMZ, TOTcovMZ, dataMZ, objMZ, Rfactmz, fitFunction, StFL, varL, name="MZ" )
modelDZ	<-mxModel(pars1, pars2, pars3, FcovDZ, SpcovDZ, TOTcovDZ, dataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )
minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cistFL	<-mxCI (c ('MZ.StandFact'))
cistFc	<-mxCI (c ('MZ.stAc[1,1]','MZ.stAc[2,1]','MZ.stAc[3,1]','MZ.stAc[2,2]','MZ.stAc[3,2]','MZ.stAc[3,3]',
				'MZ.stCc[1,1]','MZ.stCc[2,1]','MZ.stCc[3,1]','MZ.stCc[2,2]','MZ.stCc[3,2]','MZ.stCc[3,3]',
				'MZ.stEc[1,1]','MZ.stEc[2,1]','MZ.stEc[3,1]','MZ.stEc[2,2]','MZ.stEc[3,2]','MZ.stEc[3,3]') ) 	# standardized var comp from Common feactors	
cistVs	<-mxCI (c ('MZ.stAs[1,1]','MZ.stAs[2,2]','MZ.stAs[3,3]','MZ.stAs[4,4]',
				'MZ.stCs[1,1]','MZ.stCs[2,2]','MZ.stCs[3,3]','MZ.stCs[4,4]',
				'MZ.stEs[1,1]','MZ.stEs[2,2]','MZ.stEs[3,3]','MZ.stEs[4,4]') ) 	# standardized var comp from specific Factors
#cistvars	<-mxCI (c ('MZ.stAvariables','MZ.stCvariables','MZ.stEvariables'))
ACEModel	<-mxModel("ace", pars1, pars2, modelMZ, modelDZ, minus2ll, obj, cistFc, cistFL, cistVs) 

# --------------------------------------------------------------------------------------------------------------------------------
# 4 RUN ACE Factor Model: Cholesky (by Zygosity)

ACEFit		<-mxTryHardOrdinal(ACEModel, intervals=F)
#ACEFit		<-mxTryHard(ACEModel, intervals=F, bestInitsOutput=TRUE, showInits=TRUE)
(ACESumm		<-summary(ACEFit, verbose=T))

# Get some output

mxEval(MZ.Vs, ACEFit)

mxEval(MZ.FAc, ACEFit)
mxEval(MZ.FCc, ACEFit)
mxEval(MZ.FEc, ACEFit)
mxEval(MZ.FVc, ACEFit)
mxEval(MZ.PhC, ACEFit)

mxEval(MZ.stAc, ACEFit)
mxEval(MZ.stCc, ACEFit)
mxEval(MZ.stEc, ACEFit)

mxEval(MZ.stAs, ACEFit)
mxEval(MZ.stCs, ACEFit)
mxEval(MZ.stEs, ACEFit)

mxEval(MZ.StandFact, ACEFit)

mxEval(MZ.stAvariables, ACEFit)
mxEval(MZ.stCvariables, ACEFit)
mxEval(MZ.stEvariables, ACEFit)

#------------------------------------------------------------------------
# Drop Ccommon pathsfrom previous model 
# I.e. Cholesky model without causal paths but with C-spec dropped
# -----------------------------------------------------------------------
subModel1	<- mxModel(ACEFit, name="sub1")
subModel1	<- omxSetParameters(subModel1, labels=c('cs2','c11','c21','c31','c22','c32','c33'), free=FALSE, values=0)
subFit1	<- mxTryHardOrdinal(subModel1, intervals=T)
(subSum1	<- summary(subFit1))

mxCompare(ACEFit, subFit1)

mxEval(MZ.stAc, subFit1)
mxEval(MZ.stCc, subFit1)
mxEval(MZ.stEc, subFit1)


# 
#****************************************************************************************************************************
# __(Va)_____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL for SO>PD
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Factors: F1>F2>F3 & F1>F3; F1=PRS, F2=SSA, F3=PD
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
#_____________________________________________________________________________________________________________________________

nv		<- 4				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 3				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nvo 		<- 1     			# number of ordinal variables per twin
nvc 		<- nv-nvo  			# number of continuous variables per twin
poso 		<- nvo 			# position where ordinal variables start
nth		<- 4				# number of max thresholds
 nv*nv
ninc 		<- nth-1 			# number of max increments
ncovariates <- 2 				# number of covariates
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups	<- c("mz", "dz")
Vars		<- c('PRSSO','SO','Dep','Anx')
selVars	<- c('PRSSO1','SO1','Dep1','Anx1',
		     'PRSSO2','SO2','Dep2','Anx2')
useVars	<- c('PRSSO1','SO1','Dep1','Anx1',
		     'PRSSO2','SO2','Dep2','Anx2','age1','sex1','age2','sex2')

mzData		<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , useVars)
dzData		<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , useVars)

psych::describe(mzData)
psych::describe(dzData)

mzData$SO1[is.na(mzData$age1)] <- NA
mzData$SO2[is.na(mzData$age2)] <- NA
dzData$SO1[is.na(dzData$age1)] <- NA
dzData$SO2[is.na(dzData$age2)] <- NA

mzData$SO1[is.na(mzData$sex1)] <- NA
mzData$SO2[is.na(mzData$sex2)] <- NA
dzData$SO1[is.na(dzData$sex1)] <- NA
dzData$SO2[is.na(dzData$sex2)] <- NA

mzData$age1[is.na(mzData$age1)] <- 999
mzData$age2[is.na(mzData$age2)] <- 999
dzData$age1[is.na(dzData$age1)] <- 999
dzData$age2[is.na(dzData$age2)] <- 999

mzData$sex1[is.na(mzData$sex1)] <- 999
mzData$sex2[is.na(mzData$sex2)] <- 999
dzData$sex1[is.na(dzData$sex1)] <- 999
dzData$sex2[is.na(dzData$sex2)] <- 999

psych::describe(mzData)
psych::describe(dzData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- c(5,0,4.5,4.1))
(PatM		<- c(TRUE,F,TRUE,TRUE))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEs	<- c('es1','es2','es3','es3'))
(LabAs	<- c('as1','as2','as3','as3'))
(LabCs	<- c('cs1','cs2','cs3','cs3'))

PatSpe	<- c(F,F,TRUE,TRUE)
PatSpac	<- c(F,F,TRUE,TRUE)
StSpa		<- c(0,0,.5,.5)
StSpc		<- c(0,0,.5,.5)
StSpe		<- c(0,0,.5,.5)

# all 1st loadings fixed to 1
PatFl		<- c(F,F,F,F,			
		     F,F,F,F,
		     F,F,F,T)

StFl		<- c(1,0,0,0,
		     0,1,0,0,
		     0,0,1,.5)

LabFl		<- c('l1',NA,NA,NA,
	 	     NA,'l2',NA,NA,
	 	     NA,NA,'l3','l4')

PatPhC	<- c(F,T,T,
		     F,F,T,
		     F,F,F)

StPhC		<- c(0,.3,.3,
		     0,0,.3,
		     0,0,0)

LabPhC	<- c(NA,'c1on2','c1on3',
		     NA,NA,'c2on3',
		     NA,NA,NA)	 

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmean,Stmean), labels=c(mLabs,mLabs), name="expMean") 

# Threshold and covariates
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')

ThPat		<-c(T,T,T,T)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.3, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.4, labels=LabCovS, name="BsexTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2), name="expThres")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )

# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpa, labels=LabAs, name="as" )
PathsCs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpc, labels=LabCs, name="cs" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpe, labels=LabEs, name="es" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
covCs		<-mxAlgebra( expression= cs %*% t(cs), name="Cs" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Cs+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("ac22","ac32","ac33"), name="ac" ) # Component paths for factors 2 and 3
PathsCcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("cc22","cc32","cc33"), name="cc" )
PathsEcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(.8,0,.8), labels=c("ec22","ec32","ec33"), name="ec" )
PathsP11	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11", name="pc" ) # SD path for factor 3 (the PRS factor)
Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  #Padding

Ze21		<-mxMatrix(type="Zero",	nrow=2, ncol=1, free=F, name="Z21" )  #Padding
Ze12		<-mxMatrix(type="Zero",	nrow=1, ncol=2, free=F, name="Z12" )  #Padding
covAcsub	<-mxAlgebra( expression= ac %*% t(ac), name="Acsub" )
covCcsub	<-mxAlgebra( expression= cc %*% t(cc), name="Ccsub" )
covEcsub	<-mxAlgebra( expression= ec %*% t(ec), name="Ecsub" )
covPcsub	<-mxAlgebra( expression= Acsub+Ccsub+Ecsub, name="Vcsub" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)

covPc11	<-mxAlgebra( expression= pc %*% t(pc), name="Pc11" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components

covPc		<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12,Vcsub)), name="Vc") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPcMz	<-mxAlgebra(cbind(rbind(Pc11,Z21) ,rbind(Z12,Acsub+Ccsub)), name="Vcmz") #I specify the MZ between-twin covariance - excluding E parameters
covPcDz	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z21), rbind(Z12,.5%x%Acsub+Ccsub)), name="Vcdz") #I specify the DZ between-twin covariance - specifying half of A and excluding E

# Generate Covariance of Latent factor model Including Causal Paths between factors
Id3		<-mxMatrix(type="Iden",	nrow=3, ncol=3, free=F, name="I3" )
covFVc	<-mxAlgebra( expression= solve(I3-PhC) %&% Vc, name ="FVc")
covFcMz	<-mxAlgebra( expression= solve(I3-PhC) %&% Vcmz, name ="Fcmz")
covFcDz	<-mxAlgebra( expression= solve(I3-PhC) %&% Vcdz, name ="Fcdz")

# Constraint on total variance of Ordinal variable (A+C+E=1)
varL1		<- mxConstraint( expression=FVc[2,2]==1, name="L1" )
#varL2	<- mxConstraint( expression=Fcmz[1:3,1]==FVc[1:3,1], name="L2" )
#varL3	<- mxConstraint( expression=Fcdz[1,1]==.5, name="L3" )

FcovMZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcmz), cbind(Fcmz, FVc))) , name="expFCovMZ" )
FcovDZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcdz), cbind(Fcdz, FVc))) , name="expFCovDZ" )

SpcovMZ	<-mxAlgebra( expression= rbind (cbind(Vs, As+Cs), cbind(As+Cs, Vs)) , name="expSpCovMZ" )
SpcovDZ	<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As+Cs), cbind(.5%x%As+Cs, Vs)) , name="expSpCovDZ" )

TOTcovMZ	<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )
# *******************************************************************************************************
# Calculator

# Standardize the causal effects
Stcp1on2	<-mxAlgebra( expression= (PhC[2,1]* sqrt(FVc[1,1]))/sqrt(FVc[2,2]) , name="Stand_1on2" )
Stcp1on3	<-mxAlgebra( expression= (PhC[3,1]* sqrt(FVc[1,1]))/sqrt(FVc[3,3]) , name="Stand_1on3" )
Stcp2on3	<-mxAlgebra( expression= (PhC[3,2]* sqrt(FVc[2,2]))/sqrt(FVc[3,3]) , name="Stand_2on3" )

# Standardize the Total var/covariances matrices of the observed variables
Id8		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I8" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )

# Phenotypic, A, C and E correlations	
RfactAc	<-mxAlgebra( expression= solve(sqrt(I2*Acsub)) %&% Acsub, name="Ra" )
RfactCc	<-mxAlgebra( expression= solve(sqrt(I2*Ccsub)) %&% Ccsub, name="Rc" )
RfactEc	<-mxAlgebra( expression= solve(sqrt(I2*Ecsub)) %&% Ecsub, name="Re" )
RfactP	<-mxAlgebra( expression= solve(sqrt(I3*FVc)) %&% FVc, name="Rph" )

# Standardize the Common Effects
covFVc22	<-mxAlgebra( expression= FVc[2:3,2:3], name ="FVc22")
stcovAc	<-mxAlgebra( expression= Acsub/FVc22, name="stAc" )
stcovCc	<-mxAlgebra( expression= Ccsub/FVc22, name="stCc" )
stcovEc	<-mxAlgebra( expression= Ecsub/FVc22, name="stEc" )

# Standardised path estimates
StpathAc	<-mxAlgebra( expression= (sqrt(stAc)), name="stpac" )
StpathCc	<-mxAlgebra( expression= (sqrt(stCc)), name="stpcc" )
StpathEc	<-mxAlgebra( expression= (sqrt(stEc)), name="stpec" )

# Algebra to compute Rph-A, Rph-C and Rph-E
RphA		<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[2,1]*sqrt(stAc[2,2]), name = 'Rpha')
RphC		<-mxAlgebra(expression=sqrt(stCc[1,1])*Rc[2,1]*sqrt(stCc[2,2]), name = 'Rphc')
RphE		<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[2,1]*sqrt(stEc[2,2]), name = 'Rphe')

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= sqrt(As/( (FactL %&% FVc) +Vs)), name="stAs" )
stcovCs	<-mxAlgebra( expression= sqrt(Cs/( (FactL %&% FVc) +Vs)), name="stCs" )
stcovEs	<-mxAlgebra( expression= sqrt(Es/( (FactL %&% FVc) +Vs)), name="stEs" )

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% FVc / TOTexpCovMZ[1:4,1:4])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1		<-list(Means,Load,LoadTw,PhCaus,PathsAs,PathsCs,PathsEs,covAs,covCs,covEs,covPs,Id3,Id8,Id2)
pars2		<-list(PathsAcsub,PathsCcsub,PathsEcsub,PathsP11,Ze21,Ze12,Ze11,covAcsub,covCcsub,covEcsub,covPcsub,covPc11,covPc,covPcMz,covPcDz,covFVc,covFVc22,covFcMz,covFcDz)
pars3		<-list(obsage1, obsage2, obssex1, obssex2, betaA, betaS, Thr, inc, Thres)
parsst	<-list(stcovAs, stcovCs, stcovEs, stcovAc, stcovCc, stcovEc, RfactAc, RfactCc, RfactEc, RfactP,RphA,RphC,RphE,StpathAc,StpathCc,StpathEc)
parsmed	<-list(Stcp1on2, Stcp1on3, Stcp2on3)
modelMZ	<-mxModel(pars1, pars2, pars3, parsmed, FcovMZ, SpcovMZ, TOTcovMZ, dataMZ, objMZ, Rfactmz, parsst, fitFunction, StFL, varL1, name="MZ" )
modelDZ	<-mxModel(pars1, pars2, pars3, FcovDZ, SpcovDZ, TOTcovDZ, dataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )

minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cistFL	<-mxCI (c ('MZ.StandFact','MZ.Stand_1on2','MZ.Stand_1on3','MZ.Stand_2on3','MZ.PhC'))
cistVs	<-mxCI (c ('MZ.stAs[3,3]','MZ.stAs[4,4]',
				'MZ.stCs[3,3]','MZ.stCs[4,4]',
				'MZ.stEs[3,3]','MZ.stEs[4,4]') ) 	# standardized var comp from specific Factors
cistVc	<-mxCI (c ('MZ.stAc[1,1]','MZ.stAc[2,1]','MZ.stAc[2,2]','MZ.stCc[1,1]','MZ.stCc[2,1]','MZ.stCc[2,2]','MZ.stEc[1,1]','MZ.stEc[2,1]','MZ.stEc[2,2]') ) 	# standardized var comp for ACE on latent Factors
cistRc	<-mxCI (c ('MZ.Rpha','MZ.Rphc','MZ.Rphe','MZ.Ra','MZ.Rc','MZ.Re','MZ.stpac','MZ.stpcc','MZ.stpec') ) 	
ACEMsModel	<-mxModel("aceMs", pars1, pars2, modelMZ, modelDZ, minus2ll, obj, cistFL, cistVs, cistVc, cistRc) 

# --------------------------------------------------------------------------------------------------------------------------------
# 5a RUN ACEMs Factor Model with phenotypic causal mediation paths by Zygosity

ACEMsFit	<-mxTryHardOrdinal(ACEMsModel, intervals=F)
(ACEMsSumm	<-summary(ACEMsFit, verbose=F))

mxEval(MZ.Acsub, ACEMsFit)
mxEval(MZ.Ccsub, ACEMsFit)
mxEval(MZ.Ecsub, ACEMsFit)
mxEval(MZ.Vcsub, ACEMsFit)
mxEval(MZ.FVc22, ACEMsFit)

#------------------------------------------------------------------------
# Submodel 5aii: Drop correlation path c paths from previous model 
#------------------------------------------------------------------------
# Drop correlation c paths from the model
# -----------------------------------------------------------------------
AEMs2cModel		<- mxModel(ACEMsFit, name="AEMs2c")
AEMs2cModel		<- omxSetParameters(AEMs2cModel, labels=c('cc22','cc32','cc33'), free=FALSE, values=0)
AEMs2cModel		<- omxSetParameters(AEMs2cModel, labels=c('cs3'), free=FALSE, values=0)
AEMs2cFit		<- mxTryHardOrdinal(AEMs2cModel, intervals=T)
(AEMs2cSum		<- summary(AEMs2cFit))

mxCompare(ACEMsFit, AEMs2cFit)

mxEval(MZ.FactcorMZ, AEMs2cFit)
mxEval(DZ.FactcorDZ, AEMs2cFit)

mxEval(MZ.Acsub, AEMs2cFit)
mxEval(MZ.Ccsub, AEMs2cFit)
mxEval(MZ.Ecsub, AEMs2cFit)
mxEval(MZ.Vcsub, AEMs2cFit)
mxEval(MZ.FVc22, AEMs2cFit)

mxEval(MZ.FVc, AEMs2cFit)

mxEval(MZ.stAc, AEMs2cFit)
mxEval(MZ.stCc, AEMs2cFit)
mxEval(MZ.stEc, AEMs2cFit)

mxEval(MZ.stpac, AEMs2cFit)
mxEval(MZ.stpcc, AEMs2cFit)
mxEval(MZ.stpec, AEMs2cFit)

mxEval(MZ.Ra, AEMs2cFit)
mxEval(MZ.Rc, AEMs2cFit)
mxEval(MZ.Re, AEMs2cFit)
mxEval(MZ.Rph, AEMs2cFit)

mxEval(MZ.Rpha, AEMs2cFit) 
mxEval(MZ.Rphc, AEMs2cFit) 
mxEval(MZ.Rphe, AEMs2cFit) 

mxEval(MZ.StandFact, AEMs2cFit)
mxEval(MZ.Stand_1on2[1,1], AEMs2cFit)
mxEval(MZ.Stand_1on3[1,1], AEMs2cFit)
mxEval(MZ.Stand_2on3[1,1], AEMs2cFit)

# 
#****************************************************************************************************************************
# __(Vb)_____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL for PD>SO
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Phenotypic Factors: F1>F2>F3 & F1>F3; F1=PRS, F2=PD, F3=SSA
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
#_____________________________________________________________________________________________________________________________

nv		<- 5				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 3				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nvo 		<- 1     			# number of ordinal variables per twin
nvc 		<- nv-nvo  			# number of continuous variables per twin
poso 		<- nvo 			# position where ordinal variables start
nth		<- 4				# number of max thresholds
ninc 		<- nth-1 			# number of max increments
ncovariates <- 2 				# number of covariates
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups	<- c("mz", "dz")
Vars		<- c('PRSAnx','PRSDep','Dep','Anx','SO')
selVars	<- c('PRSAnx1','PRSDep1','Dep1','Anx1','SO1',
		     'PRSAnx2','PRSDep2','Dep2','Anx2','SO2')
useVars	<- c('PRSAnx1','PRSDep1','Dep1','Anx1','SO1',
		     'PRSAnx2','PRSDep2','Dep2','Anx2','SO2','age1','sex1','age2','sex2')

mzData		<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , useVars)
dzData		<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , useVars)

psych::describe(mzData)
psych::describe(dzData)

mzData$SO1[is.na(mzData$age1)] <- NA
mzData$SO2[is.na(mzData$age2)] <- NA
dzData$SO1[is.na(dzData$age1)] <- NA
dzData$SO2[is.na(dzData$age2)] <- NA

mzData$SO1[is.na(mzData$sex1)] <- NA
mzData$SO2[is.na(mzData$sex2)] <- NA
dzData$SO1[is.na(dzData$sex1)] <- NA
dzData$SO2[is.na(dzData$sex2)] <- NA

mzData$age1[is.na(mzData$age1)] <- 999
mzData$age2[is.na(mzData$age2)] <- 999
dzData$age1[is.na(dzData$age1)] <- 999
dzData$age2[is.na(dzData$age2)] <- 999

mzData$sex1[is.na(mzData$sex1)] <- 999
mzData$sex2[is.na(mzData$sex2)] <- 999
dzData$sex1[is.na(dzData$sex1)] <- 999
dzData$sex2[is.na(dzData$sex2)] <- 999

psych::describe(mzData)
psych::describe(dzData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- c(5,5,4.5,4.1,0))
(PatM		<- c(T,T,T,T,F))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEs	<- c('es1','es1','es3','es3','es5'))
(LabAs	<- c('as1','as1','as3','as3','as5'))
(LabCs	<- c('cs1','cs1','cs3','cs3','cs5'))

PatSpe	<- c(T,T,T,T,F)
PatSpac	<- c(F,F,T,T,F)
StSpa		<- c(0,0,.5,.5,0)
StSpc		<- c(0,0,.5,.5,0)
StSpe		<- c(.5,.5,.5,.5,0)

# all 1st loadings fixed to 1
PatFl		<- c(F,T,F,F,F,			
		     F,F,F,T,F,
		     F,F,F,F,F)

StFl		<- c(1,.5,0,0,0,
		     0,0,1,.5,0,
		     0,0,0,0,1)

LabFl		<- c('l1','l2',NA,NA,NA,
	 	     NA,NA,'l3','l4',NA,
	 	     NA,NA,NA,NA,'l5')

PatPhC	<- c(F,T,T,
		     F,F,T,
		     F,F,F)

StPhC		<- c(0,.3,.3,
		     0,0,.3,
		     0,0,0)

LabPhC	<- c(NA,'c1on2','c1on3',
		     NA,NA,'c2on3',
		     NA,NA,NA)	 

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmean,Stmean), labels=c(mLabs,mLabs), name="expMean") 

# Threshold and covariates
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')

ThPat		<-c(T,T,T,T)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.3, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.4, labels=LabCovS, name="BsexTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2), name="expThres")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )

# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpa, labels=LabAs, name="as" )
PathsCs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpc, labels=LabCs, name="cs" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpe, labels=LabEs, name="es" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
covCs		<-mxAlgebra( expression= cs %*% t(cs), name="Cs" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Cs+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("ac22","ac32","ac33"), name="ac" ) # Component paths for factors 2 and 3
PathsCcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("cc22","cc32","cc33"), name="cc" )
PathsEcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(.8,0,.8), labels=c("ec22","ec32","ec33"), name="ec" )
PathsP11	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11", name="pc" ) # SD path for factor 3 (the PRS factor)
Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  #Padding

Ze21		<-mxMatrix(type="Zero",	nrow=2, ncol=1, free=F, name="Z21" )  #Padding
Ze12		<-mxMatrix(type="Zero",	nrow=1, ncol=2, free=F, name="Z12" )  #Padding
covAcsub	<-mxAlgebra( expression= ac %*% t(ac), name="Acsub" )
covCcsub	<-mxAlgebra( expression= cc %*% t(cc), name="Ccsub" )
covEcsub	<-mxAlgebra( expression= ec %*% t(ec), name="Ecsub" )
covPcsub	<-mxAlgebra( expression= Acsub+Ccsub+Ecsub, name="Vcsub" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)

covPc11	<-mxAlgebra( expression= pc %*% t(pc), name="Pc11" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components

covPc		<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12,Vcsub)), name="Vc") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPcMz	<-mxAlgebra(cbind(rbind(Pc11,Z21) ,rbind(Z12,Acsub+Ccsub)), name="Vcmz") #I specify the MZ between-twin covariance - excluding E parameters
covPcDz	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z21), rbind(Z12,.5%x%Acsub+Ccsub)), name="Vcdz") #I specify the DZ between-twin covariance - specifying half of A and excluding E

# Generate Covariance of Latent factor model Including Causal Paths between factors
Id3		<-mxMatrix(type="Iden",	nrow=3, ncol=3, free=F, name="I3" )
covFVc	<-mxAlgebra( expression= solve(I3-PhC) %&% Vc, name ="FVc")
covFcMz	<-mxAlgebra( expression= solve(I3-PhC) %&% Vcmz, name ="Fcmz")
covFcDz	<-mxAlgebra( expression= solve(I3-PhC) %&% Vcdz, name ="Fcdz")

# Constraint on total variance of Ordinal variable (A+C+E=1)
varL1		<- mxConstraint( expression=FVc[3,3]==1, name="L1" )

FcovMZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcmz), cbind(Fcmz, FVc))) , name="expFCovMZ" )
FcovDZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcdz), cbind(Fcdz, FVc))) , name="expFCovDZ" )

SpcovMZ	<-mxAlgebra( expression= rbind (cbind(Vs, As+Cs), cbind(As+Cs, Vs)) , name="expSpCovMZ" )
SpcovDZ	<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As+Cs), cbind(.5%x%As+Cs, Vs)) , name="expSpCovDZ" )

TOTcovMZ	<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )
# *******************************************************************************************************
# Calculator

# Standardize the causal effects
Stcp1on2	<-mxAlgebra( expression= (PhC[2,1]* sqrt(FVc[1,1]))/sqrt(FVc[2,2]) , name="Stand_1on2" )
Stcp1on3	<-mxAlgebra( expression= (PhC[3,1]* sqrt(FVc[1,1]))/sqrt(FVc[3,3]) , name="Stand_1on3" )
Stcp2on3	<-mxAlgebra( expression= (PhC[3,2]* sqrt(FVc[2,2]))/sqrt(FVc[3,3]) , name="Stand_2on3" )

# Standardize the Total var/covariances matrices of the observed variables
Id10		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I10" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )

# Phenotypic, A, C and E correlations	
RfactAc	<-mxAlgebra( expression= solve(sqrt(I2*Acsub)) %&% Acsub, name="Ra" )
RfactCc	<-mxAlgebra( expression= solve(sqrt(I2*Ccsub)) %&% Ccsub, name="Rc" )
RfactEc	<-mxAlgebra( expression= solve(sqrt(I2*Ecsub)) %&% Ecsub, name="Re" )
RfactP	<-mxAlgebra( expression= solve(sqrt(I3*FVc)) %&% FVc, name="Rph" )

# Standardize the Common Effects
covFVc22	<-mxAlgebra( expression= FVc[2:3,2:3], name ="FVc22")
stcovAc	<-mxAlgebra( expression= Acsub/FVc22, name="stAc" )
stcovCc	<-mxAlgebra( expression= Ccsub/FVc22, name="stCc" )
stcovEc	<-mxAlgebra( expression= Ecsub/FVc22, name="stEc" )

# Standardised path estimates
StpathAc	<-mxAlgebra( expression= (sqrt(stAc)), name="stpac" )
StpathCc	<-mxAlgebra( expression= (sqrt(stCc)), name="stpcc" )
StpathEc	<-mxAlgebra( expression= (sqrt(stEc)), name="stpec" )

# Algebra to compute Rph-A, Rph-C and Rph-E
RphA		<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[2,1]*sqrt(stAc[2,2]), name = 'Rpha')
RphC		<-mxAlgebra(expression=sqrt(stCc[1,1])*Rc[2,1]*sqrt(stCc[2,2]), name = 'Rphc')
RphE		<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[2,1]*sqrt(stEc[2,2]), name = 'Rphe')

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= sqrt(As/( (FactL %&% FVc) +Vs)), name="stAs" )
stcovCs	<-mxAlgebra( expression= sqrt(Cs/( (FactL %&% FVc) +Vs)), name="stCs" )
stcovEs	<-mxAlgebra( expression= sqrt(Es/( (FactL %&% FVc) +Vs)), name="stEs" )

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% FVc / TOTexpCovMZ[1:5,1:5])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1		<-list(Means,Load,LoadTw,PhCaus,PathsAs,PathsCs,PathsEs,covAs,covCs,covEs,covPs,Id3,Id10,Id2)
pars2		<-list(PathsAcsub,PathsCcsub,PathsEcsub,PathsP11,Ze21,Ze12,Ze11,covAcsub,covCcsub,covEcsub,covPcsub,covPc11,covPc,covPcMz,covPcDz,covFVc,covFVc22,covFcMz,covFcDz)
pars3		<-list(obsage1, obsage2, obssex1, obssex2, betaA, betaS, Thr, inc, Thres)
parsst	<-list(stcovAs, stcovCs, stcovEs, stcovAc, stcovCc, stcovEc, RfactAc, RfactCc, RfactEc, RfactP,RphA,RphC,RphE,StpathAc,StpathCc,StpathEc)
parsmed	<-list(Stcp1on2, Stcp1on3, Stcp2on3)
modelMZ	<-mxModel(pars1, pars2, pars3, parsmed, FcovMZ, SpcovMZ, TOTcovMZ, dataMZ, objMZ, Rfactmz, parsst, fitFunction, StFL, varL1, name="MZ" )
modelDZ	<-mxModel(pars1, pars2, pars3, FcovDZ, SpcovDZ, TOTcovDZ, dataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )

minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cistFL	<-mxCI (c ('MZ.StandFact','MZ.Stand_1on2','MZ.Stand_1on3','MZ.Stand_2on3','MZ.PhC'))
cistVs	<-mxCI (c ('MZ.stAs[3,3]','MZ.stAs[4,4]',
				'MZ.stCs[3,3]','MZ.stCs[4,4]',
				'MZ.stEs[1,1]','MZ.stEs[2,2]','MZ.stEs[3,3]','MZ.stEs[4,4]') ) 	# standardized var comp from specific Factors
cistVc	<-mxCI (c ('MZ.stAc[1,1]','MZ.stAc[2,1]','MZ.stAc[2,2]','MZ.stCc[1,1]','MZ.stCc[2,1]','MZ.stCc[2,2]','MZ.stEc[1,1]','MZ.stEc[2,1]','MZ.stEc[2,2]') ) 	# standardized var comp for ACE on latent Factors
cistRc	<-mxCI (c ('MZ.Rpha','MZ.Rphc','MZ.Rphe','MZ.Ra','MZ.Rc','MZ.Re','MZ.stpac','MZ.stpcc','MZ.stpec') ) 	
ACEMsPDModel	<-mxModel("aceMsPD", pars1, pars2, modelMZ, modelDZ, minus2ll, obj, cistFL, cistVs, cistVc, cistRc) 

# --------------------------------------------------------------------------------------------------------------------------------
# 5c RUN ACEMsR Factor Model with phenotypic causal mediation paths by Zygosity

ACEMsPDFit	<-mxTryHardOrdinal(ACEMsPDModel, intervals=F)
(ACEMsPDSumm	<-summary(ACEMsPDFit, verbose=F))

mxEval(MZ.Acsub, ACEMsPDFit)
mxEval(MZ.Ccsub, ACEMsPDFit)
mxEval(MZ.Ecsub, ACEMsPDFit)
mxEval(MZ.Vcsub, ACEMsPDFit)
mxEval(MZ.FVc22, ACEMsPDFit)

#------------------------------------------------------------------------
# Submodel 5cii: Drop correlation path c paths from previous model 
#------------------------------------------------------------------------
# Drop correlation c paths from the model
# -----------------------------------------------------------------------
AEMs4cModel		<- mxModel(ACEMsPDFit, name="AEMs4c")
AEMs4cModel		<- omxSetParameters(AEMs4cModel, labels=c('cc22','cc32','cc33'), free=FALSE, values=0)
AEMs4cModel		<- omxSetParameters(AEMs4cModel, labels=c('cs3'), free=FALSE, values=0)
AEMs4cFit		<- mxTryHardOrdinal(AEMs4cModel, intervals=T)
(AEMs4cSum		<- summary(AEMs4cFit))

mxCompare(ACEMsPDFit, AEMs4cFit)

mxEval(MZ.FactcorMZ, AEMs4cFit)
mxEval(DZ.FactcorDZ, AEMs4cFit)

mxEval(MZ.Acsub, AEMs4cFit)
mxEval(MZ.Ccsub, AEMs4cFit)
mxEval(MZ.Ecsub, AEMs4cFit)
mxEval(MZ.Vcsub, AEMs4cFit)
mxEval(MZ.FVc22, AEMs4cFit)

mxEval(MZ.FVc, AEMs4cFit)

mxEval(MZ.stAc, AEMs4cFit)
mxEval(MZ.stCc, AEMs4cFit)
mxEval(MZ.stEc, AEMs4cFit)

mxEval(MZ.stpac, AEMs4cFit)
mxEval(MZ.stpcc, AEMs4cFit)
mxEval(MZ.stpec, AEMs4cFit)

mxEval(MZ.Ra, AEMs4cFit)
mxEval(MZ.Rc, AEMs4cFit)
mxEval(MZ.Re, AEMs4cFit)
mxEval(MZ.Rph, AEMs4cFit)

mxEval(MZ.Rpha, AEMs4cFit) 
mxEval(MZ.Rphc, AEMs4cFit) 
mxEval(MZ.Rphe, AEMs4cFit) 

mxEval(MZ.StandFact, AEMs4cFit)
mxEval(MZ.Stand_1on2[1,1], AEMs4cFit)
mxEval(MZ.Stand_1on3[1,1], AEMs4cFit)
mxEval(MZ.Stand_2on3[1,1], AEMs4cFit)

# 
#****************************************************************************************************************************
# __(VIa)_____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL for SO>PD incorporating moderation by ELA
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Phenotypic Factors: F1>F2>F3 & F1>F3;
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
#_____________________________________________________________________________________________________________________________

nv		<- 4				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 3				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nvo 		<- 1     			# number of ordinal variables per twin
nvc 		<- nv-nvo  			# number of continuous variables per twin
poso 		<- nvo 			# position where ordinal variables start
nth		<- 4				# number of max thresholds
ninc 		<- nth-1 			# number of max increments
ncovariates <- 3 				# number of covariates
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups	<- c("mz", "dz")
Vars		<- c('PRSSO','SO','Depct','Anxct')
selVars	<- c('PRSSO1','SO1','Depct1','Anxct1',
		     'PRSSO2','SO2','Depct2','Anxct2')
selVarsS1	<- c('PRSSO1','SO1','Depct1','Anxct1')
selVarsS2	<- c('PRSSO2','SO2','Depct2','Anxct2')
useVars	<- c('PRSSO1','SO1','Depct1','Anxct1',
		     'PRSSO2','SO2','Depct2','Anxct2','age1','sex1','ELA1','age2','sex2','ELA2')
useVars1	<- c('PRSSO1','SO1','Depct1','Anxct1','age1','sex1','ELA1')
useVars2	<- c('PRSSO2','SO2','Depct2','Anxct2','age2','sex2','ELA2')

mzData	<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , selVars)
dzData	<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , selVars)

mz12Data	<- TWINdata2[TWINdata2$zyg1==c(1,3) & (!is.na(TWINdata2$ELA1) & !is.na(TWINdata2$ELA2)), useVars]
dz12Data	<- TWINdata2[TWINdata2$zyg1==c(2,4,5) & (!is.na(TWINdata2$ELA1) & !is.na(TWINdata2$ELA2)), useVars]
sib1Data	<- TWINdata2[!is.na(TWINdata2$ELA1) & is.na(TWINdata2$ELA2), useVars1]
sib2Data	<- TWINdata2[is.na(TWINdata2$ELA1) & !is.na(TWINdata2$ELA2), useVars2]

psych::describe(mzData)
psych::describe(dzData)

psych::describe(mz12Data)
psych::describe(dz12Data)
psych::describe(sib1Data)
psych::describe(sib2Data)

mz12Data$SO1[is.na(mz12Data$age1)] <- NA
mz12Data$SO2[is.na(mz12Data$age2)] <- NA
dz12Data$SO1[is.na(dz12Data$age1)] <- NA
dz12Data$SO2[is.na(dz12Data$age2)] <- NA
sib1Data$SO1[is.na(sib1Data$age1)] <- NA
sib2Data$SO2[is.na(sib2Data$age2)] <- NA

mz12Data$SO1[is.na(mz12Data$sex1)] <- NA
mz12Data$SO2[is.na(mz12Data$sex2)] <- NA
dz12Data$SO1[is.na(dz12Data$sex1)] <- NA
dz12Data$SO2[is.na(dz12Data$sex2)] <- NA
sib1Data$SO1[is.na(sib1Data$sex1)] <- NA
sib2Data$SO2[is.na(sib2Data$sex2)] <- NA

mz12Data$SO1[is.na(mz12Data$ELA1)] <- NA
mz12Data$SO2[is.na(mz12Data$ELA2)] <- NA
dz12Data$SO1[is.na(dz12Data$ELA1)] <- NA
dz12Data$SO2[is.na(dz12Data$ELA2)] <- NA
sib1Data$SO1[is.na(sib1Data$ELA1)] <- NA
sib2Data$SO2[is.na(sib2Data$ELA2)] <- NA
mz12Data$Depct1[is.na(mz12Data$ELA1)] <- NA
mz12Data$Depct2[is.na(mz12Data$ELA2)] <- NA
dz12Data$Depct1[is.na(dz12Data$ELA1)] <- NA
dz12Data$Depct2[is.na(dz12Data$ELA2)] <- NA
sib1Data$Depct1[is.na(sib1Data$ELA1)] <- NA
sib2Data$Depct2[is.na(sib2Data$ELA2)] <- NA
mz12Data$Anxct1[is.na(mz12Data$ELA1)] <- NA
mz12Data$Anxct2[is.na(mz12Data$ELA2)] <- NA
dz12Data$Anxct1[is.na(dz12Data$ELA1)] <- NA
dz12Data$Anxct2[is.na(dz12Data$ELA2)] <- NA
sib1Data$Anxct1[is.na(sib1Data$ELA1)] <- NA
sib2Data$Anxct2[is.na(sib2Data$ELA2)] <- NA

mz12Data$age1[is.na(mz12Data$age1)] <- 999
mz12Data$age2[is.na(mz12Data$age2)] <- 999
dz12Data$age1[is.na(dz12Data$age1)] <- 999
dz12Data$age2[is.na(dz12Data$age2)] <- 999
sib1Data$age1[is.na(sib1Data$age1)] <- 999
sib2Data$age2[is.na(sib2Data$age2)] <- 999

mz12Data$sex1[is.na(mz12Data$sex1)] <- 999
mz12Data$sex2[is.na(mz12Data$sex2)] <- 999
dz12Data$sex1[is.na(dz12Data$sex1)] <- 999
dz12Data$sex2[is.na(dz12Data$sex2)] <- 999
sib1Data$sex1[is.na(sib1Data$sex1)] <- 999
sib2Data$sex2[is.na(sib2Data$sex2)] <- 999

mz12Data$ELA1[is.na(mz12Data$ELA1)] <- 999
mz12Data$ELA2[is.na(mz12Data$ELA2)] <- 999
dz12Data$ELA1[is.na(dz12Data$ELA1)] <- 999
dz12Data$ELA2[is.na(dz12Data$ELA2)] <- 999
sib1Data$ELA1[is.na(sib1Data$ELA1)] <- 999
sib2Data$ELA2[is.na(sib2Data$ELA2)] <- 999

psych::describe(mz12Data)
psych::describe(dz12Data)
psych::describe(sib1Data)
psych::describe(sib2Data)

# ******************************************************************************************************************************************
# (5) van der Sluis correction
# ******************************************************************************************************************************************
# I control for the cotwin moderation variable (ELA)  (van der Sluis et al, 2012) via residualisation (for SOI, SOP, Dep, And and RSB). 
# This controls for any confounding due to main effects of the co-twin/sibling's moderator variable.
# I use this approach to reduce the complexity of the model
# *********************************************************************************************

mz12Data$Depct1	<- residuals(lm(mz12Data$Depct1 ~ mz12Data$ELA2, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Depct1	<- (mz12Data$Depct1 + 3.27)
psych::describe(mz12Data$Depct1)
hist(mz12Data$Depct1)

mz12Data$Depct2	<- residuals(lm(mz12Data$Depct2 ~ mz12Data$ELA1, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Depct2	<- (mz12Data$Depct2 + 3.32)
psych::describe(mz12Data$Depct2)
hist(mz12Data$Depct2)

dz12Data$Depct1	<- residuals(lm(dz12Data$Depct1 ~ dz12Data$ELA2, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Depct1	<- (dz12Data$Depct1 + 3.30)
psych::describe(dz12Data$Depct1)
hist(dz12Data$Depct1)

dz12Data$Depct2	<- residuals(lm(dz12Data$Depct2 ~ dz12Data$ELA1, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Depct2	<- (dz12Data$Depct2 + 3.35)
psych::describe(dz12Data$Depct2)
hist(dz12Data$Depct2)

mz12Data$Anxct1	<- residuals(lm(mz12Data$Anxct1 ~ mz12Data$ELA2, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Anxct1	<- (mz12Data$Anxct1 + 8.24)
psych::describe(mz12Data$Anxct1)
hist(mz12Data$Anxct1)

mz12Data$Anxct2	<- residuals(lm(mz12Data$Anxct2 ~ mz12Data$ELA1, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Anxct2	<- (mz12Data$Anxct2 + 8.22)
psych::describe(mz12Data$Anxct2)
hist(mz12Data$Anxct2)

dz12Data$Anxct1	<- residuals(lm(dz12Data$Anxct1 ~ dz12Data$ELA2, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Anxct1	<- (dz12Data$Anxct1 + 8.26)
psych::describe(dz12Data$Anxct1)
hist(dz12Data$Anxct1)

dz12Data$Anxct2	<- residuals(lm(dz12Data$Anxct2 ~ dz12Data$ELA1, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Anxct2	<- (dz12Data$Anxct2 + 8.28)
psych::describe(dz12Data$Anxct2)
hist(dz12Data$Anxct2)

psych::describe(mz12Data)
psych::describe(dz12Data)
psych::describe(sib1Data)
psych::describe(sib2Data)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- c(5,0,4.5,4.1))
(PatM		<- c(TRUE,F,TRUE,TRUE))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEs	<- c('es1','es2','es3','es3'))
(LabAs	<- c('as1','as2','as3','as3'))

PatSpe	<- c(F,F,TRUE,TRUE)
PatSpac	<- c(F,F,TRUE,TRUE)
StSpa		<- c(0,0,.5,.5)
StSpe		<- c(0,0,.5,.5)

# all 1st loadings fixed to 1
PatFl		<- c(F,F,F,F,			
		     F,F,F,F,
		     F,F,F,T)

StFl		<- c(1,0,0,0,
		     0,1,0,0,
		     0,0,1,.5)

LabFl		<- c('l1',NA,NA,NA,
	 	     NA,'l2',NA,NA,
	 	     NA,NA,'l3','l4')

PatPhC	<- c(F,T,T,
		     F,F,T,
		     F,F,F)

StPhC		<- c(0,.1,.01,
		     0,0,.01,
		     0,0,0)

LabPhC	<- c(NA,'c1on2','c1on3',
		     NA,NA,'c2on3',
		     NA,NA,NA)	 

StPhCM	<- c(0,0,0,
		     0,0,.01,
		     0,0,0)

LabPhCM	<- c(NA,'c1on2mod','c1on3mod',
		     NA,NA,'c2on3mod',	
		     NA,NA,NA)	 

PatPhCM	<- c(F,F,F,
		     F,F,T,
		     F,F,F)

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmean,Stmean), labels=c(mLabs,mLabs), name="expMean") 
MeansS	<-mxMatrix("Full", 1, nv, free=c(PatM), values=c(Stmean), labels=c(mLabs), name="expMeanS") 

# Threshold and covariates
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')
LabCovE 	<-c('BelaThSO','BelaThSO','BelaThSO','BelaThSO')

ThPat		<-c(T,T,T,T)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

obsela1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.ELA1"), name="ELA1")
obsela2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.ELA2"), name="ELA2")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.01, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.01, labels=LabCovS, name="BsexTH" )
betaE		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.1, labels=LabCovE, name="BelaTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1 + BelaTH%x%ELA1,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2 + BelaTH%x%ELA2), name="expThres")
Thres1	<-mxAlgebra( expression= Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1 + BelaTH%x%ELA1, name="expThres1")
Thres2	<-mxAlgebra( expression= Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2 + BelaTH%x%ELA2, name="expThres2")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )
PhCausMod	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhCM, values=StPhCM, labels=LabPhCM, name="PhCM" )

# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpa, labels=LabAs, name="as" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpe, labels=LabEs, name="es" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.5,.5,.5), labels=c("ac22","ac32","ac33"), name="ac" ) # Component paths for factors 2 and 3
PathsEcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(.5,0,.5), labels=c("ec22","ec32","ec33"), name="ec" )

PathsAcMod	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,T,T), values=c(.01,.01,.1), labels=c("acMod22","acMod32","acMod33"), name="acM" ) # Component paths for factors 2 and 3
PathsEcMod	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(.1,0,.1), labels=c("ecMod22","ecMod32","ecMod33"), name="ecM" )

PathsP11	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11", name="pc" ) # SD path for factor 3 (the PRS factor)
Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  #Padding

Ze21		<-mxMatrix(type="Zero",	nrow=2, ncol=1, free=F, name="Z21" )  #Padding
Ze12		<-mxMatrix(type="Zero",	nrow=1, ncol=2, free=F, name="Z12" )  #Padding

covAcsub11	<-mxAlgebra( expression= (ac + (acM%x%ELA1))%*% t(ac + (acM%x%ELA1)), name="Acsub11" )
covAcsub22	<-mxAlgebra( expression= (ac + (acM%x%ELA2))%*% t(ac + (acM%x%ELA2)), name="Acsub22" )
covAcsub12	<-mxAlgebra( expression= (ac + (acM%x%ELA1))%*% t(ac + (acM%x%ELA2)), name="Acsub12" )
covAcsub21	<-mxAlgebra( expression= (ac + (acM%x%ELA2))%*% t(ac + (acM%x%ELA1)), name="Acsub21" )

covEcsub11	<-mxAlgebra( expression= (ec + (ecM%x%ELA1))%*% t(ec + (ecM%x%ELA1)), name="Ecsub11" )
covEcsub22	<-mxAlgebra( expression= (ec + (ecM%x%ELA2))%*% t(ec + (ecM%x%ELA2)), name="Ecsub22" )
covEcsub12	<-mxAlgebra( expression= (ec + (ecM%x%ELA1))%*% t(ec + (ecM%x%ELA2)), name="Ecsub12" )
covEcsub21	<-mxAlgebra( expression= (ec + (ecM%x%ELA2))%*% t(ec + (ecM%x%ELA1)), name="Ecsub21" )

covPcsub11	<-mxAlgebra( expression= Acsub11+Ecsub11, name="Vcsub11" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)
covPcsub22	<-mxAlgebra( expression= Acsub22+Ecsub22, name="Vcsub22" ) 

covPc11	<-mxAlgebra( expression= pc %*% t(pc), name="Pc11" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components

covPC11	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Vcsub11)), name="Vc11") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPc22	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Vcsub22)), name="Vc22") 

covPcMz12	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Acsub12)), name="Vcmz12") #I specify the MZ between-twin covariance - excluding E parameters
covPcMz21	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Acsub21)), name="Vcmz21") 

covPcDz12	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z21), rbind(Z12, .5%x%Acsub12)), name="Vcdz12") #I specify the DZ between-twin covariance - specifying half of A and excluding E
covPcDz21	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z21), rbind(Z12, .5%x%Acsub21)), name="Vcdz21")
 
# Generate Covariance of Latent factor model Including Causal Paths between factors
Id3		<-mxMatrix(type="Iden",	nrow=3, ncol=3, free=F, name="I3" )
covFVc1	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%ELA1))) %&% Vc11, name ="FVc1")
covFVc2	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%ELA2))) %&% Vc22, name ="FVc2")

covFcMz12	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%ELA1))) %*% Vcmz12 %*% solve(I3-(PhC+(PhCM%x%ELA2))), name ="Fcmz12")
covFcMz21	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%ELA2))) %*% Vcmz21 %*% solve(I3-(PhC+(PhCM%x%ELA1))), name ="Fcmz21")

covFcDz12	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%ELA1))) %*% Vcdz12 %*% solve(I3-(PhC+(PhCM%x%ELA2))), name ="Fcdz12")
covFcDz21	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%ELA2))) %*% Vcdz21 %*% solve(I3-(PhC+(PhCM%x%ELA1))), name ="Fcdz21")

# Constraint on total variance of Ordinal variable (A+C+E=1)
varL1		<- mxConstraint( expression=FVc1[2,2]==1, name="L1" )
varL2		<- mxConstraint( expression=FVc2[2,2]==1, name="L2" )

FcovMZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc1, Fcmz12), cbind(Fcmz21, FVc2))) , name="expFCovMZ" )
FcovDZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc1, Fcdz12), cbind(Fcdz21, FVc2))) , name="expFCovDZ" )
FcovSIB1	<-mxAlgebra( expression= (FactL  %&% FVc1 ), name="expFCovSIB1" )
FcovSIB2	<-mxAlgebra( expression= (FactL  %&% FVc2 ), name="expFCovSIB2" )

SpcovMZ	<-mxAlgebra( expression= rbind (cbind(Vs, As), cbind(As, Vs)) , name="expSpCovMZ" )
SpcovDZ	<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As), cbind(.5%x%As, Vs)) , name="expSpCovDZ" )

TOTcovMZ	<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )
TOTcovSIB1	<-mxAlgebra( expression= expFCovSIB1 + Vs , name="TOTexpCovSIB1" )
TOTcovSIB2	<-mxAlgebra( expression= expFCovSIB2 + Vs , name="TOTexpCovSIB2" )

# *******************************************************************************************************
# Calculator

# Standardize the Total var/covariances matrices of the observed variables
Id8		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I8" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= sqrt(As/( (FactL %&% FVc1) +Vs)), name="stAs" )
stcovEs	<-mxAlgebra( expression= sqrt(Es/( (FactL %&% FVc1) +Vs)), name="stEs" )

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% FVc1 / TOTexpCovMZ[1:4,1:4])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
ctdataMZ	<- mxData( observed=mz12Data, type="raw" )
ctdataDZ	<- mxData( observed=dz12Data, type="raw" )
ctdataSIB1	<- mxData( observed=sib1Data, type="raw" )
ctdataSIB2	<- mxData( observed=sib2Data, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
objSIB1	<- mxExpectationNormal( covariance="TOTexpCovSIB1", means="expMeanS", dimnames=selVarsS1, thresholds="expThres1", threshnames=c("SO1"))
objSIB2	<- mxExpectationNormal( covariance="TOTexpCovSIB2", means="expMeanS", dimnames=selVarsS2, thresholds="expThres2", threshnames=c("SO2"))

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1		<-list(Means,Load,LoadTw,PhCaus,PhCausMod,PathsAs,PathsEs,covAs,covEs,covPs,Id3,Id8,Id2)
pars2		<-list(PathsAcsub,PathsEcsub,PathsAcMod,PathsEcMod,PathsP11,Ze21,Ze12)
pars3		<-list(Ze11,covAcsub11,covAcsub22,covAcsub12,covAcsub21,covEcsub11,covEcsub22,covEcsub12,covEcsub21,
				covPcsub11,covPcsub22,covPc11,covPC11,covPc22,covFVc1,covFVc2)
pars4		<-list(obsage1, obsage2, obssex1, obssex2, obsela1, obsela2, betaA, betaS, betaE, Thr, inc, Thres)
parsst	<-list(stcovAs, stcovEs)
parsmz	<-list(covPcMz12,covPcMz21,covFcMz12,covFcMz21)
parsdz	<-list(covPcDz12,covPcDz21,covFcDz12,covFcDz21)

pars1a	<-list(MeansS,Load,PhCaus,PhCausMod,PathsAs,PathsEs,covAs,covEs,covPs,Id3)
pars1b	<-list(PathsAcsub,PathsEcsub,PathsAcMod,PathsEcMod,PathsP11,Ze21,Ze12)
pars1c	<-list(covAcsub11,covEcsub11)
pars1d	<-list(covPcsub11,covPc11,covPC11,covFVc1)
pars1e	<-list(obsage1, obssex1, obsela1, betaA, betaS, betaE, Thr, inc, Thres1)

pars2a	<-list(MeansS,Load,PhCaus,PhCausMod,PathsAs,PathsEs,covAs,covEs,covPs,Id3)
pars2b	<-list(PathsAcsub,PathsEcsub,PathsAcMod,PathsEcMod,PathsP11,Ze21,Ze12)
pars2c	<-list(covAcsub22,covEcsub22)
pars2d	<-list(covPcsub22,covPc11,covPc22,covFVc2)
pars2e	<-list(obsage2, obssex2, obsela2, betaA, betaS, betaE, Thr, inc, Thres2)

modelMZ	<-mxModel(parsmz, pars1, pars2, pars3, pars4, FcovMZ, SpcovMZ, TOTcovMZ, ctdataMZ, objMZ, Rfactmz, parsst, fitFunction, StFL, varL1, varL2, name="MZ" )
modelDZ	<-mxModel(parsdz, pars1, pars2, pars3, pars4, FcovDZ, SpcovDZ, TOTcovDZ, ctdataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )
modelSIB1	<-mxModel(pars1a, pars1b, pars1c, pars1d, pars1e, FcovSIB1, TOTcovSIB1, ctdataSIB1, objSIB1, fitFunction, name="SIB1" )
modelSIB2	<-mxModel(pars2a, pars2b, pars2c, pars2d, pars2e, FcovSIB2, TOTcovSIB2, ctdataSIB2, objSIB2, fitFunction, name="SIB2" )

minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective + SIB1.objective + SIB2.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
#cistFL	<-mxCI (c ('MZ.StandFact','MZ.PhCM','MZ.PhC','MZ.acM','MZ.ecM'))
#cistVs	<-mxCI (c ('MZ.stAs[3,3]','MZ.stEs[3,3]') ) 	# standardized var comp from specific Factors
cistFL	<-mxCI (c ('MZ.PhCM'))
MRDoCMod1Model	<-mxModel("MRDoCMod1", modelMZ, modelDZ, modelSIB1, modelSIB2, minus2ll, obj, cistFL) 

# --------------------------------------------------------------------------------------------------------------------------------
# 6a RUN ACEMs Factor Model with phenotypic causal mediation paths by Zygosity

MRDoCMod1Fit	<-mxTryHardOrdinal(MRDoCMod1Model, bestInitsOutput=T, intervals=T)
(MRDoCMod1Summ	<-summary(MRDoCMod1Fit, verbose=F))

mxEval(MZ.FactcorMZ, MRDoCMod1Fit)
mxEval(DZ.FactcorDZ, MRDoCMod1Fit)

mxEval(MZ.FVc1, MRDoCMod1Fit)
mxEval(MZ.FVc2, MRDoCMod1Fit)

mxEval(MZ.StandFact, MRDoCMod1Fit)

mxEval(MZ.Vs, MRDoCMod1Fit)

mxEval(MZ.Acsub11, MRDoCMod1Fit)
mxEval(MZ.Ecsub11, MRDoCMod1Fit)
mxEval(MZ.acM, MRDoCMod1Fit)
mxEval(MZ.ecM, MRDoCMod1Fit)
mxEval(MZ.PhCM, MRDoCMod1Fit)
mxEval(MZ.PhC, MRDoCMod1Fit)

mxEval(MZ.Vcsub11, MRDoCMod1Fit)
mxEval(MZ.Vc11, MRDoCMod1Fit)
mxEval(MZ.Vc22, MRDoCMod1Fit)
mxEval(MZ.expFCovMZ123, MRDoCMod1Fit)
mxEval(DZ.expFCovDZ123, MRDoCMod1Fit)

mxEval(MZ.stAs, MRDoCMod1Fit)
mxEval(MZ.stEs, MRDoCMod1Fit)


MRDoCNoMod1cModel	<- mxModel(MRDoCMod1Fit, name="MRDoCNoMod1c")
MRDoCNoMod1cModel	<- omxSetParameters(MRDoCNoMod1cModel, labels=c('c2on3mod'), free=FALSE, values=0)
MRDoCNoMod1cFit	<- mxRun(MRDoCNoMod1cModel, intervals=F)
(MRDoCNoMod1cSum	<- summary(MRDoCNoMod1cFit, verbose=F))

mxCompare(MRDoCMod1Fit, MRDoCNoMod1cFit)


# 
#****************************************************************************************************************************
# __(VIb)_____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL for SO>PD incorporating moderation by CGN4
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Phenotypic Factors: F1>F2>F3 & F1>F3;
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
#_____________________________________________________________________________________________________________________________

nv		<- 4				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 3				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nvo 		<- 1     			# number of ordinal variables per twin
nvc 		<- nv-nvo  			# number of continuous variables per twin
poso 		<- nvo 			# position where ordinal variables start
nth		<- 4				# number of max thresholds
ninc 		<- nth-1 			# number of max increments
ncovariates <- 3 				# number of covariates
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups	<- c("mz", "dz")
Vars		<- c('PRSSO','SO','Depcg4','Anxcg4')
selVars	<- c('PRSSO1','SO1','Depcg41','Anxcg41',
		     'PRSSO2','SO2','Depcg42','Anxcg42')
selVarsS1	<- c('PRSSO1','SO1','Depcg41','Anxcg41')
selVarsS2	<- c('PRSSO2','SO2','Depcg42','Anxcg42')
useVars	<- c('PRSSO1','SO1','Depcg41','Anxcg41',
		     'PRSSO2','SO2','Depcg42','Anxcg42','age1','sex1','CGN41','age2','sex2','CGN42')
useVars1	<- c('PRSSO1','SO1','Depcg41','Anxcg41','age1','sex1','CGN41')
useVars2	<- c('PRSSO2','SO2','Depcg42','Anxcg42','age2','sex2','CGN42')

mzData	<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , selVars)
dzData	<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , selVars)

mz12Data	<- TWINdata2[TWINdata2$zyg1==c(1,3) & (!is.na(TWINdata2$CGN41) & !is.na(TWINdata2$CGN42)), useVars]
dz12Data	<- TWINdata2[TWINdata2$zyg1==c(2,4,5) & (!is.na(TWINdata2$CGN41) & !is.na(TWINdata2$CGN42)), useVars]
sib1Data	<- TWINdata2[!is.na(TWINdata2$CGN41) & is.na(TWINdata2$CGN42), useVars1]
sib2Data	<- TWINdata2[is.na(TWINdata2$CGN41) & !is.na(TWINdata2$CGN42), useVars2]

psych::describe(mzData)
psych::describe(dzData)

psych::describe(mz12Data)
psych::describe(dz12Data)
psych::describe(sib1Data)
psych::describe(sib2Data)

mz12Data$SO1[is.na(mz12Data$age1)] <- NA
mz12Data$SO2[is.na(mz12Data$age2)] <- NA
dz12Data$SO1[is.na(dz12Data$age1)] <- NA
dz12Data$SO2[is.na(dz12Data$age2)] <- NA
sib1Data$SO1[is.na(sib1Data$age1)] <- NA
sib2Data$SO2[is.na(sib2Data$age2)] <- NA

mz12Data$SO1[is.na(mz12Data$sex1)] <- NA
mz12Data$SO2[is.na(mz12Data$sex2)] <- NA
dz12Data$SO1[is.na(dz12Data$sex1)] <- NA
dz12Data$SO2[is.na(dz12Data$sex2)] <- NA
sib1Data$SO1[is.na(sib1Data$sex1)] <- NA
sib2Data$SO2[is.na(sib2Data$sex2)] <- NA

mz12Data$SO1[is.na(mz12Data$CGN41)] <- NA
mz12Data$SO2[is.na(mz12Data$CGN42)] <- NA
dz12Data$SO1[is.na(dz12Data$CGN41)] <- NA
dz12Data$SO2[is.na(dz12Data$CGN42)] <- NA
sib1Data$SO1[is.na(sib1Data$CGN41)] <- NA
sib2Data$SO2[is.na(sib2Data$CGN42)] <- NA
mz12Data$Depcg41[is.na(mz12Data$CGN41)] <- NA
mz12Data$Depcg42[is.na(mz12Data$CGN42)] <- NA
dz12Data$Depcg41[is.na(dz12Data$CGN41)] <- NA
dz12Data$Depcg42[is.na(dz12Data$CGN42)] <- NA
sib1Data$Depcg41[is.na(sib1Data$CGN41)] <- NA
sib2Data$Depcg42[is.na(sib2Data$CGN42)] <- NA
mz12Data$Anxcg41[is.na(mz12Data$CGN41)] <- NA
mz12Data$Anxcg42[is.na(mz12Data$CGN42)] <- NA
dz12Data$Anxcg41[is.na(dz12Data$CGN41)] <- NA
dz12Data$Anxcg42[is.na(dz12Data$CGN42)] <- NA
sib1Data$Anxcg41[is.na(sib1Data$CGN41)] <- NA
sib2Data$Anxcg42[is.na(sib2Data$CGN42)] <- NA

mz12Data$age1[is.na(mz12Data$age1)] <- 999
mz12Data$age2[is.na(mz12Data$age2)] <- 999
dz12Data$age1[is.na(dz12Data$age1)] <- 999
dz12Data$age2[is.na(dz12Data$age2)] <- 999
sib1Data$age1[is.na(sib1Data$age1)] <- 999
sib2Data$age2[is.na(sib2Data$age2)] <- 999

mz12Data$sex1[is.na(mz12Data$sex1)] <- 999
mz12Data$sex2[is.na(mz12Data$sex2)] <- 999
dz12Data$sex1[is.na(dz12Data$sex1)] <- 999
dz12Data$sex2[is.na(dz12Data$sex2)] <- 999
sib1Data$sex1[is.na(sib1Data$sex1)] <- 999
sib2Data$sex2[is.na(sib2Data$sex2)] <- 999

mz12Data$CGN41[is.na(mz12Data$CGN41)] <- 999
mz12Data$CGN42[is.na(mz12Data$CGN42)] <- 999
dz12Data$CGN41[is.na(dz12Data$CGN41)] <- 999
dz12Data$CGN42[is.na(dz12Data$CGN42)] <- 999
sib1Data$CGN41[is.na(sib1Data$CGN41)] <- 999
sib2Data$CGN42[is.na(sib2Data$CGN42)] <- 999

psych::describe(mz12Data)
psych::describe(dz12Data)
psych::describe(sib1Data)
psych::describe(sib2Data)

# ******************************************************************************************************************************************
# (5) van der Sluis correction
# ******************************************************************************************************************************************
# I control for the cotwin and sibling moderation variable (CTQ) via residualisation (for SOI, SOP, Dep, And and RSB) (van der Sluis et al, 2012). 
# This controls for any confounding due to main effects of the co-twin/sibling's moderator variable.
# I use this approach to reduce the complexity of the model
# *********************************************************************************************

mz12Data$Depcg41	<- residuals(lm(mz12Data$Depcg41 ~ mz12Data$CGN42, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Depcg41	<- (mz12Data$Depcg41 + 6.27)
psych::describe(mz12Data$Depcg41)
hist(mz12Data$Depcg41)

mz12Data$Depcg42	<- residuals(lm(mz12Data$Depcg42 ~ mz12Data$CGN41, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Depcg42	<- (mz12Data$Depcg42 + 6.32)
psych::describe(mz12Data$Depcg42)
hist(mz12Data$Depcg42)

dz12Data$Depcg41	<- residuals(lm(dz12Data$Depcg41 ~ dz12Data$CGN42, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Depcg41	<- (dz12Data$Depcg41 + 6.30)
psych::describe(dz12Data$Depcg41)
hist(dz12Data$Depcg41)

dz12Data$Depcg42	<- residuals(lm(dz12Data$Depcg42 ~ dz12Data$CGN41, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Depcg42	<- (dz12Data$Depcg42 + 6.35)
psych::describe(dz12Data$Depcg42)
hist(dz12Data$Depcg42)

mz12Data$Anxcg41	<- residuals(lm(mz12Data$Anxcg41 ~ mz12Data$CGN42, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Anxcg41	<- (mz12Data$Anxcg41 + 4.24)
psych::describe(mz12Data$Anxcg41)
hist(mz12Data$Anxcg41)

mz12Data$Anxcg42	<- residuals(lm(mz12Data$Anxcg42 ~ mz12Data$CGN41, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Anxcg42	<- (mz12Data$Anxcg42 + 4.22)
psych::describe(mz12Data$Anxcg42)
hist(mz12Data$Anxcg42)

dz12Data$Anxcg41	<- residuals(lm(dz12Data$Anxcg41 ~ dz12Data$CGN42, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Anxcg41	<- (dz12Data$Anxcg41 + 4.26)
psych::describe(dz12Data$Anxcg41)
hist(dz12Data$Anxcg41)

dz12Data$Anxcg42	<- residuals(lm(dz12Data$Anxcg42 ~ dz12Data$CGN41, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Anxcg42	<- (dz12Data$Anxcg42 + 4.28)
psych::describe(dz12Data$Anxcg42)
hist(dz12Data$Anxcg42)

psych::describe(mz12Data)
psych::describe(dz12Data)
psych::describe(sib1Data)
psych::describe(sib2Data)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- c(5,0,6.5,4.5))
(PatM		<- c(TRUE,F,TRUE,TRUE))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEs	<- c('es1','es2','es3','es3'))
(LabAs	<- c('as1','as2','as3','as3'))

PatSpe	<- c(F,F,TRUE,TRUE)
PatSpac	<- c(F,F,TRUE,TRUE)
StSpa		<- c(0,0,.5,.5)
StSpe		<- c(0,0,.5,.5)

# all 1st loadings fixed to 1
PatFl		<- c(F,F,F,F,			
		     F,F,F,F,
		     F,F,F,T)

StFl		<- c(1,0,0,0,
		     0,1,0,0,
		     0,0,1,.5)

LabFl		<- c('l1',NA,NA,NA,
	 	     NA,'l2',NA,NA,
	 	     NA,NA,'l3','l4')

PatPhC	<- c(F,T,T,
		     F,F,T,
		     F,F,F)

StPhC		<- c(0,.05,.01,
		     0,0,.01,
		     0,0,0)

LabPhC	<- c(NA,'c1on2','c1on3',
		     NA,NA,'c2on3',
		     NA,NA,NA)	 

StPhCM	<- c(0,0,0,
		     0,0,.01,
		     0,0,0)

LabPhCM	<- c(NA,'c1on2mod','c1on3mod',
		     NA,NA,'c2on3mod',	
		     NA,NA,NA)	 

PatPhCM	<- c(F,F,F,
		     F,F,T,
		     F,F,F)

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmean,Stmean), labels=c(mLabs,mLabs), name="expMean") 
MeansS	<-mxMatrix("Full", 1, nv, free=c(PatM), values=c(Stmean), labels=c(mLabs), name="expMeanS") 

# Threshold and covariates
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')
LabCovC 	<-c('BcgnThSO','BcgnThSO','BcgnThSO','BcgnThSO')

ThPat		<-c(T,T,T,T)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

obscgn1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.CGN41"), name="CGN41")
obscgn2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.CGN42"), name="CGN42")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.01, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.01, labels=LabCovS, name="BsexTH" )
betaC		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.1, labels=LabCovC, name="BcgnTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1 + BcgnTH%x%CGN41,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2 + BcgnTH%x%CGN42), name="expThres")
Thres1	<-mxAlgebra( expression= Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1 + BcgnTH%x%CGN41, name="expThres1")
Thres2	<-mxAlgebra( expression= Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2 + BcgnTH%x%CGN42, name="expThres2")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )
PhCausMod	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhCM, values=StPhCM, labels=LabPhCM, name="PhCM" )

# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpa, labels=LabAs, name="as" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpe, labels=LabEs, name="es" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.5,.5,.5), labels=c("ac22","ac32","ac33"), name="ac" ) # Component paths for factors 2 and 3
PathsEcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(.5,0,.5), labels=c("ec22","ec32","ec33"), name="ec" )

PathsAcMod	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,T,T), values=c(.01,.01,.1), labels=c("acMod22","acMod32","acMod33"), name="acM" ) # Component paths for factors 2 and 3
PathsEcMod	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(.1,0,.1), labels=c("ecMod22","ecMod32","ecMod33"), name="ecM" )

PathsP11	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11", name="pc" ) # SD path for factor 3 (the PRS factor)
Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  #Padding

Ze21		<-mxMatrix(type="Zero",	nrow=2, ncol=1, free=F, name="Z21" )  #Padding
Ze12		<-mxMatrix(type="Zero",	nrow=1, ncol=2, free=F, name="Z12" )  #Padding

covAcsub11	<-mxAlgebra( expression= (ac + (acM%x%CGN41))%*% t(ac + (acM%x%CGN41)), name="Acsub11" )
covAcsub22	<-mxAlgebra( expression= (ac + (acM%x%CGN42))%*% t(ac + (acM%x%CGN42)), name="Acsub22" )
covAcsub12	<-mxAlgebra( expression= (ac + (acM%x%CGN41))%*% t(ac + (acM%x%CGN42)), name="Acsub12" )
covAcsub21	<-mxAlgebra( expression= (ac + (acM%x%CGN42))%*% t(ac + (acM%x%CGN41)), name="Acsub21" )

covEcsub11	<-mxAlgebra( expression= (ec + (ecM%x%CGN41))%*% t(ec + (ecM%x%CGN41)), name="Ecsub11" )
covEcsub22	<-mxAlgebra( expression= (ec + (ecM%x%CGN42))%*% t(ec + (ecM%x%CGN42)), name="Ecsub22" )
covEcsub12	<-mxAlgebra( expression= (ec + (ecM%x%CGN41))%*% t(ec + (ecM%x%CGN42)), name="Ecsub12" )
covEcsub21	<-mxAlgebra( expression= (ec + (ecM%x%CGN42))%*% t(ec + (ecM%x%CGN41)), name="Ecsub21" )

covPcsub11	<-mxAlgebra( expression= Acsub11+Ecsub11, name="Vcsub11" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)
covPcsub22	<-mxAlgebra( expression= Acsub22+Ecsub22, name="Vcsub22" ) 

covPc11	<-mxAlgebra( expression= pc %*% t(pc), name="Pc11" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components

covPC11	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Vcsub11)), name="Vc11") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPc22	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Vcsub22)), name="Vc22") 

covPcMz12	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Acsub12)), name="Vcmz12") #I specify the MZ between-twin covariance - excluding E parameters
covPcMz21	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Acsub21)), name="Vcmz21") 

covPcDz12	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z21), rbind(Z12, .5%x%Acsub12)), name="Vcdz12") #I specify the DZ between-twin covariance - specifying half of A and excluding E
covPcDz21	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z21), rbind(Z12, .5%x%Acsub21)), name="Vcdz21")
 
# Generate Covariance of Latent factor model Including Causal Paths between factors
Id3		<-mxMatrix(type="Iden",	nrow=3, ncol=3, free=F, name="I3" )
covFVc1	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%CGN41))) %&% Vc11, name ="FVc1")
covFVc2	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%CGN42))) %&% Vc22, name ="FVc2")

covFcMz12	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%CGN41))) %*% Vcmz12 %*% solve(I3-(PhC+(PhCM%x%CGN42))), name ="Fcmz12")
covFcMz21	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%CGN42))) %*% Vcmz21 %*% solve(I3-(PhC+(PhCM%x%CGN41))), name ="Fcmz21")

covFcDz12	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%CGN41))) %*% Vcdz12 %*% solve(I3-(PhC+(PhCM%x%CGN42))), name ="Fcdz12")
covFcDz21	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%CGN42))) %*% Vcdz21 %*% solve(I3-(PhC+(PhCM%x%CGN41))), name ="Fcdz21")

# Constraint on total variance of Ordinal variable (A+C+E=1)
varL1		<- mxConstraint( expression=FVc1[2,2]==1, name="L1" )
varL2		<- mxConstraint( expression=FVc2[2,2]==1, name="L2" )

FcovMZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc1, Fcmz12), cbind(Fcmz21, FVc2))) , name="expFCovMZ" )
FcovDZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc1, Fcdz12), cbind(Fcdz21, FVc2))) , name="expFCovDZ" )
FcovSIB1	<-mxAlgebra( expression= (FactL  %&% FVc1 ), name="expFCovSIB1" )
FcovSIB2	<-mxAlgebra( expression= (FactL  %&% FVc2 ), name="expFCovSIB2" )

SpcovMZ	<-mxAlgebra( expression= rbind (cbind(Vs, As), cbind(As, Vs)) , name="expSpCovMZ" )
SpcovDZ	<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As), cbind(.5%x%As, Vs)) , name="expSpCovDZ" )

TOTcovMZ	<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )
TOTcovSIB1	<-mxAlgebra( expression= expFCovSIB1 + Vs , name="TOTexpCovSIB1" )
TOTcovSIB2	<-mxAlgebra( expression= expFCovSIB2 + Vs , name="TOTexpCovSIB2" )

# *******************************************************************************************************
# Calculator

# Standardize the Total var/covariances matrices of the observed variables
Id8		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I8" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= sqrt(As/( (FactL %&% FVc1) +Vs)), name="stAs" )
stcovEs	<-mxAlgebra( expression= sqrt(Es/( (FactL %&% FVc1) +Vs)), name="stEs" )

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% FVc1 / TOTexpCovMZ[1:4,1:4])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
ctdataMZ	<- mxData( observed=mz12Data, type="raw" )
ctdataDZ	<- mxData( observed=dz12Data, type="raw" )
#ctdataSIB1	<- mxData( observed=sib1Data, type="raw" )
#ctdataSIB2	<- mxData( observed=sib2Data, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
#objSIB1	<- mxExpectationNormal( covariance="TOTexpCovSIB1", means="expMeanS", dimnames=selVarsS1, thresholds="expThres1", threshnames=c("SO1"))
#objSIB2	<- mxExpectationNormal( covariance="TOTexpCovSIB2", means="expMeanS", dimnames=selVarsS2, thresholds="expThres2", threshnames=c("SO2"))

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1		<-list(Means,Load,LoadTw,PhCaus,PhCausMod,PathsAs,PathsEs,covAs,covEs,covPs,Id3,Id8,Id2)
pars2		<-list(PathsAcsub,PathsEcsub,PathsAcMod,PathsEcMod,PathsP11,Ze21,Ze12)
pars3		<-list(Ze11,covAcsub11,covAcsub22,covAcsub12,covAcsub21,covEcsub11,covEcsub22,covEcsub12,covEcsub21,
				covPcsub11,covPcsub22,covPc11,covPC11,covPc22,covFVc1,covFVc2)
pars4		<-list(obsage1, obsage2, obssex1, obssex2, obscgn1, obscgn2, betaA, betaS, betaC, Thr, inc, Thres)
parsst	<-list(stcovAs, stcovEs)
parsmz	<-list(covPcMz12,covPcMz21,covFcMz12,covFcMz21)
parsdz	<-list(covPcDz12,covPcDz21,covFcDz12,covFcDz21)

pars1a	<-list(MeansS,Load,PhCaus,PhCausMod,PathsAs,PathsEs,covAs,covEs,covPs,Id3)
pars1b	<-list(PathsAcsub,PathsEcsub,PathsAcMod,PathsEcMod,PathsP11,Ze21,Ze12)
pars1c	<-list(covAcsub11,covEcsub11)
pars1d	<-list(covPcsub11,covPc11,covPC11,covFVc1)
pars1e	<-list(obsage1, obssex1, obscgn1, betaA, betaS, betaC, Thr, inc, Thres1)

pars2a	<-list(MeansS,Load,PhCaus,PhCausMod,PathsAs,PathsEs,covAs,covEs,covPs,Id3)
pars2b	<-list(PathsAcsub,PathsEcsub,PathsAcMod,PathsEcMod,PathsP11,Ze21,Ze12)
pars2c	<-list(covAcsub22,covEcsub22)
pars2d	<-list(covPcsub22,covPc11,covPc22,covFVc2)
pars2e	<-list(obsage2, obssex2, obscgn2, betaA, betaS, betaC, Thr, inc, Thres2)

modelMZ	<-mxModel(parsmz, pars1, pars2, pars3, pars4, FcovMZ, SpcovMZ, TOTcovMZ, ctdataMZ, objMZ, Rfactmz, parsst, fitFunction, StFL, varL1, varL2, name="MZ" )
modelDZ	<-mxModel(parsdz, pars1, pars2, pars3, pars4, FcovDZ, SpcovDZ, TOTcovDZ, ctdataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )
#modelSIB1	<-mxModel(pars1a, pars1b, pars1c, pars1d, pars1e, FcovSIB1, TOTcovSIB1, ctdataSIB1, objSIB1, fitFunction, name="SIB1" )
#modelSIB2	<-mxModel(pars2a, pars2b, pars2c, pars2d, pars2e, FcovSIB2, TOTcovSIB2, ctdataSIB2, objSIB2, fitFunction, name="SIB2" )

#minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective + SIB1.objective + SIB2.objective, name="m2LL" )
minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
#cistFL	<-mxCI (c ('MZ.StandFact','MZ.PhCM','MZ.PhC','MZ.acM','MZ.ecM'))
#cistVs	<-mxCI (c ('MZ.stAs[3,3]','MZ.stEs[3,3]') ) 	# standardized var comp from specific Factors
cistFL	<-mxCI (c ('MZ.PhCM'))
#MRDoCMod1Model	<-mxModel("MRDoCMod1", modelMZ, modelDZ, minus2ll, obj, cistFL) 
MRDoCMod1Model	<-mxModel("MRDoCMod1", modelMZ, modelDZ, minus2ll, obj, cistFL) 

# --------------------------------------------------------------------------------------------------------------------------------
# 6a RUN ACEMs Factor Model with phenotypic causal mediation paths by Zygosity

MRDoCMod1Fit	<-mxTryHardOrdinal(MRDoCMod1Model, bestInitsOutput=T, intervals=T)
(MRDoCMod1Summ	<-summary(MRDoCMod1Fit, verbose=F))

mxEval(MZ.FactcorMZ, MRDoCMod1Fit)
mxEval(DZ.FactcorDZ, MRDoCMod1Fit)

mxEval(MZ.FVc1, MRDoCMod1Fit)
mxEval(MZ.FVc2, MRDoCMod1Fit)

mxEval(MZ.StandFact, MRDoCMod1Fit)

mxEval(MZ.Vs, MRDoCMod1Fit)

mxEval(MZ.Acsub11, MRDoCMod1Fit)
mxEval(MZ.Ecsub11, MRDoCMod1Fit)
mxEval(MZ.acM, MRDoCMod1Fit)
mxEval(MZ.ecM, MRDoCMod1Fit)
mxEval(MZ.PhCM, MRDoCMod1Fit)
mxEval(MZ.PhC, MRDoCMod1Fit)

mxEval(MZ.Vcsub11, MRDoCMod1Fit)
mxEval(MZ.Vc11, MRDoCMod1Fit)
mxEval(MZ.Vc22, MRDoCMod1Fit)
mxEval(MZ.expFCovMZ123, MRDoCMod1Fit)
mxEval(DZ.expFCovDZ123, MRDoCMod1Fit)

mxEval(MZ.stAs, MRDoCMod1Fit)
mxEval(MZ.stEs, MRDoCMod1Fit)


MRDoCNoMod1cModel	<- mxModel(MRDoCMod1Fit, name="MRDoCNoMod1c")
MRDoCNoMod1cModel	<- omxSetParameters(MRDoCNoMod1cModel, labels=c('c2on3mod'), free=FALSE, values=0)
MRDoCNoMod1cFit	<- mxRun(MRDoCNoMod1cModel, intervals=F)
(MRDoCNoMod1cSum	<- summary(MRDoCNoMod1cFit, verbose=F))

mxCompare(MRDoCMod1Fit, MRDoCNoMod1cFit)

# 
#****************************************************************************************************************************
# __(VIc) MH->SO_____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL for PD>SO incorporating moderation by ELA
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Phenotypic Factors: F1>F2>F3 & F1>F3;
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
#_____________________________________________________________________________________________________________________________

nv		<- 5				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 3				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nvo 		<- 1     			# number of ordinal variables per twin
nvc 		<- nv-nvo  			# number of continuous variables per twin
poso 		<- nvo 			# position where ordinal variables start
nth		<- 4				# number of max thresholds
ninc 		<- nth-1 			# number of max increments
ncovariates <- 3 				# number of covariates
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups	<- c("mz", "dz")
Vars		<- c('PRSAnx','PRSDep','Depct','Anxct','SO')
selVars	<- c('PRSAnx1','PRSDep1','Depct1','Anxct1','SO1',
		     'PRSAnx2','PRSDep2','Depct2','Anxct2','SO2')
selVarsS1	<- c('PRSAnx1','PRSDep1','Depct1','Anxct1','SO1')
selVarsS2	<- c('PRSAnx2','PRSDep2','Depct2','Anxct2','SO2')
useVars	<- c('PRSAnx1','PRSDep1','Depct1','Anxct1','SO1',
		     'PRSAnx2','PRSDep2','Depct2','Anxct2','SO2','age1','sex1','ELA1','age2','sex2','ELA2')
useVars1	<- c('PRSAnx1','PRSDep1','Depct1','Anxct1','SO1','age1','sex1','ELA1')
useVars2	<- c('PRSAnx2','PRSDep2','Depct2','Anxct2','SO2','age2','sex2','ELA2')

mzData	<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , selVars)
dzData	<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , selVars)

mz12Data	<- TWINdata2[TWINdata2$zyg1==c(1,3) & (!is.na(TWINdata2$ELA1) & !is.na(TWINdata2$ELA2)), useVars]
dz12Data	<- TWINdata2[TWINdata2$zyg1==c(2,4,5) & (!is.na(TWINdata2$ELA1) & !is.na(TWINdata2$ELA2)), useVars]
sib1Data	<- TWINdata2[!is.na(TWINdata2$ELA1) & is.na(TWINdata2$ELA2), useVars1]
sib2Data	<- TWINdata2[is.na(TWINdata2$ELA1) & !is.na(TWINdata2$ELA2), useVars2]

psych::describe(mzData)
psych::describe(dzData)

psych::describe(mz12Data)
psych::describe(dz12Data)
psych::describe(sib1Data)
psych::describe(sib2Data)

mz12Data$SO1[is.na(mz12Data$age1)] <- NA
mz12Data$SO2[is.na(mz12Data$age2)] <- NA
dz12Data$SO1[is.na(dz12Data$age1)] <- NA
dz12Data$SO2[is.na(dz12Data$age2)] <- NA
sib1Data$SO1[is.na(sib1Data$age1)] <- NA
sib2Data$SO2[is.na(sib2Data$age2)] <- NA

mz12Data$SO1[is.na(mz12Data$sex1)] <- NA
mz12Data$SO2[is.na(mz12Data$sex2)] <- NA
dz12Data$SO1[is.na(dz12Data$sex1)] <- NA
dz12Data$SO2[is.na(dz12Data$sex2)] <- NA
sib1Data$SO1[is.na(sib1Data$sex1)] <- NA
sib2Data$SO2[is.na(sib2Data$sex2)] <- NA

mz12Data$SO1[is.na(mz12Data$ELA1)] <- NA
mz12Data$SO2[is.na(mz12Data$ELA2)] <- NA
dz12Data$SO1[is.na(dz12Data$ELA1)] <- NA
dz12Data$SO2[is.na(dz12Data$ELA2)] <- NA
sib1Data$SO1[is.na(sib1Data$ELA1)] <- NA
sib2Data$SO2[is.na(sib2Data$ELA2)] <- NA
mz12Data$Depct1[is.na(mz12Data$ELA1)] <- NA
mz12Data$Depct2[is.na(mz12Data$ELA2)] <- NA
dz12Data$Depct1[is.na(dz12Data$ELA1)] <- NA
dz12Data$Depct2[is.na(dz12Data$ELA2)] <- NA
sib1Data$Depct1[is.na(sib1Data$ELA1)] <- NA
sib2Data$Depct2[is.na(sib2Data$ELA2)] <- NA
mz12Data$Anxct1[is.na(mz12Data$ELA1)] <- NA
mz12Data$Anxct2[is.na(mz12Data$ELA2)] <- NA
dz12Data$Anxct1[is.na(dz12Data$ELA1)] <- NA
dz12Data$Anxct2[is.na(dz12Data$ELA2)] <- NA
sib1Data$Anxct1[is.na(sib1Data$ELA1)] <- NA
sib2Data$Anxct2[is.na(sib2Data$ELA2)] <- NA

mz12Data$age1[is.na(mz12Data$age1)] <- 999
mz12Data$age2[is.na(mz12Data$age2)] <- 999
dz12Data$age1[is.na(dz12Data$age1)] <- 999
dz12Data$age2[is.na(dz12Data$age2)] <- 999
sib1Data$age1[is.na(sib1Data$age1)] <- 999
sib2Data$age2[is.na(sib2Data$age2)] <- 999

mz12Data$sex1[is.na(mz12Data$sex1)] <- 999
mz12Data$sex2[is.na(mz12Data$sex2)] <- 999
dz12Data$sex1[is.na(dz12Data$sex1)] <- 999
dz12Data$sex2[is.na(dz12Data$sex2)] <- 999
sib1Data$sex1[is.na(sib1Data$sex1)] <- 999
sib2Data$sex2[is.na(sib2Data$sex2)] <- 999

mz12Data$ELA1[is.na(mz12Data$ELA1)] <- 999
mz12Data$ELA2[is.na(mz12Data$ELA2)] <- 999
dz12Data$ELA1[is.na(dz12Data$ELA1)] <- 999
dz12Data$ELA2[is.na(dz12Data$ELA2)] <- 999
sib1Data$ELA1[is.na(sib1Data$ELA1)] <- 999
sib2Data$ELA2[is.na(sib2Data$ELA2)] <- 999

psych::describe(mz12Data)
psych::describe(dz12Data)
psych::describe(sib1Data)
psych::describe(sib2Data)

# ******************************************************************************************************************************************
# (5) van der Sluis correction
# ******************************************************************************************************************************************
# I control for the cotwin and sibling moderation variable (CTQ) via residualisation (for SOI, SOP, Dep, And and RSB) (van der Sluis et al, 2012). 
# This controls for any confounding due to main effects of the co-twin/sibling's moderator variable.
# I use this approach to reduce the complexity of the model
# *********************************************************************************************

mz12Data$Depct1	<- residuals(lm(mz12Data$Depct1 ~ mz12Data$ELA2, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Depct1	<- (mz12Data$Depct1 + 3.27)
psych::describe(mz12Data$Depct1)
hist(mz12Data$Depct1)

mz12Data$Depct2	<- residuals(lm(mz12Data$Depct2 ~ mz12Data$ELA1, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Depct2	<- (mz12Data$Depct2 + 3.32)
psych::describe(mz12Data$Depct2)
hist(mz12Data$Depct2)

dz12Data$Depct1	<- residuals(lm(dz12Data$Depct1 ~ dz12Data$ELA2, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Depct1	<- (dz12Data$Depct1 + 3.30)
psych::describe(dz12Data$Depct1)
hist(dz12Data$Depct1)

dz12Data$Depct2	<- residuals(lm(dz12Data$Depct2 ~ dz12Data$ELA1, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Depct2	<- (dz12Data$Depct2 + 3.35)
psych::describe(dz12Data$Depct2)
hist(dz12Data$Depct2)

mz12Data$Anxct1	<- residuals(lm(mz12Data$Anxct1 ~ mz12Data$ELA2, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Anxct1	<- (mz12Data$Anxct1 + 8.24)
psych::describe(mz12Data$Anxct1)
hist(mz12Data$Anxct1)

mz12Data$Anxct2	<- residuals(lm(mz12Data$Anxct2 ~ mz12Data$ELA1, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Anxct2	<- (mz12Data$Anxct2 + 8.22)
psych::describe(mz12Data$Anxct2)
hist(mz12Data$Anxct2)

dz12Data$Anxct1	<- residuals(lm(dz12Data$Anxct1 ~ dz12Data$ELA2, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Anxct1	<- (dz12Data$Anxct1 + 8.26)
psych::describe(dz12Data$Anxct1)
hist(dz12Data$Anxct1)

dz12Data$Anxct2	<- residuals(lm(dz12Data$Anxct2 ~ dz12Data$ELA1, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Anxct2	<- (dz12Data$Anxct2 + 8.28)
psych::describe(dz12Data$Anxct2)
hist(dz12Data$Anxct2)

psych::describe(mz12Data)
psych::describe(dz12Data)
psych::describe(sib1Data)
psych::describe(sib2Data)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- c(5,5,4.5,4.1,0))
(PatM		<- c(TRUE,T,TRUE,TRUE,F))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEs	<- c('es1','es1','es3','es3','es5'))
(LabAs	<- c('as1','as1','as3','as3','es5'))

PatSpe	<- c(T,T,TRUE,TRUE,F)
PatSpac	<- c(F,F,TRUE,TRUE,F)
StSpe		<- c(.5,.5,.5,.5,0)
StSpa		<- c(0,0,.5,.5,0)

# all 1st loadings fixed to 1
PatFl		<- c(F,T,F,F,F,			
		     F,F,F,T,F,
		     F,F,F,F,F)

StFl		<- c(1,.5,0,0,0,
		     0,0,1,.5,0,
		     0,0,0,0,1)

LabFl		<- c('l1','l2',NA,NA,NA,
	 	     NA,NA,'l3','l4',NA,
	 	     NA,NA,NA,NA,'l5')

PatPhC	<- c(F,T,T,
		     F,F,T,
		     F,F,F)

StPhC		<- c(0,.3,.01,
		     0,0,.01,
		     0,0,0)

LabPhC	<- c(NA,'c1on2','c1on3',
		     NA,NA,'c2on3',
		     NA,NA,NA)	 

PatPhCM	<- c(F,F,F,
		     F,F,T,
		     F,F,F)

StPhCM	<- c(0,0,0,
		     0,0,.01,
		     0,0,0)

LabPhCM	<- c(NA,'c1on2mod','c1on3mod',
		     NA,NA,'c2on3mod',	
		     NA,NA,NA)	 

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmean,Stmean), labels=c(mLabs,mLabs), name="expMean") 
MeansS	<-mxMatrix("Full", 1, nv, free=c(PatM), values=c(Stmean), labels=c(mLabs), name="expMeanS") 

# Threshold and covariates
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')
LabCovE 	<-c('BelaThSO','BelaThSO','BelaThSO','BelaThSO')

ThPat		<-c(T,T,T,T)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

obsela1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.ELA1"), name="ELA1")
obsela2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.ELA2"), name="ELA2")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.1, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.1, labels=LabCovS, name="BsexTH" )
betaE		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.1, labels=LabCovE, name="BelaTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1 + BelaTH%x%ELA1,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2 + BelaTH%x%ELA2), name="expThres")
Thres1	<-mxAlgebra( expression= Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1 + BelaTH%x%ELA1, name="expThres1")
Thres2	<-mxAlgebra( expression= Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2 + BelaTH%x%ELA2, name="expThres2")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )
PhCausMod	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhCM, values=StPhCM, labels=LabPhCM, name="PhCM" )

# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpa, labels=LabAs, name="as" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpe, labels=LabEs, name="es" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.8,.1,.5), labels=c("ac22","ac32","ac33"), name="ac" ) # Component paths for factors 2 and 3
PathsEcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(.8,0,.8), labels=c("ec22","ec32","ec33"), name="ec" )

PathsAcMod	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(.01,0,.1), labels=c("acMod22","acMod32","acMod33"), name="acM" ) # Component paths for factors 2 and 3
PathsEcMod	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(0.1,0,0.1), labels=c("ecMod22","ecMod32","ecMod33"), name="ecM" )

PathsP11	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11", name="pc" ) # SD path for factor 3 (the PRS factor)
Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  #Padding

Ze21		<-mxMatrix(type="Zero",	nrow=2, ncol=1, free=F, name="Z21" )  #Padding
Ze12		<-mxMatrix(type="Zero",	nrow=1, ncol=2, free=F, name="Z12" )  #Padding

covAcsub11	<-mxAlgebra( expression= (ac + (acM%x%ELA1))%*% t(ac + (acM%x%ELA1)), name="Acsub11" )
covAcsub22	<-mxAlgebra( expression= (ac + (acM%x%ELA2))%*% t(ac + (acM%x%ELA2)), name="Acsub22" )
covAcsub12	<-mxAlgebra( expression= (ac + (acM%x%ELA1))%*% t(ac + (acM%x%ELA2)), name="Acsub12" )
covAcsub21	<-mxAlgebra( expression= (ac + (acM%x%ELA2))%*% t(ac + (acM%x%ELA1)), name="Acsub21" )

covEcsub11	<-mxAlgebra( expression= (ec + (ecM%x%ELA1))%*% t(ec + (ecM%x%ELA1)), name="Ecsub11" )
covEcsub22	<-mxAlgebra( expression= (ec + (ecM%x%ELA2))%*% t(ec + (ecM%x%ELA2)), name="Ecsub22" )
covEcsub12	<-mxAlgebra( expression= (ec + (ecM%x%ELA1))%*% t(ec + (ecM%x%ELA2)), name="Ecsub12" )
covEcsub21	<-mxAlgebra( expression= (ec + (ecM%x%ELA2))%*% t(ec + (ecM%x%ELA1)), name="Ecsub21" )

covPcsub11	<-mxAlgebra( expression= Acsub11+Ecsub11, name="Vcsub11" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)
covPcsub22	<-mxAlgebra( expression= Acsub22+Ecsub22, name="Vcsub22" ) 

covPc11	<-mxAlgebra( expression= pc %*% t(pc), name="Pc11" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components

covPC11	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Vcsub11)), name="Vc11") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPc22	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Vcsub22)), name="Vc22") 

covPcMz12	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Acsub12)), name="Vcmz12") #I specify the MZ between-twin covariance - excluding E parameters
covPcMz21	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Acsub21)), name="Vcmz21") 

covPcDz12	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z21), rbind(Z12, .5%x%Acsub12)), name="Vcdz12") #I specify the DZ between-twin covariance - specifying half of A and excluding E
covPcDz21	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z21), rbind(Z12, .5%x%Acsub21)), name="Vcdz21")
 
# Generate Covariance of Latent factor model Including Causal Paths between factors
Id3		<-mxMatrix(type="Iden",	nrow=3, ncol=3, free=F, name="I3" )
covFVc1	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%ELA1))) %&% Vc11, name ="FVc1")
covFVc2	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%ELA2))) %&% Vc22, name ="FVc2")

covFcMz12	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%ELA1))) %*% Vcmz12 %*% solve(I3-(PhC+(PhCM%x%ELA2))), name ="Fcmz12")
covFcMz21	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%ELA2))) %*% Vcmz21 %*% solve(I3-(PhC+(PhCM%x%ELA1))), name ="Fcmz21")

covFcDz12	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%ELA1))) %*% Vcdz12 %*% solve(I3-(PhC+(PhCM%x%ELA2))), name ="Fcdz12")
covFcDz21	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%ELA2))) %*% Vcdz21 %*% solve(I3-(PhC+(PhCM%x%ELA1))), name ="Fcdz21")

# Constraint on total variance of Ordinal variable (A+C+E=1)
varL1		<- mxConstraint( expression=FVc1[3,3]==1, name="L1" )
varL2		<- mxConstraint( expression=FVc2[3,3]==1, name="L2" )

FcovMZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc1, Fcmz12), cbind(Fcmz21, FVc2))) , name="expFCovMZ" )
FcovDZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc1, Fcdz12), cbind(Fcdz21, FVc2))) , name="expFCovDZ" )
FcovSIB1	<-mxAlgebra( expression= (FactL  %&% FVc1 ), name="expFCovSIB1" )
FcovSIB2	<-mxAlgebra( expression= (FactL  %&% FVc2 ), name="expFCovSIB2" )

SpcovMZ	<-mxAlgebra( expression= rbind (cbind(Vs, As), cbind(As, Vs)) , name="expSpCovMZ" )
SpcovDZ	<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As), cbind(.5%x%As, Vs)) , name="expSpCovDZ" )

TOTcovMZ	<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )
TOTcovSIB1	<-mxAlgebra( expression= expFCovSIB1 + Vs , name="TOTexpCovSIB1" )
TOTcovSIB2	<-mxAlgebra( expression= expFCovSIB2 + Vs , name="TOTexpCovSIB2" )

# *******************************************************************************************************
# Calculator

# Standardize the Total var/covariances matrices of the observed variables
Id10		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I10" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= sqrt(As/( (FactL %&% FVc1) +Vs)), name="stAs" )
stcovEs	<-mxAlgebra( expression= sqrt(Es/( (FactL %&% FVc1) +Vs)), name="stEs" )

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% FVc1 / TOTexpCovMZ[1:5,1:5])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
ctdataMZ	<- mxData( observed=mz12Data, type="raw" )
ctdataDZ	<- mxData( observed=dz12Data, type="raw" )
ctdataSIB1	<- mxData( observed=sib1Data, type="raw" )
ctdataSIB2	<- mxData( observed=sib2Data, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
objSIB1	<- mxExpectationNormal( covariance="TOTexpCovSIB1", means="expMeanS", dimnames=selVarsS1, thresholds="expThres1", threshnames=c("SO1"))
objSIB2	<- mxExpectationNormal( covariance="TOTexpCovSIB2", means="expMeanS", dimnames=selVarsS2, thresholds="expThres2", threshnames=c("SO2"))

fitFunction <- mxFitFunctionML()
 
# Combine Groups

cistFL	<-mxCI (c ('StandFact','PhCM','PhC','acM','ecM'))
cistVs	<-mxCI (c ('stAs[3,3]','stAs[4,4]','stEs[1,1]','stEs[2,2]','stEs[3,3]','stEs[4,4]') ) 	# standardized var comp from specific Factors

pars1		<-list(Means,Load,LoadTw,PhCaus,PhCausMod,PathsAs,PathsEs,covAs,covEs,covPs,Id3,Id10,Id2, cistFL, cistVs)
pars2		<-list(PathsAcsub,PathsEcsub,PathsAcMod,PathsEcMod,PathsP11,Ze21,Ze12)
pars3		<-list(Ze11,covAcsub11,covAcsub22,covAcsub12,covAcsub21,covEcsub11,covEcsub22,covEcsub12,covEcsub21,
				covPcsub11,covPcsub22,covPc11,covPC11,covPc22,covFVc1,covFVc2)
pars4		<-list(obsage1, obsage2, obssex1, obssex2, obsela1, obsela2, betaA, betaS, betaE, Thr, inc, Thres)
parsst	<-list(stcovAs, stcovEs)
parsmz	<-list(covPcMz12,covPcMz21,covFcMz12,covFcMz21)
parsdz	<-list(covPcDz12,covPcDz21,covFcDz12,covFcDz21)

pars1a	<-list(MeansS,Load,PhCaus,PhCausMod,PathsAs,PathsEs,covAs,covEs,covPs,Id3)
pars1b	<-list(PathsAcsub,PathsEcsub,PathsAcMod,PathsEcMod,PathsP11,Ze21,Ze12)
pars1c	<-list(covAcsub11,covEcsub11)
pars1d	<-list(covPcsub11,covPc11,covPC11,covFVc1)
pars1e	<-list(obsage1, obssex1, obsela1, betaA, betaS, betaE, Thr, inc, Thres1)

pars2a	<-list(MeansS,Load,PhCaus,PhCausMod,PathsAs,PathsEs,covAs,covEs,covPs,Id3)
pars2b	<-list(PathsAcsub,PathsEcsub,PathsAcMod,PathsEcMod,PathsP11,Ze21,Ze12)
pars2c	<-list(covAcsub22,covEcsub22)
pars2d	<-list(covPcsub22,covPc11,covPc22,covFVc2)
pars2e	<-list(obsage2, obssex2, obsela2, betaA, betaS, betaE, Thr, inc, Thres2)

#cistFL	<-mxCI (c ('MZ.StandFact','MZ.PhCM','MZ.PhC','MZ.acM','MZ.ecM'))
cistFL	<-mxCI (c ('MZ.PhCM'))
#cistVs	<-mxCI (c ('MZ.stAs[3,3]','MZ.stEs[1,1]','MZ.stEs[2,2]','MZ.stEs[3,3]') ) 	# standardized var comp from specific Factors

modelMZ	<-mxModel(parsmz, pars1, pars2, pars3, pars4, FcovMZ, SpcovMZ, TOTcovMZ, ctdataMZ, objMZ, Rfactmz, parsst, fitFunction, StFL, varL1, varL2, name="MZ" )
modelDZ	<-mxModel(parsdz, pars1, pars2, pars3, pars4, FcovDZ, SpcovDZ, TOTcovDZ, ctdataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )
modelSIB1	<-mxModel(pars1a, pars1b, pars1c, pars1d, pars1e, FcovSIB1, TOTcovSIB1, ctdataSIB1, objSIB1, fitFunction, name="SIB1" )
modelSIB2	<-mxModel(pars2a, pars2b, pars2c, pars2d, pars2e, FcovSIB2, TOTcovSIB2, ctdataSIB2, objSIB2, fitFunction, name="SIB2" )

minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective + SIB1.objective + SIB2.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
MRDoCMod3Model	<-mxModel("MRDoCMod3", modelMZ, modelDZ, modelSIB1, modelSIB2, minus2ll, obj, cistFL) 

# --------------------------------------------------------------------------------------------------------------------------------
# 6a RUN ACEMs Factor Model with phenotypic causal mediation paths by Zygosity
#MRDoCMod3Fit	<-mxTryHardOrdinal(MRDoCMod3Model,intervals=TRUE,checkHess=FALSE,OKstatuscodes=c(0,1,5,6))
MRDoCMod3Fit	<-mxTryHardOrdinal(MRDoCMod3Model,intervals=TRUE)
(MRDoCMod3Summ	<-summary(MRDoCMod3Fit, verbose=F))

mxEval(MZ.FactcorMZ, MRDoCMod3Fit)
mxEval(DZ.FactcorDZ, MRDoCMod3Fit)

mxEval(MZ.FVc1, MRDoCMod3Fit)
mxEval(MZ.FVc2, MRDoCMod3Fit)

mxEval(MZ.StandFact, MRDoCMod3Fit)

mxEval(MZ.Vs, MRDoCMod3Fit)

mxEval(MZ.Acsub11, MRDoCMod3Fit)
mxEval(MZ.Ecsub11, MRDoCMod3Fit)
mxEval(MZ.acM, MRDoCMod3Fit)
mxEval(MZ.ecM, MRDoCMod3Fit)
mxEval(MZ.PhCM, MRDoCMod3Fit)
mxEval(MZ.PhC, MRDoCMod3Fit)

mxEval(MZ.Vcsub11, MRDoCMod3Fit)
mxEval(MZ.Vc11, MRDoCMod3Fit)
mxEval(MZ.expFCovMZ, MRDoCMod3Fit)
mxEval(DZ.expFCovDZ, MRDoCMod3Fit)

mxEval(MZ.stAs, MRDoCMod3Fit)
mxEval(MZ.stEs, MRDoCMod3Fit)


MRDoCNoMod1cModel	<- mxModel(MRDoCMod3Fit, name="MRDoCNoMod1c")
MRDoCNoMod1cModel	<- omxSetParameters(MRDoCNoMod1cModel, labels=c('c2on3mod'), free=FALSE, values=0)
MRDoCNoMod1cFit	<- mxTryHardOrdinal(MRDoCNoMod1cModel, intervals=F)
(MRDoCNoMod1cSum	<- summary(MRDoCNoMod1cFit, verbose=F))

mxCompare(MRDoCMod3Fit, MRDoCNoMod1cFit)


# 
#****************************************************************************************************************************
# __(VId) MH>SO_____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL for PD>SO incorporating moderation by CGN
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Phenotypic Factors: F1>F2>F3 & F1>F3;
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
#_____________________________________________________________________________________________________________________________

nv		<- 5				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 3				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nvo 		<- 1     			# number of ordinal variables per twin
nvc 		<- nv-nvo  			# number of continuous variables per twin
poso 		<- nvo 			# position where ordinal variables start
nth		<- 4				# number of max thresholds
ninc 		<- nth-1 			# number of max increments
ncovariates <- 3 				# number of covariates
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups	<- c("mz", "dz")
Vars		<- c('PRSAnx','PRSDep','Depcg4','Anxcg4','SO')
selVars	<- c('PRSAnx1','PRSDep1','Depcg41','Anxcg41','SO1',
		     'PRSAnx2','PRSDep2','Depcg42','Anxcg42','SO2')
selVarsS1	<- c('PRSAnx1','PRSDep1','Depcg41','Anxcg41','SO1')
selVarsS2	<- c('PRSAnx2','PRSDep2','Depcg42','Anxcg42','SO2')
useVars	<- c('PRSAnx1','PRSDep1','Depcg41','Anxcg41','SO1',
		     'PRSAnx2','PRSDep2','Depcg42','Anxcg42','SO2','age1','sex1','CGN41','age2','sex2','CGN42')
useVars1	<- c('PRSAnx1','PRSDep1','Depcg41','Anxcg41','SO1','age1','sex1','CGN41')
useVars2	<- c('PRSAnx2','PRSDep2','Depcg42','Anxcg42','SO2','age2','sex2','CGN42')

mzData	<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , selVars)
dzData	<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , selVars)

mz12Data	<- TWINdata2[TWINdata2$zyg1==c(1,3) & (!is.na(TWINdata2$CGN41) & !is.na(TWINdata2$CGN42)), useVars]
dz12Data	<- TWINdata2[TWINdata2$zyg1==c(2,4,5) & (!is.na(TWINdata2$CGN41) & !is.na(TWINdata2$CGN42)), useVars]
sib1Data	<- TWINdata2[!is.na(TWINdata2$CGN41) & is.na(TWINdata2$CGN42), useVars1]
sib2Data	<- TWINdata2[is.na(TWINdata2$CGN41) & !is.na(TWINdata2$CGN42), useVars2]

psych::describe(mzData)
psych::describe(dzData)

psych::describe(mz12Data)
psych::describe(dz12Data)
psych::describe(sib1Data)
psych::describe(sib2Data)

mz12Data$SO1[is.na(mz12Data$age1)] <- NA
mz12Data$SO2[is.na(mz12Data$age2)] <- NA
dz12Data$SO1[is.na(dz12Data$age1)] <- NA
dz12Data$SO2[is.na(dz12Data$age2)] <- NA
sib1Data$SO1[is.na(sib1Data$age1)] <- NA
sib2Data$SO2[is.na(sib2Data$age2)] <- NA

mz12Data$SO1[is.na(mz12Data$sex1)] <- NA
mz12Data$SO2[is.na(mz12Data$sex2)] <- NA
dz12Data$SO1[is.na(dz12Data$sex1)] <- NA
dz12Data$SO2[is.na(dz12Data$sex2)] <- NA
sib1Data$SO1[is.na(sib1Data$sex1)] <- NA
sib2Data$SO2[is.na(sib2Data$sex2)] <- NA

mz12Data$SO1[is.na(mz12Data$CGN41)] <- NA
mz12Data$SO2[is.na(mz12Data$CGN42)] <- NA
dz12Data$SO1[is.na(dz12Data$CGN41)] <- NA
dz12Data$SO2[is.na(dz12Data$CGN42)] <- NA
sib1Data$SO1[is.na(sib1Data$CGN41)] <- NA
sib2Data$SO2[is.na(sib2Data$CGN42)] <- NA
mz12Data$Depcg41[is.na(mz12Data$CGN41)] <- NA
mz12Data$Depcg42[is.na(mz12Data$CGN42)] <- NA
dz12Data$Depcg41[is.na(dz12Data$CGN41)] <- NA
dz12Data$Depcg42[is.na(dz12Data$CGN42)] <- NA
sib1Data$Depcg41[is.na(sib1Data$CGN41)] <- NA
sib2Data$Depcg42[is.na(sib2Data$CGN42)] <- NA
mz12Data$Anxcg41[is.na(mz12Data$CGN41)] <- NA
mz12Data$Anxcg42[is.na(mz12Data$CGN42)] <- NA
dz12Data$Anxcg41[is.na(dz12Data$CGN41)] <- NA
dz12Data$Anxcg42[is.na(dz12Data$CGN42)] <- NA
sib1Data$Anxcg41[is.na(sib1Data$CGN41)] <- NA
sib2Data$Anxcg42[is.na(sib2Data$CGN42)] <- NA

mz12Data$age1[is.na(mz12Data$age1)] <- 999
mz12Data$age2[is.na(mz12Data$age2)] <- 999
dz12Data$age1[is.na(dz12Data$age1)] <- 999
dz12Data$age2[is.na(dz12Data$age2)] <- 999
sib1Data$age1[is.na(sib1Data$age1)] <- 999
sib2Data$age2[is.na(sib2Data$age2)] <- 999

mz12Data$sex1[is.na(mz12Data$sex1)] <- 999
mz12Data$sex2[is.na(mz12Data$sex2)] <- 999
dz12Data$sex1[is.na(dz12Data$sex1)] <- 999
dz12Data$sex2[is.na(dz12Data$sex2)] <- 999
sib1Data$sex1[is.na(sib1Data$sex1)] <- 999
sib2Data$sex2[is.na(sib2Data$sex2)] <- 999

mz12Data$CGN41[is.na(mz12Data$CGN41)] <- 999
mz12Data$CGN42[is.na(mz12Data$CGN42)] <- 999
dz12Data$CGN41[is.na(dz12Data$CGN41)] <- 999
dz12Data$CGN42[is.na(dz12Data$CGN42)] <- 999
sib1Data$CGN41[is.na(sib1Data$CGN41)] <- 999
sib2Data$CGN42[is.na(sib2Data$CGN42)] <- 999

psych::describe(mz12Data)
psych::describe(dz12Data)
psych::describe(sib1Data)
psych::describe(sib2Data)

# ******************************************************************************************************************************************
# (5) van der Sluis correction
# ******************************************************************************************************************************************
# I control for the cotwin and sibling moderation variable (CTQ) via residualisation (for SOI, SOP, Dep, And and RSB) (van der Sluis et al, 2012). 
# This controls for any confounding due to main effects of the co-twin/sibling's moderator variable.
# I use this approach to reduce the complexity of the model
# *********************************************************************************************

mz12Data$Depcg41	<- residuals(lm(mz12Data$Depcg41 ~ mz12Data$CGN42, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Depcg41	<- (mz12Data$Depcg41 + 3.27)
psych::describe(mz12Data$Depcg41)
hist(mz12Data$Depcg41)

mz12Data$Depcg42	<- residuals(lm(mz12Data$Depcg42 ~ mz12Data$CGN41, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Depcg42	<- (mz12Data$Depcg42 + 3.32)
psych::describe(mz12Data$Depcg42)
hist(mz12Data$Depcg42)

dz12Data$Depcg41	<- residuals(lm(dz12Data$Depcg41 ~ dz12Data$CGN42, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Depcg41	<- (dz12Data$Depcg41 + 3.30)
psych::describe(dz12Data$Depcg41)
hist(dz12Data$Depcg41)

dz12Data$Depcg42	<- residuals(lm(dz12Data$Depcg42 ~ dz12Data$CGN41, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Depcg42	<- (dz12Data$Depcg42 + 3.35)
psych::describe(dz12Data$Depcg42)
hist(dz12Data$Depcg42)

mz12Data$Anxcg41	<- residuals(lm(mz12Data$Anxcg41 ~ mz12Data$CGN42, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Anxcg41	<- (mz12Data$Anxcg41 + 8.24)
psych::describe(mz12Data$Anxcg41)
hist(mz12Data$Anxcg41)

mz12Data$Anxcg42	<- residuals(lm(mz12Data$Anxcg42 ~ mz12Data$CGN41, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
mz12Data$Anxcg42	<- (mz12Data$Anxcg42 + 8.22)
psych::describe(mz12Data$Anxcg42)
hist(mz12Data$Anxcg42)

dz12Data$Anxcg41	<- residuals(lm(dz12Data$Anxcg41 ~ dz12Data$CGN42, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Anxcg41	<- (dz12Data$Anxcg41 + 8.26)
psych::describe(dz12Data$Anxcg41)
hist(dz12Data$Anxcg41)

dz12Data$Anxcg42	<- residuals(lm(dz12Data$Anxcg42 ~ dz12Data$CGN41, na.action="na.exclude")) ## Repeat for each of the moderators for each of the variables
dz12Data$Anxcg42	<- (dz12Data$Anxcg42 + 8.28)
psych::describe(dz12Data$Anxcg42)
hist(dz12Data$Anxcg42)

psych::describe(mz12Data)
psych::describe(dz12Data)
psych::describe(sib1Data)
psych::describe(sib2Data)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- c(5,5,4.5,4.1,0))
(PatM		<- c(TRUE,T,TRUE,TRUE,F))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEs	<- c('es1','es1','es3','es3','es5'))
(LabAs	<- c('as1','as1','as3','as3','es5'))

PatSpe	<- c(T,T,TRUE,TRUE,F)
PatSpac	<- c(F,F,TRUE,TRUE,F)
StSpe		<- c(.5,.5,.5,.5,0)
StSpa		<- c(0,0,.5,.5,0)

# all 1st loadings fixed to 1
PatFl		<- c(F,T,F,F,F,			
		     F,F,F,T,F,
		     F,F,F,F,F)

StFl		<- c(1,.5,0,0,0,
		     0,0,1,.5,0,
		     0,0,0,0,1)

LabFl		<- c('l1','l2',NA,NA,NA,
	 	     NA,NA,'l3','l4',NA,
	 	     NA,NA,NA,NA,'l5')

PatPhC	<- c(F,T,T,
		     F,F,T,
		     F,F,F)

StPhC		<- c(0,.3,.01,
		     0,0,.1,
		     0,0,0)

LabPhC	<- c(NA,'c1on2','c1on3',
		     NA,NA,'c2on3',
		     NA,NA,NA)	 

PatPhCM	<- c(F,F,F,
		     F,F,T,
		     F,F,F)

StPhCM	<- c(0,0,0,
		     0,0,-.01,
		     0,0,0)

LabPhCM	<- c(NA,'c1on2mod','c1on3mod',
		     NA,NA,'c2on3mod',	
		     NA,NA,NA)	 

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmean,Stmean), labels=c(mLabs,mLabs), name="expMean") 
MeansS	<-mxMatrix("Full", 1, nv, free=c(PatM), values=c(Stmean), labels=c(mLabs), name="expMeanS") 

# Threshold and covariates
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')
LabCovC 	<-c('BcgnThSO','BcgnThSO','BcgnThSO','BcgnThSO')

ThPat		<-c(T,T,T,T)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

obscgn1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.CGN41"), name="CGN41")
obscgn2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.CGN42"), name="CGN42")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.1, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.1, labels=LabCovS, name="BsexTH" )
betaC		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.1, labels=LabCovC, name="BcgnTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1 + BcgnTH%x%CGN41,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2 + BcgnTH%x%CGN42), name="expThres")
Thres1	<-mxAlgebra( expression= Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1 + BcgnTH%x%CGN41, name="expThres1")
Thres2	<-mxAlgebra( expression= Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2 + BcgnTH%x%CGN42, name="expThres2")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )
PhCausMod	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhCM, values=StPhCM, labels=LabPhCM, name="PhCM" )

# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpa, labels=LabAs, name="as" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpe, labels=LabEs, name="es" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.8,.1,.5), labels=c("ac22","ac32","ac33"), name="ac" ) # Component paths for factors 2 and 3
PathsEcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(.8,0,.8), labels=c("ec22","ec32","ec33"), name="ec" )

PathsAcMod	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(.01,0,.1), labels=c("acMod22","acMod32","acMod33"), name="acM" ) # Component paths for factors 2 and 3
PathsEcMod	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=c(T,F,T), values=c(0.1,0,0.1), labels=c("ecMod22","ecMod32","ecMod33"), name="ecM" )

PathsP11	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11", name="pc" ) # SD path for factor 3 (the PRS factor)
Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  #Padding

Ze21		<-mxMatrix(type="Zero",	nrow=2, ncol=1, free=F, name="Z21" )  #Padding
Ze12		<-mxMatrix(type="Zero",	nrow=1, ncol=2, free=F, name="Z12" )  #Padding

covAcsub11	<-mxAlgebra( expression= (ac + (acM%x%CGN41))%*% t(ac + (acM%x%CGN41)), name="Acsub11" )
covAcsub22	<-mxAlgebra( expression= (ac + (acM%x%CGN42))%*% t(ac + (acM%x%CGN42)), name="Acsub22" )
covAcsub12	<-mxAlgebra( expression= (ac + (acM%x%CGN41))%*% t(ac + (acM%x%CGN42)), name="Acsub12" )
covAcsub21	<-mxAlgebra( expression= (ac + (acM%x%CGN42))%*% t(ac + (acM%x%CGN41)), name="Acsub21" )

covEcsub11	<-mxAlgebra( expression= (ec + (ecM%x%CGN41))%*% t(ec + (ecM%x%CGN41)), name="Ecsub11" )
covEcsub22	<-mxAlgebra( expression= (ec + (ecM%x%CGN42))%*% t(ec + (ecM%x%CGN42)), name="Ecsub22" )
covEcsub12	<-mxAlgebra( expression= (ec + (ecM%x%CGN41))%*% t(ec + (ecM%x%CGN42)), name="Ecsub12" )
covEcsub21	<-mxAlgebra( expression= (ec + (ecM%x%CGN42))%*% t(ec + (ecM%x%CGN41)), name="Ecsub21" )

covPcsub11	<-mxAlgebra( expression= Acsub11+Ecsub11, name="Vcsub11" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)
covPcsub22	<-mxAlgebra( expression= Acsub22+Ecsub22, name="Vcsub22" ) 

covPc11	<-mxAlgebra( expression= pc %*% t(pc), name="Pc11" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components

covPC11	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Vcsub11)), name="Vc11") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPc22	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Vcsub22)), name="Vc22") 

covPcMz12	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Acsub12)), name="Vcmz12") #I specify the MZ between-twin covariance - excluding E parameters
covPcMz21	<-mxAlgebra(cbind(rbind(Pc11,Z21), rbind(Z12, Acsub21)), name="Vcmz21") 

covPcDz12	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z21), rbind(Z12, .5%x%Acsub12)), name="Vcdz12") #I specify the DZ between-twin covariance - specifying half of A and excluding E
covPcDz21	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z21), rbind(Z12, .5%x%Acsub21)), name="Vcdz21")
 
# Generate Covariance of Latent factor model Including Causal Paths between factors
Id3		<-mxMatrix(type="Iden",	nrow=3, ncol=3, free=F, name="I3" )
covFVc1	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%CGN41))) %&% Vc11, name ="FVc1")
covFVc2	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%CGN42))) %&% Vc22, name ="FVc2")

covFcMz12	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%CGN41))) %*% Vcmz12 %*% solve(I3-(PhC+(PhCM%x%CGN42))), name ="Fcmz12")
covFcMz21	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%CGN42))) %*% Vcmz21 %*% solve(I3-(PhC+(PhCM%x%CGN41))), name ="Fcmz21")

covFcDz12	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%CGN41))) %*% Vcdz12 %*% solve(I3-(PhC+(PhCM%x%CGN42))), name ="Fcdz12")
covFcDz21	<-mxAlgebra( expression= solve(I3-(PhC+(PhCM%x%CGN42))) %*% Vcdz21 %*% solve(I3-(PhC+(PhCM%x%CGN41))), name ="Fcdz21")

# Constraint on total variance of Ordinal variable (A+C+E=1)
varL1		<- mxConstraint( expression=FVc1[3,3]==1, name="L1" )
varL2		<- mxConstraint( expression=FVc2[3,3]==1, name="L2" )

FcovMZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc1, Fcmz12), cbind(Fcmz21, FVc2))) , name="expFCovMZ" )
FcovDZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc1, Fcdz12), cbind(Fcdz21, FVc2))) , name="expFCovDZ" )
FcovSIB1	<-mxAlgebra( expression= (FactL  %&% FVc1 ), name="expFCovSIB1" )
FcovSIB2	<-mxAlgebra( expression= (FactL  %&% FVc2 ), name="expFCovSIB2" )

SpcovMZ	<-mxAlgebra( expression= rbind (cbind(Vs, As), cbind(As, Vs)) , name="expSpCovMZ" )
SpcovDZ	<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As), cbind(.5%x%As, Vs)) , name="expSpCovDZ" )

TOTcovMZ	<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )
TOTcovSIB1	<-mxAlgebra( expression= expFCovSIB1 + Vs , name="TOTexpCovSIB1" )
TOTcovSIB2	<-mxAlgebra( expression= expFCovSIB2 + Vs , name="TOTexpCovSIB2" )

# *******************************************************************************************************
# Calculator

# Standardize the Total var/covariances matrices of the observed variables
Id10		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I10" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= sqrt(As/( (FactL %&% FVc1) +Vs)), name="stAs" )
stcovEs	<-mxAlgebra( expression= sqrt(Es/( (FactL %&% FVc1) +Vs)), name="stEs" )

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% FVc1 / TOTexpCovMZ[1:5,1:5])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
cgdataMZ	<- mxData( observed=mz12Data, type="raw" )
cgdataDZ	<- mxData( observed=dz12Data, type="raw" )
#cgdataSIB1	<- mxData( observed=sib1Data, type="raw" )
#cgdataSIB2	<- mxData( observed=sib2Data, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
#objSIB1	<- mxExpectationNormal( covariance="TOTexpCovSIB1", means="expMeanS", dimnames=selVarsS1, thresholds="expThres1", threshnames=c("SO1"))
#objSIB2	<- mxExpectationNormal( covariance="TOTexpCovSIB2", means="expMeanS", dimnames=selVarsS2, thresholds="expThres2", threshnames=c("SO2"))

fitFunction <- mxFitFunctionML()
 
# Combine Groups

pars1		<-list(Means,Load,LoadTw,PhCaus,PhCausMod,PathsAs,PathsEs,covAs,covEs,covPs,Id3,Id10,Id2)
pars2		<-list(PathsAcsub,PathsEcsub,PathsAcMod,PathsEcMod,PathsP11,Ze21,Ze12)
pars3		<-list(Ze11,covAcsub11,covAcsub22,covAcsub12,covAcsub21,covEcsub11,covEcsub22,covEcsub12,covEcsub21,
				covPcsub11,covPcsub22,covPc11,covPC11,covPc22,covFVc1,covFVc2)
pars4		<-list(obsage1, obsage2, obssex1, obssex2, obscgn1, obscgn2, betaA, betaS, betaC, Thr, inc, Thres)
parsst	<-list(stcovAs, stcovEs)
parsmz	<-list(covPcMz12,covPcMz21,covFcMz12,covFcMz21)
parsdz	<-list(covPcDz12,covPcDz21,covFcDz12,covFcDz21)

pars1a	<-list(MeansS,Load,PhCaus,PhCausMod,PathsAs,PathsEs,covAs,covEs,covPs,Id3)
pars1b	<-list(PathsAcsub,PathsEcsub,PathsAcMod,PathsEcMod,PathsP11,Ze21,Ze12)
pars1c	<-list(covAcsub11,covEcsub11)
pars1d	<-list(covPcsub11,covPc11,covPC11,covFVc1)
pars1e	<-list(obsage1, obssex1, obscgn1, betaA, betaS, betaC, Thr, inc, Thres1)

pars2a	<-list(MeansS,Load,PhCaus,PhCausMod,PathsAs,PathsEs,covAs,covEs,covPs,Id3)
pars2b	<-list(PathsAcsub,PathsEcsub,PathsAcMod,PathsEcMod,PathsP11,Ze21,Ze12)
pars2c	<-list(covAcsub22,covEcsub22)
pars2d	<-list(covPcsub22,covPc11,covPc22,covFVc2)
pars2e	<-list(obsage2, obssex2, obscgn2, betaA, betaS, betaC, Thr, inc, Thres2)

modelMZ	<-mxModel(parsmz, pars1, pars2, pars3, pars4, FcovMZ, SpcovMZ, TOTcovMZ, cgdataMZ, objMZ, Rfactmz, parsst, fitFunction, StFL, varL1, varL2, name="MZ" )
modelDZ	<-mxModel(parsdz, pars1, pars2, pars3, pars4, FcovDZ, SpcovDZ, TOTcovDZ, cgdataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )
#modelSIB1	<-mxModel(pars1a, pars1b, pars1c, pars1d, pars1e, FcovSIB1, TOTcovSIB1, cgdataSIB1, objSIB1, fitFunction, name="SIB1" )
#modelSIB2	<-mxModel(pars2a, pars2b, pars2c, pars2d, pars2e, FcovSIB2, TOTcovSIB2, cgdataSIB2, objSIB2, fitFunction, name="SIB2" )

#minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective + SIB1.objective + SIB2.objective, name="m2LL" )
minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )

#cistFL	<-mxCI (c ('StandFact','PhCM','PhC','acM','ecM'))
#cistVs	<-mxCI (c ('stAs[3,3]','stAs[4,4]','stEs[1,1]','stEs[2,2]','stEs[3,3]','stEs[4,4]') ) 	# standardized var comp from specific Factors
#cistFL	<-mxCI (c ('MZ.StandFact','MZ.PhCM','MZ.PhC','MZ.acM','MZ.ecM'))
cistFL	<-mxCI (c ('MZ.PhCM'))
#cistVs	<-mxCI (c ('MZ.stAs[3,3]','MZ.stEs[1,1]','MZ.stEs[2,2]','MZ.stEs[3,3]') ) 	# standardized var comp from specific Factors

#MRDoCMod3Model	<-mxModel("MRDoCMod3", modelMZ, modelDZ, modelSIB1, modelSIB2, minus2ll, obj, cistFL) 
MRDoCMod3Model	<-mxModel("MRDoCMod3", modelMZ, modelDZ, minus2ll, obj, cistFL) 

# --------------------------------------------------------------------------------------------------------------------------------
# 6a RUN ACEMs Factor Model with phenotypic causal mediation paths by Zygosity
#MRDoCMod3Fit	<-mxTryHardOrdinal(MRDoCMod3Model,intervals=TRUE,checkHess=FALSE,OKstatuscodes=c(0,1,5,6))
MRDoCMod3Fit	<-mxTryHardOrdinal(MRDoCMod3Model,intervals=TRUE)
(MRDoCMod3Summ	<-summary(MRDoCMod3Fit, verbose=F))

mxEval(MZ.FactcorMZ, MRDoCMod3Fit)
mxEval(DZ.FactcorDZ, MRDoCMod3Fit)

mxEval(MZ.FVc1, MRDoCMod3Fit)
mxEval(MZ.FVc2, MRDoCMod3Fit)

mxEval(MZ.StandFact, MRDoCMod3Fit)

mxEval(MZ.Vs, MRDoCMod3Fit)

mxEval(MZ.Acsub11, MRDoCMod3Fit)
mxEval(MZ.Ecsub11, MRDoCMod3Fit)
mxEval(MZ.acM, MRDoCMod3Fit)
mxEval(MZ.ecM, MRDoCMod3Fit)
mxEval(MZ.PhCM, MRDoCMod3Fit)
mxEval(MZ.PhC, MRDoCMod3Fit)

mxEval(MZ.Vcsub11, MRDoCMod3Fit)
mxEval(MZ.Vc11, MRDoCMod3Fit)
mxEval(MZ.expFCovMZ, MRDoCMod3Fit)
mxEval(DZ.expFCovDZ, MRDoCMod3Fit)

mxEval(MZ.stAs, MRDoCMod3Fit)
mxEval(MZ.stEs, MRDoCMod3Fit)


MRDoCNoMod1cModel	<- mxModel(MRDoCMod3Fit, name="MRDoCNoMod1c")
MRDoCNoMod1cModel	<- omxSetParameters(MRDoCNoMod1cModel, labels=c('c2on3mod'), free=FALSE, values=0)
MRDoCNoMod1cFit	<- mxTryHardOrdinal(MRDoCNoMod1cModel, intervals=F)
(MRDoCNoMod1cSum	<- summary(MRDoCNoMod1cFit, verbose=F))

mxCompare(MRDoCMod3Fit, MRDoCNoMod1cFit)



# 
#****************************************************************************************************************************
# __(VIIa)_____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL for SO>PD with mediation by Vict
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Phenotypic Factors: F1>F2>F3 & F1>F3;
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
#_____________________________________________________________________________________________________________________________

nv		<- 5				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 4				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nvo 		<- 1     			# number of ordinal variables per twin
nvc 		<- nv-nvo  			# number of continuous variables per twin
poso 		<- nvo 			# position where ordinal variables start
nth		<- 4				# number of max thresholds
ninc 		<- nth-1 			# number of max increments
ncovariates <- 2 				# number of covariates
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups	<- c("mz", "dz")
Vars		<- c('PRSSO','SO','Vict','Dep','Anx')
selVars	<- c('PRSSO1','SO1','Vict1','Dep1','Anx1',
		     'PRSSO2','SO2','Vict2','Dep2','Anx2')
useVars	<- c('PRSSO1','SO1','Vict1','Dep1','Anx1',
		     'PRSSO2','SO2','Vict2','Dep2','Anx2','age1','sex1','age2','sex2')

mzData		<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , useVars)
dzData		<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , useVars)

psych::describe(mzData)
psych::describe(dzData)

mzData$SO1[is.na(mzData$age1)] <- NA
mzData$SO2[is.na(mzData$age2)] <- NA
dzData$SO1[is.na(dzData$age1)] <- NA
dzData$SO2[is.na(dzData$age2)] <- NA

mzData$SO1[is.na(mzData$sex1)] <- NA
mzData$SO2[is.na(mzData$sex2)] <- NA
dzData$SO1[is.na(dzData$sex1)] <- NA
dzData$SO2[is.na(dzData$sex2)] <- NA

mzData$age1[is.na(mzData$age1)] <- 999
mzData$age2[is.na(mzData$age2)] <- 999
dzData$age1[is.na(dzData$age1)] <- 999
dzData$age2[is.na(dzData$age2)] <- 999

mzData$sex1[is.na(mzData$sex1)] <- 999
mzData$sex2[is.na(mzData$sex2)] <- 999
dzData$sex1[is.na(dzData$sex1)] <- 999
dzData$sex2[is.na(dzData$sex2)] <- 999

psych::describe(mzData)
psych::describe(dzData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- c(5,0,2.5,4.5,4.1))
(PatM		<- c(TRUE,F,TRUE,TRUE,TRUE))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEs	<- c('es1','es2','es3','es4','es4'))
(LabAs	<- c('as1','as2','as3','as4','as4'))
(LabCs	<- c('cs1','cs2','cs3','cs4','cs4'))

PatSpe	<- c(F,F,F,TRUE,TRUE)
PatSpac	<- c(F,F,F,TRUE,TRUE)
StSpa		<- c(0,0,0,.5,.5)
StSpc		<- c(0,0,0,.5,.5)
StSpe		<- c(0,0,0,.5,.5)

# all 1st loadings fixed to 1
PatFl		<- c(F,F,F,F,F,			
		     F,F,F,F,F,
		     F,F,F,F,F,
		     F,F,F,F,T)

StFl		<- c(1,0,0,0,0,
		     0,1,0,0,0,
		     0,0,1,0,0,
		     0,0,0,1,.5)

LabFl		<- c('l1',NA,NA,NA,NA,
	 	     NA,'l2',NA,NA,NA,
	 	     NA,NA,'l3',NA,NA,
	 	     NA,NA,NA,'l4','l5')

PatPhC	<- c(F,T,T,T,
		     F,F,T,T,
		     F,F,F,T,
		     F,F,F,F)

StPhC		<- c(0,.3,.3,.3,
		     0,0,.3,.3,
		     0,0,0,.3,
		     0,0,0,0)

LabPhC	<- c(NA,'c1on2','c1on3','c1on4',
		     NA,NA,'c2on3','c2on4',
		     NA,NA,NA,'c3on4',
		     NA,NA,NA,NA)	 

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmean,Stmean), labels=c(mLabs,mLabs), name="expMean") 

# Threshold and covariates
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')

ThPat		<-c(T,T,T,T)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.3, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.4, labels=LabCovS, name="BsexTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2), name="expThres")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )

# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpa, labels=LabAs, name="as" )
PathsCs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpc, labels=LabCs, name="cs" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpe, labels=LabEs, name="es" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
covCs		<-mxAlgebra( expression= cs %*% t(cs), name="Cs" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Cs+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("ac22","ac32","ac42","ac33","ac43","ac44"), name="ac" ) # Component paths for factors 2 and 3
PathsCcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("cc22","cc32","cc42","cc33","cc43","cc44"), name="cc" )
PathsEcsub	<-mxMatrix(type="Diag", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("ec22","ec33","ec44"), name="ec" )
PathsP11	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11", name="pc" ) # SD path for factor 3 (the PRS factor)
Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  #Padding

Ze31		<-mxMatrix(type="Zero",	nrow=3, ncol=1, free=F, name="Z31" )  #Padding
Ze13		<-mxMatrix(type="Zero",	nrow=1, ncol=3, free=F, name="Z13" )  #Padding
covAcsub	<-mxAlgebra( expression= ac %*% t(ac), name="Acsub" )
covCcsub	<-mxAlgebra( expression= cc %*% t(cc), name="Ccsub" )
covEcsub	<-mxAlgebra( expression= ec %*% t(ec), name="Ecsub" )
covPcsub	<-mxAlgebra( expression= Acsub+Ccsub+Ecsub, name="Vcsub" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)

covPc11	<-mxAlgebra( expression= pc %*% t(pc), name="Pc11" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components

covPc		<-mxAlgebra(cbind(rbind(Pc11,Z31), rbind(Z13,Vcsub)), name="Vc") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPcMz	<-mxAlgebra(cbind(rbind(Pc11,Z31) ,rbind(Z13,Acsub+Ccsub)), name="Vcmz") #I specify the MZ between-twin covariance - excluding E parameters
covPcDz	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z31), rbind(Z13,.5%x%Acsub+Ccsub)), name="Vcdz") #I specify the DZ between-twin covariance - specifying half of A and excluding E

# Generate Covariance of Latent factor model Including Causal Paths between factors
Id4		<-mxMatrix(type="Iden",	nrow=4, ncol=4, free=F, name="I4" )
covFVc	<-mxAlgebra( expression= solve(I4-PhC) %&% Vc, name ="FVc")
covFcMz	<-mxAlgebra( expression= solve(I4-PhC) %&% Vcmz, name ="Fcmz")
covFcDz	<-mxAlgebra( expression= solve(I4-PhC) %&% Vcdz, name ="Fcdz")

# Constraint on total variance of Ordinal variable (A+C+E=1)
varL1		<- mxConstraint( expression=FVc[2,2]==1, name="L1" )

FcovMZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcmz), cbind(Fcmz, FVc))) , name="expFCovMZ" )
FcovDZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcdz), cbind(Fcdz, FVc))) , name="expFCovDZ" )

SpcovMZ	<-mxAlgebra( expression= rbind (cbind(Vs, As+Cs), cbind(As+Cs, Vs)) , name="expSpCovMZ" )
SpcovDZ	<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As+Cs), cbind(.5%x%As+Cs, Vs)) , name="expSpCovDZ" )

TOTcovMZ	<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )

# *******************************************************************************************************
# Calculator

# Standardize the causal effects
Stcp1on2	<-mxAlgebra( expression= (PhC[2,1]* sqrt(FVc[1,1]))/sqrt(FVc[2,2]) , name="Stand_1on2" )
Stcp1on3	<-mxAlgebra( expression= (PhC[3,1]* sqrt(FVc[1,1]))/sqrt(FVc[3,3]) , name="Stand_1on3" )
Stcp1on4	<-mxAlgebra( expression= (PhC[4,1]* sqrt(FVc[1,1]))/sqrt(FVc[4,4]) , name="Stand_1on4" )
Stcp2on3	<-mxAlgebra( expression= (PhC[3,2]* sqrt(FVc[2,2]))/sqrt(FVc[3,3]) , name="Stand_2on3" )
Stcp2on4	<-mxAlgebra( expression= (PhC[4,2]* sqrt(FVc[2,2]))/sqrt(FVc[4,4]) , name="Stand_2on4" )
Stcp3on4	<-mxAlgebra( expression= (PhC[4,3]* sqrt(FVc[3,3]))/sqrt(FVc[4,4]) , name="Stand_3on4" )

# Standardize the Total var/covariances matrices of the observed variables
Id10		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I10" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )

# Phenotypic, A, C and E correlations	
RfactAc	<-mxAlgebra( expression= solve(sqrt(I3*Acsub)) %&% Acsub, name="Ra" )
RfactCc	<-mxAlgebra( expression= solve(sqrt(I3*Ccsub)) %&% Ccsub, name="Rc" )
RfactEc	<-mxAlgebra( expression= solve(sqrt(I3*Ecsub)) %&% Ecsub, name="Re" )
RfactP	<-mxAlgebra( expression= solve(sqrt(I4*FVc)) %&% FVc, name="Rph" )

# Standardize the Common Effects
covFVc33	<-mxAlgebra( expression= FVc[2:4,2:4], name ="FVc33")
stcovAc	<-mxAlgebra( expression= Acsub/FVc33, name="stAc" )
stcovCc	<-mxAlgebra( expression= Ccsub/FVc33, name="stCc" )
stcovEc	<-mxAlgebra( expression= Ecsub/FVc33, name="stEc" )

# Standardised path estimates
StpathAc	<-mxAlgebra( expression= (sqrt(stAc)), name="stpac" )
StpathCc	<-mxAlgebra( expression= (sqrt(stCc)), name="stpcc" )
StpathEc	<-mxAlgebra( expression= (sqrt(stEc)), name="stpec" )

# Algebra to compute Rph-A, Rph-C and Rph-E
RphA12	<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[2,1]*sqrt(stAc[2,2]), name = 'Rpha12')
RphA13	<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[3,1]*sqrt(stAc[3,3]), name = 'Rpha13')
RphA23	<-mxAlgebra(expression=sqrt(stAc[2,2])*Ra[3,2]*sqrt(stAc[3,3]), name = 'Rpha23')
RphC12	<-mxAlgebra(expression=sqrt(stCc[1,1])*Rc[2,1]*sqrt(stCc[2,2]), name = 'Rphc12')
RphC13	<-mxAlgebra(expression=sqrt(stCc[1,1])*Rc[3,1]*sqrt(stCc[3,3]), name = 'Rphc13')
RphC23	<-mxAlgebra(expression=sqrt(stCc[2,2])*Rc[3,2]*sqrt(stCc[3,3]), name = 'Rphc23')
RphE12	<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[2,1]*sqrt(stEc[2,2]), name = 'Rphe12')
RphE13	<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[3,1]*sqrt(stEc[3,3]), name = 'Rphe13')
RphE23	<-mxAlgebra(expression=sqrt(stEc[2,2])*Re[3,2]*sqrt(stEc[3,3]), name = 'Rphe23')

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= sqrt(As/( (FactL %&% FVc) +Vs)), name="stAs" )
stcovCs	<-mxAlgebra( expression= sqrt(Cs/( (FactL %&% FVc) +Vs)), name="stCs" )
stcovEs	<-mxAlgebra( expression= sqrt(Es/( (FactL %&% FVc) +Vs)), name="stEs" )

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% FVc / TOTexpCovMZ[1:5,1:5])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1		<-list(Means,Load,LoadTw,PhCaus,PathsAs,PathsCs,PathsEs,covAs,covCs,covEs,covPs,Id2,Id3,Id4,Id10)
pars2		<-list(PathsAcsub,PathsCcsub,PathsEcsub,PathsP11,Ze31,Ze13,Ze11,covAcsub,covCcsub,covEcsub,covPcsub,covPc11,covPc,covPcMz,covPcDz,covFVc,covFVc33,covFcMz,covFcDz)
pars3		<-list(obsage1, obsage2, obssex1, obssex2, betaA, betaS, Thr, inc, Thres)
parsst	<-list(stcovAs, stcovCs, stcovEs, stcovAc, stcovCc, stcovEc, RfactAc, RfactCc, RfactEc, RfactP,RphA12,RphA13,RphA23,RphC12,RphC13,RphC23,RphE12,RphE13,RphE23,StpathAc,StpathCc,StpathEc)
parsmed	<-list(Stcp1on2, Stcp1on3, Stcp1on4, Stcp2on3, Stcp2on4, Stcp3on4)
modelMZ	<-mxModel(pars1, pars2, pars3, parsmed, FcovMZ, SpcovMZ, TOTcovMZ, dataMZ, objMZ, Rfactmz, parsst, fitFunction, StFL, varL1, name="MZ" )
modelDZ	<-mxModel(pars1, pars2, pars3, FcovDZ, SpcovDZ, TOTcovDZ, dataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )

minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cistFL	<-mxCI (c ('MZ.StandFact','MZ.Stand_1on2','MZ.Stand_1on3','MZ.Stand_1on4','MZ.Stand_2on3','MZ.Stand_2on4','MZ.Stand_3on4','MZ.PhC'))
cistVs	<-mxCI (c ('MZ.stAs[5,5]','MZ.stAs[4,4]',
				'MZ.stCs[5,5]','MZ.stCs[4,4]',
				'MZ.stEs[5,5]','MZ.stEs[4,4]') ) 	# standardized var comp from specific Factors
cistVc	<-mxCI (c ('MZ.stAc[1,1]','MZ.stAc[2,1]','MZ.stAc[3,1]','MZ.stAc[2,2]','MZ.stAc[3,2]','MZ.stAc[3,3]',
				'MZ.stCc[1,1]','MZ.stCc[2,1]','MZ.stCc[3,1]','MZ.stCc[2,2]','MZ.stCc[3,2]','MZ.stCc[3,3]',
				'MZ.stEc[1,1]','MZ.stEc[2,2]','MZ.stEc[3,3]') ) 	# standardized var comp for ACE on latent Factors
cistRc	<-mxCI (c ('MZ.Rpha12','MZ.Rpha13','MZ.Rpha23','MZ.Rphc12','MZ.Rphc13','MZ.Rphc23','MZ.Rphe12','MZ.Rphe13','MZ.Rphe23','MZ.Ra','MZ.Rc','MZ.Re','MZ.stpac','MZ.stpcc','MZ.stpec') ) 	
ACEMMs1Model	<-mxModel("aceMMs1", pars1, pars2, modelMZ, modelDZ, minus2ll, obj, cistFL, cistVs, cistVc, cistRc) 

# --------------------------------------------------------------------------------------------------------------------------------
# 5a RUN ACEMs Factor Model with phenotypic causal mediation paths by Zygosity

ACEMMs1Fit		<-mxTryHardOrdinal(ACEMMs1Model, intervals=F)
(ACEMMs1Summ	<-summary(ACEMMs1Fit, verbose=F))

mxEval(MZ.Acsub, ACEMMs1Fit)
mxEval(MZ.Ccsub, ACEMMs1Fit)
mxEval(MZ.Ecsub, ACEMMs1Fit)
mxEval(MZ.Vcsub, ACEMMs1Fit)
mxEval(MZ.FVc33, ACEMMs1Fit)

mxEval(MZ.Stand_1on2[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_1on3[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_1on4[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_2on3[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_2on4[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_3on4[1,1], ACEMMs1Fit)


#------------------------------------------------------------------------
# Submodel 5aii: Drop correlation path c paths from previous model 
#------------------------------------------------------------------------
# Drop correlation c paths from the model
# -----------------------------------------------------------------------
AEMMs1aModel	<- mxModel(ACEMMs1Fit, name="AEMMs1a")
AEMMs1aModel	<- omxSetParameters(AEMMs1aModel, labels=c('cc22','cc32','cc42','cc33','cc43','cc44'), free=FALSE, values=0)
AEMMs1aModel	<- omxSetParameters(AEMMs1aModel, labels=c('cs4'), free=FALSE, values=0)
AEMMs1aFit		<- mxTryHardOrdinal(AEMMs1aModel, intervals=T)
(AEMMs1aSum		<- summary(AEMMs1aFit))

mxCompare(ACEMMs1Fit, AEMMs1aFit)

mxEval(MZ.FactcorMZ, AEMMs1aFit)
mxEval(DZ.FactcorDZ, AEMMs1aFit)

mxEval(MZ.Acsub, AEMMs1aFit)
mxEval(MZ.Ccsub, AEMMs1aFit)
mxEval(MZ.Ecsub, AEMMs1aFit)
mxEval(MZ.Vcsub, AEMMs1aFit)
mxEval(MZ.FVc22, AEMMs1aFit)

mxEval(MZ.FVc, AEMMs1aFit)

mxEval(MZ.stAc, AEMMs1aFit)
mxEval(MZ.stCc, AEMMs1aFit)
mxEval(MZ.stEc, AEMMs1aFit)

mxEval(MZ.stpac, AEMMs1aFit)
mxEval(MZ.stpcc, AEMMs1aFit)
mxEval(MZ.stpec, AEMMs1aFit)

mxEval(MZ.Ra, AEMMs1aFit)
mxEval(MZ.Rc, AEMMs1aFit)
mxEval(MZ.Re, AEMMs1aFit)
mxEval(MZ.Rph, AEMMs1aFit)

mxEval(MZ.Rpha, AEMMs1aFit) 
mxEval(MZ.Rphc, AEMMs1aFit) 
mxEval(MZ.Rphe, AEMMs1aFit) 

mxEval(MZ.StandFact, AEMMs1aFit)
mxEval(MZ.Stand_1on2[1,1], AEMMs1aFit)
mxEval(MZ.Stand_1on3[1,1], AEMMs1aFit)
mxEval(MZ.Stand_1on4[1,1], AEMMs1aFit)
mxEval(MZ.Stand_2on3[1,1], AEMMs1aFit)
mxEval(MZ.Stand_2on4[1,1], AEMMs1aFit)
mxEval(MZ.Stand_4on4[1,1], AEMMs1aFit)


# 
#****************************************************************************************************************************
# __(VIIb)_____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL for PD>SO with mediation by Vict
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Phenotypic Factors: F1>F2>F3 & F1>F3;
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
#_____________________________________________________________________________________________________________________________

nv		<- 6				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 4				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nvo 		<- 1     			# number of ordinal variables per twin
nvc 		<- nv-nvo  			# number of continuous variables per twin
poso 		<- nvo 			# position where ordinal variables start
nth		<- 4				# number of max thresholds
ninc 		<- nth-1 			# number of max increments
ncovariates <- 2 				# number of covariates
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups	<- c("mz", "dz")
Vars		<- c('PRSAnx','PRSDep','Dep','Anx','Vict','SO')
selVars	<- c('PRSAnx1','PRSDep1','Dep1','Anx1','Vict1','SO1',
		     'PRSAnx2','PRSDep2','Dep2','Anx2','Vict2','SO2')
useVars	<- c('PRSAnx1','PRSDep1','Dep1','Anx1','Vict1','SO1',
		     'PRSAnx2','PRSDep2','Dep2','Anx2','Vict2','SO2','age1','sex1','age2','sex2')

mzData		<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , useVars)
dzData		<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , useVars)

psych::describe(mzData)
psych::describe(dzData)

mzData$SO1[is.na(mzData$age1)] <- NA
mzData$SO2[is.na(mzData$age2)] <- NA
dzData$SO1[is.na(dzData$age1)] <- NA
dzData$SO2[is.na(dzData$age2)] <- NA

mzData$SO1[is.na(mzData$sex1)] <- NA
mzData$SO2[is.na(mzData$sex2)] <- NA
dzData$SO1[is.na(dzData$sex1)] <- NA
dzData$SO2[is.na(dzData$sex2)] <- NA

mzData$age1[is.na(mzData$age1)] <- 999
mzData$age2[is.na(mzData$age2)] <- 999
dzData$age1[is.na(dzData$age1)] <- 999
dzData$age2[is.na(dzData$age2)] <- 999

mzData$sex1[is.na(mzData$sex1)] <- 999
mzData$sex2[is.na(mzData$sex2)] <- 999
dzData$sex1[is.na(dzData$sex1)] <- 999
dzData$sex2[is.na(dzData$sex2)] <- 999

psych::describe(mzData)
psych::describe(dzData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- c(5,5,4.5,4.1,2.5,0))
(PatM		<- c(TRUE,T,TRUE,TRUE,TRUE,F))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEs	<- c('es1','es1','es3','es3','es5','es6'))
(LabAs	<- c('as1','as1','as3','as3','as5','as6'))
(LabCs	<- c('cs1','cs1','cs3','cs3','cs5','cs6'))

PatSpe	<- c(T,T,T,T,F,F)
PatSpac	<- c(F,F,T,T,F,F)
StSpa		<- c(0,0,.5,.5,0,0)
StSpc		<- c(0,0,.5,.5,0,0)
StSpe		<- c(.5,.5,.5,.5,0,0)

# all 1st loadings fixed to 1
PatFl		<- c(F,T,F,F,F,F,			
		     F,F,F,T,F,F,
		     F,F,F,F,F,F,
		     F,F,F,F,F,F)

StFl		<- c(1,.5,0,0,0,0,
		     0,0,1,.5,0,0,
		     0,0,0,0,1,0,
		     0,0,0,0,0,1)

LabFl		<- c('l1','l2',NA,NA,NA,NA,
	 	     NA,NA,'l3','l4',NA,NA,
	 	     NA,NA,NA,NA,'l5',NA,
	 	     NA,NA,NA,NA,NA,'l6')

PatPhC	<- c(F,T,T,T,
		     F,F,T,T,
		     F,F,F,T,
		     F,F,F,F)

StPhC		<- c(0,.1,.1,.1,
		     0,0,.1,.1,
		     0,0,0,.1,
		     0,0,0,0)

LabPhC	<- c(NA,'c1on2','c1on3','c1on4',
		     NA,NA,'c2on3','c2on4',
		     NA,NA,NA,'c3on4',
		     NA,NA,NA,NA)	 

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmean,Stmean), labels=c(mLabs,mLabs), name="expMean") 

# Threshold and covariates
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')

ThPat		<-c(T,T,T,T)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.3, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.4, labels=LabCovS, name="BsexTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2), name="expThres")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 4 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )

# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpa, labels=LabAs, name="as" )
PathsCs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpc, labels=LabCs, name="cs" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpe, labels=LabEs, name="es" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
covCs		<-mxAlgebra( expression= cs %*% t(cs), name="Cs" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Cs+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("ac22","ac32","ac42","ac33","ac43","ac44"), name="ac" ) # Component paths for factors 2 and 3
PathsCcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("cc22","cc32","cc42","cc33","cc43","cc44"), name="cc" )
PathsEcsub	<-mxMatrix(type="Diag", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("ec22","ec33","ec44"), name="ec" )
PathsP11	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11", name="pc" ) # SD path for factor 3 (the PRS factor)
Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  #Padding

Ze31		<-mxMatrix(type="Zero",	nrow=3, ncol=1, free=F, name="Z31" )  #Padding
Ze13		<-mxMatrix(type="Zero",	nrow=1, ncol=3, free=F, name="Z13" )  #Padding
covAcsub	<-mxAlgebra( expression= ac %*% t(ac), name="Acsub" )
covCcsub	<-mxAlgebra( expression= cc %*% t(cc), name="Ccsub" )
covEcsub	<-mxAlgebra( expression= ec %*% t(ec), name="Ecsub" )
covPcsub	<-mxAlgebra( expression= Acsub+Ccsub+Ecsub, name="Vcsub" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)

covPc11	<-mxAlgebra( expression= pc %*% t(pc), name="Pc11" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components

covPc		<-mxAlgebra(cbind(rbind(Pc11,Z31), rbind(Z13,Vcsub)), name="Vc") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPcMz	<-mxAlgebra(cbind(rbind(Pc11,Z31) ,rbind(Z13,Acsub+Ccsub)), name="Vcmz") #I specify the MZ between-twin covariance - excluding E parameters
covPcDz	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z31), rbind(Z13,.5%x%Acsub+Ccsub)), name="Vcdz") #I specify the DZ between-twin covariance - specifying half of A and excluding E

# Generate Covariance of Latent factor model Including Causal Paths between factors
Id4		<-mxMatrix(type="Iden",	nrow=4, ncol=4, free=F, name="I4" )
covFVc	<-mxAlgebra( expression= solve(I4-PhC) %&% Vc, name ="FVc")
covFcMz	<-mxAlgebra( expression= solve(I4-PhC) %&% Vcmz, name ="Fcmz")
covFcDz	<-mxAlgebra( expression= solve(I4-PhC) %&% Vcdz, name ="Fcdz")

# Constraint on total variance of Ordinal variable (A+C+E=1)
varL1		<- mxConstraint( expression=FVc[4,4]==1, name="L1" )

FcovMZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcmz), cbind(Fcmz, FVc))) , name="expFCovMZ" )
FcovDZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcdz), cbind(Fcdz, FVc))) , name="expFCovDZ" )

SpcovMZ	<-mxAlgebra( expression= rbind (cbind(Vs, As+Cs), cbind(As+Cs, Vs)) , name="expSpCovMZ" )
SpcovDZ	<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As+Cs), cbind(.5%x%As+Cs, Vs)) , name="expSpCovDZ" )

TOTcovMZ	<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )

# *******************************************************************************************************
# Calculator

# Standardize the causal effects
Stcp1on2	<-mxAlgebra( expression= (PhC[2,1]* sqrt(FVc[1,1]))/sqrt(FVc[2,2]) , name="Stand_1on2" )
Stcp1on3	<-mxAlgebra( expression= (PhC[3,1]* sqrt(FVc[1,1]))/sqrt(FVc[3,3]) , name="Stand_1on3" )
Stcp1on4	<-mxAlgebra( expression= (PhC[4,1]* sqrt(FVc[1,1]))/sqrt(FVc[4,4]) , name="Stand_1on4" )
Stcp2on3	<-mxAlgebra( expression= (PhC[3,2]* sqrt(FVc[2,2]))/sqrt(FVc[3,3]) , name="Stand_2on3" )
Stcp2on4	<-mxAlgebra( expression= (PhC[4,2]* sqrt(FVc[2,2]))/sqrt(FVc[4,4]) , name="Stand_2on4" )
Stcp3on4	<-mxAlgebra( expression= (PhC[4,3]* sqrt(FVc[3,3]))/sqrt(FVc[4,4]) , name="Stand_3on4" )

# Standardize the Total var/covariances matrices of the observed variables
Id12		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I12" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I12*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I12*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )

# Phenotypic, A, C and E correlations	
Id3		<-mxMatrix(type="Iden",	nrow=3, ncol=3, free=F, name="I3" )
RfactAc	<-mxAlgebra( expression= solve(sqrt(I3*Acsub)) %&% Acsub, name="Ra" )
RfactCc	<-mxAlgebra( expression= solve(sqrt(I3*Ccsub)) %&% Ccsub, name="Rc" )
RfactEc	<-mxAlgebra( expression= solve(sqrt(I3*Ecsub)) %&% Ecsub, name="Re" )
RfactP	<-mxAlgebra( expression= solve(sqrt(I4*FVc)) %&% FVc, name="Rph" )

# Standardize the Common Effects
covFVc33	<-mxAlgebra( expression= FVc[2:4,2:4], name ="FVc33")
stcovAc	<-mxAlgebra( expression= Acsub/FVc33, name="stAc" )
stcovCc	<-mxAlgebra( expression= Ccsub/FVc33, name="stCc" )
stcovEc	<-mxAlgebra( expression= Ecsub/FVc33, name="stEc" )

# Standardised path estimates
StpathAc	<-mxAlgebra( expression= (sqrt(stAc)), name="stpac" )
StpathCc	<-mxAlgebra( expression= (sqrt(stCc)), name="stpcc" )
StpathEc	<-mxAlgebra( expression= (sqrt(stEc)), name="stpec" )

# Algebra to compute Rph-A, Rph-C and Rph-E
RphA12	<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[2,1]*sqrt(stAc[2,2]), name = 'Rpha12')
RphA13	<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[3,1]*sqrt(stAc[3,3]), name = 'Rpha13')
RphA23	<-mxAlgebra(expression=sqrt(stAc[2,2])*Ra[3,2]*sqrt(stAc[3,3]), name = 'Rpha23')
RphC12	<-mxAlgebra(expression=sqrt(stCc[1,1])*Rc[2,1]*sqrt(stCc[2,2]), name = 'Rphc12')
RphC13	<-mxAlgebra(expression=sqrt(stCc[1,1])*Rc[3,1]*sqrt(stCc[3,3]), name = 'Rphc13')
RphC23	<-mxAlgebra(expression=sqrt(stCc[2,2])*Rc[3,2]*sqrt(stCc[3,3]), name = 'Rphc23')
RphE12	<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[2,1]*sqrt(stEc[2,2]), name = 'Rphe12')
RphE13	<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[3,1]*sqrt(stEc[3,3]), name = 'Rphe13')
RphE23	<-mxAlgebra(expression=sqrt(stEc[2,2])*Re[3,2]*sqrt(stEc[3,3]), name = 'Rphe23')

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= sqrt(As/( (FactL %&% FVc) +Vs)), name="stAs" )
stcovCs	<-mxAlgebra( expression= sqrt(Cs/( (FactL %&% FVc) +Vs)), name="stCs" )
stcovEs	<-mxAlgebra( expression= sqrt(Es/( (FactL %&% FVc) +Vs)), name="stEs" )

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% FVc / TOTexpCovMZ[1:6,1:6])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1		<-list(Means,Load,LoadTw,PhCaus,PathsAs,PathsCs,PathsEs,covAs,covCs,covEs,covPs,Id2,Id3,Id4,Id12)
pars2		<-list(PathsAcsub,PathsCcsub,PathsEcsub,PathsP11,Ze31,Ze13,Ze11,covAcsub,covCcsub,covEcsub,covPcsub,covPc11,covPc,covPcMz,covPcDz,covFVc,covFVc33,covFcMz,covFcDz)
pars3		<-list(obsage1, obsage2, obssex1, obssex2, betaA, betaS, Thr, inc, Thres)
parsst	<-list(stcovAs, stcovCs, stcovEs, stcovAc, stcovCc, stcovEc, RfactAc, RfactCc, RfactEc, RfactP,RphA12,RphA13,RphA23,RphC12,RphC13,RphC23,RphE12,RphE13,RphE23,StpathAc,StpathCc,StpathEc)
parsmed	<-list(Stcp1on2, Stcp1on3, Stcp1on4, Stcp2on3, Stcp2on4, Stcp3on4)
modelMZ	<-mxModel(pars1, pars2, pars3, parsmed, FcovMZ, SpcovMZ, TOTcovMZ, dataMZ, objMZ, Rfactmz, parsst, fitFunction, StFL, varL1, name="MZ" )
modelDZ	<-mxModel(pars1, pars2, pars3, FcovDZ, SpcovDZ, TOTcovDZ, dataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )

minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cistFL	<-mxCI (c ('MZ.StandFact','MZ.Stand_1on2','MZ.Stand_1on3','MZ.Stand_1on4','MZ.Stand_2on3','MZ.Stand_2on4','MZ.Stand_3on4','MZ.PhC'))
cistVs	<-mxCI (c ('MZ.stAs[3,3]','MZ.stAs[4,4]',
				'MZ.stCs[3,3]','MZ.stCs[4,4]',
				'MZ.stEs[1,1]','MZ.stEs[2,2]','MZ.stEs[3,3]','MZ.stEs[4,4]') ) 	# standardized var comp from specific Factors
cistVc	<-mxCI (c ('MZ.stAc[1,1]','MZ.stAc[2,1]','MZ.stAc[3,1]','MZ.stAc[2,2]','MZ.stAc[3,2]','MZ.stAc[3,3]',
				'MZ.stCc[1,1]','MZ.stCc[2,1]','MZ.stCc[3,1]','MZ.stCc[2,2]','MZ.stCc[3,2]','MZ.stCc[3,3]',
				'MZ.stEc[1,1]','MZ.stEc[2,2]','MZ.stEc[3,3]') ) 	# standardized var comp for ACE on latent Factors
cistRc	<-mxCI (c ('MZ.Rpha12','MZ.Rpha13','MZ.Rpha23','MZ.Rphc12','MZ.Rphc13','MZ.Rphc23','MZ.Rphe12','MZ.Rphe13','MZ.Rphe23','MZ.Ra','MZ.Rc','MZ.Re','MZ.stpac','MZ.stpcc','MZ.stpec') ) 	
ACEMMs2Model	<-mxModel("aceMMs2", pars1, pars2, modelMZ, modelDZ, minus2ll, obj, cistFL, cistVs, cistVc, cistRc) 

# --------------------------------------------------------------------------------------------------------------------------------
# 5a RUN ACEMs Factor Model with phenotypic causal mediation paths by Zygosity

ACEMMs2Fit		<-mxTryHardOrdinal(ACEMMs2Model, intervals=F)
(ACEMMs2Summ	<-summary(ACEMMs2Fit, verbose=F))

mxEval(MZ.Acsub, ACEMMs2Fit)
mxEval(MZ.Ccsub, ACEMMs2Fit)
mxEval(MZ.Ecsub, ACEMMs2Fit)
mxEval(MZ.Vcsub, ACEMMs2Fit)
mxEval(MZ.FVc33, ACEMMs2Fit)

mxEval(MZ.Stand_1on2[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_1on3[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_1on4[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_2on3[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_2on4[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_3on4[1,1], ACEMMs2Fit)

#------------------------------------------------------------------------
# Submodel 5aii: Drop correlation path c paths from previous model 
#------------------------------------------------------------------------
# Drop correlation c paths from the model
# -----------------------------------------------------------------------
AEMMs2aModel	<- mxModel(ACEMMs2Fit, name="AEMMs2a")
AEMMs2aModel	<- omxSetParameters(AEMMs2aModel, labels=c('cc22','cc32','cc42','cc33','cc43','cc44'), free=FALSE, values=0)
AEMMs2aModel	<- omxSetParameters(AEMMs2aModel, labels=c('cs3'), free=FALSE, values=0)
AEMMs2aFit		<- mxTryHardOrdinal(AEMMs2aModel, intervals=T)
(AEMMs2aSum		<- summary(AEMMs2aFit))

mxCompare(ACEMMs2Fit, AEMMs2aFit)

mxEval(MZ.FactcorMZ, AEMMs2aFit)
mxEval(DZ.FactcorDZ, AEMMs2aFit)

mxEval(MZ.Acsub, AEMMs2aFit)
mxEval(MZ.Ccsub, AEMMs2aFit)
mxEval(MZ.Ecsub, AEMMs2aFit)
mxEval(MZ.Vcsub, AEMMs2aFit)
mxEval(MZ.FVc33, AEMMs2aFit)

mxEval(MZ.FVc, AEMMs2aFit)

mxEval(MZ.stAc, AEMMs2aFit)
mxEval(MZ.stCc, AEMMs2aFit)
mxEval(MZ.stEc, AEMMs2aFit)

mxEval(MZ.stpac, AEMMs2aFit)
mxEval(MZ.stpcc, AEMMs2aFit)
mxEval(MZ.stpec, AEMMs2aFit)

mxEval(MZ.Ra, AEMMs2aFit)
mxEval(MZ.Rc, AEMMs2aFit)
mxEval(MZ.Re, AEMMs2aFit)
mxEval(MZ.Rph, AEMMs2aFit)

mxEval(MZ.Rpha12, AEMMs2aFit) 
mxEval(MZ.Rpha13, AEMMs2aFit) 
mxEval(MZ.Rpha23, AEMMs2aFit) 

mxEval(MZ.StandFact, AEMMs2aFit)
mxEval(MZ.Stand_1on2[1,1], AEMMs2aFit)
mxEval(MZ.Stand_1on3[1,1], AEMMs2aFit)
mxEval(MZ.Stand_1on4[1,1], AEMMs2aFit)
mxEval(MZ.Stand_2on3[1,1], AEMMs2aFit)
mxEval(MZ.Stand_2on4[1,1], AEMMs2aFit)
mxEval(MZ.Stand_3on4[1,1], AEMMs2aFit)



# 
#****************************************************************************************************************************
# __(VIIIa)_____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL for SO>VICT>PD with mediation by Vict and confounding by CGN
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Phenotypic Factors: F1>F2>F3 & F1>F3;
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
#_____________________________________________________________________________________________________________________________

nv		<- 6				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 5				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nvo 		<- 1     			# number of ordinal variables per twin
nvc 		<- nv-nvo  			# number of continuous variables per twin
poso 		<- nvo 			# position where ordinal variables start
nth		<- 4				# number of max thresholds
ninc 		<- nth-1 			# number of max increments
ncovariates <- 2 				# number of covariates
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups	<- c("mz", "dz")
Vars		<- c('PRSSO','SO','Vict','Dep','Anx','CGN4')
selVars	<- c('PRSSO1','SO1','Vict1','Dep1','Anx1','CGN41',
		     'PRSSO2','SO2','Vict2','Dep2','Anx2','CGN42')
useVars	<- c('PRSSO1','SO1','Vict1','Dep1','Anx1','CGN41',
		     'PRSSO2','SO2','Vict2','Dep2','Anx2','CGN42','age1','sex1','age2','sex2')

mzData		<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , useVars)
dzData		<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , useVars)

psych::describe(mzData)
psych::describe(dzData)

mzData$SO1[is.na(mzData$age1)] <- NA
mzData$SO2[is.na(mzData$age2)] <- NA
dzData$SO1[is.na(dzData$age1)] <- NA
dzData$SO2[is.na(dzData$age2)] <- NA

mzData$SO1[is.na(mzData$sex1)] <- NA
mzData$SO2[is.na(mzData$sex2)] <- NA
dzData$SO1[is.na(dzData$sex1)] <- NA
dzData$SO2[is.na(dzData$sex2)] <- NA

mzData$age1[is.na(mzData$age1)] <- 999
mzData$age2[is.na(mzData$age2)] <- 999
dzData$age1[is.na(dzData$age1)] <- 999
dzData$age2[is.na(dzData$age2)] <- 999

mzData$sex1[is.na(mzData$sex1)] <- 999
mzData$sex2[is.na(mzData$sex2)] <- 999
dzData$sex1[is.na(dzData$sex1)] <- 999
dzData$sex2[is.na(dzData$sex2)] <- 999

psych::describe(mzData)
psych::describe(dzData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- c(5,0,2.5,4.5,4.1,4))
(PatM		<- c(TRUE,F,TRUE,TRUE,TRUE,TRUE))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEs	<- c('es1','es2','es3','es4','es4','es6'))
(LabAs	<- c('as1','as2','as3','as4','as4','as6'))
#(LabCs	<- c('cs1','cs2','cs3','cs4','cs4','cs6'))

PatSpe	<- c(F,F,F,TRUE,TRUE,F)
PatSpac	<- c(F,F,F,TRUE,TRUE,F)
StSpa		<- c(0,0,0,.5,.5,0)
#StSpc		<- c(0,0,0,.5,.5,0)
StSpe		<- c(0,0,0,.5,.5,0)

# all 1st loadings fixed to 1
PatFl		<- c(F,F,F,F,F,F,			
		     F,F,F,F,F,F,
		     F,F,F,F,F,F,
		     F,F,F,F,T,F,
		     F,F,F,F,F,F)

StFl		<- c(1,0,0,0,0,0,
		     0,1,0,0,0,0,
		     0,0,1,0,0,0,
		     0,0,0,1,.5,0,
		     0,0,0,0,0,1)

LabFl		<- c('l1',NA,NA,NA,NA,NA,
	 	     NA,'l2',NA,NA,NA,NA,
	 	     NA,NA,'l3',NA,NA,NA,
	 	     NA,NA,NA,'l4','l5',NA,
	 	     NA,NA,NA,NA,NA,'l6')

PatPhC	<- c(F,T,T,T,T,
		     F,F,T,T,F,
		     F,F,F,T,F,
		     F,F,F,F,F,
		     F,T,T,T,F)

StPhC		<- c(0,.3,.3,.3,.1,
		     0,0,.3,.3,0,
		     0,0,0,.3,0,
		     0,0,0,0,0,
		     0,.1,.1,.3,0)

LabPhC	<- c(NA,'c1on2','c1on3','c1on4','c1on5',
		     NA,NA,'c2on3','c2on4',NA,
		     NA,NA,NA,'c3on4',NA,
		     NA,NA,NA,NA,NA,
		     NA,'c5on2','c5on3','c5on4',NA)	 

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmean,Stmean), labels=c(mLabs,mLabs), name="expMean") 

# Threshold and covariates
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')

ThPat		<-c(T,T,T,T)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.3, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.4, labels=LabCovS, name="BsexTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2), name="expThres")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )

# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpa, labels=LabAs, name="as" )
#PathsCs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpc, labels=LabCs, name="cs" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpe, labels=LabEs, name="es" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
#covCs		<-mxAlgebra( expression= cs %*% t(cs), name="Cs" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("ac22","ac32","ac42","ac52","ac33","ac43","ac53","ac44","ac54","ac55"), name="ac" ) # Component paths for factors 2 and 3
#PathsCcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("cc22","cc32","cc42","cc52","cc33","cc43","cc53","cc44","cc54","cc55"), name="cc" )
PathsEcsub	<-mxMatrix(type="Diag", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("ec22","ec33","ec44","ec55"), name="ec" )
PathsP11	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11", name="pc" ) # SD path for factor 3 (the PRS factor)
Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  #Padding

Ze41		<-mxMatrix(type="Zero",	nrow=4, ncol=1, free=F, name="Z41" )  #Padding
Ze14		<-mxMatrix(type="Zero",	nrow=1, ncol=4, free=F, name="Z14" )  #Padding
covAcsub	<-mxAlgebra( expression= ac %*% t(ac), name="Acsub" )
#covCcsub	<-mxAlgebra( expression= cc %*% t(cc), name="Ccsub" )
covEcsub	<-mxAlgebra( expression= ec %*% t(ec), name="Ecsub" )
covPcsub	<-mxAlgebra( expression= Acsub+Ecsub, name="Vcsub" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)

covPc11	<-mxAlgebra( expression= pc %*% t(pc), name="Pc11" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components

covPc		<-mxAlgebra(cbind(rbind(Pc11,Z41), rbind(Z14,Vcsub)), name="Vc") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPcMz	<-mxAlgebra(cbind(rbind(Pc11,Z41) ,rbind(Z14,Acsub)), name="Vcmz") #I specify the MZ between-twin covariance - excluding E parameters
covPcDz	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z41), rbind(Z14,.5%x%Acsub)), name="Vcdz") #I specify the DZ between-twin covariance - specifying half of A and excluding E

# Generate Covariance of Latent factor model Including Causal Paths between factors
Id5		<-mxMatrix(type="Iden",	nrow=5, ncol=5, free=F, name="I5" )
covFVc	<-mxAlgebra( expression= solve(I5-PhC) %&% Vc, name ="FVc")
covFcMz	<-mxAlgebra( expression= solve(I5-PhC) %&% Vcmz, name ="Fcmz")
covFcDz	<-mxAlgebra( expression= solve(I5-PhC) %&% Vcdz, name ="Fcdz")

# Constraint on total variance of Ordinal variable (A+C+E=1)
varL1		<- mxConstraint( expression=FVc[2,2]==1, name="L1" )

FcovMZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcmz), cbind(Fcmz, FVc))) , name="expFCovMZ" )
FcovDZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcdz), cbind(Fcdz, FVc))) , name="expFCovDZ" )

SpcovMZ	<-mxAlgebra( expression= rbind (cbind(Vs, As), cbind(As, Vs)) , name="expSpCovMZ" )
SpcovDZ	<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As), cbind(.5%x%As, Vs)) , name="expSpCovDZ" )

TOTcovMZ	<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )

# *******************************************************************************************************
# Calculator

# Standardize the causal effects
Stcp1on2	<-mxAlgebra( expression= (PhC[2,1]* sqrt(FVc[1,1]))/sqrt(FVc[2,2]) , name="Stand_1on2" )
Stcp1on3	<-mxAlgebra( expression= (PhC[3,1]* sqrt(FVc[1,1]))/sqrt(FVc[3,3]) , name="Stand_1on3" )
Stcp1on4	<-mxAlgebra( expression= (PhC[4,1]* sqrt(FVc[1,1]))/sqrt(FVc[4,4]) , name="Stand_1on4" )
Stcp1on5	<-mxAlgebra( expression= (PhC[5,1]* sqrt(FVc[1,1]))/sqrt(FVc[5,5]) , name="Stand_1on5" )
Stcp2on3	<-mxAlgebra( expression= (PhC[3,2]* sqrt(FVc[2,2]))/sqrt(FVc[3,3]) , name="Stand_2on3" )
Stcp2on4	<-mxAlgebra( expression= (PhC[4,2]* sqrt(FVc[2,2]))/sqrt(FVc[4,4]) , name="Stand_2on4" )
Stcp3on4	<-mxAlgebra( expression= (PhC[4,3]* sqrt(FVc[3,3]))/sqrt(FVc[4,4]) , name="Stand_3on4" )
Stcp5on2	<-mxAlgebra( expression= (PhC[2,5]* sqrt(FVc[5,5]))/sqrt(FVc[2,2]) , name="Stand_5on2" )
Stcp5on3	<-mxAlgebra( expression= (PhC[3,5]* sqrt(FVc[5,5]))/sqrt(FVc[3,3]) , name="Stand_5on3" )
Stcp5on4	<-mxAlgebra( expression= (PhC[4,5]* sqrt(FVc[5,5]))/sqrt(FVc[4,4]) , name="Stand_5on4" )

# Standardize the Total var/covariances matrices of the observed variables
Id12		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I12" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I12*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I12*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )

# Phenotypic, A, C and E correlations	
Id4		<-mxMatrix(type="Iden",	nrow=4, ncol=4, free=F, name="I4" )
RfactAc	<-mxAlgebra( expression= solve(sqrt(I4*Acsub)) %&% Acsub, name="Ra" )
#RfactCc	<-mxAlgebra( expression= solve(sqrt(I4*Ccsub)) %&% Ccsub, name="Rc" )
RfactEc	<-mxAlgebra( expression= solve(sqrt(I4*Ecsub)) %&% Ecsub, name="Re" )
RfactP	<-mxAlgebra( expression= solve(sqrt(I5*FVc)) %&% FVc, name="Rph" )

# Standardize the Common Effects
covFVc44	<-mxAlgebra( expression= FVc[2:5,2:5], name ="FVc44")
stcovAc	<-mxAlgebra( expression= Acsub/FVc44, name="stAc" )
#stcovCc	<-mxAlgebra( expression= Ccsub/FVc44, name="stCc" )
stcovEc	<-mxAlgebra( expression= Ecsub/FVc44, name="stEc" )

# Standardised path estimates
StpathAc	<-mxAlgebra( expression= (sqrt(stAc)), name="stpac" )
#StpathCc	<-mxAlgebra( expression= (sqrt(stCc)), name="stpcc" )
StpathEc	<-mxAlgebra( expression= (sqrt(stEc)), name="stpec" )

# Algebra to compute Rph-A, Rph-C and Rph-E
RphA12	<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[2,1]*sqrt(stAc[2,2]), name = 'Rpha12')
RphA13	<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[3,1]*sqrt(stAc[3,3]), name = 'Rpha13')
RphA14	<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[4,1]*sqrt(stAc[4,4]), name = 'Rpha14')
RphA23	<-mxAlgebra(expression=sqrt(stAc[2,2])*Ra[3,2]*sqrt(stAc[3,3]), name = 'Rpha23')
RphA24	<-mxAlgebra(expression=sqrt(stAc[2,2])*Ra[4,2]*sqrt(stAc[4,4]), name = 'Rpha24')
RphA34	<-mxAlgebra(expression=sqrt(stAc[3,3])*Ra[4,3]*sqrt(stAc[4,4]), name = 'Rpha34')
#RphC12	<-mxAlgebra(expression=sqrt(stCc[1,1])*Rc[2,1]*sqrt(stCc[2,2]), name = 'Rphc12')
#RphC13	<-mxAlgebra(expression=sqrt(stCc[1,1])*Rc[3,1]*sqrt(stCc[3,3]), name = 'Rphc13')
#RphC23	<-mxAlgebra(expression=sqrt(stCc[2,2])*Rc[3,2]*sqrt(stCc[3,3]), name = 'Rphc23')
RphE12	<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[2,1]*sqrt(stEc[2,2]), name = 'Rphe12')
RphE13	<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[3,1]*sqrt(stEc[3,3]), name = 'Rphe13')
RphE14	<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[4,1]*sqrt(stEc[4,4]), name = 'Rphe14')
RphE23	<-mxAlgebra(expression=sqrt(stEc[2,2])*Re[3,2]*sqrt(stEc[3,3]), name = 'Rphe23')
RphE24	<-mxAlgebra(expression=sqrt(stEc[2,2])*Re[4,2]*sqrt(stEc[4,4]), name = 'Rphe24')
RphE34	<-mxAlgebra(expression=sqrt(stEc[3,3])*Re[4,3]*sqrt(stEc[4,4]), name = 'Rphe34')

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= sqrt(As/( (FactL %&% FVc) +Vs)), name="stAs" )
#stcovCs	<-mxAlgebra( expression= sqrt(Cs/( (FactL %&% FVc) +Vs)), name="stCs" )
stcovEs	<-mxAlgebra( expression= sqrt(Es/( (FactL %&% FVc) +Vs)), name="stEs" )
# Standardized Factor Loadings

StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% FVc / TOTexpCovMZ[1:6,1:6])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1		<-list(Means,Load,LoadTw,PhCaus,PathsAs,PathsEs,covAs,covEs,covPs,Id2,Id4,Id5,Id12)
pars2		<-list(PathsAcsub,PathsEcsub,PathsP11,Ze41,Ze14,covAcsub,covEcsub,covPcsub,covPc11,covPc,covPcMz,covPcDz,covFVc,covFVc44,covFcMz,covFcDz)
pars3		<-list(obsage1, obsage2, obssex1, obssex2, betaA, betaS, Thr, inc, Thres)
parsst	<-list(stcovAs, stcovEs, stcovAc, stcovEc, RfactAc, RfactEc, RfactP,RphA12,RphA13,RphA14,RphA23,RphA24,RphA34,RphE12,RphE13,RphE14,RphE23,RphE24,RphE34,StpathAc,StpathEc)
parsmed	<-list(Stcp1on2, Stcp1on3, Stcp1on4, Stcp1on5, Stcp2on3, Stcp2on4, Stcp3on4, Stcp5on2, Stcp5on3, Stcp5on4)
modelMZ	<-mxModel(pars1, pars2, pars3, parsmed, FcovMZ, SpcovMZ, TOTcovMZ, dataMZ, objMZ, Rfactmz, parsst, fitFunction, StFL, varL1, name="MZ" )
modelDZ	<-mxModel(pars1, pars2, pars3, FcovDZ, SpcovDZ, TOTcovDZ, dataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )

minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cistFL	<-mxCI (c ('MZ.StandFact','MZ.Stand_1on2','MZ.Stand_1on3','MZ.Stand_1on4','MZ.Stand_1on5','MZ.Stand_2on3','MZ.Stand_2on4','MZ.Stand_3on4','MZ.Stand_5on2','MZ.Stand_5on3','MZ.Stand_5on4','MZ.PhC'))
cistVs	<-mxCI (c ('MZ.stAs[3,3]','MZ.stAs[4,4]',
				'MZ.stEs[1,1]','MZ.stEs[2,2]','MZ.stEs[3,3]','MZ.stEs[4,4]') ) 	# standardized var comp from specific Factors
cistVc	<-mxCI (c ('MZ.stAc[1,1]','MZ.stAc[2,1]','MZ.stAc[3,1]','MZ.stAc[4,1]','MZ.stAc[2,2]','MZ.stAc[3,2]','MZ.stAc[4,2]','MZ.stAc[3,3]','MZ.stAc[4,3]','MZ.stAc[4,4]',
				'MZ.stEc[1,1]','MZ.stEc[2,2]','MZ.stEc[3,3]','MZ.stEc[4,4]') ) 	# standardized var comp for ACE on latent Factors
cistRc	<-mxCI (c ('MZ.Rpha12','MZ.Rpha13','MZ.Rpha14','MZ.Rpha23','MZ.Rpha24','MZ.Rpha34','MZ.Rphe12','MZ.Rphe13','MZ.Rphe14','MZ.Rphe23','MZ.Rphe24','MZ.Rphe34','MZ.Ra','MZ.Re','MZ.stpac','MZ.stpec') ) 	
ACEMMs1Model	<-mxModel("aceMMs1", pars1, pars2, modelMZ, modelDZ, minus2ll, obj, cistFL, cistVs, cistVc, cistRc) 

# --------------------------------------------------------------------------------------------------------------------------------
# 5a RUN ACEMs Factor Model with phenotypic causal mediation paths by Zygosity

ACEMMs1Fit		<-mxTryHardOrdinal(ACEMMs1Model, intervals=F)
(ACEMMs1Summ	<-summary(ACEMMs1Fit, verbose=F))

mxEval(MZ.Acsub, ACEMMs1Fit)
mxEval(MZ.Ecsub, ACEMMs1Fit)
mxEval(MZ.Vcsub, ACEMMs1Fit)
mxEval(MZ.FVc44, ACEMMs1Fit)

mxEval(MZ.Stand_1on2[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_1on3[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_1on4[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_1on5[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_2on3[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_2on4[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_3on4[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_5on2[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_5on3[1,1], ACEMMs1Fit)
mxEval(MZ.Stand_5on4[1,1], ACEMMs1Fit)

mxEval(MZ.FactcorMZ, ACEMMs1Fit)
mxEval(DZ.FactcorDZ, ACEMMs1Fit)
mxEval(MZ.FVc, ACEMMs1Fit)
mxEval(MZ.PhC, ACEMMs1Fit)

mxEval(MZ.stAc, ACEMMs1Fit)
mxEval(MZ.stCc, ACEMMs1Fit)
mxEval(MZ.stEc, ACEMMs1Fit)

mxEval(MZ.stpac, ACEMMs1Fit)
mxEval(MZ.stpcc, ACEMMs1Fit)
mxEval(MZ.stpec, ACEMMs1Fit)

mxEval(MZ.Ra, ACEMMs1Fit)
mxEval(MZ.Rc, ACEMMs1Fit)
mxEval(MZ.Re, ACEMMs1Fit)
mxEval(MZ.Rph, ACEMMs1Fit)

mxEval(MZ.Rpha12, ACEMMs1Fit) 
mxEval(MZ.Rpha13, ACEMMs1Fit) 
mxEval(MZ.Rpha23, ACEMMs1Fit) 

mxEval(MZ.Acsub, ACEMMs1Fit)
mxEval(MZ.Ecsub, ACEMMs1Fit)
mxEval(MZ.Vcsub, ACEMMs1Fit)
mxEval(MZ.FVc44, ACEMMs1Fit)

mxEval(MZ.FactcorMZ, ACEMMs1Fit)
mxEval(DZ.FactcorDZ, ACEMMs1Fit)
mxEval(MZ.FVc, ACEMMs1Fit)
mxEval(MZ.PhC, ACEMMs1Fit)

mxEval(MZ.stAc, ACEMMs1Fit)
mxEval(MZ.stCc, ACEMMs1Fit)
mxEval(MZ.stEc, ACEMMs1Fit)

mxEval(MZ.Rpha12, ACEMMs1Fit) 
mxEval(MZ.Rpha13, ACEMMs1Fit) 
mxEval(MZ.Rpha23, ACEMMs1Fit) 

mxEval(MZ.StandFact, ACEMMs1Fit)
# 
#****************************************************************************************************************************
# __(VIIb)_____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL for PD>VICT>SO with mediation by Vict and confounding by CGN
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Phenotypic Factors: F1>F2>F3 & F1>F3;
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
#_____________________________________________________________________________________________________________________________

nv		<- 7				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 5				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nvo 		<- 1     			# number of ordinal variables per twin
nvc 		<- nv-nvo  			# number of continuous variables per twin
poso 		<- nvo 			# position where ordinal variables start
nth		<- 4				# number of max thresholds
ninc 		<- nth-1 			# number of max increments
ncovariates <- 2 				# number of covariates
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups	<- c("mz", "dz")
Vars		<- c('PRSAnx','PRSDep','Dep','Anx','Vict','SO','CGN4')
selVars	<- c('PRSAnx1','PRSDep1','Dep1','Anx1','Vict1','SO1','CGN41',
		     'PRSAnx2','PRSDep2','Dep2','Anx2','Vict2','SO2','CGN42')
useVars	<- c('PRSAnx1','PRSDep1','Dep1','Anx1','Vict1','SO1','CGN41',
		     'PRSAnx2','PRSDep2','Dep2','Anx2','Vict2','SO2','CGN42','age1','sex1','age2','sex2')

mzData		<- subset(TWINdata2, zyg1%in%c(1,3)|zyg2%in%c(1,3) , useVars)
dzData		<- subset(TWINdata2, zyg1%in%c(2,4,5)|zyg2%in%c(2,4,5) , useVars)

psych::describe(mzData)
psych::describe(dzData)

mzData$SO1[is.na(mzData$age1)] <- NA
mzData$SO2[is.na(mzData$age2)] <- NA
dzData$SO1[is.na(dzData$age1)] <- NA
dzData$SO2[is.na(dzData$age2)] <- NA

mzData$SO1[is.na(mzData$sex1)] <- NA
mzData$SO2[is.na(mzData$sex2)] <- NA
dzData$SO1[is.na(dzData$sex1)] <- NA
dzData$SO2[is.na(dzData$sex2)] <- NA

mzData$age1[is.na(mzData$age1)] <- 999
mzData$age2[is.na(mzData$age2)] <- 999
dzData$age1[is.na(dzData$age1)] <- 999
dzData$age2[is.na(dzData$age2)] <- 999

mzData$sex1[is.na(mzData$sex1)] <- 999
mzData$sex2[is.na(mzData$sex2)] <- 999
dzData$sex1[is.na(dzData$sex1)] <- 999
dzData$sex2[is.na(dzData$sex2)] <- 999

psych::describe(mzData)
psych::describe(dzData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabs	<- paste("m",1:nv,sep=""))
(Stmean	<- c(5,5,4.5,4.1,2.5,0,4))
(PatM		<- c(TRUE,T,TRUE,TRUE,TRUE,F,TRUE))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEs	<- c('es1','es1','es3','es3','es5','es6','es7'))
(LabAs	<- c('as1','as1','as3','as3','as5','as6','as7'))
#(LabCs	<- c('cs1','cs1','cs3','cs3','cs5','cs6','cs7'))

PatSpe	<- c(T,T,T,T,F,F,F)
PatSpac	<- c(F,F,T,T,F,F,F)
StSpa		<- c(0,0,.5,.5,0,0,0)
#StSpc		<- c(0,0,.5,.5,0,0,0)
StSpe		<- c(.5,.5,.5,.5,0,0,0)

# all 1st loadings fixed to 1
PatFl		<- c(F,T,F,F,F,F,F,			
		     F,F,F,T,F,F,F,
		     F,F,F,F,F,F,F,
		     F,F,F,F,F,F,F,
		     F,F,F,F,F,F,F)

StFl		<- c(1,.5,0,0,0,0,0,
		     0,0,1,.5,0,0,0,
		     0,0,0,0,1,0,0,
		     0,0,0,0,0,1,0,
		     0,0,0,0,0,0,1)

LabFl		<- c('l1','l2',NA,NA,NA,NA,NA,
	 	     NA,NA,'l3','l4',NA,NA,NA,
	 	     NA,NA,NA,NA,'l5',NA,NA,
	 	     NA,NA,NA,NA,NA,'l6',NA,
	 	     NA,NA,NA,NA,NA,NA,'l7')

PatPhC	<- c(F,T,T,T,T,
		     F,F,T,T,F,
		     F,F,F,T,F,
		     F,F,F,F,F,
		     F,T,T,T,F)

StPhC		<- c(0,.1,.1,.1,.1,
		     0,0,.1,.1,0,
		     0,0,0,.1,0,
		     0,0,0,0,0,
		     0,.1,.1,.1,0)

LabPhC	<- c(NA,'c1on2','c1on3','c1on4','c1on5',
		     NA,NA,'c2on3','c2on4',NA,
		     NA,NA,NA,'c3on4',NA,
		     NA,NA,NA,NA,NA,
		     NA,'c5on2','c5on3','c5on4',NA)	 

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Means		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmean,Stmean), labels=c(mLabs,mLabs), name="expMean") 

# Threshold and covariates
LabTh		<-c('T_1','i_11','i_12','i_13')	# THs for var 1 

LabCovA	<-c('BageThSO','BageThSO','BageThSO','BageThSO')
LabCovS 	<-c('BsexThSO','BsexThSO','BsexThSO','BsexThSO')

ThPat		<-c(T,T,T,T)
StTH		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

obssex1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex1"), name="sex1")
obssex2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.sex2"), name="sex2")

# effect of age and sex on ordinal variable
betaA		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.3, labels=LabCovA, name="BageTH" )
betaS		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.4, labels=LabCovS, name="BsexTH" )
 
# thresholds
Thr		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTH, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabTh, name="Th")
inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")
Thres		<-mxAlgebra( expression= cbind(Low%*%Th + BageTH%x%age1 + BsexTH%x%sex1,
                  Low%*%Th + BageTH%x%age2 + BsexTH%x%sex2), name="expThres")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Load		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFl, name="FactL" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTw	<-mxAlgebra(I2%x%FactL, name="FactLTw")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 4 latent variables  
PhCaus	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhC, labels=LabPhC, name="PhC" )

# Define the matrix to hold the A and C effects: Specific 
PathsAs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpa, labels=LabAs, name="as" )
#PathsCs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpc, labels=LabCs, name="cs" )
PathsEs	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpe, labels=LabEs, name="es" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
covAs		<-mxAlgebra( expression= as %*% t(as), name="As" )
#covCs		<-mxAlgebra( expression= cs %*% t(cs), name="Cs" )
covEs		<-mxAlgebra( expression= es %*% t(es), name="Es" )
covPs		<-mxAlgebra( expression= As+Es, name="Vs" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("ac22","ac32","ac42","ac52","ac33","ac43","ac53","ac44","ac54","ac55"), name="ac" ) # Component paths for factors 2 and 3
#PathsCcsub	<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("cc22","cc32","cc42","cc52","cc33","cc43","cc53","cc44","cc54","cc55"), name="cc" )
PathsEcsub	<-mxMatrix(type="Diag", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=.8, labels=c("ec22","ec33","ec44","ec55"), name="ec" )
PathsP11	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11", name="pc" ) # SD path for factor 3 (the PRS factor)
Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  #Padding

Ze41		<-mxMatrix(type="Zero",	nrow=4, ncol=1, free=F, name="Z41" )  #Padding
Ze14		<-mxMatrix(type="Zero",	nrow=1, ncol=4, free=F, name="Z14" )  #Padding
covAcsub	<-mxAlgebra( expression= ac %*% t(ac), name="Acsub" )
#covCcsub	<-mxAlgebra( expression= cc %*% t(cc), name="Ccsub" )
covEcsub	<-mxAlgebra( expression= ec %*% t(ec), name="Ecsub" )
covPcsub	<-mxAlgebra( expression= Acsub+Ecsub, name="Vcsub" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)

covPc11	<-mxAlgebra( expression= pc %*% t(pc), name="Pc11" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components

covPc		<-mxAlgebra(cbind(rbind(Pc11,Z41), rbind(Z14,Vcsub)), name="Vc") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPcMz	<-mxAlgebra(cbind(rbind(Pc11,Z41) ,rbind(Z14,Acsub)), name="Vcmz") #I specify the MZ between-twin covariance - excluding E parameters
covPcDz	<-mxAlgebra(cbind(rbind(.5%x%Pc11,Z41), rbind(Z14,.5%x%Acsub)), name="Vcdz") #I specify the DZ between-twin covariance - specifying half of A and excluding E

# Generate Covariance of Latent factor model Including Causal Paths between factors
Id5		<-mxMatrix(type="Iden",	nrow=5, ncol=5, free=F, name="I5" )
covFVc	<-mxAlgebra( expression= solve(I5-PhC) %&% Vc, name ="FVc")
covFcMz	<-mxAlgebra( expression= solve(I5-PhC) %&% Vcmz, name ="Fcmz")
covFcDz	<-mxAlgebra( expression= solve(I5-PhC) %&% Vcdz, name ="Fcdz")

# Constraint on total variance of Ordinal variable (A+C+E=1)
varL1		<- mxConstraint( expression=FVc[4,4]==1, name="L1" )

FcovMZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcmz), cbind(Fcmz, FVc))) , name="expFCovMZ" )
FcovDZ	<-mxAlgebra( expression= (FactLTw  %&% rbind ( cbind(FVc, Fcdz), cbind(Fcdz, FVc))) , name="expFCovDZ" )

SpcovMZ	<-mxAlgebra( expression= rbind (cbind(Vs, As), cbind(As, Vs)) , name="expSpCovMZ" )
SpcovDZ	<-mxAlgebra( expression= rbind (cbind(Vs, .5%x%As), cbind(.5%x%As, Vs)) , name="expSpCovDZ" )

TOTcovMZ	<-mxAlgebra( expression= expFCovMZ + expSpCovMZ , name="TOTexpCovMZ" )
TOTcovDZ	<-mxAlgebra( expression= expFCovDZ + expSpCovDZ , name="TOTexpCovDZ" )

# *******************************************************************************************************
# Calculator

# Standardize the causal effects
Stcp1on2	<-mxAlgebra( expression= (PhC[2,1]* sqrt(FVc[1,1]))/sqrt(FVc[2,2]) , name="Stand_1on2" )
Stcp1on3	<-mxAlgebra( expression= (PhC[3,1]* sqrt(FVc[1,1]))/sqrt(FVc[3,3]) , name="Stand_1on3" )
Stcp1on4	<-mxAlgebra( expression= (PhC[4,1]* sqrt(FVc[1,1]))/sqrt(FVc[4,4]) , name="Stand_1on4" )
Stcp1on5	<-mxAlgebra( expression= (PhC[5,1]* sqrt(FVc[1,1]))/sqrt(FVc[5,5]) , name="Stand_1on5" )
Stcp2on3	<-mxAlgebra( expression= (PhC[3,2]* sqrt(FVc[2,2]))/sqrt(FVc[3,3]) , name="Stand_2on3" )
Stcp2on4	<-mxAlgebra( expression= (PhC[4,2]* sqrt(FVc[2,2]))/sqrt(FVc[4,4]) , name="Stand_2on4" )
Stcp3on4	<-mxAlgebra( expression= (PhC[4,3]* sqrt(FVc[3,3]))/sqrt(FVc[4,4]) , name="Stand_3on4" )
Stcp5on2	<-mxAlgebra( expression= (PhC[2,5]* sqrt(FVc[5,5]))/sqrt(FVc[2,2]) , name="Stand_5on2" )
Stcp5on3	<-mxAlgebra( expression= (PhC[3,5]* sqrt(FVc[5,5]))/sqrt(FVc[3,3]) , name="Stand_5on3" )
Stcp5on4	<-mxAlgebra( expression= (PhC[4,5]* sqrt(FVc[5,5]))/sqrt(FVc[4,4]) , name="Stand_5on4" )

# Standardize the Total var/covariances matrices of the observed variables
Id14		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I14" )
Rfactmz	<-mxAlgebra( expression= solve(sqrt(I14*TOTexpCovMZ)) %&% TOTexpCovMZ, name="FactcorMZ" )
Rfactdz	<-mxAlgebra( expression= solve(sqrt(I14*TOTexpCovDZ)) %&% TOTexpCovDZ, name="FactcorDZ" )

# Phenotypic, A, C and E correlations	
Id4		<-mxMatrix(type="Iden",	nrow=4, ncol=4, free=F, name="I4" )
RfactAc	<-mxAlgebra( expression= solve(sqrt(I4*Acsub)) %&% Acsub, name="Ra" )
#RfactCc	<-mxAlgebra( expression= solve(sqrt(I4*Ccsub)) %&% Ccsub, name="Rc" )
RfactEc	<-mxAlgebra( expression= solve(sqrt(I4*Ecsub)) %&% Ecsub, name="Re" )
RfactP	<-mxAlgebra( expression= solve(sqrt(I5*FVc)) %&% FVc, name="Rph" )

# Standardize the Common Effects
covFVc44	<-mxAlgebra( expression= FVc[2:5,2:5], name ="FVc44")
stcovAc	<-mxAlgebra( expression= Acsub/FVc44, name="stAc" )
#stcovCc	<-mxAlgebra( expression= Ccsub/FVc44, name="stCc" )
stcovEc	<-mxAlgebra( expression= Ecsub/FVc44, name="stEc" )

# Standardised path estimates
StpathAc	<-mxAlgebra( expression= (sqrt(stAc)), name="stpac" )
#StpathCc	<-mxAlgebra( expression= (sqrt(stCc)), name="stpcc" )
StpathEc	<-mxAlgebra( expression= (sqrt(stEc)), name="stpec" )

# Algebra to compute Rph-A, Rph-C and Rph-E
RphA12	<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[2,1]*sqrt(stAc[2,2]), name = 'Rpha12')
RphA13	<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[3,1]*sqrt(stAc[3,3]), name = 'Rpha13')
RphA14	<-mxAlgebra(expression=sqrt(stAc[1,1])*Ra[4,1]*sqrt(stAc[4,4]), name = 'Rpha14')
RphA23	<-mxAlgebra(expression=sqrt(stAc[2,2])*Ra[3,2]*sqrt(stAc[3,3]), name = 'Rpha23')
RphA24	<-mxAlgebra(expression=sqrt(stAc[2,2])*Ra[4,2]*sqrt(stAc[4,4]), name = 'Rpha24')
RphA34	<-mxAlgebra(expression=sqrt(stAc[3,3])*Ra[4,3]*sqrt(stAc[4,4]), name = 'Rpha34')
#RphC12	<-mxAlgebra(expression=sqrt(stCc[1,1])*Rc[2,1]*sqrt(stCc[2,2]), name = 'Rphc12')
#RphC13	<-mxAlgebra(expression=sqrt(stCc[1,1])*Rc[3,1]*sqrt(stCc[3,3]), name = 'Rphc13')
#RphC23	<-mxAlgebra(expression=sqrt(stCc[2,2])*Rc[3,2]*sqrt(stCc[3,3]), name = 'Rphc23')
RphE12	<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[2,1]*sqrt(stEc[2,2]), name = 'Rphe12')
RphE13	<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[3,1]*sqrt(stEc[3,3]), name = 'Rphe13')
RphE14	<-mxAlgebra(expression=sqrt(stEc[1,1])*Re[4,1]*sqrt(stEc[4,4]), name = 'Rphe14')
RphE23	<-mxAlgebra(expression=sqrt(stEc[2,2])*Re[3,2]*sqrt(stEc[3,3]), name = 'Rphe23')
RphE24	<-mxAlgebra(expression=sqrt(stEc[2,2])*Re[4,2]*sqrt(stEc[4,4]), name = 'Rphe24')
RphE34	<-mxAlgebra(expression=sqrt(stEc[3,3])*Re[4,3]*sqrt(stEc[4,4]), name = 'Rphe34')

# Standardize the Specific Effects
stcovAs	<-mxAlgebra( expression= sqrt(As/( (FactL %&% FVc) +Vs)), name="stAs" )
#stcovCs	<-mxAlgebra( expression= sqrt(Cs/( (FactL %&% FVc) +Vs)), name="stCs" )
stcovEs	<-mxAlgebra( expression= sqrt(Es/( (FactL %&% FVc) +Vs)), name="stEs" )

# Standardized Factor Loadings
StFL		<-mxAlgebra( expression= sqrt(diag2vec( FactL %&% FVc / TOTexpCovMZ[1:7,1:7])) , name="StandFact" )

# *******************************************************************************************************

# Data objects for Multiple Groups
dataMZ	<- mxData( observed=mzData, type="raw" )
dataDZ	<- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
objMZ		<- mxExpectationNormal( covariance="TOTexpCovMZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))
objDZ		<- mxExpectationNormal( covariance="TOTexpCovDZ", means="expMean", dimnames=selVars, thresholds="expThres", threshnames=c("SO1","SO2"))

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1		<-list(Means,Load,LoadTw,PhCaus,PathsAs,PathsEs,covAs,covEs,covPs,Id2,Id4,Id5,Id14)
pars2		<-list(PathsAcsub,PathsEcsub,PathsP11,Ze41,Ze14,covAcsub,covEcsub,covPcsub,covPc11,covPc,covPcMz,covPcDz,covFVc,covFVc44,covFcMz,covFcDz)
pars3		<-list(obsage1, obsage2, obssex1, obssex2, betaA, betaS, Thr, inc, Thres)
parsst	<-list(stcovAs, stcovEs, stcovAc, stcovEc, RfactAc, RfactEc, RfactP,RphA12,RphA13,RphA14,RphA23,RphA24,RphA34,RphE12,RphE13,RphE14,RphE23,RphE24,RphE34,StpathAc,StpathEc)
parsmed	<-list(Stcp1on2, Stcp1on3, Stcp1on4, Stcp1on5, Stcp2on3, Stcp2on4, Stcp3on4, Stcp5on2, Stcp5on3, Stcp5on4)
modelMZ	<-mxModel(pars1, pars2, pars3, parsmed, FcovMZ, SpcovMZ, TOTcovMZ, dataMZ, objMZ, Rfactmz, parsst, fitFunction, StFL, varL1, name="MZ" )
modelDZ	<-mxModel(pars1, pars2, pars3, FcovDZ, SpcovDZ, TOTcovDZ, dataDZ, objDZ, Rfactdz, fitFunction, name="DZ" )

minus2ll	<-mxAlgebra( expression=MZ.objective + DZ.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cistFL	<-mxCI (c ('MZ.StandFact','MZ.Stand_1on2','MZ.Stand_1on3','MZ.Stand_1on4','MZ.Stand_1on5','MZ.Stand_2on3','MZ.Stand_2on4','MZ.Stand_3on4','MZ.Stand_5on2','MZ.Stand_5on3','MZ.Stand_5on4','MZ.PhC'))
cistVs	<-mxCI (c ('MZ.stAs[3,3]','MZ.stAs[4,4]',
				'MZ.stEs[1,1]','MZ.stEs[2,2]','MZ.stEs[3,3]','MZ.stEs[4,4]') ) 	# standardized var comp from specific Factors
cistVc	<-mxCI (c ('MZ.stAc[1,1]','MZ.stAc[2,1]','MZ.stAc[3,1]','MZ.stAc[4,1]','MZ.stAc[2,2]','MZ.stAc[3,2]','MZ.stAc[4,2]','MZ.stAc[3,3]','MZ.stAc[4,3]','MZ.stAc[4,4]',
				'MZ.stEc[1,1]','MZ.stEc[2,2]','MZ.stEc[3,3]','MZ.stEc[4,4]') ) 	# standardized var comp for ACE on latent Factors
cistRc	<-mxCI (c ('MZ.Rpha12','MZ.Rpha13','MZ.Rpha14','MZ.Rpha23','MZ.Rpha24','MZ.Rpha34','MZ.Rphe12','MZ.Rphe13','MZ.Rphe14','MZ.Rphe23','MZ.Rphe24','MZ.Rphe34','MZ.Ra','MZ.Re','MZ.stpac','MZ.stpec') ) 	
ACEMMs2Model	<-mxModel("aceMMs2", pars1, pars2, modelMZ, modelDZ, minus2ll, obj, cistFL, cistVs, cistVc, cistRc) 

# --------------------------------------------------------------------------------------------------------------------------------
# 5a RUN ACEMs Factor Model with phenotypic causal mediation paths by Zygosity

ACEMMs2Fit		<-mxTryHardOrdinal(ACEMMs2Model, intervals=T)
(ACEMMs2Summ	<-summary(ACEMMs2Fit, verbose=F))

mxEval(MZ.Acsub, ACEMMs2Fit)
mxEval(MZ.Ecsub, ACEMMs2Fit)
mxEval(MZ.Vcsub, ACEMMs2Fit)
mxEval(MZ.FVc44, ACEMMs2Fit)

mxEval(MZ.Stand_1on2[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_1on3[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_1on4[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_1on5[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_2on3[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_2on4[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_3on4[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_5on2[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_5on3[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_5on4[1,1], ACEMMs2Fit)

mxEval(MZ.FactcorMZ, ACEMMs2Fit)
mxEval(DZ.FactcorDZ, ACEMMs2Fit)
mxEval(MZ.FVc, ACEMMs2Fit)
mxEval(MZ.PhC, ACEMMs2Fit)

mxEval(MZ.stAc, ACEMMs2Fit)
mxEval(MZ.stCc, ACEMMs2Fit)
mxEval(MZ.stEc, ACEMMs2Fit)

mxEval(MZ.stpac, ACEMMs2Fit)
mxEval(MZ.stpcc, ACEMMs2Fit)
mxEval(MZ.stpec, ACEMMs2Fit)

mxEval(MZ.Ra, ACEMMs2Fit)
mxEval(MZ.Rc, ACEMMs2Fit)
mxEval(MZ.Re, ACEMMs2Fit)
mxEval(MZ.Rph, ACEMMs2Fit)

mxEval(MZ.Rpha12, ACEMMs2Fit) 
mxEval(MZ.Rpha13, ACEMMs2Fit) 
mxEval(MZ.Rpha23, ACEMMs2Fit) 

mxEval(MZ.StandFact, ACEMMs2Fit)
mxEval(MZ.Stand_1on2[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_1on3[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_1on4[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_2on3[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_2on4[1,1], ACEMMs2Fit)
mxEval(MZ.Stand_3on4[1,1], ACEMMs2Fit)

#> (( -0.04707845)*(sqrt(0.51157371)))/sqrt(1.01976880)
#[1] -0.03334459


#------------------------------------------------------------------------
# Submodel 5aii: Drop correlation path c paths from previous model 
#------------------------------------------------------------------------
# Drop correlation c paths from the model
# -----------------------------------------------------------------------
AEMMs1aModel	<- mxModel(ACEMMs1Fit, name="AEMMs1a")
AEMMs1aModel	<- omxSetParameters(AEMMs1aModel, labels=c('cc22','cc32','cc42','cc33','cc43','cc44'), free=FALSE, values=0)
AEMMs1aModel	<- omxSetParameters(AEMMs1aModel, labels=c('cs4'), free=FALSE, values=0)
AEMMs1aFit		<- mxTryHardOrdinal(AEMMs1aModel, intervals=T)
(AEMMs1aSum		<- summary(AEMMs1aFit))

mxCompare(ACEMMs1Fit, AEMMs1aFit)

mxEval(MZ.FactcorMZ, AEMMs1aFit)
mxEval(DZ.FactcorDZ, AEMMs1aFit)

mxEval(MZ.Acsub, AEMMs1aFit)
mxEval(MZ.Ccsub, AEMMs1aFit)
mxEval(MZ.Ecsub, AEMMs1aFit)
mxEval(MZ.Vcsub, AEMMs1aFit)
mxEval(MZ.FVc33, AEMMs1aFit)

mxEval(MZ.FVc, AEMMs1aFit)

mxEval(MZ.stAc, AEMMs1aFit)
mxEval(MZ.stCc, AEMMs1aFit)
mxEval(MZ.stEc, AEMMs1aFit)

mxEval(MZ.stpac, AEMMs1aFit)
mxEval(MZ.stpcc, AEMMs1aFit)
mxEval(MZ.stpec, AEMMs1aFit)

mxEval(MZ.Ra, AEMMs1aFit)
mxEval(MZ.Rc, AEMMs1aFit)
mxEval(MZ.Re, AEMMs1aFit)
mxEval(MZ.Rph, AEMMs1aFit)

mxEval(MZ.Rpha12, AEMMs1aFit) 
mxEval(MZ.Rpha13, AEMMs1aFit) 
mxEval(MZ.Rpha23, AEMMs1aFit) 

mxEval(MZ.StandFact, AEMMs1aFit)
mxEval(MZ.Stand_1on2[1,1], AEMMs1aFit)
mxEval(MZ.Stand_1on3[1,1], AEMMs1aFit)
mxEval(MZ.Stand_1on4[1,1], AEMMs1aFit)
mxEval(MZ.Stand_2on3[1,1], AEMMs1aFit)
mxEval(MZ.Stand_2on4[1,1], AEMMs1aFit)
mxEval(MZ.Stand_3on4[1,1], AEMMs1aFit)


########################################################################################################
##SEX DIFFERENCES MODELS
########################################################################################################

#*******************************************************************************************************
# __(IIa - Sex Diff)_________________________________________________________________________________________________
# Phenotypic Covariance Model across Latent Constructs  with Vict with sex differences
# Restrictions: means and variances equated across birth-order & zygosity groups;
# One set of factor loadings; one set of correltional paths between the factors; one set of error terms
# We estimate the factor variances, giving them a scale by fixing the loading on the 1st variable to 1
# This model specifies a full var/cov structure between the latent factors for MZ and DZ twins 
#______________________________________________________________________________________________________

nv			<- 7				# number of variables for a twin = 1 in Univariate
nvo 			<- 1     			# number of ordinal variables per twin
nvc 			<- nv-nvo  			# number of continuous variables per twin
poso 			<- nvo 			# position where ordinal variables start
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nth			<- 4				# number of max thresholds
nlower		<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor			<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
ninc 			<- nth-1 			# number of max increments
ncovariates 	<- 2 				# number of covariates

nfact			<- 5				# number of Latent Factors for Mediation Model per twin
nfact2		<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nfcor			<-(nfact*(nfact+1)/2)-nfact	# number of free elements in a correlation matrix nfact*nfcat

Groups		<- c("mzm","mzf","dzm","dzf","dzo") 
#Vars			<- c('SO','Dep','Anx','Vict','PRSSO','PRSDep','PRSAnx')
Vars			<- c('SO','Dep','Anx','Vict','PRSSO','PRSDep','PRSAnx')
selVars		<- c('SO1','Dep1','Anx1','Vict1','PRSSO1','PRSDep1','PRSAnx1',
			     'SO2','Dep2','Anx2','Vict2','PRSSO2','PRSDep2','PRSAnx2')
useVars		<- c('SO1','Dep1','Anx1','Vict1','PRSSO1','PRSDep1','PRSAnx1',
			     'SO2','Dep2','Anx2','Vict2','PRSSO2','PRSDep2','PRSAnx2','age1','sex1','age2','sex2')

mzmData		<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1) , useVars)
mzfData		<- subset(TWINdata2, zyg1%in%c(3)|zyg2%in%c(3) , useVars)
dzmData		<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2) , useVars)
dzfData		<- subset(TWINdata2, zyg1%in%c(4)|zyg2%in%c(4) , useVars)
dzoData		<- subset(TWINdata2, zyg1%in%c(5)|zyg2%in%c(5) , useVars)

psych::describe(mzmData)
psych::describe(mzfData)
psych::describe(dzmData)
psych::describe(dzfData)
psych::describe(dzoData)

mzmData$SO1[is.na(mzmData$age1)] <- NA
mzmData$SO2[is.na(mzmData$age2)] <- NA
mzfData$SO1[is.na(mzfData$age1)] <- NA
mzfData$SO2[is.na(mzfData$age2)] <- NA
dzmData$SO1[is.na(dzmData$age1)] <- NA
dzmData$SO2[is.na(dzmData$age2)] <- NA
dzfData$SO1[is.na(dzfData$age1)] <- NA
dzfData$SO2[is.na(dzfData$age2)] <- NA
dzoData$SO1[is.na(dzoData$age1)] <- NA
dzoData$SO2[is.na(dzoData$age2)] <- NA

mzmData$age1[is.na(mzmData$age1)] <- 999
mzmData$age2[is.na(mzmData$age2)] <- 999
mzfData$age1[is.na(mzfData$age1)] <- 999
mzfData$age2[is.na(mzfData$age2)] <- 999
dzmData$age1[is.na(dzmData$age1)] <- 999
dzmData$age2[is.na(dzmData$age2)] <- 999
dzfData$age1[is.na(dzfData$age1)] <- 999
dzfData$age2[is.na(dzfData$age2)] <- 999
dzoData$age1[is.na(dzoData$age1)] <- 999
dzoData$age2[is.na(dzoData$age2)] <- 999

psych::describe(mzmData)
psych::describe(mzfData)
psych::describe(dzmData)
psych::describe(dzfData)
psych::describe(dzoData)


# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)

(Stmeanm	<-colMeans(mzmData[,2:nv],na.rm=TRUE))
StMm 		<-c(0, Stmeanm, 0, Stmeanm)
(LabMm	<- paste("mm",2:nv,sep=""))
MLabsm	<-c(NA,LabMm,NA,LabMm) 

(Stmeanf	<-colMeans(mzfData[,2:nv],na.rm=TRUE))
StMf 		<-c(0, Stmeanf, 0, Stmeanf)
(LabMf	<- paste("mf",2:nv,sep=""))
MLabsf	<-c(NA,LabMf,NA,LabMf) 

StMo 		<-c(0, Stmeanm, 0, Stmeanf)
MLabso	<-c(NA,LabMm,NA,LabMf) 

(LabErm	<-c("e1m","e2m","e2m","e4m","e5m","e6m","e6m"))
(LabErf	<-c("e1f","e2f","e2f","e4f","e5f","e6f","e6f"))

# Create Labels for the Factor parameters
(sdLabsm	<- paste("sdm",1:nfact,sep=""))	# SD
(rphvLabsm	<- paste("rvm",1:3,sep="")) #Within person - for the first 3 factors - different in mz
(rphPLabsm	<- paste("rPm",1,sep="")) #Within person - for the last 2 polygenic factors - same in mz
(rphvPLabsm	<- paste("rvPm",1:6,sep="")) #Within person - for the cross PRS-variable correlations - same in mz
(MZMbLabs 	<- paste("rmzm", do.call(c, sapply(seq(1, 3), function(x){ paste(x:3, x,sep="") })), sep="")) # first 2 vars in mz
(DZMbLabs 	<- paste("rdzm", do.call(c, sapply(seq(1, nfact), function(x){ paste(x:nfact, x,sep="") })), sep="")) # all vars in dz

(sdLabsf	<- paste("sdf",1:nfact,sep=""))	# SD
(rphvLabsf	<- paste("rvf",1:3,sep="")) #Within person - for the first 3 factors - different in mz
(rphPLabsf	<- paste("rPf",1,sep="")) #Within person - for the last 2 polygenic factors - same in mz
(rphvPLabsf	<- paste("rvPf",1:6,sep="")) #Within person - for the cross PRS-variable correlations - same in mz
(MZFbLabs 	<- paste("rmzf", do.call(c, sapply(seq(1, 3), function(x){ paste(x:3, x,sep="") })), sep="")) # first 2 vars in mz
(DZFbLabs 	<- paste("rdzf", do.call(c, sapply(seq(1, nfact), function(x){ paste(x:nfact, x,sep="") })), sep="")) # all vars in dz
(DZObLabs 	<- paste("rdzo", do.call(c, sapply(seq(1, nfact), function(x){ paste(x:nfact, x,sep="") })), sep="")) # all vars in dz

(StWithinpersonv 	<-c(.2))
(StWithinpersonP 	<-c(.2))
(StWithinpersonvP	<-c(.1))

(StBetweenMZM  	<-c(.5,.2,.2,.5,.2,.5))
(StBetweenDZM  	<-c(.3,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.5,.1,.5))
(StBetweenMZF  	<-c(.5,.2,.2,.5,.2,.5))
(StBetweenDZF  	<-c(.3,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.5,.1,.5))

(PatBetweenDZ  	<-c(T,T,T,T,T,T,T,T,T,T,T,T,F,T,F))

# Create Labels for the Factor Loadings (1st loadings fixed to 1)

PatFl	<- c(F,F,F,F,F,F,F,			
	     F,T,F,F,F,F,F,
	     F,F,F,F,F,F,F,
	     F,F,F,F,F,F,F,
	     F,F,F,F,F,F,T)

StFl	<- c(1,0,0,0,0,0,0,
	     0,.5,1,0,0,0,0,
	     0,0,0,1,0,0,0,
	     0,0,0,0,1,0,0,
	     0,0,0,0,0,1,.5)

LabFlm	<- c('l1m',NA,NA,NA,NA,NA,NA,
		      NA,'l2m','l3m',NA,NA,NA,NA,
		      NA,NA,NA,'l4m',NA,NA,NA,
		      NA,NA,NA,NA,'l5m',NA,NA,
		      NA,NA,NA,NA,NA,'l6m','l7m')

LabFlf	<- c('l1f',NA,NA,NA,NA,NA,NA,
		      NA,'l2f','l3f',NA,NA,NA,NA,
		      NA,NA,NA,'l4f',NA,NA,NA,
		      NA,NA,NA,NA,'l5f',NA,NA,
		      NA,NA,NA,NA,NA,'l6f','l7f')

# Free parameters
(Pat  	<- c( rep(FALSE,nvo), rep(TRUE, nvc)))

# ______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
# ______________________________________________________________________________________________________

Meanm	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(Pat, Pat), values=c(StMm), labels=c(MLabsm), name="expmm" ) #first variable is not free to be estimated
Meanf	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(Pat, Pat), values=c(StMf), labels=c(MLabsf), name="expmf" ) 
Meano	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(Pat, Pat), values=c(StMo), labels=c(MLabso), name="expmo" ) 

# I constrain the threshold and increments to be equal across birth order or zygosity and specify only one set
LabThm	<-c('T_1m','i_11m','i_12m','i_13m')	# THs for var 1 for mz twin 1 (mzm)
LabThf	<-c('T_1f','i_11f','i_12f','i_13f')	# THs for var 1 for mz twin 1 (mzf)

LabCovAm	<-c('BageThSOm','BageThSOm','BageThSOm','BageThSOm')
LabCovAf	<-c('BageThSOf','BageThSOf','BageThSOf','BageThSOf')

ThPat		<-c(TRUE,TRUE,TRUE,TRUE)
StTHm		<-c(0.84,0.57,0.14,0.20)
StTHf		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

# effect of age and sex on ordinal variable
betaAm	<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.05, labels=LabCovAm, name="BageTHm" )
betaAf	<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.05, labels=LabCovAf, name="BageTHf" )
 
# thresholds
Thrm		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTHm, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabThm, name="Thm")
incm		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=F, values=1, name="Lowm")

Thrf		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTHf, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabThf, name="Thf")
incf		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=F, values=1, name="Lowf")

Thresm	<-mxAlgebra( expression= cbind(Lowm%*%Thm + BageTHm%x%age1, Lowm%*%Thm + BageTHm%x%age2), name="expThresm")
Thresf	<-mxAlgebra( expression= cbind(Lowf%*%Thf + BageTHf%x%age1, Lowf%*%Thf + BageTHf%x%age2), name="expThresf")
Threso	<-mxAlgebra( expression= cbind(Lowm%*%Thm + BageTHm%x%age1, Lowf%*%Thf + BageTHf%x%age2), name="expThreso")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Loadm		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFlm, name="FactLm" )
Loadf		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFlf, name="FactLf" )
Ze75		<-mxMatrix("Zero", nv, nfact, free=F, name="Z75")
LoadTwm	<-mxAlgebra(rbind(cbind(FactLm,Z75), cbind(Z75, FactLm)), name="FactLTwm")
LoadTwf	<-mxAlgebra(rbind(cbind(FactLf,Z75), cbind(Z75, FactLf)), name="FactLTwf")
LoadTwo	<-mxAlgebra(rbind(cbind(FactLm,Z75), cbind(Z75, FactLf)), name="FactLTwo")

ErPathm	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=c(F,T,T,F,F,T,T), values=c(0,.5,.5,0,0,.5,.5), labels=LabErm, name="Erpm" )
ErPathf	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=c(F,T,T,F,F,T,T), values=c(0,.5,.5,0,0,.5,.5), labels=LabErf, name="Erpf" )

Erm		<-mxAlgebra(Erpm %*% t(Erpm), name="Errorm")
Erf		<-mxAlgebra(Erpf %*% t(Erpf), name="Errorf")

Ze7		<-mxMatrix("Zero", nv, nv, free=F, name="Z7")
ErTwm		<-mxAlgebra(rbind(cbind(Errorm,Z7), cbind(Z7, Errorm)), name="ErrorTwm")
ErTwf		<-mxAlgebra(rbind(cbind(Errorf,Z7), cbind(Z7, Errorf)), name="ErrorTwf")
ErTwo		<-mxAlgebra(rbind(cbind(Errorm,Z7), cbind(Z7, Errorf)), name="ErrorTwo")
 
# elements for the SD of Factors
Id2	<-mxMatrix("Iden", 2, 2, free=F, name="I2")
sdFm	<-mxMatrix("Diag", nfact, nfact, free=c(F,T,T,T,T), values=1, labels=sdLabsm, name="SDfm") 
sdFf	<-mxMatrix("Diag", nfact, nfact, free=c(F,T,T,T,T), values=1, labels=sdLabsf, name="SDff") 

sdFTwm	<-mxAlgebra(I2 %x% SDfm, name="SDftwinm")
sdFTwf	<-mxAlgebra(I2 %x% SDff, name="SDftwinf")
Ze5		<-mxMatrix("Zero", nfact, nfact, free=F, name="Z5")
sdFTwo	<-mxAlgebra(rbind(cbind(SDfm,Z5), cbind(Z5, SDff)), name="SDftwino")

# elements for the correlations of Factors
Rphvm		<-mxMatrix("Stand", 3, 3, free = TRUE, values = StWithinpersonv, labels=rphvLabsm, lbound=-.999, ubound=.999, name="withinvm") 
RphPm		<-mxMatrix("Stand", 2, 2, free = TRUE, values = StWithinpersonP, labels=rphPLabsm, lbound=-.999, ubound=.999, name="withinPm") 
RphvPm	<-mxMatrix("Full", 2, 3, free = TRUE, values = StWithinpersonvP, labels=rphvPLabsm, lbound=-.999, ubound=.999, name="withinvPm") 
RphPvm	<-mxMatrix("Full", 3, 2, free = TRUE, values = StWithinpersonvP, labels=rphvPLabsm, lbound=-.999, ubound=.999, byrow=T, name="withinPvm") 
Rphm		<-mxAlgebra(rbind(cbind(withinvm,withinPvm), cbind(withinvPm, withinPm)), name="Rwithinm")

Rphvf		<-mxMatrix("Stand", 3, 3, free = TRUE, values = StWithinpersonv, labels=rphvLabsf, lbound=-.999, ubound=.999, name="withinvf") 
RphPf		<-mxMatrix("Stand", 2, 2, free = TRUE, values = StWithinpersonP, labels=rphPLabsf, lbound=-.999, ubound=.999, name="withinPf") 
RphvPf	<-mxMatrix("Full", 2, 3, free = TRUE, values = StWithinpersonvP, labels=rphvPLabsf, lbound=-.999, ubound=.999, name="withinvPf") 
RphPvf	<-mxMatrix("Full", 3, 2, free = TRUE, values = StWithinpersonvP, labels=rphvPLabsf, lbound=-.999, ubound=.999, byrow=T, name="withinPvf") 
Rphf		<-mxAlgebra(rbind(cbind(withinvf,withinPvf), cbind(withinvPf, withinPf)), name="Rwithinf")

MZMbv		<-mxMatrix("Symm", 3, 3, free = TRUE, values = StBetweenMZM, labels=MZMbLabs, lbound=-.999, ubound=.999, name="BetweenMZMv") 
MZFbv		<-mxMatrix("Symm", 3, 3, free = TRUE, values = StBetweenMZF, labels=MZFbLabs, lbound=-.999, ubound=.999, name="BetweenMZFv") 
MZMb		<-mxAlgebra(rbind(cbind(BetweenMZMv,withinPvm), cbind(withinvPm, withinPm)), name="RbetweenMZM")
MZFb		<-mxAlgebra(rbind(cbind(BetweenMZFv,withinPvf), cbind(withinvPf, withinPf)), name="RbetweenMZF")

DZMb		<-mxMatrix("Symm", nfact, nfact, free=PatBetweenDZ, values=StBetweenDZM, labels=DZMbLabs, lbound=-.999, ubound=.999, name="RbetweenDZM") 
DZFb		<-mxMatrix("Symm", nfact, nfact, free=PatBetweenDZ, values=StBetweenDZF, labels=DZFbLabs, lbound=-.999, ubound=.999, name="RbetweenDZF") 
DZOb		<-mxMatrix("Symm", nfact, nfact, free=PatBetweenDZ, values=StBetweenDZF, labels=DZObLabs, lbound=-.999, ubound=.999, name="RbetweenDZO") 

FactCorMZM	<-mxAlgebra(rbind(cbind(Rwithinm,RbetweenMZM), cbind(RbetweenMZM, Rwithinm)), name="RMZM")
FactCorMZF	<-mxAlgebra(rbind(cbind(Rwithinf,RbetweenMZF), cbind(RbetweenMZF, Rwithinf)), name="RMZF")
FactCorDZM	<-mxAlgebra(rbind(cbind(Rwithinm,RbetweenDZM), cbind(RbetweenDZM, Rwithinm)), name="RDZM")
FactCorDZF	<-mxAlgebra(rbind(cbind(Rwithinf,RbetweenDZF), cbind(RbetweenDZF, Rwithinf)), name="RDZF")
FactCorDZO	<-mxAlgebra(rbind(cbind(Rwithinm,RbetweenDZO), cbind(RbetweenDZO, Rwithinf)), name="RDZO")

# Generate expected Covariance matrices of Factors
FactCovMZM	<-mxAlgebra(SDftwinm %&% RMZM , name="expFactCovMZM")
FactCovMZF	<-mxAlgebra(SDftwinf %&% RMZF , name="expFactCovMZF")
FactCovDZM	<-mxAlgebra(SDftwinm %&% RDZM , name="expFactCovDZM")
FactCovDZF	<-mxAlgebra(SDftwinf %&% RDZF , name="expFactCovDZF")
FactCovDZO	<-mxAlgebra(SDftwino %&% RDZO , name="expFactCovDZO")

## This second step then derives the var/cov matrix of the observed/measured variables in terms of the variance/covariances of the latent factors and the Factor Loadings
covMZM	<-mxAlgebra( expression= FactLTwm  %&% expFactCovMZM , name="ExpCovMZM" )
covMZF	<-mxAlgebra( expression= FactLTwf  %&% expFactCovMZF , name="ExpCovMZF" )
covDZM	<-mxAlgebra( expression= FactLTwm  %&% expFactCovDZM , name="ExpCovDZM" )
covDZF	<-mxAlgebra( expression= FactLTwf  %&% expFactCovDZF , name="ExpCovDZF" )
covDZO	<-mxAlgebra( expression= FactLTwo  %&% expFactCovDZO , name="ExpCovDZO" )

## Finally, we derive the total expected variance/covariances for the measured variables which go in the models
TOTcovMZM	<-mxAlgebra( expression= ExpCovMZM + ErrorTwm , name="TOTexpCovMZM" )
TOTcovMZF	<-mxAlgebra( expression= ExpCovMZF + ErrorTwf , name="TOTexpCovMZF" )
TOTcovDZM	<-mxAlgebra( expression= ExpCovDZM + ErrorTwm , name="TOTexpCovDZM" )
TOTcovDZF	<-mxAlgebra( expression= ExpCovDZF + ErrorTwf , name="TOTexpCovDZF" )
TOTcovDZO	<-mxAlgebra( expression= ExpCovDZO + ErrorTwo , name="TOTexpCovDZO" )

# Standardizing parameters **********************

# Standardized Factor Loadings
StFLm		<-mxAlgebra( expression= sqrt(diag2vec( FactLm %&% expFactCovMZM[1:5,1:5] / TOTexpCovMZM[1:7,1:7])), name="StandFactm")
StFLf		<-mxAlgebra( expression= sqrt(diag2vec( FactLf %&% expFactCovMZF[1:5,1:5] / TOTexpCovMZF[1:7,1:7])), name="StandFactf")

# Standardise error terms of the measured variables
StErm		<-mxAlgebra( expression= sqrt(diag2vec( Errorm/TOTexpCovMZM[1:7,1:7])), name="StandErm")
StErf		<-mxAlgebra( expression= sqrt(diag2vec( Errorf/TOTexpCovMZF[1:7,1:7])), name="StandErf")

# ************************************

# Data objects for Multiple Groups
dataMZM	<- mxData( observed=mzmData, type="raw" )
dataMZF	<- mxData( observed=mzfData, type="raw" )
dataDZM	<- mxData( observed=dzmData, type="raw" )
dataDZF	<- mxData( observed=dzfData, type="raw" )
dataDZO	<- mxData( observed=dzoData, type="raw" )

# Objective objects for Multiple Groups
objMZM	<- mxExpectationNormal( covariance="TOTexpCovMZM", means="expmm", dimnames=selVars, thresholds="expThresm", threshnames=c("SO1","SO2"))
objMZF	<- mxExpectationNormal( covariance="TOTexpCovMZF", means="expmf", dimnames=selVars, thresholds="expThresf", threshnames=c("SO1","SO2"))
objDZM	<- mxExpectationNormal( covariance="TOTexpCovDZM", means="expmm", dimnames=selVars, thresholds="expThresm", threshnames=c("SO1","SO2"))
objDZF	<- mxExpectationNormal( covariance="TOTexpCovDZF", means="expmf", dimnames=selVars, thresholds="expThresf", threshnames=c("SO1","SO2"))
objDZO	<- mxExpectationNormal( covariance="TOTexpCovDZO", means="expmo", dimnames=selVars, thresholds="expThreso", threshnames=c("SO1","SO2"))

fitFunction <- mxFitFunctionML()

# Combine Groups
pars1m	<-list(Meanm, Loadm, Ze75, Ze7, LoadTwm, ErPathm, Erm, ErTwm, Id2, sdFm, sdFTwm, Rphm, Rphvm, RphPm, RphvPm, RphPvm)
pars1f	<-list(Meanf, Loadf, Ze75, Ze7, LoadTwf, ErPathf, Erf, ErTwf, Id2, sdFf, sdFTwf, Rphf, Rphvf, RphPf, RphvPf, RphPvf)
pars1o	<-list(Meano, LoadTwo, ErTwo, sdFTwo)
pars2m	<-list(obsage1, obsage2, betaAm, Thrm, incm, Thresm)
pars2f	<-list(obsage1, obsage2, betaAf, Thrf, incf, Thresf)

modelMZM	<-mxModel(pars1m, pars2m, MZMb, MZMbv, FactCorMZM, FactCovMZM, covMZM, TOTcovMZM, dataMZM, objMZM, fitFunction, StFLm, StErm, name="MZM" )
modelMZF	<-mxModel(pars1f, pars2f, MZFb, MZFbv, FactCorMZF, FactCovMZF, covMZF, TOTcovMZF, dataMZF, objMZF, fitFunction, StFLf, StErf, name="MZF" )
modelDZM	<-mxModel(pars1m, pars2m, DZMb, FactCorDZM, FactCovDZM, covDZM, TOTcovDZM, dataDZM, objDZM, fitFunction, name="DZM" )
modelDZF	<-mxModel(pars1f, pars2f, DZFb, FactCorDZF, FactCovDZF, covDZF, TOTcovDZF, dataDZF, objDZF, fitFunction, name="DZF" )
modelDZO	<-mxModel(pars1m, pars1f, pars1o, pars2m, pars2f, Threso, Ze5, DZOb, FactCorDZO, FactCovDZO, covDZO, TOTcovDZO, dataDZO, objDZO, fitFunction, name="DZO" )

minus2ll	<-mxAlgebra( expression=MZM.objective + MZF.objective + DZM.objective + DZF.objective + DZO.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )

cist1m	<-mxCI (c ('MZM.Rwithinm[2,1]','MZM.Rwithinm[3,1]','MZM.Rwithinm[4,1]','MZM.Rwithinm[5,1]','MZM.Rwithinm[3,2]','MZM.Rwithinm[4,2]','MZM.Rwithinm[5,2]','MZM.Rwithinm[4,3]','MZM.Rwithinm[5,3]','MZM.Rwithinm[5,4]'))
cist2m	<-mxCI (c ('MZM.StandFactm[1,1]','MZM.StandFactm[2,1]','MZM.StandFactm[3,1]','MZM.StandFactm[4,1]','MZM.StandFactm[5,1]','MZM.StandFactm[6,1]','MZM.StandFactm[7,1]'))
cist3m	<-mxCI (c ('MZM.StandErm[1,1]','MZM.StandErm[2,1]','MZM.StandErm[3,1]','MZM.StandErm[4,1]','MZM.StandErm[5,1]','MZM.StandErm[6,1]','MZM.StandErm[7,1]'))
cist4m	<-mxCI (c ('MZM.RbetweenMZM[1,1]','MZM.RbetweenMZM[2,1]','MZM.RbetweenMZM[3,1]','MZM.RbetweenMZM[4,1]','MZM.RbetweenMZM[5,1]','MZM.RbetweenMZM[2,2]','MZM.RbetweenMZM[3,2]','MZM.RbetweenMZM[4,2]','MZM.RbetweenMZM[5,2]','MZM.RbetweenMZM[3,3]','MZM.RbetweenMZM[4,3]','MZM.RbetweenMZM[5,3]','MZM.RbetweenMZM[4,4]','MZM.RbetweenMZM[5,4]','MZM.RbetweenMZM[5,5]'))
cist5m	<-mxCI (c ('DZM.RbetweenDZM[1,1]','DZM.RbetweenDZM[2,1]','DZM.RbetweenDZM[3,1]','DZM.RbetweenDZM[4,1]','DZM.RbetweenDZM[5,1]','DZM.RbetweenDZM[2,2]','DZM.RbetweenDZM[3,2]','DZM.RbetweenDZM[4,2]','DZM.RbetweenDZM[5,2]','DZM.RbetweenDZM[3,3]','DZM.RbetweenDZM[4,3]','DZM.RbetweenDZM[5,3]','DZM.RbetweenDZM[4,4]','DZM.RbetweenDZM[5,4]','DZM.RbetweenDZM[5,5]'))

cist1f	<-mxCI (c ('MZF.Rwithinf[2,1]','MZF.Rwithinf[3,1]','MZF.Rwithinf[4,1]','MZF.Rwithinf[5,1]','MZF.Rwithinf[3,2]','MZF.Rwithinf[4,2]','MZF.Rwithinf[5,2]','MZF.Rwithinf[4,3]','MZF.Rwithinf[5,3]','MZF.Rwithinf[5,4]'))
cist2f	<-mxCI (c ('MZF.StandFactf[1,1]','MZF.StandFactf[2,1]','MZF.StandFactf[3,1]','MZF.StandFactf[4,1]','MZF.StandFactf[5,1]','MZF.StandFactf[6,1]','MZF.StandFactf[7,1]'))
cist3f	<-mxCI (c ('MZF.StandErf[1,1]','MZF.StandErf[2,1]','MZF.StandErf[3,1]','MZF.StandErf[4,1]','MZF.StandErf[5,1]','MZF.StandErf[6,1]','MZF.StandErf[7,1]'))
cist4f	<-mxCI (c ('MZF.RbetweenMZF[1,1]','MZF.RbetweenMZF[2,1]','MZF.RbetweenMZF[3,1]','MZF.RbetweenMZF[4,1]','MZF.RbetweenMZF[5,1]','MZF.RbetweenMZF[2,2]','MZF.RbetweenMZF[3,2]','MZF.RbetweenMZF[4,2]','MZF.RbetweenMZF[5,2]','MZF.RbetweenMZF[3,3]','MZF.RbetweenMZF[4,3]','MZF.RbetweenMZF[5,3]','MZF.RbetweenMZF[4,4]','MZF.RbetweenMZF[5,4]','MZF.RbetweenMZF[5,5]'))
cist5f	<-mxCI (c ('DZF.RbetweenDZF[1,1]','DZF.RbetweenDZF[2,1]','DZF.RbetweenDZF[3,1]','DZF.RbetweenDZF[4,1]','DZF.RbetweenDZF[5,1]','DZF.RbetweenDZF[2,2]','DZF.RbetweenDZF[3,2]','DZF.RbetweenDZF[4,2]','DZF.RbetweenDZF[5,2]','DZF.RbetweenDZF[3,3]','DZF.RbetweenDZF[4,3]','DZF.RbetweenDZF[5,3]','DZF.RbetweenDZF[4,4]','DZF.RbetweenDZF[5,4]','DZF.RbetweenDZF[5,5]'))
cist5o	<-mxCI (c ('DZO.RbetweenDZO[1,1]','DZO.RbetweenDZO[2,1]','DZO.RbetweenDZO[3,1]','DZO.RbetweenDZO[4,1]','DZO.RbetweenDZO[5,1]','DZO.RbetweenDZO[2,2]','DZO.RbetweenDZO[3,2]','DZO.RbetweenDZO[4,2]','DZO.RbetweenDZO[5,2]','DZO.RbetweenDZO[3,3]','DZO.RbetweenDZO[4,3]','DZO.RbetweenDZO[5,3]','DZO.RbetweenDZO[4,4]','DZO.RbetweenDZO[5,4]','DZO.RbetweenDZO[5,5]'))
HetPhCModel	<-mxModel("HetPhC", modelMZM, modelMZF, modelDZM, modelDZF, modelDZO, minus2ll, obj, cist1m, cist2m, cist3m, cist4m, cist5m, cist1f, cist2f, cist3f, cist4f, cist5f, cist5o) 

# --------------------------------------------------------------------------------------------------------------------------------
# 2 RUN Phenotypic Fact Covariance Model by Zygosity

HetPhCFit	<-mxTryHardOrdinal(HetPhCModel, intervals=T)
(HetPhCSumm	<-summary(HetPhCFit))


# 2bi - Homogeneity model for phenotypic factor correlation model 
# Purpose is to test for sex differences
# -----------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------
# Submodel 2bi - Homogeneity model
# -----------------------------------------------------------------------------------------
HomoPhCModel	<- mxModel(HetPhCFit, name="HomoPhC")
HomoPhCModel	<- omxSetParameters(HomoPhCModel, labels=sdLabsm,free=c(F,T,T,T,T), values=1, newlabels=sdLabsf)
HomoPhCModel	<- omxSetParameters(HomoPhCModel, labels=rphvLabsm,free=T, values=StWithinpersonv, newlabels=rphvLabsf)
HomoPhCModel	<- omxSetParameters(HomoPhCModel, labels=rphPLabsm,free=T, values=StWithinpersonP, newlabels=rphPLabsf)
HomoPhCModel	<- omxSetParameters(HomoPhCModel, labels=rphvPLabsm,free=T, values=StWithinpersonvP, newlabels=rphvPLabsf)
HomoPhCModel	<- omxSetParameters(HomoPhCModel, labels=MZMbLabs,free=T, values=StBetweenMZM, newlabels=MZFbLabs)
HomoPhCModel	<- omxSetParameters(HomoPhCModel, labels=DZMbLabs,free=PatBetweenDZ, values=StBetweenDZM, newlabels=DZFbLabs)
HomoPhCModel	<- omxSetParameters(HomoPhCModel, labels=DZObLabs,free=PatBetweenDZ, values=StBetweenDZM, newlabels=DZFbLabs)
HomoPhCModel	<- omxSetParameters(HomoPhCModel, labels=LabMm,free=PatBetweenDZ, values=Stmeanf, newlabels=LabMf)
HomoPhCModel	<- omxAssignFirstParameters(HomoPhCModel)

HomoPhCFit	<- mxTryHard(HomoPhCModel, intervals=F)
(HomoPhCSum	<- summary(HomoPhCFit))

mxCompare(HetPhCFit,HomoPhCFit)


# Generate confidence intervals
HetPhCCIModel	<-mxModel(HetPhCModel)
HetPhCCIFit		<-mxRun(HetPhCCIModel, intervals=TRUE)
(HetPhCCISumm	<-summary(HetPhCCIFit, verbose=T))

capture.output(print(HetPhCCISumm,row.names=F), file = "summary.HetPhCModel", append = FALSE)

mxEval(MZM.Rwithinm, HetPhCFit)
mxEval(MZM.RbetweenMZM, HetPhCFit)
mxEval(DZM.RbetweenDZM, HetPhCFit)

mxEval(MZM.expFactCovMZM, HetPhCFit)
mxEval(DZM.expFactCovDZM, HetPhCFit)

mxEval(MZM.ExpCovMZM, HetPhCFit)
mxEval(DZM.ExpCovDZM, HetPhCFit)

mxEval(MZM.TOTexpCovMZM, HetPhCFit)
mxEval(DZM.TOTexpCovDZM, HetPhCFit)

mxEval(MZM.FactLm, HetPhCFit)
mxEval(MZM.StandFactm, HetPhCFit)

mxEval(MZM.Errorm, HetPhCFit)
mxEval(MZM.StandErm, HetPhCFit)



mxEval(MZF.Rwithinf, HetPhCFit)
mxEval(MZF.RbetweenMZF, HetPhCFit)
mxEval(DZF.RbetweenDZF, HetPhCFit)

mxEval(MZF.expFactCovMZF, HetPhCFit)
mxEval(DZF.expFactCovDZF, HetPhCFit)

mxEval(MZF.ExpCovMZF, HetPhCFit)
mxEval(DZF.ExpCovDZF, HetPhCFit)

mxEval(MZF.TOTexpCovMZF, HetPhCFit)
mxEval(DZF.TOTexpCovDZF, HetPhCFit)

mxEval(MZF.FactLf, HetPhCFit)
mxEval(MZF.StandFactf, HetPhCFit)

mxEval(MZF.Errorf, HetPhCFit)
mxEval(MZF.StandErf, HetPhCFit)


mxEval(DZO.RbetweenDZO, HetPhCFit)
mxEval(DZO.expFactCovDZO, HetPhCFit)
mxEval(DZO.ExpCovDZO, HetPhCFit)
mxEval(DZO.TOTexpCovDZO, HetPhCFit)

####

#*******************************************************************************************************
# __(IIb-Sex Diff)_________________________________________________________________________________________________
# Phenotypic Covariance Model across Latent Constructs  including ELA, CGN and Vict with sex differences
# Restrictions: means and variances equated across birth-order & zygosity groups;
# One set of factor loadings; one set of correltional paths between the factors; one set of error terms
# We estimate the factor variances, giving them a scale by fixing the loading on the 1st variable to 1
# This model specifies a full var/cov structure between the latent factors for MZ and DZ twins 
#______________________________________________________________________________________________________

nv			<- 6				# number of variables for a twin = 1 in Univariate
nvo 			<- 1     			# number of ordinal variables per twin
nvc 			<- nv-nvo  			# number of continuous variables per twin
poso 			<- nvo 			# position where ordinal variables start
ntv			<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nth			<- 4				# number of max thresholds
nlower		<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor			<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
ninc 			<- nth-1 			# number of max increments
ncovariates 	<- 2 				# number of covariates

nfact			<- 5				# number of Latent Factors for Mediation Model per twin
nfact2		<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nfcor			<-(nfact*(nfact+1)/2)-nfact	# number of free elements in a correlation matrix nfact*nfcat

Groups		<- c("mzm","mzf","dzm","dzf","dzo") 
Vars			<- c('SO','CGN4','ELA','Dep','Anx','Vict')
selVars		<- c('SO1','CGN41','ELA1','Dep1','Anx1','Vict1',
			     'SO2','CGN42','ELA2','Dep2','Anx2','Vict2')
useVars		<- c('SO1','CGN41','ELA1','Dep1','Anx1','Vict1',
			     'SO2','CGN42','ELA2','Dep2','Anx2','Vict2','age1','age2')

mzmData		<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1) , useVars)
mzfData		<- subset(TWINdata2, zyg1%in%c(3)|zyg2%in%c(3) , useVars)
dzmData		<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2) , useVars)
dzfData		<- subset(TWINdata2, zyg1%in%c(4)|zyg2%in%c(4) , useVars)
dzoData		<- subset(TWINdata2, zyg1%in%c(5)|zyg2%in%c(5) , useVars)

psych::describe(mzmData)
psych::describe(mzfData)
psych::describe(dzmData)
psych::describe(dzfData)
psych::describe(dzoData)

mzmData$SO1[is.na(mzmData$age1)] <- NA
mzmData$SO2[is.na(mzmData$age2)] <- NA
mzfData$SO1[is.na(mzfData$age1)] <- NA
mzfData$SO2[is.na(mzfData$age2)] <- NA
dzmData$SO1[is.na(dzmData$age1)] <- NA
dzmData$SO2[is.na(dzmData$age2)] <- NA
dzfData$SO1[is.na(dzfData$age1)] <- NA
dzfData$SO2[is.na(dzfData$age2)] <- NA
dzoData$SO1[is.na(dzoData$age1)] <- NA
dzoData$SO2[is.na(dzoData$age2)] <- NA

mzmData$age1[is.na(mzmData$age1)] <- 999
mzmData$age2[is.na(mzmData$age2)] <- 999
mzfData$age1[is.na(mzfData$age1)] <- 999
mzfData$age2[is.na(mzfData$age2)] <- 999
dzmData$age1[is.na(dzmData$age1)] <- 999
dzmData$age2[is.na(dzmData$age2)] <- 999
dzfData$age1[is.na(dzfData$age1)] <- 999
dzfData$age2[is.na(dzfData$age2)] <- 999
dzoData$age1[is.na(dzoData$age1)] <- 999
dzoData$age2[is.na(dzoData$age2)] <- 999

psych::describe(mzmData)
psych::describe(mzfData)
psych::describe(dzmData)
psych::describe(dzfData)
psych::describe(dzoData)


# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)

(Stmeanm	<-colMeans(mzmData[,2:nv],na.rm=TRUE))
StMm 		<-c(0, Stmeanm, 0, Stmeanm)
(LabMm	<- paste("mm",2:nv,sep=""))
MLabsm	<-c(NA,LabMm,NA,LabMm) 

(Stmeanf	<-colMeans(mzfData[,2:nv],na.rm=TRUE))
StMf 		<-c(0, Stmeanf, 0, Stmeanf)
(LabMf	<- paste("mf",2:nv,sep=""))
MLabsf	<-c(NA,LabMf,NA,LabMf) 

StMo 		<-c(0, Stmeanm, 0, Stmeanf)
MLabso	<-c(NA,LabMm,NA,LabMf) 

(LabErm	<-c("e1m","e2m","e3m","e4m","e4m","e6m"))
(LabErf	<-c("e1f","e2f","e3f","e4f","e4f","e6f"))

# Create Labels for the Factor parameters
(sdLabsm	<- paste("sdm",1:nfact,sep=""))	# SD
(rphLabsm	<- paste("rvm",1:nfcor,sep="")) #Within person 
(MZMbLabs 	<- paste("rmzm", do.call(c, sapply(seq(1, nfact), function(x){ paste(x:nfact, x,sep="") })), sep="")) # all vars in mz
(DZMbLabs 	<- paste("rdzm", do.call(c, sapply(seq(1, nfact), function(x){ paste(x:nfact, x,sep="") })), sep="")) # all vars in dz

(sdLabsf	<- paste("sdf",1:nfact,sep=""))	# SD
(rphLabsf	<- paste("rvf",1:nfcor,sep="")) #Within person 
(MZFbLabs 	<- paste("rmzf", do.call(c, sapply(seq(1, nfact), function(x){ paste(x:nfact, x,sep="") })), sep="")) # within person
(DZFbLabs 	<- paste("rdzf", do.call(c, sapply(seq(1, nfact), function(x){ paste(x:nfact, x,sep="") })), sep="")) # all vars in dz
(DZObLabs 	<- paste("rdzo", do.call(c, sapply(seq(1, nfact), function(x){ paste(x:nfact, x,sep="") })), sep="")) # all vars in dz

(StWithinperson 	<-c(.2))

(StBetweenMZM  	<-c(.5,.2,.2,.2,.2,.5,.2,.2,.2,.5,.2,.2,.5,.2,.5))
(StBetweenDZM  	<-c(.3,.1,.1,.1,.1,.3,.1,.1,.1,.3,.1,.1,.3,.1,.3))
(StBetweenMZF  	<-c(.5,.2,.2,.2,.2,.5,.2,.2,.2,.5,.2,.2,.5,.2,.5))
(StBetweenDZF  	<-c(.3,.1,.1,.1,.1,.3,.1,.1,.1,.3,.1,.1,.3,.1,.3))

(PatBetweenDZ  	<-c(T))

# Create Labels for the Factor Loadings (1st loadings fixed to 1)

PatFl	<- c(F,F,F,F,F,F,			
	     F,F,F,F,F,F,
	     F,F,F,F,F,F,
	     F,F,F,T,F,F,
	     F,F,F,F,F,F)

StFl	<- c(1,0,0,0,0,0,
	     0,1,0,0,0,0,
	     0,0,1,0,0,0,
	     0,0,0,.5,1,0,
	     0,0,0,0,0,1)

LabFlm	<- c('l1m',NA,NA,NA,NA,NA,
		      NA,'l2m',NA,NA,NA,NA,
		      NA,NA,'l3m',NA,NA,NA,
		      NA,NA,NA,'l4m','l5m',NA,
		      NA,NA,NA,NA,NA,'l6m')

LabFlf	<- c('l1f',NA,NA,NA,NA,NA,
		      NA,'l2f',NA,NA,NA,NA,
		      NA,NA,'l3f',NA,NA,NA,
		      NA,NA,NA,'l4f','l5f',NA,
		      NA,NA,NA,NA,NA,'l6f')

# Free parameters
(Pat  	<- c( rep(FALSE,nvo), rep(TRUE, nvc)))

# ______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
# ______________________________________________________________________________________________________

Meanm	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(Pat, Pat), values=c(StMm), labels=c(MLabsm), name="expmm" ) #first variable is not free to be estimated
Meanf	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(Pat, Pat), values=c(StMf), labels=c(MLabsf), name="expmf" ) 
Meano	<-mxMatrix( type="Full", nrow=1, ncol=ntv, free=c(Pat, Pat), values=c(StMo), labels=c(MLabso), name="expmo" ) 

# I constrain the threshold and increments to be equal across birth order or zygosity and specify only one set
LabThm	<-c('T_1m','i_11m','i_12m','i_13m')	# THs for var 1 for mz twin 1 (mzm)
LabThf	<-c('T_1f','i_11f','i_12f','i_13f')	# THs for var 1 for mz twin 1 (mzf)

LabCovAm	<-c('BageThSOm','BageThSOm','BageThSOm','BageThSOm')
LabCovAf	<-c('BageThSOf','BageThSOf','BageThSOf','BageThSOf')

ThPat		<-c(TRUE,TRUE,TRUE,TRUE)
StTHm		<-c(0.84,0.57,0.14,0.20)
StTHf		<-c(0.84,0.57,0.14,0.20)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

# effect of age and sex on ordinal variable
betaAm	<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.05, labels=LabCovAm, name="BageTHm" )
betaAf	<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.05, labels=LabCovAf, name="BageTHf" )
 
# thresholds
Thrm		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTHm, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabThm, name="Thm")
incm		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=F, values=1, name="Lowm")

Thrf		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTHf, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabThf, name="Thf")
incf		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=F, values=1, name="Lowf")

Thresm	<-mxAlgebra( expression= cbind(Lowm%*%Thm + BageTHm%x%age1, Lowm%*%Thm + BageTHm%x%age2), name="expThresm")
Thresf	<-mxAlgebra( expression= cbind(Lowf%*%Thf + BageTHf%x%age1, Lowf%*%Thf + BageTHf%x%age2), name="expThresf")
Threso	<-mxAlgebra( expression= cbind(Lowm%*%Thm + BageTHm%x%age1, Lowf%*%Thf + BageTHf%x%age2), name="expThreso")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Loadm		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFlm, name="FactLm" )
Loadf		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFlf, name="FactLf" )
Ze65		<-mxMatrix("Zero", nv, nfact, free=F, name="Z65")
LoadTwm	<-mxAlgebra(rbind(cbind(FactLm,Z65), cbind(Z65, FactLm)), name="FactLTwm")
LoadTwf	<-mxAlgebra(rbind(cbind(FactLf,Z65), cbind(Z65, FactLf)), name="FactLTwf")
LoadTwo	<-mxAlgebra(rbind(cbind(FactLm,Z65), cbind(Z65, FactLf)), name="FactLTwo")

ErPathm	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=c(F,F,F,T,T,F), values=c(0,0,0,.5,.5,0), labels=LabErm, name="Erpm" )
ErPathf	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=c(F,F,F,T,T,F), values=c(0,0,0,.5,.5,0), labels=LabErf, name="Erpf" )

Erm		<-mxAlgebra(Erpm %*% t(Erpm), name="Errorm")
Erf		<-mxAlgebra(Erpf %*% t(Erpf), name="Errorf")

Ze6		<-mxMatrix("Zero", nv, nv, free=F, name="Z6")
ErTwm		<-mxAlgebra(rbind(cbind(Errorm,Z6), cbind(Z6, Errorm)), name="ErrorTwm")
ErTwf		<-mxAlgebra(rbind(cbind(Errorf,Z6), cbind(Z6, Errorf)), name="ErrorTwf")
ErTwo		<-mxAlgebra(rbind(cbind(Errorm,Z6), cbind(Z6, Errorf)), name="ErrorTwo")
 
# elements for the SD of Factors
Id2	<-mxMatrix("Iden", 2, 2, free=F, name="I2")
sdFm	<-mxMatrix("Diag", nfact, nfact, free=c(F,T,T,T,T), values=1, labels=sdLabsm, name="SDfm") 
sdFf	<-mxMatrix("Diag", nfact, nfact, free=c(F,T,T,T,T), values=1, labels=sdLabsf, name="SDff") 

sdFTwm	<-mxAlgebra(I2 %x% SDfm, name="SDftwinm")
sdFTwf	<-mxAlgebra(I2 %x% SDff, name="SDftwinf")
Ze5		<-mxMatrix("Zero", nfact, nfact, free=F, name="Z5")
sdFTwo	<-mxAlgebra(rbind(cbind(SDfm,Z5), cbind(Z5, SDff)), name="SDftwino")

# elements for the correlations of Factors
Rphm		<-mxMatrix("Stand", nfact, nfact, free = TRUE, values = StWithinperson, labels=rphLabsm, lbound=-.999, ubound=.999, name="Rwithinm") 
Rphf		<-mxMatrix("Stand", nfact, nfact, free = TRUE, values = StWithinperson, labels=rphLabsf, lbound=-.999, ubound=.999, name="Rwithinf") 

MZMb		<-mxMatrix("Symm", nfact, nfact, free = TRUE, values = StBetweenMZM, labels=MZMbLabs, lbound=-.999, ubound=.999, name="RbetweenMZM") 
MZFb		<-mxMatrix("Symm", nfact, nfact, free = TRUE, values = StBetweenMZF, labels=MZFbLabs, lbound=-.999, ubound=.999, name="RbetweenMZF") 

DZMb		<-mxMatrix("Symm", nfact, nfact, free=PatBetweenDZ, values=StBetweenDZM, labels=DZMbLabs, lbound=-.999, ubound=.999, name="RbetweenDZM") 
DZFb		<-mxMatrix("Symm", nfact, nfact, free=PatBetweenDZ, values=StBetweenDZF, labels=DZFbLabs, lbound=-.999, ubound=.999, name="RbetweenDZF") 
DZOb		<-mxMatrix("Symm", nfact, nfact, free=PatBetweenDZ, values=StBetweenDZF, labels=DZObLabs, lbound=-.999, ubound=.999, name="RbetweenDZO") 

FactCorMZM	<-mxAlgebra(rbind(cbind(Rwithinm,RbetweenMZM), cbind(RbetweenMZM, Rwithinm)), name="RMZM")
FactCorMZF	<-mxAlgebra(rbind(cbind(Rwithinf,RbetweenMZF), cbind(RbetweenMZF, Rwithinf)), name="RMZF")
FactCorDZM	<-mxAlgebra(rbind(cbind(Rwithinm,RbetweenDZM), cbind(RbetweenDZM, Rwithinm)), name="RDZM")
FactCorDZF	<-mxAlgebra(rbind(cbind(Rwithinf,RbetweenDZF), cbind(RbetweenDZF, Rwithinf)), name="RDZF")
FactCorDZO	<-mxAlgebra(rbind(cbind(Rwithinm,RbetweenDZO), cbind(RbetweenDZO, Rwithinf)), name="RDZO")

# Generate expected Covariance matrices of Factors
FactCovMZM	<-mxAlgebra(SDftwinm %&% RMZM , name="expFactCovMZM")
FactCovMZF	<-mxAlgebra(SDftwinf %&% RMZF , name="expFactCovMZF")
FactCovDZM	<-mxAlgebra(SDftwinm %&% RDZM , name="expFactCovDZM")
FactCovDZF	<-mxAlgebra(SDftwinf %&% RDZF , name="expFactCovDZF")
FactCovDZO	<-mxAlgebra(SDftwino %&% RDZO , name="expFactCovDZO")

## This second step then derives the var/cov matrix of the observed/measured variables in terms of the variance/covariances of the latent factors and the Factor Loadings
covMZM	<-mxAlgebra( expression= FactLTwm  %&% expFactCovMZM , name="ExpCovMZM" )
covMZF	<-mxAlgebra( expression= FactLTwf  %&% expFactCovMZF , name="ExpCovMZF" )
covDZM	<-mxAlgebra( expression= FactLTwm  %&% expFactCovDZM , name="ExpCovDZM" )
covDZF	<-mxAlgebra( expression= FactLTwf  %&% expFactCovDZF , name="ExpCovDZF" )
covDZO	<-mxAlgebra( expression= FactLTwo  %&% expFactCovDZO , name="ExpCovDZO" )

## Finally, we derive the total expected variance/covariances for the measured variables which go in the models
TOTcovMZM	<-mxAlgebra( expression= ExpCovMZM + ErrorTwm , name="TOTexpCovMZM" )
TOTcovMZF	<-mxAlgebra( expression= ExpCovMZF + ErrorTwf , name="TOTexpCovMZF" )
TOTcovDZM	<-mxAlgebra( expression= ExpCovDZM + ErrorTwm , name="TOTexpCovDZM" )
TOTcovDZF	<-mxAlgebra( expression= ExpCovDZF + ErrorTwf , name="TOTexpCovDZF" )
TOTcovDZO	<-mxAlgebra( expression= ExpCovDZO + ErrorTwo , name="TOTexpCovDZO" )

# Standardizing parameters **********************

# Standardized Factor Loadings
StFLm		<-mxAlgebra( expression= sqrt(diag2vec( FactLm %&% expFactCovMZM[1:5,1:5] / TOTexpCovMZM[1:6,1:6])), name="StandFactm")
StFLf		<-mxAlgebra( expression= sqrt(diag2vec( FactLf %&% expFactCovMZF[1:5,1:5] / TOTexpCovMZF[1:6,1:6])), name="StandFactf")

# Standardise error terms of the measured variables
StErm		<-mxAlgebra( expression= sqrt(diag2vec( Errorm/TOTexpCovMZM[1:6,1:6])), name="StandErm")
StErf		<-mxAlgebra( expression= sqrt(diag2vec( Errorf/TOTexpCovMZF[1:6,1:6])), name="StandErf")

# ************************************

# Data objects for Multiple Groups
dataMZM	<- mxData( observed=mzmData, type="raw" )
dataMZF	<- mxData( observed=mzfData, type="raw" )
dataDZM	<- mxData( observed=dzmData, type="raw" )
dataDZF	<- mxData( observed=dzfData, type="raw" )
dataDZO	<- mxData( observed=dzoData, type="raw" )

# Objective objects for Multiple Groups
objMZM	<- mxExpectationNormal( covariance="TOTexpCovMZM", means="expmm", dimnames=selVars, thresholds="expThresm", threshnames=c("SO1","SO2"))
objMZF	<- mxExpectationNormal( covariance="TOTexpCovMZF", means="expmf", dimnames=selVars, thresholds="expThresf", threshnames=c("SO1","SO2"))
objDZM	<- mxExpectationNormal( covariance="TOTexpCovDZM", means="expmm", dimnames=selVars, thresholds="expThresm", threshnames=c("SO1","SO2"))
objDZF	<- mxExpectationNormal( covariance="TOTexpCovDZF", means="expmf", dimnames=selVars, thresholds="expThresf", threshnames=c("SO1","SO2"))
objDZO	<- mxExpectationNormal( covariance="TOTexpCovDZO", means="expmo", dimnames=selVars, thresholds="expThreso", threshnames=c("SO1","SO2"))

fitFunction <- mxFitFunctionML()

# Combine Groups
pars1m	<-list(Meanm, Loadm, Ze65, Ze6, LoadTwm, ErPathm, Erm, ErTwm, Id2, sdFm, sdFTwm, Rphm)
pars1f	<-list(Meanf, Loadf, Ze65, Ze6, LoadTwf, ErPathf, Erf, ErTwf, Id2, sdFf, sdFTwf, Rphf)
pars1o	<-list(Meano, LoadTwo, ErTwo, sdFTwo)
pars2m	<-list(obsage1, obsage2, betaAm, Thrm, incm, Thresm)
pars2f	<-list(obsage1, obsage2, betaAf, Thrf, incf, Thresf)

modelMZM	<-mxModel(pars1m, pars2m, MZMb, FactCorMZM, FactCovMZM, covMZM, TOTcovMZM, dataMZM, objMZM, fitFunction, StFLm, StErm, name="MZM" )
modelMZF	<-mxModel(pars1f, pars2f, MZFb, FactCorMZF, FactCovMZF, covMZF, TOTcovMZF, dataMZF, objMZF, fitFunction, StFLf, StErf, name="MZF" )
modelDZM	<-mxModel(pars1m, pars2m, DZMb, FactCorDZM, FactCovDZM, covDZM, TOTcovDZM, dataDZM, objDZM, fitFunction, name="DZM" )
modelDZF	<-mxModel(pars1f, pars2f, DZFb, FactCorDZF, FactCovDZF, covDZF, TOTcovDZF, dataDZF, objDZF, fitFunction, name="DZF" )
modelDZO	<-mxModel(pars1m, pars1f, pars1o, pars2m, pars2f, Threso, Ze5, DZOb, FactCorDZO, FactCovDZO, covDZO, TOTcovDZO, dataDZO, objDZO, fitFunction, name="DZO" )

minus2ll	<-mxAlgebra( expression=MZM.objective + MZF.objective + DZM.objective + DZF.objective + DZO.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )

cist1m	<-mxCI (c ('MZM.Rwithinm[2,1]','MZM.Rwithinm[3,1]','MZM.Rwithinm[4,1]','MZM.Rwithinm[5,1]','MZM.Rwithinm[3,2]','MZM.Rwithinm[4,2]','MZM.Rwithinm[5,2]','MZM.Rwithinm[4,3]','MZM.Rwithinm[5,3]','MZM.Rwithinm[5,4]'))
#cist2m	<-mxCI (c ('MZM.StandFactm[1,1]','MZM.StandFactm[2,1]','MZM.StandFactm[3,1]','MZM.StandFactm[4,1]','MZM.StandFactm[5,1]','MZM.StandFactm[6,1]'))
#cist3m	<-mxCI (c ('MZM.StandErm[1,1]','MZM.StandErm[2,1]','MZM.StandErm[3,1]','MZM.StandErm[4,1]','MZM.StandErm[5,1]','MZM.StandErm[6,1]'))

cist1f	<-mxCI (c ('MZF.Rwithinf[2,1]','MZF.Rwithinf[3,1]','MZF.Rwithinf[4,1]','MZF.Rwithinf[5,1]','MZF.Rwithinf[3,2]','MZF.Rwithinf[4,2]','MZF.Rwithinf[5,2]','MZF.Rwithinf[4,3]','MZF.Rwithinf[5,3]','MZF.Rwithinf[5,4]'))
#cist2f	<-mxCI (c ('MZF.StandFactf[1,1]','MZF.StandFactf[2,1]','MZF.StandFactf[3,1]','MZF.StandFactf[4,1]','MZF.StandFactf[5,1]','MZF.StandFactf[6,1]'))
#cist3f	<-mxCI (c ('MZF.StandErf[1,1]','MZF.StandErf[2,1]','MZF.StandErf[3,1]','MZF.StandErf[4,1]','MZF.StandErf[5,1]','MZF.StandErf[6,1]'))

HetPhCbModel	<-mxModel("HetPhCb", modelMZM, modelMZF, modelDZM, modelDZF, modelDZO, minus2ll, obj, cist1m, cist1f) 

# --------------------------------------------------------------------------------------------------------------------------------
# 2b RUN Phenotypic Fact Covariance Model by Zygosity

HetPhCbFit	<-mxTryHardOrdinal(HetPhCbModel, intervals=T)
(HetPhCbSumm	<-summary(HetPhCbFit))

####

#****************************************************************************************************************************
# __(IVa-Sex diff)______________________________________________________________________________________________
# ACE Factor MODEL (including victimisation) by zygosity
# NO causal paths between Phenotypic Factors; A, C and E latent factors have Cholesky Structure
# + Asp, Csp and Esp in the bottom with constraints to identify the model on top
# Correlation between Phenotypic Factors only due to shared A, C and E influences
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1)
#_____________________________________________________________________________________________________________________________


nv		<- 4				# number of variables for a twin = 1 in Univariate
nvo 		<- 1     			# number of ordinal variables per twin
nvc 		<- nv-nvo  			# number of continuous variables per twin
poso 		<- nvo 			# position where ordinal variables start
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nth		<- 4				# number of max thresholds
nlower		<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
ninc 		<- nth-1 			# number of max increments
ncovariates 	<- 2 				# number of covariates

nfact		<- 3				# number of Latent Factors for Mediation Model per twin
nfact2		<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nfcor		<-(nfact*(nfact+1)/2)-nfact	# number of free elements in a correlation matrix nfact*nfcat

Groups		<- c("mzm", "dzm", "mzf", "dzf", "dzo")

Vars		<- c('SO','Dep','Anx','Vict')
selVars	<- c('SO1','Dep1','Anx1','Vict1',
		     'SO2','Dep2','Anx2','Vict2')
useVars	<- c('SO1','Dep1','Anx1','Vict1',
		     'SO2','Dep2','Anx2','Vict2','age1','sex1','age2','sex2')

mzmData	<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1) , useVars)
mzfData	<- subset(TWINdata2, zyg1%in%c(3)|zyg2%in%c(3) , useVars)
dzmData	<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2) , useVars)
dzfData	<- subset(TWINdata2, zyg1%in%c(4)|zyg2%in%c(4) , useVars)
dzoData	<- subset(TWINdata2, zyg1%in%c(5)|zyg2%in%c(5) , useVars)

psych::describe(mzmData)
psych::describe(mzfData)
psych::describe(dzmData)
psych::describe(dzfData)
psych::describe(dzoData)

mzmData$SO1[is.na(mzmData$age1)] <- NA
mzmData$SO2[is.na(mzmData$age2)] <- NA
mzfData$SO1[is.na(mzfData$age1)] <- NA
mzfData$SO2[is.na(mzfData$age2)] <- NA
dzmData$SO1[is.na(dzmData$age1)] <- NA
dzmData$SO2[is.na(dzmData$age2)] <- NA
dzfData$SO1[is.na(dzfData$age1)] <- NA
dzfData$SO2[is.na(dzfData$age2)] <- NA
dzoData$SO1[is.na(dzoData$age1)] <- NA
dzoData$SO2[is.na(dzoData$age2)] <- NA

mzmData$age1[is.na(mzmData$age1)] <- 999
mzmData$age2[is.na(mzmData$age2)] <- 999
mzfData$age1[is.na(mzfData$age1)] <- 999
mzfData$age2[is.na(mzfData$age2)] <- 999
dzmData$age1[is.na(dzmData$age1)] <- 999
dzmData$age2[is.na(dzmData$age2)] <- 999
dzfData$age1[is.na(dzfData$age1)] <- 999
dzfData$age2[is.na(dzfData$age2)] <- 999
dzoData$age1[is.na(dzoData$age1)] <- 999
dzoData$age2[is.na(dzoData$age2)] <- 999

psych::describe(mzmData)
psych::describe(mzfData)
psych::describe(dzmData)
psych::describe(dzfData)
psych::describe(dzoData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabsm	<- paste("mm",1:nv,sep=""))
(mLabsf	<- paste("mf",1:nv,sep=""))
(Stmeanm	<- colMeans(mzmData[,2:nv],na.rm=TRUE))
(Stmeanf	<- colMeans(mzfData[,2:nv],na.rm=TRUE))
(Stsdm 	<- sapply(mzmData[,2:nv],sd, na.rm=TRUE))
(Stsdf 	<- sapply(mzfData[,2:nv],sd, na.rm=TRUE))
(PatM		<- c(F,TRUE,TRUE,TRUE))

# Create Labels for Diagonal Matrices
# To identify this model we need to equate the sp effects of var and 2 and fix the Sp of last variable to 0)
(LabEsm	<- c(NA,'es2m','es2m',NA))
(LabAsm	<- c(NA,'as2m','as2m',NA))
(LabCsm	<- c(NA,'cs2m','cs2m',NA))

(LabEsf	<- c(NA,'es2f','es2f',NA))
(LabAsf	<- c(NA,'as2f','as2f',NA))
(LabCsf	<- c(NA,'cs2f','cs2f',NA))

PatSp		<- c(F,TRUE,TRUE,F)
StSpam		<- c(0,.5,.5,0)
StSpcm		<- c(0,.2,.2,0)
StSpem		<- c(0,.5,.5,0)

StSpaf		<- c(0,.5,.5,0)
StSpcf		<- c(0,.2,.2,0)
StSpef		<- c(0,.5,.5,0)

# all 1st loadings fixed to 1
PatFl		<- c(F,F,F,F,			
		     F,F,T,F,
		     F,F,F,F)

StFlm		<- c(1,0,0,0,
		     0,1,.5,0,
		     0,0,0,1)

StFlf		<- c(1,0,0,0,
		     0,1,.5,0,
		     0,0,0,1)

LabFlm	<- c('l1m',NA,NA,NA,
		     NA,'l2m','l3m',NA,
		     NA,NA,NA,'l4m')

LabFlf	<- c('l1f',NA,NA,NA,
		     NA,'l2f','l3f',NA,
		     NA,NA,NA,'l4f')

# ______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
# ______________________________________________________________________________________________________

Meansm		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(0,Stmeanm,0,Stmeanm), labels=c(mLabsm,mLabsm), name="expMeanm") 
Meansf		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(0,Stmeanf,0,Stmeanf), labels=c(mLabsf,mLabsf), name="expMeanf") 
Meanso		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(0,Stmeanm,0,Stmeanf), labels=c(mLabsm,mLabsf), name="expMeano") 

# Threshold and covariates
LabThm	<-c('T_1m','i_11m','i_12m','i_13m')	# THs for var 1 
LabThf	<-c('T_1f','i_11f','i_12f','i_13f')	# THs for var 1 

LabCovAm	<-c('BageThSOm','BageThSOm','BageThSOm','BageThSOm')
LabCovAf	<-c('BageThSOf','BageThSOf','BageThSOf','BageThSOf')

ThPat		<-c(T,T,T,T)
StTHm		<-c(1.4,.5,.05,.2)
StTHf		<-c(.1,1,4,.5)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

# effect of age and sex on ordinal variable
betaAm	<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.02, labels=LabCovAm, name="BageTHm" )
betaAf	<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.03, labels=LabCovAf, name="BageTHf" )
 
# thresholds
Thrm		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTHm, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabThm, name="Thm")
Thrf		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTHf, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabThf, name="Thf")

inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")

Thresm	<-mxAlgebra( expression= cbind(Low%*%Thm + BageTHm%x%age1, Low%*%Thm + BageTHm%x%age2), name="expThresm")
Thresf	<-mxAlgebra( expression= cbind(Low%*%Thf + BageTHf%x%age1, Low%*%Thf + BageTHf%x%age2), name="expThresf")
Threso	<-mxAlgebra( expression= cbind(Low%*%Thm + BageTHm%x%age1, Low%*%Thf + BageTHf%x%age2), name="expThreso")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Loadm		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFlm, labels=LabFlm, name="FactLm" )
Loadf		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFlf, labels=LabFlf, name="FactLf" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
Ze43		<-mxMatrix("Zero", nv, nfact, free=F, name="Z43")
LoadTwm	<-mxAlgebra(I2%x%FactLm, name="FactLTwm")
LoadTwf	<-mxAlgebra(I2%x%FactLf, name="FactLTwf")
LoadTwo	<-mxAlgebra(rbind(cbind(FactLm,Z43), cbind(Z43, FactLf)), name="FactLTwo")
  
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  

# Define the matrix to hold the A and C effects: Specific 
PathsAsm	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSp, values=StSpam, labels=LabAsm, name="asm" )
PathsCsm	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSp, values=StSpcm, labels=LabCsm, name="csm" )
PathsEsm	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSp, values=StSpem, labels=LabEsm, name="esm" )
PathsAsf	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSp, values=StSpaf, labels=LabAsf, name="asf" )
PathsCsf	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSp, values=StSpcf, labels=LabCsf, name="csf" )
PathsEsf	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSp, values=StSpef, labels=LabEsf, name="esf" )

covAsm		<-mxAlgebra( expression= asm %*% t(asm), name="Asm" )
covCsm		<-mxAlgebra( expression= csm %*% t(csm), name="Csm" )
covEsm		<-mxAlgebra( expression= esm %*% t(esm), name="Esm" )
covPsm		<-mxAlgebra( expression= Asm+Csm+Esm, name="Vsm" )

covAsf		<-mxAlgebra( expression= asf %*% t(asf), name="Asf" )
covCsf		<-mxAlgebra( expression= csf %*% t(csf), name="Csf" )
covEsf		<-mxAlgebra( expression= esf %*% t(esf), name="Esf" )
covPsf		<-mxAlgebra( expression= Asf+Csf+Esf, name="Vsf" )

covAsmf	<-mxAlgebra( expression= asm %*% t(asf), name="Asmf" )
covCsmf	<-mxAlgebra( expression= csm %*% t(csf), name="Csmf" )
covAsfm	<-mxAlgebra( expression= asf %*% t(asm), name="Asfm" )
covCsfm	<-mxAlgebra( expression= csf %*% t(csm), name="Csfm" )

# Define the matrices to hold the A and C effects: Common 
PathsAcm	<-mxMatrix(type="Lower", nrow=nfact, ncol=nfact, free=TRUE, values=.6, labels=c("a11m","a21m","a31m","a22m","a32m","a33m"), name="a_cm" )
PathsCcm	<-mxMatrix(type="Lower", nrow=nfact, ncol=nfact, free=TRUE, values=.3, labels=c("c11m","c21m","c31m","c22m","c32m","c33m"), name="c_cm" )
PathsEcm	<-mxMatrix(type="Lower", nrow=nfact, ncol=nfact, free=TRUE, values=.6, labels=c("e11m","e21m","e31m","e22m","e32m","e33m"), name="e_cm" )
PathsAcf	<-mxMatrix(type="Lower", nrow=nfact, ncol=nfact, free=TRUE, values=.6, labels=c("a11f","a21f","a31f","a22f","a32f","a33f"), name="a_cf" )
PathsCcf	<-mxMatrix(type="Lower", nrow=nfact, ncol=nfact, free=TRUE, values=.3, labels=c("c11f","c21f","c31f","c22f","c32f","c33f"), name="c_cf" )
PathsEcf	<-mxMatrix(type="Lower", nrow=nfact, ncol=nfact, free=TRUE, values=.6, labels=c("e11f","e21f","e31f","e22f","e32f","e33f"), name="e_cf" )

covAcm	<-mxAlgebra( expression= a_cm %*% t(a_cm), name="Acm" )
covCcm	<-mxAlgebra( expression= c_cm %*% t(c_cm), name="Ccm" )
covEcm	<-mxAlgebra( expression= e_cm %*% t(e_cm), name="Ecm" )
covPcm	<-mxAlgebra( expression= Acm+Ccm+Ecm, name="Vcm" )

covAcf	<-mxAlgebra( expression= a_cf %*% t(a_cf), name="Acf" )
covCcf	<-mxAlgebra( expression= c_cf %*% t(c_cf), name="Ccf" )
covEcf	<-mxAlgebra( expression= e_cf %*% t(e_cf), name="Ecf" )
covPcf	<-mxAlgebra( expression= Acf+Ccf+Ecf, name="Vcf" )

covAcmf	<-mxAlgebra( expression= a_cm %*% t(a_cf), name="Acmf" )
covCcmf	<-mxAlgebra( expression= c_cm %*% t(c_cf), name="Ccmf" )

covAcfm	<-mxAlgebra( expression= a_cf %*% t(a_cm), name="Acfm" )
covCcfm	<-mxAlgebra( expression= c_cf %*% t(c_cm), name="Ccfm" )

# Constraint on total variance of Ordinal variable (A+C+E=1)
varLm		<- mxConstraint( expression=Vcm[1,1]==1, name="Lm" )
varLf		<- mxConstraint( expression=Vcf[1,1]==1, name="Lf" )

# Var-Cov of measured vars in terms of latent factors and AC, Cc, and Ec
FcovMZM	<-mxAlgebra( expression= (FactLTwm  %&% rbind ( cbind(Vcm, Acm+Ccm), cbind(Acm+Ccm, Vcm) )) , name="expFCovMZM" )#This traces the path from vars to factors and back to vars
FcovDZM	<-mxAlgebra( expression= (FactLTwm  %&% rbind ( cbind(Vcm, .5%x%Acm+Ccm), cbind(.5%x%Acm+Ccm, Vcm) )) , name="expFCovDZM" )
FcovMZF	<-mxAlgebra( expression= (FactLTwf  %&% rbind ( cbind(Vcf, Acf+Ccf), cbind(Acf+Ccf, Vcf) )) , name="expFCovMZF" )#This traces the path from vars to factors and back to vars
FcovDZF	<-mxAlgebra( expression= (FactLTwf  %&% rbind ( cbind(Vcf, .5%x%Acf+Ccf), cbind(.5%x%Acf+Ccf, Vcf) )) , name="expFCovDZF" )
FcovDZO	<-mxAlgebra( expression= (FactLTwo  %&% rbind ( cbind(Vcm, .5%x%Acmf+Ccmf), cbind(.5%x%Acfm+Ccfm, Vcf) )) , name="expFCovDZO" )

SpcovMZM	<-mxAlgebra( expression= rbind (cbind(Vsm, Asm+Csm), cbind(Asm+Csm, Vsm) ) , name="expSpCovMZM" )
SpcovDZM	<-mxAlgebra( expression= rbind (cbind(Vsm, .5%x%Asm+Csm), cbind(.5%x%Asm+Csm, Vsm) ) , name="expSpCovDZM" )
SpcovMZF	<-mxAlgebra( expression= rbind (cbind(Vsf, Asf+Csf), cbind(Asf+Csf, Vsf) ) , name="expSpCovMZF" )
SpcovDZF	<-mxAlgebra( expression= rbind (cbind(Vsf, .5%x%Asf+Csf), cbind(.5%x%Asf+Csf, Vsf) ) , name="expSpCovDZF" )
SpcovDZO	<-mxAlgebra( expression= rbind (cbind(Vsm, .5%x%Asmf+Csmf), cbind(.5%x%Asfm+Csfm, Vsf) ) , name="expSpCovDZO" )

TOTcovMZM	<-mxAlgebra( expression= expFCovMZM + expSpCovMZM , name="TOTexpCovMZM" )
TOTcovDZM	<-mxAlgebra( expression= expFCovDZM + expSpCovDZM , name="TOTexpCovDZM" )
TOTcovMZF	<-mxAlgebra( expression= expFCovMZF + expSpCovMZF , name="TOTexpCovMZF" )
TOTcovDZF	<-mxAlgebra( expression= expFCovDZF + expSpCovDZF , name="TOTexpCovDZF" )
TOTcovDZO	<-mxAlgebra( expression= expFCovDZO + expSpCovDZO , name="TOTexpCovDZO" )

# *******************************************************************************************************
# Calculator

# Standardize the Total var/covariances matrices of the observed variables
Id8		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I8" )
Rfactmzm	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovMZM)) %&% TOTexpCovMZM, name="FactcorMZM" )
Rfactdzm	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovDZM)) %&% TOTexpCovDZM, name="FactcorDZM" )
Rfactmzf	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovMZF)) %&% TOTexpCovMZF, name="FactcorMZF" )
Rfactdzf	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovDZF)) %&% TOTexpCovDZF, name="FactcorDZF" )
Rfactdzo	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovDZO)) %&% TOTexpCovDZO, name="FactcorDZO" )

# Standardize the Common Effects
stcovAcm	<-mxAlgebra( expression= Acm/Vcm, name="stAcm" )
stcovCcm	<-mxAlgebra( expression= Ccm/Vcm, name="stCcm" )
stcovEcm	<-mxAlgebra( expression= Ecm/Vcm, name="stEcm" )

stcovAcf	<-mxAlgebra( expression= Acf/Vcf, name="stAcf" )
stcovCcf	<-mxAlgebra( expression= Ccf/Vcf, name="stCcf" )
stcovEcf	<-mxAlgebra( expression= Ecf/Vcf, name="stEcf" )

# Standardize the Specific Effects
stcovAsm	<-mxAlgebra( expression= Asm/( (FactLm %&% Vcm) +Vsm), name="stAsm" )
stcovCsm	<-mxAlgebra( expression= Csm/( (FactLm %&% Vcm) +Vsm), name="stCsm" )
stcovEsm	<-mxAlgebra( expression= Esm/( (FactLm %&% Vcm) +Vsm), name="stEsm" )

stcovAsf	<-mxAlgebra( expression= Asf/( (FactLf %&% Vcf) +Vsf), name="stAsf" )
stcovCsf	<-mxAlgebra( expression= Csf/( (FactLf %&% Vcf) +Vsf), name="stCsf" )
stcovEsf	<-mxAlgebra( expression= Esf/( (FactLf %&% Vcf) +Vsf), name="stEsf" )

# Standardized Effects of Individual variables from the factors (Variance components) above
#stAvar	<-mxAlgebra( expression= (FactL %&% FAc)/( (FactL %&% FVc) +Vs), name="stAvariables" )
#stCvar	<-mxAlgebra( expression= (FactL %&% FCc)/( (FactL %&% FVc) +Vs), name="stCvariables" )
#stEvar	<-mxAlgebra( expression= (FactL %&% FEc)/( (FactL %&% FVc) +Vs), name="stEvariables" )

# Standardized Factor Loadings
StFLm		<-mxAlgebra( expression= sqrt(diag2vec( FactLm %&% Vcm / TOTexpCovMZM[1:4,1:4])) , name="StandFactm" )
StFLf		<-mxAlgebra( expression= sqrt(diag2vec( FactLf %&% Vcf / TOTexpCovMZF[1:4,1:4])) , name="StandFactf" )

# *******************************************************************************************************

# Data objects for Multiple Groups
dataMZM	<- mxData( observed=mzmData, type="raw" )
dataDZM	<- mxData( observed=dzmData, type="raw" )
dataMZF	<- mxData( observed=mzfData, type="raw" )
dataDZF	<- mxData( observed=dzfData, type="raw" )
dataDZO	<- mxData( observed=dzoData, type="raw" )

# Objective objects for Multiple Groups
objMZM	<- mxExpectationNormal( covariance="TOTexpCovMZM", means="expMeanm", dimnames=selVars, thresholds="expThresm", threshnames=c("SO1","SO2"))
objMZF	<- mxExpectationNormal( covariance="TOTexpCovMZF", means="expMeanf", dimnames=selVars, thresholds="expThresf", threshnames=c("SO1","SO2"))
objDZM	<- mxExpectationNormal( covariance="TOTexpCovDZM", means="expMeanm", dimnames=selVars, thresholds="expThresm", threshnames=c("SO1","SO2"))
objDZF	<- mxExpectationNormal( covariance="TOTexpCovDZF", means="expMeanf", dimnames=selVars, thresholds="expThresf", threshnames=c("SO1","SO2"))
objDZO	<- mxExpectationNormal( covariance="TOTexpCovDZO", means="expMeano", dimnames=selVars, thresholds="expThreso", threshnames=c("SO1","SO2"))

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1m	<-list(Meansm,Loadm,LoadTwm,PathsAsm,PathsCsm,PathsEsm,covAsm,covCsm,covEsm,covPsm,Id2,Id8,Ze43)
pars1f	<-list(Meansf,Loadf,LoadTwf,PathsAsf,PathsCsf,PathsEsf,covAsf,covCsf,covEsf,covPsf,Id2,Id8,Ze43)

pars2m	<-list(PathsAcm,PathsCcm,PathsEcm,covAcm,covCcm,covEcm,covPcm,stcovAcm,stcovCcm,stcovEcm, stcovAsm, stcovCsm, stcovEsm)
pars2f	<-list(PathsAcf,PathsCcf,PathsEcf,covAcf,covCcf,covEcf,covPcf,stcovAcf,stcovCcf,stcovEcf, stcovAsf, stcovCsf, stcovEsf)

pars3m	<-list(obsage1, obsage2, betaAm, Thrm, inc, Thresm)
pars3f	<-list(obsage1, obsage2, betaAf, Thrf, inc, Thresf)

modelMZM	<-mxModel(pars1m, pars2m, pars3m, FcovMZM, SpcovMZM, TOTcovMZM, dataMZM, objMZM, Rfactmzm, fitFunction, varLm, StFLm, name="MZM" )
modelDZM	<-mxModel(pars1m, pars2m, pars3m, FcovDZM, SpcovDZM, TOTcovDZM, dataDZM, objDZM, Rfactdzm, fitFunction, name="DZM" )
modelMZF	<-mxModel(pars1f, pars2f, pars3f, FcovMZF, SpcovMZF, TOTcovMZF, dataMZF, objMZF, Rfactmzf, fitFunction, varLf, StFLf, name="MZF" )
modelDZF	<-mxModel(pars1f, pars2f, pars3f, FcovDZF, SpcovDZF, TOTcovDZF, dataDZF, objDZF, Rfactdzf, fitFunction, name="DZF" )
modelDZO	<-mxModel(pars1m, pars2m, pars3m, pars1f, pars2f, pars3f,Meanso,LoadTwo, Threso,covAcmf,covCcmf,covAcfm,covCcfm,covAsmf,covCsmf,covAsfm,covCsfm, FcovDZO, SpcovDZO, TOTcovDZO, dataDZO, objDZO, Rfactdzo, fitFunction, name="DZO" )

minus2ll	<-mxAlgebra( expression=MZM.objective + DZM.objective + MZF.objective + DZF.objective + DZO.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )
cistFL	m	<-mxCI (c ('MZM.StandFactm'))
cistFcm	<-mxCI (c ('MZM.stAcm[1,1]','MZM.stAcm[2,1]','MZM.stAcm[3,1]','MZM.stAcm[2,2]','MZM.stAcm[3,2]','MZM.stAcm[3,3]',
				'MZM.stCcm[1,1]','MZM.stCcm[2,1]','MZM.stCcm[3,1]','MZM.stCcm[2,2]','MZM.stCcm[3,2]','MZM.stCcm[3,3]',
				'MZM.stEcm[1,1]','MZM.stEcm[2,1]','MZM.stEcm[3,1]','MZM.stEcm[2,2]','MZM.stEcm[3,2]','MZM.stEcm[3,3]') ) 	# standardized var comp from Common feactors	
cistVsm	<-mxCI (c ('MZM.stAsm[1,1]','MZM.stAsm[2,2]','MZM.stAsm[3,3]','MZM.stAsm[4,4]',
				'MZM.stCsm[1,1]','MZM.stCsm[2,2]','MZM.stCsm[3,3]','MZM.stCsm[4,4]',
				'MZM.stEsm[1,1]','MZM.stEsm[2,2]','MZM.stEsm[3,3]','MZM.stEsm[4,4]') ) 	# standardized var comp from specific Factors

cistFLf	<-mxCI (c ('MZF.StandFactf'))
cistFcf	<-mxCI (c ('MZF.stAcf[1,1]','MZF.stAcf[2,1]','MZF.stAcf[3,1]','MZF.stAcf[2,2]','MZF.stAcf[3,2]','MZF.stAcf[3,3]',
				'MZF.stCcf[1,1]','MZF.stCcf[2,1]','MZF.stCcf[3,1]','MZF.stCcf[2,2]','MZF.stCcf[3,2]','MZF.stCcf[3,3]',
				'MZF.stEcf[1,1]','MZF.stEcf[2,1]','MZF.stEcf[3,1]','MZF.stEcf[2,2]','MZF.stEcf[3,2]','MZF.stEcf[3,3]') ) 	# standardized var comp from Common feactors	
cistVsf	<-mxCI (c ('MZF.stAsf[1,1]','MZF.stAsf[2,2]','MZF.stAsf[3,3]','MZF.stAsf[4,4]',
				'MZF.stCsf[1,1]','MZF.stCsf[2,2]','MZF.stCsf[3,3]','MZF.stCsf[4,4]',
				'MZF.stEsf[1,1]','MZF.stEsf[2,2]','MZF.stEsf[3,3]','MZF.stEsf[4,4]') ) 	# standardized var comp from specific Factors
HetACEModel	<-mxModel("Hetace", pars1m, pars2m, pars1f, pars2f, modelMZM, modelDZM, modelMZF, modelDZF, modelDZO, minus2ll, obj, cistFLm, cistFcm, cistVsm, cistFLf, cistFcf, cistVsf) 

# --------------------------------------------------------------------------------------------------------------------------------
# 4 RUN HetACE Factor Model: Cholesky (by Zygosity)

HetACEFit		<-mxTryHardOrdinal(HetACEModel, intervals=T)
#ACEFit		<-mxTryHard(ACEModel, intervals=F, bestInitsOutput=TRUE, showInits=TRUE)
(HetACESumm		<-summary(HetACEFit, verbose=F))

# Get some output

mxEval(MZ.Vs, ACEFit)

mxEval(MZ.FAc, ACEFit)
mxEval(MZ.FCc, ACEFit)
mxEval(MZ.FEc, ACEFit)
mxEval(MZ.FVc, ACEFit)
mxEval(MZ.PhC, ACEFit)

mxEval(MZM.stAcm, HetACEFit)
mxEval(MZM.stCcm, HetACEFit)
mxEval(MZM.stEcm, HetACEFit)

mxEval(MZF.stAcf, HetACEFit)
mxEval(MZF.stCcf, HetACEFit)
mxEval(MZF.stEcf, HetACEFit)

mxEval(MZ.stAs, ACEFit)
mxEval(MZ.stCs, ACEFit)
mxEval(MZ.stEs, ACEFit)

mxEval(MZ.StandFact, ACEFit)

mxEval(MZ.stAvariables, ACEFit)
mxEval(MZ.stCvariables, ACEFit)
mxEval(MZ.stEvariables, ACEFit)

#------------------------------------------------------------------------
# Homogeneity model
# Drop Ccommon pathsfrom previous model 
# I.e. Cholesky model without causal paths but with C-spec dropped
# -----------------------------------------------------------------------
subModel1	<- mxModel(HetACEFit, name="sub1")
subModel1	<- omxSetParameters(subModel1, labels=c('cs2m'), free=FALSE, values=c(.5), newlabels=c('cs2f'))
subModel1	<- omxSetParameters(subModel1, labels=c('c11m','c21m','c31m','c22m','c32m','c33m'), free=FALSE, values=c(.3), newlabels=c('c11f','c21f','c31f','c22f','c32f','c33f'))
subModel1	<- omxSetParameters(subModel1, labels=c('l3m'), free=FALSE, values=c(.5), newlabels=c('l4f'))
subModel1	<- omxAssignFirstParameters(subModel1)
subFit1	<- mxRun(subModel1, intervals=F)
(subSum1	<- summary(subFit1))

mxCompare(HetACEFit, subFit1)

mxEval(MZ.stAc, subFit1)
mxEval(MZ.stCc, subFit1)
mxEval(MZ.stEc, subFit1)



####

#****************************************************************************************************************************
# __(Vai-sex diff)_____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL for SO on PD with sex differences
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Phenotypic Factors: F1>F2>F3 & F1>F3;
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
# this because applying a constraint on the factor variances of 1 is problematic especially when we model the causal paths.
# To identify the model we constrain Asp, Csp and Esp variance components loading on variables 1 and 10 to zero
#_____________________________________________________________________________________________________________________________

nv		<- 4				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 3				# number of Latent Factors for Mediation Model per twin
nfact2		<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nvo 		<- 1     			# number of ordinal variables per twin
nvc 		<- nv-nvo  			# number of continuous variables per twin
poso 		<- nvo 			# position where ordinal variables start
nth		<- 4				# number of max thresholds
ninc 		<- nth-1 			# number of max increments
ncovariates 	<- 2 				# number of covariates
nlower		<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups		<- c("mzm", "dzm", "mzf", "dzf", "dzo")
Vars		<- c('PRSSO','SO','Dep','Anx')
selVars	<- c('PRSSO1','SO1','Dep1','Anx1',
		     'PRSSO2','SO2','Dep2','Anx2')
useVars	<- c('PRSSO1','SO1','Dep1','Anx1',
		     'PRSSO2','SO2','Dep2','Anx2','age1','sex1','age2','sex2')

mzmData		<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1) , useVars)
mzfData		<- subset(TWINdata2, zyg1%in%c(3)|zyg2%in%c(3) , useVars)
dzmData		<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2) , useVars)
dzfData		<- subset(TWINdata2, zyg1%in%c(4)|zyg2%in%c(4) , useVars)
dzoData		<- subset(TWINdata2, zyg1%in%c(5)|zyg2%in%c(5) , useVars)

psych::describe(mzmData)
psych::describe(mzfData)
psych::describe(dzmData)
psych::describe(dzfData)
psych::describe(dzoData)

mzmData$SO1[is.na(mzmData$age1)] <- NA
mzmData$SO2[is.na(mzmData$age2)] <- NA
mzfData$SO1[is.na(mzfData$age1)] <- NA
mzfData$SO2[is.na(mzfData$age2)] <- NA
dzmData$SO1[is.na(dzmData$age1)] <- NA
dzmData$SO2[is.na(dzmData$age2)] <- NA
dzfData$SO1[is.na(dzfData$age1)] <- NA
dzfData$SO2[is.na(dzfData$age2)] <- NA
dzoData$SO1[is.na(dzoData$age1)] <- NA
dzoData$SO2[is.na(dzoData$age2)] <- NA

mzmData$age1[is.na(mzmData$age1)] <- 999
mzmData$age2[is.na(mzmData$age2)] <- 999
mzfData$age1[is.na(mzfData$age1)] <- 999
mzfData$age2[is.na(mzfData$age2)] <- 999
dzmData$age1[is.na(dzmData$age1)] <- 999
dzmData$age2[is.na(dzmData$age2)] <- 999
dzfData$age1[is.na(dzfData$age1)] <- 999
dzfData$age2[is.na(dzfData$age2)] <- 999
dzoData$age1[is.na(dzoData$age1)] <- 999
dzoData$age2[is.na(dzoData$age2)] <- 999

psych::describe(mzmData)
psych::describe(mzfData)
psych::describe(dzmData)
psych::describe(dzfData)
psych::describe(dzoData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabsm	<- paste("mm",1:nv,sep=""))
(mLabsf	<- paste("mf",1:nv,sep=""))
(Stmeanm	<- c(5,0,4.5,4.1))
(Stmeanf	<- c(5,0,4.5,4.1))
(PatM		<- c(TRUE,F,TRUE,TRUE))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEsm	<- c('es1m','es2m','es3m','es3m'))
(LabAsm	<- c('as1m','as2m','as3m','as3m'))
(LabCsm	<- c('cs1m','cs2m','cs3m','cs3m'))

(LabEsf	<- c('es1f','es2f','es3f','es3f'))
(LabAsf	<- c('as1f','as2f','as3f','as3f'))
(LabCsf	<- c('cs1f','cs2f','cs3f','cs3f'))

PatSpe		<- c(F,F,TRUE,TRUE)
PatSpac	<- c(F,F,TRUE,TRUE)
StSpam		<- c(0,0,.2,.2)
StSpcm		<- c(0,0,.05,.05)
StSpem		<- c(0,0,.5,.5)

StSpaf		<- c(0,0,.3,.3)
StSpcf		<- c(0,0,.1,.1)
StSpef		<- c(0,0,.7,.7)

# all 1st loadings fixed to 1
PatFl		<- c(F,F,F,F,			
		     F,F,F,F,
		     F,F,F,T)

StFl		<- c(1,0,0,0,
		     0,1,0,0,
		     0,0,1,.9)

LabFlm		<- c('l1m',NA,NA,NA,
	 	     NA,'l2m',NA,NA,
	 	     NA,NA,'l3m','l4m')

LabFlf		<- c('l1f',NA,NA,NA,
	 	     NA,'l2f',NA,NA,
	 	     NA,NA,'l3f','l4f')

PatPhC		<- c(F,T,T,
		     F,F,T,
		     F,F,F)

StPhCm		<- c(0,.05,.01,
		     0,0,.05,
		     0,0,0)

StPhCf		<- c(0,.05,.01,
		     0,0,.1,
		     0,0,0)

LabPhCm	<- c(NA,'c1on2m','c1on3m',
		     NA,NA,'c2on3m',
		     NA,NA,NA)	 

LabPhCf	<- c(NA,'c1on2f','c1on3f',
		     NA,NA,'c2on3f',
		     NA,NA,NA)	 

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Meansm		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmeanm,Stmeanm), labels=c(mLabsm,mLabsm), name="expMeanm") 
Meansf		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmeanf,Stmeanf), labels=c(mLabsf,mLabsf), name="expMeanf") 
Meanso		<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmeanm,Stmeanf), labels=c(mLabsm,mLabsf), name="expMeano") 

# Threshold and covariates
LabThm		<-c('T_1m','i_11m','i_12m','i_13m')	# THs for var 1 
LabThf		<-c('T_1f','i_11f','i_12f','i_13f')	# THs for var 1 

LabCovAm	<-c('BageThSOm','BageThSOm','BageThSOm','BageThSOm')
LabCovAf	<-c('BageThSOf','BageThSOf','BageThSOf','BageThSOf')

ThPat		<-c(T,T,T,T)
StTHm		<-c(1.4,.5,.05,.2)
StTHf		<-c(.1,1,4,.5)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

# effect of age and sex on ordinal variable
betaAm		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.02, labels=LabCovAm, name="BageTHm" )
betaAf		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.03, labels=LabCovAf, name="BageTHf" )
 
# thresholds
Thrm		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTHm, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabThm, name="Thm")
Thrf		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTHf, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabThf, name="Thf")

inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")

Thresm		<-mxAlgebra( expression= cbind(Low%*%Thm + BageTHm%x%age1, Low%*%Thm + BageTHm%x%age2), name="expThresm")
Thresf		<-mxAlgebra( expression= cbind(Low%*%Thf + BageTHf%x%age1, Low%*%Thf + BageTHf%x%age2), name="expThresf")
Threso		<-mxAlgebra( expression= cbind(Low%*%Thm + BageTHm%x%age1, Low%*%Thf + BageTHf%x%age2), name="expThreso")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Loadm		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFlm, name="FactLm" )
Loadf		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFlf, name="FactLf" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTwm	<-mxAlgebra(I2%x%FactLm, name="FactLTwm")
LoadTwf	<-mxAlgebra(I2%x%FactLf, name="FactLTwf")
Ze43		<-mxMatrix("Zero", nv, nfact, free=F, name="Z43")
LoadTwo	<-mxAlgebra(rbind(cbind(FactLm,Z43), cbind(Z43, FactLf)), name="FactLTwo")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  
PhCausm	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhCm, labels=LabPhCm, name="PhCm" )
PhCausf	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhCf, labels=LabPhCf, name="PhCf" )

# Define the matrix to hold the A and C effects: Specific 
PathsAsm	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpam, labels=LabAsm, name="asm" )
PathsCsm	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpcm, labels=LabCsm, name="csm" )
PathsEsm	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpem, labels=LabEsm, name="esm" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
PathsAsf	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpaf, labels=LabAsf, name="asf" )
PathsCsf	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpcf, labels=LabCsf, name="csf" )
PathsEsf	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpef, labels=LabEsf, name="esf" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components

covAsm		<-mxAlgebra( expression= asm %*% t(asm), name="Asm" )
covCsm		<-mxAlgebra( expression= csm %*% t(csm), name="Csm" )
covEsm		<-mxAlgebra( expression= esm %*% t(esm), name="Esm" )
covAsf		<-mxAlgebra( expression= asf %*% t(asf), name="Asf" )
covCsf		<-mxAlgebra( expression= csf %*% t(csf), name="Csf" )
covEsf		<-mxAlgebra( expression= esf %*% t(esf), name="Esf" )

covPsm		<-mxAlgebra( expression= Asm+Csm+Esm, name="Vsm" )
covPsf		<-mxAlgebra( expression= Asf+Csf+Esf, name="Vsf" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsubm		<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.2,.1,.2), labels=c("ac22m","ac32m","ac33m"), name="acm" ) # Component paths for factors 2 and 3
PathsCcsubm		<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.3,.1,.01), labels=c("cc22m","cc32m","cc33m"), name="ccm" )
PathsEcsubm		<-mxMatrix(type="Diag", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.9,.3), labels=c("ec22m","ec33m"), name="ecm" )
PathsAcsubf		<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.6,.1,.3), labels=c("ac22f","ac32f","ac33f"), name="acf" ) # Component paths for factors 2 and 3
PathsCcsubf		<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.4,.1,.01), labels=c("cc22f","cc32f","cc33f"), name="ccf" )
PathsEcsubf		<-mxMatrix(type="Diag", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.7,.4), labels=c("ec22f","ec33f"), name="ecf" )

PathsP11m	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11m", name="pcm" ) # SD path for factor 3 (the PRS factor)
PathsP11f	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11f", name="pcf" ) 

Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  #Padding
Ze21		<-mxMatrix(type="Zero",	nrow=2, ncol=1, free=F, name="Z21" )  #Padding
Ze12		<-mxMatrix(type="Zero",	nrow=1, ncol=2, free=F, name="Z12" )  #Padding

covAcsubm	<-mxAlgebra( expression= acm %*% t(acm), name="Acsubm" )
covCcsubm	<-mxAlgebra( expression= ccm %*% t(ccm), name="Ccsubm" )
covEcsubm	<-mxAlgebra( expression= ecm %*% t(ecm), name="Ecsubm" )
covAcsubf	<-mxAlgebra( expression= acf %*% t(acf), name="Acsubf" )
covCcsubf	<-mxAlgebra( expression= ccf %*% t(ccf), name="Ccsubf" )
covEcsubf	<-mxAlgebra( expression= ecf %*% t(ecf), name="Ecsubf" )

covPcsubm	<-mxAlgebra( expression= Acsubm+Ccsubm+Ecsubm, name="Vcsubm" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)
covPcsubf	<-mxAlgebra( expression= Acsubf+Ccsubf+Ecsubf, name="Vcsubf" ) 

covPc11m	<-mxAlgebra( expression= pcm %*% t(pcm), name="Pc11m" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components
covPc11f	<-mxAlgebra( expression= pcf %*% t(pcf), name="Pc11f" ) 

covPcm		<-mxAlgebra(cbind(rbind(Pc11m,Z21), rbind(Z12,Vcsubm)), name="Vcm") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPcf		<-mxAlgebra(cbind(rbind(Pc11f,Z21), rbind(Z12,Vcsubf)), name="Vcf") 

covPcMzm	<-mxAlgebra(cbind(rbind(Pc11m,Z21) ,rbind(Z12,Acsubm+Ccsubm)), name="Vcmzm") #I specify the MZ between-twin covariance - excluding E parameters
covPcMzf	<-mxAlgebra(cbind(rbind(Pc11f,Z21) ,rbind(Z12,Acsubf+Ccsubf)), name="Vcmzf") 
covPcDzm	<-mxAlgebra(cbind(rbind(.5%x%Pc11m,Z21), rbind(Z12,.5%x%(acm %*% t(acm))+(ccm %*% t(ccm)))), name="Vcdzm") #I specify the DZ between-twin covariance - specifying half of A and excluding E
covPcDzf	<-mxAlgebra(cbind(rbind(.5%x%Pc11f,Z21), rbind(Z12,.5%x%(acf %*% t(acf))+(ccf %*% t(ccf)))), name="Vcdzf") 
covPcDzmf	<-mxAlgebra(cbind(rbind(.5%x%Pc11m,Z21), rbind(Z12,.5%x%(acm %*% t(acf))+(ccm %*% t(ccf)))), name="Vcdzmf") 
covPcDzfm	<-mxAlgebra(cbind(rbind(.5%x%Pc11f,Z21), rbind(Z12,.5%x%(acf %*% t(acm))+(ccf %*% t(ccm)))), name="Vcdzfm") 

# Generate Covariance of Latent factor model Including Causal Paths between factors
Id3		<-mxMatrix(type="Iden",	nrow=3, ncol=3, free=F, name="I3" )
covFVcm	<-mxAlgebra( expression= solve(I3-PhCm) %&% Vcm, name ="FVcm")
covFVcf	<-mxAlgebra( expression= solve(I3-PhCf) %&% Vcf, name ="FVcf")
covFcMzm	<-mxAlgebra( expression= solve(I3-PhCm) %&% Vcmzm, name ="Fcmzm")
covFcMzf	<-mxAlgebra( expression= solve(I3-PhCf) %&% Vcmzf, name ="Fcmzf")
covFcDzm	<-mxAlgebra( expression= solve(I3-PhCm) %&% Vcdzm, name ="Fcdzm")
covFcDzf	<-mxAlgebra( expression= solve(I3-PhCf) %&% Vcdzf, name ="Fcdzf")
covFcDzmf	<-mxAlgebra( expression= solve(I3-PhCm) %*% Vcdzmf %*% (I3-PhCf), name ="Fcdzmf")
covFcDzfm	<-mxAlgebra( expression= solve(I3-PhCf) %*% Vcdzfm %*% (I3-PhCm), name ="Fcdzfm")

# Constraint on total variance of Ordinal variable (A+C+E=1)
varL1m		<- mxConstraint( expression=FVcm[2,2]==1, name="L1m" )
varL1f		<- mxConstraint( expression=FVcf[2,2]==1, name="L1f" )

FcovMZM	<-mxAlgebra( expression= (FactLTwm  %&% rbind ( cbind(FVcm, Fcmzm), cbind(Fcmzm, FVcm))) , name="expFCovMZM" )
FcovMZF	<-mxAlgebra( expression= (FactLTwf  %&% rbind ( cbind(FVcf, Fcmzf), cbind(Fcmzf, FVcf))) , name="expFCovMZF" )
FcovDZM	<-mxAlgebra( expression= (FactLTwm  %&% rbind ( cbind(FVcm, Fcdzm), cbind(Fcdzm, FVcm))) , name="expFCovDZM" )
FcovDZF	<-mxAlgebra( expression= (FactLTwf  %&% rbind ( cbind(FVcf, Fcdzf), cbind(Fcdzf, FVcf))) , name="expFCovDZF" )
FcovDZO	<-mxAlgebra( expression= (FactLTwo  %&% rbind ( cbind(FVcm, Fcdzmf), cbind(Fcdzfm, FVcf))) , name="expFCovDZO" )

SpcovMZM	<-mxAlgebra( expression= rbind (cbind(Vsm, Asm+Csm), cbind(Asm+Csm, Vsm)) , name="expSpCovMZM" )
SpcovMZF	<-mxAlgebra( expression= rbind (cbind(Vsf, Asf+Csf), cbind(Asf+Csf, Vsf)) , name="expSpCovMZF" )
SpcovDZM	<-mxAlgebra( expression= rbind (cbind(Vsm, .5%x%Asm+Csm), cbind(.5%x%Asm+Csm, Vsm)) , name="expSpCovDZM" )
SpcovDZF	<-mxAlgebra( expression= rbind (cbind(Vsf, .5%x%Asf+Csf), cbind(.5%x%Asf+Csf, Vsf)) , name="expSpCovDZF" )
SpcovDZO	<-mxAlgebra( expression= rbind (cbind(Vsm, .5%x%(asm %*% t(asf))+(csm %*% t(csf))), cbind(.5%x%(asf %*% t(asm))+(csf %*% t(csm)), Vsf)) , name="expSpCovDZO" )

TOTcovMZM	<-mxAlgebra( expression= expFCovMZM + expSpCovMZM , name="TOTexpCovMZM" )
TOTcovMZF	<-mxAlgebra( expression= expFCovMZF + expSpCovMZF , name="TOTexpCovMZF" )
TOTcovDZM	<-mxAlgebra( expression= expFCovDZM + expSpCovDZM , name="TOTexpCovDZM" )
TOTcovDZF	<-mxAlgebra( expression= expFCovDZF + expSpCovDZF , name="TOTexpCovDZF" )
TOTcovDZO	<-mxAlgebra( expression= expFCovDZO + expSpCovDZO , name="TOTexpCovDZO" )

# *******************************************************************************************************
# Calculator

# Standardize the causal effects
Stcp1on2m	<-mxAlgebra( expression= (PhCm[2,1]* sqrt(FVcm[1,1]))/sqrt(FVcm[2,2]) , name="Stand_1on2m" )
Stcp1on3m	<-mxAlgebra( expression= (PhCm[3,1]* sqrt(FVcm[1,1]))/sqrt(FVcm[3,3]) , name="Stand_1on3m" )
Stcp2on3m	<-mxAlgebra( expression= (PhCm[3,2]* sqrt(FVcm[2,2]))/sqrt(FVcm[3,3]) , name="Stand_2on3m" )

Stcp1on2f	<-mxAlgebra( expression= (PhCf[2,1]* sqrt(FVcf[1,1]))/sqrt(FVcf[2,2]) , name="Stand_1on2f" )
Stcp1on3f	<-mxAlgebra( expression= (PhCf[3,1]* sqrt(FVcf[1,1]))/sqrt(FVcf[3,3]) , name="Stand_1on3f" )
Stcp2on3f	<-mxAlgebra( expression= (PhCf[3,2]* sqrt(FVcf[2,2]))/sqrt(FVcf[3,3]) , name="Stand_2on3f" )

# Standardize the Total var/covariances matrices of the observed variables
Id8		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I8" )
Rfactmzm	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovMZM)) %&% TOTexpCovMZM, name="FactcorMZM" )
Rfactmzf	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovMZF)) %&% TOTexpCovMZF, name="FactcorMZF" )
Rfactdzm	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovDZM)) %&% TOTexpCovDZM, name="FactcorDZM" )
Rfactdzf	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovDZF)) %&% TOTexpCovDZF, name="FactcorDZF" )
Rfactdzo	<-mxAlgebra( expression= solve(sqrt(I8*TOTexpCovDZO)) %&% TOTexpCovDZO, name="FactcorDZO" )

# Phenotypic, A, C and E correlations	
RfactAcm	<-mxAlgebra( expression= solve(sqrt(I2*Acsubm)) %&% Acsubm, name="Ram" )
RfactCcm	<-mxAlgebra( expression= solve(sqrt(I2*Ccsubm)) %&% Ccsubm, name="Rcm" )
RfactEcm	<-mxAlgebra( expression= solve(sqrt(I2*Ecsubm)) %&% Ecsubm, name="Rem" )
RfactPm	<-mxAlgebra( expression= solve(sqrt(I3*FVcm)) %&% FVcm, name="Rphm" )

RfactAcf	<-mxAlgebra( expression= solve(sqrt(I2*Acsubf)) %&% Acsubf, name="Raf" )
RfactCcf	<-mxAlgebra( expression= solve(sqrt(I2*Ccsubf)) %&% Ccsubf, name="Rcf" )
RfactEcf	<-mxAlgebra( expression= solve(sqrt(I2*Ecsubf)) %&% Ecsubf, name="Ref" )
RfactPf	<-mxAlgebra( expression= solve(sqrt(I3*FVcf)) %&% FVcf, name="Rphf" )

# Standardize the Common Effects
covFVc22m	<-mxAlgebra( expression= FVcm[2:3,2:3], name ="FVc22m")
covFVc22f	<-mxAlgebra( expression= FVcf[2:3,2:3], name ="FVc22f")

stcovAcm	<-mxAlgebra( expression= Acsubm/FVc22m, name="stAcm" )
stcovCcm	<-mxAlgebra( expression= Ccsubm/FVc22m, name="stCcm" )
stcovEcm	<-mxAlgebra( expression= Ecsubm/FVc22m, name="stEcm" )
stcovAcf	<-mxAlgebra( expression= Acsubf/FVc22f, name="stAcf" )
stcovCcf	<-mxAlgebra( expression= Ccsubf/FVc22f, name="stCcf" )
stcovEcf	<-mxAlgebra( expression= Ecsubf/FVc22f, name="stEcf" )

# Standardised path estimates
StpathAcm	<-mxAlgebra( expression= (sqrt(stAcm)), name="stpacm" )
StpathCcm	<-mxAlgebra( expression= (sqrt(stCcm)), name="stpccm" )
StpathEcm	<-mxAlgebra( expression= (sqrt(stEcm)), name="stpecm" )
StpathAcf	<-mxAlgebra( expression= (sqrt(stAcf)), name="stpacf" )
StpathCcf	<-mxAlgebra( expression= (sqrt(stCcf)), name="stpccf" )
StpathEcf	<-mxAlgebra( expression= (sqrt(stEcf)), name="stpecf" )

# Algebra to compute Rph-A, Rph-C and Rph-E
RphA12m	<-mxAlgebra(expression=sqrt(stAcm[1,1])*Ram[2,1]*sqrt(stAcm[2,2]), name = 'Rpha12m')
RphC12m	<-mxAlgebra(expression=sqrt(stCcm[1,1])*Rcm[2,1]*sqrt(stCcm[2,2]), name = 'Rphc12m')
RphE12m	<-mxAlgebra(expression=sqrt(stEcm[1,1])*Rem[2,1]*sqrt(stEcm[2,2]), name = 'Rphe12m')

RphA12f	<-mxAlgebra(expression=sqrt(stAcf[1,1])*Raf[2,1]*sqrt(stAcf[2,2]), name = 'Rpha12f')
RphC12f	<-mxAlgebra(expression=sqrt(stCcf[1,1])*Rcf[2,1]*sqrt(stCcf[2,2]), name = 'Rphc12f')
RphE12f	<-mxAlgebra(expression=sqrt(stEcf[1,1])*Ref[2,1]*sqrt(stEcf[2,2]), name = 'Rphe12f')

# Standardize the Specific Effects
stcovAsm	<-mxAlgebra( expression= sqrt(Asm/( (FactLm %&% FVcm) +Vsm)), name="stAsm" )
stcovCsm	<-mxAlgebra( expression= sqrt(Csm/( (FactLm %&% FVcm) +Vsm)), name="stCsm" )
stcovEsm	<-mxAlgebra( expression= sqrt(Esm/( (FactLm %&% FVcm) +Vsm)), name="stEsm" )

stcovAsf	<-mxAlgebra( expression= sqrt(Asf/( (FactLf %&% FVcf) +Vsf)), name="stAsf" )
stcovCsf	<-mxAlgebra( expression= sqrt(Csf/( (FactLf %&% FVcf) +Vsf)), name="stCsf" )
stcovEsf	<-mxAlgebra( expression= sqrt(Esf/( (FactLf %&% FVcf) +Vsf)), name="stEsf" )

# Standardized Factor Loadings
StFLm		<-mxAlgebra( expression= sqrt(diag2vec( FactLm %&% FVcm / TOTexpCovMZM[1:4,1:4])) , name="StandFactm" )
StFLf		<-mxAlgebra( expression= sqrt(diag2vec( FactLf %&% FVcf / TOTexpCovMZF[1:4,1:4])) , name="StandFactf" )

# *******************************************************************************************************

# Data objects for Multiple Groups
dataMZM	<- mxData( observed=mzmData, type="raw" )
dataMZF	<- mxData( observed=mzfData, type="raw" )
dataDZM	<- mxData( observed=dzmData, type="raw" )
dataDZF	<- mxData( observed=dzfData, type="raw" )
dataDZO	<- mxData( observed=dzoData, type="raw" )

# Objective objects for Multiple Groups
objMZM		<- mxExpectationNormal( covariance="TOTexpCovMZM", means="expMeanm", dimnames=selVars, thresholds="expThresm", threshnames=c("SO1","SO2"))
objMZF		<- mxExpectationNormal( covariance="TOTexpCovMZF", means="expMeanf", dimnames=selVars, thresholds="expThresf", threshnames=c("SO1","SO2"))
objDZM		<- mxExpectationNormal( covariance="TOTexpCovDZM", means="expMeanm", dimnames=selVars, thresholds="expThresm", threshnames=c("SO1","SO2"))
objDZF		<- mxExpectationNormal( covariance="TOTexpCovDZF", means="expMeanf", dimnames=selVars, thresholds="expThresf", threshnames=c("SO1","SO2"))
objDZO		<- mxExpectationNormal( covariance="TOTexpCovDZO", means="expMeano", dimnames=selVars, thresholds="expThreso", threshnames=c("SO1","SO2"))

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1m		<-list(Meansm,Loadm,LoadTwm,PhCausm,PathsAsm,PathsCsm,PathsEsm,covAsm,covCsm,covEsm,covPsm,Id2,Id3,Id8)
pars1f		<-list(Meansf,Loadf,LoadTwf,PhCausf,PathsAsf,PathsCsf,PathsEsf,covAsf,covCsf,covEsf,covPsf,Id2,Id3,Id8)
pars1o		<-list(Meanso,Ze43, LoadTwo)
pars2m		<-list(PathsAcsubm,PathsCcsubm,PathsEcsubm,PathsP11m,Ze21,Ze12,Ze11,covAcsubm,covCcsubm,covEcsubm,covPcsubm,covPc11m,covPcm,covPcMzm,covPcDzm,covFVcm,covFVc22m,covFcMzm,covFcDzm)
pars2f		<-list(PathsAcsubf,PathsCcsubf,PathsEcsubf,PathsP11f,Ze21,Ze12,Ze11,covAcsubf,covCcsubf,covEcsubf,covPcsubf,covPc11f,covPcf,covPcMzf,covPcDzf,covFVcf,covFVc22f,covFcMzf,covFcDzf)
pars2o		<-list(covPcDzmf,covPcDzfm,covFcDzmf,covFcDzfm)
pars3m		<-list(obsage1, obsage2, betaAm, Thrm, inc, Thresm)
pars3f		<-list(obsage1, obsage2, betaAf, Thrf, inc, Thresf)
parsstm	<-list(stcovAsm, stcovCsm, stcovEsm, stcovAcm, stcovCcm, stcovEcm, RfactAcm, RfactCcm, RfactEcm, RfactPm,RphA12m,RphC12m,RphE12m,StpathAcm,StpathCcm,StpathEcm)
parsstf	<-list(stcovAsf, stcovCsf, stcovEsf, stcovAcf, stcovCcf, stcovEcf, RfactAcf, RfactCcf, RfactEcf, RfactPf,RphA12f,RphC12f,RphE12f,StpathAcf,StpathCcf,StpathEcf)
parsmedm	<-list(Stcp1on2m, Stcp1on3m, Stcp2on3m)
parsmedf	<-list(Stcp1on2f, Stcp1on3f, Stcp2on3f)

modelMZM	<-mxModel(pars1m, pars2m, pars3m, parsmedm, FcovMZM, SpcovMZM, TOTcovMZM, dataMZM, objMZM, Rfactmzm, parsstm, fitFunction, StFLm, varL1m, name="MZM" )
modelMZF	<-mxModel(pars1f, pars2f, pars3f, parsmedf, FcovMZF, SpcovMZF, TOTcovMZF, dataMZF, objMZF, Rfactmzf, parsstf, fitFunction, StFLf, varL1f, name="MZF" )
modelDZM	<-mxModel(pars1m, pars2m, pars3m, FcovDZM, SpcovDZM, TOTcovDZM, dataDZM, objDZM, Rfactdzm, fitFunction, name="DZM" )
modelDZF	<-mxModel(pars1f, pars2f, pars3f, FcovDZF, SpcovDZF, TOTcovDZF, dataDZF, objDZF, Rfactdzf, fitFunction, name="DZF" )
modelDZO	<-mxModel(pars1m, pars2m, pars1f, pars2f, pars1o, pars2o, pars3m, pars3f, Threso, FcovDZO, SpcovDZO, TOTcovDZO, dataDZO, objDZO, Rfactdzo, fitFunction, name="DZO" )

minus2ll	<-mxAlgebra( expression=MZM.objective + MZF.objective + DZM.objective + DZF.objective + DZO.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )

cistFLm	<-mxCI (c ('MZM.StandFactm','MZM.Stand_1on2m','MZM.Stand_1on3m','MZM.Stand_2on3m','MZM.PhCm'))
cistFLf	<-mxCI (c ('MZF.StandFactf','MZF.Stand_1on2f','MZF.Stand_1on3f','MZF.Stand_2on3f','MZF.PhCf'))

cistVsm	<-mxCI (c ('MZM.stAsm[3,3]','MZM.stAsm[4,4]',
				'MZM.stCsm[3,3]','MZM.stCsm[4,4]',
				'MZM.stEsm[3,3]','MZM.stEsm[4,4]') ) 	# standardized var comp from specific Factors
cistVsf	<-mxCI (c ('MZF.stAsf[3,3]','MZF.stAsf[4,4]',
				'MZF.stCsf[3,3]','MZF.stCsf[4,4]',
				'MZF.stEsf[3,3]','MZF.stEsf[4,4]') ) 	# standardized var comp from specific Factors

cistVcm	<-mxCI (c ('MZM.stAcm[1,1]','MZM.stAcm[2,1]','MZM.stAcm[2,2]',
				'MZM.stCcm[1,1]','MZM.stCcm[2,1]','MZM.stCcm[2,2]',
				'MZM.stEcm[1,1]','MZM.stEcm[2,2]') ) 	# standardized var comp for ACE on latent Factors
cistVcf	<-mxCI (c ('MZF.stAcf[1,1]','MZF.stAcf[2,1]','MZF.stAcf[2,2]',
				'MZF.stCcf[1,1]','MZF.stCcf[2,1]','MZF.stCcf[2,2]',
				'MZF.stEcf[1,1]','MZF.stEcf[2,2]') ) 	# standardized var comp for ACE on latent Factors

cistRcm	<-mxCI (c ('MZM.Rpha12m','MZM.Rphc12m','MZM.Rphe12m','MZM.Ram','MZM.Rcm','MZM.Rem','MZM.stpacm','MZM.stpccm','MZM.stpecm') ) 	
cistRcf	<-mxCI (c ('MZF.Rpha12f','MZF.Rphc12f','MZF.Rphe12f','MZF.Raf','MZF.Rcf','MZF.Ref','MZF.stpacf','MZF.stpccf','MZF.stpecf') ) 	

HetACEMs1Model	<-mxModel("HetaceMs1", pars1m, pars2m, pars1f, pars2f, pars1o, pars2o, modelMZM, modelMZF, modelDZM, modelDZF, modelDZO, minus2ll, obj, cistFLm, cistVsm, cistVcm, cistRcm, cistFLf, cistVsf, cistVcf, cistRcf) 

# --------------------------------------------------------------------------------------------------------------------------------
# 5ai RUN HetACEMs Factor Model with phenotypic causal mediation paths by Zygosity

HetACEMs1Fit	<-mxTryHardOrdinal(HetACEMs1Model, intervals=F)
(HetACEMs1Summ	<-summary(HetACEMs1Fit, verbose=F))

# Generate confidence intervals
HetACEMs1CIModel	<-mxModel(HetACEMs1Model)
HetACEMs1CIFit	<-mxRun(HetACEMs1CIModel, intervals=TRUE)
(HetACEMs1CISumm	<-summary(HetACEMs1CIFit, verbose=F))

mxEval(MZM.FactcorMZM, HetACEMs1Fit)
mxEval(DZM.FactcorDZM, HetACEMs1Fit)

mxEval(MZF.FactcorMZF, HetACEMs1Fit)
mxEval(DZF.FactcorDZF, HetACEMs1Fit)

mxEval(MZM.FVcm, HetACEMs1Fit)
mxEval(MZF.FVcf, HetACEMs1Fit)

mxEval(MZM.stAcm, HetACEMs1Fit)
mxEval(MZM.stCcm, HetACEMs1Fit)
mxEval(MZM.stEcm, HetACEMs1Fit)

mxEval(MZF.stAcf, HetACEMs1Fit)
mxEval(MZF.stCcf, HetACEMs1Fit)
mxEval(MZF.stEcf, HetACEMs1Fit)

mxEval(MZM.stpacm, HetACEMs1Fit)
mxEval(MZM.stpccm, HetACEMs1Fit)
mxEval(MZM.stpecm, HetACEMs1Fit)

mxEval(MZF.stpacf, HetACEMs1Fit)
mxEval(MZF.stpccf, HetACEMs1Fit)
mxEval(MZF.stpecf, HetACEMs1Fit)

mxEval(MZM.Ram, HetACEMs1Fit)
mxEval(MZM.Rcm, HetACEMs1Fit)
mxEval(MZM.Rem, HetACEMs1Fit)
mxEval(MZM.Rphm, HetACEMs1Fit)

mxEval(MZF.Raf, HetACEMs1Fit)
mxEval(MZF.Rcf, HetACEMs1Fit)
mxEval(MZF.Ref, HetACEMs1Fit)
mxEval(MZF.Rphf, HetACEMs1Fit)

mxEval(MZM.Rpham, HetACEMs1Fit) 
mxEval(MZM.Rphcm, HetACEMs1Fit) 
mxEval(MZM.Rphem, HetACEMs1Fit) 

mxEval(MZF.Rphaf, HetACEMs1Fit) 
mxEval(MZF.Rphcf, HetACEMs1Fit) 
mxEval(MZF.Rphef, HetACEMs1Fit) 

mxEval(MZM.Acsubm, HetACEMs1Fit)
mxEval(MZM.Ccsubm, HetACEMs1Fit)

mxEval(MZM.Ecsubm, HetACEMs1Fit)
mxEval(MZM.Vcsubm, HetACEMs1Fit)
mxEval(MZM.FVc22m, HetACEMs1Fit)

mxEval(MZF.Acsubf, HetACEMs1Fit)
mxEval(MZF.Ccsubf, HetACEMs1Fit)

mxEval(MZF.Ecsubf, HetACEMs1Fit)
mxEval(MZF.Vcsubf, HetACEMs1Fit)
mxEval(MZF.FVc2f, HetACEMs1Fit)

mxEval(MZM.StandFactm, HetACEMs1Fit)
mxEval(MZM.Stand_1on2m[1,1], HetACEMs1Fit)
mxEval(MZM.Stand_1on3m[1,1], HetACEMs1Fit)
mxEval(MZM.Stand_2on3m[1,1], HetACEMs1Fit)

mxEval(MZF.StandFactf, HetACEMs1Fit)
mxEval(MZF.Stand_1on2f[1,1], HetACEMs1Fit)
mxEval(MZF.Stand_1on3f[1,1], HetACEMs1Fit)
mxEval(MZF.Stand_2on3f[1,1], HetACEMs1Fit)

#------------------------------------------------------------------------
# Submodel 5ai: Constrain causal path coefficients to be equal across sex from previous model 
#------------------------------------------------------------------------
# Drop correlation male c paths from the model
# -----------------------------------------------------------------------
HomACECausModel	<- mxModel(HetACEMs1Fit, name="HomACECaus")
HomACECausModel	<- omxSetParameters(HomACECausModel, labels=c('c2on3f'), free=T, values=.05, newlabels=c('c2on3m'))
HomACECausModel	<- omxAssignFirstParameters(HomACECausModel)
HomACECausFit		<- mxRun(HomACECausModel, intervals=F)
(HomACECausSum	<- summary(HomACECausFit))

mxCompare(HetACEMs1Fit,HomACECausFit)

#------------------------------------------------------------------------
# Submodel 5aiim: Drop correlation path male c paths from previous model 
#------------------------------------------------------------------------
# Drop correlation male c paths from the model
# -----------------------------------------------------------------------
HetAEMs1mModel	<- mxModel(HetACEMs1Fit, name="HetAEMs1m")
HetAEMs1mModel	<- omxSetParameters(HetAEMs1mModel, labels=c('ac33m','cc32m','cc33m'), free=FALSE, values=0)
HetAEMs1mModel	<- omxSetParameters(HetAEMs1mModel, labels=c('cs3m'), free=FALSE, values=0)
HetAEMs1mFit	<- mxRun(HetAEMs1mModel, intervals=F)
(HetAEMs1mSum	<- summary(HetAEMs1mFit))

mxCompare(HetACEMs1Fit, HetAEMs1mFit)

#------------------------------------------------------------------------
# Submodel 5aiif: Drop correlation path female c paths from previous model 
#------------------------------------------------------------------------
# Drop correlation female c paths from the model
# -----------------------------------------------------------------------
HetAEMs1fModel	<- mxModel(HetACEMs1Fit, name="HetAEMs1f")
HetAEMs1fModel	<- omxSetParameters(HetAEMs1fModel, labels=c('cc22f','cc32f','cc33f'), free=FALSE, values=0)
HetAEMs1fModel	<- omxSetParameters(HetAEMs1fModel, labels=c('cs3f'), free=FALSE, values=0)
HetAEMs1fFit	<- mxRun(HetAEMs1fModel, intervals=F)
(HetAEMs1fSum	<- summary(HetAEMs1fFit))

mxCompare(HetACEMs1Fit, HetAEMs1fFit)

# 
#****************************************************************************************************************************
# __(Vaii-Sex diff)_____________________________________________________________________________________________________________________
# Mendelian Randomisation Direction of Causation (MRDoC) MODEL for PD on SO with sex differences
# We specify Specific effects on the latent factors(Acsp, Ccsp and Ecsp) and add causal paths:
# Causal paths specified between Phenotypic Factors: F1>F2>F3 & F1>F3;
# Asp, Csp and Esp in the bottom with constraints to Identify the model on top
# We are estimating the variances of the factors by scaling them to the 1st indicator variable (by fixing the loading to 1), 
#_____________________________________________________________________________________________________________________________

nv		<- 5				# number of variables for a twin = 1 in Univariate
ntv		<- 2*nv			# number of variables for a pair = 2* 1 for Univariate
nfact		<- 3				# number of Latent Factors for Mediation Model per twin
nfact2	<- 2*nfact			# number of Latent Factors for Mediation Model per twin
nvo 		<- 1     			# number of ordinal variables per twin
nvc 		<- nv-nvo  			# number of continuous variables per twin
poso 		<- nvo 			# position where ordinal variables start
nth		<- 4				# number of max thresholds
ninc 		<- nth-1 			# number of max increments
ncovariates <- 2 				# number of covariates
nlower	<- nv*(nv+1)/2 		# number of free elements in a lower matrix nv*nv
ncor		<- (nv*(nv+1)/2)-nv	# number of free elements in a correlation matrix nv*nv
Groups	<- c("mzm", "dzm", "mzf", "dzf", "dzo")
Vars		<- c('PRSAnx','PRSDep','Dep','Anx','SO')
selVars	<- c('PRSAnx1','PRSDep1','Dep1','Anx1','SO1',
		     'PRSAnx2','PRSDep2','Dep2','Anx2','SO2')
useVars	<- c('PRSAnx1','PRSDep1','Dep1','Anx1','SO1',
		     'PRSAnx2','PRSDep2','Dep2','Anx2','SO2','age1','sex1','age2','sex2')

mzmData		<- subset(TWINdata2, zyg1%in%c(1)|zyg2%in%c(1) , useVars)
mzfData		<- subset(TWINdata2, zyg1%in%c(3)|zyg2%in%c(3) , useVars)
dzmData		<- subset(TWINdata2, zyg1%in%c(2)|zyg2%in%c(2) , useVars)
dzfData		<- subset(TWINdata2, zyg1%in%c(4)|zyg2%in%c(4) , useVars)
dzoData		<- subset(TWINdata2, zyg1%in%c(5)|zyg2%in%c(5) , useVars)

psych::describe(mzmData)
psych::describe(mzfData)
psych::describe(dzmData)
psych::describe(dzfData)
psych::describe(dzoData)

mzmData$SO1[is.na(mzmData$age1)] <- NA
mzmData$SO2[is.na(mzmData$age2)] <- NA
mzfData$SO1[is.na(mzfData$age1)] <- NA
mzfData$SO2[is.na(mzfData$age2)] <- NA
dzmData$SO1[is.na(dzmData$age1)] <- NA
dzmData$SO2[is.na(dzmData$age2)] <- NA
dzfData$SO1[is.na(dzfData$age1)] <- NA
dzfData$SO2[is.na(dzfData$age2)] <- NA
dzoData$SO1[is.na(dzoData$age1)] <- NA
dzoData$SO2[is.na(dzoData$age2)] <- NA

mzmData$age1[is.na(mzmData$age1)] <- 999
mzmData$age2[is.na(mzmData$age2)] <- 999
mzfData$age1[is.na(mzfData$age1)] <- 999
mzfData$age2[is.na(mzfData$age2)] <- 999
dzmData$age1[is.na(dzmData$age1)] <- 999
dzmData$age2[is.na(dzmData$age2)] <- 999
dzfData$age1[is.na(dzfData$age1)] <- 999
dzfData$age2[is.na(dzfData$age2)] <- 999
dzoData$age1[is.na(dzoData$age1)] <- 999
dzoData$age2[is.na(dzoData$age2)] <- 999

psych::describe(mzmData)
psych::describe(mzfData)
psych::describe(dzmData)
psych::describe(dzfData)
psych::describe(dzoData)

# CREATE LABELS & START VALUES as objects(to ease specification in the body of the model)
(mLabsm	<- paste("mm",1:nv,sep=""))
(mLabsf	<- paste("mf",1:nv,sep=""))
(Stmeanm	<- c(5,5,4.5,4.1,0))
(Stmeanf	<- c(5,5,4.5,4.1,0))
(PatM		<- c(TRUE,T,TRUE,TRUE,F))

# Create Labels for Diagonal Matrices
# To identify this model we equate the sp effects of the 2 indicators per factor to be equal)
(LabEsm	<- c('es1m','es1m','es3m','es3m','es5m'))
(LabAsm	<- c('as1m','as1m','as3m','as3m','as5m'))
(LabCsm	<- c('cs1m','cs1m','cs3m','cs3m','cs5m'))

(LabEsf	<- c('es1f','es1f','es3f','es3f','es5f'))
(LabAsf	<- c('as1f','as1f','as3f','as3f','as5f'))
(LabCsf	<- c('cs1f','cs1f','cs3f','cs3f','cs4f'))

PatSpe	<- c(T,T,T,T,F)
PatSpac	<- c(F,F,T,T,F)
StSpam	<- c(0,0,.5,.5,0)
StSpcm	<- c(0,0,.5,.5,0)
StSpem	<- c(.5,.5,.5,.5,0)

StSpaf	<- c(0,0,.5,.5,0)
StSpcf	<- c(0,0,.5,.5,0)
StSpef	<- c(.5,.5,.5,.5,0)

# all 1st loadings fixed to 1
PatFl		<- c(F,T,F,F,F,			
		     F,F,F,T,F,
		     F,F,F,F,F)

StFl		<- c(1,.5,0,0,0,
		     0,0,1,.5,0,
		     0,0,0,0,1)

LabFlm	<- c('l1m','l2m',NA,NA,NA,
	 	     NA,NA,'l3m','l4m',NA,
	 	     NA,NA,NA,NA,'l5m')

LabFlf	<- c('l1f','l2f',NA,NA,NA,
	 	     NA,NA,'l3f','l4f',NA,
	 	     NA,NA,NA,NA,'l5f')

PatPhC	<- c(F,T,T,
		     F,F,T,
		     F,F,F)

StPhCm	<- c(0,.05,.01,
		     0,0,.05,
		     0,0,0)

StPhCf	<- c(0,.05,.01,
		     0,0,.1,
		     0,0,0)

LabPhCm	<- c(NA,'c1on2m','c1on3m',
		     NA,NA,'c2on3m',
		     NA,NA,NA)	 

LabPhCf	<- c(NA,'c1on2f','c1on3f',
		     NA,NA,'c2on3f',
		     NA,NA,NA)	 

#______________________________________________________________________________________________________
# Define matrices to hold the Means, SD, correlations
# Use Algebra to generate expected var/cov matrices and Means
# Specify: data objects, Fitfunction, the Model, 
# Run the Model 
#______________________________________________________________________________________________________

Meansm	<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmeanm,Stmeanm), labels=c(mLabsm,mLabsm), name="expMeanm") 
Meansf	<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmeanf,Stmeanf), labels=c(mLabsf,mLabsf), name="expMeanf") 
Meanso	<-mxMatrix("Full", 1, ntv, free=c(PatM,PatM), values=c(Stmeanm,Stmeanf), labels=c(mLabsm,mLabsf), name="expMeano") 

# Threshold and covariates
LabThm	<-c('T_1m','i_11m','i_12m','i_13m')	# THs for var 1 
LabThf	<-c('T_1f','i_11f','i_12f','i_13f')	# THs for var 1 

LabCovAm	<-c('BageThSOm','BageThSOm','BageThSOm','BageThSOm')
LabCovAf	<-c('BageThSOf','BageThSOf','BageThSOf','BageThSOf')

ThPat		<-c(T,T,T,T)
StTHm		<-c(1.4,.5,.05,.2)
StTHf		<-c(.1,1,4,.5)

# Matrices to hold observed covariates (data. = definition variable)
obsage1	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="age1")
obsage2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="age2")

# effect of age and sex on ordinal variable
betaAm	<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.02, labels=LabCovAm, name="BageTHm" )
betaAf	<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=.03, labels=LabCovAf, name="BageTHf" )
 
# thresholds
Thrm		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTHm, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabThm, name="Thm")
Thrf		<-mxMatrix( type="Full", nrow=nth, ncol=nvo, free=TRUE, values=StTHf, lbound=c(-4,-4), ubound=c(4,4),
                  labels=LabThf, name="Thf")

inc		<-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=FALSE, values=1, name="Low")

Thresm	<-mxAlgebra( expression= cbind(Low%*%Thm + BageTHm%x%age1, Low%*%Thm + BageTHm%x%age2), name="expThresm")
Thresf	<-mxAlgebra( expression= cbind(Low%*%Thf + BageTHf%x%age1, Low%*%Thf + BageTHf%x%age2), name="expThresf")
Threso	<-mxAlgebra( expression= cbind(Low%*%Thm + BageTHm%x%age1, Low%*%Thf + BageTHf%x%age2), name="expThreso")

# Define matrices to specify the loadings of the dependent variables on the latent factors
Loadm		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFlm, name="FactLm" )
Loadf		<-mxMatrix(type="Full",	nrow=nv, ncol=nfact, free=PatFl, values=StFl, labels=LabFlf, name="FactLf" )
Id2		<-mxMatrix(type="Iden",	nrow=2, ncol=2, free=F, name="I2" )
LoadTwm	<-mxAlgebra(I2%x%FactLm, name="FactLTwm")
LoadTwf	<-mxAlgebra(I2%x%FactLf, name="FactLTwf")
Ze53		<-mxMatrix("Zero", nv, nfact, free=F, name="Z53")
LoadTwo	<-mxAlgebra(rbind(cbind(FactLm,Z53), cbind(Z53, FactLf)), name="FactLTwo")
 
# Define the matrix to hold the Single headed Arrows (causal paths) between the 3 latent variables  
PhCausm	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhCm, labels=LabPhCm, name="PhCm" )
PhCausf	<-mxMatrix(type="Full",	nrow=nfact, ncol=nfact, free=PatPhC, values=StPhCf, labels=LabPhCf, name="PhCf" )

# Define the matrix to hold the A and C effects: Specific 
PathsAsm	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpam, labels=LabAsm, name="asm" )
PathsCsm	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpcm, labels=LabCsm, name="csm" )
PathsEsm	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpem, labels=LabEsm, name="esm" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components
PathsAsf	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpaf, labels=LabAsf, name="asf" )
PathsCsf	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpac, values=StSpcf, labels=LabCsf, name="csf" )
PathsEsf	<-mxMatrix(type="Diag",	nrow=nv, ncol=nv, free=PatSpe, values=StSpef, labels=LabEsf, name="esf" ) # I make all the residual variances for the PRSs go into E as I am not parsing these into ACE components

covAsm	<-mxAlgebra( expression= asm %*% t(asm), name="Asm" )
covCsm	<-mxAlgebra( expression= csm %*% t(csm), name="Csm" )
covEsm	<-mxAlgebra( expression= esm %*% t(esm), name="Esm" )
covAsf	<-mxAlgebra( expression= asf %*% t(asf), name="Asf" )
covCsf	<-mxAlgebra( expression= csf %*% t(csf), name="Csf" )
covEsf	<-mxAlgebra( expression= esf %*% t(esf), name="Esf" )

covPsm	<-mxAlgebra( expression= Asm+Csm+Esm, name="Vsm" )
covPsf	<-mxAlgebra( expression= Asf+Csf+Esf, name="Vsf" )

# Define the matrices to hold the A and C effects: Common 
PathsAcsubm		<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.2,.1,.1), labels=c("ac22m","ac32m","ac33m"), name="acm" ) # Component paths for factors 2 and 3
PathsCcsubm		<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.3,.1,.01), labels=c("cc22m","cc32m","cc33m"), name="ccm" )
PathsEcsubm		<-mxMatrix(type="Diag", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.9,.3), labels=c("ec22m","ec33m"), name="ecm" )
PathsAcsubf		<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.6,.1,.3), labels=c("ac22f","ac32f","ac33f"), name="acf" ) # Component paths for factors 2 and 3
PathsCcsubf		<-mxMatrix(type="Lower", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.4,.1,.01), labels=c("cc22f","cc32f","cc33f"), name="ccf" )
PathsEcsubf		<-mxMatrix(type="Diag", nrow=nfact-1, ncol=nfact-1, free=TRUE, values=c(.7,.4), labels=c("ec22f","ec33f"), name="ecf" )

PathsP11m	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11m", name="pcm" ) # SD path for factor 3 (the PRS factor)
PathsP11f	<-mxMatrix(type="Full",  nrow=1, ncol=1, free=c(T), values=1, labels="pc11f", name="pcf" ) 

Ze11		<-mxMatrix(type="Zero",	nrow=1, ncol=1, free=F, name="Z11" )  #Padding
Ze21		<-mxMatrix(type="Zero",	nrow=2, ncol=1, free=F, name="Z21" )  #Padding
Ze12		<-mxMatrix(type="Zero",	nrow=1, ncol=2, free=F, name="Z12" )  #Padding

covAcsubm	<-mxAlgebra( expression= acm %*% t(acm), name="Acsubm" )
covCcsubm	<-mxAlgebra( expression= ccm %*% t(ccm), name="Ccsubm" )
covEcsubm	<-mxAlgebra( expression= ecm %*% t(ecm), name="Ecsubm" )
covAcsubf	<-mxAlgebra( expression= acf %*% t(acf), name="Acsubf" )
covCcsubf	<-mxAlgebra( expression= ccf %*% t(ccf), name="Ccsubf" )
covEcsubf	<-mxAlgebra( expression= ecf %*% t(ecf), name="Ecsubf" )

covPcsubm	<-mxAlgebra( expression= Acsubm+Ccsubm+Ecsubm, name="Vcsubm" ) #Matrix for the total variance of factors 2 and 3 (i.e. X and Y)
covPcsubf	<-mxAlgebra( expression= Acsubf+Ccsubf+Ecsubf, name="Vcsubf" ) 

covPc11m	<-mxAlgebra( expression= pcm %*% t(pcm), name="Pc11m" ) # variance for factor 1 (the PRS factor), I specify this separately as I do not want to resolve its variance into ACE components
covPc11f	<-mxAlgebra( expression= pcf %*% t(pcf), name="Pc11f" ) 

covPcm	<-mxAlgebra(cbind(rbind(Pc11m,Z21), rbind(Z12,Vcsubm)), name="Vcm") #I combine the PRS variance with the var-cov matrix of the other two factors.
covPcf	<-mxAlgebra(cbind(rbind(Pc11f,Z21), rbind(Z12,Vcsubf)), name="Vcf") 

covPcMzm	<-mxAlgebra(cbind(rbind(Pc11m,Z21) ,rbind(Z12,Acsubm+Ccsubm)), name="Vcmzm") #I specify the MZ between-twin covariance - excluding E parameters
covPcMzf	<-mxAlgebra(cbind(rbind(Pc11f,Z21) ,rbind(Z12,Acsubf+Ccsubf)), name="Vcmzf") 
covPcDzm	<-mxAlgebra(cbind(rbind(.5%x%Pc11m,Z21), rbind(Z12,.5%x%(acm %*% t(acm))+(ccm %*% t(ccm)))), name="Vcdzm") #I specify the DZ between-twin covariance - specifying half of A and excluding E
covPcDzf	<-mxAlgebra(cbind(rbind(.5%x%Pc11f,Z21), rbind(Z12,.5%x%(acf %*% t(acf))+(ccf %*% t(ccf)))), name="Vcdzf") 
covPcDzmf	<-mxAlgebra(cbind(rbind(.5%x%Pc11m,Z21), rbind(Z12,.5%x%(acm %*% t(acf))+(ccm %*% t(ccf)))), name="Vcdzmf") 
covPcDzfm	<-mxAlgebra(cbind(rbind(.5%x%Pc11f,Z21), rbind(Z12,.5%x%(acf %*% t(acm))+(ccf %*% t(ccm)))), name="Vcdzfm") 

# Generate Covariance of Latent factor model Including Causal Paths between factors
Id3		<-mxMatrix(type="Iden",	nrow=3, ncol=3, free=F, name="I3" )
covFVcm	<-mxAlgebra( expression= solve(I3-PhCm) %&% Vcm, name ="FVcm")
covFVcf	<-mxAlgebra( expression= solve(I3-PhCf) %&% Vcf, name ="FVcf")
covFcMzm	<-mxAlgebra( expression= solve(I3-PhCm) %&% Vcmzm, name ="Fcmzm")
covFcMzf	<-mxAlgebra( expression= solve(I3-PhCf) %&% Vcmzf, name ="Fcmzf")
covFcDzm	<-mxAlgebra( expression= solve(I3-PhCm) %&% Vcdzm, name ="Fcdzm")
covFcDzf	<-mxAlgebra( expression= solve(I3-PhCf) %&% Vcdzf, name ="Fcdzf")
covFcDzmf	<-mxAlgebra( expression= solve(I3-PhCm) %*% Vcdzmf %*% (I3-PhCf), name ="Fcdzmf")
covFcDzfm	<-mxAlgebra( expression= solve(I3-PhCf) %*% Vcdzfm %*% (I3-PhCm), name ="Fcdzfm")

# Constraint on total variance of Ordinal variable (A+C+E=1)
varL1m	<- mxConstraint( expression=FVcm[3,3]==1, name="L1m" )
varL1f	<- mxConstraint( expression=FVcf[3,3]==1, name="L1f" )

FcovMZM	<-mxAlgebra( expression= (FactLTwm  %&% rbind ( cbind(FVcm, Fcmzm), cbind(Fcmzm, FVcm))) , name="expFCovMZM" )
FcovMZF	<-mxAlgebra( expression= (FactLTwf  %&% rbind ( cbind(FVcf, Fcmzf), cbind(Fcmzf, FVcf))) , name="expFCovMZF" )
FcovDZM	<-mxAlgebra( expression= (FactLTwm  %&% rbind ( cbind(FVcm, Fcdzm), cbind(Fcdzm, FVcm))) , name="expFCovDZM" )
FcovDZF	<-mxAlgebra( expression= (FactLTwf  %&% rbind ( cbind(FVcf, Fcdzf), cbind(Fcdzf, FVcf))) , name="expFCovDZF" )
FcovDZO	<-mxAlgebra( expression= (FactLTwo  %&% rbind ( cbind(FVcm, Fcdzmf), cbind(Fcdzfm, FVcf))) , name="expFCovDZO" )

SpcovMZM	<-mxAlgebra( expression= rbind (cbind(Vsm, Asm+Csm), cbind(Asm+Csm, Vsm)) , name="expSpCovMZM" )
SpcovMZF	<-mxAlgebra( expression= rbind (cbind(Vsf, Asf+Csf), cbind(Asf+Csf, Vsf)) , name="expSpCovMZF" )
SpcovDZM	<-mxAlgebra( expression= rbind (cbind(Vsm, .5%x%Asm+Csm), cbind(.5%x%Asm+Csm, Vsm)) , name="expSpCovDZM" )
SpcovDZF	<-mxAlgebra( expression= rbind (cbind(Vsf, .5%x%Asf+Csf), cbind(.5%x%Asf+Csf, Vsf)) , name="expSpCovDZF" )
SpcovDZO	<-mxAlgebra( expression= rbind (cbind(Vsm, .5%x%(asm %*% t(asf))+(csm %*% t(csf))), cbind(.5%x%(asf %*% t(asm))+(csf %*% t(csm)), Vsf)) , name="expSpCovDZO" )

TOTcovMZM	<-mxAlgebra( expression= expFCovMZM + expSpCovMZM , name="TOTexpCovMZM" )
TOTcovMZF	<-mxAlgebra( expression= expFCovMZF + expSpCovMZF , name="TOTexpCovMZF" )
TOTcovDZM	<-mxAlgebra( expression= expFCovDZM + expSpCovDZM , name="TOTexpCovDZM" )
TOTcovDZF	<-mxAlgebra( expression= expFCovDZF + expSpCovDZF , name="TOTexpCovDZF" )
TOTcovDZO	<-mxAlgebra( expression= expFCovDZO + expSpCovDZO , name="TOTexpCovDZO" )

# *******************************************************************************************************
# Calculator

# Standardize the causal effects
Stcp1on2m	<-mxAlgebra( expression= (PhCm[2,1]* sqrt(FVcm[1,1]))/sqrt(FVcm[2,2]) , name="Stand_1on2m" )
Stcp1on3m	<-mxAlgebra( expression= (PhCm[3,1]* sqrt(FVcm[1,1]))/sqrt(FVcm[3,3]) , name="Stand_1on3m" )
Stcp2on3m	<-mxAlgebra( expression= (PhCm[3,2]* sqrt(FVcm[2,2]))/sqrt(FVcm[3,3]) , name="Stand_2on3m" )

Stcp1on2f	<-mxAlgebra( expression= (PhCf[2,1]* sqrt(FVcf[1,1]))/sqrt(FVcf[2,2]) , name="Stand_1on2f" )
Stcp1on3f	<-mxAlgebra( expression= (PhCf[3,1]* sqrt(FVcf[1,1]))/sqrt(FVcf[3,3]) , name="Stand_1on3f" )
Stcp2on3f	<-mxAlgebra( expression= (PhCf[3,2]* sqrt(FVcf[2,2]))/sqrt(FVcf[3,3]) , name="Stand_2on3f" )

# Standardize the Total var/covariances matrices of the observed variables
Id10		<-mxMatrix(type="Iden",	nrow=ntv, ncol=ntv, name="I10" )
Rfactmzm	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovMZM)) %&% TOTexpCovMZM, name="FactcorMZM" )
Rfactmzf	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovMZF)) %&% TOTexpCovMZF, name="FactcorMZF" )
Rfactdzm	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovDZM)) %&% TOTexpCovDZM, name="FactcorDZM" )
Rfactdzf	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovDZF)) %&% TOTexpCovDZF, name="FactcorDZF" )
Rfactdzo	<-mxAlgebra( expression= solve(sqrt(I10*TOTexpCovDZO)) %&% TOTexpCovDZO, name="FactcorDZO" )

# Phenotypic, A, C and E correlations	
RfactAcm	<-mxAlgebra( expression= solve(sqrt(I2*Acsubm)) %&% Acsubm, name="Ram" )
RfactCcm	<-mxAlgebra( expression= solve(sqrt(I2*Ccsubm)) %&% Ccsubm, name="Rcm" )
RfactEcm	<-mxAlgebra( expression= solve(sqrt(I2*Ecsubm)) %&% Ecsubm, name="Rem" )
RfactPm	<-mxAlgebra( expression= solve(sqrt(I3*FVcm)) %&% FVcm, name="Rphm" )

RfactAcf	<-mxAlgebra( expression= solve(sqrt(I2*Acsubf)) %&% Acsubf, name="Raf" )
RfactCcf	<-mxAlgebra( expression= solve(sqrt(I2*Ccsubf)) %&% Ccsubf, name="Rcf" )
RfactEcf	<-mxAlgebra( expression= solve(sqrt(I2*Ecsubf)) %&% Ecsubf, name="Ref" )
RfactPf	<-mxAlgebra( expression= solve(sqrt(I3*FVcf)) %&% FVcf, name="Rphf" )

# Standardize the Common Effects
covFVc22m	<-mxAlgebra( expression= FVcm[2:3,2:3], name ="FVc22m")
covFVc22f	<-mxAlgebra( expression= FVcf[2:3,2:3], name ="FVc22f")

stcovAcm	<-mxAlgebra( expression= Acsubm/FVc22m, name="stAcm" )
stcovCcm	<-mxAlgebra( expression= Ccsubm/FVc22m, name="stCcm" )
stcovEcm	<-mxAlgebra( expression= Ecsubm/FVc22m, name="stEcm" )
stcovAcf	<-mxAlgebra( expression= Acsubf/FVc22f, name="stAcf" )
stcovCcf	<-mxAlgebra( expression= Ccsubf/FVc22f, name="stCcf" )
stcovEcf	<-mxAlgebra( expression= Ecsubf/FVc22f, name="stEcf" )

# Standardised path estimates
StpathAcm	<-mxAlgebra( expression= (sqrt(stAcm)), name="stpacm" )
StpathCcm	<-mxAlgebra( expression= (sqrt(stCcm)), name="stpccm" )
StpathEcm	<-mxAlgebra( expression= (sqrt(stEcm)), name="stpecm" )
StpathAcf	<-mxAlgebra( expression= (sqrt(stAcf)), name="stpacf" )
StpathCcf	<-mxAlgebra( expression= (sqrt(stCcf)), name="stpccf" )
StpathEcf	<-mxAlgebra( expression= (sqrt(stEcf)), name="stpecf" )

# Algebra to compute Rph-A, Rph-C and Rph-E
RphA12m	<-mxAlgebra(expression=sqrt(stAcm[1,1])*Ram[2,1]*sqrt(stAcm[2,2]), name = 'Rpha12m')
RphC12m	<-mxAlgebra(expression=sqrt(stCcm[1,1])*Rcm[2,1]*sqrt(stCcm[2,2]), name = 'Rphc12m')
RphE12m	<-mxAlgebra(expression=sqrt(stEcm[1,1])*Rem[2,1]*sqrt(stEcm[2,2]), name = 'Rphe12m')

RphA12f	<-mxAlgebra(expression=sqrt(stAcf[1,1])*Raf[2,1]*sqrt(stAcf[2,2]), name = 'Rpha12f')
RphC12f	<-mxAlgebra(expression=sqrt(stCcf[1,1])*Rcf[2,1]*sqrt(stCcf[2,2]), name = 'Rphc12f')
RphE12f	<-mxAlgebra(expression=sqrt(stEcf[1,1])*Ref[2,1]*sqrt(stEcf[2,2]), name = 'Rphe12f')

# Standardize the Specific Effects
stcovAsm	<-mxAlgebra( expression= sqrt(Asm/( (FactLm %&% FVcm) +Vsm)), name="stAsm" )
stcovCsm	<-mxAlgebra( expression= sqrt(Csm/( (FactLm %&% FVcm) +Vsm)), name="stCsm" )
stcovEsm	<-mxAlgebra( expression= sqrt(Esm/( (FactLm %&% FVcm) +Vsm)), name="stEsm" )

stcovAsf	<-mxAlgebra( expression= sqrt(Asf/( (FactLf %&% FVcf) +Vsf)), name="stAsf" )
stcovCsf	<-mxAlgebra( expression= sqrt(Csf/( (FactLf %&% FVcf) +Vsf)), name="stCsf" )
stcovEsf	<-mxAlgebra( expression= sqrt(Esf/( (FactLf %&% FVcf) +Vsf)), name="stEsf" )

# Standardized Factor Loadings
StFLm		<-mxAlgebra( expression= sqrt(diag2vec( FactLm %&% FVcm / TOTexpCovMZM[1:5,1:5])) , name="StandFactm" )
StFLf		<-mxAlgebra( expression= sqrt(diag2vec( FactLf %&% FVcf / TOTexpCovMZF[1:5,1:5])) , name="StandFactf" )

# *******************************************************************************************************

# Data objects for Multiple Groups
dataMZM	<- mxData( observed=mzmData, type="raw" )
dataMZF	<- mxData( observed=mzfData, type="raw" )
dataDZM	<- mxData( observed=dzmData, type="raw" )
dataDZF	<- mxData( observed=dzfData, type="raw" )
dataDZO	<- mxData( observed=dzoData, type="raw" )

# Objective objects for Multiple Groups
objMZM	<- mxExpectationNormal( covariance="TOTexpCovMZM", means="expMeanm", dimnames=selVars, thresholds="expThresm", threshnames=c("SO1","SO2"))
objMZF	<- mxExpectationNormal( covariance="TOTexpCovMZF", means="expMeanf", dimnames=selVars, thresholds="expThresf", threshnames=c("SO1","SO2"))
objDZM	<- mxExpectationNormal( covariance="TOTexpCovDZM", means="expMeanm", dimnames=selVars, thresholds="expThresm", threshnames=c("SO1","SO2"))
objDZF	<- mxExpectationNormal( covariance="TOTexpCovDZF", means="expMeanf", dimnames=selVars, thresholds="expThresf", threshnames=c("SO1","SO2"))
objDZO	<- mxExpectationNormal( covariance="TOTexpCovDZO", means="expMeano", dimnames=selVars, thresholds="expThreso", threshnames=c("SO1","SO2"))

fitFunction <- mxFitFunctionML()
 
# Combine Groups
pars1m	<-list(Meansm,Loadm,LoadTwm,PhCausm,PathsAsm,PathsCsm,PathsEsm,covAsm,covCsm,covEsm,covPsm,Id2,Id3,Id10)
pars1f	<-list(Meansf,Loadf,LoadTwf,PhCausf,PathsAsf,PathsCsf,PathsEsf,covAsf,covCsf,covEsf,covPsf,Id2,Id3,Id10)
pars1o	<-list(Meanso,Ze53, LoadTwo)
pars2m	<-list(PathsAcsubm,PathsCcsubm,PathsEcsubm,PathsP11m,Ze21,Ze12,Ze11,covAcsubm,covCcsubm,covEcsubm,covPcsubm,covPc11m,covPcm,covPcMzm,covPcDzm,covFVcm,covFVc22m,covFcMzm,covFcDzm)
pars2f	<-list(PathsAcsubf,PathsCcsubf,PathsEcsubf,PathsP11f,Ze21,Ze12,Ze11,covAcsubf,covCcsubf,covEcsubf,covPcsubf,covPc11f,covPcf,covPcMzf,covPcDzf,covFVcf,covFVc22f,covFcMzf,covFcDzf)
pars2o	<-list(covPcDzmf,covPcDzfm,covFcDzmf,covFcDzfm)
pars3m	<-list(obsage1, obsage2, betaAm, Thrm, inc, Thresm)
pars3f	<-list(obsage1, obsage2, betaAf, Thrf, inc, Thresf)
parsstm	<-list(stcovAsm, stcovCsm, stcovEsm, stcovAcm, stcovCcm, stcovEcm, RfactAcm, RfactCcm, RfactEcm, RfactPm,RphA12m,RphC12m,RphE12m,StpathAcm,StpathCcm,StpathEcm)
parsstf	<-list(stcovAsf, stcovCsf, stcovEsf, stcovAcf, stcovCcf, stcovEcf, RfactAcf, RfactCcf, RfactEcf, RfactPf,RphA12f,RphC12f,RphE12f,StpathAcf,StpathCcf,StpathEcf)
parsmedm	<-list(Stcp1on2m, Stcp1on3m, Stcp2on3m)
parsmedf	<-list(Stcp1on2f, Stcp1on3f, Stcp2on3f)

modelMZM	<-mxModel(pars1m, pars2m, pars3m, parsmedm, FcovMZM, SpcovMZM, TOTcovMZM, dataMZM, objMZM, Rfactmzm, parsstm, fitFunction, StFLm, varL1m, name="MZM" )
modelMZF	<-mxModel(pars1f, pars2f, pars3f, parsmedf, FcovMZF, SpcovMZF, TOTcovMZF, dataMZF, objMZF, Rfactmzf, parsstf, fitFunction, StFLf, varL1f, name="MZF" )
modelDZM	<-mxModel(pars1m, pars2m, pars3m, FcovDZM, SpcovDZM, TOTcovDZM, dataDZM, objDZM, Rfactdzm, fitFunction, name="DZM" )
modelDZF	<-mxModel(pars1f, pars2f, pars3f, FcovDZF, SpcovDZF, TOTcovDZF, dataDZF, objDZF, Rfactdzf, fitFunction, name="DZF" )
modelDZO	<-mxModel(pars1m, pars2m, pars1f, pars2f, pars1o, pars2o, pars3m, pars3f, Threso, FcovDZO, SpcovDZO, TOTcovDZO, dataDZO, objDZO, Rfactdzo, fitFunction, name="DZO" )

minus2ll	<-mxAlgebra( expression=MZM.objective + MZF.objective + DZM.objective + DZF.objective + DZO.objective, name="m2LL" )
obj		<-mxFitFunctionAlgebra( "m2LL" )

cistFLm	<-mxCI (c ('MZM.StandFactm','MZM.Stand_1on2m','MZM.Stand_1on3m','MZM.Stand_2on3m','MZM.PhCm'))
cistFLf	<-mxCI (c ('MZF.StandFactf','MZF.Stand_1on2f','MZF.Stand_1on3f','MZF.Stand_2on3f','MZF.PhCf'))

cistVsm	<-mxCI (c ('MZM.stAsm[3,3]','MZM.stCsm[3,3]','MZM.stEsm[1,1]','MZM.stEsm[3,3]') ) 	# standardized var comp from specific Factors
cistVsf	<-mxCI (c ('MZF.stAsf[3,3]','MZF.stCsf[3,3]','MZF.stEsf[1,1]','MZF.stEsf[3,3]') ) 	# standardized var comp from specific Factors

cistVcm	<-mxCI (c ('MZM.stAcm[1,1]','MZM.stAcm[2,1]','MZM.stAcm[2,2]',
				'MZM.stCcm[1,1]','MZM.stCcm[2,1]','MZM.stCcm[2,2]',
				'MZM.stEcm[1,1]','MZM.stEcm[2,2]') ) 	# standardized var comp for ACE on latent Factors
cistVcf	<-mxCI (c ('MZF.stAcf[1,1]','MZF.stAcf[2,1]','MZF.stAcf[2,2]',
				'MZF.stCcf[1,1]','MZF.stCcf[2,1]','MZF.stCcf[2,2]',
				'MZF.stEcf[1,1]','MZF.stEcf[2,2]') ) 	# standardized var comp for ACE on latent Factors

cistRcm	<-mxCI (c ('MZM.Rpha12m','MZM.Rphc12m','MZM.Rphe12m','MZM.Ram','MZM.Rcm','MZM.Rem','MZM.stpacm','MZM.stpccm','MZM.stpecm') ) 	
cistRcf	<-mxCI (c ('MZF.Rpha12f','MZF.Rphc12f','MZF.Rphe12f','MZF.Raf','MZF.Rcf','MZF.Ref','MZF.stpacf','MZF.stpccf','MZF.stpecf') ) 	

HetACEMs2Model	<-mxModel("HetaceMs2", pars1m, pars2m, pars1f, pars2f, pars1o, pars2o, modelMZM, modelMZF, modelDZM, modelDZF, modelDZO, minus2ll, obj, cistFLm, cistVsm, cistVcm, cistRcm, cistFLf, cistVsf, cistVcf, cistRcf) 

# --------------------------------------------------------------------------------------------------------------------------------
# 5aii RUN HetACEMs2 Factor Model with phenotypic causal mediation paths by Zygosity

HetACEMs2Fit	<-mxTryHardOrdinal(HetACEMs2Model, intervals=T)
(HetACEMs2Summ	<-summary(HetACEMs2Fit, verbose=F))

# Generate confidence intervals
HetACEMs2CIModel	<-mxModel(HetACEMs2Model)
HetACEMs2CIFit	<-mxRun(HetACEMs2CIModel, intervals=TRUE)
(HetACEMs2CISumm	<-summary(HetACEMs2CIFit, verbose=F))

mxEval(MZM.FactcorMZM, HetACEMs2Fit)
mxEval(DZM.FactcorDZM, HetACEMs2Fit)

mxEval(MZF.FactcorMZF, HetACEMs2Fit)
mxEval(DZF.FactcorDZF, HetACEMs2Fit)

mxEval(MZM.FVcm, HetACEMs2Fit)
mxEval(MZF.FVcf, HetACEMs2Fit)

mxEval(MZM.stAcm, HetACEMs2Fit)
mxEval(MZM.stCcm, HetACEMs2Fit)
mxEval(MZM.stEcm, HetACEMs2Fit)

mxEval(MZF.stAcf, HetACEMs2Fit)
mxEval(MZF.stCcf, HetACEMs2Fit)
mxEval(MZF.stEcf, HetACEMs2Fit)

mxEval(MZM.stpacm, HetACEMs2Fit)
mxEval(MZM.stpccm, HetACEMs2Fit)
mxEval(MZM.stpecm, HetACEMs2Fit)

mxEval(MZF.stpacf, HetACEMs2Fit)
mxEval(MZF.stpccf, HetACEMs2Fit)
mxEval(MZF.stpecf, HetACEMs2Fit)

mxEval(MZM.Ram, HetACEMs2Fit)
mxEval(MZM.Rcm, HetACEMs2Fit)
mxEval(MZM.Rem, HetACEMs2Fit)
mxEval(MZM.Rphm, HetACEMs2Fit)

mxEval(MZF.Raf, HetACEMs2Fit)
mxEval(MZF.Rcf, HetACEMs2Fit)
mxEval(MZF.Ref, HetACEMs2Fit)
mxEval(MZF.Rphf, HetACEMs2Fit)

mxEval(MZM.Rpham, HetACEMs2Fit) 
mxEval(MZM.Rphcm, HetACEMs2Fit) 
mxEval(MZM.Rphem, HetACEMs2Fit) 

mxEval(MZF.Rphaf, HetACEMs2Fit) 
mxEval(MZF.Rphcf, HetACEMs2Fit) 
mxEval(MZF.Rphef, HetACEMs2Fit) 

mxEval(MZM.Acsubm, HetACEMs2Fit)
mxEval(MZM.Ccsubm, HetACEMs2Fit)
mxEval(MZM.Ecsubm, HetACEMs2Fit)
mxEval(MZM.Vcsubm, HetACEMs2Fit)
mxEval(MZM.FVc22m, HetACEMs2Fit)

mxEval(MZF.Acsubf, HetACEMs2Fit)
mxEval(MZF.Ccsubf, HetACEMs2Fit)
mxEval(MZF.Ecsubf, HetACEMs2Fit)
mxEval(MZF.Vcsubf, HetACEMs2Fit)
mxEval(MZF.FVc22f, HetACEMs2Fit)


mxEval(MZM.Stand_1on2m[1,1], HetACEMs2Fit)
mxEval(MZM.Stand_1on3m[1,1], HetACEMs2Fit)
mxEval(MZM.Stand_2on3m[1,1], HetACEMs2Fit)

mxEval(MZF.Stand_1on2f[1,1], HetACEMs2Fit)
mxEval(MZF.Stand_1on3f[1,1], HetACEMs2Fit)
mxEval(MZF.Stand_2on3f[1,1], HetACEMs2Fit)

#------------------------------------------------------------------------
# Submodel 5aiim: Drop correlation path male c paths from previous model 
#------------------------------------------------------------------------
# Drop correlation male c paths from the model
# -----------------------------------------------------------------------
HetAEMs2mModel	<- mxModel(HetACEMs2Fit, name="HetAEMs2m")
HetAEMs2mModel	<- omxSetParameters(HetAEMs2mModel, labels=c('cc22m','cc32m','cc33m'), free=FALSE, values=0)
HetAEMs2mModel	<- omxSetParameters(HetAEMs2mModel, labels=c('cs3m'), free=FALSE, values=0)
HetAEMs2mFit	<- mxRun(HetAEMs2mModel, intervals=F)
(HetAEMs2mSum	<- summary(HetAEMs2mFit))

mxCompare(HetACEMs2Fit, HetAEMs2mFit)

#------------------------------------------------------------------------
# Submodel 5aiif: Drop correlation path female c paths from previous model 
#------------------------------------------------------------------------
# Drop correlation female c paths from the model
# -----------------------------------------------------------------------
HetAEMs2fModel	<- mxModel(HetACEMs2Fit, name="HetAEMs2f")
HetAEMs2fModel	<- omxSetParameters(HetAEMs2fModel, labels=c('cc22f','cc32f','cc33f'), free=FALSE, values=0)
HetAEMs2fModel	<- omxSetParameters(HetAEMs2fModel, labels=c('cs3f'), free=FALSE, values=0)
HetAEMs2fFit	<- mxRun(HetAEMs2fModel, intervals=F)
(HetAEMs2fSum	<- summary(HetAEMs2fFit))

mxCompare(HetACEMs2Fit, HetAEMs2fFit)


#------------------------------------------------------------------------
# Submodel 5aiif: Equate causal paths in males and females 
#------------------------------------------------------------------------
# Equate causal paths in males and females
# -----------------------------------------------------------------------
HetAEMs2cModel	<- mxModel(HetACEMs2Fit, name="HetAEMs2c")
HetAEMs2cModel	<- omxSetParameters(HetAEMs2cModel, labels=c('c2on3f'), newlabels=c('c2on3m'), free=T, values=.05)
HetAEMs2cModel	<- omxAssignFirstParameters(HetAEMs2cModel)
HetAEMs2cFit	<- mxRun(HetAEMs2cModel, intervals=F)
(HetAEMs2cSum	<- summary(HetAEMs2cFit))

mxCompare(HetACEMs2Fit, HetAEMs2cFit)

