#WCOP2016 simulation PK example using PKADVAN package
	#Before running the PK example, you need to:
		#(1) Download and install the PKADVAN package from GitHib ( https://github.com/abuhelwa/PKADVAN_Rpackage )
		#(2) intall the plyr package.
		
	rm(list=ls(all=TRUE))
	graphics.off()

	library(PKADVAN)
	library(plyr)
	library(ggplot2)

	master.dir <- "D:/AAbuhelwa/Projects/6_ADVAN_Derivation_Testing/WCOP216_example"	
#--------------------------------------------------------------
#Customize ggplot2 theme -
	theme_bw2 <- theme_set(theme_bw(base_size = 20))
	theme_bw2 <- theme_update(plot.margin = unit(c(1.5,1.5,3,1.5), "lines"),
							  axis.title.x=element_text(size = 18, vjust = 0),
							  axis.title.y=element_text(size = 18, vjust = 0, angle = 90),
							  strip.text.x = element_text(size=12, margin = margin(3, 0, 3, 0)),
							  strip.text.y=element_text(size = 16, angle = 90))

#--------------------------------------------------------------
#Confidence intervals - from function utility

	CI90lo <- function(x) quantile(x, probs=0.05,na.rm=T)
	CI90hi <- function(x) quantile(x, probs=0.95,na.rm=T)

#########################################################
## Three steps for performing simulations:
# 1.  Generate a simulation data frame with the individual PK parameters and include any covariate effects on the PK parameters.
# 2.  Run the "PKADVAN" function for the selected model: the "PKADVAN" function will return the simulated amounts in each compartment of the PK system and IPREDs for the central compartment.
# 3.  Add residual un-explained variability on the IPREDs.

###################################################################################################
######## EXAMPLE 1:     Two compartment two transit absorption model       ########################
###################################################################################################
#=====================================
#   Generate simulation data frame
#=====================================
# If you have a NONMEM df ready then just read it and skip this step)
	#Set dose records:
	dosetimes <- c(seq(0,48,12))
	tlast <- 96
	
	#Now define finer sample times for after a dose to capture Cmax
	doseseq <- c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8,9,10)
	
	#Use the outer product but with addition to expand this doseseq for all dosetimes
	PKtimes <- outer(dosetimes,doseseq,"+")
	
	#set number of subjects
	nsub <- 1000
	ID <- 1:nsub
	
	#Make dataframe
	df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),PKtimes))),"AMT"=0,"MDV"=0,"DV"=NA,"CLCR"=120)
	df$CLCR[df$TIME <= 48] <- 90
	
	doserows <- subset(df, TIME%in%dosetimes)
	
	#Dose: It can be any arbitrary dose
	doserows$AMT <- 500
	doserows$MDV <- 1
	
	#Add back dose information
	df <- rbind(df,doserows)
	df <- df[order(df$ID,df$TIME,df$AMT),]       # arrange df by TIME (ascending) and by AMT (descending)
	df <- subset(df, (TIME==0 & AMT==0)==F) # remove the row that has a TIME=0 and AMT=0
	
#----------------------------------------------------------------------------------------------------
# Two compartment-2 transit first-order absorption model via PKADVAN package
#----------------------------------------------------------------------------------------------------
#Define between subject variability on PK parameters
	#BSV (Omegas as SD)
	ETA1CL  <- 0.15     
	ETA2V2  <- 0.12		
	ETA3Q   <- 0.14		
	ETA4V3  <- 0.05		
	ETA5KTR <- 0.30    

#Define residual error model
	#Residuals (Epsilons as SD)
	EPS1	<- 0.10		#Proportional residual error
	EPS2	<- 0.15		#Additive residual error

	#Use random number generator to simulate residuals from a normal distribution    
	BSVCL	<- rnorm(nsub, mean = 0, sd = ETA1CL)	#BSV on CL  
	BSVV2	<- rnorm(nsub, mean = 0, sd = ETA2V2)	#BSV on V2  
	BSVQ	<- rnorm(nsub, mean = 0, sd = ETA3Q)	#BSV on Q
	BSVV3	<- rnorm(nsub, mean = 0, sd = ETA4V3)	#BSV on V3
	BSVKTR	<- rnorm(nsub, mean = 0, sd = ETA5KTR)  #BSV on KTR

	EPS1 <- rnorm(nsub, mean = 0, sd = EPS1)	#Proportional residual error
	EPS2 <- rnorm(nsub, mean = 0, sd = EPS2)	#Additive residual error
	
#Set population PK parameters for 2-compartment First-order absorption model
	CLpop   <- 0.5      #clearance
	V2pop   <- 20       #central volume of distribution
	Qpop    <- 1        #inter-compartmental clearance
	V3pop   <- 25       #peripheral volume of distribution
	KTRpop  <- 2.05     #first-order absorption rate constant
	F1pop   <- 0.80     #Bioavailability

#Modify df for ADVAN calculation and include any covariates on PK parameters
	dfadvan <- df
	#Calculate group parameter values including any covariate effects
	dfadvan$CL  <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)   #creatinine clearance (CLCR) added as a time-changing covariate on CL
	dfadvan$V2  <- V2pop*exp(BSVV2)
	dfadvan$Q   <- Qpop*exp(BSVQ)
	dfadvan$V3  <- V3pop*exp(BSVV3)
	dfadvan$KTR <- KTRpop*exp(BSVKTR)
	dfadvan$F1  <- F1pop
	
	#Apply PKADVAN functions
	simdf <- ddply(dfadvan, .(ID), TwoCompTwoTransit)
	head(simdf)
	
#Add residual unexplained variability (within subject variability)
	#For example; Additive error model
	simdf$DV <- simdf$IPRED*(1 + EPS1) + EPS2

#Plotting	
	#Subset missing data
	simdf <- subset(simdf, MDV==0)
	simdf$SOURCE <- "PKADVAN-R package"
    
	#change working directory
	setwd(paste(master.dir,"/2comp_oral_2tranist.nm7",sep=""))
	
	#IPRED
	plotobj <- NULL
	titletext <- expression(atop("Simulated Drug Concentrations",
	                             atop(italic("Two compartment- two tranist first-order absorption model"),"Median and 90% Prediction Interval, 1000 subjects"))) 
	plotobj <- ggplot(data=simdf)
	plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
	plotobj <- plotobj + ggtitle(titletext)
	plotobj <- plotobj + scale_y_continuous("Concentration\n")
	plotobj <- plotobj + scale_x_continuous("\nTime after dose")
	plotobj <- plotobj + facet_wrap(~SOURCE)
	plotobj

    #DV
	plotobj <- NULL
	titletext <- expression(atop("Simulated Drug Concentrations",
	                             atop(italic("Two compartment- two tranist first-order absorption model"),"Median and 90% Prediction Interval, 1000 subjects"))) 
	
	plotobj <- ggplot(data=simdf)
	plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
	plotobj <- plotobj + ggtitle(titletext)
	plotobj <- plotobj + scale_y_continuous("Concentration\n")
	plotobj <- plotobj + scale_x_continuous("\nTime after dose")
	plotobj <- plotobj + facet_wrap(~SOURCE)
	plotobj

###################################################################################################
######## EXAMPLE 2:     One compartment IV bolus with 1 comp metabolite    ########################
###################################################################################################
#================================
# Generate simulation data frame
#================================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
	#Set dose records:
	dosetimes <- c(seq(0,48,24))
	tlast <- 96
	#set number of subjects
	nsub <- 1000
	ID  <- 1:nsub
#Make dataframe: CLCR: is creatinine clearance
	df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),dosetimes))),"AMT"=0,"MDV"=0,"CLCR"=90)
	df$CLCR[df$ID >= 0.5*nsub] <- 70

#Dose = 100 mg, Dose 1  at time 0
	doserows <- subset(df, TIME %in% dosetimes)
	doserows$AMT <- 100
	doserows$MDV <- 1

#Add back dose information
	df <- rbind(df,doserows)
	df <- df[order(df$ID,df$TIME,df$AMT),]       # arrange df by TIME (ascending) and by AMT (descending)
	df <- subset(df, (TIME==0 & AMT==0)==FALSE)      # remove the row that has a TIME=0 and AMT=0

#--------------------------------------------------------------
# 1 compartment-IV Bolus 1 comp metabolite via PKADVAN package
#--------------------------------------------------------------
#Define between subject variability (BSV)
	#Use random number generator to simulate residuals from a normal distribution
	BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL
	BSVV  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V

#Define residual error model
	PROP 	<- rnorm(nsub, mean = 0, sd = 0.22)

#Define population PK parameters for 1-compartment model
	CLpop <- 0.2       # clearance
	Vpop  <- 5      # central volume of distribution

#Modify df for ADVAN calculation and include any covariate effects on the PK parameters
	dfadvan <- df
	#Calculate group parameter values including any covariate effects
	dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)      #creatinine clearance (CLCR) added as a covariate on CL
	dfadvan$V  <- Vpop*exp(BSVV)

#Apply PKADVAN function
	simdf <- ddply(dfadvan, .(ID), OneCompIVbolus)

#Add residual unexplained variability (within subject variability)
	#For example: proportional error model on IPRED
	simdf$DV <- simdf$IPRED*(1+PROP)
	head(simdf)

#plotting
	#Subset MDV==0
	simdf <- subset(simdf, MDV==0)

	#IPRED
	plotobj <- NULL
	titletext <- expression(atop('Simulated Drug Concentrations',
								 atop(italic("Red line is the median. Gray band is the 90% confidence interval")
								 )))
	plotobj <- ggplot(data=simdf)
	plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
	plotobj <- plotobj + ggtitle(titletext)
	plotobj <- plotobj + scale_y_continuous("Concentration\n")
	plotobj <- plotobj + scale_x_continuous("\nTime after dose")
	plotobj

	#DV
	plotobj <- NULL
	titletext <- expression(atop('Simulated Drug Concentrations',
								 atop(italic("Red line is the median. Gray band is the 90% confidence interval")
								 )))
	plotobj <- ggplot(data=simdf)
	plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
	plotobj <- plotobj + ggtitle(titletext)
	plotobj <- plotobj + scale_y_continuous("Concentration\n")
	plotobj <- plotobj + scale_x_continuous("\nTime after dose")
	plotobj

	
#--------------------------------------------
# 2 compartment-IV Bolus via PKADVAN package
#--------------------------------------------
#Set population PK parameters for 2-compartment model
CLpop <- 2       # clearance
V1pop <- 10      # central volume of distribution
Qpop  <- 1       # inter-compartmental clearance
V2pop <- 25      # peripheral volume of distribution

#Modify df for ADVAN calculation and include any covariate effects on the PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*(dfadvan$CLCR/100)           #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V1 <- V1pop
dfadvan$Q  <- Qpop
dfadvan$V2 <- V2pop

#apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), TwoCompIVbolus)
head(simdf)

system.time(ddply(dfadvan, .(ID), TwoCompIVbolus))

#--------------------------------------------
# 3 compartment-IV Bolus via PKADVAN package
#--------------------------------------------
#Set population PK parameters for 1-compartment model
CLpop  <- 2          # clearance
V1pop  <- 10         # central volume of distribution
Q12pop <- 0.5        # inter-compartmental clearance (1)
V2pop  <- 30         # peripheral volume of distribution (1)
Q13pop <- 0.3        # inter-compartmental clearance (2)
V3pop  <- 40         # peripheral volume of distribution (2)

#Modify df for ADVAN calculation and include any covariate effects on the PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL  <- CLpop*(dfadvan$CLCR/100)    #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V1  <- V1pop
dfadvan$Q2 <- Q12pop
dfadvan$V2  <- V2pop
dfadvan$Q3 <- Q13pop
dfadvan$V3  <- V3pop

#Apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), ThreeCompIVbolus)
head(simdf)

#?processing time
system.time(ddply(dfadvan, .(ID), ThreeCompIVbolus))

#================================
# First-order absorption Models
#================================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
#Set dose records:
dosetimes <- c(seq(0,48,12))
tlast <- 72

#Now define finer sample times for after a dose to capture Cmax
doseseq <- c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8,9,10)

#Use the outer product but with addition to expand this doseseq for all dosetimes
PKtimes <- outer(dosetimes,doseseq,"+")

#set number of subjects
nsub <- 1000
ID <- 1:nsub

#Make dataframe
df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),PKtimes))),"AMT"=0,"MDV"=0,"CLCR"=90)
df$CLCR[df$ID <= 0.5*nsub] <- 70

doserows <- subset(df, TIME%in%dosetimes)

#Dose = 100 mg. It can be any arbitrary dose
doserows$AMT <- 100
doserows$MDV <- 1

#Add back dose information
df <- rbind(df,doserows)
df <- df[order(df$ID,df$TIME,df$AMT),]       # arrange df by TIME (ascending) and by AMT (descending)
df <- subset(df, (TIME==0 & AMT==0)==F) # remove the row that has a TIME=0 and AMT=0

#---------------------------------------------------------------
# 1 compartment-first order absorption via PKADVAN package
#---------------------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL
BSVF1 <- rnorm(nsub, mean = 0, sd = 0.2)  #BSV on F

#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.22)

#Define population PK parameters for 1-compartment model
CLpop <- 0.5       # clearance
Vpop  <- 10      # central volume of distribution
KApop <- 0.4     # first-order absorption rate constant
F1pop <- 1       # bioavailability

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)      #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V  <- Vpop
dfadvan$KA <- KApop
dfadvan$F1 <- F1pop*exp(BSVF1)

#Apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), OneCompFirstOrderAbs)
head(simdf)

#Add residual error model
#For example: proportional error model on IPRED
simdf$DV <- simdf$IPRED*(1+PROP)
head(simdf)

#plotting
#subset data for MDV ==0
simdf <- subset(simdf, MDV==0)

#IPRED
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#DV
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#?processing time. Depends on PK sampling duration
system.time(ddply(dfadvan, .(ID), OneCompFirstOrderAbs))

#----------------------------------------------------------------
# 2 compartment-first order absorption  via PKADVAN package
#----------------------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL
BSVF1 <- rnorm(nsub, mean = 0, sd = 0.2)  #BSV on F

#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.22)

#Set population PK parameters for 2-compartment model
CLpop  <- 0.5       # clearance
V2pop  <- 10     # central volume of distribution
Qpop  <- 1       # inter-compartmental clearance
V3pop <- 25      # peripheral volume of distribution
KApop <- 0.4     # first-order absorption rate constant
F1pop <- 1       # bioavailability

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df

#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)           #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V2 <- V2pop
dfadvan$Q  <- Qpop
dfadvan$V3 <- V3pop
dfadvan$KA <- KApop
dfadvan$F1 <- F1pop*exp(BSVF1)

#apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), TwoCompFirstOrderAbs)
head(simdf)

#Add residual error model
#For example: proportional error model on IPRED
simdf$DV <- simdf$IPRED*(1+PROP)
head(simdf)

#plotting
#subset data for MDV ==0
simdf <- subset(simdf, MDV==0)

#IPRED
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#DV
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#?processing time
system.time(ddply(dfadvan, .(ID), TwoCompFirstOrderAbs))

#---------------------------------------------------------------
# 3 compartment-first order absorption via PKADVAN package
#---------------------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL
BSVF1 <- rnorm(nsub, mean = 0, sd = 0.2)  #BSV on F

#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.22)

#Set population PK parameters for 3-compartment model
CLpop <- 0.5        # clearance
V2pop <- 20          # central volume of distribution
Q3pop <- 0.5        # inter-compartmental clearance (1)
V3pop <- 30         # peripheral volume of distribution (1)
Q4pop <- 0.3        # inter-compartmental clearance (2)
V4pop <- 40         # peripheral volume of distribution (2)

KApop <- 0.4        # first-order absorption rate constant
F1pop <- 0.5          # bioavailability

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)    #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V2 <- V2pop
dfadvan$Q3 <- Q3pop
dfadvan$V3 <- V3pop
dfadvan$Q4 <- Q4pop
dfadvan$V4 <- V4pop
dfadvan$KA <- KApop
dfadvan$F1 <- F1pop*exp(BSVF1)

#apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), ThreeCompFirstOrderAbs)
head(simdf)

#Add residual error model
#For example: proportional error model on IPRED
simdf$DV <- simdf$IPRED*(1+PROP)
head(simdf)

#plotting
#subset data for MDV ==0
simdf <- subset(simdf, MDV==0)

#IPRED
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#DV
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#?processing time
system.time(ddply(dfadvan, .(ID), ThreeCompFirstOrderAbs))

#=============================
# IV infusion Basic Models
#=============================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
#Set dose records:
#number of subjects
nsub <- 1000
ID <- 1:nsub

tlast <- 72
#Set dose records
dosetimessim <- c(0,20.5)

#Make dataframe
df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(0.5,seq(0,tlast,1),dosetimessim))),"AMT"=0,MDV=0,"RATE"=0,"CLCR"=90)

df$CLCR[df$ID>=0.5*nsub] <- 70

#subset doserows
doserowssim <- subset(df, TIME%in%dosetimessim)

#Dose & RATE can be arbitrary
doserowssim$AMT  <- 100
doserowssim$RATE <- 4
doserowssim$MDV  <- 1

#Add back dose information
df <- rbind(df,doserowssim)
df <- df[order(df$TIME,-df$AMT),]       # arrange df by TIME with the AMT from high to low
df <- subset(df, (TIME==0 & AMT==0)==F) # remove the row that has a TIME=0 and AMT=0

#----------------------------------------------------
# 1 compartment-IV infusion via PKADVAN package
#----------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL
BSVV  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V

#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.09)

#Define population PK parameters for 1-compartment model
CLpop <- 2       # clearance
Vpop  <- 10      # central volume of distribution

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)      #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V  <- Vpop*exp(BSVV)

#apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), OneCompIVinfusion)
head(simdf)

#Add residual unexplained variability (within subject variability)
#For example; additive error model
simdf$DV <- simdf$IPRED+PROP

#plotting
#subset data for MDV==0
simdf <- subset(simdf, MDV==0)

#IPRED
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#DV
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#?provcessing time
system.time(simdf <- ddply(dfadvan, .(ID), OneCompIVinfusion))

#--------------------------------------------------
# 2 compartment-IV infusion via PKADVAN package
#--------------------------------------------------
#Set population PK parameters for 1-compartment model
CLpop <- 2       # clearance
V1pop <- 10      # central volume of distribution
Qpop  <- 1       # inter-compartmental clearance
V2pop <- 25      # peripheral volume of distribution

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*(dfadvan$CLCR/100)           #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V1 <- V1pop
dfadvan$Q  <- Qpop
dfadvan$V2 <- V2pop

#Apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), TwoCompIVinfusion)

#system.time(simdf <- ddply(dfadvan, .(ID), TwoCompIVinfusion))

#-------------------------------------------------
# 3 compartment-IV infusion via PKADVAN package
#-------------------------------------------------
#Set population PK parameters for 1-compartment model
CLpop  <- 2          # clearance
V1pop  <- 10         # central volume of distribution
Q12pop <- 0.5        # inter-compartmental clearance (1)
V2pop  <- 30         # peripheral volume of distribution (1)
Q13pop <- 0.3        # inter-compartmental clearance (2)
V3pop  <- 40         # peripheral volume of distribution (2)

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
#Calculate group parameter values including any covariate effects
dfadvan$CL  <- CLpop*(dfadvan$CLCR/100)    #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V1  <- V1pop
dfadvan$Q2 <- Q12pop
dfadvan$V2  <- V2pop
dfadvan$Q3 <- Q13pop
dfadvan$V3  <- V3pop

#Apply PKADVAN function
simdf <- ddply(dfadvan, .(ID), ThreeCompIVinfusion)
#system.time(simdf <- ddply(dfadvan, .(ID), ThreeCompIVinfusion))

#=======================================================================
#   First-order absorption first-order metabolite formation models
#=======================================================================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
#Set dose records:
dosetimes <- c(seq(0,48,12))
tlast <- 96

#Now define finer sample times for after a dose to capture Cmax
doseseq <- c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8,9,10)

#Use the outer product but with addition to expand this doseseq for all dosetimes
PKtimes <- outer(dosetimes,doseseq,"+")

#set number of subjects
nsub <- 100
ID <- 1:nsub

#Make dataframe
df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),PKtimes))),"AMT"=0,"MDV"=0,"DV"=NA,"CLCR"=120)
df$CLCR[df$TIME <= 48] <- 90

doserows <- subset(df, TIME%in%dosetimes)

#Dose: It can be any arbitrary dose
doserows$AMT <- 500
doserows$MDV <- 1

#Add back dose information
df <- rbind(df,doserows)
df <- df[order(df$ID,df$TIME,df$AMT),]       # arrange df by TIME (ascending) and by AMT (descending)
df <- subset(df, (TIME==0 & AMT==0)==F) # remove the row that has a TIME=0 and AMT=0

#----------------------------------------------------------------------------------------------------
# 1 compartment-first order absorption with 1-compartment metabolite model via PKADVAN package
#----------------------------------------------------------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
#Parent BSV
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL-parent
BSVV  <- rnorm(nsub, mean = 0, sd = 0.12)  #BSV on V-parent
BSVKA <- rnorm(nsub, mean=0, sd = 0.14)
BSVF1 <- rnorm(nsub, mean=0, sd = 0.05)

#Metabolite BSV
BSVCLM <- rnorm(nsub, mean = 0, sd = 0.30)  #BSV on CL-metabolite
BSVVM  <- rnorm(nsub, mean = 0, sd = 0.13)  #BSV on VM-metabolite

#Define residual error model
RESEDUAL 	<- rnorm(nsub, mean = 0, sd = 0.15)  #parent
RESEDUALMET <- rnorm(nsub, mean = 0, sd = 0.20) #metabolite

#Set population PK parameters for 3-compartment parent model
CLpop <- 0.5          # clearance
Vpop <-  20         # central volume of distribution
KApop <- 1.36        # first-order absorption rate constant
F1pop <- 1          # bioavailability

#Set population parameters for 1-compartment metabolite model
CLpopMet <- 0.8        # Clearance of the metabolite
VpopMet  <- 7       # Volume of distribution of th metabolite

#Set fraction of parent drug converted into the metabolite
FR = 0.60

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
dfadvan$FR <- FR
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)   #creatinine clearance (CLCR) added as a time-changing covariate on CL
dfadvan$V  <- Vpop*exp(BSVV)
dfadvan$KA <- KApop*exp(BSVKA)
dfadvan$F1 <- F1pop*exp(BSVF1)
dfadvan$CLM <- CLpopMet*exp(BSVCLM)*(dfadvan$CLCR/100) #creatinine clearance (CLCR) added as a covariate on CLM
dfadvan$VM <- VpopMet*exp(BSVVM)

#Apply PKADVAN function to each ID in df
simdf <- ddply(dfadvan, .(ID), OneCompFirstOrderAbsOneCompMetab)
head(simdf)
tail(simdf)

#Add residual unexplained variability (within subject variability)
#For example; Additive error model
simdf$DVP <- simdf$IPREDP+RESEDUAL
simdf$DVM <- simdf$IPREDM+RESEDUALMET

#PLOTTING
# SUBSET mdv==0
simdf <- subset(simdf, MDV==0)

#IPRED
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),geom="ribbon", fill="red",fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#plotobj2 <- plotobj + scale_y_log10("Concentration\n")
#plotobj2


#DV
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),geom="ribbon",fill="red", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#----------------------------------------------------------------------------------------------------
# 2 compartment-first order absorption with 1-compartment metabolite model via PKADVAN package
#----------------------------------------------------------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
#Parent BSV
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL-parent
BSVV2  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V-parent

#Metabolite BSV
BSVCLM <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL-metabolite


#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.2) #parent
PROPMET <- rnorm(nsub, mean = 0, sd = 0.15) #metabolite

#Set population PK parameters for 3-compartment parent model
CLpop <- 0.5          # clearance
V2pop <- 20         # central volume of distribution
Q3pop <- 0.5        # inter-compartmental clearance (1)
V3pop <- 30         # peripheral volume of distribution (1)
KApop <- 0.4        # first-order absorption rate constant
F1pop <- 1          # bioavailability

#Set population parameters for 1-compartment metabolite model
CLpopMet <- 1        # Clearance of the metabolite
VpopMet  <- 15       # Volume of distribution of th metabolite

#Set fraction of parent drug converted into the metabolite
FR = 0.85

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
dfadvan$FR <- FR
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)   #creatinine clearance (CLCR) added as a time-changing covariate on CL
dfadvan$V2 <- V2pop*exp(BSVV2)
dfadvan$Q3 <- Q3pop
dfadvan$V3 <- V3pop
dfadvan$KA <- KApop
dfadvan$F1 <- F1pop
dfadvan$CLM <- CLpopMet**exp(BSVCLM)*(dfadvan$CLCR/100) #creatinine clearance (CLCR) added as a covariate on CLM
dfadvan$VM <- VpopMet

#Apply PKADVAN function to each ID in df
simdf <- ddply(dfadvan, .(ID), TwoCompFirstOrderAbsOneCompMetab)
head(simdf)
tail(simdf)

#Add residual unexplained variability (within subject variability)
#For example; proportional error model
simdf$DVP <- simdf$IPREDP*(1+PROP)
simdf$DVM <- simdf$IPREDM*(1+PROPMET)

#IPRED
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),geom="ribbon", fill="red",fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#DV
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),geom="ribbon",fill="red", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#----------------------------------------------------------------------------------------------------
# 3 compartment-first order absorption with 1-compartment metabolite model via PKADVAN package
#----------------------------------------------------------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
#Parent BSV
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL-parent
BSVV2  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V-parent

#Metabolite BSV
BSVCLM <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL-metabolite


#Define residual error model
PROP 	<- rnorm(nsub, mean = 0, sd = 0.2) #parent
PROPMET <- rnorm(nsub, mean = 0, sd = 0.15) #metabolite

#Set population PK parameters for 3-compartment parent model
CLpop <- 0.5          # clearance
V2pop <- 20         # central volume of distribution
Q3pop <- 0.5        # inter-compartmental clearance (1)
V3pop <- 30         # peripheral volume of distribution (1)
Q4pop <- 0.3        # inter-compartmental clearance (2)
V4pop <- 40         # peripheral volume of distribution (2)
KApop <- 0.4        # first-order absorption rate constant
F1pop <- 1          # bioavailability

#Set population parameters for 1-compartment metabolite model
CLpopMet <- 1        # Clearance of the metabolite
VpopMet  <- 15       # Volume of distribution of th metabolite

#Set fraction of parent drug converted into the metabolite
FR = 0.85

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
dfadvan$FR <- FR
#Calculate group parameter values including any covariate effects
dfadvan$CL <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)   #creatinine clearance (CLCR) added as a time-changing covariate on CL
dfadvan$V2 <- V2pop*exp(BSVV2)
dfadvan$Q3 <- Q3pop
dfadvan$V3 <- V3pop
dfadvan$Q4 <- Q4pop
dfadvan$V4 <- V4pop
dfadvan$KA <- KApop
dfadvan$F1 <- F1pop
dfadvan$CLM <- CLpopMet**exp(BSVCLM)*(dfadvan$CLCR/100) #creatinine clearance (CLCR) added as a covariate on CLM
dfadvan$VM <- VpopMet

#Apply PKADVAN function to each ID in df
simdf <- ddply(dfadvan, .(ID), ThreeCompFirstOrderAbsOneCompMetab)
head(simdf)
tail(simdf)

#Add residual unexplained variability (within subject variability)
#For example; proportional error model
simdf$DVP <- simdf$IPREDP*(1+PROP)
simdf$DVM <- simdf$IPREDM*(1+PROPMET)

#IPRED
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),geom="ribbon", fill="red",fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#DV
plotobj <- NULL
titletext <- expression(atop('Simulated Drug Concentrations',
                             atop(italic("Red line is the median. Gray band is the 90% confidence interval")
                             )))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),geom="ribbon",fill="red", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#===========================================================================================
#                     infusion first-order formation Metabolite models
#===========================================================================================
#Generate data frame: (if you havea a NONMEM df ready then just read and apply function)
#Set dose records:
#number of subjects
nsub <- 2000
ID <- 1:nsub

tlast <- 72
#Set dose records
dosetimessim <- c(0,20.5)

#Make dataframe
df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,tlast,1),dosetimessim))),"AMT"=0,"RATE"=0,"CLCR"=90,"MDV"=0)

df$CLCR[df$ID>=0.5*nsub] <- 70

#subset doserows
doserowssim <- subset(df, TIME%in%dosetimessim)

#Dose & RATE can be arbitrary
doserowssim$AMT  <- 120
doserowssim$RATE <- 4
doserowssim$MDV  <- 1

#Add back dose information
df <- rbind(df,doserowssim)
df <- df[order(df$TIME,-df$AMT),]       # arrange df by TIME with the AMT from high to low
df <- subset(df, (TIME==0 & AMT==0)==F) # remove the row that has a TIME=0 and AMT=0

#---------------------------------------------------------------
# 1 compartment-IV infusion with 1-compartment metabolite model via PKADVAN package
#---------------------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
#Parent BSV
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL-parent
BSVV  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V-parent

#Metabolite BSV
BSVCLM <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL-metabolite
BSVVM  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V-metabolite

#Define residual error model
RESEDUAL 	<- rnorm(nsub, mean = 0, sd = 0.15)  #parent
RESEDUALMET <- rnorm(nsub, mean = 0, sd = 0.05) #metabolite

# Set population PK parameters for parent 3-compartment model
CLpop  <- 1.5          # clearance
Vpop  <- 10         # central volume of distribution

#Set population parameters for 1-compartment metabolite model
CLpopMet <- 2.5        # Clearance of the metabolite
VpopMet  <- 20       # Volume of distribution of th metabolite

#Set fraction of parent drug converted into the metabolite
FR = 0.75

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
dfadvan$FR <- FR
#Calculate group parameter values including any covariate effects
dfadvan$CL  <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)      #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V   <- Vpop*exp(BSVV)
dfadvan$CLM <- CLpopMet*exp(BSVCLM)*(dfadvan$CLCR/100)
dfadvan$VM <- VpopMet*exp(BSVVM)

#Apply the ADVAN function for each ID in df
simdf <- ddply(dfadvan, .(ID), OneCompIVinfusionOneCompMetab)
head(simdf)
tail(simdf)

library(dplyr)
library(scales)
library(tidyr)
df <- simdf %>% gather(DVID, IPRED, IPREDP:IPREDM) %>% mutate(DVID = as.integer(factor(DVID, levels = c('IPREDP', 'IPREDM'))))

#Add residual unexplained variability (within subject variability)
#For example; Additive error model
simdf$DVP <- simdf$IPREDP+RESEDUAL
simdf$DVM <- simdf$IPREDM+RESEDUALMET

#plotting
#Subset MDV==0
simdf <- subset(simdf, MDV==0)

#IPRED
plotobj <- NULL
titletext <- expression(atop('One compartment multiple IV infusion administration',
                             atop(italic("with discrete changes in creatinine clearance covariate values on parent systemic clearance"),
                                  italic("One compartment first-order formation for the metabolite"))))
plotobj <- ggplot(data=df)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj <- plotobj + facet_wrap(~DVID)
plotobj


#IPRED
plotobj <- NULL
titletext <- expression(atop('One compartment multiple IV infusion administration',
                             atop(italic("with discrete changes in creatinine clearance covariate values on parent systemic clearance"),
                                  italic("One compartment first-order formation for the metabolite"))))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),geom="ribbon", fill="red",fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#DV
plotobj <- NULL
titletext <- expression(atop('One compartment multiple IV infusion administration',
                             atop(italic("with discrete changes in creatinine clearance covariate values on parent systemic clearance"),
                                  italic("One compartment first-order formation for the metabolite"))))

plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),geom="ribbon",fill="red", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#---------------------------------------------------------------
# 2 compartment-IV infusion with 1-compartment metabolite model via PKADVAN package
#---------------------------------------------------------------
#Define between subject variability (BSV)
#Use random number generator to simulate residuals from a normal distribution
#Parent BSV
BSVCL <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL-parent
BSVV  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V-parent

#Metabolite BSV
BSVCLM <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on CL-metabolite
BSVVM  <- rnorm(nsub, mean = 0, sd = 0.15)  #BSV on V-metabolite

#Define residual error model
RESEDUAL 	<- rnorm(nsub, mean = 0, sd = 0.2)  #parent
RESEDUALMET <- rnorm(nsub, mean = 0, sd = 0.10) #metabolite

# Set population PK parameters for parent 3-compartment model
CLpop  <- 1.5          # clearance
V1pop  <- 10         # central volume of distribution
Q2pop <- 0.5        # inter-compartmental clearance (1)
V2pop  <- 30         # peripheral volume of distribution (1)

#Set population parameters for 1-compartment metabolite model
CLpopMet <- 2.5        # Clearance of the metabolite
VpopMet  <- 20       # Volume of distribution of th metabolite

#Set fraction of parent drug converted into the metabolite
FR = 0.75

#Modify df for ADVAN calculation and include any covariates on PK parameters
dfadvan <- df
dfadvan$FR <- FR
#Calculate group parameter values including any covariate effects
dfadvan$CL  <- CLpop*exp(BSVCL)*(dfadvan$CLCR/100)      #creatinine clearance (CLCR) added as a covariate on CL
dfadvan$V1  <- V1pop*exp(BSVV)
dfadvan$Q2 <-  Q2pop
dfadvan$V2  <- V2pop
dfadvan$CLM <- CLpopMet*exp(BSVCLM)*(dfadvan$CLCR/100)
dfadvan$VM <- VpopMet*exp(BSVVM)

#Apply the ADVAN function for each ID in df
simdf <- ddply(dfadvan, .(ID), TwoCompIVinfusionOneCompMetab)
head(simdf)
tail(simdf)

#Add residual unexplained variability (within subject variability)
#For example; Additive error model
simdf$DVP <- simdf$IPREDP+RESEDUAL
simdf$DVM <- simdf$IPREDM+RESEDUALMET


#plotting
#Subset MDV==0
simdf <- subset(simdf, MDV==0)

#IPRED
plotobj <- NULL
titletext <- expression(atop('Two compartment multiple IV infusion administration',
                             atop(italic("with discrete changes in creatinine clearance covariate values on parent systemic clearance"),
                                  italic("One compartment first-order formation for the metabolite"))))
plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPREDM),geom="ribbon", fill="red",fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

#DV
plotobj <- NULL
titletext <- expression(atop('Two compartment multiple IV infusion administration',
                             atop(italic("with discrete changes in creatinine clearance covariate values on parent systemic clearance"),
                                  italic("One compartment first-order formation for the metabolite"))))

plotobj <- ggplot(data=simdf)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),fun.y=median, geom="line", colour="red", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVP),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),fun.y=median, geom="line", colour="black", size=1)
plotobj <- plotobj + stat_summary(aes(x=TIME, y= DVM),geom="ribbon",fill="red", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
plotobj <- plotobj + ggtitle(titletext)
plotobj <- plotobj + scale_y_continuous("Concentration\n")
plotobj <- plotobj + scale_x_continuous("\nTime after dose")
plotobj

