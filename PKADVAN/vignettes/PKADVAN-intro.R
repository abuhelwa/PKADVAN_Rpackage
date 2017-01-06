## ------------------------------------------------------------------------
rm(list=ls(all=TRUE))
	graphics.off()

#load packages
	library(PKADVAN)
	library(plyr)
	library(ggplot2)
	
#Function for 90% confidence intervals calculation
	CI90lo <- function(x) quantile(x, probs=0.05,na.rm=T)
	CI90hi <- function(x) quantile(x, probs=0.95,na.rm=T)

#-----------------------------------------------
# Generate a NONMEM-style simulation data frame
#-----------------------------------------------
	#Set dose records: It can be arbitrary.
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
	
	df$CLCR[df$TIME >= 36] <- 80  #Allow CLCR to change
	
	doserows <- subset(df, TIME%in%dosetimes)
	
	#Dose: It can be any arbitrary dose
	doserows$AMT <- 500
	doserows$MDV <- 1
	
	#Add back dose information
	df <- rbind(df,doserows)
	df <- subset(df, (TIME==0 & AMT==0)==F) # remove the row that has a TIME=0 and AMT=0
	head(df)

## ------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Two compartment-2 transit first-order absorption model via PKADVAN package
#----------------------------------------------------------------------------
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
	CLpop   <- 0.5      #central clearance
	V2pop   <- 20       #central volume of distribution
	Qpop    <- 1        #inter-compartmental clearance
	V3pop   <- 25       #peripheral volume of distribution
	KTRpop  <- 2.05     #first-order absorption rate constant
	F1pop   <- 0.80     #Bioavailability


## ------------------------------------------------------------------------
inputDataFrame <- df
#Calculate group parameter values including any covariate effects
	inputDataFrame$CL  <- CLpop*exp(BSVCL)*(inputDataFrame$CLCR/100)   #creatinine clearance added as a time-changing covariate on CL
	inputDataFrame$V2  <- V2pop*exp(BSVV2)
	inputDataFrame$Q   <- Qpop*exp(BSVQ)
	inputDataFrame$V3  <- V3pop*exp(BSVV3)
	inputDataFrame$KTR <- KTRpop*exp(BSVKTR)
	inputDataFrame$F1  <- F1pop
#This is crucial!
	inputDataFrame <- inputDataFrame[order(inputDataFrame$ID,inputDataFrame$TIME,inputDataFrame$AMT),] # Arrange df by TIME (ascending) and by AMT (descending)
	head(inputDataFrame)

## ------------------------------------------------------------------------
simdf <- ddply(inputDataFrame, .(ID), TwoCompTwoTransit)
head(simdf)

#plot IPRED
	plotobj <- NULL
	titletext <- expression(atop("Simulated Drug Concentrations",
	                             atop(italic("Two compartment- two tranist first-order absorption"),"Median and 90% Prediction Interval, 1000 subjects"))) 
	plotobj <- ggplot(data=simdf)
	plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIME, y= IPRED),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
	plotobj <- plotobj + ggtitle(titletext)
	plotobj <- plotobj + scale_y_continuous("IPREDs\n")
	plotobj <- plotobj + scale_x_continuous("\nTime after dose")
	

## ---- fig.width=7, fig.height=5------------------------------------------
  plotobj

## ------------------------------------------------------------------------
#combined additive and proportional error model
	simdf$DV <- simdf$IPRED*(1 + EPS1) + EPS2

#Plotting	
	#Subset missing data
	simdf <- subset(simdf, MDV==0)
 #DV
	plotobj <- NULL
	titletext <- expression(atop("Simulated Drug Concentrations",
	                             atop(italic("Two compartment- two tranist first-order absorption"),"Median and 90% Prediction Interval, 1000 subjects"))) 
	
	plotobj <- ggplot(data=simdf)
	plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIME, y= DV),geom="ribbon", fun.ymin="CI90lo", fun.ymax="CI90hi", alpha=0.3)
	plotobj <- plotobj + ggtitle(titletext)
	plotobj <- plotobj + scale_y_continuous("DV\n")
	plotobj <- plotobj + scale_x_continuous("\nTime after dose")


## ---- fig.width=7, fig.height=5------------------------------------------
  plotobj

