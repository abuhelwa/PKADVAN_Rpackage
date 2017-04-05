# R script for simulating a population and concentrations using the PKADVAN package
        #(2) intall ggplot2, MASS, MBESS, Rcpp and plyr packages.

## --------------------------------------------------------------------------------
rm(list=ls(all=TRUE))
	graphics.off()

#load packages
	library(PKADVAN)    #Analytical solutions package
	library(plyr)       #Split and rearrange data, ddply function
	library(ggplot2)    #plotting
	library(MASS)       #mvrnorm function
	library(MBESS)      #cor2cov function
#----------------------------------------------------------------------------------
# Use custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 16))
	theme_bw2 <- theme_update(plot.margin = unit(c(1.1, 1.1, 3, 1.1), "lines"),
	                          axis.title.x = element_text(size = 16, vjust = 0),
	                          axis.title.y = element_text(size = 16, vjust = 0, angle = 90),
	                          strip.text.x = element_text(size = 14),
	                          strip.text.y = element_text(size = 14, angle = 90))


# Function for calculating 5th and 95th percentiles for plotting concentrations
	CIlow <- function(x) quantile(x, probs = 0.05)
	CIhi <- function(x) quantile(x, probs = 0.95)
#---------------------------------------------------------------------------------

#================
# Example 1
	# 2-compartment model
	# Dosing: Oral first-order absorption
	# Population parameter variability on CL, V1, KA
	# Variance-covariance matrix for CL, V2, KA
	# Covariate effects:
	# Gender, Smoking & Creatinine Clearance on CL

#----------------------------------------------------------------------------------------------------------
# Step 1: Setup a data frame with individual PK parameters including any covariate effects on PK parameters
#----------------------------------------------------------------------------------------------------------
	# Define individual
	n <- 1000  #Number of "individuals"
	ID <- seq(from = 1, to = n, by = 1)  #Simulation ID
	WT <- 70  #Total body weight, kg
	AGE <- 60  #Age, years
	SECR <- 100  #Serum Creatinine, umol/L
	SEX <- 0  #Gender, Male (0) Female (1)
	SMOK <- 0  #Smoking Status, Not Current (0) Current (1)

	# Define parameter values
	#Thetas
	CLPOP <- 10     #Clearance, L/h
	V2POP <- 50     #Volume of central compartment, L
	QPOP <-  10     #Inter-compartmental clearance, L/h
	V3POP <- 100    #Volume of peripheral compartment, L
	KAPOP <- 0.5    #Absorption rate constant, h^-1
	F1POP <- 1      #Bioavailability

	COV1 <- 0.5   #Effect of smoking status
	COV2 <- 1.15  #Effect of creatinine clearance on clearance

	#Omegas (as SD)
	ETA1SD <- 0.16
	ETA2SD <- 0.16
	ETA3SD <- 0.16

	#Specify a correlation matrix for ETA's
	R12 <- 0.5  #Correlation coefficient for CL-V1
	R13 <- 0.7  #Correlation coefficient for CL-KA
	R23 <- 0.5  #Correlation coefficient for V1-KA

	#Epsilons (as SD)
	EPS1SD <- 0.3  #Proportional residual error
	EPS2SD <- 0  #Additional residual error (none for this model)

	# Calculate ETA values for each subject
	cor.vec <- c(
	    1, R12, R13,
	    R12, 1, R23,
	    R13, R23, 1)
	CORR <- matrix(cor.vec, 3, 3)

	# Specify the between subject variability for CL, V1, V2
	SDVAL <- c(ETA1SD, ETA2SD, ETA3SD)

	# Use this function to turn CORR and SDVAL into a covariance matrix
	OMEGA <- cor2cov(CORR, SDVAL)

	# Now use multivariate rnorm to turn the covariance matrix into ETA values
	ETAmat <- mvrnorm(n = n, mu = c(0, 0, 0), OMEGA)
	ETA1 <- ETAmat[, 1]
	ETA2 <- ETAmat[, 2]
	ETA3 <- ETAmat[, 3]

	# Define covariate effects
	SMOKCOV <- 1
	if(SMOK == 1) SMOKCOV <- SMOKCOV + COV1
	CRCL <- ((140 - AGE)*WT)/(SECR*0.815)  # Male creatinine clearance
	if(SEX == 1) CRCL <- CRCL*0.85  #Female creatinine clearance

	# Define individual parameter values
	    # Note that 2 comp oral is parameterized using parameters: CL, V2, Q, V3, KA, F1 parameters
	    # refer to function documentation for more information
	CL <- CLPOP*exp(ETA1)*((WT/70)^0.75)*SMOKCOV*((CRCL/90)^COV2)
	V2 <- V2POP*exp(ETA2)*(WT/70)
	Q  <- QPOP*(WT/70)^0.75
	V3 <- V3POP*(WT/70)
	KA <- KAPOP*exp(ETA3)
	F1 <- F1POP

	# Collect the individual parameter values in a parameter dataframe
	par.data <- data.frame(
	    ID, CL, V2, Q, V3, KA, F1,  #patient parameters
	    WT, AGE, SECR, SEX, SMOK)  #covariates
	head(par.data, 10)

#----------------------------------------------------------------------
# Step 2: Make a NONMEM-style data frame with dosing records and individual PK parameters
#---------------------------------------------------------------------
	#This includes:
	    #setting up dosing record
	    #setting up a time sequence
        #join with individual PK parameters

	#Set dose records:
	dosetimes <- c(seq(from = 0, to = 72, by = 12)) #Can be arbitrary [e,g: dosetimes <- c(0,6.5,12,48) ]

	#Make a time sequence (hours). This should include dosetimes.
	TIME <- sort(unique(c(seq(from = 0, to = 144, by = 0.25),dosetimes)))

	#generate df with the time sequence
	df <- expand.grid("ID"=ID,"TIME"=TIME,"AMT"=0,"MDV"=0)

	#Set up doses in AMT column
	doserows <- subset(df, TIME%in%dosetimes)   #subset dose rows (records at dosetimes)
	doserows$AMT <- 500                         #Doses can be arbitrary
	#doserows$AMT[doserows$TIME <= 24] <- 750
	#doserows$AMT[doserows$TIME > 24]  <- 500
	doserows$MDV <- 1

	#Add back dose information
	df <- rbind(df,doserows)

# Join par.data with the NONMEM-style data frame
	inputDataFrame <- join(df, par.data,by="ID")

#This is crucial!
	# Arrange "inputDataFrame" df by ID, TIME (ascending) and by AMT (descending)
	inputDataFrame <- inputDataFrame[order(inputDataFrame$ID,inputDataFrame$TIME,inputDataFrame$AMT),]
	# Remove extra row that has a TIME=0 and AMT=0
	inputDataFrame <- subset(inputDataFrame, (TIME==0 & AMT==0)==F)
	head(inputDataFrame, 10)
#------------------------------------------------------------
# step 3: Apply PKADVAN function of the respective PK model
#------------------------------------------------------------
	#PKADVAN functions retruns the inputDataFrame with calculated amounts in each compartment and IPREDs in the central compartment
	sim.data <- ddply(inputDataFrame, .(ID), TwoCompFirstOrderAbs)
	head(sim.data, 10)

#------------------------------------------------------------
# step 4: Add residual unexplained variability to IPREDs
#------------------------------------------------------------
# Use random number generator to simulate residuals from a normal distribution
	#no. of observations =  no. of time points
	nobs <- length(sim.data$TIME)
	EPS1 <- rnorm(nobs, mean = 0, sd = EPS1SD)  #Proportional residual error
	EPS2 <- rnorm(nobs, mean = 0, sd = EPS2SD)  #Additive residual error
	sim.data$DV <- sim.data$IPRED*(1 + EPS1) + EPS2
	head(sim.data,10)
#----------------------------------------------
# Step 5: Draw some plots of the simulated data
#----------------------------------------------
	# Factor covariate values (for example, SEX)
	#  sim.data$SEXf <- as.factor(sim.data$SEX)
	#  levels(sim.data$SEXf) <- c("Male","Female")

#subset data for MDV ==0
	sim.data <- subset(sim.data, MDV==0)

# Generate a plot of the sim.data
	plotobj <- NULL
	titletext <- expression(atop('Simulated Drug Concentrations',
	                             atop(italic("Red line is the median. Red band is the 90% CI of IPREDs, and blue band is 90% CI of Dvs ")
	                             )))
	plotobj <- ggplot(data = sim.data)
	  plotobj <- plotobj + stat_summary(aes(x = TIME, y = DV), fun.ymin = CIlow,
	    fun.ymax = CIhi, geom = "ribbon", fill = "blue", alpha = 0.2)
	plotobj <- plotobj + stat_summary(aes(x = TIME, y = IPRED), fun.ymin = CIlow,
	                                  fun.ymax = CIhi, geom = "ribbon", fill = "red", alpha = 0.3)
	plotobj <- plotobj + stat_summary(aes(x = TIME, y = IPRED),
	                                  fun.y = median, geom = "line", size = 1, colour = "red")
	plotobj <- plotobj + ggtitle(titletext)
	plotobj <- plotobj + scale_y_continuous("Concentration (mg/L) \n")
	plotobj <- plotobj + scale_x_continuous("\nTime (hours)")
	plotobj


################################################################################################
#Example 2:
	#2 compartments
	#Dosing: intravenous infusion

#-----------------------------------------
# Step 1: Define individual PK parameters
#----------------------------------------
	# Define individual parameter values: This builds on the previous example
	    # Note that 2 comp intravenous is parameterized using parameters: CL, V1, Q, V2 parameters
	    # refer to function documentation for more information
	CL <- CLPOP*exp(ETA1)*((WT/70)^0.75)*SMOKCOV*((CRCL/90)^COV2)
	V1 <- V2POP*exp(ETA2)*(WT/70)
	Q  <- QPOP*(WT/70)^0.75
	V2 <- V3POP*(WT/70)

	# Collect the individual parameter values in a parameter dataframe
	par.data <- data.frame(
	    ID, CL, V1, Q, V2, #patient parameters
	    WT, AGE, SECR, SEX, SMOK)  #covariates
	head(par.data)

#----------------------------------------------------------------------
# Step 2: Make a NONMEM-style data frame with dosing records and individual PK parameters
#---------------------------------------------------------------------
#This includes:
	#setting up dosing record
	#setting up a time sequence
	#join with individual PK parameters

	#Set dose records:
	dosetimes <-  c(0, 23.2, 47.1) #c(seq(from = 0, to = 48, by = 24)) #Can be arbitrary [e,g: dosetimes <- c(0,6,12,48) ]

	#Make a time sequence (hours). This should include dosetimes.
	TIME <- sort(unique(c(seq(from = 0, to = 120, by = 0.25),dosetimes)))

	#generate df with the time sequence
	df <- expand.grid("ID"=ID,"TIME"=TIME,"AMT"=0,"MDV"=0,"RATE"=0)

	#Set up dose in AMT column
	doserows <- subset(df, TIME%in%dosetimes)
	doserows$AMT <- 500                         #Doses can be arbitrary
	doserows$RATE <- 40
	doserows$MDV <- 1

	#Add back dose information
	df <- rbind(df,doserows)

# Join par.data with the NONMEM-style data frame
	inputDataFrame <- join(df, par.data, by="ID")

#This is crucial!
	# Arrange "inputDataFrame" df by ID, TIME (ascending) and by AMT (descending)
	inputDataFrame <- inputDataFrame[order(inputDataFrame$ID,inputDataFrame$TIME,inputDataFrame$AMT),]
	# Remove extra row that has a TIME=0 and AMT=0
	inputDataFrame <- subset(inputDataFrame, (TIME==0 & AMT==0)==F)

	head(inputDataFrame, 10)

#------------------------------------------------------------
# step 3: Apply PKADVAN function of the respective PK model
#------------------------------------------------------------
#PKADVAN functions retruns the inputDataFrame with calculated amounts in each compartment and IPREDs in the central compartment
	sim.data <- ddply(inputDataFrame, .(ID), TwoCompIVinfusion)
	head(sim.data)

#------------------------------------------------------------
# step 4: Add residual unexplained variability to IPREDs
#------------------------------------------------------------
# Use random number generator to simulate residuals from a normal distribution
	#no. of observations =  no. of time points
	nobs <- length(sim.data$TIME)
	EPS1 <- rnorm(nobs, mean = 0, sd = EPS1SD)  #Proportional residual error
	EPS2 <- rnorm(nobs, mean = 0, sd = EPS2SD)  #Additive residual error
	sim.data$DV <- sim.data$IPRED*(1 + EPS1) + EPS2

# Generate a plot of the sim.data
	plotobj <- NULL
	titletext <- expression(atop('Simulated Drug Concentrations',
	                             atop(italic("Red line is the median. Red band is the 90% CI of IPREDs, and blue band is 90% CI of Dvs ")
	                             )))
	plotobj <- ggplot(data = sim.data)
	plotobj <- plotobj + stat_summary(aes(x = TIME, y = DV), fun.ymin = CIlow,
	                                  fun.ymax = CIhi, geom = "ribbon", fill = "blue", alpha = 0.2)
	plotobj <- plotobj + stat_summary(aes(x = TIME, y = IPRED), fun.ymin = CIlow,
	                                  fun.ymax = CIhi, geom = "ribbon", fill = "red", alpha = 0.3)
	plotobj <- plotobj + stat_summary(aes(x = TIME, y = IPRED),
	                                  fun.y = median, geom = "line", size = 1, colour = "red")
	plotobj <- plotobj + ggtitle(titletext)
	plotobj <- plotobj + scale_y_continuous("Concentration (mg/L) \n")
	plotobj <- plotobj + scale_x_continuous("\nTime (hours)")
	plotobj


############################################################################
# Example 3:
	#2 compartment
	#IV bolus
	#Note than IV bolus models has the same parametrization to the IV infusion models

#----------------------------------------------------------------------
# Step 2: Make a NONMEM-style data frame with dosing records and individual PK parameters
#---------------------------------------------------------------------
	#This includes:
	#setting up dosing record
	#setting up a time sequence
	#join with individual PK parameters

	#Set dose records:
	dosetimes <-  c(seq(from = 0, to = 48, by = 8)) #Can be arbitrary [e,g: dosetimes <- c(0,6,12,48) ]

	#Make a time sequence (hours). This should include dosetimes.
	TIME <- sort(unique(c(seq(from = 0, to = 96, by = 0.25),dosetimes)))

	#generate df with the time sequence
	df <- expand.grid("ID"=ID,"TIME"=TIME,"AMT"=0,"MDV"=0)

	#Set up dose in AMT column
	doserows <- subset(df, TIME%in%dosetimes)
	doserows$AMT <- 500                         #Doses can be arbitrary
	doserows$MDV <- 1

	#Add back dose information
	df <- rbind(df,doserows)

# Join par.data with the NONMEM-style data frame
	inputDataFrame <- join(df, par.data,by="ID")

#This is crucial!
	# Arrange "inputDataFrame" df by ID, TIME (ascending) and by AMT (descending)
	inputDataFrame <- inputDataFrame[order(inputDataFrame$ID,inputDataFrame$TIME,inputDataFrame$AMT),]
	# Remove extra row that has a TIME=0 and AMT=0
	inputDataFrame <- subset(inputDataFrame, (TIME==0 & AMT==0)==F)
	head(inputDataFrame, 10)


#------------------------------------------------------------
# step 3: Apply PKADVAN function of the respective PK model
#------------------------------------------------------------
#PKADVAN functions retruns the inputDataFrame with calculated amounts in each compartment and IPREDs in the central compartment
	sim.data <- ddply(inputDataFrame, .(ID), TwoCompIVbolus)
	head(sim.data)
#------------------------------------------------------------
# step 4: Add residual unexplained variability to IPREDs
#------------------------------------------------------------
	# Use random number generator to simulate residuals from a normal distribution
	#no. of observations =  no. of time points
	nobs <- length(sim.data$TIME)
	EPS1 <- rnorm(nobs, mean = 0, sd = EPS1SD)  #Proportional residual error
	EPS2 <- rnorm(nobs, mean = 0, sd = EPS2SD)  #Additive residual error
	sim.data$DV <- sim.data$IPRED*(1 + EPS1) + EPS2

	# Generate a plot of the sim.data
	plotobj <- NULL
	titletext <- expression(atop('Simulated Drug Concentrations',
	                             atop(italic("Red line is the median. Red band is the 90% CI of IPREDs, and blue band is 90% CI of Dvs ")
	                             )))
	plotobj <- ggplot(data = sim.data)
	plotobj <- plotobj + stat_summary(aes(x = TIME, y = DV), fun.ymin = CIlow,
	                                  fun.ymax = CIhi, geom = "ribbon", fill = "blue", alpha = 0.2)
	plotobj <- plotobj + stat_summary(aes(x = TIME, y = IPRED), fun.ymin = CIlow,
	                                  fun.ymax = CIhi, geom = "ribbon", fill = "red", alpha = 0.3)
	plotobj <- plotobj + stat_summary(aes(x = TIME, y = IPRED),
	                                  fun.y = median, geom = "line", size = 1, colour = "red")
	plotobj <- plotobj + ggtitle(titletext)
	plotobj <- plotobj + scale_y_continuous("Concentration (mg/L) \n")
	plotobj <- plotobj + scale_x_continuous("\nTime (hours)")
	plotobj

##################################################################################################
# Example 4: Simulation with time changing covariates
    #2 compartment 2-transit first order absorption model
	#Creatinin clearance is a time-changing covariate on central clearance

#------------------------------------------------------
# Create a NONMEM-style data frame with dosing records
#------------------------------------------------------
	#set number of subjects
	n <- 1
	ID <- 1:n

	#Set dose records:
	dosetimes <- c(0,12)  # This can be arbitrary

	#Now define finer sample times for after a dose to capture Cmax
	doseseq <- c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8,9,10)

	#Use the outer product but with addition to expand this doseseq for all dosetimes
	PKtimes <- outer(dosetimes,doseseq, FUN="+")

	#Make dataframe
	df <- expand.grid("ID"=ID,"TIME"=sort(unique(c(seq(0,48,1),PKtimes))),"AMT"=0,"MDV"=0,"CLCR"=NA)

	#Set time-varying creatinine clearance
	df$CLCR[df$ID ==1 & df$TIME <  24 ] <- 100
	df$CLCR[df$ID ==1 & df$TIME >= 24 ] <- 30

	#Set Doserows. It can be any arbitrary dose
	doserows <- subset(df, TIME%in%dosetimes)
	doserows$AMT <- 100
	doserows$MDV <- 1

	#Add back dose information
	df <- rbind(df,doserows)

	testdf <- df[order(df$ID,df$TIME,df$AMT),] # Arrange df by ID, TIME (ascending) and by AMT (descending)
	testdf

#---------------------
# Define PK parameters
#---------------------
#Define between subject variability on PK parameters
	#BSV (Omegas as SD)
	ETA1CL  <- 0 #0.15
	ETA2V2  <- 0 #0.12
	ETA3Q   <- 0 #0.14
	ETA4V3  <- 0 #0.05
	ETA5KTR <- 0 #0.30

#Define residual error model (Epsilons as SD)
	EPS1	<- 0		#Proportional residual error
	EPS2	<- 0		#Additive residual error

#Use random number generator to simulate residuals from a normal distribution
	BSVCL	<- rnorm(n, mean = 0, sd = ETA1CL)	#BSV on CL
	BSVV2	<- rnorm(n, mean = 0, sd = ETA2V2)	#BSV on V2
	BSVQ	<- rnorm(n, mean = 0, sd = ETA3Q)	#BSV on Q
	BSVV3	<- rnorm(n, mean = 0, sd = ETA4V3)	#BSV on V3
	BSVKTR	<- rnorm(n, mean = 0, sd = ETA5KTR) #BSV on KTR

#Set population PK parameters for 2-compartment 2-transit absorption model
	CLpop   <- 0.5      #clearance
	V2pop   <- 20       #central volume of distribution
	Qpop    <- 1        #inter-compartmental clearance
	V3pop   <- 25       #peripheral volume of distribution
	KTRpop  <- 2.05     #first-order absorption rate constant
	F1pop   <- 0.80     #Bioavailability

#Modify df for PKADVAN calculations and include any covariates on PK parameters
	inputDataFrame <- df

#Calculate group parameter values including any covariate effects
	inputDataFrame$CL  <- CLpop*exp(BSVCL)*(inputDataFrame$CLCR/100)   #creatinine clearance (CLCR) added as a time-changing covariate on CL
	inputDataFrame$V2  <- V2pop*exp(BSVV2)
	inputDataFrame$Q   <- Qpop*exp(BSVQ)
	inputDataFrame$V3  <- V3pop*exp(BSVV3)
	inputDataFrame$KTR <- KTRpop*exp(BSVKTR)
	inputDataFrame$F1  <- F1pop

#This is crucial!
	# Arrange "inputDataFrame" df by ID, TIME (ascending) and by AMT (descending)
	inputDataFrame <- inputDataFrame[order(inputDataFrame$ID,inputDataFrame$TIME,inputDataFrame$AMT),]
	# Remove extra row that has a TIME=0 and AMT=0
	inputDataFrame <- subset(inputDataFrame, (TIME==0 & AMT==0)==F)
	head(inputDataFrame, 10)

#Apply PKADVAN functions
	sim.data <- ddply(inputDataFrame, .(ID), TwoCompTwoTransit)  #simulated data with CLCR as time-varying covariate
	head(sim.data)

#simulated data Without accounting for the time-varying covariate
	inputDataFrame$CL  <- CLpop*exp(BSVCL)
	sim.data2 <- ddply(inputDataFrame, .(ID), TwoCompTwoTransit) #simulated data without CLCR covariate effect
	head(sim.data2)

#Add residual unexplained variability (within subject variability)
	#nobs <- length(sim.data$TIME)
	#EPS1 <- rnorm(n, mean = 0, sd = EPS1)	#Proportional residual error
	#EPS2 <- rnorm(n, mean = 0, sd = EPS2)	#Additive residual error
	#sim.data$DV <- sim.data$IPRED*(1 + EPS1) + EPS2

#Plotting
	# Generate a plot of the sim.data
	plotobj <- NULL
	titletext <- expression(atop('Simulated Drug Amounts',
	                             atop(italic("with (red) or without (blue) the effect of time-varying covariate")
	                             )))
	#titletext <- "Simulation with and without time-changing\ncovariates"
	plotobj <- ggplot(data = sim.data)
	plotobj <- plotobj + geom_line(aes(x = TIME, y = A2), size = 1, colour = "red")     #with covariate
	plotobj <- plotobj + geom_line(aes(x = TIME, y = A2), data= sim.data2, size = 1, colour = "blue")  #Without covariate
	plotobj <- plotobj + ggtitle(titletext)
	plotobj <- plotobj + scale_y_continuous("Amount in central\n compartment\n")
	plotobj <- plotobj + scale_x_continuous("\nTime (hours)")
	plotobj


