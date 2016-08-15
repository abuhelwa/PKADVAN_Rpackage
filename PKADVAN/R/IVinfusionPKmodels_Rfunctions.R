#' Process simulations using 1-compartment IV-infusion model.
#'
#' @description    \code{OneCompIVinfusion} function accepts NONMEM-style data frame for one subject and calculates drug amount in the
#' respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, RATE, CL, V}.
#'
#' where:
#' \tabular{ll}{
#'    \code{ID}: \tab is the subject ID\cr
#'    \code{TIME}:\tab is the sampling time points\cr
#'    \code{AMT}:\tab is the dose\cr
#'    \code{RATE}:\tab is the infusion rate\cr
#'    \code{CL}:\tab is the central compartment clearance\cr
#'    \code{V}:\tab is the central volume of distribution\cr
#'    }
#'
#' @usage OneCompIVinfusion(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, RATE, CL, and V}.
#' @return The function calculates the drug amount \code{(A1)} and individual predicted concentrations \code{(IPRED)} in the central compartment and returns the output added to the \code{inputDataFrame}.
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), OneCompIVinfusion)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regimens and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior processing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{TwoCompIVinfusion}}
#' @seealso \code{\link{ThreeCompIVinfusion}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export
#-------------------------------------------------------------------
# 1 compartment IV bolus via ADVAN-style equations: RCppfunctions
#-------------------------------------------------------------------
OneCompIVinfusion <- function(inputDataFrame, A1init = 0){
#Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT, MDV, RATE, RATEALL, CL, V
#Returns a dataframe with populated columns for A1, and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k10 <- RATEALL <- TIME <- NULL

    #Sampling Times
    sampletimes <- inputDataFrame$TIME

    #Process infusion doses: This function will add end infusion time points, if they are not already there.
    inputDataFrame <- ProcessInfusionDoses(inputDataFrame)

    #Calculate micro-rate constants
    inputDataFrame$k10    <- inputDataFrame$CL/inputDataFrame$V

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- A1init    # drug amount in the central compartment at time zero.
    OneCompIVinfusionCpp( inputDataFrame )

    #Remove end infusion time points
    inputDataFrame <- subset(inputDataFrame, (TIME%in%sampletimes))

    #Calculate IPRED
    inputDataFrame$IPRED <- inputDataFrame$A1/inputDataFrame$V

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k10,RATEALL))

    #Return output
    inputDataFrame
}

#' @title Process simulations using 2-compartment IV-infusion model.
#'
#' @description \code{TwoCompIVinfusion} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, RATE, CL, V1, Q, V2}
#'
#' where:
#' \tabular{ll}{
#'    \code{ID}: \tab is the subject ID\cr
#'    \code{TIME}:\tab is the sampling time points\cr
#'    \code{AMT}:\tab is the dose\cr
#'    \code{RATE}:\tab is the infusion rate\cr
#'    \code{CL}:\tab is the central compartment clearance\cr
#'    \code{V1}:\tab is the central volume of distribution\cr
#'    \code{Q}:\tab    is the inter-compartmental clearance\cr
#'    \code{V2}:\tab is the peripheral volume of distribution\cr
#'    }
#'
#' @usage TwoCompIVinfusion(inputDataFrame)
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for: \code{ID, TIME, AMT, RATE, CL, V1, Q, V2}.\cr
#' @return The function calculates the amounts in the central (\code{A1} & individual predicted concentrations, \code{IPRED}) and peripheral compartment    \code{(A2)} and returns the output added to \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), TwoCompIVinfusion)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regimens and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior processing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{OneCompIVinfusion}}
#' @seealso \code{\link{ThreeCompIVinfusion}}
#'
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @export
#-------------------------------------------------------------------
# 2 compartment IV bolus via ADVAN-style equations: RCppfunctions
#-------------------------------------------------------------------
# the TwoCompIVinfusion function but hybrid of R and C++ code
TwoCompIVinfusion <- function(inputDataFrame, A1init = 0){
#Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV, RATE, RATEALL, CL, V1, Q, V2
#Returns a dataframe with populated columns for A1, A2, and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k10 <- k12 <- k21 <- k20    <- NULL
    RATEALL <- TIME <- NULL

    #Sampling Times
    sampletimes <- inputDataFrame$TIME

    #Process infusion doses
    inputDataFrame <- ProcessInfusionDoses(inputDataFrame)

    #Calculate micro-rate constants
    inputDataFrame$k10 <- inputDataFrame$CL/inputDataFrame$V1
    inputDataFrame$k12 <- inputDataFrame$Q/inputDataFrame$V1
    inputDataFrame$k21 <- inputDataFrame$Q/inputDataFrame$V2
    inputDataFrame$k20 <- 0

    #Adding a vector for A2 and A2 for marginal speed gain
    inputDataFrame$A1 <- 0
    inputDataFrame$A2 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- A1init	# drug amount in the central compartment at time zero.
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0		# drug amount in the peripheral compartment at time zero.

    TwoCompIVinfusionCpp( inputDataFrame )

    #Remove end infusion time points
    inputDataFrame <- subset(inputDataFrame, (TIME%in%sampletimes))

    #Calculate IPRED
    inputDataFrame$IPRED <- inputDataFrame$A1/inputDataFrame$V1

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k10,k12,k21,k20,RATEALL))

    #Return output
    inputDataFrame
}

#' @title Process simulations using 3-compartment IV-infusion model.
#' @description \code{ThreeCompIVinfusion} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, RATE, CL, V1, Q2, V2, Q3, V3}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{RATE}:\tab is the infusion rate\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V1}:\tab is the central volume of distribution\cr
#' \code{Q2}:\tab is the inter-compartmental clearance (1)\cr
#' \code{V2}:\tab is the peripheral volume of distribution (1)\cr
#' \code{Q3}:\tab is the inter-compartmental clearance (2)\cr
#' \code{V3}:\tab is the peripheral volume of distribution (2)\cr
#' }
#'
#' @usage ThreeCompIVinfusion(inputDataFrame)
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for: \code{ID, TIME, AMT, RATE, CL, V1, Q2, V2, Q3, V3}.\cr
#' @return The function calculates the amounts in the central (\code{A1} & individual predicted concentrations, \code{IPRED}) and the two peripheral compartments (\code{A2}, \code{A3}) and returns the output added to \code{inputDataFrame}
#'
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), ThreeCompIVinfusion)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regimens and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior processing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{OneCompIVinfusion}}
#' @seealso \code{\link{TwoCompIVinfusion}}
#' @seealso \code{\link{ThreeCompIVinfusionOneCompMetab}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @export

#-------------------------------------------------------------------
# 3 compartment IV bolus via ADVAN-style equations: RCppfunctions
#-------------------------------------------------------------------
# the ThreeCompIVinfusion function but hybrid of R and C++ code
ThreeCompIVinfusion <- function(inputDataFrame, A1init = 0){
#Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV,RATE, RATEALL, CL, V1, Q2, V2, Q3, V3
#Returns a dataframe with populated columns for A1, A2, A3,and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k10 <- k12 <- k21 <- k20 <- k13 <- k31 <- k30 <- NULL
    RATEALL <- TIME <- NULL

    #Sampling Times
    sampletimes <- inputDataFrame$TIME

    #Process infusion doses
    inputDataFrame <- ProcessInfusionDoses(inputDataFrame)

    #Calculate rate constants
    inputDataFrame$k10 <- inputDataFrame$CL/inputDataFrame$V1
    inputDataFrame$k12 <- inputDataFrame$Q2/inputDataFrame$V1
    inputDataFrame$k21 <- inputDataFrame$Q2/inputDataFrame$V2
    inputDataFrame$k20 <- 0
    inputDataFrame$k13 <- inputDataFrame$Q3/inputDataFrame$V1
    inputDataFrame$k31 <- inputDataFrame$Q3/inputDataFrame$V3
    inputDataFrame$k30 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- A1init     # drug amount in the central compartment at time zero.
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0     # drug amount in the 1st-peripheral compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0     # drug amount in the 2nd-peripheral compartment at time zero.

    ThreeCompIVinfusionCpp( inputDataFrame )

    #Remove end infusion time points
    inputDataFrame <- subset(inputDataFrame, (TIME%in%sampletimes))

    #Calculate IPRED
    inputDataFrame$IPRED <- inputDataFrame$A1/inputDataFrame$V1

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k10,k12,k21,k20,k13,k31,k30,RATEALL))

    #Return output
    inputDataFrame
}

#---------------------------------------------
# This function is to process infusion doses
#---------------------------------------------
ProcessInfusionDoses <- function (inputDataFrame) {

	#Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
	AMT <- DNUM <- RATEALLI <- DNUMI <- NULL

	#Calculate all amounts
	doserows <- subset(inputDataFrame, AMT!=0)
	dosecount <- nrow(doserows)    #total number of doses
	doserows$DNUM <- 1:dosecount

	#Need to add times for ending the infusions - these may not be in the database
	doserowslast <- doserows
	doserowslast$TIME <- doserowslast$TIME+doserowslast$AMT/doserowslast$RATE
	doserowslast$DNUM <- doserowslast$DNUM*(-1)

	goodcols <- c("ID","TIME","AMT","RATE","DNUM")
	badcols <- which(names(doserowslast)%in%goodcols==F)
	doserowslast[,badcols] <- NA

	#Are there any doserows without a DV value?    These need to precede the infusion change
	noDVindex <- which(doserowslast$TIME%in%inputDataFrame$TIME==F)
	doserowslastnoDV <- doserowslast[noDVindex,]
	doserowslastnoDV$AMT <- 0
	doserowslastnoDV$RATE <- 0
	doserowslastnoDV$DNUM <- NA

	#Collect the new rows
	doserows <- rbind(doserows,doserowslast,doserowslastnoDV)

	#Rewrite previous dose rows with new dose rows
	inputDataFrame$DNUM <- NA
	inputDataFrame <- rbind(inputDataFrame[inputDataFrame$AMT==0,],doserows)
	inputDataFrame <- inputDataFrame[order(inputDataFrame$ID,inputDataFrame$TIME,inputDataFrame$AMT),]

	#Set an extra last row
	lastrow <- tail(inputDataFrame, 1)
	lastrow$TIME <- lastrow$TIME+1
	inputDataFrame <- rbind(inputDataFrame,lastrow)

	#Now fill in the gaps for the covariates by locf
	for (i in badcols)
	{
     inputDataFrame[,i] <- locf(inputDataFrame[,i])
	}

#-------------------------------------------------------------------------------------
	#Process infusion doses in a loop
	inputDataFrame$RATEALL <- 0
	for (DCOUNT in 1:dosecount)
	{
		inputDataFrame$RATEALLI <- 0
		inputDataFrame$DNUMI <- inputDataFrame$DNUM
		inputDataFrame$DNUMI[abs(inputDataFrame$DNUM)!=DCOUNT] <- NA
		inputDataFrame$DNUMI <- locf(inputDataFrame$DNUMI)
		inputDataFrame$DNUMI[is.na(inputDataFrame$DNUMI)==T] <- 0
		inputDataFrame$RATEALLI[inputDataFrame$DNUMI==DCOUNT] <- inputDataFrame$RATE[which(inputDataFrame$DNUM==DCOUNT)]
		inputDataFrame$RATEALL <- inputDataFrame$RATEALL+inputDataFrame$RATEALLI
	}

	#This is crucial
	inputDataFrame$RATEALL[inputDataFrame$DNUM > 0] <- 0

	#Get rid of extra dose rows
	inputDataFrame <- subset(inputDataFrame, (DNUM > 0 | is.na(DNUM)==T))
	inputDataFrame <- subset(inputDataFrame, select = -c(DNUM,RATEALLI,DNUMI))     # no longer needed
}


#-------------------------------------------------
# Last observation carried forward function (locf)
#-------------------------------------------------

 locf <- function (x)
    # #Last observation carried forward
    # #Finds an NA and carries forward the previous value
 {
     good <- !is.na(x)
     positions <- seq(length(x))
     good.positions <- good * positions
     last.good.position <- cummax(good.positions)
     last.good.position[last.good.position == 0] <- NA
     x[last.good.position]
 }
