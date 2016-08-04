#' Process simulations using 1-compartment IV-bolus model.
#'
#' @description    \code{OneCompIVbolus} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, CL, V}.
#'
#' where:
#' \tabular{ll}{
#'    \code{ID}: \tab is the subject ID\cr
#'    \code{TIME}:\tab is the sampling time points\cr
#'    \code{AMT}:\tab is the dose\cr
#'    \code{CL}:\tab is the central compartment clearance\cr
#'    \code{V}:\tab is the central volume of distribution\cr
#'    }
#'
#' @usage OneCompIVbolus(inputDataFrame)
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, CL, and V}.
#' @return The function calculates the drug amount \code{(A1)} and individual predicted concentrations \code{(IPRED)}
#' in the central compartment and returns the output added to the \code{inputDataFrame}.
#'
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), OneCompIVbolus)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regemins and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior porocessing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{TwoCompIVbolus}}
#' @seealso \code{\link{ThreeCompIVbolus}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @useDynLib PKADVAN
#' @importFrom Rcpp evalCpp
#' @export
#-------------------------------------------------------------------
# 1 compartment IV bolus via ADVAN-style equations: RCppfunctions
#-------------------------------------------------------------------
OneCompIVbolus <- function(inputDataFrame){
    #Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT, MDV, CL, V
    #Returns a dataframe with populated columns for A1, and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k10 <- NULL

    #Calculate micro-rate constants
    inputDataFrame$k10 <- inputDataFrame$CL/inputDataFrame$V

    #Add columns for amounts for marginals speed gain!
    inputDataFrame$A1 <- 0

    #Set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]

    #process calculations using the Cpp code
    OneCompIVbolusCpp( inputDataFrame )

        #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A1/inputDataFrame$V

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k10))

    #Return output
    inputDataFrame
}


#' @title Process simulations using 2-compartment IV-bolus model.
#'
#' @description \code{TwoCompIVbolus} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, CL, V1, Q, V2}
#'
#' where:
#' \tabular{ll}{
#'    \code{ID}: \tab is the subject ID\cr
#'    \code{TIME}:\tab is the sampling time points\cr
#'    \code{AMT}:\tab is the dose\cr
#'    \code{CL}:\tab is the central compartment clearance\cr
#'    \code{V1}:\tab is the central volume of distribution\cr
#'    \code{Q}:\tab    is the inter-compartmental clearance\cr
#'    \code{V2}:\tab is the peripheral volume of distribution\cr
#'    }
#'
#' @usage TwoCompIVbolus(inputDataFrame)
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for: \code{ID, TIME, AMT, CL, V1, Q, V2}.\cr
#'
#' @return The function calculates the amounts in the central (\code{A1} & individual predicted concentrations, \code{IPRED}) and peripheral compartment    \code{(A2)} and returns the output added to \code{inputDataFrame}
#'
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), TwoCompIVbolus)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regemins and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior porocessing simulations.
#' See examples below for more details.
#'
#'
#' @seealso \code{\link{OneCompIVbolus}}
#' @seealso \code{\link{ThreeCompIVbolus}}
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
TwoCompIVbolus <- function(inputDataFrame){
    #Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV, CL, V1, Q, V2
    #Returns a dataframe with populated columns for A1, A2, and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k10 <- k12 <- k21 <- k20 <- NULL

    #Calculate micro-rate constants
    inputDataFrame$k10 <- inputDataFrame$CL/inputDataFrame$V1
    inputDataFrame$k12 <- inputDataFrame$Q/inputDataFrame$V1
    inputDataFrame$k21 <- inputDataFrame$Q/inputDataFrame$V2
    inputDataFrame$k20 <- 0

    #Add columns for amounts for marginals speed gain!
    inputDataFrame$A1 <- 0
    inputDataFrame$A2 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]    # drug amount in the central compartment at time zero.
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0    # drug amount in the peripheral compartment at time zero.

    TwoCompIVbolusCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A1/inputDataFrame$V1

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k10,k12,k21,k20))

    #Return output
    inputDataFrame
}


#' @title Process simulations using 3-compartment IV-bolus model.
#' @description \code{ThreeCompIVbolus} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, CL, V1, Q2, V2, Q3, V3}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V1}:\tab is the central volume of distribution\cr
#' \code{Q2}:\tab is the inter-compartmental clearance (1)\cr
#' \code{V2}:\tab is the peripheral volume of distribution (1)\cr
#' \code{Q3}:\tab is the inter-compartmental clearance (2)\cr
#' \code{V3}:\tab is the peripheral volume of distribution (2)\cr
#' }
#'
#' @usage ThreeCompIVbolus(inputDataFrame)
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for: \code{ID, TIME, AMT, CL, V1, Q2, V2, Q3, V3}.\cr
#' @return The function calculates the amounts in the central (\code{A1} & individual predicted concentrations, \code{IPRED}) and the two peripheral compartments (\code{A2}, \code{A3}) and returns the output added to \code{inputDataFrame}
#'
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), ThreeCompIVbolus)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regemins and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior porocessing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{OneCompIVbolus}}
#' @seealso \code{\link{TwoCompIVbolus}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @export

#-------------------------------------------------------------------
# 3 compartment IV bolus via ADVAN-style equations: RCppfunctions
#-------------------------------------------------------------------
ThreeCompIVbolus <- function(inputDataFrame){
    #Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV,CL, V1, Q12, V2, Q13, V3
    #Returns a dataframe with populated columns for A1, A2, A3,and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k10 <- k12 <- k21 <- k20 <- k13 <- k31 <- k30 <- NULL # Setting the variables to NULL first

    #Calculate micro-rate constants
    inputDataFrame$k10 <- inputDataFrame$CL/inputDataFrame$V1
    inputDataFrame$k12 <- inputDataFrame$Q2/inputDataFrame$V1
    inputDataFrame$k21 <- inputDataFrame$Q2/inputDataFrame$V2
    inputDataFrame$k20 <- 0
    inputDataFrame$k13 <- inputDataFrame$Q3/inputDataFrame$V1
    inputDataFrame$k31 <- inputDataFrame$Q3/inputDataFrame$V3
    inputDataFrame$k30 <- 0

    #Adding these into the input data for marginal speed gain
    inputDataFrame$A1 <- 0
    inputDataFrame$A2 <- 0
    inputDataFrame$A3 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]    # drug amount in the central compartment at time zero.
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0     # drug amount in the 1st-peripheral compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0     # drug amount in the 2nd-peripheral compartment at time zero.

    ThreeCompIVbolusCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A1/inputDataFrame$V1

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k10,k12,k21,k20,k13,k31,k30))

    #Return output
    inputDataFrame
}
