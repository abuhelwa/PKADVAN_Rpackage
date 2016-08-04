#' Process simulations using 1-compartment first-order absorption model.
#'
#' @description    \code{OneCompFirstOrderAbs} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, F1, KA, CL, V}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KA}:\tab is the absorption rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V}:\tab is the central volume of distribution.\cr
#' }
#'
#' @usage OneCompFirstOrderAbs(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, F1, KA, CL, V}.
#'
#' @return The function calculates the amounts in the absorption \code{(A1)} and central compartment (\code{A2} & individual predicted concentrations, \code{IPRED}) and returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), OneCompFirstOrderAbs)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regemins and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior porocessing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{TwoCompFirstOrderAbs}}, \code{\link{ThreeCompFirstOrderAbs}}
#' @seealso \code{\link{OneCompOneTransit}}, \code{\link{OneCompTwoTransit}}, \code{\link{OneCompThreeTransit}}, \code{\link{OneCompFourTransit}}
#' @seealso \code{\link{TwoCompOneTransit}}, \code{\link{TwoCompTwoTransit}}, \code{\link{TwoCompThreeTransit}}, \code{\link{TwoCompFourTransit}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export
#------------------------------------------------------------------------------
# 1 compartment-first order absorption via ADVAN-style equations: RCppfunctions
#------------------------------------------------------------------------------
OneCompFirstOrderAbs <- function(inputDataFrame){
    #Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV, CL, V, KA & F1
    #Returns a dataframe with populated columns for A1, A2 and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k10 <- NULL

    #Calculate micro-rate constants
    inputDataFrame$k10    <- inputDataFrame$CL/inputDataFrame$V

    #Add columns for amounts for marginals speed gain!
    inputDataFrame$A1 <- 0
    inputDataFrame$A2 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1] #drug amount in the absorption compartment at time zero.
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                                                #drug amount in the central compartment at time zero.

    OneCompOralCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A2/inputDataFrame$V

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k10))

    #Return output
    inputDataFrame
}


#' Process simulations using 2-compartment first-order absorption model.
#'
#' @description    \code{TwoCompFirstOrderAbs} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, F1, KA, CL, V2, Q, V3}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KA}:\tab is the absorption rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V2}:\tab is the central volume of distribution\cr
#' \code{Q}:\tab    is the inter-compartmental clearance\cr
#' \code{V3}:\tab is the peripheral volume of distribution\cr
#' }
#'
#' @usage TwoCompFirstOrderAbs(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, F1, KA, CL, V2, Q, V3}.
#'
#' @return The function calculates the amounts in the absorption \code{(A1)}, central (\code{A2} & individual predicted concentrations, \code{IPRED}) and peripheral compartment \code{(A3)} and returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), TwoCompFirstOrderAbs)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regemins and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior porocessing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{OneCompFirstOrderAbs}}, \code{\link{ThreeCompFirstOrderAbs}}
#' @seealso \code{\link{OneCompOneTransit}}, \code{\link{OneCompTwoTransit}}, \code{\link{OneCompThreeTransit}}, \code{\link{OneCompFourTransit}}
#' @seealso \code{\link{TwoCompOneTransit}}, \code{\link{TwoCompTwoTransit}}, \code{\link{TwoCompThreeTransit}}, \code{\link{TwoCompFourTransit}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export
#------------------------------------------------------------------------------
# 2 compartment-first order absorption via ADVAN-style equations: RCppfunctions
#------------------------------------------------------------------------------
TwoCompFirstOrderAbs <- function(inputDataFrame){
    #Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV, CL, V2, Q, V3, KA & F1
    #Returns a dataframe with populated columns for A1, A2, A3 and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k20 <- k23 <- k32 <- k30 <- NULL

    #Calculate micro-rate constants
    inputDataFrame$k20 <- inputDataFrame$CL/inputDataFrame$V2
    inputDataFrame$k23 <- inputDataFrame$Q/inputDataFrame$V2
    inputDataFrame$k32 <- inputDataFrame$Q/inputDataFrame$V3
    inputDataFrame$k30 <- 0

    #Add columns for amounts for marginals speed gain!
    inputDataFrame$A1 <- 0
    inputDataFrame$A2 <- 0
    inputDataFrame$A3 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1] #drug amount in the absorption compartment at time zero.
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                                                #drug amount in the central compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0

    #Process
    TwoCompOralCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A2/inputDataFrame$V2

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k20,k23,k32,k30))

    #Return output
    inputDataFrame
}

#' Process simulations using 3-compartment first-order absorption model.
#'
#' @description    \code{ThreeCompFirstOrderAbs} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, F1, KA, CL, V2, Q3, V3, Q4, V4}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KA}:\tab is the absorption rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V2}:\tab is the central volume of distribution\cr
#' \code{Q3}:\tab    is the inter-compartmental clearance (1)\cr
#' \code{V3}:\tab is the peripheral volume of distribution (1)\cr
#' \code{Q4}:\tab is the inter-compartmental clearance (2)\cr
#' \code{V4}:\tab is the peripheral volume of distribution (2)\cr
#' }
#'
#' @usage ThreeCompFirstOrderAbs(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, F1, KA, CL, V2, Q3, V3, Q4, V4}.
#'
#' @return The function calculates the amounts in the absorption \code{(A1)}, central (\code{A2} & individual predicted concentrations, \code{IPRED}) and two peripheral compartments \code{(A3, A4)} and returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), ThreeCompFirstOrderAbs)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regemins and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior porocessing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{OneCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}
#' @seealso \code{\link{OneCompOneTransit}}, \code{\link{OneCompTwoTransit}}, \code{\link{OneCompThreeTransit}}, \code{\link{OneCompFourTransit}}
#' @seealso \code{\link{TwoCompOneTransit}}, \code{\link{TwoCompTwoTransit}}, \code{\link{TwoCompThreeTransit}}, \code{\link{TwoCompFourTransit}}
#' @seealso \code{\link{ThreeCompFirstOrderAbsOneCompMetab}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export

#------------------------------------------------------------------------------
# 3 compartment-first order absorption via ADVAN-style equations: RCppfunctions
#------------------------------------------------------------------------------
ThreeCompFirstOrderAbs <- function(inputDataFrame){
    #Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV, CL, V2, Q, V3, KA & F1
    #Returns a dataframe with populated columns for A1, A2, A3 and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k20 <- k23 <- k32 <- k30 <- k24 <- k42 <- k40 <- NULL

    #Calculate micro-rate constants
    inputDataFrame$k20 <- inputDataFrame$CL/inputDataFrame$V2
    inputDataFrame$k23 <- inputDataFrame$Q3/inputDataFrame$V2
    inputDataFrame$k32 <- inputDataFrame$Q3/inputDataFrame$V3
    inputDataFrame$k30 <- 0
    inputDataFrame$k24 <- inputDataFrame$Q4/inputDataFrame$V2
    inputDataFrame$k42 <- inputDataFrame$Q4/inputDataFrame$V4
    inputDataFrame$k40 <- 0

    #Add columns for amounts for marginals speed gain!
    inputDataFrame$A1 <- 0
    inputDataFrame$A2 <- 0
    inputDataFrame$A3 <- 0
    inputDataFrame$A4 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1]        # Amount in the absorption compartment at time zero.
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                                     # Amount in the central compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0                                     # Amount in the 1st peripheral compartment at time zero.
    inputDataFrame$A4[inputDataFrame$TIME==0] <- 0                                     # Amount in the 2nd peripheral compartment at time zero.

    ThreeCompOralCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A2/inputDataFrame$V2

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k20,k23,k32,k30,k24,k42,k40))

    #Return output
    inputDataFrame
}
