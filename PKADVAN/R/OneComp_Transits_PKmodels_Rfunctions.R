#' Process simulations using 1-compartment model with 1-transit absorption compartment.
#'
#' @description    \code{OneCompOneTransit} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, F1, KTR, CL, V}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KTR}:\tab is the transit rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V}:\tab is the central volume of distribution.\cr
#' }
#'
#' @usage OneCompOneTransit(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, F1, KTR, CL, V}.
#'
#' @return The function calculates the amounts in the absorption \code{(A1)}, transit \code{(A3)} and central compartment (\code{A2} & individual predicted concentrations, \code{IPRED}) and returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), OneCompOneTransit)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regemins and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior porocessing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{OneCompTwoTransit}}, \code{\link{OneCompThreeTransit}}, \code{\link{OneCompFourTransit}}
#' @seealso \code{\link{OneCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export
#--------------------------------------------------------------------------------------------------------
# First-order absorption: 1 compartment-1transit absorption model via ADVAN-style equations: RCppfunction
#--------------------------------------------------------------------------------------------------------
OneCompOneTransit <- function(inputDataFrame){
#Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV, CL, V2, Q, V3, KTR & F1
#Returns a dataframe with populated columns for A1, A3, A2 and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k20 <- NULL

    #Calculate micro-rate constants
    inputDataFrame$k20 <- inputDataFrame$CL/inputDataFrame$V2

    #Add columns for amounts for marginal speed gain!
    inputDataFrame$A1 <- 0
    inputDataFrame$A3 <- 0
    inputDataFrame$A2 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1]     # Amount in the absorption (GUT) compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0                                     # Amount in the first transit
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                                     # Amount in the central compartment at time zero.

    #Process
    OneCompOneTransitCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A2/inputDataFrame$V2

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k20))

    #Return output
    inputDataFrame
}

#' Process simulations using 1-compartment model with 2-transit absorption compartments.
#'
#' @description    \code{OneCompTwoTransit} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, F1, KTR, CL, V}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KTR}:\tab is the transit rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V}:\tab is the central volume of distribution.\cr
#' }
#'
#' @usage OneCompTwoTransit(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, F1, KTR, CL, V}.
#'
#' @return The function calculates the amounts in the absorption \code{(A1)}, first and second transit \code{(A3, A4)}, respectively, and central compartment (\code{A2} & individual predicted concentrations, \code{IPRED}) and returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), OneCompTwoTransit)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regemins and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior porocessing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{OneCompOneTransit}}, \code{\link{OneCompThreeTransit}}, \code{\link{OneCompFourTransit}}
#' @seealso \code{\link{OneCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export
#--------------------------------------------------------------------------------------------------------
# First-order absorption: 1 compartment-2transit absorption model via ADVAN-style equations: RCppfunction
#--------------------------------------------------------------------------------------------------------
OneCompTwoTransit <- function(inputDataFrame){
#Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV, CL, V2, Q, V3, KTR & F1
#Returns a dataframe with populated columns for A1, A3, A4, A2 and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k20 <- NULL

    #Calculate micro-rate constants
    inputDataFrame$k20 <- inputDataFrame$CL/inputDataFrame$V2

    #Add columns for amounts for marginal speed gain!
    inputDataFrame$A1 <- 0
    inputDataFrame$A3 <- 0
    inputDataFrame$A4 <- 0
    inputDataFrame$A2 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1]     # Amount in the absorption (GUT) compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0                                     # Amount in the first transit
    inputDataFrame$A4[inputDataFrame$TIME==0] <- 0                                     # Amount in the 2nd transit
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                                     # Amount in the central compartment at time zero.

    #Process
    OneCompTwoTransitCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A2/inputDataFrame$V2

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k20))

    #Return output
    inputDataFrame
}

#' Process simulations using 1-compartment model with 3-transit absorption compartments.
#'
#' @description    \code{OneCompThreeTransit} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, F1, KTR, CL, V}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KTR}:\tab is the transit rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V}:\tab is the central volume of distribution.\cr
#' }
#'
#' @usage OneCompThreeTransit(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, F1, KTR, CL, V}.
#'
#' @return The function calculates the amounts in the absorption \code{(A1)}, first, second, and third transit \code{(A3, A4, A5)}, respectively, and central compartment (\code{A2} & individual predicted concentrations, \code{IPRED}) and returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), OneCompThreeTransit)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regemins and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior porocessing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{OneCompOneTransit}}, \code{\link{OneCompTwoTransit}}, \code{\link{OneCompFourTransit}}
#' @seealso \code{\link{OneCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export

#--------------------------------------------------------------------------------------------------------
# First-order absorption: 1 compartment-3transit absorption model via ADVAN-style equations: RCppfunction
#--------------------------------------------------------------------------------------------------------
OneCompThreeTransit <- function(inputDataFrame){
#Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV, CL, V2, Q, V3, KTR & F1
#Returns a dataframe with populated columns for A1, A3, A4, A5, A2 and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k20 <- NULL

    #Calculate micro-rate constants
    inputDataFrame$k20 <- inputDataFrame$CL/inputDataFrame$V2

    #Add columns for amounts for marginal speed gain!
    inputDataFrame$A1 <- 0
    inputDataFrame$A3 <- 0
    inputDataFrame$A4 <- 0
    inputDataFrame$A5 <- 0
    inputDataFrame$A2 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1]     # Amount in the absorption (GUT) compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0                                     # Amount in the 1st transit
    inputDataFrame$A4[inputDataFrame$TIME==0] <- 0                                     # Amount in the 2nd transit
    inputDataFrame$A5[inputDataFrame$TIME==0] <- 0                                     # Amount in the 3rd transit
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                                     # Amount in the central compartment at time zero.

    #Process
    OneCompThreeTransitCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A2/inputDataFrame$V2

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k20))

    #Return output
    inputDataFrame
}

#' Process simulations using 1-compartment model with 4-transit absorption compartments.
#'
#' @description    \code{OneCompFourTransit} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, F1, KTR, CL, V}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KTR}:\tab is the transit rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V}:\tab is the central volume of distribution.\cr
#' }
#'
#' @usage OneCompFourTransit(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, F1, KTR, CL, V}.
#'
#' @return The function calculates the amounts in the absorption \code{(A1)}, first, second, third, and fourth transit \code{(A3, A4, A5, A6)}, respectively, and central compartment (\code{A2} & individual predicted concentrations, \code{IPRED}) and returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), OneCompFourTransit)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regemins and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior porocessing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{OneCompOneTransit}}, \code{\link{OneCompThreeTransit}}, \code{\link{OneCompFourTransit}}
#' @seealso \code{\link{OneCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export

#--------------------------------------------------------------------------------------------------------
# First-order absorption: 1 compartment-4transit absorption model via ADVAN-style equations: RCppfunction
#--------------------------------------------------------------------------------------------------------
OneCompFourTransit <- function(inputDataFrame){
#Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV, CL, V2, Q, V3, KTR & F1
#Returns a dataframe with populated columns for A1, A3, A4, A5, A6, A2 and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k20 <- NULL

    #Calculate micro-rate constants
    inputDataFrame$k20 <- inputDataFrame$CL/inputDataFrame$V2

    #Add columns for amounts for marginal speed gain!
    inputDataFrame$A1 <- 0
    inputDataFrame$A3 <- 0
    inputDataFrame$A4 <- 0
    inputDataFrame$A5 <- 0
    inputDataFrame$A6 <- 0
    inputDataFrame$A2 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1]     # Amount in the absorption (GUT) compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0                                     # Amount in the 1st transit
    inputDataFrame$A4[inputDataFrame$TIME==0] <- 0                                     # Amount in the 2nd transit
    inputDataFrame$A5[inputDataFrame$TIME==0] <- 0                                     # Amount in the 3rd transit
    inputDataFrame$A6[inputDataFrame$TIME==0] <- 0                                     # Amount in the 4th transit
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                                     # Amount in the central compartment at time zero.

    #Process
    OneCompFourTransitCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A2/inputDataFrame$V2

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k20))

    #Return output
    inputDataFrame
}
