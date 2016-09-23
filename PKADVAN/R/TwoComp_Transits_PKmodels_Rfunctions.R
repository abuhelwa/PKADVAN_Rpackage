#' Process simulations using 2-compartment model with 1-transit absorption compartment.
#'
#' @description    \code{TwoCompOneTransit} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, F1, KTR, CL, V2, Q, V3}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KTR}:\tab is the transit rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V2}:\tab is the central volume of distribution\cr
#' \code{Q}:\tab    is the inter-compartmental clearance\cr
#' \code{V3}:\tab is the peripheral volume of distribution\cr
#' }
#'
#' @usage TwoCompOneTransit(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, F1, KTR, CL, V2, Q, V3}.
#'
#' @return The function calculates the amounts in the absorption \code{(A1)}, transit \code{(A4)}, central (\code{A2} & individual predicted concentrations, \code{IPRED}) and peripheral compartment \code{(A3)} and returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), TwoCompOneTransit)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regimens and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior processing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{TwoCompTwoTransit}}, \code{\link{TwoCompThreeTransit}}, \code{\link{TwoCompFourTransit}}
#' @seealso \code{\link{OneCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export

#--------------------------------------------------------------------------------------------------------
# First-order absorption: 2 compartment-1transit absorption model via ADVAN-style equations: RCppfunction
#--------------------------------------------------------------------------------------------------------
TwoCompOneTransit <- function(inputDataFrame){
#Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV, CL, V2, Q, V3, KTR & F1
#Returns a dataframe with populated columns for A1, A4, A2, A3 and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k20 <- k23 <- k32 <- k30 <- NULL

    #Calculate micro-rate constants
    inputDataFrame$k20 <- inputDataFrame$CL/inputDataFrame$V2
    inputDataFrame$k23 <- inputDataFrame$Q/inputDataFrame$V2
    inputDataFrame$k32 <- inputDataFrame$Q/inputDataFrame$V3
    inputDataFrame$k30 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1]     # Amount in the absorption (GUT) compartment at time zero.
    inputDataFrame$A4[inputDataFrame$TIME==0] <- 0                                     # Amount in the first transit
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                                     # Amount in the central compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0                                     # Amount in the peripheral compartment at time zero.

    #Process
    TwoCompOneTransitCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A2/inputDataFrame$V2

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k20,k23,k32,k30))

    #Return output
    inputDataFrame
}

#' Process simulations using 2-compartment model with 2-transit absorption compartment.
#'
#' @description    \code{TwoCompTwoTransit} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, F1, KTR, CL, V2, Q, V3}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KTR}:\tab is the transit rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V2}:\tab is the central volume of distribution\cr
#' \code{Q}:\tab    is the inter-compartmental clearance\cr
#' \code{V3}:\tab is the peripheral volume of distribution\cr
#' }
#'
#' @usage TwoCompTwoTransit(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, F1, KTR, CL, V2, Q, V3}.
#'
#' @return The function calculates the amounts in the absorption \code{(A1)}, first, and second transit \code{(A4, A5)}, respectively, central (\code{A2} & individual predicted concentrations, \code{IPRED}) and peripheral compartment \code{(A3)} and returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), TwoCompTwoTransit)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regimens and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior processing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{TwoCompOneTransit}}, \code{\link{TwoCompThreeTransit}}, \code{\link{TwoCompFourTransit}}
#' @seealso \code{\link{OneCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export

#--------------------------------------------------------------------------------------------------------
# First-order absorption: 2 compartment-2transit absorption model via ADVAN-style equations: RCppfunction
#--------------------------------------------------------------------------------------------------------
TwoCompTwoTransit <- function(inputDataFrame){
#Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV, CL, V2, Q, V3, KTR & F1
#Returns a dataframe with populated columns for A1, A4, A5, A2, A3 and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k20 <- k23 <- k32 <- k30 <- NULL

    #Calculate micro-rate constants
    inputDataFrame$k20 <- inputDataFrame$CL/inputDataFrame$V2
    inputDataFrame$k23 <- inputDataFrame$Q/inputDataFrame$V2
    inputDataFrame$k32 <- inputDataFrame$Q/inputDataFrame$V3
    inputDataFrame$k30 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1]     # Amount in the absorption (GUT) compartment at time zero.
    inputDataFrame$A4[inputDataFrame$TIME==0] <- 0                                     # Amount in the first transit
    inputDataFrame$A5[inputDataFrame$TIME==0] <- 0                                     # Amount in the 2nd transit
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                                     # Amount in the central compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0                                     # Amount in the peripheral compartment at time zero.

    #Process
    TwoCompTwoTransitCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A2/inputDataFrame$V2

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k20,k23,k32,k30))

    #Return output
    inputDataFrame
}

#' Process simulations using 2-compartment model with 3-transit absorption compartment.
#'
#' @description    \code{TwoCompThreeTransit} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, F1, KTR, CL, V2, Q, V3}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KTR}:\tab is the transit rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V2}:\tab is the central volume of distribution\cr
#' \code{Q}:\tab    is the inter-compartmental clearance\cr
#' \code{V3}:\tab is the peripheral volume of distribution\cr
#' }
#'
#' @usage TwoCompThreeTransit(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, F1, KTR, CL, V2, Q, V3}.
#'
#' @return The function calculates the amounts in the absorption \code{(A1)}, first, second, and third transit \code{(A4, A5, A6)}, respectively, central (\code{A2} & individual predicted concentrations, \code{IPRED}) and peripheral compartment \code{(A3)} and returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), TwoCompThreeTransit)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regimens and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior processing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{TwoCompOneTransit}}, \code{\link{TwoCompTwoTransit}}, \code{\link{TwoCompFourTransit}}
#' @seealso \code{\link{OneCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export

#--------------------------------------------------------------------------------------------------------
# First-order absorption: 2 compartment-3transit absorption model via ADVAN-style equations: RCppfunction
#--------------------------------------------------------------------------------------------------------
TwoCompThreeTransit <- function(inputDataFrame){
#Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV, CL, V2, Q, V3, KTR & F1
#Returns a dataframe with populated columns for A1, A4, A5, A6, A2, A3 and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k20 <- k23 <- k32 <- k30 <- NULL

    #Calculate micro-rate constants
    inputDataFrame$k20 <- inputDataFrame$CL/inputDataFrame$V2
    inputDataFrame$k23 <- inputDataFrame$Q/inputDataFrame$V2
    inputDataFrame$k32 <- inputDataFrame$Q/inputDataFrame$V3
    inputDataFrame$k30 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1]     # Amount in the absorption (GUT) compartment at time zero.
    inputDataFrame$A4[inputDataFrame$TIME==0] <- 0                                     # Amount in the 1st transit
    inputDataFrame$A5[inputDataFrame$TIME==0] <- 0                                     # Amount in the 2nd transit
    inputDataFrame$A6[inputDataFrame$TIME==0] <- 0                                     # Amount in the 3rd transit
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                                     # Amount in the central compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0                                     # Amount in the peripheral compartment at time zero.

    #Process
    TwoCompThreeTransitCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A2/inputDataFrame$V2

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k20,k23,k32,k30))

    #Return output
    inputDataFrame
}

#' Process simulations using 2-compartment model with 4-transit absorption compartment.
#'
#' @description    \code{TwoCompFourTransit} function accepts NONMEM-style data frame for one subject and calculates drug amount in the respective compartments. The data frame should have the following columns: \code{ID, TIME, AMT, F1, KTR, CL, V2, Q, V3}.
#'
#' where:
#' \tabular{ll}{
#' \code{ID}: \tab is the subject ID\cr
#' \code{TIME}:\tab is the sampling time points\cr
#' \code{AMT}:\tab is the dose\cr
#' \code{F1}:\tab is the bioavailability\cr
#' \code{KTR}:\tab is the transit rate contstant\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{CL}:\tab is the central compartment clearance\cr
#' \code{V2}:\tab is the central volume of distribution\cr
#' \code{Q}:\tab    is the inter-compartmental clearance\cr
#' \code{V3}:\tab is the peripheral volume of distribution\cr
#' }
#'
#' @usage TwoCompFourTransit(inputDataFrame)
#'
#'
#' @param inputDataFrame which is the NONMEM-style data frame that contains columns for \code{ID, TIME, AMT, F1, KTR, CL, V2, Q, V3}.
#'
#' @return The function calculates the amounts in the absorption \code{(A1)}, first, second, third, and fourth transit \code{(A4, A5, A6, A7)}, respectively, central (\code{A2} & individual predicted concentrations, \code{IPRED}) and peripheral compartment \code{(A3)} and returns the output added to the \code{inputDataFrame}
#' @details
#' To simulate a population (i.e. the \code{inputDataFrame} has more than one subject \code{ID}),    the function
#' has to be applied for each subject \code{ID}. One way of doing that is through using the \code{ddply} functionality
#' in the \pkg{plyr} package in R. The \code{ddply} functionality allows applying the \pkg{PKADVAN} function to each subject
#' \code{ID} and combines the results into a data frame. Please load \pkg{plyr} package in \code{R} before
#' processsing simulations.
#'
#' \code{ddply(inputDataFrame, .(ID), TwoCompFourTransit)}
#'
#' The \pkg{PKADVAN} function is capable of simulating arbitrary dosing regimens and can account for covariate structures; however,
#' covariate effects on respective parameters must be calculated prior processing simulations.
#' See examples below for more details.
#'
#' @seealso \code{\link{TwoCompOneTransit}}, \code{\link{TwoCompTwoTransit}}, \code{\link{TwoCompFourTransit}}
#' @seealso \code{\link{OneCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}
#' @author Ahmad Abuhelwa\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}
#' @importFrom Rcpp evalCpp
#' @export

#--------------------------------------------------------------------------------------------------------
# First-order absorption: 2 compartment-4transit absorption model via ADVAN-style equations: RCppfunction
#--------------------------------------------------------------------------------------------------------
TwoCompFourTransit <- function(inputDataFrame){
#Accepts a NONMEM style data frame for 1 subject with columns for TIME, AMT,MDV, CL, V2, Q, V3, KTR & F1
#Returns a dataframe with populated columns for A1, A4, A5, A6, A7, A2, A3 and IPRED

    #Setting variables to NULL first to avoid notes "no visible binding for global variable [variable name]" upon checking the package
    k20 <- k23 <- k32 <- k30 <- NULL

    #Calculate micro-rate constants
    inputDataFrame$k20 <- inputDataFrame$CL/inputDataFrame$V2
    inputDataFrame$k23 <- inputDataFrame$Q/inputDataFrame$V2
    inputDataFrame$k32 <- inputDataFrame$Q/inputDataFrame$V3
    inputDataFrame$k30 <- 0

    #set initial values in the compartments
    inputDataFrame$A1[inputDataFrame$TIME==0] <- inputDataFrame$AMT[inputDataFrame$TIME==0]*inputDataFrame$F1[1]     # Amount in the absorption (GUT) compartment at time zero.
    inputDataFrame$A4[inputDataFrame$TIME==0] <- 0                                     # Amount in the 1st transit
    inputDataFrame$A5[inputDataFrame$TIME==0] <- 0                                     # Amount in the 2nd transit
    inputDataFrame$A6[inputDataFrame$TIME==0] <- 0                                     # Amount in the 3rd transit
    inputDataFrame$A7[inputDataFrame$TIME==0] <- 0                                     # Amount in the 4th transit
    inputDataFrame$A2[inputDataFrame$TIME==0] <- 0                                     # Amount in the central compartment at time zero.
    inputDataFrame$A3[inputDataFrame$TIME==0] <- 0                                     # Amount in the peripheral compartment at time zero.

    #Process
    TwoCompFourTransitCpp( inputDataFrame )

    #Calculate IPRED for the central compartment
    inputDataFrame$IPRED <- inputDataFrame$A2/inputDataFrame$V2

    #subset extra columns
    inputDataFrame <- subset(inputDataFrame, select=-c(k20,k23,k32,k30))

    #Return output
    inputDataFrame
}
