#' PKADVAN: Faster Pharmacokinetic Simulations Using the ADVAN-style Analytical Solutions of Common Pharmacokinetic Models.
#'
#' @description The \code{PKADVAN} package presents the ADVAN-style analytical solutions, published by \emph{Abuhelwa et al.} (1),
#' for the 1, 2, and 3 compartments IV bolus, infusion, and first-order absoprtion pharmacokinetic models for
#' simulating pharmacokinetic data. The ADVAN-style analytical solutions functions simulate the time-course
#' of drug amounts in the respective compartments of a pharmacokinetic models and they \emph{ADVANCE} the solution of the model from one time point to the next, allowing for any dose or covariate factors to be accounted for.
#'
#' The \code{PKADVAN} functions presented in this package were written using an integrated R/C++ functions for better
#' performance and significantly faster computational speed compared to using R-only coded functions.
#'
#' A main application for the \code{PKADVAN} package is for simulation from stochastic population pharmacokinetic
#' models coded in R and incorporated into reactive Shiny Web applications. The speed of the integrated R/C++
#' analytical solutions is important here as they allow for simulating larger populations (e.g. 1000 subjects)
#' without compromising the speed and the reactivity benefits of the Shiny application.
#'
#' All ADVAN-style analytical functions were successfully validated against the commericially available population pharmacokinetic modelling software; NONMEM (2).
#'
#' @details \tabular{ll}{
#' Package:\tab \code{PKADVAN}\cr
#' Type:\tab Package\cr
#' Version:\tab 0.1.0\cr
#' Date:\tab 04-04-2016\cr
#' License:\tab GPL (>= 2)\cr
#' }
#'
#' The \code{PKADVAN} functions of a pharmacokinetic model accept a NONMEM-style data frame and calculate drug
#' amounts in the respective compartments.
#'
#' The data frame should have the following columns: \code{ID, TIME, AMT}, in addition to columns for the
#' disposition parameters of the respective pharmacokinetic model (e.g. \code{CL, V, Q}).
#' Please see help pages for more details!
#'
#' @section PKADVAN Functions:
#' The PKADVAN functions are capable of simulating arbitrary dosing regemins and can account for covariate
#' structures; however, covariate effects on respective parameters must be calculated prior porocessing
#' simulations.
#'
#' @seealso The following links provides more details on how to use the \code{PKADVAN} package for processing
#' simulations from various pharmacokinetic models.
#'
#' @seealso Basic Pharmacokinetic Models:
#' @seealso \code{\link{OneCompIVbolus}}, \code{\link{TwoCompIVbolus}}, \code{\link{ThreeCompIVbolus}}
#' @seealso \code{\link{OneCompIVinfusion}}, \code{\link{TwoCompIVinfusion}}, \code{\link{ThreeCompIVinfusion}}
#' @seealso \code{\link{OneCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}, \code{\link{TwoCompFirstOrderAbs}}
#'
#' @seealso Transit Models:
#' @seealso \code{\link{OneCompOneTransit}}, \code{\link{OneCompTwoTransit}}, \code{\link{OneCompThreeTransit}}, \code{\link{OneCompFourTransit}}
#' @seealso \code{\link{TwoCompOneTransit}}, \code{\link{TwoCompTwoTransit}}, \code{\link{TwoCompThreeTransit}}, \code{\link{TwoCompFourTransit}}
#'
#' @seealso Metabolite Models:
#' @seealso \code{\link{ThreeCompFirstOrderAbsOneCompMetab}}
#' @seealso \code{\link{ThreeCompIVinfusionOneCompMetab}}
#'
#'
#' @references (1) Abuhelwa, A.Y., D.J. Foster, and R.N. Upton, \emph{ADVAN-style analytical solutions for common pharmacokinetic models}. Journal of pharmacological and toxicological methods, 2015. \bold{73}: p. 42-48.
#'
#' (2) Beal, S., et al., \emph{NONMEM User's Guides, Part V}. (1989-2009), Icon Development Solutions, Ellicott City, MD, USA, 2009.

#' @name PKADVAN-package
#' @author Ahmad Abuhelwa, David Foster, Richard Upton\cr
#' Australian Center for Pharmacometrics\cr
#' School of Pharmacy and Medical Sciences\cr
#' University of South Australia\cr
#' \email{Ahmad.Abuhelwa@@myamil.unisa.edu.au}

NULL


