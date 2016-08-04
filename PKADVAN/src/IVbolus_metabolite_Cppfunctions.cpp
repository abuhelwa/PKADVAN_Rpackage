#include <Rcpp.h>
#include <math.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

//----------------------------------------------------------------------------------------------
// IV bolus- 1 compartment parent with 1 compartment first-order metabolite formation: Cpp code
//----------------------------------------------------------------------------------------------

// [[Rcpp::export]]

DataFrame OneCompIVbolusOneCompMetabCpp(DataFrame inputFrame){

  //  Create vectors of each element used in function and for constructing output dataframe
  Rcpp::DoubleVector TIME = inputFrame["TIME"];
  Rcpp::DoubleVector AMT = inputFrame["AMT"];
  Rcpp::DoubleVector k10 = inputFrame["k10"];
  Rcpp::DoubleVector kmf = inputFrame["kmf"];
  Rcpp::DoubleVector kme = inputFrame["kme"];
  Rcpp::DoubleVector A1 = inputFrame["A1"];
  Rcpp::DoubleVector AM = inputFrame["AM"];

  double currentk10, currentTIME, currentA1,currentAM, currentkmf, currentkme, E1;
  double previousA1, previousAM;

  // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
  // for counter from 1 to the number of rows in input data frame
  for(int counter = 1; counter < inputFrame.nrows(); counter++){
    // pull out all the variables that will be used for calculation
    currentk10  = k10[ counter ];
	currentkmf  = kmf[ counter ];
	currentkme  = kme[ counter ];
	
    currentTIME = TIME[ counter ] - TIME[ counter - 1];
    previousA1  = A1[ counter - 1 ];
	previousAM  = AM[ counter -1 ];
	
	//calculate E1- rate constant exiting from compartment 1
    E1          = currentk10 + currentkmf;

    // Calculate currentA1: Amount in the central compartment
    currentA1 = previousA1*exp(-currentTIME*E1);
	
	// Calculate currentAM: Amount of the metabolite
	currentAM = previousAM*exp(-currentTIME*currentkme) +
	    currentkmf*previousA1*(exp(-currentTIME*currentkme)/(E1-currentkme) + exp(-currentTIME*E1)/(currentkme-E1));

    // Fill in Amounts and check for other doses
	AM[ counter ] = currentAM;
    A1[ counter ] = currentA1 + AMT[ counter ];
	
  } // end for loop
  return(0);
}


//---------------------------------------------------------------------------------------------
// IV bolus- 2 compartment parent with 1 compartment first-order metabolite formation: Cpp code
//---------------------------------------------------------------------------------------------

// [[Rcpp::export]]
DataFrame TwoCompIVbolusOneCompMetabCpp(DataFrame inputFrame){

  //  Create vectors of each element used in function and for constructing output dataframe
  Rcpp::DoubleVector TIME = inputFrame["TIME"];
  Rcpp::DoubleVector AMT = inputFrame["AMT"];
  Rcpp::DoubleVector k10 = inputFrame["k10"];
  Rcpp::DoubleVector k12 = inputFrame["k12"];
  Rcpp::DoubleVector k21 = inputFrame["k21"];
  Rcpp::DoubleVector kmf = inputFrame["kmf"];
  Rcpp::DoubleVector kme = inputFrame["kme"];
  Rcpp::DoubleVector A1 = inputFrame["A1"];
  Rcpp::DoubleVector A2 = inputFrame["A2"];
  Rcpp::DoubleVector AM = inputFrame["AM"];

  double currentTIME, currentk10, currentk12, currentk21,currentkmf, currentkme, E1, E2, lambda1, lambda2;
  double previousA1, previousA2, previousAM, currentA1, currentA2,currentAM;

  // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
  // for counter from 1 to the number of rows in input data frame
  for(int counter = 1; counter < inputFrame.nrows(); counter++){

    // pull out all the variables that will be used for calculation
    currentk10  = k10[ counter ];
    currentk12  = k12[ counter ];
    currentk21  = k21[ counter ];
	currentkmf  = kmf[ counter ];
	currentkme  = kme[ counter ];
	
    currentTIME = TIME[ counter ] - TIME[ counter - 1];
    previousA1  = A1[ counter - 1 ];
    previousA2  = A2[ counter - 1 ];
	previousAM  = AM[ counter -1 ];


    //calculate hybrid rate constants
    E1          = currentk10 + currentk12 + currentkmf;
    E2          = currentk21 ;

    lambda1 = 0.5*((E1+E2)+sqrt(pow((E1+E2),2)-4*(E1*E2-currentk12*currentk21)));
    lambda2 = 0.5*((E1+E2)-sqrt(pow((E1+E2),2)-4*(E1*E2-currentk12*currentk21)));

    // Calculate currentA1: Amount in the central compartment
    currentA1 = (((previousA1*E2+previousA2*currentk21)-previousA1*lambda1)*exp(-currentTIME*lambda1)-((previousA1*E2+previousA2*currentk21)-previousA1*lambda2)*exp(-currentTIME*lambda2))/(lambda2-lambda1);

    //calculate currentA2: Amount in the peripheral compartment
    currentA2 = (((previousA2*E1+previousA1*currentk12)-previousA2*lambda1)*exp(-currentTIME*lambda1)-((previousA2*E1+previousA1*currentk12)-previousA2*lambda2)*exp(-currentTIME*lambda2))/(lambda2-lambda1);
    
	// Calculate currentAM: Amount of the metabolite
	currentAM = previousAM*exp(-currentTIME*currentkme) ;
	currentAM = currentAM + currentkmf*(previousA2*currentk21*(-exp(-currentTIME*lambda2)/((lambda1-lambda2)*(lambda2-currentkme)) +
                                                  exp(-currentTIME*lambda1)/((lambda1-lambda2)*(lambda1-currentkme)) -
                                                  exp(-currentTIME*currentkme)/((lambda1-currentkme)*(currentkme-lambda2))) + 
                           previousA1*((lambda2-E2)*exp(-currentTIME*lambda2)/((lambda1-lambda2)*(lambda2-currentkme)) +
                                           (currentkme-E2)*exp(-currentTIME*currentkme)/((lambda1-currentkme)*(currentkme-lambda2)) +
                                           (E2-lambda1)*exp(-currentTIME*lambda1)/((lambda1-lambda2)*(lambda1-currentkme))) ); 

	    
    // Fill in Amounts and look for other doses
    A2[ counter ] = currentA2;
    AM[ counter ] = currentAM;
    A1[ counter ] = currentA1 + AMT[ counter ] ;

  } // end for loop

  return 0;
}

//---------------------------------------------------------------------------------------------
// IV bolus- 3 compartment parent with 1 compartment first-order metabolite formation: Cpp code
//---------------------------------------------------------------------------------------------

// [[Rcpp::export]]
DataFrame ThreeCompIVbolusOneCompMetabCpp(DataFrame inputFrame){

  //  Create vectors of each element used in function and for constructing output dataframe
  Rcpp::DoubleVector TIME = inputFrame["TIME"];
  Rcpp::DoubleVector AMT = inputFrame["AMT"];
  Rcpp::DoubleVector k10 = inputFrame["k10"];
  Rcpp::DoubleVector k12 = inputFrame["k12"];
  Rcpp::DoubleVector k21 = inputFrame["k21"];
  Rcpp::DoubleVector k13 = inputFrame["k13"];
  Rcpp::DoubleVector k31 = inputFrame["k31"];
  Rcpp::DoubleVector kmf = inputFrame["kmf"];
  Rcpp::DoubleVector kme = inputFrame["kme"];  
  Rcpp::DoubleVector A1 = inputFrame["A1"];
  Rcpp::DoubleVector A2 = inputFrame["A2"];
  Rcpp::DoubleVector A3 = inputFrame["A3"];
  Rcpp::DoubleVector AM = inputFrame["AM"];
  
  double currentTIME, currentk10, currentk12, currentk21, currentk13, currentk31,currentkmf, currentkme;
  double previousA1, previousA2, previousA3,previousAM, currentA1, currentA2,currentA3, currentAM  ;
  double a, b, c, m, n, Q, alpha, beta, gamma, theta, B, C, I, J ;
  double E1, E2, E3, lambda1, lambda2, lambda3;

  // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
  // for counter from 1 to the number of rows in input data frame
  for(int counter = 1; counter < inputFrame.nrows(); counter++){

    // pull out all the variables that will be used for calculation
    currentk10  = k10[ counter ];
    currentk12  = k12[ counter ];
    currentk21  = k21[ counter ];
    currentk13  = k13[ counter ];
    currentk31  = k31[ counter ];
	currentkmf  = kmf[ counter ];
	currentkme  = kme[ counter ];
	
    currentTIME = TIME[ counter ] - TIME[ counter - 1];
    previousA1  = A1[ counter - 1 ];
    previousA2  = A2[ counter - 1 ];
    previousA3  = A3[ counter - 1 ];
    previousAM  = AM[ counter -1 ];
	
    //calculate hybrid rate constants
    E1          = currentk10 + currentk12 + currentk13 + currentkmf ;
    E2          = currentk21;
    E3          = currentk31;

    a = E1+E2+E3;
    b = E1*E2+E3*(E1+E2)-currentk12*currentk21-currentk13*currentk31;
    c = E1*E2*E3-E3*currentk12*currentk21-E2*currentk13*currentk31;

    m = (3*b - pow(a,2.0))/3;
    n = (2*pow(a,3.0) - 9*a*b + 27*c)/27;
    Q = pow(n,2.0)/4 + pow(m,3.0)/27;

    alpha = sqrt(-1*Q);
    beta  = -1*n/2;
    gamma = sqrt(pow(beta,2.0) + pow(alpha,2.0));
    theta = atan2(alpha,beta);

    lambda1 = a/3 + pow(gamma,1.0/3.0)*(cos(theta/3.0) + sqrt(3)*sin(theta/3.0));
    lambda2 = a/3 + pow(gamma,1.0/3.0)*(cos(theta/3.0) - sqrt(3)*sin(theta/3.0));
    lambda3 = a/3 -(2*pow(gamma,1.0/3.0)*cos(theta/3.0));

    B = previousA2*currentk21+previousA3*currentk31;
    C = E3*previousA2*currentk21+E2*previousA3*currentk31;
    I = previousA1*currentk12*E3-previousA2*currentk13*currentk31+previousA3*currentk12*currentk31;
    J = previousA1*currentk13*E2+previousA2*currentk13*currentk21-previousA3*currentk12*currentk21;

    // Calculate currentA1: Amount in the central compartment
    // Split equation into multiple steps to ensure arithmatic correctness
    currentA1 = previousA1*(exp(-currentTIME*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-currentTIME*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-currentTIME*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
    currentA1 = currentA1 + exp(-currentTIME*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-currentTIME*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-currentTIME*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2)) ;

    //calculate currentA2: Amount in the peripheral compartment1
    currentA2 = previousA2*(exp(-currentTIME*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-currentTIME*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-currentTIME*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
    currentA2 = currentA2 + exp(-currentTIME*lambda1)*(I-previousA1*currentk12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-currentTIME*lambda2)*(previousA1*currentk12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-currentTIME*lambda3)*(previousA1*currentk12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2));

    //calculate currentA3: Amount in the peripheral compartment2
    currentA3 = previousA3*(exp(-currentTIME*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-currentTIME*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-currentTIME*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
    currentA3 = currentA3 + exp(-currentTIME*lambda1)*(J-previousA1*currentk13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-currentTIME*lambda2)*(previousA1*currentk13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-currentTIME*lambda3)*(previousA1*currentk13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2)) ;

    // Calculate currentAM: Amount of the metabolite
    currentAM = previousAM*exp(-currentTIME*currentkme);
    currentAM = currentAM + currentkmf*(previousA1*((E3*lambda1-E2*E3-pow(lambda1,2.0)+E2*lambda1)*exp(-currentTIME*lambda1)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-currentkme)) +
                      (-E3*lambda2+E2*E3+pow(lambda2,2.0)-E2*lambda2)*exp(-currentTIME*lambda2)/((lambda1-lambda2)*(lambda2-lambda3)*(lambda2-currentkme)) +
						(-E3*lambda3+E2*E3+pow(lambda3,2.0)-E2*lambda3)*exp(-currentTIME*lambda3)/((lambda1-lambda3)*(lambda3-lambda2)*(lambda3-currentkme)) +
							(-E3*currentkme+E2*E3+pow(currentkme,2.0)-E2*currentkme)*exp(-currentTIME*currentkme)/((lambda1-currentkme)*(currentkme-lambda2)*(currentkme-lambda3))) + 	
			exp(-currentTIME*lambda1)*(B*lambda1-C)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-currentkme)) + 
				exp(-currentTIME*currentkme)*(C-B*currentkme)/((lambda1-currentkme)*(currentkme-lambda2)*(currentkme-lambda3)) +
					exp(-currentTIME*lambda2)*(C-B*lambda2)/((lambda1-lambda2)*(lambda2-lambda3)*(lambda2-currentkme)) +
						exp(-currentTIME*lambda3)*(C-B*lambda3)/((lambda1-lambda3)*(lambda3-lambda2)*(lambda3-currentkme)));	
    
    // Fill in Amounts and look for other doses
    A2[ counter ] = currentA2;
    A3[ counter ] = currentA3;
	AM[ counter ] = currentAM;
    A1[ counter ] = currentA1 + AMT[ counter ];

  } // end for loop

  return 0;
}

