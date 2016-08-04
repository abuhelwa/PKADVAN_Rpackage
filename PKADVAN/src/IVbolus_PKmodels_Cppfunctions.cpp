#include <Rcpp.h>
#include <math.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

//-------------------------------------------------------------
// 1 compartment IV bolus via ADVAN-style equations: Cpp code
//-------------------------------------------------------------

// [[Rcpp::export]]

DataFrame OneCompIVbolusCpp(DataFrame inputFrame){

    //    Create vectors of each element used in function and for constructing output dataframe
    Rcpp::DoubleVector TIME = inputFrame["TIME"];
    Rcpp::DoubleVector AMT = inputFrame["AMT"];
    Rcpp::DoubleVector k10 = inputFrame["k10"];
    Rcpp::DoubleVector A1 = inputFrame["A1"];

    double currentk10, currentTime, previousA1, currentA1;

    // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
    // for counter from 1 to the number of rows in input data frame
    for(int counter = 1; counter < inputFrame.nrows(); counter++){
        // pull out all the variables that will be used for calculation
        currentk10    = k10[ counter ];
        currentTime = TIME[ counter ] - TIME[ counter - 1];
        previousA1    = A1[ counter - 1 ];

        // Calculate currentA1
        currentA1 = previousA1*exp(-currentTime*currentk10);

        // Fill in Amounts and check for other doses
        A1[ counter ] = currentA1 + AMT[ counter ];

    } // end for loop
    return(0);
}


//-------------------------------------------------------------
// 2 compartment IV bolus via ADVAN-style equations: Cpp code
//-------------------------------------------------------------

// input Dataframe from R
// [[Rcpp::export]]
DataFrame TwoCompIVbolusCpp(DataFrame inputFrame){

    //    Create vectors of each element used in function and for constructing output dataframe
    Rcpp::DoubleVector TIME = inputFrame["TIME"];
    Rcpp::DoubleVector AMT = inputFrame["AMT"];
    Rcpp::DoubleVector k10 = inputFrame["k10"];
    Rcpp::DoubleVector k12 = inputFrame["k12"];
    Rcpp::DoubleVector k21 = inputFrame["k21"];
    Rcpp::DoubleVector A1 = inputFrame["A1"];
    Rcpp::DoubleVector A2 = inputFrame["A2"];

    double currentTIME, currentk10, currentk12, currentk21, E1, E2, lambda1, lambda2;
    double previousA1, previousA2, currentA1, currentA2;

    // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
    // for counter from 1 to the number of rows in input data frame
    for(int counter = 1; counter < inputFrame.nrows(); counter++){

        // pull out all the variables that will be used for calculation
        currentk10    = k10[ counter ];
        currentk12    = k12[ counter ];
        currentk21    = k21[ counter ];
        currentTIME = TIME[ counter ] - TIME[ counter - 1];
        previousA1    = A1[ counter - 1 ];
        previousA2    = A2[ counter - 1 ];


        //calculate hybrid rate constants
        E1                    = currentk10 + currentk12;
        E2                    = currentk21 ;

        lambda1 = 0.5*((E1+E2)+sqrt(pow((E1+E2),2)-4*(E1*E2-currentk12*currentk21)));
        lambda2 = 0.5*((E1+E2)-sqrt(pow((E1+E2),2)-4*(E1*E2-currentk12*currentk21)));

        // Calculate currentA1: Amount in the central compartment
        currentA1 = (((previousA1*E2+previousA2*currentk21)-previousA1*lambda1)*exp(-currentTIME*lambda1)-((previousA1*E2+previousA2*currentk21)-previousA1*lambda2)*exp(-currentTIME*lambda2))/(lambda2-lambda1);

        //calculate currentA2: Amount in the peripheral compartment
        currentA2 = (((previousA2*E1+previousA1*currentk12)-previousA2*lambda1)*exp(-currentTIME*lambda1)-((previousA2*E1+previousA1*currentk12)-previousA2*lambda2)*exp(-currentTIME*lambda2))/(lambda2-lambda1);

        // Fill in Amounts and look for other doses
        A2[ counter ] = currentA2;
        A1[ counter ] = currentA1 + AMT[ counter ] ;

    } // end for loop

    return 0;
}

//-------------------------------------------------------------
// 3 compartment IV bolus via ADVAN-style equations: Cpp code
//-------------------------------------------------------------

// input Dataframe from R
// [[Rcpp::export]]
DataFrame ThreeCompIVbolusCpp(DataFrame inputFrame){

    //    Create vectors of each element used in function and for constructing output dataframe
    Rcpp::DoubleVector TIME = inputFrame["TIME"];
    Rcpp::DoubleVector AMT = inputFrame["AMT"];
    Rcpp::DoubleVector k10 = inputFrame["k10"];
    Rcpp::DoubleVector k12 = inputFrame["k12"];
    Rcpp::DoubleVector k21 = inputFrame["k21"];
    Rcpp::DoubleVector k13 = inputFrame["k13"];
    Rcpp::DoubleVector k31 = inputFrame["k31"];
    Rcpp::DoubleVector A1 = inputFrame["A1"];
    Rcpp::DoubleVector A2 = inputFrame["A2"];
    Rcpp::DoubleVector A3 = inputFrame["A3"];

    double currentTIME, currentk10, currentk12, currentk21, currentk13, currentk31, previousA1, previousA2, previousA3, currentA1, currentA2,currentA3 ;
    double a, b, c, m, n, Q, alpha, beta, gamma, theta, B, C, I, J ;
    double E1, E2, E3, lambda1, lambda2, lambda3;

    // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
    // for counter from 1 to the number of rows in input data frame
    for(int counter = 1; counter < inputFrame.nrows(); counter++){

        // pull out all the variables that will be used for calculation
        currentk10    = k10[ counter ];
        currentk12    = k12[ counter ];
        currentk21    = k21[ counter ];
        currentk13    = k13[ counter ];
        currentk31    = k31[ counter ];

        currentTIME = TIME[ counter ] - TIME[ counter - 1];
        previousA1    = A1[ counter - 1 ];
        previousA2    = A2[ counter - 1 ];
        previousA3    = A3[ counter - 1 ];

        //calculate hybrid rate constants
        E1                    = currentk10 + currentk12 + currentk13;
        E2                    = currentk21;
        E3                    = currentk31;

        a = E1+E2+E3;
        b = E1*E2+E3*(E1+E2)-currentk12*currentk21-currentk13*currentk31;
        c = E1*E2*E3-E3*currentk12*currentk21-E2*currentk13*currentk31;

        m = (3*b - pow(a,2.0))/3;
        n = (2*pow(a,3.0) - 9*a*b + 27*c)/27;
        Q = pow(n,2.0)/4 + pow(m,3.0)/27;

        alpha = sqrt(-1*Q);
        beta    = -1*n/2;
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

        // Fill in Amounts and look for other doses
        A2[ counter ] = currentA2;
        A3[ counter ] = currentA3;
        A1[ counter ] = currentA1 + AMT[ counter ] ;

    } // end for loop

    return 0;
}

