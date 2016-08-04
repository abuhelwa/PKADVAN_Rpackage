#include <Rcpp.h>
#include <math.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

//-----------------------------------------------------------------------------------------------------------------
// First-order absorption: 1 compartment parent with 1 compartment first-order metabolite formation: Cpp function
//-----------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame OneCompFirstOrderAbsOneCompMetabCpp(DataFrame inputFrame){

    //Create vectors of each element used in function and for constructing output dataframe
    Rcpp::DoubleVector TIME = inputFrame["TIME"];
    Rcpp::DoubleVector AMT = inputFrame["AMT"];
    Rcpp::DoubleVector KA = inputFrame["KA"];
    Rcpp::DoubleVector F1 = inputFrame["F1"];
    Rcpp::DoubleVector k20 = inputFrame["k20"];
    Rcpp::DoubleVector kmf = inputFrame["kmf"];
    Rcpp::DoubleVector kme = inputFrame["kme"];    
    Rcpp::DoubleVector A1 = inputFrame["A1"];
    Rcpp::DoubleVector A2 = inputFrame["A2"];
    Rcpp::DoubleVector AM = inputFrame["AM"];

    double currentk20, currentKA, currentTIME, previousA1, previousA2, currentA2, currentA1;
    double currentkmf, currentkme, previousAM, currentAM, E2 ;

    // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
    // for counter from 1 to the number of rows in input data frame
    for(int counter = 1; counter < inputFrame.nrows(); counter++){
        //cout << counter << "\n";

        // pull out all the variables that will be used for calculation
        currentk20    = k20[ counter ];
        currentKA     = KA[ counter ];
		currentkmf    = kmf[ counter ];
		currentkme    = kme[ counter ];
	
        currentTIME = TIME[ counter ] - TIME[ counter - 1];
        previousA2    = A2[ counter - 1 ];
		previousAM    = AM[ counter -1 ];
		previousA1    = A1[ counter - 1 ];

		E2 = currentk20 + currentkmf ;
	
        // Calculate currentA2
        currentA2 = (previousA1 * currentKA) / (currentKA - E2) *
            (exp((-currentTIME)*E2) - exp((-currentTIME)*currentKA)) +
            (previousA2 * exp(-currentTIME*E2));
	    
		// Calculate currentAM
        currentAM = previousAM*exp(-currentTIME*currentkme);
		currentAM = currentAM + currentkmf*(previousA1*currentKA*(-exp(-currentTIME*E2)/((currentKA-E2)*(E2-currentkme)) +
				exp(-currentTIME*currentKA)/((currentKA-E2)*(currentKA-currentkme)) -
					exp(-currentTIME*currentkme)/((currentKA-currentkme)*(currentkme-E2))) + 
		previousA2*(exp(-currentTIME*E2)/(currentkme-E2) + exp(-currentTIME*currentkme)/(E2-currentkme) ));
									
        // Calculate currentA1
        currentA1 = previousA1 * exp(-1 * currentTIME * currentKA);

        // Fill in Amounts
        A2[ counter ] = currentA2;
		AM[ counter ] = currentAM;
        A1[ counter ] = currentA1 + (AMT[ counter ] * F1[ counter ]);

    } // end for loop

    return 0;
}

//-----------------------------------------------------------------------------------------------------------------
// First-order absorption: 2 compartment parent with 1 compartment first-order metabolite formation: Cpp function
//-----------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame TwoCompFirstOrderAbsOneCompMetabCpp(DataFrame inputFrame){

    //    Create vectors of each element used in function and for constructing output dataframe
    Rcpp::DoubleVector TIME = inputFrame["TIME"];
    Rcpp::DoubleVector AMT = inputFrame["AMT"];
    Rcpp::DoubleVector KA = inputFrame["KA"];
    Rcpp::DoubleVector F1 = inputFrame["F1"];
    Rcpp::DoubleVector k20 = inputFrame["k20"];
    Rcpp::DoubleVector k23 = inputFrame["k23"];
    Rcpp::DoubleVector k32 = inputFrame["k32"];
    Rcpp::DoubleVector kmf = inputFrame["kmf"];
    Rcpp::DoubleVector kme = inputFrame["kme"];    
    Rcpp::DoubleVector A1 = inputFrame["A1"];
    Rcpp::DoubleVector A2 = inputFrame["A2"];
    Rcpp::DoubleVector A3 = inputFrame["A3"];
    Rcpp::DoubleVector AM = inputFrame["AM"];    

    double currentTIME, currentKA, currentk20, currentk23, currentk32, E2, E3, lambda1, lambda2;
    double previousA1, previousA2, previousA3, currentA1, currentA2,currentA3;
    double currentkmf, currentkme, previousAM, currentAM;

    // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
    // for counter from 1 to the number of rows in input data frame
    for(int counter = 1; counter < inputFrame.nrows(); counter++){

        // pull out all the variables that will be used for calculation
        currentk20    = k20[ counter ];
        currentk23    = k23[ counter ];
        currentk32    = k32[ counter ];
        currentKA     = KA[ counter ];
		currentkmf    = kmf[ counter ];
		currentkme    = kme[ counter ];
	
        currentTIME = TIME[ counter ] - TIME[ counter - 1];
        previousA2    = A2[ counter - 1 ];
        previousA3    = A3[ counter - 1 ];
		previousAM    = AM[ counter -1 ];
        previousA1    = A1[ counter - 1 ];

        //calculate hybrid rate constants
        E2                    = currentk20 + currentk23 + currentkmf;
        E3                    = currentk32 ;

        lambda1 = 0.5*((E2+E3)+sqrt(pow((E2+E3),2)-4*(E2*E3-currentk23*currentk32)));
        lambda2 = 0.5*((E2+E3)-sqrt(pow((E2+E3),2)-4*(E2*E3-currentk23*currentk32)));

        // Calculate currentA2: Amount in the central compartment
        // Split equation into multiple steps to ensure arithmatic correctness
        currentA2 = (((previousA2*E3+previousA3*currentk32)-previousA2*lambda1)*exp(-currentTIME*lambda1)-((previousA2*E3+previousA3*currentk32)-previousA2*lambda2)*exp(-currentTIME*lambda2))/(lambda2-lambda1);
        currentA2 = currentA2 + previousA1*currentKA*(exp(-currentTIME*currentKA)*(E3-currentKA)/((lambda1-currentKA)*(lambda2-currentKA))+exp(-currentTIME*lambda1)*(E3-lambda1)/((lambda2-lambda1)*(currentKA-lambda1))+exp(-currentTIME*lambda2)*(E3-lambda2)/((lambda1-lambda2)*(currentKA-lambda2)))    ;

        //calculate currentA3: Amount in the peripheral compartment
        currentA3 = (((previousA3*E2+previousA2*currentk23)-previousA3*lambda1)*exp(-currentTIME*lambda1)-((previousA3*E2+previousA2*currentk23)-previousA3*lambda2)*exp(-currentTIME*lambda2))/(lambda2-lambda1);
        currentA3 = currentA3 + previousA1*currentKA*currentk23*(exp(-currentTIME*currentKA)/((lambda1-currentKA)*(lambda2-currentKA))+exp(-currentTIME*lambda1)/((lambda2-lambda1)*(currentKA-lambda1))+exp(-currentTIME*lambda2)/((lambda1-lambda2)*(currentKA-lambda2)));

		// Calculate currentAM
        currentAM = previousAM*exp(-currentTIME*currentkme);
		currentAM = currentAM + currentkmf*(previousA1*currentKA*(
		(lambda1-E3)*exp(-currentTIME*lambda1)/((lambda1-lambda2)*(lambda1-currentKA)*(lambda1-currentkme)) + 
        (E3-lambda2)*exp(-currentTIME*lambda2)/((lambda1-lambda2)*(lambda2-currentKA)*(lambda2-currentkme)) +
        (E3-currentKA)*exp(-currentTIME*currentKA)/((lambda1-currentKA)*(currentKA-lambda2)*(currentKA-currentkme)) + 
        (E3-currentkme)*exp(-currentTIME*currentkme)/((lambda1-currentkme)*(currentkme-lambda2)*(currentkme-currentKA)) ) +
		(previousA2*lambda2-(previousA2*E3+previousA3*currentk32))*exp(-currentTIME*lambda2)/((lambda1-lambda2)*(lambda2-currentkme)) +
		(previousA2*currentkme-(previousA2*E3+previousA3*currentk32))*exp(-currentTIME*currentkme)/((lambda1-currentkme)*(currentkme-lambda2)) + 
		((previousA2*E3+previousA3*currentk32)-previousA2*lambda1)*exp(-currentTIME*lambda1)/((lambda1-lambda2)*(lambda1-currentkme))); 
	
        // Calculate currentA1: Amount in the absorption compartment
        currentA1 = previousA1 * exp(- currentTIME * currentKA);

        // Fill in Amounts
        A2[ counter ] = currentA2;
        A3[ counter ] = currentA3;
		AM[ counter ] = currentAM;
        A1[ counter ] = currentA1 + (AMT[ counter ] * F1[ counter ]);

    } // end for loop

    return 0;
}

//-----------------------------------------------------------------------------------------------------------------
// First-order absorption: 3 compartment parent with 1 compartment first-order metabolite formation: Cpp function
//-----------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame ThreeCompFirstOrderAbsOneCompMetabCpp(DataFrame inputFrame){

    //Create vectors of each element used in function and for constructing output dataframe
    Rcpp::DoubleVector TIME = inputFrame["TIME"];
    Rcpp::DoubleVector AMT = inputFrame["AMT"];
    Rcpp::DoubleVector KA = inputFrame["KA"];
    Rcpp::DoubleVector F1 = inputFrame["F1"];
    Rcpp::DoubleVector k20 = inputFrame["k20"];
    Rcpp::DoubleVector k23 = inputFrame["k23"];
    Rcpp::DoubleVector k32 = inputFrame["k32"];
    Rcpp::DoubleVector k24 = inputFrame["k24"];
    Rcpp::DoubleVector k42 = inputFrame["k42"];
    Rcpp::DoubleVector kmf = inputFrame["kmf"];
    Rcpp::DoubleVector kme = inputFrame["kme"];    
    Rcpp::DoubleVector A1 = inputFrame["A1"];
    Rcpp::DoubleVector A2 = inputFrame["A2"];
    Rcpp::DoubleVector A3 = inputFrame["A3"];
    Rcpp::DoubleVector A4 = inputFrame["A4"];
    Rcpp::DoubleVector AM = inputFrame["AM"];
    
    double currentTIME, currentKA, currentk20, currentk23, currentk32, currentk24, currentk42, currentkmf, currentkme;
    double previousA1, previousA2, previousA3, previousA4, previousAM, currentA1, currentA2, currentA3, currentA4, currentAM ;
    double a, b, c, m, n, Q, alpha, beta, gamma, theta, B, C, I, J ;
    double E2, E3, E4, lambda1, lambda2, lambda3;

    // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
    // for counter from 1 to the number of rows in input data frame
    for(int counter = 1; counter < inputFrame.nrows(); counter++){

        // pull out all the variables that will be used for calculation
        currentk20    = k20[ counter ];
        currentk23    = k23[ counter ];
        currentk32    = k32[ counter ];
        currentk24    = k24[ counter ];
        currentk42    = k42[ counter ];
        currentKA     = KA[ counter ];
		currentkmf    = kmf[ counter ];
		currentkme    = kme[ counter ];

        currentTIME = TIME[ counter ] - TIME[ counter - 1];
        previousA2    = A2[ counter - 1 ];
        previousA3    = A3[ counter - 1 ];
        previousA4    = A4[ counter - 1 ];
		previousAM    = AM[ counter -1 ];
        previousA1    = A1[ counter - 1 ];

        //calculate hybrid rate constants
        E2                    = currentk20 + currentk23 + currentk24 + currentkmf;
        E3                    = currentk32;
        E4                    = currentk42;

        a = E2+E3+E4;
        b = E2*E3+E4*(E2+E3)-currentk23*currentk32-currentk24*currentk42;
        c = E2*E3*E4-E4*currentk23*currentk32-E3*currentk24*currentk42;

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

        B = previousA3*currentk32+previousA4*currentk42;
        C = E4*previousA3*currentk32+E3*previousA4*currentk42;
        I = previousA2*currentk23*E4-previousA3*currentk24*currentk42+previousA4*currentk23*currentk42;
        J = previousA2*currentk24*E3+previousA3*currentk24*currentk32-previousA4*currentk23*currentk32;


        // Calculate currentA2: Amount in the central compartment
        // Split equation into multiple steps to ensure arithmetic correctness
        currentA2 = previousA2*(exp(-currentTIME*lambda1)*(E3-lambda1)*(E4-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-currentTIME*lambda2)*(E3-lambda2)*(E4-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-currentTIME*lambda3)*(E3-lambda3)*(E4-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
        currentA2 = currentA2 + exp(-currentTIME*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-currentTIME*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-currentTIME*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2));
        currentA2 = currentA2 + previousA1*currentKA*(exp(-currentTIME*lambda1)*(E3-lambda1)*(E4-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(currentKA-lambda1))+exp(-currentTIME*lambda2)*(E3-lambda2)*(E4-lambda2)/((lambda1-lambda2)*(lambda3-lambda2)*(currentKA-lambda2))+exp(-currentTIME*lambda3)*(E3-lambda3)*(E4-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)*(currentKA-lambda3))+exp(-currentTIME*currentKA)*(E3-currentKA)*(E4-currentKA)/((lambda1-currentKA)*(lambda2-currentKA)*(lambda3-currentKA)));

        //calculate currentA3: Amount in the peripheral compartment1
        currentA3 = previousA3*(exp(-currentTIME*lambda1)*(E2-lambda1)*(E4-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-currentTIME*lambda2)*(E2-lambda2)*(E4-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-currentTIME*lambda3)*(E2-lambda3)*(E4-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
        currentA3 = currentA3 + exp(-currentTIME*lambda1)*(I-previousA2*currentk23*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-currentTIME*lambda2)*(previousA2*currentk23*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-currentTIME*lambda3)*(previousA2*currentk23*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2)) ;
        currentA3 = currentA3 + previousA1*currentKA*currentk23*(exp(-currentTIME*lambda1)*(E4-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(currentKA-lambda1))+exp(-currentTIME*lambda2)*(E4-lambda2)/((lambda1-lambda2)*(lambda3-lambda2)*(currentKA-lambda2))+exp(-currentTIME*lambda3)*(E4-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)*(currentKA-lambda3))+exp(-currentTIME*currentKA)*(E4-currentKA)/((lambda1-currentKA)*(lambda2-currentKA)*(lambda3-currentKA)));

        //calculate currentA4: Amount in the peripheral compartment2
        currentA4 = previousA4*(exp(-currentTIME*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-currentTIME*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-currentTIME*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
        currentA4 = currentA4 + exp(-currentTIME*lambda1)*(J-previousA2*currentk24*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-currentTIME*lambda2)*(previousA2*currentk24*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-currentTIME*lambda3)*(previousA2*currentk24*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2)); 
        currentA4 = currentA4 + previousA1*currentKA*currentk24*(exp(-currentTIME*lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(currentKA-lambda1))+exp(-currentTIME*lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2)*(currentKA-lambda2))+exp(-currentTIME*lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)*(currentKA-lambda3))+exp(-currentTIME*currentKA)*(E3-currentKA)/((lambda1-currentKA)*(lambda2-currentKA)*(lambda3-currentKA)));

		//calculate currentAM: Amount in the metabolite compartment
        currentAM = previousAM*exp(-currentTIME*currentkme)+currentkmf*previousA2*(exp(-currentTIME*lambda1)*(E3-lambda1)*(E4-lambda1)/((currentkme-lambda1)*(lambda2-lambda1)*(lambda3-lambda1))+exp(-currentTIME*lambda2)*(E3-lambda2)*(E4-lambda2)/((currentkme-lambda2)*(lambda1-lambda2)*(lambda3-lambda2))+exp(-currentTIME*lambda3)*(E3-lambda3)*(E4-lambda3)/((currentkme-lambda3)*(lambda1-lambda3)*(lambda2-lambda3))+exp(-currentTIME*currentkme)*(E3-currentkme)*(E4-currentkme)/((lambda1-currentkme)*(lambda2-currentkme)*(lambda3-currentkme)));
        currentAM = currentAM + currentkmf*(exp(-currentTIME*lambda1)*(B*lambda1-C)/((lambda1-lambda2)*(lambda1-lambda3)*(lambda1-currentkme))+exp(-currentTIME*lambda2)*(C-B*lambda2)/((lambda1-lambda2)*(lambda2-lambda3)*(lambda2-currentkme))+exp(-currentTIME*lambda3)*(C-B*lambda3)/((lambda1-lambda3)*(lambda3-lambda2)*(lambda3-currentkme))-exp(-currentTIME*currentkme)*(B*currentkme-C)/((lambda1-currentkme)*(currentkme-lambda2)*(currentkme-lambda3))) ;
        currentAM = currentAM + currentkmf*previousA1*currentKA*(exp(-currentTIME*lambda1)*(E3-lambda1)*(E4-lambda1)/((currentkme-lambda1)*(lambda2-lambda1)*(lambda3-lambda1)*(currentKA-lambda1))+exp(-currentTIME*lambda2)*(E3-lambda2)*(E4-lambda2)/((currentkme-lambda2)*(lambda1-lambda2)*(lambda3-lambda2)*(currentKA-lambda2))+exp(-currentTIME*lambda3)*(E3-lambda3)*(E4-lambda3)/((currentkme-lambda3)*(lambda1-lambda3)*(lambda2-lambda3)*(currentKA-lambda3))+exp(-currentTIME*currentKA)*(E3-currentKA)*(E4-currentKA)/((currentkme-currentKA)*(lambda1-currentKA)*(lambda2-currentKA)*(lambda3-currentKA))+exp(-currentTIME*currentkme)*(E3-currentkme)*(E4-currentkme)/((lambda1-currentkme)*(lambda2-currentkme)*(lambda3-currentkme)*(currentKA-currentkme)));
 
        // Calculate currentA1: Amount in the absorption compartment
        currentA1 = previousA1 * exp(- currentTIME * currentKA);

        // Fill in Amounts
        A2[ counter ] = currentA2;
        A3[ counter ] = currentA3;
        A4[ counter ] = currentA4;
        AM[ counter ] = currentAM;	
        A1[ counter ] = currentA1 + (AMT[ counter ] * F1[ counter ]);

    } // end for loop


    return 0;
}
