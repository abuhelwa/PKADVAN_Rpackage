#include <Rcpp.h>
#include <math.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

//--------------------------------------------------------------------------------
// 1 compartment-1transit absorption model via ADVAN-style equations: Cpp Function
//--------------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame OneCompOneTransitCpp(DataFrame inputFrame){
  
  //  Create vectors of each element used in function and for constructing output dataframe
  Rcpp::DoubleVector TIME = inputFrame["TIME"];
  Rcpp::DoubleVector AMT = inputFrame["AMT"];
  Rcpp::DoubleVector KTR = inputFrame["KTR"];
  Rcpp::DoubleVector F1 = inputFrame["F1"];
  Rcpp::DoubleVector k20 = inputFrame["k20"];
  Rcpp::DoubleVector A1 = inputFrame["A1"];		//Amount in the absorption compartment
  Rcpp::DoubleVector A3 = inputFrame["A3"];		//Amount in 1st transit 
  Rcpp::DoubleVector A2 = inputFrame["A2"];		//Amount in central compartment

  double currentTIME, currentKTR, currentk20;
  double previousA1,previousA3,  previousA2;
  double currentA1,currentA3, currentA2;
  
  // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
  // for counter from 1 to the number of rows in input data frame
  for(int counter = 1; counter < inputFrame.nrows(); counter++){
    
    // pull out all the variables that will be used for calculation
    currentk20  = k20[ counter ];
    currentKTR   = KTR[ counter ];
    currentTIME = TIME[ counter ] - TIME[ counter - 1];
    previousA2  = A2[ counter - 1 ];
	previousA3  = A3[ counter - 1 ];	
    previousA1  = A1[ counter - 1 ];

    //Transit compartments
    currentA3 = previousA3*exp(-currentTIME*currentKTR)+currentKTR*previousA1*currentTIME*exp(-currentTIME*currentKTR);
       
   // Calculate currentA2: Amount in the central compartment 
   // Split equation into multiple steps to ensure arithmetic correctness
	currentA2 = previousA2*exp(-currentTIME*currentk20);
	
	currentA2 = currentA2 + currentKTR*(previousA3*(exp(-currentTIME*currentk20)/(currentKTR-currentk20) + exp(-currentTIME*currentKTR)/(currentk20-currentKTR))
	                 + previousA1*currentKTR*(currentTIME*exp(-currentTIME*currentKTR)/(currentk20-currentKTR) + exp(-currentTIME*currentk20)/pow((currentk20-currentKTR),2.0) - exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),2.0)));
      
    // Calculate currentA1: Amount in the absorption compartment
    currentA1 = previousA1 * exp(- currentTIME * currentKTR);
    
    // Fill in Amounts and look for other doses
    A3[ counter ] = currentA3;
    A2[ counter ] = currentA2;
    A1[ counter ] = currentA1 + (AMT[ counter ] * F1[ counter ]);
    
  } // end for loop
  
  return 0;
}

//---------------------------------------------------------------------------
// 1 compartment-2transit absorption model via ADVAN-style equations: Cpp Function
//---------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame OneCompTwoTransitCpp(DataFrame inputFrame){
  
  //  Create vectors of each element used in function and for constructing output dataframe
  Rcpp::DoubleVector TIME = inputFrame["TIME"];
  Rcpp::DoubleVector AMT = inputFrame["AMT"];
  Rcpp::DoubleVector KTR = inputFrame["KTR"];
  Rcpp::DoubleVector F1 = inputFrame["F1"];
  Rcpp::DoubleVector k20 = inputFrame["k20"];
  Rcpp::DoubleVector A1 = inputFrame["A1"];		//Amount in the absorption compartment
  Rcpp::DoubleVector A3 = inputFrame["A3"];		//Amount in 1st transit 
  Rcpp::DoubleVector A4 = inputFrame["A4"];		//Amount in the 2nd transit
  Rcpp::DoubleVector A2 = inputFrame["A2"];		//Amount in central compartment

  
  double currentTIME, currentKTR, currentk20;
  double previousA1,previousA3, previousA4, previousA2;
  double currentA1,currentA3,currentA4, currentA2;
  
  // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
  // for counter from 1 to the number of rows in input data frame
  for(int counter = 1; counter < inputFrame.nrows(); counter++){
    
    // pull out all the variables that will be used for calculation
    currentk20  = k20[ counter ];
    currentKTR   = KTR[ counter ];
    currentTIME = TIME[ counter ] - TIME[ counter - 1];
    previousA2  = A2[ counter - 1 ];
	previousA3  = A3[ counter - 1 ];	
    previousA4  = A4[ counter - 1 ];
    previousA1  = A1[ counter - 1 ];

    //Transit compartments
    currentA3 = previousA3*exp(-currentTIME*currentKTR)+currentKTR*previousA1*currentTIME*exp(-currentTIME*currentKTR);
    currentA4 = previousA4*exp(-currentTIME*currentKTR)+currentKTR*previousA3*currentTIME*exp(-currentTIME*currentKTR)+0.5*pow(currentKTR,2.0)*previousA1*pow(currentTIME,2.0)*exp(-currentTIME*currentKTR);
       
   // Calculate currentA2: Amount in the central compartment 
   // Split equation into multiple steps to ensure arithmetic correctness
	currentA2 = previousA2*exp(-currentTIME*currentk20);
	
	currentA2 = currentA2 + currentKTR*(previousA4*(exp(-currentTIME*currentk20)/(currentKTR-currentk20) + exp(-currentTIME*currentKTR)/(currentk20-currentKTR))
	                 + previousA3*currentKTR*(currentTIME*exp(-currentTIME*currentKTR)/(currentk20-currentKTR) + exp(-currentTIME*currentk20)/pow((currentk20-currentKTR),2.0) - exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),2.0))
	                 + previousA1*pow(currentKTR,2.0)*(pow(currentTIME,2.0)*exp(-currentTIME*currentKTR)/(2*(currentk20-currentKTR)) - currentTIME*exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),2.0) - exp(-currentTIME*currentk20)/pow((currentk20-currentKTR),3.0) + exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),3.0)));
      
    // Calculate currentA1: Amount in the absorption compartment
    currentA1 = previousA1 * exp(- currentTIME * currentKTR);
    
    // Fill in Amounts and look for other doses
    A3[ counter ] = currentA3;
    A4[ counter ] = currentA4;
    A2[ counter ] = currentA2;
    A1[ counter ] = currentA1 + (AMT[ counter ] * F1[ counter ]);
    
  } // end for loop
  
  return 0;
}

//--------------------------------------------------------------------------------
// 1 compartment-3transit absorption model via ADVAN-style equations: Cpp Function
//--------------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame OneCompThreeTransitCpp(DataFrame inputFrame){
  
  //  Create vectors of each element used in function and for constructing output dataframe
  Rcpp::DoubleVector TIME = inputFrame["TIME"];
  Rcpp::DoubleVector AMT = inputFrame["AMT"];
  Rcpp::DoubleVector KTR = inputFrame["KTR"];
  Rcpp::DoubleVector F1 = inputFrame["F1"];
  Rcpp::DoubleVector k20 = inputFrame["k20"];
  Rcpp::DoubleVector A1 = inputFrame["A1"];		//Amount in the absorption compartment
  Rcpp::DoubleVector A3 = inputFrame["A3"];		//Amount in 1st transit 
  Rcpp::DoubleVector A4 = inputFrame["A4"];		//Amount in the 2nd transit
  Rcpp::DoubleVector A5 = inputFrame["A5"];		//Amount in the 3rd transit
  Rcpp::DoubleVector A2 = inputFrame["A2"];		//Amount in central compartment

  
  double currentTIME, currentKTR, currentk20;
  double previousA1,previousA3, previousA4,previousA5, previousA2;
  double currentA1,currentA3,currentA4,currentA5, currentA2;
  
  // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
  // for counter from 1 to the number of rows in input data frame
  for(int counter = 1; counter < inputFrame.nrows(); counter++){
    
    // pull out all the variables that will be used for calculation
    currentk20  = k20[ counter ];
    currentKTR   = KTR[ counter ];
    currentTIME = TIME[ counter ] - TIME[ counter - 1];
    previousA2  = A2[ counter - 1 ];
	previousA3  = A3[ counter - 1 ];	
    previousA4  = A4[ counter - 1 ];
    previousA5  = A5[ counter - 1 ];
    previousA1  = A1[ counter - 1 ];

    //Transit compartments
    currentA3 = previousA3*exp(-currentTIME*currentKTR)+currentKTR*previousA1*currentTIME*exp(-currentTIME*currentKTR);
    currentA4 = previousA4*exp(-currentTIME*currentKTR)+currentKTR*previousA3*currentTIME*exp(-currentTIME*currentKTR)+0.5*pow(currentKTR,2.0)*previousA1*pow(currentTIME,2.0)*exp(-currentTIME*currentKTR);
    currentA5 = previousA5*exp(-currentTIME*currentKTR)+currentKTR*previousA4*currentTIME*exp(-currentTIME*currentKTR)+0.5*pow(currentKTR,2.0)*previousA3*pow(currentTIME,2.0)*exp(-currentTIME*currentKTR)+(1.0/6.0)*pow(currentKTR,3.0)*previousA1*pow(currentTIME,3.0)*exp(-currentTIME*currentKTR);
        
   // Calculate currentA2: Amount in the central compartment 
   // Split equation into multiple steps to ensure arithmetic correctness
	currentA2 = previousA2*exp(-currentTIME*currentk20);
	
	currentA2 = currentA2 + currentKTR*(previousA5*(exp(-currentTIME*currentk20)/(currentKTR-currentk20) + exp(-currentTIME*currentKTR)/(currentk20-currentKTR))
	                 + previousA4*currentKTR*(currentTIME*exp(-currentTIME*currentKTR)/(currentk20-currentKTR) + exp(-currentTIME*currentk20)/pow((currentk20-currentKTR),2.0) - exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),2.0))
	                 + previousA3*pow(currentKTR,2.0)*(pow(currentTIME,2.0)*exp(-currentTIME*currentKTR)/(2*(currentk20-currentKTR)) - currentTIME*exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),2.0) - exp(-currentTIME*currentk20)/pow((currentk20-currentKTR),3.0) + exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),3.0))
	                 + previousA1*pow(currentKTR,3.0)*(pow(currentTIME,3.0)*exp(-currentTIME*currentKTR)/(6*(currentk20-currentKTR)) - pow(currentTIME,2.0)*exp(-currentTIME*currentKTR)/(2*pow((currentk20-currentKTR),2.0)) + currentTIME*exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),3.0)
	                                  + exp(-currentTIME*currentk20)/pow((currentk20-currentKTR),4.0) - exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),4.0) ));
      
    // Calculate currentA1: Amount in the absorption compartment
    currentA1 = previousA1 * exp(- currentTIME * currentKTR);
    
    // Fill in Amounts and look for other doses
    A3[ counter ] = currentA3;
    A4[ counter ] = currentA4;
    A5[ counter ] = currentA5;
    A2[ counter ] = currentA2;
    A1[ counter ] = currentA1 + (AMT[ counter ] * F1[ counter ]);
    
  } // end for loop
  
  return 0;
}

//--------------------------------------------------------------------------------
// 1 compartment-4transit absorption model via ADVAN-style equations: Cpp Function
//--------------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame OneCompFourTransitCpp(DataFrame inputFrame){
  
  //  Create vectors of each element used in function and for constructing output dataframe
  Rcpp::DoubleVector TIME = inputFrame["TIME"];
  Rcpp::DoubleVector AMT = inputFrame["AMT"];
  Rcpp::DoubleVector KTR = inputFrame["KTR"];
  Rcpp::DoubleVector F1 = inputFrame["F1"];
  Rcpp::DoubleVector k20 = inputFrame["k20"];
  Rcpp::DoubleVector A1 = inputFrame["A1"];		//Amount in the absorption compartment
  Rcpp::DoubleVector A3 = inputFrame["A3"];		//Amount in 1st transit 
  Rcpp::DoubleVector A4 = inputFrame["A4"];		//Amount in the 2nd transit
  Rcpp::DoubleVector A5 = inputFrame["A5"];		//Amount in the 3rd transit
  Rcpp::DoubleVector A6 = inputFrame["A6"];		//Amount in the 4th transit
  Rcpp::DoubleVector A2 = inputFrame["A2"];		//Amount in central compartment

  
  double currentTIME, currentKTR, currentk20;
  double previousA1,previousA3, previousA4,previousA5,previousA6, previousA2;
  double currentA1,currentA3,currentA4,currentA5,currentA6, currentA2;
  
  // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
  // for counter from 1 to the number of rows in input data frame
  for(int counter = 1; counter < inputFrame.nrows(); counter++){
    
    // pull out all the variables that will be used for calculation
    currentk20  = k20[ counter ];
    currentKTR   = KTR[ counter ];
    currentTIME = TIME[ counter ] - TIME[ counter - 1];
    previousA2  = A2[ counter - 1 ];
	previousA3  = A3[ counter - 1 ];	
    previousA4  = A4[ counter - 1 ];
    previousA5  = A5[ counter - 1 ];
    previousA6  = A6[ counter - 1 ];
    previousA1  = A1[ counter - 1 ];

    //Transit compartments
    currentA3 = previousA3*exp(-currentTIME*currentKTR)+currentKTR*previousA1*currentTIME*exp(-currentTIME*currentKTR);
    currentA4 = previousA4*exp(-currentTIME*currentKTR)+currentKTR*previousA3*currentTIME*exp(-currentTIME*currentKTR)+0.5*pow(currentKTR,2.0)*previousA1*pow(currentTIME,2.0)*exp(-currentTIME*currentKTR);
    currentA5 = previousA5*exp(-currentTIME*currentKTR)+currentKTR*previousA4*currentTIME*exp(-currentTIME*currentKTR)+0.5*pow(currentKTR,2.0)*previousA3*pow(currentTIME,2.0)*exp(-currentTIME*currentKTR)+(1.0/6.0)*pow(currentKTR,3.0)*previousA1*pow(currentTIME,3.0)*exp(-currentTIME*currentKTR);
    currentA6 = previousA6*exp(-currentTIME*currentKTR)+currentKTR*previousA5*currentTIME*exp(-currentTIME*currentKTR)+0.5*pow(currentKTR,2.0)*previousA4*pow(currentTIME,2.0)*exp(-currentTIME*currentKTR)+(1.0/6.0)*pow(currentKTR,3.0)*previousA3*pow(currentTIME,3.0)*exp(-currentTIME*currentKTR)+(1.0/24.0)*pow(currentKTR,4.0)*previousA1*pow(currentTIME,4.0)*exp(-currentTIME*currentKTR);
       
   // Calculate currentA2: Amount in the central compartment 
   // Split equation into multiple steps to ensure arithmetic correctness
	currentA2 = previousA2*exp(-currentTIME*currentk20);
	
	currentA2 = currentA2 + currentKTR*(previousA6*(exp(-currentTIME*currentk20)/(currentKTR-currentk20) + exp(-currentTIME*currentKTR)/(currentk20-currentKTR))
	                 + previousA5*currentKTR*(currentTIME*exp(-currentTIME*currentKTR)/(currentk20-currentKTR) + exp(-currentTIME*currentk20)/pow((currentk20-currentKTR),2.0) - exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),2.0))
	                 + previousA4*pow(currentKTR,2.0)*(pow(currentTIME,2.0)*exp(-currentTIME*currentKTR)/(2*(currentk20-currentKTR)) - currentTIME*exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),2.0) - exp(-currentTIME*currentk20)/pow((currentk20-currentKTR),3.0) + exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),3.0))
	                 + previousA3*pow(currentKTR,3.0)*(pow(currentTIME,3.0)*exp(-currentTIME*currentKTR)/(6*(currentk20-currentKTR)) - pow(currentTIME,2.0)*exp(-currentTIME*currentKTR)/(2*pow((currentk20-currentKTR),2.0)) + currentTIME*exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),3.0)
	                                  + exp(-currentTIME*currentk20)/pow((currentk20-currentKTR),4.0) - exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),4.0) )
	                 + previousA1*pow(currentKTR,4.0)*(pow(currentTIME,4.0)*exp(-currentTIME*currentKTR)/(24*(currentk20-currentKTR)) - pow(currentTIME,3.0)*exp(-currentTIME*currentKTR)/(6*pow((currentk20-currentKTR),2.0)) + pow(currentTIME,2.0)*exp(-currentTIME*currentKTR)/(2*pow((currentk20-currentKTR),3.0))
	                                  - currentTIME*exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),4.0) - exp(-currentTIME*currentk20)/pow((currentk20-currentKTR),5.0) + exp(-currentTIME*currentKTR)/pow((currentk20-currentKTR),5.0) ));
      
    // Calculate currentA1: Amount in the absorption compartment
    currentA1 = previousA1 * exp(- currentTIME * currentKTR);
    
    // Fill in Amounts and look for other doses
    A3[ counter ] = currentA3;
    A4[ counter ] = currentA4;
    A5[ counter ] = currentA5;
    A6[ counter ] = currentA6;
    A2[ counter ] = currentA2;
    A1[ counter ] = currentA1 + (AMT[ counter ] * F1[ counter ]);
    
  } // end for loop
  
  return 0;
}
