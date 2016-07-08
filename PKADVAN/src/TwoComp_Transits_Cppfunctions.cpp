#include <Rcpp.h>
#include <math.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

//--------------------------------------------------------------------------------
// 2 compartment-1transit absorption model via ADVAN-style equations: Cpp Function
//--------------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame TwoCompOneTransitCpp(DataFrame inputFrame){
  
  //  Create vectors of each element used in function and for constructing output dataframe
  Rcpp::DoubleVector TIME = inputFrame["TIME"];
  Rcpp::DoubleVector AMT = inputFrame["AMT"];
  Rcpp::DoubleVector KTR = inputFrame["KTR"];
  Rcpp::DoubleVector F1 = inputFrame["F1"];
  Rcpp::DoubleVector k20 = inputFrame["k20"];
  Rcpp::DoubleVector k23 = inputFrame["k23"];
  Rcpp::DoubleVector k32 = inputFrame["k32"];
  Rcpp::DoubleVector k30 = inputFrame["k30"];
  Rcpp::DoubleVector A1 = inputFrame["A1"];   //Amount in the absorption compartment
  Rcpp::DoubleVector A4 = inputFrame["A4"];   //Amount in the 1st transit
  Rcpp::DoubleVector A2 = inputFrame["A2"];   //Amount in central compartment
  Rcpp::DoubleVector A3 = inputFrame["A3"];   //Amount in peripheral compartment
  
  double currentTIME, currentKTR, currentk20, currentk23, currentk32,currentk30, E2, E3, lambda1, lambda2;
  double previousA1,previousA4, previousA2, previousA3;
  double currentA1,currentA4,currentA2,currentA3;
  
  // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
  // for counter from 1 to the number of rows in input data frame
  for(int counter = 1; counter < inputFrame.nrows(); counter++){
    
    // pull out all the variables that will be used for calculation
    currentk20  = k20[ counter ];
    currentk23  = k23[ counter ];
    currentk32  = k32[ counter ];
	currentk30  = k30[ counter ];
    currentKTR   = KTR[ counter ];
    currentTIME = TIME[ counter ] - TIME[ counter - 1];
    previousA2  = A2[ counter - 1 ];
    previousA4  = A4[ counter - 1 ];
    previousA3  = A3[ counter - 1 ];
    previousA1  = A1[ counter - 1 ];
    	
   //calculate hybrid rate constants
    E2          = currentk20 + currentk23;
    E3          = currentk32 + currentk30 ;
    
    lambda1 = 0.5*((E2+E3)+sqrt(pow((E2+E3),2)-4*(E2*E3-currentk23*currentk32)));
    lambda2 = 0.5*((E2+E3)-sqrt(pow((E2+E3),2)-4*(E2*E3-currentk23*currentk32)));
    
    //Transit compartments
    currentA4 = previousA4*exp(-currentTIME*currentKTR)+currentKTR*previousA1*currentTIME*exp(-currentTIME*currentKTR);
        
   // Calculate currentA2: Amount in the central compartment 
   // Split equation into multiple steps to ensure arithmatic correctness
    currentA2 = (exp(-currentTIME*lambda1)*((previousA2*E3+previousA3*currentk32)-previousA2*lambda1)-exp(-currentTIME*lambda2)*((previousA2*E3+previousA3*currentk32)-previousA2*lambda2))/(lambda2-lambda1);
    
    currentA2 = currentA2 + currentKTR*E3*(previousA4*(exp(-currentTIME*currentKTR)/((lambda1-currentKTR)*(lambda2-currentKTR))+exp(-currentTIME*lambda1)/((currentKTR-lambda1)*(lambda2-lambda1))+exp(-currentTIME*lambda2)/((currentKTR-lambda2)*(lambda1-lambda2)))
    +previousA1*currentKTR*(exp(-currentTIME*currentKTR)*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                 -exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 2.0))
                 +exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 2.0))
                 -exp(-currentTIME*currentKTR)*currentTIME/((lambda1-currentKTR)*(currentKTR-lambda2))));
    
    currentA2 = currentA2 + currentKTR*(previousA4*(exp(-currentTIME*currentKTR)*currentKTR/((lambda1-currentKTR)*(currentKTR-lambda2))+exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*(lambda2-currentKTR))-exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*(lambda1-currentKTR)))
    +previousA1*currentKTR*(exp(-currentTIME*currentKTR)*(lambda1*lambda2-pow(currentKTR,2.0))/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))+exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*pow((lambda1-currentKTR), 2.0))-exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*pow((lambda2-currentKTR), 2.0))+exp(-currentTIME*currentKTR)*currentTIME*currentKTR/((lambda1-currentKTR)*(currentKTR-lambda2))));
        
	//calculate currentA3: Amount in the peripheral compartment
  	currentA3 = (exp(-currentTIME*lambda1)*((previousA3*E2+currentk23*previousA2)-previousA3*lambda1)-exp(-currentTIME*lambda2)*((previousA3*E2+currentk23*previousA2)-previousA3*lambda2))/(lambda2-lambda1);
    
    currentA3 = currentA3 + currentKTR*currentk23*(previousA4*(exp(-currentTIME*currentKTR)/((lambda1-currentKTR)*(lambda2-currentKTR))+exp(-currentTIME*lambda1)/((currentKTR-lambda1)*(lambda2-lambda1))+exp(-currentTIME*lambda2)/((currentKTR-lambda2)*(lambda1-lambda2)))
                       +previousA1*currentKTR*(exp(-currentTIME*currentKTR)*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))-exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 2.0))+exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 2.0))-exp(-currentTIME*currentKTR)*currentTIME/((lambda1-currentKTR)*(currentKTR-lambda2))));
        
    // Calculate currentA1: Amount in the absorption compartment
    currentA1 = previousA1 * exp(- currentTIME * currentKTR);
    
    // Fill in Amounts and look for other doses
    A4[ counter ] = currentA4;
    A2[ counter ] = currentA2;
    A3[ counter ] = currentA3;
    A1[ counter ] = currentA1 + (AMT[ counter ] * F1[ counter ]);
    
  } // end for loop
  
  return 0;
}


//---------------------------------------------------------------------------
// 2 compartment-2transit absorption model via ADVAN-style equations: Cpp Function
//---------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame TwoCompTwoTransitCpp(DataFrame inputFrame){
  
  //  Create vectors of each element used in function and for constructing output dataframe
  Rcpp::DoubleVector TIME = inputFrame["TIME"];
  Rcpp::DoubleVector AMT = inputFrame["AMT"];
  Rcpp::DoubleVector KTR = inputFrame["KTR"];
  Rcpp::DoubleVector F1 = inputFrame["F1"];
  Rcpp::DoubleVector k20 = inputFrame["k20"];
  Rcpp::DoubleVector k23 = inputFrame["k23"];
  Rcpp::DoubleVector k32 = inputFrame["k32"];
  Rcpp::DoubleVector k30 = inputFrame["k30"];
  Rcpp::DoubleVector A1 = inputFrame["A1"];   //Amount in the absorption compartment
  Rcpp::DoubleVector A4 = inputFrame["A4"];   //Amount in the 1st transit
  Rcpp::DoubleVector A5 = inputFrame["A5"];   //Amount in the 2nd transit
  Rcpp::DoubleVector A2 = inputFrame["A2"];   //Amount in central compartment
  Rcpp::DoubleVector A3 = inputFrame["A3"];   //Amount in peripheral compartment
  
  double currentTIME, currentKTR, currentk20, currentk23, currentk32,currentk30, E2, E3, lambda1, lambda2;
  double previousA1,previousA4,previousA5, previousA2, previousA3;
  double currentA1,currentA4,currentA5, currentA2,currentA3;
  
  // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
  // for counter from 1 to the number of rows in input data frame
  for(int counter = 1; counter < inputFrame.nrows(); counter++){
    
    // pull out all the variables that will be used for calculation
    currentk20  = k20[ counter ];
    currentk23  = k23[ counter ];
    currentk32  = k32[ counter ];
	currentk30  = k30[ counter ];
    currentKTR   = KTR[ counter ];
    currentTIME = TIME[ counter ] - TIME[ counter - 1];
    previousA2  = A2[ counter - 1 ];
    previousA4  = A4[ counter - 1 ];
    previousA5  = A5[ counter - 1 ];
    previousA3  = A3[ counter - 1 ];
    previousA1  = A1[ counter - 1 ];
    	
   //calculate hybrid rate constants
    E2          = currentk20 + currentk23;
    E3          = currentk32 + currentk30 ;
    
    lambda1 = 0.5*((E2+E3)+sqrt(pow((E2+E3),2)-4*(E2*E3-currentk23*currentk32)));
    lambda2 = 0.5*((E2+E3)-sqrt(pow((E2+E3),2)-4*(E2*E3-currentk23*currentk32)));
    
    //Transit compartments
    currentA4 = previousA4*exp(-currentTIME*currentKTR)+currentKTR*previousA1*currentTIME*exp(-currentTIME*currentKTR);
    currentA5 = previousA5*exp(-currentTIME*currentKTR)+currentKTR*previousA4*currentTIME*exp(-currentTIME*currentKTR)+0.5*pow(currentKTR,2.0)*previousA1*pow(currentTIME,2.0)*exp(-currentTIME*currentKTR);
        
   // Calculate currentA2: Amount in the central compartment 
   // Split equation into multiple steps to ensure arithmatic correctness
    currentA2 = (exp(-currentTIME*lambda1)*((previousA2*E3+previousA3*currentk32)-previousA2*lambda1)-exp(-currentTIME*lambda2)*((previousA2*E3+previousA3*currentk32)-previousA2*lambda2))/(lambda2-lambda1);
    
    currentA2 = currentA2 + currentKTR*E3*(previousA5*(exp(-currentTIME*currentKTR)/((lambda1-currentKTR)*(lambda2-currentKTR))+exp(-currentTIME*lambda1)/((currentKTR-lambda1)*(lambda2-lambda1))+exp(-currentTIME*lambda2)/((currentKTR-lambda2)*(lambda1-lambda2)))
    +previousA4*currentKTR*(exp(-currentTIME*currentKTR)*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                 -exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 2.0))
                 +exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 2.0))
                 -exp(-currentTIME*currentKTR)*currentTIME/((lambda1-currentKTR)*(currentKTR-lambda2)))
    +previousA1*pow(currentKTR,2.0)*((exp(-currentTIME*currentKTR)*(-1*pow(lambda1,2.0)-lambda1*lambda2+3*lambda1*currentKTR-pow(lambda2, 2.0)+3*lambda2*currentKTR-3*pow(currentKTR,2.0)))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                    -exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)/(2*(lambda1-currentKTR)*(currentKTR-lambda2))+exp(-currentTIME*currentKTR)*currentTIME*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                    +exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 3.0))-exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 3.0))));
    
    currentA2 = currentA2 + currentKTR*(previousA5*(exp(-currentTIME*currentKTR)*currentKTR/((lambda1-currentKTR)*(currentKTR-lambda2))+exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*(lambda2-currentKTR))-exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*(lambda1-currentKTR)))
    +previousA4*currentKTR*(exp(-currentTIME*currentKTR)*(lambda1*lambda2-pow(currentKTR,2.0))/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))+exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*pow((lambda1-currentKTR), 2.0))-exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*pow((lambda2-currentKTR), 2.0))+exp(-currentTIME*currentKTR)*currentTIME*currentKTR/((lambda1-currentKTR)*(currentKTR-lambda2)))
    +previousA1*pow(currentKTR,2.0)*(exp(-currentTIME*currentKTR)*(pow(lambda1,2.0)*lambda2+lambda1*pow(lambda2, 2.0)-3*lambda1*lambda2*currentKTR+pow(currentKTR,3.0))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                    +exp(-currentTIME*currentKTR)*currentTIME*(lambda1*lambda2-pow(currentKTR,2.0))/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))+exp(-currentTIME*currentKTR)*currentKTR*pow(currentTIME, 2.0)/(2*(lambda1-currentKTR)*(currentKTR-lambda2))
                    -exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*pow((lambda1-currentKTR), 3.0))+exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*pow((lambda2-currentKTR), 3.0))));
        
	//calculate currentA3: Amount in the peripheral compartment
  	currentA3 = (exp(-currentTIME*lambda1)*((previousA3*E2+currentk23*previousA2)-previousA3*lambda1)-exp(-currentTIME*lambda2)*((previousA3*E2+currentk23*previousA2)-previousA3*lambda2))/(lambda2-lambda1);
    
    currentA3 = currentA3 + currentKTR*currentk23*(previousA5*(exp(-currentTIME*currentKTR)/((lambda1-currentKTR)*(lambda2-currentKTR))+exp(-currentTIME*lambda1)/((currentKTR-lambda1)*(lambda2-lambda1))+exp(-currentTIME*lambda2)/((currentKTR-lambda2)*(lambda1-lambda2)))
                       +previousA4*currentKTR*(exp(-currentTIME*currentKTR)*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))-exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 2.0))+exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 2.0))-exp(-currentTIME*currentKTR)*currentTIME/((lambda1-currentKTR)*(currentKTR-lambda2)))
                       +previousA1*pow(currentKTR,2.0)*((exp(-currentTIME*currentKTR)*(-1*pow(lambda1,2.0)-lambda1*lambda2+3*lambda1*currentKTR-pow(lambda2, 2.0)+3*lambda2*currentKTR-3*pow(currentKTR,2.0)))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                                       -exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)/(2*(lambda1-currentKTR)*(currentKTR-lambda2))+exp(-currentTIME*currentKTR)*currentTIME*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                                       +exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 3.0))-exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 3.0))));
        
    // Calculate currentA1: Amount in the absorption compartment
    currentA1 = previousA1 * exp(- currentTIME * currentKTR);
    
    // Fill in Amounts and look for other doses
    A4[ counter ] = currentA4;
    A5[ counter ] = currentA5;
    A2[ counter ] = currentA2;
    A3[ counter ] = currentA3;
    A1[ counter ] = currentA1 + (AMT[ counter ] * F1[ counter ]);
    
  } // end for loop
  
  return 0;
}

//--------------------------------------------------------------------------------
// 2 compartment-3transit absorption model via ADVAN-style equations: Cpp Function
//--------------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame TwoCompThreeTransitCpp(DataFrame inputFrame){
  
  //  Create vectors of each element used in function and for constructing output dataframe
  Rcpp::DoubleVector TIME = inputFrame["TIME"];
  Rcpp::DoubleVector AMT = inputFrame["AMT"];
  Rcpp::DoubleVector KTR = inputFrame["KTR"];
  Rcpp::DoubleVector F1 = inputFrame["F1"];
  Rcpp::DoubleVector k20 = inputFrame["k20"];
  Rcpp::DoubleVector k23 = inputFrame["k23"];
  Rcpp::DoubleVector k32 = inputFrame["k32"];
  Rcpp::DoubleVector k30 = inputFrame["k30"];
  Rcpp::DoubleVector A1 = inputFrame["A1"];		//Amount in the absorption compartment
  Rcpp::DoubleVector A4 = inputFrame["A4"];		//Amount in the 1st transit
  Rcpp::DoubleVector A5 = inputFrame["A5"];		//Amount in the 2nd transit
  Rcpp::DoubleVector A6 = inputFrame["A6"];		//Amount in the 3rd transit
  Rcpp::DoubleVector A2 = inputFrame["A2"];		//Amount in central compartment
  Rcpp::DoubleVector A3 = inputFrame["A3"];		//Amount in peripheral compartment
  
  double currentTIME, currentKTR, currentk20, currentk23, currentk32, currentk30, E2, E3, lambda1, lambda2;
  double previousA1,previousA4,previousA5,previousA6, previousA2, previousA3;
  double currentA1,currentA4,currentA5,currentA6, currentA2,currentA3;
  
  // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
  // for counter from 1 to the number of rows in input data frame
  for(int counter = 1; counter < inputFrame.nrows(); counter++){
    
    // pull out all the variables that will be used for calculation
    currentk20  = k20[ counter ];
    currentk23  = k23[ counter ];
    currentk32  = k32[ counter ];
	currentk30  = k30[ counter ];
    currentKTR   = KTR[ counter ];
    currentTIME = TIME[ counter ] - TIME[ counter - 1];
    previousA2  = A2[ counter - 1 ];
    previousA4  = A4[ counter - 1 ];
    previousA5  = A5[ counter - 1 ];
    previousA6  = A6[ counter - 1 ];
    previousA3  = A3[ counter - 1 ];
    previousA1  = A1[ counter - 1 ];
    	
   //calculate hybrid rate constants
    E2          = currentk20 + currentk23;
    E3          = currentk32 + currentk30;
    
    lambda1 = 0.5*((E2+E3)+sqrt(pow((E2+E3),2)-4*(E2*E3-currentk23*currentk32)));
    lambda2 = 0.5*((E2+E3)-sqrt(pow((E2+E3),2)-4*(E2*E3-currentk23*currentk32)));
    
    //Transit compartments
    currentA4 = previousA4*exp(-currentTIME*currentKTR)+currentKTR*previousA1*currentTIME*exp(-currentTIME*currentKTR);
    currentA5 = previousA5*exp(-currentTIME*currentKTR)+currentKTR*previousA4*currentTIME*exp(-currentTIME*currentKTR)+0.5*pow(currentKTR,2.0)*previousA1*pow(currentTIME,2.0)*exp(-currentTIME*currentKTR);
    currentA6 = previousA6*exp(-currentTIME*currentKTR)+currentKTR*previousA5*currentTIME*exp(-currentTIME*currentKTR)+0.5*pow(currentKTR,2.0)*previousA4*pow(currentTIME,2.0)*exp(-currentTIME*currentKTR)+(1.0/6.0)*pow(currentKTR,3.0)*previousA1*pow(currentTIME,3.0)*exp(-currentTIME*currentKTR);
     
   // Calculate currentA2: Amount in the central compartment 
   // Split equation into multiple steps to ensure arithmatic correctness
    currentA2 = (exp(-currentTIME*lambda1)*((previousA2*E3+previousA3*currentk32)-previousA2*lambda1)-exp(-currentTIME*lambda2)*((previousA2*E3+previousA3*currentk32)-previousA2*lambda2))/(lambda2-lambda1);
    
    currentA2 = currentA2 + currentKTR*currentk32*(previousA6*(exp(-currentTIME*currentKTR)/((lambda1-currentKTR)*(lambda2-currentKTR))+exp(-currentTIME*lambda1)/((currentKTR-lambda1)*(lambda2-lambda1))+exp(-currentTIME*lambda2)/((currentKTR-lambda2)*(lambda1-lambda2)))
    +previousA5*currentKTR*(exp(-currentTIME*currentKTR)*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                 -exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 2.0))
                 +exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 2.0))
                 -exp(-currentTIME*currentKTR)*currentTIME/((lambda1-currentKTR)*(currentKTR-lambda2)))
    +previousA4*pow(currentKTR,2.0)*((exp(-currentTIME*currentKTR)*(-1*pow(lambda1,2.0)-lambda1*lambda2+3*lambda1*currentKTR-pow(lambda2, 2.0)+3*lambda2*currentKTR-3*pow(currentKTR,2.0)))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                    -exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)/(2*(lambda1-currentKTR)*(currentKTR-lambda2))+exp(-currentTIME*currentKTR)*currentTIME*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                    +exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 3.0))-exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 3.0)))
    +previousA1*pow(currentKTR,3.0)*((exp(-currentTIME*currentKTR)*currentTIME*(-1*pow(lambda1,2.0)-lambda1*lambda2+3*lambda1*currentKTR-pow(lambda2, 2.0)+3*lambda2*currentKTR-3*pow(currentKTR,2.0)))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                    +exp(-currentTIME*currentKTR)*(-1*pow(lambda1,3.0)-pow(lambda1,2.0)*lambda2+4*pow(lambda1,2.0)*currentKTR-lambda1*pow(lambda2, 2.0)+4*lambda1*lambda2*currentKTR-6*lambda1*pow(currentKTR,2.0)-pow(lambda2, 3.0)+4*pow(lambda2, 2.0)*currentKTR-6*lambda2*pow(currentKTR,2.0)+4*pow(currentKTR,3.0))/(pow((lambda1-currentKTR), 4.0)*pow((currentKTR-lambda2), 4.0))
                    -exp(-currentTIME*currentKTR)*pow(currentTIME, 3.0)/(6*(lambda1-currentKTR)*(currentKTR-lambda2))
                    +exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)*(-lambda1-lambda2+2*currentKTR)/(2*pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))-exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 4.0))
                    +exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 4.0))));
    
    currentA2 = currentA2 + currentKTR*(previousA6*(exp(-currentTIME*currentKTR)*currentKTR/((lambda1-currentKTR)*(currentKTR-lambda2))+exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*(lambda2-currentKTR))-exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*(lambda1-currentKTR)))
    +previousA5*currentKTR*(exp(-currentTIME*currentKTR)*(lambda1*lambda2-pow(currentKTR,2.0))/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))+exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*pow((lambda1-currentKTR), 2.0))-exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*pow((lambda2-currentKTR), 2.0))+exp(-currentTIME*currentKTR)*currentTIME*currentKTR/((lambda1-currentKTR)*(currentKTR-lambda2)))
    +previousA4*pow(currentKTR,2.0)*(exp(-currentTIME*currentKTR)*(pow(lambda1,2.0)*lambda2+lambda1*pow(lambda2, 2.0)-3*lambda1*lambda2*currentKTR+pow(currentKTR,3.0))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                    +exp(-currentTIME*currentKTR)*currentTIME*(lambda1*lambda2-pow(currentKTR,2.0))/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))+exp(-currentTIME*currentKTR)*currentKTR*pow(currentTIME, 2.0)/(2*(lambda1-currentKTR)*(currentKTR-lambda2))
                    -exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*pow((lambda1-currentKTR), 3.0))+exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*pow((lambda2-currentKTR), 3.0)))
    +previousA1*pow(currentKTR,3.0)*(exp(-currentTIME*currentKTR)*currentTIME*(pow(lambda1,2.0)*lambda2+lambda1*pow(lambda2, 2.0)-3*lambda1*lambda2*currentKTR+pow(currentKTR,3.0))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                    +exp(-currentTIME*currentKTR)*(pow(lambda1,3.0)*lambda2+pow(lambda1,2.0)*pow(lambda2, 2.0)-4*pow(lambda1,2.0)*lambda2*currentKTR+lambda1*pow(lambda2, 3.0)-4*lambda1*pow(lambda2, 2.0)*currentKTR+6*lambda1*lambda2*pow(currentKTR,2.0)-pow(currentKTR,4.0))/(pow((lambda1-currentKTR), 4.0)*pow((currentKTR-lambda2), 4.0))
                    +exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)*(lambda1*lambda2-pow(currentKTR,2.0))/(2*pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                    +exp(-currentTIME*currentKTR)*currentKTR*pow(currentTIME, 3.0)/(6*(lambda1-currentKTR)*(currentKTR-lambda2))
                    +exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*pow((lambda1-currentKTR), 4.0))
                    -exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*pow((lambda2-currentKTR), 4.0))));
    
    
	//calculate currentA3: Amount in the peripheral compartment
  	currentA3 = (exp(-currentTIME*lambda1)*((previousA3*E2+currentk23*previousA2)-previousA3*lambda1)-exp(-currentTIME*lambda2)*((previousA3*E2+currentk23*previousA2)-previousA3*lambda2))/(lambda2-lambda1);
    
    currentA3 = currentA3 + currentKTR*currentk23*(previousA6*(exp(-currentTIME*currentKTR)/((lambda1-currentKTR)*(lambda2-currentKTR))+exp(-currentTIME*lambda1)/((currentKTR-lambda1)*(lambda2-lambda1))+exp(-currentTIME*lambda2)/((currentKTR-lambda2)*(lambda1-lambda2)))
                       +previousA5*currentKTR*(exp(-currentTIME*currentKTR)*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))-exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 2.0))+exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 2.0))-exp(-currentTIME*currentKTR)*currentTIME/((lambda1-currentKTR)*(currentKTR-lambda2)))
                       +previousA4*pow(currentKTR,2.0)*((exp(-currentTIME*currentKTR)*(-1*pow(lambda1,2.0)-lambda1*lambda2+3*lambda1*currentKTR-pow(lambda2, 2.0)+3*lambda2*currentKTR-3*pow(currentKTR,2.0)))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                                       -exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)/(2*(lambda1-currentKTR)*(currentKTR-lambda2))+exp(-currentTIME*currentKTR)*currentTIME*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                                       +exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 3.0))-exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 3.0)))
                       +previousA1*pow(currentKTR,3.0)*((exp(-currentTIME*currentKTR)*currentTIME*(-1*pow(lambda1,2.0)-lambda1*lambda2+3*lambda1*currentKTR-pow(lambda2, 2.0)+3*lambda2*currentKTR-3*pow(currentKTR,2.0)))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                                       +exp(-currentTIME*currentKTR)*(-1*pow(lambda1,3.0)-pow(lambda1,2.0)*lambda2+4*pow(lambda1,2.0)*currentKTR-lambda1*pow(lambda2, 2.0)+4*lambda1*lambda2*currentKTR-6*lambda1*pow(currentKTR,2.0)-pow(lambda2, 3.0)+4*pow(lambda2, 2.0)*currentKTR-6*lambda2*pow(currentKTR,2.0)+4*pow(currentKTR,3.0))/(pow((lambda1-currentKTR), 4.0)*pow((currentKTR-lambda2), 4.0))
                                       -exp(-currentTIME*currentKTR)*pow(currentTIME, 3.0)/(6*(lambda1-currentKTR)*(currentKTR-lambda2))
                                       +exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)*(-lambda1-lambda2+2*currentKTR)/(2*pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))-exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 4.0))
                                       +exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 4.0))));
        
    // Calculate currentA1: Amount in the absorption compartment
    currentA1 = previousA1 * exp(- currentTIME * currentKTR);
    
    // Fill in Amounts and look for other doses
    A4[ counter ] = currentA4;
    A5[ counter ] = currentA5;
    A6[ counter ] = currentA6;
    A2[ counter ] = currentA2;
    A3[ counter ] = currentA3;
    A1[ counter ] = currentA1 + (AMT[ counter ] * F1[ counter ]);
    
  } // end for loop
  
  return 0;
}

//--------------------------------------------------------------------------------
// 2 compartment-4transit absorption model via ADVAN-style equations: Cpp Function
//--------------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame TwoCompFourTransitCpp(DataFrame inputFrame){
  
  //  Create vectors of each element used in function and for constructing output dataframe
  Rcpp::DoubleVector TIME = inputFrame["TIME"];
  Rcpp::DoubleVector AMT = inputFrame["AMT"];
  Rcpp::DoubleVector KTR = inputFrame["KTR"];
  Rcpp::DoubleVector F1 = inputFrame["F1"];
  Rcpp::DoubleVector k20 = inputFrame["k20"];
  Rcpp::DoubleVector k23 = inputFrame["k23"];
  Rcpp::DoubleVector k32 = inputFrame["k32"];
  Rcpp::DoubleVector k30 = inputFrame["k30"];
  Rcpp::DoubleVector A1 = inputFrame["A1"];		//Amount in the absorption compartment
  Rcpp::DoubleVector A4 = inputFrame["A4"];		//Amount in the 1st transit
  Rcpp::DoubleVector A5 = inputFrame["A5"];		//Amount in the 2nd transit
  Rcpp::DoubleVector A6 = inputFrame["A6"];		//Amount in the 3rd transit
  Rcpp::DoubleVector A7 = inputFrame["A7"];		//Amount in the 4th transit
  Rcpp::DoubleVector A2 = inputFrame["A2"];		//Amount in central compartment
  Rcpp::DoubleVector A3 = inputFrame["A3"];		//Amount in peripheral compartment
  
  double currentTIME, currentKTR, currentk20, currentk23, currentk32, currentk30, E2, E3, lambda1, lambda2;
  double previousA1,previousA4,previousA5,previousA6,previousA7, previousA2, previousA3;
  double currentA1,currentA4,currentA5,currentA6,currentA7, currentA2,currentA3;
  
  // in C++ arrays start at index 0, so to start at 2nd row need to set counter to 1
  // for counter from 1 to the number of rows in input data frame
  for(int counter = 1; counter < inputFrame.nrows(); counter++){
    
    // pull out all the variables that will be used for calculation
    currentk20  = k20[ counter ];
    currentk23  = k23[ counter ];
    currentk32  = k32[ counter ];
	currentk30  = k30[ counter ];
    currentKTR   = KTR[ counter ];
    currentTIME = TIME[ counter ] - TIME[ counter - 1];
    previousA2  = A2[ counter - 1 ];
    previousA4  = A4[ counter - 1 ];
    previousA5  = A5[ counter - 1 ];
    previousA6  = A6[ counter - 1 ];
    previousA7  = A7[ counter - 1 ];
    previousA3  = A3[ counter - 1 ];
    previousA1  = A1[ counter - 1 ];
    	
   //calculate hybrid rate constants
    E2          = currentk20 + currentk23;
    E3          = currentk32 + currentk30;
    
    lambda1 = 0.5*((E2+E3)+sqrt(pow((E2+E3),2)-4*(E2*E3-currentk23*currentk32)));
    lambda2 = 0.5*((E2+E3)-sqrt(pow((E2+E3),2)-4*(E2*E3-currentk23*currentk32)));
    
    //Transit compartments
    currentA4 = previousA4*exp(-currentTIME*currentKTR)+currentKTR*previousA1*currentTIME*exp(-currentTIME*currentKTR);
    currentA5 = previousA5*exp(-currentTIME*currentKTR)+currentKTR*previousA4*currentTIME*exp(-currentTIME*currentKTR)+0.5*pow(currentKTR,2.0)*previousA1*pow(currentTIME,2.0)*exp(-currentTIME*currentKTR);
    currentA6 = previousA6*exp(-currentTIME*currentKTR)+currentKTR*previousA5*currentTIME*exp(-currentTIME*currentKTR)+0.5*pow(currentKTR,2.0)*previousA4*pow(currentTIME,2.0)*exp(-currentTIME*currentKTR)+(1.0/6.0)*pow(currentKTR,3.0)*previousA1*pow(currentTIME,3.0)*exp(-currentTIME*currentKTR);
    currentA7 = previousA7*exp(-currentTIME*currentKTR)+currentKTR*previousA6*currentTIME*exp(-currentTIME*currentKTR)+0.5*pow(currentKTR,2.0)*previousA5*pow(currentTIME,2.0)*exp(-currentTIME*currentKTR)+(1.0/6.0)*pow(currentKTR,3.0)*previousA4*pow(currentTIME,3.0)*exp(-currentTIME*currentKTR)+(1.0/24.0)*pow(currentKTR,4.0)*previousA1*pow(currentTIME,4.0)*exp(-currentTIME*currentKTR);
     
   // Calculate currentA2: Amount in the central compartment 
   // Split equation into multiple steps to ensure arithmatic correctness
    currentA2 = (exp(-currentTIME*lambda1)*((previousA2*E3+previousA3*currentk32)-previousA2*lambda1)-exp(-currentTIME*lambda2)*((previousA2*E3+previousA3*currentk32)-previousA2*lambda2))/(lambda2-lambda1);
    
    currentA2 = currentA2 + currentKTR*currentk32*(previousA7*(exp(-currentTIME*currentKTR)/((lambda1-currentKTR)*(lambda2-currentKTR))+exp(-currentTIME*lambda1)/((currentKTR-lambda1)*(lambda2-lambda1))+exp(-currentTIME*lambda2)/((currentKTR-lambda2)*(lambda1-lambda2)))
    +previousA6*currentKTR*(exp(-currentTIME*currentKTR)*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                 -exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 2.0))
                 +exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 2.0))
                 -exp(-currentTIME*currentKTR)*currentTIME/((lambda1-currentKTR)*(currentKTR-lambda2)))
    +previousA5*pow(currentKTR,2.0)*((exp(-currentTIME*currentKTR)*(-1*pow(lambda1,2.0)-lambda1*lambda2+3*lambda1*currentKTR-pow(lambda2, 2.0)+3*lambda2*currentKTR-3*pow(currentKTR,2.0)))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                    -exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)/(2*(lambda1-currentKTR)*(currentKTR-lambda2))+exp(-currentTIME*currentKTR)*currentTIME*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                    +exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 3.0))-exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 3.0)))
    +previousA4*pow(currentKTR,3.0)*((exp(-currentTIME*currentKTR)*currentTIME*(-1*pow(lambda1,2.0)-lambda1*lambda2+3*lambda1*currentKTR-pow(lambda2, 2.0)+3*lambda2*currentKTR-3*pow(currentKTR,2.0)))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                    +exp(-currentTIME*currentKTR)*(-1*pow(lambda1,3.0)-pow(lambda1,2.0)*lambda2+4*pow(lambda1,2.0)*currentKTR-lambda1*pow(lambda2, 2.0)+4*lambda1*lambda2*currentKTR-6*lambda1*pow(currentKTR,2.0)-pow(lambda2, 3.0)+4*pow(lambda2, 2.0)*currentKTR-6*lambda2*pow(currentKTR,2.0)+4*pow(currentKTR,3.0))/(pow((lambda1-currentKTR), 4.0)*pow((currentKTR-lambda2), 4.0))
                    -exp(-currentTIME*currentKTR)*pow(currentTIME, 3.0)/(6*(lambda1-currentKTR)*(currentKTR-lambda2))
                    +exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)*(-lambda1-lambda2+2*currentKTR)/(2*pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))-exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 4.0))
                    +exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 4.0)))
    +previousA1*pow(currentKTR,4.0)*((exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)*(-1*pow(lambda1,2.0)-lambda1*lambda2+3*lambda1*currentKTR-pow(lambda2, 2.0)+3*lambda2*currentKTR-3*pow(currentKTR,2.0)))/(2*pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                    +exp(-currentTIME*currentKTR)*currentTIME*(-1*pow(lambda1,3.0)-pow(lambda1,2.0)*lambda2+4*pow(lambda1,2.0)*currentKTR-lambda1*pow(lambda2, 2.0)+4*lambda1*lambda2*currentKTR-6*lambda1*pow(currentKTR,2.0)-pow(lambda2, 3.0)+4*pow(lambda2, 2.0)*currentKTR-6*lambda2*pow(currentKTR,2.0)+4*pow(currentKTR,3.0))/(pow((lambda1-currentKTR), 4.0)*pow((currentKTR-lambda2), 4.0))
                    +(exp(-currentTIME*currentKTR)/(pow((lambda1-currentKTR), 5.0)*pow((currentKTR-lambda2), 5.0)))*(-1*pow(lambda1,4.0)-pow(lambda1,3.0)*lambda2+5*pow(lambda1,3.0)*currentKTR-pow(lambda1,2.0)*pow(lambda2, 2.0)+5*pow(lambda1,2.0)*lambda2*currentKTR-10*pow(lambda1,2.0)*pow(currentKTR,2.0)-lambda1*pow(lambda2, 3.0)+5*lambda1*pow(lambda2, 2.0)*currentKTR-10*lambda1*lambda2*pow(currentKTR,2.0)+10*lambda1*pow(currentKTR,3.0)-pow(lambda2, 4.0)+5*pow(lambda2, 3.0)*currentKTR-10*pow(lambda2, 2.0)*pow(currentKTR,2.0)+10*lambda2*pow(currentKTR,3.0)-5*pow(currentKTR,4.0))
                    -exp(-currentTIME*currentKTR)*pow(currentTIME, 4.0)/(24*(lambda1-currentKTR)*(currentKTR-lambda2))+exp(-currentTIME*currentKTR)*pow(currentTIME, 3.0)*(-lambda1-lambda2+2*currentKTR)/(6*pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                    +exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 5.0))-exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 5.0))));
    
    currentA2 = currentA2 + currentKTR*(previousA7*(exp(-currentTIME*currentKTR)*currentKTR/((lambda1-currentKTR)*(currentKTR-lambda2))+exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*(lambda2-currentKTR))-exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*(lambda1-currentKTR)))
    +previousA6*currentKTR*(exp(-currentTIME*currentKTR)*(lambda1*lambda2-pow(currentKTR,2.0))/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))+exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*pow((lambda1-currentKTR), 2.0))-exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*pow((lambda2-currentKTR), 2.0))+exp(-currentTIME*currentKTR)*currentTIME*currentKTR/((lambda1-currentKTR)*(currentKTR-lambda2)))
    +previousA5*pow(currentKTR,2.0)*(exp(-currentTIME*currentKTR)*(pow(lambda1,2.0)*lambda2+lambda1*pow(lambda2, 2.0)-3*lambda1*lambda2*currentKTR+pow(currentKTR,3.0))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                    +exp(-currentTIME*currentKTR)*currentTIME*(lambda1*lambda2-pow(currentKTR,2.0))/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))+exp(-currentTIME*currentKTR)*currentKTR*pow(currentTIME, 2.0)/(2*(lambda1-currentKTR)*(currentKTR-lambda2))
                    -exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*pow((lambda1-currentKTR), 3.0))+exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*pow((lambda2-currentKTR), 3.0)))
    +previousA4*pow(currentKTR,3.0)*(exp(-currentTIME*currentKTR)*currentTIME*(pow(lambda1,2.0)*lambda2+lambda1*pow(lambda2, 2.0)-3*lambda1*lambda2*currentKTR+pow(currentKTR,3.0))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                    +exp(-currentTIME*currentKTR)*(pow(lambda1,3.0)*lambda2+pow(lambda1,2.0)*pow(lambda2, 2.0)-4*pow(lambda1,2.0)*lambda2*currentKTR+lambda1*pow(lambda2, 3.0)-4*lambda1*pow(lambda2, 2.0)*currentKTR+6*lambda1*lambda2*pow(currentKTR,2.0)-pow(currentKTR,4.0))/(pow((lambda1-currentKTR), 4.0)*pow((currentKTR-lambda2), 4.0))
                    +exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)*(lambda1*lambda2-pow(currentKTR,2.0))/(2*pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                    +exp(-currentTIME*currentKTR)*currentKTR*pow(currentTIME, 3.0)/(6*(lambda1-currentKTR)*(currentKTR-lambda2))
                    +exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*pow((lambda1-currentKTR), 4.0))
                    -exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*pow((lambda2-currentKTR), 4.0)))
    +previousA1*pow(currentKTR,4.0)*((exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)*(pow(lambda1,2.0)*lambda2+lambda1*pow(lambda2, 2.0)-3*lambda1*lambda2*currentKTR+pow(currentKTR,3.0)))/(2*pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                    +exp(-currentTIME*currentKTR)*currentTIME*(pow(lambda1,3.0)*lambda2+pow(lambda1,2.0)*pow(lambda2, 2.0)-4*pow(lambda1,2.0)*lambda2*currentKTR+lambda1*pow(lambda2, 3.0)-4*lambda1*pow(lambda2, 2.0)*currentKTR+6*lambda1*lambda2*pow(currentKTR,2.0)-pow(currentKTR,4.0))/(pow((lambda1-currentKTR), 4.0)*pow((currentKTR-lambda2), 4.0))
                    +(exp(-currentTIME*currentKTR)/(pow((lambda1-currentKTR), 5.0)*pow((currentKTR-lambda2), 5.0)))*(pow(lambda1,4.0)*lambda2+pow(lambda1,3.0)*pow(lambda2, 2.0)-5*pow(lambda1,3.0)*lambda2*currentKTR+pow(lambda1,2.0)*pow(lambda2, 3.0)-5*pow(lambda1,2.0)*pow(lambda2, 2.0)*currentKTR+10*pow(lambda1,2.0)*lambda2*pow(currentKTR,2.0)+lambda1*pow(lambda2, 4.0)-5*lambda1*pow(lambda2, 3.0)*currentKTR+10*lambda1*pow(lambda2, 2.0)*pow(currentKTR,2.0)-10*lambda1*lambda2*pow(currentKTR,3.0)+pow(currentKTR,5.0))
                    +exp(-currentTIME*currentKTR)*pow(currentTIME, 3.0)*(lambda1*lambda2-pow(currentKTR,2.0))/(6*pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                    +exp(-currentTIME*currentKTR)*currentKTR*pow(currentTIME, 4.0)/(24*(lambda1-currentKTR)*(currentKTR-lambda2))
                    -exp(-currentTIME*lambda1)*lambda1/((lambda1-lambda2)*pow((lambda1-currentKTR), 5.0))
                    +exp(-currentTIME*lambda2)*lambda2/((lambda1-lambda2)*pow((lambda2-currentKTR), 5.0))));
    
    
	//calculate currentA3: Amount in the peripheral compartment
  	currentA3 = (exp(-currentTIME*lambda1)*((previousA3*E2+currentk23*previousA2)-previousA3*lambda1)-exp(-currentTIME*lambda2)*((previousA3*E2+currentk23*previousA2)-previousA3*lambda2))/(lambda2-lambda1);
    
    currentA3 = currentA3 + currentKTR*currentk23*(previousA7*(exp(-currentTIME*currentKTR)/((lambda1-currentKTR)*(lambda2-currentKTR))+exp(-currentTIME*lambda1)/((currentKTR-lambda1)*(lambda2-lambda1))+exp(-currentTIME*lambda2)/((currentKTR-lambda2)*(lambda1-lambda2)))
                       +previousA6*currentKTR*(exp(-currentTIME*currentKTR)*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))-exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 2.0))+exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 2.0))-exp(-currentTIME*currentKTR)*currentTIME/((lambda1-currentKTR)*(currentKTR-lambda2)))
                       +previousA5*pow(currentKTR,2.0)*((exp(-currentTIME*currentKTR)*(-1*pow(lambda1,2.0)-lambda1*lambda2+3*lambda1*currentKTR-pow(lambda2, 2.0)+3*lambda2*currentKTR-3*pow(currentKTR,2.0)))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                                       -exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)/(2*(lambda1-currentKTR)*(currentKTR-lambda2))+exp(-currentTIME*currentKTR)*currentTIME*(-lambda1-lambda2+2*currentKTR)/(pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                                       +exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 3.0))-exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 3.0)))
                       +previousA4*pow(currentKTR,3.0)*((exp(-currentTIME*currentKTR)*currentTIME*(-1*pow(lambda1,2.0)-lambda1*lambda2+3*lambda1*currentKTR-pow(lambda2, 2.0)+3*lambda2*currentKTR-3*pow(currentKTR,2.0)))/(pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                                       +exp(-currentTIME*currentKTR)*(-1*pow(lambda1,3.0)-pow(lambda1,2.0)*lambda2+4*pow(lambda1,2.0)*currentKTR-lambda1*pow(lambda2, 2.0)+4*lambda1*lambda2*currentKTR-6*lambda1*pow(currentKTR,2.0)-pow(lambda2, 3.0)+4*pow(lambda2, 2.0)*currentKTR-6*lambda2*pow(currentKTR,2.0)+4*pow(currentKTR,3.0))/(pow((lambda1-currentKTR), 4.0)*pow((currentKTR-lambda2), 4.0))
                                       -exp(-currentTIME*currentKTR)*pow(currentTIME, 3.0)/(6*(lambda1-currentKTR)*(currentKTR-lambda2))
                                       +exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)*(-lambda1-lambda2+2*currentKTR)/(2*pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))-exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 4.0))
                                       +exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 4.0)))
                       +previousA1*pow(currentKTR,4.0)*((exp(-currentTIME*currentKTR)*pow(currentTIME, 2.0)*(-1*pow(lambda1,2.0)-lambda1*lambda2+3*lambda1*currentKTR-pow(lambda2, 2.0)+3*lambda2*currentKTR-3*pow(currentKTR,2.0)))/(2*pow((lambda1-currentKTR), 3.0)*pow((currentKTR-lambda2), 3.0))
                                       +exp(-currentTIME*currentKTR)*currentTIME*(-1*pow(lambda1,3.0)-pow(lambda1,2.0)*lambda2+4*pow(lambda1,2.0)*currentKTR-lambda1*pow(lambda2, 2.0)+4*lambda1*lambda2*currentKTR-6*lambda1*pow(currentKTR,2.0)-pow(lambda2, 3.0)+4*pow(lambda2, 2.0)*currentKTR-6*lambda2*pow(currentKTR,2.0)+4*pow(currentKTR,3.0))/(pow((lambda1-currentKTR), 4.0)*pow((currentKTR-lambda2), 4.0))
                                       +(exp(-currentTIME*currentKTR)/(pow((lambda1-currentKTR), 5.0)*pow((currentKTR-lambda2), 5.0)))*(-1*pow(lambda1,4.0)-pow(lambda1,3.0)*lambda2+5*pow(lambda1,3.0)*currentKTR-pow(lambda1,2.0)*pow(lambda2, 2.0)+5*pow(lambda1,2.0)*lambda2*currentKTR-10*pow(lambda1,2.0)*pow(currentKTR,2.0)-lambda1*pow(lambda2, 3.0)+5*lambda1*pow(lambda2, 2.0)*currentKTR-10*lambda1*lambda2*pow(currentKTR,2.0)+10*lambda1*pow(currentKTR,3.0)-pow(lambda2, 4.0)+5*pow(lambda2, 3.0)*currentKTR-10*pow(lambda2, 2.0)*pow(currentKTR,2.0)+10*lambda2*pow(currentKTR,3.0)-5*pow(currentKTR,4.0))
                                       -exp(-currentTIME*currentKTR)*pow(currentTIME, 4.0)/(24*(lambda1-currentKTR)*(currentKTR-lambda2))+exp(-currentTIME*currentKTR)*pow(currentTIME, 3.0)*(-lambda1-lambda2+2*currentKTR)/(6*pow((lambda1-currentKTR), 2.0)*pow((currentKTR-lambda2), 2.0))
                                       +exp(-currentTIME*lambda1)/((lambda1-lambda2)*pow((lambda1-currentKTR), 5.0))-exp(-currentTIME*lambda2)/((lambda1-lambda2)*pow((lambda2-currentKTR), 5.0))));
    
    
    
    // Calculate currentA1: Amount in the absorption compartment
    currentA1 = previousA1 * exp(- currentTIME * currentKTR);
    
    // Fill in Amounts and look for other doses
    A4[ counter ] = currentA4;
    A5[ counter ] = currentA5;
    A6[ counter ] = currentA6;
    A7[ counter ] = currentA7;
    A2[ counter ] = currentA2;
    A3[ counter ] = currentA3;
    A1[ counter ] = currentA1 + (AMT[ counter ] * F1[ counter ]);
    
  } // end for loop
  
  return 0;
}
