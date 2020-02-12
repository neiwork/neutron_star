#include "modelFunctions.h"

   // --------------------------------------------------------------------------------------------------------
    // Definicion de la funcion temperatura
    
    double Temp(double Md)
    {
        double temperature;
        
        temperature = pow(10.0, 7.887 + 0.528*(1 - exp(-0.899*(log10(Md) + 11))));
        return temperature;
    }

    //---------------------------------------------------------------------------------------------------------
    // Definicion de la funcion density 1
    
    double density1(double pos)
    {
       // double den1_a0, den1_a1, den1_a2;
        double den1;
        
        
        //den1_a0 = 1.98017e6;
       // den1_a1 = 614102;
       // den1_a2 = 12106.8;
       //  den1_a0 + den1_a1 * pos + den1_a2 * pow(pos,2);     

 //      den1 = 1.516935688136645e6 + 816817.8832152047 * pos;
 
        den1 = 1.609456e6 + 7920 * pow(pos,3);
       
       return den1;
    }
    // ---------------------------------------------------------------------------------------------------------------------------
   // Definicion de la funcion density2

//    double density2(double pos)
//    {
    
//    double den2_a0, den2_a1;
 //   double den2;
  
   
   
   
 //       den2_a0 = 6.001402381791151e7;
 //       den2_a1 = 3465.2776436715776;
        
        
  //      den2 = den2_a0 + den2_a1 * pow(pos,3.2);           
  //      return den2;
    
  //  }
    
    //------------------------------------------------------------------------------------------------------------------------------
    // Definicion de la funcion density3
    
    double density3(double pos)
    {
        double den3;

    //den3 = pow(pos,4.929545434398888);
    
      den3 = pow(pos,4.9494066);  

        return den3;
    }
    // __________________________________________________________________________________________________________________________________
    
    // Definicion de la funcion sigma_ph
    
        double sigma_ph(double den, double temp)
        {
        double  x13, x, Td, Z, A, cond_ph;       // we consider that the core is composed of iron 56
        
        Z = 26;                 // charge number
        A = 56;                 // mass number
        
        x = (Z * den)/(1.e6 *A);         
        x13 = pow(x, 1./3.);       // Fermi momentum
        Td = 2.4e6 * sqrt(2*Z/A) * pow(x13,3./2.);         // Debye temperature 

        cond_ph = 1.21e28 * (pow(x13,2) / pow(temp,2)) * sqrt(pow(temp,2) + 0.084 * pow(Td,2)); 
        return cond_ph;
        }
        
    // ________________________________________________________________________________________________________________

    // Definicion de la funcion sigma_imp

        double sigma_imp(double den)
        {
        double x, x13, Z, A, cond_imp;
        double lambda_imp = 2;
        double Q = 0.01;
        
        Z = 30;                 // charge number
        A = 130;                 // mass number
        
        x = (Z * den)/(1.e6 *A);         
        x13 = pow(x,1./3.);       // Fermi momentum
        
        cond_imp = 8.53e21 *x13* (1/lambda_imp)*(Z/Q);
        return cond_imp;
 
        }
		
// ____________________________________________________________________________________
	double conductivity(double den,double temp)	
   {     
	   double den_cr = 1.e12;					//cambio en la parametrizacion de la conductividad
		
		if (den < den_cr){
			return sigma_ph(den,temp); 
        }
     
		else {
			return sigma_imp(den);
		}
}
// ____________________________________________________________________________________

//      Definicion de la funcion de transicion 

         double sigma_trans(double den)
        {
            double cond_trans;
            
           cond_trans =-3.083026340166874e27 + 3.7654540810016265e15 * den;
           return cond_trans;
        }