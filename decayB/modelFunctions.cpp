#include "modelFunctions.h"

#include "globalParameters.h"
   // --------------------------------------------------------------------------------------------------------
    // Definicion de la funcion temperatura
    
    double Temp(double Md)
    {
        double temperature;
        
        temperature = pow(10.0, 7.887 + 0.528*(1.0 - exp(-0.899*(log10(Md) + 11.0))));
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
 
        den1 = 1.609456e6 + 7920.0 * P3(pos);
       
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
	
	////////////////////////////////////
	//parametrizacion de Geppper&Urpin 1994
	double density(double r)
	{
		double den;
		if (r_f - r <= 100.0) {
			den = density1(r_f - r);
		}
		else {
			den = density3(r_f - r);
		}
		return den;
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
		
        Td = 2.4e6 * sqrt(2.0*Z/A) * pow(x13,3./2.);         // Debye temperature 

        cond_ph = 1.21e28 * P2(x13/temp) * sqrt(P2(temp) + 0.084 * P2(Td)); 
        return cond_ph;
        }
        
 // -------------------------------------------------------------------------------------------------------------

     //   double sigma_ph(double den, double temp)
        /*{
        double  x13, x, Td, Z, A, cond_ph;       // we consider that the core is composed of iron 56
        
        Z = 26;                 // charge number
        A = 56;                 // mass number
        
        x = (Z * den)/(1.e6 *A);         
        x13 = pow(x, 1./3.);       // Fermi momentum
        Td = 2.4e6 * sqrt(2*Z/A) * pow(x13,3./2.);         // Debye temperature 

        cond_ph = 1.21e28 * (pow(x13,2) / pow(temp,2)) * sqrt(pow(temp,2) + 0.084 * pow(Td,2)); 
        return cond_ph;
        }*/
        
// ------------------------------------------------------------------------------------------------
// Parametrizacion Flor
// -----------------------------------------------------------------------------------
		
//	double sigma_ph(double den, double temp)
//	{
		// we consider that the core is composed of iron 56
		
//		double Z = 26;                 // charge number
//		double A = 56;                 // mass number
		
//		double x = pow( (Z * den/(1.e6 *A)),1./3.);   // Fermi momentum
  
		//double ni = den/(atomicMassUnit*A);
		//double a = pow(3.0/(4.0*pi*ni),1./3.);
///		double Gamma = 22.75*P2(Z)*(1.0e6/temp)*pow(den/(A*1.0e6),1./3.); //P2(Z*electronCharge)/(a*boltzmann*temp);

//		double Tm = 3.04e7*pow((Z/26),5./3.)*(170./Gamma)*x; //temperatura de fusion
		
		/*if(temp > Tm){ //dispersion por iones p
			double lambda_p = log10(Gamma); //10.0; //aprox, en realidad es log(12.0*pi*ne*Lambda_De^3/Z)
			double cond_p = 8.53e21*P3(x)/(Z*lambda_p*(1.0+P2(x)));
			return cond_p;
		}
		else{ //dispersion por phonons ph*/
//			double Td = 2.4e6 * sqrt(2.0*Z/A) * pow(x,3./2.);       // Debye temperature
//			double u = 0.45*temp/Td;
//			double cond_ph = 1.21e28 * (pow(x,4)/(2.0+P2(x))) * sqrt((P2(u) + 0.017)/(temp*u)); 
//			return cond_ph;
		//}

		
//	}
        
    // ________________________________________________________________________________________________________________

    // Definicion de la funcion sigma_imp

        double sigma_imp(double den)
        {
        double x, x13, Z, A, cond_imp;
        double lambda_imp = 2;   //esto vale para rho > 1.0e5
        
        
        Z = 30;                 // charge number
        A = 130;                 // mass number
        
        x = (Z * den)/(1.e6 *A);         
        x13 = pow(x,1./3.);       // Fermi momentum
        
        cond_imp = 8.53e21 *x13* (1.0/lambda_imp)*(Z/Q);
        return cond_imp;
 
        }
		
// ____________________________________________________________________________________
	double conductivity(double den,double temp)	
   {     
	   double den_cr = 1.e12;					//cambio en la parametrizacion de la conductividad
		
		double cond = 1.0/(1.0/sigma_ph(den,temp)+1.0/sigma_imp(den));  
		return cond;
		/*if (den < den_cr){
			return sigma_ph(den,temp); 
        }
     
		else {
			return sigma_imp(den);
		}*/
}
// ____________________________________________________________________________________

//      Definicion de la funcion de transicion 

         double sigma_trans(double den)
        {
            double cond_trans;
            
           cond_trans =-3.083026340166874e27 + 3.7654540810016265e15 * den;
           return cond_trans;
        }