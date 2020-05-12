#include "modelFunctions.h"

#include "globalParameters.h"

#include <fmath/interpolation.h>




using namespace std;

//extern "C" void CONDCONV_(double T6,double RHO,double B,double Zion,double CMI,double Zimp, // input
//				double SIGMA,double CKAPPA,double QJ,double SIGMAT,double CKAPPAT,double QJT,double SIGMAH,double CKAPPAH,double QJH);

extern "C" void condconv_(double* T6,double* RHO,double* B,double* Zion,double* CMI,double* Zimp, // input
				double* SIGMA,double* CKAPPA,double* QJ,double* SIGMAT,double* CKAPPAT,double* QJT,double* SIGMAH,double* CKAPPAH,double* QJH);



/*double eqState()
{
	
	
interpol(u, Ecreator, Ncreator, Ncreator.size() - 1);
}
*/




double atomicNumberZ(double density)
{
	double res = interpol(density, rho_Z, Z, Z.size() - 1);
	return res;
}

double massNumberA(double density)
{
	double Z = atomicNumberZ(density);
	
	double nucleonNumber = interpol(density, rho_N, N, N.size() - 1);
	double A = nucleonNumber-Z;
	return A;
}



double condPothekin(double den, double temp, double B)
	{
		
	const double CONVSIGM = 8.99e9; // ! conversion S/m (SI) to 1/s (CGS)
	const double CONVOPAC = 3.02242e-4; //conv.th.conductivity [erg cm^{-1} s^{-1} K^{-1}] to opac. [cm^2/g]

	//write(*,'('' Charge and atomic mass of ions (Z,A): ''$)')
	double Zion = atomicNumberZ(den); //30; //Z;
	double CMI = massNumberA(den); //130;//A;
	
	//write(*,'('' Impurity parameter (effective Z): ''$)')
	double Zimp = Q;

	//print*,'Note: if B is nonzero, then ee-collisions are ignored'
	//write(*,'('' Magnetic field B (in Gauss): ''$)')
	//double B12 = B/1.0e12;

	
	//write(*,'('' Temperature T6 (in 10^6 K): ''$)')
	double T6 = temp/1.0e6;

	//write(*,'('' Density of ions rho (in g/cm^3): ''$)')
	double RHO = den;
		
	double SIGMA, CKAPPA, QJ, SIGMAT, CKAPPAT, QJT, SIGMAH, CKAPPAH, QJH;
	
	condconv_(&T6,&RHO,&B,&Zion,&CMI,&Zimp, // input
				&SIGMA,&CKAPPA,&QJ,&SIGMAT,&CKAPPAT,&QJT,&SIGMAH,&CKAPPAH,&QJH);
				
	//CONDCONV_(T6,RHO,B,Zion,CMI,Zimp, // input
	//		SIGMA,CKAPPA,QJ,SIGMAT,CKAPPAT,QJT,SIGMAH,CKAPPAH,QJH);
  
/*----------------------   OUTPUT:
C%C      T=T6*1.D6
C%C      if (B.gt.0.) then
C%C         write(*,113)
C%C         write(*,111) RHO,T,SIGMA,CKAPPA,QJ,
C%C     *      SIGMAT,CKAPPAT,QJT,
C%C     *      SIGMAH,CKAPPAH,QJH
C%C      else
C%C         write(*,114)
C%C         write(*,111) RHO,T,SIGMA,CKAPPA,QJ
C%C      endif
C%C      OPAClg=dlog10(CONVOPAC/CKAPPA/RHO*T**3)
C%C      write(*,'('' log_{10}(opacity[cm^2/g]) ='',F8.3)') OPAClg
C%C      SIGMAlg=dlog10(SIGMA/CONVSIGM)
C%C      write(*,'('' log_{10}(el.conductivity[S/m] ='',F8.3)') SIGMAlg  */

	double sigma_cgs = SIGMA;              // me parece que el programa ya devuelve en cgs
	return sigma_cgs;
}

   // --------------------------------------------------------------------------------------------------------
    // Definicion de la funcion temperatura
    
    double Temp(double Md)
    {
        double temperature;
        
        temperature = pow(10.0, 7.887 + 0.528*(1.0 - exp(-0.899*(log10(Md) + 11.0))));
        return temperature;
    }

////////////////////////////////////

	//parametrizacion de Geppper&Urpin 1994
	double density(double r)
	{
		double pos = r_f - r;
		double den;
		
		if (pos <= 100.0) {
			den = 1.609456e6 + 7920.0 * P3(pos);   //den1  density1(r_f - r);
		}
		else {
			den = pow(pos,4.9494066);   //den3 density3(r_f - r);
		}
		return den;
	}
    // __________________________________________________________________________________________________________________________________
    
    // Definicion de la funcion sigma_ph
    
        /*double sigma_ph(double den, double temp)
        {
        double  x13, x, Td, Z, A, cond_ph;       // we consider that the core is composed of iron 56
        
        Z = 26;                 // charge number
        A = 56;                 // mass number
        
        x = (Z * den)/(1.e6 *A);         
        x13 = pow(x, 1./3.);       // Fermi momentum
		
        Td = 2.4e6 * sqrt(2.0*Z/A) * pow(x13,3./2.);         // Debye temperature 

        cond_ph = 1.21e28 * P2(x13/temp) * sqrt(P2(temp) + 0.084 * P2(Td)); 
        return cond_ph;
        }*/
		
	double sigma_ph(double den, double temp)
	{
		// we consider that the core is composed of iron 56
		
		double Z = 26;                 // charge number
		double A = 56;                 // mass number
		
		double x = pow( (Z * den/(1.e6 *A)),1./3.);   // Fermi momentum
  
		//double ni = den/(atomicMassUnit*A);
		//double a = pow(3.0/(4.0*pi*ni),1./3.);
		double Gamma = 22.75*P2(Z)*(1.0e6/temp)*pow(den/(A*1.0e6),1./3.); //P2(Z*electronCharge)/(a*boltzmann*temp);

		double Tm = 3.04e7*pow((Z/26),5./3.)*(170./Gamma)*x; //temperatura de fusion
		
		/*if(temp > Tm){ //dispersion por iones p
			double lambda_p = log10(Gamma); //10.0; //aprox, en realidad es log(12.0*pi*ne*Lambda_De^3/Z)
			double cond_p = 8.53e21*P3(x)/(Z*lambda_p*(1.0+P2(x)));
			return cond_p;
		}
		else{ //dispersion por phonons ph*/
			double Td = 2.4e6 * sqrt(2.0*Z/A) * pow(x,3./2.);       // Debye temperature
			double u = 0.45*temp/Td;
			double cond_ph = 1.21e28 * (pow(x,4)/(2.0+P2(x))) * sqrt((P2(u) + 0.017)/(temp*u)); 
			return cond_ph;
		//}

		
	}
        
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
		
		
		
		
		   // ---------------------------------------------------------------------------------------------------------------------------
/* 
 * // Definicion de la funcion density2

//    double density2(double pos)
//    {
    
//    double den2_a0, den2_a1;
 //   double den2;
  
   
   
   
 //       den2_a0 = 6.001402381791151e7;
 //       den2_a1 = 3465.2776436715776;
        
        
  //      den2 = den2_a0 + den2_a1 * pow(pos,3.2);           
  //      return den2;
    
  //  }
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
 
        
       
       return den1;
    }
 
    
    //------------------------------------------------------------------------------------------------------------------------------
    // Definicion de la funcion density3
    
    double density3(double pos)
    {
        double den3;

    //den3 = pow(pos,4.929545434398888);
    
      

        return den3;
    } */