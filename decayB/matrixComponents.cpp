#include "matrixComponents.h"

#include "modelFunctions.h"
#include "globalParameters.h"

#include <string>
#include <iostream>
#include <fstream>

//----------------------------------------------------------------------------------------------------------------
    // Calculo la matriz BB y AA
    // _____________________________________________________________________________________________________________
    
	
void setMatrix(Matrix& bb, Vector& cd, Vector& ad, Vector& bd)
{
	
	std::ofstream myfile3;
    std::ofstream myfile4;
    //myfile3.open ("pos-den-31.txt");
    myfile4.open ("r-den-sigma.txt");
	
    double temp = Temp(accr_rate); //ver unidades de accr_rate

    Vector rp(n_rows,0.0);
	
    
    for (int i = 0; i < n_rows; i++){           // vector para recorrer la posicion medido en metros
    
        rp[i] = r_i + delta_r * i;

 //     myfile << " pos vs densidad" << std::endl;  
 
		double den = density(rp[i]);
		
		myfile3 << rp[i] << "\t" << den <<std::endl;
    
   
		double sigma = conductivity(den,temp); sigma_ph(den,temp); 

		myfile4 << rp[i] << "\t" << den << "\t" << "sigma ph" << "\t" << sigma << std::endl;  
       

	   // myfile3 << rp[i] << "\t" << den << "\t" << sigma << std::endl;  
	   
	   // sigma = (sigma1*sigma2)/(sigma1+sigma2);
	  // myfile4 << den << "\t" << sigma << std::endl;
		
		double alpha = cLight2/(4.0*pi*sigma);  // alpha tiene unidades de cm^2/s

		//bb == es la matriz bb
		
		//  myfile << rp[i] << " " << alpha << std::endl;
		ad[0] = 1; //ver, mejor = 0?            // Defino este valor para que no invente el vector
		
		for (int j = 0; j < n_rows; j++){
            
			if (j == i){
			//bd == vector: diagonal de la matriz aa	
                
                if (i == 0) {
                    bd[i] = 1;
                    bb[i][j] = 1;
                }
                else if (i == n_rows-1){
                    bd[i] = - (1 + (delta_r/rp[i]));
                    bb[i][j] = 0;
                }
                else {
                    bd[i] = -2 *((pow(delta_r*1.e2,2)/(alpha * delta_t))+ 1 + (pow(delta_r*1.e2,2))/(pow(rp[i]*1.e2,2)));// elemento de la diagonal de la matriz AA
                    bb[i][j] = 2 *(-(pow(delta_r*1.e2,2)/(alpha * delta_t))+ 1 + (pow(delta_r*1.e2,2))/(pow(rp[i]*1.e2,2)));// elemento de la diagonal de la matriz BB
   
                    
				}
			}
            
			else if (j == i + 1){
			//cd == vector diagonal superior
                
                if (i == 0) {
                    cd[i] = 0;
                    bb[i][j] = 0;
                }
                
                else {
					cd[i] = 1;
					bb[i][j] = -1;
                }
            }
            else if (j == i - 1){
			//ad == vector diagonal inferior
			
                if (i == n_rows -1) {
                    ad[i] = 1;
                    bb[i][j] = 0;
                }
                else {
                
					ad[i] = 1;
					bb[i][j] = -1;
				}
            }
            else {
                bb[i][j] = 0;
            }
        }
     
    }

    //myfile3.close();
    myfile4.close();
	
}