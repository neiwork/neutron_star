#include "modelFunctions.h"
#include "globalParameters.h"

#include <string>
#include <iostream>
#include <fstream>



void setInitialConditions(Matrix& s, const std::string& filename)
//void setInitialConditions(Matrix& s)

{
	//int n_rows = 200;
 
	
	std::ofstream file;
	file.open(filename.c_str(),std::ios::out);
    
    // Defino el vector solucion que da las condiciones iniciales

     //double r_f = 10.6e3;                   // radio de la estrella de neutrones en metros
     //double r_i = 9.77e3;                   // radio donde comienza el crust en metros
     //double r_o = 10174.6;                                          // el tamano de la region de integracion es de 400 m 
    
	 //int n_time = 1e3;  
    
     //double s[n_rows][n_time];              // declaracion de la matriz solucion
	 matrixInit(s,n_rows,n_time,0.0);// s[n_rows][n_time];              // declaracion de la matriz solucion
     //std::vector<double> rp(n_rows);        // declaracion del vector posicion
     Vector rp(n_rows,0.0);
     
    //double delta_r = (r_f - r_i)/(n_rows-1);
    double den;

//  Condiciones iniciales
   
	for (int i = 0; i < n_rows; i++){          

		rp[i] = r_i + delta_r * i;  

		


		if (rp[i] >= r_o) {
			s[i][0] = pow(((rp[i] - r_o)/(r_f-r_o)),2);
		}
		else {
			s[i][0] = 0;
		}
		
		file  << den << "\t" << s[i][0] << std::endl;
	}

	
	file.close();
	
}