#include "modelFunctions.h"
#include "globalParameters.h"

#include <string>
#include <iostream>
#include <fstream>



void setInitialConditions(Matrix& s, const std::string& filename)
//void setInitialConditions(Matrix& s)

{

	std::ofstream file;
	file.open(filename.c_str(),std::ios::out);

	 matrixInit(s,n_rows,n_time,0.0);// s[n_rows][n_time];              // declaracion de la matriz solucion

     Vector rp(n_rows,0.0);
 

//  Condiciones iniciales
   
	for (int i = 0; i < n_rows; i++){          

		rp[i] = r_i + delta_r * i;  

		double s0 = 0.0; 

		if (rp[i] > r_o) {
			s0  = P2((rp[i]-r_o)/(r_f-r_o));//  1.0-P2((r_f-rp[i])/(r_f-r_o)); //pow(((rp[i] - r_o)/(r_f-r_o)),2);
		}
		
		s[i][0] = s0;

		double den = density(rp[i]);
		
		file  << den << "\t" << s0 << std::endl;
	}
 
	
	file.close();
	
}