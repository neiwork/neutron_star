#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>




#include "modelFunctions.h"
#include "initialConditions.h"
#include "matrixComponents.h"
#include "globalParameters.h"

#include "prueba_gauss.h"
#include <fmath/thomasMethod.h>
#include <fmath/physics.h>

int main(int argc, char **argv)
{
    std::ofstream myfile1;
    std::ofstream myfile2;

    myfile1.open ("den-s-inicial.txt");
    myfile2.open ("fr-10-6_gauss.txt");

    
    
	Matrix s;
	setInitialConditions(s, "pos-den-31.txt");
	
	Vector t(n_time,0.0); 
	t[0] = t_i;
	
	for (size_t k = 1; k < n_time; ++k){
		
		double dt = t[k-1]*(t_int - 1.0);
	
		Matrix bb;
		matrixInit(bb,n_rows,n_rows,0.0);   //double bb[n_rows][n_rows];
		
		Vector bd(n_rows,0.0);     //std::vector<double> bd(n_rows);
		Vector cd(n_rows-1,0.0);     //std::vector<double> cd(n_rows-1);
		Vector ad(n_rows,0.0);       //std::vector<double> ad(n_rows);
		Vector d(n_rows,0.0);

		setMatrix(bb, cd, ad, bd, d, dt, s, k-1);

		s[0][k] = 0.0;
			
		Vector x(n_rows,0.0);
		thomasMethod(cd, ad, bd, d, x, n_rows); //, k-1);
		
		
		//Matrix a;
		//matrixInit(a,n_rows,n_rows+1,0.0);
		//prueba(a, cd, ad, bd, d, x, s, k-1);
				
		for (size_t i = 0; i < n_rows; ++i){
			if (i == 0) {
				s[i][k] = 0.0;
			}
			else{
				s[i][k] = x[i];             
			}
		}
		
		
		t[k] = t[k-1]*t_int;
	}

    Vector rp(n_rows,0.0); 
	
	for (size_t j = 1; j < n_time; ++j){
		for (int i = 0; i < n_rows; ++i){
			
			rp[i] = r_i + delta_r * i; 
			
			double den = density(rp[i]);
			
			//myfile2 << (den) << "\t" << log10(t[j]/yr) << "\t" << (s[i][j]) << std::endl;
			myfile2 << (den) << "\t" << (s[i][j]) << std::endl; 
		}
    }
	
    myfile1.close();
    myfile2.close();

	return 0;
}

// ---------------------------------------------------------------------------------------------------------------------------------




   
  //  myfile << "Diagonal principal matriz AA" <<std::endl;
    
  //  for (int i = 0; i < n_rows; i++){           // Escribe en un archivo la diagonal de la matriz AA
  //      myfile << bd[i] << " " << std::endl;
  //     }
    
    
   // myfile << "Diagonal superior matriz AA" <<std::endl;

    
  //  for (int i = 0; i < n_rows-1; i++){           // Escribe en un archivo la diagonal superior de la matriz AA
  //          myfile << cd[i] << " " << std::endl;
       
  //  }
    
  //  myfile << "Diagonal inferior matriz AA" <<std::endl;

    
  //  for (int i = 0; i < n_rows-1; i++){           // Escribe en un archivo la diagonal inferior de la matriz AA
   //       myfile << ad[i+1] << " " << std::endl;
       
  //  }
    
 //  myfile << "Matriz BB" <<std::endl;

    
 //   for (int i = 0; i < n_rows; i++){           // Escribe en un archivo la matriz BB
 //      for (int j = 0; j < n_rows; j++){
 //         myfile << bb[i][j] << " ";
 //      }
 //     myfile << std::endl;
    
    
  //  }