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

#include <fmath/thomasMethod.h>
#include <fmath/physics.h>

int main(int argc, char **argv)
{
    std::ofstream myfile1;
    std::ofstream myfile2;

    myfile1.open ("den-s-inicial.txt");
    myfile2.open ("fr-10-6.txt");

    
    
	Matrix s;
	setInitialConditions(s, "pos-den-31.txt");
	
	Matrix bb;
	matrixInit(bb,n_rows,n_rows,0.0);   //double bb[n_rows][n_rows];
    
	Vector bd(n_rows,0.0);     //std::vector<double> bd(n_rows);
	Vector cd(n_rows-1,0.0);     //std::vector<double> cd(n_rows-1);
	Vector ad(n_rows,0.0);       //std::vector<double> ad(n_rows);


	setMatrix(bb, cd, ad, bd);

	thomasMethod(bb, cd, ad, bd, s, n_rows, n_time);


    Vector rp(n_rows,0.0); 
	
    for (int i = 0; i < n_rows; i++){
		
		rp[i] = r_i + delta_r * i; 
		
		double den;

		if (r_f - rp[i] <= 100) {

			den = density1(r_f - rp[i]);
		}

		else {
		//      myfile << r_f - rp[i] << " " << den <<std::endl;

			den = density3(r_f - rp[i]);

		}
   
		
		myfile2 << den << " " << s[i][n_time-1] << std::endl;
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