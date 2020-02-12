#include "thomasMethod.h"

 // __________________________________________________________________________________
  // __________________________________________________________________________________
 // Redefino los vectores para poder aplicar el metodo de Thomas
 // ----------------------------------------------------------------------------------
 
 //   myfile << "Diagonal superior redefinida" << std::endl;   
 
 //ad vector diagonal inferior
 //cd vector diagonal superior
 //bd diagonal de la matriz aa
 
 //bb matriz por el termino en S^k_j
 //d(i) Vector que da el termino inhomogeneo: bb[i][j] * s[j][k]
 
 
 
void thomasMethod(Matrix bb, Vector cd, Vector ad, Vector bd, Matrix& s, int rows, int ntime)
{
	
	for (int k = 0; k < ntime - 1; k++){

		Vector d(rows,0.0);
		
//____________________________________________________________________________
// Calculo el vector que da el termino inhomogeneo (vector dd)
//d(i) Vector que da el termino inhomogeneo: bb[i][j] * s[j][k]
// ___________________________________________________________________________
//   myfile << " Vector termino inhomogeneo " <<std::endl;

		for (int i = 0; i < rows; i++){
			
			double sum = 0.0;
			
			for (int j = 0; j < rows; j++){
				sum = sum + bb[i][j] * s[j][k];
            }
			d[i] = sum;    
		}


		Vector cn(rows-1,0.0);  //std::vector<double> cn(n_rows-1);
		Vector dn(rows,0.0); //ver si es rows-1 //std::vector<double> dn(n_rows-1);
		
		for (int i = 0; i < rows-1; i++){
        
			if (i == 0) {
				cn[i] = cd[i]/bd[i];
			}
			else{
				cn[i] = cd[i]/(bd[i] - ad[i] * cn[i-1]);  
			}
		}
    
   //  myfile << "Vector inhomogeneo redefinido" << std::endl;
 
		for (int i = 0; i < rows; i++){ 
			
			if (i == 0) {
				dn[i] = d[i]/bd[i];
	   //     	myfile << i <<" " << dn[i]<<std::endl;  
			}

			else{
				dn[i] = (d[i] - ad[i] * dn[i-1])/(bd[i] - ad[i] * cn[i-1]); 
		   //     myfile << i <<" " << dn[i]<<std::endl;  
			}

		}

// ------------------------------------------------------------------------------------
// Aplico el metodo de Thomas para encontrar la solucion en el paso temporal siguiente.
// _____________________________________________________________________________________   
 //   myfile2 << "Vector solucion en el paso temporal siguiente" << std::endl; 
 
		for (int i = rows-1; i >= 0; i--){
			
			if (i == rows-1){
				s[i][k+1] = dn[i];             
          //  	myfile2 << k+1 << " " << i << " " << rp[i] << " " << s[i][k+1] << std::endl;

			}
		else{
           s[i][k+1] = dn[i] - cn[i] * s[i+1][k+1];       // aqui calcula el vector solucion
        //   myfile2 << k+1 << " " << i << " " << rp[i] << " " << s[i][k+1] << std::endl;
		}

		}
	}
}