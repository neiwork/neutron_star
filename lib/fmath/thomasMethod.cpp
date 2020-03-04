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
 
 
 
void thomasMethod(Vector cd, Vector ad, Vector bd, Vector d, Vector x, int rows)
{
	

		Vector cn(rows,0.0); //rows -1 
		Vector dn(rows,0.0); 
		//Vector x(rows,0.0);  //vector solucion

        for (int i = 0; i < rows; i++){
			
			if (i == 0) {
				dn[i] = d[i]/bd[i];
				cn[i] = cd[i]/bd[i];
			}
			else if (i == rows-1) {
				dn[i] = (d[i] - ad[i] * dn[i-1])/(bd[i] - ad[i] * cn[i-1]); //cn tiene una dim menos
				}
			else{
				dn[i] = (d[i] - ad[i] * dn[i-1])/(bd[i] - ad[i] * cn[i-1]); 
				cn[i] = cd[i]/(bd[i] - ad[i] * cn[i-1]);  
			}
		}

// ------------------------------------------------------------------------------------
// Aplico el metodo de Thomas para encontrar la solucion en el paso temporal siguiente.
// _____________________________________________________________________________________   
 //   myfile2 << "Vector solucion en el paso temporal siguiente" << std::endl; 
 
		for (int i = rows-1; i >= 0; i--){
			
			if (i == rows-1){
				x[i] = dn[i];             
			}
			else{
			   x[i] = dn[i] - cn[i] * x[i+1];       // aqui calcula el vector solucion
			}

		}
	//}
}
