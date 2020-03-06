#include "prueba_gauss.h"

#include "globalParameters.h"

#include <fmath/elimiGaussiana.h>

void prueba(Matrix& a, Vector& cd, Vector& ad, Vector& bd, Vector& d, Vector& x, Matrix& s, int k)

{
	
	for (int i = 0; i < n_rows; i++){
		
		a[i][i] = bd[i];
		a[i][n_rows] = d[i];
		
		if (i == n_rows-1){
			//a[i][i+1] = cd[i];
			a[i][i-1] = ad[i];
			
		}
		else if (i == n_rows-1){
			//a[i][i+1] = cd[i];
			a[i][i-1] = ad[i];
		}
	   else{
			a[i][i+1] = cd[i];
			//a[i][i-1] = ad[i];
	   }
	}
	
	
	elimiGaussiana(n_rows, a, x);
	
	for (int i = 0; i < n_rows; ++i){

		s[i][k+1] = x[i];             
		
	}
	
}