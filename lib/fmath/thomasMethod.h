#pragma once
#include "physics.h"

/*Solves system: aa*S(k+1)_j = bb*S^k_j donde
 ad: vector diagonal inferior de aa
 cd vector diagonal superior de aa
 bd diagonal de la matriz aa
 
 bb matriz por el termino en S^k_j es el vector d
  * x es el vector solucion
  */
//void thomasMethod(Matrix bb, Vector cd, Vector ad, Vector bd, Matrix& s, int rows, int ntime);

void thomasMethod(Vector cd, Vector ad, Vector bd, Vector d, Vector x, int rows);