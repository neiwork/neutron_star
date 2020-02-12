#pragma once
#include "physics.h"

/*Solves system: aa*S(k+1)_j = bb*S^k_j donde
 ad: vector diagonal inferior de aa
 cd vector diagonal superior de aa
 bd diagonal de la matriz aa
 
 bb matriz por el termino en S^k_j
  */
void thomasMethod(Matrix bb, Vector cd, Vector ad, Vector bd, Matrix& s, int rows, int ntime);