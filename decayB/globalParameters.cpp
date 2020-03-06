#include "globalParameters.h"

#include <fmath/physics.h>
      

	const int n_rows = 10;

	const double r_f = 10.6e3;                   // radio de la estrella de neutrones en metros
	const double r_i = 9.77e3;                   // radio donde comienza el crust en metros
	const double r_o = 10329.0; //10174.6;                  // el tamano de la region de integracion es de 400 m 
	const double delta_r = (r_f - r_i)/(n_rows-1);

	const int n_time = 10; 
	const double t_i = 1.0e-3;                             // el tiempo se mide en segundos
    const double t_f = 1e6*yr;
	//const double delta_t = (t_f - t_i)/(n_time);        //CUIDADO! aqui no tengo un DO para el tiempo porque solo estoy probando como escribir AA y BB
    const double t_int = pow(t_f / t_i, (1.0 / n_time));
	
	const double accr_rate = 3.e-10;             //[solarMasses yr^-1]

	const double Q = 0.01;