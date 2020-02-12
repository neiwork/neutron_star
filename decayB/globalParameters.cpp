#include "globalParameters.h"

#include <fmath/physics.h>
      

	const int n_rows = 200;

	const double r_f = 10.6e3;                   // radio de la estrella de neutrones en metros
	const double r_i = 9.77e3;                   // radio donde comienza el crust en metros
	const double r_o = 10174.6;                  // el tamano de la region de integracion es de 400 m 
	const double delta_r = (r_f - r_i)/(n_rows-1);

	const int n_time = 1e3; 
	const double t_i = 0.0;                             // el tiempo se mide en segundos
    const double t_f = 1e6*yr;
	const double delta_t = (t_f - t_i)/(n_time);        //CUIDADO! aqui no tengo un DO para el tiempo porque solo estoy probando como escribir AA y BB
    
	const double accr_rate = 3.e-10;             //[solarMasses yr^-1]

