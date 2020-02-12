#pragma once
	
	
	extern const int n_rows;

	extern const double r_f;                   // radio de la estrella de neutrones en metros
	extern const double r_i;                   // radio donde comienza el crust en metros
	extern const double r_o;                  // el tamano de la region de integracion es de 400 m 
	extern const double delta_r;

	extern const int n_time; 
	extern const double t_i;                             // el tiempo se mide en segundos
    extern const double t_f;
	extern const double delta_t;        //CUIDADO! aqui no tengo un DO para el tiempo porque solo estoy probando como escribir AA y BB
    
	extern const double accr_rate;             //[solarMasses yr^-1]