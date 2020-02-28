#pragma once

#include <fmath/physics.h>

/*Declaracion de la funcion para la temperatura*/ 
double Temp(double Md);


/* parametrizacion de xx*/ 
double density(double r);

/*Declaracion de la funcion density1*/
//double density1(double pos);

/*Declaracion de la funcion density2*/ //por ahora no se esta usando

/*Declaracion de la funcion density3*/
//double density3(double pos);


double conductivity(double den,double temp);

/* Declaracion de la funcion Conductivity due to electron-phonon scattering, denoted sigma_ph
 * dispersion por fonones: Geppert-Urpin 1994; dispersion por iones: Potekhin 1999*/ 
double sigma_ph(double den, double temp);


/*Declaracion de la funcion Conductivity due to electron scattering on impurities, denoted sigma_imp
 * dispersion por impurezas: Geppert-Urpin 1994*/
double sigma_imp(double den);


/*Declaracion de la funcion conductividad de transicion*/
double sigma_trans(double den);

