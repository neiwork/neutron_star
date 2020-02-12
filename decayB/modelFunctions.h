#pragma once

#include <fmath/physics.h>

/*Declaracion de la funcion para la temperatura*/ 
double Temp(double Md);


/*Declaracion de la funcion density1*/
double density1(double pos);

/*Declaracion de la funcion density2*/ //por ahora no se esta usando

/*Declaracion de la funcion density3*/
double density3(double pos);

double conductivity(double den,double temp);

/* Declaracion de la funcion Conductivity due to electron-phonon scattering, denoted sigma_ph*/ 
double sigma_ph(double den, double temp);


/*Declaracion de la funcion Conductivity due to electron scattering on impurities, denoted sigma_imp*/
double sigma_imp(double den);


/*Declaracion de la funcion conductividad de transicion*/
double sigma_trans(double den);

