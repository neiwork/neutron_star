#pragma once


#include <string>
//#include <fparticle\Particle.h>
#include <iostream>
#include <fstream>
#include <fmath/physics.h>




/* Read is a function that from archive obtained logRho and y and generates the vectors x and y. */ 
void read (const std::string& archive, Vector& x, Vector& y);