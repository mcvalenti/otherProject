// Mission's Setting and Scenario selection
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <time.h>
#include "propagators.h"
#include "auxiliaries.h"
#include "LDVector.h"
#include "my_tests.h"
#include "mission.h"
#include "Orbit.h"

// Globals
double Rt=6378;  // Earth radius

int main() {

	// TESTS
	std::cout << "--------TESTs-------------\n";
	//run_tests();
	std::cout << "-------- END OF TESTs-------------\n";


	unsigned size=7;
	long double sv_m[size]={6858,0,0,0,7.7102,0,2000};
	LDVector init_sv(sv_m, 7);
	std::cout<<init_sv<<std::endl;

	long double x1=7600;
	long double x2=7900;
	long double x3=500;
	long double x4=7.1;
	long double x5=0.5;
	long double x6=0.3;
	long double mass=0;
	std::string elements="cartessian";

	StateVector sv0(x1, x2, x3, x4, x5, x6, mass, elements);
	std::cout<<sv0<<std::endl;

	Continuar con la herencia de la clase LDVector.
	Armar los constructores de StateVector a partir de LDVector.

	return 0;
}













