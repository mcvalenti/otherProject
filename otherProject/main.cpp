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

using namespace std;

// Globals
double Rt=6378;  // Earth radius

int main() {

	// TESTS
	cout << "--------TESTs-------------\n";
	run_tests();
	cout << "-------- END OF TESTs-------------\n";

	unsigned size=7;
	float mass;
	long double sv_m[size]={6858,0,0,0,7.7102,0,2000};
	LDVector init_sv(sv_m, 7);
	mass=init_sv[6];
	SpaceVehicle my_sat(init_sv, mass);
	StateVector my_sv(init_sv);


	std::cout<<my_sat<<std::endl;
	std::cout<<my_sv.x<<std::endl;

	return 0;
}













