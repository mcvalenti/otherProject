#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <time.h>
#include "propagators.h"
#include "auxiliaries.h"
#include "LDVector.h"
using namespace std;


int main() {

	unsigned size=7;
	long double sv_m[size]={6858,0,0,0,7.7102,0,2000};
	LDVector init_sv(sv_m, 7);
	double Rt=6378;
	double hp=480; // Perigee altitude
	double ha=800; // Apogee altitude
	double rp=Rt+hp; // Perigee radius
	double ra=Rt+ha; // Apogee radius
	long double semimajorAxis=(rp+ra)/2;


	//-----------------------------------------
	// Orbit 1 - Initial - Central Body
	//-----------------------------------------
	long double step=1.0;
	double period0=period(semimajorAxis);
	int total_time=ceil(period0);
	string filename0="output_files/orbit1.csv";

	// First propagation instance
	propagator cbody_prop(init_sv,total_time,step);
	cBody_param cbody;
	cbody.mu=398600.448;
	LDVector cbody_last_sv;

	// Add perturbation
	cbody_prop.addPerturbation(&central_body, &cbody); // args: function and structure
	const clock_t begin_time = clock();
	cbody_prop.propagate(filename0);
	cbody_last_sv=cbody_prop.last_sv;

	std::cout <<cbody_prop.last_sv <<std::endl;
	cbody_prop.sv2oe(cbody_prop.last_sv);
	std::cout <<"Semimajor Axis : "<<cbody_prop.a <<std::endl;

	//-----------------------------------------
	// Orbit 2 - With Continuous thrust
	//-----------------------------------------
	int thrust_time=1000;
	long double thrust_step=0.1;
	string filename1="output_files/orbit2.csv";
	propagator thrust_prop(cbody_last_sv,thrust_time,thrust_step);
	thrust_param tparam;
	tparam.isp = 300;
	tparam.thrust = 10000;
	thrust_prop.addPerturbation(&central_body, &cbody);
	thrust_prop.addPerturbation(&thrust, &tparam);
	thrust_prop.propagate(filename1);

	std::cout <<"True anomaly nu: "<<thrust_prop.a <<std::endl;

	//-----------------------------------------
	// Orbit 3 - Free propagation to Apogee
	//-----------------------------------------
	thrust_prop.sv2oe(thrust_prop.last_sv);
	double period1=period(thrust_prop.a);
	int free_prop_time=ceil(period1);
	long double free_step=1;
	string filename2="output_files/orbit3.csv";
	propagator free_prop(thrust_prop.last_sv,free_prop_time/2,free_step);
	free_prop.addPerturbation(&central_body, &cbody);
	free_prop.propagate(filename2);


	// End of Propagation
	std::cout << "Propagation time: "<<float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;

	std::cout <<free_prop.last_sv <<std::endl;
	free_prop.sv2oe(free_prop.last_sv);
	std::cout <<"Semimajor Axis : "<<free_prop.a <<std::endl;

	// End of Code
	cout<<"End of processing!";

	return 0;
}













