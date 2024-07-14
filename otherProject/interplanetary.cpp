/*
 * interplanetary.cpp
 *
 *  Created on: Jul 13, 2024
 *      Author: macec
 */
#include <iostream>
#include <cmath>
#include "interplanetary.h"

using namespace std;

void run_interplanetary(){
	double AU = 149597870.7; // 149,597,870.7 km
	double mass_sun = 1.989e30; // Masa del Sol en kg
	double mu_earth = 3.986e14; // 3.986 x 10^14 m^3 s^-2
	double G = 6.67430e-11;  // m^3 kg^-1 s^-2
	double r_earth=6371;
	double mu_sun;
	double distance_earth_sun = 1.495978707e8; // 1.495978707 x 10^8 km
	double delta_v, delta_v_final, v_park;
	double r_departure=300+r_earth; // km
	//double falconH=11.2; // km/s
	double distance_mars_sun = 1.524*AU;
	//double distance_venus_sun = 1.077504e8; // 1.077504 x 10^8 km
	double distance_mercury_sun = 5.833821e7; // 5.833821 x 10^7 km
	mu_sun=mass_sun*G;
	delta_v=delta_V_homann(mu_sun, distance_earth_sun*1000.0, 35000000*1000.0);
	v_park=parking_v(mu_earth, r_departure*1000.0);
	delta_v_final=delta_v_to_hyperbola(v_park, delta_v);
	std::cout<<"================================================="<<std::endl;
	std::cout<<"Delta V"<<std::endl;
	std::cout<<delta_v<<std::endl;
	std::cout<<"V Parking"<<std::endl;
	std::cout<<v_park<<std::endl;
	std::cout<<"Total V to departure"<<std::endl;
	std::cout<<delta_v_final<<std::endl;
	std::cout<<"================================================="<<std::endl;
}

double delta_V_homann(double mu_center, double R1, double R2){
	double delta_v;
	delta_v=sqrt(mu_center/R1)*(sqrt(2*R2/(R1+R2))-1);

	return delta_v;
}


double parking_v(double mu_center, double r_departure){
	double parking_v;
	parking_v=sqrt(mu_center/r_departure);

	return parking_v;
}

double delta_v_to_hyperbola(double parking_v, double delta_v_inf){
	double delta_v_hyp;
	delta_v_hyp=parking_v*(sqrt(2+(delta_v_inf*delta_v_inf/(parking_v*parking_v)))-1);

	return delta_v_hyp;
}




