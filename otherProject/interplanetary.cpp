/*
 * interplanetary.cpp
 *  Plantery departure: concept of Hyperbolic Excess Velocity = V_inf
 *
 *
 *
 *
 *  Created on: Jul 13, 2024
 *      Author: macec
 */
#include <iostream>
#include <cmath>
#include "interplanetary.h"
#include "orbit_dynamics.h"

using namespace std;

void run_interplanetary(){
	// variables
	double delta_v, v_elliptic, delta_v_final, v_park;
	double vel_escape, hyper_escape;
	// Constants
	double AU = 150e6; //  km 149597870.7;
	double mu_earth = 398604; // km3/s2
	double r_earth=6378;
	double h_parking=300; // km
	double r_parking=h_parking+r_earth; // km
	double mu_sun=1.33e11; // km3/s2

	// Delta V Heliocentric Homann Transfer
	double R1=AU; // km
	double R2=39.7e6; // km Decided by the project

	std::cout<<"================================================="<<std::endl;
	delta_v=delta_V_homann_heliocentric(mu_sun, R1, R2);
	std::cout<<"Delta V Homann Heliocentric - DV-H: "<<delta_v<<std::endl;
	std::cout<<"================================================="<<std::endl;
	v_elliptic=elliptic_velocity((R1+R2)/2, R1, mu_sun);
	std::cout<<"S/C velocity in Elliptic Heliocentric - DV1: "<<v_elliptic<<std::endl;
	std::cout<<"================================================="<<std::endl;
	v_park=parking_v(mu_earth, r_parking);
	delta_v_final=delta_v_to_hyperbola(v_park, delta_v);
	std::cout<<"V Parking - v_park: "<<v_park<<std::endl;
	std::cout<<"From Parking to Hiperbolic excess - Total DV: "<<delta_v_final<<std::endl;
	std::cout<<"================================================="<<std::endl;
	vel_escape= escape_vel(mu_earth, r_parking);
	std::cout<<"V escape from parking orbit - v_escape: "<<vel_escape<<std::endl;
	hyper_escape=hyperbolic_escape_velocity(vel_escape, delta_v);
	std::cout<<"Total departure velocity - sqrt(v_inf2+v_esc2)- Total V: "<<hyper_escape<<std::endl;
	std::cout<<"================================================="<<std::endl;

}

double delta_V_homann_heliocentric(double mu_sun, double R1, double R2){
	/* Planetary departure
	* Heliocentric velocity for a SC departing on a Homann trajectory,
	* to a point farther from the Sun.
	* R1: Distance of the departure planet (planet 1) to the Sun
	* R2: Distance of the arriving point/planet (planet 2) to the Sun
	*/
	double delta_v;
	delta_v=sqrt(mu_sun/R1)*(sqrt(2*R2/(R1+R2))-1);

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





// Extra data
//double falconH=11.2; // km/s
//double distance_mars_sun = 1.524*AU;
//double distance_venus_sun = 1.077504e8; // 1.077504 x 10^8 km
//double distance_mercury_sun = 5.833821e7; // 5.833821 x 10^7 km
//double mass_sun = 1.989e30; // Masa del Sol en kg
//double G = 6.67430e-11;  // m^3 kg^-1 s^-2
