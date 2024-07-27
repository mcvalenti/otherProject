/*
 * interplanetary.cpp
 *  [Ref 1] Curtis, H. D. (2020). Orbital mechanics for engineering students:
 *	 		Revised Reprint. Butterworth-Heinemann.
 *  Plantery departure: concept of Hyperbolic Excess Velocity = V_inf
 *  Planetary Rendezvous
 *
 *  Created on: Jul 13, 2024
 *      Author: macec
 */
#include <iostream>
#include <cmath>
#include "interplanetary.h"
#include "orbit_dynamics.h"

using namespace std;

// TO DO: compute phase_end according to inter_to_outer // functions

void run_interplanetary(){
	// variables
	double delta_v, v_elliptic, delta_v_final, v_park;
	double vel_escape, hyper_escape, v_eclip_norm;
	double dv_total_departure, C3_orbital;
	// Constants
	double pi = M_PI;
	double AU = 150e6; //  km 149597870.7;
	double mu_earth = 398604; // km3/s2
	double r_earth=6378;
	double h_parking=300; // km
	double r_parking=h_parking+r_earth; // km
	double mu_sun=1.33e11; // km3/s2

	// Delta V Heliocentric Homann Transfer
	double R1=AU-0.8e06; // km
	double R2=39.8e06; // km Decided by the project
	double inc=22.5; // deg inclination wrt Ecliptic plane

	std::cout<<"================================================="<<std::endl;
	delta_v=delta_V_homann_heliocentric(mu_sun, R1, R2);
	std::cout<<"Delta V Homann Heliocentric - DV-H: "<<delta_v<<std::endl;
	std::cout<<"================================================="<<std::endl;
	v_elliptic=elliptic_velocity((R1+R2)/2, R1, mu_sun);
	std::cout<<"S/C velocity in Elliptic Heliocentric - DV1: "<<v_elliptic<<std::endl;
	std::cout<<"================================================="<<std::endl;
	// Projection In the inclination plane
	v_eclip_norm=v_elliptic*tan(inc*pi/180.0);
	std::cout<<"S/C velocity - V_eclip_norm : "<<v_eclip_norm<<std::endl;
	std::cout<<"================================================="<<std::endl;
	// Total DV departure (module)
	dv_total_departure=sqrt(delta_v*delta_v+v_eclip_norm*v_eclip_norm);
	std::cout<<"Total DV departure (module): "<<dv_total_departure<<std::endl;
	C3_orbital=dv_total_departure*dv_total_departure;
	std::cout<<" C3 (orbital): "<<C3_orbital<<std::endl;
	std::cout<<"================================================="<<std::endl;
	v_park=parking_v(mu_earth, r_parking);
	delta_v_final=delta_v_to_hyperbola(v_park, delta_v);
	std::cout<<"V Parking - v_park: "<<v_park<<std::endl;
	std::cout<<"From Parking to Hyperbolic excess - Total DV: "<<delta_v_final<<std::endl;
	std::cout<<"================================================="<<std::endl;
	vel_escape= escape_vel(mu_earth, r_parking);
	std::cout<<"V escape from parking orbit - v_escape: "<<vel_escape<<std::endl;
	hyper_escape=hyperbolic_escape_velocity(vel_escape, delta_v);
	std::cout<<"Total departure velocity - sqrt(v_inf2+v_esc2)- Total V: "<<hyper_escape<<std::endl;
	std::cout<<"================================================="<<std::endl;

}

void run_trip_from_Mars_to_Earth(){
	// Computes the minimun wait time to return from Mars to Earth

	//Constants
	double R_departure = 227.9e06; // Km Mars Distance to Sun
	double R_arrival = 149.6e06;   // km Earth's Distance to Sun
	double mu_sun = 132.71e09;     // Sun mu
	double n_arrival = 0.01720;    // rad/day (2*pi/365.26)
	double n_departure = 0.0091327; // rad/day (2*pi/687.99)

	// Variables
	double t_12, t_12_days, phase_end, t_wait;
	t_12=ellipse_period(mu_sun, R_departure, R_arrival); // sec
	t_12_days = t_12/86400; // days
	phase_end=M_PI-n_arrival*t_12_days; // rad
	t_wait=time_wait(phase_end, n_arrival, n_departure); // days
	std::cout<<"================================================="<<std::endl;
	std::cout<<"Time to wait t_wait: "<<t_wait<<std::endl;
	std::cout<<"================================================="<<std::endl;
}


double synodic_period_from_velocities(double n1, double n2){
	/*
	 * Considering two orbiting bodies, the time required for the phase
	 * angle to return to its initial value value is called the Synodic Period (t_syn).
	 * If n1>n1 is the time required for phi to change  from phi0 to phi0-2*PI.
	 * If n2>n1 is the time required for phi to change from phi0 to phi0+s*PI.
	 * t_syn:[days] Earth days 24hs
	  */
	double t_syn;
	if (n1>n2){
		t_syn=2*M_PI/(n1-n2);
	} else {
		t_syn=2*M_PI/(n2-n1);
	}
	return t_syn;
}

double synodic_period_from_periods(double t1, double t2){
	/*
	 * Synodic Period (t_syn), is the orbital period of planet 2
	 * relative to planet 1.
	 * t1, t2, t_syn:[days] Earth days 24hs
	 */
	double t_syn;
	if (t1>t2){
		t_syn=t1*t2/(t1-t2);
	} else {
		t_syn=t1*t2/(t2-t1);
	}
	return t_syn;
}

void planetary_rendezvous(char inner_outer){
	// Planetary rendezvous & flyby
	//double v_inf, v_sc, v_planet;
	//seguir aca

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

double time_wait(double phase_end, double n1, double n2){
	/*
	 * n1: Arrival Traget Angular velocity
	 * n2: Departure Planet Angular velocity
	 */
	double t_wait=-1;
	int N=0;
	if (n1>n2){
		while (t_wait < 0){
			t_wait=(-2*phase_end-2*M_PI*N)/(n2-n1);
			N=N+1;
		}
	}else{
		while (t_wait < 0){
			t_wait=(-2*phase_end+2*M_PI*N)/(n2-n1);
			N=N+1;
			}
		}
	return t_wait;
};

double return_trip(double mu_center, double R_arrival, double R_departure, double n_arrival, double n_departure){
	/*
	 * Minimum time to return from p1 to p2
	 *
	 */
	double t_wait, t_12, phase_end;

	t_12=ellipse_period(mu_center, R_arrival, R_departure);
	phase_end= M_PI-n_arrival*t_12;
	t_wait=time_wait(phase_end, n_departure, n_arrival);

	return t_wait;
}



// Extra data
//double falconH=11.2; // km/s
//double distance_mars_sun = 1.524*AU;
//double distance_venus_sun = 1.077504e8; // 1.077504 x 10^8 km
//double distance_mercury_sun = 5.833821e7; // 5.833821 x 10^7 km
//double mass_sun = 1.989e30; // Masa del Sol en kg
//double G = 6.67430e-11;  // m^3 kg^-1 s^-2
// Flybys data
// Venus_SOI=6.14e05 // km
// Venus_SOI_bd=100
// Earth_SOI=9.24e05 // km
// Earth_SOI_bd=145
//te=365.26; // Earth's period days
//tm=687.99; // Mars's period days
//period t_earth=365.26;
//period t_venus=225;
//period t_mercury=88;
