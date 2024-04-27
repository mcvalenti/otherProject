/*
 * mission.cpp
 *
 *  Created on: 30 mar. 2024
 *      Author: ceci
 */
#include "LDVector.h"
#include "mission.h"
#include "Orbit.h"

/*
SpaceVehicle::SpaceVehicle()
{
	// Default constructor - LEO satellite over the Equator
	long double sv_m[]={6858,0,0,0,7.7102,0,2000};
	LDVector sv0(sv_m, 7);
	StateVector sv1(sv0);
	this->sv=sv1;
	}Ç*/

SpaceVehicle::SpaceVehicle(StateVector& sv){
	this->sv=sv;
}

/*
ostream& operator<<(ostream& os, const SpaceVehicle& sat){
	// Print SpaceVehicle attributes
	os <<"Initial Conditions \n";
    for(unsigned i=0;i<7;i++){
        os<<sat.sv<<"\t";
    }
    
    return os;
}*/


/*

SpaceVehicle::SpaceVehicle(LDVector& init_sv, float mass)
{
	this->init_sv=init_sv;
	this->mass = mass;
	
}*/
