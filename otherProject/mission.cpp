/*
 * mission.cpp
 *
 *  Created on: 30 mar. 2024
 *      Author: ceci
 */
#include "mission.h"
#include "LDVector.h"

SpaceVehicle::SpaceVehicle()
{
	// Default constructor - LEO satellite over the Equator
	long double sv_m[]={6858,0,0,0,7.7102,0,2000};
	LDVector init_sv(sv_m, 7);
	this->init_sv=init_sv;
	this->mass=init_sv[6];
	}

SpaceVehicle::SpaceVehicle(LDVector& init_sv, float mass){
	this->init_sv=init_sv;
	this->mass=init_sv[6];
}

SpaceVehicle::SpaceVehicle(StateVector& sv){
	this->sv=sv;
	this->mass=sv.mass;
}

ostream& operator<<(ostream& os, const SpaceVehicle& sat){
	// Print SpaceVehicle attributes
	os <<"Initial Conditions \n";
    for(unsigned i=0;i<7;i++){
        os<<sat.init_sv[i]<<"\t";
    }
    
    return os;
}


/*

SpaceVehicle::SpaceVehicle(LDVector& init_sv, float mass)
{
	this->init_sv=init_sv;
	this->mass = mass;
	
}*/
