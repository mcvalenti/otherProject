/*
 * StateVector.cpp
 *
 *  Created on: 4 abr. 2024
 *      Author: ceci
 */

#include "Orbit.h"
#include "LDVector.h"

StateVector::StateVector(LDVector& sv){
	this->sv0=sv;
};

ostream& operator<<(ostream& os, const StateVector& sv){
	// Print StateVector attributes
	std::cout<<sv<<std::endl;
    return os;
}

OrbitalElements::OrbitalElements(LDVector& oe){
	this->a=oe[0];
	this->e=oe[1];
	this->i=oe[2];
	this->Omega=oe[3];
	this->arg_per=oe[4];
	this->M=oe[5];
	this->mass=oe[6];
};


