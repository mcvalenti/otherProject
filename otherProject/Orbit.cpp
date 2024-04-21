/*
 * StateVector.cpp
 *
 *  Created on: 4 abr. 2024
 *      Author: ceci
 */

#include "Orbit.h"

#include "LDVector.h"

StateVector::StateVector(LDVector& sv){
	this->x=sv[0];
	this->y=sv[1];
	this->z=sv[2];
	this->vx=sv[3];
	this->vy=sv[4];
	this->vz=sv[5];
	this->mass=sv[6];
};

ostream& operator<<(ostream& os, const StateVector& sv){
	// Print StateVector attributes
	os<<sv.x<<"\t";
	os<<sv.y<<"\t";
	os<<sv.z<<"\t";
	os<<sv.vx<<"\t";
	os<<sv.vy<<"\t";
	os<<sv.vz<<"\t";
	os<<sv.mass<<"\t";

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


