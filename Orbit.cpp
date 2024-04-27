/*
 * StateVector.cpp
 *
 *  Created on: 4 abr. 2024
 *      Author: ceci
 */

#include "Orbit.h"

StateVector::StateVector()
{
	std::cout<<"StateVector constructor \n";
	this->x1=0;
	this->x2=0;
	this->x3=0;
	this->x4=0;
	this->x5=0;
	this->x6=0;
	this->mass=0;
	this->elements="notspecified";
}

StateVector::StateVector(long double x, long double y, long double z, long double vx, long double vy, long double vz,
			long double mass, std::string elements){
	std::cout<<"StateVector constructor \n";
	this->x1=x;
	this->x2=y;
	this->x3=z;
	this->x4=vx;
	this->x5=vy;
	this->x6=vz;
	this->mass=mass;
	this->elements=elements;
};


OrbitalElements::OrbitalElements(LDVector& oe){
	this->a=oe[0];
	this->e=oe[1];
	this->i=oe[2];
	this->Omega=oe[3];
	this->arg_per=oe[4];
	this->M=oe[5];
	this->mass=oe[6];
};


ostream& operator<<(ostream& os, const StateVector& vec){
	// Print LDVector
	//os << std::fixed;
	//os << std::setprecision(10);
    for(unsigned i=0;i<vec.getSize();i++){
        os<<vec[i]<<"\t";
    }
    return os;
}
