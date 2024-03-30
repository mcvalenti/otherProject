/*
 * mission.cpp
 *
 *  Created on: 30 mar. 2024
 *      Author: ceci
 */
#include "mission.h"
#include "LDVector.h"




SpaceVehicle::SpaceVehicle(LDVector& init_sv, float mass)
{
	this->init_sv=init_sv;
	this->mass = mass;
	
}