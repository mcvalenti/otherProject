/*
 * mission.h
 *
 *  Created on: 30 mar. 2024
 *      Author: ceci
 */

#ifndef MISSION_H_
#define MISSION_H_

#include "LDVector.h"
#include "Orbit.h"



class SpaceVehicle
{
	private:
		//incorporate StateVector to ARCH
		StateVector sv;
	public:
		LDVector init_sv;
		float mass;
		SpaceVehicle();
		SpaceVehicle(LDVector& init_sv, float mass);
		SpaceVehicle(StateVector& sv);
		//~SpaceVehicle();
		friend ostream& operator<<(ostream& os, const SpaceVehicle& sat);
		
};		


#endif /* MISSION_H_ */
