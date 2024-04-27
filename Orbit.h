/*
 * StateVector.h
 *
 *  Defines Mision State Vector
 *  It is a particular LDVector, with coordinates, velocities and mass
 *  Dim = 6 or 7(if mass included)
 *
 *  Created on: 4 abr. 2024
 *      Author: ceci
 */

#ifndef ORBIT_H_
#define ORBIT_H_

#include"LDVector.h"

class StateVector{
	/*
	 *  LDVector with Cartesian Coordinates of the orbit and mass
	 *
	*/
	private:
	LDVector sv0;
	public:
	StateVector()=default;
	StateVector(LDVector& sv); // TO DO check dimensions
	void toOrbitalElements(LDVector& sv, LDVector& oe);
	friend ostream& operator<<(ostream& os, const StateVector& sv);

};

class OrbitalElements{
	/*
	 *  LDVector with Keplerian elements of the orbit and mass
	 *
	*/
	public:
	long double a, e, i, Omega, arg_per, M; //TO DO, describe every element
	long double mass;
	OrbitalElements();
	OrbitalElements(LDVector& soe); // TO DO check dimensions
	friend ostream& operator<<(ostream& os, const StateVector& sv);

};


#endif /* ORBIT_H_ */
