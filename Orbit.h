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

class StateVector: public LDVector {
	/*
	 *  LDVector with Cartesian Coordinates of the orbit and mass
	 *  TO DO check dimensions
	*/
	private:
	long double x1, x2, x3, x4, x5, x6;
	long double mass;
	std::string elements;
	public:
	StateVector();
	StateVector(long double x, long double y, long double z, long double vx, long double vy, long double vz,
				long double mass, std::string elements);
	~StateVector(){std::cout<<"StateVector destructor \n";};
	//void toOrbitalElements(LDVector& sv, LDVector& oe);
	//friend ostream& operator<<(ostream& os, const StateVector& sv);
	friend ostream& operator<<(ostream& os, const StateVector& vec);

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
