/*
 * propagators.cpp
 *
 *  Created on: 28 mar. 2023
 *      Author: ceci
 */

#include <iostream>
#include <cmath>
#include <vector>
#include "propagators.h"
#include"DVector.h"
using namespace std;


void central_body(void* param){
	cBody_param* cbody=(cBody_param*)param;
	double r, r3, coef;
	r=sqrt(cbody->sv[0]*cbody->sv[0]+cbody->sv[1]*cbody->sv[1]+cbody->sv[2]*cbody->sv[2]);
	r3=r*r*r;
	cout<<"Central body radius:"<<r3<<endl;
	coef=-cbody->mu/r3;
	cbody->acc[0]=cbody->sv[3];
	cbody->acc[1]=cbody->sv[4];
	cbody->acc[2]=cbody->sv[5];
	cbody->acc[3]=cbody->sv[0]*coef;
	cbody->acc[4]=cbody->sv[1]*coef;
	cbody->acc[5]=cbody->sv[2]*coef;
}

DVector keplerian(DVector& dv, double mu){
	unsigned size=6;
	DVector acc_kep(size);
	double r, r3, coef;
	
	r=sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]);
	r3=r*r*r;
	coef=-mu/r3;
	acc_kep[0]=dv[3];
	acc_kep[1]=dv[4];
	acc_kep[2]=dv[5];
	acc_kep[3]=dv[0]*coef;
	acc_kep[4]=dv[1]*coef;
	acc_kep[5]=dv[2]*coef;
	
	return acc_kep;
}

DVector thrust(DVector& dv, double mass, double T, double isp){
/*
	    Computes thrust forces and mass consumption
	    input: stateVector_mass [array, dim:7] - The state vector considers
	            the satellite mass
	    output : acceleration vector for each components [km2/s] and
	            mass consumption
*/
	const double g0=9.81; // [km/s2]
	double vel[3]={dv[3],dv[4],dv[5]}; 
	double vel_norm;
	double coeff;
	DVector acc_thrust;
	
    vel_norm=sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);
    coeff=T/(mass*vel_norm*1000); // Thrust aligned with velocities [km/s2]
	acc_thrust[0]=0;	// np.array([stateVector[3]*coeff,stateVector[4]*coeff,stateVector[5]*coeff,-thrust/(isp*g0)])										
	acc_thrust[1]=0;
	acc_thrust[2]=0;
	acc_thrust[3]=dv[3]*coeff;									
	acc_thrust[4]=dv[4]*coeff;
	acc_thrust[5]=dv[5]*coeff;
	acc_thrust[6]=-T/(isp*g0);
	return acc_thrust;
}



DVector _deriv(DVector& dv){

	unsigned size=6;
	DVector acc0(size);
	double r, r3, coef;
	double mu=398600.448;
	r=sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]);
	r3=r*r*r;
	coef=-mu/r3;
	acc0[0]=dv[3];
	acc0[1]=dv[4];
	acc0[2]=dv[5];
	acc0[3]=dv[0]*coef;
	acc0[4]=dv[1]*coef;
	acc0[5]=dv[2]*coef;
	
	return acc0;
}



DVector RK4(double h, DVector& dv){
	
	// RK4
	DVector k1,k2,k3,k4;
	DVector y1,y2,y3;
	DVector dv1;
	
	k1=_deriv(dv)*h;
	y1=dv+k1*0.5;
	k2=_deriv(y1)*h;
	y2=dv+k2*0.5;
	k3=_deriv(y2)*h;
	y3=dv+k3;
	k4=_deriv(y3);
	dv1=dv+(k1+k2*2+k3*2+k4)*(1./6);
	
	return dv1;
}

