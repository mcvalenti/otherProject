/*
 * propagators.cpp
 *
 *  Created on: 28 mar. 2023
 *      Author: ceci
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "propagators.h"
#include "DVector.h"
using namespace std;


DVector central_body(void* param, const DVector& dv){
	/*
	 Computes Central Body accelerations.
	 input:
	 ------
	 mu param - double
	 dv state vector - DVector
	 output:
	 ------
	 cbody_acc acceleration vector - DVector
	 */

	cBody_param* cbody=(cBody_param*)param;
	DVector cbody_acc(7);
	double r, r3, coef;
	r=sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]);
	r3=r*r*r;
	coef=-cbody->mu/r3;
	cbody_acc[0]=0;
	cbody_acc[1]=0;
	cbody_acc[2]=0;
	cbody_acc[3]=dv[0]*coef;
	cbody_acc[4]=dv[1]*coef;
	cbody_acc[5]=dv[2]*coef;
	cbody_acc[6]=0;

	return cbody_acc;
}


DVector thrust(void *param, const DVector& dv){
/*
		Computes thrust forces and mass consumption
		input: thrust_param structure
		output: thrust acceleration DVector
*/
    thrust_param* thr_p = (thrust_param*)param;
	const double g0=9.81; // [km/s2]
	double vel_norm;
	double coeff;
	DVector acc_thrust(7);

    vel_norm=sqrt(dv[3]*dv[3]+dv[4]*dv[4]+dv[5]*dv[5]);
    coeff=thr_p->thrust/(dv[6]*vel_norm*1000); // Thrust aligned with velocities [km/s2]
    acc_thrust[0]=0.0;
    acc_thrust[1]=0.0;
    acc_thrust[2]=0.0;
    acc_thrust[3]= dv[3]*coeff;
    acc_thrust[4]= dv[4]*coeff;
    acc_thrust[5]= dv[5]*coeff;
    acc_thrust[6]=-thr_p->thrust/(thr_p->isp*g0);

	return acc_thrust;

}



propagator::propagator(DVector& init_sv, int total_time, double step){
	this->init_sv=init_sv;
	this->total_time=total_time;
	this->step=step;
}

DVector propagator::derivatives(DVector& dv){
	int i=0;
	double sv[7]={dv[3], dv[4], dv[5], 0.0, 0.0, 0.0, 0.0};
    DVector f_deriv(sv, 7);
    DVector accelerations(7);

	// Luego las ejecuta una a una
    for  (auto f : this->vec_functions){
    	accelerations+=f(this->vec_params[i], this->init_sv);
        i++;
    }
    f_deriv=f_deriv+accelerations;
    return f_deriv;

}


void propagator::addPerturbation(DVector (*funcptr)(void *, const DVector&),void* param)
{
        	vec_functions.push_back(funcptr);
			vec_params.push_back(param);
 }
void propagator::propagate()
{
	// Propagation with RK4
	int i;
	DVector dv1;
	vector<DVector> sv_list;
	ofstream outputFile;
	cout<<"To propagate!"<<endl;
	for (i=0; i<=this->total_time; i++){
		this->init_sv=this->RK4();
		dv1 = this->init_sv;
		sv_list.push_back(dv1);

	}
	/*outputFile.open("output.csv");
	for (DVector sv : sv_list){
		outputFile << sv<<endl;
	}
	outputFile.close();*/
}


DVector propagator::RK4(){

	// RK4
	DVector k1,k2,k3,k4;
	DVector y1,y2,y3;
	DVector dv1;
	DVector dv = this->init_sv;
	k1=derivatives(dv)*this->step;
	y1=dv+k1*0.5;
	k2=derivatives(y1)*this->step;
	y2=dv+k2*0.5;
	k3=derivatives(y2)*this->step;
	y3=dv+k3;
	k4=derivatives(y3);
	dv1=dv+(k1+k2*2+k3*2+k4)*(1./6);

	this->init_sv = dv1;
	return this->init_sv;
}

