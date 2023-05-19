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
#include "LDVector.h"
using namespace std;


LDVector thrust(void *param, const LDVector& dv){
/*
		Computes thrust forces and mass consumption
		input: thrust_param structure
		output: thrust acceleration LDVector
*/
    thrust_param* thr_p = (thrust_param*)param;
	const long double g0=9.81; // [km/s2]
	long double vel_norm;
	long double coeff;
	LDVector acc_thrust(7);

    vel_norm=sqrtl(dv[3]*dv[3]+dv[4]*dv[4]+dv[5]*dv[5]);
    coeff=thr_p->thrust/(dv[6]*vel_norm*1000); // Thrust aligned with velocities [km/s2]
    acc_thrust[0]=0.0;
    acc_thrust[1]=0.0;
    acc_thrust[2]=0.0;
    acc_thrust[3]= dv[3]*coeff;
    acc_thrust[4]= dv[4]*coeff;
    acc_thrust[5]= dv[5]*coeff;
    acc_thrust[6]=(long double)(-thr_p->thrust)/(thr_p->isp*g0);

	return acc_thrust;

}

LDVector central_body(void* param, const LDVector& dv){
	/*
	 Computes Central Body accelerations.
	 input:
	 ------
	 mu param - long double
	 dv state vector - LDVector
	 output:
	 ------
	 cbody_acc acceleration vector - LDVector
	 */

	cBody_param* cbody=(cBody_param*)param;
	LDVector cbody_acc(7);
	long double r, r3, coef;
	r=sqrtl(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]);
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






propagator::propagator(LDVector& init_sv, int total_time, long double step){
	this->init_sv=init_sv;
	this->total_time=total_time;
	this->step=step;
}

LDVector propagator::derivatives(LDVector& dv){
	int i=0;
	long double sv[7]={dv[3], dv[4], dv[5], 0.0, 0.0, 0.0, 0.0};
    LDVector f_deriv(sv, 7);
    LDVector accelerations(7);

	// Luego las ejecuta una a una
    for  (auto f : this->vec_functions){
    	accelerations+=f(this->vec_params[i], dv);
        i++;
    }
    f_deriv=f_deriv+accelerations;
    return f_deriv;

}


void propagator::addPerturbation(LDVector (*funcptr)(void *, const LDVector&),void* param)
{
    vec_functions.push_back(funcptr);
    vec_params.push_back(param);
 }
void propagator::propagate()
{
	// Propagation with RK4
	int i;
	LDVector dv1;
	vector<LDVector> sv_list;
	ofstream outputFile;
	cout<<"To propagate!"<<endl;
	for (i=0; i<=this->total_time; i++){
		this->init_sv=this->RK4();
		sv_list.push_back(this->init_sv);
	}
	outputFile.open("output.csv");
	for (LDVector sv : sv_list){
		outputFile <<sv<<endl;
	}
	outputFile.close();
}


LDVector propagator::RK4(){

	// RK4
	LDVector k1,k2,k3,k4;
	LDVector y1,y2,y3;
	LDVector dv1;
	LDVector dv = this->init_sv;

	k1=derivatives(dv)*this->step;
	y1=dv+k1*0.5;
	k2=derivatives(y1)*this->step;
	y2=dv+k2*0.5;
	k3=derivatives(y2)*this->step;
	y3=dv+k3;
	k4=derivatives(y3);
	dv1=dv+(k1+k2*2+k3*2+k4)*((long double)(1.0L/6.0L))*this->step;

	this->init_sv = dv1;
	return this->init_sv;
}

