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
	 Modifies cBody.acc acceleration parameter
	 */
	cBody_param* cbody=(cBody_param*)param;
	double r, r3, coef;
	r=sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]);
	r3=r*r*r;
	coef=-cbody->mu/r3;
	cbody->acc[0]=0;
	cbody->acc[1]=0;
	cbody->acc[2]=0;
	cbody->acc[3]=dv[0]*coef;
	cbody->acc[4]=dv[1]*coef;
	cbody->acc[5]=dv[2]*coef;
	cbody->acc[6]=0;

	return cbody->acc;
}


DVector thrust(void *param, const DVector& dv){
/*
		Computes thrust forces and mass consumption
		input: thrust_param structure
		thrust parameteres acceleration results modified
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

	// Vector of functions manipulation
	DVector acc_cb(7);
	DVector acc_thr(7);
	cBody_param cbody;

	cbody.acc=acc_cb;
	cbody.mu=398600.448;
	cbody.sv=dv;

	thrust_param tparam;
	tparam.isp = 1200;
	tparam.thrust = 10000;
	tparam.acc_thr = acc_thr;
	tparam.sv=dv;


	// Primero enlista las funciones y sus parametros
	this->addPerturbation(&central_body, &cbody); // agrega la perturbacion a la lista
	//this->addPerturbation(&thrust, &tparam);
	// Luego las ejecuta una a una
    for  (auto f : this->vec_functions){
    	accelerations+=f(this->vec_params[i], this->init_sv); // Cada Vec_param debe contener una acc que resulta de la perturbacion
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
	outputFile.open("output.csv");
	for (DVector sv : sv_list){
		outputFile << sv<<endl;
	}
	outputFile.close();
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

