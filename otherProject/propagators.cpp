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
#include "auxiliaries.h"
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


propagator::propagator(long double a, long double e, long double i, // @suppress("Class members should be properly initialized")
		long double RAAN, long double arg_per, long double nu){
	this->a=a;
	this->e=e;
	this->i=i;
	this->RAAN=RAAN;
	this->arg_per=arg_per;
	this->nu=nu;
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


void propagator::sv2oe(LDVector& init_sv){
	/*
	 -----------------------------------------------------
	 Curtis - Orbit Mechanical for Engineering Students
	 -----------------------------------------------------
	 Gets state vector (sv) values and sets
	 the keplerian orbital elements.
	 to know: a,e,i,RAAN,arg_per,nu (true anomaly)
	 Cartesian
	 pos: position vector
	 vel: velocity vector
	 Keplerian
	 a: semimajor axis [m]
	 e: eccentricity
	 i: inclination [deg]
	 RAAN: Right Ascention of Ascending Node [deg]
	 arg_per: argument of perigee [deg]
	 nu: true anomaly [deg]
	*/

	//long double mu=398600.448;
	//long double earth_radius=6378.0; //[km]
	long double GM=398600.4405; //[km3/s2]
	LDVector pos(3);
	LDVector vel(3);
	long double pos_mod;
	long double vel_mod;
	LDVector h(3);
	long double h_mod;
	LDVector h_norm(3);
	long double vr;
	LDVector e_vec(3);
	long double M_anom;
	long double E_anom;
	LDVector aux1(3);
	LDVector aux2(3);
	LDVector aux3(3);

	// Position vector
	pos[0]=init_sv[0];
	pos[1]=init_sv[1];
	pos[2]=init_sv[2];
	// Velocity vector
	vel[0]=init_sv[3];
	vel[1]=init_sv[4];
	vel[2]=init_sv[5];

	pos_mod=sqrtl(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	vel_mod=sqrtl(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);

	// Angular momentum
	h=cross(pos,vel);
	h_mod=sqrtl(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]);
	h_norm[0]=h[0]*(1.0/h_mod);
	h_norm[1]=h[1]*(1.0/h_mod);
	h_norm[2]=h[2]*(1.0/h_mod);

	// inclination
	this->i=atanl(sqrtl(h[0]*h[0]+h[1]*h[1])/h[2]);
	if (this->i < 0.0){
		this->i=M_PI+this->i;
	}
	// RAAN
	if (h[1]!=0){
		this->RAAN=atanl(-h[0]/h[1]);
		// A non-zero value (true) if the sign of x is negative; and zero (false) otherwise.
		if (signbit(-h[1])>0){
			this->RAAN=this->RAAN+M_PI;
		}
	}else{
		this->RAAN=0.0;
	}
	// Radial Velocity
	vr=pos[0]*vel[0]+pos[1]*vel[1]+pos[2]*vel[2];
	vr=vr/pos_mod;

	// eccentricity
	e_vec=(pos*(vel_mod*vel_mod-(GM/pos_mod))+vel*pos_mod*(-1)*vr)*(1/GM);
	this->e=sqrtl(e_vec[0]*e_vec[0]+e_vec[1]*e_vec[1]+e_vec[2]*e_vec[2]);

	// semimajor axis
	this->a = h_mod*h_mod/(GM*(1-this->e*this->e));

	// Mean anomaly

	E_anom=atanl(((vr*pos_mod)/(this->a*this->a*this->a))/(1-pos_mod/this->a));
	if (signbit(1-pos_mod/this->a)){E_anom=E_anom+M_PI;};
	M_anom=E_anom-this->e*sinl(E_anom);


	// Argument of Perigee
		this->arg_per=NAN; // TO DO!, not necessary now.
	// True Anomaly
	/*
	 If e*pos>0, then nu lies in I or IV quadrant.
	 If e*pos<0, then nu lies in II or III quadrant.
	 */

	aux1=e_vec*this->e;
	aux2=pos*pos_mod;
	if (vr>=0){
		this->nu=acosl(dot3(aux1,aux2));
	}
	else{
		this->nu=2*M_PI-acosl(dot3(aux1,aux2));
	}
}

