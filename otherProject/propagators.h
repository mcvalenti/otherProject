/*
 * propagators.h
 *
 *  Class and methods to collect init conditions,
 *	perturbartions and to do de integration (RK4).
 *	Structs  of the parameters that the perturbation
 *	models requires.
 *  propagator::addPerturbation will require
 *  perturbation functions to be considered.
 *
 *  LDVector: A class designed to operate with vectors efficiently.
 *
 *  Created on: 28 mar. 2023
 *      Author: ceci
 */

#ifndef PROPAGATORS_H_
#define PROPAGATORS_H_

#include <vector>
#include"LDVector.h"
#include <string>

struct cBody_param{
	long double mu;
};

struct thrust_param {
	long double isp;
	long double thrust;

};

LDVector central_body(void* param, const LDVector& dv);
LDVector thrust(void* param, const LDVector& dv);
void lagrange_coeff_from_true_anomaly(long double r, long double r0, long double h, double delta_nu,
										double& f, double& g, double &fdot, double &gdot);


class propagator
{
	private:
	LDVector total_acc;
	vector<LDVector (*)(void *, const LDVector&)> vec_functions; // vector of:(function pointer, state vector)
	vector<void *> vec_params; // vector of parameters pointers (Structures)
	int total_time;
	long double step;
	LDVector init_sv;
	LDVector RK4();
	LDVector derivatives(LDVector& dv);
	// TO DO POS and VEL modulus method.

    public:
	LDVector current_sv;
	LDVector last_sv;
	long double a,e,i,RAAN,arg_per,nu;
	propagator(LDVector& init_sv, int total_time, long double step);
	propagator(long double a, long double e, long double i,
			long double RAAN, long double arg_per, long double nu);
    void addPerturbation(LDVector (*funcptr)(void *, const LDVector&),void* param);
    void propagate(string const& filename);
    void sv2oe(LDVector& current_sv);
    LDVector sv_from_true_anomaly(LDVector& sv, double delta_nu);

};


#endif /* PROPAGATORS_H_ */
