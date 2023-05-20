/*
 * propagators.h
 *
 *  Created on: 28 mar. 2023
 *      Author: ceci
 */

#ifndef PROPAGATORS_H_
#define PROPAGATORS_H_

#include <vector>
#include"LDVector.h"

struct cBody_param{
	long double mu;
};

struct thrust_param {
	long double isp;
	long double thrust;

};

LDVector central_body(void* param, const LDVector& dv);
LDVector thrust(void* param, const LDVector& dv);


class propagator
{
	private:
	LDVector total_acc;
	vector<LDVector (*)(void *, const LDVector&)> vec_functions; // vector of function pointers
	vector<void *> vec_params; // vector of parameters pointers (Structures)
	int total_time;
	long double step;
	LDVector init_sv;
	LDVector sv;
	long double a,e,i,RAAN,arg_per,nu;
	LDVector RK4();
	LDVector derivatives(LDVector& dv);

    public:

	propagator(LDVector& init_sv, int total_time, long double step);
	propagator(long double a, long double e, long double i,
			long double RAAN, long double arg_per, long double nu);
    void addPerturbation(LDVector (*funcptr)(void *, const LDVector&),void* param);
    void propagate();
    void sv2oe(LDVector& init_sv);

};


#endif /* PROPAGATORS_H_ */
