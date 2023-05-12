/*
 * propagators.h
 *
 *  Created on: 28 mar. 2023
 *      Author: ceci
 */

#ifndef PROPAGATORS_H_
#define PROPAGATORS_H_

#include <vector>
#include"DVector.h"

struct cBody_param{
	DVector sv;
	DVector acc;
	double mu;
};

struct thrust_param {
	DVector sv;
	DVector acc_thr;
	double isp;
	double thrust;

};

DVector central_body(void* param);
DVector thrust(void *param);


class propagator
{
	private:
	DVector total_acc;
	vector<DVector (*)(void *, const DVector&)> vec_functions; // vector of function pointers
	vector<void *> vec_params; // vector of parameters pointers (Structures)
	int total_time;
	double step;
	DVector init_sv;
	DVector sv;
	DVector RK4();
	DVector derivatives(DVector& dv);

    public:

	propagator(DVector& init_sv, int total_time, double step);

    void addPerturbation(DVector (*funcptr)(void *, const DVector&),void* param);
    void propagate();

};


#endif /* PROPAGATORS_H_ */
