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

void central_body(void* param);
DVector _deriv(DVector& dv);
DVector keplerian(DVector& dv, double mu);
DVector thrust(DVector& dv, double T, double isp);
vector<DVector> collects_perturbations(DVector& dv, vector<string>& perturbations);
DVector RK4(double h, DVector& dv);



#endif /* PROPAGATORS_H_ */
