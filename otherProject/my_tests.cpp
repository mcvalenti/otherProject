/*
 * my_tests.cpp
 * =======================================================================	
 * This files will check some functions used as middle steps.
 * =======================================================================	
 *  Created on: 11 ago. 2023
 *      Author: ceci
 */

#include <iostream>
#include "LDVector.h"
#include "propagators.h"

void anomaly_propagation(){
	/* 
	 * Tests: "sv_from_true_anomaly" and including "lagrange_coeff_from_true_anomaly"
	 */
	long double sv_vector[6]={8182.4,-6865.9,0,0.47572,8.8116,0};
	long double sv_result[6]={1454.9864223053, 8251.4609450803,0.0000000000,-8.1323880267,5.6785406852,-0.0000000000};
	double tolerance=1e-10;
	LDVector sv_0(sv_vector,6);
	LDVector result(sv_result,6);
	LDVector sv_1(6);
	double delta_nu=120.0;
	sv_1=sv_from_true_anomaly(sv_0,delta_nu);
	if (sv_1.areEqual(result,tolerance)){
		cout << "TEST OK" << endl;
	}else{
		cout << "TEST NOK" << endl;
	}
	
}


