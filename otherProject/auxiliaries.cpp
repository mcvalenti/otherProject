/*
 * auxiliaries.cpp
 *
 *  Created on: 19 may. 2023
 *      Author: ceci
 */

#ifndef AUXILIARIES_CPP_
#define AUXILIARIES_CPP_

#include "LDVector.h"

LDVector cross(LDVector& a, LDVector& b){
	LDVector cross_prod(3);
	cross_prod[0]=a[1]*b[2]-b[1]*a[2];
	cross_prod[1]=-(a[0]*b[2]-b[0]*a[2]);
	cross_prod[2]=a[0]*b[1]-b[0]*a[1];
	return cross_prod;	
}




#endif /* AUXILIARIES_CPP_ */
