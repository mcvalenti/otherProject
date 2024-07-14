/*
 * interplanetary.h
 *
 *  Created on: Jul 13, 2024
 *      Author: macec
 */

#ifndef INTERPLANETARY_H_
#define INTERPLANETARY_H_

void run_interplanetary();
double delta_V_homann(double mu_center, double R1, double R2);
double departure_v(double v_escape, double mu_center, double r_departure);
double parking_v(double mu_center, double r_departure);




#endif /* INTERPLANETARY_H_ */
