/*
 * scenario1.cpp
 *
 *  Created on: 1 abr. 2024
 *      Author: ceci
 */

//double hp=480; // Perigee altitude
//double ha=800; // Apogee altitude
//double rp=Rt+hp; // Perigee radius
//double ra=Rt+ha; // Apogee radius
//long double semimajorAxis=(rp+ra)/2;
//double aux=get_max_absolute(init_sv); ????????????????????????????????
//cout << aux << endl;

void scenario01(){

	//-----------------------------------------
	// Orbit 1 - Initial - Central Body
	//-----------------------------------------
/*	std::cout <<std::endl;
	std::cout <<std::endl;
	std::cout <<"---------------------"<<std::endl;
	std::cout <<"- CENTRAL BODY "<<std::endl;
	std::cout <<"---------------------"<<std::endl;
	float step=1.0;
	double period0=period(semimajorAxis);
	int total_time=ceil(period0);
	string filename0="output_files/orbit1.csv";

	// First propagation instance
	propagator cbody_prop(init_sv,total_time,step);
	cBody_param cbody;
	cbody.mu=398600.448;
	LDVector cbody_last_sv;

	// Add perturbation
	cbody_prop.addPerturbation(&central_body, &cbody); // args: function and structure
	const clock_t begin_time = clock();
	cbody_prop.propagate(filename0);
	cbody_last_sv=cbody_prop.last_sv;

	std::cout <<cbody_prop.last_sv <<std::endl;
	cbody_prop.sv2oe(cbody_prop.last_sv);
	std::cout <<"Semimajor Axis : "<<cbody_prop.a <<std::endl;*/

}


