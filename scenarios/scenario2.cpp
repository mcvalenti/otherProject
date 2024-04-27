/*
 * scenario02.cpp
 *
 *  Created on: 4 abr. 2024
 *      Author: ceci
 */

	//-----------------------------------------
	// Orbit 2 - With Continuous thrust
	//-----------------------------------------
	/*
	std::cout <<"---------------------"<<std::endl;
	std::cout <<"- CONTINUOUS THRUST "<<std::endl;
	std::cout <<"---------------------"<<std::endl;
	float thrust_time=261.1127;
	float thrust_step=0.1;
	string filename1="output_files/orbit2.csv";
	cBody_param cbody;
	cbody.mu=398600.448;
	propagator thrust_prop(init_sv,thrust_time,thrust_step);
	thrust_param tparam;
	tparam.isp = 300.0;
	tparam.thrust = 10000.0;
	thrust_prop.addPerturbation(&central_body, &cbody);
	thrust_prop.addPerturbation(&thrust, &tparam);
	thrust_prop.propagate(filename1);

	cout<<" End-of-burn state vector: "<<thrust_prop.last_sv<<endl;
	thrust_prop.sv2oe(thrust_prop.last_sv);
	std::cout <<"True anomaly nu: "<<thrust_prop.nu <<std::endl;
	double delta_nu=(M_PI-thrust_prop.nu)*180.0/M_PI; // computes delta true anomaly to apogee
	std::cout <<"Delta anomaly nu: "<<delta_nu <<std::endl;
	LDVector final = sv_from_true_anomaly(thrust_prop.last_sv,delta_nu); //to compute radius at final point
	cout << "At apogee: " << final << endl;*/



	/*




	//-----------------------------------------
	// Orbit 2 - With Continuous thrust
	//-----------------------------------------
	std::cout <<"---------------------"<<std::endl;
	std::cout <<"- CONTINUOUS THRUST "<<std::endl;
	std::cout <<"---------------------"<<std::endl;
	float thrust_time=261.1127;
	float thrust_step=0.1;
	string filename1="output_files/orbit2.csv";
	propagator thrust_prop(init_sv,thrust_time,thrust_step);
	thrust_param tparam;
	tparam.isp = 300.0;
	tparam.thrust = 10000.0;
	thrust_prop.addPerturbation(&central_body, &cbody);
	thrust_prop.addPerturbation(&thrust, &tparam);
	thrust_prop.propagate(filename1);


	cout<<" End-of-burn state vector: "<<thrust_prop.last_sv<<endl;
	thrust_prop.sv2oe(thrust_prop.last_sv);
	std::cout <<"True anomaly nu: "<<thrust_prop.nu <<std::endl;
	double delta_nu=(M_PI-thrust_prop.nu)*180.0/M_PI; // computes delta true anomaly to apogee
	std::cout <<"Delta anomaly nu: "<<delta_nu <<std::endl;
	LDVector final = sv_from_true_anomaly(thrust_prop.last_sv,delta_nu); //to compute radius at final point
	cout << "At apogee: " << final << endl;

	// if not ra=22378 --> increase thrust_time

	//Estimation
	long double semimajor_axis_transfer=14618.0;
	double period2=period(semimajor_axis_transfer);
	// transfer Free propagation
	propagator free_prop(thrust_prop.last_sv,period2-thrust_time,step);
	cBody_param free_body;
	free_body.mu=398600.448;
	LDVector free_body_last_sv;

	// Add perturbation
	free_prop.addPerturbation(&central_body, &free_body); // args: function and structure
	string filename3="output_files/orbit3.csv";
	free_prop.propagate(filename3);
	free_prop.sv2oe(free_prop.last_sv);
	std::cout <<"Semimajor Axis : "<<free_prop.a <<std::endl;

	// End of Propagation
	std::cout << "Propagation time: "<<float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;





	// End of Code
	std::cout << "Propagation time: "<<float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
	cout << endl;
	cout<<"End of processing!";
	*/



