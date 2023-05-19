#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "propagators.h"
#include "LDVector.h"
using namespace std;


int main() {

	unsigned size=7;
	long double sv_m[size]={6858,0,0,0,7.7102,0,2000};
	LDVector init_sv(sv_m, 7);


	// Propagator as a Class
	cout<<"============================="<<endl;
	cout<<"PROPAGATOR CLASS"<<endl;
	cout<<"============================="<<endl;
	long double step=1.0;
	int total_time=5851.01245965262;
	
	propagator myprop(init_sv,total_time,step);

	cBody_param cbody;
	thrust_param tparam;
	cbody.mu=398600.448;
	tparam.isp = 300;
	tparam.thrust = 10000;


	// Agrego las perturbaciones que van a ser consideradas
	myprop.addPerturbation(&central_body, &cbody);
	//myprop.addPerturbation(&thrust, &tparam);
	const clock_t begin_time = clock();
	myprop.propagate();
	std::cout << "Excecution time!"<<float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
;
	// End of Code
	cout<<"End of processing!";
			
	return 0;
}
