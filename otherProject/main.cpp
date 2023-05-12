#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "propagators.h"
#include "DVector.h"
using namespace std;


int main() {
	//double mu=398600.448;
	unsigned size=7;
	double sv_m[size]={42164.0, 0.0, 0.0, 0.0, 3.0, 0.0,1200.0};
	DVector init_sv(sv_m, 7);


	// Propagator as a Class
	cout<<"============================="<<endl;
	cout<<"PROPAGATOR CLASS"<<endl;
	cout<<"============================="<<endl;
	double step=0.5;
	int total_time=100;
	
	propagator myprop(init_sv,total_time,step);
	myprop.propagate();


	// End of Code
	cout<<"End of processing!";
			
	return 0;
}
