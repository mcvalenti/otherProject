#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "propagators.h"
#include "DVector.h"
using namespace std;

struct thrust_param {
	double isp;
	double thrust;
};

void sum(void *param){
    thrust_param* thr_p = (thrust_param*)param;
    cout<<"Sum of thrust parameters: "<<thr_p->isp+thr_p->thrust<<endl;
}

void add_perturbation(void (*funcptr)(void *),void* param){
	vector<void (*)(void *)> vec_functions; // vector of function pointers
	vector<void *> vec_params; // vector of parameters pointers (Structures)
	int i;
	vec_functions.push_back(funcptr);
	vec_params.push_back(param);

    i=0;
    for (auto f : vec_functions){
         f(vec_params[i]);
         i++;
    }
}

class propagator
{
	private:
	vector<void (*)(void *)> vec_functions; // vector of function pointers
	vector<void *> vec_params; // vector of parameters pointers (Structures)
	// DVector de acceleraciones.

    public:
    void addPerturbation(void (*funcptr)(void *),void* param);
    void propagate();
};
void propagator::addPerturbation(void (*funcptr)(void *),void* param)
{
        	vec_functions.push_back(funcptr);
			vec_params.push_back(param);
 }
void propagator::propagate()
{	// Debe devolver el DVector de acceleraciones acumulado
	// Se puede?
	int i;
    i=0;
    for  (auto f : this->vec_functions){
         f(this->vec_params[i]); // Cada f debe devolver un DVector de acceleraciones
         i++;
    }
}

int main() {
	//double mu=398600.448;
	unsigned size=6;
	double sv[size]={42164.0, 0.0, 0.0, 0.0, 3.0, 0.0};
	DVector dv(sv, size);
	DVector dv1;
	double h=0.5;
	int i;
	int time=86400;
	vector<DVector> sv_list;
	ofstream outputFile;
	
	// Propagation
	for (i=0; i<=time; i++){
		dv1=RK4(h,dv);
		sv_list.push_back(dv1);
		dv=dv1;
	}
	
	outputFile.open("output.csv");
	for (DVector sv : sv_list)
		outputFile << sv<<endl;
	outputFile.close();
	
	cout<<"Keplerian propagation DONE!"<<endl;
	
	// Vector of functions manipulation
	DVector acc_cb(size);
	cBody_param cbody;
	cbody.sv=dv;
	cbody.acc=acc_cb;
	cbody.mu=398600.448;
	thrust_param tparam;
	tparam.isp = 1200;
	tparam.thrust = 10000;

	// Add both functions to vec_functions
	add_perturbation(&central_body, &cbody);
	add_perturbation(&sum, &tparam);


	// Propagator as a Class
	cout<<"============================="<<endl;
	cout<<"PROPAGATOR CLASS"<<endl;
	cout<<"============================="<<endl;
	propagator myprop;

	myprop.addPerturbation(&central_body, &cbody);
	myprop.addPerturbation(&sum, &tparam);
	myprop.propagate();


	// End of Code
	cout<<"End of processing!";
			
	return 0;
}
