#include <iostream>
#include <cmath>
#include <vector>
#include "statistics.h"
#include "vmc.h"

using namespace std;


VMC::VMC(double a, int dim, double (*fun)(std::vector<double>), void(RWalk::*type_of_rw)()) : Metropolis(a, dim , fun, type_of_rw)
{

}


VMC::VMC(double a, int dim, double (*fun)(vector <double>), void(RWalk::*type_of_rw)(), double (*cost_fun)(vector <double>)) : Metropolis(a, dim , fun, type_of_rw)
{
	cost=cost_fun;
}


VMC::~VMC()
{
	cost=nullptr;
}


double VMC::sample_cost() //puoi anche evitarlo forse,m nel senso usi sample point e quadno serve applichi cost senza avere un metodo
{
	this->sample_point(); //execute metropolis step
	return cost(this->get_coords()); //evaluate the cost function at the current RWalk configuration
}


vector<double> VMC::get_cost_vector(int Nblocks, int Nstep) 
{
	vector <double> A(Nblocks);
	
	for (int i=0; i<Nblocks; i++)
	{
		for(int j=0; j<Nstep; j++)
		{
			A[i]+=this->sample_cost();
		}
		A[i]/=Nstep;
	}
	
	return A;	
}


void VMC::add_block(vector<double> & vec, int Nstep) 
{
	double appo=0;
		
	for(int j=0; j<Nstep; j++)
		appo+=this->sample_cost();
	
	appo/=Nstep;
	vec.push_back(appo);	
}


vector<double> VMC::get_cost_value(int Nblocks, int Nstep) 
{
	vector<double> A(Nblocks);
	vector<double> result(2);

	A=this->get_cost_vector(Nblocks, Nstep);
	result[0]=mean(A);
	result[1]=sd_mean(A);

	return result;	
}


vector<double> VMC::get_cost_value(vector<double> A) 
{
	vector<double> result(2);

	result[0]=mean(A);
	result[1]=sd_mean(A);

	return result;	
}

