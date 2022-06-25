#include <iostream>
#include <cmath>
#include <vector>
#include "metropolis.h"
#include "random.h"
#include "random_walk.h"

using namespace std;

Metropolis::Metropolis(double a, int dim, void(RWalk::*type_of_rw)()) : RWalk(a, dim) 
{	
	Tyx=type_of_rw; 
	mypdf=nullptr; 
	accepted=0;
	rejected=0;
}

Metropolis::Metropolis(double a, int dim, double (*fun)(vector <double>), void(RWalk::*type_of_rw)() ) : RWalk(a, dim) 
{	
	Tyx=type_of_rw; 
	mypdf=fun; 
	accepted=0;
	rejected=0;
}

Metropolis::~Metropolis()
{
Tyx=nullptr;
mypdf=nullptr;
}


void Metropolis::set_pdf(double (*fun)(vector <double>))
{
	mypdf=fun;
}


void Metropolis::set_Tyx(void(RWalk::*type_of_rw)()) //return_type(Class::*name)(type_arg1, ..., type_argn)
{
	Tyx=type_of_rw; 
}


void Metropolis::set_initial_point(int i, double r) 
{
	if(i > this->get_dim()) {cout<<"Error"<<endl;}
	this->set_i_coord(i, r); 
}


void Metropolis::equilibrate(int Nequil) //sarebbe bello farlo in modo più intelligente, 
{
	for(int i=0; i<Nequil; i++)	
		this->sample_point(); 
}


void Metropolis::find_step(double step_in, int Nequil) //choose wisely the intial step anyway, to avoid being inefficent
{
	double step=step_in;
 	double rate=0;
	double rate_min=0.47;
	double rate_max=0.52;
	int counter=0;
	
	this->set_step_length(step);
	do
	{
			accepted=0;
			rejected=0;
			counter++;
			this->equilibrate(Nequil); 
			rate=(double)accepted/(accepted+rejected);
			//cout<<"trial rate	"<<rate<<endl;	
			if(this->get_step_length()<0.00001) //it seems it solved the problem, but it is obvious it is a ad hoc situation
			{
				this->set_step_length(1);
				this->equilibrate(10);
			}
			
			if(rate<rate_min)
			{
				this->set_step_length(step-step/10.); //subtract ten percent
				step=step-step/10.;
			}
			else if(rate>rate_max)
			{
				this->set_step_length(step+step/10.);
				step=step+step/10.;
			}
			else
			{
				this->set_step_length(step);
			}
			if(counter>1000) {break;}
	}
	while(!(rate_min<rate and rate<rate_max));
	//cout<<"acceptance rate	from find_step() call"<<(double)accepted/(accepted+rejected)<<endl;
	accepted=0;
	rejected=0;
}


void Metropolis::sample_point() 
{
	double x;
	double alpha=0;
	
	x=mypdf(this->get_coords());//save pdf value at previous point
	//this->*Tyx --> this->data_member and then we add (). it is a bit confusing is like the Tyx variable kepsp the function
	(this->*Tyx)(); //make the random walk move
	

	alpha=mypdf(this->get_coords())/x;
	
	if(alpha>1)
	{
		accepted++;
	}
	else if((this->get_gen()).Rannyu()<=alpha)  //get_gen return the generator by reference
	{
		accepted++;
	}
	else
	{
		this->move_back(); 
		rejected++;
	}	
}


bool Metropolis::sample_point(double xnew, double xold) 
{
	double alpha=0;
	
	alpha=xnew/xold;
	
	if (alpha>1)
	{
		accepted++;
		return true;
	}
	else if((this->get_gen()).Rannyu()<=alpha)  //get_gen return the random generator by reference
	{
		accepted++;
		return true;
	}
	else
	{
		this->move_back(); 
		rejected++;
		return false;
	}	
}

