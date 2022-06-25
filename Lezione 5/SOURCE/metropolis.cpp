#include <iostream>
#include <cmath>
#include <vector>
#include "metropolis.h"
#include "random.h"
#include "random_walk.h"

using namespace std;

Metropolis::Metropolis(double a, int dim, double (*fun)(vector <double>), void(RWalk::*type_of_rw)() ) : RWalk(a, dim) 
{	
	Tyx=type_of_rw; //set the tyoe of random walk used for the transition probability
	mypdf=fun; //set the pdf we want to sample
	accepted=0;
	rejected=0;
}

//aggiungi un distruttore che mette a nullptr i puntatori;

void Metropolis::set_pdf(double (*fun)(vector<double>))
{
	mypdf=fun;
}

void Metropolis::set_Tyx( void(RWalk::*type_of_rw)() ) //return_type(Class::*name)(type_arg1, ..., type_argn)
{
	Tyx=type_of_rw; 
}

void Metropolis::set_initial_point(int i, double r) //it wasn't really necessary to implement this function
{
	if(i > this->get_dim()) {cout<<"Error"<<endl;}
		this->set_i_coord(i, r); 
}

vector<double> Metropolis::equilibrate(int Nequil) //it isn't really general
{
	vector<double> Obs; //Obs observable choosen to check equilibration
	for(int i=0; i<Nequil; i++)	
	{
		this->sample_point(); 
		Obs.push_back(sqrt(this->get_norm_squared()));
	}	
	return Obs;
}

void Metropolis::tuning() //I wanted to do an lagorithm to select the metropolis step to have 50% acceptance rate, but here I didn't implement it
{

	cout<<(double)accepted/(accepted+rejected)<<endl;
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
	
	if (alpha>1)
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

