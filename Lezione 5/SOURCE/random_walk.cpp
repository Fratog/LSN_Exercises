#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "random.h"
#include "random_walk.h"

using namespace std;

//construnctor
RWalk::RWalk(double a, int d)
{
	dim=d;
	step.push_back(-a); 
	step.push_back(a);
	coords.resize(d);
	last_coords=coords;
}

//set get
void RWalk::set_dim(int d) {dim=d;}

void RWalk::set_step_length (double a) 
{
	step[0]=-a;
	step[1]=a;
}

int RWalk::get_dim() 
{
	return dim;
}

 Random & RWalk::get_gen() 
{
	return rnd;
}

vector<double> RWalk::get_coords()
{
	return coords;
}

void RWalk::set_origin(vector<double> origin)
{
	coords=origin;
}

void RWalk::set_i_coord(int i, double x) 
{
	coords[i]=x;
}

void RWalk::print_coords()
{
	for(auto i : this->get_coords()) {cout<<i<<"\t";}
}

double RWalk::get_step_length()
{
	return abs(step[0]);
}

//move functions

void RWalk::move_to_origin()
{
	fill(coords.begin(),coords.end(),0);
}

void RWalk::move_discrete()
{
	int dir=int(rnd.Rannyu(0,dim)); 
	int move=step[(int(rnd.Rannyu(0, step.size())))];

	last_coords=coords;
	coords[dir]+=move;

}


void RWalk::move_discrete(int nstep)
{
	for(int i=0; i<nstep; i++)
		this->move_discrete();
}

void RWalk::move_continuos() //works only in 3dim obviously
{
	double theta=0;
	double phi=0;
	double a=abs(step[0]);

	theta=this->rnd.Theta();
	phi=this->rnd.Rannyu(0,2*M_PI);

	if( theta < 0 or theta > M_PI) 
		cout<<"error, theta :"<<theta<<endl;

	if( phi < 0 or  phi > 2*M_PI) 
		cout<<"error, phi :"<<phi<<endl;

	last_coords=coords;

	coords[0]+=a*sin(theta)*cos(phi);
	coords[1]+=a*sin(theta)*sin(phi);
	coords[2]+=a*cos(theta);

}


void RWalk::move_continuos(int nstep)
{
	for(int i=0; i<nstep; i++)
		this->move_continuos();
}


void RWalk::move_back()
{
	coords=last_coords;
}

void RWalk::move_uniform() 
{
	double a=abs(step[0]);

	last_coords=coords;
	for (int i=0; i<dim; i++)
		coords[i]+=rnd.Rannyu(-a,a);
}

void RWalk::move_uniform2() 
{
	double a=abs(step[0]);
	int dir=int(rnd.Rannyu(0,dim));
	
	last_coords=coords;
	coords[dir]+=rnd.Rannyu(0,a);
}


void RWalk::move_gauss() 
{
	double a=abs(step[0]);
	
	last_coords=coords;
	for (int i=0; i<dim; i++)
		coords[i]+=rnd.Gauss(0,a);	
}

double RWalk::get_norm_squared()
{
	return inner_product(coords.begin(), coords.end(), coords.begin(), 0.0L); 
}


