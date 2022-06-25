#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "random.h"
#include "random_walk.h"

using namespace std;

RWalk::RWalk(double a, int d)
{
	dim=d;
	step.push_back(-a); //postresti generalizzare a un vettore di step piu grande ma per ora fa niente;
	step.push_back(a);
	coords.resize(d);
	last_coords.resize(d);
}

vector <double> RWalk::get_coords()
{
	return coords;
}

void RWalk::move_to_origin()
{
	fill(coords.begin(),coords.end(),0);
}

vector <double> RWalk::get_step()
{
	return step;
}

void RWalk::RW_move_discrete()
{
	int dir=int(rnd.Rannyu(0,dim));
	int move=step[(int(rnd.Rannyu(0, step.size())))];
	
	last_coords=coords;
	coords[dir]+=move;
}


void RWalk::RW_move_discrete(int nstep)
{
	for (int i=0; i<nstep; i++)
		this->RW_move_discrete();
}

void RWalk::RW_move_back()
{
	coords=last_coords;
}


void RWalk::RW_move_continuos()///potresti usare coordinate sferiche in n diensioni e generalizzare?? non so se è cos' scontato, bisogna capire come campionare le variabili
//credo si possa fare girare il programma con dim>3, la parte del programma che lavora sul continuo non da problemi, semplicemente usa le prime tre celle e rida quindi i risultati di un caso dim=3. ma non dovreebbe dare errori dato che il resto del vector rimane sempre zero
{
	double theta=0;
	double phi=0;
	double a=abs(step[0]);
	
	last_coords=coords;
	//rnd genera teta e phi
	theta=this->rnd.Theta();
	phi=this->rnd.Rannyu(0,2*M_PI);

	if( theta < 0 or theta > M_PI) 
		cout<<"error, theta :"<<theta<<endl;

	if( phi < 0 or  phi > 2*M_PI) 
		cout<<"error, phi :"<<phi<<endl;

	coords[0]+=a*sin(theta)*cos(phi);
	coords[1]+=a*sin(theta)*sin(phi);
	coords[2]+=a*cos(theta);
	
}


void RWalk::RW_move_continuos(int nstep)
{
	for (int i=0; i<nstep; i++)
		this->RW_move_continuos();
}


double RWalk::get_norm_squared()
{
	double norm=0;
	
	for (auto i : coords)
		norm+=i*i;

	return norm;
//return sqrt(inner_product(coords.begin(), coords.end(), coords.begin(), 0.0L)); 
}


