#include <iostream>
#include	<cmath>
//#include "random.h"
#include "cities.h"

using namespace std;

Cities::Cities(int N, int dim, int primes_row )
{
Ndim=dim;
Ncities=N;

	if(primes_row>31500)
	{
		cout<<"primes row too large"<<endl;
		rnd.starting(1);
	}		
	else	
	{	
		rnd.starting(primes_row);
	}

	for(int i=0; i<Ncities; i++)
	{
		vector<double> v;
		for(int j=0; j<Ndim; j++)
		{
			v.push_back(0);
		}
	coords.push_back(v);	
	}

//vector<double> colonna(dim);
//	for(int i=0; i<Ncities; i++)
//		coords.push_back(colonna);
		
}

void Cities::cities_on_hypercube(double L) 
{
	for(int i=0; i<Ncities; i++)
	{
		for(int j=0; j<Ndim; j++)
		{
			coords[i][j]=rnd.Rannyu()*L;
		}
	}
}

void Cities::cities_on_circle(double R)
{
double x=0;
	for(int i=0; i<Ncities; i++)
	{
			x=rnd.Rannyu();
			coords[i][0]=R + R*cos(x*2*M_PI);
			coords[i][1]=R + R*sin(x*2*M_PI);
	}
}

void Cities::print_coords() const
{
	for(int i=0; i<Ncities; i++)
	{
	cout<<endl;
		for(int j=0; j<Ndim; j++)
		{
			cout<<coords[i][j]<<" ";
		}
	}
	cout<<endl;
}

vector<vector<double>> Cities::get_coords() 
{
	return coords;
}

double Cities::get_coords(int i, int j) 
{
	return coords[i][j];
}
