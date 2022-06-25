#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm> //for_each
#include <numeric> //accumulate, partial_sum, reduce
#include <vector>
#include "statistics.h"

using namespace std;


double mean(vector<double> vec)
{
	return accumulate(vec.begin(), vec.end(), 0.0)/vec.size(); 
}


//https://docs.microsoft.com/it-it/cpp/cpp/lambda-expressions-in-cpp?view=msvc-170
double variance(vector<double> vec) //vector variance
{
	double m=0;
	double sum=0;
	m=mean(vec);

	for_each(vec.begin(), vec.end(), [&] (const double d) {   //& lmambd,body gets acces to values out of his scope, by reference
		  sum += (d - m) * (d - m);
			});

	return sum/vec.size(); 
}

double sd(vector<double> vec)
{
	return sqrt(variance(vec));
}

double sd_mean(vector<double> vec)
{
	return sqrt(variance(vec)/vec.size());
}

void autocorrelation(vector<double> A, int Nstep, string filename)
{
int dim=A.size();
double aa, a1, a2, var;

ofstream out;
	out.open(filename);
	
	var=variance(A);
	for (int i=0; i<Nstep; i++)
	{
		aa=0;
		a1=0;
		a2=0;
		for (int j=0; j<(dim-i); j++)
		{
			aa+=A[j]*A[j+i];
			a1+=A[j];
			a2+=A[j+i];
		}
	  aa/=(dim-i);
    a1/=(dim-i);
    a2/=(dim-i);
		out<<(aa-a1*a2)/(var)<<endl;
	}
	
	out.close();
}

vector<double> blocking(vector <double> ist, int Nblocks, int Nstep_x_block) //it expects the vector of instant values
{
	vector<double> blocks(Nblocks);
	
	for(int i=0; i<Nblocks; i++)
	{
		auto it1 = ist.begin()+i*Nstep_x_block;
			//auto it2 = it1 + Nsteps;
		blocks[i]=accumulate(it1, it1+Nstep_x_block, 0.0);
		blocks[i]/=Nstep_x_block;
	}	
return blocks;
}

vector<double> blocking_lento(vector <double> ist, int Nblocks, int Nstep_x_block) //it expexts the vector of instant values
{
	vector<double> blocks(Nblocks);
	
	for(int i=0; i<Nblocks; i++)
	{
		for(int j=0; j<Nstep_x_block; j++)
		{
			blocks[i]+=ist[j+i*Nblocks];
		}
		blocks[i]/=Nstep_x_block;
	}	
return blocks;
}

vector<double> running_av(vector <double> A) //it expexts a vector of Nblock with entries already divided by Nstep x block
{
	vector<double> running_av(A.size());
	partial_sum(A.begin(), A.end(), running_av.begin(), [appo = 1] (double x, double y) mutable {
																								appo+=1; 
																								return((x*(appo-1)+y))/appo;
																								});
	return running_av; 
} 


vector<double> running_err(vector <double> A) 
{
	vector<double> running_av_A(A.size());
	vector<double> running_err(A.size());

	running_av_A=running_av(A);
	//here i get the partial sum vector of the squared elemnts of A, the first elemnt of this vector won't be actually x^2, but doesn't mind since the first element of running_err will be 0
	
	//running_err[0]=pow(A[0],2);
	for (auto & i : A) { i=i*i; }//i am channging A, thought it is passed by value to the function, so it is not a problem
	partial_sum(A.begin(), A.end(), running_err.begin(), [appo = 1.] (double x, double y) mutable {
																									appo+=1; 
																									return((x*(appo-1)+y))/appo;
																									});											

	running_err[0]=0;
	for (unsigned int i=1; i<A.size(); i++)
		running_err[i]=sqrt((running_err[i]-pow(running_av_A[i],2))/i);											

	return running_err; 
}



vector<double> bm_running_err(vector<double> A) 
{
	vector <double> A2(A.size());
	vector <double> running_av_A(A.size());
	vector <double> running_err(A.size());

	running_av_A=running_av(A);

	for(unsigned int i=0; i<A.size(); i++)
		A2[i]=pow(A[i],2);

	partial_sum(A2.begin(), A2.end(), running_err.begin(), [appo = 1] (double x, double y) mutable {
																									appo+=1; 
																									return((x*(appo-1)+y))/appo;
																									});											

	running_err[0]=0;
	for (unsigned int i=1; i<A.size(); i++)
		running_err[i]=sqrt((running_err[i]-pow(running_av_A[i],2))/i);											

	return running_err; 
}

bool compatible(vector<double> x2, vector<double> x1) //should return true if compatible
{
	return (abs(x2[0]-x1[0])<(x2[1]+x1[1]));
}
