#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include "block_method.h"

using namespace std;

Block_method::Block_method() 
{
	count=0;
	Nblocks=0;
}

Block_method::Block_method(int i, int j)
{
	count=i;
	Nblocks=j;
}

void Block_method::set_count(int i) {count=i;} 
void Block_method::set_Nblocks(int i) {Nblocks=i;} 

const int Block_method::get_count() {return count;}
const int Block_method::get_Nblocks() {return Nblocks;}

vector<double> Block_method::get_running_av(vector <double> A) 
{
	vector<double> running_av(A.size());

	partial_sum(A.begin(), A.end(), running_av.begin(), [appo = 1] (double x, double y) mutable {
																								appo+=1; 
																								return((x*(appo-1)+y))/appo;
																								});
	return running_av; 
}

vector<double> Block_method::get_running_err(vector <double> A) 
{
	vector<double> A2(A.size());
	vector<double> running_av_A(A.size());
	vector<double> running_err(A.size());

	running_av_A=get_running_av(A);

	for(unsigned int i=0; i<A.size(); i++)
		A2[i]=pow(A[i],2);

	partial_sum(A2.begin(), A2.end(), running_err.begin(), [appo = 1] (double x, double y) mutable {
																									appo+=1; 
																									return((x*(appo-1)+y))/appo;
																									});											

	running_err[0]=0;
	for(unsigned int i=1; i<A.size(); i++)
		running_err[i]=sqrt((running_err[i]-pow(running_av_A[i],2))/i);											

	return running_err; 
}
											 

double Block_method::get_mean(vector<double> A) 
{
	return (accumulate(A.begin(), A.end(),0.0)/A.size()); 
}

double Block_method::get_error(vector<double> A) //A: vettore delle quantità stimate in N blocchi, che vanno gia divise per Nthrows per block!!!
{
	int dim=A.size();
	double error=0;
	double AV_A2=0;
	vector<double> A2(dim);

	for (int i=0; i<dim; i++)
		A2[i]=pow(A[i],2);

	AV_A2=this->get_mean(A2);
	error=sqrt((AV_A2-pow(this->get_mean(A),2))/(dim-1)); 

	return error;
}


