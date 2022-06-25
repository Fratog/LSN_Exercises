#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include "block_method.h"

using namespace std;

//Attention: I am aware this class doesn't make sense. In the sense that is just a collection of public functions, so it could have been just a collection of functions. Later on during the cours I fixed this, in the sense that I made a file statistic.cpp. When I made these class I had something in mind (which i can't remeber now) that in the end i didn't carry on. Though since the program works using this functions, i won't change it to avoid problems.

Block_method::Block_method(){}
//Block_method::Block_method(int i=1, int j =100) : count(i), Nblocks(j) {}

void Block_method::set_count(int i) {count=i;}

const int Block_method::get_count() {return count;}

double Block_method::Raverage_predicate(double x, double y)
{ 
	count++; //si asepetta che alla prima chiamata count sia uguale a 1 
	return((x*(count-1)+y))/count;
}


double Block_method::get_mean(vector <double> A) 
{
	return (accumulate(A.begin(), A.end(),0.0)/A.size()); 
}

double Block_method::get_error(vector <double> A) //A: vettore delle quantità stimate in N blocchi, che vanno gia divise per Nthrows per block!!!
{
int dim=A.size();
double error=0;
double AV_A2=0;
	vector <double> A2(dim);

	for (int i=0; i<dim; i++)
		A2[i]=pow(A[i],2);

	AV_A2=this->get_mean(A2);
	error=sqrt((AV_A2-pow(this->get_mean(A),2))/(dim-1)); 

	return error;
}


