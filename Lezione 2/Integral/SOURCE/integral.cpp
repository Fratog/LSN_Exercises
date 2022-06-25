#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "random.h"

using namespace std;
 
//functions and global variables
const double pi = 2*acos(0.0); 
double integrand(double x) {return (pi/2)*cos((pi*x)/2);} 

double myfun(double x, float n) {return ((n+1)/n)*(1-pow(x,n));}

int k=1; 
double Raverage_predicate(double x, double y) //this is a function we use to get the running averages, we use it as a predicate of partial_sum
{
	k++;
	//cout<<"k "<<k<<endl;
	//cout<<"x y "<<x<<"	"<<y<<"  e poi  ((x+y))/k ="<<((x+y))/k<<endl; //mi serviva per assicurarmi partial_sum facesse quello che desidero
	return((x*(k-1)+y))/k;
}
 

int main (int argc, char *argv[])
{

// setting  up the random generator
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
   	for(int i=0; i<34; i++)
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


///////////////////////////////Evaluation of the integral from uniform distribution
int Nblocks=100;
int Nthrows=1000000;
int L=Nthrows/Nblocks;

vector<double> I(Nblocks);
vector<double> sigmaI(Nblocks);
vector<double> runningAv(Nblocks);
vector<double> runningAv2(Nblocks);
vector<double> runningErr(Nblocks);

	for (int i=0; i<Nblocks; i++)
	{
		for (int j=0; j<L; j++)
			I[i]+=integrand(rnd.Rannyu());

		I[i]/=L;	
		sigmaI[i]=pow(I[i],2);
	}

	partial_sum(I.begin(), I.end(), runningAv.begin(), Raverage_predicate); //https://en.cppreference.com/w/cpp/algorithm/partial_sum
	k=1;
	partial_sum(sigmaI.begin(), sigmaI.end(), runningAv2.begin(), Raverage_predicate);
	k=1;

ofstream fout;
fout.open("I_uniform.txt");
 
	fout<<"here there are the progressive averages and errors for the MC calculation of the integral sampling a uniform distribution"<<endl;

	runningErr[0]=0;
	fout<<runningAv[0]<<"	"<<runningErr[0]<<endl;
	for (int i=1; i<Nblocks; i++)
	{
		runningErr[i]=sqrt((runningAv2[i]-pow(runningAv[i],2))/i);
		fout<<runningAv[i]<<"	"<<runningErr[i]<<endl;
	}
fout.close();

///////////////////////////////////////Evaluation of the integral with importance sampling
//integral you want
//find a suitbale pdf for importance sampling
//implement accept reject since you d'ont know the cumulative of the new pdf
//use imp samp for calculating the integral



double xi=0;//unifor random number
int count=0;//a counter
float n=2;//power index for myfun function
	
	fill(I.begin(), I.end(), 0);
	for (int i=0; i<Nblocks; i++)
	{
		count=0;
		while(count!=L)//L is the nuber of throws
		{
			xi=rnd.Rannyu(); 
			if(rnd.Rannyu()<myfun(xi,n)/myfun(0,n)) //accept reject method to sample myfun , mypdf(0) is the max value of my pdf since is a 1-x^2 functon
			{
				I[i]+=integrand(xi)/myfun(xi,n);//importance sampling method
				count++;
			}
		}
		I[i]/=L;	
		sigmaI[i]=pow(I[i],2);
	}

	//for (const auto & i : I)  {cout<<i<<endl;} 

	fill(runningAv.begin(), runningAv.end(), 0);
	fill(runningAv2.begin(), runningAv2.end(), 0);
	partial_sum(I.begin(), I.end(), runningAv.begin(), Raverage_predicate); //https://en.cppreference.com/w/cpp/algorithm/partial_sum
	k=1;
	partial_sum(sigmaI.begin(), sigmaI.end(), runningAv2.begin(), Raverage_predicate);
	k=1;

//for (const auto & i: runningAv )  {cout<<i<<endl;}

fout.open ("I_imp_samp.txt");
	fout<<"here there are the progressive averages and errors for the MC calculation of the integral using importance sampling"<<endl;

	runningErr[0]=0;
	fout<<runningAv[0]<<"	"<<runningErr[0]<<endl;
	for (int i=1; i<Nblocks; i++)
	{
		runningErr[i]=sqrt((runningAv2[i]-pow(runningAv[i],2))/i);
		fout<<runningAv[i]<<"	"<<runningErr[i]<<endl;
	}

fout.close();

rnd.SaveSeed();
return 0;

}
