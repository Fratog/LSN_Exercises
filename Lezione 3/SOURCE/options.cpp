#include <iostream>
#include <fstream>
//#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "random.h"
#include "block_method.h"

using namespace std;

int main(int argcv, char *argv[])
{
//parameters
double t=0;
double T=1;
double K=100;
double r=0.1;//fa il ruolo di r, è il parametro che caratterizza il drift, la crescita determinisitica costante
double sigma=0.25;
double S0=100;


//block method parameters
int Nthrows=10000;
int Nblocks=100;
int N=Nthrows/Nblocks;

//objects
Random rnd;
Block_method bm;


vector <double> S_block(Nblocks);//asset price per block //in realtà questo non ti serve , non è necessario salvare i valori di S

vector <double> C_block(Nblocks);//call price per block
vector <double> Call(Nblocks); //running averages for call prices
vector <double> Call_error(Nblocks);//running errors for call prices

vector <double> P_block(Nblocks); //analogous for put options
vector <double> Put(Nblocks);
vector <double> Put_error(Nblocks);


/////////////////////////full time evolution i.e. we sample directly the price of the asset a the final T 

	for (int i=0; i<Nblocks; i++)
	{
		for (int j=0; j<N; j++)
		{
			S_block[i]=S0*exp((r-0.5*sigma*sigma)*T+sigma*rnd.Gauss(0,T)); //sample asset price at final T //change S from a vector to a double
			if(S_block[i]-K>0) //check wheter the holder uses the option or not, both in call and put case
				C_block[i]+=exp(-r*T)*(S_block[i]-K);
			else
				P_block[i]+=-exp(-r*T)*(S_block[i]-K);
		}
		C_block[i]/=N;
		P_block[i]/=N;
	}

	Call=bm.get_running_av(C_block); //use block_method class to get the running averages and errors
	Call_error=bm.get_running_err(C_block);
//for (auto i : Call_error) {cout<<i<<endl;}

ofstream fout;
fout.open("call_fulltime.txt");

	for(int i=0; i<Nblocks; i++)
		fout<<Call[i]<<"	"<<Call_error[i]<<endl;

fout.close();

	Put=bm.get_running_av(P_block);
	Put_error=bm.get_running_err(P_block);

fout.open("put_fulltime.txt");

	for(int i=0; i<Nblocks; i++)
		fout<<Put[i]<<"	"<<Put_error[i]<<endl;

fout.close();





/////////////////discrete time evolution  
double passo=(T-t)/100.;  //define the length of time step
double S=0; //variable used to keep the value of the asset price at every instant of the time series

	fill(S_block.begin(), S_block.end(), 0.0);
	fill(C_block.begin(), C_block.end(), 0.0);
	fill(P_block.begin(), P_block.end(), 0.0);

	T=1;
	//cout<<passo<<endl;

	for (int i=0; i<Nblocks; i++)
	{

		for (int j=0; j<N; j++)
		{
			t=0;
			S=S0; //set the asset price to the intial value
			while (t<T) //evolve the asset price step by step
			{
				t+=passo;
				S=S*exp((r-0.5*sigma*sigma)*passo+sigma*rnd.Gauss(0,1)*sqrt(passo));
			}
			if(S-K>0) //check wheter the holder uses the option or not, both in call and put case
				C_block[i]+=exp(-r*T)*(S-K);
			else
				P_block[i]+=-exp(-r*T)*(S-K);
		}
		C_block[i]/=N;
		P_block[i]/=N;
	}

	Call=bm.get_running_av(C_block);
	Call_error=bm.get_running_err(C_block);


fout.open("call_discrete.txt");

	for (int i=0; i<Nblocks; i++)
		fout<<Call[i]<<"	"<<Call_error[i]<<endl;

fout.close();


	Put=bm.get_running_av(P_block);
	Put_error=bm.get_running_err(P_block);

fout.open("put_discrete.txt");

	for (int i=0; i<Nblocks; i++)
		fout<<Put[i]<<"	"<<Put_error[i]<<endl;

fout.close();

return 0;
}


