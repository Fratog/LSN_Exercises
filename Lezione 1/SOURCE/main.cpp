/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "random.h"

using namespace std;
 
const double pi = 2*acos(0.0);// use M_PI instead
double chii(int Oi){return (pow(Oi-100,2))/100;}  //not really elegant


int main (int argc, char *argv[]){

Random rnd;
   
///////////////////////////////////////////////ESERCIZIO 1.1

int M=1000000; //Number of random numbers to sample 
int N=100; //NBlocks
int L=int(M/N); //Nthrows per block
double sum=0; 

vector<double> ave(N);
//cout<<ave.size()<<endl;
//for (auto i : ave) {cout<<i<<endl;}
vector<double> ave2(N);
vector<double> sumprog(N);
vector<double> sum2prog(N);
vector<double> err_prog(N);

	for(int i=0; i<N; i++)
	{
		sum=0;
		for(int j=0; j<L; j++)	
		{
			sum+=rnd.Rannyu();
		}
		ave[i]=sum/L;
		ave2[i]=pow(ave[i],2);		
	}

	partial_sum(ave.begin(), ave.end(), sumprog.begin()); //fill sumprog with the partial sums terms of the series of elements in ave
	partial_sum(ave2.begin(), ave2.end(), sum2prog.begin()); 

	err_prog[0]=0;
	for(int i=1; i<N; i++)
	{
		sumprog[i]/=(i+1);
		sum2prog[i]/=(i+1);
		err_prog[i]=sqrt((sum2prog[i]-pow(sumprog[i],2))/i);
	}

ofstream fout;
fout.open ("averages.txt");
  
	fout<<"here there are the progressive averages and errors"<<endl;
	for (int i=0; i<N; i++ ) 
		fout<<sumprog[i]<<"	"<<err_prog[i]<<endl;
  
fout.close();

//now we do the same for the variance of the uniform distribution, I could have done both mean and variance in a single for cycle to reduce the execution time, though it was easier to just paste the previous code here and avoid creating other vectors to do both simultaneously.
	for(int i=0; i<N; i++)
	{
		sum=0;
		for(int j=0; j<L; j++)	
		{
			sum+=pow(rnd.Rannyu()-0.5,2); 
		}
		
		ave[i]=sum/L;
		ave2[i]=pow(sum/L,2);		
	}

	partial_sum(ave.begin(), ave.end(), sumprog.begin());
	partial_sum(ave2.begin(), ave2.end(), sum2prog.begin()); 

	//for (auto i : sumprog) cout<<i<<endl;
	err_prog[0]=0;
	for(int i=1; i<N; i++)
	{
		sumprog[i]/=(i+1);
		sum2prog[i]/=(i+1);
		err_prog[i]=sqrt((sum2prog[i]-pow(sumprog[i],2))/i);
		//cout<<err_prog[i]<<endl;	
	}

fout.open ("averages2.txt");
  
	fout<<"here there are the progressive averages and errors"<<endl;
	for (int i=0; i<N; i++ ) 
		fout<<sumprog[i]<<"	"<<err_prog[i]<<endl;
  
fout.close();
/////////////////////////////////////////////chi quadro test
  

//cout<<"CHI QUAD TEST"<<endl;
int Nchi=1000;// number of Chi to sample
int Nint=100;  
int Nthrows=10000;
int E=Nthrows/Nint;
vector<double> Chi(Nchi);
vector<int> ni(Nint);
vector<double> appo(Nint);

fout.open ("chiquad.txt");

	for(int i=0; i<Nchi; i++)
	{
		fill(ni.begin(), ni.end(), 0);

			for(int j=0; j<Nthrows; j++)
				ni[int(rnd.Rannyu()*Nint)]+=1; //build up the histogram

		transform(ni.begin(), ni.end(), appo.begin() , chii); //transform the histogram values to get the summands of the chi_qud value
		Chi[i]=accumulate(appo.begin(), appo.end(), 0.0); //get the sum
		fout<<Chi[i]<<endl;//print to file
	}

fout.close();

////////////////////////////////////////////////////ESERCIZIO 01.2

ofstream fout1;
ofstream fout2;
ofstream fout3;

fout1.open ("randu.txt");
	fout1<<"SN for randu with N=1,2,10,100"<<endl;
	
fout2.open ("exponential.txt");
	fout2<<"SN for exponential with N=1,2,10,100"<<endl;
	
fout3.open ("cauchy.txt");
	fout3<<"SN for randu with N=1,2,10,100"<<endl;

Nthrows=10000;
vector <int> Nind={1,2,10,100};

double SNu=0;
double SNexp=0;
double SNcauchy=0;
//for(auto n : Nind) {cout<<n<<endl;}

	for(auto n : Nind) 
	{
	//fout1<<Nthrows<<"valori di SN con N = "<<n<<endl<<endl;
	//fout2<<Nthrows<<"valori di SN con N = "<<n<<endl<<endl;
	//fout3<<Nthrows<<"valori di SN con N = "<<n<<endl<<endl;

		for(int i=0; i<Nthrows; i++)
		{
			SNu=0;
			SNexp=0;
			SNcauchy=0;
			for (int j=0; j<n; j++)
			{
				SNu+=rnd.Rannyu();
				SNexp+=rnd.Exponential(1);
				SNcauchy+=rnd.Cauchy(0,1);
			}
			
			fout1<<SNu/n<<endl;
			fout2<<SNexp/n<<endl;
			fout3<<SNcauchy/n<<endl;
		}	
	}

fout1.close();
fout2.close();
fout3.close();


////////////////////////////////////////////Buffon

double l=4.; //needle length
double d=8.; //lines space

int Nblocks=100;  
Nthrows=100000*Nblocks;
E=Nthrows/Nblocks;//throws per block
int Nhit=0;

vector <double> pi_block(Nblocks);//averages
vector <double> pi_block_2(Nblocks);//squared averages

double x_med=0;//we can characterize a needle with it's barycenter and an angle theta
double theta=0;//angle
double proj=0;//projection of needle on axis

	//we sample theta withotuh using pi. We are sampling theta in [0, Pi], though, due the model symmetry, is equivalent to sampling theta in [0,2*Pi]. 	
	for(int i=0; i<Nblocks; i++)
	{
		Nhit=0;
		for(int j=0; j<E; j++)
		{
			x_med=rnd.Rannyu(0,d);//seampling the barycenter

			theta=atan(abs(rnd.Rannyu()-rnd.Rannyu())/abs((rnd.Rannyu()-rnd.Rannyu())));//sampling theta
			proj=l*sin(theta)/2.;//projection of needle on x axis

			if((x_med + proj)>d or (x_med - proj)<0) 
			{
				Nhit++;
			}
		}
		pi_block[i]=(2*l*E)/(d*Nhit);
		pi_block_2[i]=pow(pi_block[i],2);
	}

	//reset alrady existing vectors
	sumprog.resize(Nblocks);
	sum2prog.resize(Nblocks);
	err_prog.resize(Nblocks);
	fill(sumprog.begin(), sumprog.end(), 0.0);
	fill(sum2prog.begin(), sum2prog.end(), 0.0);

	//staring to calculate running averages and errors
	partial_sum(pi_block.begin(), pi_block.end(), sumprog.begin());
	partial_sum(pi_block_2.begin(), pi_block_2.end(), sum2prog.begin()); 

	err_prog[0]=0;
	for(int i=1; i<Nblocks; i++)
	{
		sumprog[i]/=(i+1);
		sum2prog[i]/=(i+1);
		err_prog[i]=sqrt((sum2prog[i]-pow(sumprog[i],2))/i);
	}

fout.open("Buffon.txt");

	for (int i=0; i<Nblocks; i++ ) 
		fout<<sumprog[i]<<"	"<<err_prog[i]<<endl;

fout.close();

rnd.SaveSeed();
   
   return 0;
}



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
