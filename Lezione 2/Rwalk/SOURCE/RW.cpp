#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "random.h"
#include "block_method.h"
#include "random_walk.h"

using namespace std;


//////////////////////////////MAIN

int main(int argcv, char *argv[])
{

int Nstep=atoi(argv[1]);//number of rwalk steps
int Nsim=atoi(argv[2]);//number of rw lengths to sample for each step
int Nblocks=atoi(argv[3]);//number of blocks for data blocking
int dim=atoi(argv[4]);//dimnesion of space in which rwalk lives

int Nthrows=Nsim/Nblocks;
vector<double> Rn(Nstep);//average distance vector
vector<double> RnErr(Nstep);//
vector<double> Rn_per_block(Nblocks);

RWalk rw(1,dim);
Block_method bm;

cout<<dim<<endl<<endl;

////////////////////////////////////DISCRETE RWALK  
for(int i=0; i<Nstep; i++)
{
	fill(Rn_per_block.begin(), Rn_per_block.end(),0);
	for(int j=0; j<Nblocks; j++)
	{
		for(int l=0; l<Nthrows; l++)
		{			
			rw.RW_move_discrete(i);
			Rn_per_block[j]+=rw.get_norm_squared();
			rw.move_to_origin();
		}

		Rn_per_block[j]/=Nthrows;	
	}

	Rn[i]=sqrt(bm.get_mean(Rn_per_block));	
	RnErr[i]=bm.get_error(Rn_per_block);
	RnErr[i]*=0.5/sqrt(Rn[i]);//error propagation
}	

ofstream fout;
string filename="Rw_discrete";
fout.open (string(filename) + to_string(dim) + ".txt");

	fout<<"here there are the progressive distance avaraged over "<<Nsim<<" Random walks and the relative errors for Nstep= "<<Nstep<<" in "<<dim<<" dimensions. DISCRETE CASE"<<endl; 


	RnErr[0]=0;
	fout<<Rn[0]<<"\t"<<RnErr[0]<<endl;//perchè l'ho messo fuori?
		for (int i=1; i<Nstep; i++)
	fout<<Rn[i]<<"\t"<<RnErr[i]<<endl;

fout.close();
////////////////////////////CONTINUOUS RWALK 

if(dim==3)
{
	for (int i=0; i<Nstep; i++)
	{
		fill(Rn_per_block.begin(), Rn_per_block.end(),0);
		for (int j=0; j<Nblocks; j++)
		{
			for(int l=0; l<Nthrows; l++)
			{			
				rw.RW_move_continuos(i); 
				Rn_per_block[j]+=rw.get_norm_squared();
				rw.move_to_origin();	
			}
			Rn_per_block[j]/=Nthrows;	
		}
		
		Rn[i]=sqrt(bm.get_mean(Rn_per_block));		
		RnErr[i]=bm.get_error(Rn_per_block);
		RnErr[i]*=0.5/sqrt(Rn[i]); //error propagation
	}	

	fout.open ("Rw_continuos.txt");
		fout<<"here there are the progressive distance avaraged over "<<Nsim<<" Random walks and the relative errors for Nstep= "<<Nstep<<" in "<<dim<<" dimensions. CONTINUOS CASE"<<endl; 


		RnErr[0]=0;
		fout<<Rn[0]<<"\t"<<RnErr[0]<<endl;
		for (int i=1; i<Nstep; i++)
			fout<<Rn[i]<<"\t"<<RnErr[i]<<endl;

	fout.close();
}
return 0;
}


