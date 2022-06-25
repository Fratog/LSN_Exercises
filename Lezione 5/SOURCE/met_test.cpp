#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <numeric>
#include "block_method.h"
#include "metropolis.h"

using namespace std;

double psi100(vector<double> r_vec) 
{
	return exp(-2*sqrt(inner_product(r_vec.begin(), r_vec.end(), r_vec.begin(), 0.0L)));
}


double psi210(vector<double> r_vec) 
{
	double r=sqrt(inner_product(r_vec.begin(), r_vec.end(), r_vec.begin(), 0.0L));
	return pow(r*exp(-r/2.)*r_vec[r_vec.size()-1]/r,2);
} 

int main(int argcv, char *argv[])
{

///////////////////////////////variables

//data blocking variables
int Nthrows=1000000;
int Nblocks=50;
int Ntxb=Nthrows/Nblocks;
Block_method bm;
vector<double> running_r(Nblocks);
vector<double> running_err_r(Nblocks);

//output files variables
ofstream pos; //points in space
ofstream rout; //r distances

//observables vector
vector<double> r_unif(Nblocks);
vector<double> r_gauss(Nblocks);

//metropolis variables
float psi100_max=100;
float step_unif=1.2;  
float step_gauss=0.8;
int dim=3;
Metropolis m_unif(step_unif, dim, psi100, &RWalk::move_uniform);//create metropolis objects 
Metropolis m_gauss(step_gauss, dim, psi100, &RWalk::move_gauss);

//equilibration variables
int Nequil=5000;
vector<double> r_equil(Nequil);

///////////////////////////////////////////100 orbital

////////////////////////////100 orbital equilibration 
	m_unif.set_initial_point(dim-1, psi100_max); //the first argument select the coordinate i of the random walk vector
	m_gauss.set_initial_point(dim-1, psi100_max);

	r_equil=m_unif.equilibrate(Nequil); //get a vector of Nequil values of r

	rout.open("equil_unif_psi100.txt");
	for(int i=0; i<Nequil; i++)
		rout<<r_equil[i]<<endl;

	rout.close();

	r_equil=m_gauss.equilibrate(Nequil); //get a vector of Nequil values of r

	rout.open("equil_gauss_psi100.txt");
	for(int i=0; i<Nequil; i++)
		rout<<r_equil[i]<<endl;

	rout.close();

///////////////////////////////100 get results
	pos.open("unif_pos100.txt");
	for(int i=0; i<Nblocks; i++)
	{
		for(int j=0; j<Ntxb; j++)
		{
			m_unif.sample_point();
			m_gauss.sample_point();
			
			if(j%100==0)
			{
				for(const auto & i : m_unif.get_coords()) 
					pos<<i<<"\t";
				pos<<endl;
			}
			
			r_unif[i]+=sqrt(m_unif.get_norm_squared());
			r_gauss[i]+=sqrt(m_gauss.get_norm_squared());
		}
		r_unif[i]/=Ntxb;	
		r_gauss[i]/=Ntxb;
	}
	pos.close();
	
	cout<<"uniform move, acceptance rate : "<<(double)m_unif.get_accepted()/(m_unif.get_accepted()+m_unif.get_rejected())<<endl;
	cout<<"gauss move, acceptance rate : "<<(double)m_gauss.get_accepted()/(m_gauss.get_accepted()+m_gauss.get_rejected())<<endl;

	running_r=bm.get_running_av(r_unif);//get running averages
	running_err_r=bm.get_running_err(r_unif); //get running errors

	rout.open("unif_psi100.txt");
	for (int i=0; i<Nblocks; i++)
		rout<<running_r[i]<<"\t"<<running_err_r[i]<<endl;

	rout.close();

	running_r=bm.get_running_av(r_gauss);
	running_err_r=bm.get_running_err(r_gauss);

	rout.open("gauss_psi100.txt");
	for (int i=0; i<Nblocks; i++)
		rout<<running_r[i]<<"\t"<<running_err_r[i]<<endl;

	rout.close();

/////////////////////////////////////210 orbital

/////////////resetting the parameters
double psi210_max=sqrt(2); 
	
	step_unif=2.9;
	step_gauss=1.8;
	Nequil=1000;

	m_unif.set_step_length(step_unif);//change step length to rwalk of metropolis object
	m_unif.move_to_origin();//reset to zero the rwalk coordinates vector
	m_unif.set_initial_point(dim-1, psi210_max);//choose initial position
	m_unif.set_pdf(psi210);//change pdf to metropolis object

	m_gauss.set_step_length(step_gauss);
	m_gauss.move_to_origin();
	m_gauss.set_initial_point(dim-1, psi210_max);
	m_gauss.set_pdf(psi210);

//////////////////////////////210 equilibration
	//here we have chosen the correct initial point, though i'll make it move for Nequil moves just to be sure.
	//r_equil=m_unif.equilibrate(Nequil); //get a vector of Nequil values of r

///////////////////////////////210 get results
	fill(r_unif.begin(), r_unif.end(), 0.0);
	fill(r_gauss.begin(), r_gauss.end(), 0.0);

	pos.open("unif_pos210.txt");
	for(int i=0; i<Nblocks; i++)
	{
		for(int j=0; j<Ntxb; j++)
		{
			m_unif.sample_point();//oss ho fatto l'algoritmo di metroplis in modo che restituisca il punto, ma questo non serve!!
			m_gauss.sample_point();
			
			if(j%100==0)
			{
				for(const auto & i : m_unif.get_coords()) 
					pos<<i<<"\t";
				pos<<endl;
			}
			
			r_unif[i]+=sqrt(m_unif.get_norm_squared());
			r_gauss[i]+=sqrt(m_gauss.get_norm_squared());
		}
		r_unif[i]/=Ntxb;	
		r_gauss[i]/=Ntxb;
	}
	pos.close();
	
	cout<<"uniform move, acceptance rate : "<<(double)m_unif.get_accepted()/(m_unif.get_accepted()+m_unif.get_rejected())<<endl;
	cout<<"gauss move, acceptance rate : "<<(double)m_gauss.get_accepted()/(m_gauss.get_accepted()+m_gauss.get_rejected())<<endl;

	running_r=bm.get_running_av(r_unif);//get running averages
	running_err_r=bm.get_running_err(r_unif); //get running errors

	rout.open("unif_psi210.txt");
	for (int i=0; i<Nblocks; i++)
		rout<<running_r[i]<<"\t"<<running_err_r[i]<<endl;

	rout.close();

	running_r=bm.get_running_av(r_gauss);
	running_err_r=bm.get_running_err(r_gauss);

	rout.open("gauss_psi210.txt");
	for (int i=0; i<Nblocks; i++)
		rout<<running_r[i]<<"\t"<<running_err_r[i]<<endl;

	rout.close();

return 0;
}
