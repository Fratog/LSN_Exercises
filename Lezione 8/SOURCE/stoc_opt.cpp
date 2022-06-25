/***************************
With this program we want to find the values mu sigma which characterize the ground state wavefunction of our hamiltonian.
We approach the problem with simulated annealing.
The program works using two Metropilis algorithms, one in the VMC object, used to sample the trial wavefunction. Another one in the SA part, applied to the trial energies obtained by the VMC code, while moving mu sigma with a random walk.
We will refer to the diffrent metropolis as VMC Metropolis and SA metropolis.
***************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include "random.h"
#include "statistics.h"
#include "random_walk.h"
#include "vmc.h"

#include <thread>
#include <chrono>

using namespace std;

///choose the starting point in a clever way
double mu=1; //1
double sigma=1; //1

double psiT2(vector<double> r_vec); //we are in 1-dim case
double HT_psiT(vector<double> r_vec); 

double oscill(vector <double> r_vec); //mu =0 sigma =1
double HT_oscill(vector <double> r_vec); 


int main(int argcv, char *argv[])
{

////////////////////////////VMC, here we build the part of the code that allows to get estimates of the energy. We test it on the quantum harmonic oscillator

cout<<"--------------VMC testing start-----------"<<endl;

//VMC class variables
double a=1; // VMC Mettropolis step, is the step used by the rwalk needed to sample Ht
int dim=1; //dimension of VMC Rwalk
int Nequil=1000; //numebr of metropolis step to find acceptance rate with find delta


//we change mu and sigma to work on harmonic oscillator
mu=0;
sigma=1;
VMC testing(a, dim, oscill, &RWalk::move_uniform, HT_oscill);


//statistics variables, autocorrelation and data blocking
int Nist_val=10000; //number of istantaneus vales for autocorrrelation calculation
int Nac_step=1000; //number of autocorellation steps (i.e 1000 values of Ac)
int Ht_Nblock=100; //number of blocks for plotting Ht whith mu=0 sigma=1
int Ht_Nstep=100000; //it is Nstep per block
vector<double> ist_values(Nist_val); //istant values vector

//observbales variables
vector<double> Ht(Ht_Nblock); //to save the vector with the block averages
vector<double> running_Ht(Ht_Nblock);//tu save the running average of the previous vector
vector<double> running_err_Ht(Ht_Nblock);//the same for the errors

ofstream out;

	testing.set_initial_point(0, mu);//start the VMC random walk from a clever point
	//opt.equilibrate(Nequil); //not really necessary
	cout<<"step before find_step call "<<testing.get_step_length()<<endl;
	testing.find_step(1, Nequil); //automatically adjust the step, at the given mu sigma
	cout<<"step after find_step call "<<testing.get_step_length()<<endl;

	//testing if the metropolis samples correctly the squared modulus of the HO ground state
	out.open("VMC_test_sampling.txt");
	for(int i=0; i<10000; i++)
	{
	testing.sample_point();
			for(const auto & j : testing.get_coords())
			{
				out<<j<<"\t";
			}
			out<<endl;
	}	
	out.close();	

  //here i study the autocorrelation,
	for (int i=0; i<Nist_val; i++)
	{
		ist_values[i]=testing.sample_cost();
	}				
	
	autocorrelation(ist_values, Nac_step, "VMC_E_AutoCorr.txt"); //correalation goes to zero very fast, so no problem with data blocking.
	
	//now i write a part of the program to get the cumulative value of the energy to check wether the variational monte carlo works
	Ht=testing.get_cost_vector(Ht_Nblock, Ht_Nstep);//get the blocking averages vector, (cost is meant as cost function, not a great name choice)
	running_Ht=running_av(Ht);
	running_err_Ht=bm_running_err(Ht);
 
	out.open("VMC_test_E0.txt");//to save the running averages and errors
 
	for(int i=0; i<Ht_Nblock; i++)
		out<<i<<"\t"<<running_Ht[i]<<"\t"<<running_err_Ht[i]<<endl;
	out.close();
	
cout<<"--------VMC testing end-------------"<<endl;

/////////////////////////////////Simulated Annealing, here we implement the other part of the exercisation

cout<<"--------SA starting-------------"<<endl;

//reseting the VMC object 
a=1;
mu=1;
sigma=1;
VMC opt(a, dim, psiT2, &RWalk::move_uniform, HT_psiT);

opt.set_initial_point(0, mu);
cout<<"step before find_step call "<<opt.get_step_length()<<endl;
testing.find_step(1, Nequil); //automatically adjust the step, at the given mu sigma
cout<<"step after find_step call "<<opt.get_step_length()<<endl;


//simulated annealing variables
double T=50.;
double beta=1/T;//intial temperatur, oss T intial must be st: T>>H_typical
int NSA_step=1000;// Number of SA metropolis steps

//statistics variables
int Nblocks=20; //each measure of Ht to pass to the annealing metropolis will be taken with this Nblocks and Nstep per block
int Nsteps=10000; //it is a nstep per block

//rwalk variabeles for the SA metropolis
int Nparam=2; //mu sigma
vector<double> delta={2.5, 1.7, 1.3, 0.7, 0.5, 0.35, 0.18, 0.1, 0.05, 0.01};//those are the values for the steps of SA random walk. i.e mu and sigma random walk. 
//I am choosing the delta values by direct observation. There isn't  an algorithm to adapt automatically the step for  the SA metropolis in order to satisfy the 50% rule. While there is one for the VMC metropolis, implemented in metropolis class.

double rate=0;//to save acceptance rate of SA metrpolis
RWalk rw(delta[0], Nparam);//mu sigma random walk
rw.set_i_coord(0,mu);//start parameters random walk from the initital point choosen in the global variables
rw.set_i_coord(1,sigma);

//other variables for the SA metropolis 
Random rnd; //SA metropolis random generator
vector<double> H_old(Nblocks); //it will be a vector of nblocks which containes the averages per block
vector<double> H_new(Nblocks);
vector<double> res_new(2);//it containts from H_new the estiamte of the energy and the error (with data blocking), this res_new/old will be the SA metropolis variables
vector<double> res_old(2);//the same for H_old

int accepted=0;
int rejected=0;
double alpha=0;
int counter=0;

//observable variables
vector<double> mu_vec;//history of mu parameter evolution
vector<double> sigma_vec;//history of sigma parameter evolution
vector<double> H_vec; 
vector<double> H_err_vec;

ofstream min;
ofstream param;
ofstream t_law;

	min.open("SA_Ht_istant.txt");
	param.open("SA_parameters.txt");
	t_law.open("SA_geometric_T.txt");

	opt.find_step(opt.get_step_length(), Nequil); 
	//cout<<opt.get_step_length()<<endl;	

	H_old=opt.get_cost_vector(Nblocks, Nsteps);//saves a vector of Nblocks
	res_old=opt.get_cost_value(H_old);// gets the mean and standard deviation of the mean of the passed vector

//Simulated annealing algorithm, I have choosen to not equilibrate the system at each temperature with N_Met steps, to avoid having to understand how many metropolis steps are necessary to equilibrate at varying temperature.  Instead i vary the temperature slowly, so that the system is always equilibrated by continuing from the previous SA step. And we perform many NSA_step to reach T->infinite even with our slow varying law for the temperature. 
	for(int i=0; i<NSA_step; i++)
	{
		cout<<"-------------------------------------------------"<<endl;
		cout<<"Sim Annealing step number "<<i<<endl;
		
		//mu and sigma rwalk
		rw.move_uniform(); //move mu sigma random walks
		mu=rw.get_coords()[0];//save new parameters, they are global variables
		sigma=rw.get_coords()[1];//save new parameters
		param<<mu<<"	"<<sigma<<endl;//write their values on file
		mu_vec.push_back(mu);//for later use
		sigma_vec.push_back(sigma);
		
		//if(i<500)//when the parameters change slowly it is better to keep it fixed 
		opt.find_step(opt.get_step_length(), 1000);//as we move mu and sigma, we need to adapt the step of the VMC random walk in orde to keep 50% rate
		cout<<"vmc metropolis step length	"<<opt.get_step_length()<<endl;
		
		
		H_new=opt.get_cost_vector(Nblocks, Nsteps);//get new data blocking vector for the energy estimate
		counter=0;
		do 
		{
			res_new=opt.get_cost_value(H_new);//get the energy estimate and it's error starting from H_new vector
			counter++;
			
			if(counter>50) {break;}
			if(compatible(res_new, res_old))//check if the measures are compatible
				opt.add_block(H_new, Nsteps);//this modifies the length of H_new to add a block and lower the error but then we must do a resize of H_new
			
			//cout<<"counter	"<<counter<<endl;
			//cout<<"H values	"<<res_new[0]<<"	"<<res_old[0]<<endl;
			//cout<<"check discrepancy"<<res_new[0]-res_old[0]<<"	"<<res_new[1]+res_old[1]<<endl;
			//cout<<"compatible??	"<<compatible(res_new, res_old)<<endl;
		}
		while(compatible(res_new, res_old));
		
		cout<<"counter	"<<counter<<endl;//checking when break enters
		H_new.resize(Nblocks);
		
		//Simulated annealing metropolis applicated on res_new/old
		alpha=exp(-beta*(res_new[0]-res_old[0]));
		if(alpha>1)
		{
			cout<<"beta	"<<beta<<"	Hnew	"<<res_new[0]<<endl;
			res_old=res_new;
			min<<res_new[0]<<"	"<<res_new[1]<<endl;
			accepted++;
		}
		else if(rnd.Rannyu()<=alpha)  
		{
			cout<<"beta	"<<beta<<"	Hnew	"<<res_new[0]<<endl;
			res_old=res_new;
			min<<res_new[0]<<"	"<<res_new[1]<<endl;
			accepted++;
		}
		else
		{
			cout<<"rejeceted"<<endl;
			rw.move_back();
			min<<res_old[0]<<"	"<<res_old[1]<<endl;
			//res_sold remains the same
			rejected++;
		}
		
		H_vec.push_back(res_old[0]);//saving the energy found
		H_err_vec.push_back(res_old[1]);

		t_law<<beta<<endl;
		beta/=0.99;	//temperature update
		cout<<"beta-> infinite	"<<beta<<endl;
	
		//Here we adapt the SA random walk step_length by using the delta vector, which must be created directly checking the acceptance rate printed out here in a trial run.
		if(i!=0 and i%100==0)
		{
			cout<<"*******************************************"<<endl;
			cout<<"here is the delta used for the SA metropolis for the previous 100 SA steps, and the relative acceptance rate"<<endl;
			cout<<int(i*10/NSA_step)<<endl;//delta index
			cout<<"delta	"<<delta[int(i*10/NSA_step-1)]<<endl;	//value of the delta
			rate=(double)accepted/(accepted+rejected);//acceptance rate of the related delta
			cout<<"rate	"<<rate<<endl;
			cout<<"*******************************************"<<endl;
			
			rw.set_step_length(delta[int(i*10/NSA_step)]);//modify the step with new delta
			std::this_thread::sleep_for(3s);
			accepted=0;//reset
			rejected=0;
		
		}
		
	}

t_law.close();
min.close();
param.close();
cout<<"--------SA end-------------"<<endl;

//question1: Show a picture of ⟨𝐻̂ ⟩𝑇 (with statistical uncertainties) as a function of the SA steps of the algorithm

	//resetting data blocking variables
	Ht_Nblock=50;
	Ht_Nstep=int(NSA_step/Ht_Nblock);
	Ht.resize(Ht_Nblock);
	running_Ht.resize(Ht_Nblock);
	running_err_Ht.resize(Ht_Nblock);

	//make data blocking starting from istant values, using blocking from statistics
	Ht=blocking(H_vec, Ht_Nblock, Ht_Nstep);//get blocking estimates from instantenous values 
	running_Ht=running_av(Ht);
	running_err_Ht=running_err(Ht);

	
out.open("SA_Ht_blocking.txt");
	
	for (int i=0; i<Ht_Nblock; i++)
		out <<i<<"\t"<<running_Ht[i]<<"\t"<<running_err_Ht[i]<<endl;

out.close();
	
//question 3 show a picture of the estimation of ⟨𝐻̂ ⟩𝑇 and its statistical uncertainty as a function of the number of blocks/MC steps for the set of parameters which minimize ⟨𝐻̂ ⟩𝑇
//now we should have our minimum parameters
int min_index=0;
double min_mu, min_sigma;

	min_index=distance(H_vec.begin(), min_element(H_vec.begin(),H_vec.end()));//find index of energy min
	min_mu=mu_vec[min_index];
	min_sigma=sigma_vec[min_index];

	cout<<"Hmin"<<*min_element(H_vec.begin(),H_vec.end())<<endl;
	cout<<"min index	"<<min_index<<"	mu_min	"<<min_mu<<"	sigma_min	"<<min_sigma<<endl;
	cout<<"min index	"<<min_index<<"	H_min	"<<H_vec[min_index]<<"	error	"<<H_err_vec[min_index]<<endl;

	mu=abs(min_mu);//i setted the new parameters value. oss the function is simmetric for parity transformations on the paramters
	sigma=abs(min_sigma);

//now with the new parameters setted we can get the running average for the energy, just as at the start

	//resetting data blocking variables
	Ht_Nblock=100;
	Ht_Nstep=100000;
	Ht.resize(Ht_Nblock);
	running_Ht.resize(Ht_Nblock);
	running_err_Ht.resize(Ht_Nblock);
	
	Ht=opt.get_cost_vector(Ht_Nblock, Ht_Nstep);
	//for(auto i : Ht ) cout<<i<<endl;
	running_Ht=running_av(Ht);
	running_err_Ht=bm_running_err(Ht);
 
 	

	out.open("SA_Hmin_blocking.txt");
	for (int i=0; i<Ht_Nblock; i++)
		out <<i<<"\t"<<running_Ht[i]<<"\t"<<running_err_Ht[i]<<endl;

	out.close();
	
//question 4 show also a picture of the sampled |Ψ𝑇(𝑥)|2 by filling a histogram with the sampled configurations, moreover compare it with the analytic curve of |Ψ𝑇(𝑥)|2 and with the numerical solution obtained by transforming the Schrodinger equation into a matrix equation 
	int Npoints=100000;

	out.open("ground_state_pdf.txt");
	for(int i=0; i<Npoints; i++)//remember you have to do this with mu=min_mu...
	{	
		opt.sample_point();
		out<<opt.get_coords()[0]<<endl;//we are one dim, so opt.get_coords()[0]
	}	
	out.close();

return 0;
}



double psiT2(vector <double> r_vec) //pdf to sample //we are in 1-dim case
{
	return pow(exp(-pow((r_vec[0]-mu)/sigma,2)/2) + exp(-pow((r_vec[0]+mu)/sigma,2)/2),2);
}

double HT_psiT(vector <double> r_vec) //this is the function we evaluate to get the trial energy estimates
{
	double x=r_vec[0];
	double x1=pow(x-mu,2);
	double x2=pow(x+mu,2);
	double e1=exp(-x1/(2*sigma*sigma));
	double e2=exp(-x2/(2*sigma*sigma));
	
	return (-0.5*(-(e1+e2)/(sigma*sigma) + 1/pow(sigma,4)*(x1*e1 + x2*e2)) + (pow(x,4) - 5./2*pow(x,2))*(e1+e2))/sqrt(psiT2(r_vec));
}

double oscill(vector <double> r_vec) //mu =0 sigma =1
{
	return pow(exp(-pow((r_vec[0]-mu)/sigma,2)/2),2);
}

double HT_oscill(vector <double> r_vec) 
{
	double x=r_vec[0];
	double e1=exp(-pow(x,2)/2.);
	
	return (-0.5*(e1*(pow(x,2)-1)) + 0.5*pow(x,2)*e1)/e1;
}

