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
#include <ostream>
#include <cmath>
#include <numeric>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(int argcv, char *argv[])
{ 
	Input(argcv, argv); //Inizialization
	for(int i=0; i<2000; i++) //equilibration
		Move(metro);
  
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);  //metropolis move
      Measure();		//measure on new state
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(int argcv, char *argv[]) 
{
  ifstream ReadConf, Seed;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;



//Read input informations
  //ReadInput.open("input.dat");


  //ReadInput >> temp;
  temp=atof(argv[1]);
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  //ReadInput >> nspin;
  nspin=atoi(argv[2]);
  cout << "Number of spins = " << nspin << endl;

  //ReadInput >> J;
  J=atof(argv[3]);
  cout << "Exchange interaction = " << J << endl;

  //ReadInput >> h;
  h=atof(argv[4]);
  cout << "External field = " << h << endl << endl;
    
  //ReadInput >> metro; // if=1 Metropolis else Gibbs
	metro=atoi(argv[5]);
  if(metro==1)
		algo="M";
  else
  	algo="G";
  	
  //ReadInput >> nblk;
	nblk=atoi(argv[6]);
	
  //ReadInput >> nstep;
	nstep=atoi(argv[7]);
	
	//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

	restart=atoi(argv[8]);
  if(restart) Seed.open("seed.out");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();
	
  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  //ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

	//output file names relevant
	std::stringstream ss;
	ss << std::fixed << std::setprecision(2) << temp;
	std::string t = ss.str();

//initial configuration

  if(restart)
  {
    ReadConf.open("config.final_"+ t + "_" + algo);
    for(int i=0; i<m_spin; ++i) ReadConf >> s[i];
  }
  else
  {
		for(int i=0; i<nspin; ++i)
		{
		  if(rnd.Rannyu() >= 0.5) s[i] = 1;
		  else s[i] = -1;
		}
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;

}


void Move(int metro)
{
  int o; //randomly selected spin
  double p, energy_old, energy_new, sm; //met variables 
  double energy_up, energy_down; //gibbs variables

  for(int i=0; i<nspin; ++i) //one moves= one try on every spin, it is not better to select randomly one at a time???
  {

    o = (int)(rnd.Rannyu()*nspin); //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
		attempted++;
		
    if(metro==1) //Metropolis
    {
			energy_old=Boltzmann(s[o], o); //old energy
			sm=-1*s[o]; //flip the spin
			energy_new=Boltzmann(sm, o); //new energy
			p=exp(-beta*(energy_new-energy_old)); //calculate transition porbability
				
			if(p>1) 
			{
				s[o]*=-1;		//accept spin flip							
				accepted++;
			}
			else if(rnd.Rannyu()<=p)
			{
				s[o]*=-1;		//accept spin filp							
				accepted++;
			}
			else ;
			
    }
    else //Gibbs sampling
    {
    	energy_up=Boltzmann(1, o);
    	energy_down=Boltzmann(-1, o);
    	p=1/(1+exp(-beta*(energy_up-energy_down)));//prob up, probdown=1-p
    	if(rnd.Rannyu()<=p)
    	{
    		s[o]=-1;
    		accepted++; //gibbs always accepts, it is trivial
    	}
    	else 
    	{
    		s[o]=+1;	
    		accepted++;
    	}	
    }
  }
}

double Boltzmann(int sm, int ip) 
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
 //s[Pbc(ip-1)] + s[Pbc(ip+1)] are the nerest neighbours to the selected spin indexed by ip, which will correspond to o in move i guess. sm is the
 // sm is the value of the spin o (the one selected randomly). so ene gives the energy part that will differ from the flipped state.
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0;
	
	//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
		u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
		m += s[i];		
  }
  walker[iu] = u; //it doesn't make sense to divide now for nspin, we can do it only one time later
  walker[ic] = u*u; 
  
  walker[im] = m;
  walker[ix] = m*m; // oss at h=0 we must not put the <m> since it shoud be exactly =0. (no spontaneus magnetization)

}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)// if it is the first block, make sure they are 0
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i) //this needs to be 0 at every new block
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0; //imp: should i have 50 % at every block?? yes is normal, since at every block i want to have an estimate of the quantitiy to calculate with the Metropolis
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
		ofstream Ene, Heat, Mag, Chi, file_names;
		const int wd=12;
		string fname;
		cout<<temp<<endl;

  //avoid trailing zeroes that to_string in double case would gvive 		 
		std::stringstream ss;
		ss << std::fixed << std::setprecision(2) << temp;
		std::string t = ss.str();

		
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    if(h==0)
    {
			fname="Energy_" + t + "_" + algo + ".txt";
			Ene.open(fname,ios::app);
			stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
			glob_av[iu]  += stima_u;
			glob_av2[iu] += stima_u*stima_u;
			err_u=Error(glob_av[iu],glob_av2[iu],iblk);
			Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
			Ene.close();
			
			fname="Heat_" + t + "_" + algo + ".txt";
			Heat.open(fname,ios::app);
			stima_c = pow(beta,2)*(blk_av[ic]/blk_norm-pow((stima_u*((double)nspin)),2))/(double)nspin; 
			glob_av[ic] += stima_c;
			glob_av2[ic] += stima_c*stima_c;
			err_c=Error(glob_av[ic], glob_av2[ic], iblk);
			Heat<< setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
			Heat.close();
			

			fname="Chi_" + t + "_" + algo + ".txt";
			Chi.open(fname,ios::app);
			stima_x = (beta)*blk_av[ix]/blk_norm/(double)nspin; 
			glob_av[ix] += stima_x;
			glob_av2[ix] += stima_x*stima_x;
			err_x=Error(glob_av[ix], glob_av2[ix],iblk);
			Chi<< setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
			Chi.close();
		}	

		if(h!=0) 
		{
		  fname="Mag_" + t + "_0.02_" + algo + ".txt"; //i'll use only h=0.02 or 0 for now
		  Mag.open(fname,ios::app);
		  stima_m = blk_av[im]/blk_norm/(double)nspin; //Energy
		  glob_av[im] += stima_m;
		  glob_av2[im] += stima_m*stima_m;
		  err_m=Error(glob_av[im],glob_av2[im],iblk);
		  Mag<< setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
		  Mag.close();
    }
    
    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;
		//output file names
	std::stringstream ss;
	ss << std::fixed << std::setprecision(2) << temp;
	std::string t = ss.str();

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final_"+ t + "_" + algo);
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed(); //imp
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
