/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __ISING__
#define __ISING__

//Random numbers
//#include <sstream>
#include <iomanip>
#include <string>
#include "random.h"
#include <vector>

int seed[4];
int restart;
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props,iu,ic,im,ix,ig; //ig???
double nbins;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_u,stima_c,stima_m,stima_x,stima_g;
double err_u,err_c,err_m,err_x,err_g;

//configuration
std::string algo;
const int m_spin=50;
std::vector<double> s(m_spin);  

// thermodynamical state
int nspin;
double beta, temp , J , h;

// simulation
int nstep, nblk, metro;



//functions
void Input(int argcv, char *argv[]);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(int);
void ConfFinal(void);
void Measure(void);
double Boltzmann(int, int);
int Pbc(int);
double Error(double,double,int);
void Equilibrate(int Nstep);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
