#ifndef __VMC__
#define __VMC__

#include "metropolis.h"

class VMC : public Metropolis
{

public:
  // constructors
	VMC(double a, int dim, double (*fun)(std::vector<double>), void(RWalk::*type_of_rw)());  
  VMC(double a, int dim, double (*fun)(std::vector<double>), void(RWalk::*type_of_rw)(), double (*cost_fun)(std::vector<double>));
  // destructor //lo devi implementare per mettere a 0 i puntatori, non puoi lasciare quello implicito
  ~VMC();
  
  // methods
  double sample_cost(); 
  std::vector <double> get_cost_vector(int Nblocks, int Nstep);
	void add_block(std::vector<double> & A, int Nstep);
	std::vector<double> get_cost_value(int Nblocks, int Nstep);  
  std::vector<double> get_cost_value(std::vector<double> A);
  
  //set get
	

private:

 //cost function
 double (*cost)(std::vector<double>); //IMP!!//for us this will be the Ht*psi/psi function
};

#endif // __VMC__

