#ifndef __Metropolis__
#define __Metropolis__

#include "random.h"
#include "random_walk.h"

class Metropolis : public RWalk
{

public:
  // constructors
  Metropolis(double a, int dim, void(RWalk::*type_of_rw)());
  Metropolis(double a, int dim, double (*fun)(std::vector <double>), void(RWalk::*type_of_rw)());
  
  // destructor 
  ~Metropolis();
  
  // methods
  void equilibrate(int Nequil);
  void find_step(double step_in, int Nequil);
  void sample_point();
  bool sample_point(double xnew, double xold);
  
  //set get
  int get_accepted() {return accepted;}
  int get_rejected() {return rejected;}
  void set_initial_point(int i, double r);
  void set_pdf(double (*fun)(std::vector <double>));
  void set_Tyx(void(RWalk::*type_of_rw)());
  void set_step(double step);
	
protected:
	void (RWalk::*Tyx)();
	
private:
	double (*mypdf)(std::vector <double>); 
	int accepted, rejected;
};

#endif // __Metropolis__

