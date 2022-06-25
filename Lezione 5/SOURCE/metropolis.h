#ifndef __Metropolis__
#define __Metropolis__

#include "random.h"
#include "random_walk.h"

class Metropolis : public RWalk
{

public:
  // constructors
  Metropolis(double a, int dim, double (*fun)(std::vector<double>), void(RWalk::*type_of_rw)());
  
  // destructor
  
  // methods
  std::vector<double> equilibrate(int Nequil);
  void tuning();
  void sample_point();
  
  //set get
  int get_accepted() {return accepted;}
  int get_rejected() {return rejected;}
  void set_initial_point(int i, double r);
  void set_pdf(double (*fun)(std::vector<double>));
  void set_Tyx( void(RWalk::*type_of_rw)() );

private:
	void (RWalk::*Tyx)(); // pointer to Rwalk data member of type void that takes nothing as input 
	double (*mypdf)(std::vector<double>); //pointer to function that take a vector and returns a double
	int accepted, rejected;
};

#endif // __Metropolis__

