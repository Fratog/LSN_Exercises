#ifndef __Cities__
#define __Cities__

#include<armadillo>
#include "random.h"

class Cities 
{

public:
  // constructors
  Cities(int, int, int);
  
  // destructor
  //~Cities();
  
  // methods
	void cities_on_hypercube(double);
	void cities_on_circle(double);
	void print_coords() const;
	std::vector<std::vector<double>> get_coords();
	double get_coords(int i, int j);
	
	
private:
	int Ndim; 
	int Ncities;
	std::vector<std::vector<double>> coords; //ncities rows, ndim columns
	Random rnd;
	//arma::Mat<double> coords;//, arma::fill::zeros);
};

#endif // __Cities__
