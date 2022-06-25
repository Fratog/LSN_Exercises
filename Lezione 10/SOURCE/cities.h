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
	void cities_from_file(std::string filename);
	void print_coords() const;
	std::vector<std::vector<double>> get_coords();
	double get_coords(int i, int j);
	
	
private:
	int Ndim; 
	int Ncities;
	std::vector<std::vector<double>> coords; //ncities rows, ndim columns
	Random rnd;
};

#endif // __Cities__
