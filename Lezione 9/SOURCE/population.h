#ifndef __Population__
#define __Population__

//#include<armadillo>
//#include <set>
#include "chromosome.h"
#include "cities.h"

class Population 
{

public:
  // constructors
	Population(int N, int start); 
  Population(int N, int start, const std::vector<std::vector<double>> &, int primes_row);
  
  // destructor
  ~Population();

	//GA algoritm methods and operators
  bool check_population() const;
  void sort_pop(); //on a fitness basis
  void sort_fits();
  Chromosome selection1();
	Chromosome selection2(double p); 
  void crossover(Chromosome &, Chromosome &);
  void new_generation(const std::vector<std::vector<double>> &, int);
  
  //set get print methods
  std::vector<Chromosome> get_pop_fail() const;
  Chromosome get_chrom_i(int) const;
  void print_fits() const;
  void set_fit(int);
  std::vector<double> get_fits() const;
  Random & get_gen();
	
private:
	
	int Nchrom;//number of chromosomes
	int start;//cities to start with. Not really necessary
	std::vector<Chromosome> Pop;
	std::vector<double> fits; //vector of fits, associated to the population. I don't like this
	Random rnd;
	
	Chromosome & get_chrom_i_ref(int i);
	std::vector<Chromosome> & get_pop(); 
};

#endif // __Population__int N, int start
