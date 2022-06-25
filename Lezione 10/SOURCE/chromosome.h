#ifndef __Chromosome__
#define __Chromosome__

//#include<armadillo>
#include "random.h"
#include <vector>


class Chromosome { //should i do it : public cities, so i have the cities vector??

public:
  // constructors
  Chromosome();
  Chromosome(int N, int start, Random & rnd);
  
  // destructor
  ~Chromosome();
  
  //copy and move constructors/assginment op
  Chromosome(const Chromosome &);// = default;
  Chromosome & operator=(const Chromosome &); // = default;
  Chromosome(Chromosome &&);// = default;
  Chromosome & operator=(Chromosome &&);// = default;
  
  //operator overloading-to order on a fitness basis
	bool operator<(const Chromosome &) const;//to sort population by fit

  //methods

  //set get print
  std::vector<int> & get_chromosome_ref();
  std::vector<int> get_chromosome() const;
  int & get_gene_i_ref(int);
  void print() const;
  double get_fit() const;
	
  //others
  void calculate_fit(const std::vector<std::vector<double>> &);
  bool check() const;
  bool check_2() const; 
  void print_config(std::string , const std::vector<std::vector<double>> &);
  
  //mutations 
  void mutation1(Random &);
  void mutation2(Random &);
  void mutation3(Random &);
  void mutation4(Random &);
  void mutation5(Random &);
	
	//crossover
	void crossover(Chromosome & Chrom1, Random & rnd);

private:	
	
	int Ngenes; 
	std::vector<int> Chrom;
	double fit;
};

#endif // __Chromosome__


