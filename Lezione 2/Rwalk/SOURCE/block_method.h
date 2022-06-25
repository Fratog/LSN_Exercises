#ifndef __Block_method__
#define __Block_method__

class Block_method
{
	public:
	Block_method();

	void set_count(int i); 
	const int get_count(); 
	double Raverage_predicate(double x, double y); 
	double get_mean(std::vector <double> A); 
	double get_error(std::vector <double> A); 

	private:
	int count;
	int Nblocks;
};


#endif // __Block_method__
