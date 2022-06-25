#ifndef __Block_method__
#define __Block_method__




class Block_method
{
public:
//constrcutor
	Block_method(); 	
	Block_method(int i, int j);

//methods
	void set_count(int i); 
	void set_Nblocks(int i);
	const int get_count();
	const int get_Nblocks(); 
	std::vector<double> get_running_av(std::vector <double> A);
	std::vector<double> get_running_err(std::vector <double> A); 
	double get_mean(std::vector<double> A); 
	double get_error(std::vector<double> A); 


private:
	int count;
	int Nblocks;
};

#endif // __Block_method__
