#ifndef __RWalk__
#define __RWalk__



class RWalk
{
public:

	RWalk(double a, int d);

	void set_dim(int d);
	void set_step_length(double a);
	void set_origin(std::vector<double>);
	void set_i_coord(int i, double x);

	int get_dim();
	double get_step_length();
	std::vector<double> get_coords();
	Random & get_gen(); 

	void print_coords();

	void move_to_origin();
	void move_discrete();
	void move_discrete(int nstep);
	void move_continuos();
	void move_continuos(int nstep);
	
	void move_uniform();
	void move_gauss();
	void move_back(); 
	
	double get_norm_squared();
	 

private:

	int dim;
	std::vector<double> step;
	std::vector<double> coords;
	Random rnd;
	std::vector<double> last_coords;
};

#endif // __RWalk__
