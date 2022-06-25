#ifndef __RWalk__
#define __RWalk__



class RWalk
{
	public:

	RWalk(double a, int d);

	std::vector<double> get_coords();
	std::vector<double> get_step();
	double get_norm_squared();

	void move_to_origin();
	void RW_move_discrete();
	void RW_move_discrete(int nstep);
	void RW_move_continuos(); //should be used only for dim=3 case, it is not general
	void RW_move_continuos(int nstep);
	void RW_move_back(); 


	private:
	int dim;
	std::vector<double> step;
	std::vector<double> coords;
	std::vector<double> last_coords;
	Random rnd;
	
};

#endif // __RWalk__
