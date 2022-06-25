#include <iostream>
#include <vector>
#include <string>
//#include <ranges>
#include <numeric>
#include <algorithm>
#include "random.h"
//#include "statistics.h"
#include "chromosome.h"
#include "population.h"
#include "cities.h"
#include <unordered_set>
#include <chrono>
//#include "random.h"

using namespace std;


int main(int argc, char **argv)
{

int row1=atoi(argv[1]);// row of primes for Cities random generator
int row2=atoi(argv[2]);//row of primes for Population random generator
int Ncities=atoi(argv[3]);//number of cities
int L=atoi(argv[4]);	//length of the square
int select=0;//index for selection operator, // 1 for SA, 2 for standard selctio operator
double av=0;
string type_op;

int start=0;
int Npop=200;
int dim=2;
int N_gen=700;


	//choose selection operator, //standard selction operator is the default choiche if not specified
	if(argc==5)
	{
		select=2; 
		type_op="standard";
	}
	else if(argc==6)
	{
		select=atoi(argv[5]);
		type_op="SA";
	}
	
	//code for TSP inside a square, we already checked it works on circle, and now we switch pernamently to the square.	
	Cities Europe(Ncities, dim, row1);
	Europe.cities_on_hypercube(L);

	Population pop(Npop, start, Europe.get_coords(), row2);
	
	ofstream out;	
	out.open("fit_square.txt");
	for(int i=0; i<N_gen; i++)
	{
		pop.new_generation(Europe.get_coords(), select);//the flag select determines which selection operator will be used
		
		auto fits_vec=pop.get_fits();
		av=reduce(fits_vec.begin(), fits_vec.begin()+fits_vec.size()/2., 0.0)/(fits_vec.size()/2.);
		out<<fits_vec[0]<<"		"<<av<<endl;
		
		if(i%10==0)
		{
		//given that Pop is sorted, the first element correspond to minimum fit
		//IF YOU WANT TO PRINT CONFIGURATIONS UNCOMMENT THE FOLLOWING LINE
		//pop.get_chrom_i(0).print_config("DATA/SQUARE/CONFIG/Config_" + to_string(i) + ".txt", Europe.get_coords()); 
		}
		
	}
	out.close();
	
	//attention make sure that before using ./run.sh inputfile there isn't already a bound_fixedL.txt/ bound_fixedL_SA.txt file in the same folder of the main.x.
	if(select==1)
	{
		cout<<"you are using the SA selection operator"<<endl;
		out.open("bound_fixedL_SA.txt", ios::app);
	}
	else if(select==2)
	{
		cout<<"you are using the standard selection operator"<<endl;
		out.open("bound_fixedL.txt", ios::app);
	}
	else
	{
		cout<<"selection operator index is wrong"<<endl;
	}
	
	
	out<<pop.get_fits()[0]<<endl;//we save the best fit of the whole process, for different executions of the program with different seeds, to get our estimate of the expected value of the distnace for N randomly distributed points inside a square with L=1. This way we can confront our results with the theoretical and experimental bound.
	out.close();
	


return 0;
}
