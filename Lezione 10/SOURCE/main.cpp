#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include "random.h"
#include "chromosome.h"
#include "population.h"
#include "cities.h"
#include "mpi.h"

//#include "random.h"

using namespace std;

int main(int argc, char **argv)
{

int size, rank;
MPI_Init(&argc,&argv); //Initialize mpi
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//GA variables
int row1=atoi(argv[2]);//save integer for setting cities random generator
int row2=atoi(argv[3]);//save integer for setting populatin random generators, (there will be even a rank dependece)
const int Ncities=atoi(argv[4]);//number of cities
int select=atoi(argv[5]);//index to choose among standard selection operator (2), and SA selection operator (1)
int dim=2;
int Npop=200;
int start=0;//begin of cities path
int N_gen=700;
double L=1;

//parallelization variables
int Nmigr=20; //migration every 20 new generations
vector<int> exchange_index={1};//prima facevo 1,2,3,5,10 ma forse è troppo cosi non sto scmabiando informazione, ma la sto spostando e basta
int index1=0;//index of first core selected randomly for information exchange
int index2=0;//second index

int itag1=1;
int itag2=2;
MPI_Status stat1, stat2;
MPI_Request req;


Cities USA(Ncities, dim, row1);
	USA.cities_on_hypercube(L);
	//USA.cities_from_file("usa_only_coords");
	//USA.print_coords();
	
Population pop(Npop, start, USA.get_coords(), row2+10*rank);//every rank has it's own pop with different seeds, as rank changes which row of primes passed to the generator;
	

	for(int i=0; i<N_gen; ++i)
	{
		pop.new_generation(USA.get_coords(), select);		
		if(atoi(argv[1])) //1 parallelization with comuication, 0 trivial parallelization
		{
			for(const auto & j : exchange_index) 
			{			
				if(i%Nmigr==0 and i!=0)//every Nmigr new generation do the exchange
				{
					//if(rank==0)
						//cout<<"in Nmigr if"<<i<<endl;
					
					//if i use the random gen in pop, all processes will have different values of index1 index2 which won't work if you look how send and recevie are used. index1 and index 2 must be common to all processes. a simple but non elegant way to solve this would be to use the USA random gen, since the seed it's common to all processes will give same indexes. the probably more correct way is to use collective comuniactions
					if(rank==0) //rank =0 slectes the cores who will exchange information
					{
						index1=int(pop.get_gen().Rannyu(0,size));
						do
						{
							index2=int(pop.get_gen().Rannyu(0,size));
						}
						while(index2==index1);
					}

					MPI_Bcast(&index1,1,MPI_INTEGER,0, MPI_COMM_WORLD);//rank 0 broadcast the index results to all the others
					MPI_Bcast(&index2,1,MPI_INTEGER,0, MPI_COMM_WORLD);
					
					//cout<<"rank	"<<rank<<"	index1	"<<index1<<endl;
					//cout<<"rank	"<<rank<<"  index2	"<<index2<<endl;
									
					//pop.get_chrom_i_ref(0).check();//not necessary
					auto msg1=&(pop.get_chrom_i_ref(j).get_gene_i_ref(0));//saving the message to exchange
					auto msg2=&(pop.get_chrom_i_ref(j).get_gene_i_ref(0));
					
					//int * msg1=&(pop.get_chrom_i_ref(0).get_gene_i_ref(0));
					//int * msg2=&(pop.get_chrom_i_ref(0).get_gene_i_ref(0));
					
					///////////////OSS IF I USE MPI TO EXCHANGE ONLY THE PATH AND NOT THE CHROMOSOME, THEN I HAVE TO UPDATE THE FIT OF THE RELATIVE CHROMOSOME AND RUN A SORT FOR ONLY FEW ELEMENTS
					//it would be nice to create our own data type chromosme to use send and receive directly on those obejcts
					if(rank==index1)
					{
						MPI_Isend(msg1, Ncities, MPI_INTEGER, index2, itag1 ,MPI_COMM_WORLD, &req);
						//MPI_Send(msg1, Ncities, MPI_INTEGER, index2, itag1 ,MPI_COMM_WORLD);
						MPI_Recv(msg2, Ncities, MPI_INTEGER, index2, itag2, MPI_COMM_WORLD, &stat2);
						pop.get_chrom_i_ref(j).calculate_fit(USA.get_coords());
						pop.sort_pop();//maybe if I use the SA selection this isn't necessary
					}
					else if(rank==index2)
					{
						MPI_Send(msg2, Ncities, MPI_INTEGER, index1, itag2, MPI_COMM_WORLD);
						MPI_Recv(msg1, Ncities, MPI_INTEGER, index1, itag1, MPI_COMM_WORLD, &stat1);
						pop.get_chrom_i_ref(j).calculate_fit(USA.get_coords());
						pop.sort_pop();
					}
				}	
			}		
		}
	}

	pop.sort_pop(); //probably not necessary, but i want to be sure
	
ofstream out;
double min = 0;
//double red_buff=pop.get_fits()[0];
double red_buff=pop.get_chrom_i_ref(0).get_fit();

//	cout<<"rank "<<rank<<"	fit	"<<red_buff<<endl;    
	MPI_Reduce(&red_buff, &min, 1 , MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);//find the minimum fit among all the cores

	if(rank==0)
	{
		
		if(select==1)
			//out.open("fit_8cores_SA.txt", ios::app);
			out.open("fit_varying_ncores_" + to_string(Ncities) + "_SA_nonpar.txt", ios::app);
		else if(select==2)
			//out.open("fit_8cores_Standard.txt", ios::app);
			out.open("fit_varying_ncores_" + to_string(Ncities) + "_Standard_nonpar.txt", ios::app);
		else
			cout<<"error in opening output file"<<endl;
			
		out<<(size)<<"\t"<<min<<endl;
		cout<<(size)<<"\t"<<min<<endl;
		out.close();	
		
		if(select==1)
			cout<<"you are using the SA selection operator"<<endl;
		else if(select==2)
			cout<<"you are using the standard selection operator"<<endl;
		else
			cout<<"selection operator index is wrong"<<endl;
	}


	//cout<<"node number	"<<rank<<"	best fit	"<<pop.get_fits()[0]<<endl;	
MPI_Finalize();

return 0;
}
