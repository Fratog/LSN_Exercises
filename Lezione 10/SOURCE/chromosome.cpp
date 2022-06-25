#include <iostream>
#include <fstream>
#include <cmath>
//#include <armadillo>
#include <numeric>
#include <algorithm>
#include <unordered_set>
//#include "random.h"
#include "chromosome.h"

using namespace std;

Chromosome::Chromosome()
{
	//cout<<"CTOR EMPTY"<<endl;
}

Chromosome::Chromosome(int N, int start, Random & rnd) : Ngenes(N) 
{
	vector<int> allele;
	int appo=0;
	
	Chrom.resize(Ngenes);//NCity 
	Chrom.front()=start; //first city fixed as home
	
	for(int i=0; i<Ngenes; i++) 
		allele.push_back(i); //Cities index, so 0,1,2,3,4,5,....N
	
	allele.erase(allele.begin()+start); //erase home index from possible values
	
	for(auto it=Chrom.begin()+1; it<Chrom.end(); ++it) 
	{
		appo=int(rnd.Rannyu()*allele.size());//extract int from 0 to allele.size() 
		*it=allele[appo]; 
		allele.erase(allele.begin()+appo);
	}
	
	this->fit=0;
	this->check();//every time a chromosome is built check will be applied
	//cout<<"CTOR"<<endl;
}

//having just a vector as data member which gets acces to a resource, it wasn't necessary to implement destructor, copy, move... there would be the default ones working fine, but I wanted to check how they worked.
Chromosome::~Chromosome()//vector works on heap memory by itself, i don't have to control it
{
	//cout<<"DTOR"<<endl;
}

Chromosome::Chromosome(const Chromosome & Ch1) : Ngenes(Ch1.Ngenes), Chrom(Ch1.Chrom), fit(Ch1.fit)
{
	//cout<<"COPY"<<endl;
}

Chromosome& Chromosome::operator=(const Chromosome & Ch1) 
{
	//cout<<"COPY ASSIGNMENT"<<endl;
	this->Ngenes=Ch1.Ngenes;
	this->Chrom=Ch1.Chrom;
	this->fit=Ch1.fit;
	return *this;
}

Chromosome::Chromosome(Chromosome && Ch_old) : Ngenes(Ch_old.Ngenes), Chrom(move(Ch_old.Chrom)), fit(Ch_old.fit)
{
//invalidate Ch_old. DO i have to do it, or probably Chrom(std::move(Ch_old.Chrom)) will call the move constructor of vector class which does everything for me
	//cout<<"MOVE"<<endl;
}


Chromosome& Chromosome::operator=(Chromosome && Ch_old) 
{
	//cout<<"MOVE ASSIGNEMENT"<<endl;
	this->Ngenes=Ch_old.Ngenes; //it is copied, move is not made for double int... but for objects that deal with resources like memory...
	this->fit=Ch_old.fit;
	this->Chrom=Ch_old.Chrom; //i suppose it uses the move assignement of the vector's class
	return *this;
}

bool Chromosome::operator<(const Chromosome & Ch) const
{
	return this->get_fit()<Ch.get_fit(); 
}


vector<int> & Chromosome::get_chromosome_ref() 
{
	return Chrom;
}

vector<int> Chromosome::get_chromosome() const
{
	return Chrom;
}

int & Chromosome::get_gene_i_ref(int i)
{
return Chrom[i];
}

void Chromosome::print() const
{

	for(const auto & i : Chrom)
		cout<<i<<" ";
	cout<<endl;
	
}

void Chromosome::calculate_fit(const vector<vector<double>> & cities) //matrix of citites coordinates
{
	this->fit=0;
	Chrom.push_back(*(Chrom.begin())); //just temporaneily add home at end of path
	
	for(auto it=Chrom.begin(); it<Chrom.end()-1; it++) 
	{
		fit+=sqrt(transform_reduce(cities[*it].begin(), cities[*it].end(), cities[*(it+1)].begin(), 0.0 , plus<>(), [](double x, double y){ 
																																																											return pow(x-y,2);}));
	}
	Chrom.erase(Chrom.end()-1); //eliminate home end that was added
}

double Chromosome::get_fit() const
{
	return fit;
}

bool Chromosome::check() const //we need it to be efficent since is a method that will eb called continously
{
	unordered_set<int> s(Chrom.begin(), Chrom.end()); 
	if(s.size()!=Chrom.size())
		cout<<"error: breaking bound"<<endl;
		
	return (s.size()==Chrom.size()); 
}

bool Chromosome::check_2() const//we need it to be efficent since is a method that will eb called continously
{

	bool flag=true;
	
	for(auto it=Chrom.begin(); it<Chrom.end(); ++it) 
	{
		if(count(Chrom.begin(), Chrom.end(), *it)!=1) 
		{
			flag=false;
			cout<<"error: breaking bound"<<endl;
		}
	}
	
	return flag;
}

void Chromosome::print_config(string filename, const vector<vector<double>> & Coords)
{
	ofstream out;
	
	out.open(filename);	
	for(auto it=Chrom.begin(); it<Chrom.end(); it++)
	{
		for(unsigned int i=0; i<Coords[Chrom[0]].size(); i++)
		{
			out<<Coords[*it][i]<<"\t";
		}
	out<<endl;	
	}
	
	for(unsigned int i=0; i<Coords[Chrom[0]].size(); i++)
	{
		out<<Coords[*Chrom.begin()][i]<<"\t"; // reprint the starting point coords
	}

	out.close();
}

///////////////////mutations

void Chromosome::mutation1(Random & rnd)//i should be between 1 and Ngenes-1?
{
	int i=rnd.Rannyu(1, Ngenes-1); //gives [1, Ngenes-2]. Chrom[1].... Chrom[Ngenes-1] accesibale values
	std::swap(Chrom[i], Chrom[i+1]); //std::swap != std::vector::swap 
}

void Chromosome::mutation2(Random & rnd)
{
	int i=rnd.Rannyu(1, Ngenes-2);//[1, Ngens-3] 
	int m=rnd.Rannyu(2, Ngenes-i);//[2, Ngenes-i-1]
	int j=rnd.Rannyu(i+m+1, Ngenes+1);//[i+m+1, Ngenes] v.begin()+v.size()-1 -> last element iterator

  rotate(Chrom.begin() + i, Chrom.begin() + i + m, Chrom.begin() + j);//	
}


//m : length of block to be permuted, i starting point (include), j starting point (include) of other block
void Chromosome::mutation3(Random & rnd) 
{

	//to get rand unif int numbers in [int1,int2] (with int1<int2) pass to Rannyu() int1 and int2+1 (this refets to INT!!!! not double)
	int m=rnd.Rannyu(2, Ngenes/2.);//m<Ngenes/2 //ma i casi m=0 e m=1 li voglio?? forse no, m=0 non fa niente, m=1 permuta le città i j, forse è un po brutale, meglio togleire m=0 e m=1 -> [2, int(Ngenes/2)-1] (IMP 2. THE . IS FUNDAMNETAL TO MAKE TI WORK BOTH IN Negenes EVEN AND ODD)
	int i=rnd.Rannyu(1, Ngenes-2*m+1); //1<=i<=Ngenes-2m 
	int j=rnd.Rannyu(i+m, Ngenes-m+1);//i+m<=j<=N-m
	
	std::swap_ranges(Chrom.begin()+i, Chrom.begin()+i+m, Chrom.begin()+j);
}

void Chromosome::mutation4(Random & rnd) 
{
	int i=rnd.Rannyu()*Ngenes+1; //[1, Ngenes] 
	int j=rnd.Rannyu()*Ngenes+1;
	if(i<j)	
		reverse(Chrom.begin()	+	i, Chrom.begin() + j);
	else
		reverse(Chrom.begin()	+	j, Chrom.begin() + i);
}

/*
void Chromosome::mutation5(Random & rnd)
{
	//imagine Ngenes is Ngenes -1 since the first miust be fixed
	int m=rnd.Rannyu(1, Ngenes+1); //1<=m<=Ngenes//point of initial shift
	//int n=rnd.Rannyu(m, Ngenes+1); //n>=m and n<=Ngenes number of position to be shifted. !!!!!!!!!!!!!!non sto usando n
	//cout<<"m	"<<m<<endl;
	rotate(Chrom.begin()+1, Chrom.begin()+m, Chrom.end()); //non va bene ?
}
*/


/////////////crossover ----------------use the one in population, this isn't tested

void Chromosome::crossover(Chromosome & Chrom2, Random & rnd) 
{

	vector<int> & C1=this->get_chromosome_ref(); //non modificarli se non per copia solo alla fine, usali per lettura, modifica appo1 e appo2
	vector<int> & C2=Chrom2.get_chromosome_ref();
	int i=rnd.Rannyu(1, C1.size()); //index da 1 a Ngenes-1

	vector<int> appo1=C1; //ma una copia di una ref è una ref?no ho controllato nel main
	vector<int> appo2=C2;
	
	int count1=i;
	int count2=i;
	cout<<i<<endl;
	for(unsigned int k=1; k<C2.size(); k++)
	{	
		//if(any_of(C1.begin()+i, C1.end(), [&](int x){return x==C2[k];})) //any of was better i think
		if(C2[k]==*(find(C1.begin()+i, C1.end(), C2[k])))
		{
		appo1[count1]=C2[k];
		count1++;
		}
		
		//if(any_of(C2.begin()+i, C2.end(), [&](int x){return x==C1[k];}))
		if(C1[k]==*(find(C2.begin()+i, C2.end(), C1[k])))
		{
		appo2[count2]=C1[k];
		count2++;
		}
	}
	
	C1=move(appo1);//TO BE CHECKED
	C2=move(appo2);

}
