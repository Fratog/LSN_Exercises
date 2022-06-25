#include <iostream>
#include <algorithm>
#include <unordered_set>
#include "population.h"

using namespace std;

Population::Population(int N, int start)
{
	this->start=start;
	this->Nchrom=N;
	fits.resize(Nchrom);
	Pop.resize(Nchrom);
} 

Population::Population(int N, int start, const vector<vector<double>> & Coords, int primes_row) 
{

	if(primes_row>31500)
	{
		cout<<"primes row too large"<<endl;
		rnd.starting(1);
	}		
	else	
	{	
		rnd.starting(primes_row);
	}
	
	this->start=start;//non credo gli serva salvarselo
	this->Nchrom=N;
	Pop.resize(Nchrom);
	
	for(auto it=Pop.begin(); it<Pop.end(); it++) 
	{
		(*it)=(Chromosome(Coords.size(), start, this->rnd));
		(*it).calculate_fit(Coords);
		fits.push_back((*it).get_fit());
		//cout<<" read fit	"<<(*it).get_fit()<<endl;;
	}
	
	sort(fits.begin(), fits.end());
	this->check_population();
	this->sort_pop(); //on a fitness basis thanks to < overload in chromosome
}

Population::~Population()
{
}

vector<Chromosome> & Population::get_pop() 
{
	return Pop;
}

vector<Chromosome> Population::get_pop_fail() const
{
	return Pop; 
}

Chromosome & Population::get_chrom_i_ref(int i) 
{
	return Pop[i]; 
}

Chromosome Population::get_chrom_i(int i) const
{
	return *(Pop.begin()+i); //rimetti Pop[i] dopo
}

Random & Population::get_gen()
{
	return this->rnd;
}

bool Population::check_population() const
{
	bool flag=true;
	
	for(const auto & i : Pop)
		flag=i.check();

	return flag;			
}

void Population::print_fits() const
{
	for (const auto & i : fits)
		cout<<i<<" ";
	cout<<endl;		
}

void Population::set_fit(int i)
{
	fits[i]=get_chrom_i(i).get_fit();
	//fits.push_back(get_chrom_i(i).get_fit());
}

vector<double> Population::get_fits() const
{
	return fits;
}
  

void Population::sort_pop()
{
	sort(Pop.begin(), Pop.end());
}

void Population::sort_fits()
{
	sort(fits.begin(), fits.end());
}

Chromosome Population::selection1()	//extracts chromosomes, from a pdf based on their fitness (or better 1-fitness)? 
{
	
	double T=50.;
	static double beta =1/T;
	int i=rnd.Rannyu(0, Nchrom/2.);//truncated population
	
	static int index=Nchrom-1;//starting point of metropolis SA
	static int old=fits[Nchrom-1];//starting value of metroplis SA
	static int accepted=0;
	static int rejected=0;
	static int counter=0;
	
	++counter;
	
	if(this->rnd.Rannyu()<=exp(-beta*(fits[i]-old)))
	{
		index=i;
		old=fits[i];
		accepted++;
	}	
	else
	{
	//old remains the same
		rejected++;
	}		
/*
	if(counter%Nchrom==0)
	{

	cout<<"----------------"<<endl;
	cout<<"actaual best	"<<fits[0]<<endl;
	cout<<"index	"<<index<<endl;
	cout<<"beta	"<<beta<<endl;
	cout<<"rate "<<(double)(accepted)/(accepted+rejected)<<endl;
	cout<<"------------------"<<endl;
	//cout<<"beta	"<<beta<<endl;
	//cout<<"rate "<<(double)(accepted)/(accepted+rejected)<<endl;;
	
	}
	*/
	if(counter%(Nchrom)==0)//update temperaature, on time each new generaation.
	{
		beta/=0.99;
	}
	return Pop[index];
	
}

Chromosome Population::selection2(double p) //useses the ordering of the fintess. I guess with varying p this can be both democratic or elitist
{
	return Pop[int(Nchrom*pow(rnd.Rannyu(), p))]; 
}

void Population::crossover(Chromosome & Chrom1, Chromosome & Chrom2) 
{
	vector<int> & C1=Chrom1.get_chromosome_ref(); //non modificarli se non per copia solo alla fine, usali per lettura, modifica appo1 e appo2
	vector<int> & C2=Chrom2.get_chromosome_ref();
	int i=rnd.Rannyu(1, C1.size()); //index da 1 a Ngenes-1

	vector<int> appo1=C1; 
	vector<int> appo2=C2;
	
	int count1=i;
	int count2=i;

	for(unsigned int k=1; k<C2.size(); k++)
	{	

		if(any_of(C1.begin()+i, C1.end(), [=](int x){return x==C2[k];}))
		{
		appo1[count1]=C2[k];
		count1++;
		}

		if(any_of(C2.begin()+i, C2.end(), [=](int x){return x==C1[k];}))
		{
		appo2[count2]=C1[k];
		count2++;
		}
	}

	C1=move(appo1);
	C2=move(appo2);
}

void Population::new_generation(const vector<vector<double>> & Coords, int flag)
{
	
	ofstream out_fit_av;
	
	double p1, p2, p3, p4, p_select, p_cross;
	int count=0;
	static int pippo=0;
	Chromosome appo1, appo2;
	Population new_pop(this->Nchrom, this->start);
	
	p1=0.08;
	p2=0.08;
	p3=0.08;
	p4=0.08;
	p_select=2;
	p_cross=0.65;
	
	while(count!=Nchrom)
	{
	//slect 1 2
	//crossover
	//mutatoin
	//add to new gen
		if(flag==1)
		{
			appo1=this->selection1(); //there should be RV0 optimization here. or a move. try to add std::move to see if it gives warning	
			appo2=this->selection1();
		}
		else if(flag==2)	
		{	
			appo1=this->selection2(p_select); //there should be RV0 optimization here. or a move. try to add std::move to see if it gives warning	
			appo2=this->selection2(p_select);
		}
		else
		{	
			cout<<"error in choosing selection operator"<<endl;	
		}	
	
		if(this->rnd.Rannyu()<p_cross)
			this->crossover(appo1 ,appo2);  
		
		if(rnd.Rannyu()<p1)
			appo1.mutation1(this->get_gen());
		if(rnd.Rannyu()<p2)
			appo1.mutation2(this->get_gen());
		if(rnd.Rannyu()<p3)
			appo1.mutation3(this->get_gen());
		if(rnd.Rannyu()<p4)
			appo1.mutation4(this->get_gen());
	
		
		if(rnd.Rannyu()<p1)
			appo2.mutation1(this->get_gen());
		if(rnd.Rannyu()<p2)
			appo2.mutation2(this->get_gen());
		if(rnd.Rannyu()<p3)
			appo2.mutation3(this->get_gen());
		if(rnd.Rannyu()<p4)
			appo2.mutation4(this->get_gen());
		
		appo1.check();
		appo2.check();
		
		appo1.calculate_fit(Coords);//they should have new values for the fit as mutations changed these
		appo2.calculate_fit(Coords);

		new_pop.get_chrom_i_ref(count)=move(appo1); //maybe force a move with std::move	
		new_pop.set_fit(count); //it isn't
		
		count++;
	
		if(count==Nchrom) //if Nchrom odd we risk adding one more than necessary
			break;

		new_pop.get_chrom_i_ref(count)=move(appo2); //maybe force a move with std::move
		new_pop.set_fit(count);
		
		count++;
	}	

	//this->rnd=new_pop.get_gen(); //i don't think i want to do this, it is okay to keep always the same generator
	this->fits=new_pop.get_fits();
	this->Pop=new_pop.get_pop(); //updating the popolation to the one just created
	
	
	//this->check_population(); //probably redundant since all appo are already checked
	this->sort_fits();
	this->sort_pop();

	//cout<<"new generation best fit: "<<this->fits[0]<<endl;//minimum fit of the new (now current) population.
	

	/*
	out_fit_av.open("fit_square.txt", ios::app);

	if(pippo==0)
		out_fit_av<<"best fit of each generation"<<"		"<<"avarage value of the fits on best half of the poulation"<<"	"<<"errors"<<endl;
	pippo=1;
	
	double av=0;
	double variance=0;
	av=reduce(fits.begin(), fits.begin()+fits.size()/2., 0.0)/(fits.size()/2.);
	
	for_each(fits.begin(), fits.begin()+fits.size()/2., [&] (const double d) {   //& lmambd,body gets acces to values out of his scope, by reference
		  variance += (d - av) * (d - av);
			});
			
	variance=variance/(fits.size()/2.);
	
	//out_fit_av<<fits[0]<<"	"<<reduce(fits.begin(), fits.begin()+fits.size()/2.)/(fits.size()/2.)<<"	"<<sqrt(variance/fits.size()/2.)<<endl;	
	//out_fit_av<<fits[0]<<"	"<<av<<"	"<<sqrt(variance/fits.size()/2.)<<endl;	
	out_fit_av<<fits[0]<<"	"<<av<<endl;	
	out_fit_av.close();	
	*/
	
}

