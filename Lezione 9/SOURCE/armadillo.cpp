#include<iostream>
#include<armadillo>

using namespace std;

class appo
{
public:
void print() {cout<<"hello"<<endl;}
private:
int a,b;
};

int main()
{
	arma::mat A(2,3, arma::fill::zeros);
    
  // .n_rows and .n_cols are read only
  cout << "A.n_rows: " << A.n_rows << endl;
  cout << "A.n_cols: " << A.n_cols << endl;
     
  A(1,2) = 456.0;  // access an element (indexing starts at 0)
  A.print("A:");

appo pippo;
/////////////////7
arma::field<appo> kappa(10);

//for (auto i : kappa) {cout<<(i.print())<<endl;}

for(int i=0; i<10; i++)
kappa[i].print();
return 0;
}
