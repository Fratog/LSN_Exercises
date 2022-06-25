#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>									//mpicxx -std=c++17 SOURCE/main.cpp //mpiexec -np 4 a.out
#include "mpi.h"

using namespace std;

//https://en.wikipedia.org/wiki/Message_Passing_Interface#Communicator

const int n = 100;

int main(int argc, char *argv[])
{
//mio: cosa stiamo facendo con queste MPI
//1) la tua architteura consisite di un NUMA NODE, a 16 core (di cui me en fa utilizzare 8 ) quindi magari non so bene cosa intenda con 16 core... vomunque per vedere le informazioni usa lscpu 
//2) io penso che bisogni distingure tra architteura fisica (numa node-> distributed memory with non uniform memory acces and shared addres space ) e tra il modello di parallel computing che usiamo: es shared memory, message passing, thread... hybrid... //Magari non posso direlo con generalità assoluta, ma penso che anche su un computer come il mio con un numa node, puoi usare delle librerie per fare un parallel computing basato su shared memory.
//3) cosa sto facendo qua?, secondo me sto facendo parallel computing con message passing model (and so obviously i think separated memeory and addres space). so this is the framework in which i'll think. and SPMD (single program multple data) so imagini n cores identified by rank, which ran the same executable



//////////
/*
int size, rank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

vector<int> v(100);

for(int i=0; i<v.size(); i++) {v[i]=1;}
cout<<"rank	"<<rank<<endl;

for(auto i : v) {cout<<i<<"	";}
cout<<endl<<endl;

cout<<"rank	"<<rank<<endl;

int i=0;
double a=0;
double b=0;

if(rank==0)
	a=accumulate(v.begin(), v.begin() + int(v.size()/2.), 0.0);
if(rank==1)
	b=accumulate(v.begin()+int(v.size()/2.), v.end(), 0.0);

 

cout<<" Sono il nodo "<<rank<<" dei "<<size<<" che haiutilizzato!"<<endl;
MPI_Finalize();



cout<<"rank	"<<rank<<"	result	"<<"	a	"<<a<<" b	"<<b<<"	a+b	"<<a+b<<endl;
*/
/*
int size, rank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int my_values[3];
	for(int i=0;i<3;i++)
		if(rank==1) 
			my_values[i]=i+1;
		else 
			my_values[i]=0;
			
	cout<< "Prima: "<< my_values[0]<< " "<< my_values[1]<<" "<< my_values[2]<< " per il processo "<< rank<< endl;
	//http://mpi.deino.net/mpi_functions/MPI_Bcast.html : buffer, entry in buffer, type, rank of root, coomunicotr (mio:i.e family of processes involved)
	MPI_Bcast(my_values,3,MPI_INTEGER,1, MPI_COMM_WORLD); //root is the starting core/process from wich the message will be broadcasted
	cout<< "Dopo: "<< my_values[0]<< " "<< my_values[1]<<" "<< my_values[2]<< " per il processo "<< rank<< endl;


MPI_Finalize();
*/
/*
///////////////////////////////////SCATTER https://www.open-mpi.org/doc/v3.1/man3/MPI_Scatter.3.php
int size, rank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//vector<int> appo(size);
	int isend[size];
	int irecv=0; // 0 for every rank
	
	for(int i=0;i<size;i++) 
		if(rank==0)
			isend[i]=i; //what we want to scatter
		else
			isend[i]=0;
//int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,MPI_Comm comm)		
MPI_Scatter(&isend, 1, MPI_INTEGER, &irecv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
cout<<"rank "<<rank<<"	received	"<<irecv<<endl;


//rank 6	received	6
//rank 7	received	7
//rank 0	received	0
//rank 1	received	1
//rank 2	received	2
//rank 3	received	3
//rank 4	received	4
//rank 5	received	5


MPI_Finalize();
*/
////////////////////GATHER
/*
int size, rank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(size>3)
	{
		cout<<"Hai scelto troppi processi"<<endl;
		return 1;
	}

	int irecv[3];
	for(int i=0;i<3;i++) 
		irecv[i]=0; //empty array
	
	int isend = rank + 1; //every process will have a isend differnt, with the next  lines we want to gather them in the root process
	
//int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
	MPI_Gather(&isend,1,MPI_INTEGER, irecv,1,MPI_INTEGER,0, MPI_COMM_WORLD); //why irecv doesn't have &?
	if(rank==0) 
		cout<< "irecv: " <<irecv[0] <<" "<<irecv[1] <<" " <<irecv[2] <<endl;//irecv: 1 2 3
		
MPI_Finalize();
*/

//////////////////////reduce
/*
int size, rank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int isend[2], sum, prod;
	for(int i=0;i<2;i++) 
		isend[i]=rank+i+1;

//int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
	MPI_Reduce(&isend[0], &sum, 1, MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD); //sum the values of isend[0] of every process-> 1+2+3+..
	MPI_Reduce(&isend[1], &prod, 1, MPI_INTEGER, MPI_PROD,0,MPI_COMM_WORLD); //prod

	if(rank==0)//the values are sent to the root caller i guess
	cout<<"irecv: "<<sum<<" "<<prod<<endl;
	
MPI_Finalize();
*/

//////////////////////SPLIT //https://www.open-mpi.org/doc/v4.1/man3/MPI_Comm_split.3.php
/*

comm
    Communicator (handle). 
color
    Control of subset assignment (nonnegative integer). 
key
    Control of rank assignment (integer). ->>>>>>>>>>>>>>???
*/
/*
int size, rank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(size!=4)
	{
		cout<<"Servono 4 processi, non"<<size<<"!!"<<endl; 
		return 1;
		}

	int icolor, ikey;
	if(rank==0){icolor=1;ikey=2;}// i think icolor identify the one that will form a new substet, ikey will be the new rank in the given subset
	if(rank==1){icolor=1;ikey=1;}//so rank==1 becomes 0 di 1
	if(rank==2){icolor=2;ikey=1;}
	if(rank==3){icolor=2;ikey=2;}
	
	//int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
	MPI_Comm nuovocom;
	MPI_Comm_split(MPI_COMM_WORLD,icolor,ikey,&nuovocom);
	
	int newsize,newrank;
	MPI_Comm_size(nuovocom, &newsize);
	MPI_Comm_rank(nuovocom, &newrank);
	cout<<"Ero: "<<rank<<" di "<<size<<" ... e adesso sono: "<<newrank<<" di "<<newsize<<endl;
MPI_Finalize();
*/

/////////////////////////Blocking unidirectional send and receive //https://www.open-mpi.org/doc/v4.0/man3/MPI_Recv.3.php
/* send arguments
buf
    Initial address of send buffer (choice). 
count
    Number of elements send (nonnegative integer). 
datatype
    Datatype of each send buffer element (handle). 
dest
    Rank of destination (integer). 
tag
    Message tag (integer). 
comm
    Communicator (handle)

*/
/*
int size, rank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Status stat;
//MPI_Status stat2;//i was trying to send the same message imesg ti differrent cores, but maybe this can't be done here since this should be a unidirectional message? i don't know if what i am saying is tr

//https://www.rookiehpc.com/mpi/docs/mpi_status.php
//MPI_Status represents the status of a reception operation, returned by receive operations (MPI_Recv), non-blocking operation wait (MPI_Wait, MPI_Waitall, MPI_Waitany and MPI_Waitsome) or test (MPI_Test, MPI_Testall, MPI_Testany and MPI_Testsome). In C the MPI_Status is a structure that contains at least 3 attributes: MPI_SOURCE, MPI_TAG and MPI_ERROR.


	int itag=1;//non mi è chiaro a cosa serva itag, da una tag al messagio e questa tag viene automaticamente slavata in stat inseieme alla source -> a rank 1????????'
	int imesg = rank;
	int itag2=2;
	int imesg2= rank;
	if(rank==1)
		//int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
		MPI_Send(&imesg,1,MPI_INTEGER,0,itag,MPI_COMM_WORLD);
	else if(rank==0)
		//int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
		MPI_Recv(&imesg,1,MPI_INTEGER,1,itag,MPI_COMM_WORLD, &stat);
	else if(rank==2)
		MPI_Recv(&imesg2,1,MPI_INTEGER,1,itag2, MPI_COMM_WORLD, &stat2);
	if(rank!=1)
		cout<<"this rank:	"<<rank<<" received messaggio = "<<imesg<<endl;


MPI_Finalize();
*/

//another send receive to try underand it a bit better

//mpiexec -np 3 a.out               3 or greater, but doesn't work
/*
int size, rank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Status stat;

	int itag=1;
	int imesg = rank;
	
	if(rank==0)
	{	//int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
		for(int i=1; i<size; i++)
		{
			MPI_Send(&imesg,1,MPI_INTEGER, i ,itag, MPI_COMM_WORLD);
		}
		cout<<rank<<"	sending	"<<endl;
	}
	else if(rank!=0)
	{	
		MPI_Recv(&imesg,1,MPI_INTEGER, 0 ,itag ,MPI_COMM_WORLD , &stat);
		cout<<rank<<"	receiving"<<endl;
	}
	
	if(rank!=0)
		cout<<"this rank:	"<<rank<<" received messaggio = "<<imesg<<endl;


MPI_Finalize();
*/
//BIDIRECTIONAL MESSAGE, ALWAYS USING THE SAME SEND AND RECEIVE

int size, rank;
MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Status stat1, stat2;

	int* imesg = new int[n]; 
	int* imesg2 = new int[n];
	int itag=1; 
	int itag2=2;
	
	for(int i=0;i<n;i++)
	{
		imesg[i]=rank; 
		imesg2[i]=rank+1;
	}

	if(rank==1)
	{
		cout<<"rank	"<<rank<<"	was	"<<imesg2[0]<<endl;//expect 2
		MPI_Send(&imesg[0],n, MPI_INTEGER,0,itag,MPI_COMM_WORLD);
		MPI_Recv(&imesg2[0],n,MPI_INTEGER,0,itag2, MPI_COMM_WORLD, &stat2);
		cout<<"rank	"<<rank<<" received messaggio = "<<imesg2[67]<<endl;// mu aspetto sia 1 , quando prima di ricevere era 2 suppongo
	}
	else if(rank==0)
	{
		cout<<"rank	"<<rank<<"	was	"<<imesg[0]<<endl; //expext 0
		MPI_Send(&imesg2[0],n, MPI_INTEGER,1,itag2, MPI_COMM_WORLD);
		MPI_Recv(&imesg[0],n,MPI_INTEGER,1, itag, MPI_COMM_WORLD,&stat1);
		cout<<"rank	"<<rank<<" received messaggio = "<<imesg[45]<<endl; //mi aspetto si 1, qunado prima di ricevere era 0
	}
	
MPI_Finalize();

/*
rank	1was	2
rank	1 received messaggio = 1
rank	0was	0
rank	0 received messaggio = 1
*/
return 0;
}
