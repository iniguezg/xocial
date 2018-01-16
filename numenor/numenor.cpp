/* Numenor v1.0
 * A simple model to describe opinion formation on coevolving social networks
 * Copyright (C) 2018 Gerardo Iñiguez
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * If you use this program for your research, please cite any relevant from the following articles:
 *
 * G. Iñiguez, J. Kertész, K. Kaski, R. A. Barrio
 * Opinion and community formation in coevolving networks
 * Physical Review E 80, 066119 (2009)
 * DOI: 10.1103/PhysRevE.80.066119
 * arXiv: 0908.1068v2
 *
 * G. Iñiguez, R. A. Barrio, J. Kertész, K. Kaski
 * Modelling opinion formation driven communities in social networks
 * Computer Physics Communications 182, 1866–1869 (2011)
 * DOI: 10.1016/j.cpc.2010.11.020
 * arXiv: 1007.4177
 *
 * G. Iñiguez, J. Kertész, K. Kaski, R. A. Barrio
 * Phase change in an opinion-dynamics model with separation of time scales
 * Physical Review E 83, 016111 (2011)
 * DOI: 10.1103/PhysRevE.83.016111
 * arXiv: 1009.2643
 *
 * G. Iñiguez, J. Tagüeña-Martínez, K. Kaski, R. A. Barrio
 * Are opinions based on science: Modelling social response to scientific facts
 * PLoS ONE 7, e42122 (2012)
 * DOI: 10.1371/journal.pone.0042122
 * arXiv: 1109.1488
 */


#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include "./randomc.h" //include for random library

//number of agents in the network
int N;

//number of entries in connectivity matrix, R = N^2
int R;

//some variables for node rewiring and y-rewiring
int inc; //number of agents willing to be cut from agent i
int inw; //number of agents willing to be rewired to agent i
int inwY; //number of agents willing to be y-rewired to agent i

//external field in the equations of motion
double h;

//internal parameter used to control the relative amounts of rewiring and y-rewiring
double y;

//rewiring window
double win;

//bifurcation parameter for the terms in the equations of motion
int mcut;

//generation time
int tgen;

//mean connectivity for the initial random network
int z;

//structure for the output variables
struct out
{
	double nb, ng, nw;			//number of black, gray and white agents
	double nclusters;			//number of clusters in network
	double cmax;				//maximum size of cluster in network
	double avclu;				//average cluster size in network
	double suscep;				//susceptibility of the network
	double meanz;				//mean connectivity per site
	double ccoef;				//clustering coefficient of the network
	double plen;				//characteristic path length in the network
	double neigh;				//mean number of nth-order neighbours in the network
	double Fplus, Fminus, fplus, fminus;	//F,f rewiring variables
	double rewnodes;			//number of rewiring nodes

};

//maximum diameter for an arbitrary network
int dmax;

//absolut limit and cut-off parameters used to rewire agents
double lim;
double cutoff;

//number of rewirings per generation
int rewcheck;

//number of times each dynamics is run
int ntimes;
double dntimes; //the double version

//pointer to the x array, which stand for the agents' dynamical variables in the equations of motion
double * arrx;

//pointer to the buffer of the previous array
double * buffx;

//pointer to the array of parameters for the equations of motion
double * alpha;

//pointer to the connectivity array, and its buffer
int * arrA;
int * buffA;

//pointers to the rewire arrays
int * arrP, * arrQ, * arrR;

//pointer to the shortest path array
int * arrsp;

//pointer to the x check array, used to check if agent i is still alive
int * xcheck;

//pointer to the array of c which holds the average number of clusters of a given size, and the counter for number of elements in the average
double * valc;
int ccount;

//pointers to the arrays of z which hold, for each i, <znn> and the number of elements in that average
double * valznn;
int * znncount;

//structs for the average time series, and variables for the size and number of time steps, and the average counter
out Tone, Tmean;
out * Tacc;
int twrite, twmax, wcount;

//pointers to the arrays where the cut and rewired agents of agent i are stored
int * incut;
int * inwire;
int * inwireY;

//pointers for the arrays of exposed and live vertices, counted respectively with Xi and Yi
int * Li;
int * Oi;

//pointer to the c check array, used to see if agent i has already been placed in a cluster
int * ccheck;

//pointer to the array of cluster sizes
int * clust;

//file descriptor for the clusters' distribution
FILE * cdib;
char cdibf2w[255];

//file descriptor for the z correlation output
FILE * zcor;
char zcorf2w[255];

//file descriptor for the time series output
FILE * tser;
char tserf2w[255];

//file descriptors for the himmeli files
FILE * config;        //himmeli configuration file
FILE * edges;         //himmeli edges file
FILE * vertices;      //himmeli vertices file
char configf2w[255];
char edgesf2w[255];
char verticesf2w[255];

//this will come in handy for the function declarations
typedef double * pdouble;
typedef int * pint;

//random number generator library
TRanrotBGenerator * randing;

//random number generator initializer
void RandInit()
{
	long timech;
	timech = time((time_t)NULL);
	randing = new TRanrotBGenerator(timech);
	randing->RandomInit(timech);
}

//we characterise the place of a pair of agents in A by linear (i=1,..,R) and matrix (r;c=1,..,N;1,..,N) position
//the next two little functions give us the direct and inverse transformations between these two positions

int Amat2lin(int ra, int ca)
{
	int pos;
	pos = (ra - 1)*N + ca;
	return pos;
}

void Alin2mat(int ia, int& ra, int& ca)
{
	double dia = (double) ia;
	double dN = (double) N;
	double iaN = (double) ( dia/dN );
	if (ia % N)
	{
		ra = (int) ( floor(iaN) + 1.0 );
		ca = ia % N;
	}
	else
	{
		ra = (int) iaN;
		ca = N;
	}
}

//these tiny functions compute the probabilities pij, qij, y-qij associated with agents i, j

double pcut(int ia, int ja)
{
	double xi; //x for agents i,j
	double xj;
	double pija; //probability of breaking the bond between them

	xi = arrx[ia - 1];
	xj = arrx[ja - 1];

	pija = fabs( xi - xj )/(2.0*lim);

	return pija;
}

double qwire(int ia, int ja)
{
	double xi; //x for agents i,j
	double xj;
	double qija; //probability of creating a bond between them

	xi = arrx[ia - 1];
	xj = arrx[ja - 1];

	qija = fabs( xi + xj )/(2.0*lim);

	return qija;
}

double qwireY(int ia, int ja)
{
	double xi; //x for agents i,j
	double xj;
	double qija; //probability of creating a bond between them

	xi = arrx[ia - 1];
	xj = arrx[ja - 1];

	qija = 1.0 - ( fabs( xi - xj )/(2.0*lim) );

	return qija;
}

//the next two small functions calculate the initial conditions for the x array, taking samples either from a uniform
//distribution in xlimit*[-1,1], or from a normal distribution specified by mean and var

void uniXarray(double xlimit)
{
	double rand; //random number in [0,1] interval

	for (int i = 1; i <= N; i++)
	{
		rand = randing->Random();
		arrx[i - 1] = xlimit*( 2.0*rand - 1.0 ); //we normalize to xlimit*[-1,1] interval

		xcheck[i - 1] = 1; //all agents are initially alive
	}
}

void normXarray(double mean, double var)
{
	double x1, x2; //uniform deviates in [0,1] interval
	double v1, v2; //x transformation to [-1,1] interval
	double Rsquare; //R^2 = v1^2 + v2^2
	double fac; //Box-Muller transformation
	double n1, n2; //normal deviates, chosen in lim*[-1,1] interval

	int halfN, M; //integer half of N, computed with the help of M
	int limitpass; //used to check if n1,n2 are inside the lim*[-1,1] interval
	int circlepass; //used to check if (v1,v2) it's inside the unit circle

	if ( (N % 2) == 0 ) M = N - 1; else M = N;
	halfN = (int) ceil( ( (double)M )/2.0 ); //we set the integer half of N

	for (int i = 1; i <= halfN; i++)
	{
		limitpass = 1; //we initialize the limit checker
		while (limitpass)
		{
			circlepass = 1; //we initialize the circle checker
			while (circlepass)
			{
				x1 = randing->Random(); x2 = randing->Random(); //find the x's
				v1 = 2.0*x1 - 1.0; v2 = 2.0*x2 - 1.0; //get the v's
				Rsquare = pow(v1, 2.0) + pow(v2, 2.0); //get R^2

				if ( (Rsquare < 1.0) && (Rsquare != 0.0) ) //if we're inside the unit circle
					circlepass = 0; //get out!
			}

			fac = sqrt( -2.0*log(Rsquare)/Rsquare ); //get fac

			n1 = sqrt(var)*(v1*fac) + mean; //we set the two normal deviates, scaling with mean and var
			n2 = sqrt(var)*(v2*fac) + mean;

			if ( (fabs(n1) < lim) && (fabs(n2) < lim) ) //if the normal deviates are in the required interval
				limitpass = 0; //get out!
		}

		//finally, set the two normal deviates in the x array!
		arrx[i - 1] = n1;	arrx[N - i] = n2;

		xcheck[i - 1] = 1;	xcheck[N - i] = 1; //all agents are initially alive
	}
}

//this little function computes the alpha constants for the equations of motion
void theAlphas()
{
	double rand; //random number in [0,1] interval

	for (int i = 1; i <= N; i++)
	{
		rand = randing->Random();
		alpha[i - 1] = 2.0*rand - 1.0; //we normalize to [-1,1] interval
	}
}

//this mighty little function computes the alpha constants in the interval [amin,amax]
void mightyAlphas(double amin, double amax)
{
	double rand; //random number in [0,1] interval

	for (int i = 1; i <= N; i++)
	{
		rand = randing->Random();
		alpha[i - 1] = (amax - amin)*rand + amin; //we normalize to [-1,1] interval
	}
}

//this little function calculates the number of non-zero elements in the array arr. It is used to compute Xi, Yi from Li, Oi.
int count(pint arr)
{
	int acc = 0; //counter
	for (int i = 1; i <= N; i++)
	{
		if ( arr[i - 1] != 0 ) //sum over the non-zero elements
			acc++;
	}
	return acc;
}

//the next two little functions sum the elements (and squares) of array arr. It is used with the array clust for suscep
int sum(pint arr)
{
	int el; //array element
	int acc = 0; //counter
	for (int i = 1; i <= N; i++)
	{
		el = arr[i - 1];
		acc += el;
	}
	return acc;
}

double sum2(pint arr)
{
	int el; //array element
	double acc = 0.0; //counter
	for (int i = 1; i <= N; i++)
	{
		el = arr[i - 1];
		acc += pow( (double)el, 2.0 );
	}
	return acc;
}

//this little function calculates the standard deviation of the system's dynamical variables
double squares()
{
	double mean = 0.0; //std dev counter, initialized to 0
	double xdiff; //difference in x per site per time step

	for (int i = 1; i <= N; i++)
	{
		xdiff = buffx[i - 1] - arrx[i - 1]; //we set xdiff
		mean += fabs(xdiff); //we add the squares
	}

	return mean;
}

//this baby function allocates memory space for n double arrays of size len pointed by arr
void Dalloc(pdouble& arr, int len, int n)
{
	int nelm; //number of elements in the array

	nelm = n*len; //len elements for each one of the n arrays in arr

	arr = (double *) calloc(nelm, sizeof(double)); //we allocate memory

	if (arr == NULL) //if something goes wrong...
	{
		printf ("error allocating memory\n");
		//close all open files and free any allocated memory
		exit (1);
	}
}

//this baby function allocates memory space for n int arrays of size len pointed by arr
void Ialloc(pint& arr, int len, int n)
{
	int nelm; //number of elements in the array

	nelm = n*len; //len elements for each one of the n arrays in arr

	arr = (int *) calloc(nelm, sizeof(int)); //we allocate memory

	if (arr == NULL) //if something goes wrong...
	{
		printf ("error allocating memory\n");
		//close all open files and free any allocated memory
		exit (1);
	}
}

//here come the prototypes of the functions called along the program

//this function calculates m powers of the int array arr after each rewiring process
void IpowersA(pint& arr, int m);

//this function constructs the connectivity matrix for a random network with mean connectivity z
void ranAmatrix(int z);

//this function finds the shortest path between all (i,j) pairs, and fills the sp array with them
void spath();

//this function computes the f1,f2 terms for agent i
void f1f2(int ia, double& f1, double& f2);

//this function computes the x dynamical term for agent i
double dynam(int ia);

//this function calculates nb, ng and nw
void ncount(out& one);

//this function computes the characteristic path length and the mean number of nth-order neighbours for the network
void length(out& one, int n);

//this function fills arrays incut and inwire with respect to agent i
void incount(int ia);

//this function fills arrays incut and inwireY with respect to agent i
void incountY(int ia);

//this function performs the network rewiring due to agent i over the A buffer
void rewire(int ia);

//this function performs the network y-rewiring due to agent i over the A buffer
void rewireY(int ia);

//this function computes the F,f average variables after each rewiring
void theFvars(out& one);

//this function percolates the network over a certain threshold thres dependent on the dynamical variables
void perc(double thres);

//this function searches the network for clusters and its sizes
void cluster();

//this function computes the c,z variables in the out struct, namely nclusters, cmax, avclu, suscep, meanz and ccoef, and all the correlations
void CZvars(out& one);

//this function constructs the c distribution graph and writes it down in cdib according to
//stage = 1,2,3, which corresponds respectively to setup, writing of data and construction of graph
void cgraph(int stage);

//this function constructs the z correlation graph that determines the assortativity of the network and writes it down in zcor according to
//stage = 1,2,3, which corresponds respectively to setup, writing of data and construction of graph
void zcorr(int stage);

//this function computes the average time series according to
//stage = 1,2,3, which corresponds respectively to setup, writing of data and construction of graph
void tseries(int stage);

//for himmeli and beyond!
void himmeli(char * str);

//this function initializes n double arrays arr with name str and size len
void DarrayInit(pdouble& arr, int len, int n, char * str);

//this function initializes n int arrays arr with name str and size len
void IarrayInit(pint& arr, int len, int n, char * str);

//this function does the memory mmaping of n double arrays arr with name str, size len and file descriptor filedesc.
void Dmmaparray(pdouble& arr, int len, int n, char * str, int& filedesc);

//this function does the memory mmaping of n int arrays arr with name str, size len and file descriptor filedesc.
void Immaparray(pint& arr, int len, int n, char * str, int& filedesc);


//Behold the main function
int main (int argc, char ** argv)
{
	//time variable
	int t;

	//time step
	double dt;

	//starting time
	int tstart;

	//...the beginning and the end...
	bool life;

	//parameters for the initial distributions in the x array
	double xlimit; //limit for the uniform distribution
	double mean; //mean and variance of the normal distribution
	double var;

	//variables needed to specify the parameter loop
	int initp; //initial value
	int finp; //final value
	int dp; //step

	//counter for the averages loop
	int nt;

	//standard deviation of the set of dynamical variables, and its buffer
	double stdev;
	double stdevbuff;

	//precision variables, used to compare stdev and stdevbuff, and limit rewcheck
	double zero;
	double rewzero;

	//x variable for agent i
	double xi;

	//dynamic term for agent i
	double dynxi;

	//structs for the output variables, their accumulators and average values
	out one;
	out acc;
	out Mean;

	//random number in [0,1] interval
	double rand;

	//percolation threshold
	double thres;

	//all the switches used in the process
	int mmapswitch;  //the basic arrays (x,alpha,A,sp) are mmaped (1) or allocated (0)
	int distswitch;  //the c,z distros are to be (1) or not to be (0) calculated
	int timeswitch;  //the time series are to be (1) or not to be (0) calculated
	int initswitch;  //the initial conditions for the basic arrays are read from files (1) or set manually (0)
	int lifeswitch;  //life is allowed to endure (1) or not (0)
	int percswitch;  //the network is percolated (1) or not (0)
	int himmswitch;  //the himmeli graph is drawn (1) or not (0)
	int finalswitch; //the end is set by rewiring (1) or by dynamics (0)

	//mean value of the terms added in the euler method sum, its buffer and counter
	double eul;
	double eulbuff;
	int eulcount;

	//some variables used in the dynamics
	int arrcheck; //used to check if the arrays have to be recalculated

	//variables used to mmap the basic arrays
	char xstr[] = "arrx";
	int xfiledesc;
	char alstr[] = "alpha";
	int alfiledesc;
	char Astr[] = "arrA";
	int Afiledesc;
	char spstr[] = "arrsp";
	int spfiledesc;

	//file descriptor for the output array
	FILE * gnu;
	char gnuf2w[255];

	//himmeli chars
	char inistr[] = "ini";
	char intstr[] = "inter";
	char finstr[] = "fin";

	//reading the command line arguments
	if(argc == 1)
	{
		N = 100;
		initp = 1;
		finp = 10;
		dp = 1;
	}
	else if(argc == 2)
	{
		N = atoi(&argv[1][0]);
		initp = 1;
		finp = 10;
		dp = 1;
	}
	else if(argc == 3)
	{
		N = atoi(&argv[1][0]);
		initp = atoi(&argv[2][0]);
		finp = 10;
		dp = 1;
	}
	else if (argc == 4)
	{
		N = atoi(&argv[1][0]);
		initp = atoi(&argv[2][0]);
		finp = atoi(&argv[3][0]);
		dp = 1;
	}
	else
	{
		N = atoi(&argv[1][0]);
		initp = atoi(&argv[2][0]);
		finp = atoi(&argv[3][0]);
		dp = atoi(&argv[4][0]);
	}

	//topological parameters
	R = (int) pow( (double)N, 2.0 );
	dmax = N;

	//dynamic parameters
	mcut = 1;
	h = 0.0;

	//rewiring
	lim = 1.0;
	win = 1.0;
	y = 0.0;

	//initial conditions
	z = 4;
	xlimit = 1.0*lim;
	mean = 0.0;
	var = 1.0;

	//times
	tstart = 0;
	dt = 0.0001;
//	tgen = 100;
	ntimes = 1;
	dntimes = (double) ntimes;

	//additional parameters
	zero = 0.001*dt;
	rewzero = 0.0*(double)N;
	cutoff = 0.0;
	twmax = 12000;
	thres = 0.5*lim;

	//switchs
	mmapswitch = 0;
	distswitch = 0;
	timeswitch = 0;
	initswitch = 0;
	lifeswitch = 1;
	percswitch = 0;
	himmswitch = 1;
	finalswitch = 0;

	//we initialize the random library
	RandInit();

	//we open the file descriptor
	sprintf(gnuf2w, "./numenor-N%i-initp%i-finp%i-dp%i.txt", N, initp, finp, dp);
	if((gnu = fopen(gnuf2w, "w+")) == (FILE *) NULL )
	{
		perror("\ncant open file to write\n");
		exit(1);
	}

	//we still are the programmers who say kni!

	//kni!

	//memory mmaping of loads of arrays

	//we get memory for the basic arrays according to their switch
	if ( mmapswitch )
	{
		Dmmaparray(arrx, N, 1, xstr, xfiledesc);
		Dmmaparray(alpha, N, 1, alstr, alfiledesc);

		Immaparray(arrA, R, 3, Astr, Afiledesc);
		Immaparray(arrsp, R, 1, spstr, spfiledesc);
	}
	else
	{
		Dalloc(arrx, N, 1);
		Dalloc(alpha, N, 1);

		Ialloc(arrA, R, 3);
		Ialloc(arrsp, R, 1);
	}

	//and we allocate the rest of them
	Dalloc(buffx, N, 1);

	Ialloc(buffA, R, 1);
	Ialloc(arrP, R, 1);
	Ialloc(arrQ, R, 1);
	Ialloc(arrR, R, 1);
	Ialloc(xcheck, N, 1);
	Ialloc(incut, N, 1);
	Ialloc(inwire, N, 1);
	Ialloc(inwireY, N, 1);
	Ialloc(Li, N, 1);
	Ialloc(Oi, N, 1);
	Ialloc(ccheck, N, 1);
	Ialloc(clust, N, 1);

	if ( distswitch ) //the distros' according to their switch
	{
		Dalloc(valc, N, 1);
		Dalloc(valznn, N, 1);

		Ialloc(znncount, N, 1);
	}

	//here we come! the dynamics are almost here

	for (tgen = initp; tgen <= finp; tgen += dp) //parameter loop, this time over the generation time tgen
	{
		//we initialize the accumulators
		acc.nb = 0.0; acc.ng = 0.0; acc.nw = 0.0;
		acc.nclusters = 0.0;
		acc.cmax = 0.0;
		acc.avclu = 0.0;
		acc.suscep = 0.0;
		acc.meanz = 0.0;
		acc.ccoef = 0.0;
		acc.plen = 0.0;
		acc.neigh = 0.0;

		//we write the distros according to their stage and switch
		if ( distswitch )
		{
			cgraph(1);
			zcorr(1);
		}

		//we write the time series according to its stage and switch
		if ( timeswitch )
			tseries(1);

		nt = ntimes; //we set the counter
		while ( nt ) //nt loop, used to obtain average measurements
		{
			//initial conditions for the basic arrays

			//we read them or set them manually according to their flags
			if ( (mmapswitch == 1) && (initswitch == 1) )
			{
				DarrayInit(arrx, N, 1, xstr);
				/*DarrayInit(alpha, N, 1, alstr);*/	mightyAlphas(-2.5, -2.0); //mighty alphas!
				IarrayInit(arrA, R, 3, Astr);
				IarrayInit(arrsp, R, 1, spstr);

				for (int i = 1; i <= N; i++)
					xcheck[i - 1] = 1; //all agents are initially alive
			}
			else if ( initswitch == 0 )
			{
//				uniXarray(xlimit); //uniformly distributed random initial conditions for the x array
				normXarray(mean, var); //normally distributed random initial conditions for the x array

				theAlphas(); //random alpha constants for the equations of motion
				ranAmatrix(z); //connectivity matrix of a random network

				IpowersA(arrA, 3); //we get A^2 and A^3
				spath(); //we ready the sp array
			}

			//we fill the initial A buffer and initialize the rewire arrays
			for (int i = 1; i <= R; i++)
			{
				buffA[i - 1] = arrA[i - 1];

				arrP[i - 1] = 0;	arrQ[i - 1] = 0;
			}

			//we begin time, life... and other stuff as well
			t = tstart;

			if ( lifeswitch )	life = true;
			else			life = false;

			twrite = tgen;
			wcount = twmax + 1;
			stdev = 1.0;
			arrcheck = 0;
			eulbuff = 0.0;	eulcount = 0;

			//we write the time series according to its stage and switch
			if ( timeswitch )
				tseries(2);

			//we construct the himmeli graph of the initial state of the network
			if ( himmswitch )
				himmeli(inistr);


			//Behold thou fools! For the dynamics art before thee!!

			while( life ) //while the agents' dynamics are still working
			{
				t++; //we increase the time counter

//				if ( (t % (10*tgen)) == 0 ) //to know where we stand...
//				{
//					ncount(one); //we get nb, ng, nw
//					std::cout << " at gen = " << t/tgen << " < nb ng nw > = < " << one.nb << "  " << one.ng << "  " << one.nw << " > " << std::endl;
//				}

				//the World is in motion
				for (int i = 1; i <= N; i++)
				{
					if ( xcheck[i - 1] ) //for the live agents
					{
						xi = arrx[i - 1]; //x variable for agent i

						dynxi = dynam(i); //dynamical term for the wave equation

						xi = xi + dt*dynxi; //equation of motion

						eulbuff += fabs(dt*dynxi); //we add terms in the euler check
						eulcount++;

						buffx[i - 1] = xi;  //we fill the buffer
					}
				}

				//array update
				for (int i = 1; i <= N; i++ )
				{
					if ( xcheck[i - 1] ) //we update the x of live agents
						arrx[i - 1] = buffx[i - 1];
				}

				//we fix the limit opinions
				for (int i = 1; i <= N; i++)
				{
					if ( xcheck[i - 1] ) //for the live agents, since we just have to fix each opinion once
					{
						xi = arrx[i - 1]; //we get the x

						if ( fabs(xi) > lim ) //if agent has reached its maximum opinion
						{
							//we fix it
							if ( xi > lim )
								arrx[i - 1] = lim;
							else // ( xi < -lim )
								arrx[i - 1] = -lim;

							xcheck[i - 1] = 0; //the agent's dead
						}
					}
				}

				//generation rewiring
				if ( (t % tgen) == 0 )
				{
					if (arrcheck) //we reset the rewire arrays from previous rewirings
					{
						for (int i = 1; i <= R; i++)
						{
							arrP[i - 1] = 0;	arrQ[i - 1] = 0;
						}
					}

					rewcheck = 0; 	//we initialize all counters
					arrcheck = 0;

					for (int i = 1; i <= N; i++) //we rewire the network
					{
						xi = arrx[i - 1]; //we get the x

						if ( fabs(xi) < win ) //if the agent is still inside the rewiring window
						{
							rand = randing->Random();

							if ( rand < y ) //with probability y we y-rewire
							{
								rewireY(i); //we y-rewire! cautiously over the A buffer

								if ( (inc != 0) && (inwY != 0) ) //if node i is actually y-rewiring
								{
									arrcheck = 1; //we need to recalculate the arrays
									rewcheck++; //we increase the rew counter
								}
							}
							else //with probability 1-y we rewire
							{
								rewire(i); //we rewire! cautiously over the A buffer

								if ( (inc != 0) && (inw != 0) ) //if node i is actually rewiring
								{
									arrcheck = 1; //we need to recalculate the arrays
									rewcheck++; //we increase the rew counter
								}
							}
						}
					}

					if ( arrcheck ) //we recalculate the arrays
					{
						for (int i = 1; i <= R; i++) //we first actualize the A array
							arrA[i - 1] = buffA[i - 1];

						IpowersA(arrA, 3); //we then get A^2 and A^3
						spath(); //and finally we get sp
					}

					//show the actual number of rewirings per generation, once in a while
					if ( (t % (10*tgen)) == 0 )
						std::cout << " \t\t< gen rewcheck > " << t/tgen << ", " << rewcheck << std::endl;

				}

//				//we construct the himmeli graph of an intermediate state of the network
//				if ( ( himmswitch == 1 ) && ( ( t % (2500000) ) == 0 ) )
//					himmeli(intstr);

				//we write the time series according to its stage and switch
				if ( (timeswitch == 1) && (wcount > 0) && ((t % twrite) == 0) )
					tseries(2);

				//now! regarding the end of the dynamics...


				if ( timeswitch ) //if we are writing time series...
				{
					if ( t == (twmax*twrite) )
						life = false;
				}
				else //or not
				{
					//each generation we rewcheck and calculate the standard deviation, as to review our Work
					if ( (t % tgen) == 0 )
					{
						if ( finalswitch ) //if the end is determined by rewiring dynamics
						{
							if ( (double)rewcheck <= rewzero ) //a few rewiring nodes
								life = false;
						}
						else //or instead by opinion dynamics
						{
							stdevbuff = squares();

							if ( fabs(stdevbuff - stdev) <= zero )
								life = false; //life has endeth...

							stdev = stdevbuff; //we actualize stdev
						}

						//if everything else should fail...
						if ( t == (12000*tgen) )
							life = false;
					}
				}

			} //life loop

			if ( percswitch ) //we percolate the network according to the corresponding switch
				perc(thres);

			//we get the final results for the present value of the parameter and nt
			ncount(one); //we get nb, ng, nw
			CZvars(one); //we get nclusters, cmax, avclu, suscep, meanz and ccoef
			length(one, 2); //we get plen and neigh for n=2

			//we print the final state of the system
			std::cout << std::endl << " final at | p nt | : < nb ng nw nclu cmax meanz ccoef suscep plen neigh avclu > " << "\n\t" << " | " << tgen << "  " << nt << " | : < " << one.nb << "  " << one.ng << "  " << one.nw << "  " << one.nclusters << "  " << one.cmax << "  " << one.meanz << "  " << one.ccoef << "  " << one.suscep << "  " << one.plen << "  " << one.neigh << "  " << one.avclu << " >" << std::endl;

			//and the euler check
			if ( eulcount )
			{
				eul = eulbuff/( (double)eulcount ); //we get the mean value of the added term in the euler method for this run
				std::cout << " \t\t< eul > " << eul << std::endl; //and show it
			}

			//we fill the accumulators
			acc.nb += one.nb;	acc.ng += one.ng;	acc.nw += one.nw;
			acc.nclusters += one.nclusters;
			acc.cmax += one.cmax;
			acc.avclu += one.avclu;
			acc.suscep += one.suscep;
			acc.meanz += one.meanz;
			acc.ccoef += one.ccoef;
			acc.plen += one.plen;
			acc.neigh += one.neigh;

			//we write the distros according to their stage and switch
			if ( distswitch )
			{
				cgraph(2);
				zcorr(2);
			}

			nt--; //we decrease the counter
		} //nt loop

		//we construct the himmeli graph of the final state of the network
		if ( himmswitch )
			himmeli(finstr);

		//now we get the stuffed accumulators, make averages
		Mean.nb = (acc.nb)/dntimes; Mean.ng = (acc.ng)/dntimes; Mean.nw = (acc.nw)/dntimes;
		Mean.nclusters = (acc.nclusters)/dntimes;
		Mean.cmax = (acc.cmax)/dntimes;
		Mean.avclu = (acc.avclu)/dntimes;
		Mean.suscep = (acc.suscep)/dntimes;
		Mean.meanz = (acc.meanz)/dntimes;
		Mean.ccoef = (acc.ccoef)/dntimes;
		Mean.plen = (acc.plen)/dntimes;
		Mean.neigh = (acc.neigh)/dntimes;

		//and write down in gnu
		fprintf(gnu, "%i %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", tgen, Mean.nb, Mean.ng, Mean.nw, Mean.nclusters, Mean.cmax, Mean.meanz, Mean.ccoef, Mean.suscep, Mean.plen, Mean.neigh, Mean.avclu);

		//we write the distros according to their stage and switch
		if ( distswitch )
		{
			cgraph(3);
			zcorr(3);
		}

		//we write the time series according to its stage and switch
		if ( timeswitch )
			tseries(3);
	} //parameter loop

	//we close the file descriptor
	fclose(gnu);

	//we update the init files according to their switchs
	if ( (mmapswitch == 1) && (initswitch == 0) && (lifeswitch == 0) )
	{
		char cpx[255];
		char cpa[255];
		char cpA[255];
		char cpsp[255];

		sprintf(cpx, "cp -p ./numenorN%iarrx ./numenorN%iarrxinit", N, N);
		sprintf(cpa, "cp -p ./numenorN%ialpha ./numenorN%ialphainit", N, N);
		sprintf(cpA, "cp -p ./numenorN%iarrA ./numenorN%iarrAinit", N, N);
		sprintf(cpsp, "cp -p ./numenorN%iarrsp ./numenorN%iarrspinit", N, N);

		system(cpx);
		system(cpa);
		system(cpA);
		system(cpsp);
	}

	//alas! the end of main
	return 0;
}

//last but no least, all the secondary functions used in main

//this function calculates m powers of the int array arr after each rewiring process
void IpowersA(pint& arr, int m)
{
	int posAord; //linear position for A^ord = A * A^n, n = ord-1
	int posA;    //linear position for A
	int posAn;   //linear position of A^n

	int Aordij;  //element of A^ord
	int Aiind;   //element of A
	int Anindj;  //element of A^n
	for (int ord = 2; ord <= m; ord++)
	{
		for (int i = 1; i <= N; i++)
		{
			for (int j = i; j <= N; j++)
			{
				posAord = Amat2lin(i,j); //we set the ord position
				Aordij = 0; //we initialize the ord element value to 0

				for (int ind = 1; ind <= N; ind++)
				{
					posA = Amat2lin(i, ind); //A and An positions
					posAn = Amat2lin(ind, j);

					Aiind = arr[posA - 1]; //A and An elements
					Anindj = arr[ (ord - 2)*R + posAn - 1];

					Aordij = Aordij + Aiind*Anindj; //we make the product

				}

				arr[ (ord - 1)*R + posAord - 1 ] = Aordij; //we fill the array

				if (i != j) //if we're outside the diagonal...
				{
					posAord = Amat2lin(j, i);
					arr[ (ord - 1)*R + posAord - 1 ] = Aordij; //we symmetrize
				}
			}
		}
	}
}

//this function constructs the connectivity matrix for a random network with mean connectivity z
void ranAmatrix(int z)
{
	double zN = (double) (z*N); //evidently, z*N
	int Nbond = (int) ( zN/2.0 ); //number of bonds in network
	int acc = 0; //counter, initialized to 0

	for (int i = 1; i <= R; i++ ) //just in case we first initialize the A array to 0
		arrA[i - 1] = 0;

	while ( acc < Nbond ) //while we still have bonds to allocate...
	{
		int m = randing->IRandom(1, N); //we choose random sites
		int n = randing->IRandom(1, N);
		if (m != n) //that are different
		{
			int pos = Amat2lin(m, n);
			if( arrA[pos - 1] == 0 ) //and then
			{
				arrA[pos - 1] = 1; //connecting people!

				pos = Amat2lin(n, m); //and symmetrize!
				arrA[pos - 1] = 1;

				acc++; //finally we increase the counter
			}
		}
	}
}

//this function finds the shortest path between all (i,j) pairs, and fills the sp array with them
void spath()
{
	int i; //nodes i,j,k
	int j;
	int k;
	int s; //trial distance for nodes j,k
	int posji; //matrix positions
	int posik;
	int posjk;

	for (int ind = 1; ind <= R; ind++) //we initialize the sp array
	{
		if ( arrA[ind - 1] == 1 )
			arrsp[ind - 1] = 1;
		else
			arrsp[ind - 1] = dmax;
	}

	for (i = 1; i <= N; i++)
	{
//		std::cout << "spath i = " << i << std::endl;

		for (j = 1; j <= N; j++)
		{
			posji = Amat2lin(j,i);

			if ( arrsp[posji - 1] < dmax )
			{
				for (k = 1; k <= N; k++)
				{
					posik = Amat2lin(i, k);

					if ( arrsp[posik - 1] < dmax )
					{
						posjk = Amat2lin(j, k);
						s = arrsp[posji - 1] + arrsp[posik - 1];

						if ( s < arrsp[posjk - 1] )
							arrsp[posjk - 1] = s;
					}
				}
			}
		}
	}
}

//this function computes the f1,f2 terms for agent i
void f1f2(int ia, double& f1, double& f2)
{
	int j; //agent j
	double xj; //x variable for agent j
	int n; //neighbouring order for i
	int posij; //matrix position
	double fk; //normalizing factor

	f1 = 0.0; //sum of x values (weighted by fk) in the m-vicinity of node i, initialized to 0
	f2 = 0.0; //sum of x values (weighted by fk) outside the m-vicinity of node i, initialized to 0

	for (j = 1; j <= N; j++) //we construct the f1,f2 terms
	{
		if ( j != ia ) //we review all nodes except i
		{
			xj = arrx[j - 1]; //we set the x for node j
			posij = Amat2lin(ia, j); //we set the matrix position
			n = arrsp[posij - 1]; //we get the neighbouring order for node j
			fk = 1.0/((double)n); //we get fk

			if ( n < dmax ) //check if nodes i,j are actually connected by a path
			{
				if ( n <= mcut ) //if we are in the m-vicinity of node i
					f1 += fk*xj; //we add terms to f1
				else //if we are in the rest of the network around node i
					f2 += fk*xj; //we add terms to f2
			}
		}
	}
}

//converting vegetarians!! into the midnight!

//this function computes the x dynamical term for agent i
double dynam(int ia)
{
	double xi; //x variable for agents i
	double signi; //sign of xi
	double f1; //sum of x values (weighted by fk) in the m-vicinity of node i
	double f2; //sum of x values (weighted by fk) outside the m-vicinity of node i
	double alphai; //alpha constant for agent i
	double dterm; //dynamical term

	xi = arrx[ia - 1]; //we set the x for node i

	if ( xi > 0.0 ) //we get the sign
		signi = 1.0;
	else if (xi == 0.0)
		signi = 0.0;
	else
		signi = -1.0;

	f1f2(ia, f1, f2); //we get the f1,f2 terms for node i
	alphai = alpha[ia - 1]; //we get the alpha constant for node i

	dterm = signi*f1*xi + alphai*f2 + h; //the dynamical term

	return dterm;
}

//this function calculates nb, ng and nw
void ncount(out& one)
{
	double xi; //x for agent i
	one.nb = 0.0; //we initialize all values to 0
	one.ng = 0.0;
	one.nw = 0.0;

	for (int i = 1; i <= N; i++)
	{
		xi = arrx[i - 1];
		if ( xi <= -lim ) //black!
			one.nb++;
		else if ( (xi > -lim) && (xi < lim) ) //gray!
			one.ng++;
		else if ( xi >= lim ) //white!
			one.nw++;
	}
}

//this function computes the characteristic path length and the mean number of nth-order neighbours for the network
void length(out& one, int n)
{
	int posij; //matrix position
	int spentry; //element in the sp array

	int sumlen = 0; //sum of the individual lengths, initialized to 0
	int pcounter = 0; //length counter, initialized to 0
	int ncounter = 0; //neighbour counter, initialized to 0

	//we sail through the whole sp matrix
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			if ( i != j ) //we keep away from the diagonal
			{
				posij = Amat2lin(i, j); //we set the position in the sp array
				spentry = arrsp[posij - 1]; //we get the sp element

				if ( spentry < dmax ) //if nodes i,j are connected by a path...
				{
					sumlen += spentry; //we sum the path lengths
					pcounter++; //we increase the counter
				}

				if ( spentry == n ) //if nodes i,j are nth-order neighbours...
					ncounter++; //we increase the counter
			}
		}
	}

	//we make averages, and we're done
	one.plen = ( (double)sumlen )/( (double)pcounter ); //the characteristic path length
	one.neigh = ( (double)ncounter )/( (double)N ); //the mean number of nth-order neighbours
}

//this function fills arrays incut and inwire with respect to agent i
void incount(int ia)
{
	int posij; //linear position for pair i,j in A
	double pij; //probabilities of breaking and rewiring with agent i
	double qij;

	//we initialize the counters
	inc = 0;
	inw = 0;

	//let us get on with it
	for (int j = 1; j <= N; j++)
	{
		incut[j - 1] = 0; //we initialize the arrays
		inwire[j - 1] = 0;

		posij = Amat2lin(ia, j);
		pij = pcut(ia, j); //we cut from similar people,
		qij = qwire(ia, j); //and rewire with opposing people to try and convince them

		if ( j != ia ) //to be sure we are not choosing the same agent...
		{
			if ( arrA[posij - 1] == 1 ) //check for Aij (connection)
			{
//				if ( arrA[R + posij - 1] == 0 ) //check for 1 - Tij (lack of triangles)
//				{
					if ( cutoff < pij ) //we cut over cutoff
					{
						incut[inc] = j; //we wish to cut bond with agent j;
						inc++; //we increase the counter
					}
//				}
			}
			else //check for 1 - Aij (disconnection)
			{
				if ( arrA[R + posij - 1] != 0 ) //check for Tij (presence of triangles)
				{
					if ( cutoff < qij ) //we rewire over cutoff
					{
						inwire[inw] = j; //we wish to rewire bond with agent j;
						inw++; //we increase the counter
					}
				}
			}
		}
	}
}

//this function fills arrays incut and inwireY with respect to agent i
void incountY(int ia)
{
	int posij; //linear position for pair i,j in A
	double pij; //probabilities of breaking and y-rewiring with agent i
	double qij;
	int n; //neighbouring order for i

	//we initialize the counters
	inc = 0;
	inwY = 0;

	//let us get on with it
	for (int j = 1; j <= N; j++)
	{
		incut[j - 1] = 0; //we initialize the arrays
		inwireY[j - 1] = 0;

		posij = Amat2lin(ia, j);
		pij = pcut(ia, j); //we cut from similar people,
		qij = qwireY(ia, j); //and y-rewire with opposing people to try and convince them
		n = arrsp[posij - 1]; //we get the neighbouring order for node j

		if ( j != ia ) //to be sure we are not choosing the same agent...
		{
			if ( arrA[posij - 1] == 1 ) //check for Aij (connection)
			{
//				if ( arrA[R + posij - 1] == 0 ) //check for 1 - Tij (lack of triangles)
//				{
					if ( cutoff < pij ) //we cut over cutoff
					{
						incut[inc] = j; //we wish to cut bond with agent j;
						inc++; //we increase the counter
					}
//				}
			}
			else //check for 1 - Aij (disconnection)
			{
				if ( n > 2 ) //check for longer shortest paths than 2
				{
					if ( cutoff < qij ) //we y-rewire over cutoff
					{
						inwireY[inwY] = j; //we wish to y-rewire bond with agent j;
						inwY++; //we increase the counter
					}
				}
			}
		}
	}
}

//this function performs the network rewiring due to agent i over the A buffer
void rewire(int ia)
{
	int j; //agents j,k
	int k;
	int posij; //linear position of i,j pair in A
	int posji; //transpose
	int posik; //linear position of i,k pair in A
	int poski; //transpose
	double pi; //probabilities of breaking and rewiring with agent i
	double qi;
	int nchange; //number of changes to be made (the minimum of inw and inc)

	double storedq; //various variables needed for storage
	double storedp;
	int storedk;
	int storedind;

	incount(ia); //we fill the incut and inwire arrays

	//rewire the society!

	if ( inw > inc ) //more wiring than breaking
	{
		nchange = inc; //we set the counter

		//first we break!
		for (int ind = 1; ind <= N; ind++)
		{
			if ( incut[ind - 1] )
			{
				j = incut[ind - 1];
				posij = Amat2lin(ia, j);
				posji = Amat2lin(j, ia);

				buffA[posij - 1] = 0; //break!
				buffA[posji - 1] = 0; //symmetrize!

				arrP[posij - 1] = -1; //we fill the cut matrix P
			}
		}

		//and then we rewire!
		while (nchange)
		{
			storedq = 0.0;
			for (int ind = 1; ind <= N; ind++)
			{
				if ( inwire[ind - 1] )
				{
					k = inwire[ind - 1];
					qi = qwire(ia, k);
					if ( qi > storedq ) //we pick the agent k with greatest q
					{
						storedq = qi; //we update the stored value of q
						storedk = k; //we store the value of k and its place in inwire
						storedind = ind;
					}
				}
			}
			k = storedk;
			posik = Amat2lin(ia, k);
			poski = Amat2lin(k, ia);

			buffA[posik - 1] = 1; //rewire!
			buffA[poski - 1] = 1; //symmetrize!

			arrQ[posik - 1] = 1; //we fill the wire matrix Q

			inwire[storedind - 1] = 0; //to make sure we are not repeating agents
			nchange--; //we decrease the counter
		}
	}
	else if ( inw == inc ) //same wiring than breaking
	{
		//we break and rewire!
		for (int ind = 1; ind <= N; ind++)
		{
			if ( incut[ind - 1] )
			{
				j = incut[ind - 1];
				posij = Amat2lin(ia, j);
				posji = Amat2lin(j, ia);

				buffA[posij - 1] = 0; //break!
				buffA[posji - 1] = 0; //symmetrize!

				arrP[posij - 1] = -1; //we fill the cut matrix P
			}

			if ( inwire[ind - 1] )
			{
				j = inwire[ind - 1];
				posij = Amat2lin(ia, j);
				posji = Amat2lin(j, ia);

				buffA[posij - 1] = 1; //rewire!
				buffA[posji - 1] = 1; //symmetrize!

				arrQ[posij - 1] = 1; //we fill the wire matrix Q
			}
		}
	}
	else // ( inw < inc ) less wiring than breaking
	{
		nchange = inw; //we set the counter

		//we break!
		while (nchange)
		{
			storedp = 0.0;
			for (int ind = 1; ind <= N; ind++)
			{
				if ( incut[ind - 1] )
				{
					k = incut[ind - 1];
					pi = pcut(ia, k);
					if ( pi > storedp ) //we pick the agent k with greatest p
					{
						storedp = pi; //we update the stored value of p
						storedk = k; //we store the value of k and its place in incut
						storedind = ind;
					}
				}
			}
			k = storedk;
			posik = Amat2lin(ia, k);
			poski = Amat2lin(k, ia);

			buffA[posik - 1] = 0; //break!
			buffA[poski - 1] = 0; //symmetrize!

			arrP[posik - 1] = -1; //we fill the cut matrix P

			incut[storedind - 1] = 0; //to make sure we are not repeating agents
			nchange--; //we decrease the counter
		}

		//and then we rewire!
		for (int ind = 1; ind <= N; ind++)
		{
			if ( inwire[ind - 1] )
			{
				j = inwire[ind - 1];
				posij = Amat2lin(ia, j);
				posji = Amat2lin(j, ia);

				buffA[posij - 1] = 1; //rewire!
				buffA[posji - 1] = 1; //symmetrize!

				arrQ[posij - 1] = 1; //we fill the wire matrix Q
			}
		}
	}
}

//this function performs the network y-rewiring due to agent i over the A buffer
void rewireY(int ia)
{
	int j; //agents j,k
	int k;
	int posij; //linear position of i,j pair in A
	int posji; //transpose
	int posik; //linear position of i,k pair in A
	int poski; //transpose
	double pi; //probabilities of breaking and y-rewiring with agent i
	double qi;
	int nchange; //number of changes to be made (the minimum of inwY and inc)

	double storedq; //various variables needed for storage
	double storedp;
	int storedk;
	int storedind;

	incountY(ia); //we fill the incut and inwireY arrays

	//y-rewire the society!

	if ( inwY > inc ) //more y-wiring than breaking
	{
		nchange = inc; //we set the counter

		//first we break!
		for (int ind = 1; ind <= N; ind++)
		{
			if ( incut[ind - 1] )
			{
				j = incut[ind - 1];
				posij = Amat2lin(ia, j);
				posji = Amat2lin(j, ia);

				buffA[posij - 1] = 0; //break!
				buffA[posji - 1] = 0; //symmetrize!

				arrP[posij - 1] = -1; //we fill the cut matrix P
			}
		}

		//and then we y-rewire!
		while (nchange)
		{
			storedq = 0.0;
			for (int ind = 1; ind <= N; ind++)
			{
				if ( inwireY[ind - 1] )
				{
					k = inwireY[ind - 1];
					qi = qwireY(ia, k);
					if ( qi > storedq ) //we pick the agent k with greatest q
					{
						storedq = qi; //we update the stored value of q
						storedk = k; //we store the value of k and its place in inwire
						storedind = ind;
					}
				}
			}
			k = storedk;
			posik = Amat2lin(ia, k);
			poski = Amat2lin(k, ia);

			buffA[posik - 1] = 1; //y-rewire!
			buffA[poski - 1] = 1; //symmetrize!

			arrQ[posik - 1] = 1; //we fill the wire matrix Q

			inwireY[storedind - 1] = 0; //to make sure we are not repeating agents
			nchange--; //we decrease the counter
		}
	}
	else if ( inwY == inc ) //same y-wiring than breaking
	{
		//we break and y-rewire!
		for (int ind = 1; ind <= N; ind++)
		{
			if ( incut[ind - 1] )
			{
				j = incut[ind - 1];
				posij = Amat2lin(ia, j);
				posji = Amat2lin(j, ia);

				buffA[posij - 1] = 0; //break!
				buffA[posji - 1] = 0; //symmetrize!

				arrP[posij - 1] = -1; //we fill the cut matrix P
			}

			if ( inwireY[ind - 1] )
			{
				j = inwireY[ind - 1];
				posij = Amat2lin(ia, j);
				posji = Amat2lin(j, ia);

				buffA[posij - 1] = 1; //y-rewire!
				buffA[posji - 1] = 1; //symmetrize!

				arrQ[posij - 1] = 1; //we fill the wire matrix Q
			}
		}
	}
	else // ( inwY < inc ) less y-wiring than breaking
	{
		nchange = inwY; //we set the counter

		//we break!
		while (nchange)
		{
			storedp = 0.0;
			for (int ind = 1; ind <= N; ind++)
			{
				if ( incut[ind - 1] )
				{
					k = incut[ind - 1];
					pi = pcut(ia, k);
					if ( pi > storedp ) //we pick the agent k with greatest p
					{
						storedp = pi; //we update the stored value of p
						storedk = k; //we store the value of k and its place in incut
						storedind = ind;
					}
				}
			}
			k = storedk;
			posik = Amat2lin(ia, k);
			poski = Amat2lin(k, ia);

			buffA[posik - 1] = 0; //break!
			buffA[poski - 1] = 0; //symmetrize!

			arrP[posik - 1] = -1; //we fill the cut matrix P

			incut[storedind - 1] = 0; //to make sure we are not repeating agents
			nchange--; //we decrease the counter
		}

		//and then we y-rewire!
		for (int ind = 1; ind <= N; ind++)
		{
			if ( inwireY[ind - 1] )
			{
				j = inwireY[ind - 1];
				posij = Amat2lin(ia, j);
				posji = Amat2lin(j, ia);

				buffA[posij - 1] = 1; //y-rewire!
				buffA[posji - 1] = 1; //symmetrize!

				arrQ[posij - 1] = 1; //we fill the wire matrix Q
			}
		}
	}
}

//this function computes the F,f average variables after each rewiring
void theFvars(out& one)
{
	int plusacc, minusacc, totacc; //the accumulators for the F,f variables
	int posij, posji; //matrix positions

	//we firstly set the rewire matrix R with the filled P,Q arrays
	for (int i = 1; i <= R; i++)
		arrR[i - 1] = arrP[i - 1] + arrQ[i - 1];

	//we initialize the F,f variables
	one.Fplus = 0.0;	one.Fminus = 0.0;
	one.fplus = 0.0;	one.fminus = 0.0;

	//we review each potential node in the rewiring
	for (int i = 1; i <= N; i++)
	{
		plusacc = 0; //we initialize the accumulators for each rewiring node i
		minusacc = 0;
		totacc = 0;

		//we find the diagonal elements of the squared rewire matrices
		for (int j = 1; j <= N; j++)
		{
			if (i != j) //choose different agents
			{
				posij = Amat2lin(i, j); //we get the matrix positions
				posji = Amat2lin(j, i);

				if ( (arrP[posij - 1] != 0) && (arrP[posji - 1] != 0) ) //if the diagonal element's component of P^2 is non-zero
					plusacc++;

				if ( (arrQ[posij - 1] != 0) && (arrQ[posji - 1] != 0) ) //if the diagonal element's component of Q^2 is non-zero
					minusacc++;

				if ( (arrR[posij - 1] != 0) && (arrR[posji - 1] != 0) ) //if the diagonal element's component of R^2 is non-zero
					totacc++;
			}
		}

		if ( totacc ) //if there is something to add to the F,f variables...
		{
			one.Fplus += (double)plusacc;				one.Fminus += (double)minusacc;
			one.fplus += ( (double)plusacc/(double)totacc );	one.fminus += ( (double)minusacc/(double)totacc );
		}
		else //and if not...
		{
			one.fplus += 0.5;	one.fminus += 0.5;
		}
	}

	//finally, we average the F,f variables over this rewiring
	one.Fplus /= (double)N;	one.Fminus /= (double)N;
	one.fplus /= (double)N;	one.fminus /= (double)N;
}

//this function percolates the network over a certain threshold thres dependent on the dynamical variables
void perc(double thres)
{
	int posij, posji; //linear positions for pair i,j in A
	double pij; //probability of percolating nodes i,j
	int acc = 0; //counter, initialized to 0

	//we review the upper triangle of the A matrix
	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j <= N; j++)
		{
			posij = Amat2lin(i, j);
			posji = Amat2lin(j, i);
			pij = pcut(i, j); //we compute the percolation probability

			if ( (arrA[posij - 1] == 1) && (pij >= thres) ) //if nodes i,j are connected and the threshold is surpassed...
			{
				arrA[posij - 1] = 0; //we percolate
				arrA[posji - 1] = 0; //we symmetrize

				acc++; //we increase the counter
			}
		}
	}

	if ( acc ) //if there are changes we have to recalculate arrays
	{
		IpowersA(arrA, 3); //we get A^2 and A^3
		spath(); //and we get sp
	}
}

//this function searches the network for clusters and its sizes
void cluster()
{
	int Xi; //number of non-zero elements in the Li, Oi arrays
	int Yi;
	int nclu = 0; //we start with no clusters counted
	int u; //random agent to be exposed
	int pos; //matrix position

	for (int i = 1; i <= N; i++)
	{
		clust[i - 1] = 0;
		ccheck[i - 1] = 1; //at the beginning we can pick any node
	}

	for (int v = 1; v <= N; v++)
	{
		if ( ccheck[v - 1] != 0 ) //if we haven't picked this agent yet
		{
			for (int i = 1; i <= N; i++)
			{
				Li[i - 1] = 0; //we initialize everything to 0
				Oi[i - 1] = 0;
			}

			Li[v - 1] = v; //the initial V set
			ccheck[v - 1] = 0; //node v has been checked

			Oi[v - 1] = v; //the initial O set

			Xi = count(Li); //the initial cardinality of the V,O sets, namely 1
			Yi = count(Oi);

			while ( Yi ) //until there is nothing left to be exposed
			{
				u = randing->IRandom(1, N); //choose randomly
				if ( Oi[u - 1] != 0 ) //check if u is a live vertex
				{
					Oi[u - 1] = 0; //we kill u!
					for (int n = 1; n <= N; n++) //we sail through the neighbours of u
					{
						pos = Amat2lin(u, n);
						if ( arrA[pos - 1] != 0 )
						{
							if (Li[n - 1] == 0) //we mark n live if it is not in Li
								Oi[n - 1] = n; //we actualize the O set

							Li[n - 1] = n; //we actualize the L set
							ccheck[n - 1] = 0; //node n has been checked
						}
					}

					Xi = count(Li); //the compute the new cardinality of the V,O sets
					Yi = count(Oi);
				}
			}

			clust[nclu] = Xi; //the size of the cluster
			nclu++; //one more cluster
		}
	}
}

//this function computes the c,z variables in the out struct, namely nclusters, cmax, avclu, suscep, meanz and ccoef, and all the correlations
void CZvars(out& one)
{
	double dN = (double) N; //double version of N

	int nclu; //number of clusters
	int csize; //cluster size
	int csmax; //maximum cluster size
	double csmax2; //square of the LCC size, that is to say csmax^2

	int posii; //matrix position
	int zi; //connectivity of node i
	int mzi; //mean z counter for node i
	double ccimax; //(double of) maximum number of bonds between neighbours of node i, ki(ki-1) with ki=(A^2)ii
	double ni; //(double of) number of bonds between neighbours of node i, (A^3)ii
	double cci; //clustering coefficient counter for node i, ni/ccimax

	//first the c variables

	cluster(); //we find the clusters and fill the clust array
	nclu = count(clust); //we get nclusters

	for (int c = 1; c <= nclu; c++) //we search for cmax
	{
		csize = clust[c - 1]; //we get the corresponding size element

		if ( c == 1 ) csmax = csize; //we initialize cmax
		else
		{
			if (csize > csmax)
				csmax = csize; //we eventually get cmax
		}
	}

	one.nclusters = (double) nclu; //the number of clusters
	one.cmax = (double) csmax; //the maximum cluster size
	one.avclu = ( (double)sum(clust) )/one.nclusters; //the average cluster size

	if ( nclu > 1 ) //if there is something apart from the LCC (largest connected cluster)
	{
		csmax2 = pow( one.cmax, 2.0 );
		//the susceptibility excludes the LCC...
		one.suscep = ( sum2(clust) - csmax2 )/( (double)(sum(clust) - csmax) );
	}
	else
		one.suscep = 0.0; //we define the undefined susceptibility as zero (go figure...)


	//then the z variables and the correlations

	mzi = 0; //we initialize the counters
	cci = 0.0;

	for (int i = 1; i <= N; i++)
	{
		posii = Amat2lin(i, i);
		zi = arrA[R + posii - 1]; //the connectivity of agent i is element (A^2)ii

		mzi += zi; //we fill the mz counter

		if ( zi > 1 ) //if there is something to find, we get the clustering coefficient for agent i
		{
			ccimax = (double) ( zi*(zi - 1) );
			ni = (double) arrA[2*R + posii - 1];
			cci += ni/ccimax; //we fill the cc counter
		}
	}

	one.meanz = ( (double)mzi )/dN; //we get the mean z
	one.ccoef = cci/dN; //we get the mean clustering coefficient for the network
}

//this function constructs the c distribution graph and writes it down in cdib according to
//stage = 1,2,3, which corresponds respectively to setup, writing of data and construction of graph
void cgraph(int stage)
{
	int nclu; //number of clusters in network
	int i; //node i
	int ci; //cluster i
	int si; //size of cluster i

	if ( stage == 1 ) //setup
	{
		//we open the file descriptor for this value of the parameter
		sprintf(cdibf2w, "./cdibN-N%i-p%i.txt", N, z);
		if((cdib = fopen(cdibf2w, "w+")) == (FILE *) NULL )
		{
			perror("\ncant open file to write\n");
			exit(1);
		}

		for (i = 1; i <= N; i++)
			valc[i - 1] = 0.0; //we initialize the array for the average number of clusters with size i

		ccount = 0; //we initialize the c counter
	}
	else if ( stage == 2 ) //writing of data
	{
		cluster(); //we find the clusters and fill the clust array
		nclu = count(clust); //we initialize the number of clusters

		for (ci = 1; ci <= nclu; ci++)
		{
			si = clust[ci - 1]; //we get the size of cluster i
			valc[si - 1] += 1.0; //we increase the terms in the c array
		}

		ccount++; //we increase the c counter
	}
	else // ( stage == 3 ) // construction of graph
	{
		if ( ccount ) //if there is something to do...
		{
			for (si = 1; si <= N; si++)
			{
				valc[si - 1] /= ( (double)ccount ); //we normalize the values in the c array

				if ( valc[si - 1] >= 0.0001 ) //print relevant results
					fprintf(cdib, "%i %.4f\n", si, valc[si - 1]); //we write the array elements in cdib
			}
		}

		//we close the file descriptor for this value of the parameter
		fclose(cdib);
	}
}

//this function constructs the z correlation graph that determines the assortativity of the network and writes it down in zcor according to
//stage = 1,2,3, which corresponds respectively to setup, writing of data and construction of graph
void zcorr(int stage)
{
	int i, j; //nodes i,j
	int posii, posij, posjj; //matrix positions
	int zi; //connectivity of node i
	int mzj; //mean z counter for node j
	int zcount; //counter for the neighbours of node i
	double znn; //mean z for the neighbours of node i, mzj/zcount

	if ( stage == 1 ) //setup
	{
		//we open the file descriptor for this value of the parameter
		sprintf(zcorf2w, "./zcorN-N%i-p%i.txt", N, z);
		if((zcor = fopen(zcorf2w, "w+")) == (FILE *) NULL )
		{
			perror("\ncant open file to write\n");
			exit(1);
		}

		for (i = 1; i <= N; i++)
		{
			valznn[i - 1] = 0.0; //we initialize the arrays of values of znn, and the znn counter array
			znncount[i - 1] = 0;
		}
	}
	else if ( stage == 2 ) //writing of data
	{
		for (i = 1; i <= N; i++)
		{
			posii = Amat2lin(i, i);
			zi = arrA[R + posii - 1]; //the connectivity of agent i is element (A^2)ii

			if ( zi ) //if there is something to do
			{
				mzj = 0; //we initialize the counters
				zcount = 0;

				for (j = 1; j <= N; j++) //we review the neighbours...
				{
					posij = Amat2lin(i, j);

					if ( arrA[posij - 1] ) //if nodes i,j are actually neighbours
					{
						posjj = Amat2lin(j, j);
						mzj += arrA[R + posjj - 1]; //we fill the counters
						zcount++;
					}
				}

				znn = ( (double)mzj )/( (double)zcount ); //we get the mean z of the neighbours of node i

				valznn[zi - 1] += znn; //finally, we increase the terms in the the znn array
				znncount[zi - 1]++; //and increase the counter
			}
		}
	}
	else // ( stage == 3 ) // construction of graph
	{
		for (zi = 1; zi <= N; zi++)
		{
			if ( znncount[zi - 1] ) //if there are elements to sort out...
			{
				valznn[zi - 1] /= ( (double)znncount[zi - 1] ); //we normalize the values in the znn array

				fprintf(zcor, "%i %.4f\n", zi, valznn[zi - 1]); //we write the array elements in zcor
			}
		}

		//we close the file descriptor for this value of the parameter
		fclose(zcor);
	}
}

//this function computes the average time series according to
//stage = 1,2,3, which corresponds respectively to setup, writing of data and construction of graph
void tseries(int stage)
{
	int gen; //generation

	if ( stage == 1 ) //setup
	{
		//we open the file descriptor for this value of the parameter
		sprintf(tserf2w, "./tser-N%i-p%i.txt", N, z);
		if((tser = fopen(tserf2w, "w+")) == (FILE *) NULL )
		{
			perror("\ncant open file to write\n");
			exit(1);
		}

		//we allocate memory for the accumulator array, a struct of type out for each time step 0,...,twmax
		Tacc = (out *) calloc(twmax + 1, sizeof(out));

		if (Tacc == NULL) //if something goes wrong...
		{
			printf ("error allocating memory\n");
			//close all open files and free any allocated memory
			exit (1);
		}
	}
	else if ( stage == 2 ) //writing of data
	{
		wcount--; //we decrease the write counter

		gen = twmax - wcount; //we set the generation

		//we fill the Tone struct
		ncount(Tone); //we get nb, ng, nw
		CZvars(Tone); //we get nclusters, cmax, avclu, suscep, meanz and ccoef
		length(Tone, 2); //we get plen and neigh for n=2
		theFvars(Tone); //we get the F,f rewiring variables

		//we print the current state of the system
//		std::cout << std::endl << "at gen < nb ng nw nclu cmax meanz ccoef suscep plen neigh avclu Fplus Fminus fplus fminus rewcheck > " << "\n\t" << gen << " < " << Tone.nb << "  " << Tone.ng << "  " << Tone.nw << "  " << Tone.nclusters << "  " << Tone.cmax << "  " << Tone.meanz << "  " << Tone.ccoef << "  " << Tone.suscep << "  " << Tone.plen << "  " << Tone.neigh << "  " << Tone.avclu << "  " << Tone.Fplus << "  " << Tone.Fminus << "  " << Tone.fplus << "  " << Tone.fminus << "  " << rewcheck << " >" << std::endl;

		//and fill the entry gen in the accumulator array
		Tacc[gen].nb += Tone.nb;			Tacc[gen].ng += Tone.ng;		Tacc[gen].nw += Tone.nw;
		Tacc[gen].nclusters += Tone.nclusters;
		Tacc[gen].cmax += Tone.cmax;
		Tacc[gen].avclu += Tone.avclu;
		Tacc[gen].suscep += Tone.suscep;
		Tacc[gen].meanz += Tone.meanz;
		Tacc[gen].ccoef += Tone.ccoef;
		Tacc[gen].plen += Tone.plen;
		Tacc[gen].neigh += Tone.neigh;
		Tacc[gen].Fplus += Tone.Fplus;			Tacc[gen].Fminus += Tone.Fminus;
		Tacc[gen].fplus += Tone.fplus;			Tacc[gen].fminus += Tone.fminus;
		Tacc[gen].rewnodes += rewcheck;
	}
	else // ( stage == 3 ) //construction of graph
	{
		//we get the stuffed accumulators and make averages for all generations
		for (gen = 0; gen <= twmax; gen++)
		{

			Tmean.nb = (Tacc[gen].nb)/dntimes;
			Tmean.ng = (Tacc[gen].ng)/dntimes;
			Tmean.nw = (Tacc[gen].nw)/dntimes;
			Tmean.nclusters = (Tacc[gen].nclusters)/dntimes;
			Tmean.cmax = (Tacc[gen].cmax)/dntimes;
			Tmean.avclu = (Tacc[gen].avclu)/dntimes;
			Tmean.suscep = (Tacc[gen].suscep)/dntimes;
			Tmean.meanz = (Tacc[gen].meanz)/dntimes;
			Tmean.ccoef = (Tacc[gen].ccoef)/dntimes;
			Tmean.plen = (Tacc[gen].plen)/dntimes;
			Tmean.neigh = (Tacc[gen].neigh)/dntimes;
			Tmean.Fplus = (Tacc[gen].Fplus)/dntimes;	Tmean.Fminus = (Tacc[gen].Fminus)/dntimes;
			Tmean.fplus = (Tacc[gen].fplus)/dntimes;	Tmean.fminus = (Tacc[gen].fminus)/dntimes;
			Tmean.rewnodes = (Tacc[gen].rewnodes)/dntimes;

			//and finally, we write down in tser
			fprintf(tser, "%i %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", gen, Tmean.nb, Tmean.ng, Tmean.nw, Tmean.nclusters, Tmean.cmax, Tmean.meanz, Tmean.ccoef, Tmean.suscep, Tmean.plen, Tmean.neigh, Tmean.avclu, Tmean.Fplus, Tmean.Fminus, Tmean.fplus, Tmean.fminus, Tmean.rewnodes);
		}

		//we close the file descriptor for this value of the parameter
		fclose(tser);
	}
}

//for himmeli and beyond!
void himmeli(char * str)
{
	//we open the 3 himmeli files
	sprintf(configf2w, "./config.txt");

	if((config = fopen(configf2w, "w+")) == (FILE *) NULL )
	{
		perror("\ncant open config file to write\n");
		exit(1);
	}

	sprintf(edgesf2w, "./edges.txt");

	if((edges = fopen(edgesf2w, "w+")) == (FILE *) NULL )
	{
		perror("\ncant open edges file to write\n");
		exit(1);
	}

	sprintf(verticesf2w, "./vertices.txt");

	if((vertices = fopen(verticesf2w, "w+")) == (FILE *) NULL )
	{
		perror("\ncant open vertices file to write\n");
		exit(1);
	}

	//we write down the config file
	fprintf(config, "GraphName                graph-N%i-p%i%s\n", N, tgen, str);
	fprintf(config, "EdgeFile                 edges.txt\n");
	fprintf(config, "VertexFile               vertices.txt\n");
	fprintf(config, "EdgeHeadVariable         HEAD\n");
	fprintf(config, "EdgeTailVariable         TAIL\n");
	fprintf(config, "VertexNameVariable       NAME\n");
	fprintf(config, "VertexColorVariable      COLOR\n");
	fprintf(config, "VertexShapeVariable      SHAPE\n");
	fprintf(config, "VertexColorInfo          white node             999999\n");
	fprintf(config, "VertexColorInfo          black node             000000\n");
	fprintf(config, "VertexShapeInfo          decided agent          circle\n");
	fprintf(config, "VertexShapeInfo          undecided agent        square\n");
	fprintf(config, "LabelMode                off\n");
	fprintf(config, "DistanceUnit             0.9\n");

	//we write down the header of the edges and vertices files
	fprintf(edges, "HEAD\t\tTAIL\n");
	fprintf(vertices, "NAME\t\tCOLOR\t\tSHAPE\n");

	int posij; //linear position of pair i,j in A
	double xi; //x for agent i
	double dobcol; //colour double version
	int intcol; //colour int version

	for (int i = 1; i <= N; i++)
	{
		//first the edges
		if ( i < N )
		{
			for (int j = i + 1; j <= N; j++)
			{
				posij = Amat2lin(i, j);
				if ( arrA[posij - 1] != 0 )
					fprintf(edges, "v%i\t\tv%i\n", i, j);
			}
		}

		//and then the vertices
		xi = arrx[i - 1]; //we set xi

		if ( xi == lim ) //white agents -> white!
			fprintf(vertices, "v%i\t\t999999\t\tcircle\n", i);
		else if ( xi == -lim ) //black agents -> black!
			fprintf(vertices, "v%i\t\t000000\t\tcircle\n", i);
		else //the rest of them
		{
			dobcol = 50*(xi + 1); //renormalise opinion
			intcol = (int) floor(dobcol); //get colour coding

			if ( intcol < 10 )
				fprintf(vertices, "v%i\t\t0%i0%i0%i\t\tsquare\n", i, intcol, intcol, intcol);
			else
				fprintf(vertices, "v%i\t\t%i%i%i\t\tsquare\n", i, intcol, intcol, intcol);
		}

//ok, if you really want discrete colouring...
//		if ( xi >= lim ) //white agents -> white!
//			fprintf(vertices, "v%i\t\t999999\n", i);
//		else if ( xi <= -lim ) //black agents -> black!
//			fprintf(vertices, "v%i\t\t000000\n", i);
//		else if ( (xi > -lim) && (xi <= -0.1) ) //black gray agents -> blue!
//			fprintf(vertices, "v%i\t\t000099\n", i);
//		else if ( (xi > -0.1) && (xi < 0.1) ) //the ones in the middle -> grey!
//			fprintf(vertices, "v%i\t\t505050\n", i);
//		else  //( (xi >= 0.1) && (xi < lim) ) // white gray agents -> red!
//			fprintf(vertices, "v%i\t\t990000\n", i);

	}

	//we close the 3 himmeli files
	fclose(config);
	fclose(edges);
	fclose(vertices);

	//we run himmeli
	char himm[255];

	sprintf(himm, "himmeli config.txt");

	system(himm);
}

//this function initializes n double arrays arr with name str and size len
void DarrayInit(pdouble& arr, int len, int n, char * str)
{
	//we open the file to be read
	char path[255];

	sprintf(path, "./numenorN%i%sinit", N, str);

	//joder! (clap)(clap)

	int initfd;

	initfd = open(path, O_RDONLY);

	if(initfd == -1)
	{
		perror("\ncant open init file\n");
		exit(1);
	}

	//we read from it and write it down in the array, one double at a time
	double buffer[1];
	int rlen = n*len;

	for (int i = 0; i < rlen; i++)
	{
		int ret;
		ret = read(initfd, (void *)buffer, sizeof(double));

		if(ret == -1)
		{
			perror("\ncant read init file\n");
			exit(1);
		}

		arr[i] = buffer[0];
	}

	//and finally we close the file
	int ret;

	ret = close(initfd);

	if(ret == -1)
		{
			perror("\ncant close init file\n");
			exit(1);
		}
}

//this function initializes n int arrays arr with name str and size len
void IarrayInit(pint& arr, int len, int n, char * str)
{
	//we open the file to be read
	char path[255];

	sprintf(path, "./numenorN%i%sinit", N, str);

	//joder! (clap)(clap)

	int initfd;

	initfd = open(path, O_RDONLY);

	if(initfd == -1)
	{
		perror("\ncant open init file\n");
		exit(1);
	}

	//we read from it and write it down in the array, one int at a time
	int buffer[1];
	int rlen = n*len;

	for (int i = 0; i < rlen; i++)
	{
		int ret;
		ret = read(initfd, (void *)buffer, sizeof(int));

		if(ret == -1)
		{
			perror("\ncant read init file\n");
			exit(1);
		}

		arr[i] = buffer[0];
	}

	//and finally we close the file
	int ret;

	ret = close(initfd);

	if(ret == -1)
		{
			perror("\ncant close init file\n");
			exit(1);
		}
}

//this function does the memory mmaping of n double arrays arr with name str, size len and file descriptor filedesc.
void Dmmaparray(pdouble& arr, int len, int n, char * str, int& filedesc)
{
	char path[255];

	sprintf(path, "./numenorN%i%s", N, str);

	//joder! (clap)(clap)

	filedesc = open(path, O_CREAT|O_RDWR);

	if(filedesc == -1)
	{
		perror("\ncant open file to mmap the array\n");
		exit(1);
	}

	//we change the permissions of the file
	char chmod[255];

	sprintf(chmod, "chmod 0666 %s", path);

	system(chmod);

	//we stuff the array with 0's
	unsigned char stuffing[1] = {0x00}; //0000 0000
	int rlen = n*len*sizeof(double);

	for(int i = 0; i < rlen; i++)
	{
		int ret;
		ret = write(filedesc, (void *)stuffing, 1);

		if(ret == -1)
		{
			perror("\ncant write to mmap file\n");
			exit(1);
		}
	}

	int ret;
	ret = lseek(filedesc, SEEK_SET, 0);

	if(ret == -1)
	{
		perror("\ncant go to the beggining of the mmap file\n");
		exit(1);
	}

	//dany dice que soy un mal companero de cama
	//lo soy?
	//joder! (clap)(clap)

	//memory mapping of the array

	arr = (pdouble) mmap((void *)0, rlen, PROT_READ|PROT_WRITE, MAP_SHARED, filedesc, 0);

	if(arr == (void *)-1)
	{
		perror("\ncant mmap array\n");
		exit(1);
	}
}

//this function does the memory mmaping of n int arrays arr with name str, size len and file descriptor filedesc.
void Immaparray(pint& arr, int len, int n, char * str, int& filedesc)
{
	char path[255];

	sprintf(path, "./numenorN%i%s", N, str);

	//joder! (clap)(clap)

	filedesc = open(path, O_CREAT|O_RDWR);

	if(filedesc == -1)
	{
		perror("\ncant open file to mmap the array\n");
		exit(1);
	}

	//we change the permissions of the file
	char chmod[255];

	sprintf(chmod, "chmod 0666 %s", path);

	system(chmod);

	//we stuff the array with 0's
	unsigned char stuffing[1] = {0x00}; //0000 0000
	int rlen = n*len*sizeof(int);

	for(int i = 0; i < rlen; i++)
	{
		int ret;
		ret = write(filedesc, (void *)stuffing, 1);

		if(ret == -1)
		{
			perror("\ncant write to mmap file\n");
			exit(1);
		}
	}

	int ret;
	ret = lseek(filedesc, SEEK_SET, 0);

	if(ret == -1)
	{
		perror("\ncant go to the beggining of the mmap file\n");
		exit(1);
	}

	//dany dice que soy un mal companero de cama
	//lo soy?
	//joder! (clap)(clap)

	//memory mapping of the array

	arr = (pint) mmap((void *)0, rlen, PROT_READ|PROT_WRITE, MAP_SHARED, filedesc, 0);

	if(arr == (void *)-1)
	{
		perror("\ncant mmap array\n");
		exit(1);
	}
}
