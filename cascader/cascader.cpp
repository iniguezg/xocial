/*	Cascader v1.0
 *	A simple model to describe complex contagion on techno-social networks
 *	Copyright (C) 2018 Gerardo Iñiguez
 *
 *	This program is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * If you use this program for your research, please cite any relevant from the following articles:
 *
 * M. Karsai, G. Iñiguez, K. Kaski, J. Kertész
 * Complex contagion process in spreading of online innovation
 * Journal of the Royal Society Interface 11, 20140694 (2014)
 * DOI: 10.1098/rsif.2014.0694
 * arXiv: 1405.6879

 * Z. Ruan, G. Iñiguez, M. Karsai, J. Kertész
 * Kinetics of social contagion
 * Physical Review Letters 115, 218702 (2015)
 * DOI: 10.1103/PhysRevLett.115.218702
 * arXiv: 1506.00251

 * M. Karsai, G. Iñiguez, R. Kikas, K. Kaski, J. Kertész
 * Local cascades induced global contagion: How heterogeneous thresholds, exogenous effects, and unconcerned behaviour govern online adoption spreading
 * Scientific Reports 6, 27178 (2016)
 * DOI: 10.1038/srep27178
 * arXiv: 1601.07995
 *
 * SUMMARY
 *
 * Implementation of a model for global cascades in complex networks: Threshold model in a uniform random graph with single threshold, rate of innovators and immune nodes
 *
 * Details:
 *
 * - Network has size N and average degree z, and all agents have threshold phi
 * - Innovators are selected randomly either as an initial multiseed (ini_innov), or uniformly in time with rate pnew
 * - There is a fraction r of immunes nodes
 * - {vulnerable, stable} nodes have int_thres = { 1, >1 }
 *
 * Conditions:
 *
 * - Static multiseed model: pnew = 0, ini_innov >  0
 * - Dynamic multiseed model: pnew > 0, ini_innov >= 0
 * - In general: ini_innov <= N * (1 - r)
 */

#include <iostream>
#include <exception>
#include <algorithm>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_rng.h>


using namespace std;

//constants for strings and files
#define STR_MAX 255  //max string/line length
#define NUM_VAR 21   //number of variables
#define VNAME_MAX 20 //max name length
#define BIG_NETS 0   //using small (0) or BIG (1) networks

//types for function declarations
typedef bool * pbool;
typedef int * pint;
typedef double * pdouble;
typedef FILE * pfile;

//types for data storage
#if BIG_NETS==0
typedef int lint;             //size for microtime of adoption
typedef short int sint;       //size for adoption timeseries
#else
typedef long long int lint; //or for BIG networks...
typedef int sint;
#endif
typedef lint * plint;         //and a pointer

//classy struct for node properties
struct node {
	//dynamical properties
	bool immune;      //immunity state
	bool state;       //adoption state
	double threshold; //adoption threshold

	//structural properties
	int degree;      //degree in network
	pint neighbours; //pointer to neighbours array

	//auxiliary
	int int_thres;    //int version of threshold, ceil(phi*k)
	int con_thres;    //conditional threshold, {innov, vuln/stab} = {0, int_thres}
	double eff_thres; //efective threshold, determined at adoption
	lint time_adopt;  //time of adoption (microtime units)
	int adopt_neighs; //number of adopting neighbours at adoption
	int pos;          //position in neighbours array

	node(); //constructor
};
node::node() {
	immune = false;	state = false;	threshold = 0.0; //initialise all
	degree = 0;	neighbours = 0;
	int_thres = 0;	con_thres = 0;	eff_thres = 0.0;	time_adopt = 0;	adopt_neighs = 0;	pos = 0;
}
typedef node * pnode; //define type for struct pointer

//classy struct for loops in average degree and threshold
struct loop {
	double ini, fin, inc; //start, end and increment

	loop(); //constructor
};
loop::loop() {
	ini = 0.0;	fin = 0.0;	inc = 0.0; //initialise all
}

//classy struct for output
struct output {
	int num_adopt; //number of adopters (0,...,N),
	int num_innov; //innovators,
	int num_vuln;  //vulnerable,
	int num_stab;  //and stable nodes!

	int neff_innov; //effective innov, vuln and stable nodes!
	int neff_vuln;
	int neff_stab;

	int relax_time; //time to steady state (MC units)

	output(); //constructor
};
output::output() {
	num_adopt = 0;	num_innov = 0;	num_vuln = 0;	num_stab = 0; //initialise all
	neff_innov = 0;	neff_vuln = 0;	neff_stab = 0;
	relax_time = 0;
}
typedef output * poutput; //define type for struct pointer


//some global variables

int N;         //total number of agents
int ini_innov; //number of initial innovators
double pnew;   //rate of innovators
double phi;    //adoption threshold
double z;      //average degree in random network
double r;      //fraction of immune nodes

pnode nodesAll;     //pointer to nodes array
poutput tseriesOut; //pointer to timeseries array

int tmax; //maximum time for dynamics (MC units)
int nt;   //loop counter for realisations

gsl_rng * randGene; //random number generator

//config file variables
pfile conf;
char conff2w[STR_MAX];
char varnames[][VNAME_MAX] = { "N", "ini_innov", "pnew", "r", "tmax", "tstart", "tsnaps", "ntimes", "randSeed", "dataswitch", "clusswitch", "timeswitch", "tendswitch", "distswitch", "coutswitch", "avdeg.ini", "avdeg.fin", "avdeg.inc", "thres.ini", "thres.fin", "thres.inc" };
char varvalues[NUM_VAR][VNAME_MAX];


//random number generator initialiser
void RandInit(long randSeed)
{
	const gsl_rng_type * randType; // type of random generator

	//specify type and seed, if necessary
	randType = gsl_rng_mt19937;
	if ( randSeed == 0 )	randSeed = time( (time_t)NULL );

	randGene = gsl_rng_alloc(randType); //create instance of random generator
	gsl_rng_set(randGene, randSeed); //and seed it
}

//allocate memory for (bool, int, lint, node or output) array arr of size len
void alloc_mem(pbool& arr, int len)
{
	try {	arr = new bool[len];	} //allocate bool array of size len
	catch (exception& e) //catch sneaky exceptions!
	{
		cout << "Standard exception: " << e.what() << endl;
		exit(1);
	}
}
void alloc_mem(pint& arr, int len)
{
	try {	arr = new int[len];	} //allocate int array of size len
	catch (exception& e) //catch sneaky exceptions!
	{
		cout << "Standard exception: " << e.what() << endl;
		exit(1);
	}
}
#if BIG_NETS>0 //only needed for BIG networks...
void alloc_mem(plint& arr, int len)
{
	try {	arr = new lint[len];	} //allocate lint array of size len
	catch (exception& e) //catch sneaky exceptions!
	{
		cout << "Standard exception: " << e.what() << endl;
		exit(1);
	}
}
#endif
void alloc_mem(pnode& arr, int len)
{
	try {	arr = new node[len];	} //allocate node array of size len
	catch (exception& e) //catch sneaky exceptions!
	{
		cout << "Standard exception: " << e.what() << endl;
		exit(1);
	}
}
void alloc_mem(poutput& arr, int len)
{
	try {	arr = new output[len];	} //allocate output array of size len
	catch (exception& e) //catch sneaky exceptions!
	{
		cout << "Standard exception: " << e.what() << endl;
		exit(1);
	}
}

//calculate number of true elements in (bool or int) array with given size
int count(pbool arr, int size)
{
	int acc = 0; //counter
	for (int i = 0; i < size; i++)
		if ( arr[i] ) //sum over the non-zero elements
			acc++;
	return acc;
}
int count(pint arr, int size)
{
	int acc = 0; //counter
	for (int i = 0; i < size; i++)
		if ( arr[i] ) //sum over the non-zero elements
			acc++;
	return acc;
}


//and there come Our Prototypes!

//this function writes both the edge file for a Poisson random network and the threshold file
//according to stage: edges (0) and thresholds (1)
void write_iniFiles(int stage);

//this function reads the edge and threshold files and fills the nodes array
//according to stage: edges (0) and thresholds (1)
void read_iniFiles(int stage);

//this function gets initial innovators randomly
void start_innov();

//this function writes some node properties compactly: degree, int_thres, con_thres, time_adopt and adopt_neighs
void write_data();

//this function searches a network of allowed nodes for clusters and their sizes
void clust_finder( pbool allowed_nodes, pint& cluster_sizes, int maxsize );

//this function gets cluster sizes for different networks in all realisations and stores them in output files
//network clusters go by stage: adoption cascades (0)
void write_clusters(int stage);

//this function stores the output timeseries compactly
void write_tseries();

//this function stores the output effective timeseries compactly
void write_eff_tseries();

//this functions writes the distributions of effective thresholds and waiting times
void write_distributions();

//this function reads the contents of the configuration file
void readCfile();


//there goes thy Main!
int main (int argc, char ** argv)
{
	//define Our Variables!

	//...for the beginning and the end...

	int t; 	    //time (MC units)
	int tstart; //starting time (MC units)
	int tsnaps; //length of time snapshot (MC units)
	int microt; //update counter in a MC time step

	int ntimes; //number of realisations
	loop avdeg, thres; //loop structs for phi and z

	int dataswitch; //the node data is not stored (0) or it is (1)
	int clusswitch; //the cluster finding is off (0) or on (1)
	int timeswitch; //the time series calculation is off (0) or on (1)
	int tendswitch; //the end is defined by max time (0) or steady state (1)
	int distswitch; //the distribution calculation is off (0) or on (1)
	int coutswitch; //the screen output is off (0) or on (1)

	pnode pnodei, pnodej; //node pointers
	int pos; //position

	int nadopt_buff; //buffer for number of adopters

	long randSeed; //seed of random generator

	pfile oFinal;            //file descriptor for final-state output
	char oFinalf2w[STR_MAX];
	output out;		 //current output struct

	//read configuration file name
	if ( argc == 1 )
		sprintf(conff2w, "conf.txt");
	else
		sprintf(conff2w, "%s", argv[1]);

	//fill varvalues array from conf file...
	readCfile();

	//...and set values with it

	//the whole system
	N = atoi(&varvalues[0][0]);
	ini_innov = atoi(&varvalues[1][0]);
	pnew = atof(&varvalues[2][0]);
	r = atof(&varvalues[3][0]);

	//dynamics
	tmax = atoi(&varvalues[4][0]);
	tstart = atoi(&varvalues[5][0]);
	tsnaps = atoi(&varvalues[6][0]);

	//auxiliary
	ntimes = atoi(&varvalues[7][0]);
	randSeed = atoi(&varvalues[8][0]);

	//switchs
	dataswitch = atoi(&varvalues[9][0]);
	clusswitch = atoi(&varvalues[10][0]);
	timeswitch = atoi(&varvalues[11][0]);
	tendswitch = atoi(&varvalues[12][0]);
	distswitch = atoi(&varvalues[13][0]);
	coutswitch = atoi(&varvalues[14][0]);

	//average degree loop
	avdeg.ini = atof(&varvalues[15][0]);
	avdeg.fin = atof(&varvalues[16][0]);
	avdeg.inc = atof(&varvalues[17][0]);

	//adoption threshold loop
	thres.ini = atof(&varvalues[18][0]);
	thres.fin = atof(&varvalues[19][0]);
	thres.inc = atof(&varvalues[20][0]);


	//allocate memory!
	alloc_mem(nodesAll, N); //nodes array
	if ( timeswitch )
		alloc_mem(tseriesOut, tmax+1); //timeseries array

	//initialize random number generator
	RandInit(randSeed);


	//we print a welcome!
	cout << endl << "A threshold model of global cascades on a uniform random network with single threshold" << endl;
	cout << endl << "network size N = " << N << endl;
	cout << endl << "number of initial innovators ini_innov = " << ini_innov << endl;
	cout << endl << "rate of innovators p_new = " << pnew << endl;
	cout << endl << "fraction of immune nodes r = " << r << endl;

	//we live in a dynamical World!

	for (z = avdeg.ini; z <= avdeg.fin; z += avdeg.inc) //average degree loop
	{
		//print parameters
		cout << endl << "average degree z = " << z << endl;

		for (nt = 0; nt < ntimes; nt++) //nt loop, used to obtain average measurements
		{
			//initialise structural/auxiliary variables for all nodes
			for (pos = 0; pos < N; pos++)
			{
				pnodei = &nodesAll[pos]; //get node pointer
				pnodei->degree = 0;	pnodei->neighbours = 0;	pnodei->pos = 0; //initialise all
			}

			//initialise edge file
			write_iniFiles(0); //Poisson random graph (p=z/N)

			//read edge file
			read_iniFiles(0);


			for (phi = thres.ini; phi <= thres.fin; phi += thres.inc) //adoption threshold loop
			{
				//initialise dynamical/auxiliary variables for all nodes
				for (pos = 0; pos < N; pos++)
				{
					pnodei = &nodesAll[pos]; //get node pointer
					pnodei->immune = false;		pnodei->state = false;		pnodei->threshold = 0.0; //initialise all
					pnodei->int_thres = 0;		pnodei->con_thres = 0;		pnodei->eff_thres = 0.0;
					pnodei->adopt_neighs = -1;	pnodei->time_adopt = lint(N)*tmax + 1;
				}

				//initialise theshold file
				write_iniFiles(1); //single threshold phi

				//read threshold file
				read_iniFiles(1);

				//get initial innovators, if any
				if ( ini_innov )
					start_innov();


				//print parameters, in the beginning
				if ( nt == 0 )
					cout << endl << "adoption threshold phi = " << phi << endl;

				//we print the initial state of the system, according to flag
				if ( coutswitch )
				{
					cout << "\tsystem state for realisation nt = " << nt << ", every time snapshot of size " << tsnaps << endl;
					cout << "\t\t{ time  adopt }\t{ innov  vuln  stab }\t{ innov_eff  vuln_eff  stab_eff }" << endl;
					cout << "\t\t{ t  A(t) }\t\t{ I(t)  V(t)  S(t) }\t\t{ I_eff(t)  V_eff(t)  S_eff(t) }" << endl;
					cout << "\t\t{ " << tstart << "  " << ini_innov << " }\t\t{ " << ini_innov << "  0  0 }\t\t{ " << ini_innov << "  0  0 }" << endl; //initial state
				}


				//initialise current output
				out.num_adopt = ini_innov; //initial number of adopters,
				out.num_innov = ini_innov;	out.num_vuln = 0;	out.num_stab = 0; //innov, vuln and stable nodes!
				out.neff_innov = ini_innov;	out.neff_vuln = 0;	out.neff_stab = 0; //and their effective values!
				out.relax_time = tmax; //maximum time
				nadopt_buff = out.num_adopt; //buffer for number of adopters

				//initialise output timeseries to current output
				if ( timeswitch )
					for (pos = 0; pos <= tmax; pos++)
						tseriesOut[pos] = out;

				//initialise output file, just in the beginning
				sprintf(oFinalf2w, "./oFinal-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.txt", N, ini_innov, pnew, r, phi, z);
				if ( nt == 0 )
				{
					if ( ( oFinal = fopen(oFinalf2w, "w") ) == (FILE *) NULL )
					{	perror("\ncant open output file to write\n");		exit(1);	}
					fclose(oFinal);
				}


				//Behold thou fools! For the dynamics art before thee!!

				for (t = tstart+1; t <= tmax; t++) //while there's still some life left in the bones
				{
					for (microt = 0; microt < N; microt++) //N updates per MC time step!
					{
						pnodei = &nodesAll[ gsl_rng_uniform_int(randGene, N) ]; //get random node pointer

						//if node is still susceptible (i.e. no immunity and no adoption)...
						if ( (pnodei->immune == false) && (pnodei->state == false) )
						{
							//first, update condition threshold for adoption
							if ( gsl_rng_uniform_pos(randGene) < pnew ) //if we have an innovator with rate pnew...
								pnodei->con_thres = 0;
							//else we have by default con_thres = int_thres

							//then calculate number of adopting neighbours
							pnodei->adopt_neighs = 0; //initialise count
							for (pos = 0; pos < pnodei->degree; pos++) //loop through neighbours array
							{
								pnodej = &nodesAll[ pnodei->neighbours[pos] ]; //get neighbour node pointer
								if ( pnodej->state ) //if neighbour is an adopter...
									pnodei->adopt_neighs++; //increase counter!
							}

							//if at least a fraction phi out of k neighs are adopters... [ int_thres = ceil( phi*k ) ]
							//(innovators always fulfill this condition)
							if ( pnodei->adopt_neighs >= pnodei->con_thres )
							{
								pnodei->state = true; //adopt the hell out of it!
								pnodei->time_adopt = lint(N)*(t - 1) + microt + 1; //and record microtime

								//if there are too many adopting neighbours, compute effective threshold!
								//(nodes with k=0 never enter here)
								if ( pnodei->adopt_neighs > pnodei->con_thres )
									pnodei->eff_thres = pnodei->adopt_neighs / double( pnodei->degree );

								//finally, update normal counters appropriately...
								out.num_adopt++; //total number of adopters
								if      ( pnodei->con_thres == 0 )	out.num_innov++; //one more innovator,
								else if ( pnodei->con_thres == 1 )	out.num_vuln++; //or vulnerable,
								else					out.num_stab++; //or stable node!

								//and effective counters as well
								if      ( pnodei->adopt_neighs == 0 )	out.neff_innov++; //one more innovator,
								else if ( pnodei->adopt_neighs == 1 )	out.neff_vuln++; //or vulnerable,
								else					out.neff_stab++; //or stable node!
							}
						}
					}  //microt loop


				//after The Dynamics, let us store and show!!

					//show the system's current state, once in a while, according to flag
					if ( (coutswitch == 1) && ( (t % tsnaps) == 0 ) )
						cout << "\t\t{ " << t << "  " << out.num_adopt << " }\t\t{ " << out.num_innov << "  " << out.num_vuln << "  " << out.num_stab << " }\t\t{ " << out.neff_innov << "  " << out.neff_vuln << "  " << out.neff_stab << " }" << endl;

					//store timeseries output
					if ( timeswitch )
						tseriesOut[t] = out; //current output (relax_time is not used)

					//when The End is nigh...
					if ( (t % tsnaps) == 0 ) //every time snapshot we check our creation
					{
						if ( out.num_adopt > nadopt_buff ) //if there's still hope for life, stand up!
							nadopt_buff = out.num_adopt;
						else //we have reach thy Stationary State!
						{
							out.relax_time = t; //update out struct
							if ( tendswitch ) //and if we choose to end with a steady state...
								break;
						}
					}
				} //time loop

				//store node data (int_thres, time_adopt and adopt_neighs)
				if ( dataswitch )	write_data();

				//find adoption cascades
				if ( clusswitch )	write_clusters(0);

				//write timeseries output (normal and effective)
				if ( timeswitch )	{ write_tseries();	write_eff_tseries(); }

				//calculate threshold/waiting-time distributions
				if ( distswitch )	write_distributions();


				//write final-state output
				if( ( oFinal = fopen(oFinalf2w, "a") ) == (FILE *) NULL ) //open file descriptor
				{	perror("\ncant open output file to write\n");		exit(1);	}

				fprintf(oFinal, "%i %i %i %i %i %i %i %i\n", out.relax_time, out.num_adopt, out.num_innov, out.num_vuln, out.num_stab, out.neff_innov, out.neff_vuln, out.neff_stab); //print out struct
				fclose(oFinal); //close file descriptor

			} //adoption threshold loop

			//free memory for neighbours arrays
			for (pos = 0; pos < N; pos++)
				delete [] nodesAll[pos].neighbours;

		} //nt averages loop

	} //average degree loop

	//free memory from arrays
	delete [] nodesAll;
	if ( timeswitch )	delete [] tseriesOut;

	gsl_rng_free(randGene); //free random number generator

	//beware of thy end of main
	return 0;
}


//Behold! For here are all the secondary functions thou needst

//this function writes both the edge file for a Poisson random network and the threshold file
//according to stage: edges (0) and thresholds (1)
void write_iniFiles(int stage)
{
	double p; //edge probability in random network
	pfile edgeF;	char edgeF2W[STR_MAX]; //file descriptor for edges/thresholds
	pfile thresF;	char thresF2W[STR_MAX];
	int i, j; //nodes
	double rand;   //random variable

	if (stage == 0) //first the edge file!
	{
		p = z / double(N); //edge probability

		//write filename and open
		sprintf(edgeF2W, "./edgefile-N%i-z%.2f.txt", N, z);
		if ((edgeF = fopen(edgeF2W, "w")) == (FILE *) NULL ) //open file to write
		{	perror("\ncant open edge file to write\n");	exit(1);	}

		//find random edges!
		for (i = 0; i < N-1; i++)
		{
			for (j = i+1; j < N; j++)
			{
				rand = gsl_rng_uniform_pos(randGene); //get double rand in (0,1)

				if ( rand < p )
					fprintf(edgeF, "%i\t%i\n", i, j); //store edge!
			}
		}
		fclose(edgeF); //close file descriptor! we're A done!
	}
	else if (stage == 1) //then the threshold file!
	{
		//write filename and open
		sprintf(thresF2W, "./thresholds-N%i-phi%.2f.txt", N, phi);
		if ((thresF = fopen(thresF2W, "w")) == (FILE *) NULL ) //open file to write
		{	perror("\ncant open threshold file to write\n");	exit(1);	}

		//go for thresholds!
		for (i = 0; i < N; i++)
			fprintf(thresF, "%.4f\n", phi); //all nodes

		fclose(thresF); //close file descriptor! we're A done!
	}
}

//this function reads the edge and threshold files and fills the nodes array
//according to stage: edges (0) and thresholds (1)
void read_iniFiles(int stage)
{
	pfile edgeF;	char edgeF2W[STR_MAX]; //file descriptor for edges/thresholds
	pfile thresF;	char thresF2W[STR_MAX];
	char line[STR_MAX]; //line in file
	int linenum;        //line number
	int i, j; //nodes
	pnode pnodei, pnodej; //node pointers
	float thres; //threshold
	double min_thres = 0.00001; //min thres to separate zero/nonzero thresholds

	if (stage == 0) //first the edge file!
	{
		//write filename and open
		sprintf(edgeF2W, "./edgefile-N%i-z%.2f.txt", N, z);
		if ((edgeF = fopen(edgeF2W, "r")) == (FILE *) NULL ) //open file to read
		{	perror("\ncant open edge file to read\n");	exit(1);	}

		//get line by line
		linenum = 0;
		while ( fgets(line, STR_MAX, edgeF) != NULL )
		{
			linenum++; //increase line number and skip comments
			if (line[0] == '#') continue;

			//check for correct syntax and store line in integers
			if(sscanf(line, "%i %i", &i, &j) != 2)
			{
				fprintf(stderr, "\nsyntax error in edge file, line %d\n", linenum);
				exit(1);
			}

			nodesAll[i].degree++;	nodesAll[j].degree++; //increase degrees in nodes array
		}
		fclose(edgeF); //close file descriptor! we're A done!

		//allocate memory for neighbours array, for all nodes
		for (i = 0; i < N; i++)
			alloc_mem( nodesAll[i].neighbours , nodesAll[i].degree );

		//open file to read (again)
		if ((edgeF = fopen(edgeF2W, "r")) == (FILE *) NULL )
		{	perror("\ncant open edge file to read\n");	exit(1);	}

		//get line by line
		linenum = 0;
		while ( fgets(line, STR_MAX, edgeF) != NULL )
		{
			linenum++; //increase line number and skip comments
			if (line[0] == '#') continue;

			//check for correct syntax and store line in integers
			if(sscanf(line, "%i %i", &i, &j) != 2)
			{
				fprintf(stderr, "\nsyntax error in edge file, line %d\n", linenum);
				exit(1);
			}

			pnodei = &nodesAll[i];	pnodej = &nodesAll[j]; //get node pointers
			pnodei->neighbours[ pnodei->pos ] = j; //fill up neighbours arrays
			pnodej->neighbours[ pnodej->pos ] = i;
			pnodei->pos++;	pnodej->pos++; //increase positions for next neighbour
		}
		fclose(edgeF); //close file descriptor! we're A done!
		for (i = 0; i < N; i++)	nodesAll[i].pos = 0; //reset position in neighbours array
	}
	else if (stage == 1) //then the threshold file!
	{
		//write filename and open
		sprintf(thresF2W, "./thresholds-N%i-phi%.2f.txt", N, phi);
		if ((thresF = fopen(thresF2W, "r")) == (FILE *) NULL ) //open file to read
		{	perror("\ncant open threshold file to read\n");	exit(1);	}

		//get line by line
		linenum = 0;
		while ( fgets(line, STR_MAX, thresF) != NULL )
		{
			linenum++; //increase line number and skip comments
			if (line[0] == '#') continue;

			//check for correct syntax and store line in double
			if(sscanf( line, "%f", &thres ) != 1)
			{
				fprintf(stderr, "\nsyntax error in threshold file, line %d\n", linenum);
				exit(1);
			}

			nodesAll[linenum-1].threshold = double(thres); //store threshold in nodes array
			nodesAll[linenum-1].eff_thres = double(thres); //and copy to effective threshold
		}
		fclose(thresF); //close file descriptor! we're A done!

		//get immunity state and set/fix versions of adoption thresholds, for all nodes
		for (i = 0; i < N; i++)
		{
			pnodei = &nodesAll[i]; //get node pointer

			if ( i < int(N*r) )
				pnodei->immune = true; //first rN nodes are immune to adoption!

			//first innovators (if phi=0, for example)
			if ( pnodei->threshold < min_thres )
				pnodei->int_thres = 0; //zeros for the bastards!
			else //then everyone else!
			{
				if (pnodei->degree) //for friendly nodes (k > 0)
					//these are vulnerable or stable nodes
					pnodei->int_thres = int( ceil( (pnodei->degree) * (pnodei->threshold) ) ); //up to ceiling int
				else //or for loners...
					pnodei->int_thres = 1; //all vulnerable! if they were connected :(
			}

			//then, initialise condition threshold!
			pnodei->con_thres = pnodei->int_thres;
		}
	}
}

//this function gets initial innovators randomly
void start_innov()
{
	int nInnov; //innovator counter
	double pInnov; //probability of choosing innovator
	int pos; //position
	pnode pnodei; //node pointer

	//initialise all
	nInnov = ini_innov; //all initial innovators left to find
	pInnov = ini_innov / double(N); //get probability
	pos = int(N*r); //and position (ignoring the immune nodes)

	while (nInnov)
	{
		if ( pos == N )	pos = int(N*r); //go back to the start, if necessary
		pnodei = &nodesAll[pos]; //get node pointer

		if ( pnodei->state == false ) //if the state is still off...
			if ( gsl_rng_uniform_pos(randGene) < pInnov ) //and with prob pInnov...
			{
				pnodei->state = true; //switch on initial innovator!
				pnodei->time_adopt = 0; //that adopts at microtime 0!
				pnodei->eff_thres = 0.0; //and has zero effective threshold!

				nInnov--; //one less to the count
			}

		pos++; //go to next position
	}
}

//this function writes some node properties compactly: degree, int_thres, con_thres, time_adopt and adopt_neighs
void write_data()
{
	int pos; //array position
	pnode pnodei; //node pointer

	pfile degreeF;		char degreeF2W[STR_MAX]; //file descriptors
	pfile int_thresF;	char int_thresF2W[STR_MAX];
	pfile con_thresF;	char con_thresF2W[STR_MAX];
	pfile time_adoptF;	char time_adoptF2W[STR_MAX];
	pfile adopt_neighsF;	char adopt_neighsF2W[STR_MAX];
	char degree_mode[STR_MAX]; //modes to open files
	char int_thres_mode[STR_MAX];
	char con_thres_mode[STR_MAX];
	char time_adopt_mode[STR_MAX];
	char adopt_neighs_mode[STR_MAX];

	//prepare files

	sprintf(degreeF2W, "./degree-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(degree_mode, "w"); //first erase previous contents
	else		sprintf(degree_mode, "a"); //and then append each realisation
	if ( ( degreeF = fopen(degreeF2W, degree_mode) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");		exit(1);	}

	sprintf(int_thresF2W, "./int_thres-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(int_thres_mode, "w"); //first erase previous contents
	else		sprintf(int_thres_mode, "a"); //and then append each realisation
	if ( ( int_thresF = fopen(int_thresF2W, int_thres_mode) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");		exit(1);	}

	sprintf(con_thresF2W, "./con_thres-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(con_thres_mode, "w"); //first erase previous contents
	else		sprintf(con_thres_mode, "a"); //and then append each realisation
	if ( ( con_thresF = fopen(con_thresF2W, con_thres_mode) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");		exit(1);	}

	sprintf(time_adoptF2W, "./time_adopt-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(time_adopt_mode, "w"); //first erase previous contents
	else		sprintf(time_adopt_mode, "a"); //and then append each realisation
	if ( ( time_adoptF = fopen(time_adoptF2W, time_adopt_mode) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");		exit(1);	}

	sprintf(adopt_neighsF2W, "./adopt_neighs-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(adopt_neighs_mode, "w"); //first erase previous contents
	else		sprintf(adopt_neighs_mode, "a"); //and then append each realisation
	if ( ( adopt_neighsF = fopen(adopt_neighsF2W, adopt_neighs_mode) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");		exit(1);	}

	//write all data compactly!
	for (pos = 0; pos < N; pos++)
	{
		pnodei = &nodesAll[pos]; //get node pointer

		fwrite( &(pnodei->degree), sizeof(int), 1, degreeF ); //and write in files!
		fwrite( &(pnodei->int_thres), sizeof(int), 1, int_thresF );
		fwrite( &(pnodei->con_thres), sizeof(int), 1, con_thresF );
		fwrite( &(pnodei->time_adopt), sizeof(lint), 1, time_adoptF );
		fwrite( &(pnodei->adopt_neighs), sizeof(int), 1, adopt_neighsF );
	}

	//close file descriptors! we're A done!
	fclose(degreeF);
	fclose(int_thresF);
	fclose(con_thresF);
	fclose(time_adoptF);
	fclose(adopt_neighsF);
}

//this function searches a network of allowed nodes for clusters and their sizes
void clust_finder( pbool allowed_nodes, pint& cluster_sizes, int maxsize )
{
	pint nodeATpos; //pointer to array of nodes at given positions
	pint posOFnode; //pointer to array of positions of given nodes
	pbool ccheck; //pointer to array of checked node positions
	pbool Li, Oi; //position sets for cluster detection
	int posi, posv, posu, posn; //positions
	int i, u, n; //nodes
	int nclu; //number of clusters
	int Xi, Yi; //number of nonzero elements in the Li, Oi arrays
	pnode pnodeu; //node pointer
	int neigh; //position in neighbours array

	//allocate arrays
	alloc_mem( nodeATpos, maxsize ); //nodes at positions
	alloc_mem( posOFnode, N ); //positions of nodes
	alloc_mem( ccheck, maxsize ); //checked node positions
	alloc_mem( Li, maxsize ); //V and O position sets
	alloc_mem( Oi, maxsize );

	//initialise stuff
	posi = 0; //start position of node i
	for (i = 0; i < N; i++)
	{
		if ( allowed_nodes[i] )
		{
			nodeATpos[posi] = i; //store allowed node i at position posi
			posOFnode[i] = posi; //and conversely
			posi++; //and go to next one
		}
		else
			posOFnode[i] = -1; //the rest of the nodes don't have positions
	}
	for (posi = 0; posi < maxsize; posi++)
	{
		cluster_sizes[posi] = 0; //initialise cluster sizes
		ccheck[posi] = true; //at the beginning we can pick any node
	}
	nclu = 0; //we start with no clusters counted

	//finding clusters!
	for (posv = 0; posv < maxsize; posv++) //loop through all node positions
		if ( ccheck[posv] ) //if we haven't picked this node position yet
		{
			for (posi = 0; posi < maxsize; posi++)
			{
				Li[posi] = false; //we initialize everything to false
				Oi[posi] = false;
			}

			Li[posv] = true; //the initial V set
			ccheck[posv] = false; //position posv has been checked

			Oi[posv] = true; //the initial O set

			Xi = count(Li, maxsize); //the initial cardinality of the V,O sets, namely 1
			Yi = count(Oi, maxsize);

			while ( Yi ) //until there is nothing left to be exposed
			{
				posu = gsl_rng_uniform_int(randGene, maxsize); //choose random position
				if ( Oi[posu] ) //check if posu is a live position
				{
					Oi[posu] = false; //we kill posu!

					u = nodeATpos[posu]; //get node u at position posu
					pnodeu = &nodesAll[u]; //get node pointer

					for (neigh = 0; neigh < pnodeu->degree; neigh++) //loop through neighbours array
					{
						n = pnodeu->neighbours[neigh]; //get neighbour of node u

						if ( allowed_nodes[n] ) //if the neighbour is allowed as well...
						{
							posn = posOFnode[n]; //get position of node n

							if (Li[posn] == false) //we mark posn live if it is not in Li
								Oi[posn] = true; //we actualize the O set

							Li[posn] = true; //we actualize the L set
							ccheck[posn] = false; //position posn has been checked
						}
					}

					Xi = count(Li, maxsize); //compute the new cardinality of the V,O sets
					Yi = count(Oi, maxsize);
				}
			}

			cluster_sizes[nclu] = Xi; //the size of the cluster
			nclu++; //one more cluster
		}

	//clean up!
	delete [] nodeATpos;	delete [] posOFnode;
	delete [] ccheck;	delete [] Li;	delete [] Oi;
}

//this function gets cluster sizes for different networks in all realisations and stores them in output files
//network clusters go by stage: adoption cascades (0)
void write_clusters(int stage)
{
	pbool allowed_nodes; //pointer to array of allowed nodes
	pint cluster_sizes; //pointer to array of cluster sizes
	int maxsize; //max cluster size, or max number of clusters
	int i; //node
	pnode pnodei; //node pointer
	pfile fobj;	char fobjF2W[STR_MAX]; //file descriptor
	char mode[STR_MAX]; //mode to open files

	//allocate allowed nodes array
	alloc_mem( allowed_nodes, N );

	if (stage == 0) //first the adoption cascades!
	{
		//get allowed nodes
		for (i = 0; i < N; i++)
		{
			pnodei = &nodesAll[i]; //get node pointer

			if (pnodei->state)	allowed_nodes[i] = true; //all adopters are allowed
			else			allowed_nodes[i] = false; //but not the rest
		}

		//initialise output
		sprintf(fobjF2W, "./sizes_adoptClust-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.txt", N, ini_innov, pnew, r, phi, z);
	}

	//initialise output
	if ( nt == 0 )	sprintf(mode, "w"); //first erase previous contents
	else		sprintf(mode, "a"); //and then append each realisation
	if ( ( fobj = fopen(fobjF2W, mode) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");		exit(1);	}

	//look for clusters!
	maxsize = count(allowed_nodes, N); //max size is number of allowed nodes
	if (maxsize) //if there is anything to do...
	{
		alloc_mem( cluster_sizes, maxsize ); //allocate cluster sizes array
		clust_finder( allowed_nodes, cluster_sizes, maxsize ); //get clusters!

		//write output
		for (i = 0; i < maxsize; i++) //loop through all possible clusters
			if ( cluster_sizes[i] ) //if there is actually one...
				fprintf( fobj, "%i ", cluster_sizes[i] ); //write down cluster size

		delete [] cluster_sizes; //free memory from array
	}
	else //in the case of no clusters...
		fprintf( fobj, "0 " );
	fprintf(fobj, "\n"); //go to next line (i.e. realisation)

	//clean up!
	delete [] allowed_nodes; //free memory from array
	fclose(fobj); //close file descriptor
}

//this function stores the output timeseries compactly
void write_tseries()
{
	sint tiny_innov, tiny_vuln, tiny_stab; //tiny variables for node numbers
	int pos; //array position

	pfile fobj_innov;	char fobjF2W_innov[STR_MAX]; //file descriptors
	pfile fobj_vuln;	char fobjF2W_vuln[STR_MAX];
	pfile fobj_stab;	char fobjF2W_stab[STR_MAX];
	char mode_innov[STR_MAX]; //modes to open files
	char mode_vuln[STR_MAX];
	char mode_stab[STR_MAX];

	//prepare file for innovators
	sprintf(fobjF2W_innov, "./tseries_innov-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(mode_innov, "w"); //first erase previous contents
	else		sprintf(mode_innov, "a"); //and then append each realisation
	if ( ( fobj_innov = fopen(fobjF2W_innov, mode_innov) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");		exit(1);	}

	//then for vulnerable nodes
	sprintf(fobjF2W_vuln, "./tseries_vuln-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(mode_vuln, "w"); //first erase previous contents
	else		sprintf(mode_vuln, "a"); //and then append each realisation
	if ( ( fobj_vuln = fopen(fobjF2W_vuln, mode_vuln) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");		exit(1);	}

	//and finally for stable nodes
	sprintf(fobjF2W_stab, "./tseries_stab-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(mode_stab, "w"); //first erase previous contents
	else		sprintf(mode_stab, "a"); //and then append each realisation
	if ( ( fobj_stab = fopen(fobjF2W_stab, mode_stab) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");		exit(1);	}

	//write output timeseries compactly!
	for (pos = 0; pos <= tmax; pos++)
	{
		tiny_innov = sint( tseriesOut[pos].num_innov ); //transform values to tiny data
		tiny_vuln = sint( tseriesOut[pos].num_vuln );
		tiny_stab = sint( tseriesOut[pos].num_stab );

		fwrite( &(tiny_innov), sizeof(sint), 1, fobj_innov ); //and write in files!
		fwrite( &(tiny_vuln), sizeof(sint), 1, fobj_vuln );
		fwrite( &(tiny_stab), sizeof(sint), 1, fobj_stab );
	}

	//close file descriptors! we're A done!
	fclose(fobj_innov);	fclose(fobj_vuln); fclose(fobj_stab);
}

//this function stores the output effective timeseries compactly
void write_eff_tseries()
{
	sint tiny_innov, tiny_vuln, tiny_stab; //tiny variables for node numbers
	int pos; //array position

	pfile fobj_innov;	char fobjF2W_innov[STR_MAX]; //file descriptors
	pfile fobj_vuln;	char fobjF2W_vuln[STR_MAX];
	pfile fobj_stab;	char fobjF2W_stab[STR_MAX];
	char mode_innov[STR_MAX]; //modes to open files
	char mode_vuln[STR_MAX];
	char mode_stab[STR_MAX];

	//prepare file for innovators
	sprintf(fobjF2W_innov, "./tseries_eff_innov-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(mode_innov, "w"); //first erase previous contents
	else		sprintf(mode_innov, "a"); //and then append each realisation
	if ( ( fobj_innov = fopen(fobjF2W_innov, mode_innov) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");		exit(1);	}

	//then for vulnerable nodes
	sprintf(fobjF2W_vuln, "./tseries_eff_vuln-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(mode_vuln, "w"); //first erase previous contents
	else		sprintf(mode_vuln, "a"); //and then append each realisation
	if ( ( fobj_vuln = fopen(fobjF2W_vuln, mode_vuln) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");		exit(1);	}

	//and finally for stable nodes
	sprintf(fobjF2W_stab, "./tseries_eff_stab-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(mode_stab, "w"); //first erase previous contents
	else		sprintf(mode_stab, "a"); //and then append each realisation
	if ( ( fobj_stab = fopen(fobjF2W_stab, mode_stab) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");		exit(1);	}

	//write output timeseries compactly!
	for (pos = 0; pos <= tmax; pos++)
	{
		tiny_innov = sint( tseriesOut[pos].neff_innov ); //transform values to tiny data
		tiny_vuln = sint( tseriesOut[pos].neff_vuln );
		tiny_stab = sint( tseriesOut[pos].neff_stab );

		fwrite( &(tiny_innov), sizeof(sint), 1, fobj_innov ); //and write in files!
		fwrite( &(tiny_vuln), sizeof(sint), 1, fobj_vuln );
		fwrite( &(tiny_stab), sizeof(sint), 1, fobj_stab );
	}

	//close file descriptors! we're A done!
	fclose(fobj_innov);	fclose(fobj_vuln); fclose(fobj_stab);
}

//this functions writes the distributions of effective thresholds and waiting times
void write_distributions()
{
	int posi, posj; //array positions
	pnode pnodei, pnodej; //node pointers
	plint wTimes, eff_wTimes; //pointer to waiting times arrays
	plint neigh_times; //pointer to array of neighbours' adoption times
	lint wtime, eff_wtime; //effective waiting time for i
	int early_neighs; //number of neighbours adopting before i
	lint deltat; //time difference between i,j's adoptions

	pfile thresF;		char thresF2W[STR_MAX]; //file descriptors
	pfile wtimeF;		char wtimeF2W[STR_MAX];
	pfile eff_wtimeF;	char eff_wtimeF2W[STR_MAX];
	char thres_mode[STR_MAX];     //modes to open files
	char wtime_mode[STR_MAX];
	char eff_wtime_mode[STR_MAX];

//first go with the distribution of effective thresholds

	//prepare file
	sprintf(thresF2W, "./eff_thresholds-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(thres_mode, "w"); //first erase previous contents
	else		sprintf(thres_mode, "a"); //and then append each realisation
	if ( ( thresF = fopen(thresF2W, thres_mode) ) == (FILE *) NULL )
	{	perror("\ncant open threshold file to write\n");	exit(1);	}

	//write effective thresholds compactly!
	for (posi = 0; posi < N; posi++)
	{
		pnodei = &nodesAll[posi]; //get node pointer
		fwrite( &(pnodei->eff_thres), sizeof(double), 1, thresF ); //and write in file!
	}

	fclose(thresF); //close file descriptor! we're A done!


//then deal with the normal/effective waiting time distribution

	//get waiting times

	alloc_mem( wTimes, N ); //allocate arrays of waiting times
	alloc_mem( eff_wTimes, N );

	for (posi = 0; posi < N; posi++) //loop through nodes
	{
		pnodei = &nodesAll[posi]; //get node pointer

		if ( pnodei->state == false ) //if node has not adopted...
		{
			wTimes[posi] = pnodei->time_adopt; //waiting times are maximum -> N*tmax + 1
			eff_wTimes[posi] = pnodei->time_adopt;
		}
		else
		{
			if ( pnodei->degree == 0 ) //for loners the thres condition is satisfied at t=0
			{
				wTimes[posi] = pnodei->time_adopt; //so waiting time is adoption time!
				eff_wTimes[posi] = pnodei->time_adopt;
			}
			else
			{
				alloc_mem( neigh_times, pnodei->degree ); //allocate array of neigh adoption times

				for (posj = 0; posj < pnodei->degree; posj++) //loop through neighbours array
				{
					pnodej = &nodesAll[ pnodei->neighbours[posj] ]; //get neighbour node pointer
					neigh_times[posj] = pnodej->time_adopt; //and its adoption time
				}

				//create vector from neigh times array, and sort the hell out of it!!
				vector<lint> times_vec( neigh_times, neigh_times+pnodei->degree );
				sort( times_vec.begin(), times_vec.end() );

				wtime = pnodei->time_adopt; //initialise all
				eff_wtime = pnodei->time_adopt;
				early_neighs = 0;

				//browse through ordered adoption times
				for ( vector<lint>::iterator iter = times_vec.begin(); iter != times_vec.end(); ++iter )
				{
					deltat = pnodei->time_adopt - *iter; //get diff in adoption times
					early_neighs++; //one more early neighbour

					if ( deltat > 0 ) //if we have an early adopting neighbour...
					{
						eff_wtime = deltat; //update eff waiting time!

						if ( early_neighs == pnodei->int_thres ) //and when int thres is reached...
							wtime = deltat; //update waiting time!
					}
					else
						break; //go away when late adopters arrive
				}

				wTimes[posi] = wtime; //at last, store waiting times!
				eff_wTimes[posi] = eff_wtime;

				delete [] neigh_times; //free array
			}
		}
	}

	//prepare files

	sprintf(wtimeF2W, "./waiting_times-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(wtime_mode, "w"); //first erase previous contents
	else		sprintf(wtime_mode, "a"); //and then append each realisation
	if ( ( wtimeF = fopen(wtimeF2W, wtime_mode) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");	exit(1);	}

	sprintf(eff_wtimeF2W, "./eff_waiting_times-N%i-ini_innov%i-pnew%.5f-r%.2f-phi%.2f-z%.2f.bin", N, ini_innov, pnew, r, phi, z);
	if ( nt == 0 )	sprintf(eff_wtime_mode, "w"); //first erase previous contents
	else		sprintf(eff_wtime_mode, "a"); //and then append each realisation
	if ( ( eff_wtimeF = fopen(eff_wtimeF2W, eff_wtime_mode) ) == (FILE *) NULL )
	{	perror("\ncant open output file to write\n");	exit(1);	}

	//write waiting times compactly!
	for (posi = 0; posi < N; posi++)
	{
		pnodei = &nodesAll[posi]; //get node pointer
		fwrite( &(wTimes[posi]), sizeof(lint), 1, wtimeF ); //and write in file!
		fwrite( &(eff_wTimes[posi]), sizeof(lint), 1, eff_wtimeF );
	}

	delete [] wTimes; //clean up!
	delete [] eff_wTimes;
	fclose(wtimeF); //close file descriptors! we're A done!
	fclose(eff_wtimeF);
}

//this function reads the contents of the configuration file
void readCfile()
{
	char line[STR_MAX]; // line in conf file
	int linenum;        // line number
	bool check;         // check counter
	int cvar;           // checked variable
	int svar[NUM_VAR];  // success array
	int stotal;         // total of succesful reads

	//we open the file descriptor
	if( ( conf = fopen(conff2w, "r") ) == (FILE *) NULL )
	{
		perror("\ncant open file to read\n");
		exit(1);
	}

	//reading the configuration file

	//initialize all
	linenum = 0;
	for (cvar = 0; cvar < NUM_VAR; cvar++)	svar[cvar] = 0;

	//get line by line
	while ( fgets(line, STR_MAX, conf) != NULL )
	{
		char name[VNAME_MAX], value[VNAME_MAX]; //name and value strings

		//increase line number and skip comments
		linenum++;
		if(line[0] == '#') continue;

		//check for correct syntax and store line in strings
		if(sscanf(line, "%s %s", name, value) != 2)
		{
			fprintf(stderr, "\nsyntax error in config file, line %d\n", linenum);
			exit(1);
		}

		//initialize all
		check = true;	cvar = 0;

		//store value in proper place in varvalues array, according to name
		while (check)
		{
			//compare name with each entry in varnames array
			if ( strncmp(varnames[cvar], name, VNAME_MAX) == 0 )
			{
				strcpy( varvalues[cvar], value ); //store value
				svar[cvar] = 1; //add one succesful variable
				check = false; //go out!
			}
			else
			{
				cvar++; //increase cvar
				if ( cvar == NUM_VAR ) //if nothing fits..
				{
					fprintf(stderr, "\nerror in config file, variable %s doesnt exist\n", name);
					exit(1);
				}
			}
		}
	}

	//final check

	stotal = 0; //initialize stotal

	for (cvar = 0; cvar < NUM_VAR; cvar++) //check success reads
	{
		if ( svar[cvar] == 1 )	stotal++;
	}

	if ( stotal == NUM_VAR )
		fclose(conf); //we close the file descriptor
	else
	{
		perror("\nerror in config file, not enough variables to read\n");
		exit(1);
	}
}
