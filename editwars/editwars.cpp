/*	EditWars v1.0
 *	A simple model to describe opinion formation and conflict resolution
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
 * J. Török, G. Iñiguez, T. Yasseri, M. San Miguel, K. Kaski, J. Kertész
 * Opinions, conflicts, and consensus: Modeling social dynamics in a collaborative environment
 * Physical Review Letters 110, 088701 (2013)
 * DOI: 10.1103/PhysRevLett.110.088701
 * arXiv: 1207.4914
 *
 * G. Iñiguez, J. Török, T. Yasseri, K. Kaski, J. Kertész
 * Modeling social dynamics in a collaborative environment
 * EPJ Data Science 3, 7 (2014)
 * DOI: 10.1140/epjds/s13688-014-0007-z
 * arXiv: 1403.3568
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>

using namespace std;

//constants for strings and files
#define STR_MAX 255	// max string/line length
#define NUM_VAR 26	// number of variables
#define VNAME_MAX 20	// max name length

//types for function declarations
typedef double * pdouble;
typedef int * pint;
typedef FILE * pfile;


int N; //total number of agents

double epsT, muT; //tolerance-convergence of opinion dynamics
double epsA, muA; //tolerance-convergence of article editing

double eta; //fraction of noisy agents substituted each time step

pdouble xAll;  //pointer to the opinion array for article (0) and agents (1...N)
pint xNum;     //pointer to the number distribution array for opinions

pdouble xArt; //pointer to the storage array of article values

int tmax;          //maximum time for dynamics
int ntimes;        //number of realisations
int tstart, tlen;  //start and number of time steps for storage
int lagmax;	   //maximum lag to calculate autocorrelation
double xmin, xmax; //min-max value of opinion/article value
int nboxes;        //number of boxes in opinion distribution
double dx;         //size of opinion box

gsl_rng * randGene; //random number generator


//loop structure for tolerance, convergence and noise
struct loop
{
	double ini, fin, inc; //start, end and increment
} tol, conv, noi;

//config file variables
pfile conf;
char conff2w[STR_MAX];
char varnames[][VNAME_MAX] = { "N" , "xmin" , "xmax", "adelta", "epsT", "muT", "tmax", "tstart", "tsnaps", "ntimes", "randSeed", "lagmax", "nboxes", "loopswitch", "timeswitch", "fastswitch", "coutswitch", "tol.ini", "tol.fin", "tol.inc", "conv.ini", "conv.fin", "conv.inc", "noi.ini", "noi.fin", "noi.inc"};
char varvalues[NUM_VAR][VNAME_MAX];


//this micro function computes the box in which xvalue resides, from 0 to (xmax-xmin)/dx - 1 = nboxes - 1
int xbox(double xvalue)
{
	int ibox = int ( floor( (xvalue - xmin)/dx ) );
	return ibox;
}

//this micro function prints the fast-time/final-state output to out at time ti according to flag
void printbox(pfile& out, int ti, int flag)
{
	fprintf(out, "%i %.4f",ti , xAll[0]); //time and article value

	if (flag) //if we want the opinion distribution too...
	{
		for (int box = 0; box < nboxes-1; box++)
			fprintf(out, " %i", xNum[box]); //distribution of opinions
		fprintf(out, " %i\n", xNum[nboxes - 1] ); //and the last box
	}
	else	fprintf(out, "\n"); //go to next line
}

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

//this baby function allocates memory space for n double arrays of size len pointed by arr
void Dalloc(pdouble& arr, int len, int n)
{
	int nelm; //number of elements in the array

	nelm = n*len; //len elements for each one of the n arrays in arr

	arr = (pdouble) calloc(nelm, sizeof(double)); //we allocate memory

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

	arr = (pint) calloc(nelm, sizeof(int)); //we allocate memory

	if (arr == NULL) //if something goes wrong...
	{
		printf ("error allocating memory\n");
		//close all open files and free any allocated memory
		exit (1);
	}
}


//and there were Our Prototypes!

//this function stores the article time series compactly
void Xseries();

//this function computes the statistics for the article time series distribution
void Xstats();

//this function reads the contents of the configuration file
void readCfile();


//there goes thy Main!
int main (int argc, char ** argv)
{
	//define Our Variables!

	//...for the beginning and the end...

	int t; 	    //time (in MC units)
	int tsnaps; //length of time snapshot
	int microt; //update counter in a MC time step

	int nt;       //loop counter for realisations
	double thres; //precision threshold

	int loopswitch; //coupling article-opinions is off (0) or on (1)
	int timeswitch; //the time series calculation is off (0) or on (1)
	int fastswitch; //fast-time output writing is off (0) or on (1)
	int coutswitch; //the screen output is off (0) or on (1)

	double adelta, amin, amax; //half-width and limits for initial article value

	double xmean, xvar; //mean and variance of article/opinions

	int i, j;          //pair of agents
	double xi, xj, xa; //opinions and article

	int box, pos;	   //boxes and positions

	double dynTerm, dynBuff;  //dynamical term for article and opinions, and its buffer

	long randSeed; //seed of random generator
	double rand;   //random variable

	pfile oTime;           //file descriptor for fast-time output
	char oTimef2w[STR_MAX];
	pfile oFinal;           //file descriptor for final-state output
	char oFinalf2w[STR_MAX];


	//read configuration file name
	if ( argc == 1 )
		sprintf(conff2w, "conf.txt");
	else
		sprintf(conff2w, "%s", argv[1]);
//	if ( argc == 1 )
//	{
//		sprintf(conff2w, "conf.txt");
//		epsA = 0.1;			muA = 0.5;
//	}
//	else if ( argc == 2 )
//	{
//		sprintf(conff2w, "%s", argv[1]);
//		epsA = 0.1;			muA = 0.5;
//	}
//	else if ( argc == 3 )
//	{
//		sprintf(conff2w, "%s", argv[1]);
//		epsA = atof(&argv[2][0]);	muA = 0.5;
//	}
//	else
//	{
//		sprintf(conff2w, "%s", argv[1]);
//		epsA = atof(&argv[2][0]);	muA = atof(&argv[3][0]);
//	}

	//fill varvalues array from conf file...
	readCfile();

	//...and set values with it

	//the whole system
	N = atoi(&varvalues[0][0]);
	xmin = atof(&varvalues[1][0]);
	xmax = atof(&varvalues[2][0]);
	adelta = atof(&varvalues[3][0]);

	//discussion and editing
	epsT = atof(&varvalues[4][0]);
	muT = atof(&varvalues[5][0]);

	//dynamics
	tmax = atoi(&varvalues[6][0]);
	tstart = atoi(&varvalues[7][0]);
	tsnaps = atoi(&varvalues[8][0]);

	//auxiliary
	ntimes = atoi(&varvalues[9][0]);
	randSeed = atoi(&varvalues[10][0]);
	lagmax = atoi(&varvalues[11][0]);
	nboxes = atoi(&varvalues[12][0]);

	//switchs
	loopswitch = atoi(&varvalues[13][0]);
	timeswitch = atoi(&varvalues[14][0]);
	fastswitch = atoi(&varvalues[15][0]);
	coutswitch = atoi(&varvalues[16][0]);

	//article tolerance loop
	tol.ini = atof(&varvalues[17][0]);
	tol.fin = atof(&varvalues[18][0]);
	tol.inc = atof(&varvalues[19][0]);

	//article convergence loop
	conv.ini = atof(&varvalues[20][0]);
	conv.fin = atof(&varvalues[21][0]);
	conv.inc = atof(&varvalues[22][0]);

	//noise rate loop
	noi.ini = atof(&varvalues[23][0]);
	noi.fin = atof(&varvalues[24][0]);
	noi.inc = atof(&varvalues[25][0]);

	//what remains out of the conf file...
	thres = 0.0001;
	tlen = tmax - tstart + 1;
	dx = (xmax - xmin)/(double)nboxes;
	amin = (xmin + xmax)/2.0 - adelta;	amax = (xmin + xmax)/2.0 + adelta;

	//allocate memory for arrays
	Dalloc(xAll, N+1, 1);
	Ialloc(xNum, nboxes, 1);

	if ( timeswitch )
		Dalloc(xArt, ntimes*tlen, 1);


	//we live in a dynamical World!

	for (epsA = tol.ini; epsA <= tol.fin; epsA += tol.inc) //article tolerance loop
	{
		for (muA = conv.ini; muA <= conv.fin; muA += conv.inc) //article convergence loop
		{
			for (eta = noi.ini; eta <= noi.fin; eta += noi.inc) //noise rate loop
			{
				//initialize random number generator
				RandInit(randSeed);

				//initialise article storage array, if needed
				if ( timeswitch )
				{
					for (pos = 0; pos < ntimes*tlen; pos++)
						xArt[pos] = 0.0;
				}

				//we initialise the final-state output according to switch
				if ( timeswitch == 0 )
				{
					sprintf(oFinalf2w, "./oFinal-N%i-epsT%.6f-muT%.2f-epsA%.6f-muA%.2f-eta%.2f.txt", N, epsT, muT, epsA, muA, eta);
					if ( ( oFinal = fopen(oFinalf2w, "w") ) == (FILE *) NULL )
					{
						perror("\ncant open file to write\n");		exit(1);
					}
					fclose(oFinal);
				}

				//we print relevant parameters
				cout << endl << " initial at N = " << N << " and { epsT  muT  epsA  muA  eta } = { " << epsT << "  " << muT << "  " << epsA << "  " << muA << "  " << eta << " }" << endl;

				for (nt = 0; nt < ntimes; nt++) //nt loop, used to obtain average measurements
				{
					//initialize arrays
					for (box = 0; box < nboxes; box++)
						xNum[box] = 0;

					//initial conditions for the article/opinion array

					xAll[0] = amin + (amax - amin)*gsl_rng_uniform_pos(randGene); //first the article!

					for (i = 1; i <= N; i++) //then the agents
					{
						rand = gsl_rng_uniform_pos(randGene); //get double rand in (0,1)
						xi = xmin + (xmax - xmin)*rand;   //renormalise and get opinions

						xAll[i] = xi;	xNum[ xbox(xi) ]++; //store in arrays
					}

					//for time and life... initialise other stuff as well
					dynBuff = 0.0; //dynamics buffer
					if ( ( timeswitch == 1 ) && ( tstart == 0 ) )
						xArt[nt*tlen] = xAll[0]; //store initial article in xArt, if needed


					//we initialise the fast-time output according to switch
					if ( (fastswitch == 1) && (nt == 0) )
					{
						sprintf(oTimef2w, "./oTime-N%i-epsT%.6f-muT%.2f-epsA%.6f-muA%.2f-eta%.2f.txt", N, epsT, muT, epsA, muA, eta);
						if ( ( oTime = fopen(oTimef2w, "w") ) == (FILE *) NULL )
						{
							perror("\ncant open file to write\n");		exit(1);
						}

						printbox(oTime, 0, 1); //output initial state (with flag)
					}

					//we print the initial state of the system, according to flag
					if ( coutswitch )
					{
						cout << "\tinitial article value a(0) = " << xAll[0] << endl;

						xmean = gsl_stats_mean(xAll, 1, N+1); //get opinion statistics
						xvar = gsl_stats_variance(xAll, 1, N+1);
						cout << "\t\tinitial opinion statistics { E[x(0)]  Var[x(0)] } = { " << xmean << "  " << xvar << " }" << endl;

						cout << "\tarticle evolution for nt = " << nt << ", every snapshot of size " << tsnaps << endl;
						cout << "\t\t{ t  a(t) }" << endl;
					}


					//Behold thou fools! For the dynamics art before thee!!

					for (t = 1; t <= tmax; t++) //while there's still some life left in the bones
					{
					//let us carry on with the discussing, editing and renewing business...
						for (microt = 0; microt < N; microt++) //N updates per MC time step!
						{
						//first we discuss in the talk page!

							i = gsl_rng_uniform_int(randGene, N) + 1; //Life is fun in pairs!
							j = gsl_rng_uniform_int(randGene, N) + 1;

							xi = xAll[i];	xj = xAll[j]; //get opinions
							dynTerm = muT*( xj - xi ); //and fill dynamical term

							if ( fabs(dynTerm) < (muT*epsT) ) //if agents have similar opinions... uuhh
							{
								xAll[i] += dynTerm;	xAll[j] -= dynTerm; //opinion changing

								xNum[ xbox(xi) ]--;		xNum[ xbox(xj) ]--;
								xNum[ xbox(xi + dynTerm) ]++;	xNum[ xbox(xj - dynTerm) ]++; //update distribution

								dynBuff += fabs( dynTerm ); //increase buffer
							}

						//then we edit the hell out of the article!!

							i = gsl_rng_uniform_int(randGene, N) + 1; //random editor!

							xi = xAll[i];	xa = xAll[0]; //get opinions
							dynTerm = muA*( xi - xa ); //fill dynamical term

							if ( fabs(dynTerm) > (muA*epsA) ) //when the editor disagrees with the article... double uuhh
							{
								xAll[0] += dynTerm;  //article changing

								dynBuff += fabs( dynTerm ); //increase buffer
							}
							else if ( loopswitch == 1 ) //when the editor agrees with the article... triple uuhh!
							{
								xAll[i] -= dynTerm; //opinion changing

								xNum[ xbox(xi) ]--;	xNum[ xbox(xi - dynTerm) ]++; //update distribution

								dynBuff += fabs( dynTerm ); //increase buffer
							}

						//give us some new agents!

							rand = gsl_rng_uniform_pos(randGene); //get double rand in (0,1)

							if ( rand < eta ) //if there is trouble coming
							{
								i = gsl_rng_uniform_int(randGene, N) + 1; //choose a random troublemaker

								xi = xAll[i];		xNum[ xbox(xi) ]--; //previous opinion
								xi = xmin + (xmax - xmin)*gsl_rng_uniform_pos(randGene);
								xAll[i] = xi;		xNum[ xbox(xi) ]++; //new random opinion!
							}
						} //microt loop

					//after The Dynamics, let us store and show!!

						//evolve xArt! if needed
						if ( ( timeswitch == 1 ) && ( t >= tstart ) ) //after tstart, of course
							xArt[nt*tlen + t - tstart] = xAll[0];

						//we write the fast-time output according to switch
						if ( (fastswitch == 1) && (nt == 0) )
							printbox(oTime, t, 1); //output intermediate state (with flag)

						//show the system's current state, once in a while, according to flag
						if ( (coutswitch == 1) && ( (t % tsnaps) == 0 ) )
							cout << "\t\t" << t << "  " << xAll[0] << endl;

						//when The End is nigh...
						if ( timeswitch == 0 ) //we ain't need no time series!
						{
							if ( (t % tsnaps) == 0 ) //every time snapshot we check our creation
							{
								if ( dynBuff > thres ) //if there's still hope for life, stand up!
									dynBuff = 0.0;
								else //we have reach thy Stationary State!
									break;
							}
						}
					} //time loop

					//we close the fast-time output according to switch
					if ( (fastswitch == 1) && (nt == 0) )
						fclose(oTime); //close file descriptor

					//we write the final-state output according to switch
					if ( timeswitch == 0 )
					{
						if( ( oFinal = fopen(oFinalf2w, "a") ) == (FILE *) NULL ) //open file descriptor
						{
							perror("\ncant open file to write\n");		exit(1);
						}

						printbox(oFinal, t, 0); //output final state (no flag)

						fclose(oFinal); //close file descriptor
					}

					//we print the final state of the system, according to flag
					if ( coutswitch )
					{
						cout << "\tfinal article value a(t) = " << xAll[0] << endl;

						xmean = gsl_stats_mean(xAll, 1, N+1); //get opinion statistics
						xvar = gsl_stats_variance(xAll, 1, N+1);
						cout << "\t\tfinal opinion statistics { E[x(t)]  Var[x(t)] } = { " << xmean << "  " << xvar << " }" << endl;
					}
//					cout << endl << " final for nt = " << nt << endl;
				} //nt averages loop

				//we write the time series according to switch
				if ( timeswitch )
					Xseries();

				//free random number generator
				gsl_rng_free(randGene);
			} //noise rate loop
		} //article convergence loop
	} //article tolerance loop


	//beware of thy end of main
	return 0;
}


//Behold! For here are all the secondary functions thou needst

//this function stores the article time series compactly
void Xseries()
{
	typedef unsigned short int ushort; //tiny variable definition
	ushort t_art; //tiny version of article value

	pfile artF;	char artF2W[STR_MAX]; //file descriptor for the time series output
	int pos; //array position
	double prec = 10000.0; //precision (e.g. 10^-4) to store tiny data

	//write file name and open
	sprintf(artF2W, "./art-N%i-epsT%.6f-muT%.2f-epsA%.6f-muA%.2f-eta%.2f.bin", N, epsT, muT, epsA, muA, eta);
	if((artF = fopen(artF2W, "w")) == (FILE *) NULL ) //open file to write
	{	perror("\ncant open file to write\n");	exit(1);	}

	//write article time series compactly!
	for (pos = 0; pos < ntimes*tlen; pos++)
	{
		t_art = (ushort) floor( xArt[pos]*prec ); //transform article value to tiny data
		fwrite( &(t_art), sizeof(ushort), 1, artF ); //and write in file!
	}
	fclose(artF); //close file descriptor! we're A done!
}


//this function computes the statistics for the article time series distribution
void Xstats()
{
	//absolute change, average value and autocorrelation
	short t_chan, t_aval, t_acor; //tiny version
	double d_chan, d_aval, d_acor; //double version

	pfile chanF;	char chanF2W[STR_MAX]; //file descriptors for the time series output
	pfile avalF;	char avalF2W[STR_MAX];
	pfile acorF;	char acorF2W[STR_MAX];

	int t, lag; //time and lag
	int nt;     //particular realisation
	pdouble mean; //pointer to the mean array
	double yOne, yTwo; //no-mean article values for autocorrelation
	double dntimes = (double)ntimes; //double version of ntimes
	double prec = 10000.0; //precision (e.g. 10^-4) to store tiny data

	//write file names
	sprintf(chanF2W, "./chan-N%i-epsT%.6f-muT%.2f-epsA%.6f-muA%.2f-eta%.2f.bin", N, epsT, muT, epsA, muA, eta);
	sprintf(avalF2W, "./aval-N%i-epsT%.6f-muT%.2f-epsA%.6f-muA%.2f-eta%.2f.bin", N, epsT, muT, epsA, muA, eta);
	sprintf(acorF2W, "./acor-N%i-epsT%.6f-muT%.2f-epsA%.6f-muA%.2f-eta%.2f.bin", N, epsT, muT, epsA, muA, eta);

	Dalloc(mean, tlen, 1); //allocate mean array


	//first the absolute change!

	if((chanF = fopen(chanF2W, "w")) == (FILE *) NULL ) //open file to write
	{	perror("\ncant open file to write\n");	exit(1);	}
	for (t = tstart; t <= tmax; t++) //go on with the calculating! one time at a time
	{
		d_chan = 0.0; //initialise
		if (t > tstart)
		{
			for (nt = 0; nt < ntimes; nt++)
				d_chan += fabs( xArt[nt*tlen + t - tstart] - xArt[nt*tlen + t-1 - tstart] ); //accumulate
			d_chan /= dntimes; //and normalise!
		}
		t_chan = (short) floor( d_chan*prec ); //transform to tiny data
		fwrite( &(t_chan), sizeof(short), 1, chanF ); //write in file!
	}
	fclose(chanF); //close file descriptor! we're A done!


	//then the average value!

	if((avalF = fopen(avalF2W, "w")) == (FILE *) NULL ) //open file to write
	{	perror("\ncant open file to write\n");	exit(1);	}
	for (t = tstart; t <= tmax; t++) //go on with the calculating! one time at a time
	{
		d_aval = 0.0; //initialise!
		for (nt = 0; nt < ntimes; nt++)
			d_aval += xArt[nt*tlen + t - tstart]; //accumulate
		d_aval /= dntimes; //and normalise!

		mean[t - tstart] = d_aval; //store in mean array
		t_aval = (short) floor( d_aval*prec ); //transform to tiny data
		fwrite( &(t_aval), sizeof(short), 1, avalF ); //write in file!
	}
	fclose(avalF); //close file descriptor! we're A done!


	//finally the autocorrelation!

	if((acorF = fopen(acorF2W, "w")) == (FILE *) NULL ) //open file to write
	{	perror("\ncant open file to write\n");	exit(1);	}
	for (lag = 0; lag <= lagmax; lag++) //for each lag... until lagmax
	{
		for (t = tstart; t <= (tmax-lag); t++) //and one time at a time...
		{
			d_acor = 0.0; //initialise
			for (nt = 0; nt < ntimes; nt++)
			{
				yOne = xArt[nt*tlen + t - tstart] - mean[t - tstart]; //substract the meany means!
				yTwo = xArt[nt*tlen + t+lag - tstart] - mean[t+lag - tstart];
				d_acor += yOne*yTwo; //accumulate
			}
			d_acor /= dntimes; //and normalise!

			t_acor = (short) floor( d_acor*prec ); //transform to tiny data
			fwrite( &(t_acor), sizeof(short), 1, acorF ); //write in file!
		}
	}
	fclose(acorF); //close file descriptor! we're A done!
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
