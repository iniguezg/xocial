=========================================================
numenor - Opinion formation on coevolving social networks
=========================================================

This program implements a simple model to describe opinion formation on coevolving social networks. The description, analysis, and results of this model can be found in the following publications:

G. Iñiguez, J. Kertész, K. Kaski, R. A. Barrio
Opinion and community formation in coevolving networks
Physical Review E 80, 066119 (2009)
DOI: 10.1103/PhysRevE.80.066119
arXiv: 0908.1068v2

G. Iñiguez, R. A. Barrio, J. Kertész, K. Kaski
Modelling opinion formation driven communities in social networks
Computer Physics Communications 182, 1866–1869 (2011)
DOI: 10.1016/j.cpc.2010.11.020
arXiv: 1007.4177

G. Iñiguez, J. Kertész, K. Kaski, R. A. Barrio
Phase change in an opinion-dynamics model with separation of time scales
Physical Review E 83, 016111 (2011)
DOI: 10.1103/PhysRevE.83.016111
arXiv: 1009.2643

G. Iñiguez, J. Tagüeña-Martínez, K. Kaski, R. A. Barrio
Are opinions based on science: Modelling social response to scientific facts
PLoS ONE 7, e42122 (2012)
DOI: 10.1371/journal.pone.0042122
arXiv: 1109.1488

================
INITIAL CONTENTS
================

- README.* (this file)
- numenor.cpp (source code)
- compile (compilation executable)
- randomc.h (header for random number generator)
- ranrotb.cpp (source code for random number generator)

========================================
REQUIREMENTS (for network visualisation)
========================================

- himmeli executable (download and documentation available here: http://www.finndiane.fi/software/himmeli/)
Note: Executable should be available globally (i.e. '$himmeli config.txt' works in terminal)

======================
INSTRUCTIONS FOR USAGE
======================

COMPILING SOURCE:

- Compile by running executable 'compile':
  $./compile
  to obtain executable 'numenor'.


PARAMETER VALUES:

Parameter values may be introduced as command line arguments or directly into the code. A brief description is as follows:

Command line arguments:

- N - number of agents in network (system size)
- initp, finp, dp - initial value, final value, and step size for the loop over the generation time 'tgen'

Other parameters (lines 654-691):

Dynamics:

- mcut - maximum layer of neighbours considered for short-interaction term
- h - external field

Rewiring:

- lim - normalization parameter for rewiring probabilities (1 for the original formulation of the model)
- win - maximum opinion magnitude to allow for rewiring (1 for the original formulation of the model)
- y - probability of global rewiring (as opposed to local rewiring)

Initial conditions:

- z - average degree of initial ER network
- xlimit - maximum magnitude for uniform distribution of initial opinions (not used in the original formulation of the model)
- mean, var - mean and variance of the Gaussian distribution of initial opinions

Times:

- tstart - initial time for dynamics (usually 0)
- dt - time step for numerical integration of opinion equation
- ntimes - number of realisations of dynamics

Additional parameters:

- zero - precision parameter used to stop dynamics when opinions are not changing anymore (for finalswitch == 0)
- rewzero - precision parameter used to stop dynamics when edges are not being rewired anymore (for finalswitch == 1)
- cutoff - minimum cutoff value for rewiring probabilities
- twmax - maximum number of generations in dynamics (if timeswitch == 1)
- thres - threshold to percolate network (i.e. cut links) depending on opinion differences (for percswitch == 1)

Switchs:

- mmapswitch - basic arrays ('x', 'alpha', 'A', 'sp') are memory-maped (1) or allocated (0)
- distswitch - cluster and assortativity information are (1) or not (0) calculated
- timeswitch - time series are (1) or not (0) calculated
- initswitch - initial conditions for basic arrays are read from files (1) or set manually (0)
- lifeswitch - dynamics is calculated (1) or not (0)
- percswitch - network is percolated (1) or not (0)
- himmswitch - himmeli graph is drawn (1) or not (0)
- finalswitch - end of dynamics is set by lack of change in rewiring (1) or in opinions (0)


RUNNING SOURCE:

- Run executable 'numenor':
	$./numenor <N> <initp> <finp> <dp>


TYPICAL USE:

Minimum screen output (with lifeswitch == 1 and the rest of switches zero) is:

 < gen rewcheck > (counter of rewirings per generation, every 10 generations)
 < gen rewcheck > ...

 final at | p nt | : < nb ng nw nclu cmax meanz ccoef suscep plen neigh avclu > (properties for final state of dynamics in each realisation, with averages written in output file numenor-*)

 < eul > (mean value of added terms in Euler integration method for opinion equation)

Minimum file output is:

numenor-* -	average properties for final state of dynamics. Quantities correspond in order to (nb ng nw nclu cmax meanz ccoef suscep plen neigh avclu):

- nb, ng, nw - number of black (x == -1), gray (|x| < 1) and white (x == 1) agents
- nclu - number of clusters in network
- cmax - maximum size of cluster in network
- meanz - average degree in network
- ccoef - average clustering coefficient in network
- suscep - susceptibility of network (according to percolation theory)
- plen - average shortest path length in network
- neigh - average number of nth-order neighbours in network (set to n == 2)
- avclu - average cluster size in network

NETWORK VISUALISATION:

Use himmswitch == 1 to visualise initial and final networks with Himmeli.

Additional screen output: Himmeli screen output

Additional file output (files get rewritten for each visualised network):

- config.txt - configuration file for Himmeli
- edges.txt - edge list for network
- vertices.txt - node attributes for network

Himmeli file output:

- graph-*ini*.txt, graph-*fin*.txt - Himmeli configuration, edge and node files for initial/final networks
- graph-*ini_*.svg, graph-*fin_*.svg - visualisation for initial/final networks
- graph-*ini.info.svg, graph-*ini.info.svg - information on node color and shape
