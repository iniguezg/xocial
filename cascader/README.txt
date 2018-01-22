======================================================
cascader - Complex contagion on techno-social networks
======================================================

This program implements a simple model to describe complex contagion on techno-social networks. The description, analysis, and results of this model can be found in the following publications:

M. Karsai, G. Iñiguez, K. Kaski, J. Kertész
Complex contagion process in spreading of online innovation
Journal of the Royal Society Interface 11, 20140694 (2014)
DOI: 10.1098/rsif.2014.0694
arXiv: 1405.6879

Z. Ruan, G. Iñiguez, M. Karsai, J. Kertész
Kinetics of social contagion
Physical Review Letters 115, 218702 (2015)
DOI: 10.1103/PhysRevLett.115.218702
arXiv: 1506.00251

M. Karsai, G. Iñiguez, R. Kikas, K. Kaski, J. Kertész
Local cascades induced global contagion: How heterogeneous thresholds, exogenous effects, and unconcerned behaviour govern online adoption spreading
Scientific Reports 6, 27178 (2016)
DOI: 10.1038/srep27178
arXiv: 1601.07995

=======
SUMMARY
=======

Implementation of a model for global cascades in complex networks: Threshold model in a uniform random graph with single threshold, rate of innovators and immune nodes.

Details:

- Network has size N and average degree z, and all agents have threshold phi
- Innovators are selected randomly either as an initial multiseed (ini_innov), or uniformly in time with rate pnew
- There is a fraction r of immunes nodes
- {vulnerable, stable} nodes have int_thres = { 1, >1 }

Conditions:

- Static multiseed model:  pnew = 0, ini_innov >  0
- Dynamic multiseed model: pnew > 0, ini_innov >= 0
- In general: ini_innov <= N * (1 - r)

================
INITIAL CONTENTS
================

- README.* (this file)
- cascader.cpp (source code)
- compile (compilation executable)
- conf.txt (configuration file)
- cascader.py (Python script for plotting)

======================
INSTRUCTIONS FOR USAGE
======================

COMPILING SOURCE:

- Compile by running executable 'compile':
	$./compile
  to obtain executable 'cascader'.


USING CONFIGURATION FILE:

The configuration file 'conf.txt' stores tunable parameters for the code. A brief description is as follows:

The whole system:

- N - system size
- ini_innov - number of initial innovators
- pnew - rate of innovators
- r - fraction of immune nodes

Dynamics:

- tmax - maximum time of dynamics (Monte Carlo [MC] units)
- tstart - starting time (MC units)
- tsnaps - time snapshot to check final state (not used if tendswitch == 0)

Auxiliary:

- ntimes - number of realisations
- randSeed - seed of random generator (0 to use time since Epoch as seed)

Switchs:

- dataswitch - node data is not stored (0) or it is (1). Output files: degree*, int_thres*, con_thres*, time_adopt*, and adopt_neighs*
- clusswitch - cluster finding is off (0) or on (1). Output files: sizes_adoptClust*
- timeswitch - time series calculation is off (0) or on (1). Output files: tseries_innov*, tseries_vuln*, tseries_stab*, tseries_eff_innov*, tseries_eff_vuln*, and tseries_eff_stab*
- tendswitch - end of dynamics is defined by maximum time 'tmax' (0) or steady state (1)
- distswitch - distribution calculation is off (0) or on (1). Output files: eff_thresholds*, waiting_times*, and eff_waiting_times*
- coutswitch - screen output is off (0) or on (1)

Average degree loop:

- avdeg.ini - initial value
- avdeg.fin - final value
- avdeg.inc - increment

Threshold loop:

- thres.ini - initial value
- thres.fin - final value
- thres.inc - increment


RUNNING SOURCE:

Run executable 'cascader':
	$./cascader


TYPICAL USE:

Minimum screen output (with coutswitch == 0) is:

A threshold model of global cascades on a uniform random network with single threshold

network size N = <N>

number of initial innovators ini_innov = <ini_innov>

rate of innovators p_new = <p_new>

fraction of immune nodes r = <r>

average degree z = <z>

adoption threshold phi = <phi>


Minimum file output (with timeswitch == 1 and other switchs in 0) is:

Main output:

oFinal* - Each row is a realisation of the final state of the dynamics (at t == tmax for tendswitch == 0). Quantities correspond in order to (relax_time, num_adopt, num_innov, num_vuln, num_stab, neff_innov, neff_vuln, neff_stab):

- relax_time - time of final state (MC units)
- num_adopt - number of adopters (0,...,N)
- num_innov, num_vuln, num_stab - number of innovators, vulnerable, and stable adopters
- neff_innov, neff_vuln, neff_stab - number of effective innovators, vulnerable, and stable adopters

Network and thresholds:

- edgefile* - edge file for ER network (node labels in [0, N-1])
- thresholds* - threshold file (row number is node label; for single threshold all values are equal)

Time series:

- tseries_innov*, tseries_vuln*, and tseries_stab*
- tseries_eff_innov*, tseries_eff_vuln*, and tseries_eff_stab*

Binary files where each number corresponds to a MC time step (tstart,...,tmax, for example) and gives the number of adopters (innovators, vulnerable, stable, etc.) as a short integer. Time series for each realisation appear consecutively.


PLOTTING RESULTS:

The Python script 'cascader.py' plots time series of adoption for varying 'pnew'.

Required files: tseries_innov*, tseries_vuln*, and tseries_stab* for all considered parameter values (in the same directory).

Run script with Python (First time use loadflag == n to create numpy array data of adoption time series; later on loadflag == y may be used to slightly improve performance):
  $python cascader.py <loadflag>

Argument and parameter values (equivalent to configuration file in source code):

- loadflag - load numpy array data (y) or not (n)

Simulation parameters:

- z - average degree
- phi - single threshold
- r - fraction of immune nodes

Other parameters:

- N - system size
- ini_innov - number of initial innovators
- tmax - maximum time for dynamics
- ntimes - number of realisations

Arrays:

- pn_vals - array with values of innovator rate 'pnew' (called 'pn' in Python script)

Output files:

- all_tseries_pn* - Numpy array with time series of adoption for all 'pn' values. Array shape: [pn] [innov, vuln, stab] [0, tmax]
- fig_tseries* - Image with plot of time series of adoption for varying 'pn'
