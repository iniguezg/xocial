# editwars - Opinion formation and conflict resolution

This program implements a simple model to describe opinion formation and conflict resolution. The description, analysis, and results of this model can be found in the following publications:

J. Török, G. Iñiguez, T. Yasseri, M. San Miguel, K. Kaski, J. Kertész  
*Opinions, conflicts, and consensus: Modeling social dynamics in a collaborative environment*  
Physical Review Letters **110**, 088701 (2013)  
[DOI: 10.1103/PhysRevLett.110.088701](http://link.aps.org/doi/10.1103/PhysRevLett.110.088701)  
[arXiv: 1207.4914](https://arxiv.org/abs/1207.4914)

G. Iñiguez, J. Török, T. Yasseri, K. Kaski, J. Kertész  
*Modeling social dynamics in a collaborative environment*  
EPJ Data Science **3**, 7 (2014)  
[DOI: 10.1140/epjds/s13688-014-0007-z](http://dx.doi.org/10.1140/epjds/s13688-014-0007-z)  
[arXiv: 1403.3568](https://arxiv.org/abs/1403.3568)

## INITIAL CONTENTS

- `README.*` (this file)
- `editwars.cpp` (source code)
- `compile` (compilation executable)
- `conf.txt` (configuration file)

## REQUIREMENTS (in i.e. Ubuntu)

- `libgsl2` - GNU Scientific Library (GSL) - library package
- `libgsl-dev` - GNU Scientific Library (GSL) - development package

## INSTRUCTIONS FOR USAGE

COMPILING SOURCE:

- Compile by running executable `compile`:  
	`$./compile`  
to obtain executable `editwars`.

USING CONFIGURATION FILE:

The configuration file `conf.txt` stores tunable parameters for the code. A brief description is as follows:

The whole system:

- `N` - number of agents
- `xmin` - min value of opinion/article value
- `xmax` - max value of opinion/article value
- `adelta` - half-width for initial article value (uniformly chosen at random)

Discussion and editing:

- `epsT` - tolerance of opinion dynamics
- `muT` - convergence of opinion dynamics

Dynamics:

- `tmax` - maximum time of dynamics (Monte Carlo [MC] units)
- `tstart` - starting time (MC units)
- `tsnaps` - time snapshot to check final state (not used if `timeswitch == 1`)

Auxiliary:

- `ntimes` - number of realisations
- `randSeed` - seed of random generator (0 to use time since Epoch as seed)
- `lagmax` - maximum lag to calculate autocorrelation in article time series (not used currently)
- `nboxes` - number of intervals to store opinions as discrete distribution

Switchs:

- `loopswitch` - coupling article-opinions is off (0) or on (1). Leave always on for standard version of model
- `timeswitch` - the article time series storage is off (0) or on (1). Saves all realisations in binary file `art*`
- `fastswitch` - the time / article / opinion discrete distribution storage is off (0) or on (1). Saves first realisation only in file `oTime*`
- `coutswitch` - detailed screen output is off (0) or on (1)

Article tolerance loop:

- `tol.ini` - initial value
- `tol.fin` - final value
- `tol.inc` - increment

Article convergence loop:

- `conv.ini` - initial value
- `conv.fin` - final value
- `conv.inc` - increment

Noise rate loop:

- `noi.ini` - initial value
- `noi.fin` - final value
- `noi.inc` - increment

RUNNING SOURCE:

- Run executable `editwars`:  
  `$./editwars`

TYPICAL USE:

Screen output (with `coutswitch == 1`) is:
```
initial at N = <N> and { epsT  muT  epsA  muA  eta } = { <epsT>  <muT>  <epsA>  <muA>  <eta> }
initial article value a(0) = <a(0)>
initial opinion statistics { E[x(0)]  Var[x(0)] } = { <E[x(0)]>  <Var[x(0)]> }
article evolution for nt = <nt>, every snapshot of size <tsnaps>
{ t  a(t) }
<t>  <a(t)>
...

final article value a(t) = <a(t)>
final opinion statistics { E[x(t)]  Var[x(t)] } = { <E[x(t)]>  <Var[x(t)]> }
```

File output is:

Final state output (if `timeswitch == 0`):

- `oFinal*` - Each row is a realisation of the final state of the dynamics. Quantities correspond in order to: final time, final article value

Article time series output (if `timeswitch == 1`):

- `art*` - Article values, stored in time order, one realisation after another, in binary format (to save space)

Time / article / opinion discrete distribution output (if `fastswitch == 1`):

- `oTime*` - Time, and corresponding values of article and opinion distribution (i.e. number of agents with certain opinion, accumulated over discrete intervals), for the first realisation only
