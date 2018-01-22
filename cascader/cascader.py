#! /usr/bin/env python

#	Cascader v1.0
#	A simple model to describe complex contagion on techno-social networks
#	Visualisation script
#	Copyright (C) 2018 Gerardo Iniguez
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# If you use this program for your research, please cite any relevant from the following articles:
#
# M. Karsai, G. Iniguez, K. Kaski, J. Kertesz
# Complex contagion process in spreading of online innovation
# Journal of the Royal Society Interface 11, 20140694 (2014)
# DOI: 10.1098/rsif.2014.0694
# arXiv: 1405.6879
#
# Z. Ruan, G. Iniguez, M. Karsai, J. Kertesz
# Kinetics of social contagion
# Physical Review Letters 115, 218702 (2015)
# DOI: 10.1103/PhysRevLett.115.218702
# arXiv: 1506.00251
#
# M. Karsai, G. Iniguez, R. Kikas, K. Kaski, J. Kertesz
# Local cascades induced global contagion: How heterogeneous thresholds, exogenous effects, and unconcerned behaviour govern online adoption spreading
# Scientific Reports 6, 27178 (2016)
# DOI: 10.1038/srep27178
# arXiv: 1601.07995

#import modules
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from brewer2mpl import qualitative


#function to get mean time series from numerical simulations
def meanTseries(len_time, ntimes, N, fnames):
	"""Mean time series from numerical simulations"""

	mean_tseries = np.zeros(( 3, len_time )) #initialise mean timeseries: [innov, vuln, stab] [0, tmax]

	for type in range(3): #loop through node types
		allruns = np.fromfile( fnames[type], dtype=np.short ) #load all runs

		for pos in range( ntimes ): #loop through realisations
			mean_tseries[type, :] += allruns[ pos*len_time : (pos+1)*len_time ] #accumulate

	mean_tseries /= N*ntimes #get mean timeseries!

	return mean_tseries


#small function to format strings with numbers in exp notation
def eformat(f, prec, exp_digits):
	"""Format strings with numbers in exp notation"""
	s = "%.*e"%(prec, f)
	mantissa, exp = s.split('e')
	# add 1 to digits as 1 is taken by sign +/-
	return "%se%+0*d"%(mantissa, exp_digits+1, int(exp))


### SETUP ###

loadflag = sys.argv[1] #load numpy array data (y) or not (n)

#simulation parameters
z = 7.0   #average degree
phi = 0.2 #single threshold
r = 0.0   #fraction of immune nodes

#other parameters
N = 10000 #system size
ini_innov = 0 #number of initial innovators
tmax = 10000 #maximum time for dynamics
ntimes = 100 #number of realisations

#arrays
pn_vals = np.array([ 2e-5, 2e-4, 1e-3, 5e-3 ]) #values of innovator rate
time = np.linspace(0, tmax, tmax+1) #time array [0, tmax] (integer t)


### DATA ###

fname_pn = '_z{0:.2f}_phi{1:.2f}_r{2:.2f}'.format(z, phi, r) #filename

if loadflag == 'y': #load numpy array data

	#data for varying pn (and fixed r)
	all_tseries_pn = np.load( open( 'all_tseries_pn'+fname_pn+'.npy', 'r' ) )

else: #or compute it!

	#data for varying pn (and fixed r)

	all_tseries_pn = np.zeros(( len(pn_vals), 3, len(time) )) #initialise all timeseries: [pn] [innov, vuln, stab] [0, tmax]
	for pospn, pn in enumerate(pn_vals): #loop through pn values

		Pstring = '-N{0}-ini_innov{1}-pnew{2:.5f}-r{3:.2f}-phi{4:.2f}-z{5:.2f}.bin'.format(N, ini_innov, pn, r, phi, z) #parameter string
		fnames = ( 'tseries_innov'+Pstring, 'tseries_vuln'+Pstring, 'tseries_stab'+Pstring ) #filenames to load

		mean_tseries = meanTseries( tmax+1, ntimes, N, fnames ) #compute mean timeseries [innov, vuln, stab]
		all_tseries_pn[pospn, :, :] = mean_tseries #store num timeseries

	np.save( 'all_tseries_pn'+fname_pn, all_tseries_pn ) #save numpy array data


## 1 ## Figure: Numerical simulation of cascade model

#plot variables
fig_num = 1 #figure properties
fig_size = (10, 7)
aspect_ratio = (1, 1) #grid properties
grid_params = dict( left=0.075, bottom=0.09, right=0.975, top=0.98 )
savename = 'fig_tseries_z{0:.2f}_phi{1:.2f}_r{2:.2f}'.format(z, phi, r) #filename to save
bmap = qualitative.Paired[ len(pn_vals) ] #colormap (note: there is a hard limit for number of pn values; check docs)

#initialise plot
plt.figure( fig_num, figsize=fig_size )
plt.clf()
grid = gridspec.GridSpec( *aspect_ratio )
grid.update(**grid_params)

#sizes/widths/coords
xylabel = 22
ticklabel = 18
textsize = 22
linewidth = 3
tickwidth = 1.5
legend_prop = { 'size':18 }
legend_hlen = 1
legend_cspa = 0.8
frame_ec = '0.7'
frame_lw = 2

## A ## Plot timeseries of rho for varying innovator rate pn

#subplot variables
bmap_pos = 0 #position in bmap

#initialise subplot
axA = plt.subplot( grid[0, 0] )
plt.xlabel(r'$t$', size=xylabel, labelpad=2)
plt.ylabel(r'$\rho$', size=xylabel, labelpad=5)

#subplot plot
for pospn, pn in enumerate(pn_vals): #loop through pn values

	#data plot
	xplot = time #get xplot
	yplot = np.sum( all_tseries_pn[pospn, :, :], axis=0 ) #sum innov, vuln and stab: [0, tmax]
	plt.plot( xplot, yplot, lw=linewidth, ls='-', c=bmap.mpl_colors[bmap_pos], label=eformat(pn, 0, 1) ) #subplot plot!

	bmap_pos += 1 #go to next color

#ticks/labels
axA.tick_params('both', width=tickwidth, which='major')
axA.tick_params('both', width=tickwidth, which='minor')
plt.text( 0.04, 0.95, '$r=$ {:.0f}'.format(r), fontsize=textsize, va='top', ha='left', transform=axA.transAxes )

#finalise subplot
axA.set_xscale('log')
plt.axis([ 1e0, 1e4, -0.05, 1.05 ])
plt.xticks( np.logspace( 0, 4, 5 ), size=ticklabel )
plt.yticks( np.linspace( 0, 1, 6 ), size=ticklabel )

#legend
leg = plt.legend( loc='lower right', bbox_to_anchor=(0.98, 0.03), title='$p=$', prop=legend_prop, handlelength=legend_hlen, fancybox=True )
leg.get_title().set_fontsize( str(textsize) )
frame = leg.get_frame()
frame.set_edgecolor(frame_ec)
frame.set_linewidth(frame_lw)


#finalise plot
#plt.show()

#save figure
if savename != '':
	plt.savefig(savename+'.pdf', format='pdf')
	#plt.savefig(savename+'.jpg', format='jpg')
	#plt.savefig(savename+'.png', format='png')
	#plt.savefig(savename+'.eps', format='eps')
