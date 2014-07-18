import os
import pylab as pl
import re
import numpy as np

### MATPLOTLIB RC PARAMS
params = {'legend.fontsize' : 20,
	'legend.linewidth': 2,
	'axes.linewidth'  : 3.5,
	'axes.labelsize'  : 30,
	'xtick.major.width' : 3,
	'xtick.major.pad'   : 20,
	'ytick.major.width' : 3,
	'xtick.labelsize'    : 30,
	'ytick.labelsize'    : 30,
	'text.usetex'      : True }
#          'mathtext.default'  : 'serif:bold'}
pl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
pl.rcParams.update(params)

# load the data from the file
X = np.loadtxt("gse_ifo_b.dat")

# plot the data
fig = pl.figure()
fig.subplots_adjust(left=0.16,bottom=0.16)
ax = fig.add_subplot(111)
for i in range(1,len(X[0])):
	pl.plot(X[:,0],X[:,i],lw=3)

ax.set_xlabel(r"$b$")
ax.set_ylabel(r"$E_{\textrm{g.s}}$")

pl.savefig("gse_ifo_b.pdf")
pl.show()
