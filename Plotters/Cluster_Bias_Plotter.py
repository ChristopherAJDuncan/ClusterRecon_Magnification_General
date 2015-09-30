import numpy as np
import pylab as pl
import math
import sys

Cluster_Colour = ['b','g','r','c','m','y','k']
STAGES_Cluster_Label = ['A901a', 'A901b', 'A902', 'SW']

######READ in PARAMETER FILE#######

STAGES_Bias_File = sys.argv[1]

#Standard output is "Input, Output, Variance, Input_Ap_Mass, Output_Ap_Mass, Variance" with a space between Clusters.#
STAGES_Bias = np.genfromtxt(STAGES_Bias_File)
nCluster = 4

if( ((1.*STAGES_Bias.shape[0]/nCluster).is_integer() == False) \
    or (1.*STAGES_Bias.shape[0]/nCluster) != 1):
    print 'STAGES input file is not of the expected shape, exiting'
    exit()

f = pl.figure(1, (6,2*6))
ax = f.add_subplot(211)
ax2 = f.add_subplot(212)

#Plot 1-to-1 line#
for i in range(nCluster):
    ax.errorbar(STAGES_Bias[i,0], STAGES_Bias[i,1], yerr = STAGES_Bias[i,2], color = Cluster_Colour[i], label = STAGES_Cluster_Label[i], linewidth = 3., marker = 'x')
    ax2.errorbar(STAGES_Bias[i,0], (STAGES_Bias[i,1]-STAGES_Bias[i,0])/STAGES_Bias[i,0], yerr = STAGES_Bias[i,2]/STAGES_Bias[i,0], color = Cluster_Colour[i], label = STAGES_Cluster_Label[i], linewidth = 3., marker = 'x')    
    
x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
ax.plot(x,x, linestyle = '--')

ax.set_xlabel(r'Input Virial Radius ($r_{200}\; \left[\frac{\rm Mpc}{h}\right]$)'); ax.set_ylabel(r'Output Virial Radius ($r_{200}\; \left[\frac{\rm Mpc}{h}\right]$)')
ax2.set_xlabel(r'Input Virial Radius ($r_{200}\; \left[\frac{\rm Mpc}{h}\right]$)'); ax2.set_ylabel(r'Fractional Variation in Output Radius ($r_{200}\; \left[\frac{\rm Mpc}{h}\right]$)')
ax.legend(loc = 2)
ax2.legend(loc = 3)
f.subplots_adjust(left = 0.18)
pl.show()

if(len(sys.argv) == 2):
    nCluster = 1
    R200_Bias_File = sys.argv[2]

    R200_Bias = np.genfromtxt(R200_Bias_File)

    if( ((1.*R200_Bias.shape[0]/nCluster).is_integer() == False) ):
        print 'R200 input file is not of the expected shape, exiting'
        exit()
