#Produces ML-point-estimate plots of posteriors for 4 experiments, for two masses. Uses the same output as the Non-Biased, Single Cluster Plot

import numpy as np
import pylab as pl
import math
import sys

Experiment_Colour = ['r','b']
#STAGES_Cluster_Label = ['A901a', 'A901b', 'A902', 'SW']

do_Logscale = False

######READ in PARAMETER FILE#######
nMass = 2
Bias_Directory_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/'
Experiment_Directory = ['SizeMag/SingleA901a_Bias/COMBO/Bias/', 'SizeOnly/SingleA901a_Bias/COMBO/Bias/', 'SizeMag/SingleA901a_Bias/STAGES+COMBO/Bias/', 'SizeOnly/SingleA901a_Bias/STAGES+COMBO/Bias/']
Mass_Directory = ['r200_0.8/', 'r200_1.2/']

Bias_File = 'BiasesAndErrorVariance.dat'

Input_Mass_Col = 5
Mass_Bias_Col = 6
#Antisymmetric error
Mass_Error_Cols = [7,8]

f = pl.figure()
ax = f.add_subplot(111)
#axStN = f.add_subplot(212)

##Specify X positions for Experiment labels
xPos = [1,2,3,4]
#Pad firs element as problem setting xtick label
Experiment_Labels = [' ','COMBO, SizeMag', 'COMBO, SizeOnly', 'STAGES, SizeMag', 'STAGES, SizeOnly']

##Read in Biases
Ex11_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[0]+Mass_Directory[0]+Bias_File)
Ex12_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[0]+Mass_Directory[1]+Bias_File)

Ex21_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[1]+Mass_Directory[0]+Bias_File)
Ex22_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[1]+Mass_Directory[1]+Bias_File)

Ex31_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[2]+Mass_Directory[0]+Bias_File)
Ex32_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[2]+Mass_Directory[1]+Bias_File)

Ex41_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[3]+Mass_Directory[0]+Bias_File)
Ex42_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[3]+Mass_Directory[1]+Bias_File)


ax.errorbar(xPos[0], Ex11_Bias[Input_Mass_Col], yerr = [[Ex11_Bias[Mass_Error_Cols[1]]],[Ex11_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 3., marker = 'x')
ax.errorbar(xPos[0], Ex12_Bias[Input_Mass_Col], yerr = [[Ex12_Bias[Mass_Error_Cols[1]]],[Ex12_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 3., marker = 'x')

ax.errorbar(xPos[1], Ex21_Bias[Input_Mass_Col], yerr = [[Ex21_Bias[Mass_Error_Cols[1]]],[Ex21_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 3., marker = 'x')
ax.errorbar(xPos[1], Ex22_Bias[Input_Mass_Col], yerr = [[Ex22_Bias[Mass_Error_Cols[1]]],[Ex22_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 3., marker = 'x')

ax.errorbar(xPos[2], Ex31_Bias[Input_Mass_Col], yerr = [[Ex31_Bias[Mass_Error_Cols[1]]],[Ex31_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 3., marker = 'x')
ax.errorbar(xPos[2], Ex32_Bias[Input_Mass_Col], yerr = [[Ex32_Bias[Mass_Error_Cols[1]]],[Ex32_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 3., marker = 'x')

ax.errorbar(xPos[3], Ex41_Bias[Input_Mass_Col], yerr = [[Ex41_Bias[Mass_Error_Cols[1]]],[Ex41_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 3., marker = 'x')
ax.errorbar(xPos[3], Ex42_Bias[Input_Mass_Col], yerr = [[Ex42_Bias[Mass_Error_Cols[1]]],[Ex42_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 3., marker = 'x')


ax.margins(0.1)
pl.subplots_adjust(bottom = 0.25)

#Plot Input Mass line
xmin, xmax = pl.xlim()
ax.plot((xmin, xmax), (Ex11_Bias[Input_Mass_Col]-Ex11_Bias[Mass_Bias_Col],Ex11_Bias[Input_Mass_Col]-Ex11_Bias[Mass_Bias_Col]), linestyle = '--', color = 'k')
ax.plot((xmin, xmax), (Ex12_Bias[Input_Mass_Col]-Ex12_Bias[Mass_Bias_Col],Ex12_Bias[Input_Mass_Col]-Ex12_Bias[Mass_Bias_Col]), linestyle = '--', color = 'k')


ax.set_xticklabels(Experiment_Labels, rotation = 45.)
ax.set_ylabel(r'Mass $\left(M_{200}\; \left[\frac{\rm M_\odot}{h}\right]\right)$')

if(do_Logscale):
    pl.yscale('log')

pl.show()
        
#Standard output is "Input, Output, Variance, Input_Ap_Mass, Output_Ap_Mass, Variance" with a space between Clusters.#
'''
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
    ax.errorbar(STAGES_Bias[i,0], (STAGES_Bias[i,1]-STAGES_Bias[i,0])/STAGES_Bias[i,0], yerr = STAGES_Bias[i,2]/STAGES_Bias[i,0], color = Experiment_Colour[1], linewidth = 3., marker = 'x')    
    
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
'''
