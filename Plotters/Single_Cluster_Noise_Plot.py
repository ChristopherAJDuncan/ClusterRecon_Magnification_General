#Produces ML-point-estimate plots of posteriors of a Combined run, with the aim of showing no bias for a given experiment

import numpy as np
import pylab as pl
import math
import sys
import os

getPDF = True

Experiment_Colour = ['r']

######READ in PARAMETER FILE#######
Bias_Directory_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Results/3Jul2014/SizeMag/SingleCluster_Bias/Experiment_Comparison/SL_Unbiased/r200_1.2/Bias_Error_Run_CatID45/'
nRealisation = 10
Bias_File = 'Mass_Estimates.dat'#'Mass_Estimates.dat'

Experiment_Label = ['STAGES, SizeMag']
'''
Recovered_Mass_Col = 2
Input_Mass_Col = 1
#Antisymmetric error
Mass_Error_Cols = [3,4]
'''
#Masses
Recovered_Mass_Col = 6
Input_Mass_Col = 5
#Antisymmetric error
Mass_Error_Cols = [7,8]


f = pl.figure()
ax = f.add_subplot(111)

xPos = range(1, nRealisation+1)

for c in range(0, nRealisation):
    Bias_Input = np.genfromtxt(Bias_Directory_Head+str(c+1)+'/'+Bias_File)


    ax.errorbar(xPos[c], Bias_Input[Recovered_Mass_Col], yerr = [[Bias_Input[Mass_Error_Cols[1]]],[Bias_Input[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')

ax.set_xlim(0, nRealisation+1)
#Plot Input Mass Line
ax.plot(ax.get_xlim(), (Bias_Input[Input_Mass_Col],Bias_Input[Input_Mass_Col]), linestyle = '--')

ax.set_xlabel(r'Mock Realisation '); ax.set_ylabel(r'Recovered Mass $\left(M_{200}\; \left[\frac{\rm M_\odot}{h}\right]\right)$')


if(getPDF):
    output_name = 'Single_Cluster_Noise.pdf'
    pl.savefig(output_name,format = 'pdf',bbox_inches = 'tight')
    print 'Produced plot output to '+output_name
    os.system('evince '+output_name+' &')
else:
    pl.show()

                                  
