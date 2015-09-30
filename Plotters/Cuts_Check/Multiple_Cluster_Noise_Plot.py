#Produces ML-point-estimate plots of posteriors of a Combined run, with the aim of showing no bias for a given experiment

import numpy as np
import pylab as pl
import math
import sys
import os

getPDF = True

Experiment_Colour = ['r', 'b', 'g', 'c']
nCluster = 4

######READ in PARAMETER FILE#######
Bias_Directory_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/STAGES_Clusters_Bias/STAGES+COMBO/WLWL/SM+Mag_P3.3_D6.5/Bias_Error_Run_CatID45/'
nRealisation = 20
Bias_File = 'Mass_Estimates.dat'

Experiment_Label = ['STAGES, SizeMag']

Recovered_Mass_Col = 2
Input_Mass_Col = 1
#Antisymmetric error
Mass_Error_Cols = [3,4]

f = pl.figure()

xPos = range(1, nRealisation+1)

for a in range(0, nCluster):
    ax = f.add_subplot(nCluster,1,a+1)

    if(ax != nCluster-1):
        pl.setp(ax.get_xticklabels(), visible = False)
    for c in range(0, nRealisation):
        Bias_Input = np.genfromtxt(Bias_Directory_Head+str(c+1)+'/'+Bias_File)

        print 'Doing ', Bias_Directory_Head+str(c+1)+'/'+Bias_File

        if(c==nRealisation-1):
            ax.set_xlim(0, nRealisation+1)
            ax.plot(ax.get_xlim(), (Bias_Input[a,Input_Mass_Col],Bias_Input[a,Input_Mass_Col]), linestyle = '--')
            
        ax.errorbar(xPos[c], Bias_Input[a,Recovered_Mass_Col], yerr = [[Bias_Input[a,Mass_Error_Cols[1]]],[Bias_Input[a,Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')



#Plot Input Mass Line


ax.set_xlabel(r'Mock Realisation ');


if(getPDF):
    output_name = 'Multiple_Cluster_Noise.pdf'
    pl.savefig(output_name,format = 'pdf',bbox_inches = 'tight')
    print 'Produced plot output to '+output_name
    os.system('evince '+output_name+' &')
else:
    pl.show()

                                  
