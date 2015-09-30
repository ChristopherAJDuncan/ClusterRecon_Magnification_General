#Produces ML-point-estimate plots of posteriors of a Combined run, with the aim of showing no bias for a given experiment

import numpy as np
import pylab as pl
import math
import sys
import os

getPDF = True

Experiment_Colour = ['r','b']
#STAGES_Cluster_Label = ['A901a', 'A901b', 'A902', 'SW']

######READ in PARAMETER FILE#######

Bias_Directory_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Results/3Jul2014/SizeMag/SingleCluster_Bias/STAGES+COMBO/Pzm_Bias/'
Mass_Directories = ['r200_0.4/', 'r200_0.8/', 'r200_1.2/', 'r200_1.6/']
Bias_File = 'BiasesAndErrorVariance.dat'
#Bias_File = 'Mass_Conversion_Errors.dat'

Experiment_Label = ['STAGES, SizeMag']

Recovered_Mass_Col = 5
Mass_Bias_Col = 6
#Antisymmetric error
Mass_Error_Cols = [7,8]

f = pl.figure()
ax = f.add_subplot(111)

for c in range(0, len(Mass_Directories)):
    Bias_Input = np.genfromtxt(Bias_Directory_Head+Mass_Directories[c]+Bias_File)

    ax.errorbar(Bias_Input[Recovered_Mass_Col]-Bias_Input[Mass_Bias_Col], Bias_Input[Mass_Bias_Col]/(Bias_Input[Recovered_Mass_Col]-Bias_Input[Mass_Bias_Col]), yerr = [[Bias_Input[Mass_Error_Cols[1]]/(Bias_Input[Recovered_Mass_Col]-Bias_Input[Mass_Bias_Col])],[Bias_Input[Mass_Error_Cols[0]]/(Bias_Input[Recovered_Mass_Col]-Bias_Input[Mass_Bias_Col])]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')

#Plot zero line
ax.plot(ax.get_xlim(), (0.,0.), linestyle = '--')

ax.legend()
ax.set_xlabel(r'Input Virial Mass $\left(M_{200}\; \left[\frac{\rm M_\odot}{h}\right]\right)$'); ax.set_ylabel(r'Fractional Bias in Recovered Mass')

if(getPDF):
    output_name = 'Single_Cluster_noBias_Plot.pdf'
    pl.savefig(output_name,format = 'pdf',bbox_inches = 'tight')
    print 'Produced plot output to '+output_name
    os.system('evince '+output_name+' &')
else:
    pl.show()

        
