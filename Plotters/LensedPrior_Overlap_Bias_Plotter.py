#Produces ML-point-estimate plots of posteriors of a Combined run, with the aim of showing the bias coming from Overlap and lensed prior

import numpy as np
import pylab as pl
import math
import sys
import os

getPDF = True
Experiment_Colour = ['r','b']
Plot_Type = 1 #0:Bias, 1:Fractional Bias
linewid  = 1.
#STAGES_Cluster_Label = ['A901a', 'A901b', 'A902', 'SW']

######READ in PARAMETER FILE#######
Input_Value = [0.4, 0.8, 1.2, 1.6]
Lensed_Bias_Directory_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Results/3Jul2014/SizeMag/SingleCluster_Bias/STAGES+COMBO/LensedPrior/'
Lensed_Bias_File = 'BiasesAndErrorVariance.dat'

Overlap_Bias_Directory_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Results/3Jul2014/SizeMag/Overlap_Bias/STAGES+COMBO/'
Overlap_Bias_File = 'BiasesAndErrorVariance.dat' #Lensed_Bias_File
Overlap_Cluster = 2
##Convert to python labelling
Overlap_Cluster -= 1

'''
Recovered_Mass_Col = 1
Mass_Bias_Col = 2
#Antisymmetric error
Mass_Error_Cols = [3,4]
'''


Recovered_Mass_Col = 5
Mass_Bias_Col = 6
#Antisymmetric error
Mass_Error_Cols = [7,8]


f = pl.figure()
ax = f.add_subplot(111)

for c in range(1, len(Input_Value)+1):
    Lensed_Bias_Input = np.genfromtxt(Lensed_Bias_Directory_Head+'r200_'+str(Input_Value[c-1])+'/'+Lensed_Bias_File)
    Overlap_Bias_Input = np.genfromtxt(Overlap_Bias_Directory_Head+'r200_'+str(Input_Value[c-1])+'/'+Overlap_Bias_File)

    Lensed_Input_Mass = Lensed_Bias_Input[Recovered_Mass_Col]-Lensed_Bias_Input[Mass_Bias_Col]
    Overlap_Input_Mass = Overlap_Bias_Input[Overlap_Cluster, Recovered_Mass_Col]-Overlap_Bias_Input[Overlap_Cluster, Mass_Bias_Col]

    print 'Lensed Input:', Lensed_Input_Mass, Lensed_Bias_Input[Recovered_Mass_Col], Lensed_Bias_Input[Mass_Bias_Col]

    ##Points to Plot#
    if(Plot_Type == 0):
        Lensed_Bias_Point =  Lensed_Bias_Input[Mass_Bias_Col]; Lensed_Bias_Error = [[Lensed_Bias_Input[Mass_Error_Cols[1]]],[Lensed_Bias_Input[Mass_Error_Cols[0]]]]
        Overlap_Bias_Point =  Overlap_Bias_Input[Overlap_Cluster,Mass_Bias_Col]; Overlap_Bias_Error = [[Overlap_Bias_Input[Overlap_Cluster,Mass_Error_Cols[1]]],[Overlap_Bias_Input[Overlap_Cluster,Mass_Error_Cols[0]]]]  
    elif(Plot_Type ==1):
        Lensed_Bias_Point =  Lensed_Bias_Input[Mass_Bias_Col]/Lensed_Input_Mass; Lensed_Bias_Error = [[Lensed_Bias_Input[Mass_Error_Cols[1]]/Lensed_Input_Mass],[Lensed_Bias_Input[Mass_Error_Cols[0]]/Lensed_Input_Mass]]
        Overlap_Bias_Point =  Overlap_Bias_Input[Overlap_Cluster,Mass_Bias_Col]/Overlap_Input_Mass; Overlap_Bias_Error = [[Overlap_Bias_Input[Overlap_Cluster,Mass_Error_Cols[1]]/Overlap_Input_Mass],[Overlap_Bias_Input[Overlap_Cluster,Mass_Error_Cols[0]]/Overlap_Input_Mass]]  
    else:
        print 'Error Setting Plot Type'
        exit()


    if(c==1):
        ax.errorbar(Lensed_Input_Mass, Lensed_Bias_Point, yerr = Lensed_Bias_Error, color = Experiment_Colour[0], linewidth = linewid, marker = 'x', label = 'Lensed Prior')
        ax.errorbar(Overlap_Input_Mass, Overlap_Bias_Point, yerr = Overlap_Bias_Error, color = Experiment_Colour[1], linewidth = linewid, marker = 'o', label = 'Overlap')
    else:
        ax.errorbar(Lensed_Input_Mass, Lensed_Bias_Point, yerr = Lensed_Bias_Error, color = Experiment_Colour[0], linewidth = linewid, marker = 'x')
        ax.errorbar(Overlap_Input_Mass, Overlap_Bias_Point, yerr = Overlap_Bias_Error, color = Experiment_Colour[1], linewidth = linewid, marker = 'o')

        
#Plot zero line
ax.plot(ax.get_xlim(), (0.,0.), linestyle = '--')

ax.legend(loc = 1)
ax.set_xlabel(r'Input Virial Mass ($M_{200}\; \left[\frac{\rm M_\odot}{h}\right]$)');
if(Plot_Type ==0):
    ax.set_ylabel(r'Bias in Recovered Mass ($M_{200}\; \left[\frac{\rm M_\odot}{h}\right]$)')
elif(Plot_Type == 1):
    ax.set_ylabel(r'Fractional Bias in Recovered Mass')


if(getPDF):
    output_name = 'ClusterBias_PointPlot.pdf'
    pl.savefig(output_name,format = 'pdf',bbox_inches = 'tight')
    print 'Produced plot output to '+output_name
    os.system('evince -w '+output_name+' &')
else:
    pl.show()

        
#Standard output is "Input, Output, Variance, Input_Ap_Mass, Output_Ap_Mass, Variance" with a space between Clusters.#
