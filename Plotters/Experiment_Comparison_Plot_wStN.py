#Produces ML-point-estimate plots of posteriors for 4 experiments, for two masses. Uses the same output as the Non-Biased, Single Cluster Plot
##Files need edited to plot a single run

import numpy as np
from scipy import stats
import pylab as pl
import math
import sys
import os

getPDF = True

Experiment_Colour = ['r','b']
#STAGES_Cluster_Label = ['A901a', 'A901b', 'A902', 'SW']

do_Logscale = False

######READ in PARAMETER FILE#######
nMass = 2
Bias_Directory_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Results/3Jul2014/SizeMag/SingleCluster_Bias/Experiment_Comparison/'
Experiment_Directory = ['COMBO_Size-Mag/', 'STAGES_Mag-Only/', 'STAGES_Size-Only/', 'SL_Unbiased/', 'Master_Cat/SM/', 'Master_Cat/SM+M/']
Mass_Directory = ['r200_1.2/', 'r200_1.6/']
Run_Directory = ['Bias_Error_Run_CatID5/', 'Bias_Error_Run_CatID45/', 'Bias_Error_Run_CatID45/', 'Bias_Error_Run_CatID45/', 'Bias_Error_Run_CatID45/', 'Bias_Error_Run_CatID45/']

Run_Plot = [4,3,1,4,1,1]

#Average_Decalrations
Average_Mass_Col = 0
Average_nRun = 10

Bias_File = 'Mass_Estimates.dat'

#Masses
Input_Mass_Col = 5
Recovered_Mass_Col = 6
#Antisymmetric error
Mass_Error_Cols = [7,8]


'''
#r200
Input_Mass_Col = 1
Recovered_Mass_Col = 2
#Antisymmetric error
Mass_Error_Cols = [3,4]   
'''

def Average_StN(Dir_head, Bias_File, nRun, Mass_Col, Mass_Err_Col):
    #Ex12_Bias[Input_Mass_Col]/(0.5*(Ex12_Bias[Mass_Error_Cols[0]]+Ex12_Bias[Mass_Error_Cols[1]]))
    StN = np.zeros(nRun); ErrorWidth = np.zeros(nRun)
    for i in range(1, nRun+1):
        Input = np.genfromtxt(Dir_head+str(i)+'/'+Bias_File)
        ErrorWidth[i-1] = (0.5*(Input[Mass_Error_Cols[0]]+Input[Mass_Error_Cols[1]]))
        StN[i-1] = Input[Mass_Col]/ErrorWidth[i-1]


    print 'StN:', stats.nanmean(StN)
    print 'Average Error width:', stats.nanmean(ErrorWidth)
    print ' '
    return stats.nanmean(StN)

f = pl.figure()

ax = pl.subplot2grid( (6,6), (0,0), colspan = 6, rowspan = 3)
axStN = pl.subplot2grid( (6,6), (3,0), colspan = 6, rowspan = 3) 

#ax = f.add_subplot(211)
#axStN = f.add_subplot(212)

##Specify X positions for Experiment labels
xPos = [1,2,3,4]
#Pad firs element as problem setting xtick label
Experiment_Labels = ['COMBO, SizeMag', 'STAGES, Mag-Only', 'STAGES, Size-Only', 'STAGES, Size-Mag']

##Read in Biases
#ij = Experiment, Mass
Ex11_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[0]+Mass_Directory[0]+Run_Directory[0]+str(Run_Plot[0])+'/'+Bias_File)
Ex12_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[0]+Mass_Directory[1]+Run_Directory[0]+str(Run_Plot[0])+'/'+Bias_File)

Ex21_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[1]+Mass_Directory[0]+Run_Directory[1]+str(Run_Plot[1])+'/'+Bias_File)
Ex22_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[1]+Mass_Directory[1]+Run_Directory[1]+str(Run_Plot[1])+'/'+Bias_File)

Ex31_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[2]+Mass_Directory[0]+Run_Directory[2]+str(Run_Plot[2])+'/'+Bias_File)
Ex32_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[2]+Mass_Directory[1]+Run_Directory[2]+str(Run_Plot[2])+'/'+Bias_File)

Ex41_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[3]+Mass_Directory[0]+Run_Directory[3]+str(Run_Plot[3])+'/'+Bias_File)
Ex42_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[3]+Mass_Directory[1]+Run_Directory[3]+str(Run_Plot[3])+'/'+Bias_File)

Ex51_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[4]+Mass_Directory[0]+Run_Directory[3]+str(Run_Plot[3])+'/'+Bias_File)
Ex52_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[4]+Mass_Directory[1]+Run_Directory[3]+str(Run_Plot[3])+'/'+Bias_File)

Ex61_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[5]+Mass_Directory[0]+Run_Directory[3]+str(Run_Plot[3])+'/'+Bias_File)
Ex62_Bias = np.genfromtxt(Bias_Directory_Head+Experiment_Directory[5]+Mass_Directory[1]+Run_Directory[3]+str(Run_Plot[3])+'/'+Bias_File)


ax.errorbar(xPos[0], Ex11_Bias[Recovered_Mass_Col], yerr = [[Ex11_Bias[Mass_Error_Cols[1]]],[Ex11_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')
ax.errorbar(xPos[0], Ex12_Bias[Recovered_Mass_Col], yerr = [[Ex12_Bias[Mass_Error_Cols[1]]],[Ex12_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')

ax.errorbar(xPos[1], Ex21_Bias[Recovered_Mass_Col], yerr = [[Ex21_Bias[Mass_Error_Cols[1]]],[Ex21_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')
ax.errorbar(xPos[1], Ex22_Bias[Recovered_Mass_Col], yerr = [[Ex22_Bias[Mass_Error_Cols[1]]],[Ex22_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')

ax.errorbar(xPos[2], Ex31_Bias[Recovered_Mass_Col], yerr = [[Ex31_Bias[Mass_Error_Cols[1]]],[Ex31_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')
ax.errorbar(xPos[2], Ex32_Bias[Recovered_Mass_Col], yerr = [[Ex32_Bias[Mass_Error_Cols[1]]],[Ex32_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')

ax.errorbar(xPos[3], Ex41_Bias[Recovered_Mass_Col], yerr = [[Ex41_Bias[Mass_Error_Cols[1]]],[Ex41_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')
ax.errorbar(xPos[3], Ex42_Bias[Recovered_Mass_Col], yerr = [[Ex42_Bias[Mass_Error_Cols[1]]],[Ex42_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')

'''
ax.errorbar(xPos[3], Ex51_Bias[Recovered_Mass_Col], yerr = [[Ex51_Bias[Mass_Error_Cols[1]]],[Ex51_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')
ax.errorbar(xPos[3], Ex52_Bias[Recovered_Mass_Col], yerr = [[Ex52_Bias[Mass_Error_Cols[1]]],[Ex52_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')

ax.errorbar(xPos[3], Ex61_Bias[Recovered_Mass_Col], yerr = [[Ex61_Bias[Mass_Error_Cols[1]]],[Ex61_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')
ax.errorbar(xPos[3], Ex62_Bias[Recovered_Mass_Col], yerr = [[Ex62_Bias[Mass_Error_Cols[1]]],[Ex62_Bias[Mass_Error_Cols[0]]]], color = Experiment_Colour[0], linewidth = 1.5, marker = 'x')
'''

##Plot StN for Mass 2 (Larger Mass)
#axStN.plot(xPos[0],  Ex12_Bias[Input_Mass_Col]/(0.5*(Ex12_Bias[Mass_Error_Cols[0]]+Ex12_Bias[Mass_Error_Cols[1]])), 'ro')
#axStN.plot(xPos[1],  Ex22_Bias[Input_Mass_Col]/(0.5*(Ex22_Bias[Mass_Error_Cols[0]]+Ex22_Bias[Mass_Error_Cols[1]])), 'ro')
#axStN.plot(xPos[2],  Ex32_Bias[Input_Mass_Col]/(0.5*(Ex32_Bias[Mass_Error_Cols[0]]+Ex32_Bias[Mass_Error_Cols[1]])), 'ro')
#axStN.plot(xPos[3],  Ex42_Bias[Input_Mass_Col]/(0.5*(Ex42_Bias[Mass_Error_Cols[0]]+Ex42_Bias[Mass_Error_Cols[1]])), 'ro')

axStN.plot(xPos[0],  Average_StN(Bias_Directory_Head+Experiment_Directory[0]+Mass_Directory[Average_Mass_Col]+Run_Directory[0], Bias_File, Average_nRun, Input_Mass_Col, Mass_Error_Cols), 'rx')
axStN.plot(xPos[1],  Average_StN(Bias_Directory_Head+Experiment_Directory[1]+Mass_Directory[Average_Mass_Col]+Run_Directory[1], Bias_File, Average_nRun, Input_Mass_Col, Mass_Error_Cols), 'rx')
axStN.plot(xPos[2],  Average_StN(Bias_Directory_Head+Experiment_Directory[2]+Mass_Directory[Average_Mass_Col]+Run_Directory[2], Bias_File, Average_nRun, Input_Mass_Col, Mass_Error_Cols), 'rx')
axStN.plot(xPos[3],  Average_StN(Bias_Directory_Head+Experiment_Directory[3]+Mass_Directory[Average_Mass_Col]+Run_Directory[3], Bias_File, Average_nRun, Input_Mass_Col, Mass_Error_Cols), 'rx', label = 'All sizes')
axStN.plot(xPos[3],  Average_StN(Bias_Directory_Head+Experiment_Directory[5]+Mass_Directory[Average_Mass_Col]+Run_Directory[5], Bias_File, Average_nRun, Input_Mass_Col, Mass_Error_Cols), 'b^', label = 'GALFIT sample, SM+M')
axStN.plot(xPos[3],  Average_StN(Bias_Directory_Head+Experiment_Directory[4]+Mass_Directory[Average_Mass_Col]+Run_Directory[4], Bias_File, Average_nRun, Input_Mass_Col, Mass_Error_Cols), 'ko', label = 'GALFIT sample')
axStN.legend(loc = 2, prop = {'size':11}, numpoints = 1)



axStN.set_ylim(axStN.get_ylim()[0]-(0.1*np.diff(axStN.get_ylim())), axStN.get_ylim()[1]+(0.1*np.diff(axStN.get_ylim())))

pl.subplots_adjust(bottom = 0.25)

#Plot Input Mass line
xmin, xmax = xPos[0], xPos[-1]
ax.plot((xmin, xmax), (Ex11_Bias[Input_Mass_Col],Ex11_Bias[Input_Mass_Col]), linestyle = '--', color = 'k')
ax.plot((xmin, xmax), (Ex12_Bias[Input_Mass_Col],Ex12_Bias[Input_Mass_Col]), linestyle = '--', color = 'k')


ax.set_xticks(xPos); axStN.set_xticks(xPos)
ax.set_ylabel(r'Mass $\left(M_{200}\; \left[\frac{\rm M_\odot}{h}\right]\right)$', fontsize = 14.)
pl.setp(ax.get_xticklabels(),visible = False)

axStN.set_xticklabels(Experiment_Labels, rotation = 45., fontsize = 12.)
axStN.tick_params(axis = 'y', which = 'major', labelsize = 12.)
ax.tick_params(axis = 'y', which = 'major', labelsize = 12.)
axStN.set_ylabel('Signal to Noise', fontsize = 14.)

axStN.set_yticks(axStN.get_yticks()[::2])

ax.margins(0.1); axStN.set_xlim(ax.get_xlim())

if(do_Logscale):
    pl.yscale('log')

if(getPDF):
    output_name = 'Experiment_Comparison_wStN.pdf'
    pl.savefig(output_name,format = 'pdf',bbox_inches = 'tight')
    print 'Produced plot output to '+output_name
    os.system('evince -w '+output_name+' &')
else:
    pl.show()

        
