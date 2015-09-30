

import numpy as np
import pylab as pl
import math
import sys

Directory = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/Data/STAGES/Calibrated_TrJ/Single_Cluster/Split_Catalogue/TrJ/Mag-Only_SNR30+_wCOMBO_zLimit_DELETE/'
#'/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Data-Run_31Jul/SizeMag/Data/STAGES+COMBO/P0.0_D2.2/'
File = 'Foreground_Contamination_NumberDensity.dat'

xLimit = [0., 4.]
ls = '--'#'none'

print 'REMEMBER to enter Source Number!!'

nSources = 74477 #74421
UnmaskedArea = 822.56657778 #1089. for Mocks #822.566577778 #Arcminutes Sq# Default = 822.566577778 for the STAGES field, accounting for masking

## 6-Cluster
#Cluster_Labels = ['A901a', 'A901b', 'A902', 'CB1', 'SWa', 'SWb']

##4 Cluster
Cluster_Labels = ['A901a', 'A901b', 'A902', 'SW']
Cluster_Colors = ['b','g','r','c', 'k', 'y']

Input = np.genfromtxt(Directory+File)
print 'Got Input'

def Get_Global_Average(In,xPoints, limit):
    Aver = np.zeros(In.shape[1]-2)
    nPoint = np.zeros(In.shape[1]-2)
    for Ap in range(0, In.shape[1]-2):
        for i in range(0, In.shape[0]):
            if(xPoints[i] >= limit):
                Aver[Ap] += In[i,Ap+2]
                nPoint[Ap] += 1
        Aver[Ap] = Aver[Ap]/nPoint[Ap]
        print nPoint[Ap]
    return Aver
                                     

Area = np.zeros(Input.shape[0])
Errors = np.zeros((Input.shape[0],Input.shape[1]-2))
for i in range (0, Area.shape[0]):
    Area[i] = math.pi*((Input[i,1]**2.)-(Input[i,0]**2))
    for Ap in range (0, Errors.shape[1]):
        Errors[i, Ap] = math.sqrt(Input[i,Ap+2]*Area[i])/Area[i]

#Plot

f = pl.figure()


xPoints = np.zeros(Errors.shape[0])
for i in range(xPoints.shape[0]):
    xPoints[i] = 0.5*(Input[i,0]+Input[i,1])

#Get Averages:
Global_Average = nSources/UnmaskedArea
#Global_Averages = Get_Global_Average(Input, xPoints, 1.5)

for Ap in range (0,Errors.shape[1]):
    ax = f.add_subplot(2,2,Ap+1)
    ax.errorbar(xPoints, Input[:,Ap+2], yerr = Errors[:,Ap], label = Cluster_Labels[Ap], color = Cluster_Colors[Ap], ls = ls, marker = 'o', mec = 'white')

    ax.axhline(Global_Average, linestyle = '-', color = 'k')
#    ax.axhline(Global_Averages[Ap], linestyle = '--', color = 'k')

#    ax.set_xlim(0.,1.5)
    ax.legend(loc = 1)

    ax.set_xlim(xLimit)

    if(ax.is_last_row()):
    #if(Ap == 2 or Ap == 3 or Ap == 4 or Ap == 5):
        print ' '
        #ax.set_yscale('log')
        ax.set_xlim(0.,2.0)
    if(Ap == 0 or Ap == 2 or Ap == 4):
        ax.set_ylabel(r'Number Density [${\rm arcminute}^{-2}$]', fontsize = 16)
    if(Ap == 3 or Ap == 2):
    #if(Ap == 4 or Ap == 5):
        ax.set_xlabel('Annulus Radius [Arcminute]', fontsize = 16)

    pl.setp(ax.get_xticklabels(), rotation = 45, fontsize = 14)
    pl.setp(ax.get_yticklabels(), fontsize = 14)



pl.show()
