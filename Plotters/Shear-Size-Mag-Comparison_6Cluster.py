import numpy as np
import pylab as pl

Plot_Cluster = [1,0,1,1,1,1]
Cl_Marker = ['x','x','x','x', 'x', 'x']
Cl_Color =  ['r', 'g', 'b', 'k', 'c', 'm']
Cl_Label = ['A901a', 'A901b', 'A902', 'CB1', 'SWa', 'SWb']

###Shear Results
#Taken from Table 1 of Heymans 2008, one halo fits
Shear_r200 = [1.194, 1.180, 0.766, 0.608, 0.689, 0.742]
Shear_r200_Stat_Error = [[0.086, 0.101], [0.088, 0.104], [0.137, 0.140], [0.127, 0.138], [0.141,0.192], [0.123, 0.153]]


### Magnification Results
Mag_Dir = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Data/SizeMag+Mag/SingleClusterFit/RRG+GALFIT/SM+M/'
#Mag_Dir = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Data/SizeMag+Mag/SingleClusterFit/GALFIT/GALIT_MAG/wGEMS/M_PriorCutTEST2/'
Mag_Filename = 'Mass_Estimates.dat'
Title = ''

## Input Column Information
Rec_Col = 2
Input__Col = 1
Err_Col = [3,4]

Mag = np.genfromtxt(Mag_Dir+Mag_Filename)

f = pl.figure()
ax = f.add_subplot(1,1,1)

    #Plot r200 Mag (x) vs r200 Shear (y)

for i in range(0, len(Plot_Cluster)):
    if(Plot_Cluster[i] == 0):
        continue
    
    ax.errorbar(Mag[i,Rec_Col], Shear_r200[i], xerr = [[Mag[i,Err_Col[1]]],[Mag[i,Err_Col[0]]]], yerr = [[Shear_r200_Stat_Error[i][1]],[Shear_r200_Stat_Error[i][0]]], marker = Cl_Marker[i], color = Cl_Color[i], label = Cl_Label[i])#, xerr = (Mag[i,Err_Col[1]],Mag[i,Err_Col[0]]))#, yerr = (Shear_r200_Stat_Error[i][1], Shear_r200_Stat_Error[i][0]), 

ax.set_xlim(0.0, 1.5)

yeqx = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 50)
ax.plot(yeqx, yeqx, linestyle = '--', color = 'k')

ax.set_xlim(yeqx[0], yeqx[-1]);ax.set_ylim(yeqx[0], yeqx[-1])

ax.set_title(Title, fontsize = 16)
ax.set_xlabel(r'$r_{200}\; [h^{-1} {\rm Mpc}] \; ({\rm Magnification})$', fontsize = 16)
ax.set_ylabel(r'$r_{200}\; [h^{-1} {\rm Mpc}] \; ({\rm Shear})$', fontsize = 16)

ax.legend(loc = 4)
pl.show()
