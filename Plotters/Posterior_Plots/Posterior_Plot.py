
#Can also be used to Plot Combined Posteriors

import pylab as pl
import numpy as np

#Plotting Declarations
lwidth = 2.
fsize = 20.

Input_Value = [1.194, 1.180, 0.799, 0.894]
Cluster_Labels = ['A901a', 'A901b', 'A902', 'SW']

Posterior_Filename = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/MagCutsTesting_Data/SingleClusterFit/4Cluster/RRG+STAGES/SM+M/M_less_27p5/Run2/RRGMag/Posterior_per_Aperture.dat'
#'/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Data-Run_4Jul/SizeMag/Data/STAGES+COMBO/detJ_P0.0_D2.2/Posterior_per_Aperture.dat'
#

#Posterior_Filename = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/Cluster_Bias/COMBO/STAGES_Bias/Bias_Combined_Posterior.dat'

#Posterior_Filename = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/COMBO/Bias/Bias_Error_Run_CatID5/1/Posterior_per_Aperture.dat'

Posterior_Label = 'r_{200}'
Data_Label = Posterior_Label+'|m' #'|\\theta, m'

Cluster_Colours = ['r', 'g', 'b', 'k']

Posterior = np.genfromtxt(Posterior_Filename)


f = pl.figure()
ax = f.add_subplot(111)
for i in range(1, Posterior.shape[1]):
    ax.plot(Posterior[:,0], Posterior[:,i], label = Cluster_Labels[i-1], color = Cluster_Colours[i-1], linewidth = lwidth)

for i in range(1, Posterior.shape[1]):
    ax.plot((Input_Value[i-1],Input_Value[i-1]), ax.get_ylim(), color = Cluster_Colours[i-1], linestyle = '--', linewidth = lwidth)

ax.set_xlim((0., 2.))
#ax.set_xlim((np.amin(Input_Value)- 0.2, np.amax(Input_Value)+ 0.3))
ax.set_xlabel(r'$r_{200} \; [h^{-1}{\rm Mpc}]$', fontsize = fsize)
ax.set_ylabel(r'$'+('p('+Data_Label+')')+'$', fontsize = fsize)
#ax.set_ylabel('Posterior', fontsize = fsize)

pl.legend(loc = 1)#1:UR, 2:UL, 9:UC
#pl.show()

output_name = 'Posterior_Plot.pdf'
pl.savefig(output_name,format = 'pdf',bbox_inches = 'tight', pad_inches = 0.6)
print 'Produced plot output to '+output_name
