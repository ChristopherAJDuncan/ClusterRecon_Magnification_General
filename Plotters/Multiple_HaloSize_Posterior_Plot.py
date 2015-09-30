#For use with routines that output multiple mock realisation combined posteriors to multiple different directories (e.g. Overlap/Cluster Bias routines)
#Can also be used to Plot Combined Posteriors

import pylab as pl
import numpy as np

#Plotting Declarations
lwidth = 2.
fsize = 20.

Input_Value = [0.7, 1.1, 1.5, 1.9]

Posterior_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/SingleA901a_Bias/STAGES+COMBO/StrongLensingTests2/StrongLensing_Nopzm_unEvalPrior_CoreCut0p014/'
Halo_Head = 'r200_' #Virial radius size will be appended onto this
Posterior_Filename = 'Bias_Combined_Posterior.dat'


Posterior_Label = 'r_{200}'
Data_Label = Posterior_Label+'|\\theta, m' #'\theta'

Cluster_Colours = ['r', 'g', 'b', 'k']

f = pl.figure()
ax = f.add_subplot(111)
for i in range(0, len(Input_Value)):
    Posterior = np.genfromtxt(Posterior_Head+Halo_Head+str(Input_Value[i])+'/'+Posterior_Filename)
    ax.plot(Posterior[:,0], Posterior[:,1], label = str(Input_Value[i]), color = Cluster_Colours[i], linewidth = lwidth)

for i in range(0, len(Input_Value)):
    ax.plot((Input_Value[i],Input_Value[i]), ax.get_ylim(), color = Cluster_Colours[i], linestyle = '--', linewidth = lwidth)

ax.set_xlim((np.amin(Input_Value)- 0.2, np.amax(Input_Value)+ 0.4))
ax.set_xlabel(r'$'+Posterior_Label+'$', fontsize = fsize)
ax.set_ylabel(r'$'+('p('+Data_Label+')')+'$', fontsize = fsize)
#ax.set_ylabel('Posterior', fontsize = fsize)

#pl.legend(loc = 1)#1:UR, 2:UL, 9:UC
#pl.show()

output_name = 'Posterior_Plot.pdf'
pl.savefig(output_name,format = 'pdf',bbox_inches = 'tight', pad_inches = 0.6)
print 'Produced plot output to '+output_name
