
#Can also be used to Plot Combined Posteriors

import pylab as pl
import numpy as np

#Plotting Declarations
lwidth = 2.
fsize = 20.

Input_Value = [1.194, 1.180, 0.799, 0.894]
Cluster_Labels = ['A901a', 'A901b', 'A902', 'SW']

Do_Aperture = [0,0,1,0]
nRun = [11,15]

Posterior_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/STAGES+COMBO/Bias_Error_Run_CatID45/'
Posterior_Filename = 'Posterior_per_Aperture.dat'

#Posterior_Filename = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeOnly/COMBO/Bias/Bias_Combined_Posterior.dat'

#Posterior_Filename = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/COMBO/Bias/Bias_Error_Run_CatID5/1/Posterior_per_Aperture.dat'

Posterior_Label = 'r_{200}'

Cluster_Colours = ['r', 'g', 'b', 'k']

f = pl.figure()
ax = f.add_subplot(111)

for i in range(nRun[0], nRun[1]+1):
    Posterior = np.genfromtxt(Posterior_Head+str(i)+'/'+Posterior_Filename)

    for j in range(1,len(Do_Aperture)+1):
        if(Do_Aperture[j-1] == 0):
            continue
        ax.plot(Posterior[:,0], Posterior[:,j], label = Cluster_Labels[j-1]+'/'+str(i), linewidth = lwidth)

        ax.plot((Input_Value[j-1],Input_Value[j-1]), ax.get_ylim(), linestyle = '--', linewidth = lwidth)

ax.set_xlim((0.6, 1.5))
#ax.set_xlim((np.amin(Input_Value)- 0.2, np.amax(Input_Value)+ 0.3))
ax.set_xlabel(r'$'+Posterior_Label+'$', fontsize = fsize)
ax.set_ylabel(r'$p('+Posterior_Label+')$', fontsize = fsize)
#ax.set_ylabel('Posterior', fontsize = fsize)

pl.legend(loc = 1)#1:UR, 2:UL, 9:UC
#pl.show()

output_name = 'Mock_Single_Posterior_Plot.pdf'
pl.savefig(output_name,format = 'pdf',bbox_inches = 'tight', pad_inches = 0.6)
print 'Produced plot output to '+output_name
