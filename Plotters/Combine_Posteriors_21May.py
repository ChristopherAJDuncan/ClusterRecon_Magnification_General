
#Can also be used to Plot Combined Posteriors

import pylab as pl
import numpy as np

#Plotting Declarations
lwidth = 2.
fsize = 20.

Input_Value = [1.194, 1.180, 0.799, 0.894]
Cluster_Labels = ['A901a', 'A901b', 'A902', 'SW']

Do_Aperture = [0,0,0,1]
nRun = [1,20]
nIgnore = []

Posterior_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/COMBO_Revert_6_KappaRenormON/Bias_Error_Run_CatID5/'
#Posterior_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/COMBO_Revert_5_CombinedPosterior_Lookup2/Bias_Error_Run_CatID5/'
Posterior_Filename = 'Posterior_per_Aperture.dat'

#Posterior_Filename = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeOnly/COMBO/Bias/Bias_Combined_Posterior.dat'

#Posterior_Filename = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/COMBO/Bias/Bias_Error_Run_CatID5/1/Posterior_per_Aperture.dat'

Posterior_Label = 'r_{200}'

Cluster_Colours = ['r', 'g', 'b', 'k']

#def get_Outlier_Removal_Indexes(Post

def Combine(Combined, Post):
    NCombined = np.zeros(Post.shape[0])
    for i in range(0, Post.shape[0]):
        if(Post[i] == 0.):
            NCombined[i] = Combined[i] - 100.
        else:
            NCombined[i] = Combined[i] + np.log(Post[i])

    return NCombined
            
def Renorm_lnLikelihood(L):
    for j in range(0, L.shape[1]):
        L[:,j] = L[:,j] - np.amax(L[:,j])
    return L
    

f = pl.figure()
ax = f.add_subplot(211)
ax_C = f.add_subplot(212)

for i in range(nRun[0], nRun[1]+1):
    if(i in nIgnore):
        print 'Skipping Run:', i
        continue
    Posterior = np.genfromtxt(Posterior_Head+str(i)+'/'+Posterior_Filename)

    if(i==nRun[0]):
        #Grid, Aperture
        Combined_Posterior = np.zeros(( Posterior.shape[0], Posterior.shape[1]-1))
        Combined_Posterior[:,:] = -100.

    for j in range(1,len(Do_Aperture)+1):
        Combined_Posterior[:,j-1] = np.copy(Combine(Combined_Posterior[:,j-1], Posterior[:,j]))
        
        if(Do_Aperture[j-1] == 0):
            continue
        ax.plot(Posterior[:,0], Posterior[:,j], label = Cluster_Labels[j-1]+'/'+str(i), linewidth = lwidth)

        ax.plot((Input_Value[j-1],Input_Value[j-1]), ax.get_ylim(), linestyle = '--', linewidth = lwidth)

Combined_Posterior = Renorm_lnLikelihood(Combined_Posterior)
Combined_Posterior = np.exp(Combined_Posterior)

for j in range(0,len(Do_Aperture)):
    ax_C.plot(Posterior[:,0],Combined_Posterior[:,j], label = Cluster_Labels[j], color = Cluster_Colours[j])
    ax_C.plot((Input_Value[j],Input_Value[j]), ax_C.get_ylim(), linestyle = '--', linewidth = lwidth, color = Cluster_Colours[j])

ax.set_xlim((0.6, 1.5))
ax_C.set_xlim(ax.get_xlim())
#ax.set_xlim((np.amin(Input_Value)- 0.2, np.amax(Input_Value)+ 0.3))
ax.set_xlabel(r'$'+Posterior_Label+'$', fontsize = fsize)
ax.set_ylabel(r'$p('+Posterior_Label+')$', fontsize = fsize)
#ax.set_ylabel('Posterior', fontsize = fsize)


pl.legend(loc = 1)#1:UR, 2:UL, 9:UC
#pl.show()

output_name = 'Combined_Posterior_Plot.pdf'
pl.savefig(output_name,format = 'pdf',bbox_inches = 'tight', pad_inches = 0.6)
print 'Produced plot output to '+output_name
