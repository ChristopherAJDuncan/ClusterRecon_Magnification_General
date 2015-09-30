
#Can also be used to Plot Combined Posteriors

import pylab as pl
import numpy as np

#Plotting Declarations
lwidth = 2.
fsize = 20.
getPDF = False

Input_Value = [1.6]
#Input_Value = [1.194, 1.180, 0.799, 0.894]
Cluster_Labels = ['A901a', 'A901b', 'A902', 'SW']


Do_Aperture = [1]
#Do_Aperture = [1,0,0,0]
nRun = [0,15]
nIgnore = []
Remove_Outliers = 0

Posterior_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/BiasCheck/CutsBias/Known_Redshift/PzTest/SM/SizeCut_7.5+_8.5-_NoMagCut/STAGES+COMBO/r200_1.6/Bias_Error_Run_CatID45/'
#Posterior_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/STAGES_Revert_6_KappaRenormON/Bias_Error_Run_CatID4/'
#Posterior_Head = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/COMBO_Revert_5_CombinedPosterior_Lookup2/Bias_Error_Run_CatID5/'
Posterior_Filename = 'Posterior_per_Aperture.dat'

#Posterior_Filename = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeOnly/COMBO/Bias/Bias_Combined_Posterior.dat'

#Posterior_Filename = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/COMBO/Bias/Bias_Error_Run_CatID5/1/Posterior_per_Aperture.dat'

Posterior_Label = 'r_{200}'

Cluster_Colours = ['r', 'g', 'b', 'k']

def Average_StN(Dir_head, Bias_File, nRun, Mass_Col, Mass_Err_Col):
    #Signal-To-Noise defined as the ratio of the *Input Mass* divided by half the width of the error bars
    StN = np.zeros(nRun); ErrorWidth = np.zeros(nRun)
    for i in range(1, nRun+1):
        Input = np.genfromtxt(Dir_head+str(i)+'/'+Bias_File)
        ErrorWidth[i-1] = (0.5*(Input[Mass_Error_Cols[0]]+Input[Mass_Error_Cols[1]]))
        StN[i-1] = Input[Mass_Col]/ErrorWidth[i-1]
    
    ##Account for the fact that when the StN is a NaN, the posterior has returned a value consistent with zero such that the prior dominates.
    ##This can be switched off at will
    print 'Signal-to-Noise accounts for unconstrained data, removing',  np.sum(np.isnan(StN)), ' NaNs...'
    StN[np.isnan(StN)] = 0.0

    print 'NaNs in array:', np.sum(np.isnan(StN))
    print 'StN:', stats.nanmean(StN), '+-',  stats.nanstd(StN)/np.sum(~np.isnan(StN))
    print 'Average Error width:', stats.nanmean(ErrorWidth)
    print ' '
    return stats.nanmean(StN), stats.nanstd(StN)/np.sum(~np.isnan(StN))


def get_Outlier_Removal_Indexes_PerAperture(Post, nRemove, Ap):
    ##Returns the indexes of the Runs with the nRun largest and smallest mode values
    ##Run once per aperture
    if(nRemove == 0):
        return [[],[]]
    else:
        Mode_Indexs = np.zeros(nRun[1]-nRun[0]+1)
        counter = 0
        for i in range(nRun[0], nRun[1]+1):
            Mode_Indexs[counter] = np.argmax(Post[i-1,:, Ap])
            counter+=1

        ## Get the nRemove Smallest and nRemove largest Mode_Indexs values
        Index_Sort = np.argsort(Mode_Indexs) #Returns the indexs (in the original array) of the smallest to largest
        
        Outliers = np.array((Index_Sort[:nRemove],Index_Sort[-nRemove:]))
        return Outliers+1

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

#Read in Posteriors


#Posterior is nRun, AlphaGrid, Aperture
##Originally nRun[0], nRun[1]+1
for i in range(0, nRun[1]):
    PosteriorIn = np.genfromtxt(Posterior_Head+str(i+1)+'/'+Posterior_Filename)
    if(i==0):
        #Construct Posteriors
        Posterior = np.zeros((nRun[1]+1, PosteriorIn.shape[0], PosteriorIn.shape[1]))
    Posterior[i,:,:] = np.copy(PosteriorIn)
    

    #Grid, Aperture
Combined_Posterior = np.zeros(( Posterior.shape[1], Posterior.shape[2]-1))
Combined_Posterior[:,:] = -100.

#Loop over all apertures
for j in range(1,len(Do_Aperture)+1):
    #Get Outliers for this aperture

    nIgnore_Outliers = get_Outlier_Removal_Indexes_PerAperture(Posterior, Remove_Outliers, j)
    print 'Removing Outliers:', nIgnore_Outliers, ' for aperture:', j

    for i in range(nRun[0], nRun[1]+1):
        if(i in nIgnore):
            print 'Skipping Run:', i
            continue
        if(i in nIgnore_Outliers):
            print 'Skipping Run (outlier):', i
            continue
                                
        
#Add in Combined Posterior one at a time
    
        Combined_Posterior[:,j-1] = np.copy(Combine(Combined_Posterior[:,j-1], Posterior[i-1,:,j]))
        
        if(Do_Aperture[j-1] == 0):
            continue
        ax.plot(Posterior[i-1,:,0], Posterior[i-1,:,j], label = Cluster_Labels[j-1]+'/'+str(i), linewidth = lwidth)

        ax.plot((Input_Value[j-1],Input_Value[j-1]), ax.get_ylim(), linestyle = '--', linewidth = lwidth)

#    ax.legend(loc = 2)
    Combined_Posterior = Renorm_lnLikelihood(Combined_Posterior)
    Combined_Posterior = np.exp(Combined_Posterior)

    ax_C.plot(Posterior[0,:,0],Combined_Posterior[:,j-1], label = Cluster_Labels[j-1], color = Cluster_Colours[j-1])
    ax_C.plot((Input_Value[j-1],Input_Value[j-1]), ax_C.get_ylim(), linestyle = '--', linewidth = lwidth, color = Cluster_Colours[j-1])

#ax.set_xlim((0.6, 2.0))
ax_C.set_xlim(ax.get_xlim())
#ax.set_xlim((np.amin(Input_Value)- 0.2, np.amax(Input_Value)+ 0.3))
ax_C.set_xlabel(r'$'+Posterior_Label+'$', fontsize = fsize)
ax.set_ylabel(r'$p('+Posterior_Label+')$', fontsize = fsize)
ax_C.set_ylabel(r'$p('+Posterior_Label+')$', fontsize = fsize)
#ax.set_ylabel('Posterior', fontsize = fsize)


pl.legend(loc = 1)#1:UR, 2:UL, 9:UC
#pl.show()

if(getPDF):
    output_name = 'Combined_Posterior_Plot.pdf'
    pl.savefig(output_name,format = 'pdf',bbox_inches = 'tight', pad_inches = 0.6)
    print 'Produced plot output to '+output_name
else:
    pl.show()
