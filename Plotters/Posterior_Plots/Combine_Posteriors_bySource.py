'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
Author: cajd (@roe.ac.uk)
Date: 20 Mar 2015
Purpose: Combines individual PDFs as independent measures. Plots output.
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

import numpy as np
import pylab as pl

#FGrid = 'Likelihood_RVirGrid_all_galaxies_corrected_shear'
#FLike = 'Likelihood_perGalaxy_all_galaxies_corrected_shear'

FLike = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/MagCutsTesting_Data/SingleClusterFit/4Cluster/RRG+STAGES/SM+M/M_less_27p5/Run2/SExMasterMag/Aperture_3__Posterior_per_Source_SM.dat'
#Grid = np.genfromtxt(FGrid)
Like = np.genfromtxt(FLike)

nGalMin = 0
nGalMax = Like.shape[1]

debug = 0

lnCombine = True ##Must be true (alternative not coded). Taken into account later


Grid = Like[:,0]
CLike  = Like[:,1]
Like = Like[:,2:] #First Grid, 2nd Total Likelihood


print 'Likelihood read in has ', Like.shape[1],' contributions'



## Use this to remove certain galaxies form sample (if their likelihood is suspect.) Be careful not to pick galaxies to get the result you want (Confirmation Bias)
##Note to Juliet: I've probably removed some here that don't need to be in my serch for the problem one, (i.e. 248-254)
remove_index = [] #A901a[62, 832, 982]

 #Use this to plot problem (removed) posteriors
if(debug):
    f = pl.figure()
    ax = f.add_subplot(111)
    
    for i in range(len(remove_index)):
        print 'Removed Gal:', remove_index[i], ':',  Like[:,remove_index[i]]
        ax.plot(Grid, Like[:,remove_index[i]], label = str(i))
    #ax.scatter(Grid, CL)
    ax.set_title('Removed Sources')
    ax.legend()
    raw_input('Check')

Like = np.delete(Like, remove_index, 1)



Like = Like[:,nGalMin:nGalMax+1]

print 'Like Readin:', Like

if(debug):
    # use this to plot individual for debugging
    f = pl.figure()
    ax = f.add_subplot(111)
    
    for i in range(Like.shape[1]):
        ax.plot(Grid, Like[:,i], label = str(i))
    #ax.scatter(Grid, CL)
    ax.legend()


if(lnCombine or True): #"or True" takes into accont alternative not coded
    T = Like <= 0.
    print T.sum(),  'contibutions to the likelihood are negative or zero', T.shape, Like.shape

    print 'Like neg or zero:', Like[:,Like<=0]

    Like[T] = -100.
    Like[~T] = np.log(Like[~T]) + 1.
    #--Like is now ln(Likelihood). +1 as a some general form of renormalisation to try and avoid numerical error.

    print 'lnLikelihood:', Like

    CL = Like.sum(axis = 1)

    print 'CL sum:', CL
    
    CL = CL - np.nanmax(CL)

    CL = np.exp(CL)


##Checks
print CL.shape, np.isnan(CL).sum()

print 'CL:', CL

###------- Plotting
f = pl.figure()
ax = f.add_subplot(111)

print 'CLike:', CLike
CLike = np.exp(CLike- np.nanmax(CLike))
print CLike

ax.plot(Grid, CL)
#ax.scatter(Grid, CL) # Use this to plot sampling points.
#ax.plot(np.ones(2)*(1.2/0.7), ax.get_ylim(), linestyle = '--', color = 'g')
ax.plot(Grid, CLike, label = 'Total Likelihood')

ax.set_xlabel(r'Virial Radius [Mpc]')
ax.set_ylabel(r'Likelihood')

pl.show()


