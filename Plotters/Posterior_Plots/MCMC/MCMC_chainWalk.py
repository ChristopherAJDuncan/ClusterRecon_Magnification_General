## Plots the variation of the ln-likelihood (always assumed to be the last column of the input) as a function of chanin link. Useful to see the burn-in period.
import numpy as np
import pylab as pl
import os

Directory = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/Data/STAGES/Calibrated_6ClusterMCMC_SSNRCalib2/MassOrdering_PzLookup/'
File_Prefix = 'Group1_MCMC_Chain_'
File_Suffix = '.dat'
nChain = 8
minLink = 0000
maxLink = 100000

f = pl.figure(0)
ax = f.add_subplot(111)

for C in range(nChain):
    Chain = np.genfromtxt(os.path.join(Directory, File_Prefix)+str(C+1)+File_Suffix)

    Chain = Chain[np.logical_and(Chain[:,0]<= maxLink, Chain[:,0]>=minLink),:]

    ax.plot(Chain[:,0], Chain[:,-1], label = 'Chain '+str(C+1))

ax.legend(loc = 4)
ax.set_xlabel('Chain Link')
ax.set_ylabel(r'$\ln L$')

pl.show()
