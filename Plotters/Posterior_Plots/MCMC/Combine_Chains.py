##Program to read in N independent MCMC chains and combine them into one long chain
import numpy as np


Directory = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/25Aug2015/MCMC/4Cluster/SM+M/Mass-Only/'
#'/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/Data/STAGES/Calibrated_6ClusterMCMC/'
#'/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/Data/STAGES/Uncalibrated_SingleAperture_MCMC_4Cluster/'
#'/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/Data/STAGES/Calibrated_4ClusterMCMC/SizeMag+Mag/'

#Mag-Only Data Result'/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SimultaneousFitting/MCMC/r200+Centroid/'
File_Head = 'Group1_MCMC_Chain_'
File_Tail = '.dat'
nChain = 6
nBurnin = 0 #1000

Output_File = 'Group1_MCMC_CombinedChain.dat'

##Read in first seperately
Combined = np.genfromtxt(Directory+File_Head+str(1)+File_Tail)[nBurnin:]

print 'Chain 1 contains:', Combined.shape[0]

for i in range(2,nChain+1):
    Combined = np.append(Combined, np.genfromtxt(Directory+File_Head+str(i)+File_Tail)[nBurnin:], axis = 0)


if Combined.shape[0] == 0:
    print 'Error - Combined chain has no elements!!'

print 'Combined chain has ', Combined.shape[0], ' links'

np.savetxt(Directory+Output_File, Combined)
