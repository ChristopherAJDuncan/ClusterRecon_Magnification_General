import numpy as np
import pylab as pl

Chain_Input = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/Overlap_Bias_2ClusterMCMC/r200+Centroid/r200_1.2/Bias_Error_Run_CatID45/1/Group1_MCMC_CombinedChain.dat'
Parameter_Col = 8
nBurnin = 0

Chain = np.genfromtxt(Chain_Input)

f = pl.figure()
ax = f.add_subplot(1,1,1)

Parameter_Chain =  np.copy(Chain[nBurnin:,Parameter_Col])

ax.hist(Parameter_Chain[::1], 20)

pl.show()
