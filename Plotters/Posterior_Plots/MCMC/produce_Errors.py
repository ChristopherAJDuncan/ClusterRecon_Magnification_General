## Code to read in an MCMC chain, produce a finely binned distribution, find the mode and 1 sigma variance of each side of the mode
import numpy as np
import pylab as pl
from scipy.stats import gaussian_kde

Chain_Dir = '/disk1/ps1/cajd/STAGES_Bayesian/Bayesian_DM_Profile_Constraints_Output/3Sept2015/MCMC/4Cluster/RRG_RRG/NoSNRCut_RelaxedCore/'
Chain_File = 'Group1_MCMC_CombinedChain.dat'
doIndex = [1,5,9,13] #If -1, then do all indexes


Chain = np.genfromtxt(Chain_Dir+Chain_File)


if(doIndex <= 0):
    idoIndex = range(1,Chain.shape[0])
else:
    idoIndex = doIndex


def get_Mode_and_Error(Chain, Index):

    #Produce Histogram (Fine)
    '''
    H,ed = np.histogram(Chain[:,Index], bins = 128, normed = True)
    modeIndex = np.argmax(H)
    mode = ed[np.argmax(H)]; modeval = H[np.argmax(H)]
    print 'Mode at:', mode
    X = ed[:-1] + 0.5*np.diff(ed)
    '''

    ## Get KDE Smoothed version
    density = gaussian_kde(Chain[:,Index])
    density.covariance_factor = lambda : Chain[:,Index].std()*1.0
    density._compute_covariance()

    start = np.nanmin(Chain[:,Index]); stop = np.nanmax(Chain[:,Index]); step = (stop-start)/1000.
    X = np.arange(start, stop, step)
    H = density(X)

    modeIndex = np.argmax(H)
    mode = X[modeIndex]; modeval = H[modeIndex]

    '''
    f = pl.figure()
    ax = f.add_subplot(111)

    ax.plot(X,H)
    pl.show()
    '''

    #Produce cumulative histogram (No Bin)
    '''
    N = 1.0*Chain.shape[0]
    X = np.sort(Chain[:,Index])

    print X.shape[0], N, Chain.shape
    
    C = np.array(range(int(N)))/N

    print 'C,Y:', C.shape, N, X.shape
    '''
    ## Produce Cumulative, from histogram
    dx = X[1]-X[0]
    #dx = ed[1]-ed[0]
    C = np.cumsum(H)*dx

    '''
    f = pl.figure()
    ax = f.add_subplot(111)

    ax.plot(X,C)
    pl.show()
    '''

    C_at_mu = np.interp(mode, X, C)


    ##Seperate Distribution about mode
    DivideArea = np.nanmax(C)-C_at_mu


    SD = [1.]*2
    ## Calculate the StD below the mode point, by looking for 32 percent (1-0.68 to account for the fact that we are summing from zero here, not from the mode) of the cumulative up to the mode point
    DivideArea = C_at_mu
    SD[0] = mode - np.interp(0.32*DivideArea, C[:modeIndex], X[:modeIndex])

    ## Calculate the StD above the mode point, by lookig for the 68% partto the right of the mode
    DivideArea = np.nanmax(C)-C_at_mu
    SD[1] = np.interp(0.68*DivideArea, C[modeIndex:]-C_at_mu, X[modeIndex:]) - mode

    #print 'Mode, -+ Var, -+ StD:', mode, np.power(SD,2.), SD

    return mode, SD

CMode = [0.]*len(doIndex); CStD = [0.]*len(CMode)
for II,I in enumerate(doIndex):
    CMode[II], CStD[II] = get_Mode_and_Error(Chain, I)
    print 'Cluster:', II, ' has mode, StD:', CMode[II], CStD[II]
