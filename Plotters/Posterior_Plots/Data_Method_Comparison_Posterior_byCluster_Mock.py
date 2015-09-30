#Produces a series of plots, one for each aperture, which compares posteriors

import numpy as np
import pylab as pl

Dir_Head = '/disk2/cajd/STAGES_Bayesian/Bayesian_DM_Profile_Constraints_Output/BiasCheck/FaintCut/KnownZ_PriorReEval/M/'
Directories = [Dir_Head+'m_26/STAGES+COMBO/', Dir_Head+'m_26.8/STAGES+COMBO/', Dir_Head+'m_27.5/STAGES+COMBO/']
Experiment_Label = ["m<26","m<26.8", "m<27.5"]
r200 = [0.4, 0.8, 1.2, 1.6]
#Directories = ['/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Data/SizeMag+Mag/SingleClusterFit/RRG+GALFIT/M/','/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Data/SizeMag+Mag/SingleClusterFit/RRG+GALFIT/SM/', '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Data/SizeMag+Mag/SingleClusterFit/RRG+GALFIT/SM+M/']
#Experiment_Label = ['Magnitude-Only', 'Size-Mag, SNR 30+', 'Size-Mag+Mag, SNR 30+']
Aperture_Label = [' ', 'A901a-like']#[' ','A901a', 'A901b', 'A902', 'CB1', 'SWa', 'SWb'] #Always leave first empty, or discardable
Colours = ['k', 'b', 'r', 'c']
Marker = ['o','x','v', 's']

#For Plot_Type == 1
Posterior_File = 'Posterior_per_Aperture.dat'
#For Plot_Type == 2
Estimator_File = 'BiasesAndErrorVariance.dat'
Alpha_Col = 2
Alpha_Err_Col = [3,4]


Plot_Type = 2 #Posterior Plot. 2: Point Plot (only Plot 2 coded for now [5th June])


if(Plot_Type == 1):

    ## Read in Posterior File
    Posteriors = []
    for i in range(len(Directories)):
        Posteriors.append(np.genfromtxt(Directories[i]+Posterior_File))
    
    f =pl.figure()
    for i in range(1,Posteriors[0].shape[1]):
        ax = f.add_subplot(2,3,i) ##Manually Change
        
        ax.set_title(Aperture_Label[i-1])
        if(i==1 or i==4):
            ax.set_ylabel(r'$p(r_{200})$', fontsize = 14)
        if(i==4 or i==5 or i==6):
            ax.set_xlabel(r'$r_{200}$', fontsize = 14)
            
        ax.set_xlim((0.,2.0))
    
        for ex in range(len(Posteriors)):
            ax.plot(Posteriors[ex][:,0], Posteriors[ex][:,i], label = Experiment_Label[ex])
            
        if(i==2):
            ax.legend(prop = {'size':11})

if(Plot_Type == 2):
    f =pl.figure()
    ax = f.add_subplot(111)

    Offset = 0.1/float(len(Directories))

    Alpha = np.zeros((len(r200), len(Directories))); Err = np.zeros((len(r200), len(Directories), 2))
    for i in range(len(Directories)):
        Offset_Index = i - int(len(Directories)/2.)

        for al in range(len(r200)):
            Input = np.genfromtxt(Directories[i]+'r200_'+str(r200[al])+'/'+Estimator_File)

            print 'Input:', Input.shape, Input
            Alpha[al,i] = Input[Alpha_Col] ##Assume only one cluster

        
            if(len(Alpha_Err_Col) == 1):
                Err[al,i,0] = Input[Alpha_Err_Col]; Err[al,i,1] = Input[Alpha_Err_Col]
            elif(len(Alpha_Err_Col) ==2):
                Err[al,i,1] = Input[Alpha_Err_Col[0]]; Err[al,i,0] = Input[Alpha_Err_Col[1]]
            else:
                print 'PROBLEM'

            ##Plot
            if(al == 0):
                Label = Experiment_Label[i]
            else:
                Label = None
            ax.errorbar(np.array(r200[al])+Offset*Offset_Index, Alpha[al,i], yerr = Err[al,i,:].reshape(2,1), color = Colours[i], label = Label, marker = Marker[i])


    ##Plot Horizontal
    ax.plot(ax.get_xlim(), (0., 0.), '--')

    ax.set_ylim((ax.get_ylim()[0], 1.2*ax.get_ylim()[1]))

    ax.set_ylabel(r'$r^{\rm rec}_{200}-r^{\rm in}_{200}$ [Mpc/h]')
    ax.set_xlabel(r'$r^{\rm in}_{200}$ [Mpc/h]')
    
    ax.legend(loc = 3)

pl.show()
