# Produces contours for fractional bias for three analysis cases:
#1: WL analysis, WL Sims
#2: WL Analysis, SL Sims
#3: SL Analysis, SL Sims

import numpy as np
import pylab as pl

ParamType = 0 #0:Mass, 1:Alpha

Plot_Ex = [1,1,1] #WLWL, WLSL, SLSL

Offset = 0.08e14 #(in Mass)

WLWLFilename = 'BiasesAndErrorVariance.dat'
WLWLHead = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Results/3Jul2014/SizeMag/SingleCluster_Bias/STAGES+COMBO/WL_Unbiased/'

WLSLFilename = 'BiasesAndErrorVariance.dat'
WLSLHead = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Results/3Jul2014/SizeMag/SingleCluster_Bias/STAGES+COMBO/SLWL_Bias/'
#'/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Results/29Jun2014/Mocks/SLWL_Comparison/StrongLensing_WeakLensingSims/'

SLSLFilename = 'BiasesAndErrorVariance.dat'
SLSLHead = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Results/3Jul2014/SizeMag/SingleCluster_Bias/Experiment_Comparison/SL_Unbiased/'


ExColor = ['r','b','g']
Marks = ['o','x','s']

Input_Value = [0.4, 0.8, 1.2, 1.6]
Cluster_Head = 'r200_'


if(ParamType ==0):
    Recovered_Param_Col = 5
    Param_Bias_Col = 6
    Param_Error_Cols = [7,8]

    Param_Label = 'Input Virial Mass'

elif(ParamType ==1):
    Recovered_Param_Col = 1
    Param_Bias_Col = 2
    Param_Error_Cols = [3,4]

    Param_Label = '$r_{200}$'

f = pl.figure()
ax = f.add_subplot(111)

def Return_Fractional_Bias(Input, Row, Bias_Col, Recovered_Col, Error_Cols):

    Input_Mass = Input[Recovered_Col] - Input[Bias_Col]
    FracBias = Input[Bias_Col]/Input_Mass
    FracError = [[Input[Error_Cols[1]]/Input_Mass],[Input[Error_Cols[0]]/Input_Mass]]

    '''
    Input_Mass = Input[Row,Recovered_Col] + Input[Row,Bias_Col]
    FracBias = Input[Row,Bias_Col]/Input_Mass
    FracError = [[Input[Row,Error_Cols[1]]/Input_Mass],[Input[Row,Error_Cols[0]]/Input_Mass]]
    '''
    return Input_Mass, FracBias, FracError

for i in range(0, len(Input_Value)):
    if Plot_Ex[0] == 1:
        WLWL = np.genfromtxt(WLWLHead+Cluster_Head+str(Input_Value[i])+'/'+WLWLFilename)
    if Plot_Ex[1] == 1:   
        WLSL = np.genfromtxt(WLSLHead+Cluster_Head+str(Input_Value[i])+'/'+WLSLFilename)
    if Plot_Ex[2] == 1: 
        SLSL = np.genfromtxt(SLSLHead+Cluster_Head+str(Input_Value[i])+'/'+SLSLFilename)

    #WL#
    if(i == 0):
        if Plot_Ex[0] == 1:
            InputMass, Frac, FracError = Return_Fractional_Bias(WLWL, 0, Param_Bias_Col, Recovered_Param_Col, Param_Error_Cols)
            ax.errorbar(InputMass-Offset, Frac, yerr = FracError, linewidth = 1.5, color =  ExColor[0], marker = Marks[0], label = 'WL')

    #WLSL#
        InputMass, Frac, FracError = Return_Fractional_Bias(WLSL, 0, Param_Bias_Col, Recovered_Param_Col, Param_Error_Cols)
        ax.errorbar(InputMass+Offset, Frac, yerr = FracError, linewidth = 1.5, color =  ExColor[1], marker = Marks[1], label = 'ML Mocks, WL Pipeline')

    #SLSL#
        InputMass, Frac, FracError = Return_Fractional_Bias(SLSL, 0, Param_Bias_Col, Recovered_Param_Col, Param_Error_Cols)
        ax.errorbar(InputMass, Frac, yerr = FracError, linewidth = 1.5, color =  ExColor[2], marker = Marks[2], label = 'ML')

    else:
        if Plot_Ex[0] == 1:
            InputMass, Frac, FracError = Return_Fractional_Bias(WLWL, 0, Param_Bias_Col, Recovered_Param_Col, Param_Error_Cols)
            ax.errorbar(InputMass-Offset, Frac, yerr = FracError, linewidth = 1.5, color =  ExColor[0], marker = Marks[0])

    #WLSL#
        InputMass, Frac, FracError = Return_Fractional_Bias(WLSL, 0, Param_Bias_Col, Recovered_Param_Col, Param_Error_Cols)
        ax.errorbar(InputMass+Offset, Frac, yerr = FracError, linewidth = 1.5, color =  ExColor[1], marker = Marks[1])

    #SLSL#
        InputMass, Frac, FracError = Return_Fractional_Bias(SLSL, 0, Param_Bias_Col, Recovered_Param_Col, Param_Error_Cols)
        ax.errorbar(InputMass, Frac, yerr = FracError, linewidth = 1.5, color =  ExColor[2], marker = Marks[2])

        
ax.set_xlim(ax.get_xlim()[0]-0.1*(np.diff(ax.get_xlim())), ax.get_xlim()[1]+0.1*(np.diff(ax.get_xlim())))

ax.plot( ax.get_xlim(), (0.,0.), linestyle = '--', color = 'k')

#ax.set_yscale('log')

ax.legend(loc = 2)
ax.set_xlabel(r'Input Virial Mass $\left(M_{200}\; \left[\frac{\rm M_\odot}{h}\right]\right)$', fontsize = 16)
ax.set_ylabel('Fractional Bias in Recovered Mass', fontsize = 16)

#Legend#
pl.show()
