import pylab as pl
import sys
import numpy as np
import os

def Renormalise(Lim, Num):
    Sum = 0.0
    for i in range(0,Num.shape[0]):
        Sum += (Lim[i,1]-Lim[i,0])*Num[i]
    return Sum


def Mean_of_Histogram(Lim, Num):
    Sum = 0.; Sum_NX = 0.
    for i in range(0, Num.shape[0]):
        Sum_NX += Num[i]*0.5*(Lim[i,1]+Lim[i,0])
        Sum += Num[i]

    return Sum_NX/Sum

do_Renormalise = True

Input_Filename = sys.argv[1]
Output_Directory = os.path.dirname(Input_Filename)
Output_Filename = Output_Directory+'/'+os.path.basename(Input_Filename).split('.')[0]+'.pdf'

#Read in Bar plot 1#
Input = np.genfromtxt(Input_Filename)
Limits = np.zeros( (Input.shape[0],2) )
Number = np.zeros( Input.shape[0] )
Limits[:,0] = Input[:,0]
Limits[:,1] = Input[:,1]
Number =  1.0*Input[:,2]


Normalise = 1.0
if(do_Renormalise):
    Normalise = Renormalise(Limits, Number)
    #np.divide(Number,  Normalise)
    Number /= Normalise


xTicks = []
for i in range(0, Limits.shape[0]):
    xTicks.append((Limits[i,0],Limits[i,1]))

fig, ax = pl.subplots()
rects = ax.bar(Limits[:,0], Number, np.subtract(Limits[:,1], Limits[:,0]))

Global_Mean = Mean_of_Histogram(Limits, Number)

ax.set_title(Input_Filename)
if(len(sys.argv) == 3): #Plot 2#
    Input_Filename2 = sys.argv[2]
    Input2 = np.genfromtxt(Input_Filename2)
    Limits2 = np.zeros( (Input2.shape[0],2) )
    Number2 = np.zeros( Input2.shape[0] )
    Limits2[:,0] = Input2[:,0]
    Limits2[:,1] = Input2[:,1]
    Number2 =  Input2[:,2]
    # )
    Normalise = 1.0
    if(do_Renormalise):
        Normalise = Renormalise(Limits2, Number2)
        Number2 /= Normalise

    
    rects2 = ax.bar(Limits2[:,0], Number2, np.subtract(Limits2[:,1], Limits2[:,0]), color =  'green')
    ax.set_title(Input_Filename2)
    Output_Filename = Output_Directory+'/'+os.path.basename(Input_Filename2).split('.')[0]+'.pdf'

    ##Add Shift Information##
    fig.text(0.9, 0.8, 'Shift = '+str(Mean_of_Histogram(Limits2, Number2)), ha = 'right')
    fig.text(0.9, 0.75, 'Kappa = '+str((Mean_of_Histogram(Limits2, Number2)/Global_Mean)-1.), ha = 'right')

if(do_Renormalise):
    ax.set_ylabel('Renormalised Number')
else:
    ax.set_ylabel('Number')
ax.set_xlabel('Size')
#ax.set_title(os.path.basename(Input_Filename))


pl.savefig(Output_Filename,format = 'pdf',bbox_inches = 'tight')
print 'Produced plot output to '+Output_Filename
os.system('gv '+Output_Filename+' &')
