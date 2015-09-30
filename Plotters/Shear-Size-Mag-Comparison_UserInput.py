import numpy as np
import pylab as pl

Plot_Cluster = [1,1,1,1]
Cl_Marker = ['x','o','s','^']
Cl_Color =  ['r', 'g', 'b', 'k']
Cl_Label = ['A901a', 'A901b', 'A902', 'SW']

###Shear Results
#Taken from Table 1 of Heymans 2008, one halo fits
Shear_r200 = [1.194, 1.180, 0.799, 0.894]
Shear_r200_Stat_Error = [[0.086, 0.101], [0.088, 0.104], [0.122,0.125], [0.106, 0.117]]

'''
### Magnification Results
Mag_Dir = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Data/SizeMag+Mag/SingleClusterFit/RRG+GALFIT/SM+M/'
Mag_Filename = 'Mass_Estimates.dat'

## Input Column Information
Rec_Col = 2
Err_Col = [3,4]

Mag = np.genfromtxt(Mag_Dir+Mag_Filename)
'''

Mag_r200 = [1.107, 1.258, 0.631, 0.782]
Mag_r200_Stat_Error = [[0.072, 0.119],[0.187, 0.235], [0.095, 0.182], [0.108, 0.153]]



fig_width_pt = 246.0*2 # Get this from LaTex using \the\columnwidth
inches_per_pt = 1.0/72.27
golden_mean = 2.5*(np.sqrt(5)-1.0)/2.0
fig_width  = fig_width_pt*inches_per_pt # width in inches
fig_height = fig_width*golden_mean # height in inches
fig_size = [fig_width, fig_height]
params = {'axes.labelsize':16,
          'text.fontsize':16,
          'legend.fontsize':14,
          'xtick.labelsize':14,
          'ytick.labelsize':14,
          #'figure.figsize':fig_size,
          'font.family': 'serif'}

pl.rcParams.update(params)
#pl.clf()
f = pl.figure()
#ax = f.add_subplot(2,1,1)                                    
ax = pl.subplot2grid((3,3),(0,0), colspan = 3, rowspan = 2)

    #Plot r200 Mag (x) vs r200 Shear (y)

for i in range(0, len(Plot_Cluster)):
    if(Plot_Cluster[i] == 0):
        continue
    
    ax.errorbar(Mag_r200[i], Shear_r200[i], xerr = [[Mag_r200_Stat_Error[i][1]],[Mag_r200_Stat_Error[i][0]]], yerr = [[Shear_r200_Stat_Error[i][1]],[Shear_r200_Stat_Error[i][0]]], marker = Cl_Marker[i], color = Cl_Color[i], label = Cl_Label[i])#, xerr = (Mag[i,Err_Col[1]],Mag[i,Err_Col[0]]))#, yerr = (Shear_r200_Stat_Error[i][1], Shear_r200_Stat_Error[i][0]), 

yeqx = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 50)
ax.plot(yeqx, yeqx, linestyle = '--', color = 'k')

ax.set_xlim(yeqx[0], yeqx[-1]);ax.set_ylim(yeqx[0], yeqx[-1])

ax.set_xlabel(r'$r_{200}\; [h^{-1} {\rm Mpc}] \; ({\rm Magnification})$')
ax.set_ylabel(r'$r_{200}\; [h^{-1} {\rm Mpc}] \; ({\rm Shear})$')

ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

ax.legend(loc = 4, numpoints = 1)


#ax = f.add_subplot(2,1,2)
ax = pl.subplot2grid((3,3),(2,0), colspan = 3, rowspan = 1)

xTicks = range(1, len(Cl_Label)+1)

for C in range(len(Cl_Label)):
    ax.plot(xTicks[C], sum(Mag_r200_Stat_Error[C])/sum(Shear_r200_Stat_Error[C]), 'kx', markersize = 8)
ax.set_xticklabels(Cl_Label, rotation = 45.)
ax.set_xticks(ax.get_xticks()[::2])
ax.set_ylabel(r'$\langle \sigma_{\rm Mag} \rangle / \langle \sigma_{\rm Shear} \rangle $')
ax.margins(0.1)

pl.show()
