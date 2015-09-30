import numpy as np
from matplotlib import cm
import pylab as pl
from mpl_toolkits.mplot3d.axes3d import Axes3D

animate = False

Filename = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/Data-Run_4Jul/MagOnly/Data/STAGES+COMBO/MagSize_Distribution_KDE.dat'

Input = np.genfromtxt(Filename)
#Row Column# Mag is along Column

print np.amax(Input), np.amin(Input)

MagGrid = Input[0,1:]
SizeGrid = Input[1:,0]
Dist = Input[1:,1:]


for i in range(0, Dist.shape[0]):
    for j in range(Dist.shape[1]):
        if(np.isnan(Dist[i,j])):
            print 'NaN found at:', SizeGrid[i-1], MagGrid[j-1]

#Remove NaNs
Dist[np.isnan(Dist)] = 0.

print np.amax(Dist), np.amin(Dist)

##Clip data#
SizeGrid[SizeGrid > 20] = np.nan
Dist[SizeGrid>20,:] = np.nan

MagGrid , SizeGrid = np.meshgrid(MagGrid, SizeGrid)

if(animate):
    pl.ion()

fig = pl.figure()

ax = fig.add_subplot(1,1,1, projection = '3d')

surf = ax.plot_surface(SizeGrid, MagGrid, Dist, cmap = cm.Reds, rstride = 3, cstride = 3)


#ax.set_xlim3d(0., 20.)

ax.set_xlabel(r'$\theta$')
ax.set_ylabel(r'$m$')
ax.set_zlabel(r'$p(\theta,m)$')

fig.colorbar(surf, shrink = 0.5, aspect = 10)

if(animate):
    for angle in range(0, 360):
        print angle
        ax.view_init(30, angle)
        pl.draw()
else:
    ax.view_init(30, 60)
    pl.show()

