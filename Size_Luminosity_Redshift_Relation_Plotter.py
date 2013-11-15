import numpy as np
import pylab as pl
import os
import sys


Input_Filename_SLR = sys.argv[1]
Output_Directory = os.path.dirname(Input_Filename_SLR)

zLabel = sys.argv[2]

SLR_Input = np.genfromtxt(Input_Filename_SLR)

SLR_M = SLR_Input[1:,0]
SLR_Z = SLR_Input[0,1:]
SLR = SLR_Input[1:SLR_Input.shape[0]-1,1:SLR_Input.shape[1]-1]


'''
for i in range(0, SLR.shape[1]):
    SLR[:,i] = SLR[:,i]/(1+SLR_Z[i])
'''
#SLR /= (1+SLR_Z)

ax = pl.subplot(1,1,1)
im = pl.pcolormesh(SLR_Z, SLR_M, SLR)
cbar = pl.colorbar()

ax.set_xlabel(r"Redshift")
ax.set_ylabel(r"Magnitude")
cbar.ax.set_ylabel(zLabel)

#f = pl.figure()
#ax = f.add_subplot(111)
#im = ax.imshow(SLR, interpolation = 'None', aspect = 'equal')#, extent = (SLR_M[0], SLR_M[SLR_M.shape[0]-1], SLR_Z[0], SLR_Z[SLR_Z.shape[0]-1]), aspect = 'equal', interpolation = None)
#f.colorbar(im)


Output_Filename = Output_Directory+'/'+os.path.basename(Input_Filename_SLR).split('.')[0]+'.pdf'
#'Size_Magnitude_Redshift.pdf'
pl.savefig(Output_Filename,format = 'pdf',bbox_inches = 'tight')
print 'Produced plot output to '+Output_Filename
os.system('evince '+Output_Filename+' &')
           
