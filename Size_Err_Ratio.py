import numpy as np
import pylab as pl


Kappa_In = np.genfromtxt('Plotters/Output/Average_Size_KappaEst_in_RA_Dec_Grid.dat')
Kappa = Kappa_In[1:Kappa_In.shape[0]-1, 1:Kappa_In.shape[1]-1]

KappaErr_In = np.genfromtxt('Plotters/Output/Average_Size_KappaError_in_RA_Dec_Grid.dat')
KappaErr = KappaErr_In[1:KappaErr_In.shape[0]-1, 1:KappaErr_In.shape[1]-1]

RatioK = np.zeros(Kappa.shape)

for i in range(1,Kappa.shape[0]-1):
    for j in range(1,Kappa.shape[1]-1): 
        if(KappaErr[i,j] != 0.):
            RatioK[i,j] = Kappa[i,j]/KappaErr[i,j]
        else:
            RatioK[i,j] = 0.

f = pl.figure()
ax = f.add_subplot(1,1,1)
im = ax.imshow(RatioK)
f.colorbar(im)

ax.set_title(r'$\kappa_{Size}/\sigma_{\kappa}$')
pl.show()


Size_In = np.genfromtxt('Plotters/Output/Average_Size_in_RA_Dec_Grid.dat')
Size = Size_In[1:Size_In.shape[0]-1, 1:Size_In.shape[1]-1]

SizeErr_In = np.genfromtxt('Plotters/Output/Average_Size_Error_in_RA_Dec_Grid.dat')
SizeErr = SizeErr_In[1:SizeErr_In.shape[0]-1, 1:SizeErr_In.shape[1]-1]

Ratio = np.zeros(Size.shape)

for i in range(1,Size.shape[0]-1):
    for j in range(1,Size.shape[1]-1): 
        if(SizeErr[i,j] != 0.):
            Ratio[i,j] = Size[i,j]/SizeErr[i,j]
        else:
            Ratio[i,j] = 0.

f = pl.figure()
ax = f.add_subplot(1,1,1)
im = ax.imshow(Ratio)
f.colorbar(im)

ax.set_title(r'$R/\sigma_{R}$')
pl.show()


RatioErr = np.zeros(SizeErr.shape)
for i in range(1,Size.shape[0]-1):
    for j in range(1,Size.shape[1]-1): 
        if(SizeErr[i,j] != 0.):
            RatioErr[i,j] = KappaErr[i,j]/SizeErr[i,j]
        else:
            RatioErr[i,j] = 0.

f = pl.figure()
ax = f.add_subplot(1,1,1)
im = ax.imshow(RatioErr)
f.colorbar(im)

ax.set_title(r'$\sigma_{\kappa}/\sigma_{R}$')
pl.show()

