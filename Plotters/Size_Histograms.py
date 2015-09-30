import numpy as np
import pylab as pl
import os

Catalogue_Dir = '../Catalogues/'
Catalogue_Flnm = 'STAGES_shear.cat'
Output_Dir = 'Output/'
getPDF = True

Cat = np.genfromtxt(Catalogue_Dir+Catalogue_Flnm)

FWHM = np.zeros(Cat.shape[0])
FWHM = np.copy(Cat[:,9])
FR = np.zeros(Cat.shape[0])
FR = np.copy(Cat[:,10])
KSB = np.zeros(Cat.shape[0])
KSB = np.copy(Cat[:,11])


##Note that Cols 9-FWHM, 10-FR, 11-KSB##
f = pl.figure()
ax_FWHM = f.add_subplot(3,1,1)
ax_FR = f.add_subplot(3,1,2)
ax_KSB = f.add_subplot(3,1,3)

ax_FR.set_ylabel('N')
ax_KSB.set_xlabel('Size (Pixel)')

FWHM_range = (0.,70.)
ax_FWHM.hist(FWHM, bins=100, range=FWHM_range, normed=False, weights=None,\
       cumulative=False, bottom=None, histtype='bar', align='mid',\
       orientation='vertical', rwidth=None, log=False,\
       color=None, label='FWHM')
ax_FWHM.legend()

FR_range = (0.,70.)
ax_FR.hist(FR, bins=100, range=FR_range, normed=False, weights=None,\
       cumulative=False, bottom=None, histtype='bar', align='mid',\
       orientation='vertical', rwidth=None, log=False,\
       color=None, label='FR')
ax_FR.legend()

KSB_range = (0.,70.)
ax_KSB.hist(KSB, bins=100, range=KSB_range, normed=False, weights=None,\
       cumulative=False, bottom=None, histtype='bar', align='mid',\
       orientation='vertical', rwidth=None, log=False,\
       color=None, label='KSB')
ax_KSB.legend()

print 'Maximum/Minimum:'
print 'FWHM:', np.amax(FWHM), np.amin(FWHM)
print 'FR:', np.amax(FR), np.amin(FR)
print 'KSB:', np.amax(KSB), np.amin(KSB)

print 'Mean of FWHM is:', np.mean(FWHM), '; Variance:', np.var(FWHM)
#print 'Nignored: Too Small       Too Large         Total'
#print 'FWHM:', len(FWHM[FWHM<FWHM[0]]), '       ', len(FWHM[FWHM>FWHM[1]]), '        ', FWHM.shape[0]
print 'Mean of FR is:', np.mean(FR), '; Variance:', np.var(FR)
print 'Mean of KSB is:', np.mean(KSB), '; Variance:', np.var(KSB)

print 
print 'Nignored: Too Small       Too Large         Total'
print 'FWHM:    ', len(FWHM[FWHM<FWHM_range[0]]), '       ', len(FWHM[FWHM>FWHM_range[1]]), '        ', FWHM.shape[0]
print 'FR:     ', len(FR[FR<FR_range[0]]), '       ', len(FR[FR>FR_range[1]]), '        ', FR.shape[0]
print 'KSB:       ', len(KSB[KSB<KSB_range[0]]), '       ', len(KSB[KSB>KSB_range[1]]
), '        ', KSB.shape[0]


if(getPDF):
    output_name = 'Catalogue_Size_Histograms.pdf'
    pl.savefig(Output_Dir+output_name,format = 'pdf',bbox_inches = 'tight')
    print 'Produced plot output to '+Output_Dir+output_name
    os.system('evince '+Output_Dir+output_name+' &')
else:
    pl.show()


