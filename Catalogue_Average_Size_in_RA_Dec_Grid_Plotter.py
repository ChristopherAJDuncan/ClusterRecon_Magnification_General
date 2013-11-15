import numpy as np
from scipy import signal
from scipy import ndimage
from scipy import interpolate
import pylab as pl
import os
import sys
import pyfits as pf
import math
import time

Input_Filename_AvSize = sys.argv[1]
Output_Directory = os.path.dirname(Input_Filename_AvSize)
Input_Filename_KappaSize = sys.argv[2]
Input_Filename_ErrorKappaSize = sys.argv[3]
Input_Filename_OccGrid =  sys.argv[4]
Input_Filename_Smoothed_OccGrid = Output_Directory + '/' +os.path.basename(Input_Filename_OccGrid).split('.')[0]+'_Smoothed.dat'

#Directory = re.match('(.*)(/)',Input_Filename).group()

#Basename is the filename, dirname is the directory part. Split s into a list of two parts - up to the . and after -> thus Output_List[1] contains the file suffix only

Output_List = os.path.basename(Input_Filename_AvSize).split('.')
Output_Filename_avSize = Output_Directory+'/'+Output_List[0]+'.ps'
'''
Output_List = os.path.basename(Input_Filename_KappaSize).split('.')
Output_Filename_kappaScatter = Output_Directory+'/'+Output_List[0]+'.ps'
'''
Output_Filename_kappaScatter = Output_Directory + '/' + 'KappaSize_Shape_Scatter.ps'

Output_List = os.path.basename(Input_Filename_OccGrid).split('.') 
Output_Filename_OccGrid = Output_Directory + '/' + Output_List[0]+'.ps'

Output_Filename_Noise_Map = Output_Directory + '/' + 'KappaNoiseMap.ps'
Output_Filename_StNoise = Output_Directory + '/' + 'Kappa_SN_Map.ps'

'''Produce Plots?'''
scatter_Plot = True
plot_Occupation_Grid = True
plot_Noise_Maps = True
plot_Kappa_Size_Map = True
plot_Signal_to_Noise = True
####################

plot_KappaSize_Contours = True
plot_Size_Contours = True

'''STAGES DM MAP'''
plotSTAGES_DM_Map = True
DM_fits_file = './Catalogues/STAGES_DMmap.fits'
DM_SN_fits_file = './Catalogues/STAGES_SN_DMmap.fits'

'''STRUCTURES'''
labelStructures = True
A901a = [149.1099,-9.9561]
A901b = [148.9889,-9.9841]
A902 = [149.1424,-10.1666]
SW = [148.9101,-10.1719]

def Smooth(Data, x, y):
    #Smooths on A Gaussian with width (sigma) = Guass_Width
    #xgrid is assumed to be on ROWS (or 1st index), y on COLS (or 2nd index)
    if(Data.shape[0] != x.shape[0]):
        print "ERROR in SMOOTHING - Shapes of x and data are not the same"
        exit()
    if(Data.shape[1] != y.shape[0]):
        print "ERROR in SMOOTHING - Shapes of x and data are not the same"
        exit()
    Gauss_Width = 0.75/60. #Converted from arcminutes to degrees
    print "Gauss Width Pixels:", Gauss_Width/(x[1]-x[0])
    print "Number of Pixels (x,y):", x.shape[0], y.shape[0]
    return ndimage.filters.gaussian_filter(Data, (Gauss_Width/(x[1]-x[0]), Gauss_Width/(y[1]-y[0])), mode='constant', cval = 0.0)


def read_in_DM_Map_ShearSTAGES(Filename):
    Map = pf.getdata(Filename)
    #Fits file is 256x256 pixels - each pixel has a scale of 7.6 arcseconds#                         
    RA_DM = np.zeros(256)
    Dec_DM = np.zeros(256)
    #Minvalue and dn have been taken from fits file header                                           
    for i in range(0,Dec_DM.shape[0]):
        RA_DM[i] = 149.3524 + i*(-0.002120962)
#RA_DM_Map[i] = 148.85 + i*(149.31-148.85)/256                                                       
        Dec_DM[i] = -10.29099 + i*(0.002120962)
    return Map, RA_DM, Dec_DM 


Input_AvSize = np.genfromtxt(Input_Filename_AvSize)
#Seperate out into RA, Dec and Average_Size. Note that binned such that RA(i) gives the lower limit to bin i, and RA(i+1) the upper limit to bin i+1 such that RA(i)<RA<=RA(i+1), also for Dec. Av_Size should be bounded by zeros in lowest and highest row/col.
###INPUT###
RA_In = Input_AvSize[1:,0]
Dec_In = Input_AvSize[0,1:]
Av_Size_In = Input_AvSize[1:,1:]

Input_Kappa = np.genfromtxt(Input_Filename_KappaSize)
KappaSize = Input_Kappa[1:Input_Kappa.shape[0]-1,1:Input_Kappa.shape[1]-1]

Input_ErrorKappa = np.genfromtxt(Input_Filename_ErrorKappaSize)
KappaSize_Noise = Input_ErrorKappa[1:Input_ErrorKappa.shape[0]-1,1:Input_ErrorKappa.shape[1]-1]

Input_OccGrid = np.genfromtxt(Input_Filename_OccGrid)
OccGrid = Input_OccGrid[1:Input_OccGrid.shape[0]-1, 1:Input_OccGrid.shape[1]-1]

Input_OccGrid_Sm = np.genfromtxt(Input_Filename_Smoothed_OccGrid)
SmOccGrid = Input_OccGrid_Sm[1:Input_OccGrid_Sm.shape[0]-1, 1:Input_OccGrid_Sm.shape[1]-1]

DM_Map, RA_DM_Map, Dec_DM_Map = read_in_DM_Map_ShearSTAGES(DM_fits_file)
DM_SN_Map, RA_Disc, Dec_Disc = read_in_DM_Map_ShearSTAGES(DM_SN_fits_file)
DM_Noise_Map = np.zeros(DM_Map.shape)
for i in range(0,DM_SN_Map.shape[0]):
    for j in range(0,DM_SN_Map.shape[1]):
        if(DM_SN_Map[i,j] != 0.):
            DM_Noise_Map[i,j]  = DM_Map[i,j]/DM_SN_Map[i,j]

Kappa_Size_SN = np.zeros(KappaSize.shape)
for i in range(0,KappaSize.shape[0]):
    for j in range(0,KappaSize.shape[1]):
        if(KappaSize_Noise[i,j] != 0.):
            Kappa_Size_SN[i,j]  = KappaSize[i,j]/KappaSize_Noise[i,j]

#DM_Noise_Map = np.divide(DM_Map,DM_SN_Map)
####END INPUT####

#Rescale RA/Dec such that they contain the mean value for each bin
RA = np.zeros(RA_In.shape[0]-1)
Dec = np.zeros(Dec_In.shape[0]-1)
Av_Size = np.zeros((Av_Size_In.shape[0]-1,Av_Size_In.shape[1]-1))
for i in range(0,RA.shape[0]):
    RA[i] = (RA_In[i]+RA_In[i+1])/2.
    for j in range(0,Dec.shape[0]):
        if(i==0):
            Dec[j] = (Dec_In[j]+Dec_In[j+1])/2.
        Av_Size[i,j] = Av_Size_In[i,j]

if KappaSize.shape[0] != RA.shape[0]:
    print 'WARNING - Kappa entered not conformal with RA', KappaSize.shape[0], RA.shape[0]
if KappaSize.shape[1] != Dec.shape[0]:
    print 'WARNING - Kappa entered not conformal with Dec', KappaSize.shape[1], Dec.shape[0]

'''
if(True):
    print 'WARNING: Smoothing in python!!!!'
    if(RA.shape[0] != 256):
        print "WARNING - RA not of the same pixel size as FITS file"
    if(Dec.shape[0] != 256):
        print "WARNING - Dec not of the same pixel size as FITS file"
    AvSize_Copy = np.copy(Av_Size)
#    Av_Size = 0.
    Av_Size = Smooth(AvSize_Copy, RA, Dec)
    Kappa_Copy = np.copy(KappaSize)
    KappaSize = Smooth(KappaSize,RA,Dec)
    print 'Smoothing Done'
'''

if(plot_Kappa_Size_Map):
    f = pl.figure()
    ax = f.add_subplot(1,1,1)
    '''Plot Structures'''
    if(labelStructures):
        ax.plot(A901a[0], A901a[1], marker='o', markerfacecolor='red', markersize=3.5)
        ax.plot(A901b[0], A901b[1], marker='o', markerfacecolor='red', markersize=3.5)
        ax.plot(A902[0], A902[1], marker='o', markerfacecolor='red', markersize=3.5)
        ax.plot(SW[0], SW[1], marker='o', markerfacecolor='red', markersize=3.5)
    
        
    contour_levels = np.zeros(10)
    for i in range(0,contour_levels.shape[0]):
        contour_levels[i] = np.nanmin(Av_Size) + i*(np.nanmax(Av_Size)-np.nanmin(Av_Size))/(contour_levels.shape[0]-1)

#Use transpose as we want RA as x coord, but input file has it as increasing by row.
#CS = ax.contourf(RA,Dec,np.transpose(Av_Size), levels = contour_levels, cmap = pl.cm.bone, origin = None) #100 after av_size
    CS = ax.contourf(RA,Dec,np.transpose(Av_Size), 25, cmap = pl.cm.bone, origin = None) #100 after av_size
    cbar = pl.colorbar(CS)
    cbar.ax.set_ylabel('Average Size')
    ax.invert_xaxis()

#im = ax.imshow(np.transpose(Av_Size), interpolation = None)
    if(plot_Size_Contours):
        CS2 = ax.contour(RA,Dec,np.transpose(Av_Size),levels = contour_levels[-3:],origin = None, hold = 'on', colors = ('r'))
    if(plot_KappaSize_Contours):
        CS_Kappa = ax.contour(RA,Dec,np.transpose(KappaSize), origin = None, hold = 'on', colors = 'k')
#pl.clabel(CS2, fmt = '%2.1f', colors = 'k', fontsize=3)
    cbar.add_lines(CS2)

    pl.xlabel('RA')
    pl.ylabel('Dec')
    ax.set_title('Red(&Filled): Av Size. Black: Kappa_Size. White: Kappa_Shape')

    ax.xaxis.set_major_formatter(pl.FormatStrFormatter('%4.1f'))

    if(plotSTAGES_DM_Map):
        DM_Map_Plot = ax.contour(RA_DM_Map,Dec_DM_Map,DM_Map,origin = None, hold = 'on', colors = 'w')

        '''OUTPUT'''
    pl.savefig(Output_Filename_avSize,format = 'pdf',bbox_inches = 'tight')
    print 'Produced plot output to '+Output_Filename_avSize
    os.system('evince '+Output_Filename_avSize+' &')


if(scatter_Plot):
    #Produces FITS kappa versus size scatter plot
    if(plotSTAGES_DM_Map == False):
        '''
        DM_Map = pf.getdata(DM_fits_file)
    #Fits file is 256x256 pixels - each pixel has a scale of 7.6 arcseconds#                         
        RA_DM_Map = np.zeros(256)
        Dec_DM_Map = np.zeros(256)
    #Minvalue and dn have been taken from fits file header                                           
        for i in range(0,Dec_DM_Map.shape[0]):
            RA_DM_Map[i] = 149.3524 + i*(-0.002120962)
#RA_DM_Map[i] = 148.85 + i*(149.31-148.85)/256                                                       
            Dec_DM_Map[i] = -10.29099 + i*(0.002120962)
        '''
        DM_Map, RA_DM_Map, Dec_DM_Map = read_in_DM_Map_ShearSTAGES(DM_fits_file)
        

    #THIS INTERPOLATION DOES NOT WORK YET - OUT OF RANGE#
    '''
    if((RA_DM_Map.shape[0] != RA.shape[0]) or (Dec_DM_Map.shape[0] != Dec.shape[0])):                
        print "ERROR - KAPPA MAP and AvSize are not conformal in RA or Dec"
        exit()
    '''
    fig = pl.figure()
    axes = fig.add_subplot(1,1,1)

    counter = 0
    for i in range(0, RA.shape[0],10): #Skip 10 at a time to remove correlations due to smoothing                              
        for j in range(0,Dec.shape[0],10):
            if(OccGrid[i,j] >= 1):
                if(math.isnan(DM_Noise_Map[i,j])):
                    continue
                counter += 1

    xSc = np.zeros(counter)#RA.shape[0]*Dec.shape[0])
    xError = np.zeros(xSc.shape[0])
    ySc = np.zeros(xSc.shape[0])
    yError = np.zeros(ySc.shape[0])
    counter = 0
    iSucc = 0
    jSucc = 0
    nNaN = 0
    for i in range(0, RA.shape[0],10): #Skip 10 at a time to remove correlations due to smoothing
        for j in range(0,Dec.shape[0],10):
            if(OccGrid[i,j] >= 1):
#                if( math.sqrt(float(i-iSucc)**2+float(j-jSucc)**2) >= 10. ): #This only ensures that it is 10 pix from the last successful pixel - whihc means there may still be correlations between other successful pixels
                if(math.isnan(DM_Noise_Map[i,j])):
                    nNaN += 1
                    continue
                iSucc = i
                jSucc = j
                xSc[counter] = DM_Map[i,j]
                xError[counter] = DM_Noise_Map[i,j]
                ySc[counter] = KappaSize[i,j]
                yError[counter] = KappaSize_Noise[i,j]
                counter += 1

    print 'Number of successful scatter plot points:', counter, ' with ', nNaN, " NaN's ignored"
    axes.errorbar(xSc, ySc, xerr = xError, yerr = yError, linewidth = 0.1, mfc = None, fmt = 'x')
    axes.scatter(xSc,ySc, s = 0.3)
    axes.plot(axes.get_ylim(), axes.get_ylim())
    axes.set_xlabel(r'$\kappa_{\epsilon}$')
    axes.set_ylabel(r'$\kappa_{\rm Size}$')

    #Define new statistics which defines the differnece between the two estimates of kappa#
    print 'Difference between convergence from size and shape:'
    DiffKappa = np.zeros(xSc.shape[0])
    DiffKappa = np.subtract(ySc, xSc)
    print 'Mean:', np.sum(DiffKappa)/DiffKappa.shape[0]

    DiffKappa_MeanVar = 0.
    for i in range(0, DiffKappa.shape[0]):
        DiffKappa_MeanVar += math.pow(xError[i],2.) + math.pow(yError[i],2.)
    DiffKappa_MeanVar = math.sqrt(DiffKappa_MeanVar)/DiffKappa.shape[0]
    print 'Variance of mean:', DiffKappa_MeanVar
    

    print 'Outputting Scatter Plot'
    pl.savefig(Output_Filename_kappaScatter,format = 'pdf',bbox_inches = 'tight')
    print 'Produced plot output to '+Output_Filename_kappaScatter
    os.system('evince '+Output_Filename_kappaScatter+' &')

if(plot_Occupation_Grid):
    fig = pl.figure(figsize = (8.,12.))
    axes = fig.add_subplot(2,1,1)
    
    #-Plot unsmoothed occupation grid-#
    CS = axes.imshow(np.transpose(OccGrid), extent = (np.nanmin(RA), np.nanmax(RA), np.nanmin(Dec), np.nanmax(Dec)), origin = 'lower',cmap = pl.cm.bone, interpolation = None)
    cbar = pl.colorbar(CS, ax = axes)
    cbar.ax.set_ylabel('Occupation')

    axes.invert_xaxis()

    '''

    CS = axes.contourf(RA,Dec,np.transpose(OccGrid), 25, cmap = pl.cm.bone) #100 after av_size
    cbar = pl.colorbar(CS)
    cbar.ax.set_ylabel('Occupation Grid')
    axes.invert_xaxis()
    '''
    axes.set_title('Unsmoothed')
    axes.set_xlabel('RA')
    axes.set_ylabel('Dec')
    
    axes.xaxis.set_major_formatter(pl.FormatStrFormatter('%4.1f'))
    
    axes = fig.add_subplot(2,1,2, aspect = 'equal')
    CS = axes.contourf(RA,Dec,np.transpose(SmOccGrid), 25, cmap = pl.cm.bone)
    axes.invert_xaxis()

    CS2 = axes.contour(RA,Dec,np.transpose(SmOccGrid), hold = 'on', colors = 'k')
    cbar = pl.colorbar(CS, ax = axes)
    cbar.ax.set_ylabel('Occupation')

    axes.set_title('Smoothed')
    axes.set_xlabel('RA')
    axes.set_ylabel('Dec')

    axes.xaxis.set_major_formatter(pl.FormatStrFormatter('%4.1f'))

    if(labelStructures):
        axes.plot(A901a[0], A901a[1], marker='o', markerfacecolor='red', markersize=3.5)
        axes.plot(A901b[0], A901b[1], marker='o', markerfacecolor='red', markersize=3.5)
        axes.plot(A902[0], A902[1], marker='o', markerfacecolor='red', markersize=3.5)
        axes.plot(SW[0], SW[1], marker='o', markerfacecolor='red', markersize=3.5)



    pl.savefig(Output_Filename_OccGrid,format = 'pdf',bbox_inches = 'tight')
    print 'Produced plot output to '+Output_Filename_OccGrid
    os.system('evince '+Output_Filename_OccGrid+' &')
              
if(plot_Noise_Maps):
    fig = pl.figure()
        
    axes = fig.add_subplot(1,1,1)# aspect = 'equal')
    if(labelStructures):
        axes.plot(A901a[0], A901a[1], marker='o', markerfacecolor='red', markersize=3.5)
        axes.plot(A901b[0], A901b[1], marker='o', markerfacecolor='red', markersize=3.5)
        axes.plot(A902[0], A902[1], marker='o', markerfacecolor='red', markersize=3.5)
        axes.plot(SW[0], SW[1], marker='o', markerfacecolor='red', markersize=3.5)
    

    
    CS = axes.contourf(RA,Dec,np.transpose(KappaSize_Noise), 25, cmap = pl.cm.bone)
    axes.invert_xaxis()

    CS2 = axes.contour(RA,Dec,np.transpose(KappaSize_Noise), hold = 'on', colors = 'k')
    cbar = pl.colorbar(CS, ax = axes)
    cbar.ax.set_ylabel('Error')

    axes.set_title(r'Error on Smoothed $\kappa_{\rm Size}$ Field')
    axes.set_xlabel('RA')
    axes.set_ylabel('Dec')

    #CS_Shape = axes.contourf(RA_DM_Map,Dec_DM_Map,np.transpose(DM_Noise_Map), hold = 'on', cmap = pl.cm.bone)# colors = 'r')

    axes.xaxis.set_major_formatter(pl.FormatStrFormatter('%4.1f'))


    pl.savefig(Output_Filename_Noise_Map,format = 'pdf',bbox_inches = 'tight')
    print 'Produced plot output to '+Output_Filename_Noise_Map
    os.system('evince '+Output_Filename_Noise_Map+' &')

if(plot_Signal_to_Noise):
    f = pl.figure()
    ax = f.add_subplot(1,1,1)
    '''Plot Structures'''
    if(labelStructures):
        ax.plot(A901a[0], A901a[1], marker='o', markerfacecolor='red', markersize=3.5)
        ax.plot(A901b[0], A901b[1], marker='o', markerfacecolor='red', markersize=3.5)
        ax.plot(A902[0], A902[1], marker='o', markerfacecolor='red', markersize=3.5)
        ax.plot(SW[0], SW[1], marker='o', markerfacecolor='red', markersize=3.5)
    
    '''    
    contour_levels = np.zeros(10)
    for i in range(0,contour_levels.shape[0]):
        contour_levels[i] = np.nanmin(Av_Size) + i*(np.nanmax(Av_Size)-np.nanmin(Av_Size))/(contour_levels.shape[0]-1)
    '''
#Use transpose as we want RA as x coord, but input file has it as increasing by row.
#CS = ax.contourf(RA,Dec,np.transpose(Av_Size), levels = contour_levels, cmap = pl.cm.bone, origin = None) #100 after av_size
    CS = ax.contourf(RA,Dec,np.absolute(np.transpose(Kappa_Size_SN)), 25, cmap = pl.cm.bone, origin = None) #100 after av_size
    cbar = pl.colorbar(CS)
    cbar.ax.set_ylabel('(Absolute) Signal to Noise')
    ax.invert_xaxis()

    CS_Kappa = ax.contour(RA,Dec,np.absolute(np.transpose(Kappa_Size_SN)), origin = None, hold = 'on', colors = 'k')
#pl.clabel(CS2, fmt = '%2.1f', colors = 'k', fontsize=3)

    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')
    ax.set_title('White: Shape. Black (and filled): Size')

    ax.xaxis.set_major_formatter(pl.FormatStrFormatter('%4.1f'))

    if(plotSTAGES_DM_Map):
        DM_Map_Plot = ax.contour(RA_DM_Map,Dec_DM_Map,np.absolute(DM_SN_Map),origin = None, hold = 'on', colors = 'w')


    '''OUTPUT'''
    pl.savefig(Output_Filename_StNoise,format = 'pdf',bbox_inches = 'tight')
    print 'Produced plot output to '+Output_Filename_StNoise
    os.system('evince '+Output_Filename_StNoise+' &')
