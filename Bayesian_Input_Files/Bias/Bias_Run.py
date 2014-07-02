#Produces the Bias for STAGES clusters, including overlap, and intrinsic bias in the method (e.g. CHo8 redsift distribution), but not the bias that results from using a lensed distribution

import os
import subprocess

#To Do:
# ~~ Delete old contents of directory

CWD = os.getcwd()
Program_Pathway = '/disk1/cajd/Size_Magnification/'

# Where Do Run == 0, this run is ignored #
Do_Run = [1,0,0,0,0,0]

## THERE MUST NOT BE A TRAILING '/' IN THESE DIRECTORY DECLARATIONS
Output_Paths = [Program_Pathway+"Bayesian_DM_Profile_Constraints_Output/SizeMag/STAGES_Clusters_Bias/STAGES+COMBO/WLWL/SM+Mag_P3.3_D6.5", \
                Program_Pathway+"Bayesian_DM_Profile_Constraints_Output/SizeOnly/STAGES+COMBO", \
                Program_Pathway+'Bayesian_DM_Profile_Constraints_Output/SizeMag/COMBO', \
                Program_Pathway+'Bayesian_DM_Profile_Constraints_Output/SizeMag/STAGES/StrongLensingSTAGES_Check', \
                Program_Pathway+'Bayesian_DM_Profile_Constraints_Output/SizeOnly/COMBO', \
                Program_Pathway+'Bayesian_DM_Profile_Constraints_Output/SizeOnly/STAGES']

Cluster_Pathway = Program_Pathway+'Cluster_Parameters'
Cluster_Filename = 'STAGES.ini'

Ini_Pathway = CWD+'/'

Ini_Files = ['Mock_STAGES+COMBO_Bias_SizeMag.ini',\
             'Mock_STAGES+COMBO_Bias_SizeOnly.ini',\
             'Mock_COMBO_Bias_SizeMag.ini', \
             'Mock_STAGES_Bias_SizeMag.ini', \
             'Mock_COMBO_Bias_SizeOnly.ini', \
             'Mock_STAGES_Bias_SizeOnly.ini' ]

Program_Name = './Bayesian_DM_Profile_Constraints.exe'


def Total_Split(Pathway):
    Path = []
    head = Pathway
    while 1:
        head, tail = os.path.split(head)
        if tail != "":
            Path.append(tail)
        else:
            if head!="":
                Path.append(head)
            break
    return Path[::-1]

def regEx_Path(Pathway):
    Split = Total_Split(Pathway)
    regEx = '\/'
    for i in range(1, len(Split)):
        regEx += Split[i]+'\/'

    return regEx

def set_OutputDirectory_IniFile(ODir, IniFile):
    regDir = regEx_Path(ODir)
    com = 's/'+'Output_Directory =.*'+'/'+'Output_Directory = \"'+regDir+'\" /g'
    com = 'sed \''+com+'\' <'+IniFile+'> '+Ini_Pathway+'new.ini'
    subprocess.call(com, shell = True)
    os.system('mv '+Ini_Pathway+'new.ini '+IniFile)

def set_ClusterDirectory_IniFile(CDir, CFile, IniFile):
    regDir = regEx_Path(CDir)
    com = 's/'+'Cluster_Filename =.*'+'/'+'Cluster_Filename = \"'+regDir+CFile+'\" /g'
    com = 'sed \''+com+'\' <'+IniFile+'> '+Ini_Pathway+'new.ini'
    subprocess.call(com, shell = True)
    os.system('mv '+Ini_Pathway+'new.ini '+IniFile)

def create_Directory(Dir):
    #Unpacks Directory and attempts to create a directory level by level if it doesn't already exist
    Split = Total_Split(Dir)
    TestDirectory = '/'
    for i in range(0, len(Split)):
        TestDirectory += Split[i]+'/'
        if(os.path.exists(TestDirectory) == False):
            os.system('mkdir '+TestDirectory)
                                            


##Change Directory to the Program##
os.chdir(Program_Pathway)

print 'Run:'
for i in range(len(Ini_Files)):
    if(Do_Run[i] == 0):
        continue
    set_OutputDirectory_IniFile(Output_Paths[i], Ini_Pathway+Ini_Files[i])
    set_ClusterDirectory_IniFile(Cluster_Pathway, Cluster_Filename, Ini_Pathway+Ini_Files[i])
    create_Directory(Output_Paths[i])
    if(os.path.exists(Ini_Pathway+Ini_Files[i]) == False):
        print 'ERROR - '+Ini_Pathway+Ini_Files[i]+' does not exist'
        continue
     
    os.system(Program_Name+' '+Ini_Pathway+Ini_Files[i]+ ' > '+Output_Paths[i]+'/BiasLog &')

    print Program_Name+' '+Ini_Pathway+Ini_Files[i]+ ' > '+Output_Paths[i]+ '/BiasLog &'
    print ' '

os.chdir(CWD)
