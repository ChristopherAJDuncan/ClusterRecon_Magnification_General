import os
import subprocess

#To Do:
# ~~ Delete old contents of directory

CWD = os.getcwd()

# Where Do Run == 0, this run is ignored #
Do_Run = [1,0,0,0]

## OUTPUT DIRECTORIES MUST HAVE EACH DIRECTORY LEVEL SEPERATED WITH '\/' NOT JUST '/' TO ENSURE SED WILL WORK ##
Output_Paths = ['\/disk1\/cajd\/Size_Magnification\/Bayesian_DM_Profile_Constraints_Output\/SizeMag\/COMBO\/Bias\/', \
                '\/disk1\/cajd\/Size_Magnification\/Bayesian_DM_Profile_Constraints_Output\/SizeMag\/STAGES\/Bias\/', \
                '\/disk1\/cajd\/Size_Magnification\/Bayesian_DM_Profile_Constraints_Output\/SizeOnly\/COMBO\/Bias\/', \
                '\/disk1\/cajd\/Size_Magnification\/Bayesian_DM_Profile_Constraints_Output\/SizeOnly\/STAGES\/Bias\/']

Ini_Pathway = CWD+'/'

Ini_Files = ['Mock_COMBO_Bias_SizeMag.ini', \
             'Mock_STAGES_Bias_SizeMag.ini', \
             'Mock_COMBO_Bias_SizeOnly.ini', \
             'Mock_STAGES_Bias_SizeOnly.ini' ]

Program_Pathway = '/disk1/cajd/Size_Magnification/'
Program_Name = './Bayesian_DM_Profile_Constraints.exe'

def set_OutputDirectory_IniFile(ODir, IniFile):
    # NOT WORKING YET - Invalid SED option...#
    com = 's/'+'Output_Directory =.*'+'/'+'Output_Directory = \"'+ODir+'\" /g'
    com = 'sed \''+com+'\' <'+IniFile+'> '+Ini_Pathway+'new.ini'
    subprocess.call(com, shell = True)
    os.system('mv '+Ini_Pathway+'new.ini '+IniFile)

##Change Directory to the Program##
os.chdir(Program_Pathway)

print 'Run:'
for i in range(len(Ini_Files)):
    if(Do_Run[i] == 0):
        continue
    set_OutputDirectory_IniFile(Output_Paths[i], Ini_Pathway+Ini_Files[i])
    if(os.path.exists(Ini_Pathway+Ini_Files[i]) == False):
        print 'ERROR - '+Ini_Pathway+Ini_Files[i]+' does not exist'
        continue
#    if(os.path.exists(Output_Paths[i]) == False):
#        print 'ERROR - '+Output_Paths[i]+' does not exist'
#        exit()
#        pass
     
    os.system(Program_Name+' '+Ini_Pathway+Ini_Files[i]+ ' > '+Output_Paths[i]+'BiasLog &')

    print Program_Name+' '+Ini_Pathway+Ini_Files[i]+ ' > '+Output_Paths[i]+ 'BiasLog &'
    print ' '

os.chdir(CWD)
