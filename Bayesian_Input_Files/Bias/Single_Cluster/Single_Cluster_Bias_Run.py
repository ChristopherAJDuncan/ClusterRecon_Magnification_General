import os

#To Do:
# ~~ Delete old contents of directory

CWD = os.getcwd()

Output_Paths = ['/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeOnly/SingleA901a_Bias/COMBO/', \
                '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeOnly/SingleA901a_Bias/STAGES/', \
                '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/SingleA901a_Bias/COMBO/', \
                '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/SizeMag/SingleA901a_Bias/STAGES/']

Ini_Pathway = CWD+'/'

Ini_Files = ['Single_A901a_COMBO_SizeOnly.ini', \
             'Single_A901a_STAGES_SizeOnly.ini', \
             'Single_A901a_COMBO_SizeMag.ini', \
             'Single_A901a_STAGES_SizeMag.ini' ]

Program_Pathway = '/disk1/cajd/Size_Magnification/'
Program_Name = './Bayesian_DM_Profile_Constraints.exe'

def set_OutputDirectory_IniFile(ODir, IniFile):
    # NOT WORKING YET - Invalid SED option...#
    com = 's/'+'Output_Directory =.*'+'/'+'Output_Directory = \"'+ODir+'\" /g'
    com = 'sed \''+com+'\' <'+IniFile+'> new.ini'
    os.system(com)
    print 'Check Output Directory for Input File..', IniFile
    exit()
    os.system('mv new.ini '+IniFile)

##Change Directory to the Program##
os.chdir(Program_Pathway)

print 'Run:'
for i in range(len(Ini_Files)):
#    set_OutputDirectory_IniFile(Output_Paths[i], Ini_Pathway+Ini_Files[i])
    if(os.path.exists(Ini_Pathway+Ini_Files[i]) == False):
        print 'ERROR - '+Ini_Pathway+Ini_Files[i]+' does not exist'
        pass
    if(os.path.exists(Output_Paths[i]) == False):
        print 'ERROR - '+Output_Paths[i]+' does not exist'
        exit()
        pass
     
    os.system(Program_Name+' '+Ini_Pathway+Ini_Files[i]+ ' > '+Output_Paths[i]+'SingleClusterLog &')

    print Program_Name+' '+Ini_Pathway+Ini_Files[i]+ ' > '+Output_Paths[i]+ 'SingleClusterLog &'
    print ' '

