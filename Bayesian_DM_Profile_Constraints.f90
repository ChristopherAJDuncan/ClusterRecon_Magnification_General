program Bayesian_DM_Profile_Constraints
  !--Uses Bayesian Means to Calculate Posterior for the Dark Matter Profile variables--!                                                                   
  use Param_Types; use Catalogues; use Bayesian_Routines; use Foreground_Clusters
  implicit none

  character(200), parameter::Output_Directory_Default = '/disk1/cajd/Size_Magnification/Bayesian_DM_Profile_Constraints_Output/DELETE/'
   character(200):: Output_Directory = Output_Directory_Default
   character(200)::Cluster_Filename = '/disk1/cajd/Size_Magnification/Cluster_Parameters/STAGES.ini'
   integer,allocatable:: Catalogue_Identifier(:)
   integer,allocatable::Blank_Field_Catalogue_Identifier(:)
   character(200),allocatable:: Bias_Output_Directory(:)
   integer:: Bin_By_Magnitude_Type = 2 !-1:Absolute Magnitude, 2: Apparent Magnitude-!
   integer:: Run_Type = 1 !-1:SingleRun, 2:Bias And Error-!
   logical:: ReRun_Mocks = .false.

   type(Foreground):: Clusters
!   real(double)::Lens_Redshift = 0.165e0_double

    !----Mock Parameters-------!
    integer:: nSources = 70000
    real(double)::frac_z = 0.1e0_double

   !--Command Line Argument Entry--!                
   integer::narg, i
   character(120)::arg
   character(200):: Input_Filename

   real(double),allocatable::Cluster_Posteriors(:,:,:)

   real(double),allocatable:: Cluster_Aperture_Radius(:)

   logical:: Here

   !--------------------------------------------Read in by Command Line-------------------------------------------!              
   !##############################################################################################################!                                         
   allocate(Catalogue_Identifier(1)); Catalogue_Identifier = 5
   allocate(Blank_Field_Catalogue_Identifier(1)); Blank_Field_Catalogue_Identifier = -5
   Analyse_with_Physical_Sizes = .false.
   narg = iargc()
   if(narg == 0) then
      !-Defaults-!                                                                                                                                          
      Output_Directory = trim(Output_Directory_Default)
      write(*,'(A)') 'No arguments entered, enter "Catalogue_Identifier"; "Output_Directory"; "Use Physical Sizes"; "Run_Type"; "Blank_Field_Catalogue_Identifier", "Cluster Parameter FileName (longhand!)" in that order, or an Input File ONLY'
      write(*,'(A)') 'Note:: Catalogue_Identifier:: 1:STAGES_shear, 2: COMBO matched, 3:RRG, 4:STAGES-Mock, 5:COMBO-Mock; -1:Mock_COMBO and Mock_STAGES (only for Run_Type = 2)'
      write(*,'(A)') '    :: Run_Type:: 1: Single_Run, 2:Bias and Error (MultipleRun), 3: Bias and Error (STAGES Clusters, Single Cluster w Varying Mass)'
      write(*,'(A)') ' '
      write(*,'(A)') 'Press <ENTER> to use defaults:'
      print *, 'Catalogue Identifier:', Catalogue_Identifier, ';'
      print *, 'Output_Directory: ', trim(Output_Directory)
      print *, 'Use Physical Sizes: ', Analyse_with_Physical_Sizes
      print *, 'Run Type: ', Run_Type
      print *, 'Blank Field Identifier:', Blank_Field_Catalogue_Identifier
      print *, 'Cluster Filename:', Cluster_Filename
      write(*,'(A)') ' all other declarations will use their hardwired defaults'
      read(*,*)
   elseif(narg == 1) then !-Input File-!
      call get_command_argument(1,arg)
      Input_Filename = trim(adjustl(arg))
      call Parameter_Input(Input_Filename)
   else
      print *, 'Number of arguements:', narg
      do i =1, narg
         call get_command_argument(i,arg)
         select case(i)
         case(1) !-Catalogue_Identifier-!                                                                                                                
            read(arg, *) Catalogue_Identifier(1)
            Catalogue_Identifier = Catalogue_Identifier(1)
        case(2)
           Output_Directory = trim(adjustl(arg))
        case(3)
           read(arg, *) Analyse_with_Physical_Sizes
        case(4)
           read(arg, '(I1)') Run_Type
        case(5)
            read(arg, '(I3)') Blank_Field_Catalogue_Identifier(1)
            Blank_Field_Catalogue_Identifier = Blank_Field_Catalogue_Identifier(1)
         case(6)
            Cluster_Filename = trim(adjustl(arg))
        case default
           STOP 'INCORRECT NUMBER OF ARGUEMENTS ENTERED'
        end select
     end do
     if(Cluster_Filename(1:1) /= '/') STOP 'Input Cluster Filename MUST be referenced using an absolute address, stopping...'
     if(Output_Directory(1:1) /= '/') STOP 'Output Directory MUST be referenced using an absolute address, stopping...'
  end if
  
  !--Set double Identifier if zero--!
  if(Catalogue_Identifier(1) == 0) then
     deallocate(Catalogue_Identifier); allocate(Catalogue_Identifier(2));
     Catalogue_Identifier = (/4,5/)
  end if
  if(Blank_Field_Catalogue_Identifier(1) == 0) then
     deallocate(Blank_Field_Catalogue_Identifier); allocate(Blank_Field_Catalogue_Identifier(2));
     Blank_Field_Catalogue_Identifier = (/-4,-5/)
  end if


   call get_Clusters(Clusters, Cluster_Filename)
   allocate(Cluster_Aperture_Radius(size(Clusters%Position,1))); Cluster_Aperture_Radius = 1.e0_double/60.e0_double !-Degrees-! Default - 1/60
   write(*,'(A)') 'Aperture Radius HAS BEEN SET OT 1 ARCMINUTE'
   call distance_between_Clusters(Clusters%Position, Clusters%Redshift(1))



  !--Set the Blank Field Catalogue to the Data Catalogue if invalid value (>0) entered--!
!!$  where((Blank_Field_Catalogue_Identifier > 0))
!!$     Blank_Field_Catalogue_Identifier = Catalogue_Identifier
!!$  end where

  if(Run_Type ==3) then
     !--Set Defaults for this Run Type--!
     Blank_Field_Catalogue_Identifier = Catalogue_Identifier
     !--This routine tests the effect of not having a blank feild (i.e. knowledge of the iontrinsic distributions, only taken from lensed fields.) Therefore the intrinsic distributions should come from the lensed data--!
  end if

  if(Run_Type == 1) then
     print *, 'Using Catalogue_Identifier: ', Catalogue_Identifier(1)
  else
     print *, 'Using Catalogue_Identifier: ', Catalogue_Identifier
  end if
  print *, 'and Blank Field Catalogue:', Blank_Field_Catalogue_Identifier
  print *, 'and Output_Directory: ', trim(Output_Directory)
  print *, 'Using Physical Sizes Data-Side (with longer run-time):', Analyse_with_Physical_Sizes
  print *, 'Run_Type:', Run_Type
  print *, 'Cluster Filename:', Cluster_Filename
  print *, ' '

  inquire(directory = trim(Output_Directory), exist = here)
  if(here == .false.) call system('mkdir '//trim(Output_Directory))

  select case(Run_Type)
  case(1)
     print *, '!---- Producing a single Mass Estimate:'
     call Mass_Estimate_Single_Run(Output_Directory, Cluster_Posteriors, Clusters, Cluster_Aperture_Radius, Dist_Directory = Output_Directory, Cat_Ident = Catalogue_Identifier(1), Blank_Field_Cat_Ident = Blank_Field_Catalogue_Identifier(1))
  case(2)
     !--Set Output Directory Names--!
     allocate(Bias_Output_Directory(size(Catalogue_Identifier))); Bias_Output_Directory = ' '
     do i = 1, size(Catalogue_Identifier)
        Bias_Output_Directory(i) = trim(Output_Directory)
     end do
     call Posterior_Maximum_Likelihood_Bias_Error(Catalogue_Identifier, Bias_Output_Directory, Blank_Field_Catalogue_Identifier, Cluster_Filename, Cluster_Aperture_Radius)
  case(3)
     allocate(Bias_Output_Directory(size(Catalogue_Identifier))); Bias_Output_Directory = ' '
     do i = 1, size(Catalogue_Identifier)
        write(Bias_Output_Directory(i), '(I1)') Catalogue_Identifier(i)
        Bias_Output_Directory(i) = trim(Output_Directory)
     end do     
     call Cluster_Bias_Outputs(Catalogue_Identifier, Bias_Output_Directory, Blank_Field_Catalogue_Identifier)
  case default
     STOP 'Run_type Entered not supported'
  end select

contains

  subroutine Parameter_Input(Input_File)
    use Bayesian_Routines, only: Surface_Mass_Profile, Posterior_Method, use_KDE_Smoothed_Distributions, Survey_Magnitude_Limits, Survey_Size_Limits
    character(*), intent(in):: Input_File

    !Global Declarations! (see above)
    !--Cluster_Filename, Run_Type, Output_Directory, Blank_Field_Catalogue_Identifier, Catalogue_Identifier--!

    logical:: here

    namelist/Run/Surface_Mass_Profile, Posterior_Method, use_KDE_Smoothed_Distributions, Run_Type, Cluster_Filename, Output_Directory, Blank_Field_Catalogue_Identifier, Catalogue_Identifier, Survey_Size_Limits, Survey_Magnitude_Limits
    namelist/Mocks/ nSources, frac_z, ReRun_Mocks

    !--Check for existence of Input File
    inquire(file = Input_File, exist = Here)
    if(here == .false.) then 
       print *, 'Input_File:', Input_File
       STOP 'Cluster Input File not present'
    end if

    open(unit = 92, file = Input_File)

    read(92, nml = Run)
    read(92, nml = Mocks)

    close(92)

    !--Copy Input File to output directory--!
    call system('cp '//trim(Input_File)//' '//trim(Output_Directory)//'Input_File.ini')

  end subroutine Parameter_Input

  subroutine Cluster_Bias_Outputs(Cat_Ident, Directory, Blank_Field_Cat_Ident)
    !----- Description:
    !~~~~~ Compares the input cluster parameters (as run through Mock Catalogue) to those returned by the Pipeline. Does thi in two ways:
    !~~~~~~~~~~~~~ 1. Takes (1) cluster (editable), and varies the parameter values through user-defined inputs. For each Parameter Value (e.g. r200 for NFW) it returns a mode point of multiple posterior runs, with the variance of that mode point
    !~~~~~~~~~~~~ 2. Takes the four STAGES clusters, and their expected Parameter Value (from CH08), and retunrs the mode point and vauraicne of multiple runs for this input. Does not loop over parameter values. Alos will contain the effects of any overlap in mass profiles.
    !~~~~ Mode values and variance are taken from calls to Posterior_Maximum_Likelihood_Bias_Error
 
    use Mass_Profiles, only:integrated_mass_within_radius
    !~ToDo:
    
    
    integer, intent(in)::Cat_Ident(:)
    character(*), intent(in)::Directory(:)
    integer, intent(in)::Blank_Field_Cat_Ident(:)
!    integer,intent(in)::Type !3: STAGES Clusters 4:Single Cluster

    !--Cluster Parameter Declarations--!
    integer::nClusters
    real(double),dimension(:,:),allocatable::ParameterValues
    integer::SMD_Type = 3
    real(double),dimension(:),allocatable::Position !--Move To Centre?--!
    character(200):: Mock_Cluster_Filename = '/disk1/cajd/Size_Magnification/Cluster_Parameters/Cluster_Bias_Parameter.ini'
    real(double)::Redshift = 0.165e0_double
    real(double),allocatable:: iCluster_Aperture_Radius(:)

    !--Free Parameter Loop--!
    !--Used only in the Single Cluster Multiple Run Case--!
    integer::nParam = 6
    real(double):: Param_Min = 0.5, Param_Max = 2.5

    integer:: nRun

    character(120):: Run_Directory(size(Directory))

    real(double),allocatable:: Bias_Mode(:,:), Mode_Error(:,:) !-Catalogue, Aperture-!
    real(double),allocatable:: Param_Bias_Mode(:,:,:), Param_Mode_Error(:,:,:) !-Parameter, Catalogue, Aperture-!

    real(double)::Aperture_Mass, eAperture_Mass(1), Aperture_Mass_Input, Discardable(1)

    type(Foreground)::Cluster
    character(4)::Profile_String
    character(5)::fmtString
    logical::here

    !--Loops--!
    integer::i, cid, j, Ap

    Run_Directory = 'Cluster_Bias_Output/'
    do i =1, size(Directory)
       Run_Directory(i) = trim(Directory(i))//trim(Run_Directory(i))
       inquire(directory = Run_Directory(i), exist = here)
       if(here == .false.) call system('mkdir '//trim(Run_Directory(i)))
    end do

    select case(SMD_Type)
    case(1)
       Profile_String = 'Flat'
    case(2)
       Profile_String = 'SIS'
    case(3)
       Profile_String = 'NFW'
    end select

    print *, 'Cluster Bias output to:', Directory
    print *, 'For catalogues:', Cat_Ident
    print *, 'And Blank Field Catalogues:', Blank_Field_Cat_Ident

    !--Get The STAGES Clusters Biases--!
    !--Read in STAGES lcuster file to get positions an ddefault parameter values--!
    call get_Clusters(Cluster, '/disk1/cajd/Size_Magnification/Cluster_Parameters/STAGES.ini')
    nParam = 1; nClusters = size(Cluster%Position,1)
    allocate(Param_Bias_Mode(nParam, size(Cat_Ident), nClusters)); Param_Bias_Mode =0.e0_double
    allocate(Param_Mode_Error(nParam, size(Cat_Ident), nClusters)); Param_Mode_Error =0.e0_double
    allocate(ParameterValues(nParam,nClusters)); ParameterValues = 0.e0_double

    allocate(Position(2*nClusters))
    do i =1, size(Cluster%Position,1)
       Position(2*i-1) = Cluster%Position(i,1); Position(2*i) = Cluster%Position(i,2)
    end do
    ParameterValues(1,:) = Cluster%DM_Profile_Parameter

    !---Create the Parameter File----!
    call create_ClusterFile(Mock_Cluster_Filename, nClusters, SMD_Type, Position, ParameterValues(1,:), (/Redshift/))

    !--Read in Cluster--!
    call get_Clusters(Cluster, Mock_Cluster_Filename)
    allocate(ICluster_Aperture_Radius(size(Cluster%Position,1))); ICluster_Aperture_Radius = Cluster_Aperture_Radius(1) !-Degrees-! Default - 1/60

    !--Run the Bias Routine---!
    !--Bias_Mode and Mode_Variacne contain the ML point and s.d. (error) in the DM model free parameter--!

    call Posterior_Maximum_Likelihood_Bias_Error(Cat_Ident, Run_Directory, Blank_Field_Cat_Ident, Mock_Cluster_Filename, iCluster_Aperture_Radius, Bias_Mode, Mode_Error)
    
    if( (size(Param_Bias_Mode(1,:,:)) /= size(Bias_Mode)) ) STOP 'Cluster_Bias_Outputs - Error in Bias Mode size'
    if( (size(Param_Mode_Error(1,:,:)) /= size(Mode_Error)) ) STOP 'Cluster_Bias_Outputs - Error in Mode_Error size'
    Param_Bias_Mode(1,:,:) = Bias_Mode; Param_Mode_Error(1,:,:) = Mode_Error
    deallocate(Bias_Mode, Mode_Error)
       
    !---Output Result (Input, Output, Variance)---!
    !--Input Parameter, Output Parameter, Variance, Input Ap_Mass, Output Ap_Mass, Variance--! 
    !--Space between each aperture--!
    do cid = 1, size(Cat_Ident)
       open(36, file = trim(Run_Directory(cid))//'Mass_Bias_duetoForeground_STAGES.dat')
       !-Header-!
       write(36, '(A)') '## Input_Parameter, Output_Parameter, Variance, Input Aperture Mass, Output Aperture Mass, Varaince'
       write(36, '(A)') '## Space between apertures'
       write(fmtString, '(I1)') size(ICluster_Aperture_Radius); write(36, '(A,'//trim(fmtString)//'(e12.5,x))') '## Aperture (Degrees):', ICluster_Aperture_Radius !--Edit fmt to array
       write(36, '(A,A)') '## Mass Profile:', Profile_String !--Edit fmt to array  
     
       do Ap = 1, nClusters
          do j =1, nParam
             call Integrated_Mass_Within_Radius(SMD_Type, ICluster_Aperture_Radius(Ap), ParameterValues(j,Ap), (/0.e0_double/), Aperture_Mass_Input, Discardable, Redshift)
             call Integrated_Mass_Within_Radius(SMD_Type, ICluster_Aperture_Radius(Ap), Param_Bias_Mode(j,cid,Ap), (/Param_Mode_Error(j,cid,Ap)/), Aperture_Mass, eAperture_Mass, Redshift)

             write(36, '(6(e12.5,x))') ParameterValues(j,Ap), Param_Bias_Mode(j,cid,Ap), Param_Mode_Error(j,cid,Ap), Aperture_Mass_Input, Aperture_Mass,  eAperture_Mass 
          end do
          write(36,*)
       end do
    end do
    close(36)

    deallocate(Param_Bias_Mode, Param_Mode_Error, Position, ParameterValues, iCluster_Aperture_Radius)
    print *, 'Finished Mock Cluster Bias Run STAGES normally'

    !--nClusters of equal mass, varying r200--!
    nClusters = 1
    allocate(Position(2*nClusters))
    Position = (/ 149.1099,-9.9561 /)
    allocate(Param_Bias_Mode(nParam, size(Cat_Ident), nClusters)); Param_Bias_Mode =0.e0_double
    allocate(Param_Mode_Error(nParam, size(Cat_Ident), nClusters)); Param_Mode_Error =0.e0_double

    allocate(ParameterValues(nParam,nClusters)); ParameterValues = 0.e0_double
    allocate(ICluster_Aperture_Radius(nClusters)); ICluster_Aperture_Radius = 1.e0_double/60.e0_double !-Degrees-! Default - 1/60
    do i =1, nParam
       !---Create the Parameter File----!
       ParameterValues(i,:) = Param_Min + (i-1)*((Param_Max-Param_min)/(nParam-1))
       call create_ClusterFile(Mock_Cluster_Filename, nClusters, SMD_Type, Position, ParameterValues(i,:), (/Redshift/))
       print *, 'Cluster File Created'
       !--Run the Bias Routine---!
       !--Bias_Mode and Mode_Variacne contain the variance in the DM model free parameter--!
       call Posterior_Maximum_Likelihood_Bias_Error(Cat_Ident, Run_Directory, Blank_Field_Cat_Ident, Mock_Cluster_Filename, iCluster_Aperture_Radius, Bias_Mode, Mode_Error)

       if( (size(Param_Bias_Mode(1,:,:)) /= size(Bias_Mode)) ) STOP 'Cluster_Bias_Outputs - Error in Bias Mode size'
       if( (size(Param_Mode_Error(1,:,:)) /= size(Mode_Error)) ) STOP 'Cluster_Bias_Outputs - Error in Mode_Error size'
       Param_Bias_Mode(i,:,:) = Bias_Mode; Param_Mode_Error(i,:,:) = Mode_Error
       deallocate(Bias_Mode, Mode_Error)
    end do
       
    !---Output Result (Input, Output, Variance)---!
    !--Input Parameter, Output Parameter, Variance, Input Ap_Mass, Output Ap_Mass, Variance--! 
    !--Space between each apeture--!
    do cid = 1, size(Run_Directory)!size(Cat_Ident)
       open(36, file = trim(Run_Directory(cid))//'Mass_Bias_duetoForeground.dat')
       !-Header-!
       write(36, '(A)') '## Input_Parameter, Output_Parameter, Error (s.d.), Input Aperture Mass, Output Aperture Mass, Error (s.d)'
       write(36, '(A)') '## Space between apertures'
       write(fmtString, '(I1)') size(ICluster_Aperture_Radius); write(36, '(A,'//trim(fmtString)//'(e12.5,x))') '## Aperture (Degrees):', ICluster_Aperture_Radius !--Edit fmt to array
       write(36, '(A,A)') '## Mass Profile:', Profile_String !--Edit fmt to array  
     
       do Ap = 1, nClusters
          do j =1, nParam
             !--WRONG - Radius needs to be in Mpc/h, not degrees!---!
             call Integrated_Mass_Within_Radius(SMD_Type, ICluster_Aperture_Radius(Ap), ParameterValues(j,Ap), (/0.e0_double/), Aperture_Mass_Input, Discardable, Redshift)
             call Integrated_Mass_Within_Radius(SMD_Type, ICluster_Aperture_Radius(Ap), Param_Bias_Mode(j,cid,Ap), (/Param_Mode_Error(j,cid,Ap)/), Aperture_Mass, eAperture_Mass, Redshift)

             write(36, '(6(e12.5,x))') ParameterValues(j,Ap), Param_Bias_Mode(j,cid,Ap), Param_Mode_Error(j,cid,Ap), Aperture_Mass_Input, Aperture_Mass,  eAperture_Mass 
          end do
          write(36,*)
       end do
    end do
    close(36)
    
    deallocate(Param_Bias_Mode, Param_Mode_Error, Position, ICluster_Aperture_Radius)
    print *, 'Finished Mock Cluster Bias Run'

    !--Plot Result--!
    do cid = 1, size(Run_Directory)
!       call system('python Cluster_Bias_Plotter.py '//trim(Run_Directory(cid))//'Mass_Bias_duetoForeground_STAGES.dat'//' '//trim(Run_Directory(cid))//'Mass_Bias_duetoForeground.dat')
    end do

  end subroutine Cluster_Bias_Outputs

  subroutine Posterior_Maximum_Likelihood_Bias_Error(Cat_Ident, Directory, Blank_Field_Cat_Ident, Input_Clusters_Filename, Aperture_Radius, Bias_Mode_Out, Mode_Error_Out)
    use Cosmology; use Statistics, only:variance_discrete; use Matrix_methods, only:Subset; use Bayesian_Routines, only:Surface_Mass_Profile; use Mass_Profiles, only:integrated_mass_within_radius
    !--Produces multiple Posteriors for many different mock catalogue realisations, to calculate the ML-point bias and variance--!
    !--Unlensed distribution is calculated from the Blank_Field_Catlaogue for the first run, and then read in on consecutive runs--!
    integer, intent(in)::Cat_Ident(:)
    character(*), intent(in)::Directory(:)
    integer, intent(in)::Blank_Field_Cat_Ident(:)
    character(*),intent(in):: Input_Clusters_Filename
    real(double), intent(in)::Aperture_Radius(:)
    real(double),allocatable,optional::Bias_Mode_Out(:,:), Mode_Error_Out(:,:) !-Catalogue, Ap-!

    real(double),allocatable::Bias_Mode(:), Mode_Error(:)

    character(200)::Run_Output_Directory, Run_Parent_Directory
    integer::nR
    integer::nRun = 20

    integer::Ap, nAp, Id

    type(Foreground)::iClusters
    real(double),allocatable::Single_Run_Posterior(:,:,:) !-Aperture, Grid/Posterior, Value-!
    real(double),allocatable:: Posteriors(:,:,:,:,:) !-CatalogueID, Run, Aperture, Grid/Posterior, Value-!
    real(double),allocatable:: Combined_Posterior(:,:,:) !-Aperture, Grid/Posterior, Value-!
    real(double),allocatable:: ML_Point(:,:) !-Ap, Point-!

    !--Conversion Declarations--!
    real(double)::D_l, Area
    character(12)::Error_String, Bias_String

    character(7)::fmtstring

    logical::here

    logical::tPhysical_Size

    type(Catalogue)::Cat, BFCat
    character(500)::Catalogue_Directory, Catalogue_Filename
    integer,dimension(:),allocatable::Catalogue_Cols

    real(double)::eCombined_Posterior(2)
    real(double):: Virial_Mass_Input, Discardable(2), Virial_Mass, eVirial_Mass(2)

    !----Mock Parameters-------!
    !--Directory formats must be referenced with respect to teh run directory, and these catalogues must be stored with the same level structure wrt to program--!
    character(500)::Mock_Output = 'Output/', Mock_Input, Mock_Input_Parent = 'Catalogues/Mock_Catalogues/'


    !--Temporary Storage of User-entered variable--!
    tPhysical_Size = Analyse_with_Physical_Sizes

!    if(all( (/4,5/) /= Cat_Ident)) STOP 'Posterior_Maximum_Likelihood_Bias_Error - Invalid Catalogue Indenifier entered, it must correspond to a Mock Catalogue'
    if(SubSet(Cat_Ident, (/4,5, 45/)) ==.false.)STOP' Posterior_Maximum_Likelihood_Bias_Error - Catalogue Identiers entered are not valid'

    if(size(Directory) /= size(Cat_Ident)) STOP' Posterior_Maximum_Likelihood_Bias_Error - Directory and Cat Identifier not of the same size, exiting'

    !--Create Mock Output Directory, which is just used to Temporarily store the mock output--!
    if(ReRun_Mocks) then
       !--Create Mock Output Directory--!
       Mock_Output = trim(Directory(1))//'Mock_Catalogues/'
       write(*,'(2(A))') 'Mock catalogues will be recreated and stored in:', trim(Mock_Output)
       inquire(directory = Mock_Output, exist = here)
       if(here == .false.) call system('mkdir '//trim(Mock_Output))
    else
       write(*,'(2(A))') 'Mock catalogues read in from:', trim(Mock_Input_Parent)
    end if

    do nR = 1, nRun
       write(*,'(A)') '!##################################################################!'
       write(*,'(A, I3)') '!------------------------------------------------------Run', nR

       if(ReRun_Mocks) then
          !- Run Mock Catalogue Production Script -!
          call run_Mock_Production_Script(Mock_Output, nSources, frac_z, Input_Clusters_Filename)
          Mock_Input = Mock_Output
       else
          if(nR < 10) then
             fmtstring = '(I1)'
          elseif(nR< 100) then
             fmtstring = '(I2)'
          else
             STOP 'Posterior_Maximum_Likelihood_Bias_Error - Number of runs set too hight for fmtstring'
          end if
          write(Mock_Input, fmtstring) nR
          Mock_Input = trim(Mock_Input_Parent)//trim(Mock_Input)//'/'
       end if

       

       call get_Clusters(iClusters, Input_Clusters_Filename)
       nAp = size(iClusters%Position,1)

       !--Initialise Output--!
       if(nR==1) then
          if(present(Mode_Error_Out)) then
             if(allocated(Mode_Error_Out)) deallocate(Mode_Error_Out)
             allocate(Mode_Error_Out(size(Cat_Ident), nAp)); Mode_Error_Out = 0.e0_double
          end if
          if(present(Bias_Mode_Out)) then
             if(allocated(Bias_Mode_Out)) deallocate(Bias_Mode_Out)
             allocate(Bias_Mode_Out(size(Cat_Ident), nAp)); Bias_Mode_Out = 0.e0_double
          end if
       end if
       
       do Id = 1, size(Cat_Ident)
          !-Check for directory existence-!
          if(nR==1) then
             if(dabs(1.e0_double*Cat_Ident(ID)) < 10) then
                write(Run_Output_Directory, '(I1)') Cat_Ident(ID)
             else
                write(Run_Output_Directory, '(I2)') Cat_Ident(ID)
             end if
!             Run_Output_Directory = trim(Output_Directory)
             Run_Parent_Directory = trim(Directory(ID))//'Bias_Error_Run_CatID'//trim(Run_Output_Directory)//'/'
             inquire(directory = Run_Parent_Directory, exist = here)
             if(here == .false.) call system('mkdir '//trim(Run_Parent_Directory))
          end if

          if(nR<10) then
             write(Run_Output_Directory, '(I1)') nR
          else
             write(Run_Output_Directory, '(I2)') nR
          end if
          Run_Output_Directory = trim(Run_Parent_Directory)//trim(Run_Output_Directory)//'/'

          !--Read in Mocks--!
          call common_Catalogue_directories(Cat_Ident(ID), Catalogue_Directory, Catalogue_Filename, Catalogue_Cols)
          Catalogue_Directory = Mock_Input
          call catalogue_readin(Cat, trim(adjustl(Catalogue_Directory))//trim(adjustl(Catalogue_Filename)), 'Tr(J)', Catalogue_Cols)
          if(Blank_Field_Cat_Ident(ID) == Cat_Ident(ID)) then
             BFCat = Cat
          else
             call common_Catalogue_directories(Blank_Field_Cat_Ident(ID), Catalogue_Directory, Catalogue_Filename, Catalogue_Cols)
             Catalogue_Directory = Mock_Input
             call catalogue_readin(BFCat, trim(adjustl(Catalogue_Directory))//trim(adjustl(Catalogue_Filename)), 'Tr(J)', Catalogue_Cols)
          end if

!!$          if(nR==1) then
!!$             !--Create Distribution--!
             call Mass_Estimate_Single_Run(Run_Output_Directory, Single_Run_Posterior, iClusters, Aperture_Radius, Dist_Directory = Directory(ID), inCat = Cat, inBFCat = BFCat)
!!$          else
             !--Distribution read in--!
!!$             call Mass_Estimate_Single_Run(Run_Output_Directory, Single_Run_Posterior, iClusters, Aperture_Radius, Dist_Directory = Directory(ID), inCat = Cat)
!!$          end if
!          call Mass_Estimate_Single_Run(Cat_Ident(ID), Run_Output_Directory, Single_Run_Posterior, Blank_Field_Cat_Ident(ID), iClusters, Aperture_Radius)
          call Catalogue_Destruct(Cat); call Catalogue_Destruct(BFCat)

          Analyse_with_Physical_Sizes = tPhysical_Size

          if(nR ==1 .and. ID == 1) then
             allocate(Posteriors(size(Cat_Ident),nRun, size(Single_Run_Posterior,1), size(Single_Run_Posterior,2), size(Single_Run_Posterior,3))); Posteriors = 0.e0_double
          end if
          Posteriors(ID,nR,:,:,:) = Single_Run_Posterior
          
          deallocate(Single_Run_Posterior)
       end do
    end do

    !--Combine Posteriors and Get ML Point for Combined_Posteriors--!
    do ID = 1, size(Cat_Ident)
       allocate(ML_Point(nAp, nRun)); ML_Point = 0.e0_double
       allocate(Combined_Posterior(nAp, 2, size(Posteriors,5))); Combined_Posterior = 0.e0_double

       Combined_Posterior(:,1,:) = Posteriors(1,1,:,1,:)
       do Ap = 1, nAp
          do nR = 1, nRun
             call Posterior_Statistics(Posteriors(ID,nR,Ap,1,:), Posteriors(ID,nR,Ap,2,:), ModeVal = ML_Point(Ap,nR))
          end do
          call Combine_Posteriors(Posteriors(ID,1,Ap,1,:), Posteriors(ID,:,Ap,2,:), .true., Combined_Posterior(Ap,2,:))
       end do

       !---Output---!
       open(37, file = trim(Directory(ID))//'Bias_Combined_Posterior.dat')
       write(fmtstring,'(I1)') nAp+1
       do i =1, size(Combined_Posterior,3)
          write(37,  '('//trim(adjustl(fmtstring))//'(e14.7,x))') Combined_Posterior(1,1,i), Combined_Posterior(:,2,i)
       end do
       close(37)
       print *, 'Output Combined Posterior for Bias to: ', trim(Directory(ID))//'Bias_Combined_Posterior.dat' 

       allocate(Bias_Mode(nAp)); Bias_Mode = 0.e0_double
       allocate(Mode_Error(nAp)); Mode_Error = 0.e0_double
       open(unit = 51, file = trim(Directory(ID))//'BiasesAndErrorVariance.dat')
       write(*,'(2(A))') '!------ Details of Biases output to ', trim(Directory(ID))//'BiasesAndErrorVariance.dat'
       write(51, '(A)') '#Following is Mode and Error around Mode of Combined Posterior:'
       write(51, '(A)') '#Aperture, Combined_Mode, Bias, Error(+ve/-ve), Input_Virial_Mass, Virial_Mass_Bias, Error_Virial_Mass(+ve/-ve)'
       write(51, '(A)') '#Difference between mode of Combined to known True Value gives indication of Bias'!       write(51, '(A)') '#Error on ML points gives indication on the relyability of the Errors from a single-run Bayesian Analysis'
!       write(51, '(A)') '# Combined ML Point, Bias, S.D. (error) of ML points'
       
       write(*,'(A)') '!##################################################################!'
       write(*,'(A)') '!##################################################################!'
       write(*,'(A)') ' Finished Combination, following is Mode (shift from known is bias) and variance (Error on Single Run Posterior Errors) of DM profile free parameter:'
       do Ap = 1, nAp
          D_l =  angular_diameter_distance_fromRedshift(0.e0_double, iClusters%Redshift(Ap))
          Area = 3.142e0_double*(D_l*((3.142e0_double*Cluster_Aperture_Radius(Ap))/180.e0_double))**2.e0_double
       
          call Posterior_Statistics(Combined_Posterior(Ap,1,:), Combined_Posterior(Ap,2,:), ModeVal = Bias_Mode(Ap), AntiSymm_Error = eCombined_Posterior)

          Mode_Error(Ap) = variance_discrete(ML_Point(Ap,:), ML_Point(Ap,:))
          if(Mode_Error(Ap) < 0.e0_double) then
             print *, '**** Aperture:', Ap, ' has invalid error on combined posterior ******'
             Mode_Error(Ap) = 1.e0_double/0.e0_double
          else
             Mode_Error(Ap) = dsqrt(Mode_Error(Ap))/(1.e0_double*size(ML_Point,2))
          end if

!          call Posterior_Statistics(Combined_Posterior(Ap,1,:), Combined_Posterior(Ap,2,:), AntiSymm_Error = eCombined_Posterior)
!          print *, 'Got Error on Combined'

          if(present(Bias_Mode_Out)) Bias_Mode_Out(ID, Ap) = Bias_Mode(Ap)
          if(present(Mode_Error_Out)) Mode_Error_Out(ID, Ap) = Mode_Error(Ap)

          !--Mass of Input--!
          call Integrated_Mass_Within_Radius(Surface_Mass_Profile, -1.e0_double, Clusters%DM_Profile_Parameter(Ap), (/0.e0_double, 0.e0_double/), Virial_Mass_Input, Discardable, iClusters%Redshift(Ap))
          !--Mass of Output--!
          call Integrated_Mass_Within_Radius(Surface_Mass_Profile, -1.e0_double, Bias_Mode(Ap), eCombined_Posterior, Virial_Mass, eVirial_Mass, iClusters%Redshift(Ap))

          write(51,'(I2,x,8(e10.4,x))') Ap, Bias_Mode(Ap), Bias_Mode(Ap) - Clusters%DM_Profile_Parameter(Ap), eCombined_Posterior(2), eCombined_Posterior(1), Virial_Mass, Virial_Mass-Virial_Mass_Input, eVirial_Mass(2), eVirial_Mass(1)

          write(Bias_String, '(e10.4)') Bias_Mode(Ap)
          write(Error_String, '(e10.4)') Mode_Error(Ap)
          write(*,'(A,I2,A,A,A,e10.4,A,e10.4)') 'Cluster ', Ap, ' has Mode:', Bias_String, ' + ',eCOmbined_Posterior(2), ' - ', eCombined_Posterior(1)
          write(*,'(A,A)') ' and ML-Error:', Error_String

          !--Calculate and output the parameter bias--!
          write(Bias_String, '(e10.4)') Bias_Mode(Ap) - Clusters%DM_Profile_Parameter(Ap)
          write(*,'(3A,e9.3)') '     giving Bias:', Bias_String, ' B/N:', (Bias_Mode(Ap) - Clusters%DM_Profile_Parameter(Ap))/minval(eCombined_Posterior)

          write(*,'(A,e10.4,A,e10.4,A,e10.4)') 'and Virial Mass:', Virial_Mass, ' + ', eVirial_Mass(2), ' - ', eVirial_Mass(1)
          call Integrated_Mass_Within_Radius(Surface_Mass_Profile, -1.e0_double, Bias_Mode(Ap), (/Mode_Error(Ap),Mode_Error(Ap)/), Virial_Mass, eVirial_Mass, iClusters%Redshift(Ap))
          write(*,'(A, e10.4)') '    with ML-Error:', eVirial_Mass(1)


       end do
       write(*,'(A)') '!##################################################################!'
       write(*,'(A)') '!##################################################################!'
       close(51)
       deallocate(Bias_mode, Combined_Posterior, ML_Point, Mode_Error)
    end do
    deallocate(Posteriors)
    call Foreground_Destruct(iClusters)

  end subroutine Posterior_Maximum_Likelihood_Bias_Error


  subroutine Mass_Estimate_Single_Run(run_Output_Dir, returned_Cluster_Posteriors, Clusters_In, Aperture_Radius, Dist_Directory, Cat_Ident, Blank_Field_Cat_Ident, inCat, inBFCat)
    !--inBFCat or Blank_Field_Cat_Ident should contain/point ot catalogues from which the intrinsic size-magnitude distribution can be obtained. If neither of these are entered, then the code will attempt to read in the distribution from the run_Output_Dir--!
    use MC_Redshift_Sampling; use Bayesian_Routines 
    real(double),allocatable,intent(out):: returned_Cluster_Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-!
    type(Foreground):: Clusters_In
    real(double),intent(in):: Aperture_Radius(:)
    character(*), intent(in):: Dist_Directory

    integer,intent(in),optional::Cat_Ident, Blank_Field_Cat_Ident
    type(Catalogue), intent(in),optional:: inCat, inBFCat

    character(*), intent(in):: run_Output_Dir

    type(Catalogue)::Catt, BFCatt
    character(120)::Catalogue_Directory, Catalogue_Filename
    integer,dimension(:),allocatable::Catalogue_Cols
    logical::here

    INTERFACE
       subroutine Mass_Estimate_Single_Run(run_Output_Dir, returned_Cluster_Posteriors, Clusters_In, Aperture_Radius, Dist_Directory, Cat_Ident, Blank_Field_Cat_Ident, inCat, inBFCat)
         use Param_Types; use Foreground_Clusters; use Catalogues
         real(double),allocatable,intent(out):: returned_Cluster_Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-!
         type(Foreground):: Clusters_In
         real(double),intent(in):: Aperture_Radius(:)
         character(*), intent(in):: run_Output_Dir
         character(*), intent(in):: Dist_Directory

         integer,intent(in),optional::Cat_Ident, Blank_Field_Cat_Ident
         type(Catalogue), intent(in),optional:: inCat, inBFCat
       END subroutine Mass_Estimate_Single_Run
    END INTERFACE

    inquire(directory = trim(run_Output_Dir), exist = here)
    if(here == .false.) call system('mkdir '//trim(adjustl(run_Output_Dir)))
    Bayesian_Routines_Output_Directory = trim(run_Output_Dir)

    !## 1: STAGES shear, 2: COMBO17, 3:RRG, 4: Mocks_STAGES; 5:Mocks_COMBO!          
    If(present(inCat)) then
       Catt = inCat
    elseif(present(Cat_Ident)) then 
       call common_Catalogue_directories(Cat_Ident, Catalogue_Directory, Catalogue_Filename, Catalogue_Cols)
       call catalogue_readin(Catt, trim(adjustl(Catalogue_Directory))//trim(adjustl(Catalogue_Filename)), 'Tr(J)', Catalogue_Cols)
    else
       STOP 'Single_Run - Either Cat_Ident of inCat must be specified'
    end If

    if(present(inBFCat) .or. present(Blank_Field_Cat_Ident)) then

       !--Readin--!
       if(present(inBFCat)) then
          BFCatt = inBFCat
       elseif(present(Blank_Field_Cat_Ident)) then
          if(Blank_Field_Cat_Ident == Cat_Ident) then
             BFCatt = Catt
          else
             call common_Catalogue_directories(Blank_Field_Cat_Ident, Catalogue_Directory, Catalogue_Filename, Catalogue_Cols)
             call catalogue_readin(BFCatt, trim(adjustl(Catalogue_Directory))//trim(adjustl(Catalogue_Filename)), 'Tr(J)', Catalogue_Cols)
          end if
       else
          STOP 'Error - Single_Run'
       end if

       !--Cuts on data catalogue--!
       call Cut_by_Magnitude(Catt, 23.e0_double) !-Taken from CH08 P1435-!   
       if(Analyse_with_Physical_Sizes) then
          call Monte_Carlo_Redshift_Sampling_Catalogue(Catt)
       end if
       call Cut_By_PhotoMetricRedshift(Catt, 0.21e0_double) !--Cut out foreground--!                                                                            
!    call Cut_By_PixelSize(Catt, 0.e0_double, 25.e0_double) !!!!!!!!!!!!!!!!!!!!!!!

       
       !--Cuts on Catalogue--!
!!$       call Cut_by_Magnitude(BFCatt, 23.e0_double)
!!$       if(Analyse_with_Physical_Sizes) then
!!$          call Monte_Carlo_Redshift_Sampling_Catalogue(BFCatt)
!!$       end if
!!$       call Cut_By_PhotoMetricRedshift(BFCatt, 0.21e0_double) !--Cut out foreground-
!!$       call Cut_By_PixelSize(BFCatt, 0.e0_double, 25.e0_double) !!!!!!!!!!!!!!!!!!!!!

       call DM_Profile_Variable_Posteriors_CircularAperture(Catt, Clusters_In%Position, Aperture_Radius, returned_Cluster_Posteriors, Distribution_Directory = Dist_Directory, Blank_Field_Catalogue = BFCatt)
    else
       !--If no Blank Field Information, then attempt a read in--!
       print *, 'Bayesian_Routines_Output_Directory'
       call DM_Profile_Variable_Posteriors_CircularAperture(Catt, Clusters_In%Position, Aperture_Radius, returned_Cluster_Posteriors, Distribution_Directory = Dist_Directory)
    end if

    call catalogue_destruct(Catt); call catalogue_destruct(BFCatt)

    print *, 'Finished Single Run Normally'

  end subroutine Mass_Estimate_Single_Run

  !######################################################!
  !-------Mock Catalogue Multiple Run Routines-----------!
  !######################################################!
  subroutine run_Mock_Production_Script(Mock_Output, nSources, frac_z, Cluster_Parameter_Filename)
    character(*), intent(in):: Mock_Output, Cluster_Parameter_Filename
    integer,intent(in):: nSources
    real(double),intent(in)::frac_z

    character(120)::Script_Name = './Mock_Catalogue_Production.sh'
    character(500)::Arguement_String
    integer::callcount

    write(*,'(A)')       '###############################################################################################################################################################################'
    write(*, '(A,I3,A)') '##################################################### Running Mock for the ', callcount, ' time ###############################################################################'
    write(*,'(A)')       '###############################################################################################################################################################################'

    write(Arguement_String, '(A,x,I6,x,e12.5,x,A)') trim(adjustl(Mock_Output)), nSources, frac_z, trim(adjustl(Cluster_Parameter_Filename))
    print *, 'Running Mock:', trim(Script_Name)//' '//trim(Arguement_String)
    call system(trim(Script_Name)//' '//trim(Arguement_String))

    write(*,'(A)')       '###############################################################################################################################################################################'
    write(*, '(A,I3,A)') '################################################################### End of ', callcount, ' Mock ###############################################################################'
    write(*,'(A)')       '###############################################################################################################################################################################'


  end subroutine run_Mock_Production_Script

  subroutine create_ClusterFile(File, nCluster, SMD_Profile, Position, FreeParameters, Redshift)
    character(*),intent(in)::File
    integer,intent(in)::nCluster, SMD_Profile
    real(double),intent(in)::Position(:), FreeParameters(:), Redshift(:)

    character(20):: fmtString

    !--Error Catching--!
    if(size(Position) /= 2*nCluster) STOP 'create_ClusterFile - Position Entered is not the correct size'
    if(size(FreeParameters) /= nCluster) STOP 'create_ClusterFile - FreeParameters Entered is not the correct size'
    if( ((size(Redshift) == 1) .or. (size(Redshift) == nCluster)) == .false. ) STOP 'create_ClusterFile - Redshift Entered is not the correct size' 

    open(unit =  99, file = trim(adjustl(File)))
    !--Header--!
    write(99, '(A)') '# Cluster Parameter file constructed automatically #'
    !--1st Namelist--!
    write(99, '(A)') '&Number'
    write(99,'(A,I1)') 'nClusters = ',nCluster
    write(99,'(A,I1)') 'SMD_Profile_Type = ',SMD_Profile
    write(99, '(A)') '/'
    write(99,*) 
    !--2nd Namelist--!
    write(99, '(A)') '&Cluster_Variables'
    write(fmtString, '(I2)') size(Position)
    write(99,'(A,'//trim(fmtString)//'(e14.7,x))') 'Positions = ',Position
    write(fmtString, '(I2)') size(FreeParameters)
    write(99,'(A,'//trim(fmtString)//'(e14.7,x))') 'DM_Parameters = ',FreeParameters
    write(fmtString, '(I2)') size(Redshift)
    write(99,'(A,'//trim(fmtString)//'(e14.7,x))') 'Redshift = ',Redshift
    write(99, '(A)') '/'

    close(99)

  end subroutine create_ClusterFile

  subroutine edit_File_String(File, SearchString, SetString)
    character(*), intent(in):: File, SearchString, SetString

    character(500)::command

    !--Copied From Run.py in FM, Untested--!

!!$    command = 's/'//trim(adjustl(SearchString))//' = .*'//'/'//trim(adjustl(SearchString))//' = '//trim(adjustl(SetString))//'/g'
!!$    command = 'sed \''//trim(adjustl(command))//'\' <'//trim(adjustl(File))//'> new.ini'
!!$    print *, 'Calling Command (read):', command
!!$    call system(command); call system('mv new.ini '//trim(adjustl(File)))

  end subroutine edit_File_String


end program Bayesian_DM_Profile_Constraints
