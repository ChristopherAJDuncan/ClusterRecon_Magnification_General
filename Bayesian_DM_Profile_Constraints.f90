program Bayesian_DM_Profile_Constraints
  !--Uses Bayesian Means to Calculate Posterior for the Dark Matter Profile variables--!                                                                   
  use Param_Types; use Catalogues; use Bayesian_Routines; use Foreground_Clusters
  implicit none

  character(200), parameter::Output_Directory_Default = 'Bayesian_DM_Profile_Constraints_Output/'
   character(200):: Output_Directory = Output_Directory_Default
   character(200)::Cluster_Filename = 'Cluster_Parameters/STAGES.ini'
   integer,allocatable:: Catalogue_Identifier(:)
   integer,allocatable::Blank_Field_Catalogue_Identifier(:)
   character(200),allocatable:: Bias_Output_Directory(:)
   integer:: Bin_By_Magnitude_Type = 2 !-1:Absolute Magnitude, 2: Apparent Magnitude-!
   integer:: Run_Type = 1 !-1:SingleRun, 2:Bias And Error-!

   type(Foreground):: Clusters
!   real(double)::Lens_Redshift = 0.165e0_double

   !--Command Line Argument Entry--!                
   integer::narg, i
   character(120)::arg

   real(double),allocatable::Cluster_Posteriors(:,:,:)

   real(double),allocatable:: Cluster_Aperture_Radius(:)

   logical:: Here

   !--------------------------------------------Read in by Command Line-------------------------------------------!              
   !##############################################################################################################!                                         
   allocate(Catalogue_Identifier(1)); Catalogue_Identifier = 5
   allocate(Blank_Field_Catalogue_Identifier(1)); Blank_Field_Catalogue_Identifier = -5
   Analyse_with_Physical_Sizes = .true.
   narg = iargc()
   if(narg == 0) then
      !-Defaults-!                                                                                                                                          
      Output_Directory = trim(Output_Directory_Default)
      write(*,'(A)') 'No arguments entered, enter "Catalogue_Identifier"; "Output_Directory"; "Use Physical Sizes"; "Run_Type"; "Blank_Field_Catalogue_Identifier", "Cluster Parameter FileName (longhand!)" in that order'
      write(*,'(A)') 'Note:: Catalogue_Identifier:: 1:STAGES_shear, 2: COMBO matched, 3:RRG, 4:STAGES-Mock, 5:COMBO-Mock; -1:Mock_COMBO and Mock_STAGES (only for Run_Type = 2)'
      write(*,'(A)') '    :: Run_Type:: 1: Single_Run, 2:Bias and Error (MultipleRun), 3: Bias and Error (STAGES Clusters), 4: Bias and Error (1 Cluster, Varying Mass)'
      write(*,'(A)') ' '
      write(*,'(A)') 'Press <ENTER> to use defaults:'
      print *, 'Catalogue Identifier:', Catalogue_Identifier, ';'
      print *, 'Output_Directory: ', trim(Output_Directory)
      print *, 'Use Physical Sizes: ', Analyse_with_Physical_Sizes
      print *, 'Run Type: ', Run_Type
      print *, 'Blank Field Identifier:', Blank_Field_Catalogue_Identifier
      print *, 'Cluster Filename:', Cluster_Filename
      read(*,*)
   else
      print *, 'Number of arguements:', narg
      do i =1, narg
         call get_command_argument(i,arg)
         select case(i)
         case(1) !-Catalogue_Identifier-!                                                                                                                
            read(arg, *) Catalogue_Identifier(1)
            Catalogue_Identifier = Catalogue_Identifier(1)
            if(Catalogue_Identifier(1) == 0) then
               deallocate(Catalogue_Identifier); allocate(Catalogue_Identifier(2));
               Catalogue_Identifier = (/4,5/)
            end if
        case(2)
           Output_Directory = trim(adjustl(arg))
        case(3)
           read(arg, *) Analyse_with_Physical_Sizes
        case(4)
           read(arg, '(I1)') Run_Type
        case(5)
            read(arg, '(I3)') Blank_Field_Catalogue_Identifier(1)
            Blank_Field_Catalogue_Identifier = Blank_Field_Catalogue_Identifier(1)
            if(Blank_Field_Catalogue_Identifier(1) == 0) then
               deallocate(Blank_Field_Catalogue_Identifier); allocate(Blank_Field_Catalogue_Identifier(2));
               Blank_Field_Catalogue_Identifier = (/-4,-5/)
            end if
         case(6)
            Cluster_Filename = trim(adjustl(arg))
        case default
           STOP 'INCORRECT NUMBER OF ARGUEMENTS ENTERED'
        end select
     end do
     if(Cluster_Filename(1:1) /= '/') STOP 'Input Cluster Filename MUST be referenced using an absolute address, stopping...'
  end if



   call get_Clusters(Clusters, Cluster_Filename)
   allocate(Cluster_Aperture_Radius(size(Clusters%Position,1))); Cluster_Aperture_Radius = 1.e0_double/60.e0_double !-Degrees-! Default - 1/60
   call distance_between_Clusters(Clusters%Position, Clusters%Redshift(1))



  !--Set the Blank Field Catalogue to the Data Catalogue if invalid value (>0) entered--!
  where((Blank_Field_Catalogue_Identifier > 0))
     Blank_Field_Catalogue_Identifier = Catalogue_Identifier
  end where

  if(Run_Type ==3) then
     !--Set Defaults for this Run Type--!
     Blank_Field_Catalogue_Identifier = Catalogue_Identifier
     !--This routine tests the affect of not having a blank feild (i.e. knowledge of the iontrinsic distributions, only taken from lensed fields.) Therefore the intrinsic distributions should come from the lensed data--!
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
     call Mass_Estimate_Single_Run(Catalogue_Identifier(1), Output_Directory, Cluster_Posteriors, Blank_Field_Catalogue_Identifier(1),Clusters, Cluster_Aperture_Radius)
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
    integer::nParam = 5
    real(double):: Param_Min = 1.0, Param_Max = 4.0

    integer:: nRun

    character(120):: Run_Directory(size(Directory))

    real(double),allocatable:: Bias_Mode(:,:), Mode_Variance(:,:) !-Catalogue, Aperture-!
    real(double),allocatable:: Param_Bias_Mode(:,:,:), Param_Mode_Variance(:,:,:) !-Parameter, Catalogue, Aperture-!

    real(double)::Aperture_Mass, eAperture_Mass, Aperture_Mass_Input, Discardable

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

    !--Get The STAGES Clsuters Biases--!
    call get_Clusters(Cluster, '/disk1/cajd/Size_Magnification/Cluster_Parameters/STAGES.ini')
    nParam = 1; nClusters = size(Cluster%Position,1)
    allocate(Position(2*nClusters))
    do i =1, size(Clusters%Position,1)
       Position(2*i-1) = Clusters%Position(i,1); Position(2*i) = Clusters%Position(i,2)
    end do
!    Position = (/ 149.1099,-9.9561, 148.9889,-9.9841, 149.1424,-10.1666, 148.9101,-10.1719 /)
    allocate(Param_Bias_Mode(nParam, size(Cat_Ident), nClusters)); Param_Bias_Mode =0.e0_double
    allocate(Param_Mode_Variance(nParam, size(Cat_Ident), nClusters)); Param_Mode_Variance =0.e0_double
    allocate(ParameterValues(nParam,nClusters)); ParameterValues = 0.e0_double
    !---Create the Parameter File----!
    ParameterValues(1,:) = Clusters%DM_Profile_Parameter
    call create_ClusterFile(Mock_Cluster_Filename, nClusters, SMD_Type, Position, ParameterValues(1,:), (/Redshift/))

    !--Read in Cluster--!
    call get_Clusters(Cluster, Mock_Cluster_Filename)
    allocate(ICluster_Aperture_Radius(size(Cluster%Position,1))); ICluster_Aperture_Radius = 1.e0_double/60.e0_double !-Degrees-! Default - 1/60

    !--Run the Bias Routine---!
    !--Bias_Mode and Mode_Variacne contain the variance in the DM model free parameter--!

    call Posterior_Maximum_Likelihood_Bias_Error(Cat_Ident, Run_Directory, Blank_Field_Cat_Ident, Mock_Cluster_Filename, iCluster_Aperture_Radius, Bias_Mode, Mode_Variance)
    
    if( (size(Param_Bias_Mode(1,:,:)) /= size(Bias_Mode)) ) STOP 'Cluster_Bias_Outputs - Error in Bias Mode size'
    if( (size(Param_Mode_Variance(1,:,:)) /= size(Mode_Variance)) ) STOP 'Cluster_Bias_Outputs - Error in Mode_Variance size'
    Param_Bias_Mode(1,:,:) = Bias_Mode; Param_Mode_Variance(1,:,:) = Mode_Variance
    deallocate(Bias_Mode, Mode_Variance)
       
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
             call Integrated_Mass_Within_Radius(SMD_Type, ICluster_Aperture_Radius(Ap), ParameterValues(j,Ap), 0.e0_double, Aperture_Mass_Input, Discardable, Redshift)
             call Integrated_Mass_Within_Radius(SMD_Type, ICluster_Aperture_Radius(Ap), Param_Bias_Mode(j,cid,Ap), Param_Mode_Variance(j,cid,Ap), Aperture_Mass, eAperture_Mass, Redshift)

             write(36, '(6(e12.5,x))') ParameterValues(j,Ap), Param_Bias_Mode(j,cid,Ap), Param_Mode_Variance(j,cid,Ap), Aperture_Mass_Input, Aperture_Mass,  eAperture_Mass 
          end do
          write(36,*)
       end do
    end do
    close(36)

    !--Plot Result--!
    
    deallocate(Param_Bias_Mode, Param_Mode_Variance, Position, ParameterValues, iCluster_Aperture_Radius)
    print *, 'Finished Mock Cluster Bias Run STAGES normally'

    !--nClusters of equal mass, varying r200--!
    nClusters = 1
    allocate(Position(2*nClusters))
    Position = (/ 149.1099,-9.9561 /)
    allocate(Param_Bias_Mode(nParam, size(Cat_Ident), nClusters)); Param_Bias_Mode =0.e0_double
    allocate(Param_Mode_Variance(nParam, size(Cat_Ident), nClusters)); Param_Mode_Variance =0.e0_double

    allocate(ParameterValues(nParam,nClusters)); ParameterValues = 0.e0_double
    allocate(ICluster_Aperture_Radius(nClusters)); ICluster_Aperture_Radius = 1.e0_double/60.e0_double !-Degrees-! Default - 1/60
    do i =1, nParam
       !---Create the Parameter File----!
       ParameterValues(i,:) = Param_Min + (i-1)*((Param_Max-Param_min)/(nParam-1))
       call create_ClusterFile(Mock_Cluster_Filename, nClusters, SMD_Type, Position, ParameterValues(i,:), (/Redshift/))
       !--Run the Bias Routine---!
       !--Bias_Mode and Mode_Variacne contain the variance in the DM model free parameter--!
       print *, 'Calling Bias routine'
       call Posterior_Maximum_Likelihood_Bias_Error(Cat_Ident, Run_Directory, Blank_Field_Cat_Ident, Mock_Cluster_Filename, iCluster_Aperture_Radius, Bias_Mode, Mode_Variance)
       print *, 'Called Bias Routine'

       if( (size(Param_Bias_Mode(1,:,:)) /= size(Bias_Mode)) ) STOP 'Cluster_Bias_Outputs - Error in Bias Mode size'
       if( (size(Param_Mode_Variance(1,:,:)) /= size(Mode_Variance)) ) STOP 'Cluster_Bias_Outputs - Error in Mode_Variance size'
       Param_Bias_Mode(i,:,:) = Bias_Mode; Param_Mode_Variance(i,:,:) = Mode_Variance
       deallocate(Bias_Mode, Mode_Variance)
    end do
       
    !---Output Result (Input, Output, Variance)---!
    !--Input Parameter, Output Parameter, Variance, Input Ap_Mass, Output Ap_Mass, Variance--! 
    !--Space between each apeture--!
    do cid = 1, size(Run_Directory)!size(Cat_Ident)
       open(36, file = trim(Run_Directory(cid))//'Mass_Bias_duetoForeground.dat')
       !-Header-!
       write(36, '(A)') '## Input_Parameter, Output_Parameter, Variance, Input Aperture Mass, Output Aperture Mass, Varaince'
       write(36, '(A)') '## Space between apertures'
       write(fmtString, '(I1)') size(ICluster_Aperture_Radius); write(36, '(A,'//trim(fmtString)//'(e12.5,x))') '## Aperture (Degrees):', ICluster_Aperture_Radius !--Edit fmt to array
       write(36, '(A,A)') '## Mass Profile:', Profile_String !--Edit fmt to array  
     
       do Ap = 1, nClusters
          do j =1, nParam
             call Integrated_Mass_Within_Radius(SMD_Type, ICluster_Aperture_Radius(Ap), ParameterValues(j,Ap), 0.e0_double, Aperture_Mass_Input, Discardable, Redshift)
             call Integrated_Mass_Within_Radius(SMD_Type, ICluster_Aperture_Radius(Ap), Param_Bias_Mode(j,cid,Ap), Param_Mode_Variance(j,cid,Ap), Aperture_Mass, eAperture_Mass, Redshift)

             write(36, '(6(e12.5,x))') ParameterValues(j,Ap), Param_Bias_Mode(j,cid,Ap), Param_Mode_Variance(j,cid,Ap), Aperture_Mass_Input, Aperture_Mass,  eAperture_Mass 
          end do
          write(36,*)
       end do
    end do
    close(36)
    
    deallocate(Param_Bias_Mode, Param_Mode_Variance, Position, ICluster_Aperture_Radius)
    print *, 'Finished Mock Cluster Bias Run'

    !--Plot Result--!
    do cid = 1, size(Run_Directory)
       call system('python Cluster_Bias_Plotter.py '//trim(Run_Directory(cid))//'Mass_Bias_duetoForeground_STAGES.dat'//' '//trim(Run_Directory(cid))//'Mass_Bias_duetoForeground.dat')
    end do

  end subroutine Cluster_Bias_Outputs

  subroutine Posterior_Maximum_Likelihood_Bias_Error(Cat_Ident, Directory, Blank_Field_Cat_Ident, Input_Clusters_Filename, Aperture_Radius, Bias_Mode_Out, Mode_Variance_Out)
    use Cosmology; use Statistics, only:variance_discrete; use Matrix_methods, only:Subset
    !--Produces multiple Posteriors for many different mock catalogue realisations, to calculate the ML-point bias and variance--!
    integer, intent(in)::Cat_Ident(:)
    character(*), intent(in)::Directory(:)
    integer, intent(in)::Blank_Field_Cat_Ident(:)
    character(*),intent(in):: Input_Clusters_Filename
    real(double), intent(in)::Aperture_Radius(:)
    real(double),allocatable,optional::Bias_Mode_Out(:,:), Mode_Variance_Out(:,:) !-Catalogue, Ap-!

    real(double),allocatable::Bias_Mode(:), Mode_Variance(:)

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

    !----Mock Parameters-------!
    integer:: nSources = 70000
    real(double)::frac_z = 1.e0_double
    character(120)::Mock_Output = 'Output/' !-Generalise?-!


    !--Temporary Storage of User-entered variable--!
    tPhysical_Size = Analyse_with_Physical_Sizes

!    if(all( (/4,5/) /= Cat_Ident)) STOP 'Posterior_Maximum_Likelihood_Bias_Error - Invalid Catalogue Indenifier entered, it must correspond to a Mock Catalogue'
    if(SubSet(Cat_Ident, (/4,5/)) ==.false.)STOP' Posterior_Maximum_Likelihood_Bias_Error - Catalogue Identiers entered are not valid'

    if(size(Directory) /= size(Cat_Ident)) STOP' Posterior_Maximum_Likelihood_Bias_Error - Directory and Cat Identifier not of the same size, exiting'

    

    do nR = 1, nRun
       write(*,'(A)') '!##################################################################!'
       write(*,'(A, I3)') '!------------------------------------------------------Run', nR

       !- Run Mock Catalogue Production Script -!
       call run_Mock_Production_Script(Mock_Output, nSources, frac_z, Input_Clusters_Filename)

       call get_Clusters(iClusters, Input_Clusters_Filename)
       nAp = size(iClusters%Position,1)

       !--Initialise Output--!
       if(nR==1) then
          if(present(Mode_Variance_Out)) then
             if(allocated(Mode_Variance_Out)) deallocate(Mode_Variance_Out)
             allocate(Mode_Variance_Out(size(Cat_Ident), nAp)); Mode_Variance_Out = 0.e0_double
          end if
          if(present(Bias_Mode_Out)) then
             if(allocated(Bias_Mode_Out)) deallocate(Bias_Mode_Out)
             allocate(Bias_Mode_Out(size(Cat_Ident), nAp)); Bias_Mode_Out = 0.e0_double
          end if
       end if
       
       do Id = 1, size(Cat_Ident)
          !-Check for directory existence-!
          if(nR==1) then
             write(Run_Output_Directory, '(I1)') Cat_Ident(ID)
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
          

          if(Cat_Ident(ID)==5) Analyse_with_Physical_Sizes = .true. !-Ensure Redshift Information used at Data end if available for COMBO-!
          call Mass_Estimate_Single_Run(Cat_Ident(ID), Run_Output_Directory, Single_Run_Posterior, Blank_Field_Cat_Ident(ID), iClusters, Aperture_Radius)
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
          call Combine_Posteriors(Posteriors(ID,1,Ap,1,:), Posteriors(ID,:,Ap,2,:), .false., Combined_Posterior(Ap,2,:))
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
       allocate(Mode_Variance(nAp)); Mode_Variance = 0.e0_double
       open(unit = 51, file = trim(Directory(ID))//'BiasesAndErrorVariance.dat')
       write(51, '(A)') '#Following is Mode of Combined Posterior and Variance of ML points of each run'
       write(51, '(A)') '#Difference between mode of Combined to known True Value gives indication of Bias'
       write(51, '(A)') '#Variance of ML points gives indication on the relyability of the Errors from a single-run Bayesian Analysis'
       
       write(*,'(A)') '!##################################################################!'
       write(*,'(A)') '!##################################################################!'
       write(*,'(A)') ' Finished Combination, following is Mode (shift from known is bias) and variance (Error on Single Run Posterior Errors) of DM profile free parameter:'
       do Ap = 1, nAp
          D_l =  angular_diameter_distance_fromRedshift(0.e0_double, iClusters%Redshift(Ap))
          Area = 3.142e0_double*(D_l*((3.142e0_double*Cluster_Aperture_Radius(Ap))/180.e0_double))**2.e0_double
       
          call Posterior_Statistics(Combined_Posterior(Ap,1,:), Combined_Posterior(Ap,2,:), ModeVal = Bias_Mode(Ap))
          Mode_Variance(Ap) = variance_discrete(ML_Point(Ap,:), ML_Point(Ap,:))

          write(Bias_String, '(e10.4)') Bias_Mode(Ap)
          write(Error_String, '(e10.4)') Mode_Variance(Ap)

          print *, 'Setting Bias Mode Out'
          if(present(Bias_Mode_Out)) Bias_Mode_Out(ID, Ap) = Bias_Mode(Ap)
          if(present(Mode_Variance_Out)) Mode_Variance_Out(ID, Ap) = Mode_Variance(Ap)

          write(51,*) 'Cluster ', Ap, ' has Mode:', Bias_String
          write(51,*) '     and variance:', Error_String
          
          write(*,*) 'Cluster ', Ap, ' has Mode:', Bias_String
          write(*,*) '     and variance:', Error_String
       end do
       write(*,'(A)') '!##################################################################!'
       write(*,'(A)') '!##################################################################!'
       close(51)
       deallocate(Bias_mode, Combined_Posterior, ML_Point, Mode_Variance)
    end do
    deallocate(Posteriors)
    call Foreground_Destruct(iClusters)

  end subroutine Posterior_Maximum_Likelihood_Bias_Error


  subroutine Mass_Estimate_Single_Run(Cat_Ident, run_Output_Dir, returned_Cluster_Posteriors, Blank_Field_Cat_Ident, Clusters_In, Aperture_Radius)
    use MC_Redshift_Sampling; use Bayesian_Routines
    integer,intent(in)::Cat_Ident
    real(double),allocatable,intent(out):: returned_Cluster_Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-!
    integer,intent(in)::Blank_Field_Cat_Ident
    type(Foreground):: Clusters_In
    real(double),intent(in):: Aperture_Radius(:)

    character(*), intent(in):: run_Output_Dir

    type(Catalogue)::Catt, BFCatt
    character(120)::Catalogue_Directory, Catalogue_Filename
    integer,dimension(:),allocatable::Catalogue_Cols
    logical::here

    inquire(directory = trim(run_Output_Dir), exist = here)
    if(here == .false.) call system('mkdir '//trim(adjustl(run_Output_Dir)))
    Bayesian_Routines_Output_Directory = trim(run_Output_Dir)



    !## 1: STAGES shear, 2: COMBO17, 3:RRG, 4: Mocks_STAGES; 5:Mocks_COMBO!          
    call common_Catalogue_directories(Cat_Ident, Catalogue_Directory, Catalogue_Filename, Catalogue_Cols)
    call catalogue_readin(Catt, trim(adjustl(Catalogue_Directory))//trim(adjustl(Catalogue_Filename)), 'Tr(J)', Catalogue_Cols)

    call Cut_by_Magnitude(Catt, 23.e0_double) !-Taken from CH08 P1435-!   
    if(Analyse_with_Physical_Sizes) then
       call Monte_Carlo_Redshift_Sampling_Catalogue(Catt)
    end if
    call Cut_By_PhotoMetricRedshift(Catt, 0.21e0_double) !--Cut out foreground--!                                                                            
    if(Blank_Field_Cat_Ident == Cat_Ident) then
       BFCatt = Catt 
    else
       call common_Catalogue_directories(Blank_Field_Cat_Ident, Catalogue_Directory, Catalogue_Filename, Catalogue_Cols)
       call catalogue_readin(BFCatt, trim(adjustl(Catalogue_Directory))//trim(adjustl(Catalogue_Filename)), 'Tr(J)', Catalogue_Cols)

       call Cut_by_Magnitude(BFCatt, 23.e0_double)
       if(Analyse_with_Physical_Sizes) then
          call Monte_Carlo_Redshift_Sampling_Catalogue(BFCatt)
       end if
       call Cut_By_PhotoMetricRedshift(BFCatt, 0.21e0_double) !--Cut out foreground-
    end if

    call DM_Profile_Variable_Posteriors_CircularAperture(Catt, Clusters_In%Position, Aperture_Radius, returned_Cluster_Posteriors, Blank_Field_Catalogue = BFCatt)

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
