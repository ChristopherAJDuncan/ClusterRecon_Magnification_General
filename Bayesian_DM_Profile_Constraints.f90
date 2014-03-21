program Bayesian_DM_Profile_Constraints
  !--Uses Bayesian Means to Calculate Posterior for the Dark Matter Profile variables--!                                                                   
  use Param_Types; use Catalogues; use Bayesian_Routines; use Clusters
  implicit none

  character(200), parameter::Output_Directory_Default = 'Bayesian_DM_Profile_Constraints_Output/'
   character(200):: Output_Directory = Output_Directory_Default
   integer,allocatable:: Catalogue_Identifier(:)
   integer,allocatable::Blank_Field_Catalogue_Identifier(:)
   character(200),allocatable:: Bias_Output_Directory(:)
   integer:: Bin_By_Magnitude_Type = 2 !-1:Absolute Magnitude, 2: Apparent Magnitude-!
   integer:: Run_Type = 1 !-1:SingleRun, 2:Bias And Error-!

   real(double)::Lens_Redshift = 0.165e0_double

   !--Command Line Argument Entry--!                
   integer::narg, i
   character(120)::arg

   real(double),allocatable::Cluster_Posteriors(:,:,:)

   real(double),allocatable::Cluster_Pos(:,:), Cluster_Aperture_Radius(:)

   logical:: Here

   allocate(Cluster_Pos(4,2)); Cluster_Pos = 0.e0_double
   Cluster_Pos(1,:) = (/149.1099e0_double,-9.9561e0_double/)
   Cluster_Pos(2,:) = (/148.9889e0_double,-9.9841e0_double/)
   Cluster_Pos(3,:) = (/149.1424e0_double,-10.1666e0_double/)
   Cluster_Pos(4,:) = (/148.9101e0_double,-10.1719e0_double/)
   allocate(Cluster_Aperture_Radius(4)); Cluster_Aperture_Radius = 1.e0_double/60.e0_double !-Degrees-! Default - 1/60

   call distance_between_Clusters(Cluster_Pos, Lens_Redshift)

   !--------------------------------------------Read in by Command Line-------------------------------------------!              
   !##############################################################################################################!                                         
   allocate(Catalogue_Identifier(1)); Catalogue_Identifier = 5
   allocate(Blank_Field_Catalogue_Identifier(1)); Blank_Field_Catalogue_Identifier = 5
   Analyse_with_Physical_Sizes = .true.
   narg = iargc()
   if(narg == 0) then
      !-Defaults-!                                                                                                                                          
      Output_Directory = trim(Output_Directory_Default)
      write(*,'(A)') 'No arguments entered, enter "Catalogue_Identifier"; "Output_Directory"; "Use Physical Sizes"; "Run_Type"; "Blank_Field_Catalogue_Identifier" in that order'
      write(*,'(A)') 'Note:: Catalogue_Identifier:: 1:STAGES_shear, 2: COMBO matched, 3:RRG, 4:STAGES-Mock, 5:COMBO-Mock; -1:Mock_COMBO and Mock_STAGES (only for Run_Type = 2)'
      write(*,'(A)') '    :: Run_Type:: 1: Single_Run, 2:Bias and Error (MultipleRun), 3: Bias and Error (STAGES Clusters), 4: Bias and Error (1 Cluster, Varying Mass)'
      write(*,'(A)') ' '
      write(*,'(A)') 'Press <ENTER> to use defaults:'
      print *, 'Catalogue Identifier:', Catalogue_Identifier, ';'
      print *, 'Output_Directory: ', trim(Output_Directory)
      print *, 'Use Physical Sizes: ', Analyse_with_Physical_Sizes
      print *, 'Run Type: ', Run_Type
      print *, 'Blank Field Identifier:', Blank_Field_Catalogue_Identifier
      read(*,*)
   else
      print *, 'Number of arguements:', narg
      do i =1, narg
         call get_command_argument(i,arg)
         select case(i)
         case(1) !-Catalogue_Identifier-!                                                                                                                
            read(arg, *) Catalogue_Identifier(1)
            Catalogue_Identifier = Catalogue_Identifier(1)
            if(Catalogue_Identifier(1) == -1) then
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
            read(arg, *) Blank_Field_Catalogue_Identifier(1)
            Blank_Field_Catalogue_Identifier = Blank_Field_Catalogue_Identifier(1)
            if(Blank_Field_Catalogue_Identifier(1) == -1) then
               deallocate(Blank_Field_Catalogue_Identifier); allocate(Blank_Field_Catalogue_Identifier(2));
               Blank_Field_Catalogue_Identifier = (/-4,-5/)
            end if
            read(arg, '(I1)') Blank_Field_Catalogue_Identifier
        case default
           STOP 'INCORRECT NUMBER OF ARGUEMENTS ENTERED'
        end select
     end do
  end if
  !--Set the Blank Field Catalogue to the Data Catalogue if invalid value (>0) entered--!
  where((Blank_Field_Catalogue_Identifier > 0))
     Blank_Field_Catalogue_Identifier = Catalogue_Identifier
  end where

  if(RunType ==3) then
     !--Set Defaults for this Run Type--!
     Blank_Field_Catalogue_Identifier = -Catalogue_Identifier
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
  print *, ' '

  inquire(directory = trim(Output_Directory), exist = here)
  if(here == .false.) call system('mkdir '//trim(Output_Directory))

  select case(Run_Type)
  case(1)
     print *, '!---- Producing a single Mass Estimate:'
     call Mass_Estimate_Single_Run(Catalogue_Identifier(1), Output_Directory, Cluster_Posteriors, Blank_Field_Catalogue_Identifier(1))
  case(2)
     !--Set Output Directory Names--!
     allocate(Bias_Output_Directory(size(Catalogue_Identifier))); Bias_Output_Directory = ' '
     do i = 1, size(Catalogue_Identifier)
        write(Bias_Output_Directory(i), '(I1)') Catalogue_Identifier(i)
        Bias_Output_Directory(i) = trim(Output_Directory)//'BiasRun_ID'//trim(Bias_Output_Directory(i))//'/'
     end do
     call Posterior_Maximum_Likelihood_Bias_Error(Catalogue_Identifier, Bias_Output_Directory, Blank_Field_Catalogue_Identifier)
  case default
     STOP 'Run_type Entered not supported'
  end select

contains

  subroutine Posterior_Maximum_Likelihood_Bias_Error(Cat_Ident, Directory, Blank_Field_Cat_Ident)
    use Cosmology; use Statistics, only:variance_discrete; use Matrix_methods, only:Subset
    !--Produces multiple Posteriors for many different mock catalogue realisations, to calculate the ML-point bias and variance--!
    integer, intent(in)::Cat_Ident(:)
    character(*), intent(in)::Directory(:)
    integer, intent(in)::Blank_Field_Cat_Ident(:)

    character(200)::Run_Output_Directory
    integer::nR
    integer::nRun = 20

    integer::Ap, nAp, Id

    real(double),allocatable::Single_Run_Posterior(:,:,:) !-Aperture, Grid/Posterior, Value-!
    real(double),allocatable:: Posteriors(:,:,:,:,:) !-CatalogueID, Run, Aperture, Grid/Posterior, Value-!
    real(double),allocatable:: Combined_Posterior(:,:,:) !-Aperture, Grid/Posterior, Value-!
    real(double),allocatable:: ML_Point(:,:) !-Ap, Point-!

    !--Conversion Declarations--!
    real(double)::D_l, Area
    character(12)::Error_String, Bias_String
    real(double),allocatable::Bias_Mode(:)

    character(7)::fmtstring

    logical::here

    logical::tPhysical_Size

    !--Temporary Storage of User-entered variable--!
    tPhysical_Size = Analyse_with_Physical_Sizes

!    if(all( (/4,5/) /= Cat_Ident)) STOP 'Posterior_Maximum_Likelihood_Bias_Error - Invalid Catalogue Indenifier entered, it must correspond to a Mock Catalogue'
    if(SubSet(Cat_Ident, (/4,5/)) ==.false.)STOP' Posterior_Maximum_Likelihood_Bias_Error - Catalogue Identiers entered are not valid'

    if(size(Directory) /= size(Cat_Ident)) STOP' Posterior_Maximum_Likelihood_Bias_Error - Directory and Cat Identifier not of the same size, exiting'

    nAp = size(Cluster_Pos,1)

    do nR = 1, nRun
       write(*,'(A)') '!##################################################################!'
       write(*,'(A, I3)') '!------------------------------------------------------Run', nR

       !- Run Mock Catalogue Production Script -!
       call system('./Mock_Catalogue_Production.sh')

       do Id = 1, size(Cat_Ident)
          !-Check for directory existence-!
          inquire(directory = trim(Directory(ID)), exist = here)
          if(here == .false.) call system('mkdir '//trim(Directory(ID)))

          if(nR<10) then
             write(Run_Output_Directory, '(I1)') nR
          else
             write(Run_Output_Directory, '(I2)') nR
          end if
          Run_Output_Directory = trim(Directory(ID))//trim(Run_Output_Directory)//'/'
          

          if(Cat_Ident(ID)==5) Analyse_with_Physical_Sizes = .true. !-Ensure Redshift Information used at Data end if available for COMBO-!
          call Mass_Estimate_Single_Run(Cat_Ident(ID), Run_Output_Directory, Single_Run_Posterior, Blank_Field_Cat_Ident(ID))
          Analyse_with_Physical_Sizes = tPhysical_Size

          if(nR ==1 .and. ID == 1) then
             allocate(Posteriors(size(Cat_Ident),nRun, size(Single_Run_Posterior,1), size(Single_Run_Posterior,2), size(Single_Run_Posterior,3))); Posteriors = 0.e0_double
          end if
          Posteriors(ID,nR,:,:,:) = Single_Run_Posterior
          
          deallocate(Single_Run_Posterior)
       end do
    end do

!!!!             allocate(ML_Point(size(Single_Run_Posterior,1), nRun)); ML_Point = 0.e0_double
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
       open(unit = 51, file = trim(Directory(ID))//'BiasesAndErrorVariance.dat')
       write(51, '(A)') '#Following is Mode of Combined Posterior and Variance of ML points of each run'
       write(51, '(A)') '#Difference between mode of Combined to known True Value gives indication of Bias'
       write(51, '(A)') '#Variance of ML points gives indication on the relyability of the Errors from a single-run Bayesian Analysis'
       
       write(*,'(A)') '!##################################################################!'
       write(*,'(A)') '!##################################################################!'
       write(*,'(A)') ' Finished Combination, following is Mode (shift from known is bias) and variance (Error on Single Run Posterior Errors) of DM profile free parameter:'
       do Ap = 1, nAp
          D_l =  angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
          Area = 3.142e0_double*(D_l*((3.142e0_double*Cluster_Aperture_Radius(Ap))/180.e0_double))**2.e0_double
       
          call Posterior_Statistics(Combined_Posterior(Ap,1,:), Combined_Posterior(Ap,2,:), ModeVal = Bias_Mode(Ap))
          write(Bias_String, '(e10.4)') Bias_Mode(Ap)
          write(Error_String, '(e10.4)') variance_discrete(ML_Point(Ap,:), ML_Point(Ap,:))
          
          write(51,*) 'Cluster ', Ap, ' has Mode:', Bias_String
          write(51,*) '     and variance:', Error_String
          
          write(*,*) 'Cluster ', Ap, ' has Mode:', Bias_String
          write(*,*) '     and variance:', Error_String
       end do
       write(*,'(A)') '!##################################################################!'
       write(*,'(A)') '!##################################################################!'
       close(51)
       deallocate(Bias_mode, Combined_Posterior, ML_Point)
    end do
    deallocate(Posteriors)

  end subroutine Posterior_Maximum_Likelihood_Bias_Error


  subroutine Mass_Estimate_Single_Run(Cat_Ident, run_Output_Dir, returned_Cluster_Posteriors, Blank_Field_Cat_Ident)
    use MC_Redshift_Sampling; use Bayesian_Routines
    integer,intent(in)::Cat_Ident
    real(double),allocatable,intent(out):: returned_Cluster_Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-!
    integer,intent(in)::Blank_Field_Cat_Ident

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

    call DM_Profile_Variable_Posteriors_CircularAperture(Catt, Cluster_Pos, Cluster_Aperture_Radius, returned_Cluster_Posteriors, Blank_Field_Catalogue = BFCatt)

  end subroutine Mass_Estimate_Single_Run

  !######################################################!
  !-------Mock Catalogue Multiple Run Routines-----------!
  !######################################################!


end program Bayesian_DM_Profile_Constraints
