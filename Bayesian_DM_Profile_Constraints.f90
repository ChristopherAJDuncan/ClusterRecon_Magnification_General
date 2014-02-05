program Bayesian_DM_Profile_Constraints
  !--Uses Bayesian Means to Calculate Posterior for the Dark Matter Profile variables--!                                                                   
  use Param_Types; use Catalogues; use Bayesian_Routines
  implicit none

    character(200):: Output_Directory, Output_Directory_Default = 'Bayesian_DM_Profile_Constraints_Output/'
  integer:: Catalogue_Identifier = 4
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
  allocate(Cluster_Aperture_Radius(4)); Cluster_Aperture_Radius = 1.e0_double/60.e0_double !-Degrees-!

  !--------------------------------------------Read in by Command Line-------------------------------------------!              
  !##############################################################################################################!                                         
  narg = iargc()
  if(narg == 0) then
     !-Defaults-!                                                                                                                                          
     Output_Directory = trim(Output_Directory_Default)
     write(*,'(A)') 'No arguments entered, enter "Catalogue_Identifier"; "Output_Directory"; "Use Physical Sizes"; "Run_Type" in that order'
     write(*,'(A)') 'Note:: Catalogue_Identifier:: 1:STAGES_shear, 2: COMBO matched, 3:RRG, 4:STAGES-Mock, 5:COMBO-Mock'
     write(*,'(A)') '    :: Run_Type:: 1: Single_Run, 2:Bias and Error (MultipleRun)'
     write(*,'(A)') ' '
     write(*,'(A)') 'Press <ENTER> to use defaults:'
     print *, 'Catalogue Identifier:', Catalogue_Identifier, ';'
     print *, 'Output_Directory: ', trim(Output_Directory)
     print *, 'Use Physical Sizes: ', Analyse_with_Physical_Sizes
     print *, 'Run Type: ', Run_Type
     read(*,*)
  else
     print *, 'Number of arguements:', narg
     do i =1, narg
        call get_command_argument(i,arg)
        select case(i)
        case(1) !-Catalogue_Identifier-!                                                                                                                
           read(arg, '(I1)') Catalogue_Identifier
        case(2)
           Output_Directory = trim(adjustl(arg))
        case(3)
           read(arg, *) Analyse_with_Physical_Sizes
        case(4)
           read(arg, '(I1)') Run_Type
        case default
           STOP 'INCORRECT NUMBER OF ARGUEMENTS ENTERED'
        end select
     end do
  end if
  print *, 'Using Catalogue_Identifier: ', Catalogue_Identifier
  print *, 'and Output_Directory: ', trim(Output_Directory)
  print *, 'Using Physical Sizes Data-Side (with longer run-time):', Analyse_with_Physical_Sizes
  print *, 'Run_Type:', Run_Type
  print *, ' '

  select case(Run_Type)
  case(1)
     print *, '!---- Producing a single Mass Estimate:'
     call Mass_Estimate_Single_Run(Catalogue_Identifier, Output_Directory, Cluster_Posteriors)
  case(2)
     call Posterior_Maximum_Likelihood_Bias_Error(Catalogue_Identifier)
  case default
     STOP 'Run_type Entered not supported'
  end select

contains

  subroutine Posterior_Maximum_Likelihood_Bias_Error(Cat_Ident)
    use Cosmology; use Statistics, only:variance_discrete
    !--Produces multiple Posteriors for many different mock catalogue realisations, to calculate the ML-point bias and variance--!
    integer, intent(in)::Cat_Ident

    character(200)::Run_Output_Directory
    integer::nR
    integer::nRun = 20

    integer::Ap, nAp

    real(double),allocatable::Single_Run_Posterior(:,:,:) !-Aperture, Grid/Posterior, Value-!
    real(double),allocatable:: Posteriors(:,:,:,:) !-Run, Aperture, Grid/Posterior, Value-!
    real(double),allocatable:: Combined_Posterior(:,:,:) !-Aperture, Grid/Posterior, Value-!
    real(double),allocatable:: ML_Point(:,:) !-Ap, Point-!

    !--Conversion Declarations--!
    real(double)::D_l, Area
    character(12)::Error_String, Bias_String
    real(double),allocatable::Bias_Mode(:)

    character(7)::fmtstring

    if(all( (/4,5/) /= Cat_Ident)) STOP 'Posterior_Maximum_Likelihood_Bias_Error - Invalid Catalogue Indenifier entered, it must correspond to a Mock Catalogue'
    
    nAp = size(Cluster_Pos,1)

    do nR = 1, nRun
       write(*,'(A)') '!##################################################################!'
       write(*,'(A, I3)') '!------------------------------------------------------Run', nR

       !- Run Mock Catalogue Production Script -!
       call system('./Mock_Catalogue_Production.sh')

       if(nR<10) then
          write(Run_Output_Directory, '(I1)') nR
       else
          write(Run_Output_Directory, '(I2)') nR
       end if
       Run_Output_Directory = trim(Output_Directory)//trim(Run_Output_Directory)//'/'

       call Mass_Estimate_Single_Run(Cat_Ident, Run_Output_Directory, Single_Run_Posterior)

       if(nR ==1) then
          allocate(Posteriors(nRun, size(Single_Run_Posterior,1), size(Single_Run_Posterior,2), size(Single_Run_Posterior,3))); Posteriors = 0.e0_double
          allocate(ML_Point(size(Single_Run_Posterior,1), nRun)); ML_Point = 0.e0_double
       end if
       Posteriors(nR,:,:,:) = Single_Run_Posterior
   
       deallocate(Single_Run_Posterior)
    end do

    !--Combine Posteriors and Get ML Point for Combined_Posteriors--!
    
    allocate(Combined_Posterior(nAp, 2, size(Posteriors,4))); Combined_Posterior = 0.e0_double
    Combined_Posterior(:,1,:) = Posteriors(1,:,1,:)
    do Ap = 1, nAp
       do nR = 1, nRun
          call Posterior_Statistics(Posteriors(nR,Ap,1,:), Posteriors(nR,Ap,2,:), ModeVal = ML_Point(Ap,nR))
       end do
       call Combine_Posteriors(Posteriors(1,Ap,1,:), Posteriors(:,Ap,2,:), .false., Combined_Posterior(Ap,2,:))
    end do
    deallocate(Posteriors)

    !---Output---!
    open(37, file = trim(Output_Directory)//'Bias_Combined_Posterior.dat')
    write(fmtstring,'(I1)') nAp+1
    do i =1, size(Combined_Posterior,3)
       write(37,  '('//trim(adjustl(fmtstring))//'(e14.7,x))') Combined_Posterior(1,1,i), Combined_Posterior(:,2,i)
    end do
    close(37)
    print *, 'Output Combined Posterior for Bias to: ', trim(Output_Directory)//'Bias_Combined_Posterior.dat' 

    allocate(Bias_Mode(nAp)); Bias_Mode = 0.e0_double
    open(unit = 51, file = trim(Output_Directory)//'BiasesAndErrorVariance.dat')
    write(51, '(A)') '#Following is Mode of Combined Posterior and Varaince of ML points of each run'
    write(51, '(A)') '#Difference between mode of Combined to known True Value gives indication of Bias'
    write(51, '(A)') '#Variance of ML points gives indication on the relyability of the Errors from a single-run Bayesian Analysis'

    write(*,'(A)') '!##################################################################!'
    write(*,'(A)') '!##################################################################!'
    write(*,'(A)') ' Finished Combination, following is Mode (shift from known is bias) and variance (Error on Single Run Posterior Errors):'
    do Ap = 1, nAp
       D_l =  angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
       Area = 3.142e0_double*(D_l*((3.142e0_double*Cluster_Aperture_Radius(Ap))/180.e0_double))**2.e0_double
       
       print *, 'Getting Stats for Ap:', Ap
       call Posterior_Statistics(Combined_Posterior(Ap,1,:), Combined_Posterior(Ap,2,:), ModeVal = Bias_Mode(Ap))
       write(Bias_String, '(e10.4)') Area*Bias_Mode(Ap)*1.e18_double
       write(Error_String, '(e10.4)') Area*1.e18_double*variance_discrete(ML_Point(Ap,:), ML_Point(Ap,:))

       write(51,*) 'Cluster ', Ap, ' has Mode:', Bias_String
       write(51,*) '     and variance:', Error_String

       write(*,*) 'Cluster ', Ap, ' has Mode:', Bias_String
       write(*,*) '     and variance:', Error_String
    end do
    write(*,'(A)') '!##################################################################!'
    write(*,'(A)') '!##################################################################!'
    close(51)

    deallocate(Bias_mode, Combined_Posterior)

  end subroutine Posterior_Maximum_Likelihood_Bias_Error


  subroutine Mass_Estimate_Single_Run(Cat_Ident, run_Output_Dir, returned_Cluster_Posteriors)
    use MC_Redshift_Sampling; use Bayesian_Routines
    integer,intent(in)::Cat_Ident
    real(double),allocatable,intent(out):: returned_Cluster_Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-!

    character(*), intent(in):: run_Output_Dir

    type(Catalogue)::Catt
    character(120)::Catalogue_Directory, Catalogue_Filename
    integer,dimension(:),allocatable::Catalogue_Cols
    logical::here

    inquire(directory = trim(run_Output_Dir), exist = here)
    if(here == .false.) call system('mkdir '//trim(adjustl(run_Output_Dir)))
    Bayesian_Routines_Output_Directory = trim(run_Output_Dir)

    !## 1: STAGES shear, 2: COMBO17, 3:RRG, 4: Mocks_STAGES; 5:Mocks_COMBO!          
    call common_Catalogue_directories(Cat_Ident, Catalogue_Directory, Catalogue_Filename, Catalogue_Cols)
    call catalogue_readin(Catt, trim(adjustl(Catalogue_Directory))//trim(adjustl(Catalogue_Filename)), 'Tr(J)', Catalogue_Cols)
    
!   call Clip_Sizes(Catt, (/0.e0_double, 20.e0_double/) )
    !call convert_Size_from_Pixel_to_Physical(Catt)                                                                                                          
    call Cut_by_Magnitude(Catt, 23.e0_double) !-Taken from CH08 P1435-!                                                                                      
    if(Analyse_with_Physical_Sizes) call Monte_Carlo_Redshift_Sampling_Catalogue(Catt)
    call Cut_By_PhotoMetricRedshift(Catt, 0.21e0_double) !--Cut out foreground--!                                                                            
    call DM_Profile_Variable_Posteriors_CircularAperture(Catt, Cluster_Pos, Cluster_Aperture_Radius, returned_Cluster_Posteriors)

  end subroutine Mass_Estimate_Single_Run

  !######################################################!
  !-------Mock Catalogue Multiple Run Routines-----------!
  !######################################################!


end program Bayesian_DM_Profile_Constraints
