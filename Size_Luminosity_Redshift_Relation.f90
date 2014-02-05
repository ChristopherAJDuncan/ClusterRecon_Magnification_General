program Size_Luminosity_Redshift_Relation
  !--Returns the averages size in a number of bins in redshift and luminosity. Magnitude entered MUST be rest frame absolute magnitude.--!
  use Catalogues; use Param_Types; use RunTime_Input, only:Verbose; use Convergence_Estimation, only: Average_Size_in_Luminosity_Redshift_Bins, Average_Size_in_Luminosity_Redshift_Errors, calculate_Beta
  implicit none


  character(120)::Catalogue_Filename = 'Catalogues/STAGES_COMBO17_gsMag_size_matched.pzcat'
  integer,dimension(13)::Catalogue_Cols = (/-1,10,11,-1,-1,8,-1,-1,-1,14,16,17,5/) !-ntile, RA, Dec, xpos, ypos, Mag, MagErr, Flux, FluxErr, Size, g1, g2, redshift-!   
  logical::here

  character(120)::Output_Directory = 'Size_Luminosity_Redshift_Relation_Output/'
  character(120)::Output_DirFilename_SLR

  type(Catalogue)::Cat
  integer::nBinZ = 4, nBinM = 5,i
  real(double), dimension(:,:), allocatable::Redshift_Limits, Magnitude_Limits 
  real(double),allocatable::SizeMagRelation(:,:), SizeMagRelation_Error(:,:), Beta(:,:)

  !--Readin In Catalogue--!
  inquire(file = Catalogue_Filename, exist = here)
  if(here == .false.) then
     print *, 'Catalogue:', trim(adjustl(Catalogue_Filename)), ' does not exist, stopping..'
     STOP
  end if
  call catalogue_readin(Cat, Catalogue_Filename, 'FR', Catalogue_Cols)
  call PSF_Correction(Cat, 1)
  call Cut_By_PhotoMetricRedshift(Cat, 0.21e0_double) !--Cut out foreground--!       


  !--Get Bin Limits. This gives **approximately the same number in each bin**--!
  call Calculate_Bin_Limits_by_equalNumber(Cat%Redshift, nBinZ, Redshift_Limits)
  call Calculate_Bin_Limits_by_equalNumber(Cat%Absolute_Magnitude, nBinM, Magnitude_Limits)

  allocate(SizeMagRelation(size(Magnitude_Limits,1), size(Redshift_Limits,1))); SizeMagRelation = 0.e0_double
  allocate(SizeMagRelation_Error(size(Magnitude_Limits,1), size(Redshift_Limits,1))); SizeMagRelation_Error = 0.e0_double
  allocate(Beta(size(Magnitude_Limits,1), size(Redshift_Limits,1))); Beta = 0.e0_double
  !--Do the actual work--!
  call Average_Size_in_Luminosity_Redshift_Bins(Cat, Magnitude_Limits,Redshift_Limits, SizeMagRelation)
  call Average_Size_in_Luminosity_Redshift_Errors(Cat, Magnitude_Limits,Redshift_Limits, SizeMagRelation_Error)
  call calculate_Beta(SizeMagRelation, SizeMagRelation_Error, Magnitude_Limits, Beta)
  ! call fit_functional(

  !--ERRORS?--!

  !--Output--!
  Output_DirFilename_SLR = trim(adjustl(Output_Directory))//'Size_Magnitude_Redshift_Relation.dat'
  call Size_Luminosity_Output(Magnitude_Limits, Redshift_Limits, SizeMagRelation, Output_DirFilename_SLR)
 
  !--Plot--!
  call Luminosity_Redshift_Plotter(Output_DirFilename_SLR, 'SMR')

  Output_DirFilename_SLR = trim(adjustl(Output_Directory))//'Size_Magnitude_Redshift_Relation_Error.dat'
  call Size_Luminosity_Output(Magnitude_Limits, Redshift_Limits, SizeMagRelation_Error, Output_DirFilename_SLR)

  call Luminosity_Redshift_Plotter(Output_DirFilename_SLR, 'SMR_Error')

  Output_DirFilename_SLR = trim(adjustl(Output_Directory))//'Size_Magnitude_Redshift_Relation_beta.dat'
  call Size_Luminosity_Output(Magnitude_Limits, Redshift_Limits, Beta, Output_DirFilename_SLR)

  call Luminosity_Redshift_Plotter(Output_DirFilename_SLR, 'Beta')                

contains

!!$  subroutine Average_Size_in_Luminosity_Redshift(Cat, Lum_Limits, Red_Limits, SLR)
!!$    type(Catalogue), intent(in)::Cat
!!$    real(double),dimension(:,:)::Red_Limits, Lum_Limits
!!$    real(double),intent(out)::SLR(:,:)
!!$
!!$    type(Binned_Catalogue)::BCat_Red
!!$    type(Binned_Catalogue),allocatable::BCat_Lum_Red(:) !-One for each redshift Bin-!
!!$    integer::BinL, BinZ
!!$
!!$    character(10)::by_Size_Type = 'Physical'
!!$
!!$    if(Verbose) print *, 'Calculating Average Size in Redshift and Luminosity Bins....'
!!$
!!$    !-Check Sizes-!
!!$    if(size(SLR,1) /= size(Lum_Limits,1)) stop 'Average_Size_in_Luminosity_Redshift - SLR and Luminosity Bins not conformal...'
!!$    if(size(SLR,2) /= size(Red_Limits,1)) stop 'Average_Size_in_Luminosity_Redshift - SLR and Redshift Bins not conformal...'
!!$
!!$    !--Bin first by redshift using pre-defined Redshift--!
!!$    call bin_catalogue_by_redshift(Cat, Red_Limits, BCat_Red)
!!$
!!$    allocate(BCat_Lum_Red(size(Red_Limits,1)))
!!$    do BinZ = 1, size(BCat_Lum_Red)
!!$       call Bin_Catalogue_by_magnitude(BCat_Red%Cat(BinZ), Lum_Limits, BCat_Lum_Red(BinZ))
!!$    end do
!!$
!!$    if(Verbose) print *, 'Finished Binning'
!!$
!!$    do BinZ = 1, size(Red_Limits,1)
!!$       do BinL = 1, size(Lum_Limits,1)
!!$          SLR(BinL,BinZ) = global_mean_size(BCat_Lum_Red(BinZ)%Cat(BinL), trim(adjustl(by_Size_Type)))
!!$          if(isNAN(SLR(BinL,BinZ))) SLR(BinL,BinZ) = 0.e0_double
!!$       end do
!!$    end do
!!$
!!$    if(Verbose) print *, 'Done.'
!!$
!!$  end subroutine Average_Size_in_Luminosity_Redshift
!!$
!!$  subroutine Average_Size_in_Luminosity_Redshift_Errors(Cat, Lum_Limits, Red_Limits, SLR_ERROR)
!!$    type(Catalogue), intent(in)::Cat
!!$    real(double),dimension(:,:)::Red_Limits, Lum_Limits
!!$    real(double),intent(out)::SLR_ERROR(:,:)
!!$                                                                                            
!!$    !--Random Catalogue Declarations--!
!!$    real(double),dimension(:,:),allocatable::Ran
!!$    Integer,allocatable::seed(:)
!!$    integer::Nseed, Clock
!!$    
!!$    type(Catalogue)::Randomised_Cat
!!$    integer::NBoot
!!$    integer::n,i,j,k
!!$    integer::Boot_Loop
!!$    real(double),allocatable,dimension(:,:,:)::Boot_AvSize
!!$    real(double),allocatable,dimension(:,:)::tBoot_AvSize
!!$    real(double),allocatable::mean_boot_avsize(:,:)
!!$
!!$    !--Used to test for convergence of BootStrap by Varying BootStrap
!!$    logical::Test_Convergence = .true.
!!$    integer::Convergence_Loop, Convergence_Loop_Maximum
!!$    integer::NBoot_Convergence_Lower
!!$    real(double),allocatable:: Previous_Error(:,:)
!!$    real(double)::Convergence_Tolerance = 1.e-6_double
!!$
!!$    if(Verbose) print *, 'Calculating Bootstrap Error Average Size in Redshift and Luminosity Bins....'
!!$
!!$    if(size(SLR_ERROR,1) /= size(Lum_Limits,1)) stop 'Average_Size_in_Luminosity_Redshift - SLR_ERROR and Luminosity Bins not conformal...'
!!$    if(size(SLR_ERROR,2) /= size(Red_Limits,1)) stop 'Average_Size_in_Luminosity_Redshift - SLR_ERROR and Redshift Bins not conformal...'
!!$
!!$    if(Test_Convergence) then
!!$       Convergence_Loop_Maximum = 20
!!$       allocate(Previous_Error(size(SLR_Error,1), size(SLR_Error,2))); Previous_Error = 0.e0_double
!!$       print *, 'Testing for Convergence:'
!!$       Verbose = .false.
!!$    else
!!$       Convergence_Loop_Maximum = 1
!!$    end if
!!$
!!$    NBoot = 1024;
!!$
!!$    do Convergence_Loop = 1, Convergence_Loop_Maximum
!!$       if(Test_Convergence) NBoot = 64*(2**Convergence_loop)
!!$
!!$       allocate(Ran(NBoot, size(Cat%Sizes))); Ran = 0.e0_double
!!$       call RANDOM_SEED(size = NSeed)
!!$       allocate(Seed(NSeed))
!!$       call SYSTEM_CLOCK(COUNT = Clock)
!!$       seed = Clock + (/ (i-1,i=1,NSeed) /)
!!$       call RANDOM_SEED(PUT = seed)
!!$       deallocate(Seed); NSeed = 0; Clock = 0
!!$       
!!$       call RANDOM_NUMBER(Ran)
!!$       !-Convert the random numbers from range 0-1 to 1-Ngal                                                                                                                          
!!$       Ran = nint(Ran*(size(Cat%Sizes)-1) + 1)
!!$       
!!$       allocate(tBoot_AvSize(Size(Lum_Limits,1), size(Red_Limits,1)))
!!$       allocate(Boot_AvSize(NBoot, size(tBoot_AvSize,1), size(tBoot_AvSize,2))); Boot_AvSize = 0.e0_double
!!$       do Boot_loop = 1, NBoot
!!$          !--Set up Randomised Catalogue--!
!!$          Randomised_Cat = Cat; Randomised_Cat%Sizes = 0.e0_double
!!$          do k = 1, size(Cat%Sizes)
!!$             Randomised_Cat%Sizes(k) = Cat%Sizes(Ran(Boot_Loop,k))
!!$          end do
!!$          call convert_Size_from_Pixel_to_Physical(Randomised_Cat)
!!$          
!!$          tBoot_AvSize = 0.e0_double
!!$          call Average_Size_in_Luminosity_Redshift(Randomised_Cat, Lum_Limits, Red_Limits, tBoot_AvSize) 
!!$          Boot_AvSize(Boot_Loop,:,:) = tBoot_AvSize; tBoot_AvSize = 0.e0_double
!!$          
!!$          call Catalogue_Destruct(Randomised_Cat)
!!$       end do
!!$       deallocate(tBoot_AvSize)
!!$       
!!$       allocate(Mean_Boot_AvSize(size(Boot_AvSize,2), size(Boot_AvSize,3))); Mean_Boot_AvSize = 0.e0_double
!!$       do Boot_Loop = 1, NBoot
!!$          Mean_Boot_AvSize = Mean_Boot_AvSize(:,:) + Boot_AvSize(Boot_Loop,:,:)
!!$       end do
!!$       Mean_Boot_AvSize = Mean_Boot_AvSize/nBoot
!!$       
!!$       !--Calculate Variance--!
!!$       do Boot_Loop = 1, NBoot
!!$          SLR_Error = SLR_Error + (Boot_AvSize(Boot_Loop,:,:) - Mean_Boot_AvSize)**2.e0_double
!!$       end do
!!$       SLR_Error = dsqrt(SLR_Error/(1.e0_double*(nBoot-1)))
!!$
!!$       !--Allocate large error for bins with zero galaxies--!
!!$       !--Assumes zero average size is applicable only for bins with no information--!
!!$       do i = 1, size(Boot_AvSize,2)
!!$          do j = 1, size(Boot_AvSize,3)
!!$!             if(all(Boot_AvSize(:,i,j) == 0.e0_double)) SLR_Error(i,j) = 1.e30_double
!!$          end do
!!$       end do
!!$
!!$       deallocate(Ran, Boot_AvSize, Mean_Boot_AvSize)
!!$
!!$       if(Test_Convergence) then
!!$          if(Convergence_Loop == 1) then
!!$             Previous_Error = SLR_Error
!!$             cycle
!!$          end if
!!$
!!$          if(all(SLR_Error-Previous_Error <= Convergence_Tolerance)) then
!!$             !--Success - Converged for all bins--!
!!$             print *, 'Errors have converged to tolerance:', Convergence_Tolerance, ' with nBoot = ', nBoot/2
!!$             print *, 'Press Enter to continue:'
!!$             read(*,*)
!!$             Verbose = .true.
!!$             exit
!!$          end if
!!$
!!$          if(Convergence_Loop == Convergence_Loop_Maximum) then
!!$             print *, 'Failed to find convergence by nBoot = ', nBoot
!!$             STOP
!!$          END if
!!$
!!$          print *, count(SLR_Error-Previous_Error <= Convergence_Tolerance), ' of ', size(SLR_Error), ' bins have converged for nBoot = ', nBoot/2,'. Relooping'
!!$          Previous_Error = SLR_Error
!!$       end if
!!$
!!$    end do
!!$
!!$    if(Verbose) print *, 'Done.'
!!$
!!$  end subroutine Average_Size_in_Luminosity_Redshift_Errors
!!$
!!$
!!$  subroutine calculate_Beta(SLR, SLR_Error, Mag_Limits, Beta_Array)
!!$    !--Assumes result input is binned in terms of ABSOLUTE MAGNITUDE-!
!!$    use nr, only: fit
!!$    real(double), intent(in), dimension(:,:)::SLR, SLR_Error, Mag_Limits
!!$    
!!$    real(double), intent(out), dimension(:,:):: Beta_Array
!!$
!!$
!!$    integer::Method = 1
!!$    real(double)::fit_a, fit_b, fit_siga, fit_sigb, fit_chi2, fit_q
!!$    integer::BinZ, BinL, counter
!!$
!!$    real(double),allocatable::Fit_Mag(:),Fit_Size(:), Fit_Error(:) 
!!$    !--Two Methods 1: Fit using data on eaither side; 2: SubDivide each luminosity bin and fit to subdivisions within bin--!
!!$
!!$    if(Verbose) print *, 'Finding Beta'
!!$
!!$    select case(Method)
!!$    case(1) !-Fit using data on either side-!
!!$       !-Doesn't yet account for grid points with n = 0-!
!!$       do BinZ = 1, size(SLR,2)
!!$          do BinL = 1, size(SLR,1)
!!$             if(BinL ==1) then
!!$                allocate(Fit_Size(2)); Fit_Size = 0.e0_double
!!$                allocate(Fit_Mag(2)); Fit_Mag = 0.e0_double
!!$                allocate(Fit_Error(2)); Fit_Error = 0.e0_double
!!$
!!$                counter = 1
!!$                do i = BinL, BinL+1
!!$                   if(SLR(i,BinZ) == 0.e0_double .or. isNaN(SLR_Error(i,BinZ))) then
!!$                      Fit_Size(counter) = 0.e0_double
!!$                      Fit_Error(counter) = 1.e30_double !-Set arbirarily large as no information on this scale-!
!!$                   else
!!$                      Fit_Size(counter) = dlog10(SLR(i, BinZ))
!!$                      Fit_Error(counter) = 0.434e0_double*SLR_Error(i,BinZ)/SLR(i,BinZ)
!!$                   end if
!!$                   Fit_Mag(counter) = sum(Mag_Limits(i,:))/2.e0_double
!!$                   counter = counter + 1
!!$                end do
!!$
!!$
!!$             elseif(BinL==size(SLR,1)) then
!!$                allocate(Fit_Size(2)); Fit_Size = 0.e0_double
!!$                allocate(Fit_Mag(2)); Fit_Mag = 0.e0_double
!!$                allocate(Fit_Error(2)); Fit_Error = 0.e0_double
!!$
!!$                counter = 1
!!$                do i = BinL-1, BinL
!!$                   if(SLR(i,BinZ) == 0.e0_double .or. isNaN(SLR_Error(i,BinZ))) then
!!$                      Fit_Size(counter) = 0.e0_double
!!$                      Fit_Error(counter) = 1.e30_double
!!$                   else
!!$                      Fit_Size(counter) = dlog10(SLR(i, BinZ))
!!$                      Fit_Error(counter) = 0.434e0_double*SLR_Error(i,BinZ)/SLR(i,BinZ)
!!$                   end if
!!$                   Fit_Mag(counter) = sum(Mag_Limits(i,:))/2.e0_double
!!$                   counter = counter + 1
!!$                end do
!!$             
!!$             else
!!$                allocate(Fit_Size(3)); Fit_Size = 0.e0_double
!!$                allocate(Fit_Mag(3)); Fit_Mag = 0.e0_double
!!$                allocate(Fit_Error(3)); Fit_Error = 0.e0_double
!!$
!!$                counter = 1
!!$                do i = BinL-1, BinL+1
!!$                   if(SLR(i,BinZ) == 0.e0_double .or. isNaN(SLR_Error(i,BinZ))) then
!!$                      Fit_Size(counter) = 0.e0_double
!!$                      Fit_Error(counter) = 1.e30_double
!!$                   else
!!$                      Fit_Size(counter) = dlog10(SLR(i, BinZ))
!!$                      Fit_Error(counter) = 0.434e0_double*SLR_Error(i,BinZ)/SLR(i,BinZ)
!!$                   end if
!!$                   Fit_Mag(counter) = sum(Mag_Limits(i,:))/2.e0_double
!!$                   counter = counter + 1
!!$                end do
!!$
!!$             end if
!!$
!!$             call fit(Fit_Mag, Fit_Size, fit_a, fit_b, fit_siga, fit_sigb, fit_chi2, fit_q, fit_Error)
!!$             Beta_Array(BinL, BinZ) = -2.5e0_double*fit_b
!!$
!!$             deallocate(Fit_Size, Fit_Mag, Fit_Error);
!!$          end do
!!$       end do
!!$
!!$    case default
!!$       print *, 'calculate_Beta -incorrect method chosen, returning'
!!$       return
!!$    end select
!!$
!!$    if(Verbose) print *, 'Found Beta'
!!$
!!$  end subroutine calculate_Beta

  subroutine Size_Luminosity_Output(Lum_Limits, Red_Limits, SLR, Filename)
    real(double),intent(in)::Lum_Limits(:,:), Red_Limits(:,:), SLR(:,:)
    character(*), intent(in)::Filename

    integer::i
    character(5)::fmtstring

    open(59, file = filename)
    write(fmtstring,'(I5)') size(Red_Limits,1)+2
    write(59, '('//trim(adjustl(fmtstring))//'(e14.7,x))') 0.e0_double, Red_Limits(1,1), Red_Limits(:,2)
    do i = 1, size(Lum_Limits,1)
       write(59, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Lum_Limits(i,1), SLR(i,:), 0.e0_double
    end do
    write(59, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Lum_Limits(size(Lum_Limits,1), 2), (/(0.e0_double, I = 1, size(SLR,2)+1)/)
    close(59)
    
    print *, 'File output to:', trim(adjustl(Filename))

  end subroutine Size_Luminosity_Output

  subroutine Luminosity_Redshift_Plotter(Filename, Label)
    character(*),intent(in)::Filename, label
    
    character(120)::Plotter_Name  = 'Size_Luminosity_Redshift_Relation_Plotter.py'
    logical::here

    inquire(file = Filename, exist = here)
    if(here==.false.) print *, 'Plotter - Filename does not exist'
    inquire(file = Plotter_Name, exist = here)
    if(here==.false.) print *, 'Plotter - Plotter does not exist'

    if(here) then
       print *, 'Calling Plotter:'
       call system('python '//trim(adjustl(Plotter_Name))//' '//trim(adjustl(Filename))//' '//trim(adjustl(Label)))
    end if

  end subroutine Luminosity_Redshift_Plotter


end program Size_Luminosity_Redshift_Relation
