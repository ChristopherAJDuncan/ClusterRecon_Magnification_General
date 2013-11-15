program size_luminosity_relation
  use Param_Types; use Catalogues
  implicit none

  character(120):: Catalogue_Filename = 'Catalogues/STAGES_COMBO17_gsMag_size_matched.pzcat'  
  character(120)::Output_Dir = 'Size_Luminosity_Relation_Output/'
  integer,dimension(13)::Catalogue_Cols = (/-1,10,11,-1,-1,8,-1,-1,-1,14,16,17,5/) !-6th Element denotes the magnitude to be used-!!(/-1,10,11,-1,-1,8,-1,-1,-1,14,16,17,5/) !-6th Element denotes the magnitude to be used-!
!-ntile, RA, Dec, xpos, ypos, Mag, MagErr, Flux, FluxErr, Size, g1, g2, redshift-!
  call Size_Luminosity()

  contains

    subroutine size_luminosity()
      type(Binned_Catalogue)::BCat
      type(Catalogue)::Cat

      real(double),allocatable::Size_byMag(:,:), Error_AvSize(:,:) !-RedshiftBin, Mag-!
      integer,allocatable::Occupation(:,:)
      real(double),allocatable::MagBin(:,:), RedshiftBin(:,:)
      real(double)::Mag_low, Mag_high, z_low, z_High
      integer::nMag, nZ

      real(double)::Magnitude_Cut_Lower = 23.e0_double
      logical::Cut_in_PhotZ =.false.
      real(double)::PhotZ_Lower_Cut = 0.21e0_double

      real(double)::fit_lim(2)

      real(double),allocatable::fit_toSize(:,:)

      integer::i

      character(200)::Output_Filename
      character(15)::fmtstring

      !-Read in Catalogue-!
      call catalogue_readin(Cat, Catalogue_Filename, 'FR', Catalogue_Cols)
      call PSF_Correction(Cat, 1)
!     call Cut_by_Magnitude(Cat, Magnitude_Cut_Lower)

      if(Cut_In_PhotZ) call Cut_By_PhotometricRedshift(Cat, PhotZ_Lower_Cut)

      print *, 'Max/Min Sizes:', maxval(Cat%Sizes), minval(Cat%Sizes), ' mean:', sum(Cat%Sizes)/size(Cat%Sizes)

      !--Binning Declarations--!

      Mag_low = minval(Cat%Mag)
      !-Determine Mag_low only from sample which is greater than zero-!
!!$      Mag_low = 100.e0_double
!!$      do i =1, size(Cat%Mag)
!!$         if( (Cat%Mag(i) > 0.e0_double) .and. (Cat%Mag(i) < Mag_Low) ) Mag_low = Cat%Mag(i)
!!$      end do
      Mag_high = maxval(Cat%Mag)
      nMag = 50
      
      z_low = minval(Cat%Redshift)
      z_High = maxval(Cat%Redshift)
      nz = 2
      !-------------------------!
      allocate(MagBin(nMag,2)); MagBin = -1.e0_double
      allocate(RedshiftBin(nZ,2)); RedshiftBin = -1.e0_double

      RedshiftBin(1,:) = (/ 0.e0_double, 0.5e0_double/)
      RedshiftBin(2,:) = (/ 0.5e0_double, 4.5e0_double/)

      print *, 'Magnitude Binning Between limits:', Mag_Low, Mag_High
      

      do i = 1, maxval((/nZ,nMag/))
!!$         if(i <= nz) then
!!$            RedshiftBin(i,1)= z_low + (i-1)*(z_High-z_low)/float(nZ)
!!$            RedshiftBin(i,2) = RedshiftBin(i,1)+(z_High-z_low)/float(nZ)
!!$         end if
         if(i <=nMag) then
            MagBin(i,1) = Mag_low + (i-1)*(Mag_High-Mag_low)/float(nMag)
            MagBin(i,2) = MagBin(i,1)+(Mag_High-Mag_low)/float(nMag)
         end if
      end do
      
      call Bin_catalogue_by_redshift(Cat, RedshiftBin, BCat)

      allocate(Size_byMag(nZ,nMag)); Size_byMag = 0.e0_double
      allocate(Error_AvSize(nZ, nMag)); Error_AvSize = 0.e0_double
      allocate(Occupation(nZ, nMag)); Occupation = 0
      do i = 1, nZ
         call average_size_by_magnitudebin(BCat%Cat(i), MagBin, Size_byMag(i,:), Occupation(i,:))
         call Average_Size_in_Magnitude_Grid_Errors(BCat%Cat(i), MagBin, Error_AvSize(i,:))
      end do

      !-OutputResult-!
      Output_Filename = 'Size_magnitude_relation.dat'
      Output_Filename = trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename))
      write(fmtstring,'(I2)') 2*nZ+1
      fmtstring = '('//trim(adjustl(fmtstring))//'(e14.7,x))'
      print *, 'Format String:', fmtstring
      open(65, file = Output_Filename)
      do i = 1, nMag
         write(65, fmtstring) (MagBin(i,1)+MagBin(i,2))/2.e0_double, Size_byMag(:,i), Error_avSize(:,i) 
      end do
      close(65)
      print *, 'Output to file:', Output_Filename

      Output_Filename = 'Size_magnitude_relation_Occupation.dat'
      Output_Filename = trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename))
      write(fmtstring,'(I2)') nZ+1
      fmtstring = '(e14.7,'//trim(adjustl(fmtstring))//'(I4,x))'
      print *, 'Format String:', fmtstring
      open(65, file = Output_Filename)
      do i = 1, nMag
         write(65, fmtstring) (MagBin(i,1)+MagBin(i,2))/2.e0_double, Occupation(:,i) 
      end do
      close(65)
      print *, 'Output to file:', Output_Filename


      !--Fit--!
      allocate(fit_toSize(size(Size_byMag,1), size(Size_byMag,2))); fit_toSize = 0.e0_double
      fit_lim = (/MagBin(1,1), MagBin(size(MagBin,1),2)/)
      do i = 1, Nz
         print *, 'Fitting Redshift Bin:', i
         call fit_power_law(Size_byMag(i,:), MagBin, Fit_lim, fit_toSize(i,:), Error_avSize(i,:))
      end do

      Output_Filename = 'Size_magnitude_relation_fit.dat'
      Output_Filename = trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename))
      write(fmtstring,'(I2)') nZ+1
      fmtstring = '('//trim(adjustl(fmtstring))//'(e14.7,x))'
      print *, 'Format String:', fmtstring
      open(65, file = Output_Filename)
      do i = 1, nMag
         write(65, fmtstring) (MagBin(i,1)+MagBin(i,2))/2.e0_double, fit_toSize(:,i) 
      end do
      close(65)
      print *, 'Output to file:', Output_Filename


      print *, 'Finished Normally'

    end subroutine size_luminosity

    subroutine average_size_by_magnitudebin(Cat, Bin_Limits, Res, Occupation)
      !--RETURNS THE AVERAGE SIZE IN A MAGNITUDE BIN--!
      real(double),intent(out)::Res(:)
      real(double),intent(in)::Bin_limits(:,:)
      type(Catalogue)::Cat
      integer,optional::Occupation(:)

      integer::i,j,c
      integer,allocatable::iOccupation(:)

      if(size(Res) /= size(Bin_Limits,1)) STOP 'average_size_by_magnitudebin - Result to be output is not conformal with the entered bin limits, stopping ....'
      allocate(iOccupation(size(Res))); iOccupation = 0

      do c = 1, size(Cat%Mag)
         do i = 1, size(Bin_Limits,1)
            if((Cat%Mag(c) > Bin_Limits(i,1)) .and. (Cat%Mag(c) <= Bin_Limits(i,2))) then
               if(Cat%Physical_Sizes(c) < 0.e0_double) cycle
               iOccupation(i) = iOccupation(i) + 1
               Res(i) = Res(i) + Cat%Physical_Sizes(c)
               exit
            end if
         end do
      end do
      where(iOccupation >= 1)
         Res = Res/iOccupation
      end where

      if(any(Res < 0.e0_double)) then
         print *, 'average_size_by_magnitudebin - FATAL ERROR - Some averages sizes are being returned as negative'
         print *, 'Number of failures:', count(Res < 0.e0_double)
      end if

      if(present(Occupation)) then
         if(size(Occupation) /= size(Res)) STOP 'average_size_by_magnitudebin - Occupation is not of the correct size, stopping..' 
         Occupation = -1.e0
         Occupation = iOccupation
         if(any(Occupation < 0.e0_double)) STOP 'average_size_by_magnitudebin - Occupation contains negatives, stopping..'
      end if


    end subroutine average_size_by_magnitudebin
      
    subroutine fit_power_law(Size_Mag, Mag_limits, Fit_limits, fit_functional, Errors)
      !-Fits R = A m^{\beta} by chi^2 minimisation. Returns A and Beta
      !-R = A m^B -> log(R) = log(A) + B log(m) -> fit a straight line to the logarithm
      !--Magnitude used here is ASSUMED to be the apparent magnitude
      use nr; use nrtype
      real(double),intent(in)::Size_Mag(:), Errors(:)
      real(double), intent(in)::Mag_limits(:,:) !-Binning info-!
      real(double),intent(in)::Fit_Limits(:)
      real(double),intent(out)::fit_functional(:)

      real(double),allocatable::Apparent_Luminosity(:)
      real(SP),allocatable::Sub_Size_Mag(:), Sub_Mag(:), fit_Error(:)
      integer::counter, i, Mag_Index_lower, Mag_Index_Higher

      real(double)::fit_a, fit_b, fit_siga, fit_sigb, fit_chi2, fit_q

      if(size(Fit_Limits) /= 2) STOP 'fit_power_law - fit limits not of the correct dimension, stopping...'
      if(size(fit_functional) /= size(Size_Mag)) STOP 'fit_power_law - fit functional is not the right size'
      
      print *, 'Fitting Average Size Luminosity Realtion across limits:', fit_limits

      !-Identify the subset of the Size and Magnitude arrays that falls within the magnitude limits entered-!
      Mag_Index_Lower = 1; Mag_Index_Higher = size(Mag_Limits,1)
      do i = 1, size(Mag_limits,1)
        if( (Mag_Limits(i,1) < Fit_Limits(1)) .and. (Mag_Limits(i,2) >= Fit_limits(1)) ) then
           Mag_Index_lower = i
        end if
        if( (Mag_Limits(i,1) < Fit_Limits(2)) .and. (Mag_Limits(i,2) >= Fit_limits(2)) ) then 
           Mag_Index_higher = i
        end if
     end do
     if(Fit_limits(1) <= minval(Mag_Limits)) Mag_Index_lower = 1

     print *, 'Fit_Limits index:', Mag_Index_lower, Mag_Index_higher
     print *, count(Size_mag == 0.e0_double)

     counter = 0
     do i = Mag_Index_lower, Mag_Index_higher
        if( (Size_Mag(i) > 0.e0_double) ) then
           counter = counter + 1
        end if
     end do

     print *, counter

     !--Set up luminosity--!
!!$     allocate(Apparent_Luminosity(size(Mag_Limits,1))); Apparent_Luminosity = 0.e0_double
!!$     do i = 1, size(Apparent_Luminosity)
!!$        Apparent_Luminosity = 

     allocate(Sub_size_Mag(counter)); Sub_size_Mag = 0.e0_double!Mag_Index_higher-Mag_Index_lower+1)); Sub_Size_mag = 0.e0_double
     allocate(Sub_mag(size(Sub_size_Mag))); Sub_mag = 0.e0_double
     allocate(fit_Error(size(Sub_Size_Mag))); fit_Error = 1.e0_double

     !-Ignores magnitude bins where Average_size == 0.e0_double, i.e. where there is no information-!
     counter = 1
     do i = Mag_Index_lower, Mag_Index_higher
        if( Size_Mag(i) > 0.e0_double) then
           Sub_Mag(counter) = sum(Mag_limits(i,:))/2.e0_double
           Sub_size_mag(counter) = dlog10(Size_Mag(i))
           fit_Error(counter) = (0.434e0_double*Errors(i))/Size_Mag(i)
           counter  = counter + 1
        end if
     end do

     !-Convert to Logarithms-!
     !-Cut out zeros as we don't want to fit to them?-!
     !-Convert to Fluxes, and consequently luminosities-!

     print *, 'Calling Fit'
     !-Fit a line to this-!
     call fit(Sub_Mag, Sub_Size_Mag, fit_a, fit_b, fit_siga, fit_sigb, fit_chi2, fit_q, fit_Error) !,fit_sigma !--Errors on y--!

     print *, 'Power Law fit (log(R) = A +  b*m) with:'
     print *, 'A = ', fit_a, ' ::: B = ', fit_b
     print *, 'Giving fit to R  = AL^{beta} with:'
     print *, 'beta  = ', -2.5e0_double*fit_b

     fit_functional = dexp(fit_a)*( 10.e0_double**(((Mag_Limits(:,1)+Mag_Limits(:,2))/2.e0_double))**fit_b)

     !--Testing output--!
     open(54, file = trim(adjustl(Output_Dir))//trim(adjustl('Fit_functional_logR_m.dat')))
     do i= 1, size(Sub_Mag)
        write(54,'(4(e14.7,x))') Sub_Mag(i), Sub_Size_Mag(i), fit_Error(i), fit_a+fit_b*sub_mag(i)
     end do
     close(54)
     print *, 'Output to:', trim(adjustl(Output_Dir))//trim(adjustl('Fit_functional_logR_m.dat'))

     deallocate(Sub_Size_Mag, Sub_Mag, fit_Error)

    end subroutine fit_power_law

!!$    subroutine Power_Law_ForFitting(M,Size)
!!$      real(double)::M, Size


    !--Produce Error Bars via Bootstrapping Method--!
  subroutine Average_Size_in_Magnitude_Grid_Errors(Cat_In, Mag_Limits, Error_AvSize)
    use Catalogues, only: Catalogue_Destruct
    !-Calculates the errors per pixel by randomly assinging galaxy sizes to galaxy positions-!
    !-Errors output as sigma, NOT sigma^2
    type(Catalogue),intent(in)::Cat_In
    real(double),intent(in)::Mag_Limits(:,:)
    real(double),dimension(:)::Error_AvSize

    type(Catalogue)::Randomised_Cat
    integer::NBoot
    integer::n,i,j,k
    real(double),allocatable,dimension(:,:)::Boot_AvSize
    real(double),allocatable,dimension(:)::tBoot_AvSize
    integer,allocatable:: tBoot_nGrid(:)
    real(double),allocatable,dimension(:)::Mean_Boot_AvSize
    integer,allocatable::Boot_nGrid(:,:)

    !-Random Number Declarations
    real(double),dimension(:,:),allocatable::Ran
    Integer,allocatable::seed(:)
    integer::Nseed, Clock

    print *, 'Bootstrapping errors...'

    NBoot = 100
    !-Calculate array or random numbers, consisting of NBoot sets of random samplings (with replacement) from the entered galaxy size distribution
    allocate(Ran(NBoot, size(Cat_In%Sizes))); Ran = 0.e0_double
    call RANDOM_SEED(size = NSeed)
    allocate(Seed(NSeed))
    call SYSTEM_CLOCK(COUNT = Clock)
    seed = Clock + (/ (i-1,i=1,NSeed) /)
    call RANDOM_SEED(PUT = seed)
    deallocate(Seed); NSeed = 0; Clock = 0
    
    call RANDOM_NUMBER(Ran)
    !-Convert the random numbers from range 0-1 to 1-Ngal
    Ran = nint(Ran*(size(Cat_In%Sizes)-1) + 1)

    !-Construct Average Size grid and Kappa for each bootstrap realisation
    allocate(Boot_AvSize(NBoot, size(Mag_Limits,1))); Boot_AvSize = 0.e0_double
    allocate(Boot_nGrid(size(Boot_AvSize,1), size(Boot_AvSize,2))); Boot_nGrid = 0

    allocate(tBoot_AvSize(size(Mag_Limits,1))); tBoot_AvSize = 0.e0_double
    allocate(tBoot_nGrid(size(Mag_Limits,1))); tBoot_nGrid = 0.e0_double
    do n = 1, NBoot
       !--Set up the Randomised Catalogue--!
       Randomised_Cat = Cat_In; Randomised_Cat%Sizes = 0.e0_double
       do k = 1, size(Cat_In%Sizes)
          Randomised_Cat%Sizes(k) = Cat_In%Sizes(Ran(n,k))
       end do
       !--Calculate Physical Sizes--!
       call convert_Size_from_Pixel_to_Physical(Randomised_Cat)
       tBoot_AvSize = 0.e0_double; tBoot_nGrid = 0.e0_double
       call average_size_by_magnitudebin(Randomised_Cat, Mag_Limits, tBoot_AvSize, tBoot_nGrid)
       Boot_AvSize(n,:) = tBoot_AvSize;
       Boot_nGrid(n,:) = tBoot_nGrid;
       call Catalogue_Destruct(Randomised_Cat)
    end do
    deallocate(tBoot_AvSize,tBoot_nGrid)

    !-Bootstrap sampling are not complete - calculate variance of each.
    allocate(Mean_Boot_AvSize(size(Boot_AvSize,2))); Mean_Boot_AvSize = 0.e0_double

    DO n =1, nBoot
       Mean_Boot_AvSize = Mean_Boot_AvSize + Boot_AvSize(n,:)
    end DO
    Mean_Boot_AvSize = Mean_Boot_AvSize/nBoot

    if(size(Error_AvSize) /= size(MAg_Limits,1)) STOP 'Average_Size_in_Magnitude_Grid_Errors - Error_avSize is not the correct size, stopping...'
    do n = 1, nBoot
       !-Sum_i^nBoot (x^i - <x>^i)**2 over 2 dimensions simultaneously
       Error_AvSize = Error_AvSize + (Boot_AvSize(n,:) - Mean_Boot_AvSize)**2.
    end do

    !-Error output as sigma, NOT sigma^2-!
    Error_AvSize = dsqrt(Error_AvSize/(1.e0_double*(nBoot-1)))

    deallocate(Ran, Boot_nGrid, Boot_AvSize, Mean_Boot_AvSize)

    print *, 'Finished Bootstrap'

  end subroutine Average_Size_in_Magnitude_Grid_Errors
      


end program size_luminosity_relation
