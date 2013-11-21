module Catalogues_Declarations
  !--Modules encompasses declarations needed for referencesing in Interfaces--!
  use Param_Types; use RunTime_Input, only:verbose

end module Catalogues_Declarations

module Catalogues
  use Param_Types; use Catalogues_Declarations
  implicit none

  !--Directory and ReadIn Declarations--!
  character(50),private,parameter::Cat_Dir = 'Catalogues/'
!  character(60),private::Catalogue_Name = 'STAGES_shear.cat'

  !--Catalogue Derived Type--!

  type Catalogue
     character(60)::Name
     character(len = 10)::Mag_Label
     real(double),dimension(:),allocatable::flux, fluxerr, mag, magerr
     integer,dimension(:),allocatable::ntile
     real(double),dimension(:),allocatable::xpos, ypos
     real(double),dimension(:),allocatable::RA,Dec
     character(len = 4)::Sizes_Label
     logical::log_sizes
     real(double),dimension(:),allocatable::Sizes, Physical_Sizes!_FWHM, Size_FR, Size_KSB
     real(double),dimension(:), allocatable::Redshift
     real(double),dimension(:),allocatable:: g1,g2 !-These will contain the shear components.!
     !-Missing are observed ellipticity, anisotropy corrected ellipticity and shear responsivity correction-!
  end type Catalogue

  !-- Binned Catlogue contains a seperate Catalogue derived type for each reshift bin, which contains the subset of galaxies in that Bin--!
  type Binned_Catalogue
     character(5)::Label
     real(double),dimension(:,:),allocatable::Bin_Limits
     integer,dimension(:),allocatable::Occupation
     type(Catalogue),dimension(:),allocatable::Cat
  end type Binned_Catalogue


  !--Global Parameters--!
  real(double)::ACSPixel_Size = 0.049 !-in arcseconds-!


  interface Catalogue_Assign_byGalaxy
     module procedure Catalogue_Assign_byGalaxy_byIndex, Catalogue_Assign_byGalaxy_byCatalogue 
  end interface Catalogue_Assign_byGalaxy

  contains

    !---Catalogue Binning Routines----!
    subroutine Calculate_Bin_Limits_by_equalNumber(Array, nBin, Limits)
      !--Returns the Bin Limits requires for nBin bins on Array, where the Bin is set so that there are an equal number of elements of Array in each bin (therefore suited to an array that contains information for a set of galaxies, and each element is a different galaxy--!
      use nr, only:sort
      real(double),intent(in)::Array(:)
      integer,intent(in)::nBin
      real(double), allocatable,intent(out):: Limits(:,:)

      integer:: nPerBin, i
      real(double),dimension(size(Array))::Sorted_Array

      nPerBin = int(size(Array)/nBin + 1)

      if(allocated(Limits)) deallocate(Limits)
      allocate(Limits(nBin,2)); Limits = 0.e0_double

      Sorted_Array = Array
      !--Sort array--!
      call sort(Sorted_Array)
      
      Limits(1,:) = (/Sorted_Array(1), Sorted_Array(minval( (/nPerBin,size(Sorted_Array)/) ))/)
      do i = 2, nBin
         Limits(i,1) = Limits(i-1,2)
         Limits(i,2) = Sorted_Array(minval((/i*nPerBin,size(Sorted_Array)/)))
      end do
      Limits(nBin,2) = Sorted_Array(size(Sorted_Array)) 

    end subroutine Calculate_Bin_Limits_by_equalNumber

    subroutine bin_catalogue_by_magnitude(Cat,Limits,BCat)
      !--Bins the catalogue by magnitude, returning BCat wof type Binned_Catalogue, which contains a seperate catalogue for each magnitude bin entered in Limits
      type(Catalogue),intent(in)::Cat
      type(Binned_Catalogue)::BCat
      real(double), dimension(:,:),intent(in)::Limits
      
      integer::c,i
      
      type(Catalogue)::tCat !-temporary allocation-!                                                                                                                                                                                                                          
      integer,allocatable::Expected_Occupation(:), Occupation(:)
      
      if(Verbose) print *, 'Binning Catalogue by Magnitude...'
      
      allocate(Expected_Occupation(size(Limits,1))); Expected_Occupation = 0
      allocate(Occupation(size(Limits,1))); Occupation = 0
      
      call Binned_Catalogue_Construct(BCat, 'By Magnitude', size(Limits,1))
      BCat%Bin_Limits = Limits
      
      !-Construct Expected Occupation-!                                                                                                                                                                                    
      do i = 1, size(Limits,1)
         Expected_Occupation(i) = count(Cat%Mag <= Limits(i,2)) - count(Cat%Mag <= Limits(1,1)) - sum(Expected_Occupation(:i))
         Expected_Occupation(size(Expected_Occupation)) = Expected_Occupation(size(Expected_Occupation)) + count(Cat%Mag==Limits(size(Limits,1),2))
         call Catalogue_Construct(BCat%Cat(i), Expected_Occupation(i))
      end do
      
      do c = 1, size(Cat%Sizes)
         do i =1, size(Limits,1)
            if( (Cat%Mag(c) > Limits(i,1)) .and. (Cat%Mag(c) <= Limits(i,2)) ) then
               Occupation(i) = Occupation(i) + 1
               call Catalogue_Assign_byGalaxy_byCatalogue(BCat%Cat(i), Occupation(i), Cat, c)
            end if
            if(i == size(Limits,1) .and. Cat%Mag(c) == Limits(i,2)) then
               Occupation(i) = Occupation(i) + 1
               call Catalogue_Assign_byGalaxy_byCatalogue(BCat%Cat(i), Occupation(i), Cat, c)
            end if
         end do
      end do

      do i = 1, size(Limits,1)
         if(Occupation(i) /= Expected_Occupation(i)) print *, 'Expected Error (FATAL): bin_catalogue_by_redshift: Occupation of a bin is not equal to the expected occupation', Occupation(i), Expected_Occupation(i)
      end do
      
      if(Verbose) then
         print *, 'Sample of ', size(Cat%Sizes), ' galaxies split as:'
         do i =1, size(Limits,1)
            print *, 'Bin:', i, ' has ', size(BCat%Cat(i)%Sizes), ' galaxies. Bin is:', Limits(i,:)
         end do
      end if

      BCat%Occupation = Occupation
      
      if(Verbose) print *, 'Done.'
      
    end subroutine bin_catalogue_by_magnitude

    subroutine bin_catalogue_by_redshift(Cat,Limits,BCat)
      !--Bins the catalogue by redshift, returning BCat wof type Binned_Catalogue, which contains a seperate catalogue for each redshift bin entered in Limits
      type(Catalogue),intent(in)::Cat
      type(Binned_Catalogue)::BCat
      real(double), dimension(:,:),intent(in)::Limits
      
      integer::c,i
      
      type(Catalogue)::tCat !-temporary allocation-!                                                                                                                                                                                                                          
      integer,allocatable::Expected_Occupation(:), Occupation(:)
      
      if(Verbose) print *, 'Binning Catalogue by Redshift...'
      
      allocate(Expected_Occupation(size(Limits,1))); Expected_Occupation = 0
      allocate(Occupation(size(Limits,1))); Occupation = 0
      
      call Binned_Catalogue_Construct(BCat, 'By Redshift', size(Limits,1))
      BCat%Bin_Limits = Limits
      
      !-Construct Expected Occupation-!                                                                                                                                                                                                                                       
      do i = 1, size(Limits,1)
         Expected_Occupation(i) = count(Cat%Redshift <= Limits(i,2)) - count(Cat%Redshift <= Limits(1,1))  - sum(Expected_Occupation(:i))
         !--Edit--!
         Expected_Occupation(size(Expected_Occupation)) = Expected_Occupation(size(Expected_Occupation)) + count(Cat%Redshift==Limits(size(Limits,1),2))
         !--------!
         call Catalogue_Construct(BCat%Cat(i), Expected_Occupation(i))
      end do
      
      do c = 1, size(Cat%Sizes)
         do i =1, size(Limits,1)
            if( (Cat%Redshift(c) > Limits(i,1)) .and. (Cat%Redshift(c) <= Limits(i,2)) ) then
               Occupation(i) = Occupation(i) + 1
               call Catalogue_Assign_byGalaxy_byCatalogue(BCat%Cat(i), Occupation(i), Cat, c)
            end if
            if(i == size(Limits,1) .and. Cat%Redshift(c) == Limits(i,2)) then
               Occupation(i) = Occupation(i) + 1
               call Catalogue_Assign_byGalaxy_byCatalogue(BCat%Cat(i), Occupation(i), Cat, c)
            end if
         end do
      end do
      !--Assign all the galaxies at the upper boundary to the highest redshift bin--!

      do i = 1, size(Limits,1)
         if(Occupation(i) /= Expected_Occupation(i)) print *, 'Expected Error (FATAL): bin_catalogue_by_redshift: Occupation of a bin is not equal to the expected occupation', Occupation(i), Expected_Occupation(i)
      end do
      
      if(Verbose) then
         print *, 'Sample of ', size(Cat%Sizes), ' galaxies split as:'
         do i =1, size(Limits,1)
            print *, 'Bin:', i, ' has ', size(BCat%Cat(i)%Sizes), ' galaxies. Bin is:', Limits(i,:)
         end do
      end if

      BCat%Occupation = Occupation
      
      if(Verbose) print *, 'Done.'
      
    end subroutine bin_catalogue_by_redshift

    subroutine bin_catalogue_by_Redshift_and_Magnitude(Cat, RedLimits, MagLimits, BBCat, Bin_Occupation)
      type(Binned_Catalogue),allocatable,intent(out)::BBCat(:) !-Size nZ, i.e. one Binned_Catalogue for each redshift bin-!

      real(double),intent(in)::MagLimits(:,:), RedLimits(:,:)
      type(Catalogue),intent(in)::Cat
      integer,allocatable,intent(out),optional::Bin_Occupation(:,:)

      type(Binned_Catalogue)::BCat
      integer::BinZ, BinM

      !--Bin first by redshift using pre-defined Redshift--!                                                                                                                                                                                 
      call bin_catalogue_by_redshift(Cat, RedLimits, BCat)

      allocate(BBCat(size(RedLimits,1)))
      do BinZ = 1, size(BBCat)
         call Bin_Catalogue_by_magnitude(BCat%Cat(BinZ), MagLimits, BBCat(BinZ))
         if(present(Bin_Occupation)) then
            allocate(Bin_Occupation(size(RedLimits,1), size(MagLimits,1))); Bin_Occupation = 0.e0_double
            do BinM = 1, size(MagLimits,1)
               Bin_Occupation(BInZ, BinM) = BBCat(BinZ)%Occupation(BinM)
            end do
         end if
      end do


    end subroutine bin_catalogue_by_Redshift_and_Magnitude

    !------------------END BINNING ROUTINES-------------------------!

    subroutine convert_Size_from_Pixel_to_Physical(Cat)
      use Cosmology
      !-Converts from a pixel size to physical size. Can only be done if there is redshift information for each galaxy!!!-!
      !--PHysical Size in Mpc/h, and is in proper size coo-ords (ie extra factor of a)
      type(Catalogue)::Cat

      integer::g

      real(double)::Pixel_Size_Radians

      if(Verbose) print *, 'Calculating Physical Size using Pixel Size = ', ACSPixel_Size, ' arcseconds'

      if(all(Cat%Redshift < 0.e0_double)) STOP 'convert_Size_from_Pixel_to_Physical - All Galaxies in catalogue have not been assigned a redshift, exiting...'

    if(allocated(Cat%Physical_Sizes)) deallocate(Cat%Physical_Sizes)
    allocate(Cat%Physical_Sizes(size(Cat%Sizes))); Cat%Physical_Sizes = -1.e0_double

    Pixel_Size_Radians = (ACSPixel_Size*3.141592654e0_double)/(180.e0_double*3600.e0_double)
!!$    do g = 1, size(Cat%Sizes)
!!$       if(Cat%Redshift(g) > 0.e0_double) Cat%Physical_Sizes(g) = Pixel_Size_Radians*Cat%Sizes(g)*angular_diameter_distance_fromRedshift(Cat%Redshift(g))
!!$    end do

    where(Cat%Redshift > 0.e0_double)
       Cat%Physical_Sizes = (Pixel_Size_Radians*Cat%Sizes*angular_diameter_distance_fromRedshift(Cat%Redshift))/(1.e0_double+Cat%Redshift)
    end where

    if(Verbose) print *, 'Done'

    end subroutine convert_Size_from_Pixel_to_Physical

    subroutine match_Redshift_byCatalogue_byPosition(Cat, Ref_Cat)
      !-For each galaxy in Cat, finds the *exact* match (within small tolerance for rounding errors) in RA and Dec in teh Ref_Cat, and sets the Redshift of that galaxy to the redshift stored in ref Cat
      type(Catalogue)::Cat, Ref_Cat

      integer::nFail, nPass
      real(double)::tol = 1.e-6_double

      integer::i,j

      nPass = 0.e0_double; nFail = 0.e0_double
      do i = 1, size(Cat%Redshift)
         do j = 1, size(Ref_Cat%Redshift)
            if( (dabs(Cat%RA(i) - Ref_Cat%RA(j)) <= tol) .and. (dabs(Cat%Dec(i) - Ref_Cat%Dec(j)) <= tol) ) then
               Cat%Redshift(i) = Ref_Cat%Redshift(j)
               nPass = nPass + 1
               exit
            end if
            if(j==size(Ref_Cat%Redshift)) nFail = nFail + 1
         end do
      end do
      if(nPass + nFail /= size(Cat%Redshift)) print *, 'match_Redshift_byCatalogue_byPosition - Error in assigning redshifts: nPass and nFail do not add up to original length'
      
      if(Verbose) print *, 'Matching Redshifts:', nPass, ' galaxies out of', size(Cat%Redshift),' were matched'

    end subroutine match_Redshift_byCatalogue_byPosition

    subroutine Cut_By_PhotometricRedshift(Cat, LowerCut, UpperCut, Discard_Gals_wout_Redshift)
      !-Galaxies with Redshifts can be isolated by setting LowerCut <= 0.e0_double and Discard_Gals_wout_Redshift = .true.-!
      !-if Discard_Gals_wout_Redshift is present and true, then only galaxies with redhift information (defined as having redshift > 0) are returned-!
      type(Catalogue)::Cat
      real(double), intent(in),optional::LowerCut,UpperCut
      logical, intent(in),optional::Discard_Gals_wout_Redshift

      real(double)::ilower, iupper
      type(catalogue)::Temp_Cat
      logical::Discard_Gals
      integer::nPass,i

!!$      INTERFACE
!!$         subroutine Cut_By_PhotometricRedshift(Cat, LowerCut, UpperCut, Discard_Gals_wout_Redshift)
!!$           use Param_Types; use Catalogues_Declarations
!!$           type(Catalogue)::Cat
!!$           
!!$           real(double), intent(in),optional::LowerCut,UpperCut
!!$           logical, intent(in),optional::Discard_Gals_wout_Redshift
!!$         end subroutine Cut_By_PhotometricRedshift
!!$      END INTERFACE
      
      if((present(lowerCut)==.false.) .and. (present(upperCut)==.false.)) then
         print *, 'Error - Cut_By_PhotometricRedshift - Either a lower or upper cut needs to be entered, returning without cutting'
         return
      end if

      Discard_Gals = .false.
      if(present(Discard_Gals_wout_Redshift)) Discard_Gals = Discard_Gals_wout_Redshift

      if(present(lowerCut)) then
         ilower =lowerCut
      else
         ilower = minval(Cat%Redshift)
      end if
      ilower = maxval((/0.e0_double, ilower/)) !-Ensures lower never goes below zero-!
      if(present(upperCut)) then
         iupper =upperCut
      else
         iupper = maxval(Cat%Redshift)
      end if
  
      if(Verbose) print *, 'Number of assigned galaxy redshifts:', count(Cat%Redshift >= 0.e0_double)
      if(Verbose) print *, 'Cutting Catalogue to Phot-Z between limits:', ilower, iupper

      if(ilower >= iupper) then
         STOP 'FATAL ERROR - Cut_By_PhotometricRedshift - lower cut larger than upper cut'
      end if
      
      Temp_Cat = Cat

      call Catalogue_Destruct(Cat)

      if(Discard_Gals) then
         nPass = count(Temp_Cat%Redshift <= iupper) - count(Temp_Cat%Redshift < iLower) 
      else
         nPass = count(Temp_Cat%Redshift <= iupper) - count(Temp_Cat%Redshift < iLower) + count(Temp_Cat%Redshift < 0.e0_double)
      end if

      call Catalogue_Construct(Cat, nPass)

      nPass = 0
      if(Discard_Gals) then
         do i = 1, size(Temp_Cat%Redshift)
            if((Temp_Cat%Redshift(i) >= ilower) .and.(Temp_Cat%Redshift(i) <= iupper)) then
               nPass = nPass + 1
               call Catalogue_Assign_byGalaxy_byCatalogue(Cat, nPass, Temp_Cat, i)
            end if
         end do
      else
         do i = 1, size(Temp_Cat%Redshift)
            if( (Temp_Cat%Redshift(i)< 0.e0_double) .or. ((Temp_Cat%Redshift(i) >= ilower) .and.(Temp_Cat%Redshift(i) <= iupper)) ) then
               nPass = nPass + 1
               call Catalogue_Assign_byGalaxy_byCatalogue(Cat, nPass, Temp_Cat, i)
            end if
         end do
      end if
      if(nPass /= size(Cat%Mag)) then
         print *, 'Error - Cut_By_PixelSize - Error in Assigning galxies that pass cuts', nPass, size(Cat%Mag)
         STOP
      end if

      if(Verbose) print *, 'Out of:', size(Temp_Cat%Redshift),' orginal galaxies, ', size(Cat%Redshift), ' galaxies passed the redshift cut' 

      call Catalogue_Destruct(Temp_Cat)
            
    end subroutine Cut_By_PhotometricRedshift


    subroutine Cut_By_Magnitude(Cat, lowerCut, upperCut)
      type(Catalogue)::Cat
      real(double),intent(in),optional::lowerCut, upperCut

      real(double)::ilower, iupper
      type(Catalogue)::Temp_Cat
      integer::i
      integer::nPass

      if((present(lowerCut)==.false.) .and. (present(upperCut)==.false.)) then
         print *, 'Error - Cut_By_PixelSize - Either a lower or upper cut needs to be entered, returning without cutting'
         return
      end if

      if(present(lowerCut)) then
         ilower =lowerCut
      else
         ilower = minval(Cat%Mag)
      end if
      if(present(upperCut)) then
         iupper =upperCut
      else
         iupper = maxval(Cat%Mag)
      end if

      if(ilower >= iupper) then
         print*,  'FATAL ERROR - Cut_By_Magnitude - lower cut larger than upper cut'
         print *, ilower, iupper
         STOP
      end if

      print *, 'Cutting Catalogue by magnitude between limits:', ilower, iupper
      
      Temp_Cat = Cat

      call Catalogue_Destruct(Cat)

!!$      nPass = 0
!!$      do i = 1, size(Temp_Cat%Mag)
!!$         if((Temp_Cat%Mag(i) >= ilower) .and. (Temp_Cat%Mag(i) <= iupper)) nPass = nPass + 1
!!$      end do

      nPass = count(Temp_cat%Mag <= iupper) + count(Temp_Cat%Mag >= ilower) - size(Temp_Cat%Mag)


      call Catalogue_Construct(Cat, nPass)

      nPass = 0
      do i = 1, size(Temp_Cat%Mag)
         if((Temp_Cat%Mag(i) >= ilower) .and.(Temp_Cat%Mag(i) <= iupper)) then
            nPass = nPass + 1
            call Catalogue_Assign_byGalaxy_byCatalogue(Cat, nPass, Temp_Cat, i)
         end if
      end do
      if(nPass /= size(Cat%Mag)) then
         print *, 'Error - Cut_By_PixelSize - Error in Assigning galxies that pass cuts', nPass, size(Cat%Mag)
         STOP
      end if

      call Catalogue_Destruct(Temp_Cat)

    end subroutine Cut_By_Magnitude


    subroutine Cut_By_PixelSize(Cat, lowerCut, upperCut)
      type(Catalogue)::Cat
      real(double),intent(in),optional::lowerCut, upperCut

      real(double)::ilower, iupper
      type(Catalogue)::Temp_Cat
      integer::i
      integer::nPass

      if((present(lowerCut)==.false.) .and. (present(upperCut)==.false.)) then
         print *, 'Error - Cut_By_PixelSize - Either a lower or upper cut needs to be entered, returning without cutting'
         return
      end if

      if(present(lowerCut)) then
         ilower =lowerCut
      else
         ilower = minval(Cat%Sizes)
      end if
      if(present(upperCut)) then
         iupper =upperCut
      else
         iupper = maxval(Cat%Sizes)
      end if

      if(ilower >= iupper) then
         STOP 'FATAL ERROR - Cut_By_PixelSize - lower cut larger than upper cut'
      end if
      
      Temp_Cat = Cat

      call Catalogue_Destruct(Cat)

!!$      nPass = 0
!!$      do i = 1, size(Temp_Cat%Sizes)
!!$         if((Temp_Cat%Sizes(i) >= ilower) .and. (Temp_Cat%Sizes(i) <= iupper)) nPass = nPass + 1
!!$      end do

      nPass = count(Temp_cat%Sizes <= iupper) + count(Temp_Cat%Sizes >= ilower) - size(Temp_Cat%Sizes)


      call Catalogue_Construct(Cat, nPass)

      nPass = 0
      do i = 1, size(Temp_Cat%Sizes)
         if((Temp_Cat%Sizes(i) >= ilower) .and.(Temp_Cat%Sizes(i) <= iupper)) then
            nPass = nPass + 1
            call Catalogue_Assign_byGalaxy_byCatalogue(Cat, nPass, Temp_Cat, i)
         end if
      end do
      if(nPass /= size(Cat%Sizes)) then
         print *, 'Error - Cut_By_PixelSize - Error in Assigning galxies that pass cuts', nPass, size(Cat%Sizes)
         STOP
      end if

      call Catalogue_Destruct(Temp_Cat)

    end subroutine Cut_By_PixelSize
    
    subroutine PSF_Correction(Cat, Correction_Type)
      !--Applies size correction by correction type--!
      type(Catalogue)::Cat
      integer::Correction_Type !-0:None, 1:Pixel Size Cut, 2:Subtract PSF size in quadrature-!

      real(double)::Pixel_Size_Cut(2)
      type(Catalogue)::Temporary_Cat
      integer::ii, i

      if(Correction_Type == 0) return
      if(any((/0,1,2/) == Correction_Type)==.false.) print *, 'WARNING - PSF_Correction - Correction Type not supported' 
      
      !--Main body, applies correction--!
      select case(Correction_Type)
      case(1) !-Cut by Pixel_Size_Cut-!
         !-Set pixel size cut, lower to reduce PSF Corrections, Upper to remove outliers-!
         Pixel_Size_Cut = (/4.e0_double, 20.e0_double/)
         call Clip_Sizes(Cat, Pixel_Size_Cut)
      end select

    end subroutine PSF_Correction

    subroutine Clip_Sizes(Cat, SizeClip)
      type(catalogue),intent(in)::Cat
      real(double),intent(in),dimension(2)::SizeClip !-Lower, Upper-!
      
      print *, 'Clipping Size by:', SizeClip

      call Cut_By_PixelSize(Cat, SizeClip(1), SizeClip(2))
    end subroutine Clip_Sizes

    function global_mean_size(Cat, by_Size_Type)
      !-Returns the mean size of galaxies in the catalogue, as 1/N*sum(size)-!
      type(Catalogue)::Cat
      real(double)::global_mean_size
      character(*), intent(in)::by_Size_Type

      if(by_Size_Type == 'Physical' .or. by_Size_Type == 'physical' .or. by_Size_Type == 'Phys'.or. by_Size_Type == 'phys') then
         global_mean_size = sum(Cat%Physical_Sizes)/size(Cat%Physical_Sizes)
      elseif(by_Size_Type == 'Pixel' .or. by_Size_Type == 'pixel' .or. by_Size_Type == 'Pix'.or. by_Size_Type == 'pix') then
         global_mean_size = sum(Cat%Sizes)/size(Cat%Sizes)
      else
         STOP 'global_mean_size - invalid size type entered'
      end if

    end function global_mean_size

    function size_variance(Cat, by_Size_Type, global_mean)
      !-Returns the mean size of galaxies in the catalogue, as 1/N*sum(size)-!                                                                                                                                                               
      use Statistics, only:get_variance
      type(Catalogue)::Cat
      real(double)::size_variance
      character(*), intent(in)::by_Size_Type
      real(double), intent(in),optional::global_mean

      real(double):: mean_size_global

      if(present(global_mean)) then
         mean_size_global = global_mean
      else
!         print *, 'Size Varaince: Calculating mean from Catalogue entered'
         mean_size_global =  global_mean_size(Cat, by_Size_Type)
      end if

      if(by_Size_Type == 'Physical' .or. by_Size_Type == 'physical' .or. by_Size_Type == 'Phys'.or. by_Size_Type == 'phys') then
         size_variance =  get_variance(Cat%Physical_Sizes, Cat%Physical_Sizes, mean_size_global, mean_size_global)
      elseif(by_Size_Type == 'Pixel' .or. by_Size_Type == 'pixel' .or. by_Size_Type == 'Pix'.or. by_Size_Type == 'pix') then
         size_variance =  get_variance(Cat%Sizes, Cat%Sizes, mean_size_global, mean_size_global)
      else
         STOP 'global_mean_size - invalid size type entered'
      end if

    end function size_variance

    !--Catalogue IO----!

    subroutine Catalogue_ReadIn(Cat, Cat_Filename, Size_label, Cols)
      use IO, only: ReadIn, IO_Verbose, skip_header; use Strings_lib, only:remove_FileExtension_fromString
      !-Direct ReadIn of Catalogue File-!
      type(Catalogue)::Cat
      character(120),intent(in)::Cat_Filename

      real(double),allocatable::Cat_2D(:,:)

      integer,intent(in)::Cols(:) !-ntile, RA, Dec, xpos, ypos, Mag, MagErr, Flux, FluxErr, Size, g1, g2, redshift-!
      integer::Cols_Length = 13
      integer, allocatable::iCols(:)
      integer::Column_Index
      character(200)::header_Filename
      logical::here

      character(*)::Size_Label

      !--Cols describes where the code is to find the relelvant information--!
      if( all(Cols == 0) ) then
         !--Attempt read in from header--!
         header_Filename = remove_FileExtension_fromString(Cat_Filename)
         if(len(header_filename) == 1) STOP 'Catalogue_ReadIn - Error with constructing header filename, length of 1'
         header_Filename = trim(adjustl(header_Filename(1:len(trim(adjustl(header_filename)))-1)))//'_header.dat'
         inquire(file = header_filename, exist = here)
         if(here) then
            !-Read in cols-!
            open(37, file = header_filename)
            call skip_header(37, '#')
            allocate(iCols(Cols_Length)); iCols = -100
            read(37, *) iCols
            if(any(iCols == -100)) STOP 'Catalogue_ReadIn - Error  in iCols header readin, not all columns are accounted for'
         else
            print *, 'Header_Filename:', Header_Filename
            STOP 'Catalogue_ReadIn - Header_Filename does not exist'
         end if
      else
         if(size(Cols) /= Cols_Length) then 
            print *, 'Catalogue Read In - Cols Input is not of the correct size, please input (any < 0 will be ignored):'
            print *, 'ntile, RA, Dec, xpos, ypos, Mag, MagErr, Flux, FluxErr, Size, g1, g2'
            allocate(iCols(Cols_Length)); iCols = -1
            read(*,*) iCols
         else
            allocate(iCols(size(Cols))); iCols = Cols
         end if
      end if

      !--Read into 2D format--!
      call ReadIn(Cat_2D, filename  = trim(adjustl(Cat_Filename)), tabbed = .false., header_label = '#')
      !--This outputs as Cat2D(cols, rows), so 2nd dimension lists all galaxies in catalogue--!
      print *, 'Catalogue read in from:', trim(adjustl(Cat_Filename))

      print *, 'Constructing Catalogue with:', size(Cat_2d,2), ' entries'

      !-Split 2D format into permanent Catalogue Storage--!
      call Catalogue_Construct(Cat, size(Cat_2D,2))
      Cat%Name = trim(adjustl(Cat_Filename))

      do Column_Index = 1, size(iCols)
         if(iCols(Column_Index) < 0) cycle
         select case (Column_Index)
         case(1) !-ntile-!
            Cat%ntile = Cat_2D(iCols(Column_Index),:)
         case(2) !-RA-!
            Cat%RA = Cat_2D(iCols(Column_Index),:)
         case(3) !-Dec-!
            Cat%Dec = Cat_2D(iCols(Column_Index),:)
         case(4) !-xpos-!
            Cat%xpos = Cat_2D(iCols(Column_Index),:)
         case(5) !-ypos-!
            Cat%ypos = Cat_2D(iCols(Column_Index),:)
         case(6) !-Mag-!
            Cat%Mag = Cat_2D(iCols(Column_Index),:)
         case(7) !-MagErr-!
            Cat%MagErr = Cat_2D(iCols(Column_Index),:)
         case(8) !-Flux-!
            Cat%Flux = Cat_2D(iCols(Column_Index),:)
         case(9) !-FluxErr-!
            Cat%FluxErr = Cat_2D(iCols(Column_Index),:)
         case(10) !-Size-!
            Cat%Sizes_Label= Size_Label
            Cat%Sizes = Cat_2D(iCols(Column_Index),:)
         case(11) !-g1-!
            Cat%g1 = Cat_2D(iCols(Column_Index),:)
         case(12) !-g2-!
            Cat%g2 = Cat_2D(iCols(Column_Index),:)
         case(13)
            Cat%Redshift = Cat_2D(iCols(Column_Index),:)
         case default
            STOP 'iCols constructed to a size which I cannot deal with yet'
         end select
      end do

      !--Calculate Physical_sizes for the Catalogue--!
      call convert_Size_from_Pixel_to_Physical(Cat)

      Cat_2D = 0.e0_double; deallocate(Cat_2D)

    end subroutine Catalogue_ReadIn

    subroutine Catalogue_Output(Cat, Output_Filename, Pad)
      type(Catalogue),intent(in)::Cat
      character(len = *),intent(in)::Output_filename
      logical,intent(in),optional::Pad

      logical::iPad
      integer::i

      if(present(Pad)) then
         iPad = Pad
      else
         iPad  = .false.
      end if

      open(unit = 30, file =  trim(adjustl(Output_Filename)))
      do i = 1, size(Cat%Sizes)
         if(iPad) then
            if(trim(adjustl(Cat%Sizes_Label))=='FWHM') then
               write(30, '(I2, x, 20(e14.7, x))') Cat%ntile(i), Cat%RA(i), Cat%Dec(i), Cat%xpos(i), Cat%ypos(i), Cat%Mag(i), Cat%MagErr(i), Cat%Flux(i), Cat%FluxErr(i), Cat%Sizes(i), 0., 0., 0., 0., 0., 0., 0., 0., Cat%g1(i), Cat%g2(i)
            elseif(trim(adjustl(Cat%Sizes_Label))=='FR') then
               write(30, '(I2, x, 20(e14.7,x))') Cat%ntile(i), Cat%RA(i), Cat%Dec(i), Cat%xpos(i), Cat%ypos(i), Cat%Mag(i), Cat%MagErr(i), Cat%Flux(i), Cat%FluxErr(i),0., Cat%Sizes(i),  0., 0., 0., 0., 0., 0., 0., Cat%g1(i), Cat%g2(i)
            elseif(trim(adjustl(Cat%Sizes_Label))=='KSB') then
               write(30, '(I2, x, 20(e14.7,x))') Cat%ntile(i), Cat%RA(i), Cat%Dec(i), Cat%xpos(i), Cat%ypos(i), Cat%Mag(i), Cat%MagErr(i), Cat%Flux(i), Cat%FluxErr(i), 0., 0., Cat%Sizes(i), 0., 0., 0., 0., 0., 0., Cat%g1(i), Cat%g2(i)
            end if
         else
            write(30, '(I2, x, 11(e14.7,x))') Cat%ntile(i), Cat%RA(i), Cat%Dec(i), Cat%xpos(i), Cat%ypos(i), Cat%Mag(i), Cat%MagErr(i), Cat%Flux(i), Cat%FluxErr(i), Cat%Sizes(i), Cat%g1(i), Cat%g2(i)
         end if
      end do
      close(30)

    end subroutine Catalogue_Output

    !----------------CATALOGUE DERIVED TYPE SUB/FUNC----------------------------!
    subroutine Binned_Catalogue_Construct(BCat, Label, nBin)
      type(Binned_Catalogue),intent(out)::BCat
      character(*),intent(in)::Label
      integer, intent(in)::nBin!#Bins

      integer::i

      BCat%Label = Label
      allocate(BCat%Bin_limits(nBin,2)); BCat%Bin_Limits = 0.e0_double
      allocate(BCat%Occupation(nBin)); BCat%Occupation = 0
      allocate(BCat%Cat(nBin))

    end subroutine Binned_Catalogue_Construct

    subroutine Binned_Catalogue_Destruct(BCat)
      type(Binned_Catalogue)::BCat

      integer::i

      BCat%Label = ' '
      if(allocated(BCat%Bin_Limits)) deallocate(BCat%Bin_Limits)
      if(allocated(BCat%Cat)) then
         do i = 1, size(BCat%Cat)
            call Catalogue_Destruct(BCat%Cat(i))
         end do
         deallocate(BCat%Cat)
      end if

    end subroutine Binned_Catalogue_Destruct


    subroutine Catalogue_Destruct(Cat)
      type(Catalogue)::Cat

      Cat%Name = ' '

      if(allocated(Cat%ntile)) deallocate(Cat%ntile)
      if(allocated(Cat%flux)) deallocate(Cat%flux)
      if(allocated(Cat%fluxerr)) deallocate(Cat%fluxerr)
      if(allocated(Cat%mag)) deallocate(Cat%mag)
      if(allocated(Cat%magerr)) deallocate(Cat%magerr)
      if(allocated(Cat%xpos)) deallocate(Cat%xpos)
      if(allocated(Cat%ypos)) deallocate(Cat%ypos)
      if(allocated(Cat%RA)) deallocate(Cat%RA)
      if(allocated(Cat%Dec)) deallocate(Cat%Dec)
      if(allocated(Cat%Sizes)) deallocate(Cat%Sizes)
      if(allocated(Cat%Physical_Sizes)) deallocate(Cat%Physical_Sizes)
      if(allocated(Cat%Redshift)) deallocate(Cat%Redshift)
      if(allocated(Cat%g1)) deallocate(Cat%g1)
      if(allocated(Cat%g2)) deallocate(Cat%g2)
!      if(allocated(Cat%Size_FR)) deallocate(Cat%Size_KSB)
!      if(allocated(Cat%Size_KSB)) deallocate(Cat%Size_KSB)

    end subroutine Catalogue_Destruct

    subroutine Catalogue_Construct(Cat, nObj)
      type(Catalogue)::Cat
      integer,intent(in)::nObj
       
      !-If Catalogue contains information then destroy-!
      if(Catalogue_Constructed(Cat)) call Catalogue_Destruct(Cat)

      Cat%Name = ' '

      allocate(Cat%flux(nObj)); Cat%flux = 0.e0_double
      Cat%Mag_Label = ' '
      allocate(Cat%mag(nObj)); Cat%mag = 0.e0_double
      allocate(Cat%xpos(nObj)); Cat%xpos = 0.e0_double
      allocate(Cat%ypos(nObj)); Cat%ypos = 0.e0_double 
      Cat%Sizes_Label = ' '
      allocate(Cat%Sizes(nObj)); Cat%Sizes = 0.e0_double
      allocate(Cat%Physical_Sizes(nObj)); Cat%Physical_Sizes = 0.e0_double
 
      allocate(Cat%Fluxerr(nObj)); Cat%Fluxerr = 0.e0_double
      allocate(Cat%Magerr(nObj)); Cat%Magerr = 0.e0_double
      allocate(Cat%ntile(nObj)); Cat%ntile = 0
      allocate(Cat%RA(nObj)); Cat%RA = 0.e0_double
      allocate(Cat%Dec(nObj)); Cat%Dec = 0.e0_double
      allocate(Cat%g1(nObj)); Cat%g1 = 0.e0_double
      allocate(Cat%g2(nObj)); Cat%g2 = 0.e0_double
      allocate(Cat%Redshift(nObj)); Cat%Redshift = -100.e0_double
!      allocate(Cat%Size_FR(nObj)); Cat%Size_FR = 0.e0_double
!      allocate(Cat%Size_KSB(nObj));  Cat%Size_KSB = 0.e0_double

    end subroutine Catalogue_Construct

    logical function Catalogue_Constructed(Cat)
      type(Catalogue),intent(in)::Cat

      if(Cat%Name /= ' ') Catalogue_Constructed = .true.

      Catalogue_Constructed = .false.
      if(allocated(Cat%flux)) Catalogue_Constructed = .true.
      if(allocated(Cat%mag)) Catalogue_Constructed = .true.
      if(allocated(Cat%xpos)) Catalogue_Constructed = .true.
      if(allocated(Cat%ypos)) Catalogue_Constructed = .true.
      if(allocated(Cat%Sizes)) Catalogue_Constructed = .true.
      if(allocated(Cat%Physical_Sizes)) Catalogue_Constructed = .true.
!      if(allocated(Cat%Size_FR)) Catalogue_Constructed = .true.
!      if(allocated(Cat%Size_KSB)) Catalogue_Constructed = .true.
      if(allocated(Cat%ntile)) Catalogue_Constructed = .true.
      if(allocated(Cat%fluxerr)) Catalogue_Constructed = .true.
      if(allocated(Cat%magerr)) Catalogue_Constructed = .true.
      if(allocated(Cat%RA)) Catalogue_Constructed = .true.
      if(allocated(Cat%Dec)) Catalogue_Constructed = .true.
      if(allocated(Cat%g1)) Catalogue_Constructed = .true.
      if(allocated(Cat%g2)) Catalogue_Constructed = .true.
      if(allocated(Cat%Redshift)) Catalogue_Constructed = .true.

    end function Catalogue_Constructed

    subroutine Catalogue_Assign_byGalaxy_byIndex(Cat, Galaxy_Index, flux, mag, xpos, ypos, Sizes,Status)
      !-Assigns a single value by entered index. Function is overloaded with equivalent which assigns by Catalogue ()-!
      type(Catalogue)::Cat
      integer,intent(in)::Galaxy_Index
      real(double),intent(in),optional::flux, mag, xpos,ypos, Sizes
      integer,optional::Status ! "-1 Failed", "0 Unresolved", "1 Success" !

      integer::iStatus

!!$  INTERFACE
!!$     subroutine Catalogue_Assign_byGalaxy_byIndex(Cat, Galaxy_Index, flux, mag, xpos,ypos, Size_FWHM, Size_FR, Size_KSB, Status)
!!$       use Param_Types
!!$       type(Catalogue)::Cat
!!$       integer,intent(in)::Galaxy_Index
!!$       
!!$       real(double),intent(in),optional::flux, mag, xpos,ypos, Size_FWHM, Size_FR, Size_KSB
!!$       integer,optional::Status
!!$     end subroutine Catalogue_Assign_byGalaxy_byIndex
!!$  END INTERFACE


      if(Galaxy_Index > size(Cat%mag)) then; print *, 'Catalogue_Assign_byGalaxy - Index Entered not valid - Too Large'; iStatus = -1; return;
      elseif(Galaxy_Index <= 0) then; print *, 'Catalogue_Assign_byGalaxy - Index Entered not valid - Too Small'; iStatus = -1; return;
      end if

      if(present(Status)) Status = iStatus
      if(iStatus < 0) return

      if(present(flux)) Cat%flux(Galaxy_Index) = flux
      if(present(mag)) Cat%mag(Galaxy_Index) = mag
      if(present(xpos)) Cat%xpos(Galaxy_Index) =xpos
      if(present(ypos)) Cat%ypos(Galaxy_Index) =ypos
      if(present(Sizes)) Cat%Sizes(Galaxy_Index) =Sizes
      !if(present(flux)) Cat%Size_FR(Galaxy_Index) =Size_FR
      !if(present(flux)) Cat%Size_KSB(Galaxy_Index) =Size_KSB

    end subroutine Catalogue_Assign_byGalaxy_byIndex

    subroutine Catalogue_Assign_byGalaxy_byCatalogue(Cat, Index, Cat_Ref, Index_Ref)
      type(Catalogue), intent(inout)::Cat
      integer,intent(in)::Index
      type(Catalogue), intent(in)::Cat_Ref
      integer,intent(in)::Index_Ref
      
      integer::iStatus

      if(Index > size(Cat%RA)) then; print *, 'Catalogue_Assign_byGalaxy_byCatalogue - Index Entered not valid - Too Large', Index, size(Cat%RA); iStatus = -1; return;
      elseif(Index <= 0) then; print *, 'Catalogue_Assign_byGalaxy - Index Entered not valid - Too Small'; iStatus = -1; return;
      end if

      if(Index_Ref > size(Cat_Ref%RA)) then; print *, 'Catalogue_Assign_byGalaxy - Index_Ref Entered not valid - Too Large'; iStatus = -1; return;
      elseif(Index_Ref <= 0) then; print *, 'Catalogue_Assign_byGalaxy - Index_Ref Entered not valid - Too Small'; iStatus = -1; return;
      end if

      Cat%flux(Index) = Cat_Ref%flux(Index_Ref)
      Cat%fluxerr(Index) = Cat_Ref%fluxerr(Index_Ref)
      Cat%mag(Index) = Cat_Ref%mag(Index_Ref)
      Cat%magerr(Index) = Cat_Ref%magerr(Index_Ref)
      Cat%xpos(Index) = Cat_Ref%xpos(Index_Ref)
      Cat%ypos(Index) = Cat_Ref%ypos(Index_Ref)
      Cat%RA(Index) = Cat_Ref%RA(Index_Ref)
      Cat%Dec(Index) = Cat_Ref%Dec(Index_Ref)
      Cat%Sizes(Index) = Cat_Ref%Sizes(Index_Ref)
      Cat%Physical_Sizes(Index) = Cat_Ref%Physical_Sizes(Index_Ref)
      !Cat%Size_FR(Index) = Cat_Ref%Size_FR(Index_Ref)
      !Cat%Size_KSB(Index) = Cat_Ref%Size_KSB(Index_Ref)
      Cat%g1(Index) = Cat_Ref%g1(Index_Ref)
      Cat%g2(Index) = Cat_Ref%g2(Index_Ref)
      Cat%Redshift(Index) = Cat_Ref%Redshift(Index_Ref)

    end subroutine Catalogue_Assign_byGalaxy_byCatalogue

    function get_Total_Corrected_Ellipticity_g(Cat)
      type(Catalogue)::Cat
      
      real(double),dimension(:),allocatable::get_Total_Corrected_Ellipticity_g

      integer::i

      if(Catalogue_Constructed(Cat)==.false.) STOP 'FATAL ERROR - et_Total_Corrected_Ellipticity_g - Catalogue entered not constructed'

      allocate(get_Total_Corrected_Ellipticity_g(size(Cat%g1))); get_Total_Corrected_Ellipticity_g = 0.e0_double

      do i = 1, size(Cat%g1)
         get_Total_Corrected_Ellipticity_g(i) =dsqrt( Cat%g1(i)*Cat%g1(i) + Cat%g2(i)*Cat%g2(i) )
      end do
    
    end function get_Total_Corrected_Ellipticity_g

    !"---------------------------------------DEPRECATED CODE------------------------------------------------------------------------!


    !--COMBO17 Redshift Routines---!

    subroutine match_COMBO17_Redshifts(Cat, return_only_Matches)
      !-Matches redshift information from COMBO17 file to those galaxies contained in the STAGES catalogue-!
      !-30 Aug 13 -- this version uses Catherine's pre-matched file. Deprecated version which uses the COMBO17 data in Deprecated_Code File
      !-if "return_only_Matches" is present and true, the catalogue is reduced to only those glaxies with redshift information-!
      type(Catalogue)::Cat
      logical,optional::return_only_Matches

      logical::here
      character(120)::COMBO17_Filename = trim(adjustl(Cat_Dir))//'STAGES_shear_pz_matched.dat'
      real(double)::COMBO17_Input(10788, 15)

      integer::nMatch_Single, nMatch_Total
      integer::n_Multiple_Matches, nFailed_Matches
      integer::i,j, COMBO_Index, STAGES_Index, counter
      integer::COMBO_z_Index = 6
      character(200)::COMBO17_Line

      inquire(file = COMBO17_Filename, exist = here)
      if(here == .false.) then
         print *, 'Filename:', COMBO17_Filename
         STOP 'match_COMBO17_Redshifts - COMBO17 redshift file not found - check existence, stopping'
      end if

      !--Read in the COMBO17 file--!
      print *, 'Reading in COMBO17 file...'
      open(unit = 35, file = trim(adjustl(COMBO17_Filename)))
      !-Read Line by Line-!
      i = 0; j = 0; counter = 0
      do while (j < size(COMBO17_Input,1))
         counter = counter + 1
         if(counter <= 16) then
            read(35, *) COMBO17_Line
            cycle
         end if
         j = j+1
         read(35, *, IOSTAT=i) COMBO17_Input(j,:)
         if(i < 0) exit
      end do
      if(j /= size(COMBO17_Input,1)) STOP 'ERROR - COMBO17 input - Input File Length not as expected'
      i = 0; j = 0
      close(35)
      print *, 'Done.'

      if(allocated(Cat%Redshift)) deallocate(Cat%Redshift)
      allocate(Cat%Redshift(size(Cat%Sizes))); Cat%Redshift = -1.e0_double !-Default < 0, meaning unassigned-!
      nMatch_Total = 0; n_Multiple_Matches = 0; nFailed_Matches =  0
      do COMBO_Index = 1, size(COMBO17_Input,1) !-Loop through all COMBO17 Galaxies-!
         nMatch_Single = 0
         do STAGES_Index = 1, size(Cat%RA)
            !-Match on both RA and Dec to within a tolerance-!
            if( (dabs(Cat%RA(STAGES_Index)-COMBO17_Input(COMBO_Index,11)) <= 1.e-8_double) .and. (dabs(Cat%Dec(STAGES_Index)-COMBO17_Input(COMBO_Index,12)) <= 1.e-8_double)) then
               nMatch_Single = nMatch_Single + 1
               if(Cat%Redshift(STAGES_Index) < 0.e0_double) nMatch_Total = nMatch_Total + 1
               Cat%Redshift(STAGES_Index) = COMBO17_Input(COMBO_Index,COMBO_z_Index)
               
               !--Test shear measurements to ensure they match those held in my code--!
               if(dabs(Cat%g1(STAGES_Index) - COMBO17_Input(COMBO_Index,14)) > 1.e-0_double) PRINT *, 'WARNING - SHEARS (g1) DO NOT MATCH FOR MATHCED GALAXY (combo_index, stages_index):', COMBO_Index, STAGES_Index, Cat%g1(STAGES_Index), COMBO17_Input(COMBO_Index,14)
!               if(dabs(Cat%g2(STAGES_Index) - COMBO17_Input(COMBO_Index,15)) > 1.e-0_double) PRINT *,  'WARNING - SHEARS (g2) DO NOT MATCH FOR MATHCED GALAXY (combo_index, stages_index):', COMBO_Index, STAGES_Index
               exit
            end if
         end do
         if(nMatch_Single == 0) then
!            print *, 'WARNING - A COMBO17 Galaxy has not been matched to a STAGES galaxy'
            nFailed_Matches = nFailed_Matches + 1
            print *, 'Catherines mathced galaxy, line:', COMBO_Index+16, ' failed to find a match', COMBO17_Input(COMBO_Index,:)
         end if
         if(nMatch_Single > 1) then
!            print *, 'WARNING - A COMBO17 Galaxy has been matched to MULTIPLE STAGES galaxies'
            n_Multiple_Matches = n_Multiple_Matches + 1
         end if
      end do
      if(nFailed_Matches > 0) print *, 'WARNING:', nFailed_Matches,' COMBO17 galaxies failed to find a match'
      if(n_Multiple_Matches > 0) print *, 'WARNING:', n_Multiple_Matches,' COMBO17 galaxies found multiple matches'
      print *, nMatch_Total, ' of ', size(COMBO17_Input,1), ' galaxies where matched to STAGES galaxies'
      STOP 'STOPPING FOR TESTING'

      if((present(return_only_Matches)) .and. (return_only_Matches == .true.)) call reduce_Catalogue_toGalaxieswithRedshift(Cat)

    end subroutine match_COMBO17_Redshifts

    subroutine reduce_Catalogue_toGalaxieswithRedshift(Cat)
      !--Reduces the input catalogue to only those with redshift information. Matching is already assumed, and those with Redshift > 0 are considered as having redshift information--!
      type(Catalogue)::Cat
    

    end subroutine reduce_Catalogue_toGalaxieswithRedshift




end module Catalogues
