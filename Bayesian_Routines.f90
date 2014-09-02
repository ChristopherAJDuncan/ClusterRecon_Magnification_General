!~~To DO: Adaptive method ofr finding mode of posterior, moving towards peak, checking gradient on wither side until peak is found to within a ceertain tolerance. Would require that each cluster is evaluated on a different grid range (of alpha), and possbile complications in setting the order fo the grid
!~~This would need to address the issue on whether renormalisation of the posterior for each galaxy is important, aseach point would need combined individually

module  Bayesian_Routines
  !--Uses Bayesian Means to Calculate Posterior for the Dark Matter Profile variables--!
  !--Method must have redshifts for every galaxy (20Jan2014)--!
  use Param_Types; use Catalogues
  implicit none

  character(500):: Bayesian_Routines_Output_Directory
  
  real(double)::Default_Source_Redshift = 1.4e0_double
  logical:: Analyse_with_Physical_Sizes = .false. !#Default Value#

  logical::Combine_log_Posteriors = .true. !-If False, then combined by multiplication-!

  logical,private::Debug_Mode = .true.

  integer::Surface_Mass_Profile = 3 !-1:Flat, 2:SIS, 3:NFW-!

  !--Method: 1: Size-Only, 2: Size-Magnitude, 3:Magnitude Only, 4: SizeMag - MagOnly Combination (SizeLimits)--!
  integer:: Posterior_Method = 2

  logical:: Enforce_Weak_Lensing = .false.
  real(double):: Core_Cut_Radius(4) = 0.5e0_double !--Default
!  real(double):: Core_Cut_Radius(4) = (/1.50e0_double,1.50e0_double,0.5e0_double, 0.9e0_double/) !-In Arcminutes-! Data

  logical:: use_KDE_Smoothed_Distributions = .true., KDE_onTheFly = .false., allow_KDE_Extrapolation = .false.
  logical::use_lnSize_Prior = .false.
  logical:: Cuts_Renormalise_Likelihood = .true.
  real(double),dimension(2):: Survey_Magnitude_Limits = (/23.e0_double, 27.5e0_double/), Survey_Size_Limits = (/0.e0_double, 100.e0_double/)    
  real(double),dimension(2):: Prior_Magnitude_Limits = (/23.e0_double, 27.5e0_double/), Prior_Size_Limits = (/0.0e0_double, 100.e0_double/) !-3.3


contains

  subroutine Combine_Posteriors(PosteriorGrid, Posteriors, Combine_by_ln, Combined_Posterior)
    !--Combines posteriors by looping over the first dimension--!
    real(double), intent(in):: PosteriorGrid(:),Posteriors(:,:)
    real(double), intent(out):: Combined_Posterior(:)
    logical, intent(in):: Combine_by_ln

    integer::c, j
    real(double)::Renorm, Combination_Normalisation

    integer::nPosteriorsSkipped

    INTERFACE
       subroutine Combine_Posteriors(PosteriorGrid, Posteriors, Combine_by_ln, Combined_Posterior)
         use Param_Types
         real(double), intent(in):: PosteriorGrid(:), Posteriors(:,:)
         real(double), intent(out):: Combined_Posterior(:)
         logical, intent(in):: Combine_by_ln
       END subroutine Combine_Posteriors
    END INTERFACE

    nPosteriorsSkipped = 0
    if(Combine_by_ln == .false.) then
       Combination_Normalisation = 1.e0_double/maxval(Posteriors)!or 1.e0_double/(0.5e0_double*maxval(Posteriors(c,:))) within loop
       Combined_Posterior = 1.e0_double
       do c = 1, size(Posteriors,1) !-Loop over galaxies-!   
          if(all(Posteriors(c,:) == 0.e0_double) .or. all(isNaN(Posteriors(c,:)))) then
             nPosteriorsSkipped = nPosteriorsSkipped + 1
             cycle
          end if

          if(all(Posteriors(c,:)*Combination_Normalisation == 0.e0_double)) then
             print *, 'Invalid Posterior for galaxy:', c, ' (==0) press [ENTER] to output an stop..'; READ(*,*)
             print *, Posteriors(c,:)
             STOP
          end if
          if(any(Posteriors(c,:)*Combination_Normalisation < 0.e0_double)) then
             print *, 'Invalid Posterior for galaxy:', c, ' (<0) presS [ENTER] to output an stop..'; READ(*,*)
             print *, Posteriors(c,:)
             STOP
          end if
            
          Combined_Posterior(:) = Combined_Posterior(:)*(Posteriors(c,:)*Combination_Normalisation)

          if(all(Combined_Posterior == 0.e0_double)) then
             print *, 'Invalid CPosterior for galaxy:', c, ' press (==0) [ENTER] to output an stop..'; READ(*,*)
             print *, Combination_Normalisation
             read(*,*)
             print *, Combined_Posterior(:)
             read(*,*)
             print *, Posteriors(c,:)
             STOP
          end if
          if(anY(Combined_Posterior < 0.e0_double)) then
             print *, 'Invalid CPosterior for galaxy:', c, ' press (<0) [ENTER] to output an stop..'; READ(*,*)
             print *, Posteriors(c,:)
             STOP
          end if
          
       end do
    else
       print *, 'Combining using logs'
       Combined_Posterior = -100_double
       do c = 1, size(Posteriors,1) !-Loop over galaxies-!
          !-Error Catching--!
          if(all(Posteriors(c,:) == 1.e-75_double) .or. all(isNaN(Posteriors(c,:)))) then
             nPosteriorsSkipped = nPosteriorsSkipped + 1
             cycle
          end if
          !--Sum log posteriors--!
          where(Posteriors(c,:) == 0.e0_double)
             Combined_Posterior = Combined_Posterior - 100.e0_double
          elsewhere
             Combined_Posterior = Combined_Posterior + dlog(Posteriors(c,:))
          end where

             if(any(isNAN(Combined_Posterior(:)))) then
                print *, 'Any NaNs in Combined Posterior?, galaxy:',c, any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
                STOP
             end if
       end do
       !--Convert to PDF, not ln(PDF)--!
       Combined_Posterior = dexp(Combined_Posterior - maxval(Combined_Posterior))
    end if
 
    if(any(isNAN(Combined_Posterior(:)))) then
       print *, 'Any NaNs in Combined Posterior?', any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
       print *, 'Stopping'
       STOP
    END if

    if(nPosteriorsSkipped > 0) write(*,'(A,I4,A,I4,A)') '################### Combine_Posteriors - ', nPosteriorsSkipped, ' of ', size(Posteriors,1), ' posteriors were skipped as they we invlaid - NaNs or 0 ###########################'
    if((1.e0_double*nPosteriorsSkipped)/size(Posteriors,1) > 0.1) STOP 'Combine_Posteriors - number of skipped posteriors too large, stopping!'

    !--Renormalise--!                                                                                                                                     
    Renorm = 0.e0_double
    do j = 1, size(Combined_Posterior)-1
       Renorm = Renorm + 0.5e0_double*(Combined_Posterior(j) + Combined_Posterior(j+1))*(PosteriorGrid(j+1)-PosteriorGrid(j))
    end do
    if(Renorm <= 0.e0_double) then
       print *, 'Renormalisation:', Renorm
       STOP 'Combine_Posteriors - Invalid Renormalisation for combined Posterior'
    end if
    Combined_Posterior(:) = Combined_Posterior(:)/Renorm

    if(any(isNAN(Combined_Posterior(:)))) then
       print *, 'Any NaNs in Renormalised Combined Posterior?', any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
       print *, 'Stopping'
       STOP
    END if

  end subroutine Combine_Posteriors

  subroutine Posterior_Statistics(PosteriorGrid, Posterior, MeanVal, ModeVal, Error, AntiSymm_Error)
    use Statistics, only: mean, mode_distribution, variance_distribution, Antisymmetric_Variance_Distribution
    real(double),intent(in)::Posterior(:), PosteriorGrid(:)
    real(double), intent(out),optional::MeanVal, ModeVal, Error, AntiSymm_Error(2)

    INTERFACE
       subroutine Posterior_Statistics(PosteriorGrid, Posterior, MeanVal, ModeVal, Error, AntiSymm_Error)
         use Param_Types
         real(double),intent(in)::Posterior(:), PosteriorGrid(:)
         
         real(double), intent(out),optional::MeanVal, ModeVal, Error
         real(double), intent(out),optional:: AntiSymm_Error(2)
       end subroutine Posterior_Statistics
    END INTERFACE

    if(present(ModeVal)) ModeVal = mode_distribution(PosteriorGrid, Posterior)
    if(present(MeanVal)) MeanVal = mean(Posterior, PosteriorGrid)
    if(present(Error)) Error = dsqrt(variance_distribution(PosteriorGrid, Posterior))
    if(present(AntiSymm_Error)) then
       if(present(ModeVal)) then
          !--Take error about mode as a preference--!
          AntiSymm_Error = Antisymmetric_Variance_Distribution(PosteriorGrid, Posterior, ModeVal)
       elseif(present(Meanval)) then
          !--Both the following use the mean value, however re-use if already calculated--!
          AntiSymm_Error = Antisymmetric_Variance_Distribution(PosteriorGrid, Posterior, MeanVal)
       else
          AntiSymm_Error = Antisymmetric_Variance_Distribution(PosteriorGrid, Posterior)
       end if
    end if

  end subroutine Posterior_Statistics

  subroutine Maximise_Convergence_byShifts_inAperture(Cat, BFCat, Ap_Pos, Ap_Radius, Core_Radius)
    use Statistics, only: mean_discrete
    type(Catalogue), intent(in)::Cat, BFCat
    real(double),intent(in)::Ap_Pos(:), Ap_Radius
    real(double), intent(in), Optional:: Core_Radius
    
    type(Catalogue)::Ap_Cat
    integer:: i,j
    integer,parameter:: nSide = 2500
    real(double)::Start_Pos(2)
    real(double):: Try_Pos(2), dTheta
    real(double):: Global_Mean_Size, Aperture_Mean_Size, Global_Mean_Mag, Aperture_Mean_Mag
    real(double), dimension(nSide, nSide, 2):: Convergence

    integer:: Maximum_Position(2)

    !--Cuts on Prior--!

    Global_Mean_Size = mean_discrete(BFCat%Sizes)
    Global_Mean_Mag =  mean_discrete(BFCat%MF606W)

    Start_Pos = Ap_Pos - Ap_Radius

    dTheta = Ap_Radius/(nSide-1)

    do i = 1, nSide
       Try_Pos(1) = Start_Pos(1) + (i-1)*dTheta
       do j = 1, nSide
          Try_Pos(2) = Start_Pos(2) + (j-1)*dTheta

          call Identify_Galaxys_in_Circular_Aperture(Cat, Try_Pos, Ap_Radius, Ap_Cat)

          Aperture_Mean_Size = mean_discrete(Ap_Cat%Sizes)
          Aperture_Mean_Mag = mean_discrete(Ap_Cat%MF606W)

          call Catalogue_Destruct(Ap_Cat)

          Convergence(i,j,1) = (Aperture_Mean_Size/Global_Mean_Size) - 1.e0_double !-Size-!
          Convergence(i,j,2) = (Global_Mean_Mag-Aperture_Mean_Mag)/2.17e0_double !-Magnitude-!
       end do
    end do
    !--In below, if not dabs, then maximum positive convergence is returned--!
    print *, 'Start Pos for Aperture:', Ap_Pos
    !--Output Maximum Pos for Mag--!
    Maximum_Position = maxloc(Convergence(:,:,2))
    print *, '** Maximum convergence for magnitude at:', Start_Pos(1)+(Maximum_Position(1)-1)*dTheta, Start_Pos(2)+(Maximum_Position(2)-1)*dTheta
    print *, '   which has value:', Convergence(Maximum_Position(1), Maximum_Position(2), 2)

    !--Output Maximum Position for Size--!
    Maximum_Position = maxloc(Convergence(:,:,1))
    print *, '** Maximum convergence for magnitude at:', Start_Pos(1)+(Maximum_Position(1)-1)*dTheta, Start_Pos(2)+(Maximum_Position(2)-1)*dTheta
    print *, '   which has value:', Convergence(Maximum_Position(1), Maximum_Position(2), 1)

    !--Output Maximum Position for both--!

  end subroutine Maximise_Convergence_byShifts_inAperture



  subroutine Convert_Alpha_Posteriors_to_VirialMass(SMD_Profile, Posteriors, Mass_Posteriors, Lens_Redshift, Output_Label)
    !--Converts from posteriors on DM free parameter (alpha) to posteriors on virial mass using the conservation of probability
    use Mass_Profiles, only: Halo_Mass, virial_Radius_from_ProfileFreeParameter; use Integration, only: Integrate
    
    integer,intent(in):: SMD_Profile
    real(double), intent(in):: Posteriors(:,:) !-Posteriors(1,:) assumed to be grid of DM Free parameters, (2,:) the actual posterior
    real(double), intent(out):: Mass_Posteriors(:,:)
    real(double):: Lens_Redshift !-Used only in NFW-!
    character(*), intent(in), optional:: Output_Label

    integer:: i
    real(double):: Discardable(1)
    character(20)::fmt
    real(double):: Redshift

    !-Posteriors entered must be in the form of the output mass posterior: that is, Dimension 1 labels grid/cluster, so that element A of Dim. 1 is Grid if A==1, or Cluster A-1 if A>1

    if(size(Posteriors,1) < 2) STOP 'Convert_VirialRadius_Posteriors_to_VirialMass - Error in input posterior - no grid/value'
    if(size(Mass_Posteriors) /= size(Posteriors)) STOP 'Convert_VirialRadius_Posteriors_to_VirialMass - Error in size of mass posterior, not equal to Original posterior' !Can be deleted in the case of interpolation

    print *, 'Converting alpha posterior to Mass Posterior'

!!$    if(size(Lens_Redshift) /= 1 .and. size(Lens_Redshift) /= size(Redshift)) STOP 'Convert_Alpha_Posteriors_to_VirialMass - Either enter 1 redshift or a redshift for each cluster'
!!$    if(size(Lens_Redshift) == 1) then
!!$       Redshift = Lens_Redshift(1)
!!$    else
!!$       Redshift = Lens_Redshift
!!$    end if
    Redshift = Lens_Redshift


    !--Assumes that all the clusters entered are at teh same redshift--!
    do i = 1, size(Posteriors,2)
       call Halo_Mass(SMD_Profile, Posteriors(1,i), (/0.e0_double/), Mass_Posteriors(1,i), Discardable, Redshift)
       
       !--p(M)dM = p(r)dr -> p(M) \propto p(r)/r^2 -- Proportionality just affects renormalisation
!       Mass_Posteriors(2:,i) = Posteriors(2:,i)/(virial_Radius_from_ProfileFreeParameter(SMD_Profile, Posteriors(1,i))**2.e0_double)
       Mass_Posteriors(2:,i) = Posteriors(2:,i)/(Mass_Posteriors(1,i)**(2.e0_double/3.e0_double))
    end do

    print *, '**Mass Posterior output on grid of 10^14 Msun/h'
    Mass_Posteriors(1,:) = Mass_Posteriors(1,:)/1.e14_double

    !--Renormalise
    do i = 2, size(Posteriors,1)
       Mass_Posteriors(i,:) = Mass_Posteriors(i,:)/Integrate(Mass_Posteriors(1,:), Mass_Posteriors(i,:), 2, lim = (/minval(Mass_Posteriors(1,:)), maxval(Mass_Posteriors(1,:))/))
    end do

    if(present(Output_Label)) then
       open(unit = 78, file = trim(Output_Label)//'VirialMass_Posterior.dat')
       write(fmt, *) size(Mass_Posteriors,1)
       do i = 1, size(Mass_Posteriors,2)
          write(78, '('//trim(fmt)//'(e9.3,x))') Mass_Posteriors(:,i)
       end do
       write(*,'(2A)') '* Output file to: ', trim(Output_Label)//'VirialMass_Posterior.dat'
       close(78)
    end if
    print *, '** Virial Mass Posterior output to: ', trim(Output_Label)//'VirialMass_Posterior.dat'

    print *, '---Finished Conversion'

  end subroutine Convert_Alpha_Posteriors_to_VirialMass

  subroutine DM_Profile_Variable_Posteriors_CircularAperture(Cat, Ap_Pos, Ap_Radius, Posteriors, Distribution_Directory, reproduce_Prior, Blank_Field_Catalogue)
    use Statistics, only: mean_discrete, mode_distribution, variance_distribution; use Mass_profiles
    use Distributions; use Cosmology, only: angular_diameter_distance_fromRedshift;
    !--Main routine that returns the Posteriors on DM free parameter (alpha) over all apertures.
    !--Ap_Radius in DEGREES
    !--If Blank_Field_Catalogue is entered, then intrinsic distributions are produced form this Catalogue
    !--To Do: 
    !~~Recent edits to use a joint size magnitude have ignored magnitude binning


    type(Catalogue), intent(in)::Cat
    real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:)
    real(double),intent(out),allocatable::Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-! - Second dimension allows for a different grid for each aperture
    character(*), intent(in):: Distribution_Directory
    logical, intent(in):: reproduce_Prior
    type(Catalogue), intent(in),optional::Blank_Field_Catalogue


    real(double),dimension(Size(Ap_Pos,1))::iAp_Radius

    type(catalogue),dimension(size(Ap_pos,1))::Ap_Cats
    integer::i, ap, m, j

    !--Size Grid Declarations--!     
    integer::nSizeGrid = 100        
!    real(double)::SizeGrid_Lower = 0.e0_double, SizeGrid_Upper = 70.e0_double, dSizeGrid !-Pixel-! 
    real(double)::SizeGrid_Lower = 0.e0_double, SizeGrid_Upper = 0.0035e0_double, dSizeGrid
    real(double),allocatable::SizeGrid(:)

    real(double),allocatable:: Joint_Size_Magnitude_Distribution(:,:), Magnitude_Distribution(:)

    real(double),allocatable::SizePrior_byMag(:,:) !-MagBin, GridValue-!
    real(double),allocatable::Aperture_Posterior_byMag(:,:,:) !-MagBin, Grid/Posterior, Value-! 

    type(Binned_Catalogue)::BCat
    real(double),allocatable::Posterior_Single(:,:) !-Grid/Posterior, Value-!
    real(double),allocatable::Prior_Single(:,:) !!-Grid/Prior, Value-! 

    real(double)::Lens_Redshift = 0.165e0_double
    real(double)::Renorm

    character(2)::fmtString, apString

    !--Conversion to Mass--!
    real(double)::D_l, Area
    real(double),allocatable::Cluster_Mean(:), Cluster_Variance(:), Cluster_Mode(:), AntiSymm_Variance(:,:)
    real(double):: Scale_Mass, Scale_Mass_Error(2)
    character(10):: Mass_String, Error_Mass_String_Positive, Error_Mass_String_Negative

    !--MAgnitude Binning (by Absolute Magnitude)--!
    integer::nMag = 1
    real(double),allocatable::MagBins(:,:)
    integer::Magnitude_Binning_Type = 1 !-1:Absolute, 2:Apparent (MF606)-!

    character(500):: Output_File_Prefix

    real(double),allocatable::MagGrid(:) !--Discardable for now as marginalised over--!

    type(catalogue)::TCat

    INTERFACE
         subroutine DM_Profile_Variable_Posteriors_CircularAperture(Cat, Ap_Pos, Ap_Radius, Posteriors, Distribution_Directory, reproduce_Prior, Blank_Field_Catalogue)
           use Param_Types; use Catalogues
           !--Main routine that returns the Posteriors over all apertures.
           !--Ap_Radius in DEGREES
           type(Catalogue), intent(in)::Cat
           real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:)
           real(double),intent(out),allocatable::Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-! - Second dimension allows for a different grid for each aperture
           character(*), intent(in):: Distribution_Directory
           logical, intent(in):: Reproduce_Prior

           type(Catalogue), intent(in),optional::Blank_Field_Catalogue
         end subroutine DM_Profile_Variable_Posteriors_CircularAperture
      end INTERFACE

    if(size(Ap_Radius)==1) then
       iAp_Radius = Ap_Radius(1)
    else
       iAp_Radius = Ap_Radius
    end if

    !--Get Mean of size and magnitude distributions in masked apertures
    call Mask_Circular_Aperture(TCat, Ap_Pos, Core_Cut_Radius/60.e0_double)
    TCat = Cat

    print *, '**Global catalogue has mean [Core mask accounted for] (Size, Mag):', mean_discrete(TCat%Sizes), mean_discrete(TCat%MF606W)
    call Catalogue_Destruct(TCat)

    !--Identify Reduced Catalogue for each aperture--!
    do i =1, size(Ap_Cats)
       print *, 'Searching for position of maxima in Size-Magnitude Shifts for Cluster:', i
!       call Maximise_Convergence_byShifts_inAperture(Cat, Blank_Field_Catalogue, Ap_Pos(i,:), Ap_Radius(i))

       print *, 'Cutting on a core radius of:', Core_Cut_Radius(i), ' arcminutes for aperture:', i
       call Identify_Galaxys_in_Circular_Aperture(Cat, Ap_Pos(i,:), iAp_Radius(i), Ap_Cats(i), Core_Radius = Core_Cut_Radius(i)/60.e0_double)
       
       print *, '* Ap ', i ,' has mean (Size,Mag):',  mean_discrete(Ap_Cats(i)%Sizes), mean_discrete(Ap_Cats(i)%MF606W)
    end do

    do i =1, size(Ap_Cats)
       print *, 'Aperture:', i, ' contains:', size(Ap_Cats(I)%RA), ' galaxies'
    end do

    print *, count(Cat%Redshift >= 0.e0_double), ' of ',size(Cat%Redshift), ' galaxies have redshift information'

    !--If a catalogue for the blank field is passed in, then use this catalogue to get the intrinsic distributions--!
    if(KDE_OnTheFly == .false.) then
       if(reproduce_Prior) then
          if(present(Blank_Field_Catalogue) == .false.) STOP 'DM_Profile_Variable_Posteriors_CircularAperture - Blankf Field Catalogue must be entered to allow for the production of the prior on a grid'
          write(*,'(A)') 'Producing Distribution from Catalogue'
          call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, Distribution_Directory, Blank_Field_Catalogue)
       else
          write(*,'(A)') 'Reading in distribution from:', Distribution_Directory
          call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, Distribution_Directory)
          print *, 'Success:', Distribution_Directory
       end if
    else
       print *, ' '
       print *, 'PRIOR CONSTRUCTION TAKEN OUT FOR ON-THE-FLY'
       print *, ' '
    end if

    if(size(magBins,1) > 1) STOP 'I have had to disable Magnitude Binning for now, youll have to edit the code to get this to work, stopping'

    Do Ap = 1, size(Ap_Cats)
       write(Output_File_Prefix,'(I2)') Ap
       Output_File_Prefix = trim(adjustl(Bayesian_Routines_Output_Directory))//'Aperture_'//trim(adjustl(Output_File_Prefix))//'_'


!  subroutine DM_Profile_Variable_Posterior(Cat, Mass_Profile, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior, PriorCatalogue, PriorMagGrid, PriorSizeGrid, Prior)
       if(KDE_OnTheFly) then
          if(present(Blank_Field_Catalogue) == .false.) STOP 'KDE_OnTheFly can only be done if a prior catalogue is present'
          print *, 'Calling DM_Profile_Variable_Posterior', present(Blank_Field_Catalogue)
          call DM_Profile_Variable_Posterior(Ap_Cats(Ap), Surface_Mass_Profile, Lens_Redshift, Ap_Pos(Ap,:), Posterior_Single, Output_File_Prefix, use_lnSize_Prior, PriorCatalogue = Blank_Field_Catalogue)  
       else
          if(present(Blank_Field_Catalogue)) then
             call DM_Profile_Variable_Posterior(Ap_Cats(Ap), Surface_Mass_Profile, Lens_Redshift, Ap_Pos(Ap,:), Posterior_Single, Output_File_Prefix, use_lnSize_Prior, PriorCatalogue = Blank_Field_Catalogue, PriorMagGrid = MagGrid, PriorSizeGrid = SizeGrid, Prior = Joint_Size_Magnitude_Distribution, MagPrior = Magnitude_Distribution)  
          else
             call DM_Profile_Variable_Posterior(Ap_Cats(Ap), Surface_Mass_Profile, Lens_Redshift, Ap_Pos(Ap,:), Posterior_Single, Output_File_Prefix, use_lnSize_Prior, PriorMagGrid = MagGrid, PriorSizeGrid = SizeGrid, Prior = Joint_Size_Magnitude_Distribution, MagPrior = Magnitude_Distribution)  
          end if
       end if

       if(Ap == 1) allocate(Posteriors(size(Ap_Cats), 2, size(Posterior_Single,2)))
       Posteriors(Ap,:,:) = Posterior_Single
       deallocate(Posterior_Single)
    end Do
    if(allocated(SizeGrid)) deallocate(SizeGrid)
    if(allocated(MagGrid)) deallocate(MagGrid)
    if(allocated(Joint_Size_Magnitude_Distribution)) deallocate(Joint_Size_Magnitude_Distribution)

    !--Get Prior by getting refernce for each bin, and binning each reduced cat using the same definition (mag)--!
    !--Produce Posterior for each mag bin--!

!!    !$OMP PARALLEL DEFAULT(PRIVATE), SHARED(Ap_Cats, MagBins, Surface_Mass_Profile, MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Lens_Redshift,Ap_Pos, use_lnSize_Prior, Posteriors, Combine_log_Posteriors)
!!    !$OMP DO
!!$       !--Convert Posteriors in NFW from Virial Radius to VirialMass--!
!       if(Surface_Mass_Profile == 3) then
!          print *, '**Converting NFW posteriors from Virial radius to Virial Mass:....'
!          do m =1, size(Aperture_Posterior_byMag,1)
!             Aperture_Posterior_byMag(m,1,:) = get_NFW_VirialMass_from_VirialRadius(Lens_Redshift, Aperture_Posterior_byMag(m,1,:))
!          end do
!       end if
!!    !$OMP END PARALLEL
!!end parallel

    !--Output Posterior--!
    open(unit = 51, file = trim(Bayesian_Routines_Output_Directory)//'Posterior_per_Aperture.dat')
    write(fmtstring,'(I2)') size(Posteriors,1)+1 !-Assumes all posteriors described by the same posterior grid-!
    do j = 1, size(Posteriors,3)
       write(51, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posteriors(1,1,j), Posteriors(:,2,j)
    end do
    close(51)
    print *,'Output file to: ', trim(Bayesian_Routines_Output_Directory)//'Posterior_per_Aperture.dat'


  end subroutine DM_Profile_Variable_Posteriors_CircularAperture

  subroutine DM_Profile_Variable_Posterior(Cat, Mass_Profile, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior, PriorCatalogue, PriorMagGrid, PriorSizeGrid, Prior, MagPrior)
    use cosmology, only:angular_diameter_distance_fromRedshift; use MC_Redshift_Sampling, only: Monte_Carlo_Redshift_Sampling_SigmaCritical; use Mass_Profiles; use Distributions, only: ch08_redshift_distribution_Array, CH08_Redshift_Distribution_Scalar; use Interpolaters, only: Linear_Interp
    use Integration, only:TrapInt, Integrate; use Matrix_methods, only: Determinant, Matrix_Invert; use Smoothing, only: KDE_BiVariate_Gaussian_Scalar, KDE_UniVariate_Gaussian; use Statistics, only:Discrete_Covariance, mean_discrete
    !--Returns the Bayesian posterior for the DM Profile Variables, as defined by Mass_Profile
    !-For Mass_Profile = 1 (Flat): Sigma_0
    !-                   2 (SIS) : Velocity_Dispersion**2
    !-                   3 (NFW) : Virial Radius (r200)
    !--Prior is the prior distribution of p(m,R), or possibly p(m,lnR), m is *apparent magnitude*, R is *apparent size*
    !--Cat entered already assumed to be only the galaxies that are to be considered (e.g. in aperture), and well-described by posterior (e.g. same mag bin) 

    !---TO DO:
    !--Conversion from posterior to surface mass density (possibly not here, but later)

    type(Catalogue)::Cat
    integer,intent(in)::Mass_Profile !-1:Flat, 2:SIS, 3: NFW-!
    real(double),intent(in)::Lens_Redshift
    real(double),intent(in)::Lens_Position(2)
    real(double),allocatable,intent(Out)::Posterior(:,:) !-Grid/Posterior, Value-!
    character(*), intent(in)::Output_Prefix
    logical,intent(in)::lnSize_Prior
    type(Catalogue), intent(in), optional:: PriorCatalogue
    real(double),intent(in),optional:: PriorSizeGrid(:), PriorMagGrid(:), MagPrior(:), Prior(:,:) !-Magnitude, Size-!

    integer::i, c, j, z, m
    integer:: Galaxy_Posterior_Method

    !--Variable Grid Declarations-!
    integer::nGrid = 1000
    real(double)::VGrid_Lower, VGrid_Higher

    !-Size_Only_Prior contains the priors which depends only on size, which is evalutaed for each redshift and integrates over all magnitudes
    real(double),dimension(:),allocatable:: Size_Only_Prior, Mag_Prior, Mag_Only_Prior
    !--Size_Only_Mag_Prior contains the prior for size which evaluated over the mag grid, which will be integrated over. Contains p_[theta_0, m_0|z]*p[z|m_0] for all m in grid
    real(double),dimension(:,:),allocatable::  Kappa_Renormalised_Prior, Survey_Renormalised_Prior, Size_Only_Mag_Prior
    real(double),dimension(:),allocatable::Kappa_Renormalised_MagPrior, Survey_Renormalised_MagPrior
    real(double),allocatable::Posterior_perGalaxy(:,:) !-Galaxy, Posterior-! - Uses Same Grid as overall Posterior -
    real(double),allocatable:: Effective_Magnification(:,:) !-Galaxy, DM Free Parameter-! - Only Model Dependant

    logical:: MC_Sigma_Critical = .false.
    logical:: Marginalise_Redshift_Distribution = .true.

    logical::Known_Redshift

    !--KDE_Smoothing Declarations--!
    real(double),allocatable:: KDE_Gaussian_Covariance(:,:), Data_Vectors(:,:), KDE_Covariance_Inverse(:,:)
    real(double):: KDE_Covariance_Determinant, KDE_Gaussian_Covariance_Reduction = 0.01e0_double !-How much is sig^2 which give KDE width reduced from measured covariance?--!
    logical:: do_KDE_Extrapolation, do_KDE_OnTheFly
    logical:: need_Extrapolate

    real(double),allocatable::Sigma_Crit(:,:), Sigma_Crit_MC(:)
    real(double)::D_l, D_s, D_ls
    real(double)::Distance_from_Mass_Center

    real(double)::Renorm

    logical:: Output_Posterior_Per_Galaxy = .true.
    character(7)::fmtstring
    character(500)::Filename
    logical::here

    !--Redshift Distribution Declarations--!
    integer,parameter:: nRedshift_Sampling = 50
    real(double),parameter::Redshift_Lower = 0.156e0_double, Redshift_Higher = 4.e0_double !!!Edit to Lens_Redshift
    real(double), dimension(nRedshift_Sampling):: RedshiftGrid
    real(double),allocatable::RedshiftPDF(:)
    real(double),allocatable:: Posterior_perGalaxy_Redshift(:,:,:) !-Galaxy, Posterior, Redshift-!

    !--Kappa dependant renormalisation--!
!!    real(double),dimension(2):: Survey_Magnitude_Limits, Survey_Size_Limits !--User defined needs edit to below, these are set by default from distribution-!
    real(double),dimension(2):: Renormalisation_Magnitude_Limits, Renormalisation_Size_Limits
    real(double),allocatable:: MagnificationGrid(:), Renormalisation_by_Magnification(:), Convergence_Renorm_PerGalaxy(:,:),  MagOnly_Renormalisation_by_Magnification(:)
    !vv Must be set here vv!
    integer:: nMagnificationGrid = 2000
    real(double):: MagFactorGridLower = 1.e0_double, MagFactorGridHigher = 65.e0_double
    integer:: IntegrationFlag = -1000

    integer:: nMagPosterior, nSizePosterior, nSizeMagPosterior

    !--Testing Declarations--!
    real(double),dimension(3,size(Cat%RA)):: Convergence_per_Cluster
    integer:: n_Default_Source_Redshift_Used, nGal_Ignored_MagLimits, nGal_Ignored_SizeLimits, nGal_Ignored_NaN
    real(double), allocatable:: Aperture_Smoothed_Size_PDF(:), Aperture_Smoothed_Mag_PDF(:)
    logical:: Produce_Shift_inAperture = .false.
    real::Time1, Time2, Time3, Time4
    integer::test

    INTERFACE
       subroutine DM_Profile_Variable_Posterior(Cat, Mass_Profile, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior, PriorCatalogue, PriorMagGrid, PriorSizeGrid, Prior, MagPrior)
         use Param_Types; use Catalogues
         type(Catalogue)::Cat
         integer,intent(in)::Mass_Profile !-1:Flat, 2:SIS, 3:NFW-!
         real(double),intent(in)::Lens_Redshift
         real(double),intent(in)::Lens_Position(2)
         real(double),allocatable,intent(Out)::Posterior(:,:) !-Grid/Posterior, Value-!
         character(*), intent(in)::Output_Prefix
         logical,intent(in)::lnSize_Prior

         type(Catalogue), intent(in), optional:: PriorCatalogue
         real(double),intent(in),optional:: PriorSizeGrid(:), PriorMagGrid(:), MagPrior(:), Prior(:,:) !-Magnitude, Size-!
       END subroutine DM_Profile_Variable_Posterior
    END INTERFACE

    print *, 'Called DM_Profile_Variable_Posterior'

    if(Analyse_with_Physical_Sizes) STOP 'DM_Profile_Variable_Posterior - I HAVE DISABLED THE ABILITY TO USE PHYISCAL SIZES AS UNNECESSARY, code still to be edited'

    do_KDE_Extrapolation = .false.; do_KDE_OnTheFly = .false.
    if(present(Prior) .or. present(MagPrior)) then
       if((present(PriorMagGrid) == .false.) .or. (present(PriorSizeGrid)== .false.)) STOP 'DM_Profile_Variable_Posterior - Prior must be accompanied by grids'
       allocate(Kappa_Renormalised_Prior(size(Prior,1), size(Prior,2))); Kappa_Renormalised_Prior = 0.e0_double
       allocate(Survey_Renormalised_Prior(size(Prior,1),size(Prior,2))); Survey_Renormalised_Prior = 0.e0_double
       if(present(MagPrior)) then
          allocate(Survey_Renormalised_MagPrior(size(PriorMagGrid))); Survey_Renormalised_MagPrior = 0.e0_double
          allocate(Kappa_Renormalised_MagPrior(size(PriorMagGrid))); Kappa_Renormalised_MagPrior =0.e0_double
       end if
       allocate(Size_Only_Mag_Prior(size(Prior,1),size(Prior,2))); Size_Only_Mag_Prior = 0.e0_double
       if(allow_KDE_Extrapolation .and. present(PriorCatalogue)) then
          do_KDE_Extrapolation = .true.
          print *, ' '
          print *, 'Attempting Prior Interpolation with KDE Extrapolation'
          print *, ' '
       else
          print *, ' '
          print *, 'Attempting Prior Interpolation without Extrapolation'
          print *, ' '
       end if
    elseif(present(PriorCatalogue)) then
       do_KDE_OnTheFly = .true.
       print *, ' '
       print *, 'Attempting KDE on the Fly - NOTE This does not include Kappa-Renormalisation (Seg fault straight after this..)'
       print *, ' '
    end if

    !-Set Up Posterior Grid-!
    select case(Mass_Profile)
    case(1) !-Flat-!
       nGrid = 100000
       VGrid_Lower = -5.e-3_double; VGrid_Higher = 1.e-2_double !-- SMD ~ Masses 10^12 -> 10^15 Msun/h, in Units of 10^18 Msun/h  
    case(2) !-SIS-!
       nGrid = 10000 !--This needs to be sigma_v^2--!
       VGrid_Lower= -2.5e5_double; VGrid_Higher = 9.e6_double    !--------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case(3)
       nGrid = 500
       VGrid_Lower= 0.05e0_double; VGrid_Higher = 3.e0_double
    case default
       STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
    end select
    allocate(Posterior(2, nGrid)); Posterior = 1.e0_double !*!
    do i =1, nGrid
       Posterior(1,i) = VGrid_Lower + (i-1)*((VGrid_Higher-VGrid_Lower)/(nGrid-1))
    end do
    !----------------------!

    if(allow_KDE_Extrapolation .or. KDE_OnTheFly) then
       if(present(PriorCatalogue) == .false.) STOP 'DM_Profile_Variable_Posterior - Prior Catalogue needs to be entered to allow KDE Extrapolation'
       !--Construct the Covariance that will be used for the KDE Smoothing--!
       allocate(Data_Vectors(2,size(PriorCatalogue%Sizes))); Data_Vectors(1,:) = PriorCatalogue%MF606W; Data_Vectors(2,:) = PriorCatalogue%Sizes
       call Discrete_Covariance(Data_Vectors, KDE_Gaussian_Covariance)
       KDE_Gaussian_Covariance = KDE_Gaussian_Covariance_Reduction*KDE_Gaussian_Covariance
       call Matrix_Invert(KDE_Gaussian_Covariance, KDE_Covariance_Inverse, 'S')
       KDE_Covariance_Determinant = Determinant(KDE_Gaussian_Covariance)
       deallocate(KDE_Gaussian_Covariance)
    end if

    if(produce_Shift_inAperture) then
       !--Testing: output smoothed size and mag distributions within the aperture to see if there is any noticable shift
       allocate(Data_Vectors(2,size(PriorCatalogue%Sizes))); Data_Vectors(1,:) = PriorCatalogue%MF606W; Data_Vectors(2,:) = PriorCatalogue%Sizes
       call Discrete_Covariance(Data_Vectors, KDE_Gaussian_Covariance)
       KDE_Gaussian_Covariance = KDE_Gaussian_Covariance_Reduction*KDE_Gaussian_Covariance
       allocate(Aperture_Smoothed_Size_PDF(size(PriorSizeGrid))); Aperture_Smoothed_Size_PDF = 0.e0_double
       allocate(Aperture_Smoothed_Mag_PDF(size(PriorMagGrid))); Aperture_Smoothed_Mag_PDF = 0.e0_double
       call KDE_Univariate_Gaussian(Cat%Sizes, dsqrt(KDE_Gaussian_Covariance(2,2)), PriorSizeGrid, Aperture_Smoothed_Size_PDF)
       call KDE_Univariate_Gaussian(Cat%MF606W, dsqrt(KDE_Gaussian_Covariance(1,1)), PriorMagGrid, Aperture_Smoothed_Mag_PDF)
       open(unit = 31, file = trim(adjustl(Output_Prefix))//'KDE_Distributions_in_Aperture_Size.dat')
       do i =1, size(PriorSizeGrid)
          write(31, *) PriorSizeGrid(i), Aperture_Smoothed_Size_PDF(i)
       end do
       close(31)
       open(unit = 31, file = trim(adjustl(Output_Prefix))//'KDE_Distributions_in_Aperture_Mag.dat')
       do i =1, size(PriorSizeGrid)
          write(31, *) PriorMagGrid(i), Aperture_Smoothed_Mag_PDF(i)
       end do
       close(31)
       write(*, '(4A)') '**Output distributions in aperture to: ', trim(adjustl(Output_Prefix))//'KDE_Distributions_in_Aperture_Size.dat', ' : ', trim(adjustl(Output_Prefix))//'KDE_Distributions_in_Aperture_Mag.dat'
       print *, '**with mean (Size, Mag):', mean_discrete(Cat%Sizes), mean_discrete(Cat%MF606W) 
       deallocate(Aperture_Smoothed_Size_PDF, Aperture_Smoothed_Mag_PDF)
    !---End of distributions in Aperture
    end if

    D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
    !--SET UP REDSHIFT GRID--!
    if(Marginalise_Redshift_Distribution) then
       do z = 1, nRedshift_Sampling
          RedshiftGrid(z) = Redshift_Lower + (z-1)*((Redshift_Higher-Redshift_Lower)/(nRedshift_Sampling-1))
       end do
    end if


    !--Renormalise the prior within these size and magnitude limits--!
    print *, 'Renormalisation of the intrinsic distribution:', Integrate(PriorMagGrid, PriorSizeGrid, Prior, 2, lim1 = Survey_Magnitude_Limits, lim2 = Survey_Size_Limits)
    Survey_Renormalised_Prior = Prior/Integrate(PriorMagGrid, PriorSizeGrid, Prior, 2, lim1 = Survey_Magnitude_Limits, lim2 = Survey_Size_Limits)

    if(present(MagPrior)) then
       !--Allow for seperate renormalisation of the magnitude prior. This is required if the magnitude prior is constructed from galaxies which are excluded from the joint size-magnitude analysis, e.g. due to size cuts--!
       print *, 'Renormalisation of the intrinsic magnitude distribution:', Integrate(PriorMagGrid, MagPrior, 2, lim = Survey_Magnitude_Limits)
       Survey_Renormalised_MagPrior = MagPrior/Integrate(PriorMagGrid, MagPrior, 2, lim = Survey_Magnitude_Limits)
    end if

    !------This section renormalises the likelihhod, taking into account kappa-dependent mag and size cuts, however this need to be implemented in the prior-----!
    allocate(MagnificationGrid(nMagnificationGrid)); MagnificationGrid = 0.e0_double
    allocate(Renormalisation_by_Magnification(size(MagnificationGrid))); Renormalisation_by_Magnification = 0.e0_double
    if(present(MagPrior)) then
       allocate(MagOnly_Renormalisation_by_Magnification(size(MagnificationGrid))); MagOnly_Renormalisation_by_Magnification = 0.e0_double
    end if
    allocate(Size_Only_Prior(size(Prior,2))); Size_Only_Prior = 0.e0_double
    !Choose MagGrid lower to be the point where the full magnitude range is swept out:
    MagFactorGridHigher = 1.05e0_double*(10.e0_double**((Survey_Magnitude_Limits(2)-Survey_Magnitude_Limits(1))/2.5e0_double)) !1.05 gives lee-way!
    do i = 1, size(MagnificationGrid)
       MagnificationGrid(i) = MagFactorGridLower + (i-1)*((MagFactorGridHigher- MagFactorGridLower)/(size(MagnificationGrid)-1))
       Renormalisation_Size_Limits = Survey_Size_Limits/dsqrt(MagnificationGrid(i))
       Renormalisation_Magnitude_Limits = Survey_Magnitude_Limits + 2.5e0_double*dlog10(MagnificationGrid(i))

       Renormalisation_by_Magnification(i) = Integrate(PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, 2, lim1 = Renormalisation_Magnitude_Limits, lim2 = Renormalisation_Size_Limits)
       if(present(MagPrior)) then

          !--Edit this code to calculate the magnitude renormalisation over the whole size grid if Posterior_Method == 3, and over (0, Survey_Size_Limit(1)) if Posterior_Method == 4
          !-- vv This is only true if there are no cuts on the size-mag prior vv
!MagRenorm          MagOnly_Renormalisation_by_Magnification(i) =  Integrate(PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, 2, lim1 = Renormalisation_Magnitude_Limits, lim2 = Renormalisation_Size_Limits)
          MagOnly_Renormalisation_by_Magnification(i) = Integrate(PriorMagGrid, Survey_Renormalised_MagPrior, 2, lim = Renormalisation_Magnitude_Limits)
       end if
    end do
    deallocate(Size_Only_Prior)
       
    open(unit = 53, file = trim(adjustl(Output_Prefix))//'Renormalisation_by_Magnification.dat')
    do i = 1, size(Renormalisation_by_Magnification)
       write(53, *) MagnificationGrid(i), Renormalisation_by_Magnification(i)
    end do
    write(*,'(2(A))') 'File output: ', trim(adjustl(Output_Prefix))//'Renormalisation_by_Magnification.dat'

    !--Start of Posterior Routines--!
    allocate(Posterior_perGalaxy(size(Cat%RA), size(Posterior,2))); Posterior_perGalaxy = 1.e-100_double
    allocate(Effective_Magnification(size(Cat%RA), size(Posterior,2))); Effective_Magnification = 0.e0_double
    Convergence_Per_Cluster = 0.e0_double


    write(*,'(A)') '--------------------------------------------------------------------------------------------------'
    write(*,'(A)',advance = 'no') 'Getting Posterior for Cluster '
    if(Posterior_Method == 1) write(*,'(A)') 'using Sizes Only'
    if(Posterior_Method == 2) write(*,'(A)') 'using Sizes and Magnitudes'
    if(Posterior_Method == 3) write(*,'(A)') 'using Magnitudes Only'
    if(Posterior_Method == 4) write(*,'(A)') 'using Size and Magnitudes, and Magnitudes Only below the size limit'
    if(Enforce_Weak_Lensing) print *, '*** Weak lensing assumptions have been enforced'

    n_Default_Source_Redshift_Used = 0; nGal_Ignored_MagLimits = 0; nGal_Ignored_SizeLimits = 0; nGal_Ignored_NaN = 0
    nMagPosterior = 0; nSizePosterior = 0; nSizeMagPosterior = 0
    Time1 = 0.; Time2 = 0.
    do c = 1, size(Effective_Magnification,1) !-Loop over galaxies-!
       !--Select the method of posterior reconstruction for that point based in input method

!!$       Time1 = Time2
!!$       call CPU_TIME(Time2)
!!$
!!$       print *, 'Galaxy:', c, Time2-Time1

       select case(Posterior_Method)
       case(1) !--SizeOnly--!
          if( (Cat%Sizes(c) > Survey_Size_Limits(2)) .or. (Cat%Sizes(c) < Survey_Size_Limits(1))) then
             nGal_Ignored_SizeLimits = nGal_Ignored_SizeLimits + 1
             cycle
          end if
!MagRenorm          iSurvey_size_Limits = Survey_Size_Limits

          nSizePosterior = nSizePosterior + 1
          Galaxy_Posterior_Method = 1
       case(2)!--SizeMag--!
          if( (Cat%Sizes(c) > Survey_Size_Limits(2)) .or. (Cat%Sizes(c) < Survey_Size_Limits(1))) then
             nGal_Ignored_SizeLimits = nGal_Ignored_SizeLimits + 1
             cycle
          end if

!MagRenorm          iSurvey_size_Limits = Survey_Size_Limits

          nSizeMagPosterior = nSizeMagPosterior + 1
          Galaxy_Posterior_Method = 2
       case(3) !-Magnitude Only--!
          nMagPosterior = nMagPosterior + 1

!MagRenorm          iSurvey_size_Limits = (/0.e0_double, 1.e30_double/) !--Should encompase the whole data set

          if(present(MagPrior) == .false.) STOP 'Attempting to produce posterior using mag only, but no magnitde prior entered, stopping'
          Galaxy_Posterior_Method = 3
       case(4) !-Size mag above size data limit, Mag-Only below size data limit-!
          if(present(MagPrior) == .false.) STOP 'Attempting to produce posterior using mag only under size limit, but no magnitde prior entered, stopping'
          STOP 'Posterior Method 4 has been disabled as it needs further thought into its construction. In particular, with respect to the construction of the prior from the size-mag distriution, the effect of prior cuts, and correct renormalisation. The skeleton code has been added for this, but disabled for now. See commented code labeled MagRenorm'

          !--Note, in this case the data-renormalisation should take inot account that the mag only is constructed from a sample which takes small, faint galaxies
          if(Cat%Sizes(c) < Survey_Size_Limits(1)) then
!MagRenorm             iSurvey_Size_Limits = (/0.e0_double, Survey_Size_Limits(1)/)

             Galaxy_Posterior_Method = 3
             nMagPosterior = nMagPosterior + 1
          else
!MagRenorm             iSurvey_Size_Limits = Survey_Size_Limits

             Galaxy_Posterior_Method = 2
             nSizeMagPosterior = nSizeMagPosterior + 1
          end if
       end select

       if( (Cat%MF606W(c) > Survey_Magnitude_Limits(2)) .or. (Cat%MF606W(c) < Survey_Magnitude_Limits(1))) then
          nGal_Ignored_MagLimits = nGal_Ignored_MagLimits + 1
          cycle
       end if

       if(isNaN(Cat%Sizes(c)) .or. isNaN(Cat%MF606W(c))) then
          nGal_Ignored_NaN = nGal_Ignored_NaN + 1
          cycle
       end if

       if(c == size(Effective_Magnification,1)/2) print *, 'Halfway done for this Aperture..'

       Known_Redshift = .false.
       if(Cat%Redshift(c) >= 0.e0_double) then
          !--Ignore Galaxies with redshift less than the foreground-!
          if(Cat%Redshift(c) < Lens_Redshift) cycle
          !--Use Galaxy Redshift if available--!

          allocate(Posterior_perGalaxy_Redshift(size(Posterior_perGalaxy,1), size(Posterior_perGalaxy,2), 1)); Posterior_perGalaxy_Redshift = 0.e0_double
          allocate(Sigma_Crit(size(Cat%RA),1)); Sigma_Crit = -1.e0_double

          D_s = angular_diameter_distance_fromRedshift(0.e0_double, Cat%Redshift(c))
          D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, Cat%Redshift(c))

          Sigma_Crit(c,1) = 1.66492e0_double*(D_s/(D_l*D_ls))

          allocate(RedshiftPDF(1)); RedshiftPDF = -1.
          Known_Redshift = .true.
          
       elseif(MC_Sigma_Critical == .false. .and. Marginalise_Redshift_Distribution == .false.) then
          !--Use Default Redshift if everything else fails--!
          allocate(Posterior_perGalaxy_Redshift(size(Posterior_perGalaxy,1), size(Posterior_perGalaxy,2), 1)); Posterior_perGalaxy_Redshift = 0.e0_double
          allocate(Sigma_Crit(size(Cat%RA),1)); Sigma_Crit = -1.e0_double

          n_Default_Source_Redshift_Used = n_Default_Source_Redshift_Used + 1
          D_s = angular_diameter_distance_fromRedshift(0.e0_double, Default_Source_Redshift)
          D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, Default_Source_Redshift)
          
          Sigma_Crit(c,1) = 1.66492e0_double*(D_s/(D_l*D_ls)) !-(10^18 M_Sun/h)-!

          allocate(RedshiftPDF(1)); RedshiftPDF= -1.
       elseif(MC_Sigma_Critical .and. Marginalise_Redshift_Distribution == .false.) then
          !--Not really used anymore--!
          allocate(Posterior_perGalaxy_Redshift(size(Posterior_perGalaxy,1), size(Posterior_perGalaxy,2), 1)); Posterior_perGalaxy_Redshift = 0.e0_double
          allocate(Sigma_Crit(size(Cat%RA),1)); Sigma_Crit = -1.e0_double

          call Monte_Carlo_Redshift_Sampling_SigmaCritical(Cat, Lens_Redshift, Sigma_Crit_MC)
          Sigma_Crit(:,1) = Sigma_Crit_MC
          deallocate(Sigma_Crit_MC)

       elseif(Marginalise_Redshift_Distribution .and. MC_Sigma_Critical == .false.) then
          !--Marginalising over the redshift distribution for that galaxy--!
          allocate(Posterior_perGalaxy_Redshift(size(Posterior_perGalaxy,1), size(Posterior_perGalaxy,2), nRedshift_Sampling)); Posterior_perGalaxy_Redshift = 0.e0_double
          allocate(Sigma_Crit(size(Cat%RA),nRedshift_Sampling)); Sigma_Crit = -1.e0_double

          !--Get Sigma_Critical for each point on the Redshift PDF grid--!
          do z = 1, nRedshift_Sampling
             D_s = angular_diameter_distance_fromRedshift(0.e0_double, RedshiftGrid(z))
             D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, RedshiftGrid(z))
             Sigma_Crit(c,z) = 1.66492e0_double*(D_s/(D_l*D_ls)) !-(10^18 M_Sun/h)-!
          end do

          !--Get Redshift PDF for this galaxy--!
          allocate(RedshiftPDF(size(RedshiftGrid))); RedshiftPDF = 0.e0_double
       else
          STOP 'DM_Profile_Variable_Posterior - Both MC and Redshift Distribution methods set - this cannot be'
       end if
       if(any(Sigma_Crit(c,:) < 0.e0_double)) STOP 'DM_Profile_Variable_Posterior - Invalid Sigma Critical Entered, negative'

       Distance_from_Mass_Center = dsqrt( (Cat%RA(c)-Lens_Position(1))**2.e0_double + (Cat%Dec(c)-Lens_Position(2))**2.e0_double ) !-in Degrees-!
       Distance_from_Mass_Center = (D_l*Distance_from_Mass_Center*(3.142e0_double/(180.e0_double))) !-in Mpc/h-! 

!!$       if(c==1) write(*,'(A)', advance = 'no') 'Galaxy:'
!!$       write(*,'(I3,x,A)',advance = 'no') c, ','


       do i = 1, size(Effective_Magnification,2) !--Loop over posterior values--!

          do z = 1, size(RedshiftPDF)

             select case(Mass_Profile)
             case(1) !-Flat-!
                STOP 'The Strong Lensing approximation has not been implemented for this profile yet'
                Effective_Magnification(c,i) = 1.e0_double+2.e0_double*(Posterior(1,i)/Sigma_Crit(c,z))
             case(2) !-SIS-!
                STOP 'The Strong Lensing approximation has not been implemented for this profile yet'
                Effective_Magnification(c,i) = 1.e0_double+2.e0_double*(SMD_SIS(Posterior(1,i), Distance_From_Mass_Center)/(Sigma_Crit(c,z)*1.e18_double))
                Convergence_Per_Cluster(:,c)  = (/Distance_From_Mass_Center,  SMD_SIS(1073.e0_double*1073.e0_double, Distance_From_Mass_Center), Sigma_Crit(c,z)/) !---TESTING--!
             case(3) !-NFW-!
                if(Enforce_Weak_Lensing) then
                   Effective_Magnification(c,i) = 1.e0_double + 2.e0_double*(SMD_NFW(Distance_From_Mass_Center, Lens_Redshift, Posterior(1,i))/(Sigma_Crit(c,z)*1.e18_double))
                else
                   Effective_Magnification(c,i) = Magnification_Factor(3, Distance_From_Mass_Center, Posterior(1,i),  Lens_Redshift, Sigma_Crit(c,z)*1.e18_double)
                end if
             case default
                STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
             end select

             if(isNaN(Effective_Magnification(c,i))) then
                print *, 'Magnification Factor is a NaN:', Distance_From_Mass_Center, Posterior(1,i), Lens_Redshift, Sigma_Crit(c,z)*1.e18_double
                STOP
             end if


!             if(Effective_Magnification(c,i) < 1.e0_double) print *, '** DE-AMPLIFICATION!', Effective_Magnification(c,i), SMD_NFW_Scalar(Distance_From_Mass_Center, Lens_Redshift, Posterior(1,i)), Differential_SMD_Scalar(Distance_From_Mass_Center, Lens_Redshift, Posterior(1,i)), Sigma_Crit(c,z)*1.e18_double, Posterior(1,i), Distance_From_Mass_Center
             if(Effective_Magnification(c,i) < minval(MagnificationGrid) .or. Effective_Magnification(c,i) > maxval(MagnificationGrid)) then
                !--Skipping as outside limits on which magnification was evaluated--!
!                print *, 'Skipping since mu lie outside renormalisation limits (including < 1)', Effective_Magnification(c,i)
                Posterior_perGalaxy_Redshift(c,i,z) = 1.e-100_double
                cycle
             end if

             if(Effective_Magnification(c,i) <= 0.e0_double) then
!                print *, 'Skipping since Mu is negative'
!                Posterior_perGalaxy_Redshift(c,i,z) = 1.e-100_double
!                cycle
                !--This won't be picked up unless the above cycle is relaxed--!
                print *, Effective_Magnification(c,i), SMD_NFW_Scalar(Distance_From_Mass_Center, Lens_Redshift, Posterior(1,i)), Differential_SMD_Scalar(Distance_From_Mass_Center, Lens_Redshift, Posterior(1,i)), Sigma_Crit(c,z)*1.e18_double, Posterior(1,i), Distance_From_Mass_Center
                STOP 'Effective Magnification invalid - negative. Stopping'
             end if


             !--Construct Joint Size-Magnitde Prior which is renormalised according to the convergence value--!

             Renormalisation_Size_Limits = Survey_Size_Limits/dsqrt(Effective_Magnification(c,i))
             Renormalisation_Magnitude_Limits = Survey_Magnitude_Limits + 2.5e0_double*dlog10(Effective_Magnification(c,i))

             if(Cuts_Renormalise_Likelihood) then
                if(i == 1 .and. z ==1 .and. c == 1) print *, 'KAPPA RENORMALISATION TURNED ON'
                Renorm = Linear_Interp(Effective_Magnification(c,i), MagnificationGrid, Renormalisation_by_Magnification, ExValue = 1.e30_double)
                if(Renorm == 0) then
                   Kappa_Renormalised_Prior = 0.e0_double
                else
                   Kappa_Renormalised_Prior = Survey_Renormalised_Prior/Renorm
                end if
                Renorm = 0.e0_double
             else
                if(i == 1 .and. z ==1 .and. c == 1) print *, 'KAPPA RENORMALISATION TURNED OFF'
                Kappa_Renormalised_Prior =  Survey_Renormalised_Prior
                Kappa_Renormalised_MagPrior = Survey_Renormalised_MagPrior
             end if

             
             select case (Galaxy_Posterior_Method)
             case(1) !--Size Only--!

                !--Evaluate p_{theta_0, m_0}*p_{z|m_0} for the whole magnitude grid
                do m =  1, size(PriorMagGrid)
                   !--m_0, theta_0--!
                   if((PriorMagGrid(m) < Renormalisation_Magnitude_Limits(1)) .or. (PriorMagGrid(m) > Renormalisation_Magnitude_Limits(2))) then
                      !-Since Integrand does not extend over this region anyway-!
                      Size_Only_Mag_Prior(m,:) = 1.e-100_double
                   else
                      if(Known_Redshift) then
                         Size_Only_Mag_Prior(m,:) =  Kappa_Renormalised_Prior(m,:)
                      else
                         Size_Only_Mag_Prior(m,:) =  Kappa_Renormalised_Prior(m,:)*CH08_Redshift_Distribution_Scalar(PriorMagGrid(m), RedshiftGrid(z))
                      end if
                   end if
                end do


                if(lnSize_Prior) then
                   STOP 'ln Size prior has been switched off since the conversion to strong lensing - check the theory carefully before use'
                   if(dlog(Cat%Sizes(c)/dsqrt(Effective_Magnification(c,i))) > maxval(PriorSizeGrid) .or. (dlog(Cat%Sizes(c)/dsqrt(Effective_Magnification(c,i))) < minval(PriorSizeGrid))) then
                      !--Extrapolation--!
                      Posterior_perGalaxy_Redshift(c,i,z) = 0.e0_double
                      cycle
                   end if

                   !--Find Boundary Indexs for the Size at which the Prior will be evaluated - Aims to improve run-tim by minimising the number of integrations required
                   allocate(Size_Only_Prior(2)); Size_Only_Prior = 0.e0_double
                   do j = 1, size(PriorSizeGrid)-1
                      if( ( PriorSizeGrid(j)<= dlog(Cat%Sizes(c)/dsqrt(Effective_Magnification(c,i)))) .and. (  PriorSizeGrid(j+1) > dlog(Cat%Sizes(c)/dsqrt(Effective_Magnification(c,i))))) then
!                      if( (PriorSizeGrid(j)<= dlog(Cat%Sizes(c))-Effective_Convergence(c,i)) .and. ( PriorSizeGrid(j+1) > dlog(Cat%Sizes(c))-Effective_Convergence(c,i)) ) then
                         Size_Only_Prior(1) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
                         Size_Only_Prior(2) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j+1), 2, lim = Renormalisation_Magnitude_Limits)
                         exit
                      end if
                   end do

!!$                   Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(dlog(Cat%Sizes(c))-Effective_Convergence(c,i), PriorSizeGrid(j:j+1), Size_Only_Prior(:), ExValue = 0.e0_double) !-Reinstate with SLensing--!
                else
                   if((Cat%Sizes(c)/dsqrt(Effective_Magnification(c,i)) > maxval(PriorSizeGrid)) .or. (Cat%Sizes(c)/dsqrt(Effective_Magnification(c,i)) < minval(PriorSizeGrid))) then
!                   if(Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)) > maxval(PriorSizeGrid) .or. (Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)) < minval(PriorSizeGrid))) then !-WL-!
                      !--Extrapolation--!
                      Posterior_perGalaxy_Redshift(c,i,z) = 1.e-100_double
                      cycle
                   else
                      
                      !--Find Boundary Indexs for the Size at which the Prior will be evaluated - Aims to improve run-time by minimising the number of integrations required
                      allocate(Size_Only_Prior(2)); Size_Only_Prior = 0.e0_double
                      do j = 1, size(PriorSizeGrid)-1
                         if( (PriorSizeGrid(j)<= Cat%Sizes(c)/dsqrt(Effective_Magnification(c,i))) .and. ( PriorSizeGrid(j+1) > Cat%Sizes(c)/dsqrt(Effective_Magnification(c,i))) ) then
                            Size_Only_Prior(1) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
                            Size_Only_Prior(2) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j+1), 2, lim = Renormalisation_Magnitude_Limits)
                            exit
                         end if
                      end do

                      
                      Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(Cat%Sizes(c)/dsqrt(Effective_Magnification(c,i)), PriorSizeGrid(j:j+1), Size_Only_Prior(:), ExValue =  1.e-100_double)*(1.e0_double/dsqrt(Effective_Magnification(c,i)))

                   end if
                end if
                deallocate(Size_Only_Prior)
             case(2)!-Size and Magnitude-!

                if(Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification(c,i)) < 21.e0_double) then
                   print *, 'Possible problem with de-lensed magnitude - falls outwith bright limit of Scrabback Fit (21)'
                   print *, Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification(c,i))
                end if
                if(Known_Redshift) then
                   RedshiftPDF(z) = 1.e0_double
                else
                   RedshiftPDF(z) = CH08_Redshift_Distribution_Scalar(Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification(c,i)), RedshiftGrid(z))
                end if
                if(lnSize_Prior) then
                   STOP 'DM_Profile_Variable_Posterior - lnSize with Size_Magnitude Method - I cannae do that captain!'
                else
                   !--Uses distributions of the apparent size--!
                   !--Test for need for extrapolation---!
                   need_Extrapolate = .false.
                   if(( (Cat%MF606W(c)+2.5e0_double*dlog10(Effective_Magnification(c,i)) > maxval(PriorMagGrid)) .or. (Cat%MF606W(c)+2.5e0_double*dlog10(Effective_Magnification(c,i)) < minval(PriorMagGrid)) ) .or. ((Cat%Sizes(c)/dsqrt(Effective_Magnification(c,i)) > maxval(PriorSizeGrid)) .or. (Cat%Sizes(c)/dsqrt(Effective_Magnification(c,i)) > maxval(PriorSizeGrid)) )) then
                      need_Extrapolate = .true.
                   else
                      need_Extrapolate = .false.
                   end if
      
                   if((need_Extrapolate .and. do_KDE_Extrapolation) .or. do_KDE_OnTheFly) then !-.or. KDE_OnTheFly
                      !--KDE_Extrapolation / KDE_OnTheFly(? - What about entry of prior?)
                      Posterior_perGalaxy_Redshift(c,i,z) = KDE_BiVariate_Gaussian_Scalar(PriorCatalogue%MF606W, PriorCatalogue%Sizes, Cat%MF606W(c)+2.5e0_double*dlog10(Effective_Magnification(c,i)), Cat%Sizes(c)/dsqrt(Effective_Magnification(c,i)), Inverse_Covar = KDE_Covariance_Inverse, Det_Covar = KDE_Covariance_Determinant)*(1.e0_double/dsqrt(Effective_Magnification(c,i)))*RedshiftPDF(z)                         
                   elseif(need_Extrapolate .and. (do_KDE_Extrapolation == .false.)) then
                      !--Extrapolation, set to default, (zero)
                      Posterior_perGalaxy_Redshift(c,i,z) = 1.e-100_double
                   else
                      !--Interpolation--!
                      Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(Cat%MF606W(c)+2.5e0_double*dlog10(Effective_Magnification(c,i)), Cat%Sizes(c)/dsqrt(Effective_Magnification(c,i)), PriorMagGrid, PriorSizeGrid, Kappa_Renormalised_Prior, ExValue = 1.e-100_double)*(1.e0_double/dsqrt(Effective_Magnification(c,i)))*RedshiftPDF(z)
                   end if

                      !--If no KDE Extrapolation, then this will set to a default value (effectively zero) outside the prior grid range
                   end if
                case(3) !--Magnitude Only--!
                   !--Renoramlise Magnitude Prior Distribution--!

                   if(Cuts_Renormalise_Likelihood) then
                      Renorm = Linear_Interp(Effective_Magnification(c,i), MagnificationGrid, MagOnly_Renormalisation_by_Magnification, ExValue = 1.e30_double)                      
                      if(Renorm == 0) then
                         Kappa_Renormalised_MagPrior = 0.e0_double
                      else
                         Kappa_Renormalised_MagPrior = Survey_Renormalised_MagPrior/Renorm
                      end if
                   end if

                   if(Known_Redshift) then
                      RedshiftPDF(z) = 1.e0_double
                   else
                      RedshiftPDF(z) = CH08_Redshift_Distribution_Scalar(Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification(c,i)), RedshiftGrid(z))
                   end if
                   !--Could be edited for KDE Extrapolation, not done yet--!

                   !--Construct p_[m_0] as the magnitude distribution for a sample of galaxies between cuts-corrected survey size limits 
                   !-MagRenorm
!!$                   allocate(Mag_Only_Prior(2)); Mag_Only_Prior = 0.e0_double
!!$                   !--Find point where de-lensed magnitude lie on intrinsic distribution - could also use locate.
!!$                   do j = 1, size(PriorMagGrid)-1
!!$                      if( ( PriorMagGrid(j)<=  Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification(c,i)) ) .and. (  PriorMagGrid(j+1) > Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification(c,i))) ) then
!!$                         Mag_Only_Prior(1) = Integrate(PriorSizeGrid, Kappa_Renormalised_Prior(j,:), 2, lim = Renormalisation_Size_Limits)
!!$                         Mag_Only_Prior(2) = Integrate(PriorSizeGrid, Kappa_Renormalised_Prior(j+1,:), 2, lim = Renormalisation_Size_Limits)
!!$                         exit
!!$                      end if
!!$                   end do
!!$                   Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification(c,i)), PriorMagGrid(j:j+1), Mag_Only_Prior, ExValue =  1.e-100_double)*RedshiftPDF(z)
!!$                   deallocate(Mag_Only_Prior)
                   !---------------------------------------------------------------------------------------------------------------------

                   !--The following gives unbiased results if no size cuts are used, however is known to be biased in the presence of size cuts, I believe that this is due to the fact that the prior itself should depend on the size cuts in the presence of a size-magnitude correlation.
                   Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification(c,i)), PriorMagGrid, Kappa_Renormalised_MagPrior, ExValue =  1.e-100_double)*RedshiftPDF(z)
             case default
                STOP 'DM_Profile_Variable_Posterior - Error in choosing method of finding posterior'
             end select
          end do
       end do

       !--Get Posterior on Size for each galaxy by marginalising over the redshift distribution--!
       if(size(RedshiftPDF) == 1) then
          Posterior_perGalaxy(c,:) = Posterior_perGalaxy_Redshift(c,:,1)
       else
          !--Integrate over the redshift Information--!
          do i = 1, size(Posterior_perGalaxy,2) !--Loop over SMD alpha--!
             Posterior_perGalaxy(c,i) = TrapInt(RedshiftGrid, Posterior_perGalaxy_Redshift(c,i,:))
          end do
          if(all(Posterior_perGalaxy(c,:) == 0.e0_double)) then
             print *,  'Posterior per galxy is zero! for galaxy no:', c
             STOP
          END if
       end if
       deallocate(Posterior_perGalaxy_Redshift, RedshiftPDF, Sigma_Crit)

       IF(all(Posterior_perGalaxy(C,:) == 0.E0_DOUBLE)) then
          print *, 'Galxy ', c, ' has zero posterior'
          read(*,*)
       end IF

       !--Renormalise the posterior per galaxy - could reasonably be removed.
       Renorm = 0.e0_double
       do j = 1, size(Posterior_perGalaxy,2)-1
          Renorm = Renorm + 0.5e0_double*(Posterior_perGalaxy(c,j) + Posterior_perGalaxy(c,j+1))*(Posterior(1,j+1)-Posterior(1,j))
       end do
       if(Renorm == 0.e0_double) THEN
          PRINT *, 'DM_Profile_Variable_Posterior - Error Renormalising Galaxy posteriors - Posterior Empty', C
          !STOP
       END if

       Posterior_perGalaxy(c,:) = Posterior_perGalaxy(c,:)/Renorm
       Renorm = 0.e0_double
    end do !--End of galaxy loop-!
    if(n_Default_Source_Redshift_Used > 0) write(*,'(A,I3,A)') '****** Used the default redshift for:', n_Default_Source_Redshift_Used, ' galaxies ********'
    if(nGal_Ignored_SizeLimits > 0) write(*,'(A,I3)') '****** Number of galxies ignored as they fell outside the survey size limits:', nGal_Ignored_SizeLimits
    if(nGal_Ignored_NaN > 0) write(*,'(A,I3)') '****** Number of galaxies ignored as they were NaNs:', nGal_Ignored_NaN


    print *, '-------------------------------------------------------------'
    print *, 'Constructed ', nSizePosterior, ' size-only posteriors'
    print *, 'Constructed ', nSizeMagPosterior, ' size-magnitude posteriors'
    print *, 'Constructed ', nMagPosterior, ' magnitude-only posteriors'
    print *, '-------------------------------------------------------------'

    if(Debug_Mode .and. Output_Posterior_Per_Galaxy) then
!!$       Filename = trim(adjustl(Output_Prefix))//'Convergence_Per_Cluster.dat'
!!$       open(unit = 18, file = Filename)
!!$       write(fmtstring,'(I1)') 3
!!$       do i =1, size(Cat%RA)
!!$          write(18, '(I5, x, '//trim(adjustl(fmtstring))//'(e14.7,x))') Cat%Galaxy_Number(i), Convergence_Per_Cluster(:,i)
!!$       end do
!!$      close(18)
!!$      print *, 'Output Convergence per Cluster to: ', trim(adjustl(Filename))

       Filename = trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat'
       open(unit = 17, file = Filename)
       write(fmtstring,'(I7)') size(Posterior_perGalaxy,1)+1
       do i =1, size(Posterior,2)
          write(17, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posterior(1,i), Posterior_perGalaxy(:,i)
       end do
      close(17)
      print *, 'Output Posterior per galaxy to: ', trim(adjustl(Filename))
   end if

   call Combine_Posteriors(Posterior(1,:), Posterior_perGalaxy(:,:), Combine_log_Posteriors, Posterior(2,:))

    Filename= trim(adjustl(Output_Prefix))//'Posterior_Combined_Renormalised.dat'
    open(unit = 82, file = Filename)
    write(fmtstring,'(I1)') 2
    do i =1, size(Posterior,2)
       write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posterior(1,i), Posterior(2,i)
    end do
    close(82)
    print *, 'Output Combined Posterior to: ', trim(adjustl(Filename))

    if(any(isNAN(Posterior(2,:)))) then
       print *, 'Any NaNs in aperture posterior?', any(isNAN(Posterior(2,:))), count(isNAN(Posterior(2,:)) == .true.)
       print *, 'Stopping'
       STOP
    END if

    !--Remove any information less than a tolerance - Numerical Error--!
    where(Posterior(2,:) < 1.e-12_double)
       Posterior(2,:) = 0.e0_double
    end where

    deallocate(Effective_Magnification, Posterior_perGalaxy)
    if(allocated(Kappa_Renormalised_Prior)) deallocate(Kappa_Renormalised_Prior)
    if(allocated(Survey_Renormalised_Prior)) deallocate(Survey_Renormalised_Prior)
    if(allocated(Size_Only_Mag_Prior)) deallocate(Size_Only_Mag_Prior)

    !--On Successful Completion delete Poster per galaxy as large file--!
!!$    inquire(file = trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat', exist = here)
!!$    if(here == .false.) STOP "Posterior per galaxy doesn't exist, stopping before accidental deletion"
!!$    call system('rm '//trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat')

  end subroutine DM_Profile_Variable_Posterior


  subroutine return_Prior_Distributions(MagGrid, SizeGrid, Dist, MagDist, Dir, BFCat, do_KDE)
    use IO, only: readin; use Distributions, only: produce_Joint_Size_Magnitude_Distribution, produce_Magnitude_Distribution; use Integration; use Interpolaters, only: linear_Interp
    !--Returns the joint size magnitude distribution--!
    !-- If BFCat entered, then the Distribution is calculated from the catalogue, otherwise a read in is attempted --!
    !--Dist is Size-Magnitude Distribution. MagDist is magnitude-only distribution, which may differ from the projected Dist due to size cuts on the prior. Magnitude_Distribution is constructed from a catalogue which only cuts between entered prior magnitude limits--!
    real(double),intent(out),allocatable:: MagGrid(:), SizeGrid(:), Dist(:,:), MagDist(:)
    character(*),intent(in):: Dir
    type(Catalogue),intent(in),optional:: BFCat
    logical,optional:: do_KDE

    character(500):: Filename, Output_Filename, Input_Filename, MagFilename
    real(double),allocatable:: Input_Array(:,:)

    character(5)::fmtstring
    integer::i,j
    logical:: iDo_KDE

    real(double),allocatable:: tMagGrid(:), tMagDist(:)

    !--Cuts on the prior distribution--!
    real(double):: Redshift_Cuts(2)  = (/0.22e0_double, 100.e0_double/) !-- If no cuts, then lower should still be < -1 to ensure only galaxies without redshift information are cut
    type(Catalogue)::Cut_Catalogue

    !--TESTING--!
    real(double),allocatable::Size_Only(:), Mag_Only(:)
    logical:: here
    
    if(present(do_KDE)) then
       iDo_KDE = do_KDE
    else
       iDO_KDE = use_KDE_Smoothed_Distributions
    end if


 
    if(iDO_KDE) then
       Filename = 'MagSize_Distribution_KDE.dat'
       MagFilename = 'Magnitude_Distribution_KDE.dat'
    else
       Filename= 'MagSize_Distribution_Histogram.dat'
       MagFilename = 'Magnitude_Distribution_Histogram.dat'
    end if


    if(present(BFCat)) then
       print *, '--Producing Joint Size Magnitude Distribution from Catalogue.....'
       
!!$    print *, '**Applying Masks to Prior Catalogue:'                                                                                                                                                            
       Cut_Catalogue = BFCat

       print *, '------Applying Cuts on the Catalogue from which the prior is constructed---------'
!       call Mask_Circular_Aperture(Cut_Catalogue, Clusters_In%Position, (/2.e0_double, 2.e0_double, 2.e0_double, 2.e0_double/)/60.e0_double)
       call Cut_By_Magnitude(Cut_Catalogue, Prior_Magnitude_Limits(1), Prior_Magnitude_Limits(2))
       call Cut_By_PixelSize(Cut_Catalogue, Prior_Size_Limits(1), Prior_Size_Limits(2))
       call Cut_by_PhotometricRedshift(Cut_Catalogue, Redshift_Cuts(1), Redshift_Cuts(2))
       print *, '----- Finished Cuts on BF catalogue----------------------------------------------'

       call produce_Joint_Size_Magnitude_Distribution(SizeGrid, MagGrid, Dist, Cut_Catalogue, use_Physical_Sizes = Analyse_with_Physical_Sizes, Magnitude_Type = 2, Output_Dir = trim(Dir), ln_size_Distribution = use_lnSize_Prior, KDE_Smooth = ido_KDE)
       !--Output (matches the method of output used above)--!
       Output_Filename = trim(Dir)//trim(Filename)
       open(unit = 45, file = Output_Filename)
       
       write(fmtstring, '(I5)') size(MagGrid) + 1
       write(45, '('//trim(fmtstring)//'(e14.7,x))') 0.0e0_double, MagGrid
       do i= 1, size(SizeGrid)
          write(45, '('//trim(fmtstring)//'(e14.7,x))') SizeGrid(i), Dist(:,i)
       end do
       close(45)
       write(*,'(2(A))') 'Output to: ', trim(Output_Filename)
    
       call Catalogue_Destruct(Cut_Catalogue)

       print *, 'Producing magnitude Distribution using Integration Method'
       allocate(MagDist(size(MagGrid))); MagDist = 0.e0_double
       do i= 1, size(SizeGrid)
          MagDist(i) = Integrate(SizeGrid, Dist(i,:), 2, lim = (/minval(SizeGrid), maxval(SizeGrid)/))
       end do


!!$       print *, '--Producing Magnitude Distribution from Catalogue.....'
!!$       Cut_Catalogue = BFCat
!!$       
!!$       print *, '------Applying Cuts on the Catalogue from which the prior is constructed---------'
!!$       call Cut_By_Magnitude(Cut_Catalogue, Prior_Magnitude_Limits(1), Prior_Magnitude_Limits(2))
!!$       call Cut_by_PhotometricRedshift(Cut_Catalogue, Redshift_Cuts(1), Redshift_Cuts(2))
!!$       print *, '----- Finished Cuts on BF catalogue----------------------------------------------'
!!$
!!$       call produce_Magnitude_Distribution(MagGrid, MagDist, Cut_Catalogue, iDO_KDE)

       !--Output (matches the method of output used above)--!
       Output_Filename = trim(Dir)//trim(MagFilename)
       open(unit = 46, file = Output_Filename)
       do i= 1, size(magGrid)
          write(46, '(2(e14.7,x))') magGrid(i), MagDist(i)
       end do
       close(46)
       write(*,'(2(A))') 'Output to: ', trim(Output_Filename)
       

    else
       Input_Filename = trim(Dir)//trim(Filename)

       write(*,'(2A)') 'Reading In Size-Magnitude Distribution from:', trim(Input_Filename)

       inquire(file = Input_Filename, exist = here)
       if(here == .false.) STOP 'return_Prior_Distributions - Size-Magnitude Distribution File does not exist, stopping...'

       call ReadIn(Input_Array, filename  = trim(adjustl(Input_Filename)), tabbed = .false., header_label = '#')
       !--Output of ReadIn is (Col, Row)

       print *, 'Input_Array has:', size(Input_Array,1), size(Input_Array,2)

       allocate(MagGrid(size(Input_Array,1)-1)); MagGrid = 0.e0_double
       allocate(SizeGrid(size(Input_Array,2)-1)); SizeGrid = 0.e0_double
       allocate(Dist(size(MagGrid), size(SizeGrid))); Dist = 0.e0_double

       MagGrid = Input_Array(2:,1)
       SizeGrid = Input_Array(1,2:)
       Dist = Input_Array(2:,2:)

       deallocate(Input_Array)

       !--Magnitude Distribution-!
       print *, 'Producing magnitude Distribution using Integration Method'
       allocate(MagDist(size(MagGrid))); MagDist = 0.e0_double
       do i= 1, size(SizeGrid)
          MagDist(i) = Integrate(SizeGrid, Dist(i,:), 2, lim = (/minval(SizeGrid), maxval(SizeGrid)/))
       end do

!!$
!!$       Input_Filename = trim(Dir)//trim(MagFilename)
!!$
!!$       write(*,'(2A)') 'Reading In Magnitude Distribution from:', trim(Input_Filename)
!!$       inquire(file = Input_Filename, exist = here)
!!$       if(here == .false.) STOP 'return_Prior_Distributions - Magnitude Distribution File does not exist, stopping...'
!!$
!!$       call ReadIn(Input_Array, filename  = trim(adjustl(Input_Filename)), tabbed = .false., header_label = '#')
!!$       !--Output of ReadIn is (Col, Row)  
!!$
!!$       print *, 'Input_Array has:', size(Input_Array,1), size(Input_Array,2)
!!$
!!$       allocate(tMagGrid(size(Input_Array,2))); tMagGrid = 0.e0_double
!!$       allocate(tMagDist(size(tMagGrid))); MagDist = 0.e0_double
!!$
!!$       tMagGrid = Input_Array(1,:)
!!$       tMagDist = Input_Array(2,:)
!!$
!!$       deallocate(Input_Array)
!!$       !--Linearly Interpolate onto Joint Size-Magnitude Distribution (this shouldn't be a problem since the output is the same as the input)--!
!!$       allocate(MagDist(size(MagGrid))); MagDist = 0.e0_double
!!$       call Linear_Interp(MagGrid, MagDist, tMagGrid, tMagDist)

       !--Test output--!
       Output_Filename = trim(Dir)//'Magnitude_Distribution_Prior_ReadinTest_DELETE.dat'
       open(unit = 46, file = Output_Filename)
       do i= 1, size(magGrid)
          write(46, '(2(e14.7,x))') MagGrid(i), MagDist(i)
       end do
       close(46)
       write(*,'(2(A))') 'Output to: ', trim(Output_Filename)
       
    end if

!    print *, 'Testing of Distribution Input:'
    allocate(Size_Only(size(SizeGrid))); Size_Only = 0.e0_double
    allocate(Mag_Only(size(MagGrid))); Mag_Only = 0.e0_double
    open(unit = 23, file = trim(Dir)//'SizeOnlyDist_Test_DELETE.dat')
    open(unit =24, file = trim(Dir)//'MagOnlyDist_Test_DELETE.dat')
    do i= 1, size(Size_Only)
       Size_Only(i) = Integrate(MagGrid, Dist(:,i),2, lim = (/minval(MagGrid), maxval(MagGrid)/))
       write(23,*) SizeGrid(i), Size_Only(i)
    end do
    do i= 1, size(Mag_Only)
       Mag_Only(i) = Integrate(SizeGrid, Dist(i,:),2, lim = (/minval(SizeGrid), maxval(SizeGrid)/))
       write(24,*) MagGrid(i),Mag_Only(i)
    end do
    close(23)
    close(24)
!!$    print *, 'Output to ', trim(Dir), 'SizeOnlyDist_Test_DELETE.dat, MagOnlyDist_Test_DELETE.dat'
!!$    read(*,*)

  end subroutine return_Prior_Distributions


  integer function Nearest_Neighbour_Index(Grid, Val)
    !--Returns the index of the point of the grid which is closest to the value - this requires the grid to be well sampled WHICH CAN'T BE TESTED--!
    !--Assumes ORDERED GRID--!
    real(double), intent(in):: Grid(:), Val

    integer:: i, index
    real(double)::Dist(2)

    !--Find indexes between which the value lies
    do i =1, size(Grid)
       if( (Grid(i) <= Val) .and. (Grid(i+1) > Val) ) then
          index = i
          exit
       end if
    end do

    !--Calculate Dist
    Dist = (/Val-Grid(index),Grid(index+1)-Val/)

    if(Dist(1) < Dist(2)) then
       Nearest_Neighbour_Index = index
    else
       Nearest_Neighbour_Index = index + 1
    end if

  end function Nearest_Neighbour_Index
    

  real(double) function Linear_Interpolation(Ref_Grid, Ref_Value, Grid_Wanted)
    real(double), intent(in)::Ref_Grid(:), Ref_Value(:)
    real(double),intent(in)::Grid_Wanted

    integer::i

    if(Grid_Wanted > maxval(Ref_Grid)) then
!       print *, 'Linear_Interpolation - Fatal Error - Grid value wanted is larger than refernce grid:', Grid_Wanted, maxval(Ref_Grid)
       Linear_Interpolation = 0.e0_double
!       STOP
    end if
    if(Grid_Wanted < minval(Ref_Grid)) then
!       print *,'Linear_Interpolation - Fatal Error - Grid value wanted is smaller than refernce grid, extrapolation needed:',Grid_Wanted, minval(Ref_Grid)
       Linear_Interpolation = 0.e0_double
!       STOP
    end if

    do i =1 , size(Ref_Grid)-1
       if(Ref_Grid(i) == Grid_Wanted) then
          Linear_Interpolation = Ref_Value(i)
          exit
       elseif(Ref_Grid(i+1) == Grid_Wanted) then
          Linear_Interpolation = Ref_Value(i+1)
          exit
       elseif( (Ref_Grid(i) < Grid_Wanted) .and. (Ref_Grid(i+1) > Grid_Wanted) ) then 
          Linear_Interpolation = Ref_Value(i) + ((Ref_Value(i+1)-Ref_Value(i))/(Ref_Grid(i+1)-Ref_Grid(i)))*(Grid_Wanted-Ref_Grid(i))
          exit
       end if
    end do

  end function Linear_Interpolation

end module Bayesian_Routines
