module  Bayesian_Routines
  !--Uses Bayesian Means to Calculate Posterior for the Dark Matter Profile variables--!
  !--Method must have redshifts for every galaxy (20Jan2014)--!
  use Param_Types; use Catalogues
  implicit none

  character(200):: Bayesian_Routines_Output_Directory
  
  real(double)::Default_Source_Redshift = 1.4e0_double
  logical:: Analyse_with_Physical_Sizes = .false. !#Default Value#

  logical::Combine_log_Posteriors = .false. !-If False, then combined by multiplication-!

  logical,private::Debug_Mode = .true.

  integer::Surface_Mass_Profile = 3 !-1:Flat, 2:SIS, 3:NFW-!

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
       Combined_Posterior = 1.0e-100_double
       Combined_Posterior = log(Combined_Posterior)
       do c = 1, size(Posteriors,1) !-Loop over galaxies-!
          Combination_Normalisation = dlog(maxval(Posteriors(c,:)))!1.e0_double/size(Posteriors,1)
          if(all(Posteriors(c,:) == 0.e0_double) .or. all(isNaN(Posteriors(c,:)))) then
             nPosteriorsSkipped = nPosteriorsSkipped + 1
             cycle
          end if
!          where(Posteriors(c,:) /= 0.e0_double)
             Combined_Posterior = Combined_Posterior + dlog(Posteriors(c,:)) - Combination_Normalisation
!          end where
       end do
       where(Combined_Posterior /= 0.e0_double)
          Combined_Posterior = dexp(Combined_Posterior)
       end where
    end if
 
    if(any(isNAN(Combined_Posterior(:)))) then
       print *, 'Any NaNs in Combined Posterior?', any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
       print *, 'Stopping'
       STOP
    END if

    if(nPosteriorsSkipped > 0) write(*,'(A,I4,A,I4,A)') '################### Combine_Posteriors - ', nPosteriorsSkipped, ' of ', size(Posteriors,1), ' posteriors were skipped as they we invlaid - NaNs or 0 ###########################'
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
    if(present(Error)) Error = dsqrt(variance_distribution(PosteriorGrid, Posterior))
    if(present(AntiSymm_Error)) then
       print *, 'Getting AntiSymmetric Error:'
       if(present(ModeVal)) then
          !--Take error about mode--!
          AntiSymm_Error = Antisymmetric_Variance_Distribution(PosteriorGrid, Posterior, ModeVal)
       else
          AntiSymm_Error = Antisymmetric_Variance_Distribution(PosteriorGrid, Posterior)
       end if
       print *, 'Got Antisymmetric Error'
    end if
    if(present(MeanVal)) MeanVal = mean(Posterior, PosteriorGrid)

  end subroutine Posterior_Statistics

  subroutine Identify_Galaxys_in_Circular_Aperture(Cat, Ap_Pos, Ap_Radius, Ap_Cat)
    !--Returns a reduced catalogue containing only the galaxies in the aperture
    !--Ap_Pos and Ap_Radius should be in DEGREES
    type(Catalogue), intent(in)::Cat
    real(double),intent(in)::Ap_Pos(:), Ap_Radius
    type(Catalogue), intent(out)::Ap_Cat

    integer::c, ac
    integer:: Expected_Number_In_Aperture

    Expected_Number_in_Aperture = count( dsqrt( (Cat%RA-Ap_Pos(1))**2.e0_double +(Cat%Dec-Ap_Pos(2))**2.e0_double ) <= Ap_Radius)

    if(Expected_Number_in_Aperture == 0) STOP 'Identify_Galaxys_in_Circular_Aperture - There are no galaxies within the aperture, suggest increasing aperture size'

    call Catalogue_Construct(Ap_Cat, Expected_Number_in_Aperture)

    ac = 0
    do c = 1, size(Cat%RA)
       if(dsqrt( (Cat%RA(c)-Ap_Pos(1))**2.e0_double +(Cat%Dec(c)-Ap_Pos(2))**2.e0_double ) <= Ap_Radius ) then
          if(ac > Expected_Number_in_Aperture) STOP 'Identify_Galaxys_in_Circular_Aperture - Error in assigning aperture galaxy - GALAXY ASSINGATION IS LARGER THAN EXPECTED, stopping..'
          ac = ac + 1
          call Catalogue_Assign_byGalaxy_byCatalogue(Ap_Cat, ac, Cat, c)
       end if
    end do

  end subroutine Identify_Galaxys_in_Circular_Aperture

  subroutine DM_Profile_Variable_Posteriors_CircularAperture(Cat, Ap_Pos, Ap_Radius, Posteriors, Blank_Field_Catalogue)
    use Statistics, only: mean, mode_distribution, variance_distribution; use Mass_profiles
    use Distributions; use Cosmology, only: angular_diameter_distance_fromRedshift;
    !--Main routine that returns the Posteriors over all apertures.
    !--Ap_Radius in DEGREES
    !--To Do: 
    !~~Recent edits to use a joint size magnitude have ignored magnitude binning - either remove this option, or edit to allow


    type(Catalogue), intent(in)::Cat
    real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:)
    real(double),intent(out),allocatable::Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-! - Second dimension allows for a different grid for each aperture
    type(Catalogue), intent(in),optional::Blank_Field_Catalogue

    real(double),dimension(Size(Ap_Pos,1))::iAp_Radius

    type(catalogue),dimension(size(Ap_pos,1))::Ap_Cats
    integer::i, ap, m, j

    !--Size Grid Declarations--!     
    integer::nSizeGrid = 100        
!    real(double)::SizeGrid_Lower = 0.e0_double, SizeGrid_Upper = 70.e0_double, dSizeGrid !-Pixel-! 
    real(double)::SizeGrid_Lower = 0.e0_double, SizeGrid_Upper = 0.0035e0_double, dSizeGrid
    real(double),allocatable::SizeGrid(:)

    real(double),allocatable:: Joint_Size_Magnitude_Distribution(:,:)
    real(double),allocatable:: Size_Only_Distribution(:)

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

    character(200):: Output_File_Prefix

    logical::use_lnSize_Prior = .false.
    logical:: use_KDE_Smoothed_Distributions = .true.
    real(double),allocatable::MagGrid(:) !--Discardable for now as marginalised over--!

    INTERFACE
         subroutine DM_Profile_Variable_Posteriors_CircularAperture(Cat, Ap_Pos, Ap_Radius, Posteriors, Blank_Field_Catalogue)
           use Param_Types; use Catalogues
           !--Main routine that returns the Posteriors over all apertures.
           !--Ap_Radius in DEGREES
           type(Catalogue), intent(in)::Cat
           real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:)
           real(double),intent(out),allocatable::Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-! - Second dimension allows for a different grid for each aperture
           type(Catalogue), intent(in),optional::Blank_Field_Catalogue
         end subroutine DM_Profile_Variable_Posteriors_CircularAperture
      end INTERFACE

    if(size(Ap_Radius)==1) then
       iAp_Radius = Ap_Radius(1)
    else
       iAp_Radius = Ap_Radius
    end if

    !--Identify Reduced Catalogue for each aperture--!
    do i =1, size(Ap_Cats)
       call Identify_Galaxys_in_Circular_Aperture(Cat, Ap_Pos(i,:), iAp_Radius(i), Ap_Cats(i))
    end do

    !--Get Prior Distributions--!
    print *, count(Cat%Redshift >= 0.e0_double), ' of ',size(Cat%Redshift), ' galaxies have redshift information'

    !--Not used anymore---!
    if(Magnitude_Binning_Type == 1) then
       call Calculate_Bin_Limits_by_equalNumber(Cat%Absolute_Magnitude, nMag, MagBins)
    elseif(Magnitude_Binning_Type == 2) then
       call Calculate_Bin_Limits_by_equalNumber(Cat%MF606W, nMag, MagBins)
    else
       STOP 'DM_Profile_Variable_Posteriors_CircularAperture - Invalid Magnitude Type Specified stopping...'
    END if
    !^^^^^^^^^^^^^^^^^^^!

    !--If a catalogue for the blank field is passed in, then use this catalogue to get the intrinsic distributions--!
    if(present(Blank_Field_Catalogue)) then
       call get_Joint_Size_Magnitude_Distribution(SizeGrid, MagGrid, Joint_Size_Magnitude_Distribution, Blank_Field_Catalogue, use_Physical_sizes = Analyse_with_Physical_Sizes, Magnitude_Type = 2, ln_size_Distribution = use_lnSize_Prior, KDE_Smooth = use_KDE_Smoothed_Distributions)
    else
       call get_Joint_Size_Magnitude_Distribution(SizeGrid, MagGrid, Joint_Size_Magnitude_Distribution, Cat, use_Physical_sizes = Analyse_with_Physical_Sizes, Magnitude_Type = 2, ln_size_Distribution = use_lnSize_Prior, KDE_Smooth = use_KDE_Smoothed_Distributions)
    end if

!!$    Size_Only_Distribution =0.e0_double
!!$    do j = 1, size(SizeGrid)
!!$       Size_Only_Distribution(j) = Integrate(MagGrid, Joint_Size_Magnitude_Distribution(:,j), 2)
!!$    end do
!!$    allocate(SizePrior_byMag(1, size(SizeGrid))); SizePrior_byMag(1,:) = Size_Only_Distribution !-Should be correctly renormalised if p(theta,m) is

    if(size(magBins,1) > 1) STOP 'I have had to disable Magnitude Binning for now, youll have to edit the code to get this to work, stopping'

    !--Get Prior by getting refernce for each bin, and binning each reduced cat using the same definition (mag)--!
    !--Produce Posterior for each mag bin--!

    do ap = 1, size(Ap_Cats)
       call bin_catalogue_by_magnitude(Ap_Cats(Ap), MagBins, BCat)

       print *, 'Aperture Catalogue binned in the same way with magnitudes: Aperture: ', Ap, ' has ', size(Ap_Cats(Ap)%RA), ' galaxies binned accordjing to:', BCat%Occupation

       do m = 1, size(MagBins,1)
          !--Set up Prior_Single for that Aperture and Mag Bin-!
!DELETE (or Edit to Joint)          allocate(Prior_Single(2, size(SizePrior_byMag,2))); Prior_Single(1,:) = SizeGrid; Prior_Single(2,:) = SizePrior_byMag(m,:)

          !-Get Posterior for that Magnitude-!
          write(Output_File_Prefix,'(I2)') Ap
          Output_File_Prefix = trim(adjustl(Bayesian_Routines_Output_Directory))//'Aperture_'//trim(adjustl(Output_File_Prefix))//'_'
          call DM_Profile_Variable_Posterior(BCat%Cat(m), Surface_Mass_Profile, MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Lens_Redshift, Ap_Pos(Ap,:), Posterior_Single, Output_File_Prefix, use_lnSize_Prior)


          if(any(isNAN(Posterior_Single(2,:)))) then
             print *, 'Any NaNs in returned posterior?', any(isNAN(Posterior_Single(2,:))), count(isNAN(Posterior_Single(2,:)) == .true.)
             print *, 'Stopping'
             STOP
          END if


          if(m==1) then
             allocate(Aperture_Posterior_byMag(size(MagBins,1), 2, size(Posterior_Single,2))); Aperture_Posterior_byMag = 0.e0_double !-Assumes posterior on the same grid as the prior
          end if
          Aperture_Posterior_byMag(m,1,:) = Posterior_Single(1,:); Aperture_Posterior_byMag(m,2,:) = Posterior_Single(2,:)
          

          deallocate(Posterior_Single)
       end do
       call Binned_Catalogue_Destruct(BCat)

       !--Convert Posteriors in NFW from Virial Radius to VirialMass--!
!!$       if(Surface_Mass_Profile == 3) then
!!$          print *, '**Converting NFW posteriors from Virial radius to Virial Mass:....'
!!$          do m =1, size(Aperture_Posterior_byMag,1)
!!$             Aperture_Posterior_byMag(m,1,:) = get_NFW_VirialMass_from_VirialRadius(Lens_Redshift, Aperture_Posterior_byMag(m,1,:))
!!$          end do
!!$       end if
 
       !--Recombine Posterior for that Aperture, over all mag bins - Assumes all on the same grid-!
       if(Ap==1) then
          allocate(Posteriors(size(Ap_Cats),2, size(Aperture_Posterior_byMag,3))); Posteriors = 1.e0_double
       end if
       Posteriors(Ap,1,:) = Aperture_Posterior_byMag(1,1,:);

       !--Combine as posteriors--!
       call  Combine_Posteriors(Posteriors(Ap,1,:), Aperture_Posterior_byMag(:,2,:), Combine_log_Posteriors, Posteriors(Ap,2,:))

       !--Output Posterior per Mag Bin--!
       write(apString, '(I1)') Ap
       open(38, file =trim(Bayesian_Routines_Output_Directory)//'Posterior_Aperture_'//trim(apString)//'_byMagBin.dat')
       !--Header--!!!!!!!!!!!!!!!!!
       do j = 1, size(MagBins,1)
          write(38, '(A1, 2(e14.7,x))') '#', MagBins(j,:)
       end do
       write(fmtstring,'(I2)') size(Aperture_Posterior_byMag,1)+1 !--Assumes all on the same grid--!
       do j =1, size(Aperture_Posterior_byMag,3)
          write(38, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Aperture_Posterior_byMag(1,1,j), Aperture_Posterior_byMag(:,2,j)
       end do
       close(38)
       print *, 'Output file to: ',trim(Bayesian_Routines_Output_Directory)//'Posterior_Aperture_'//trim(apString)//'_byMagBin.dat' 

       deallocate(Aperture_Posterior_byMag)
    end do


    !--Output Posterior--!
    open(unit = 51, file = trim(Bayesian_Routines_Output_Directory)//'Posterior_per_Aperture.dat')
    write(fmtstring,'(I2)') size(Posteriors,1)+1 !-Assumes all posteriors described by the same posterior grid-!
    do j = 1, size(Posteriors,3)
       write(51, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posteriors(1,1,j), Posteriors(:,2,j)
    end do
    close(51)
    print *,'Output file to: ', trim(Bayesian_Routines_Output_Directory)//'Posterior_per_Aperture.dat'

    !--Convert to a Mass as A*Sigma(mean)--!
    allocate(Cluster_Mean(size(Posteriors,1))); Cluster_Mean = 0.e0_double
    allocate(Cluster_Mode(size(Posteriors,1))); Cluster_Mode = 0.e0_double
    allocate(Cluster_Variance(size(Posteriors,1))); Cluster_Variance = 0.e0_double
    allocate(AntiSymm_Variance(size(Posteriors,1),2)); AntiSymm_Variance = 0.e0_double

    open(51, file = trim(Bayesian_Routines_Output_Directory)//'Mass_Estimates.dat')
    write(51, '(A)') '# Aperture, Mode, Mode_Error'

    print *, '!----------------------------------------------------------------------------------!'
    do Ap = 1, size(Posteriors,1)
       !--Get Mean, Mode, Variance--!
       call Posterior_Statistics(Posteriors(Ap,1,:), Posteriors(Ap,2,:), Cluster_Mean(Ap), Cluster_Mode(Ap), Cluster_Variance(Ap), AntiSymm_Variance(Ap,:))
       if(Surface_Mass_Profile == 1) then
          !--Convert into the correct units (i.e. Msun/h)--!
          Cluster_Mean(Ap) = 1.e18_double*Cluster_Mean(Ap); Cluster_Mode(Ap) =  1.e18_double*Cluster_Mode(Ap); Cluster_Variance(Ap) = 1.e18_double*Cluster_Variance(Ap); AntiSymm_Variance(Ap,:) = 1.e18_double*AntiSymm_Variance(Ap,:)
       end if

       D_l =  angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
       Area = 3.142e0_double*(D_l*((3.142e0_double*Ap_Radius(Ap))/180.e0_double))**2.e0_double
!       print *, 'Typical Conversion from Sigma to Mass (x10^18 Msun/h):', Area
       write(Error_Mass_String_Positive, '(e8.2)') AntiSymm_Variance(Ap,2); write(Error_Mass_String_Negative, '(e8.2)') AntiSymm_Variance(Ap,1)
       write(Mass_String, '(e8.2)') Cluster_Mean(Ap)
       if(Surface_Mass_Profile == 1) write(*,'(A,I2,6(A))') 'Mean Sigma_0 of Cluster :', Ap, ', is (Mean): ', trim(Mass_String), '   +', trim(Error_Mass_String_Positive), ' -', Error_Mass_String_Negative
       if(Surface_Mass_Profile == 2) write(*,'(A,I2,6(A))') 'Mean Velocity Dispersion of Cluster :', Ap, ', is (Mean): ', trim(Mass_String), '   +', trim(Error_Mass_String_Positive), ' -', Error_Mass_String_Negative
       if(Surface_Mass_Profile == 3) write(*,'(A,I2,6(A))') 'Mean Virial Radius of Cluster :', Ap, ', is (Mean): ', trim(Mass_String), '   +', trim(Error_Mass_String_Positive), ' -', Error_Mass_String_Negative
       
       write(Mass_String, '(e8.2)') Cluster_Mode(Ap)
       write(*,'(6(A))') '                                   and (Mode): ', trim(Mass_String), '   +', trim(Error_Mass_String_Positive), ' -', Error_Mass_String_Negative

       !--Get Mean Mass--!
       call  Integrated_Mass_Within_Radius(Surface_Mass_Profile, D_l*iAp_Radius(Ap)*(3.142e0_double/180.e0_double), Cluster_Mean(Ap), AntiSymm_Variance(Ap,:), Scale_Mass, Scale_Mass_Error, Redshift = Lens_Redshift)
       write(Mass_String, '(e8.2)') Scale_Mass; write(Error_Mass_String_Positive, '(e8.2)') Scale_Mass_Error(2); write(Error_Mass_String_Negative, '(e8.2)') Scale_Mass_Error(1)
       write(*,'(6A)') "Mass of Cluster within 1' is (Mean):", Mass_String, '   +', trim(Error_Mass_String_Positive), ' -', Error_Mass_String_Negative
       !--Get Mode Mass--!
       call  Integrated_Mass_Within_Radius(Surface_Mass_Profile, D_l*iAp_Radius(Ap)*(3.142e0_double/180.e0_double), Cluster_Mode(Ap), AntiSymm_Variance(Ap,:), Scale_Mass, Scale_Mass_Error, Redshift = Lens_Redshift)
       write(Mass_String, '(e8.2)') Scale_Mass; write(Error_Mass_String_Positive, '(e8.2)') Scale_Mass_Error(2); write(Error_Mass_String_Negative, '(e8.2)') Scale_Mass_Error(1)
       write(*,'(6A)') "Mass of Cluster within 1' is (Mode):", Mass_String, '   +', trim(Error_Mass_String_Positive), ' -', Error_Mass_String_Negative
       print *, ' '

       write(51, '(I2,3(e9.3,x))') Ap, Scale_Mass,  Scale_Mass_Error
    end do
    print *, '!----------------------------------------------------------------------------------!'

  end subroutine DM_Profile_Variable_Posteriors_CircularAperture



  subroutine DM_Profile_Variable_Posterior(Cat, Mass_Profile, PriorMagGrid, PriorSizeGrid, Prior, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior)
    use cosmology, only:angular_diameter_distance_fromRedshift; use MC_Redshift_Sampling, only: Monte_Carlo_Redshift_Sampling_SigmaCritical; use Mass_Profiles, only:SMD_SIS, SMD_NFW; use Distributions, only: ch08_redshift_distribution; use Interpolaters, only: Linear_Interp
    use Integration, only:TrapInt, Integrate
    !--Returns the Bayesian posterior for the DM Profile Variables, as defined by Mass_Profile
    !-For Mass_Profile = 1 (Flat): Sigma_0
    !-                   2 (SIS) : Velocity_Dispersion**2
    !-                   3 (NFW) : Virial Radius (r200)
    !--Prior is the prior distribution of p(m,R), or possibly p(m,lnR), m is *apparent magnitude*, R is *apparent size*
    !--Cat entered already assumed to be only the galaxies that are to be considered (e.g. in aperture), and well-described by posterior (e.g. same mag bin) 

    !---TO DO:
    !--Conversion from posterior to surface mass density (possibly not here, but later)

    type(Catalogue)::Cat
    integer,intent(in)::Mass_Profile !-1:Flat, 2:SIS-!
    real(double),intent(in)::Lens_Redshift
    real(double),intent(in)::Lens_Position(2)
    real(double),intent(in):: PriorSizeGrid(:), PriorMagGrid(:), Prior(:,:) !-Magnitude, Size-!
    real(double),allocatable,intent(Out)::Posterior(:,:) !-Grid/Posterior, Value-!
    character(*), intent(in)::Output_Prefix
    logical,intent(in)::lnSize_Prior

    !--Method: 1: Size-Only, 2: Size-Magnitude--!
    integer:: Posterior_Method  = 1 !--Eventually pass in--!

    integer::i, c, j, z

    !--Variable Grid Declarations-!
    integer::nGrid = 1000
    real(double)::VGrid_Lower, VGrid_Higher
    
    real(double),dimension(size(Prior,2)):: Size_Only_Prior
    real(double),dimension(size(Prior,1)):: Mag_Only_Prior
    real(double),dimension(size(Prior,1),size(Prior,2)):: Size_Given_Mag_Prior
    real(double),allocatable::Posterior_perGalaxy(:,:) !-Galaxy, Posterior-! - Uses Same Grid as overall Posterior -
    real(double),allocatable:: Effective_Convergence(:,:) !-Galaxy, Posterior-! -Same Size as Posterior Grid, may need reevaluated per galaxy

    logical:: MC_Sigma_Critical = .false.
    logical:: Marginalise_Redshift_Distribution = .true.

    real(double),allocatable::Sigma_Crit(:,:), Sigma_Crit_MC(:)
    real(double)::D_l, D_s, D_ls
    real(double)::Distance_from_Mass_Center

    real(double)::Renorm

    logical:: Output_Posterior_Per_Galaxy = .true.
    character(7)::fmtstring
    character(200)::Filename
    logical::here

    !--Redshift Distribution Declarations--!
    integer,parameter:: nRedshift_Sampling = 50
    real(double),parameter::Redshift_Lower = 0.156e0_double, Redshift_Higher = 3.e0_double !!!Edit to Lens_Redshift
    real(double), dimension(nRedshift_Sampling):: RedshiftGrid
    real(double),allocatable::RedshiftPDF(:)
    real(double),allocatable:: Posterior_perGalaxy_Redshift(:,:,:) !-Galaxy, Posterior, Redshift-!

    !--Kappa dependant renormalisation--!
    real(double),dimension(2):: Survey_Magnitude_Limits, Survey_Size_Limits !--User defined needs edit to below, these are set by default from distribution-!
    real(double),dimension(2):: Renormalisation_Magnitude_Limits, Renormalisation_Size_Limits
    real(double),allocatable:: ConvergenceGrid(:), Renormalisation_by_Convergence(:) 
    !vv Must be set here vv!
    integer:: nConvergenceGrid = 1000
    real(double):: KappaGridLower = 0.e0_double, KappaGridHigher = 4.0e0_double

    !--Testing Declarations--!
    real(double),dimension(3,size(Cat%RA)):: Convergence_per_Cluster
    integer:: n_Default_Source_Redshift_Used

    if(Analyse_with_Physical_Sizes) STOP 'DM_Profile_Variable_Posterior - I HAVE DISABLED THE ABILITY TO USE PHYISCAL SIZES AS UNNECESSARY, code still to be edited'


    !-Set Up Posterior Grid-!
    select case(Mass_Profile)
    case(1) !-Flat-!
       nGrid = 100000
       VGrid_Lower = -5.e-3_double; VGrid_Higher = 1.e-2_double !-- SMD ~ Masses 10^12 -> 10^15 Msun/h, in Units of 10^18 Msun/h  
    case(2) !-SIS-!
       nGrid = 10000 !--This needs to be sigma_v^2--!
       VGrid_Lower= -2.5e5_double; VGrid_Higher = 9.e6_double    !--------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case(3)
       nGrid = 10000
       VGrid_Lower= 0.1e0_double; VGrid_Higher = 6.e0_double
    case default
       STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
    end select
    allocate(Posterior(2, nGrid)); Posterior = 1.e0_double !*!
    do i =1, nGrid
       Posterior(1,i) = VGrid_Lower + (i-1)*((VGrid_Higher-VGrid_Lower)/(nGrid-1))
    end do
    !----------------------!

    D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
    !--SET UP REDSHIFT GRID--!
    if(Marginalise_Redshift_Distribution) then
       do z = 1, nRedshift_Sampling
          RedshiftGrid(z) = Redshift_Lower + (z-1)*((Redshift_Higher-Redshift_Lower)/(nRedshift_Sampling-1))
       end do
    end if


!!$    !-----Construct the renormalisation as a function of convergence-----------!
!!$    !--Assumes p(theta,m) is correctly renormalised from the routine that spawns it--!
!!$    Size_Only_Prior = 0.e0_double
!!$
!!$    !--Limits may need edited --!!!
!!$    Survey_Magnitude_Limits = (/minval(PriorMagGrid), maxval(PriorMagGrid)/)
!!$    Survey_Size_Limits = (/minval(PriorSizeGrid), maxval(PriorSizeGrid)/)
!!$
!!$    allocate(ConvergenceGrid(nConvergenceGrid)); ConvergenceGrid = 0.e0_double
!!$    allocate(Renormalisation_by_Convergence(size(ConvergenceGrid))); Renormalisation_by_Convergence = 0.e0_double
!!$    do i = 1, size(ConvergenceGrid)
!!$       ConvergenceGrid(i) = KappaGridLower + (i-1)*((KappaGridHigher- KappaGridLower)/(size(ConvergenceGrid)-1))
!!$       Renormalisation_Size_Limits = Survey_Size_Limits/(1.e0_double+ConvergenceGrid(i))
!!$       Renormalisation_Magnitude_Limits = Survey_Magnitude_Limits + 2.17e0_double*ConvergenceGrid(i)
!!$       
!!$       !--Construct Size only PDF--!
!!$       Size_Only_Prior = 0.e0_double
!!$       do j = 1, size(PriorSizeGrid)
!!$          Size_Only_Prior(j) = Integrate(PriorMagGrid, Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
!!$       end do
!!$
!!$       Renormalisation_by_Convergence(i) = Integrate(PriorSizeGrid, Size_Only_Prior, 2, lim = Renormalisation_Size_Limits)
!!$    end do

    !--Construct Size Only Prior that will be used in the Size Only Routine, p(theta)--!
    do j = 1, size(PriorSizeGrid)
       Size_Only_Prior(j) = Integrate(PriorMagGrid, Prior(:,j), 2, lim = (/minval(PriorMagGrid),maxval(PriorMagGrid)/))
    end do

    !--Construct Magntide Only Prior p(m)--!
    do j = 1, size(PriorMagGrid)
       Mag_Only_Prior(j) = Integrate(PriorSizeGrid, Prior(j,:), 2, lim = (/minval(PriorSizeGrid),maxval(PriorSizeGrid)/))
    end do

    !--Construct Size Given Mag Prior, p(theta|m)
    !~~~~Found from p(m) = int dtheta p(theta,m) so p(theta|m) = p(theta,m)/p(m)
    !~~~~ In pratice, take p(theta,m), for each m evaluation renormalise by the integral over theta, or p(m) [calc'd above] for that m
    Size_Given_Mag_Prior = 0.e0_double
    do j = 1, size(PriorMagGrid)
       Size_Given_Mag_Prior(j,:) = Prior(j,:)/Mag_Only_Prior(j)
    end do


    !--Start of Posterior Routines--!
    allocate(Posterior_perGalaxy(size(Cat%RA), size(Posterior,2))); Posterior_perGalaxy = 1.e-100_double
    allocate(Effective_Convergence(size(Cat%RA), size(Posterior,2))); Effective_Convergence = 0.e0_double
    Convergence_Per_Cluster = 0.e0_double

    print *, 'Getting Posterior for Cluster...'
    n_Default_Source_Redshift_Used = 0
    do c = 1, size(Effective_Convergence,1) !-Loop over galaxies-!
       
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
          allocate(Posterior_perGalaxy_Redshift(size(Posterior_perGalaxy,1), size(Posterior_perGalaxy,2), 1)); Posterior_perGalaxy_Redshift = 0.e0_double
          allocate(Sigma_Crit(size(Cat%RA),1)); Sigma_Crit = -1.e0_double

          call Monte_Carlo_Redshift_Sampling_SigmaCritical(Cat, Lens_Redshift, Sigma_Crit_MC)
          Sigma_Crit(:,1) = Sigma_Crit_MC
          deallocate(Sigma_Crit)

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
          call CH08_redshift_distribution(Cat%MF606W(c), RedshiftGrid, RedshiftPDF)

       else
          STOP 'DM_Profile_Variable_Posterior - Both MC and Redshift Distribution methods set - this cannot be'
       end if
       if(any(Sigma_Crit(c,:) < 0.e0_double)) STOP 'DM_Profile_Variable_Posterior - Invalid Sigma Critical Entered, negative'

       Distance_from_Mass_Center = dsqrt( (Cat%RA(c)-Lens_Position(1))**2.e0_double + (Cat%Dec(c)-Lens_Position(2))**2.e0_double ) !-in Degrees-!
       Distance_from_Mass_Center = (D_l*Distance_from_Mass_Center*(3.142e0_double/(180.e0_double))) !-in Mpc/h-! 


       do i = 1, size(Effective_Convergence,2) !--Loop over posterior values--!
          do z = 1, size(RedshiftPDF)
             select case(Mass_Profile)
             case(1) !-Flat-!
                Effective_Convergence(c,i) = Posterior(1,i)/Sigma_Crit(c,z)
             case(2) !-SIS-!
                Effective_Convergence(c,i) = SMD_SIS(Posterior(1,i), Distance_From_Mass_Center)/(Sigma_Crit(c,z)*1.e18_double)
                Convergence_Per_Cluster(:,c)  = (/Distance_From_Mass_Center,  SMD_SIS(1073.e0_double*1073.e0_double, Distance_From_Mass_Center), Sigma_Crit(c,z)/) !---TESTING--!
             case(3)
                Effective_Convergence(c,i) = SMD_NFW(Distance_From_Mass_Center, Lens_Redshift, Posterior(1,i))/(Sigma_Crit(c,z)*1.e18_double)
             case default
                STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
             end select
             !--Linearly Interpolate--!
             select case (Posterior_Method)
             case(1) !--Size Only--!
                if(lnSize_Prior) then
                   if(Analyse_with_Physical_Sizes) then
                      !DELETE                   Posterior_perGalaxy(c,i) = Linear_Interpolation(Prior(1,:), Prior(2,:), dlog(Cat%Physical_Sizes(c))-Effective_Convergence(c,i))
                   else
                      Posterior_perGalaxy(c,i) = Linear_Interpolation(PriorSizeGrid(:), Size_Only_Prior(:), dlog(Cat%Sizes(c))-Effective_Convergence(c,i))
                   end if
                else
                   if(Analyse_with_Physical_Sizes) then
                      !DELETE                   Posterior_perGalaxy(c,i) = Linear_Interpolation(Prior(1,:), Prior(2,:), Cat%Physical_Sizes(c)/(1.e0_double+Effective_Convergence(c,i)))*(1.e0_double/(1+Effective_Convergence(c,i)))
                   else
                      Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interpolation(PriorSizeGrid(:), Size_Only_Prior(:), Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)))*(1.e0_double/(1+Effective_Convergence(c,i)))
                   end if
                end if
             case(2)!-Size and Magnitude-!
                if(lnSize_Prior) then
                   STOP 'DM_Profile_Variable_Posterior - lnSize with Size_Magnitude Method - I cannae do that captain!'
                else
                   !--Uses distributions of the apparent size--!
                   if(Cat%MF606W(c)+2.17e0_double*Effective_Convergence(c,i) > maxval(PriorMagGrid)) then
                      !--If outside the range in which p(theta|m) is evaluated, then set to zero--!
                      Posterior_perGalaxy_Redshift(c,i,z) = 0.e0_double
                   else
                      Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interpolation(PriorSizeGrid, Size_Given_Mag_Prior(Nearest_Neighbour_Index(PriorMagGrid,Cat%MF606W(c)+2.17e0_double*Effective_Convergence(c,i)),:), Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)))*Linear_Interpolation(PriorMagGrid, Mag_Only_Prior, Cat%MF606W(c)+2.17e0_double*Effective_Convergence(c,i))*(1.e0_double/(1+Effective_Convergence(c,i)))
                   end if
                end if
             case default
                STOP 'DM_Profile_Variable_Posterior - Error in choosing method of finding posterior'
             end select
             !--Convergence Dependant Renormalisation--!
             !-Need to be careful about renormalisation being zero, esp. as a result of falling outside interpolation limits
!             Posterior_perGalaxy_Redshift(c,i,z) = Posterior_perGalaxy_Redshift(c,i,z)/Linear_Interp(Effective_Convergence(c,i), ConvergenceGrid, Renormalisation_by_Convergence)
          end do
       end do

       !--Get Posterior on Size for each galaxy by marginalising over the redshift distribution--!
       if(size(RedshiftPDF) == 1) then
          Posterior_perGalaxy(c,:) = Posterior_perGalaxy_Redshift(c,:,1)
       else
          !--Integrate over the redshift Information--!
          do i = 1, size(Posterior_perGalaxy,2) !--Loop over redshift information--!
!             if(all(Posterior_perGalaxy_Redshift(c,i,:)*RedshiftPDF(:) == 0.e0_double)) STOP 'Product of redshift dependent posterior and PDF is zero'
             Posterior_perGalaxy(c,i) = TrapInt(RedshiftGrid, Posterior_perGalaxy_Redshift(c,i,:)*RedshiftPDF(:))
          end do
          if(all(Posterior_perGalaxy == 0.e0_double)) STOP 'Posterior per galxy is zero!'
       end if
       deallocate(Posterior_perGalaxy_Redshift, RedshiftPDF, Sigma_Crit)

       !--Renomralised to the posterior per galaxy
       Renorm = 0.e0_double
       do j = 1, size(Posterior_perGalaxy,2)-1
          Renorm = Renorm + 0.5e0_double*(Posterior_perGalaxy(c,j) + Posterior_perGalaxy(c,j+1))*(Posterior(1,j+1)-Posterior(1,j))
       end do
!       if(Renorm == 0.e0_double) STOP 'DM_Profile_Variable_Posterior - Error Renormalising Galaxy posteriors - Posterior Empty'
       Posterior_perGalaxy(c,:) = Posterior_perGalaxy(c,:)/Renorm
       Renorm = 0.e0_double
    end do
    if(n_Default_Source_Redshift_Used > 0) print *, '****** Used the default redshift for:', n_Default_Source_Redshift_Used, ' galaxies ********'

    if(Debug_Mode .and. Output_Posterior_Per_Galaxy) then
       Filename = trim(adjustl(Output_Prefix))//'Convergence_Per_Cluster.dat'
       open(unit = 18, file = Filename)
       write(fmtstring,'(I1)') 3
       do i =1, size(Cat%RA)
          write(18, '(I5, x, '//trim(adjustl(fmtstring))//'(e14.7,x))') Cat%Galaxy_Number(i), Convergence_Per_Cluster(:,i)
       end do
      close(18)
      print *, 'Output Convergence per Cluster to: ', trim(adjustl(Filename))

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
!!$    do c = 1, size(Effective_Convergence,1) !-Loop over galaxies-!
!!$       Posterior(2,:) = Posterior(2,:)*(Posterior_perGalaxy(c,:)/(0.5e0_double*maxval(Posterior_perGalaxy(c,:))))
!!$    end do
!!$
!!$    Filename= trim(adjustl(Output_Prefix))//'Posterior_Combined.dat'
!!$    open(unit = 82, file = Filename)
!!$    write(fmtstring,'(I1)') 2
!!$    do i =1, size(Posterior,2)
!!$       write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posterior(1,i), Posterior(2,i)
!!$    end do
!!$    close(82)
!!$    print *, 'Output Combined Posterior to: ', trim(adjustl(Filename))
!!$
!!$!    do j = 1, size(Posterior,2)
!!$!       if(Posterior(2,j) /= 0.e0_double) Posterior(2,j) = dexp(Posterior(2,j))
!!$!    end do
!!$
!!$

!!$!    where(Posterior(2,:) /= 0.e0_double)
!!$!       Posterior(2,:) = dexp(Posterior(2,:))
!!$!    end where

!!$
!!$    if(any(isNAN(Posterior(2,:)))) then
!!$       print *, 'Any NaNs in aperture posterior?', any(isNAN(Posterior(2,:))), count(isNAN(Posterior(2,:)) == .true.)
!!$       print *, 'Stopping'
!!$       STOP
!!$    END if
!!$
!!$
!!$    !--Renormalise--!
!!$    Renorm = 0.e0_double
!!$    do j = 1, size(Posterior,2)-1
!!$       Renorm = Renorm + 0.5e0_double*(Posterior(2,j) + Posterior(2,j+1))*(Posterior(1,j+1)-Posterior(1,j))
!!$    end do
!!$    Posterior(2,:) = Posterior(2,:)/Renorm

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



    deallocate(Effective_Convergence, Posterior_perGalaxy)

    !--On Successful Completion delete Poster per galaxy as large file--!
!!$    inquire(file = trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat', exist = here)
!!$    if(here == .false.) STOP "Posterior per galaxy doesn't exist, stopping before accidental deletion"
!!$    call system('rm '//trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat')

  end subroutine DM_Profile_Variable_Posterior

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
