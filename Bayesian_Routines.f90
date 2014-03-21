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

  subroutine Posterior_Statistics(PosteriorGrid, Posterior, MeanVal, ModeVal, Error)
    use Statistics, only: mean, mode_distribution, variance_distribution
    real(double),intent(in)::Posterior(:), PosteriorGrid(:)
    real(double), intent(out),optional::MeanVal, ModeVal, Error

    INTERFACE
       subroutine Posterior_Statistics(PosteriorGrid, Posterior, MeanVal, ModeVal, Error)
         use Param_Types
         real(double),intent(in)::Posterior(:), PosteriorGrid(:)
         
         real(double), intent(out),optional::MeanVal, ModeVal, Error
       end subroutine Posterior_Statistics
    END INTERFACE

    if(present(ModeVal)) ModeVal = mode_distribution(PosteriorGrid, Posterior)
    if(present(Error)) Error = dsqrt(variance_distribution(PosteriorGrid, Posterior))
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
    !--Main routine that returns the Posteriors over all apertures.
    !--Ap_Radius in DEGREES
    !--To Do:

    use Distributions; use Cosmology, only: angular_diameter_distance_fromRedshift
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
    real(double),allocatable::Cluster_Mean(:), Cluster_Variance(:), Cluster_Mode(:)
    real(double):: Scale_Mass, Scale_Mass_Error
    character(10):: Mass_String, Error_Mass_String

    !--MAgnitude Binning (by Absolute Magnitude)--!
    integer::nMag = 1
    real(double),allocatable::MagBins(:,:)
    integer::Magnitude_Binning_Type = 1 !-1:Absolute, 2:Apparent (MF606)-!

    character(200):: Output_File_Prefix

    logical::use_lnSize_Prior = .false.

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
    !-Set up size grid-!
!!$    allocate(SizeGrid(nSizeGrid)); SizeGrid = 0.e0_double
!!$    dSizeGrid = (SizeGrid_Upper-SizeGrid_Lower)/(1.e0_double*nSizeGrid)
!!$    do i = 1, nSizeGrid
!!$       SizeGrid(i) = SizeGrid_Lower + i*dSizeGrid
!!$    end do

    if(Magnitude_Binning_Type == 1) then
       call Calculate_Bin_Limits_by_equalNumber(Cat%Absolute_Magnitude, nMag, MagBins)
    elseif(Magnitude_Binning_Type == 2) then
       call Calculate_Bin_Limits_by_equalNumber(Cat%MF606W, nMag, MagBins)
    else
       STOP 'DM_Profile_Variable_Posteriors_CircularAperture - Invalid Magnitude Type Specified stopping...'
    END if

    !--If a catalogue for the blank field is passed in, then use this catalogue to get the intrinsic distributions--!
    if(present(Blank_Field_Catalogue)) then
       call get_Size_Distribution_MagnitudeBinning_byCatalogue(MagBins, SizeGrid, SizePrior_byMag, Blank_Field_Catalogue, use_Physical_sizes = Analyse_with_Physical_Sizes, Magnitude_Type = 1, ln_size_Distribution = use_lnSize_Prior)
    else
       call get_Size_Distribution_MagnitudeBinning_byCatalogue(MagBins, SizeGrid, SizePrior_byMag, Cat, use_Physical_sizes = Analyse_with_Physical_Sizes, Magnitude_Type = 1, ln_size_Distribution = use_lnSize_Prior)
    end if

    !--Get Prior by getting refernce for each bin, and binning each reduced cat using the same definition (mag)--!
    !--Produce Posterior for each mag bin--!

    do ap = 1, size(Ap_Cats)
       call bin_catalogue_by_magnitude(Ap_Cats(Ap), MagBins, BCat)

       print *, 'Aperture Catalogue binned in the same way with magnitudes: Aperture: ', Ap, ' has ', size(Ap_Cats(Ap)%RA), ' galaxies binned accordjing to:', BCat%Occupation

       do m = 1, size(MagBins,1)
          !--Set up Prior_Single for that Aperture and Mag Bin-!
          allocate(Prior_Single(2, size(SizePrior_byMag,2))); Prior_Single(1,:) = SizeGrid; Prior_Single(2,:) = SizePrior_byMag(m,:)

          !-Get Posterior for that Magnitude-!
          write(Output_File_Prefix,'(I2)') Ap
          Output_File_Prefix = trim(adjustl(Bayesian_Routines_Output_Directory))//'Aperture_'//trim(adjustl(Output_File_Prefix))//'_'
          call DM_Profile_Variable_Posterior(BCat%Cat(m), Surface_Mass_Profile, Prior_Single, Lens_Redshift, Ap_Pos(Ap,:), Posterior_Single, Output_File_Prefix, use_lnSize_Prior)


          if(any(isNAN(Posterior_Single(2,:)))) then
             print *, 'Any NaNs in returned posterior?', any(isNAN(Posterior_Single(2,:))), count(isNAN(Posterior_Single(2,:)) == .true.)
             print *, 'Stopping'
             STOP
          END if


          if(m==1) then
             allocate(Aperture_Posterior_byMag(size(MagBins,1), 2, size(Posterior_Single,2))); Aperture_Posterior_byMag = 0.e0_double !-Assumes posterior on the same grid as the prior
          end if
          Aperture_Posterior_byMag(m,1,:) = Posterior_Single(1,:); Aperture_Posterior_byMag(m,2,:) = Posterior_Single(2,:)
          

          deallocate(Posterior_Single, Prior_Single)
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

    open(51, file = trim(Bayesian_Routines_Output_Directory)//'Mass_Estimates.dat')
    write(51, '(A)') '# Aperture, Mode, Mode_Error, Mean, Mean_Error'

    print *, '!----------------------------------------------------------------------------------!'
    do Ap = 1, size(Posteriors,1)
       !--Get Mean, Mode, Variance--!
       call Posterior_Statistics(Posteriors(Ap,1,:), Posteriors(Ap,2,:), Cluster_Mean(Ap), Cluster_Mode(Ap), Cluster_Variance(Ap))
       if(Surface_Mass_Profile == 1) then
          !--Convert into the correct units (i.e. Msun/h)--!
          Cluster_Mean(Ap) = 1.e18_double*Cluster_Mean(Ap); Cluster_Mode(Ap) =  1.e18_double*Cluster_Mode(Ap); Cluster_Variance(Ap) = 1.e18_double*Cluster_Variance(Ap)
       end if

       D_l =  angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
       Area = 3.142e0_double*(D_l*((3.142e0_double*Ap_Radius(Ap))/180.e0_double))**2.e0_double
!       print *, 'Typical Conversion from Sigma to Mass (x10^18 Msun/h):', Area
       write(Error_Mass_String, '(e8.2)') Cluster_Variance(Ap)
       write(Mass_String, '(e8.2)') Cluster_Mean(Ap)
       if(Surface_Mass_Profile == 1) write(*,'(A,I2,4(A))') 'Mean Sigma_0 of Cluster :', Ap, ', is (Mean): ', trim(Mass_String), ' +- ', Error_Mass_String
       if(Surface_Mass_Profile == 2) write(*,'(A,I2,4(A))') 'Mean Velocity Dispersion of Cluster :', Ap, ', is (Mean): ', trim(Mass_String), ' +- ', Error_Mass_String
       if(Surface_Mass_Profile == 3) write(*,'(A,I2,4(A))') 'Mean Virial Radius of Cluster :', Ap, ', is (Mean): ', trim(Mass_String), ' +- ', Error_Mass_String
       
       write(Mass_String, '(e8.2)') Cluster_Mode(Ap)
       write(*,'(4(A))') '                                   and (Mode): ', trim(Mass_String), ' +- ', Error_Mass_String

       if(Surface_Mass_Profile == 3) then
          Scale_Mass =  NFW_Mass_withinRadius(D_l*iAp_Radius(Ap)*(3.142e0_double/180.e0_double), Lens_Redshift, r_200 = Cluster_Mean(Ap))
          Scale_Mass_Error = Error_NFW_Mass_withinRadius(D_l*iAp_Radius(Ap)*(3.142e0_double/180.e0_double), Cluster_Mean(Ap), Cluster_Variance(Ap), Lens_Redshift)
          !!Error?!!
       else       
          call  Integrated_Mass_Within_Radius(Surface_Mass_Profile, D_l*iAp_Radius(Ap)*(3.142e0_double/180.e0_double), Cluster_Mean(Ap), Cluster_Variance(Ap), Scale_Mass, Scale_Mass_Error)
       end if
       write(Mass_String, '(e8.2)') Scale_Mass; write(Error_Mass_String, '(e8.2)') Scale_Mass_Error
       print *, "Mass of Cluster within 1' is (Mean):", Mass_String, ' +- ', Error_Mass_String !-Error?-!
       if(Surface_Mass_Profile == 3) then
          Scale_Mass =  NFW_Mass_withinRadius(D_l*iAp_Radius(Ap)*(3.142e0_double/180.e0_double), Lens_Redshift, r_200 = Cluster_Mode(Ap))
          Scale_Mass_Error = Error_NFW_Mass_withinRadius(D_l*iAp_Radius(Ap)*(3.142e0_double/180.e0_double), Cluster_Mode(Ap), Cluster_Variance(Ap), Lens_Redshift)
       else
          call  Integrated_Mass_Within_Radius(Surface_Mass_Profile, D_l*iAp_Radius(Ap)*(3.142e0_double/180.e0_double), Cluster_Mode(Ap), Cluster_Variance(Ap), Scale_Mass, Scale_Mass_Error)
       end if
       write(Mass_String, '(e8.2)') Scale_Mass; write(Error_Mass_String, '(e8.2)') Scale_Mass_Error
       print *, "Mass of Cluster within 1' is (Mode):", Mass_String, ' +- ', Error_Mass_String !-Error?-!
       print *, ' '

       write(51, '(I2,4(e9.3,x))') Ap, Area*Cluster_Mode(Ap)*1.e18_double,  Area*1.e18_double*Cluster_Variance(Ap), Area*Cluster_Mean(Ap)*1.e18_double, Area*1.e18_double*Cluster_Variance(Ap)
    end do
    print *, '!----------------------------------------------------------------------------------!'

  end subroutine DM_Profile_Variable_Posteriors_CircularAperture



  subroutine DM_Profile_Variable_Posterior(Cat, Mass_Profile, Prior, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior)
    use cosmology, only:angular_diameter_distance_fromRedshift; use MC_Redshift_Sampling, only: Monte_Carlo_Redshift_Sampling_SigmaCritical; use Mass_Profiles, only:SMD_SIS, SMD_NFW
    !--Returns the Bayesian posterior for the DM Profile Variables, as defined by Mass_Profile
    !-For Mass_Profile = 1 (Flat): Sigma_0
    !-                   2 (SIS) : Velocity_Dispersion**2
    !-                   3 (NFW) : Virial Radius (r200)
    !--Prior is the prior distribution of p(R), or possibly p(lnR)[not yet implemented]
    !--Cat entered already assumed to be only the galaxies that are to be considered (e.g. in aperture), and well-described by posterior (e.g. same mag bin) 

    !---TO DO:
    !--Conversion from posterior to surface mass density (possibly not here, but later)

    type(Catalogue)::Cat
    integer,intent(in)::Mass_Profile !-1:Flat, 2:SIS-!
    real(double),intent(in)::Lens_Redshift
    real(double),intent(in)::Lens_Position(2)
    real(double),intent(inout)::Prior(:,:) 
    real(double),allocatable,intent(Out)::Posterior(:,:) !-Grid/Posterior, Value-!
    character(*), intent(in)::Output_Prefix
    logical,intent(in)::lnSize_Prior

    integer::i, c, j

    !--Variable Grid Declarations-!
    integer::nGrid = 1000
    real(double)::VGrid_Lower, VGrid_Higher
    
    real(double),allocatable::Posterior_perGalaxy(:,:) !-Galaxy, Value-! - Uses Same Grid as overall Posterior -
    real(double),allocatable:: Effective_Convergence(:,:) !-Galaxy, Value-! -Same Size as Posterior Grid, may need reevaluated per galaxy

    logical:: MC_Sigma_Critical = .false.
    real(double),allocatable::Sigma_Crit(:)
    real(double)::D_l, D_s, D_ls
    real(double)::Distance_from_Mass_Center

    real(double)::Renorm

    logical:: Output_Posterior_Per_Galaxy = .true.
    character(7)::fmtstring
    character(200)::Filename
    logical::here

    !--Testing Declarations--!
    real(double),dimension(3,size(Cat%RA)):: Convergence_per_Cluster
    integer:: n_Default_Source_Redshift_Used

    !--Renomralise the Priors her to ensure that they are renormalised and to avoid issues later:--!
    Renorm = 0.e0_double
    do j = 1, size(Prior,2)-1
       Renorm = Renorm + 0.5e0_double*(Prior(2,j) + Prior(2,j+1))*(Prior(1,j+1)-Prior(1,j))
    end do
    Prior(2,:) = Prior(2,:)/Renorm
    Renorm = 0.e0_double
    

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

    if(MC_Sigma_Critical) then
       call Monte_Carlo_Redshift_Sampling_SigmaCritical(Cat, Lens_Redshift, Sigma_Crit)
    else
       D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
       allocate(Sigma_Crit(size(Cat%RA))); Sigma_Crit = -1.e0_double
    end if
    !--
    allocate(Posterior_perGalaxy(size(Cat%RA), size(Posterior,2))); Posterior_perGalaxy = 1.e-100_double
    allocate(Effective_Convergence(size(Cat%RA), size(Posterior,2))); Effective_Convergence = 0.e0_double

    Convergence_Per_Cluster = 0.e0_double
    print *, 'Getting Posterior for Cluster...'
    do c = 1, size(Effective_Convergence,1) !-Loop over galaxies-!
       !--Ignore Galaxies with redshift less than the foreground-!
       if(Cat%Redshift(c) < Lens_Redshift) cycle
       if(MC_Sigma_Critical == .false.) then
          if(Cat%Redshift(c) >= 0.e0_double) then
             D_s = angular_diameter_distance_fromRedshift(0.e0_double, Cat%Redshift(c))
             D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, Cat%Redshift(c))
          else
             n_Default_Source_Redshift_Used = n_Default_Source_Redshift_Used + 1
             D_s = angular_diameter_distance_fromRedshift(0.e0_double, Default_Source_Redshift)
             D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, Default_Source_Redshift)
             !--Check number of times used--!
          end if
          Sigma_Crit(c) = 1.66492e0_double*(D_s/(D_l*D_ls)) !-(10^18 M_Sun/h)-!
       end if
       if(Sigma_Crit(c) < 0.e0_double) STOP 'DM_Profile_Variable_Posterior - Invalid Sigma Critical Entered, negative'

       Distance_from_Mass_Center = dsqrt( (Cat%RA(c)-Lens_Position(1))**2.e0_double + (Cat%Dec(c)-Lens_Position(2))**2.e0_double ) !-in Degrees-!
       Distance_from_Mass_Center = (D_l*Distance_from_Mass_Center*(3.142e0_double/(180.e0_double))) !-in Mpc/h-! 


       do i = 1, size(Effective_Convergence,2) !--Loop over posterior values--!
          select case(Mass_Profile)
          case(1) !-Flat-!
             Effective_Convergence(c,i) = Posterior(1,i)/Sigma_Crit(c)
          case(2) !-SIS-!
             Effective_Convergence(c,i) = SMD_SIS(Posterior(1,i), Distance_From_Mass_Center)/(Sigma_Crit(c)*1.e18_double)
             Convergence_Per_Cluster(:,c)  = (/Distance_From_Mass_Center,  SMD_SIS(1073.e0_double*1073.e0_double, Distance_From_Mass_Center), Sigma_Crit(c)/) !---TESTING--!
          case(3)
             Effective_Convergence(c,i) = SMD_NFW(Distance_From_Mass_Center, Lens_Redshift, Posterior(1,i))/(Sigma_Crit(c)*1.e18_double)
          case default
             STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
          end select
          !--Linearly Interpolate--!
          if(lnSize_Prior) then
             if(Analyse_with_Physical_Sizes) then
                Posterior_perGalaxy(c,i) = Linear_Interpolation(Prior(1,:), Prior(2,:), dlog(Cat%Physical_Sizes(c))-Effective_Convergence(c,i))
             else
                Posterior_perGalaxy(c,i) = Linear_Interpolation(Prior(1,:), Prior(2,:), dlog(Cat%Sizes(c))-Effective_Convergence(c,i))
             end if
          else
             if(Analyse_with_Physical_Sizes) then
                Posterior_perGalaxy(c,i) = Linear_Interpolation(Prior(1,:), Prior(2,:), Cat%Physical_Sizes(c)/(1.e0_double+Effective_Convergence(c,i)))*(1.e0_double/(1+Effective_Convergence(c,i)))
             else
                Posterior_perGalaxy(c,i) = Linear_Interpolation(Prior(1,:), Prior(2,:), Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)))*(1.e0_double/(1+Effective_Convergence(c,i)))
             end if
          end if
       end do
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
