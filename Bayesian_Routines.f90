module  Bayesian_Routines
  !--Uses Bayesian Means to Calculate Posterior for the Dark Matter Profile variables--!
  !--Method must have redshifts for every galaxy (20Jan2014)--!
  use Param_Types; use Catalogues
  implicit none

  character(200):: Bayesian_Routines_Output_Directory
  
  real(double)::Default_Source_Redshift = 1.4e0_double
  logical:: Analyse_with_Physical_Sizes = .false. !#Default Value#

  logical::Combine_log_Posteriors = .true. !-If False, then combined by multiplication-!

  logical,private::Debug_Mode = .true.

  integer::Surface_Mass_Profile = 3 !-1:Flat, 2:SIS, 3:NFW-!

  !--Method: 1: Size-Only, 2: Size-Magnitude--!
  integer:: Posterior_Method = 2 !--Eventually pass in--!
  logical:: use_KDE_Smoothed_Distributions = .true., KDE_onTheFly = .false.
  logical::use_lnSize_Prior = .false.
  real(double),dimension(2):: Survey_Magnitude_Limits = (/23.e0_double, 27.5e0_double/), Survey_Size_Limits = (/0.e0_double, 100.e0_double/)    

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

  subroutine DM_Profile_Variable_Posteriors_CircularAperture(Cat, Ap_Pos, Ap_Radius, Posteriors, Distribution_Directory, Blank_Field_Catalogue)
    use Statistics, only: mean, mode_distribution, variance_distribution; use Mass_profiles
    use Distributions; use Cosmology, only: angular_diameter_distance_fromRedshift;
    !--Main routine that returns the Posteriors over all apertures.
    !--Ap_Radius in DEGREES
    !--If Blank_Field_Catalogue is entered, then intrinsic distributions are produced form this Catalogue
    !--To Do: 
    !~~Recent edits to use a joint size magnitude have ignored magnitude binning


    type(Catalogue), intent(in)::Cat
    real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:)
    real(double),intent(out),allocatable::Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-! - Second dimension allows for a different grid for each aperture
    character(*), intent(in):: Distribution_Directory
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

    real(double),allocatable::MagGrid(:) !--Discardable for now as marginalised over--!

    INTERFACE
         subroutine DM_Profile_Variable_Posteriors_CircularAperture(Cat, Ap_Pos, Ap_Radius, Posteriors, Distribution_Directory, Blank_Field_Catalogue)
           use Param_Types; use Catalogues
           !--Main routine that returns the Posteriors over all apertures.
           !--Ap_Radius in DEGREES
           type(Catalogue), intent(in)::Cat
           real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:)
           real(double),intent(out),allocatable::Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-! - Second dimension allows for a different grid for each aperture
           character(*), intent(in):: Distribution_Directory

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

    do i =1, size(Ap_Cats)
       print *, 'Aperture:', i, ' contains:', size(Ap_Cats(I)%RA), ' galaxies'
    end do

    print *, count(Cat%Redshift >= 0.e0_double), ' of ',size(Cat%Redshift), ' galaxies have redshift information'

    !--If a catalogue for the blank field is passed in, then use this catalogue to get the intrinsic distributions--!

    if(KDE_OnTheFly == .false.) then
       if(present(Blank_Field_Catalogue)) then
          write(*,'(A)') 'Producing Distribution from Catalogue'
          call return_Size_Magnitude_Distribution(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Distribution_Directory, Blank_Field_Catalogue)
       else
          write(*,'(A)') 'Reading in distribution from:', Distribution_Directory
          call return_Size_Magnitude_Distribution(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Distribution_Directory)
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

       if(KDE_OnTheFly) then
          call DM_Profile_Variable_Posterior_Catalogue(Ap_Cats(Ap), Surface_Mass_Profile, Blank_Field_Catalogue, Lens_Redshift, Ap_Pos(Ap,:), Posterior_Single, Output_File_Prefix, use_lnSize_Prior) 
       else
          call DM_Profile_Variable_Posterior(Ap_Cats(Ap), Surface_Mass_Profile, MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Lens_Redshift, Ap_Pos(Ap,:), Posterior_Single, Output_File_Prefix, use_lnSize_Prior)  
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
       write(Mass_String, '(e9.3)') Cluster_Mean(Ap)
       if(Surface_Mass_Profile == 1) write(*,'(A,I2,6(A))') 'Sigma_0 of Cluster :', Ap, ', is (Mean): ', trim(Mass_String), '   +', trim(Error_Mass_String_Positive), ' -', Error_Mass_String_Negative
       if(Surface_Mass_Profile == 2) write(*,'(A,I2,6(A))') 'Velocity Dispersion of Cluster :', Ap, ', is (Mean): ', trim(Mass_String), '   +', trim(Error_Mass_String_Positive), ' -', Error_Mass_String_Negative
       if(Surface_Mass_Profile == 3) write(*,'(A,I2,6(A))') 'Virial Radius of Cluster :', Ap, ', is (Mean): ', trim(Mass_String), '   +', trim(Error_Mass_String_Positive), ' -', Error_Mass_String_Negative
       
       write(Mass_String, '(e9.3)') Cluster_Mode(Ap)
       write(*,'(6(A))') '                                   and (Mode): ', trim(Mass_String), '   +', trim(Error_Mass_String_Positive), ' -', Error_Mass_String_Negative

       !--Get Mean Mass--!
       call  Integrated_Mass_Within_Radius(Surface_Mass_Profile, D_l*iAp_Radius(Ap)*(3.142e0_double/180.e0_double), Cluster_Mean(Ap), AntiSymm_Variance(Ap,:), Scale_Mass, Scale_Mass_Error, Redshift = Lens_Redshift)
       write(Mass_String, '(e9.3)') Scale_Mass; write(Error_Mass_String_Positive, '(e8.2)') Scale_Mass_Error(2); write(Error_Mass_String_Negative, '(e8.2)') Scale_Mass_Error(1)
       write(*,'(6A)') "Mass of Cluster within 1' is (Mean):", Mass_String, '   +', trim(Error_Mass_String_Positive), ' -', Error_Mass_String_Negative
       !--Get Mode Mass--!
       call  Integrated_Mass_Within_Radius(Surface_Mass_Profile, D_l*iAp_Radius(Ap)*(3.142e0_double/180.e0_double), Cluster_Mode(Ap), AntiSymm_Variance(Ap,:), Scale_Mass, Scale_Mass_Error, Redshift = Lens_Redshift)
       write(Mass_String, '(e9.3)') Scale_Mass; write(Error_Mass_String_Positive, '(e8.2)') Scale_Mass_Error(2); write(Error_Mass_String_Negative, '(e8.2)') Scale_Mass_Error(1)
       write(*,'(6A)') "Mass of Cluster within 1' is (Mode):", Mass_String, '   +', trim(Error_Mass_String_Positive), ' -', Error_Mass_String_Negative
       print *, ' '

       write(51, '(I2,3(e9.3,x))') Ap, Scale_Mass,  Scale_Mass_Error
    end do
    print *, '!----------------------------------------------------------------------------------!'
    close(51)


  end subroutine DM_Profile_Variable_Posteriors_CircularAperture

  subroutine DM_Profile_Variable_Posterior(Cat, Mass_Profile, PriorMagGrid, PriorSizeGrid, Prior, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior)
    use cosmology, only:angular_diameter_distance_fromRedshift; use MC_Redshift_Sampling, only: Monte_Carlo_Redshift_Sampling_SigmaCritical; use Mass_Profiles, only:SMD_SIS, SMD_NFW; use Distributions, only: ch08_redshift_distribution_Array, CH08_Redshift_Distribution_Scalar; use Interpolaters, only: Linear_Interp
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

    integer::i, c, j, z, m

    !--Variable Grid Declarations-!
    integer::nGrid = 1000
    real(double)::VGrid_Lower, VGrid_Higher

    !-Size_Only_Prior contains the priors which depends only on size, which is evalutaed for each redshift and integrates over all magnitudes
    real(double),dimension(:),allocatable:: Size_Only_Prior
    !--Size_Only_Mag_Prior contains the prior for size which evaluated over the mag grid, which will be integrated over. Contains p_[theta_0, m_0|z]*p[z|m_0] for all m in grid
    real(double),dimension(size(Prior,1),size(Prior,2))::Size_Only_Mag_Prior 
    real(double),dimension(size(Prior,1)):: Mag_Only_Prior
    real(double),dimension(size(Prior,1),size(Prior,2)):: Size_Given_Mag_Prior, Kappa_Renormalised_Prior, Survey_Renormalised_Prior
    real(double),allocatable::Posterior_perGalaxy(:,:) !-Galaxy, Posterior-! - Uses Same Grid as overall Posterior -
    real(double),allocatable:: Effective_Convergence(:,:) !-Galaxy, Posterior-! -Same Size as Posterior Grid, may need reevaluated per galaxy

    logical:: MC_Sigma_Critical = .false.
    logical:: Marginalise_Redshift_Distribution = .true.

    logical::Known_Redshift

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
    real(double),parameter::Redshift_Lower = 0.156e0_double, Redshift_Higher = 4.e0_double !!!Edit to Lens_Redshift
    real(double), dimension(nRedshift_Sampling):: RedshiftGrid
    real(double),allocatable::RedshiftPDF(:)
    real(double),allocatable:: Posterior_perGalaxy_Redshift(:,:,:) !-Galaxy, Posterior, Redshift-!

    !--Kappa dependant renormalisation--!
!!    real(double),dimension(2):: Survey_Magnitude_Limits, Survey_Size_Limits !--User defined needs edit to below, these are set by default from distribution-!
    real(double),dimension(2):: Renormalisation_Magnitude_Limits, Renormalisation_Size_Limits
    real(double),allocatable:: ConvergenceGrid(:), Renormalisation_by_Convergence(:), Convergence_Renorm_PerGalaxy(:,:)
    !vv Must be set here vv!
    integer:: nConvergenceGrid = 2000
    real(double):: KappaGridLower = 0e0_double, KappaGridHigher = 2.5e0_double
    integer:: IntegrationFlag = -1000

    !--Testing Declarations--!
    real(double),dimension(3,size(Cat%RA)):: Convergence_per_Cluster
    integer:: n_Default_Source_Redshift_Used, nGal_Ignored_MagLimits, nGal_Ignored_SizeLimits

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

    D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
    !--SET UP REDSHIFT GRID--!
    if(Marginalise_Redshift_Distribution) then
       do z = 1, nRedshift_Sampling
          RedshiftGrid(z) = Redshift_Lower + (z-1)*((Redshift_Higher-Redshift_Lower)/(nRedshift_Sampling-1))
       end do
    end if

    !-----Construct the renormalisation as a function of convergence-----------!
    !--Assumes p(theta,m) is correctly renormalised from the routine that spawns it--!

    !--Limits may need edited --!!!
!    Survey_Magnitude_Limits = (/minval(PriorMagGrid), maxval(PriorMagGrid)/)
    Survey_Magnitude_Limits = (/23.e0_double, 27.5e0_double/) !27.5
!    Survey_Size_Limits = (/minval(PriorSizeGrid), maxval(PriorSizeGrid)/)
    Survey_Size_Limits = (/0.e0_double, 100.e0_double/) 


    !--Renormalise the prior within these size and magnitude limits--!
    print *, 'Renormalisation of the intrinsic distribution:', Integrate(PriorMagGrid, PriorSizeGrid, Prior, 2, lim1 = Survey_Magnitude_Limits, lim2 = Survey_Size_Limits)
    Survey_Renormalised_Prior = Prior/Integrate(PriorMagGrid, PriorSizeGrid, Prior, 2, lim1 = Survey_Magnitude_Limits, lim2 = Survey_Size_Limits)

    !--Remove any probability outside the survey limits - Do we really want to do that? We would surely want to know what galaxies would be brought in outside the limits?--!
!!$    do i = 1, size(PriorMagGrid)
!!$       do j = 1, size(PriorSizeGrid)
!!$          if( (PriorMagGrid(i) < Survey_Magnitude_Limits(1)) .or. (PriorMagGrid(i) > Survey_Magnitude_Limits(2)) .or. (PriorSizeGrid(j) < Survey_Size_Limits(1)) .or. (PriorSizeGrid(j) > Survey_Size_Limits(2)) ) Survey_Renormalised_Prior(i,j) = 0.e0_double
!!$       end do
!!$    end do

    !------This section renormalises the likelihhod, taking into account kappa-dependent mag and size cuts, however this need to be implemented in the prior-----!
    allocate(ConvergenceGrid(nConvergenceGrid)); ConvergenceGrid = 0.e0_double
    allocate(Renormalisation_by_Convergence(size(ConvergenceGrid))); Renormalisation_by_Convergence = 0.e0_double
    allocate(Size_Only_Prior(size(Prior,2))); Size_Only_Prior = 0.e0_double
    do i = 1, size(ConvergenceGrid)
       ConvergenceGrid(i) = KappaGridLower + (i-1)*((KappaGridHigher- KappaGridLower)/(size(ConvergenceGrid)-1))
       Renormalisation_Size_Limits = Survey_Size_Limits/(1.e0_double+ConvergenceGrid(i))
       Renormalisation_Magnitude_Limits = Survey_Magnitude_Limits + 2.17e0_double*ConvergenceGrid(i)

       !--This may not be the way to renormalise this distribution--!
!!$       do j = 1, size(PriorSizeGrid)
!!$          Size_Only_Prior(j) = Integrate(PriorMagGrid, Survey_Renormalised_Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
!!$       end do
!!$       Renormalisation_by_Convergence(i) = Integrate(PriorSizeGrid, (1.e0_double+ConvergenceGrid(i))*Size_Only_Prior, 2, lim = Renormalisation_Size_Limits)
!!$       
       Renormalisation_by_Convergence(i) = Integrate(PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, 2, lim1 = Renormalisation_Magnitude_Limits, lim2 = Renormalisation_Size_Limits)
    end do
    deallocate(Size_Only_Prior)
       
    open(unit = 53, file = trim(adjustl(Output_Prefix))//'Renormalisation_by_Convergence.dat')
    do i = 1, size(Renormalisation_by_Convergence)
       write(53, *) ConvergenceGrid(i), Renormalisation_by_Convergence(i)
    end do
    write(*,'(2(A))') 'File output: ', trim(adjustl(Output_Prefix))//'Renormalisation_by_Convergence.dat'


    !--Construct Size Given Mag Prior, p(theta|m)
    !~~~~Found from p(m) = int dtheta p(theta,m) so p(theta|m) = p(theta,m)/p(m)
    !~~~~ In pratice, take p(theta,m), for each m evaluation renormalise by the integral over theta, or p(m) [calc'd above] for that m
!!$    Size_Given_Mag_Prior = 0.e0_double
!!$    do j = 1, size(PriorMagGrid)
!!$       Size_Given_Mag_Prior(j,:) = Prior(j,:)/Mag_Only_Prior(j)
!!$    end do

    !--Start of Posterior Routines--!
    allocate(Posterior_perGalaxy(size(Cat%RA), size(Posterior,2))); Posterior_perGalaxy = 1.e-100_double
    allocate(Effective_Convergence(size(Cat%RA), size(Posterior,2))); Effective_Convergence = 0.e0_double
    Convergence_Per_Cluster = 0.e0_double

!    allocate(Convergence_Renorm_perGalaxy(size(Cat%RA), size(Posterior,2))); Convergence_Renorm_perGalaxy = 0.e0_double

    write(*,'(A)',advance = 'no') 'Getting Posterior for Cluster '
    if(Posterior_Method == 1) write(*,'(A)') 'using Sizes Only'
    if(Posterior_Method == 2) write(*,'(A)') 'using Sizes and Magnitudes'
    n_Default_Source_Redshift_Used = 0; nGal_Ignored_MagLimits = 0; nGal_Ignored_SizeLimits = 0
    do c = 1, size(Effective_Convergence,1) !-Loop over galaxies-!
       if( (Cat%MF606W(c) > Survey_Magnitude_Limits(2)) .or. (Cat%MF606W(c) < Survey_Magnitude_Limits(1))) then
          nGal_Ignored_MagLimits = nGal_Ignored_MagLimits + 1
          cycle
       end if
       if( (Cat%Sizes(c) > Survey_Size_Limits(2)) .or. (Cat%Sizes(c) < Survey_Size_Limits(1))) then
          nGal_Ignored_SizeLimits = nGal_Ignored_SizeLimits + 1
          cycle
       end if


       if(c == size(Effective_Convergence,1)/2) print *, 'Halfway done for this Aperture..'

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
             
             !--Construct Joint Size-Magnitde Prior which is renormalised according to the convergence value--!
             !-if(i == 1 .and. z ==1 .and. c == 1) print *, 'KAPPA RENORMALISATION TURNED OFF'
             if(i == 1 .and. z ==1 .and. c == 1) print *, 'KAPPA RENORMALISATION TURNED ON'
             Renormalisation_Size_Limits = Survey_Size_Limits/(1.e0_double+Effective_Convergence(c,i))
             Renormalisation_Magnitude_Limits = Survey_Magnitude_Limits + 2.17e0_double*Effective_Convergence(c,i)

             Renorm = Linear_Interp(Effective_Convergence(c,i), ConvergenceGrid, Renormalisation_by_Convergence, ExValue = 1.e30_double)
             if(Renorm == 0) then
                Kappa_Renormalised_Prior = 0.e0_double
             else
                Kappa_Renormalised_Prior = Survey_Renormalised_Prior/Renorm
             end if
             Renorm = 0.e0_double

             
             select case (Posterior_Method)
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
                   if(dlog(Cat%Sizes(c))-Effective_Convergence(c,i) > maxval(PriorSizeGrid) .or. (dlog(Cat%Sizes(c))-Effective_Convergence(c,i) < minval(PriorSizeGrid))) then
                      !--Extrapolation--!
                      Posterior_perGalaxy_Redshift(c,i,z) = 0.e0_double
                      cycle
                   end if

                   !--Find Boundary Indexs for the Size at which the Prior will be evaluated - Aims to improve run-tim by minimising the number of integrations required
                   allocate(Size_Only_Prior(2)); Size_Only_Prior = 0.e0_double
                   do j = 1, size(PriorSizeGrid)-1
                      if( (PriorSizeGrid(j)<= dlog(Cat%Sizes(c))-Effective_Convergence(c,i)) .and. ( PriorSizeGrid(j+1) > dlog(Cat%Sizes(c))-Effective_Convergence(c,i)) ) then
                         Size_Only_Prior(1) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
                         Size_Only_Prior(2) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j+1), 2, lim = Renormalisation_Magnitude_Limits)
                         exit
                      end if
                   end do

                   Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(dlog(Cat%Sizes(c))-Effective_Convergence(c,i), PriorSizeGrid(j:j+1), Size_Only_Prior(:), ExValue = 0.e0_double)
!OBSOLETE (REDUCED GRID)                      Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(dlog(Cat%Sizes(c))-Effective_Convergence(c,i), PriorSizeGrid(:), Size_Only_Prior(:), ExValue = 0.e0_double)
!OBSOLETE                      Posterior_perGalaxy(c,i) = Linear_Interpolation(PriorSizeGrid(:), Size_Only_Prior(:), dlog(Cat%Sizes(c))-Effective_Convergence(c,i))*RedshiftPDF(z)
                else

                   !--Find Boundary Indexs for the Size at which the Prior will be evaluated - Aims to improve run-time by minimising the number of integrations required
                   if(Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)) > maxval(PriorSizeGrid) .or. (Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)) < minval(PriorSizeGrid))) then
                      !--Extrapolation--!
                      Posterior_perGalaxy_Redshift(c,i,z) = 0.e0_double
                      cycle
                   end if

                   allocate(Size_Only_Prior(2)); Size_Only_Prior = 0.e0_double
                   do j = 1, size(PriorSizeGrid)-1
                      if( (PriorSizeGrid(j)<= Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i))) .and. ( PriorSizeGrid(j+1) > Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i))) ) then
                         Size_Only_Prior(1) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
                         Size_Only_Prior(2) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j+1), 2, lim = Renormalisation_Magnitude_Limits)
                         exit
                      end if
                   end do

                   Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)), PriorSizeGrid(j:j+1), Size_Only_Prior(:), ExValue =  0.e0_double)*(1.e0_double/(1+Effective_Convergence(c,i)))

!!$                   if(Posterior_perGalaxy_Redshift(c,i,z) == 0.e0_double) then
!!$                      print *, 'Poster per gal zero, check for extrapolation:', PriorSizeGrid(j:j+1), Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i))
!!$                      print *, Linear_Interp(Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)), PriorSizeGrid(j:j+1), Size_Only_Prior(:), ExValue =  0.e0_double)
!!$                      print *, Size_Only_Prior
!!$                      print *, Renormalisation_Magnitude_Limits, Effective_Convergence(c,i)
!!$                      print *, i, Cat%Redshift(c)
!!$                  !    read(*,*)
!!$                   end if

!OBSOLETE (REDUCEDGRID)                      Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)), PriorSizeGrid(:), Size_Only_Prior(:), ExValue =  0.e0_double)*(1.e0_double/(1+Effective_Convergence(c,i)))
                end if
                deallocate(Size_Only_Prior)
             case(2)!-Size and Magnitude-!
                if(Known_Redshift) then
                   RedshiftPDF(z) = 1.e0_double
                else
                   RedshiftPDF(z) = CH08_Redshift_Distribution_Scalar(Cat%MF606W(c) + 2.17e0_double*Effective_Convergence(c,i), RedshiftGrid(z))
                end if
                if(lnSize_Prior) then
                   STOP 'DM_Profile_Variable_Posterior - lnSize with Size_Magnitude Method - I cannae do that captain!'
                else
                   !--Uses distributions of the apparent size--!
                   if(Cat%MF606W(c)+2.17e0_double*Effective_Convergence(c,i) > maxval(PriorMagGrid)) then
                      !--If outside the range in which p(theta|m) is evaluated, then set to zero--!
                      Posterior_perGalaxy_Redshift(c,i,z) = 0.e0_double
                   else
                      Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(Cat%MF606W(c)+2.17e0_double*Effective_Convergence(c,i), Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)), PriorMagGrid, PriorSizeGrid, Kappa_Renormalised_Prior, ExValue = 0.e0_double)*(1.e0_double/(1+Effective_Convergence(c,i)))*RedshiftPDF(z)
                   end if
                end if
             case default
                STOP 'DM_Profile_Variable_Posterior - Error in choosing method of finding posterior'
             end select
             !--Convergence Dependant Renormalisation--!
             !-Need to be careful about renormalisation being zero. If Convergence is outside the limits, the renormalisation is set to a large value to force zero probability
!!$             Renorm = Linear_Interp(Effective_Convergence(c,i), ConvergenceGrid, Renormalisation_by_Convergence, ExValue = 1.e30_double)
!!$             Convergence_Renorm_PerGalaxy(c, i) = Renorm
!!$             if(Renorm < 0.e0_double) then
!!$                STOP 'DM_Profile_Variable_Posterior - Convergence Dependant Renormalisation is negative, stopping'
!!$             elseif( Renorm == 0.e0_double) then
!!$                !--Likelihood in size mag plane has gone outside the limits here--!
!!$                Posterior_perGalaxy_Redshift(c,i,z) = 0.e0_double
!!$             else
!!$                Posterior_perGalaxy_Redshift(c,i,z) = (Posterior_perGalaxy_Redshift(c,i,z)*(1.e4_double))/(Renorm*1.e4_double)
!!$             end if
!!$             Renorm = 0.e0_double
!             if(Effective_Convergence(c,i) > maxval(ConvergenceGrid)) print *, c, i, Cat%Redshift(c), Effective_Convergence(c,i), Linear_Interp(Effective_Convergence(c,i), ConvergenceGrid, Renormalisation_by_Convergence, ExValue = 1.e30_double),  Posterior(1,i), Posterior_perGalaxy_Redshift(c,i,z)
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

       !--Renormalise the posterior per galaxy
       Renorm = 0.e0_double
       do j = 1, size(Posterior_perGalaxy,2)-1
          Renorm = Renorm + 0.5e0_double*(Posterior_perGalaxy(c,j) + Posterior_perGalaxy(c,j+1))*(Posterior(1,j+1)-Posterior(1,j))
       end do
       if(Renorm == 0.e0_double) THEN
          PRINT *, 'DM_Profile_Variable_Posterior - Error Renormalising Galaxy posteriors - Posterior Empty', C
          !STOP
       END if

!!!!!!TESTING !!!!!!
!!$       if(c == 1) print *, 'Renorm for galxy 1:', Renorm
!!$       open(unit = 45, file = trim(adjustl(Output_Prefix))//'Posterior_Single.dat')
!!$       do j = 1, size(Posterior_perGalaxy,2)
!!$          write(45, *) Posterior(1,j), Posterior_perGalaxy(c,j)
!!$       end do
!!$       print *, 'Output to:', trim(adjustl(Output_Prefix))//'Posterior_Single.dat'
!!$       close(45)
!!$
!!$       if(c == 1) print *, 'Renorm for galxy 1:', Renorm
!!$       open(unit = 45, file = trim(adjustl(Output_Prefix))//'Kappa_Renorm_Single.dat')
!!$       do j = 1, size(Posterior_perGalaxy,2)
!!$          write(45, *) Posterior(1,j), Convergence_Renorm_perGalaxy(c,j)
!!$       end do
!!$       print *, 'Output to:', trim(adjustl(Output_Prefix))//'Kappa_Renorm_Single.dat'
!!$       close(45)
!!$
!!$
!!!! TESTING !!!!!

       Posterior_perGalaxy(c,:) = Posterior_perGalaxy(c,:)/Renorm
       Renorm = 0.e0_double
    end do !--End of galaxy loop-!
    if(n_Default_Source_Redshift_Used > 0) write(*,'(A,I3,A)') '****** Used the default redshift for:', n_Default_Source_Redshift_Used, ' galaxies ********'
    if(nGal_Ignored_SizeLimits > 0) write(*,'(A,I3)') '****** Number of galxies ignored as they fell outside the survey size limits:', nGal_Ignored_SizeLimits
    if(nGal_Ignored_MagLimits > 0) write(*,'(A,I3)') '****** Number of galxies ignored as they fell outside the survey magnitude limits:', nGal_Ignored_MagLimits

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


  subroutine DM_Profile_Variable_Posterior_Catalogue(Cat, Mass_Profile, PriorCatalogue, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior)
    !--Work in progress - attempts to do KDE smoothing on each galaxy in sample--!
    use cosmology, only:angular_diameter_distance_fromRedshift; use MC_Redshift_Sampling, only: Monte_Carlo_Redshift_Sampling_SigmaCritical; use Mass_Profiles, only:SMD_SIS, SMD_NFW; use Distributions, only: ch08_redshift_distribution_Array, CH08_Redshift_Distribution_Scalar; use Interpolaters, only: Linear_Interp; use Smoothing, only: KDE_BiVariate_Gaussian_Scalar; use Statistics, only:Discrete_Covariance
    use Integration, only:TrapInt, Integrate; use Matrix_Methods, only: Matrix_Invert, Determinant
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
    type(Catalogue), intent(in):: PriorCatalogue
    real(double),allocatable,intent(Out)::Posterior(:,:) !-Grid/Posterior, Value-!
    character(*), intent(in)::Output_Prefix
    logical,intent(in)::lnSize_Prior

    integer::i, c, j, z, m

    !--KDE_Smoothing Declarations--!
    real(double),allocatable:: KDE_Gaussian_Covariance(:,:), Data_Vectors(:,:), KDE_Covariance_Inverse(:,:)
    real(double):: KDE_Covariance_Determinant, KDE_Gaussian_Covariance_Reduction = 0.01e0_double !-How much is sig^2 which give KDE width reduced from measured covariance?--!

    !--Variable Grid Declarations-!
    integer::nGrid = 1000
    real(double)::VGrid_Lower, VGrid_Higher

    !-Size_Only_Prior contains the priors which depends only on size, which is evalutaed for each redshift and integrates over all magnitudes
    real(double),dimension(:),allocatable:: Size_Only_Prior
    !--Size_Only_Mag_Prior contains the prior for size which evaluated over the mag grid, which will be integrated over. Contains p_[theta_0, m_0|z]*p[z|m_0] for all m in grid
    real(double),allocatable::Posterior_perGalaxy(:,:) !-Galaxy, Posterior-! - Uses Same Grid as overall Posterior -
    real(double),allocatable:: Effective_Convergence(:,:) !-Galaxy, Posterior-! -Same Size as Posterior Grid, may need reevaluated per galaxy

    logical:: MC_Sigma_Critical = .false.
    logical:: Marginalise_Redshift_Distribution = .true.

    logical::Known_Redshift

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
!!    real(double),dimension(2):: Survey_Magnitude_Limits, Survey_Size_Limits !--User defined needs edit to below, these are set by default from distribution-!
    real(double),dimension(2):: Renormalisation_Magnitude_Limits, Renormalisation_Size_Limits
    real(double),allocatable:: ConvergenceGrid(:), Renormalisation_by_Convergence(:), Convergence_Renorm_PerGalaxy(:,:)
    !vv Must be set here vv!
    integer:: nConvergenceGrid = 2000
    real(double):: KappaGridLower = 0e0_double, KappaGridHigher = 2.5e0_double
    integer:: IntegrationFlag = -1000

    !--Testing Declarations--!
    real(double),dimension(3,size(Cat%RA)):: Convergence_per_Cluster
    integer:: n_Default_Source_Redshift_Used, nGal_Ignored_MagLimits, nGal_Ignored_SizeLimits

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

    D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
    !--SET UP REDSHIFT GRID--!
    if(Marginalise_Redshift_Distribution) then
       do z = 1, nRedshift_Sampling
          RedshiftGrid(z) = Redshift_Lower + (z-1)*((Redshift_Higher-Redshift_Lower)/(nRedshift_Sampling-1))
       end do
    end if

    !--Construct the Covariance that will be used for the KDE Smoothing--!
    allocate(Data_Vectors(2,size(PriorCatalogue%Sizes))); Data_Vectors(1,:) = PriorCatalogue%MF606W; Data_Vectors(2,:) = PriorCatalogue%Sizes
    call Discrete_Covariance(Data_Vectors, KDE_Gaussian_Covariance)
    KDE_Gaussian_Covariance = KDE_Gaussian_Covariance_Reduction*KDE_Gaussian_Covariance
    call Matrix_Invert(KDE_Gaussian_Covariance, KDE_Covariance_Inverse, 'S')
    KDE_Covariance_Determinant = Determinant(KDE_Gaussian_Covariance)
    deallocate(KDE_Gaussian_Covariance)

    !-----Construct the renormalisation as a function of convergence-----------!
    !--Assumes p(theta,m) is correctly renormalised from the routine that spawns it--!

    !--Limits may need edited --!!!
!    Survey_Magnitude_Limits = (/minval(PriorMagGrid), maxval(PriorMagGrid)/)
    Survey_Magnitude_Limits = (/23.e0_double, 27.5e0_double/) !27.5
!    Survey_Size_Limits = (/minval(PriorSizeGrid), maxval(PriorSizeGrid)/)
    Survey_Size_Limits = (/0.e0_double, 100.e0_double/) 


    !--Renormalise the prior within these size and magnitude limits--!
!!$    print *, 'Renormalisation of the intrinsic distribution:', Integrate(PriorMagGrid, PriorSizeGrid, Prior, 2, lim1 = Survey_Magnitude_Limits, lim2 = Survey_Size_Limits)
!!$    Survey_Renormalised_Prior = Prior/Integrate(PriorMagGrid, PriorSizeGrid, Prior, 2, lim1 = Survey_Magnitude_Limits, lim2 = Survey_Size_Limits)

    !------This section renormalises the likelihhod, taking into account kappa-dependent mag and size cuts, however this need to be implemented in the prior-----!
!!$    allocate(ConvergenceGrid(nConvergenceGrid)); ConvergenceGrid = 0.e0_double
!!$    allocate(Renormalisation_by_Convergence(size(ConvergenceGrid))); Renormalisation_by_Convergence = 0.e0_double
!!$    allocate(Size_Only_Prior(size(Prior,2))); Size_Only_Prior = 0.e0_double
!!$    do i = 1, size(ConvergenceGrid)
!!$       ConvergenceGrid(i) = KappaGridLower + (i-1)*((KappaGridHigher- KappaGridLower)/(size(ConvergenceGrid)-1))
!!$       Renormalisation_Size_Limits = Survey_Size_Limits/(1.e0_double+ConvergenceGrid(i))
!!$       Renormalisation_Magnitude_Limits = Survey_Magnitude_Limits + 2.17e0_double*ConvergenceGrid(i)
!!$
!!$       Renormalisation_by_Convergence(i) = Integrate(PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, 2, lim1 = Renormalisation_Magnitude_Limits, lim2 = Renormalisation_Size_Limits)
!!$    end do
!!$    deallocate(Size_Only_Prior)
       
!!$    open(unit = 53, file = trim(adjustl(Output_Prefix))//'Renormalisation_by_Convergence.dat')
!!$    do i = 1, size(Renormalisation_by_Convergence)
!!$       write(53, *) ConvergenceGrid(i), Renormalisation_by_Convergence(i)
!!$    end do
!!$    write(*,'(2(A))') 'File output: ', trim(adjustl(Output_Prefix))//'Renormalisation_by_Convergence.dat'


    !--Start of Posterior Routines--!
    allocate(Posterior_perGalaxy(size(Cat%RA), size(Posterior,2))); Posterior_perGalaxy = 1.e-100_double
    allocate(Effective_Convergence(size(Cat%RA), size(Posterior,2))); Effective_Convergence = 0.e0_double
    Convergence_Per_Cluster = 0.e0_double

!    allocate(Convergence_Renorm_perGalaxy(size(Cat%RA), size(Posterior,2))); Convergence_Renorm_perGalaxy = 0.e0_double

    write(*,'(A)',advance = 'no') 'Getting Posterior for Cluster '
    if(Posterior_Method == 1) write(*,'(A)') 'using Sizes Only'
    if(Posterior_Method == 2) write(*,'(A)') 'using Sizes and Magnitudes'
    n_Default_Source_Redshift_Used = 0; nGal_Ignored_MagLimits = 0; nGal_Ignored_SizeLimits = 0
    do c = 1, size(Effective_Convergence,1) !-Loop over galaxies-!
       if( (Cat%MF606W(c) > Survey_Magnitude_Limits(2)) .or. (Cat%MF606W(c) < Survey_Magnitude_Limits(1))) then
          nGal_Ignored_MagLimits = nGal_Ignored_MagLimits + 1
          cycle
       end if
       if( (Cat%Sizes(c) > Survey_Size_Limits(2)) .or. (Cat%Sizes(c) < Survey_Size_Limits(1))) then
          nGal_Ignored_SizeLimits = nGal_Ignored_SizeLimits + 1
          cycle
       end if


       if(c == size(Effective_Convergence,1)/2) print *, 'Halfway done for this Aperture..'

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
             
             !--Construct Joint Size-Magnitde Prior which is renormalised according to the convergence value--!
             Renormalisation_Size_Limits = Survey_Size_Limits!/(1.e0_double+Effective_Convergence(c,i))
             Renormalisation_Magnitude_Limits = Survey_Magnitude_Limits! + 2.17e0_double*Effective_Convergence(c,i)

!!$             Renorm = Linear_Interp(Effective_Convergence(c,i), ConvergenceGrid, Renormalisation_by_Convergence, ExValue = 1.e30_double)
!!$             if(Renorm == 0) then
!!$                Kappa_Renormalised_Prior = 0.e0_double
!!$             else
!                Kappa_Renormalised_Prior = Survey_Renormalised_Prior!/Renorm
!!$             end if
!!$             Renorm = 0.e0_double

             
             select case (Posterior_Method)
             case(1) !--Size Only--! THIS HAS NOT BEEN DONE YET
                STOP 'I havent coded up the method to do this yet'
                !--Evaluate p_{theta_0, m_0}*p_{z|m_0} for the whole magnitude grid
!!$                do m =  1, size(PriorMagGrid)
!!$                   !--m_0, theta_0--!
!!$                   if(Known_Redshift) then
!!$                      Size_Only_Mag_Prior(m,:) =  Kappa_Renormalised_Prior(m,:)
!!$                   else
!!$                      Size_Only_Mag_Prior(m,:) =  Kappa_Renormalised_Prior(m,:)*CH08_Redshift_Distribution_Scalar(PriorMagGrid(m), RedshiftGrid(z))
!!$                   end if
!!$                end do
!!$
!!$                if(lnSize_Prior) then
!!$                   if(dlog(Cat%Sizes(c))-Effective_Convergence(c,i) > maxval(PriorSizeGrid) .or. (dlog(Cat%Sizes(c))-Effective_Convergence(c,i) < minval(PriorSizeGrid))) then
!!$                      !--Extrapolation--!
!!$                      Posterior_perGalaxy_Redshift(c,i,z) = 0.e0_double
!!$                      cycle
!!$                   end if
!!$
!!$                   !--Find Boundary Indexs for the Size at which the Prior will be evaluated - Aims to improve run-tim by minimising the number of integrations required
!!$                   allocate(Size_Only_Prior(2)); Size_Only_Prior = 0.e0_double
!!$                   do j = 1, size(PriorSizeGrid)-1
!!$                      if( (PriorSizeGrid(j)<= dlog(Cat%Sizes(c))-Effective_Convergence(c,i)) .and. ( PriorSizeGrid(j+1) > dlog(Cat%Sizes(c))-Effective_Convergence(c,i)) ) then
!!$                         Size_Only_Prior(1) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
!!$                         Size_Only_Prior(2) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j+1), 2, lim = Renormalisation_Magnitude_Limits)
!!$                         exit
!!$                      end if
!!$                   end do
!!$
!!$                   Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(dlog(Cat%Sizes(c))-Effective_Convergence(c,i), PriorSizeGrid(j:j+1), Size_Only_Prior(:), ExValue = 0.e0_double)
!!$                else
!!$
!!$                   !--Find Boundary Indexs for the Size at which the Prior will be evaluated - Aims to improve run-time by minimising the number of integrations required
!!$                   if(Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)) > maxval(PriorSizeGrid) .or. (Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)) < minval(PriorSizeGrid))) then
!!$                      !--Extrapolation--!
!!$                      Posterior_perGalaxy_Redshift(c,i,z) = 0.e0_double
!!$                      cycle
!!$                   end if
!!$
!!$                   allocate(Size_Only_Prior(2)); Size_Only_Prior = 0.e0_double
!!$                   do j = 1, size(PriorSizeGrid)-1
!!$                      if( (PriorSizeGrid(j)<= Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i))) .and. ( PriorSizeGrid(j+1) > Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i))) ) then
!!$                         Size_Only_Prior(1) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
!!$                         Size_Only_Prior(2) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j+1), 2, lim = Renormalisation_Magnitude_Limits)
!!$                         exit
!!$                      end if
!!$                   end do
!!$
!!$                   Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)), PriorSizeGrid(j:j+1), Size_Only_Prior(:), ExValue =  0.e0_double)*(1.e0_double/(1+Effective_Convergence(c,i)))
!!$
!!$                end if
!!$                deallocate(Size_Only_Prior)
             case(2)!-Size and Magnitude-!
                if(Known_Redshift) then
                   RedshiftPDF(z) = 1.e0_double
                else
                   RedshiftPDF(z) = CH08_Redshift_Distribution_Scalar(Cat%MF606W(c) + 2.17e0_double*Effective_Convergence(c,i), RedshiftGrid(z))
                end if
                if(lnSize_Prior) then
                   STOP 'DM_Profile_Variable_Posterior - lnSize with Size_Magnitude Method - I cannae do that captain!'
                else
                   !--Uses distributions of the apparent size--!
                      Posterior_perGalaxy_Redshift(c,i,z) = KDE_BiVariate_Gaussian_Scalar(PriorCatalogue%MF606W, PriorCatalogue%Sizes, Cat%MF606W(c)+2.17e0_double*Effective_Convergence(c,i), Cat%Sizes(c)/(1.e0_double+Effective_Convergence(c,i)), Inverse_Covar = KDE_Covariance_Inverse, Det_Covar = KDE_Covariance_Determinant)*(1.e0_double/(1+Effective_Convergence(c,i)))*RedshiftPDF(z)
                end if
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
          print *, 'Galaxy ', c, ' has zero posterior'
          read(*,*)
       end IF

       !--Renormalise the posterior per galaxy
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
    if(nGal_Ignored_MagLimits > 0) write(*,'(A,I3)') '****** Number of galxies ignored as they fell outside the survey magnitude limits:', nGal_Ignored_MagLimits

    if(Debug_Mode .and. Output_Posterior_Per_Galaxy) then
!!$       Filename = trim(adjustl(Output_Prefix))//'Convergence_Per_Cluster.dat'
!!$       open(unit = 18, file = Filename)
!!$       write(fmtstring,'(I1)') 3
!!$       do i =1, size(Cat%RA)
!!$          write(18, '(I5, x, '//trim(adjustl(fmtstring))//'(e14.7,x))') Cat%Galaxy_Number(i), Convergence_Per_Cluster(:,i)
!!$       end do
!!$      close(18)
!!$      print *, 'Output Convergence per Cluster to: ', trim(adjustl(Filename))

!!$       Filename = trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat'
!!$       open(unit = 17, file = Filename)
!!$       write(fmtstring,'(I7)') size(Posterior_perGalaxy,1)+1
!!$       do i =1, size(Posterior,2)
!!$          write(17, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posterior(1,i), Posterior_perGalaxy(:,i)
!!$       end do
!!$      close(17)
!!$      print *, 'Output Posterior per galaxy to: ', trim(adjustl(Filename))
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

  end subroutine DM_Profile_Variable_Posterior_Catalogue



  subroutine return_Size_Magnitude_Distribution(MagGrid, SizeGrid, Dist, Dir, BFCat, do_KDE)
    use IO, only: readin; use Distributions, only: produce_Joint_Size_Magnitude_Distribution; use Integration
    !--Returns the joint size magnitude distribution--!
    !-- If BFCat entered, then the Distribution is calculated from the catalogue, otherwise a read in is attempted --!
    real(double),intent(out),allocatable:: MagGrid(:), SizeGrid(:), Dist(:,:)
    character(*),intent(in):: Dir
    type(Catalogue),intent(in),optional:: BFCat
    logical,optional:: do_KDE

    character(200):: Filename, Output_Filename, Input_Filename
    real(double),allocatable:: Input_Array(:,:)

    character(5)::fmtstring
    integer::i,j
    logical:: iDo_KDE

    !--TESTING--!
    real(double),allocatable::Size_Only(:), Mag_Only(:)
    
    if(present(do_KDE)) then
       iDo_KDE = do_KDE
    else
       iDO_KDE = use_KDE_Smoothed_Distributions
    end if

 
    if(iDO_KDE) then
       Filename = 'MagSize_Distribution_KDE.dat'
    else
       Filename= 'MagSize_Distribution_Histogram.dat'
    end if

    if(present(BFCat)) then
       print *, 'Producing Joint Size Magnitude Distribution from Catalogue'

       call produce_Joint_Size_Magnitude_Distribution(SizeGrid, MagGrid, Dist, BFCat, use_Physical_Sizes = Analyse_with_Physical_Sizes, Magnitude_Type = 2, Output_Dir = trim(Dir), ln_size_Distribution = use_lnSize_Prior, KDE_Smooth = ido_KDE)
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
    
    else
       Input_Filename = trim(Dir)//trim(Filename)

       print *, 'Reading In Distribution from:', trim(Input_Filename)

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

  end subroutine return_Size_Magnitude_Distribution


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
