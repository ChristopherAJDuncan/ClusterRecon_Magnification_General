program  Bayesian_DM_Profile_Constraints
  !--Uses Bayesian Means to Calculate Posterior for the Dark Matter Profile variables--!
  !--Method must have redshifts for every galaxy (20Jan2014)--!
  use Param_Types; use Catalogues
  implicit none

  character(200)::Output_Directory = 'Bayesian_DM_Profile_Constriants_Output/'

  real(double),allocatable::Cluster_Posteriors(:,:,:)
  type(Catalogue)::Catt
  
  character(120)::Catalogue_Directory, Catalogue_Filename
  integer,dimension(:),allocatable::Catalogue_Cols

  real(double),allocatable::Cluster_Pos(:,:), Cluster_Aperture_Radius(:)

  real(double)::Default_Source_Redshift = 1.4e0_double

  allocate(Cluster_Pos(4,2)); Cluster_Pos = 0.e0_double
  Cluster_Pos(1,:) = (/149.1099e0_double,-9.9561e0_double/)
  Cluster_Pos(2,:) = (/148.9889e0_double,-9.9841e0_double/)
  Cluster_Pos(3,:) = (/149.1424e0_double,-10.1666e0_double/)
  Cluster_Pos(4,:) = (/148.9101e0_double,-10.1719e0_double/)
  allocate(Cluster_Aperture_Radius(4)); Cluster_Aperture_Radius = 1.e0_double/60.e0_double !-Degrees-!

  !--Read in the Catalogue with reshifts--!                                                                                                     
  !## 1: STAGES shear, 2: COMBO17, 3:RRG, 4: Mocks_STAGES; 5:Mocks_COMBO!
  call common_Catalogue_directories(5, Catalogue_Directory, Catalogue_Filename, Catalogue_Cols)
  call catalogue_readin(Catt, trim(adjustl(Catalogue_Directory))//trim(adjustl(Catalogue_Filename)), 'Tr(J)', Catalogue_Cols)

  call Clip_Sizes(Catt, (/0.e0_double, 20.e0_double/) )
  !call convert_Size_from_Pixel_to_Physical(Catt)
  call Cut_by_Magnitude(Catt, 23.e0_double) !-Taken from CH08 P1435-!
  call Monte_Carlo_Redshift_Sampling(Catt)
  call Cut_By_PhotoMetricRedshift(Catt, 0.21e0_double) !--Cut out foreground--!

  call DM_Profile_Variable_Posteriors_CircularAperture(Catt, Cluster_Pos, Cluster_Aperture_Radius, Cluster_Posteriors)

contains

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

  subroutine DM_Profile_Variable_Posteriors_CircularAperture(Cat, Ap_Pos, Ap_Radius, Posteriors)
    !--Main routine that returns the Posteriors over all apertures.
    !--To Do:
    !-------: Edit to produce magnitude binning at this level (Control over use of priors)
    !-------: Produce priors using data fed in from catalogue (e.g. when using mock, produce the prior from the mock rather than the reference catalogue)

    use Statistics, only: mean, variance_Distribution, mode_distribution; use Distributions; use Cosmology, only: angular_diameter_distance_fromRedshift
    type(Catalogue), intent(in)::Cat
    real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:)
    real(double),intent(out),allocatable::Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-! - Second dimension allows for a different grid for each aperture

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
    character(10):: Mass_String, Error_Mass_String

    !--MAgnitude Binning (by Absolute Magnitude)--!
    integer::nMag = 1
    real(double),allocatable::MagBins(:,:)

    character(200):: Output_File_Prefix

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

    call Calculate_Bin_Limits_by_equalNumber(Cat%Absolute_Magnitude, nMag, MagBins)

    call get_Size_Distribution_MagnitudeBinning_byCatalogue(MagBins, SizeGrid, SizePrior_byMag, Cat, .true.)
    
    !--Get Prior by getting refernce for each bin, and binning each reduced cat using the same definition (mag)--!
    !--Produce Posterior for each mag bin--!

    do ap = 1, size(Ap_Cats)
       call bin_catalogue_by_magnitude(Ap_Cats(Ap), MagBins, BCat)

       print *, 'Aperture Catalogue binned in the same way with magnitudes: Aperture:', Ap, ' has ', size(Ap_Cats(Ap)%RA), ' galaxies binned accordjing to:', BCat%Occupation

       do m = 1, size(MagBins,1)
          !--Set up Prior_Single for that Aperture and Mag Bin-!
          allocate(Prior_Single(2, size(SizePrior_byMag,2))); Prior_Single(1,:) = SizeGrid; Prior_Single(2,:) = SizePrior_byMag(m,:)

          !-Get Posterior for that Magnitude-!
          write(Output_File_Prefix,'(I2)') Ap
          Output_File_Prefix = 'Aperture_'//trim(adjustl(Output_File_Prefix))//'_'
          call DM_Profile_Variable_Posterior(BCat%Cat(m), 1, Prior_Single, Lens_Redshift, Posterior_Single, Output_File_Prefix)

          print *, 'Any NaNs in returned posterior?', any(isNAN(Posterior_Single(2,:))), count(isNAN(Posterior_Single(2,:)) == .true.)

          if(m==1) then
             print *, 'Allocating Posterior by mag:', size(MagBins,1), 2, size(Posterior_Single,2)
             allocate(Aperture_Posterior_byMag(size(MagBins,1), 2, size(Posterior_Single,2))); Aperture_Posterior_byMag = 0.e0_double !-Assumes posterior on the same grid as the prior
          end if
          Aperture_Posterior_byMag(m,1,:) = Posterior_Single(1,:); Aperture_Posterior_byMag(m,2,:) = Posterior_Single(2,:)
          

          deallocate(Posterior_Single, Prior_Single)
       end do
       call Binned_Catalogue_Destruct(BCat)

       !--Recombine Posterior for that Aperture, over all mag bins - Assumes all on the same grid-!
       if(Ap==1) then
          allocate(Posteriors(size(Ap_Cats),2, size(Aperture_Posterior_byMag,3))); Posteriors = 0.e0_double
       end if
       Posteriors(Ap,1,:) = Aperture_Posterior_byMag(1,1,:);
       !--Combine as log posteriors--!
       do m = 1, size(MagBins,1)
          where(Aperture_Posterior_byMag(m,2,:) /= 0.e0_double)
             Posteriors(Ap,2,:) = Posteriors(Ap,2,:) + log(Aperture_Posterior_byMag(m,2,:))
          end where
!          Posteriors(Ap,2,:) = Posteriors(Ap,2,:)*Aperture_Posterior_byMag(m,2,:)
       end do
       

       !--Output Posterior per Mag Bin--!
       write(apString, '(I1)') Ap
       open(38, file =trim(Output_Directory)//'lnPosterior_Aperture_'//trim(apString)//'_byMagBin.dat')
       !--Header--!!!!!!!!!!!!!!!!!
       do j = 1, size(MagBins,1)
          write(38, '(A1, 2(e14.7,x))') '#', MagBins(j,:)
       end do
       write(fmtstring,'(I2)') size(Aperture_Posterior_byMag,1)+1 !--Assumes all on the same grid--!
       do j =1, size(Aperture_Posterior_byMag,3)
          write(38, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Aperture_Posterior_byMag(1,1,j), Aperture_Posterior_byMag(:,2,j)
       end do
       close(38)
       print *, 'Output file to: ',trim(Output_Directory)//'lnPosterior_Aperture_'//trim(apString)//'_byMagBin.dat' 

       deallocate(Aperture_Posterior_byMag)
    end do

    where(Posteriors(:,2,:) /= 0.e0_double)
       Posteriors(:,2,:) = dexp(Posteriors(:,2,:))
    end where

    !--Renormalising--!
    do Ap = 1, size(Posteriors,1)
       Renorm = 0.e0_double
       do j = 1, size(Posteriors,3)-1
          Renorm = Renorm + 0.5e0_double*(Posteriors(Ap,2,j) + Posteriors(Ap,2,j+1))*(Posteriors(Ap,1,j+1)-Posteriors(Ap,1,j))
       end do
       Posteriors(Ap,2,:) = Posteriors(Ap,2,:)/Renorm
    end do

    !--Output Posterior--!
    open(unit = 51, file = trim(Output_Directory)//'Posterior_per_Aperture.dat')
    write(fmtstring,'(I2)') size(Posteriors,1)+1 !-Assumes all posteriors described by the same posterior grid-!
    do j = 1, size(Posteriors,3)
       write(51, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posteriors(1,1,j), Posteriors(:,2,j)
    end do
    close(51)
    print *,'Output file to: ', trim(Output_Directory)//'Posterior_per_Aperture.dat'

    !--Convert to a Mass as A*Sigma(mean)--!
    allocate(Cluster_Mean(size(Posteriors,1))); Cluster_Mean = 0.e0_double
    allocate(Cluster_Mode(size(Posteriors,1))); Cluster_Mode = 0.e0_double
    allocate(Cluster_Variance(size(Posteriors,1))); Cluster_Variance = 0.e0_double
    print *, '!----------------------------------------------------------------------------------!'
    do Ap = 1, size(Posteriors,1)
       !--Get Mean, Mode, Variance--!
       Cluster_Mean(Ap) = mean(Posteriors(Ap,2,:), Posteriors(Ap,1,:))
       Cluster_Mode(Ap) = mode_distribution(Posteriors(Ap,1,:), Posteriors(Ap,2,:))
       Cluster_Variance(Ap) = dsqrt(variance_distribution(Posteriors(Ap,1,:), Posteriors(Ap,2,:)))

       D_l =  angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
       Area = 3.142e0_double*(D_l*((3.142e0_double*Ap_Radius(Ap))/180.e0_double))**2.e0_double


       write(Error_Mass_String, '(e8.2)') Area*1.e18_double*Cluster_Variance(Ap)
       write(Mass_String, '(e8.2)') Area*Cluster_Mean(Ap)*1.e18_double
       print *, 'Mass of Cluster :', Ap, ', is: ', trim(Mass_String), ' +- ', Error_Mass_String
       write(Mass_String, '(e8.2)') Area*Cluster_Mode(Ap)*1.e18_double
       print *, '                                    and ', trim(Mass_String), ' +- ', Error_Mass_String
    end do
    print *, '!----------------------------------------------------------------------------------!'

  end subroutine DM_Profile_Variable_Posteriors_CircularAperture



  subroutine DM_Profile_Variable_Posterior(Cat, Mass_Profile, Prior, Lens_Redshift, Posterior, Output_Prefix)
    use cosmology, only:angular_diameter_distance_fromRedshift
    !--Returns the Bayesian posterior for the DM Profile Variables, as defined by Mass_Profile
    !-For Mass_Profile = 1 (Flat): Sigma_0
    !-                   2 (SIS) : Velocity_Dispersion**2
    !--Prior is the prior distribution of p(R), or possibly p(lnR)[not yet implemented]
    !--Cat entered already assumed to be only the galaxies that are to be considered (e.g. in aperture), and well-described by posterior (e.g. same mag bin) 

    !---TO DO:
    !----Variable Grid boundary values
    !----Radial Dpendance in Eff_Convergence for SIS

    type(Catalogue)::Cat
    integer,intent(in)::Mass_Profile !-1:Flat, 2:SIS-!
    real(double),intent(in)::Lens_Redshift
    real(double),intent(inout)::Prior(:,:) 
    real(double),allocatable,intent(Out)::Posterior(:,:) !-Grid/Posterior, Value-!
    character(*), intent(in)::Output_Prefix

    integer::i, c, j

    !--Variable Grid Declarations-!
    integer::nGrid = 1000
    real(double)::VGrid_Lower, VGrid_Higher
    
    real(double),allocatable::Posterior_perGalaxy(:,:) !-Galaxy, Value-! - Uses Same Grid as overall Posterior -
    real(double),allocatable:: Effective_Convergence(:,:) !-Galaxy, Value-! -Same Size as Posterior Grid, may need reevaluated per galaxy

    real(double)::Sigma_Crit
    real(double)::D_l, D_s, D_ls
    
    real(double)::Renorm

    logical:: Output_Posterior_Per_Galaxy = .true.
    character(7)::fmtstring
    character(200)::Filename

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
       VGrid_Lower = -5.e-2_double; VGrid_Higher = 1.e-1_double !-- SMD ~ Masses 10^12 -> 10^15 Msun/h, in Units of 10^18 Msun/h  
    case(2) !-SIS-!
       VGrid_Lower= 0.e0_double; VGrid_Higher = 0.e0_double    !--------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case default
       STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
    end select
    allocate(Posterior(2, nGrid)); Posterior = 0.e0_double
    do i =1, nGrid
       Posterior(1,i) = VGrid_Lower + (i-1)*((VGrid_Higher-VGrid_Lower)/(nGrid-1))
    end do
    !----------------------!
   
    D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
    !--
    allocate(Posterior_perGalaxy(size(Cat%RA), size(Posterior,2))); Posterior_perGalaxy = 0.e0_double
    allocate(Effective_Convergence(size(Cat%RA), size(Posterior,2))); Effective_Convergence = 0.e0_double
    do c = 1, size(Effective_Convergence,1) !-Loop over galaxies-!
       if(Cat%Redshift(c) >= 0.e0_double) then
          D_s = angular_diameter_distance_fromRedshift(0.e0_double, Cat%Redshift(c))
          D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, Cat%Redshift(c))
       else
          D_s = angular_diameter_distance_fromRedshift(0.e0_double, Default_Source_Redshift)
          D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, Default_Source_Redshift)
          !--Check number of times used--!
       end if
       Sigma_Crit = 1.66492e0_double*(D_s/(D_l*D_ls)) !-(10^18 M_Sun/h)-!
       do i = 1, size(Effective_Convergence,2)
          select case(Mass_Profile)
          case(1) !-Flat-!
             Effective_Convergence(c,i) = Posterior(1,i)/Sigma_Crit
          case(2) !-SIS-!
             Effective_Convergence(c,i) = Posterior(1,i)/Sigma_Crit !all multiplied by A(r) for each galaxy for an SIS profile [A(r) = 1/(2Gr)] Careful about units, and radius definition (r is distance, not angle) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          case default
             STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
          end select
          !--This bit is true when the prior is in terms of p(R), would need modified here if in  terms of p(lnR)
          !--Linearly Interpolate--!
          Posterior_perGalaxy(c,i) = Linear_Interpolation(Prior(1,:), Prior(2,:), Cat%Physical_Sizes(c)/(1.e0_double+Effective_Convergence(c,i)))*(1.e0_double/(1+Effective_Convergence(c,i)))
       end do
       !--Renomralised to the posterior per galaxy
       Renorm = 0.e0_double
       do j = 1, size(Posterior_perGalaxy,2)-1
          Renorm = Renorm + 0.5e0_double*(Posterior_perGalaxy(c,j) + Posterior_perGalaxy(c,j+1))*(Posterior(1,j+1)-Posterior(1,j))
       end do
       Posterior_perGalaxy(c,:) = Posterior_perGalaxy(c,:)/Renorm
       Renorm = 0.e0_double


       do j = 1, size(Posterior,2)
          if(Posterior_perGalaxy(c,j) /= 0.e0_double) Posterior(2,j) = Posterior(2,j)+log(Posterior_perGalaxy(c,j)) !-Assumes they are on the same grid. Alternatively, ln(Posterior) could be summed
       end do
    end do

    !--Renormalised ln(P) by an arbirtary factor chosen to ensure ln(P) does not get so large as to introduce numerical errors. This will be accounted for when P = exp(lnP) is renormalised--!
    Posterior(2,:) = Posterior(2,:)/maxval(Posterior(2,:))



    if(Output_Posterior_Per_Galaxy) then
       Filename = trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat'
       open(unit = 17, file = Filename)
       write(fmtstring,'(I7)') size(Posterior_perGalaxy,1)+1
       do i =1, size(Posterior,2)
          write(17, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posterior(1,i), Posterior_perGalaxy(:,i)
       end do
      close(17)
      print *, 'Output Posterior per galaxy to:', trim(adjustl(Filename))
   end if

   print *, 'Any NaNs in aperture ln posterior?', any(isNAN(Posterior(2,:))), count(isNAN(Posterior(2,:)) == .true.)

    Filename= trim(adjustl(Output_Prefix))//'lnPosterior_Combined.dat'
    open(unit = 82, file = Filename)
    write(fmtstring,'(I1)') 2
    do i =1, size(Posterior,2)
       write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posterior(1,i), Posterior(2,i)
    end do
    close(82)
    print *, 'Output Combined lnPosterior to: ', trim(adjustl(Filename))

    do j = 1, size(Posterior,2)
       if(Posterior(2,j) /= 0.e0_double) Posterior(2,j) = dexp(Posterior(2,j))
    end do


!!$    where(Posterior(2,:) /= 0.e0_double)
!!$       Posterior(2,:) = dexp(Posterior(2,:))
!!$    end where

    print *, 'Any NaNs in aperture posterior?', any(isNAN(Posterior(2,:))), count(isNAN(Posterior(2,:)) == .true.)

    Filename= trim(adjustl(Output_Prefix))//'Posterior_Combined.dat'
    open(unit = 82, file = Filename)
    write(fmtstring,'(I1)') 2
    do i =1, size(Posterior,2)
       write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posterior(1,i), Posterior(2,i)
    end do
    close(82)
    print *, 'Output Combined Posterior to: ', trim(adjustl(Filename))

    !--Renormalise--!
    Renorm = 0.e0_double
    do j = 1, size(Posterior,2)-1
       Renorm = Renorm + 0.5e0_double*(Posterior(2,j) + Posterior(2,j+1))*(Posterior(1,j+1)-Posterior(1,j))
    end do
    Posterior(2,:) = Posterior(2,:)/Renorm


    deallocate(Effective_Convergence, Posterior_perGalaxy)

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

end program Bayesian_DM_Profile_Constraints
