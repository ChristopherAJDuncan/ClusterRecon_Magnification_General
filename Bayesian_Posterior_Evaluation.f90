module Bayesian_Posterior_Evaluation
  !--Contains the methods to evaluate the Bayesian Posterior for an input of free parameters
  use Param_Types

  !--THIS HAS NOT BEEN IMPLEMENTED
  
  implicit none
  
  !--Function Overload-------------------------------------
  INTERFACE Likelihood_Evaluation_atVirialRadius_perSource
     module procedure Likelihood_atVirialRadius_SingleCluster_perSource, Likelihood_atVirialRadius_MultipleCluster_perSource
  END INTERFACE Likelihood_Evaluation_atVirialRadius_perSource

  INTERFACE lnLikelihood_Evaluation_atVirialRadius_perSourceSample
     module procedure lnLikelihood_atVirialRadius_SingleCluster_perSourceSample, lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample
  END INTERFACE lnLikelihood_Evaluation_atVirialRadius_perSourceSample

  interface Combine_Posteriors
     module procedure Combine_Posteriors_Scalar, Combine_Posteriors_Vector
  end interface Combine_Posteriors

contains  

  real(double) function lnLikelihood_atVirialRadius_SingleCluster_perSourceSample(Alpha, Method, Profile, Cluster_Position, Lens_Redshift, Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, Sigma_Crit)
    !--Sigma_Crit must be entered in units (10^18 MSun/h), and the first dimension should hold sigma critical for all lenses considered (not on a lens redshift grid)
    use Cosmology, only: angular_diameter_distance_fromRedshift
    use Mass_Profiles; use Distributions, only: CH08_Redshift_Distribution_Scalar
    use Interpolaters, only: Linear_Interp; use Integration
    real(double),intent(in):: Alpha
    !--Mass Model Declarations
    integer:: Profile
    real(double), intent(in)::  Lens_Redshift, Cluster_Position(:) !- RA/Dec-!
    !--Galaxy Data Input
    real(double), intent(in):: Theta(:), Magnitude(:), Source_Redshift(:), Position(:,:) !-RA,Dec-!
    !--Supplementary
    integer:: Method
    real(double), intent(in):: Pr_MagGrid(:), Pr_SizeGrid(:), SM_Prior(:,:), M_Prior(:)
    real(double), intent(in):: Size_Limits(:), Magnitude_Limits(:)
    real(double), intent(in):: MagnificationGrid(:), M_Pr_Renormalisation(:), SM_Pr_Renormalisation(:)
    real(double), intent(in):: RedshiftGrid(:), Sigma_Crit(:) !-Source_Redshift-!

    real(double),allocatable:: Posterior_perSource(:)

    !--Internal_Declarations--!
    real(double):: tCluster_Position(1,2), tSigma_Crit(1,size(Sigma_Crit))

    tCluster_Position(1,:) = Cluster_Position; tSigma_Crit(1,:) = Sigma_Crit

    lnLikelihood_atVirialRadius_SingleCluster_perSourceSample = lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample((/Alpha/), Method, Profile, tCluster_Position, (/Lens_Redshift/), Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, tSigma_Crit)

  end function lnLikelihood_atVirialRadius_SingleCluster_perSourceSample

  
  real(double) function lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample(Alpha, Method, Profile, Cluster_Position, Lens_Redshift, Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, Sigma_Crit)
    !--Sigma_Crit must be entered in units (10^18 MSun/h), and the first dimension should hold sigma critical for all lenses considered (not on a lens redshift grid)
    use Cosmology, only: angular_diameter_distance_fromRedshift
    use Mass_Profiles; use Distributions, only: CH08_Redshift_Distribution_Scalar
    use Interpolaters, only: Linear_Interp; use Integration
    real(double),intent(in):: Alpha(:)
    !--Mass Model Declarations
    integer, intent(in):: Profile
    real(double), intent(in)::  Lens_Redshift(:), Cluster_Position(:,:) !-Cluster, RA/Dec-!
    !--Galaxy Data Input
    real(double), intent(in):: Theta(:), Magnitude(:), Source_Redshift(:), Position(:,:) !-RA,Dec-!
    !--Supplementary
    integer, intent(in):: Method
    real(double), intent(in):: Pr_MagGrid(:), Pr_SizeGrid(:), SM_Prior(:,:), M_Prior(:)
    real(double), intent(in):: Size_Limits(:), Magnitude_Limits(:)
    real(double), intent(in):: MagnificationGrid(:), M_Pr_Renormalisation(:), SM_Pr_Renormalisation(:)
    real(double), intent(in):: RedshiftGrid(:), Sigma_Crit(:,:) !-Lens, Source_Redshift-!

    real(double),allocatable:: Posterior_perSource(:)

    integer:: nGal_Ignored_MagLimits, nGal_Ignored_NaN

    integer::c

    real(double):: Res

    allocate(Posterior_perSource(size(Theta))); Posterior_perSource = dsqrt(-1.e0_double)

    do c = 1, size(Theta)
       Posterior_perSource(c) = 1.e0_double

       !--Skip Evaluation if Galaxy Input falls outside limits 
       if( (Magnitude(c) > Magnitude_Limits(2)) .or. (Magnitude(c) < Magnitude_Limits(1))) then
          nGal_Ignored_MagLimits = nGal_Ignored_MagLimits + 1
          cycle
       end if
       if(isNaN(Theta(c)) .or. isNaN(Magnitude(c))) then
          nGal_Ignored_NaN = nGal_Ignored_NaN + 1
          cycle
       end if

       Posterior_perSource(c) = Likelihood_atVirialRadius_MultipleCluster_perSource(Alpha, Method, Profile, Cluster_Position, Lens_Redshift, Theta(c), Magnitude(c), Source_Redshift(c), Position(c,:), Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, Sigma_Crit)
    end do

    !--Recombine
    call Combine_Posteriors_Scalar(Alpha(1), Posterior_perSource(:), .true., Renormalise = .false., Return_lnP = .true., Combined_Posterior = Res)

    lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample = Res

  end function lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample

  !--Evalution of the posterior on a single virial radius grid point, per source
  real(double) function Likelihood_atVirialRadius_SingleCluster_perSource(Alpha, Method, Profile, Cluster_Position, Lens_Redshift, Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, Sigma_Crit)
    !~~~~~DESCRIPTION: Wrapper routine which returns the posterior evulated at a single point, //where only a single cluster is modelled//. Uses a call to the standard multi-cluster routine to get it's result.
    !--Sigma_Crit must be entered in units (10^18 MSun/h), and the first dimension should hold sigma critical for all lenses considered (not on a lens redshift grid)
    real(double),intent(in):: Alpha
    !--Mass Model Declarations
    integer:: Profile
    real(double), intent(in)::  Lens_Redshift, Cluster_Position(:) !-Cluster, RA/Dec-!
    !--Galaxy Data Input
    real(double), intent(in):: Theta, Magnitude, Source_Redshift, Position(2) !-RA,Dec-!
    !--Supplementary
    integer:: Method
    real(double), intent(in):: Pr_MagGrid(:), Pr_SizeGrid(:), SM_Prior(:,:), M_Prior(:)
    real(double), intent(in):: Size_Limits(:), Magnitude_Limits(:)
    real(double), intent(in):: MagnificationGrid(:), M_Pr_Renormalisation(:), SM_Pr_Renormalisation(:)
    real(double), intent(in):: RedshiftGrid(:), Sigma_Crit(:)
    
    !--Internal_Declarations--!
    real(double):: tCluster_Position(1,2), tSigma_Crit(1,size(Sigma_Crit))

    tCluster_Position(1,:) = Cluster_Position; tSigma_Crit(1,:) = Sigma_Crit

    Likelihood_atVirialRadius_SingleCluster_perSource = Likelihood_atVirialRadius_MultipleCluster_perSource((/Alpha/), Method, Profile, tCluster_Position, (/Lens_Redshift/), Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, tSigma_Crit)

  end function Likelihood_atVirialRadius_SingleCluster_perSource

  real(double) function Likelihood_atVirialRadius_MultipleCluster_perSource(Alpha, Method, Profile, Cluster_Position, Lens_Redshift, Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, Sigma_Crit)
    !-This is the workhorse of the routine as it stands. This constructs the likelihood of the data (theta and magnitude) given a set of cluster mass profile free parameters (alpha). The input priors MUST be normalised to unity in the size and magnitude range of the survey/source data set. Nearly everything reuiwred here is discernable fromt eh catalogue itself, or can be obtained using a call to get_Likelihood_Evaluation_Precursors.
    !--Sigma_Crit must be entered in units (10^18 MSun/h), and the first dimension should hold sigma critical for all lenses considered (not on a lens redshift grid)

    use Cosmology, only: angular_diameter_distance_fromRedshift
    use Mass_Profiles; use Distributions, only: CH08_Redshift_Distribution_Scalar
    use Interpolaters, only: Linear_Interp; use Integration
    real(double),intent(in):: Alpha(:)
    !--Mass Model Declarations
    integer:: Profile
    real(double), intent(in)::  Lens_Redshift(:), Cluster_Position(:,:) !-Cluster, RA/Dec-!
    !--Galaxy Data Input
    real(double), intent(in):: Theta, Magnitude, Source_Redshift, Position(2) !-RA,Dec-!
    !--Supplementary
    integer:: Method
    real(double), intent(in):: Pr_MagGrid(:), Pr_SizeGrid(:), SM_Prior(:,:), M_Prior(:)
    real(double), intent(in):: Size_Limits(:), Magnitude_Limits(:)
    real(double), intent(in):: MagnificationGrid(:), M_Pr_Renormalisation(:), SM_Pr_Renormalisation(:)
    real(double), intent(in):: RedshiftGrid(:), Sigma_Crit(:,:) !-Lens, Source_Redshift-!
    
    !--Internal Declarations
    real(double), allocatable:: Posterior_perGalaxy_Redshift(:)
    real(double):: Effective_Magnification
    logical:: Known_Redshift
    integer:: z, m, j, l
    real(double):: Galaxy_Sigma_Critical(size(Sigma_Crit,1))
    real(double):: RedshiftPDF, Renorm
    real(double), dimension(size(Alpha)):: Distance_From_Mass_Center
    real(double):: D_l

    real(double), dimension(2):: Renormalisation_Size_Limits, Renormalisation_Magnitude_Limits
    real(double), dimension(size(SM_Prior,1), size(SM_Prior,2)):: Kappa_Renormalised_Prior, Size_Only_Mag_Prior
    real(double), dimension(size(M_Prior)):: Kappa_Renormalised_MagPrior
    real(double),allocatable:: Size_Only_Prior(:)

    integer:: Galaxy_Posterior_Method !-Converts from inoput posterior method to the individual method, based on the Method input
    integer:: nMagPosterior, nSizePosterior, nSizeMagPosterior
    integer:: nGal_Ignored_SizeLimits


    integer::nZ

    !--Internal Decalrations that could eventually be set/passed in
    logical:: Enforce_Weak_Lensing = .false., Cuts_Renormalise_Likelihood = .true.

    !--Error Catching on Input---!
    if(size(Alpha) /= size(Distance_From_Mass_Center) .or. (size(Alpha) /= size(Lens_Redshift))) then
       print *,  'Likelihood_Evaluation_atVirialRadius - Error on input of Mass Profile variables:'
       print *, 'Sizes (Alpha, Distance, Lens_Redshift):', size(Alpha), size(Distance_From_Mass_Center), size(lens_Redshift)
       STOP
    END if
    
    !--Redshift Knowledge--!
    if(Source_Redshift > 0.e0_double) then
       Known_Redshift = .true.
       nZ = 1
    else
       Known_Redshift = .false.
       nZ = size(RedshiftGrid)
    end if
    
    !--Get Distance of source from each Cluster
    do l = 1, size(Alpha)
       D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift(l))
       Distance_From_Mass_Center(l) = dsqrt( (dabs(Position(1)-Cluster_Position(l, 1))**2.e0_double) + (dabs(Position(2)-Cluster_Position(l, 2))**2.e0_double) ) !-in Degrees-!
       Distance_From_Mass_Center(l) = D_l*Distance_From_Mass_Center(l)*(3.142e0_double/180.e0_double) !-In Mpc/h-!
    end do

    Likelihood_atVirialRadius_MultipleCluster_perSource = 1.e0_double

    !--Set which method will be used
    select case(Method)
    case(1) !--SizeOnly--!                                                                                                                                                                              
       if( (Theta > Size_Limits(2)) .or. (Theta < Size_Limits(1))) then
          nGal_Ignored_SizeLimits = nGal_Ignored_SizeLimits + 1
          return
       end if
       !MagRenorm          iSize_Limits = Size_Limits                                                                                                                                         
       nSizePosterior = nSizePosterior + 1
       Galaxy_Posterior_Method = 1
    case(2)!--SizeMag--!                                                                                                                                                                                    
       if( (Theta > Size_Limits(2)) .or. (Theta < Size_Limits(1))) then
          nGal_Ignored_SizeLimits = nGal_Ignored_SizeLimits + 1
          return
       end if
       !MagRenorm          iSize_Limits = Size_Limits
       nSizeMagPosterior = nSizeMagPosterior + 1
       Galaxy_Posterior_Method = 2
    case(3) !-Magnitude Only--!                                                                                                                                                                             
       nMagPosterior = nMagPosterior + 1
       !MagRenorm          iSize_Limits = (/0.e0_double, 1.e30_double/) !--Should encompase the whole data set                                                                                       
!DELETE       if(present(MagPrior) == .false.) STOP 'Attempting to produce posterior using mag only, but no magnitde prior entered, stopping'
       Galaxy_Posterior_Method = 3
    case(4) !-Size mag above size data limit, Mag-Only below size data limit-!                                                                                                                              
!DELETE       if(present(MagPrior) == .false.) STOP 'Attempting to produce posterior using mag only under size limit, but no magnitde prior entered, stopping'
       STOP 'Posterior Method 4 has been disabled as it needs further thought into its construction. In particular, with respect to the construction of the prior from the size-mag distribution, the effect of prior cuts, and correct renormalisation. The skeleton code has been added for this, but disabled for now. See commented code labeled MagRenorm'
       !--Note, in this case the data-renormalisation should take inot account that the mag only is constructed from a sample which takes small, faint galaxies                                             
       if(Theta < Size_Limits(1)) then
          !--Use Magnitude information only--!
          !MagRenorm             iSize_Limits = (/0.e0_double, Size_Limits(1)/)                                                                                                               
          Galaxy_Posterior_Method = 3
          nMagPosterior = nMagPosterior + 1
       else
          !--Use full magnitude-size information--!
          !MagRenorm             iSize_Limits = Size_Limits
          Galaxy_Posterior_Method = 2
          nSizeMagPosterior = nSizeMagPosterior + 1
       end if
    end select

    allocate(Posterior_perGalaxy_Redshift(nZ)); Posterior_perGalaxy_Redshift = 0.e0_double
    
    do z = 1, nZ !-Vary source redshift
       !--Set Sigma Critical for that galaxy, across all lenses
       if(Known_Redshift) then
          Do l = 1, size(Sigma_Crit,1)
             Galaxy_Sigma_Critical(l) = Linear_Interp(Source_Redshift, RedshiftGrid, Sigma_Crit(l,:))
          end Do
       else
          !--As RedshiftPDF and Sigma_Crit are on the same grid, take value directly
          Galaxy_Sigma_Critical = Sigma_Crit(:,z)
       END if

       if(any(isNaN(Galaxy_Sigma_Critical)) .or. any((Galaxy_Sigma_Critical > huge(1.e0_double)))) STOP 'DM_Profile_Variable_Posterior - FATAL - Galaxy Sigma Critical not set correctly'
       
       !~~~~Ths section will need editing to account for multiple clusters in flat and SIS case
       select case(Profile)
       case(1) !-Flat-!                                                                                                                                                                                      
          STOP 'The Strong Lensing approximation has not been implemented for this profile yet'
          Effective_Magnification = 1.e0_double+2.e0_double*(Alpha(1)/Galaxy_Sigma_Critical(1))
       case(2) !-SIS-!                                                                                                                                                                                       
          STOP 'The Strong Lensing approximation has not been implemented for this profile yet'
          Effective_Magnification = 1.e0_double+2.e0_double*(SMD_SIS(Alpha(1), Distance_From_Mass_Center(1))/(Galaxy_Sigma_Critical(1)*1.e18_double))
       case(3) !-NFW-!                                                                                                                                                                                       
          if(Enforce_Weak_Lensing) then
             STOP 'Weak Lensing disabled for NFW for now, until the equivalent for multiple clusters is coded [dont worry, its actually easy]'
             Effective_Magnification = 1.e0_double + 2.e0_double*(SMD_NFW(Distance_From_Mass_Center(1), Lens_Redshift(1), Alpha(1))/(Galaxy_Sigma_Critical(1)*1.e18_double))
          else
             Effective_Magnification = Total_MagnificationFactor_MultipleClusters(3, Position, Cluster_Position, Alpha, Lens_Redshift, Sigma_Crit = Galaxy_Sigma_Critical*1.e18_double)
!              Effective_Magnification = Magnification_Factor(3, Distance_From_Mass_Center(1), Alpha(1),  Lens_Redshift(1), Galaxy_Sigma_Critical*1.e18_double)
          end if
       case default
          STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
       end select

       
       if(isNaN(Effective_Magnification)) then
          print *, 'Magnification Factor is a NaN:', Distance_From_Mass_Center, Alpha, Lens_Redshift, Galaxy_Sigma_Critical*1.e18_double
          STOP
       end if
       
       if(Effective_Magnification < minval(MagnificationGrid) .or. Effective_Magnification > maxval(MagnificationGrid)) then
          !--Skipping as outside limits on which magnification was evaluated--!                                                                                                                              
          Posterior_perGalaxy_Redshift(z) = 1.e-100_double
          cycle
       end if
       
       !---Set RedshiftPDF, to unity if redshift known, and to a point in a particualr distribution if it is not
       if(Known_Redshift) then
          RedshiftPDF = 1.e0_double
       else
          RedshiftPDF = CH08_Redshift_Distribution_Scalar(Magnitude + 2.5e0_double*dlog10(Effective_Magnification), RedshiftGrid(z))
       end if
       

       !--Construct the correctly renormalised prior
       Renormalisation_Size_Limits = Size_Limits/dsqrt(Effective_Magnification)
       Renormalisation_Magnitude_Limits = Magnitude_Limits + 2.5e0_double*dlog10(Effective_Magnification)
       
       if(Cuts_Renormalise_Likelihood) then
!          if(z ==1) print *, 'KAPPA RENORMALISATION TURNED ON'
          Renorm = Linear_Interp(Effective_Magnification, MagnificationGrid, SM_Pr_Renormalisation, ExValue = 1.e30_double)
          if(Renorm == 0) then
             Kappa_Renormalised_Prior = 0.e0_double
          else
             Kappa_Renormalised_Prior = SM_Prior/Renorm
          end if
          Renorm = 0.e0_double
       else
 !         if(i == 1 .and. z ==1 .and. c == 1) print *, 'KAPPA RENORMALISATION TURNED OFF'
          Kappa_Renormalised_Prior =  SM_Prior
          Kappa_Renormalised_MagPrior = M_Prior
       end if
       

       select case(Galaxy_Posterior_Method)
       case(1) !--Size-Only--!
          
          !--Evaluate p_{theta_0, m_0}*p_{z|m_0} for the whole magnitude grid                                                                                                                                
          do m =  1, size(Pr_MagGrid)
             !--m_0, theta_0--!                                                                                                                                                                              
             if((Pr_MagGrid(m) < Renormalisation_Magnitude_Limits(1)) .or. (Pr_MagGrid(m) > Renormalisation_Magnitude_Limits(2))) then
                !-Since Integrand does not extend over this region anyway-!                                                                                                                                  
                Size_Only_Mag_Prior(m,:) = 1.e-100_double
             else
                if(Known_Redshift) then
                   Size_Only_Mag_Prior(m,:) =  Kappa_Renormalised_Prior(m,:)
                else
                   Size_Only_Mag_Prior(m,:) =  Kappa_Renormalised_Prior(m,:)*CH08_Redshift_Distribution_Scalar(Pr_MagGrid(m), RedshiftGrid(z))
                end if
             end if
          end do
          
          if((Theta/dsqrt(Effective_Magnification) > maxval(PR_SizeGrid)) .or. (Theta/dsqrt(Effective_Magnification) < minval(PR_SizeGrid))) then
             !--Extrapolation--!                                                                                                                                                                          
             Posterior_perGalaxy_Redshift(z) = 1.e-100_double
             cycle
          else
             
             !--Find Boundary Indexs for the Size at which the Prior will be evaluated - Aims to improve run-time by minimising the number of integrations required                                       
             allocate(Size_Only_Prior(2)); Size_Only_Prior = 0.e0_double
             do j = 1, size(PR_SizeGrid)-1
                if( (PR_SizeGrid(j)<= Theta/dsqrt(Effective_Magnification)) .and. ( PR_SizeGrid(j+1) > Theta/dsqrt(Effective_Magnification)) ) then
                   print *, 'Calling Integrate'
                   Size_Only_Prior(1) = Integrate(Pr_MagGrid, Size_Only_Mag_Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
                   Size_Only_Prior(2) = Integrate(Pr_MagGrid, Size_Only_Mag_Prior(:,j+1), 2, lim = Renormalisation_Magnitude_Limits)
                   exit
                end if
             end do
             
             !--EVALUATE SIZE-ONLY LIKELIHOOD
             Posterior_perGalaxy_Redshift(z) = Linear_Interp(Theta/dsqrt(Effective_Magnification), PR_SizeGrid(j:j+1), Size_Only_Prior(:), ExValue =  1.e-100_double)*(1.e0_double/dsqrt(Effective_Magnification))
          end if
          deallocate(Size_Only_Prior)
          
       CASE(2) !--Joint Size-Magnitude Analysis--!
          !--Error Catching
          if(Magnitude + 2.5e0_double*dlog10(Effective_Magnification) < 21.e0_double) then
             print *, 'Possible problem with de-lensed magnitude - falls outwith bright limit of Scrabback Fit (21)'
             print *, Magnitude + 2.5e0_double*dlog10(Effective_Magnification)
          end if
          
          !--EVALUATE SIZE-MAGNITUDE LIKELIHOOD
          Posterior_perGalaxy_Redshift(z) = Linear_Interp(Magnitude+2.5e0_double*dlog10(Effective_Magnification), Theta/dsqrt(Effective_Magnification), Pr_MagGrid, PR_SizeGrid, Kappa_Renormalised_Prior, ExValue = 1.e-100_double)*(1.e0_double/dsqrt(Effective_Magnification))*RedshiftPDF
          
       case(3) !--Magnitude-only
          
          !--Renormalise Magnitude Prior Distribution--!                                                                                                                                                     
          if(Cuts_Renormalise_Likelihood) then
             Renorm = Linear_Interp(Effective_Magnification, MagnificationGrid,  M_Pr_Renormalisation, ExValue = 1.e30_double)
             if(Renorm == 0) then
                Kappa_Renormalised_MagPrior = 0.e0_double
             else
                Kappa_Renormalised_MagPrior = M_Prior/Renorm
             end if
          end if
          
          
          !--Construct p_[m_0] as the magnitude distribution for a sample of galaxies between cuts-corrected survey size limits                                                                              
          !-MagRenorm                                                                                                                                                                                        
!!$                   allocate(Mag_Only_Prior(2)); Mag_Only_Prior = 0.e0_double                                                                                                                                   !!$                   !--Find point where de-lensed magnitude lie on intrinsic distribution - could also use locate.                                                                                              !!$                   do j = 1, size(Pr_MagGrid)-1                                                                                                                                                               
!!$                      if( ( Pr_MagGrid(j)<=  Magnitude + 2.5e0_double*dlog10(Effective_Magnification) ) .and. (  Pr_MagGrid(j+1) > Magnitude + 2.5e0_double*dlog10(Effective_Magnification)) ) then
!!$                         Mag_Only_Prior(1) = Integrate(PR_SizeGrid, Kappa_Renormalised_Prior(j,:), 2, lim = Renormalisation_Size_Limits)                                                                      
!!$                         Mag_Only_Prior(2) = Integrate(PR_SizeGrid, Kappa_Renormalised_Prior(j+1,:), 2, lim = Renormalisation_Size_Limits)                                                                    
!!$                         exit                                                                                                                                                                                  !!$                      end if                                                                                                                                                                                   !!$                   end do                                                                                                                                                                                      !!$                   Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(Magnitude + 2.5e0_double*dlog10(Effective_Magnification), Pr_MagGrid(j:j+1), Mag_Only_Prior, ExValue =  1.e-100_double)*RedshiftPDF(z)  !!$                   deallocate(Mag_Only_Prior)                                                                                                                                                                   
          !---------------------------------------------------------------------------------------------------------------------                                                                             
          
          !--The following gives unbiased results if no size cuts are used, however is known to be biased in the presence of size cuts, I believe that this is due to the fact that the prior itself should depend on the size cuts in the presence of a size-magnitude correlation.                                                                                                                                           
          
          !---EVALUATE MAGNITUDE-ONLY LIKELIHOOD
          Posterior_perGalaxy_Redshift(z) = Linear_Interp(Magnitude + 2.5e0_double*dlog10(Effective_Magnification), Pr_MagGrid, Kappa_Renormalised_MagPrior, ExValue =  1.e-100_double)*RedshiftPDF
          
       case default
          STOP 'DM_Profile_Variable_Posterior - Error in choosing method of finding posterior'
       END select
       
    END do !--End of Redshift Loop
    
    !--Evaluate the Posterior for the free parameter as a function of the input parameter values only
    !-Integrate over redshift
    if(size(Posterior_perGalaxy_Redshift,1) == 1) then
       !--Redshift was taken to be exact                                                                                                                                                                     
       Likelihood_atVirialRadius_MultipleCluster_perSource = Posterior_perGalaxy_Redshift(1)
    else
       !--Integrate over the redshift Information--!                                                                                                                                                         
       Likelihood_atVirialRadius_MultipleCluster_perSource = TrapInt(RedshiftGrid, Posterior_perGalaxy_Redshift(:))
    end if

    !--Deallocate Internal Declarations
    deallocate(Posterior_perGalaxy_Redshift)
    
    
  end function Likelihood_atVirialRadius_MultipleCluster_perSource

  subroutine get_Likelihood_Evaluation_Precursors(Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, RedshiftGrid, Sigma_Crit, MagnificationGrid,Renormalisation_by_Magnification, MagOnly_Renormalisation_by_Magnification, PriorSizeGrid, PriorMagGrid, MagPrior, Prior, Survey_Magnitude_Limits, Survey_Size_Limits, Lens_Redshift, Redshift_Limit, Output_Prefix)
    use Integration; use Cosmology, only: angular_diameter_distance_fromRedshift
    !--Routine that:
    !---Normalises the Prior Distributions to integrate to one within the survey limits
    !---Sets up Sigma_Critical evalaution for all lenses considered
    !---Sets up a redshift grid on which the redshiftPDF (marginalisation) and Sigma_Critical will be evaluated
    !---Determines the normalisation of the prior as a function of magnification factor in the presence of cuts.

    !---This should *always* be called before any of the evaluation routines above

    !--Input
    real(double):: PriorSizeGrid(:), PriorMagGrid(:), MagPrior(:), Prior(:,:) !-Magnitude, Size-!
    real(double):: Survey_Magnitude_Limits(:), Survey_Size_Limits(:), Redshift_Limit
    real(double):: Lens_Redshift(:)
    character(*):: Output_Prefix

    !--Output
    real(double),allocatable:: Survey_Renormalised_Prior(:,:), Survey_Renormalised_MagPrior(:)
    real(double),allocatable:: RedshiftGrid(:)
    real(double),allocatable:: Sigma_Crit(:,:), MagnificationGrid(:),Renormalisation_by_Magnification(:), MagOnly_Renormalisation_by_Magnification(:)

    !--Internal
    real(double),dimension(2):: Renormalisation_Magnitude_Limits, Renormalisation_Size_Limits
    integer:: z, C, i
    real(double):: D_l, D_ls, D_s
    
    integer,parameter:: nRedshift_Sampling = 50
    real(double)::Redshift_Lower, Redshift_Higher = 4.e0_double

    integer:: nMagnificationGrid = 2000
    real(double):: MagFactorGridLower = 1.e0_double, MagFactorGridHigher = 65.e0_double


    allocate(Survey_Renormalised_Prior(size(Prior,1),size(Prior,2))); Survey_Renormalised_Prior = 0.e0_double
    allocate(Survey_Renormalised_MagPrior(size(PriorMagGrid))); Survey_Renormalised_MagPrior = 0.e0_double
    !--Renormalise the prior within these size and magnitude limits--!                                                                                                           
    print *, 'Renormalisation of the intrinsic distribution:', Integrate(PriorMagGrid, PriorSizeGrid, Prior, 2, lim1 = Survey_Magnitude_Limits, lim2 = Survey_Size_Limits)
    Survey_Renormalised_Prior = Prior/Integrate(PriorMagGrid, PriorSizeGrid, Prior, 2, lim1 = Survey_Magnitude_Limits, lim2 = Survey_Size_Limits)
    !--Allow for seperate renormalisation of the magnitude prior. This is required if the magnitude prior is constructed from galaxies which are excluded from the joint size-magnitude analysis, e.g. due to size cuts--!                                                                                                                                                                                                       
    print *, 'Renormalisation of the intrinsic magnitude distribution:', Integrate(PriorMagGrid, MagPrior, 2, lim = Survey_Magnitude_Limits)
    Survey_Renormalised_MagPrior = MagPrior/Integrate(PriorMagGrid, MagPrior, 2, lim = Survey_Magnitude_Limits)
    

    !--Set up the redshift grid--!                                                                                                                                                                                
    Redshift_Lower = Redshift_Limit
    allocate(RedshiftGrid(nRedshift_Sampling)); RedshiftGrid = 0.e0_double
    do z = 1, nRedshift_Sampling
       RedshiftGrid(z) = Redshift_Lower + (z-1)*((Redshift_Higher-Redshift_Lower)/(nRedshift_Sampling-1))
    end do

    !--Get Sigma_Critical for each point on the Redshift PDF grid--!                                                                                                                     
    allocate(Sigma_Crit(size(Lens_Redshift),size(RedshiftGrid))); Sigma_Crit = 0.e0_double
    do C = 1, size(Lens_Redshift)
       do z = 1, size(RedshiftGrid)
          D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift(C))
          D_s = angular_diameter_distance_fromRedshift(0.e0_double, RedshiftGrid(z))
          D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift(C), RedshiftGrid(z))
          Sigma_Crit(C,z) = 1.66492e0_double*(D_s/(D_l*D_ls)) !-(10^18 M_Sun/h)-!
       end do
    end do

    allocate(MagnificationGrid(nMagnificationGrid)); MagnificationGrid = 0.e0_double
    allocate(Renormalisation_by_Magnification(size(MagnificationGrid))); Renormalisation_by_Magnification = 0.e0_double
    allocate(MagOnly_Renormalisation_by_Magnification(size(MagnificationGrid))); MagOnly_Renormalisation_by_Magnification = 0.e0_double
    !Choose MagGrid lower to be the point where the full magnitude range is swept out:                                                                                                                             
    MagFactorGridHigher = 1.05e0_double*(10.e0_double**((Survey_Magnitude_Limits(2)-Survey_Magnitude_Limits(1))/2.5e0_double)) !1.05 gives lee-way!                                                               
    do i = 1, size(MagnificationGrid)
       MagnificationGrid(i) = MagFactorGridLower + (i-1)*((MagFactorGridHigher- MagFactorGridLower)/(size(MagnificationGrid)-1))
       Renormalisation_Size_Limits = Survey_Size_Limits/dsqrt(MagnificationGrid(i))
       Renormalisation_Magnitude_Limits = Survey_Magnitude_Limits + 2.5e0_double*dlog10(MagnificationGrid(i))
       Renormalisation_by_Magnification(i) = Integrate(PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, 2, lim1 = Renormalisation_Magnitude_Limits, lim2 = Renormalisation_Size_Limits)
       !--Edit this code to calculate the magnitude renormalisation over the whole size grid if Posterior_Method == 3, and over (0, Survey_Size_Limit(1)) if Posterior_Method == 4                             
       !-- vv This is only true if there are no cuts on the size-mag prior vv                                                                                                                                  
       !MagRenorm          MagOnly_Renormalisation_by_Magnification(i) =  Integrate(PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, 2, lim1 = Renormalisation_Magnitude_Limits, lim2 = Renormalisation_Size_Limits)                                                                                                                                                                                                                 
       MagOnly_Renormalisation_by_Magnification(i) = Integrate(PriorMagGrid, Survey_Renormalised_MagPrior, 2, lim = Renormalisation_Magnitude_Limits)
    end do

    open(unit = 53, file = trim(adjustl(Output_Prefix))//'Renormalisation_by_Magnification.dat')
    do i = 1, size(Renormalisation_by_Magnification)
       write(53, *) MagnificationGrid(i), Renormalisation_by_Magnification(i)
    end do
    write(*,'(2(A))') 'File output: ', trim(adjustl(Output_Prefix))//'Renormalisation_by_Magnification.dat'


  end subroutine get_Likelihood_Evaluation_Precursors

 subroutine Combine_Posteriors_Scalar(GridValue, Posteriors,  Combine_by_ln, Renormalise, Return_lnP, Combined_Posterior)
    !--Combines the posterior on a single grid value (free parameter alpha), using a call to the "normal" (vector) combined posterior subroutine                                                                  
     real(double), intent(in):: GridValue,Posteriors(:) !-Galaxy-!                                                                                                                                                
     real(double), intent(out):: Combined_Posterior
     logical, intent(in):: Combine_by_ln, Renormalise, Return_lnP

     !--Internal Declarations                                                                                                                                                                                     
     real(double),dimension(1):: tGrid, tCombinedPosterior
     real(double), dimension(size(Posteriors),1):: tPosteriors

     !--Set up internals                                                                                                                                                                                          
     tGrid = GridValue
     tPosteriors(:,1) = Posteriors

     call Combine_Posteriors(tGrid, tPosteriors, Combine_by_ln, .false., Return_lnP, tCombinedPosterior)
     Combined_Posterior = tCombinedPosterior(1)

   end subroutine Combine_Posteriors_Scalar

  subroutine Combine_Posteriors_Vector(PosteriorGrid, Posteriors, Combine_by_ln, Renormalise, Return_lnP, Combined_Posterior)
    !--Combines posteriors by looping over the first dimension--!                                                                                                                                                 
    real(double), intent(in):: PosteriorGrid(:),Posteriors(:,:) !-Galaxy, Grid/Value-!                                                                                                                            
    real(double), intent(out):: Combined_Posterior(:)
    logical, intent(in):: Combine_by_ln, Renormalise, Return_lnP

    integer::c, j
    real(double)::Renorm, Combination_Normalisation
    logical::iDoRenormalise

    integer::nPosteriorsSkipped

    INTERFACE
       subroutine Combine_Posteriors_Vector(PosteriorGrid, Posteriors, Combine_by_ln, Renormalise, Return_lnP, Combined_Posterior)
         use Param_Types
         real(double), intent(in):: PosteriorGrid(:), Posteriors(:,:)
         real(double), intent(out):: Combined_Posterior(:)
         logical, intent(in):: Combine_by_ln, Renormalise, Return_lnP
       END subroutine Combine_Posteriors_Vector
    END INTERFACE
    
    !--Set renormalisation by input. If a single alpha value entered, then do not renormalise                                                                                                                     
    iDoRenormalise = Renormalise
    if(size(PosteriorGrid) == 1) iDoRenormalise = .false.
    
    nPosteriorsSkipped = 0
    if(Combine_by_ln == .false.) then
       Combination_Normalisation = 1.e0_double/maxval(Posteriors)!or 1.e0_double/(0.5e0_double*maxval(Posteriors(c,:))) within loop                                                                               
       Combined_Posterior = 1.e0_double
       do c = 1, size(Posteriors,1) !-Loop over galaxies-!                                                                                                                                                        
          !--Skip if zero (lnP not defined then) or NaN                                                                                                                                                           
          if(all(Posteriors(c,:) == 0.e0_double) .or. all(isNaN(Posteriors(c,:)))) then
             nPosteriorsSkipped = nPosteriorsSkipped + 1
             cycle
          end if

          !--Skip when `renormalised' posteriors are invalid (zero/negative)                                                                                                                                      
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
       if(return_lnP) iDoRenormalise = .false.
       !       print *, 'Combining using logs'                                                                                                                                                                   
       !-Set lnP to a large negative value as default, equivalent to P ~ 0                                                                                                                                        
       Combined_Posterior = -1000.0_double
       Combination_Normalisation = size(Posteriors,1)
       do c = 1, size(Posteriors,1) !-Loop over galaxies-!                                                                                                                                                         
          !-Error Catching--!                        
          if(all(Posteriors(c,:) == 1.e-75_double) .or. all(isNaN(Posteriors(c,:)))) then
             nPosteriorsSkipped = nPosteriorsSkipped + 1
             cycle
          end if
          !--Sum log posteriors--!                                                                                                                                                                                
!!$          print *, 'CombinePosterior, loop:', c                                                                                                                                                                
!!$          print *, Combined_Posterior, dlog(Posteriors(c,:))                                                                                                                                                   
          where(Posteriors(c,:) <= 1.e-100_double)
             Combined_Posterior = Combined_Posterior - 1000.e0_double
          elsewhere
             Combined_Posterior = Combined_Posterior + dlog(Posteriors(c,:)) + 1.e0_double
          end where


          if(any(isNAN(Combined_Posterior(:)))) then
             print *, 'Any NaNs in Combined Posterior?, galaxy:',c, any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
             STOP
          end if
       end do
       !--Convert to PDF, not ln(PDF)--!                                                                                                                                                                          
       if(Return_lnP) then
          return
       else
          if(size(Combined_Posterior)/= 1) then
             Combined_Posterior = dexp(Combined_Posterior - maxval(Combined_Posterior))
          else
             Combined_Posterior = dexp(Combined_Posterior)
          end if
       end if
    end if

    if(any(isNAN(Combined_Posterior(:)))) then
       print *, 'Any NaNs in Combined Posterior?', any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
       print *, 'Stopping'
       STOP
    END if

    if(nPosteriorsSkipped > 0) write(*,'(A,I4,A,I4,A)') '################### Combine_Posteriors - ', nPosteriorsSkipped, ' of ', size(Posteriors,1), ' posteriors were skipped as they we invlaid - NaNs or 0 ##########################'
    if((1.e0_double*nPosteriorsSkipped)/size(Posteriors,1) > 0.1) STOP 'Combine_Posteriors - number of skipped posteriors too large, stopping!'

    !--Renormalise--!                                                                                                                                                                                             
    if(iDoRenormalise) then
       Renorm = 0.e0_double
       do j = 1, size(Combined_Posterior)-1
          Renorm = Renorm + 0.5e0_double*(Combined_Posterior(j) + Combined_Posterior(j+1))*(PosteriorGrid(j+1)-PosteriorGrid(j))
       end do
       if(Renorm <= 0.e0_double) then
          print *, 'Renormalisation:', Renorm
          STOP 'Combine_Posteriors - Invalid Renormalisation for combined Posterior'
       end if
       Combined_Posterior(:) = Combined_Posterior(:)/Renorm
    end if

    if(any(isNAN(Combined_Posterior(:)))) then
       print *, 'Any NaNs in Renormalised Combined Posterior?', any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
       print *, 'Stopping'
       STOP
    END if

  end subroutine Combine_Posteriors_Vector


end module Bayesian_Posterior_Evaluation
