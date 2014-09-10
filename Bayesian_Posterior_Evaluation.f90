module Bayesian_Posterior_Evaluation
  !--Contains the methods to evaluate the Bayesian Posterior for an input of free parameters
  use Param_Types

  !--THIS HAS NOT BEEN IMPLEMENTED
  
  implicit none
  
  !--Function Overload-------------------------------------
  INTERFACE Likelihood_Evaluation_atVirialRadius
     module procedure Likelihood_Evaluation_atVirialRadius_SingleCluster, Likelihood_Evaluation_atVirialRadius_MultipleCluster
  END INTERFACE Likelihood_Evaluation_atVirialRadius

  
contains  

  
  !--Evalution of the posterior on a single virial radius grid point.
  real(double) function Likelihood_Evaluation_atVirialRadius_SingleCluster(Alpha, Method, Profile, Cluster_Position, Lens_Redshift, Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, Sigma_Crit)
    real(double),intent(in):: Alpha !-could easily be extended to multiple alpha
    !--Mass Model Declarations --> these may need editing to pass in multiple values
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
    real(double):: tCluster_Position(1,2)

    tCluster_Position(1,:) = Cluster_Position

    Likelihood_Evaluation_atVirialRadius_SingleCluster = Likelihood_Evaluation_atVirialRadius_MultipleCluster((/Alpha/), Method, Profile, tCluster_Position, (/Lens_Redshift/), Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, Sigma_Crit)

  end function Likelihood_Evaluation_atVirialRadius_SingleCluster

  real(double) function Likelihood_Evaluation_atVirialRadius_MultipleCluster(Alpha, Method, Profile, Cluster_Position, Lens_Redshift, Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, Sigma_Crit)
    use Cosmology, only: angular_diameter_distance_fromRedshift
    use Mass_Profiles; use Distributions, only: CH08_Redshift_Distribution_Scalar
    use Interpolaters, only: Linear_Interp; use Integration
    real(double),intent(in):: Alpha(:) !-could easily be extended to multiple alpha
    !--Mass Model Declarations --> these may need editing to pass in multiple values
    integer:: Profile
    real(double), intent(in)::  Lens_Redshift(:), Cluster_Position(:,:) !-Cluster, RA/Dec-!
    !--Galaxy Data Input
    real(double), intent(in):: Theta, Magnitude, Source_Redshift, Position(2) !-RA,Dec-!
    !--Supplementary
    integer:: Method
    real(double), intent(in):: Pr_MagGrid(:), Pr_SizeGrid(:), SM_Prior(:,:), M_Prior(:)
    real(double), intent(in):: Size_Limits(:), Magnitude_Limits(:)
    real(double), intent(in):: MagnificationGrid(:), M_Pr_Renormalisation(:), SM_Pr_Renormalisation(:)
    real(double), intent(in):: RedshiftGrid(:), Sigma_Crit(:)
    
    !--Internal Declarations
    real(double), allocatable:: Posterior_perGalaxy_Redshift(:)
    real(double):: Effective_Magnification
    logical:: Known_Redshift
    integer:: z, m, j, l
    real(double):: Galaxy_Sigma_Critical
    real(double):: RedshiftPDF, Renorm
    real(double), dimension(size(Alpha)):: Distance_From_Mass_Center
    real(double):: D_l

    real(double), dimension(2):: Renormalisation_Size_Limits, Renormalisation_Magnitude_Limits
    real(double), dimension(size(SM_Prior,1), size(SM_Prior,2)):: Kappa_Renormalised_Prior, Size_Only_Mag_Prior
    real(double), dimension(size(M_Prior)):: Kappa_Renormalised_MagPrior
    real(double),allocatable:: Size_Only_Prior(:)

    integer::nZ

    !--Internal Decalrations that could eventually be set/passed in
    logical:: Enforce_Weak_Lensing = .false., Cuts_Renormalise_Likelihood = .true.

    !--Error Catching on Input---!
    if(size(Alpha) /= size(Distance_From_Mass_Center) .or. (size(Alpha) /= size(Lens_Redshift))) STOP 'Likelihood_Evaluation_atVirialRadius - Error on input of Mass Profile variables'
    
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

    allocate(Posterior_perGalaxy_Redshift(nZ)); Posterior_perGalaxy_Redshift = 0.e0_double
    
    do z = 1, nZ
       !--Set Simga Critical for that galaxy
       if(Known_Redshift) then
          Galaxy_Sigma_Critical = Linear_Interp(Source_Redshift, RedshiftGrid, Sigma_Crit)
       else
          Galaxy_Sigma_Critical = Sigma_Crit(z)
       END if
       if(isNaN(Galaxy_Sigma_Critical) .or. (Galaxy_Sigma_Critical > huge(1.e0_double))) STOP 'DM_Profile_Variable_Posterior - FATAL - Galaxy Sigma Critical not set correctly'
       
       !~~~~Ths section will need editing to account for multiple clusters in flat and SIS case
       select case(Profile)
       case(1) !-Flat-!                                                                                                                                                                                      
          STOP 'The Strong Lensing approximation has not been implemented for this profile yet'
          Effective_Magnification = 1.e0_double+2.e0_double*(Alpha(1)/Galaxy_Sigma_Critical)
       case(2) !-SIS-!                                                                                                                                                                                       
          STOP 'The Strong Lensing approximation has not been implemented for this profile yet'
          Effective_Magnification = 1.e0_double+2.e0_double*(SMD_SIS(Alpha(1), Distance_From_Mass_Center(1))/(Galaxy_Sigma_Critical*1.e18_double))
       case(3) !-NFW-!                                                                                                                                                                                       
          if(Enforce_Weak_Lensing) then
             STOP 'Weak Lensing disabled for NFW for now, until the equivalent for multiple clusters is coded [dont worry, its actually easy]'
             Effective_Magnification = 1.e0_double + 2.e0_double*(SMD_NFW(Distance_From_Mass_Center(1), Lens_Redshift(1), Alpha(1))/(Galaxy_Sigma_Critical*1.e18_double))
          else
!             Effective_Magnification = Total_MagnificationFactor_MultipleClusters(3, Position, Source_Redshift, Cluster_Position, Alpha, Lens_Redshift)
              Effective_Magnification = Magnification_Factor(3, Distance_From_Mass_Center(1), Alpha(1),  Lens_Redshift(1), Galaxy_Sigma_Critical*1.e18_double)
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
       
       select case(Method)
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
!!$                   allocate(Mag_Only_Prior(2)); Mag_Only_Prior = 0.e0_double                                                                                                                                    
!!$                   !--Find point where de-lensed magnitude lie on intrinsic distribution - could also use locate.                                                                                               
!!$                   do j = 1, size(Pr_MagGrid)-1                                                                                                                                                               
!!$                      if( ( Pr_MagGrid(j)<=  Magnitude + 2.5e0_double*dlog10(Effective_Magnification) ) .and. (  Pr_MagGrid(j+1) > Magnitude + 2.5e0_double*dlog10(Effective_Magnification)) ) then
!!$                         Mag_Only_Prior(1) = Integrate(PR_SizeGrid, Kappa_Renormalised_Prior(j,:), 2, lim = Renormalisation_Size_Limits)                                                                      
!!$                         Mag_Only_Prior(2) = Integrate(PR_SizeGrid, Kappa_Renormalised_Prior(j+1,:), 2, lim = Renormalisation_Size_Limits)                                                                    
!!$                         exit                                                                                                                                                                                   
!!$                      end if                                                                                                                                                                                    
!!$                   end do                                                                                                                                                                                       
!!$                   Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(Magnitude + 2.5e0_double*dlog10(Effective_Magnification), Pr_MagGrid(j:j+1), Mag_Only_Prior, ExValue =  1.e-100_double)*RedshiftPDF(z)                                                                                                                                                                                                               
!!$                   deallocate(Mag_Only_Prior)                                                                                                                                                                   
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
       Likelihood_Evaluation_atVirialRadius_MultipleCluster = Posterior_perGalaxy_Redshift(1)
    else
       !--Integrate over the redshift Information--!                                                                                                                                                         
       Likelihood_Evaluation_atVirialRadius_MultipleCluster = TrapInt(RedshiftGrid, Posterior_perGalaxy_Redshift(:))
    end if
    deallocate(Posterior_perGalaxy_Redshift)
    
    !--Deallocate Internal Declarations
    
  end function Likelihood_Evaluation_atVirialRadius_MultipleCluster

end module Bayesian_Posterior_Evaluation
