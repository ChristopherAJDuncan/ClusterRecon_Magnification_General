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
  real(double),allocatable:: Core_Cut_Position(:,:), Core_Cut_Radius(:)

  logical:: use_KDE_Smoothed_Distributions = .true., KDE_onTheFly = .false., allow_KDE_Extrapolation = .false.
  logical::use_lnSize_Prior = .false.
  logical:: Cuts_Renormalise_Likelihood = .true.
  real(double),dimension(2):: Survey_Magnitude_Limits = (/23.e0_double, 27.5e0_double/), Survey_Size_Limits = (/0.e0_double, 100.e0_double/)    
  real(double),dimension(2):: Prior_Magnitude_Limits = (/23.e0_double, 27.5e0_double/), Prior_Size_Limits = (/0.0e0_double, 100.e0_double/) !-3.3
  real(double),parameter:: Lower_Redshift_Cut = 0.21

  !--Parallisation Decalaration (To Be Passed In)
  integer:: nOMPThread = 1

  !----MCMC Declarations-----!
  integer:: nBurnin = 50
  integer:: nChains = 6, nChainOut = 6, nMinChain = 1050
  logical:: tune_MCMC = .false.
  logical:: allow_Automatic_Exit_MCMC = .false.
  logical:: output_Burnin = .true.
  integer:: fit_Parameter(3) = (/1,0,0/) !-r_200, concentration, position
  real(double):: centroid_Prior_Width = 0.03e0_double
  real(double), dimension(4):: MCMC_Proposal_Width = (/0.04e0_double, 0.0e0_double, 0.0015e0_double, 0.0015e0_double/) !-r200, c, RA, Dec-!

  !--Overload function for Combine Posteriors
!!$  interface Combine_Posteriors
!!$     module procedure Combine_Posteriors_Scalar, Combine_Posteriors_Vector
!!$  end interface Combine_Posteriors


contains

!!$  subroutine Combine_Posteriors_Scalar(GridValue, Posteriors,  Combine_by_ln, Renormalise, Return_lnP, Combined_Posterior)
!!$    !--Combines the posterior on a single grid value (free parameter alpha), using a call to the "normal" (vector) combined posterior subroutine
!!$     real(double), intent(in):: GridValue,Posteriors(:) !-Galaxy-!  
!!$     real(double), intent(out):: Combined_Posterior
!!$     logical, intent(in):: Combine_by_ln, Renormalise, Return_lnP
!!$
!!$     !--Internal Declarations
!!$     real(double),dimension(1):: tGrid, tCombinedPosterior
!!$     real(double), dimension(size(Posteriors),1):: tPosteriors
!!$
!!$     !--Set up internals
!!$     tGrid = GridValue
!!$     tPosteriors(:,1) = Posteriors
!!$
!!$     call Combine_Posteriors(tGrid, tPosteriors, Combine_by_ln, .false., Return_lnP, tCombinedPosterior)
!!$     Combined_Posterior = tCombinedPosterior(1)
!!$
!!$   end subroutine Combine_Posteriors_Scalar
!!$
!!$  subroutine Combine_Posteriors_Vector(PosteriorGrid, Posteriors, Combine_by_ln, Renormalise, Return_lnP, Combined_Posterior)
!!$    !--Combines posteriors by looping over the first dimension--!
!!$    real(double), intent(in):: PosteriorGrid(:),Posteriors(:,:) !-Galaxy, Grid/Value-!
!!$    real(double), intent(out):: Combined_Posterior(:)
!!$    logical, intent(in):: Combine_by_ln, Renormalise, Return_lnP
!!$
!!$    integer::c, j
!!$    real(double)::Renorm, Combination_Normalisation
!!$    logical::iDoRenormalise
!!$
!!$    integer::nPosteriorsSkipped
!!$
!!$    INTERFACE
!!$       subroutine Combine_Posteriors_Vector(PosteriorGrid, Posteriors, Combine_by_ln, Renormalise, Return_lnP, Combined_Posterior)
!!$         use Param_Types
!!$         real(double), intent(in):: PosteriorGrid(:), Posteriors(:,:)
!!$         real(double), intent(out):: Combined_Posterior(:)
!!$         logical, intent(in):: Combine_by_ln, Renormalise, Return_lnP
!!$       END subroutine Combine_Posteriors_Vector
!!$    END INTERFACE
!!$
!!$    !--Set renormalisation by input. If a single alpha value entered, then do not renormalise
!!$    iDoRenormalise = Renormalise
!!$    if(size(PosteriorGrid) == 1) iDoRenormalise = .false.
!!$
!!$    nPosteriorsSkipped = 0
!!$    if(Combine_by_ln == .false.) then
!!$       Combination_Normalisation = 1.e0_double/maxval(Posteriors)!or 1.e0_double/(0.5e0_double*maxval(Posteriors(c,:))) within loop
!!$       Combined_Posterior = 1.e0_double
!!$       do c = 1, size(Posteriors,1) !-Loop over galaxies-!   
!!$          !--Skip if zero (lnP not defined then) or NaN
!!$          if(all(Posteriors(c,:) == 0.e0_double) .or. all(isNaN(Posteriors(c,:)))) then
!!$             nPosteriorsSkipped = nPosteriorsSkipped + 1
!!$             cycle
!!$          end if
!!$
!!$          !--Skip when `renormalised' posteriors are invalid (zero/negative)
!!$          if(all(Posteriors(c,:)*Combination_Normalisation == 0.e0_double)) then
!!$             print *, 'Invalid Posterior for galaxy:', c, ' (==0) press [ENTER] to output an stop..'; READ(*,*)
!!$             print *, Posteriors(c,:)
!!$             STOP
!!$          end if
!!$          if(any(Posteriors(c,:)*Combination_Normalisation < 0.e0_double)) then
!!$             print *, 'Invalid Posterior for galaxy:', c, ' (<0) presS [ENTER] to output an stop..'; READ(*,*)
!!$             print *, Posteriors(c,:)
!!$             STOP
!!$          end if
!!$            
!!$          Combined_Posterior(:) = Combined_Posterior(:)*(Posteriors(c,:)*Combination_Normalisation)
!!$
!!$          if(all(Combined_Posterior == 0.e0_double)) then
!!$             print *, 'Invalid CPosterior for galaxy:', c, ' press (==0) [ENTER] to output an stop..'; READ(*,*)
!!$             print *, Combination_Normalisation
!!$             read(*,*)
!!$             print *, Combined_Posterior(:)
!!$             read(*,*)
!!$             print *, Posteriors(c,:)
!!$             STOP
!!$          end if
!!$          if(anY(Combined_Posterior < 0.e0_double)) then
!!$             print *, 'Invalid CPosterior for galaxy:', c, ' press (<0) [ENTER] to output an stop..'; READ(*,*)
!!$             print *, Posteriors(c,:)
!!$             STOP
!!$          end if
!!$          
!!$       end do
!!$    else
!!$       if(return_lnP) iDoRenormalise = .false.
!!$!       print *, 'Combining using logs'
!!$       !-Set lnP to a large negative value as default, equivalent to P ~ 0
!!$       Combined_Posterior = -100_double
!!$       Combination_Normalisation = size(Posteriors,1)
!!$       do c = 1, size(Posteriors,1) !-Loop over galaxies-!
!!$          !-Error Catching--!
!!$          if(all(Posteriors(c,:) == 1.e-75_double) .or. all(isNaN(Posteriors(c,:)))) then
!!$             nPosteriorsSkipped = nPosteriorsSkipped + 1
!!$             cycle
!!$          end if
!!$          !--Sum log posteriors--!
!!$!!!$          print *, 'CombinePosterior, loop:', c
!!$!!!$          print *, Combined_Posterior, dlog(Posteriors(c,:))
!!$
!!$          where(Posteriors(c,:) == 0.e0_double)
!!$             Combined_Posterior = Combined_Posterior - 100.e0_double
!!$          elsewhere
!!$             Combined_Posterior = Combined_Posterior + dlog(Posteriors(c,:)) + 1.e0_double
!!$          end where
!!$
!!$
!!$          if(any(isNAN(Combined_Posterior(:)))) then
!!$             print *, 'Any NaNs in Combined Posterior?, galaxy:',c, any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
!!$             STOP
!!$          end if
!!$       end do
!!$       !--Convert to PDF, not ln(PDF)--!
!!$       if(Return_lnP) then
!!$          return
!!$       else
!!$          if(size(Combined_Posterior)/= 1) then
!!$             Combined_Posterior = dexp(Combined_Posterior - maxval(Combined_Posterior))
!!$          else
!!$             Combined_Posterior = dexp(Combined_Posterior)
!!$          end if
!!$       end if
!!$    end if
!!$ 
!!$    
!!$
!!$    if(any(isNAN(Combined_Posterior(:)))) then
!!$       print *, 'Any NaNs in Combined Posterior?', any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
!!$       print *, 'Stopping'
!!$       STOP
!!$    END if
!!$
!!$    if(nPosteriorsSkipped > 0) write(*,'(A,I4,A,I4,A)') '################### Combine_Posteriors - ', nPosteriorsSkipped, ' of ', size(Posteriors,1), ' posteriors were skipped as they we invlaid - NaNs or 0 ###########################'
!!$    if((1.e0_double*nPosteriorsSkipped)/size(Posteriors,1) > 0.1) STOP 'Combine_Posteriors - number of skipped posteriors too large, stopping!'
!!$
!!$    !--Renormalise--!
!!$    if(iDoRenormalise) then
!!$       Renorm = 0.e0_double
!!$       do j = 1, size(Combined_Posterior)-1
!!$          Renorm = Renorm + 0.5e0_double*(Combined_Posterior(j) + Combined_Posterior(j+1))*(PosteriorGrid(j+1)-PosteriorGrid(j))
!!$       end do
!!$       if(Renorm <= 0.e0_double) then
!!$          print *, 'Renormalisation:', Renorm
!!$          STOP 'Combine_Posteriors - Invalid Renormalisation for combined Posterior'
!!$       end if
!!$       Combined_Posterior(:) = Combined_Posterior(:)/Renorm
!!$    end if
!!$
!!$    if(any(isNAN(Combined_Posterior(:)))) then
!!$       print *, 'Any NaNs in Renormalised Combined Posterior?', any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
!!$       print *, 'Stopping'
!!$       STOP
!!$    END if
!!$
!!$  end subroutine Combine_Posteriors_Vector

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

    if(any(isNaN(PosteriorGrid)) .or. any(PosteriorGrid > huge(1.e0_double))) STOP 'Posterior_Statistics: NaNs or Inf in Posterior Grid'

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
    use Mass_Profiles, only: Halo_Mass, virial_Radius_from_ProfileFreeParameter; use Integration, only: Integrate; use gridintervals, only: equalscale
    
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


  function Select_Source_Sample(Cat, Ap_Pos, Ap_Radius)
    !--Ap Radius must be in degrees
    type(Catalogue), intent(in):: Cat
    real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:)
    type(Catalogue):: Select_Source_Sample

    integer:: C
    type(Catalogue):: tCat, tSource_Catalogue

    !--Identify the source sample as the combination of sources in each aperture
    tCat = Cat
    do C = 1, size(Ap_Pos,1)
       !          print *, 'Cutting on a core radius of:', Core_Cut_Radius(Group_Index(C)), ' arcminutes for aperture:', Group_Index(C)
       call Identify_Galaxys_in_Circular_Aperture(tCat, Ap_Pos(C,:), Ap_Radius(C), TSource_Catalogue)!, Core_Radius = Core_Cut_Radius(Group_Index(C))/60.e0_double)
       call Concatonate_Catalogues(Select_Source_Sample, TSource_Catalogue) 
       call Catalogue_Destruct(TSource_Catalogue)
       !--Mask Source Galaxies for that aperture to ensure no double counting
       call Mask_Circular_Aperture(tCat, Ap_Pos(C,:),  Ap_Radius(C))          
    end do
    call Catalogue_Destruct(tCat)
    
  end function Select_Source_Sample

  !----------------------------------------------------------POSTERIOR PRODUCTION ROUTINES-------------------------------------------------------------------------------------------------------------!

  subroutine DM_Profile_Fitting_Simultaneous_MCMC(Cat, Ap_Pos, Ap_Radius, Marginalised_Posteriors, Distribution_Directory, reproduce_Prior, Fit_Group, Output_Prefix, Blank_Field_Catalogue)
    use Bayesian_Posterior_Evaluation, only: lnLikelihood_Evaluation_atVirialRadius_perSourceSample, get_likelihood_evaluation_precursors
    use Interpolaters, only: Linear_Interp; use MCMC; use Statistics, only: create_Histogram; use Common_Functions, only: Pyth_distance_between_Points
    !-Routine that produces posteriors on dark matter profile free parameters, by simultaneously fitting to all clusters that belong to the same 'Fit_Group'.
    !--Fit_Group should be in assending order starting from 1, with no numerical gaps, and should hold a value for each aperture considered. The number of free parameters for each fit is therefore determined by the number of clusters considered in each group. 

    type(Catalogue), intent(in)::Cat
    real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:)
    real(double),intent(out),allocatable::Marginalised_Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-! 
    character(*), intent(in):: Distribution_Directory
    logical, intent(in):: reproduce_Prior
    type(Catalogue), intent(in),optional::Blank_Field_Catalogue
    integer:: Fit_Group(:)
    character(*),intent(in):: Output_Prefix

    type(Catalogue):: tSource_Catalogue, tCat
    type(Catalogue):: Group_Cat
    real(double), dimension(Size(Ap_Pos,1))::iAp_Radius
    integer, allocatable:: Group_Index(:)

    integer:: C,G,i,j

    character(500):: Group_Output_Prefix

    real(double),allocatable:: Source_Positions(:,:)

    !--Result Output Decalrations
    integer::nAlpha_Out = 500
    real(double):: AlphaLimit_Out(2) = (/0.05e0_double, 3.e0_double/)

    character(500):: Combined_Chain_Output, Convergence_Test_Output, Acceptance_Rate_Output
    character(25):: fmt, conv_fmt, acc_fmt
    character(100):: Filename


    !--Distribution Declarations
    real(double),allocatable:: Joint_Size_Magnitude_Distribution(:,:), Magnitude_Distribution(:)
    real(double),allocatable::SizeGrid(:), MagGrid(:)

    !--Temporary Allocations--!
    real(double), allocatable:: tSigma_Crit(:,:)

    !--Likelihood Evaluation Precursor Declarations
    real(double),dimension(:,:),allocatable:: Survey_Renormalised_Prior
    real(double),dimension(:),allocatable:: Survey_Renormalised_MagPrior
    real(double),allocatable:: RedshiftGrid(:)
    real(double),allocatable::Sigma_Crit(:,:)    
    real(double),allocatable:: MagnificationGrid(:), Renormalisation_by_Magnification(:), MagOnly_Renormalisation_by_Magnification(:)
    real(double)::Lens_Redshift = 0.165e0_double    

    !--MCMC routines
    type(Array_of_MCMCChains),allocatable::ChainArray(:)
    integer:: M, out, nParameter_per_Cluster
    real(double),allocatable:: ChainWidths(:)
    real(double),allocatable:: Parameter_Start_Limits(:,:)
    character(300),allocatable:: ChainsOutput_Filename(:)
    integer:: chainCount
    real(double),allocatable:: Likelihood(:,:)
    real(double),allocatable:: Acceptance_Rate(:)
    logical:: Chains_Converged
    real(double),allocatable:: tAlpha(:), tAp_Pos(:,:), centroid_Prior_Offset(:)
    logical:: Acceptance
    integer:: nConvergenceTestPoint = 10
    real(double),allocatable:: ConvergenceStatistic(:)
    integer:: nMaxChain = 1000000
    logical, dimension(:),allocatable:: ifree_Parameter !used as a handle to tell convergence routine which parameter to check

    !--Temporary Posterior construction
    real(double),allocatable:: tMarginalised_Posterior_Grid(:), tMarginalised_Posterior(:)
    real(double), allocatable:: Combined_Chain(:,:), Combined_Likelihood(:)

    !--Paralleisation Declarations
    real(double), dimension(:), allocatable:: tSizes, tMF606W, tRedshift
    integer:: nOpenThreads, OMP_GET_NUM_THREADS

    !--Testing Declarations:
    real(double):: Time1, Time2
    
    if(size(Ap_Radius)==1) then
       iAp_Radius = Ap_Radius(1)
    else
       iAp_Radius = Ap_Radius
    end if


    !--Error Catching on Fit Group
    if(size(Fit_Group) /= size(Ap_Pos,1)) STOP 'DM_Profile_Fitting_Simultaneous - Fit_Group entered is not of the correct size.'
    if(minval(Fit_Group) /= 1) STOP 'DM_Profile_Fitting_Simultaneous - Fit_Group entered does not satisfy minimum value conditions (should be 1)'

    !_____________________________________________________________________________Process of Posterior Evaluation_____________________________________________________________________________!
    !--Construct/Read in prior distribution
    !~~~Repeat for all groups in Fit_Group:
    !-Process: For each group identify source sample: Should be all galaxies in aperture radius entered. If apertures do not overlap, then output an error message but continue with evaluation
    !- Set up new 'propose' point on chain for all free parameters
    !- Evaluate Posterior on chain.
    !- Accept/Reject point on chain
    !- Output Joint Posterior for each group to file, taken as the histogram of each column in the chain
    !- Assign Marginalised Posterior output through linear interpolation
    !~~~~ Repeat for next group
    !_____________________________________________________________________________Process of Posterior Evaluation_____________________________________________________________________________!


    !--Get Prior Distribution
    !--If a catalogue for the blank field is passed in, then use this catalogue to get the intrinsic distributions--!
    if(reproduce_Prior) then
       if(present(Blank_Field_Catalogue) == .false.) STOP 'DM_Profile_Variable_Posteriors_CircularAperture - Blankf Field Catalogue must be entered to allow for the production of the prior on a grid'
       write(*,'(A)') 'Producing Distribution from Catalogue'
       call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, Distribution_Directory, Blank_Field_Catalogue)
    else
       write(*,'(A)') 'Reading in distribution from:', Distribution_Directory
       call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, Distribution_Directory)
       print *, 'Success:', Distribution_Directory
    end if

    !--Get Precursor
    call get_Likelihood_Evaluation_Precursors(Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, RedshiftGrid, Sigma_Crit, MagnificationGrid, Renormalisation_by_Magnification, MagOnly_Renormalisation_by_Magnification, SizeGrid, MagGrid, Magnitude_Distribution, Joint_Size_Magnitude_Distribution, Survey_Magnitude_Limits, Survey_Size_Limits, (/(Lens_Redshift,i=1,size(Fit_Group))/), Lower_Redshift_Cut, Output_Prefix)

    deallocate(Joint_Size_Magnitude_Distribution, Magnitude_Distribution)

    !--Set up output grid, on which the marginalised posteriors will be output, as the linear interpolation of the constructed marginalised posterior below
    allocate(Marginalised_Posteriors(size(Ap_Pos,1), 2, nAlpha_Out));
    do i =1, nAlpha_Out
       Marginalised_Posteriors(:,1,i) = (/(AlphaLimit_Out(1) + (i-1)*((AlphaLimit_Out(2)-AlphaLimit_Out(1))/(nAlpha_Out-1)),j=1,size(Marginalised_Posteriors,1))/)
    end do

    do G = 1, 100 !-Assume no more than 100 group will be present
       !--Find all clusters with that grouping
       !---Group Index contains the index of the clusters within that group, and can be used to access the correct Aperture Location and Radius. Size of Group_Index is the number of free parameters being used
       allocate(Group_Index(count(Fit_Group == G)));
       !--Exit when all groups exhausted
       if(size(Group_Index) == 0) then
          if(maxval(Fit_Group) > G) then
             deallocate(Group_Index)
             cycle
          else
             exit
          end if
       end if

       print *, ' '
       print *, 'Joint Fitting Cluster Group ', G, ' of ', maxval(Fit_Group),'....'
       
       !--Store the index of the clusters being fit, for easy identification in the following
       i = 0
       do C = 1, size(Fit_Group)
          if(Fit_Group(C) == G) then
             i = i + 1
             Group_Index(i) = C
          end if
       end do

       !--Identify the source sample as the combination of sources in each aperture
       tCat = Cat
       do C = 1, size(Group_Index)
!          print *, 'Cutting on a core radius of:', Core_Cut_Radius(Group_Index(C)), ' arcminutes for aperture:', Group_Index(C)
          call Identify_Galaxys_in_Circular_Aperture(tCat, Ap_Pos(Group_Index(C),:), iAp_Radius(Group_Index(C)), TSource_Catalogue)!, Core_Radius = Core_Cut_Radius(Group_Index(C))/60.e0_double)
          call Concatonate_Catalogues(Group_Cat, TSource_Catalogue) 
          call Catalogue_Destruct(TSource_Catalogue)
          !--Mask Source Galaxies for that aperture to ensure no double counting
          call Mask_Circular_Aperture(tCat, Ap_Pos(Group_Index(C),:),  iAp_Radius(Group_Index(C)))          
       end do
       call Catalogue_Destruct(tCat)

       print *, 'Source Catalogue contains:', size(Group_Cat%RA), ' of ', size(Cat%RA), 'galaxies'

       if(allocated(Source_Positions)) deallocate(Source_Positions)
       allocate(Source_Positions(size(Group_Cat%RA),2));
       Source_Positions(:,1) = Group_Cat%RA; Source_Positions(:,2) = Group_Cat%Dec

       !--Set up temporary arrays of Sigma_Critical (usually aperture position)
       allocate(tSigma_Crit(size(Group_Index),size(Sigma_Crit,2))); !tSigma_Crit(1,:) = Sigma_Crit(Group_index(1),:); tSigma_Crit(2,:) = Sigma_Crit(Group_index(2),:) 
       do C = 1, size(Group_index)
          tSigma_Crit(C,:) = Sigma_Crit(Group_Index(C),:)
       end do

       !--Check for overlap between clusters

       !--Set up original chains (nChain in total, used to calculate convergence)
       allocate(ChainArray(nChains))
       allocate(Acceptance_Rate(size(ChainArray))); Acceptance_Rate = 0

       !--Set output filename
       allocate(ChainsOutput_Filename(size(ChainArray)))
       write(ChainsOutput_Filename(1), '(I3)') G
       Combined_Chain_Output = trim(adjustl(Output_Prefix))//'Group'//trim(adjustl(ChainsOutput_Filename(1)))//'_MCMC_CombinedChain.dat'
       Convergence_Test_Output = trim(adjustl(Output_Prefix))//'Group'//trim(adjustl(ChainsOutput_Filename(1)))//'_MCMC_ConvergenceTest_R.dat'
       Acceptance_Rate_Output = trim(adjustl(Output_Prefix))//'Group'//trim(adjustl(ChainsOutput_Filename(1)))//'_MCMC_AcceptanceRate.dat'
       ChainsOutput_Filename = trim(adjustl(Output_Prefix))//'Group'//trim(adjustl(ChainsOutput_Filename(1)))//'_MCMC_Chain_'
       do M = 1, size(ChainArray)
          write(ChainsOutput_Filename(M),'(A,I1,A)') trim(adjustl(ChainsOutput_Filename(M))),M,'.dat'
       end do


       !-Set Parameter_Start_Limits, which also sets which parameters are being varied
       !--Chain is set up to include all parameters by default, however if they are not marginalised over then the proposal distribution has width zero in that parameters direction
       !--Chain Width specifies the typical width of the proposal distribution in that axis of parameter space
       !--ifree_Parameter labels whether all parameters in chain are free (true) or fixed (false). When fixed, chain width is set to zero, and further in code, starting point is set to the fixed value
       nParameter_per_Cluster = size(fit_Parameter)+1 !-Increment as position requires two parameters
       allocate(Parameter_Start_Limits(size(Group_Index)*nParameter_per_Cluster, 2)); Parameter_Start_Limits = 0.e0_double
       allocate(ChainWidths(size(Group_Index)*nParameter_per_Cluster)); ChainWidths = 0.e0_double
       allocate(ifree_Parameter(size(Group_Index)*nParameter_per_Cluster)); ifree_Parameter = .true.
       do C = 1, size(Group_Index)
          !--Set initial values
          Parameter_Start_Limits((C-1)*nParameter_per_Cluster+1,:) = (/0.05e0_double, 2.5e0_double/) !-r200
          Parameter_Start_Limits((C-1)*nParameter_per_Cluster+2,:) = (/0.0e0_double, 0.0e0_double/) !-Concentration
          Parameter_Start_Limits((C-1)*nParameter_per_Cluster+3,:) = (/maxval((/Ap_Pos(Group_Index(C),1)-centroid_Prior_Width, 148.7707e0_double/)), minval((/Ap_Pos(Group_Index(C),1)+centroid_Prior_Width, 149.3524e0_double/))/) !-RA
          Parameter_Start_Limits((C-1)*nParameter_per_Cluster+4,:) = (/maxval((/Ap_Pos(Group_Index(C),2)-centroid_Prior_Width, -10.291e0_double/)), minval((/Ap_Pos(Group_Index(C),2)+centroid_Prior_Width, -9.748e0_double/))/) !-Dec

          !-0.1 works well for 6 Cluster (r200 only).
          ChainWidths((C-1)*nParameter_per_Cluster+1) = MCMC_Proposal_Width(1) !-r200
          ChainWidths((C-1)*nParameter_per_Cluster+2) = MCMC_Proposal_Width(2) !-Concentration (not coded up yet)
          ChainWidths((C-1)*nParameter_per_Cluster+3) = MCMC_Proposal_Width(3) !-RA
          ChainWidths((C-1)*nParameter_per_Cluster+4) = MCMC_Proposal_Width(4) !-Dec

          !--Set values when not free
          if(fit_Parameter(1)==0) then !--r200
             ifree_Parameter((C-1)*nParameter_per_Cluster+1) = .false.
          end if
          if(fit_Parameter(2)==0) then !--Concentration
             ifree_Parameter((C-1)*nParameter_per_Cluster+2) = .false.
          end if
          if(fit_Parameter(3)==0) then !--Position
             ifree_Parameter((C-1)*nParameter_per_Cluster+3:(C-1)*nParameter_per_Cluster+4) = .false.
             ChainWidths((C-1)*nParameter_per_Cluster+3:(C-1)*nParameter_per_Cluster+4) = 0.e0_double
          end if
       end do


       !--Set starting link for chain, by randomly placing within set parameter limits
       do M = 1,size(ChainArray)
          call MCMC_StartChain_Random(Parameter_Start_Limits, ChainArray(M)%Chain)
          !--Start chain postions seperately, since prior is defined on an aperture (Reset over those positions set in previous step)
          print *, 'WARNING, starting position set very small indeed to test acceptance rate, THIS NEEDS FIXED'
          if(fit_Parameter(3) == 1) then
             do C = 1, size(Group_Index)
                ChainArray(M)%Chain(1,(C-1)*nParameter_per_Cluster+3:(C-1)*nParameter_per_Cluster+4) = Start_MCMC_Chain_TopHatAperture(Ap_Pos(Group_Index(C),:), minval( (/0.25e0_double,0.25e0_double*centroid_Prior_Width/))) !-maximum width of 0.25 ensures that arbitrarily large  centroid_Prior can be used without placing apeture too far from center, prefactor on centroid prior notes that usually the entered position is well known
                !--Ensure that centroid is not positioned outside survey limits
                !-(Ignored for now, as probably shouldn't affec the result, however this would need relaxed if prior on survey RA/Dec limits imposed)
                if(Pyth_distance_between_Points(ChainArray(M)%Chain(1,(C-1)*nParameter_per_Cluster+3:(C-1)*nParameter_per_Cluster+4), Ap_Pos(Group_Index(C),:)) > centroid_Prior_Width) STOP 'Error in setting initial starting position for cluster'
             end do
          end if
       end do

       deallocate(Parameter_Start_Limits)

       !--If a certain parameter is not being fit, set to input values (never true for r200, but possibly true for position and concentration)
       if(fit_Parameter(1) == 0) then
          !-r200
          print *, 'Default r200 value has been set by hand, THIS SHOULD BE EDITED'
          do C = 1, size(Group_Index)
             ChainWidths((C-1)*nParameter_per_Cluster+1) = 0.e0_double
             do M = 1, size(ChainArray)
                ChainArray(M)%Chain(1,(C-1)*nParameter_per_Cluster+1) = 1.2e0_double
             END do
          END do
          STOP 'Whilst r200 may not necessarily need to be fit, I cant see why one wouldnt. Therefore, I am stopping'
       end if
       if(fit_Parameter(2) == 0) then
          !-Concentration
       end if
       if(fit_Parameter(3) == 0) then
          !--Position
          do C = 1, size(Group_Index)
             do M = 1, size(ChainArray)
                ChainArray(M)%Chain(1,(C-1)*nParameter_per_Cluster+3:(C-1)*nParameter_per_Cluster+4) = Ap_pos(Group_Index(C),:)
             end do
          end do
       end if

       !--Testing
       print *, '------------------------------------------------------------------'
       print *, 'Chains starting at:'
       do M = 1, size(ChainArray)
          print *, M, ':', ChainArray(M)%Chain(1,:)
       end do
       print *, '------------------------------------------------------------------'

       print *, 'Chain width used:'
       print *, ChainWidths
       print *, ' '

       !_____________________________________________________Proceed with chain_______________________________________________________________!

       !----Open output files
       !------Ensure no attempt is made to output more chains than exists
       nChainOut = minval((/nChains, nChainOut/))
       do out = 1, nChainOut
          open(unit = 30+out, file = ChainsOutput_Filename(out))
          write(*,'(2(A))') '---Outputting Chain to file: ', ChainsOutput_Filename(out)
          !--Header--!
          write(30+out,'(A,I3,A)') '# Output is: ChainLink (r200, c, RA, DEC)x', size(Group_Index), ' clusters. Final column is ln-likelihood (not renormalised).'
          if(output_Burnin) write(30+out, '(A)') '# Burnin is included.'
       end do
       write(fmt, '(I3)') size(ChainArray(1)%Chain,2) + 1
       fmt = '(I7,x,'//trim(adjustl(fmt))//'(e16.9,x))'
       
       !-Convergence_Test
       open(unit = 28, file = Convergence_Test_Output)
       !--Header--!
       write(28, '(A)') '# Included is R (Gelman-Rubin) across parameters'
       write(conv_fmt, '(I3)') size(ChainArray(1)%Chain,2)
       conv_fmt = '('//trim(adjustl(conv_fmt))//'(e12.5,x))'

       open(unit = 27, file = Acceptance_Rate_Output)
       write(27, '(A)') '# Chain Link; Acceptance Rate by chain'
       write(acc_fmt,'(I4)') size(Acceptance_Rate)
       acc_fmt = '(I5,x,'//trim(adjustl(acc_fmt))//'(e10.3,x))'

       !--Set minimum chain length so that it is at least 20 times the burnin
       nMinChain = maxval((/nMinChain, 20*nBurnin/))

       write(*,'(A,x,I4,x,A)') 'A minimum of:', nMinChain, ' chain links will be evaluated.'

       !--Set Paralellisation definitions
       if(nOMPThread <= 0) then
          nOpenThreads = size(ChainArray)
          !nOpenThreads = minval((/OMP_GET_NUM_THREADS(), size(ChainArray)/))
       else
          nOpenThreads = nOMPThread
       end if
       write(*,'(A,I3,A)') '!---- Running Parallelisation over: ', nOpenThreads, ' threads.'

       !--Initialise counters
       chainCount = 0; allocate(Likelihood(size(ChainArray),nMaxChain)); Acceptance_Rate = 0.e0_double
       chains_Converged = .false.
       do
          !--Exit Conditions
          if(Chains_Converged .and. chainCount >= nMinChain .and. allow_Automatic_Exit_MCMC) exit
          chainCount = chainCount + 1


          if(chainCount > nMaxChain) then
             print *, 'WARNING: Reached the maximum number of allowed chains without convergence, exiting ~~~~~~~~~'
             exit
          end if
          
          !--Reset Acceptance Rate after burnin
          if(chainCount == nBurnin) Acceptance_Rate = 0.e0_double

          !_______________________Propose New Point for each chain_____________________________________________________!
          if(chainCount > 1) then
             do M = 1, size(ChainArray)
                !--Take new point (skip on first run to allow for evaluation of first point
                call MCMC_Propose_Point_TopHat(ChainArray(M)%Chain, ChainWidths)
             end do
          end if

          !______________________EVALUATE LIKELIHOOD AT NEW CHAIN POINT________________________________________________!
          !--Set up Temporary Source Sample array (Parallelisation does not like derived types)
          allocate(tMF606W(size(Group_Cat%MF606W))); tMF606W = Group_Cat%MF606W
          allocate(tRedshift(size(Group_Cat%Redshift))); tRedshift = Group_Cat%Redshift
          allocate(tSizes(size(Group_Cat%Sizes))); tSizes = Group_Cat%Sizes

          call OMP_SET_NUM_THREADS(nOpenThreads)
          !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(M, tAlpha, tAp_Pos, centroid_Prior_Offset)
          !$OMP DO
          do M = 1, size(ChainArray)
             if(allocated(tAlpha)) deallocate(tAlpha); if(allocated(tAp_Pos)) deallocate(tAp_Pos)

             !--Set up temporary arrays to pass in free parameters
             allocate(tAlpha(size(Group_Index))); allocate(tAp_Pos(size(Group_Index),2))
             do C = 1, size(Group_Index)
                tAlpha(C) = ChainArray(M)%Chain(chainCount, (C-1)*nParameter_per_Cluster+1)
                tAp_Pos(C,:) = ChainArray(M)%Chain(chainCount, (C-1)*nParameter_per_Cluster+3:(C-1)*nParameter_per_Cluster+4)
             end do

             !__________________________________________Evaluate Posterior
             !--Priors on Parameters
             !---If using position, to avoid degeneracy need alpha1>alpha2>...., and then 1 not longer necessarily labels the 1st cluster, but the largest (identifiable by location)
             if(any(tAlpha <= 0.05)) then
                !--Put a low probability on such low values of alpha so they are never accepted
                Likelihood(M,chainCount) = -1.e0_ldp*huge(1.e0_ldp)!-100000
                cycle
             end if
             if(fit_Parameter(3) == 1) then
                !--If fitting centroid position....
                                
                !--Apply prior on position of centroid. Note: Care must be taken here when two clusters can share the same parameter space, as this will lead to degeneracies which should be accounted for by setting an ordering prior on alpha
                if(allocated(centroid_Prior_Offset)) deallocate(centroid_Prior_Offset)
                allocate(centroid_Prior_Offset(size(Group_Index))); centroid_Prior_Offset = 0.e0_double
                centroid_Prior_Offset = (/ (Pyth_distance_between_Points(tAp_Pos(C,:), Ap_Pos(Group_Index(C),:)), C = 1, size(Group_Index)) /)
                if(any(centroid_Prior_Offset > centroid_Prior_Width)) then
                   !--If any of the cntroids fall outside a circular top hat centered on the entered aperture position, with width given by centroid_Prior_Width...
                   Likelihood(M,chainCount) = -1.e0_ldp*huge(1.e0_ldp)
                   cycle
                end if
                
!!$                if(real_Ordered_Array_Decreasing(tAlpha) == .false.) then
!!$                   !--If fitting position, and Cluster N is larger than N-1 then set probability to zero
!!$                   !----This allows parameters to be fit over the whole field space, where labelling of cluster is not conserved (i.e. returned cluster 1 is the most massive)
!!$                   Likelihood(M,chainCount) = -1.e0_ldp*huge(1.e0_ldp)
!!$                   cycle
!!$                end if
             end if

             !--Likelihood evaluation provided priors have been passed
             Likelihood(M,chainCount) = lnLikelihood_Evaluation_atVirialRadius_perSourceSample(tAlpha, Posterior_Method, Surface_Mass_Profile, tAp_Pos, (/(Lens_Redshift, i = 1, size(tAlpha))/), tSizes, tMF606W, tRedshift, Source_Positions, MagGrid, SizeGrid, Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, Survey_Size_Limits, Survey_Magnitude_Limits, MagnificationGrid, MagOnly_Renormalisation_by_Magnification, Renormalisation_by_Magnification, RedshiftGrid, tSigma_Crit)

             if(allocated(centroid_Prior_Offset)) deallocate(centroid_Prior_Offset)
             if(allocated(tAlpha)) deallocate(tAlpha); if(allocated(tAp_Pos)) deallocate(tAp_Pos)
          end do
          !$OMP END PARALLEL
          deallocate(tMF606W, tRedshift, tSizes)

          !--Check that initial value of the chain has not been placed outside priors
          if(chainCount == 1 .and. any(Likelihood(:,chainCount) <= -(huge(1.e0_double)+100.e0_double))) then
             print *, 'MCMC - Possible error with initial positioning of chain, outside prior limits. FATAL.'
             do M = 1, size(Likelihood,1)
                if(Likelihood(M,1) <= (huge(1.e0_double)+100.e0_double)) then
                   print *, '!--Chain:', M, ChainArray(M)%Chain(1,:)
                   print *, ' '
                end if
             end do
             read(*,*)
          end if

          !___________________________________ACCEPT/REJECT NEW POINT ON EACH CHAIN_________________________________________________!
          do M = 1, size(ChainArray)
             !--Test for acceptance
             if(chainCount > 1) then
                call MCMC_accept_reject_Point_Metropolis(ChainArray(M)%Chain(chainCount-1:chainCount,:), Likelihood(M,chainCount-1:chainCount), acceptance, .true., Acceptance_Rate(M), chainCount-1)

                !--Test Behaviour of priors
                if(any( (/(ChainArray(M)%Chain(chainCount, (C-1)*nParameter_per_Cluster+1), C = 1, size(Group_Index))/) < 0.e0_double) .and. Likelihood(M,chainCount) > -100000) then
                   print *, 'Group, Chain, chainLink:', G, M, chainCount
                   STOP 'MCMC- Fatal- A negative Value of Alpha has been accepted!'
                end if
             end if

             if(chainCount > 1 .and. tune_MCMC) then
                write(*,'(A,x,I2,x,I5,x,e12.5,x,L)') 'Chain, ChainLink, Acceptance Rate: ', M, chainCount, Acceptance_Rate(M), acceptance
                if(M == size(ChainArray)) print *, ' '
             end if

          end do

          !--Output acceptance rate to file
          if(chainCount > 1) write(27, acc_fmt) chainCount, Acceptance_Rate

          !--Output Point to file
          if(output_Burnin .or. chainCount > nBurnin) then
             do out = 1, nChainOut
                write(30+out, fmt) chainCount, ChainArray(out)%Chain(chainCount, :), Likelihood(out,chainCount) 
             end do
          end if

          !--Cycle on first point in chain
          if(chainCount < 2) cycle

          !--Every nConvergenceTestPoint test for convergence across the chains
          if(chainCount > nBurnin .and. mod(chainCount-nBurnin, nConvergenceTestPoint) == 0) then
             Chains_Converged = convergence_Test_GelmanRubin(ChainArray, (/nBurnin/), ConvergenceStatistic, ifree_Parameter, 1.03e0_double)
             write(28,conv_fmt) ConvergenceStatistic
          end if
          if(allocated(ConvergenceStatistic)) deallocate(ConvergenceStatistic)

       end do !--End of chain loop
       print *, 'Finished Chain' !--Could output maximum value using maxloc, maybe average across chains
       do out = 1, nChainOut
          close(unit = 30+out)
       end do
       close(28)
       close(27)

       !--Combine all chains onto on large chain
       allocate(Combined_Chain(size(ChainArray)*size(ChainArray(1)%Chain(nBurnin+1:,1),1), size(ChainArray(1)%Chain,2))); Combined_Chain = dsqrt(-1.e0_double)
       allocate(Combined_Likelihood(size(Combined_Chain,1))); Combined_Likelihood = dsqrt(-1.e0_double)
       do M = 1, size(ChainArray)
          Combined_Chain((M-1)*size(ChainArray(1)%Chain(nBurnin+1:,1),1)+1:M*size(ChainArray(1)%Chain(nBurnin+1:,1),1),:) = ChainArray(M)%Chain(nBurnin+1:,:)
          Combined_Likelihood((M-1)*size(ChainArray(1)%Chain(nBurnin+1:,1),1)+1:M*size(ChainArray(1)%Chain(nBurnin+1:,1),1)) = Likelihood(M,nBurnin+1:)
       end do
       if(any(isNAN(Combined_Chain)) .or. any(isNAN(Combined_Likelihood))) STOP'FATAL: NaNs in combined chain (MCMC)'

       open(unit = 29, file = Combined_Chain_Output)
       write(29,'(A,I3,A)') '# Output is: ChainLink (r200, c, RA, DEC)x', size(Group_Index), ' clusters. Final column is likelihood.'
       write(29, '(A)') '# Burnin is NOT included.'
       do i = 1, size(Combined_Chain,1)
          write(29, fmt) i, Combined_Chain(i, :), Combined_Likelihood(i)
       end do
       close(29)

       !--Construct Marginalised Posteriors from one chain
       print *, 'I am only constructing Marginalised Posterior of Alpha......' !--For other variables, perhaps consider looping over P, and check that parent routine loops over this value
       !--Could use nChainLinks/n for second arguement (nBins), where the ''n'' here sets the number of chain links in each bin
       do C = 1, size(Group_Index)
          call create_Histogram(Combined_Chain(:,(C-1)*nParameter_per_Cluster+1), 40, tMarginalised_Posterior_Grid, tMarginalised_Posterior)
          Marginalised_Posteriors(Group_Index(C), 1+1, :) = Linear_Interp(Marginalised_Posteriors(Group_Index(C),1,:), tMarginalised_Posterior_Grid, tMarginalised_Posterior, ExValue =1.e-100_double)
          deallocate(tMarginalised_Posterior_Grid, tMarginalised_Posterior)
       end do

       deallocate(Combined_Chain, Combined_Likelihood)

       deallocate(ChainsOutput_Filename, Likelihood, Source_Positions, tSigma_Crit, Group_Index, ChainWidths)
       call Catalogue_Destruct(Group_Cat)

       do M = 1, size(ChainArray)
          deallocate(ChainArray(M)%Chain)
       end do
       deallocate(ChainArray)
       deallocate(Acceptance_Rate)

    end do !-End of group loop

    deallocate(Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, RedshiftGrid, Sigma_Crit, MagnificationGrid, MagOnly_Renormalisation_by_Magnification, Renormalisation_by_Magnification, MagGrid, SizeGrid)



  end subroutine DM_Profile_Fitting_Simultaneous_MCMC

  function Start_MCMC_Chain_TopHatAperture(Ap_Center, Width)
    use Common_Functions, only: return_Random_Set
    real(double), intent(in):: Ap_Center(:)
    real(double), intent(in):: Width

    real(double), dimension(size(Ap_Center)):: Start_MCMC_Chain_TopHatAperture
    real(double),dimension(size(Ap_Center)):: Ran
    real(double):: Sum2, Pen_Width
    integer::i

    Ran = return_Random_Set(size(Ran))
    !--Shift random onto interval [-1,1]
    Ran = 2.e0_double*(Ran-0.5e0_double)

    !--Randomly set 1st N-1 parameters
    Sum2 = 0.e0_double
    do i = 1, size(Ap_Center)-1
       Start_MCMC_Chain_TopHatAperture(i) = Ap_Center(i) + Ran(i)*Width
       Sum2 = Sum2 + (Start_MCMC_Chain_TopHatAperture(i)-Ap_Center(i))**2.e0_double
    end do

    !--Set width of penultimate co-ordinate as y = sqrt(Width - sum(x^2))
    Pen_Width = dsqrt(Width*Width - Sum2)
    i = size(Start_MCMC_Chain_TopHatAperture)
    Start_MCMC_Chain_TopHatAperture(i) = Ap_Center(i) + Ran(i)*Pen_Width

  end function Start_MCMC_Chain_TopHatAperture

  subroutine DM_Profile_Fitting_Simultaneous_2Cluster(Cat, Ap_Pos, Ap_Radius, Marginalised_Posteriors, Distribution_Directory, reproduce_Prior, Fit_Group, Alpha_Limit, Grid_Tolerance, Alpha_Tolerance, Output_Prefix, Blank_Field_Catalogue)
    use Bayesian_Posterior_Evaluation, only: lnLikelihood_Evaluation_atVirialRadius_perSourceSample, get_likelihood_evaluation_precursors, lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample
    use gridintervals, only: equalscale; use Integration, only: Integrate; use Interpolaters, only: Linear_Interp
    !-Routine that produces posteriors on dark matter profile free parameters, by simultaneously fitting to all clusters that belong to the same 'Fit_Group'.
    !--Fit_Group should be in assending order starting from 1, with no numerical gaps, and should hold a value for each aperture considered. The number of free parameters for each fit is therefore determined by the number of clusters considered in each group.
    !--Alpha_Limit should be a one-dimensional array, size of 2*nCluster. Each pair or elements labels the alpha limits on which the posterior should be evaluated. Outwith these limits, the posterior is set to zero, so these limits CAN BE THOUGHT OF AS THE LIMITS OF A FLAT PRIOR ON ALPHA.
    !--Grid_tolerance is a nCl size 1D array with sets the grid size on which the posterior is evaluated
    !--Alpha Tolerance sets the tolerance to which the maximum of the grid is found. If Alpha_Tolerance < Grid_Tolerance, then a recursive call evaluates the grid on a finer grid around the location of the maximum (either in terms of the marginalised posteriors, or the likelihood)
    !--First implementation will evaluate on a basic grid. Extensions will search for the mode value by: Evaluating on a finer grid around maxloc; Numerical Recipes Simplex method around mode;
    !---Further extension could implement MCMC if more than two free parameters are to be evaluated
    !-- Output is marginalised posterior interpolated for true output. Difficulty in passing out joint posterior results from the unknown size of the joint distribution (rank), and how many need passed out. Instead, they are written to file.

    type(Catalogue), intent(in)::Cat
    real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:)
    real(double),intent(out),allocatable::Marginalised_Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-! 
    character(*), intent(in):: Distribution_Directory
    logical, intent(in):: reproduce_Prior
    type(Catalogue), intent(in),optional::Blank_Field_Catalogue
    integer:: Fit_Group(:)
    real(double):: Alpha_Limit(:), Grid_tolerance(:), Alpha_Tolerance
    character(*),intent(in):: Output_Prefix

    
    real(double),dimension(2*size(Ap_Pos,1)):: iAlpha_Limit
    type(Catalogue):: tSource_Catalogue, tCat
    type(Catalogue):: Group_Cat
    real(double), dimension(Size(Ap_Pos,1))::iAp_Radius
    integer, allocatable:: Group_Index(:)

    integer:: C,G,i,j

    character(500):: Group_Output_Prefix

    real(double),allocatable:: Source_Positions(:,:)

    !--Result Output Decalrations
    integer::nAlpha_Out = 500
    real(double):: AlphaLimit_Out(2) = (/0.05e0_double, 3.e0_double/)

    character(15):: fmt
    character(100):: Filename

    !--Distribution Declarations
    real(double),allocatable:: Joint_Size_Magnitude_Distribution(:,:), Magnitude_Distribution(:)
    real(double),allocatable::SizeGrid(:), MagGrid(:)

    !--Posterior Grid Declarations
    integer, allocatable:: nAlpha(:)
    real(double),allocatable:: AlphaGrid1(:), AlphaGrid2(:)

    real(double),allocatable:: Likelihood(:,:)
    real(double),allocatable:: Marginalised_Cl1(:), Marginalised_Cl2(:)

    !--Temporary Allocations--!
    real(double), allocatable:: tSigma_Crit(:,:), tAp_Pos(:,:)
    !----As Part of OPENMP
    real(double), allocatable:: tMF606W(:), tRedshift(:), tSizes(:)

    !--Likelihood Evaluation Precursor Declarations
    real(double),dimension(:,:),allocatable:: Survey_Renormalised_Prior
    real(double),dimension(:),allocatable:: Survey_Renormalised_MagPrior
    real(double),allocatable:: RedshiftGrid(:)
    real(double),allocatable::Sigma_Crit(:,:)    
    real(double),allocatable:: MagnificationGrid(:), Renormalisation_by_Magnification(:), MagOnly_Renormalisation_by_Magnification(:)
    real(double)::Lens_Redshift = 0.165e0_double    

    !--Fine Grid declaration
    integer,dimension(2):: Max_Point
    real(double),allocatable:: Max_Marginalised_Posteriors(:,:,:)
    real(double), dimension(2*size(Ap_Pos,1)):: Max_Limit

    !--Single Cluster Fitting Declarations
    real(double),allocatable:: Posterior_Single(:,:)

    !-- Parallelisation
    integer:: OMP_GET_NUM_THREADS

    if(size(Ap_Radius)==1) then
       iAp_Radius = Ap_Radius(1)
    else
       iAp_Radius = Ap_Radius
    end if

    if(size(Alpha_Limit) == 2) then
       do i = 1, size(Ap_Pos,1)
          iAlpha_Limit(2*i-1:2*i) = Alpha_Limit
       end do
    elseif(size(Alpha_Limit) /= size(Ap_Pos,1)) then
       STOP '2Cluster - Alpha_limit Entered is not of the correct size'
    else
       iAlpha_Limit = Alpha_Limit
    end if

!DELETE    Output_Prefix = trim(adjustl(Bayesian_Routines_Output_Directory))

    !--Error Catching on Fit Group
    if(size(Fit_Group) /= size(Ap_Pos,1)) STOP 'DM_Profile_Fitting_Simultaneous - Fit_Group entered is not of the correct size.'
    if(minval(Fit_Group) /= 1) STOP 'DM_Profile_Fitting_Simultaneous - Fit_Group entered does not satisfy minimum value conditions (should be 1)'

    !_____________________________________________________________________________Process of Posterior Evaluation_____________________________________________________________________________!
    !--Construct/Read in prior distribution
    !~~~Repeat for all groups in Fit_Group
    !-Process: For each group identify source sample: Should be all galaxies in aperture radius entered. If apertures do not overlap, then output an error message but continue with evaluation
    !- Set up coarse grid for all free parameters
    !- Evaluate Posterior on coarse grid.
    !- Output Joint Posterior for each group to file
    !- Calculate and output marginalised posterior for each group.
    !- Assign Marginalised Posterior output through linear interpolation
    !~~~~ Repeat for next group
    !_____________________________________________________________________________Process of Posterior Evaluation_____________________________________________________________________________!


    !--Get Prior Distribution
    !--If a catalogue for the blank field is passed in, then use this catalogue to get the intrinsic distributions--!
    if(reproduce_Prior) then
       if(present(Blank_Field_Catalogue) == .false.) STOP 'DM_Profile_Variable_Posteriors_CircularAperture - Blank Field Catalogue must be entered to allow for the production of the prior on a grid'
       write(*,'(A)') 'Producing Distribution from Catalogue'
       call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, Distribution_Directory, Blank_Field_Catalogue)
    else
       write(*,'(A)') 'Reading in distribution from:', Distribution_Directory
       call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, Distribution_Directory)
       print *, 'Success:', Distribution_Directory
    end if

    !--Get Precursor
    call get_Likelihood_Evaluation_Precursors(Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, RedshiftGrid, Sigma_Crit, MagnificationGrid, Renormalisation_by_Magnification, MagOnly_Renormalisation_by_Magnification, SizeGrid, MagGrid, Magnitude_Distribution, Joint_Size_Magnitude_Distribution, Survey_Magnitude_Limits, Survey_Size_Limits, (/(Lens_Redshift,i=1,size(Fit_Group))/), Lower_Redshift_Cut, Output_Prefix)


    !--Set up output grid, on which the marginalised posteriors will be output, as the linear interpolation of the constructed marginalised posterior below
    allocate(Marginalised_Posteriors(size(Ap_Pos,1), 2, nAlpha_Out));
    do i =1, nAlpha_Out
       Marginalised_Posteriors(:,1,i) = (/(AlphaLimit_Out(1) + (i-1)*((AlphaLimit_Out(2)-AlphaLimit_Out(1))/(nAlpha_Out-1)),j=1,size(Marginalised_Posteriors,1))/)
    end do

    do G = 1, 100 !-Assume no more than 100 group will be present
       !--Find all clusters with that grouping
       !---Group Index contains the index of the clusters within that group, and can be used to access the correct Aperture Location and Radius. Size of Group_Index is the number of free parameters being used
       allocate(Group_Index(count(Fit_Group == G)));
       !--Exit when all groups exhausted
       if(size(Group_Index) == 0) then
          if(maxval(Fit_Group) > G) then
             deallocate(Group_Index)
             cycle
          else
             exit
          end if
       end if
       !--Exit if the grouping is too large. This is limited mainly by the inability to write a simultaneous fitting routine for a generally sized array - BOLLOCKS TO FORTRAN
       if(size(Group_Index) > 2) STOP 'DM_Profile_Fitting_Simultaneous - I have not been written to deal with the simultaneous fitting of more than two clusters at once, stopping.'

       print *, ' '
       print *, 'Joint Fitting Cluster Group ', G, ' of ', maxval(Fit_Group),'....'
       
       i = 0
       do C = 1, size(Fit_Group)
          if(Fit_Group(C) == G) then
             i = i + 1
             Group_Index(i) = C
          end if
       end do

       tCat = Cat
       do C = 1, size(Group_Index)
!          print *, 'Cutting on a core radius of:', Core_Cut_Radius(Group_Index(C)), ' arcminutes for aperture:', Group_Index(C)
          call Identify_Galaxys_in_Circular_Aperture(tCat, Ap_Pos(Group_Index(C),:), iAp_Radius(Group_Index(C)), TSource_Catalogue)!, Core_Radius = Core_Cut_Radius(Group_Index(C))/60.e0_double)
          call Concatonate_Catalogues(Group_Cat, TSource_Catalogue) 
          call Catalogue_Destruct(TSource_Catalogue)
          !--Mask Source Galaxies for that aperture to ensure no double counting
          call Mask_Circular_Aperture(tCat, Ap_Pos(Group_Index(C),:),  iAp_Radius(Group_Index(C)))          
       end do
       call Catalogue_Destruct(tCat)
       print *, 'Got source sample'

       if(allocated(Source_Positions)) deallocate(Source_Positions)
       allocate(Source_Positions(size(Group_Cat%RA),2));
       Source_Positions(:,1) = Group_Cat%RA; Source_Positions(:,2) = Group_Cat%Dec

       !--Check for overlap between clusters
       
       !--Set up a coarse grid on which the posterior will be evaluated
       !---For this point on the code is limited to two dimensions
       if(allocated(nAlpha)) deallocate(nAlpha)
       allocate(nAlpha(2)); nAlpha = 1
       do C= 1, size(Group_Index)
          nAlpha(C) = nint((IAlpha_Limit(2*Group_index(C)) - IAlpha_Limit(2*Group_index(C)-1))/Grid_tolerance(Group_index(C))+1)
       end do
       print *, 'Simultaneous fit for Clusters:', Group_Index, ' will take ', nAlpha(1)*nAlpha(2), 'grid points'
       call equalscale(IAlpha_Limit(2*Group_index(1)-1), IAlpha_Limit(2*Group_index(1)), nAlpha(1), AlphaGrid1)
       
       if(size(Group_Index) == 1) then
          !--Do single fitting of cluster
          print *, 'Just attempting a single fit...'
          allocate(Posterior_Single(2,size(AlphaGrid1))); Posterior_Single = dsqrt(-1.e0_double); Posterior_Single(1,:) = AlphaGrid1
          call DM_Profile_Variable_Posterior_SingleFit(Group_Cat, Surface_Mass_Profile, Lens_Redshift, Ap_Pos(Group_Index(1),:), Posterior_Single, Output_Prefix, .false., PriorMagGrid = MagGrid, PriorSizeGrid = SizeGrid, Prior = Joint_Size_Magnitude_Distribution, MagPrior = Magnitude_Distribution)
          Marginalised_Posteriors(Group_Index(1),2,:) = Linear_Interp(Marginalised_Posteriors(Group_Index(1),1,:), Posterior_Single(1,:), Posterior_Single(2,:), ExValue =1.e-100_double)
          deallocate(Posterior_Single, AlphaGrid1, Group_Index, Source_Positions)
          call Catalogue_Destruct(Group_Cat)
          print *, 'Got single fit of cluster'
          cycle
       end if

       call equalscale(IAlpha_Limit(2*Group_index(2)-1), IAlpha_Limit(2*Group_index(2)), nAlpha(2), AlphaGrid2)
       
       print *, ' '
       write(*,'(A,2(e9.2,x),A,I4,A,I1,A)') '--Alpha is being evaluated between:', minval(AlphaGrid1), maxval(AlphaGrid1), ', in ', nAlpha(1), ' grid points (Cluster ', Group_index(1), ')'
       write(*,'(A,2(e9.2,x),A,I4,A,I1,A)') '--Alpha is being evaluated between:', minval(AlphaGrid2), maxval(AlphaGrid2), ', in ', nAlpha(2), ' grid points (Cluster ', Group_index(2), ')'
       print *, ' '

       allocate(Likelihood(nAlpha(1), nAlpha(2))); Likelihood = 0.e0_double
       
       allocate(tSigma_Crit(2,size(Sigma_Crit,2))); tSigma_Crit(1,:) = Sigma_Crit(Group_index(1),:); tSigma_Crit(2,:) = Sigma_Crit(Group_index(2),:) 
       allocate(tAp_Pos(2, 2)); tAp_pos(1,:) = Ap_Pos(Group_index(1),:); tAp_Pos(2,:) = Ap_Pos(Group_index(2),:)

       !-----Evaluate Posterior - Parallelisation used, however can be switched off by setting the number of threads to 1. All varabiables passed to likelihood evaluation need to be intent(in) only, and be careful of library calls which depend on global declarations.
       allocate(tMF606W(size(Group_Cat%MF606W))); tMF606W = Group_Cat%MF606W
       allocate(tRedshift(size(Group_Cat%Redshift))); tRedshift = Group_Cat%Redshift
       allocate(tSizes(size(Group_Cat%Sizes))); tSizes = Group_Cat%Sizes
       !--Evaluate Posterior

       call OMP_SET_NUM_THREADS(nOMPThread)
!       call OMP_SET_NUM_THREADS(minval((/OMP_GET_NUM_THREADS(), nOMPThread/)))
       write(*,'(A,I2,A)') '---Paralellising over:', nOMPThread, ' threads.'
       !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(i,j)
       !--OMP DO parallelises the *outer loop only*- is nalpha(1) == nalpha(2), this could be collapsed using collapse(2)
       !$OMP DO
       do i = 1, nAlpha(1)
          do j = 1, nAlpha(2)
             Likelihood(i,j) = lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample((/AlphaGrid1(i),AlphaGrid2(j)/), Posterior_Method, Surface_Mass_Profile, tAp_Pos, (/Lens_Redshift,Lens_Redshift/), tSizes, tMF606W, tRedshift, Source_Positions, MagGrid, SizeGrid, Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, Survey_Size_Limits, Survey_Magnitude_Limits, MagnificationGrid, MagOnly_Renormalisation_by_Magnification, Renormalisation_by_Magnification, RedshiftGrid, tSigma_Crit)
             print *, 'Done:', i, j,  Likelihood(i,j)
          end do
       end do
       !$OMP END PARALLEL
       deallocate(tRedshift, tMF606W, tSizes)

       !--Convert from lnP to P
       Likelihood = dexp(Likelihood - maxval(Likelihood))
       
       !--Renormalise Posterior
       !-(Flat prior assumed at this point)
       !--This may not be strictly necessary
       Likelihood = Likelihood/Integrate(AlphaGrid1, AlphaGrid2, Likelihood, 2, lim1 = (/minval(AlphaGrid1), maxval(AlphaGrid1)/), lim2 = (/minval(AlphaGrid2), maxval(AlphaGrid2)/))
       
       deallocate(tSigma_Crit, Source_Positions, tAp_Pos)
       call Catalogue_Destruct(Group_Cat)
       
       !--Output Full Likelihood (or posterior, assuming flat prior bounded by grid)
       write(Filename, '(I1)') G
       open(unit = 21, file = trim(adjustl(Output_Prefix))//'_SimultaneousFitting_Likelihood_Group'//trim(Filename)//'.dat')
       write(21,'(2(A,x,I2,x))') '# Column 1 contains Virial Radius for cluster:', Group_Index(1), ' whilst Row 1 constains cluster:', Group_Index(2)

       !--Set Output Format
       write(fmt, '(I3)') size(Likelihood,2)+1
       fmt = '('//trim(fmt)//'(e12.5,x))'
       !--Write first row containing grid along second cluster
       write(21, fmt) 0.e0_double, AlphaGrid2
       do i = 1, size(Likelihood,1)
          write(21, fmt) AlphaGrid1(i), Likelihood(i,:)
       end do
       print *, '---Likelihood for group:', G, ' output to ', trim(adjustl(Output_Prefix))//'_SimultaneousFitting_Likelihood_Group'//trim(Filename)//'.dat'
       close(21)
       
       !--Get Marginalised Distributions
       !(At this point, a flat prior is assumed)
       allocate(Marginalised_Cl1(nAlpha(1))); Marginalised_Cl1 = 0.e0_double
       do i = 1, size(Marginalised_Cl1)
          Marginalised_Cl1(i) = Integrate(AlphaGrid2, Likelihood(i,:), 2, lim = (/minval(AlphaGrid2), maxval(AlphaGrid2)/))
       end do
       
       write(Filename, '(I1)') Group_Index(1)
       open(unit = 22, file = trim(adjustl(Output_Prefix))//'_SimultaneousFitting_MarginalisedLikelihood_Cluster'//trim(Filename)//'.dat')
       do i= 1, size(AlphaGrid1)
          write(22, '(2(e12.5,x))') AlphaGrid1(i), Marginalised_Cl1(i)
       end do
       print *, '--- Marginalised Likehood for Cluster:', Group_Index(1),' of group:', G, ' output to ', trim(adjustl(Output_Prefix))//'_SimultaneousFitting_MarginalisedLikelihood_Cluster'//trim(Filename)//'.dat'
       close(22)

       allocate(Marginalised_Cl2(nAlpha(2))); Marginalised_Cl2 = 0.e0_double
       do i = 1, size(Marginalised_Cl2)
          Marginalised_Cl2(i) = Integrate(AlphaGrid1, Likelihood(:,i), 2, lim = (/minval(AlphaGrid1), maxval(AlphaGrid1)/))
       end do

       write(Filename, '(I1)') Group_Index(2)
       open(unit = 22, file = trim(adjustl(Output_Prefix))//'_SimultaneousFitting_MarginalisedLikelihood_Cluster'//trim(Filename)//'.dat')
       do i= 1, size(AlphaGrid2)
          write(22, '(2(e12.5,x))') AlphaGrid2(i), Marginalised_Cl2(i)
       end do
       print *, '--- Marginalised Likehood for Cluster:', Group_Index(2),' of group:', G, ' output to ', trim(adjustl(Output_Prefix))//'_SimultaneousFitting_MarginalisedLikelihood_Cluster'//trim(Filename)//'.dat'
       close(22)
       
       !--Find Maximum of Likihood
       Max_Point=  maxloc(likelihood)
       !--Use of Group Tolerance ensures that the true maximum is bounded within the points considered
       Max_Limit(2*Group_Index(1)-1:2*Group_Index(1)) = (/AlphaGrid1(Max_Point(1))-Grid_Tolerance(Group_Index(1)),AlphaGrid1(Max_Point(1))+Grid_Tolerance(Group_Index(1))/)
       Max_Limit(2*Group_Index(2)-1:2*Group_Index(2)) = (/AlphaGrid1(Max_Point(2))-Grid_Tolerance(Group_Index(2)),AlphaGrid1(Max_Point(2))+Grid_Tolerance(Group_Index(2))/)
       
       !--Output Marginalised Posteriors
       !-Cluster, Grid/Value, Value
       Marginalised_Posteriors(Group_Index(1),2,:) = Linear_Interp(Marginalised_Posteriors(Group_Index(1),1,:), AlphaGrid1, Marginalised_Cl1, ExValue = 1.e-100_double)
       Marginalised_Posteriors(Group_Index(2),2,:) = Linear_Interp(Marginalised_Posteriors(Group_Index(2),1,:), AlphaGrid2, Marginalised_Cl2, ExValue = 1.e-100_double)

       !-Reset for next group
       deallocate(Marginalised_Cl1, Marginalised_Cl2, Likelihood, AlphaGrid1, AlphaGrid2)
       deallocate(Group_Index)
    end do

    deallocate(Joint_Size_Magnitude_Distribution, Magnitude_Distribution, SizeGrid, MagGrid)
    deallocate(Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, RedshiftGrid, Sigma_Crit, MagnificationGrid, MagOnly_Renormalisation_by_Magnification, Renormalisation_by_Magnification)
    
    !--Recursive call to find maximum of likelihood to a given tolerance
    !--Switched off for now to test ability of code to simultaneously fit two clusters
!!$    if(any(Grid_Tolerance > Alpha_Tolerance)) then
!!$       !--This could be edited so that the prior does not need to be read in again, or indeed that the precursors do not need to be read in
!!$       !--Pass new Alpha tolerance ten times greater to enusre recursive loop does not continue more than once
!!$       print *, 'Max limit is:', Max_Limit
!!$
!!$       call DM_Profile_Fitting_Simultaneous_2Cluster(Cat, Ap_Pos, Ap_Radius, Max_Marginalised_Posteriors, Distribution_Directory, .false., Fit_Group, Max_Limit, (/(Alpha_Tolerance, i =1, size(Ap_Pos,1))/), Alpha_Tolerance*10.e0_double, trim(adjustl(Output_Prefix))//'_FineGrid_AroundMaximum')
!!$       
!!$       deallocate(Max_Marginalised_Posteriors, Max_Limit)
!!$    end if
    
    
  end subroutine DM_Profile_Fitting_Simultaneous_2Cluster

  

  subroutine DM_Profile_Variable_Posteriors_CircularAperture(Cat, Ap_Pos, Ap_Radius, Posteriors, Distribution_Directory, reproduce_Prior, Blank_Field_Catalogue)
    use Statistics, only: mean_discrete, mode_distribution, variance_distribution; use Mass_profiles
    use Distributions; use Cosmology, only: angular_diameter_distance_fromRedshift; use Interpolaters, only: Linear_Interp; use gridintervals, only: equalscale
    !--Main routine that returns the Posteriors on DM free parameter (alpha) over all apertures.
    !--If the posterior construction is such that grid may vary, this interpolates onto a finer grid (linear interpolation). Thus returned posteriors on on the same grid
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

    !--Coarse Grid declarations
    integer::nCoarseGrid
    real(double):: Alpha_min, Alpha_Max
    real(double),allocatable:: Coarse_LikelihoodGrid(:)

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
!       call Maximise_Convergenceb_yShifts_inAperture(Cat, Blank_Field_Catalogue, Ap_Pos(i,:), Ap_Radius(i))

!       print *, 'Cutting on a core radius of:', Core_Cut_Radius(i), ' arcminutes for aperture:', i
       call Identify_Galaxys_in_Circular_Aperture(Cat, Ap_Pos(i,:), iAp_Radius(i), Ap_Cats(i))!, Core_Radius = Core_Cut_Radius(i)/60.e0_double)
       
       print *, '* Ap ', i ,' has mean (Size,Mag):',  mean_discrete(Ap_Cats(i)%Sizes), mean_discrete(Ap_Cats(i)%MF606W)
    end do

    do i =1, size(Ap_Cats)
       print *, 'Aperture:', i, ' contains:', size(Ap_Cats(I)%RA), ' galaxies'
    end do

    print *, count(Cat%Redshift >= 0.e0_double), ' of ',size(Cat%Redshift), ' galaxies have redshift information'

    !--If a catalogue for the blank field is passed in, then use this catalogue to get the intrinsic distributions--!
    if(reproduce_Prior) then
       if(present(Blank_Field_Catalogue) == .false.) STOP 'DM_Profile_Variable_Posteriors_CircularAperture - Blankf Field Catalogue must be entered to allow for the production of the prior on a grid'
       write(*,'(A)') 'Producing Distribution from Catalogue'
       call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, Distribution_Directory, Blank_Field_Catalogue)
    else
       write(*,'(A)') 'Reading in distribution from:', Distribution_Directory
       call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, Distribution_Directory)
       print *, 'Success:', Distribution_Directory
    end if

    if(size(magBins,1) > 1) STOP 'I have had to disable Magnitude Binning for now, youll have to edit the code to get this to work, stopping'

    !--Declarations for the use of a coarse likelihood grid
    !--Note that the Fitting Process used may include a search for the maximum, adding points to this grid
    select case(Surface_Mass_Profile)
    case(1) !-Flat-!
       nCoarseGrid = 100000
       Alpha_Min = -5.e-3_double; Alpha_Max = 1.e-2_double !-- SMD ~ Masses 10^12 -> 10^15 Msun/h, in Units of 10^18 Msun/h  
    case(2) !-SIS-!
       nCoarseGrid = 10000 !--This needs to be sigma_v^2--!
       Alpha_Min= -2.5e5_double; Alpha_Max = 9.e6_double    !--------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case(3) !-NFW-!
       nCoarseGrid = 100
       Alpha_Min= 0.05e0_double; Alpha_Max = 3.e0_double
    case default
       STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
    end select
    call equalscale(Alpha_min, Alpha_Max, nCoarseGrid, Coarse_LikelihoodGrid)

    Do Ap = 1, size(Ap_Cats)
       write(Output_File_Prefix,'(I2)') Ap
       Output_File_Prefix = trim(adjustl(Bayesian_Routines_Output_Directory))//'Aperture_'//trim(adjustl(Output_File_Prefix))//'_'

       !--Set input Posterior Grid
       allocate(Posterior_Single(2,size(Coarse_LikelihoodGrid))); Posterior_Single(1,:) = Coarse_LikelihoodGrid

       !--Evaluate posterior on grid
       call DM_Profile_Variable_Posterior_SingleFit(Ap_Cats(Ap), Surface_Mass_Profile, Lens_Redshift, Ap_Pos(Ap,:), Posterior_Single, Output_File_Prefix, use_lnSize_Prior, PriorMagGrid = MagGrid, PriorSizeGrid = SizeGrid, Prior = Joint_Size_Magnitude_Distribution, MagPrior = Magnitude_Distribution)  

       !--Allow the finalised posterior to be determined on a seperate grid to the input posterior, to account for the fact that the posterior routine may have added points as part of the search for the maximum of the posterior
       if(Ap == 1) then
          allocate(Posteriors(size(Ap_Cats), 2, maxval((/500, 2*size(Posterior_Single,2)/)))); Posteriors = 0.e0_double
          if(size(Posteriors,3) < size(Posterior_Single,2)) print *, 'WARNING - DM_Profile_Variable_Posteriors_CircularAperture - Posterior interpolated on COARSER grid than output'
          !--Set up Posterior Grid
          do i = 1, size(Posteriors,3)
             Posteriors(:,1,i) = minval(Posterior_Single(1,:)) + (i-1)*((maxval(Posterior_Single(1,:))-minval(Posterior_Single(1,:)))/(size(Posteriors,3)-1))
          end do
       end if
       Posteriors(Ap,2,:) = Linear_Interp(Posteriors(Ap,1,:), Posterior_Single(1,:), Posterior_Single(2,:))

       deallocate(Posterior_Single)
    end Do
    if(allocated(SizeGrid)) deallocate(SizeGrid)
    if(allocated(MagGrid)) deallocate(MagGrid)
    if(allocated(Joint_Size_Magnitude_Distribution)) deallocate(Joint_Size_Magnitude_Distribution)


    !--Output Posterior--!
    !--Moved to Single Run Routine
!!$    open(unit = 51, file = trim(Bayesian_Routines_Output_Directory)//'Posterior_per_Aperture.dat')
!!$    write(fmtstring,'(I2)') size(Posteriors,1)+1 !-Assumes all posteriors described by the same posterior grid-!
!!$    do j = 1, size(Posteriors,3)
!!$       write(51, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posteriors(1,1,j), Posteriors(:,2,j)
!!$    end do
!!$    close(51)
!!$    print *,'Output file to: ', trim(Bayesian_Routines_Output_Directory)//'Posterior_per_Aperture.dat'


  end subroutine DM_Profile_Variable_Posteriors_CircularAperture


  subroutine find_Maximum_by_Bisection(Likelihood, New_Point, tol, Keep_Going)
    use derivatives_lib, only:derivative_2pt; use Common_Functions, only: setNaN
    !--Assumes that the likelihood is well behaved and has a single maximum
    !--Returns the Likelihood with the addition of the bisected point, and the index of that point
    real(double), intent(inout),allocatable:: Likelihood(:,:)
    integer, intent(out):: New_Point
    real(double), intent(in)::tol
    logical, intent(out):: Keep_Going

    integer:: Maximum_Index
    integer:: i

    integer, save:: callcount, calltolerance = 100
    integer:: Bracket(3)
    logical::Bracket_Found

    real(double),dimension(size(Likelihood,1), size(Likelihood,2)+1):: tLikelihood
    real(double)::NewGridPoint

    !--Attemptd means using derivatives, but I am not convinced this works, therefroe obsolete
!!$    iMaximum_Boundary = (/1, size(Likelihood,2)/)
!!$    if(present(Maximum_Boundary)) then
!!$       iMaximum_Boundary = Maximum_Boundary
!!$    end if
!!$
!!$    !--Find maximum between boundary points - possibly do this using maxloc?
!!$    Maximum_Index = -1
!!$    do i = Maximum_Boundary(1), Maximum_Boundary(2)
!!$       if((derivative_2pt(Likelihood(2,i-1), Likelihood(2,i), Likelihood(1,i)-Likelihood(1,i-1)) > 0.e0_double) .and. (derivative_2pt(Likelihood(2,i+1), Likelihood(2,i+2), Likelihood(1,i+2)-Likelihood(1,i+1))< 0.e0_double)) Maximum_Index = i
!!$       !-Maximum_Index marks the lower boundary of the interval in which the maximum is exptect to occur
!!$    end do
!!$    if(Maximum_Index < 0) STOP 'find_Maximum_by_Bisection - FAILED TO FIND MAXIMUM'
!!$ 
!!$    print *, 'Maximum found to be between:', Likelihood(1,i), Likelihood(1,i+1)

    !--If call count too high, output warning, assign flag and output (flag can be used toeither stop the code or use the result as it stands)
    if(callcount > calltolerance) then
    print *, 'WARNING: CALL TOLERANCE FOR MAXIMUM FINDING IS TOO HIGH, EXITING'
       CALLCOUNT = 0
       Keep_Going = .false.
       return
    end if

    !--Find Bracket of maximum on original course grid
    !----This could potentially be sped up by keeping the previous braceting interval and only updating within some region of it
    Bracket_Found = .false.; Bracket = 0.
    do i = 2, size(Likelihood,2)-1
       if(Likelihood(2,i-1)<Likelihood(2,i) .and. Likelihood(2,i)>=Likelihood(2,i+1)) then
          if(Bracket_Found) then
             !--Compare the brackets to find the smallest (overcoming local minima
             print *, 'Multiple Mimina found!', Likelihood(:,i), ':', Likelihood(:,Bracket(2))
             if(Likelihood(2,i) > Likelihood(2, Bracket(2))) then
                cycle
             end if
             print *, 'NEW MAXIMA SET!'
          end if

          !--Bracket Found-!
          Bracket = (/i-1, i, i+1/) !--Middle of bracket is best guess at the minimum
          Bracket_Found = .true.
       end if
    end do

    if(Bracket(1) == 0 .or. Bracket(3) == size(Likelihood,2)) then
       print *, 'Maximum found to be intolerably close to the limits of the original grid, exiting without result'
       callcount = 0
       Keep_Going = .false.
       return
    end if

    !--Found within tolerance
    if( maxval((/Likelihood(1,Bracket(2))-Likelihood(1,Bracket(1)),Likelihood(1,Bracket(3))-Likelihood(1,Bracket(2))/)) <= tol ) then
       print *, 'Found Maximum to within a tolerance, expected to be:', Likelihood(1,Bracket(2))
       callcount = 0
       Keep_Going = .false.
       return
    end if

    !--Bisect the Grid and add this point onto the posterior 
!    NewGridPoint = 0.5e0_double*(Likelihood(1,Bracket(2))+Likelihood(1,Bracket(2+((-1)**callcount)))) !--1 allows for alternate sampling to left and right of entered bracket
 
    !--Add new point using golden ratio
    !--(which sets the new point at 0.38197 into the larger of the intervals)
    if( (Likelihood(1,Bracket(2))-Likelihood(1,Bracket(1))) >= (Likelihood(1,Bracket(3))-Likelihood(1,Bracket(2))) ) then
       NewGridPoint = Likelihood(1,Bracket(2)) - 0.38197e0_double*(Likelihood(1,Bracket(2))-Likelihood(1,Bracket(1)))
       Maximum_Index = Bracket(1)
    else
       NewGridPoint = Likelihood(1,Bracket(2)) + 0.38197e0_double*(Likelihood(1,Bracket(3))-Likelihood(1,Bracket(2)))
       Maximum_Index = Bracket(2)
    end if

!-- Use this if not using golden ratio--!
!!$    if( (-1)**callcount > 0) then
!!$       Maximum_Index = Bracket(2)
!!$    else
!!$       Maximum_Index = Bracket(1)
!!$    END if
    !--Assign new point
    tLikelihood(:,:Maximum_Index) = Likelihood(:,:Maximum_Index)
    tLikelihood(1,Maximum_Index+1) = NewGridPoint; tLikelihood(2,Maximum_Index+1) = setNaN()
    tLikelihood(:,Maximum_Index+2:) = Likelihood(:,Maximum_Index+1:)

    New_Point = Maximum_Index+1


    deallocate(Likelihood)
    allocate(Likelihood(size(tLikelihood,1), size(tLikelihood,2))); Likelihood = tLikelihood

  end subroutine find_Maximum_by_Bisection

  subroutine DM_Profile_Variable_Posterior_SingleFit(Cat, Mass_Profile, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior, PriorMagGrid, PriorSizeGrid, Prior, MagPrior)
    use cosmology, only:angular_diameter_distance_fromRedshift; use MC_Redshift_Sampling, only: Monte_Carlo_Redshift_Sampling_SigmaCritical; use Mass_Profiles; use Distributions, only: ch08_redshift_distribution_Array, CH08_Redshift_Distribution_Scalar; use Interpolaters, only: Linear_Interp; use Bayesian_Posterior_Evaluation, only: lnLikelihood_Evaluation_atVirialRadius_perSourceSample, get_Likelihood_Evaluation_Precursors
    use Integration, only:TrapInt, Integrate; use Matrix_methods, only: Determinant, Matrix_Invert; use Smoothing, only: KDE_BiVariate_Gaussian_Scalar, KDE_UniVariate_Gaussian; use Statistics, only:Discrete_Covariance, mean_discrete; use Common_Functions, only: setNaN
    !--Returns the Bayesian posterior for the DM Profile Variables, as defined by Mass_Profile
    !-For Mass_Profile = 1 (Flat): Sigma_0
    !-                   2 (SIS) : Velocity_Dispersion**2
    !-                   3 (NFW) : Virial Radius (r200)
    !--Prior is the prior distribution of p(m,R), or possibly p(m,lnR), m is *apparent magnitude*, R is *apparent size*
    !--Cat entered already assumed to be only the galaxies that are to be considered (e.g. in aperture), and well-described by posterior (e.g. same mag bin) 
    !--Mag-Prior entered seperately as in use we would require that no size cuts are used in the construction of this prior (including with KDE smoothing).

    !---TO DO:

    type(Catalogue)::Cat
    integer,intent(in)::Mass_Profile !-1:Flat, 2:SIS, 3: NFW-!
    real(double),intent(in)::Lens_Redshift
    real(double),intent(in)::Lens_Position(2)
    !--Posterior(1,:) MUST contain the grid on which the posterior is to be evaluated
    real(double),allocatable,intent(InOut)::Posterior(:,:) !-Grid/Posterior, Value-!
    character(*), intent(in)::Output_Prefix
    logical,intent(in)::lnSize_Prior
    real(double),intent(in),optional:: PriorSizeGrid(:), PriorMagGrid(:), MagPrior(:), Prior(:,:) !-Magnitude, Size-!

    !---INTERNAL DECLARATIONS
    integer::i, c, j, z, m
    integer:: Galaxy_Posterior_Method

    real(double),dimension(:,:),allocatable:: Survey_Renormalised_Prior
    real(double),dimension(:),allocatable:: Survey_Renormalised_MagPrior
    real(double),allocatable::Posterior_perGalaxy(:) !-Galaxy, Posterior-! - Uses Same Grid as overall Posterior -

    !--Cluster Model Delcarations
    real(double),allocatable::Sigma_Crit(:,:)
    real(double)::D_l, D_s, D_ls
    real(double),allocatable::Distance_from_Mass_Center(:) !-Galaxy-!

    logical:: Output_Posterior_Per_Galaxy = .true.
    character(7)::fmtstring
    character(500)::Filename
    logical::here

    !--Redshift Distribution Declarations--!
!!$    integer,parameter:: nRedshift_Sampling = 50
!!$    real(double),parameter::Redshift_Lower = Lower_Redshift_Cut, Redshift_Higher = 4.e0_double 
    real(double), allocatable:: RedshiftGrid(:)

    real(double),dimension(size(Cat%RA),2):: Source_Positions

    !--Kappa dependant renormalisation--!
!DELETE    real(double),dimension(2):: Renormalisation_Magnitude_Limits, Renormalisation_Size_Limits
    real(double),allocatable:: MagnificationGrid(:), Renormalisation_by_Magnification(:), Convergence_Renorm_PerGalaxy(:,:),  MagOnly_Renormalisation_by_Magnification(:)
    !vv Must be set here vv!
    integer:: nMagnificationGrid = 2000
    real(double):: MagFactorGridLower = 1.e0_double, MagFactorGridHigher = 65.e0_double
    integer:: IntegrationFlag = -1000

!    integer:: nMagPosterior, nSizePosterior, nSizeMagPosterior

    !--Maximum Search Options
    logical:: Continue_To_Evaluate
    integer:: Alpha_Loop
    logical:: Search_For_Maximum = .true.
    !~~ Initial Grid Size and Tolerance to which the maximum are found set the error on confidence limits, and mode respectively, as well as setting minimum time for run
    integer:: Search__Coarse_Grid_Size = 100
    real(double):: Find_Maximum_Tolerance = 1.e-2_double

    !--Testing Declarations--!
    integer:: n_Default_Source_Redshift_Used, nGal_Ignored_MagLimits, nGal_Ignored_SizeLimits, nGal_Ignored_NaN

    INTERFACE
       subroutine DM_Profile_Variable_Posterior_SingleFit(Cat, Mass_Profile, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior, PriorMagGrid, PriorSizeGrid, Prior, MagPrior)
         use Param_Types; use Catalogues
         type(Catalogue)::Cat
         integer,intent(in)::Mass_Profile !-1:Flat, 2:SIS, 3:NFW-!
         real(double),intent(in)::Lens_Redshift
         real(double),intent(in)::Lens_Position(2)
         real(double),allocatable,intent(Out)::Posterior(:,:) !-Grid/Posterior, Value-!
         character(*), intent(in)::Output_Prefix
         logical,intent(in)::lnSize_Prior

         real(double),intent(in),optional:: PriorSizeGrid(:), PriorMagGrid(:), MagPrior(:), Prior(:,:) !-Magnitude, Size-!
       END subroutine DM_Profile_Variable_Posterior_SingleFit
    END INTERFACE
    
    print *, 'Called DM_Profile_Variable_Posterior'
    
    if(Analyse_with_Physical_Sizes) STOP 'DM_Profile_Variable_Posterior - I HAVE DISABLED THE ABILITY TO USE PHYISCAL SIZES AS UNNECESSARY, code still to be edited'
    
!!$    if((present(PriorMagGrid) == .false.) .or. (present(PriorSizeGrid)== .false.)) STOP 'DM_Profile_Variable_Posterior - Prior must be accompanied by grids'
!!$    allocate(Survey_Renormalised_Prior(size(Prior,1),size(Prior,2))); Survey_Renormalised_Prior = 0.e0_double
!!$    if(present(MagPrior)) then
!!$       allocate(Survey_Renormalised_MagPrior(size(PriorMagGrid))); Survey_Renormalised_MagPrior = 0.e0_double
!!$    end if
    print *, ' '
    print *, 'Attempting Prior Interpolation without Extrapolation'
    print *, ' '
    
    !--Check posterior grid set up correctly
    if(allocated(Posterior) == .false.) STOP 'DM_Profile_Variable_Posterior - Posterior Grid MUST be entered.'

    !--Start of Posterior Routines--!
    allocate(Posterior_perGalaxy(size(Cat%RA))); Posterior_perGalaxy = 1.e-100_double

    write(*,'(A)') '--------------------------------------------------------------------------------------------------'
    write(*,'(A)',advance = 'no') 'Getting Posterior for Cluster '
    if(Posterior_Method == 1) write(*,'(A)') 'using Sizes Only'
    if(Posterior_Method == 2) write(*,'(A)') 'using Sizes and Magnitudes'
    if(Posterior_Method == 3) write(*,'(A)') 'using Magnitudes Only'
    if(Posterior_Method == 4) write(*,'(A)') 'using Size and Magnitudes, and Magnitudes Only below the size limit'
    if(Enforce_Weak_Lensing) print *, '*** Weak lensing assumptions have been enforced'

    !--Set defaults for supplementary delarations
    n_Default_Source_Redshift_Used = 0; nGal_Ignored_MagLimits = 0; nGal_Ignored_SizeLimits = 0; nGal_Ignored_NaN = 0
!!DELETE    nMagPosterior = 0; nSizePosterior = 0; nSizeMagPosterior = 0

    !--Set up precursors
    call get_Likelihood_Evaluation_Precursors(Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, RedshiftGrid, Sigma_Crit, MagnificationGrid,Renormalisation_by_Magnification, MagOnly_Renormalisation_by_Magnification, PriorSizeGrid, PriorMagGrid, MagPrior, Prior, Survey_Magnitude_Limits, Survey_Size_Limits, (/Lens_Redshift/), Lower_Redshift_Cut, Output_Prefix)

    Source_Positions(:,1) = Cat%RA; Source_Positions(:,2) = Cat%Dec

    i = 0; Alpha_Loop = 0
    Continue_To_Evaluate = .true.
    do while (Continue_To_Evaluate)
       Posterior_perGalaxy = setNaN()
       !--Alpha_Loop counts the number of times the Alpha value has been looped over
       Alpha_Loop = Alpha_Loop+1
       i = i+1

       !--Set exit strategy in normal case
       if(Search_For_Maximum == .false. .and. Alpha_Loop > size(Posterior,2)) Continue_To_Evaluate = .false.
       if(Search_For_Maximum .and. Alpha_Loop>size(Posterior,2)) then
          !--Find the Maximum. In this case, i will be set to the new grid point (added to grid), and exit case will be set when either tolerance or call limits are met 
          call  find_Maximum_by_Bisection(Posterior, i, Find_Maximum_Tolerance, Continue_To_Evaluate)
       end if

       !--Exit here to remove need to evaluate an extra point when maximum has been found
       if(Continue_to_Evaluate == .false.) exit

       Posterior(2,i) = 1.e0_double
       if(Alpha_Loop == Size(Posterior,2)/2) print *, 'Approximately halfway done for this Aperture..'          
!!$       do c = 1, size(Posterior_perGalaxy,1) !-Loop over galaxies-!
          
          !~~Select the method of posterior reconstruction for that point based in input method   
!!$          select case(Posterior_Method)
!!$          case(1) !--SizeOnly--!
!!$             if( (Cat%Sizes(c) > Survey_Size_Limits(2)) .or. (Cat%Sizes(c) < Survey_Size_Limits(1))) then
!!$                nGal_Ignored_SizeLimits = nGal_Ignored_SizeLimits + 1
!!$                cycle
!!$             end if
!!$             !MagRenorm          iSurvey_size_Limits = Survey_Size_Limits
!!$             
!!$             nSizePosterior = nSizePosterior + 1
!!$             Galaxy_Posterior_Method = 1
!!$          case(2)!--SizeMag--!
!!$             if( (Cat%Sizes(c) > Survey_Size_Limits(2)) .or. (Cat%Sizes(c) < Survey_Size_Limits(1))) then
!!$                nGal_Ignored_SizeLimits = nGal_Ignored_SizeLimits + 1
!!$                cycle
!!$             end if
!!$             
!!$             !MagRenorm          iSurvey_size_Limits = Survey_Size_Limits
!!$             
!!$             nSizeMagPosterior = nSizeMagPosterior + 1
!!$             Galaxy_Posterior_Method = 2
!!$          case(3) !-Magnitude Only--!
!!$             nMagPosterior = nMagPosterior + 1
!!$             
!!$             !MagRenorm          iSurvey_size_Limits = (/0.e0_double, 1.e30_double/) !--Should encompase the whole data set
!!$             
!!$             if(present(MagPrior) == .false.) STOP 'Attempting to produce posterior using mag only, but no magnitde prior entered, stopping'
!!$             Galaxy_Posterior_Method = 3
!!$          case(4) !-Size mag above size data limit, Mag-Only below size data limit-!
!!$             if(present(MagPrior) == .false.) STOP 'Attempting to produce posterior using mag only under size limit, but no magnitde prior entered, stopping'
!!$             STOP 'Posterior Method 4 has been disabled as it needs further thought into its construction. In particular, with respect to the construction of the prior from the size-mag distriution, the effect of prior cuts, and correct renormalisation. The skeleton code has been added for this, but disabled for now. See commented code labeled MagRenorm'
!!$             
!!$             !--Note, in this case the data-renormalisation should take inot account that the mag only is constructed from a sample which takes small, faint galaxies
!!$             if(Cat%Sizes(c) < Survey_Size_Limits(1)) then
!!$                !MagRenorm             iSurvey_Size_Limits = (/0.e0_double, Survey_Size_Limits(1)/)
!!$                
!!$                Galaxy_Posterior_Method = 3
!!$                nMagPosterior = nMagPosterior + 1
!!$             else
!!$                !MagRenorm             iSurvey_Size_Limits = Survey_Size_Limits
!!$                
!!$                Galaxy_Posterior_Method = 2
!!$                nSizeMagPosterior = nSizeMagPosterior + 1
!!$             end if
!!$          end select
          
          !~~~ Cycle if the observed magnitude falls outside the survey limits
!!$          if( (Cat%MF606W(c) > Survey_Magnitude_Limits(2)) .or. (Cat%MF606W(c) < Survey_Magnitude_Limits(1))) then
!!$             nGal_Ignored_MagLimits = nGal_Ignored_MagLimits + 1
!!$             cycle
!!$          end if
!!$          if(isNaN(Cat%Sizes(c)) .or. isNaN(Cat%MF606W(c))) then
!!$             nGal_Ignored_NaN = nGal_Ignored_NaN + 1
!!$             cycle
!!$          end if
!!$
!!$
!!$          !--Evaluate the Posterior 
!!$          Posterior_perGalaxy(c) = Likelihood_Evaluation_atVirialRadius_SingleCluster(Posterior(1,i), Galaxy_Posterior_Method, Mass_Profile, Lens_Position, Lens_Redshift, Cat%Sizes(c), Cat%MF606W(c), Cat%Redshift(c), (/Cat%RA(c),Cat%Dec(c)/), PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, Survey_Size_Limits, Survey_Magnitude_Limits, MagnificationGrid, MagOnly_Renormalisation_by_Magnification, Renormalisation_by_Magnification, RedshiftGrid, Sigma_Crit(1,:))
!!$                    
!!$       end do !--End of galaxy loop       
!!$
!!$       !~~Return lnP to ensure that we can renormalise in an alpha-independent way in the combined posterior, and to ensure that the PDF is correctly recovered even with ronding error (large negative lnP across all alphas). The renormalisation is done after the outer loop is finished
!!$       call Combine_Posteriors(Posterior(1,i), Posterior_perGalaxy(:), Combine_log_Posteriors, Renormalise = .false., Return_lnP = Combine_log_Posteriors, Combined_Posterior = Posterior(2,i))

          Posterior(2,i) = lnLikelihood_Evaluation_atVirialRadius_perSourceSample(Posterior(1,i), Posterior_Method, Mass_Profile, Lens_Position, Lens_Redshift, Cat%Sizes, Cat%MF606W, Cat%Redshift, Source_Positions, PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, Survey_Size_Limits, Survey_Magnitude_Limits, MagnificationGrid, MagOnly_Renormalisation_by_Magnification, Renormalisation_by_Magnification, RedshiftGrid, Sigma_Crit(1,:))

       end do !--End of Posterior Loop
!    deallocate(Posterior_perGalaxy)

    if(any(dabs(Posterior) > huge(1.e0_double)) .or. any(isNaN(Posterior))) then
       print *, 'NaNs or infinities found in posterior:'
       do i =1, size(Posterior,2)
          print *, Posterior(:,i)
       end do
       STOP
    end if

    !--Convert from lnP to P for posterior, renormalise so that max(P) = 1, to avoid rounding errors
    Posterior(2,:) = dexp(Posterior(2,:) - maxval(Posterior(2,:)))

    if(n_Default_Source_Redshift_Used > 0) write(*,'(A,I3,A)') '****** Used the default redshift for:', n_Default_Source_Redshift_Used, ' galaxies ********'
    if(nGal_Ignored_SizeLimits > 0) write(*,'(A,I3)') '****** Number of galxies ignored as they fell outside the survey size limits:', nGal_Ignored_SizeLimits
    if(nGal_Ignored_NaN > 0) write(*,'(A,I3)') '****** Number of galaxies ignored as they were NaNs:', nGal_Ignored_NaN
    
    
!!DELETE$    print *, '-------------------------------------------------------------'
!!$    print *, 'Constructed ', nSizePosterior, ' size-only posteriors'
!!$    print *, 'Constructed ', nSizeMagPosterior, ' size-magnitude posteriors'
!!$    print *, 'Constructed ', nMagPosterior, ' magnitude-only posteriors'
!!$    print *, '-------------------------------------------------------------'

    !---DO NOT DELETE THIS. Instead, come up with a way to store the posterior per galaxy on a general grid for output
!!!$    if(Debug_Mode .and. Output_Posterior_Per_Galaxy) then
!!!$       
!!!$       Filename = trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat'
!!!$       open(unit = 17, file = Filename)
!!!$       write(fmtstring,'(I7)') size(Posterior_perGalaxy,1)+1
!!!$       do i =1, size(Posterior,2)
!!!$          write(17, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posterior(1,i), Posterior_perGalaxy(:,i)
!!!$       end do
!!!$       close(17)
!!!$       print *, 'Output Posterior per galaxy to: ', trim(adjustl(Filename))
!!!$    end if
    
    
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
    
    if(allocated(Survey_Renormalised_Prior)) deallocate(Survey_Renormalised_Prior)
    if(allocated(Survey_Renormalised_MagPrior)) deallocate(Survey_Renormalised_MagPrior)
    
    !--On Successful Completion delete Poster per galaxy as large file--!
!!$    inquire(file = trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat', exist = here)
!!$    if(here == .false.) STOP "Posterior per galaxy doesn't exist, stopping before accidental deletion"
!!$    call system('rm '//trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat')
    

  end subroutine DM_Profile_Variable_Posterior_SingleFit
  

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
    real(double):: Redshift_Cuts(2)  = (/Lower_Redshift_Cut, 100.e0_double/) !-- If no cuts, then lower should still be < -1 to ensure only galaxies without redshift information are cut
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


!---------------POTENTIALLY OBSOLETE CODE---------------------------------------------------------------------------------------------------------------------------------------------------------------------!

!--10 Sept 2014 
!----DM_Profile_Variable_Posterior - Copy kept with sets up a grid on which the posterior is to be evaluated. Next step will edit to allow for a grid to be passed in, and will delete unecessary code, INCLUDING ANY KDE COVARAINCE

!!$  subroutine DM_Profile_Variable_Posterior(Cat, Mass_Profile, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior, PriorCatalogue, PriorMagGrid, PriorSizeGrid, Prior, MagPrior)
!!$    use cosmology, only:angular_diameter_distance_fromRedshift; use MC_Redshift_Sampling, only: Monte_Carlo_Redshift_Sampling_SigmaCritical; use Mass_Profiles; use Distributions, only: ch08_redshift_distribution_Array, CH08_Redshift_Distribution_Scalar; use Interpolaters, only: Linear_Interp
!!$    use Integration, only:TrapInt, Integrate; use Matrix_methods, only: Determinant, Matrix_Invert; use Smoothing, only: KDE_BiVariate_Gaussian_Scalar, KDE_UniVariate_Gaussian; use Statistics, only:Discrete_Covariance, mean_discrete; use Common_Functions, only: setNaN
!!$    !--Returns the Bayesian posterior for the DM Profile Variables, as defined by Mass_Profile
!!$    !-For Mass_Profile = 1 (Flat): Sigma_0
!!$    !-                   2 (SIS) : Velocity_Dispersion**2
!!$    !-                   3 (NFW) : Virial Radius (r200)
!!$    !--Prior is the prior distribution of p(m,R), or possibly p(m,lnR), m is *apparent magnitude*, R is *apparent size*
!!$    !--Cat entered already assumed to be only the galaxies that are to be considered (e.g. in aperture), and well-described by posterior (e.g. same mag bin) 
!!$    !--Mag-Prior entered seperately as in use we would require that no size cuts are used in the construction of this prior (including with KDE smoothing).
!!$
!!$    !---TO DO:
!!$
!!$    type(Catalogue)::Cat
!!$    integer,intent(in)::Mass_Profile !-1:Flat, 2:SIS, 3: NFW-!
!!$    real(double),intent(in)::Lens_Redshift
!!$    real(double),intent(in)::Lens_Position(2)
!!$    !--Posterior(1,:) MUST contain the grid on which the posterior is to be evaluated
!!$    real(double),allocatable,intent(InOut)::Posterior(:,:) !-Grid/Posterior, Value-!
!!$    character(*), intent(in)::Output_Prefix
!!$    logical,intent(in)::lnSize_Prior
!!$    type(Catalogue), intent(in), optional:: PriorCatalogue
!!$    real(double),intent(in),optional:: PriorSizeGrid(:), PriorMagGrid(:), MagPrior(:), Prior(:,:) !-Magnitude, Size-!
!!$
!!$    !---INTERNAL DECLARATIONS
!!$    integer::i, c, j, z, m
!!$    integer:: Galaxy_Posterior_Method
!!$
!!$    !-Size_Only_Prior contains the priors which depends only on size, which is evalutaed for each redshift and integrates over all magnitudes
!!$    real(double),dimension(:),allocatable:: Size_Only_Prior, Mag_Prior, Mag_Only_Prior
!!$    !--Size_Only_Mag_Prior contains the prior for size which evaluated over the mag grid, which will be integrated over. Contains p_[theta_0, m_0|z]*p[z|m_0] for all m in grid
!!$    real(double),dimension(:,:),allocatable::  Kappa_Renormalised_Prior, Survey_Renormalised_Prior, Size_Only_Mag_Prior
!!$    real(double),dimension(:),allocatable::Kappa_Renormalised_MagPrior, Survey_Renormalised_MagPrior
!!$    real(double),allocatable::Posterior_perGalaxy(:) !-Galaxy, Posterior-! - Uses Same Grid as overall Posterior -
!!$    real(double):: Effective_Magnification !- Only Model Dependant
!!$
!!$    logical:: Marginalise_Redshift_Distribution = .true.
!!$    logical::Known_Redshift
!!$
!!$    !--KDE_Smoothing Declarations--!
!!$    real(double),allocatable:: KDE_Gaussian_Covariance(:,:), Data_Vectors(:,:), KDE_Covariance_Inverse(:,:)
!!$    real(double):: KDE_Covariance_Determinant, KDE_Gaussian_Covariance_Reduction = 0.01e0_double !-How much is sig^2 which give KDE width reduced from measured covariance?--!
!!$    logical:: do_KDE_Extrapolation, do_KDE_OnTheFly
!!$    logical:: need_Extrapolate
!!$
!!$    real(double),allocatable::Sigma_Crit(:), Sigma_Crit_MC(:)
!!$    real(double):: Galaxy_Sigma_Critical
!!$    real(double)::D_l, D_s, D_ls
!!$    real(double),allocatable::Distance_from_Mass_Center(:) !-Galaxy-!
!!$
!!$    real(double)::Renorm
!!$
!!$    logical:: Output_Posterior_Per_Galaxy = .true.
!!$    character(7)::fmtstring
!!$    character(500)::Filename
!!$    logical::here
!!$
!!$    !--Redshift Distribution Declarations--!
!!$    integer,parameter:: nRedshift_Sampling = 50
!!$    real(double),parameter::Redshift_Lower = Lower_Redshift_Cut, Redshift_Higher = 4.e0_double !!!Edit to Lens_Redshift
!!$    real(double), dimension(nRedshift_Sampling):: RedshiftGrid
!!$    real(double)::RedshiftPDF
!!$    real(double),allocatable:: Posterior_perGalaxy_Redshift(:) !-Galaxy, Posterior, Redshift-!
!!$
!!$    !--Kappa dependant renormalisation--!
!!$!!    real(double),dimension(2):: Survey_Magnitude_Limits, Survey_Size_Limits !--User defined needs edit to below, these are set by default from distribution-!
!!$    real(double),dimension(2):: Renormalisation_Magnitude_Limits, Renormalisation_Size_Limits
!!$    real(double),allocatable:: MagnificationGrid(:), Renormalisation_by_Magnification(:), Convergence_Renorm_PerGalaxy(:,:),  MagOnly_Renormalisation_by_Magnification(:)
!!$    !vv Must be set here vv!
!!$    integer:: nMagnificationGrid = 2000
!!$    real(double):: MagFactorGridLower = 1.e0_double, MagFactorGridHigher = 65.e0_double
!!$    integer:: IntegrationFlag = -1000
!!$
!!$    integer:: nMagPosterior, nSizePosterior, nSizeMagPosterior
!!$
!!$    !--Maximum Search Options
!!$    logical:: Continue_To_Evaluate
!!$    integer:: Alpha_Loop
!!$    logical:: Search_For_Maximum = .true.
!!$    !~~ Initial Grid Size and Tolerance to which the maximum are found set the error on confidence limits, and mode respectively, as well as setting minimum time for run
!!$    integer:: Search__Coarse_Grid_Size = 100
!!$    real(double):: Find_Maximum_Tolerance = 1.e-2_double
!!$
!!$    !--Testing Declarations--!
!!$    real(double),dimension(3,size(Cat%RA)):: Convergence_per_Cluster
!!$    integer:: n_Default_Source_Redshift_Used, nGal_Ignored_MagLimits, nGal_Ignored_SizeLimits, nGal_Ignored_NaN
!!$    real(double), allocatable:: Aperture_Smoothed_Size_PDF(:), Aperture_Smoothed_Mag_PDF(:)
!!$    logical:: Produce_Shift_inAperture = .false.
!!$    real::Time1, Time2, Time3, Time4
!!$    integer::test
!!$
!!$    INTERFACE
!!$       subroutine DM_Profile_Variable_Posterior(Cat, Mass_Profile, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior, PriorCatalogue, PriorMagGrid, PriorSizeGrid, Prior, MagPrior)
!!$         use Param_Types; use Catalogues
!!$         type(Catalogue)::Cat
!!$         integer,intent(in)::Mass_Profile !-1:Flat, 2:SIS, 3:NFW-!
!!$         real(double),intent(in)::Lens_Redshift
!!$         real(double),intent(in)::Lens_Position(2)
!!$         real(double),allocatable,intent(Out)::Posterior(:,:) !-Grid/Posterior, Value-!
!!$         character(*), intent(in)::Output_Prefix
!!$         logical,intent(in)::lnSize_Prior
!!$
!!$         type(Catalogue), intent(in), optional:: PriorCatalogue
!!$         real(double),intent(in),optional:: PriorSizeGrid(:), PriorMagGrid(:), MagPrior(:), Prior(:,:) !-Magnitude, Size-!
!!$       END subroutine DM_Profile_Variable_Posterior
!!$    END INTERFACE
!!$    
!!$    print *, 'Called DM_Profile_Variable_Posterior'
!!$    
!!$    if(Analyse_with_Physical_Sizes) STOP 'DM_Profile_Variable_Posterior - I HAVE DISABLED THE ABILITY TO USE PHYISCAL SIZES AS UNNECESSARY, code still to be edited'
!!$    
!!$    do_KDE_Extrapolation = .false.; do_KDE_OnTheFly = .false.
!!$    if(present(Prior) .or. present(MagPrior)) then
!!$       if((present(PriorMagGrid) == .false.) .or. (present(PriorSizeGrid)== .false.)) STOP 'DM_Profile_Variable_Posterior - Prior must be accompanied by grids'
!!$       allocate(Kappa_Renormalised_Prior(size(Prior,1), size(Prior,2))); Kappa_Renormalised_Prior = 0.e0_double
!!$       allocate(Survey_Renormalised_Prior(size(Prior,1),size(Prior,2))); Survey_Renormalised_Prior = 0.e0_double
!!$       if(present(MagPrior)) then
!!$          allocate(Survey_Renormalised_MagPrior(size(PriorMagGrid))); Survey_Renormalised_MagPrior = 0.e0_double
!!$          allocate(Kappa_Renormalised_MagPrior(size(PriorMagGrid))); Kappa_Renormalised_MagPrior =0.e0_double
!!$       end if
!!$       allocate(Size_Only_Mag_Prior(size(Prior,1),size(Prior,2))); Size_Only_Mag_Prior = 0.e0_double
!!$       if(allow_KDE_Extrapolation .and. present(PriorCatalogue)) then
!!$          do_KDE_Extrapolation = .true.
!!$          print *, ' '
!!$          print *, 'Attempting Prior Interpolation with KDE Extrapolation'
!!$          print *, ' '
!!$       else
!!$          print *, ' '
!!$          print *, 'Attempting Prior Interpolation without Extrapolation'
!!$          print *, ' '
!!$       end if
!!$    elseif(present(PriorCatalogue)) then
!!$       do_KDE_OnTheFly = .true.
!!$       print *, ' '
!!$       print *, 'Attempting KDE on the Fly - NOTE This does not include Kappa-Renormalisation (Seg fault straight after this..)'
!!$       print *, ' '
!!$    end if
!!$
!!$    !-Set Up Posterior Grid-!
!!$!!!$    select case(Mass_Profile)
!!$!!!$    case(1) !-Flat-!
!!$!!!$       nGrid = 100000
!!$!!!$       VGrid_Lower = -5.e-3_double; VGrid_Higher = 1.e-2_double !-- SMD ~ Masses 10^12 -> 10^15 Msun/h, in Units of 10^18 Msun/h  
!!$!!!$    case(2) !-SIS-!
!!$!!!$       nGrid = 10000 !--This needs to be sigma_v^2--!
!!$!!!$       VGrid_Lower= -2.5e5_double; VGrid_Higher = 9.e6_double    !--------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!$    case(3)
!!$!!!$       if(Search_for_Maximum) then
!!$!!!$          !--Set initial grid to a smaller value since the search for the maximum will reduce the number of points needed
!!$!!!$          nGrid = Search__Coarse_Grid_Size
!!$!!!$       else
!!$!!!$          nGrid = 500
!!$!!!$       end if
!!$!!!$       VGrid_Lower= 0.05e0_double; VGrid_Higher = 3.e0_double
!!$!!!$    case default
!!$!!!$       STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
!!$!!!$    end select
!!$!!!$    allocate(Posterior(2, nGrid)); Posterior = setNaN() !*!
!!$!!!$    do i =1, nGrid
!!$!!!$       Posterior(1,i) = VGrid_Lower + (i-1)*((VGrid_Higher-VGrid_Lower)/(nGrid-1))
!!$!!!$    end do
!!$    !----------------------!
!!$
!!$    if(allow_KDE_Extrapolation .or. KDE_OnTheFly) then
!!$       if(present(PriorCatalogue) == .false.) STOP 'DM_Profile_Variable_Posterior - Prior Catalogue needs to be entered to allow KDE Extrapolation'
!!$       !--Construct the Covariance that will be used for the KDE Smoothing--!
!!$       allocate(Data_Vectors(2,size(PriorCatalogue%Sizes))); Data_Vectors(1,:) = PriorCatalogue%MF606W; Data_Vectors(2,:) = PriorCatalogue%Sizes
!!$       call Discrete_Covariance(Data_Vectors, KDE_Gaussian_Covariance)
!!$       KDE_Gaussian_Covariance = KDE_Gaussian_Covariance_Reduction*KDE_Gaussian_Covariance
!!$       call Matrix_Invert(KDE_Gaussian_Covariance, KDE_Covariance_Inverse, 'S')
!!$       KDE_Covariance_Determinant = Determinant(KDE_Gaussian_Covariance)
!!$       deallocate(KDE_Gaussian_Covariance)
!!$    end if
!!$
!!$    if(produce_Shift_inAperture) then
!!$       !--Testing: output smoothed size and mag distributions within the aperture to see if there is any noticable shift
!!$       allocate(Data_Vectors(2,size(PriorCatalogue%Sizes))); Data_Vectors(1,:) = PriorCatalogue%MF606W; Data_Vectors(2,:) = PriorCatalogue%Sizes
!!$       call Discrete_Covariance(Data_Vectors, KDE_Gaussian_Covariance)
!!$       KDE_Gaussian_Covariance = KDE_Gaussian_Covariance_Reduction*KDE_Gaussian_Covariance
!!$       allocate(Aperture_Smoothed_Size_PDF(size(PriorSizeGrid))); Aperture_Smoothed_Size_PDF = 0.e0_double
!!$       allocate(Aperture_Smoothed_Mag_PDF(size(PriorMagGrid))); Aperture_Smoothed_Mag_PDF = 0.e0_double
!!$       call KDE_Univariate_Gaussian(Cat%Sizes, dsqrt(KDE_Gaussian_Covariance(2,2)), PriorSizeGrid, Aperture_Smoothed_Size_PDF)
!!$       call KDE_Univariate_Gaussian(Cat%MF606W, dsqrt(KDE_Gaussian_Covariance(1,1)), PriorMagGrid, Aperture_Smoothed_Mag_PDF)
!!$       open(unit = 31, file = trim(adjustl(Output_Prefix))//'KDE_Distributions_in_Aperture_Size.dat')
!!$       do i =1, size(PriorSizeGrid)
!!$          write(31, *) PriorSizeGrid(i), Aperture_Smoothed_Size_PDF(i)
!!$       end do
!!$       close(31)
!!$       open(unit = 31, file = trim(adjustl(Output_Prefix))//'KDE_Distributions_in_Aperture_Mag.dat')
!!$       do i =1, size(PriorSizeGrid)
!!$          write(31, *) PriorMagGrid(i), Aperture_Smoothed_Mag_PDF(i)
!!$       end do
!!$       close(31)
!!$       write(*, '(4A)') '**Output distributions in aperture to: ', trim(adjustl(Output_Prefix))//'KDE_Distributions_in_Aperture_Size.dat', ' : ', trim(adjustl(Output_Prefix))//'KDE_Distributions_in_Aperture_Mag.dat'
!!$       print *, '**with mean (Size, Mag):', mean_discrete(Cat%Sizes), mean_discrete(Cat%MF606W) 
!!$       deallocate(Aperture_Smoothed_Size_PDF, Aperture_Smoothed_Mag_PDF)
!!$    !---End of distributions in Aperture
!!$    end if
!!$
!!$    D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
!!$    !--Set up the redshift grid--!
!!$    if(Marginalise_Redshift_Distribution) then
!!$       do z = 1, nRedshift_Sampling
!!$          RedshiftGrid(z) = Redshift_Lower + (z-1)*((Redshift_Higher-Redshift_Lower)/(nRedshift_Sampling-1))
!!$       end do
!!$    end if
!!$
!!$    !--Get Sigma_Critical for each point on the Redshift PDF grid--!
!!$    allocate(Sigma_Crit(size(RedshiftGrid))); Sigma_Crit = 0.e0_double
!!$    do z = 1, size(RedshiftGrid)
!!$       D_s = angular_diameter_distance_fromRedshift(0.e0_double, RedshiftGrid(z))
!!$       D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, RedshiftGrid(z))
!!$       Sigma_Crit(z) = 1.66492e0_double*(D_s/(D_l*D_ls)) !-(10^18 M_Sun/h)-!
!!$    end do
!!$    
!!$    !--Renormalise the prior within these size and magnitude limits--!
!!$    print *, 'Renormalisation of the intrinsic distribution:', Integrate(PriorMagGrid, PriorSizeGrid, Prior, 2, lim1 = Survey_Magnitude_Limits, lim2 = Survey_Size_Limits)
!!$    Survey_Renormalised_Prior = Prior/Integrate(PriorMagGrid, PriorSizeGrid, Prior, 2, lim1 = Survey_Magnitude_Limits, lim2 = Survey_Size_Limits)
!!$
!!$    if(present(MagPrior)) then
!!$       !--Allow for seperate renormalisation of the magnitude prior. This is required if the magnitude prior is constructed from galaxies which are excluded from the joint size-magnitude analysis, e.g. due to size cuts--!
!!$       print *, 'Renormalisation of the intrinsic magnitude distribution:', Integrate(PriorMagGrid, MagPrior, 2, lim = Survey_Magnitude_Limits)
!!$       Survey_Renormalised_MagPrior = MagPrior/Integrate(PriorMagGrid, MagPrior, 2, lim = Survey_Magnitude_Limits)
!!$    end if
!!$    
!!$    !------This section renormalises the likelihood, taking into account kappa-dependent mag and size cuts, however this need to be implemented in the prior-----!
!!$    allocate(MagnificationGrid(nMagnificationGrid)); MagnificationGrid = 0.e0_double
!!$    allocate(Renormalisation_by_Magnification(size(MagnificationGrid))); Renormalisation_by_Magnification = 0.e0_double
!!$    if(present(MagPrior)) then
!!$       allocate(MagOnly_Renormalisation_by_Magnification(size(MagnificationGrid))); MagOnly_Renormalisation_by_Magnification = 0.e0_double
!!$    end if
!!$    allocate(Size_Only_Prior(size(Prior,2))); Size_Only_Prior = 0.e0_double
!!$    !Choose MagGrid lower to be the point where the full magnitude range is swept out:
!!$    MagFactorGridHigher = 1.05e0_double*(10.e0_double**((Survey_Magnitude_Limits(2)-Survey_Magnitude_Limits(1))/2.5e0_double)) !1.05 gives lee-way!
!!$    do i = 1, size(MagnificationGrid)
!!$       MagnificationGrid(i) = MagFactorGridLower + (i-1)*((MagFactorGridHigher- MagFactorGridLower)/(size(MagnificationGrid)-1))
!!$       Renormalisation_Size_Limits = Survey_Size_Limits/dsqrt(MagnificationGrid(i))
!!$       Renormalisation_Magnitude_Limits = Survey_Magnitude_Limits + 2.5e0_double*dlog10(MagnificationGrid(i))
!!$
!!$       Renormalisation_by_Magnification(i) = Integrate(PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, 2, lim1 = Renormalisation_Magnitude_Limits, lim2 = Renormalisation_Size_Limits)
!!$       if(present(MagPrior)) then
!!$
!!$          !--Edit this code to calculate the magnitude renormalisation over the whole size grid if Posterior_Method == 3, and over (0, Survey_Size_Limit(1)) if Posterior_Method == 4
!!$          !-- vv This is only true if there are no cuts on the size-mag prior vv
!!$!MagRenorm          MagOnly_Renormalisation_by_Magnification(i) =  Integrate(PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, 2, lim1 = Renormalisation_Magnitude_Limits, lim2 = Renormalisation_Size_Limits)
!!$          MagOnly_Renormalisation_by_Magnification(i) = Integrate(PriorMagGrid, Survey_Renormalised_MagPrior, 2, lim = Renormalisation_Magnitude_Limits)
!!$       end if
!!$    end do
!!$    deallocate(Size_Only_Prior)
!!$       
!!$    open(unit = 53, file = trim(adjustl(Output_Prefix))//'Renormalisation_by_Magnification.dat')
!!$    do i = 1, size(Renormalisation_by_Magnification)
!!$       write(53, *) MagnificationGrid(i), Renormalisation_by_Magnification(i)
!!$    end do
!!$    write(*,'(2(A))') 'File output: ', trim(adjustl(Output_Prefix))//'Renormalisation_by_Magnification.dat'
!!$
!!$    !--Start of Posterior Routines--!
!!$    allocate(Posterior_perGalaxy(size(Cat%RA))); Posterior_perGalaxy = 1.e-100_double
!!$    Effective_Magnification = 0.e0_double
!!$    Convergence_Per_Cluster = 0.e0_double
!!$
!!$    write(*,'(A)') '--------------------------------------------------------------------------------------------------'
!!$    write(*,'(A)',advance = 'no') 'Getting Posterior for Cluster '
!!$    if(Posterior_Method == 1) write(*,'(A)') 'using Sizes Only'
!!$    if(Posterior_Method == 2) write(*,'(A)') 'using Sizes and Magnitudes'
!!$    if(Posterior_Method == 3) write(*,'(A)') 'using Magnitudes Only'
!!$    if(Posterior_Method == 4) write(*,'(A)') 'using Size and Magnitudes, and Magnitudes Only below the size limit'
!!$    if(Enforce_Weak_Lensing) print *, '*** Weak lensing assumptions have been enforced'
!!$
!!$
!!$    !--Get the distance from the mass centre for each galaxy
!!$    allocate(Distance_from_Mass_Center(size(Cat%RA))); Distance_from_Mass_Center = 0.e0_double
!!$    Distance_from_Mass_Center = dsqrt( (Cat%RA(:)-Lens_Position(1))**2.e0_double + (Cat%Dec(:)-Lens_Position(2))**2.e0_double ) !-in Degrees-!
!!$    Distance_from_Mass_Center = (D_l*Distance_from_Mass_Center*(3.142e0_double/(180.e0_double))) !-in Mpc/h-! 
!!$
!!$
!!$    n_Default_Source_Redshift_Used = 0; nGal_Ignored_MagLimits = 0; nGal_Ignored_SizeLimits = 0; nGal_Ignored_NaN = 0
!!$    nMagPosterior = 0; nSizePosterior = 0; nSizeMagPosterior = 0
!!$    Time1 = 0.; Time2 = 0.
!!$
!!$    i = 0; Alpha_Loop = 0
!!$    Continue_To_Evaluate = .true.
!!$    do while (Continue_To_Evaluate)
!!$       Posterior_perGalaxy = setNaN()
!!$       !--Alpha_Loop counts the number of times the Alpha value has been looped over
!!$       Alpha_Loop = Alpha_Loop+1
!!$       i = i+1
!!$
!!$       !--Set exit strategy in normal case
!!$       if(Search_For_Maximum == .false. .and. Alpha_Loop > size(Posterior,2)) Continue_To_Evaluate = .false.
!!$       if(Search_For_Maximum .and. Alpha_Loop>size(Posterior,2)) then
!!$          !--Find the Maximum. In this case, i will be set to the new grid point (added to grid), and exit case will be set when either tolerance or call limits are met 
!!$          call  find_Maximum_by_Bisection(Posterior, i, Find_Maximum_Tolerance, Continue_To_Evaluate)
!!$       end if
!!$
!!$       !--Exit here to remove need to evaluate an extra point when maximum has been found
!!$       if(Continue_to_Evaluate == .false.) exit
!!$
!!$       Posterior(2,i) = 1.e0_double
!!$       if(Alpha_Loop == Size(Posterior,2)/2) print *, 'Approximately halfway done for this Aperture..'          
!!$       do c = 1, size(Posterior_perGalaxy,1) !-Loop over galaxies-!
!!$          !             print *, 'Doing Galaxy Loop:', c, size(Posterior_pergalaxy,1)
!!$          
!!$          !~~Select the method of posterior reconstruction for that point based in input method   
!!$          select case(Posterior_Method)
!!$          case(1) !--SizeOnly--!
!!$             if( (Cat%Sizes(c) > Survey_Size_Limits(2)) .or. (Cat%Sizes(c) < Survey_Size_Limits(1))) then
!!$                nGal_Ignored_SizeLimits = nGal_Ignored_SizeLimits + 1
!!$                cycle
!!$             end if
!!$             !MagRenorm          iSurvey_size_Limits = Survey_Size_Limits
!!$             
!!$             nSizePosterior = nSizePosterior + 1
!!$             Galaxy_Posterior_Method = 1
!!$          case(2)!--SizeMag--!
!!$             if( (Cat%Sizes(c) > Survey_Size_Limits(2)) .or. (Cat%Sizes(c) < Survey_Size_Limits(1))) then
!!$                nGal_Ignored_SizeLimits = nGal_Ignored_SizeLimits + 1
!!$                cycle
!!$             end if
!!$             
!!$             !MagRenorm          iSurvey_size_Limits = Survey_Size_Limits
!!$             
!!$             nSizeMagPosterior = nSizeMagPosterior + 1
!!$             Galaxy_Posterior_Method = 2
!!$          case(3) !-Magnitude Only--!
!!$             nMagPosterior = nMagPosterior + 1
!!$             
!!$             !MagRenorm          iSurvey_size_Limits = (/0.e0_double, 1.e30_double/) !--Should encompase the whole data set
!!$             
!!$             if(present(MagPrior) == .false.) STOP 'Attempting to produce posterior using mag only, but no magnitde prior entered, stopping'
!!$             Galaxy_Posterior_Method = 3
!!$          case(4) !-Size mag above size data limit, Mag-Only below size data limit-!
!!$             if(present(MagPrior) == .false.) STOP 'Attempting to produce posterior using mag only under size limit, but no magnitde prior entered, stopping'
!!$             STOP 'Posterior Method 4 has been disabled as it needs further thought into its construction. In particular, with respect to the construction of the prior from the size-mag distriution, the effect of prior cuts, and correct renormalisation. The skeleton code has been added for this, but disabled for now. See commented code labeled MagRenorm'
!!$             
!!$             !--Note, in this case the data-renormalisation should take inot account that the mag only is constructed from a sample which takes small, faint galaxies
!!$             if(Cat%Sizes(c) < Survey_Size_Limits(1)) then
!!$                !MagRenorm             iSurvey_Size_Limits = (/0.e0_double, Survey_Size_Limits(1)/)
!!$                
!!$                Galaxy_Posterior_Method = 3
!!$                nMagPosterior = nMagPosterior + 1
!!$             else
!!$                !MagRenorm             iSurvey_Size_Limits = Survey_Size_Limits
!!$                
!!$                Galaxy_Posterior_Method = 2
!!$                nSizeMagPosterior = nSizeMagPosterior + 1
!!$             end if
!!$          end select
!!$          
!!$          !~~~ Cycle if the observed magnitude falls outside the survey limits
!!$          if( (Cat%MF606W(c) > Survey_Magnitude_Limits(2)) .or. (Cat%MF606W(c) < Survey_Magnitude_Limits(1))) then
!!$             nGal_Ignored_MagLimits = nGal_Ignored_MagLimits + 1
!!$             cycle
!!$          end if
!!$          if(isNaN(Cat%Sizes(c)) .or. isNaN(Cat%MF606W(c))) then
!!$             nGal_Ignored_NaN = nGal_Ignored_NaN + 1
!!$             cycle
!!$          end if
!!$          
!!$          !~~~Determine method of setting source redshift: Either Marginalise over distribution; Use Redshift of source if known; Use a default redshift
!!$          
!!$          !~~~Set Redshift PDF. If redshift is known, this is not needed. Alternatively, this could be determined externally and interpolated
!!$          Known_Redshift = .false.
!!$          !-Set to NaN as default (indicates not properly set
!!$          Galaxy_Sigma_Critical = dsqrt(-1.e0_double) 
!!$          if(Cat%Redshift(c) >= 0.e0_double) then
!!$             !--Use Galaxy Redshift if available--!
!!$             !--Ignore Galaxies with redshift less than the foreground-!
!!$             if(Cat%Redshift(c) < Lens_Redshift) cycle
!!$             
!!$             Galaxy_Sigma_Critical = Linear_Interp(Cat%Redshift(c), RedshiftGrid, Sigma_Crit)!1.66492e0_double*(D_s/(D_l*D_ls))
!!$             Known_Redshift = .true.
!!$             
!!$          elseif(Marginalise_Redshift_Distribution == .false.) then
!!$             !--Use Default Redshift if everything else fails (usually this should not be considered)--!
!!$             n_Default_Source_Redshift_Used = n_Default_Source_Redshift_Used + 1
!!$             
!!$             Galaxy_Sigma_Critical = Linear_Interp(Default_Source_Redshift, RedshiftGrid, Sigma_Crit)
!!$             Known_Redshift = .true.
!!$          elseif(Marginalise_Redshift_Distribution) then
!!$             !--Marginalising over the redshift distribution for that galaxy--!
!!$          else
!!$             STOP 'DM_Profile_Variable_Posterior - Both MC and Redshift Distribution methods set - this cannot be'
!!$          end if
!!$          if(Galaxy_Sigma_Critical < 0.e0_double) STOP 'DM_Profile_Variable_Posterior - Invalid Sigma Critical Entered, negative'
!!$          
!!$          if(Known_Redshift) then
!!$             allocate(Posterior_perGalaxy_Redshift(1)); Posterior_perGalaxy_Redshift = 0.e0_double
!!$          else
!!$             allocate(Posterior_perGalaxy_Redshift(nRedshift_Sampling)); Posterior_perGalaxy_Redshift = 0.e0_double
!!$          end if
!!$
!!$          do z = 1, size(Posterior_perGalaxy_Redshift,1)             
!!$
!!$             !--There must be a better way to do this, perhaps using interpolation? (i.e. it is very seperated from other methods
!!$             if(Marginalise_Redshift_Distribution .and. (Known_Redshift == .false.)) Galaxy_Sigma_Critical = Sigma_Crit(z)
!!$             
!!$             
!!$             if(isNaN(Galaxy_Sigma_Critical) .or. (Galaxy_Sigma_Critical > huge(1.e0_double))) STOP 'DM_Profile_Variable_Posterior - FATAL - Galaxy Sigma Critical not set correctly'
!!$             select case(Mass_Profile)
!!$             case(1) !-Flat-!
!!$                STOP 'The Strong Lensing approximation has not been implemented for this profile yet'
!!$                Effective_Magnification = 1.e0_double+2.e0_double*(Posterior(1,i)/Galaxy_Sigma_Critical)
!!$             case(2) !-SIS-!
!!$                STOP 'The Strong Lensing approximation has not been implemented for this profile yet'
!!$                Effective_Magnification = 1.e0_double+2.e0_double*(SMD_SIS(Posterior(1,i), Distance_From_Mass_Center(c))/(Galaxy_Sigma_Critical*1.e18_double))
!!$             case(3) !-NFW-!
!!$                if(Enforce_Weak_Lensing) then
!!$                   Effective_Magnification = 1.e0_double + 2.e0_double*(SMD_NFW(Distance_From_Mass_Center(c), Lens_Redshift, Posterior(1,i))/(Galaxy_Sigma_Critical*1.e18_double))
!!$                else
!!$                   Effective_Magnification = Magnification_Factor(3, Distance_From_Mass_Center(c), Posterior(1,i),  Lens_Redshift, Galaxy_Sigma_Critical*1.e18_double)
!!$                end if
!!$             case default
!!$                STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
!!$             end select
!!$
!!$             if(isNaN(Effective_Magnification)) then
!!$                print *, 'Magnification Factor is a NaN:', Distance_From_Mass_Center(c), Posterior(1,i), Lens_Redshift, Galaxy_Sigma_Critical*1.e18_double
!!$                print *, 'Loops, posterior, galaxy, redshift:', i, c, z
!!$                STOP
!!$             end if
!!$
!!$             if(Effective_Magnification < minval(MagnificationGrid) .or. Effective_Magnification > maxval(MagnificationGrid)) then
!!$                !--Skipping as outside limits on which magnification was evaluated--!
!!$                Posterior_perGalaxy_Redshift(z) = 1.e-100_double
!!$                cycle
!!$             end if
!!$
!!$             if(Effective_Magnification <= 0.e0_double) then
!!$                !--This won't be picked up unless the above cycle is relaxed--!
!!$                print *, Effective_Magnification, SMD_NFW_Scalar(Distance_From_Mass_Center(c), Lens_Redshift, Posterior(1,i)), Differential_SMD_Scalar(Distance_From_Mass_Center(c), Lens_Redshift, Posterior(1,i)), Galaxy_Sigma_Critical*1.e18_double, Posterior(1,i), Distance_From_Mass_Center(c)
!!$                STOP 'Effective Magnification invalid - negative. Stopping'
!!$             end if
!!$             
!!$             !--Get RedshiftPDF which depends on the unlensed magnitude - this could possibly be interpolated, which would hopefully be faster
!!$             if(Known_Redshift) then
!!$                RedshiftPDF = 1.e0_double
!!$             else
!!$                RedshiftPDF = CH08_Redshift_Distribution_Scalar(Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification), RedshiftGrid(z))
!!$             end if
!!$
!!$             !--Construct Joint Size-Magnitde Prior which is renormalised according to the convergence value--!
!!$
!!$             Renormalisation_Size_Limits = Survey_Size_Limits/dsqrt(Effective_Magnification)
!!$             Renormalisation_Magnitude_Limits = Survey_Magnitude_Limits + 2.5e0_double*dlog10(Effective_Magnification)
!!$
!!$             if(Cuts_Renormalise_Likelihood) then
!!$                if(i == 1 .and. z ==1 .and. c == 1) print *, 'KAPPA RENORMALISATION TURNED ON'
!!$                Renorm = Linear_Interp(Effective_Magnification, MagnificationGrid, Renormalisation_by_Magnification, ExValue = 1.e30_double)
!!$                if(Renorm == 0) then
!!$                   Kappa_Renormalised_Prior = 0.e0_double
!!$                else
!!$                   Kappa_Renormalised_Prior = Survey_Renormalised_Prior/Renorm
!!$                end if
!!$                Renorm = 0.e0_double
!!$             else
!!$                if(i == 1 .and. z ==1 .and. c == 1) print *, 'KAPPA RENORMALISATION TURNED OFF'
!!$                Kappa_Renormalised_Prior =  Survey_Renormalised_Prior
!!$                Kappa_Renormalised_MagPrior = Survey_Renormalised_MagPrior
!!$             end if
!!$
!!$             select case (Galaxy_Posterior_Method)
!!$             case(1) !--Size Only--!
!!$                
!!$                !--Evaluate p_{theta_0, m_0}*p_{z|m_0} for the whole magnitude grid
!!$                do m =  1, size(PriorMagGrid)
!!$                   !--m_0, theta_0--!
!!$                   if((PriorMagGrid(m) < Renormalisation_Magnitude_Limits(1)) .or. (PriorMagGrid(m) > Renormalisation_Magnitude_Limits(2))) then
!!$                      !-Since Integrand does not extend over this region anyway-!
!!$                      Size_Only_Mag_Prior(m,:) = 1.e-100_double
!!$                   else
!!$                      if(Known_Redshift) then
!!$                         Size_Only_Mag_Prior(m,:) =  Kappa_Renormalised_Prior(m,:)
!!$                      else
!!$                         Size_Only_Mag_Prior(m,:) =  Kappa_Renormalised_Prior(m,:)*CH08_Redshift_Distribution_Scalar(PriorMagGrid(m), RedshiftGrid(z))
!!$                      end if
!!$                   end if
!!$                end do
!!$
!!$
!!$                if(lnSize_Prior) then
!!$                   STOP 'ln Size prior has been switched off since the conversion to strong lensing - check the theory carefully before use'
!!$                   if(dlog(Cat%Sizes(c)/dsqrt(Effective_Magnification)) > maxval(PriorSizeGrid) .or. (dlog(Cat%Sizes(c)/dsqrt(Effective_Magnification)) < minval(PriorSizeGrid))) then
!!$                      !--Extrapolation--!
!!$                      Posterior_perGalaxy_Redshift(z) = 0.e0_double
!!$                      cycle
!!$                   end if
!!$
!!$                   !--Find Boundary Indexs for the Size at which the Prior will be evaluated - Aims to improve run-tim by minimising the number of integrations required
!!$                   allocate(Size_Only_Prior(2)); Size_Only_Prior = 0.e0_double
!!$                   do j = 1, size(PriorSizeGrid)-1
!!$                      if( ( PriorSizeGrid(j)<= dlog(Cat%Sizes(c)/dsqrt(Effective_Magnification))) .and. (  PriorSizeGrid(j+1) > dlog(Cat%Sizes(c)/dsqrt(Effective_Magnification)))) then
!!$!                      if( (PriorSizeGrid(j)<= dlog(Cat%Sizes(c))-Effective_Convergence(c,i)) .and. ( PriorSizeGrid(j+1) > dlog(Cat%Sizes(c))-Effective_Convergence(c,i)) ) then
!!$                         Size_Only_Prior(1) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
!!$                         Size_Only_Prior(2) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j+1), 2, lim = Renormalisation_Magnitude_Limits)
!!$                         exit
!!$                      end if
!!$                   end do
!!$
!!$!!!$                   Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(dlog(Cat%Sizes(c))-Effective_Convergence(c,i), PriorSizeGrid(j:j+1), Size_Only_Prior(:), ExValue = 0.e0_double) !-Reinstate with SLensing--!
!!$                else
!!$                   if((Cat%Sizes(c)/dsqrt(Effective_Magnification) > maxval(PriorSizeGrid)) .or. (Cat%Sizes(c)/dsqrt(Effective_Magnification) < minval(PriorSizeGrid))) then
!!$                      !--Extrapolation--!
!!$                      Posterior_perGalaxy_Redshift(z) = 1.e-100_double
!!$                      cycle
!!$                   else
!!$                      
!!$                      !--Find Boundary Indexs for the Size at which the Prior will be evaluated - Aims to improve run-time by minimising the number of integrations required
!!$                      allocate(Size_Only_Prior(2)); Size_Only_Prior = 0.e0_double
!!$                      do j = 1, size(PriorSizeGrid)-1
!!$                         if( (PriorSizeGrid(j)<= Cat%Sizes(c)/dsqrt(Effective_Magnification)) .and. ( PriorSizeGrid(j+1) > Cat%Sizes(c)/dsqrt(Effective_Magnification)) ) then
!!$                            Size_Only_Prior(1) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
!!$                            Size_Only_Prior(2) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j+1), 2, lim = Renormalisation_Magnitude_Limits)
!!$                            exit
!!$                         end if
!!$                      end do
!!$
!!$                      
!!$                      Posterior_perGalaxy_Redshift(z) = Linear_Interp(Cat%Sizes(c)/dsqrt(Effective_Magnification), PriorSizeGrid(j:j+1), Size_Only_Prior(:), ExValue =  1.e-100_double)*(1.e0_double/dsqrt(Effective_Magnification))
!!$
!!$                   end if
!!$                end if
!!$                deallocate(Size_Only_Prior)
!!$             case(2)!-Size and Magnitude-!
!!$                
!!$                if(Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification) < 21.e0_double) then
!!$                   print *, 'Possible problem with de-lensed magnitude - falls outwith bright limit of Scrabback Fit (21)'
!!$                   print *, Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification)
!!$                end if
!!$                if(lnSize_Prior) then
!!$                   STOP 'DM_Profile_Variable_Posterior - lnSize with Size_Magnitude Method - I cannae do that captain!'
!!$                else
!!$                   !--Uses distributions of the apparent size--!
!!$                   !--Test for need for extrapolation---!
!!$                   need_Extrapolate = .false.
!!$                   if(( (Cat%MF606W(c)+2.5e0_double*dlog10(Effective_Magnification) > maxval(PriorMagGrid)) .or. (Cat%MF606W(c)+2.5e0_double*dlog10(Effective_Magnification) < minval(PriorMagGrid)) ) .or. ((Cat%Sizes(c)/dsqrt(Effective_Magnification) > maxval(PriorSizeGrid)) .or. (Cat%Sizes(c)/dsqrt(Effective_Magnification) > maxval(PriorSizeGrid)) )) then
!!$                      need_Extrapolate = .true.
!!$                   else
!!$                      need_Extrapolate = .false.
!!$                   end if
!!$      
!!$                   if((need_Extrapolate .and. do_KDE_Extrapolation) .or. do_KDE_OnTheFly) then !-.or. KDE_OnTheFly
!!$                      !--KDE_Extrapolation / KDE_OnTheFly(? - What about entry of prior?)
!!$                      Posterior_perGalaxy_Redshift(z) = KDE_BiVariate_Gaussian_Scalar(PriorCatalogue%MF606W, PriorCatalogue%Sizes, Cat%MF606W(c)+2.5e0_double*dlog10(Effective_Magnification), Cat%Sizes(c)/dsqrt(Effective_Magnification), Inverse_Covar = KDE_Covariance_Inverse, Det_Covar = KDE_Covariance_Determinant)*(1.e0_double/dsqrt(Effective_Magnification))*RedshiftPDF
!!$                   elseif(need_Extrapolate .and. (do_KDE_Extrapolation == .false.)) then
!!$                      !--Extrapolation, set to default, (zero)
!!$                      Posterior_perGalaxy_Redshift(z) = 1.e-100_double
!!$                   else
!!$                      !--Interpolation--!
!!$                      Posterior_perGalaxy_Redshift(z) = Linear_Interp(Cat%MF606W(c)+2.5e0_double*dlog10(Effective_Magnification), Cat%Sizes(c)/dsqrt(Effective_Magnification), PriorMagGrid, PriorSizeGrid, Kappa_Renormalised_Prior, ExValue = 1.e-100_double)*(1.e0_double/dsqrt(Effective_Magnification))*RedshiftPDF
!!$                   end if
!!$                   
!!$                   !--If no KDE Extrapolation, then this will set to a default value (effectively zero) outside the prior grid range
!!$                end if
!!$             case(3) !--Magnitude Only--!
!!$                !--Renormalise Magnitude Prior Distribution--!
!!$                
!!$                if(Cuts_Renormalise_Likelihood) then
!!$                   Renorm = Linear_Interp(Effective_Magnification, MagnificationGrid, MagOnly_Renormalisation_by_Magnification, ExValue = 1.e30_double)                      
!!$                   if(Renorm == 0) then
!!$                      Kappa_Renormalised_MagPrior = 0.e0_double
!!$                   else
!!$                      Kappa_Renormalised_MagPrior = Survey_Renormalised_MagPrior/Renorm
!!$                   end if
!!$                end if
!!$                
!!$                !--Could be edited for KDE Extrapolation, not done yet--!
!!$                
!!$                !--Construct p_[m_0] as the magnitude distribution for a sample of galaxies between cuts-corrected survey size limits 
!!$                !-MagRenorm
!!$!!!!$                   allocate(Mag_Only_Prior(2)); Mag_Only_Prior = 0.e0_double
!!$!!!$                   !--Find point where de-lensed magnitude lie on intrinsic distribution - could also use locate.
!!$!!!$                   do j = 1, size(PriorMagGrid)-1
!!$!!!$                      if( ( PriorMagGrid(j)<=  Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification) ) .and. (  PriorMagGrid(j+1) > Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification)) ) then
!!$!!!$                         Mag_Only_Prior(1) = Integrate(PriorSizeGrid, Kappa_Renormalised_Prior(j,:), 2, lim = Renormalisation_Size_Limits)
!!$!!!$                         Mag_Only_Prior(2) = Integrate(PriorSizeGrid, Kappa_Renormalised_Prior(j+1,:), 2, lim = Renormalisation_Size_Limits)
!!$!!!$                         exit
!!$!!!$                      end if
!!$!!!$                   end do
!!$!!!$                   Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification), PriorMagGrid(j:j+1), Mag_Only_Prior, ExValue =  1.e-100_double)*RedshiftPDF(z)
!!$!!!$                   deallocate(Mag_Only_Prior)
!!$                !---------------------------------------------------------------------------------------------------------------------
!!$                
!!$                !--The following gives unbiased results if no size cuts are used, however is known to be biased in the presence of size cuts, I believe that this is due to the fact that the prior itself should depend on the size cuts in the presence of a size-magnitude correlation.
!!$                
!!$                Posterior_perGalaxy_Redshift(z) = Linear_Interp(Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification), PriorMagGrid, Kappa_Renormalised_MagPrior, ExValue =  1.e-100_double)*RedshiftPDF
!!$
!!$             case default
!!$                STOP 'DM_Profile_Variable_Posterior - Error in choosing method of finding posterior'
!!$             end select
!!$          end do !--End of Redshift Loop
!!$          
!!$          if(size(Posterior_perGalaxy_Redshift,1) == 1) then
!!$             !--Redshift was taken to be exact
!!$             Posterior_perGalaxy(c) = Posterior_perGalaxy_Redshift(1)
!!$          else
!!$             !--Integrate over the redshift Information--!
!!$             Posterior_perGalaxy(c) = TrapInt(RedshiftGrid, Posterior_perGalaxy_Redshift(:))
!!$          end if
!!$          deallocate(Posterior_perGalaxy_Redshift)
!!$          
!!$       end do !--End of galaxy loop       
!!$
!!$       !~~Return lnP to ensure that we can renormalise in an alpha-independent way in the combined posterior, and to ensure that the PDF is correctly recovered even with ronding error (large negative lnP across all alphas). The renormalisation is done after the outer loop is finished
!!$       call Combine_Posteriors(Posterior(1,i), Posterior_perGalaxy(:), Combine_log_Posteriors, Renormalise = .false., Return_lnP = Combine_log_Posteriors, Combined_Posterior = Posterior(2,i))
!!$
!!$    end do !--End of Posterior Loop
!!$    deallocate(Posterior_perGalaxy)
!!$
!!$    if(any(dabs(Posterior) > huge(1.e0_double)) .or. any(isNaN(Posterior))) then
!!$       print *, 'NaNs or infinities found in posterior:'
!!$       do i =1, size(Posterior,2)
!!$          print *, Posterior(:,i)
!!$       end do
!!$       STOP
!!$    end if
!!$
!!$    !--Convert from lnP to P for posterior, renormalise so that max(P) = 1, to avoid rounding errors
!!$    if(Combine_log_Posteriors) Posterior(2,:) = dexp(Posterior(2,:) - maxval(Posterior(2,:)))
!!$
!!$    if(n_Default_Source_Redshift_Used > 0) write(*,'(A,I3,A)') '****** Used the default redshift for:', n_Default_Source_Redshift_Used, ' galaxies ********'
!!$    if(nGal_Ignored_SizeLimits > 0) write(*,'(A,I3)') '****** Number of galxies ignored as they fell outside the survey size limits:', nGal_Ignored_SizeLimits
!!$    if(nGal_Ignored_NaN > 0) write(*,'(A,I3)') '****** Number of galaxies ignored as they were NaNs:', nGal_Ignored_NaN
!!$    
!!$    
!!$    print *, '-------------------------------------------------------------'
!!$    print *, 'Constructed ', nSizePosterior, ' size-only posteriors'
!!$    print *, 'Constructed ', nSizeMagPosterior, ' size-magnitude posteriors'
!!$    print *, 'Constructed ', nMagPosterior, ' magnitude-only posteriors'
!!$    print *, '-------------------------------------------------------------'
!!$
!!$    !---DO NOT DELETE THIS. Instead, come up with a way to store the posterior per galaxy on a general grid for output
!!$!!!$    if(Debug_Mode .and. Output_Posterior_Per_Galaxy) then
!!$!!!$       
!!$!!!$       Filename = trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat'
!!$!!!$       open(unit = 17, file = Filename)
!!$!!!$       write(fmtstring,'(I7)') size(Posterior_perGalaxy,1)+1
!!$!!!$       do i =1, size(Posterior,2)
!!$!!!$          write(17, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posterior(1,i), Posterior_perGalaxy(:,i)
!!$!!!$       end do
!!$!!!$       close(17)
!!$!!!$       print *, 'Output Posterior per galaxy to: ', trim(adjustl(Filename))
!!$!!!$    end if
!!$    
!!$    
!!$    Filename= trim(adjustl(Output_Prefix))//'Posterior_Combined_Renormalised.dat'
!!$    open(unit = 82, file = Filename)
!!$    write(fmtstring,'(I1)') 2
!!$    do i =1, size(Posterior,2)
!!$       write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posterior(1,i), Posterior(2,i)
!!$    end do
!!$    close(82)
!!$    print *, 'Output Combined Posterior to: ', trim(adjustl(Filename))
!!$    
!!$
!!$    if(any(isNAN(Posterior(2,:)))) then
!!$       print *, 'Any NaNs in aperture posterior?', any(isNAN(Posterior(2,:))), count(isNAN(Posterior(2,:)) == .true.)
!!$       print *, 'Stopping'
!!$       STOP
!!$    END if
!!$    
!!$    !--Remove any information less than a tolerance - Numerical Error--!
!!$    where(Posterior(2,:) < 1.e-12_double)
!!$       Posterior(2,:) = 0.e0_double
!!$    end where
!!$    
!!$    deallocate(Distance_From_Mass_Center)
!!$    if(allocated(Kappa_Renormalised_Prior)) deallocate(Kappa_Renormalised_Prior)
!!$    if(allocated(Survey_Renormalised_Prior)) deallocate(Survey_Renormalised_Prior)
!!$    if(allocated(Size_Only_Mag_Prior)) deallocate(Size_Only_Mag_Prior)
!!$    
!!$    !--On Successful Completion delete Poster per galaxy as large file--!
!!$!!!$    inquire(file = trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat', exist = here)
!!$!!!$    if(here == .false.) STOP "Posterior per galaxy doesn't exist, stopping before accidental deletion"
!!$!!!$    call system('rm '//trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat')
!!$    
!!$  end subroutine DM_Profile_Variable_Posterior

