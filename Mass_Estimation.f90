module Mass_Estimation
  !--Contains Subroutines that return mass estimates--!
  use Param_Types
  implicit none

contains

  subroutine Mass_Estimate_Circular_Aperture_Catalogue(Cat, Profile, Aperture_Positions, Aperture_Radiuses, global_average_size)
    use Catalogues; use Statistics, only: mean_discrete, get_variance; use Cosmology, only:angular_diameter_distance_fromRedshift
    !--Aperture Positions must be in (/RA, Dec/), and Aperture Radiuses in DEGREES--!
    !--Each Galaxy must have a redshift assigned, and physical size associated with it--!
    type(Catalogue),intent(in)::Cat
    integer::Profile !-0:Flat, 1:SIS, 2:NFW-!
    real(double), intent(in)::Aperture_Positions(:,:), Aperture_Radiuses(:)
    real(double),intent(in),optional::Global_Average_Size

    !--Internal Declarations--!
    real(double),allocatable::iAperture_Radiuses(:)
    real(double)::field_mean_size
    integer:: Ap, Gal
    real(double),allocatable::Galaxy_Convergence(:)
    real(double)::Convergence_Error, Mean_Convergence

    real(double)::Projected_Mass(size(Aperture_Positions,1)), Error_Projected_Mass(size(Aperture_Positions,1))

    real(double),allocatable::Weight(:,:), Sigma_Crit(:), Area(:), Renormalisation(:,:)
    real(double)::D_l, D_ls , D_s, Lens_Redshift = 0.165e0_double, Default_Source_Redshift = 1.e0_double, Source_Redshift

    !---TESTING DECLARATIONS--!
    integer::nNoRedshift

    allocate(iAperture_Radiuses(size(Aperture_Positions,1)))
    if(size(Aperture_Radiuses) == 1) then
       iAperture_Radiuses = Aperture_Radiuses(1)
    elseif(size(Aperture_Radiuses) /= size(iAperture_Radiuses)) then
       STOP 'Mass_Estimate_Circular_Aperture_Catalogue - FATAL ERROR - Error assigning Aperture_Radiuses internal'
    else
       iAperture_Radiuses = Aperture_Radiuses
    end if

    if(any( (/0,1,2/) == Profile)==.false.) STOP 'Mass_Estimate_Circular_Aperture_Catalogue - FATAL ERROR - Supported profiles are 0:Flat, 1:SIS, 2:NFW.'
    
    if(present(global_average_size)) then
       field_mean_size = global_average_size
    else
       print *, 'Calculating the mass in circular aperture using mean physical size of the full catalogue'
       field_mean_size = mean_discrete(Cat%Physical_Sizes)
    end if
    
    !--Get Convergence for each galaxy--!
    allocate(Galaxy_Convergence(size(Cat%Physical_Sizes))); Galaxy_Convergence = 0.e0_double
    do Gal = 1, size(Cat%Physical_Sizes)
       Galaxy_Convergence(Gal) = (Cat%Physical_Sizes(Gal)/field_mean_size) - 1.e0_double
    end do

    Mean_Convergence = mean_discrete(Galaxy_Convergence)
    Convergence_Error = get_variance(Galaxy_Convergence, Galaxy_Convergence, Mean_Convergence, Mean_Convergence)

    print *, 'Mean/Varaince [Convergence]:', Mean_Convergence, dsqrt(Convergence_Error)

    !--Calculate Sigma_Critical for each galaxy, and area in each aperture--!
    nNoRedshift = 0
    allocate(Sigma_Crit(size(Cat%Redshift)))
    D_l = angular_diameter_distance_fromRedshift(Lens_Redshift) !-In units of Mpc/h                                                                                                                                                         
    do Gal = 1, size(Cat%Redshift) 
       Source_Redshift = Cat%Redshift(Gal)
       if(Cat%Redshift(Gal) < 0.e0_double) then
          nNoRedshift = nNoRedshift + 1
          Source_Redshift = Default_Source_Redshift
       end if
          
       D_s = angular_diameter_distance_fromRedshift(Source_Redshift)

       D_ls = D_s-D_l
       Sigma_Crit(Gal) = 1.66492e0_double*(D_s/(D_l*D_ls)) !-(10^18 M_Sun/h)-!      
    end do
    allocate(Area(size(Aperture_Positions,1))); Area = 0.e0_double
    Area = 3.142e0_double*( (D_l*1.746e-2_double*iAperture_Radiuses)**2.e0_double )
    
    print *, 'Mass_Estimate_Circular_Aperture_Catalogue - ', nNoRedshift, ' galaxies were assigned the default redshift information of ', Default_Source_Redshift, ' as no redshift information was input'

    allocate(Weight(size(Aperture_Positions,1),size(Cat%Physical_Sizes))); Weight = 0.e0_double
    allocate(Renormalisation(size(Aperture_Positions,1),size(Cat%Physical_Sizes))); Renormalisation = 0.e0_double
    Projected_Mass = 0.e0_double
    do Ap = 1, size(Aperture_Positions,1)
       do Gal = 1, size(Cat%Physical_Sizes)
          if( distance_between_points( (/Cat%RA(Gal), Cat%Dec(Gal)/), Aperture_Positions(Ap,:) ) > iAperture_Radiuses(Ap) ) cycle
          select case(Profile)
          case(0)
             !--Flat Mass Profile--! !---CHECK THIS MATHS--!
             Projected_Mass(Ap) = Projected_Mass(Ap) + (Galaxy_Convergence(Gal)/Sigma_Crit(Gal))
             Renormalisation(Ap,Gal) = 1.e0_double/(Sigma_Crit(Gal)*Sigma_Crit(Gal))
             Weight(Ap,Gal) = 1.e0_double/(Sigma_Crit(Gal))
          case default
             STOP 'Mass_Estimate_Circular_Aperture_Catalogue - FATAL ERROR - Supported profiles are 0:Flat'
          end select
       end do

       print *, 'Done Aperture:', Ap, ', of', size(Aperture_Positions,1)
    end do
    
    !--Renormalise--!
    select case(Profile)
    case(0)
       do Ap = 1, size(Projected_Mass)
          if(sum(Renormalisation(Ap,:)) /= 0.e0_double) then
             Projected_Mass(Ap) =  Projected_Mass(Ap)/sum(Renormalisation(Ap,:))
             Error_Projected_Mass(Ap) = Area(Ap)*dsqrt(Convergence_Error)/sum(Renormalisation(Ap,:))
          else
             STOP 'Mass_Estimate_Circular_Aperture_Catalogue - Renormalisation (and Error) for Aperture', Ap,' Renormalisation is zero'
          end if
       end do
    end select

    !--Convert from projected Surface Mass Density to Mass--!
    print *, 'Multiplying by Area:'
    Projected_Mass = Area*Projected_Mass

    do Ap = 1, size(Projected_Mass)
       print *, 'Mass for cluster:', Ap, ' is:', Projected_Mass(Ap), '+-', Error_Projected_Mass(Ap)
    end do

  end subroutine Mass_Estimate_Circular_Aperture_Catalogue
    

  subroutine Mass_Estimate_CircularAperture(Convergence_Map, Error_Convergence, x, y, Aperture_Positions, Aperture_Radiuses, Masses, Error_Masses, Source_Redshift, Convergence_Map_Occupation)
    use cosmology
    !-Returns the mass estimate with a circular Aperture--!
    !-Masses are calculated in a model independant way, as detailed in Heymasn et al 2008. Note that in that anaylsis, radius of 0.75' is used-!
    !-Convergence_Map, x, y are considered to be the same struture as those defined in Convergence_Estimates - that is x(i) defines the lower x values for bin i, x(i+1) the upper x value for bin i. ALSO ASSUMED TO BE EQUALLY BINNED

    !--x and y are in DEGREES. Radius should also be consistent (in DEGRESS).  **** dx = dy not necessary, provided it is accounted for on Apix  ****

    !-NOTE IF APERTUE GOES OUTSIDE X/Y RANGE, RESULT WILL NOT BE TRUSTWORTHY-!
    !-Source and Lens Redshift are hardwired. This would need generalised for code which has a different source redshift for each aperature-!
    !-Error_Convergence is used to measure the error on the mass estimates, and must be entered as an RMS value (i.e. sigma, not sigma^2)-!

    !---NOTE: THERE MAY BE A BUG HERE, AS TAKING THE AREA OUT OF THE SUM AND MULTIPLYING BY THE SUMMED TOTAL AREA GIVES DIFFERENT RESULTS. 15 NOV 2013

    real(double),intent(in)::Convergence_Map(:,:), Error_Convergence(:,:), x(:), y(:) !-x is 1st dimension, y is 2nd dimension-!
    real(double), intent(in)::Aperture_Positions(:,:) !-Labels centers of apertures, D1: Center Label, D2: 1: x, 2:y--!
    real(double), intent(in)::Aperture_Radiuses(:) !-Must have same size of D1 as Positions or 1, if 1 then same radius is applied for all apertures-!
    real(double), intent(out),allocatable::Masses(:), Error_Masses(:)
    real(double),intent(in),optional::Source_Redshift
    integer, intent(in),optional::Convergence_Map_Occupation(:,:)

    integer::i, j, x_index, y_index, x_rad_index, y_rad_index, x_index_range, y_index_range
    real(double)::Sigma_Crit = 1.e0_double
    real(double)::iSource_Redshift, Lens_Redshift = 0.165e0_double
    real(double)::D_l, D_s, D_ls
    real(double)::A_pix

    real(double),dimension(size(Aperture_Positions,1))::iAperture_Radiuses
    real(double)::nRandomPoints = 1000, NoRandomPoints_Success
    integer:: nPixel, nFullPixel, nEmptyPixel, nPartialPixel !-Testing - contains information on the number of pixels fully in aperture

    real(double)::Total_Area_Enclosed(size(Aperture_Positions,1))

    !-Random Number declarations-!
    real(double),dimension(:),allocatable::Ran(:,:)
    Integer,allocatable::seed(:)
    integer::Nseed, Clock

    real(double),dimension(size(Aperture_Positions,1)):: Sum_Convergence_in_Aperture, NPix_in_Aperture !-Used for testing only. Sums the convergence within the aperture, and the total fraction of pixels in aperture-!

    logical::Inverse_Variance_Weight = .false.
    real(double)::Sum_Weight
    real(double)::Sum_Convergence

    INTERFACE
       subroutine Mass_Estimate_CircularAperture(Convergence_Map, Error_Convergence, x, y, Aperture_Positions, Aperture_Radiuses, Masses, Error_Masses, Source_Redshift, Convergence_Map_Occupation)
         use Param_Types
         real(double),intent(in)::Convergence_Map(:,:), Error_Convergence(:,:), x(:), y(:) !-x is 1st dimension, y is 2nd dimension-!
         real(double), intent(in)::Aperture_Positions(:,:) !-Labels centers of apertures, D1: Center Label, D2: 1: x, 2:y--!              
         real(double), intent(in)::Aperture_Radiuses(:) !-Must have same size of D1 as Positions or 1, if 1 then same radius is applied for all apertures-!
         real(double), intent(out),allocatable::Masses(:), Error_Masses(:)
         
         real(double),intent(in),optional::Source_Redshift
         integer, intent(in),optional::Convergence_Map_Occupation(:,:)
       end subroutine Mass_Estimate_CircularAperture
    END INTERFACE

    print *, 'Estimating Masses:....'

    if(any(isNaN(Convergence_Map))) STOP 'Mass_Estimate_CircularAperture - FATAL ERROR - Convergence Map contains NaNs'
    
    iSource_Redshift = 1.4e0_double
    if(present(Source_Redshift)) iSource_Redshift = Source_Redshift
    print *, 'Estimating Mass with sources at z =', Source_Redshift, 'and lens at z = ', Lens_Redshift

    !--Input Error Catching--!
    if(size(x) /= size(Convergence_Map,1)) then
       print *, size(x), size(Convergence_Map,1)
       STOP 'Error - Mass_Estimate_CircularAperture - x and Kappa Map are not conformal, stopping..'
    end if
    if(size(y) /= size(Convergence_Map,2)) STOP'Error - Mass_Estimate_CircularAperture - x and Kappa Map are not conformal, stopping..'
    if((size(Aperture_Radiuses) /= 1) .and. (size(Aperture_Radiuses) /= size(Aperture_Positions,1))) then
       print *, 'Error - Mass_Estimate_CircularAperture - Aperture_Radiuses is not of the correct size, stopping..'
       print *, size(Aperture_Radiuses), size(Aperture_Positions,1)
       stop
    end if
    if(size(Aperture_Positions,2) /= 2) STOP 'Error - Mass_Estimate_CircularAperture - Aperture_Positions must contain only 2 positions (x,y) for each aperture, stopping..'


    if(size(Aperture_Radiuses)==1) then
       iAperture_Radiuses = Aperture_Radiuses(1)
    else
       iAperture_Radiuses = Aperture_Radiuses
    end if

    !-Determine Sigma_Critical-!
    D_l = angular_diameter_distance_fromRedshift(Lens_Redshift) !-In units of Mpc/h
    D_s = angular_diameter_distance_fromRedshift(iSource_Redshift)
    !-Convert from comoving distances to physical distances-!
    !D_l = D_l/(1.e0_double+Lens_Redshift)
    !D_s = D_s/(1.e0_double+Source_Redshift)

    D_ls = D_s-D_l
    Sigma_Crit = 1.66492e0_double*(D_s/(D_l*D_ls)) !-(10^18 M_Sun/h)-!

    print *, 'Calculating Masses with D_l = ', D_l
    print *, 'and Sigma_Crit = ', Sigma_Crit

    if(Verbose) print *, 'Sigma_Critical = ', Sigma_Crit
!!$    print *, 'Dl = ', D_l
!!$read(*,*)

    allocate(Masses(size(Aperture_Positions,1))); Masses = 0.e0_double
    allocate(Error_Masses(size(Aperture_Positions,1))); Error_Masses = 0.e0_double

    Sum_Convergence_in_Aperture = 0.e0_double; NPix_in_Aperture = 0.e0_double; Total_Area_Enclosed = 0.e0_double
    do i =1, size(Masses)
       if(size(Convergence_Map) == 1) then
          !-Single convergence estimate-!
          Masses(i) = 3.142e0_double*( (D_l*1.746e-2_double*Aperture_Radiuses(i))**2.e0_double )* Sigma_Crit * Convergence_Map(1,1)
          cycle
       end if

       Sum_Weight = 0.e0_double
       if( (Aperture_Positions(i,1)-iAperture_Radiuses(i) < x(1)) .or. (Aperture_Positions(i,1)+iAperture_Radiuses(i) > x(size(x))) ) print *, 'Warning - Aperture:', i, ' falls outside x range', iAperture_Radiuses(i), x(1), Aperture_Positions(i,1)-iAperture_Radiuses(i),' :',  x(size(x)), Aperture_Positions(i,1)+iAperture_Radiuses(i)
       if( (Aperture_Positions(i,2)-iAperture_Radiuses(i)< y(1)) .or. (Aperture_Positions(i,2)+iAperture_Radiuses(i) > y(size(y))) ) print *, 'Warning - Aperture:', i, ' falls outside y range'

       !--Find x and y positions of the center of the aperture--!
       do x_index = 1, size(x)-1
          if( (x(x_index) <= Aperture_Positions(i,1)) .and. (x(x_index+1) > Aperture_Positions(i,1)) ) then
             exit
          end if
       end do
       do y_index = 1, size(y)-1
          if( (y(y_index) <= Aperture_Positions(i,2)) .and. (y(y_index+1) > Aperture_Positions(i,2)) ) then
             exit
          end if
       end do

       !--Find the number of index points that x and y must cover (defines the boundary of the aperture)--!
       !--ASSUMES EQUAL BINNING IN X, AND ALSO IN Y--!
       x_index_range = int( (iAperture_Radiuses(i)/(x(2)-x(1)) )+2)
       y_index_range = int( (iAperture_Radiuses(i)/(y(2)-y(1)) )+2) !-+2 takes into account that x_index and y_index may not lie on apex

       !--Loop over all pixels within box defined by x/y_index_range--!
       nPixel = 0; nFullPixel = 0; nPartialPixel = 0; nEmptyPixel = 0
       do x_rad_index = maxval((/x_index-x_index_range,1/)), minval((/x_index+x_index_range,size(x)-1/)),1
          do y_rad_index = maxval((/y_index-y_index_range,1/)),  minval((/y_index+y_index_range,size(y)-1/)),1

             nPixel = nPixel + 1
             !--x/y_rad_index now loop over indexs of points within box defined by radius of aperture--!
             !-Calculate A_pix - MUST BE IN CORRECT UNITS, AND ACCOUNT FOR THE UNITS IN WHICH X AND Y ARE ENTERED. WHEN ==1, IT IS IN A PIXEL SCALE-!
             A_pix = (D_l*D_l*(x(x_rad_index+1)-x(x_rad_index))*(y(y_rad_index+1)-y(y_rad_index))*(6.283185e0_double/(360.e0_double))*(6.283185e0_double/(360.e0_double))) !-in (MPC/h)^2. Uses the fact that x and y are in DEGREES, and the small angle formula to convert into distance using steradians


             !--Throw down a set of random points within the pixel, test how many are in the aperture to  estimate the amount of the pixel to include in the mass estimate
             !----Construct Random Numbers----!
             allocate(Ran(2,nint(NRandomPoints))); Ran = 0.e0_double
             call RANDOM_SEED(size = NSeed)
             allocate(Seed(NSeed))
             call SYSTEM_CLOCK(COUNT = Clock)
             seed = Clock + (/ (i-1,i=1,NSeed) /)
             call RANDOM_SEED(PUT = seed)
             deallocate(Seed); NSeed = 0; Clock = 0
             
             call RANDOM_NUMBER(Ran)

             !-Convert Ran to  random positions within pixel-!
             Ran(1,:) = Ran(1,:)*(x(x_rad_index+1)-x(x_rad_index)) + x(x_rad_index) !-x-!
             Ran(2,:) = Ran(2,:)*(y(y_rad_index+1)-y(y_rad_index)) + y(y_rad_index) !-y-!

             NoRandomPoints_Success = 0
             do j = 1, nRandomPoints
                if(distance_between_points(Ran(:,j), Aperture_Positions(i,:)) <= iAperture_Radiuses(i)) then
                   NoRandomPoints_Success = NoRandomPoints_Success + 1
                end if
             end do
             if(NoRandomPoints_Success == 0) then 
                nEmptyPixel = nEmptyPixel + 1
             elseif(NoRandomPoints_Success==nRandomPoints) then
                nFullPixel = nFullPixel + 1
             else
                nPartialPixel = nPartialPixel + 1
             end if

             Total_Area_Enclosed(i) =  Total_Area_Enclosed(i) + A_pix*(1.e0_double*NoRandomPoints_Success/nRandomPoints)
             if(Inverse_Variance_Weight) then
                   !--Using the occupation numbers as a weight
                if(present(Convergence_Map_Occupation)) then
                   if(i == 1 .and. x_rad_index == maxval((/x_index-x_index_range,1/)) .and. y_rad_index == maxval((/y_index-y_index_range,1/)) ) print *, 'Getting masses by summing using occupation of pixel as weight'
                   if(Convergence_Map_Occupation(x_rad_index,y_rad_index) > 0) print *, Convergence_Map(x_rad_index,y_rad_index), Convergence_Map_Occupation(x_rad_index,y_rad_index)

                   Masses(i) = Masses(i) + (A_pix*(1.e0_double*NoRandomPoints_Success/nRandomPoints)*Sigma_Crit)*(Convergence_Map(x_rad_index,y_rad_index)*Convergence_Map_Occupation(x_rad_index,y_rad_index))
                   Sum_Convergence_in_Aperture(i) = Sum_Convergence_in_Aperture(i) + (Convergence_Map(x_rad_index,y_rad_index)*Convergence_Map_Occupation(x_rad_index,y_rad_index))
                   Sum_Weight = Sum_Weight + Convergence_Map_Occupation(x_rad_index,y_rad_index)
                else
                   if(i == 1 .and. x_rad_index == maxval((/x_index-x_index_range,1/)) .and. y_rad_index == maxval((/y_index-y_index_range,1/)) ) print *, 'Getting masses by summing using inverse error as weight'
                   if(Error_Convergence(x_rad_index,y_rad_index) <= 0.e0_double) cycle
                   !--This bit is inverse variance weighting using the error bars. This is unsuccessful (possibly due to error bars being incorrect?).
                   Masses(i) = Masses(i) + (A_pix*(1.e0_double*NoRandomPoints_Success/nRandomPoints)*Sigma_Crit)*(Convergence_Map(x_rad_index,y_rad_index)/(Error_Convergence(x_rad_index,y_rad_index)*Error_Convergence(x_rad_index,y_rad_index)))
                   Sum_Convergence_in_Aperture(i) = Sum_Convergence_in_Aperture(i) + (Convergence_Map(x_rad_index,y_rad_index)/(Error_Convergence(x_rad_index,y_rad_index)*Error_Convergence(x_rad_index,y_rad_index)))
                   Sum_Weight = Sum_Weight + 1.e0_double/((Error_Convergence(x_rad_index,y_rad_index)*Error_Convergence(x_rad_index,y_rad_index)))
                end if
                !--ERROR??--!
                if(i == 1 .and. x_rad_index == maxval((/x_index-x_index_range,1/)) .and. y_rad_index == maxval((/y_index-y_index_range,1/)) ) print *, 'I am not calculating errors on inverse varaince weighting'
             else
                !Masses(i) = Masses(i) + A_pix*(1.e0_double*NoRandomPoints_Success/nRandomPoints)*Sigma_Crit*Convergence_Map(x_rad_index,y_rad_index)
                !if(Convergence_Map_Occupation(x_rad_index,y_rad_index) > 0) 
                print *, Convergence_Map(x_rad_index,y_rad_index), Convergence_Map_Occupation(x_rad_index,y_rad_index), A_pix
                Masses(i) = Masses(i) +  A_pix*((1.e0_double*NoRandomPoints_Success)/nRandomPoints)*Sigma_Crit*Convergence_Map(x_rad_index,y_rad_index)
                Sum_Convergence_in_Aperture(i) = Sum_Convergence_in_Aperture(i) + Convergence_Map(x_rad_index,y_rad_index); NPix_in_Aperture(i) = NPix_in_Aperture(i) + (1.e0_double*NoRandomPoints_Success/nRandomPoints)
                !-Add Errors in Quadrature - Error_Masses is Sig^2 until sqrt taken later-!
                Error_Masses(i) = Error_Masses(i) + (A_pix*(1.e0_double*NoRandomPoints_Success/nRandomPoints)*Sigma_Crit*Error_Convergence(x_rad_index,y_rad_index))**2.e0_double
             end if
             deallocate(Ran)
          END do
       END do
!       Masses(i) = Masses(i) * Total_Area_Enclosed(i)
       print *, 'nFullPixel:', nFullPixel, nEmptyPixel, nPartialPixel, nPixel

       if(Inverse_Variance_Weight) then 
          Masses(i) = Masses(i) / Sum_Weight
          Sum_Convergence_in_Aperture(i) = Sum_Convergence_in_Aperture(i)/Sum_Weight
       end if

       print *, 'Done Mass for:', i, ' reading'
       print *, 'Summed Convergence = ', Sum_Convergence_in_Aperture(i)
       print *, 'Expected Mass:', ( 3.142e0_double*(D_l*1.746e-2_double*iAperture_Radiuses(i))**2.e0_double )*Sigma_Crit* Sum_Convergence_in_Aperture(i), ' Measured:',  Masses(i), ' Ratio:', ( 3.142e0_double*(D_l*1.746e-2_double*iAperture_Radiuses(i))**2.e0_double )*Sigma_Crit* Sum_Convergence_in_Aperture(i)/Masses(i)
       read(*,*)
    end do

    if(verbose) then
       print *, 'Total Area enclosed in Apertures is:', Total_Area_Enclosed
       print *, 'Which is total fraction of the expected:', Total_Area_Enclosed/( 3.142e0_double*(D_l*1.746e-2_double*iAperture_Radiuses)**2.e0_double )
    end if


!!$    print *, 'Sum of Convergence within ap is:', Sum_Convergence_in_Aperture
!!$    print *, 'With NPix in each aperture of:', NPix_in_Aperture
!!$    read(*,*)

    !--Convert Error_Masses into RMS (sigma)--!
    if(any(Error_Masses < 0.e0_double)) STOP 'Mass_Estimate_CircularAperture - Error on Masses is imaginary, stopping..'
    Error_Masses = dsqrt(Error_Masses)

    print *, 'Finished.'

  end subroutine Mass_Estimate_CircularAperture

  function distance_between_points( C1, C2 )
    !-Calculates the straight line distance between grid points, were C1 and C2 are (x,y) pairs-!
    real(double),intent(in)::C1(2), C2(2)

    real(double)::distance_between_points

    distance_between_points = dsqrt( ((C1(1)-C2(1))**2.) + ((C1(2)-C2(2))**2.) )

  end function distance_between_points


end module Mass_Estimation
