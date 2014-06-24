!--#Version that uses r200 defined as the radius on which density contained in halo is 200rho_crit - See Wright and Brainerd
module Mass_Profiles
  use Param_Types
  implicit none

  INTERFACE SMD_SIS
     module procedure SMD_SIS_Array, SMD_SIS_Scalar
  END INTERFACE SMD_SIS
  INTERFACE SMD_NFW
     module procedure SMD_NFW_Array, SMD_NFW_Scalar
  END INTERFACE SMD_NFW
  
  INTERFACE get_NFW_VirialRadius_from_VirialMass
     module procedure get_NFW_VirialRadius_from_VirialMass_Scalar
  END INTERFACE

  contains

    !--------------------------------GENERAL ROUTINES---------------------------------------------------------!
    real(double) function Magnification_Factor(Profile, Radius, Param, Redshift, Sigma_Critical)
      integer, intent(in):: Profile
      real(double), intent(in):: Param
      real(double), intent(in):: Radius, Redshift, Sigma_Critical

      real(double):: Gamma, Kappa

      select case(Profile)
      case(1)
         STOP 'I cannot calculate the Magnification Factor for this profile yet'
      case(2)
         STOP 'I cannot calculate the Magnification Factor for this profile yet'
      case(3) !-NFW-!
         Kappa = SMD_NFW_Scalar(Radius, Redshift, Param)/Sigma_Critical
         Gamma = Differential_SMD_Scalar(Radius, Redshift, Param)/Sigma_Critical
      end select

      Magnification_Factor = 1.e0_double/( ((1.e0_double-Kappa)**2.e0_double) - Gamma*Gamma)

    end function Magnification_Factor

    subroutine Halo_Mass(Profile, Param, Param_Error, Integrated_Mass, Mass_Error, Redshift)
      !--Gets Virial Radius for SMD type and returns Halo Mass--!
      real(double), intent(in):: Param, Param_Error(:)
      real(double), intent(out):: Integrated_Mass
      real(double),intent(out):: Mass_Error(:)
      real(double), intent(in),optional::Redshift
      integer::Profile
      
      real(double):: Virial_Radius

      Virial_Radius = virial_Radius_from_ProfileFreeParameter(Profile, Param)
      
      if(present(Redshift)) then
         call Integrated_Mass_Within_Radius(Profile, Virial_Radius, Param, Param_Error, Integrated_Mass, Mass_Error, Redshift)
      else
         call Integrated_Mass_Within_Radius(Profile, Virial_Radius, Param, Param_Error, Integrated_Mass, Mass_Error)
      end if

    end subroutine Halo_Mass

    real(double) function virial_Radius_from_ProfileFreeParameter(Profile, Param)
      integer, intent(in):: Profile
      real(double), intent(in):: Param

      select case(Profile)
      case(1) !-Flat-!
         STOP 'virial_Radius_from_ProfileFreeParameter - I cannot do FLAT yet'
      case(2) !--SIS--!
         STOP 'virial_Radius_from_ProfileFreeParameter - I cannot do SIS yet'
      case(3) !--NFW, r200 only--!
         virial_Radius_from_ProfileFreeParameter = Param
      case default
         STOP 'virial_Radius_from_ProfileFreeParameter - Incorrect profile entered'
      end select

    end function virial_Radius_from_ProfileFreeParameter

    subroutine Integrated_Mass_Within_Radius(Profile, Scale, Param, Param_Error, Integrated_Mass, Mass_Error, Redshift)
      !-Param is the free parameter for the profile: Flat->Sigma_0; SIS->Sigma_v^2; 3: NFW -> r_200, (c fixed)
      !--Scale must be in Mpc/h --!                                                                        
      !-Assumes Spherically Symmetric-!
      real(double), intent(in)::Scale, Param, Param_Error(:)
      real(double), intent(out):: Integrated_Mass
      real(double),intent(out):: Mass_Error(:)
      real(double), intent(in),optional::Redshift
      integer::Profile

      integer::i
      real(double)::iScale

      if(size(Param_Error) /= size(Mass_Error)) STOP 'Integrated_Mass_Within_Radius - Param Error and Mass Error MUST be the same size'

      if(Scale <= 0.0e0_double) then
         !--Use Virial Radius--!
         iScale = virial_Radius_from_ProfileFreeParameter(Profile, Param)
      else
         iScale = Scale
      end if

      select case(Profile)
      case(1) !-Flat-!
         Integrated_Mass = Flat_Mass_withinRadius(iScale, Param)
         do i = 1, size(Param_Error)
            Mass_Error(i) = Error_Flat_Mass_withinRadius(iScale, Param_Error(i))
         end do
      case(2) !-SIS-!
         Integrated_Mass = SIS_Mass_withinRadius(iScale, Param)
         do i = 1, size(Param_Error)
            Mass_Error(i) = SIS_Mass_withinRadius(iScale, Param_Error(i))
         end do
      case(3) !-NFW-!
         if(present(Redshift) == .false.) STOP 'Integrated_Mass_Within_Radius - Cluster Redshift MUST be present'
         Integrated_Mass = NFW_Mass_withinRadius(iScale, Redshift,r_200 = Param)
         do i = 1, size(Param_Error)
            Mass_Error(i) = Error_NFW_Mass_withinRadius(iScale, Param, Param_Error(i), Redshift)
         end do
      case default 
         STOP 'Integrated_Mass_Within_Radius - Incorrect profile entered'
      end select

    end subroutine Integrated_Mass_Within_Radius



    !------------------------------------NFW------------------------------------------------------------------!
    real(double) function Error_NFW_Mass_withinRadius(R, r200, Er200, z, r_s)
      !--Returns the error on the integrated NFW mass within a radius--!
      !--Uses an analytic for for dDelcdr200 given by Wolfram Alpha, but could alos be done analytically--!
      !--Er200 is the measured error on r200-!
      !--TESTED with values 18Mar2014--!
      real(double), intent(in):: R, r200, Er200, z
      real(double), intent(in), optional::r_s

      real(double)::rho_c, M200, delta_c, c, rs

      real(double):: dDelcdr200, dMdr200, MR

      if(present(r_s)) then
         STOP 'Error_NFW_Mass_withinRadius - I cannot deal with rs yet'
      end if

      call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c)
      rs = r200/c

      !--Following is analytic form for the derivative, taken from Wolfram Alpha--!
      dDelcdr200 = (200.e0_double/3.e0_double)*((c*c)/rs)*( (3.e0_double*((1.e0_double+c)**2.e0_double)*dlog(1.e0_double+c) - c*(3.e0_double+4.e0_double*c))/((c-(1.e0_double+c)*dlog(1.e0_double+c))**2.e0_double) )

      MR = NFW_Mass_withinRadius(R, z, r_200 = r200)
      dMdr200 = (MR/delta_c)*dDelcdr200 !--Isolates the parts of the integrated mass which do not depend on delta_c, and therefore r200--!

      Error_NFW_Mass_withinRadius = dMdr200*Er200
           
  end function Error_NFW_Mass_withinRadius

    real(double) function NFW_Mass_withinRadius(R, z,r_200, M_200, r_s)
      !--Returns the mass contained within radius R--!
      !--TESTED with values 18Mar2014--!
      !--Can be used to obtain M200 and eM200 directly by entering R = r200--!
      real(double), intent(in):: R, z
      real(double), intent(in),optional:: M_200, r_200
      real(double), intent(in), optional::r_s

      real(double)::rho_c, M200, delta_c, c, rs, r200

      INTERFACE
         real(double) function NFW_Mass_withinRadius(R, z,r_200, M_200, r_s)
           use Param_Types
           real(double), intent(in):: R, z

           real(double), intent(in),optional:: M_200, r_200
           real(double), intent(in), optional::r_s
         end function NFW_Mass_withinRadius
      END INTERFACE

      if(present(r_200) .and. present(M_200)) then
         STOP 'NFW_Mass_withinRadius - ONLY ONE OF virial radius or virial mass must be entered, stopping...'
      elseif(present(r_200)) then
         r200 = r_200
      elseif(present(m_200)) then
         r200 = get_NFW_VirialRadius_from_VirialMass(z, m_200)
      else
         STOP 'NFW_Mass_withinRadius - Either virial radius or virial mass must be entered, stopping...'
      END if

      if(r200 == 0.e0_double) then
         NFW_Mass_withinRadius = 0.e0_double
         return
      end if

      if(present(r_s)) then
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c, r_s)
         rs = r_s
      else
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c)
         rs = r200/c
      end if
      !Check entered M200 against returned?!

      NFW_Mass_withinRadius = 4.e0_double*3.142e0_double*delta_c*rho_c*rs*rs*rs*(dlog(1.e0_double +(r/rs)) - ((r/rs)/(1.e0_double+(r/rs))) )


    end function NFW_Mass_withinRadius

    subroutine rescale_NFW_byMass(M_scale, r_scale, z, r200)
      !-M_Scale is the mass wanted within the scale r_scale
      !--Returns r200 so that the mass contained within a scale r_scale is equal to M_scale
      !--Assumes r200 is the only free parameter, and concentration is set by mass (fit) so that rs is set--!
      use Interpolaters, only:Linear_Interp
      real(double), intent(in)::M_scale, r_scale, z
      real(double),intent(out):: r200
      
      integer::i
      real(double), dimension(1000)::Mass_Storage, r200Grid
      real(double)::r200_l = 0.1e0_double, r200_u = 20.e0_double, dr200

      real(double)::rho_c, M200, delta_c, c, rs

      !--Loop over reasonable r200 values until the mass within the scale is as expected--!
      dr200 = (r200_u-r200_l)/size(Mass_Storage)
      do i = 1, size(Mass_Storage)
         r200Grid(i) = r200_l+(i-1)*dr200
         call get_NFW_Parameters(z, r200Grid(i), rho_c, M200, delta_c, c)
         rs = r200Grid(i)/c
         Mass_Storage(i) = NFW_Mass_withinRadius(R_Scale, Z, r_200 = r200Grid(i), r_s = rs)
      end do      

      if(M_Scale > maxval(Mass_Storage)) then
         print *, 'Mass_Scale, Max Tabulated:', M_Scale, maxval(Mass_Storage)
         STOP 'rescale_NFW_byMass - Mass Scale too large for tabulated values'
      end if
      if(M_Scale < minval(Mass_Storage)) STOP 'rescale_NFW_byMass - Mass Scale too small for tabulated values'

      r200 = Linear_Interp(M_Scale, Mass_Storage,r200Grid)

      !--Output--!
!!$      open(unit = 62, file = 'NFW_rescale_byr200.dat')
!!$      do i = 1, size(Mass_Storage)
!!$         write(62, *) r200Grid(i), Mass_Storage(i)
!!$      end do
!!$      close(62)

    end subroutine rescale_NFW_byMass

    function SMD_NFW_Array(r, z, r200, r_s)
      real(double),intent(in)::r(:), r200, z
      real(double), intent(in),optional::r_s
      
      real(double),dimension(size(r))::SMD_NFW_Array
      
      integer::i

      do i = 1, size(r)
         if(present(r_s)) then
            SMD_NFW_Array(i) = SMD_NFW_Scalar(r(i), z, r200, r_s)
         else
            SMD_NFW_Array(i) = SMD_NFW_Scalar(r(i), z, r200)
         end if
      end do

    end function SMD_NFW_Array

    real(double) function SMD_NFW_Scalar(r, z, r200, r_s)
      real(double),intent(in)::r, r200, z
      real(double), intent(in),optional::r_s

      real(double)::delta_c, c, rs, M200
      real(double)::omega_matter !-Evaluated at the redshift of the Halo-!
      real(double)::rho_c
      real(double)::x

      if(present(r_s)) then
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c,  r_s)
      else
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c)
      end if

      if(r200 == 0.e0_double) then
         SMD_NFW_Scalar = 0.e0_double
      end if

      rs = r200/c

      x = r/rs
      !--Following from Wright, Brainerd--!
      if(x < 1) then
         SMD_NFW_Scalar = ( (2.e0_double*rs*delta_c*rho_c)/(x*x - 1.e0_double) )*(1.e0_double - (2.e0_double/dsqrt(1.e0_double-x*x))*atanh(dsqrt((1.e0_double-x)/(1.e0_double+x))) )
      elseif(x == 1) then
         SMD_NFW_Scalar = (2.e0_double*rs*delta_c*rho_c)/3.e0_double
      else
         SMD_NFW_Scalar = ( (2.e0_double*rs*delta_c*rho_c)/(x*x - 1.e0_double) )*(1.e0_double - (2.e0_double/dsqrt(x*x-1.e0_double))*atan(dsqrt((x-1.e0_double)/(1.e0_double+x))) )
      end if

    end function SMD_NFW_Scalar

    function Differential_SMD_Array(r, z, r200, r_s)
      real(double),intent(in)::r(:), r200, z
      real(double), intent(in),optional::r_s

      real(double), dimension(size(r)):: Sigma, SigmaBar, Differential_SMD_Array

      if(present(r_s)) then
         Sigma = SMD_NFW_Array(r,z,r200,r_s)
         SigmaBar = Mean_SMD_NFW_Array(r, z, r200, r_s)
      else
         Sigma = SMD_NFW_Array(r,z,r200)
         SigmaBar = Mean_SMD_NFW_Array(r, z, r200)
      end if

      Differential_SMD_Array = SigmaBar-Sigma

    end function Differential_SMD_Array

    real(double) function Differential_SMD_Scalar(r, z, r200, r_s)
      !--Can be converted into a tangential shear by dividing by 1/Sigma_{Crit}
      real(double),intent(in)::r, r200, z
      real(double), intent(in),optional::r_s

      real(double):: Sigma, SigmaBar

      if(present(r_s)) then
         Sigma = SMD_NFW_Scalar(r,z,r200,r_s)
         SigmaBar = Mean_SMD_NFW_Scalar(r, z, r200, r_s)
      else
         Sigma = SMD_NFW_Scalar(r,z,r200)
         SigmaBar = Mean_SMD_NFW_Scalar(r, z, r200)
      end if

      Differential_SMD_Scalar = SigmaBar-Sigma

    end function Differential_SMD_Scalar

    real(double) function Shear_NFW(r,z,r200, Sigma_Critical, r_s)
      !--Eqn 14 Wright and Brainerd. Tested against Differential_SMD/Sigma_Critical and works
      real(double),intent(in)::r, r200, z, Sigma_Critical
      real(double), intent(in),optional::r_s

      real(double)::delta_c, c, rs, M200
      real(double)::omega_matter !-Evaluated at the redshift of the Halo-!
                                                                                                                                                                                                                  
      real(double)::rho_c, g
      real(double)::x

      if(present(r_s)) then
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c,  r_s)
      else
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c)
      end if

      rs = r200/c

      x = r/rs
      !--Following from Wright, Brainerd--!                                                                                                                                                                       
      if(x < 1) then
         g = ( (8.e0_double*atanh(dsqrt((1.e0_double-x)/(1.e0_double+x))))/((x*x)*dsqrt(1.e0_double-(x*x)))  ) + ((4.e0_double/(x*x))*dlog(0.5e0_double*x)) - (2.e0_double/((x*x)-1.e0_double)) + ( (4.e0_double*atanh(dsqrt((1.e0_double-x)/(1.e0_double+x))))/( ((x*x)-1.e0_double)*((1.e0_double-(x*x))**0.5e0_double)  )  )
         Shear_NFW = (rs*delta_c*rho_c*g)/Sigma_Critical
      elseif(x==1) then
         Shear_NFW = ((rs*delta_c*rho_c)/Sigma_Critical)*( (10.e0_double/3.e0_double) + 4.e0_double*dlog(0.5e0_double) )
      else
         g = ( (8.e0_double*atan(dsqrt((x-1.e0_double)/(1.e0_double+x))))/((x*x)*dsqrt((x*x)-1.e0_double))  ) + ((4.e0_double/(x*x))*dlog(0.5e0_double*x)) - (2.e0_double/((x*x)-1.e0_double)) + ((4.e0_double*atan(dsqrt((x-1.e0_double)/(1.e0_double+x))))/(((x*x)-1.e0_double)**1.5e0_double)   )
         Shear_NFW = (rs*delta_c*rho_c*g)/Sigma_Critical
      end if

    end function Shear_NFW

    function Mean_SMD_NFW_Array(r, z, r200, r_s)
      real(double),intent(in)::r(:), r200, z
      real(double), intent(in),optional::r_s
      
      real(double),dimension(size(r))::Mean_SMD_NFW_Array

      integer::i

      do i = 1, size(r)
         if(present(r_s)) then
            Mean_SMD_NFW_Array(i) = Mean_SMD_NFW_Scalar(r(i), z, r200, r_s)
         else
            Mean_SMD_NFW_Array(i) = Mean_SMD_NFW_Scalar(r(i), z, r200)
         end if
      end do

    end function Mean_SMD_NFW_Array

    real(double) function Mean_SMD_NFW_Scalar(r, z, r200, r_s)
      real(double),intent(in)::r, r200, z
      real(double), intent(in),optional::r_s
      
      real(double)::delta_c, c, rs, M200
      real(double)::omega_matter !-Evaluated at the redshift of the Halo-!                                                                                                                                         
      real(double)::rho_c
      real(double)::x

      if(present(r_s)) then
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c,  r_s)
      else
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c)
      end if

      rs = r200/c

      x = r/rs
      !--Following from Wright, Brainerd--!
      if(x < 1) then
         Mean_SMD_NFW_Scalar = (4.e0_double/(x*x))*(rs*delta_c*rho_c)*( (2.e0_double/dsqrt(1.e0_double-(x*x)))*atanh(dsqrt((1.e0_double-x)/(1.e0_double+x))) + dlog(0.5e0_double*x) )
!         Mean_SMD_NFW_Scalar = (4.e0_double/(x*x))*rs*delta_c*rho_c*( (2.e0_double/dsqrt(1.e0_double-x*x))*atanh(dsqrt((1.e0_double-x)/(1.e0_double+x))) +dlog(0.5e0_double*x) )
      elseif(x == 1) then
         Mean_SMD_NFW_Scalar = 4.e0_double*rs*delta_c*rho_c*(1.e0_double + dlog(0.5e0_double))
!         Mean_SMD_NFW_Scalar =  4.e0_double*rs*delta_c*rho_c*(1.e0_double+dlog(0.5e0_double))
      else
         Mean_SMD_NFW_Scalar = (4.e0_double/(x*x))*(rs*delta_c*rho_c)*( (2.e0_double/dsqrt((x*x)-1.e0_double))*atan(dsqrt((x-1.e0_double)/(1.e0_double+x))) + dlog(0.5e0_double*x) )
!         Mean_SMD_NFW_Scalar = (4.e0_double/(x*x))*rs*delta_c*rho_c*( (2.e0_double/dsqrt(x*x-1.e0_double))*atan(dsqrt((x-1.e0_double)/(1.e0_double+x))) + dlog(0.5e0_double*x) )
      end if

    end function Mean_SMD_NFW_Scalar

    function get_NFW_VirialMass_from_VirialRadius(z, VirialRadius)
      use Cosmology, only: Normalised_Hubble_Parameter
      real(double), intent(in)::VirialRadius(:), z
      real(double),dimension(size(VirialRadius))::get_NFW_VirialMass_from_VirialRadius

      real(double)::rho_c

      rho_c = (3.e0_double/(8.e0_double*3.142e0_double*4.3017e-13_double))*Normalised_Hubble_Parameter(z)

      get_NFW_VirialMass_from_VirialRadius = ( (800e0_double*3.142e0_double)/(3.e0_double) )*(VirialRadius**3.e0_double) * rho_c

    end function get_NFW_VirialMass_from_VirialRadius

    function get_NFW_VirialRadius_from_VirialMass_Scalar(z,VirialMass)
      use Cosmology, only: Normalised_Hubble_Parameter
      real(double), intent(in)::VirialMass, z
      real(double)::get_NFW_VirialRadius_from_VirialMass_Scalar

      real(double)::rho_c

      rho_c = (3.e0_double/(8.e0_double*3.142e0_double*4.3017e-13_double))*Normalised_Hubble_Parameter(z)

      get_NFW_VirialRadius_from_VirialMass_Scalar = ( (VirialMass*3.e0_double)/(800e0_double*3.142e0_double*rho_c) )**(1.e0_double/3.e0_double)

    end function get_NFW_VirialRadius_from_VirialMass_Scalar


    subroutine get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c, r_s)
      !--TESTED with values 18Mar2014--!
      use Cosmology, only: Normalised_Hubble_Parameter
      
      real(double),intent(in)::r200, z
      real(double), intent(in),optional::r_s
      real(double),intent(out)::delta_c, c, M200, rho_c

      real(double)::omega_matter,rs !-Evaluated at the redshift of the Halo-!                                                                                                                                        

      !--Get M200 for entered r200--!
      rho_c = (3.e0_double/(8.e0_double*3.142e0_double*4.3017e-13_double))*Normalised_Hubble_Parameter(z) !-in Units [M_Sun/h][h/Mpc]^3. Uses G/H_0 = 4.3017e-13 (Mpc/h)^3(h/M_sun)-!

      M200 = ( (800e0_double*3.142e0_double)/(3.e0_double) )*(r200**3.e0_double) * rho_c !-Could also use rho_bar = rho_crit * Omega_matter(z) (eg CH08)


      !--Get Concentration/rs
      if(present(r_s)) then
         rs = r_s      
         c = r200/rs
      else
         c = NFW_get_concentration_fit(z, M200)
         rs = r200/c
      end if

      Omega_matter = 1.e0_double!0.3e0_double*((1+z)**3.e0_double) -Incorrect, as rho_crit needs to vary also
      delta_c =  ((200.e0_double*Omega_Matter*c*c*c)/3.e0_double)/(dlog(1.e0_double+c) - c/(1.e0_double+c))
      
    end subroutine get_NFW_Parameters

    real(double) function NFW_get_Concentration_fit(z, M)
      !--Uses the fit of Dolag et al 2004 (see CH08) to get the concentration from the halo virial mass (M, usually M200) and the redshift of the halo--! 
      !--M must be in M_Sun/h--!
      real(double), intent(in)::z, M

      !--Fit Parameters--!
      real(double), parameter:: alpha = -0.102e0_double, c_0 = 9.59e0_double

      NFW_get_Concentration_fit = (c_0/(1.e0_double+z))*((M/10.e14_double)**alpha)

    end function NFW_get_Concentration_fit

    !--------------------------_End NFW Routines_-------------------------------------------------------------!

    !--------------------------------- Singular Isothermal Sphere -----------------------------------------------!
    real(double) function SIS_Mass_withinRadius(R, Sigv2)
      !--Free Parameter here is velocity disperion squared--! 
      real(double), intent(in)::R, Sigv2

      SIS_Mass_withinRadius = (6.2832e0_double*Sigv2*R)/(2.e0_double*4.3017e-9_double) !-(2*Pi*Sigv^2*R)/(2G)-!

    end function SIS_Mass_withinRadius

    real(double) function Error_SIS_Mass_withinRadius(R, eSigv2)
      !--Free Parameter here is velocity dispersion squared--!
      real(double), intent(in)::R, eSigv2

      Error_SIS_Mass_withinRadius = (6.2832e0_double*eSigv2*R)/(2.e0_double*4.3017e-9_double)

    end function Error_SIS_Mass_withinRadius

    function SMD_SIS_Array(velocity_dispersion2, radius, core, truncation)
      real(double), intent(in):: velocity_dispersion2, radius(:)
      real(double), intent(in),optional:: core, truncation
      
      real(double), dimension(size(radius)):: SMD_SIS_Array

      real(double)::Core_Radius, Truncation_Radius

      integer::i

      Core_Radius = 0.e0_double; Truncation_Radius = 1.e50_double
      if(present(Core)) Core_Radius = Core
      if(present(Truncation)) Truncation_Radius = Truncation

      do i = 1, size(radius)
         SMD_SIS_Array(i) = SMD_SIS_Scalar(velocity_dispersion2, radius(i), Core_Radius, Truncation_Radius)
      end do

    end function SMD_SIS_Array


    real(double) function SMD_SIS_Scalar(velocity_dispersion2, radius, core, truncation)
      real(double), intent(in):: velocity_dispersion2, radius
      real(double), intent(in),optional:: core, truncation

      real(double)::Core_Radius, Truncation_Radius

      Core_Radius = 0.e0_double; Truncation_Radius = 1.e50_double
      if(present(Core)) Core_Radius = Core
      if(present(Truncation)) Truncation_Radius = Truncation

      SMD_SIS_Scalar =  ( (Velocity_Dispersion2)/(2.e0_double*4.3017e-9_double) )*( 1.e0_double/(dsqrt(Radius*Radius + Core_Radius*Core_Radius)) - 1.e0_double/(dsqrt(Radius*Radius + Truncation_Radius*Truncation_Radius)) )


    end function SMD_SIS_Scalar

    real(double) function SIS_velocity_dispersion_byMass(Mass, Mass_Scale)
      !--Does not work if there is a core radius-!
      !--Mass_Scale needs to be in Mpc/h--!
      !--Returns Km/s--!
      real(double), intent(in):: Mass, Mass_Scale

      SIS_velocity_dispersion_byMass = dsqrt((4.3017e-9_double*Mass)/(3.142e0_double*Mass_Scale))

    end function SIS_velocity_dispersion_byMass

    !--------Flat Profile---------------------------!

    real(double) function Flat_Mass_withinRadius(R, Sigma_0)
      !--Scale must be in Mpc/h --!
      !-Assumes Spherically Symmetric-!
      real(double), intent(in)::R, Sigma_0

      Flat_Mass_withinRadius = 3.142e0_double*R*R*Sigma_0

    end function Flat_Mass_withinRadius


    real(double) function Error_Flat_Mass_withinRadius(R, eSigma_0)
      !--Scale must be in Mpc/h --!
      !-Assumes Spherically Symmetric-!
      real(double), intent(in)::R, eSigma_0

      Error_Flat_Mass_withinRadius = 3.142e0_double*R*R*eSigma_0

    end function Error_Flat_Mass_withinRadius


end module Mass_Profiles
