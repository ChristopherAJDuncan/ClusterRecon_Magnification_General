module cosmology
  use Param_Types; 
  implicit none
  
  character(50),private::Rad_Dir
  logical::Verbose


  type Cosmology_Parameters
     real(double)::Om_CDM, Om_B, Om_tot, Om_l, w, h_100
  end type Cosmology_Parameters

  !-Globals Used for Integration Routine-!
  type(Cosmology_Parameters),private::Global_Cosmology_Params


  INTERFACE angular_diameter_distance_fromRedshift
     MODULE PROCEDURE angular_diameter_distance_fromRedshift_array, angular_diameter_distance_fromRedshift_scalar
  END INTERFACE angular_diameter_distance_fromRedshift

  INTERFACE luminosity_distance
     MODULE PROCEDURE luminosity_distance_array, luminosity_distance_scalar
  END INTERFACE luminosity_distance
  
contains

  function LCDM
    type(Cosmology_Parameters)::LCDM
    
    LCDM%Om_CDM = 0.2544e0_double
    LCDM%Om_B = 0.0456e0_double
    LCDM%Om_l = 0.7e0_double
    LCDM%Om_tot = LCDM%Om_CDM+LCDM%Om_B+LCDM%Om_l
    LCDM%w = -1.e0_double
    LCDM%h_100 = 0.7e0_double

  end function LCDM

  function luminosity_distance_array(z1,z2,Cos)
        !-Returned in units of Mpc/h                                                                                                                                                                                                         
    type(Cosmology_Parameters),optional::cos
    real(double),intent(in)::z1
    real(double),intent(in)::z2(:)

    real(double),dimension(size(z2))::luminosity_distance_array

    real(double):: angular_diameter_distance(size(z2))

    type(Cosmology_Parameters)::Internal_Cos

    if(present(cos)) then
       Internal_Cos = Cos
    else
       Internal_Cos = LCDM()
    end if

    angular_diameter_distance =angular_diameter_distance_fromRedshift_array(z1,z2, Internal_Cos)
    luminosity_distance_array = ( ((1+z2)/(1+z1))**2.e0_double )*angular_diameter_distance

  end function luminosity_distance_array

  function luminosity_distance_scalar(z1,z2,Cos)
    !-Returned in units of Mpc/h                                                                                                                                                                                                         
    type(Cosmology_Parameters),optional::cos
    real(double),intent(in)::z1
    real(double),intent(in)::z2

    real(double)::luminosity_distance_scalar

    real(double)::lum_distance(1)

    lum_distance = luminosity_distance_array(z1,(/z2/),Cos)
    luminosity_distance_scalar = lum_distance(1)

  end function luminosity_distance_scalar

  function angular_diameter_distance_fromRedshift_array(z1,z2, Cos)
       !-Returned in units of Mpc/h                                                                                     
    type(Cosmology_Parameters),optional::cos
    real(double),intent(in)::z1
    real(double),intent(in)::z2(:)

    real(double),dimension(size(z2))::angular_diameter_distance_fromRedshift_array

    type(Cosmology_Parameters)::Internal_Cos
    real(double),allocatable::Radius(:)

    if(present(cos)) then
       Internal_Cos = Cos
    else
       if(Verbose) print *, 'Using LCDM to get Distances'
       Internal_Cos = LCDM()
    end if

    call getrarray(z1, z2, Radius, Internal_Cos)
    
    !--This is flat case only--!
    !--Extra h_100 converts from Mpc to Mpc/h--!
    angular_diameter_distance_fromRedshift_array = (Radius*Internal_Cos%h_100)/(1.e0_double+z2)

    deallocate(Radius)


  end function angular_diameter_distance_fromRedshift_array

 function angular_diameter_distance_fromRedshift_scalar(z1, z2, cos)
   !-Returned in units of Mpc/h
   !--This is returned as a proper distance--!
    type(Cosmology_Parameters),optional::cos
    real(double),intent(in)::z1
    real(double),intent(in)::z2

    real(double)::angular_diameter_distance_fromRedshift_scalar

    type(Cosmology_Parameters)::Internal_Cos
    real(double),allocatable::Radius(:)

    if(present(cos)) then
       Internal_Cos = Cos
    else
       if(Verbose) print *, 'Using LCDM to get Distances'
       Internal_Cos = LCDM()
    end if

    call getrarray(z1, (/z2/), Radius, Internal_Cos)
    
    !--Extra h_100 converts from Mpc to Mpc/h--!
    angular_diameter_distance_fromRedshift_scalar = (Radius(1)*Internal_Cos%h_100)/(1.e0_double+z2)

    deallocate(Radius)

  end function angular_diameter_distance_fromRedshift_scalar


  function zradiusintegrand(z)
    use nrtype;
    real(SP),dimension(:),intent(in)::z
    real(SP),dimension(size(z))::zradiusintegrand
    !needs generalised to account for how neutrinos vary with redshift
    

    Global_Cosmology_Params%Om_tot = Global_Cosmology_Params%Om_l + Global_Cosmology_Params%Om_B + Global_Cosmology_Params%Om_CDM !+ Global_Cosmology_Parametersn

    if(Global_Cosmology_Params%Om_tot .eq. 0) then
       print *, 'zradiusintegrand - Global_Cosmology_Params%Om_tot equal to zero'
    end if
    if (Global_Cosmology_Params%h_100 .eq. 0) then
       print *, 'zradiusintegrand - h not set'
       STOP
    end if

    if(any(z .lt. 0)) then
       print *, 'zradiusintegrand - redshift arg less than zero'
       STOP
    END if
    
    !does this need generalised to account for omega_neutrinos
    zradiusintegrand = 3000.e0_SP/(Global_Cosmology_Params%h_100*dsqrt(Omega_lambda(z) +(Global_Cosmology_Params%Om_CDM+Global_Cosmology_Params%Om_B)*((1+z)**3.e0_SP) + (1-Global_Cosmology_Params%Om_tot)*((1+z)**2.e0_SP)))

    if(any(zradiusintegrand .lt. 0)) then
       print *, 'zradiusintegrand - less than zero'
       STOP
    end if
    

  end function zradiusintegrand 


  function Omega_lambda(z)
    use nrtype;
    real(SP),dimension(:),intent(in)::z
    real(SP),dimension(size(z))::Omega_lambda
    real(SP), dimension(size(z))::w
    real(SP)::int
    integer::i
    integer::count = 0
    real(double)::wa = 0.e0_double

    count = count + 1

    Omega_lambda = Global_Cosmology_Params%Om_l*((1.e0_SP+z)**(3.e0_SP*(1+Global_Cosmology_Params%w)))
    if(wa .ne. 0.e0_double) then !this section needs tested
       STOP 'USING WA /= 0!!!'
       do i = 1, size(z)
          Omega_lambda(i) = Omega_lambda(i)*((1.e0_double+z(i))**(3.e0_double*wa))*dexp((z(i)*wa)/(1.e0_double+z(i)))
       end do
    end if

  end function Omega_lambda

  function Omega_lambda_int(z)
    !This isn't used yet, what is the integral wrt? Probably a or redshift - would not be used unless wa \= 0?
    use nrtype;
    real(SP),dimension(:),intent(in)::z
    real(SP),dimension(size(z))::Omega_lambda_int
    real(SP), dimension(size(z))::w
    real(SP)::wa = 0.e0_SP !in generalisation will be made a global parameter

    w = Global_Cosmology_Params%w + (z/(1.e0_SP+z))*wa
    Omega_lambda_int = (-1.e0_SP)*(1.e0_SP + w)/(1 + z)

  end function Omega_lambda_int

  subroutine getrarray(z1,z2,r,cos)
    !calculates radius values for z values in zarray, and places them in appropriate slot of rarray
    !- Returned in units of Mpc-!
    !-Calculates CO-MOVING RADIUS [dChi or R_0dr]--!
    use nr;
    type(cosmology_parameters),optional::cos
    real(double), dimension(:),allocatable::r
    real(double),intent(in)::z1
    real(double)::z2(:)
    integer::i 
    real(double)::ztemp
    character(len=200)::filename

    if(allocated(r)) then
       deallocate(r)
    end if
    allocate(r(size(z2))); r = -1.e0_double

    !set cosmology params that integrand will use. If coss present, these will be set, otherwise we assume global cosmology is set appropriately
    if(present(cos)) then
       Global_Cosmology_Params = Cos
    else
       Global_Cosmology_Params = LCDM()
    end if
    
    !get radius for all z values in z, construct r
    do i=1, size(z2)
       if(z2(i) >= z1) then
          r(i) = qromo(zradiusintegrand,z1, z2(i), midpnt, 1.e-6_double)
       else
          r(i) = 0.e0_double
       end if
    end do
    
    if(any(r<0.e0_double)) STOP 'Error determining distances from redshift - negatives still exist'

  end subroutine getrarray


end module cosmology
