module Distributions
  use Param_Types
implicit none

character(200),private::Reference_Catalogue = 'Catalogues/STAGES_COMBO17_gsMag_size_matched.pzcat'
integer::Reference_Catalogue_Columns(13) = (/-1,10,11,-1,-1,8,-1,-1,-1,14,16,17,5/)

contains


  subroutine randomly_Sample_From_Distribution(Grid, PDF, Res)
    !--Populates the array Res by Randomly Sampling from the PDF described by Grid and PDF--!
    !--Assumes that the grid is equally binned, and the PDF is renormalised--!
    real(double), intent(in)::Grid(:), PDF(:)
    real(double), intent(out):: Res(:)

    real(double)::dGrid

    !--Random Number Generator Declarations--!
    real(double),dimension(:),allocatable::Ran
    Integer,allocatable::seed(:)
    integer::Nseed, Clock
    integer::NRandom 

    integer::i, F

    !--Intrinsic PDF declarations--!
    real(double), allocatable::CumulativePDF(:)
    real(double), allocatable::Renormalised_PDF(:)
    real(double)::AreaUnderCurve

    if(size(Grid) /= size(PDF)) STOP 'randomly_Sample_From_Distribution - Grid and PDF entered not conformal'

    NRandom = size(Res) !-1 random number for each galaxy-!
    allocate(Ran(NRandom)); Ran = 0.e0_double
    call RANDOM_SEED(size = NSeed)
    allocate(Seed(NSeed))
    call SYSTEM_CLOCK(COUNT = Clock)
    seed = Clock + (/ (i-1,i=1,NSeed) /)
    call RANDOM_SEED(PUT = seed)
    deallocate(Seed); NSeed = 0; Clock = 0
    
    call RANDOM_NUMBER(Ran)

    if(size(Grid) < 4) STOP 'randomly_Sample_From_Distribution - Grid is too small, something is possibly up'
    !-_Assumes that the grid is equally binned-!
    dGrid =  Grid(2)-Grid(1)!0.5e0_double*(Grid(4)-Grid(2))

    !--Renormalise PDF such that CumulativePDF always reaches 1--!
    AreaUnderCurve = dGrid*PDF(1)
    do i =2, size(PDF)
       dGrid = (Grid(i)-Grid(i-1)) !-0.5e0_double*(Grid(i+1)-Grid(i-1))-!
       AreaUnderCurve = AreaUnderCurve + dGrid*PDF(i)
    end do
    allocate(Renormalised_PDF(size(PDF))); Renormalised_PDF = 0.e0_double
    Renormalised_PDF = PDF/AreaUnderCurve

    !-Construct Cumulative PDF-!
    dGrid =Grid(2)-Grid(1)
    allocate(CumulativePDF(size(PDF))); CumulativePDF = 0.e0_double
    CumulativePDF(1) = Renormalised_PDF(1)*dGrid
    do i = 2, size(Renormalised_PDF)
       dGrid = (Grid(i)-Grid(i-1))
       CumulativePDF(i) = CumulativePDF(i-1)+(Renormalised_PDF(i)*dGrid)
       if(CumulativePDF(i)-1.e0_double > 1.e-6_double) then
          print *, CumulativePDF(i)-1.e0_double
          STOP 'randomly_Sample_From_Distribution - Cumulative PDF is getting too large, larger than rounding tolerance of 1'
       end if
    end do
    !--If by integration error it is within a rounding tolerance of 1, then the cumualtive PDF should be set to one at it's max value

    if(dabs(CumulativePDF(size(CUmulativePDF))-1.e0_double) <= 1.e-6_double) then
       CumulativePDF(size(CUmulativePDF)) = 1.e0_double
!!$       if(CumulativePDF(size(CUmulativePDF)-1)-CumulativePDF(size(CUmulativePDF)) > 0.e0_double) then
!!$          print *, CumulativePDF(size(CUmulativePDF)-1), CumulativePDF(size(CUmulativePDF)), CumulativePDF(size(CUmulativePDF)-1)-CumulativePDF(size(CUmulativePDF))
!!$          STOP 'randomly_Sample_From_Distribution - Error in adjusting for rounding error - Cumulative PDF not increasing near upper limit...'
!!$       end if
    end if

    if(any(Ran < 0.e0_double) .or. any(Ran > 1.e0_double)) STOP 'randomly_Sample_From_Distribution - FATAL ERROR - Random array values lies outside bounds (0,1)'
    if(any(CumulativePDF < 0.e0_double) .or. any(CumulativePDF > 1.e0_double+1.e-3_double)) then
       print *, 'Cumulative PDF:', CumulativePDF
       print *, any(CumulativePDF < 0.e0_double), any(CumulativePDF > 1.e0_double)
       STOP 'randomly_Sample_From_Distribution - FATAL ERROR - Cumulative PDF array values lies outside bounds (0,1)'
    end if

    if(maxval(CumulativePDF) < 1.e0_double) then
       print *, maxval(CumulativePDF)
       STOP 'randomly_Sample_From_Distribution - Cumulatie PDF has not reached appropriate precision, and has not reached one, stopping'
    end if

    !--Randomly Sample from PDF-!
    !--Assumes that the PDF is NORMALISED--!
    do i= 1, size(Ran)
       do F = 1, size(CumulativePDF)-1
          !--First two cases account for fintie size of Cumulative PDF and Rounding Errors--!
          if(Ran(i) < CumulativePDF(1)) then
             Res(i) = Grid(1)
          elseif(Ran(i) > maxval(CumulativePDF)) then
             print *, Ran(i), maxval(CumulativePDF)
             STOP 'rANDOM NUMBER IS LARGER THAN CUMULATIVE PDF'
             Res(i) = Grid(size(Grid))
          elseif(CumulativePDF(F) == Ran(i)) then
             Res(i) = Grid(F)
             exit
          elseif( (CumulativePDF(F) < Ran(i)) .and. (CumulativePDF(F+1) > Ran(i)) ) then
             !-Linearly Interpolate-!
             Res(i) = Grid(F) + ((Grid(F+1)-Grid(F))/(CumulativePDF(F+1)-CumulativePDF(F)))*(Ran(i)-CumulativePDF(F))
             exit
          elseif(CumulativePDF(F+1) == Ran(i)) then
             Res(i) = Grid(F+1)
             exit
          elseif(F == size(CumulativePDF)-1) then
             print *, 'Loop:', i
             print *, 'Max Cumulative, Random:', CumulativePDF(size(CUmulativePDF)), Ran(i)
             print *, 'Min Cumulative, Random:', CumulativePDF(1), Ran(i)
             STOP 'Assign_Intrinsic_Sizes - FATAL ERROR - Cannot assign magnitude'
          end if
       end do
    end do 

  end subroutine randomly_Sample_From_Distribution


  !---------------REDSHIFT DISTRIBUTIONS--------------------------------------------------------!

    subroutine Analytic_Source_Redshift_Distribution(alpha, beta, z_med, Redshift, pdf, Output_Directory)
    real(double), intent(in)::Redshift(:)
    real(double), intent(out),allocatable::pdf(:)
    real(double),intent(in):: alpha, beta, z_med
    character(*),optional::Output_Directory

    integer::Method = 1 !-1:Smail-!                                                                                                                                                                                                           
    integer::i

    real(double)::AreaUnderCurve, dRedshift

    !--Smail-!      
    real(double):: z_0

    allocate(PDF(size(Redshift))); PDF = -1.e0_double

    if (z_med <= 0.e0_double) STOP 'Analytic_Source_Redshift_Distribution - MEDIAN REDSHIFT NEGATIVE OR ZERO, FATAL'

    AreaUnderCurve = 0.e0_double; dRedshift = Redshift(2)-Redshift(1)
    select case(Method)
    case(1) !-Smail-!
       !print *, 'Producing a Smail et al source redshift distribution with alpha = ', alpha, '; beta = ', beta, '; z_med  = ', z_med 
       z_0 = z_med/1.412e0_double
       do i = 1, size(Redshift)
          PDF(i) = (Redshift(i)/z_0)**2.e0_double*dexp(-(Redshift(i)/z_0)**1.5e0_double)
          if(i>1) AreaUnderCurve = AreaUnderCurve + 0.5e0_double*(Redshift(i)-Redshift(i-1))*(PDF(i)+PDF(i-1))
       end do
       !-Renormalise-!
       PDF = PDF/AreaUnderCurve
    end select

    if(any(PDF < 0.e0_double)) STOP 'Source_Redshift_Distribution - FALTAL ERROR - PDF contains negatives'

   !-Output-!                                                                                                                                                                                                                                
    if(present(Output_Directory)) then
       open(32, file = trim(adjustl(Output_Directory))//'Source_Redshift_Distribution_Smail.dat')
       !print *, 'Source Redshift Distribution output to:' ,trim(adjustl(Output_Directory))//'Source_Redshift_Distribution_Smail.dat'                                                                                                         
    else
       open(32, file = 'Distributions/Source_Redshift_Distribution_Smail.dat')
       !print *, 'Source Redshift Distribution output to: Distributions/Source_Redshift_Distribution_Smail.dat'   
    end if
    do i = 1, size(PDF)
       write(32, *) Redshift(i), PDF(i)
    end do
    close(32)

  end subroutine Analytic_Source_Redshift_Distribution



  subroutine get_Redshift_Distribution_MagnitudeBinning_byCatalogue(MagBins, Redshifts, PDF, RefCat, Output_Dir)
    use Catalogues; use Statistics, only: variance_discrete, median
    real(double), intent(out),allocatable:: Redshifts(:)
    real(double), intent(out),allocatable::pdf(:,:) !-Redshift Bin, GridValue-!
    real(double),dimension(:,:),allocatable,intent(inout)::MagBins
    character(*), intent(in), optional:: Output_Dir
    type(Catalogue), intent(in)::RefCat

    character(200)::Catalogue_Filename
    integer::Catalogue_Columns(13)

    type(Catalogue)::Cat
    
    integer::nMags 
    real(double)::Renormalisation

    integer,parameter::nZ = 20
    real(double)::Z_lower, Z_Higher, dZ
    real(double),dimension(:,:),allocatable::ZBins
    
    integer::i,j,c
    
    character(5)::fmtstring
    
    type(Binned_Catalogue)::BCat
    real(double),allocatable::Temporary_Z_Array(:)

    !--Fit declarations [ to z_m =  a*m + b]--!
    real(double),allocatable::median_redshift(:)
    real(double),allocatable::mean_Mags(:)
    real(double)::fita, fitb, siga, sigb, chi2, fitq

    Catalogue_Filename = trim(Reference_Catalogue)
    Catalogue_Columns = Reference_Catalogue_Columns
    
    PRINT *, 'Reconstructing p(z,m) from Reference Catalogue:'
    Cat = RefCat
!!$    
!!$    PRINT *, 'Reconstructing p(z,m) distribution from the Catalogue:', trim(Catalogue_Filename),'....'
!!$    !--Read in Reference Catalogue--!
!!$    call Catalogue_ReadIn(Cat, trim(Catalogue_Filename), 'COMBO', Catalogue_Columns)

    !--Test for redshifts for all galaxies in reference catalogue--!

    Z_Higher = maxval(Cat%redshift); Z_Lower = 0.e0_double
    allocate(Redshifts(nZ)); Redshifts = 0.e0_double
    allocate(ZBins(nZ,2)); ZBins = 0.e0_double
    dZ = (( Z_Higher- Z_Lower )/(1.e0_double*(nZ-1)) )
    do i = 1, nZ
       !--Use i-2 so that no galaxies fall into the first bin--!                                                                                                                                                                         
       ZBins(i,1) = Z_Lower + (i-2)*dZ
       ZBins(i,2) = ZBins(i,1) + dZ
       if(i==nZ) ZBins(i,2) = ZBins(i,2) + 1.e-3_double*dZ
       Redshifts(i) = 0.5e0_double*(ZBins(i,1)+ZBins(i,2))
    end do
    Redshifts(1) = ZBins(1,2)

    !----START HERE----!
    if(allocated(MagBins) == .false.) then
       nMags = 5
       print *, 'Magnitude Limits not passed so getting Mag Limitsf ro nBin', nMags
       call Calculate_Bin_Limits_by_equalNumber(Cat%Absolute_Magnitude, nMags, MagBins)
    end if
    nMags = size(MagBins,1)
    !--Bin Catalogue by Magnitude--!
    call bin_catalogue_by_magnitude(Cat,MagBins,BCat)

    do i =1 , nMags
       print *, 'Mag Bin ', i, ' : ', MagBins(i,:)
    end do
    
    allocate(PDF(nMags, nZ)); PDF = 0.e0_double
    do i = 1, nMags
       allocate(Temporary_Z_Array(size(BCat%Cat(i)%Redshift))); Temporary_Z_Array = BCat%Cat(i)%Redshift
       do c = 1, size(Temporary_Z_Array)
          do j = 1, nZ
             if( (Temporary_Z_Array(c) > ZBins(j,1)) .and. (Temporary_Z_Array(c) <= ZBins(j,2)) ) then
                PDF(i,j) = PDF(i,j) + 1
                exit
             end if
          end do
       end do
       deallocate(Temporary_Z_Array)
    end do

    !--Renormalise for each Magnitude Bin--!                                                                                                                                                                                                 
    do i =1, nMags
       Renormalisation = 0.e0_double
       do j = 1, size(PDF,2)
          Renormalisation = Renormalisation + PDF(i,j)*(ZBins(i,2)-ZBins(i,1))
       end do
       if(Renormalisation> 0.e0_double) PDF(i,:) = PDF(i,:)/Renormalisation
    end do

    !--Output PDFS--!                                                                                                                                                                                                                   
    if(present(output_dir)) then
       open(unit = 49, file = trim(Output_Dir)//'Redshift_Distribution_MagnitudeBinning_Catalogue.dat')
           print *, 'Redshift Distribution, by Magnitude Bin, output to '//trim(Output_Dir)//'Redshift_Distribution_MagnitudeBinning_Catalogue.dat'
    else
       open(unit = 49, file = 'Distributions/Redshift_Distribution_MagnitudeBinning_Catalogue.dat')
       print *, 'Size Distribution, by Magnitude Bin, output to Distributions/Redshift_Distribution_MagnitudeBinning_Catalogue.dat'
    end if
    !--Write Header--!
    do j = 1, size(MagBins,1)
       write(49, '(A1, 2(e14.7,x))') '#', MagBins(j,:)
    end do
    write(49,'(A)')
    write(fmtstring, '(I5)') size(PDF,1)+1
    do j = 1, size(PDF,2)
       write(49, '('//trim(fmtstring)//'(e14.7,x))') Redshifts(j), PDF(:,j)
    end do
    close(49)

!!$    print *, 'Getting fit of median redshift to *absolute magnitude*:'
!!$    allocate(median_redshift(nMags)); allocate(mean_Mags(nMags));
!!$    do i = 1, nMags
!!$       mean_Mags(i) = 0.5e0_double*(sum(MagBins(i,:)))
!!$       median_redshift(i) = median(pdf(i,:), redshifts)
!!$    end do
!!$    call fit(mean_Mags, median_redshift, fita, fitb, siga, sigb, chi2, fitq) !--Doesn't use error for median redshift--!
!!$    print *, 'Redshift distribution fit to: z_m = ',fita, '*m +', fitb

    deallocate(median_redshift, mean_Mags)

  end subroutine get_Redshift_Distribution_MagnitudeBinning_byCatalogue
  
  !----------------SIZE DISTRIBUTIONS-----------------------------------------------------------!
  subroutine get_Size_Distribution_MagnitudeBinning_byCatalogue(MagBins, Sizes, PDF, RefCat, use_Physical_sizes,  Magnitude_Type, Output_Dir)
    !--Edit to also return from file input
    use Catalogues; use Statistics, only: variance_discrete
       real(double), intent(out),allocatable:: Sizes(:)
       real(double), intent(out),allocatable::pdf(:,:) !-Mag Bin, GridValue-!
       real(double),dimension(:,:),allocatable,intent(inout)::MagBins
       character(*), intent(in), optional:: Output_Dir
       type(Catalogue),intent(in)::RefCat
       logical,intent(in)::use_Physical_Sizes
       integer, intent(in)::Magnitude_Type

       character(200)::Catalogue_Filename
       integer::Catalogue_Columns(13)

       type(Catalogue)::Cat
       
       integer::nMags = 5
       real(double)::Renormalisation
       
       !--Cuts--!  May not be implemented as want the result to resemble the catalogue as much as possible                                                                                                                                  
       integer,parameter::nSizes = 20
       logical::Apply_Size_Cuts = .false.
       real(double):: Size_Cuts(2) = (/0.,20./)
       real(double)::Size_lower, Size_Higher
       !    integer:: Catalogue_Size_Column = 14 !!!!!                                                                                                                                                                                               
       real(double)::dSize
       real(double),dimension(:,:),allocatable::SizeBins
       
       integer::i,j,c
       
       character(5)::fmtstring
       
       type(Binned_Catalogue)::BCat
       real(double),allocatable::Temporary_Sizes_Array(:)


       INTERFACE
          subroutine get_Size_Distribution_MagnitudeBinning_byCatalogue(MagBins, Sizes, PDF, RefCat, use_Physical_sizes, Magnitude_Type, Output_Dir)
            use Param_Types; use Catalogues
            real(double), intent(out),allocatable:: Sizes(:)
            real(double), intent(out),allocatable::pdf(:,:) !-Mag Bin, GridValue-!      
            real(double),dimension(:,:),allocatable,intent(inout)::MagBins
            type(Catalogue),intent(in)::RefCat
            logical,intent(in)::use_Physical_Sizes
            integer, intent(in)::Magnitude_Type

            character(*), intent(in), optional:: Output_Dir
          end subroutine get_Size_Distribution_MagnitudeBinning_byCatalogue
       END INTERFACE

!!$       Catalogue_Filename = trim(Reference_Catalogue)
!!$       Catalogue_Columns = Reference_Catalogue_Columns
!!$       
!!$       PRINT *, 'Reconstructing magnitude distribution from the Catalogue:', trim(Catalogue_Filename), ' to get size distributions'
!!$       !--Read in Catalogue--!                                                                                                                                                                                                                   
!!$       call Catalogue_ReadIn(Cat, trim(Catalogue_Filename), 'COMBO', Catalogue_Columns)
       print *, 'Reconstructing size-magnitude distribution from the reference catalogue'
       Cat = RefCat


       !-Use the MC redshift sampling to assign physical sizes-!
!       call Monte_Carlo_Redshift_Sampling(RefCat)
!!$       call convert_Size_from_Pixel_to_Physical(Cat)
       
       if(use_Physical_Sizes) then
          Size_Higher = maxval(Cat%Physical_Sizes); Size_Lower = 0.e0_double!minval(Cat%Physical_Sizes)
       else
          Size_Higher = maxval(Cat%Sizes); Size_Lower = 0.e0_double!minval(Cat%Sizes);
       end if
       if(Apply_Size_Cuts) then
          Size_Lower = maxval((/Size_Cuts(1), Size_Lower/)); Size_Higher = minval( (/Size_Cuts(2), Size_Higher /) )
          print *, 'Applying Size Cuts of:', Size_Lower, Size_Higher
       end if

       print *, 'Size Dist: Limits:', Size_Lower, Size_Higher

       allocate(Sizes(nSizes)); Sizes = 0.e0_double
       allocate(SizeBins(nSizes,2)); SizeBins = 0.e0_double
       dSize = (( Size_Higher- Size_Lower )/(1.e0_double*(nSizes-1)) )
       if(dSize == 0.e0_double) STOP 'get_Size_Distribution_MagnitudeBinning_byCatalogue - Error setting up Size Grid, dSize = 0'
       do i = 1, nSizes
          !--Use i-2 so that no galaxies fall into the first bin--!                                                                                                                                                                         
          SizeBins(i,1) = Size_Lower + (i-2)*dSize
          SizeBins(i,2) = SizeBins(i,1) + dSize
          if(i==nSizes) SizeBins(i,2) = SizeBins(i,2) + 1.e-3_double*dSize
          Sizes(i) = 0.5e0_double*(SizeBins(i,1)+SizeBins(i,2))
       end do
       Sizes(1) = SizeBins(1,2)       

       !--Bin by equal number density in magnitude bins--!                                                                                                                                                                
       if(allocated(MagBins) == .false.) then
          print *, 'Size Distribution, Getting Mag Limits', nMags
          if(all(Cat%Absolute_Magnitude >= 0.e0_double)) STOP 'size_Disribution by magnitude:, Absolute magnitude not set'
          if(Magnitude_Type == 1) then !-Absolute_Magnitude-!
             call Calculate_Bin_Limits_by_equalNumber(Cat%Absolute_Magnitude, nMags, MagBins)
          elseif(Magnitude_Type == 2) then
             call Calculate_Bin_Limits_by_equalNumber(Cat%MF606W, nMags, MagBins)
          else
             STOP 'get_Size_Distribution_MagnitudeBinning_byCatalogue - Invalid Magnitude Type entered'
          end if
       end if
       !--Bin Catalogue by Magnitude--!
       nMags = size(MagBins,1)
       call bin_catalogue_by_magnitude(Cat,MagBins,BCat, Magnitude_Type)

       print *, 'Contructing from total of:', size(Cat%RA), ' galaxies, binned according to:', BCat%Occupation

       if(any(BCat%Occupation == 0)) STOP 'get_Size_Distribution_MagnitudeBinning_byCatalogue - Binning FATAL ERROR - Some bins contain no galaxies, check bin limits and type, stopping..'

       do i =1 , nMags
          print *, 'Mag Bin ', i, ' : ', MagBins(i,:), ' Variance in Physical Sizes:', dsqrt(variance_discrete(BCat%Cat(i)%Physical_Sizes, BCat%Cat(i)%Physical_Sizes))
       end do

           allocate(PDF(nMags, nSizes)); PDF = 0.e0_double
    do i = 1, nMags
       allocate(Temporary_Sizes_Array(size(BCat%Cat(i)%Sizes)));
       if(use_Physical_Sizes) then
          Temporary_Sizes_Array = BCat%Cat(i)%Physical_Sizes
       else
          Temporary_Sizes_Array= BCat%Cat(i)%Sizes
       end if
       
       
       do c = 1, size(Temporary_Sizes_Array)
          do j = 1, nSizes
             if( (Temporary_Sizes_Array(c) > SizeBins(j,1)) .and. (Temporary_Sizes_Array(c) <= SizeBins(j,2)) ) then
                PDF(i,j) = PDF(i,j) + 1
                exit
             end if
          end do
       end do
       deallocate(Temporary_Sizes_Array)
    end do

    !--Renormalise for each Magnitude Bin--!                                                                                                                                                                                                 
    do i =1, nMags
       Renormalisation = 0.e0_double
       do j = 1, size(PDF,2)
          Renormalisation = Renormalisation + PDF(i,j)*(SizeBins(i,2)-SizeBins(i,1))
       end do
       if(Renormalisation> 0.e0_double) PDF(i,:) = PDF(i,:)/Renormalisation
    end do

    !--Output Size PDFS--!                                                                                                                                                                                                                   
    if(present(output_dir)) then
       open(unit = 49, file = trim(Output_Dir)//'Size_Distribution_MagnitudeBinning_Catalogue.dat')
           print *, 'Size Distribution, by Magnitude Bin, output to '//trim(Output_Dir)//'Size_Distribution_MagnitudeBinning_Catalogue.dat'
    else
       open(unit = 49, file = 'Distributions/Size_Distribution_MagnitudeBinning_Catalogue.dat')
       print *, 'Size Distribution, by Magnitude Bin, output to Distributions/Size_Distribution_MagnitudeBinning_Catalogue.dat'
    end if
    !--Write Header--!
    do j = 1, size(MagBins,1)
       write(49, '(A1, 2(e14.7,x))') '#', MagBins(j,:)
    end do
    write(49,'(A)')
    write(fmtstring, '(I5)') size(PDF,1)+1
    do j = 1, size(PDF,2)
       write(49, '('//trim(fmtstring)//'(e14.7,x))') Sizes(j), PDF(:,j)
    end do
    close(49)

  end subroutine get_Size_Distribution_MagnitudeBinning_byCatalogue


end module Distributions
