module Distributions
  use Param_Types
implicit none

character(200),private::Reference_Catalogue = 'Catalogues/STAGES_COMBO17_gsMag_size_matched.pzcat'
integer::Reference_Catalogue_Columns(13) = (/-1,10,11,-1,-1,8,-1,-1,-1,14,16,17,5/)

contains
  !---------------REDSHIFT DISTRIBUTIONS--------------------------------------------------------!

  subroutine get_Redshift_Distribution_MagnitudeBinning_byCatalogue(MagBins, Redshifts, PDF, Output_Dir)
    use Catalogues; use Statistics, only: variance_discrete, median
    real(double), intent(out),allocatable:: Redshifts(:)
    real(double), intent(out),allocatable::pdf(:,:) !-Redshift Bin, GridValue-!
    real(double),dimension(:,:),allocatable,intent(inout)::MagBins
    character(*), intent(in), optional:: Output_Dir

    character(200)::Catalogue_Filename
    integer::Catalogue_Columns(13)

    type(Catalogue)::Cat
    
    integer::nMags = 5
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
    
    PRINT *, 'Reconstructing p(z,m) distribution from the Catalogue:', trim(Catalogue_Filename),'....'
    !--Read in Reference Catalogue--!
    call Catalogue_ReadIn(Cat, trim(Catalogue_Filename), 'COMBO', Catalogue_Columns)

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
       print *, 'Getting Mag Limits', nMags
       call Calculate_Bin_Limits_by_equalNumber(Cat%Mag, nMags, MagBins)
    end if
    nMags = size(MagBins,1)
    !--Bin Catalogue by Magnitude--!
    print *, 'Binning By Magnitude:'
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

    print *, 'Getting fit of median redshift to *absolute magnitude*:'
    allocate(median_redshift(nMags)); allocate(mean_Mags(nMags));
    do i = 1, nMags
       mean_Mags(i) = 0.5e0_double*(sum(MagBins(i,:)))
       median_redshift(i) = median(pdf(i,:), redshifts)
    end do
    call fit(mean_Mags, median_redshift, fita, fitb, siga, sigb, chi2, fitq) !--Doesn't use error for median redshift--!
    print *, 'Redshift distribution fit to: z_m = ',fita, '*m +', fitb

    deallocate(median_redshift, mean_Mags)

  end subroutine get_Redshift_Distribution_MagnitudeBinning_byCatalogue
  
  !----------------SIZE DISTRIBUTIONS-----------------------------------------------------------!
  subroutine get_Size_Distribution_MagnitudeBinning_byCatalogue(MagBins, Sizes, PDF, Output_Dir)
    !--Edit to also return from file input
    use Catalogues; use Statistics, only: variance_discrete
       real(double), intent(out),allocatable:: Sizes(:)
       real(double), intent(out),allocatable::pdf(:,:) !-Mag Bin, GridValue-!
       real(double),dimension(:,:),allocatable,intent(inout)::MagBins
       character(*), intent(in), optional:: Output_Dir

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
       
       logical::use_Physical_Sizes = .true.


       Catalogue_Filename = trim(Reference_Catalogue)
       Catalogue_Columns = Reference_Catalogue_Columns
       
       PRINT *, 'Reconstructing magnitude distribution from the Catalogue:', trim(Catalogue_Filename), ' to get size distributions'
       !--Read in Catalogue--!                                                                                                                                                                                                                   
       call Catalogue_ReadIn(Cat, trim(Catalogue_Filename), 'COMBO', Catalogue_Columns)
       call convert_Size_from_Pixel_to_Physical(Cat)
       
       if(use_Physical_Sizes) then
          Size_Higher = maxval(Cat%Physical_Sizes); Size_Lower = 0.e0_double!minval(Cat%Physical_Sizes)
       else
          Size_Higher = maxval(Cat%Sizes); Size_Lower = 0.e0_double!minval(Cat%Sizes);
       end if
       if(Apply_Size_Cuts) then
          Size_Lower = maxval((/Size_Cuts(1), Size_Lower/)); Size_Higher = minval( (/Size_Cuts(2), Size_Higher /) )
          print *, 'Applying Size Cuts of:', Size_Lower, Size_Higher
       end if

       allocate(Sizes(nSizes)); Sizes = 0.e0_double
       allocate(SizeBins(nSizes,2)); SizeBins = 0.e0_double
       dSize = (( Size_Higher- Size_Lower )/(1.e0_double*(nSizes-1)) )
       do i = 1, nSizes
          !--Use i-2 so that no galaxies fall into the first bin--!                                                                                                                                                                         
          SizeBins(i,1) = Size_Lower + (i-2)*dSize
          SizeBins(i,2) = SizeBins(i,1) + dSize
          if(i==nSizes) SizeBins(i,2) = SizeBins(i,2) + 1.e-3_double*dSize
          Sizes(i) = 0.5e0_double*(SizeBins(i,1)+SizeBins(i,2))
       end do
       Sizes(1) = SizeBins(1,2)
       print *,  'Got SizesGrid'
       

       !--Bin by equal number density in magnitude bins--!                                                                                                                                                                
       print *, 'Getting Mag Limits', nMags
       call Calculate_Bin_Limits_by_equalNumber(Cat%Mag, nMags, MagBins)
       !--Bin Catalogue by Magnitude--!
       print *, 'Binning By Magnitude:'
       call bin_catalogue_by_magnitude(Cat,MagBins,BCat)

       do i =1 , nMags
          print *, 'Mag Bin ', i, ' : ', MagBins(i,:), ' Variance:', dsqrt(variance_discrete(BCat%Cat(i)%Physical_Sizes, BCat%Cat(i)%Physical_Sizes))
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
