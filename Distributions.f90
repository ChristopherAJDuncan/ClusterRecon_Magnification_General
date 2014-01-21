module Distributions
  use Param_Types
implicit none

character(200),private::Reference_Catalogue = 'Catalogues/STAGES_COMBO17_gsMag_size_matched.pzcat'
integer::Reference_Catalogue_Columns(13) = (/-1,10,11,-1,-1,8,-1,-1,-1,14,16,17,5/)

contains

  subroutine get_Size_Distribution_MagnitudeBinning_byCatalogue(MagBins, Sizes, PDF, Output_Dir)
    !--Edit to also return from file input
    use Catalogues; use Statistics, only: variance_discrete
       real(double), intent(out),allocatable:: Sizes(:)
       real(double), intent(out),allocatable::pdf(:,:) !-Mag Bin, GridValue-!
       real(double),dimension(:,:),allocatable,intent(inout)::MagBins
       character(200), intent(in), optional:: Output_Dir

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
