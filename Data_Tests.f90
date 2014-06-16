module Data_Tests
  !---Contains routines that test the sample and data, which will be used in a given analysis---!
  use Catalogues; use Param_Types
  implicit none
  
  
contains

  subroutine Foreground_Contamination_NumberDensity(Cat, Ap_Pos, Output_Directory)
    type(Catalogue), intent(in):: Cat
    real(double), intent(in):: Ap_Pos(:,:) !--in DEGREES; Apeture, RA/Dec-!
    character(*), intent(in):: Output_Directory

    integer:: Ap, i

    real(double):: Core_Radius = 0.e0_double
    real(double), allocatable:: NumberDensity(:,:) !-Ap, Seperation-!
    real(double), allocatable:: SeperationGrid(:), SeperationBins(:,:)
    integer:: nSep = 20
    real(double),parameter:: Seperation_Limits(2) = (/0.e0_double, 3.e0_double/) !--In ArcMinutes--!
    real(double),allocatable:: Annulus_Area(:)

    type(Catalogue):: Annulus_Cat

    character(3)::fmtstring

    !--Set up number density grid--!
    allocate(SeperationGrid(nSep)); SeperationGrid = 0.e0_double
    allocate(SeperationBins(nSep, 2)); SeperationBins = 0.e0_double
    allocate(Annulus_Area(size(SeperationGrid))); Annulus_Area = 0.e0_double
    do i =1 , nSep
       if(i == 1) then
          SeperationBins(i,1) = Seperation_Limits(1)
       else
          SeperationBins(i,1) = SeperationBins(i-1,2)
       end if
       SeperationBins(i,2) = SeperationBins(i,1) + ((Seperation_Limits(2)-Seperation_Limits(1))/nSep)

       SeperationGrid(i) = 0.5e0_double*( SeperationBins(i,2) + SeperationBins(i,1) )
       Annulus_Area(i) = 3.14159e0_double*( (SeperationBins(i,2)**2.e0_double) - (SeperationBins(i,1)**2.e0_double) ) !--Pi*r^2--!
    end do

    allocate(NumberDensity(size(Ap_Pos,1), size(SeperationGrid))); NumberDensity = 0.e0_double
    do Ap = 1, size(Ap_Pos,1)
       do i = 1, size(SeperationGrid)
          call Identify_Galaxys_in_Circular_Aperture(Cat, Ap_Pos(Ap,:), SeperationBins(i,2)/60.e0_double, Annulus_Cat, Core_Radius = SeperationBins(i,1)/60.e0_double)

          NumberDensity(Ap,i) = (1.e0_double*size(Annulus_Cat%RA))/Annulus_Area(i)

          call Catalogue_Destruct(Annulus_Cat)
       end do
    end do

    !---Output---!
    open(unit = 31, file = trim(Output_Directory)//'Foreground_Contamination_NumberDensity.dat')
    !--Header--!
    write(31, '(A)') '#1,2:Seperation Bin [Annulus Limits], 3+: Number Density of galaxies in annulus'
    write(fmtstring, '(I2)') size(Ap_Pos)+2
    do i = 1, size(SeperationGrid)
       write(31, '('//trim(fmtstring)//'(e9.3,x))') SeperationBins(i,:), NumberDensity(:,i)
    end do
    close(31)

    write(*,'(2A)') '*** Cluster Contamination Test output to: ', trim(Output_Directory)//'Foreground_Contamination_NumberDensity.dat'
    
    deallocate(Annulus_Area,SeperationGrid,SeperationBins,NumberDensity)

  end subroutine Foreground_Contamination_NumberDensity


end module Data_Tests
