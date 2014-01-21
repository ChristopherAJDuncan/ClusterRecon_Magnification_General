program Catalogue_Average_Size_in_RA_Dec_Grid
  use Catalogues; use Param_Types;
  use Statistics, only: Variance_Discrete
  use Convergence_Estimation
  implicit none

  character(120):: Catalogue_Filename = 'Catalogues/STAGES_shear.cat'!'Catalogues/STAGES_Size_pz_matched.dat'
  character(120)::Output_Dir = 'Catalogue_Average_Size_in_RA_Dec_Grid_Output/Input_Change_TEST2/'
  character(120)::Output_Filename_AvSize = 'Average_Size_in_RA_Dec_Grid.dat'
  character(120)::Output_Filename_KappaEst = 'Average_Size_KappaEst_in_RA_Dec_Grid.dat'
  character(120)::Output_Filename_KappaError = 'Average_Size_KappaError_in_RA_Dec_Grid.dat'
  character(120)::Output_Filename_SizeError = 'Average_Size_Error_in_RA_Dec_Grid.dat'
  character(120)::Output_Filename_OccGrid = 'nGrid.dat'
  character(5)::fmtstring


  call work()


contains

  subroutine Work()
    use Mass_Estimation
  type(Catalogue):: Cat

  integer::Cut_Type = 1 !-0:None, 1:By Minimum Pixel Size, 2:Subtract By Quadrature-!       
  integer::nRA, nDec = 256
  real(double),allocatable::RAGrid(:), DecGrid(:), AvSize(:,:), Kappa(:,:)
  integer,allocatable::nGrid(:,:)
  real(double),allocatable::Smoothed_nGrid(:,:)
  real(double),allocatable::Error_AvSize(:,:), Error_Kappa(:,:)

  real(double)::Magnitude_Cut_lower = 23.e0_double
  real(double)::mean_size_global
  logical::output_size_fluctuation = .false.
  logical::Cut_in_PhotZ = .true.
  character(120)::PhotZ_Filename = 'Catalogues/STAGES_Size_pz_matched.dat'
  type(Catalogue)::PhotZ_Cat
  real(double)::PhotZ_Lower_Cut = 2.1e0_double

  integer::i

  integer::narg
  character(1)::Recalc_String
  logical::Recalculate

  real(double),allocatable::Cluster_Pos(:,:) !-Defines Center of Apertures for Mass Estimation-!
  real(double),allocatable::Cluster_Masses(:), Cluster_Masses_Error(:)

  !-Cluster Positions taken from Heymans 08, (/RA,Dec/)-!
  allocate(Cluster_Pos(4,2)); Cluster_Pos = 0.e0_double
  Cluster_Pos(1,:) = (/149.1099e0_double,-9.9561e0_double/)
  Cluster_Pos(2,:) = (/148.9889e0_double,-9.9841e0_double/)
  Cluster_Pos(3,:) = (/149.1424e0_double,-10.1666e0_double/)
  Cluster_Pos(4,:) = (/148.9101e0_double,-10.1719e0_double/)
!!$


  narg = iargc()
  if(narg == 0) Recalculate = .True.

  call get_command_argument(1, Recalc_String)
  Recalculate = .TRUE.
  if(Recalc_String == 'F' .or. Recalc_String == 'f') Recalculate = .FALSE.

  if(Recalculate) then
     nRA = nDec
     
     call Catalogue_Readin(Cat, Catalogue_Filename, 'FR', (/0/))
     call PSF_Correction(Cat, Cut_Type)
     call Cut_by_Magnitude(Cat, Magnitude_Cut_Lower)
     
     if(Cut_in_PhotZ) then !-Matches to COMBO17 (or other file), assigns redshifts by matching EXACTLY in RA and DEC, and cuts by a lower limit (to remove cluster galaxies)-!
        call Catalogue_Readin(PhotZ_Cat, PhotZ_Filename, 'FR', (/0/)) !-Assumes header file which is correct-!
        call match_Redshift_byCatalogue_byPosition(Cat, PhotZ_Cat)
        call Cut_By_PhotometricRedshift(Cat, LowerCut = PhotZ_Lower_Cut)
     end if

     print *, size(Cat%Redshift),' galaxies passed cuts'

     mean_size_global = global_mean_size(Cat, "Pixel")
     print *, 'Global mean size of catalogue is:', mean_size_global

     print *, 'with variance in size (per galaxy) of:', variance_discrete(Cat%Sizes, Cat%Sizes, mean_size_global, mean_size_global)
     
     call construct_RA_Dec_Gridding(Cat, nRA, nDec, RAGrid, DecGrid)
     call Average_Size_in_RA_Dec_Grid(cat, RAGrid, DecGrid, AvSize, OccupationGrid = nGrid, KappaEst = Kappa, Smoothed_OccupationGrid = Smoothed_nGrid)

     call Average_Size_in_RA_Dec_Grid_Errors(Cat, mean_size_global, RAGrid ,DecGrid, Error_AvSize, Error_Kappa)     

     
     !-Note Smoothing exists in the routine which does the binning also--!
!!$  call smooth_gaussiankernal(AvSize, Sigma = (/0.75e0_double/60., 0.75e0_double/60./), Box_Length = (/RAGrid(size(RAGrid))-RAGrid(1), DecGrid(size(DecGrid))-DecGrid(1)/))
!!$  call smooth_gaussiankernal(Kappa, Sigma = (/0.75e0_double/60., 0.75e0_double/60./), Box_Length = (/RAGrid(size(RAGrid))-RAGrid(1), DecGrid(size(DecGrid))-DecGrid(1)/))       

    !-Output OccupationGrid-!          
    open(82,file = trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_OccGrid)))
    write(fmtstring,'(I5)') size(DecGrid)+1
    write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') 0.e0_double, DecGrid(:)
    write(fmtstring,'(I5)') size(NGrid,2)+1
    do i = 2, size(NGrid,1) + 1
       write(82, '('//trim(adjustl(fmtstring))//'(I5,x))') RAGrid(i-1), NGrid(i-1,:)
    end do
    close(82)

    print *, 'Occupation:'
    print *, ' #(n>1):', count(NGrid>1), ' #(n==1):', count(NGrid==1), ' #(n<1)[unoccupied]:', count(NGrid<1)
    print *, ' #(n==2):', count(NGrid==2), ' #(n==3):', count(NGrid==3), ' #(n==5):', count(NGrid==5)

    
    print *, 'Output to file: ', trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_OccGrid))

    
    !--Ouput Smoothed nGrid--!
    open(82,file = trim(adjustl(Output_Dir))//'nGrid_Smoothed.dat')
    write(fmtstring,'(I5)') size(DecGrid)+1
    write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') 0.e0_double, DecGrid(:)
    write(fmtstring,'(I5)') size(Smoothed_nGrid,2)+1
    do i = 2, size(NGrid,1) + 1
       write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') RAGrid(i-1), Smoothed_nGrid(i-1,:)
    end do
    close(82)
    print *, 'Output to file: ', trim(adjustl(Output_Dir))//'nGrid_Smoothed.dat'

     
     !-Output Average Size-!
     open(82,file = trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_AvSize)))
     write(fmtstring,'(I5)') size(DecGrid)+1
     write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') 0.e0_double, DecGrid(:)
     write(fmtstring,'(I5)') size(AvSize,2)+1
     do i = 2, size(AvSize,1) + 1
        write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') RAGrid(i-1), AvSize(i-1,:)
     end do
     close(82)
     print *, 'Output to file: ', trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_AvSize))
     
     print *, 'Maxval of Kappa before output:', maxval(Kappa)
     
     !-Output Kappa Estimator--!
     open(56, file = trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_KappaEst)))
     write(fmtstring,'(I5)') size(DecGrid)+1
     write(56, '('//trim(adjustl(fmtstring))//'(e14.7,x))') 0.e0_double, DecGrid(:)
     write(fmtstring,'(I5)') size(Kappa,2)+1
     do i = 2, size(Kappa,1) + 1
        write(56, '('//trim(adjustl(fmtstring))//'(e14.7,x))') RAGrid(i-1), Kappa(i-1,:)
     end do
     close(56)
     
     print *, 'Output to file: ', trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_KappaEst))
     
     !-Output Kappa Error--!
     open(87, file = trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_KappaError)))
     write(fmtstring,'(I5)') size(DecGrid)+1
     write(87, '('//trim(adjustl(fmtstring))//'(e14.7,x))') 0.e0_double, DecGrid(:)
     write(fmtstring,'(I5)') size(Error_Kappa,2)+1
     do i = 2, size(Kappa,1) + 1
        write(87, '('//trim(adjustl(fmtstring))//'(e14.7,x))') RAGrid(i-1), Error_Kappa(i-1,:)
     end do
     close(87)
     
     print *, 'Output to file: ', trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_KappaError))


     !--Output Average Size Error--!
     open(87, file = trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_SizeError)))
     write(fmtstring,'(I5)') size(DecGrid)+1
     write(87, '('//trim(adjustl(fmtstring))//'(e14.7,x))') 0.e0_double, DecGrid(:)
     write(fmtstring,'(I5)') size(Error_Kappa,2)+1
     do i = 2, size(Kappa,1) + 1
        write(87, '('//trim(adjustl(fmtstring))//'(e14.7,x))') RAGrid(i-1), Error_AvSize(i-1,:)
     end do
     close(87)
     
     print *, 'Output to file: ', trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_SizeError))

  end if

  !-Get Mass Estimates-!
  call Mass_Estimate_CircularAperture_byPixel(Kappa, Error_Kappa, RAGrid, DecGrid, Cluster_Pos, (/1.e0_double/60.e0_double/), Cluster_Masses, Cluster_Masses_Error)

  print *, '---------------------------------------------------------------------------------------------------------------------'
  print *, 'Cluster Masses are (in 10^18 M_Sun/h):'
  print *, 'A901a:', Cluster_Masses(1), ' +- ', Cluster_Masses_Error(1)
  print *, 'A901b:', Cluster_Masses(2),' +- ',Cluster_Masses_Error(2)
  print *, 'A902:', Cluster_Masses(3),' +- ',Cluster_Masses_Error(3)
  print *, 'SW Group:', Cluster_Masses(4),' +- ',Cluster_Masses_Error(4)
  print *, '---------------------------------------------------------------------------------------------------------------------'


 
  call Run_Plotting_Routine(trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_AvSize)), trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_KappaEst)), trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_KappaError)), trim(adjustl(Output_Dir))//trim(adjustl(Output_Filename_OccGrid)))
    

  end subroutine Work

  subroutine Run_Plotting_Routine(Filename1, Filename2, Error_Filename2, OccFilename)
    character(*)::Filename1, Filename2, Error_Filename2, OccFilename

    character(50)::Plotter_Name = 'Catalogue_Average_Size_in_RA_Dec_Grid_Plotter.py' !!!
    logical::Here

    inquire(file = Filename1, exist = here)
    if(here == .false.) print *, 'Error in calling plotting routine: Filename 1 does not exist'
    inquire(file = Filename2, exist = here)
    if(here == .false.)print *, 'Error in calling plotting routine: Filename 2 does not exist'

    if(here) then
       print *, 'Running Plotting routine'
       call system('python '//trim(adjustl(Plotter_Name))//' '//trim(adjustl(Filename1))//' '//trim(adjustl(Filename2))//' '//trim(adjustl(Error_Filename2))//' '//trim(adjustl(OccFilename)))
    end if

  end subroutine Run_Plotting_Routine

end program Catalogue_Average_Size_in_RA_Dec_Grid
