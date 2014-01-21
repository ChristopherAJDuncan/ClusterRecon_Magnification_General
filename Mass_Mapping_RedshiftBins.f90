program Mass_Mapping_RedshiftBins
  use Catalogues; use Param_Types; use Mass_Estimation; use Statistics, only: Variance_Discrete, Mean_Discrete; use Convergence_Estimation
  implicit none

  character(120)::Catalogue_Filename = 'Catalogues/STAGES_COMBO17_gsMag_size_matched.pzcat'
  integer,dimension(13)::Catalogue_Cols = (/-1,10,11,-1,-1,3,-1,-1,-1,14,16,17,5/) !-ntile, RA, Dec, xpos, ypos, Mag, MagErr, Flux, FluxErr, Size, g1, g2, redshift-!  

  character(120)::Output_Directory = 'Mass_Mapping_RedshiftBins_Output/'
  character(120)::Output_Filename_AvSize = 'Average_Size_in_RA_Dec_Grid.dat'
  character(120)::Output_Filename_KappaEst = 'Average_Size_KappaEst_in_RA_Dec_Grid.dat'
  character(120)::Output_Filename_KappaError = 'Average_Size_KappaError_in_RA_Dec_Grid.dat'
  character(120)::Output_Filename_SizeError = 'Average_Size_Error_in_RA_Dec_Grid.dat'
  character(120)::Output_Filename_OccGrid = 'nGrid.dat'
  character(120)::Output_Filename_SmoothedOccGrid = 'nGrid_Smoothed.dat' 
  character(120)::Output_FilenameDir_AvSize, Output_FilenameDir_KappaEst, Output_FilenameDir_KappaError, Output_FilenameDir_SizeError, Output_FilenameDir_OccGrid, Output_FilenameDir_SmoothedOccGrid
  character(10)::bin_by_sIze_Type = 'Pixel'


  character(120)::Bin_Output_Directory
  logical::here

  type(Catalogue)::Cat
  
  integer::nRA = 256, nDec = 256

  !--Redshift Binning Decalrations--!
  real(double)::lz = 0.21e0_double, hz = 2.0e0_double
  type(Binned_Catalogue):: BCat
  integer, parameter::nRedshiftBin = 1
  real(double)::RedshiftBinLimits(nRedshiftBin, 2)
  character::BinString
  !----------------------------------!

  !--Average Size Binning---!
  real(double)::mean_size_global
  real(double),allocatable::RAGrid(:), DecGrid(:), AvSize(:,:), Kappa(:,:)
  integer,allocatable::nGrid(:,:)
  real(double),allocatable::Smoothed_nGrid(:,:)
  real(double),allocatable::Error_AvSize(:,:), Error_Kappa(:,:)

  !--Cluster Information--!
  real(double),allocatable::Cluster_Pos(:,:) !-Defines Center of Apertures for Mass Estimation-!                      
  real(double),allocatable::Cluster_Masses(:), Cluster_Masses_Error(:) !-Cluster-!

  integer::Bin_Loop
  
  allocate(Cluster_Pos(4,2)); Cluster_Pos = 0.e0_double
  Cluster_Pos(1,:) = (/149.1099e0_double,-9.9561e0_double/)
  Cluster_Pos(2,:) = (/148.9889e0_double,-9.9841e0_double/)
  Cluster_Pos(3,:) = (/149.1424e0_double,-10.1666e0_double/)
  Cluster_Pos(4,:) = (/148.9101e0_double,-10.1719e0_double/)

  inquire(directory = Output_Directory, exist = here)
  if(here==.false.) then 
     print *, 'Output Directory does not exist, please create:', Output_Directory
     STOP
  END if

  inquire(file = Catalogue_Filename, exist = here)
  if(here == .false.) then
     print *, 'Catalogue:', trim(adjustl(Catalogue_Filename)), ' does not exist, stopping..'
     STOP
  end if
  !--Read in the Catalogue with reshifts--!
  call catalogue_readin(Cat, Catalogue_Filename, 'FR', Catalogue_Cols)
  call PSF_Correction(Cat, 1)
  call Cut_By_PhotoMetricRedshift(Cat, lz) !--Cut out foreground--!

  !-Construct Redshift Binning--!
  if(hz > maxval(Cat%Redshift)) hz = maxval(Cat%Redshift)
  do Bin_Loop = 1, size(RedshiftBinLimits,1)
     RedshiftBinLimits(Bin_Loop,1) = lz + (Bin_Loop-1)*(hz-lz)/real(nRedshiftBin)
     RedshiftBinLimits(Bin_Loop,2) = lz + (Bin_Loop)*(hz-lz)/real(nRedshiftBin)
     print *, 'Bin:', Bin_Loop, ' taken between limits:', RedshiftBinLimits(Bin_Loop,:)
  end do
  call bin_catalogue_by_redshift(Cat, RedshiftBinLimits, BCat)

  call construct_RA_Dec_Gridding(Cat, nRA, nDec, RAGrid, DecGrid)
  do Bin_Loop = 1, size(RedshiftBinLimits,1)
     mean_size_global = global_mean_size(BCat%Cat(Bin_Loop), trim(adjustl(Bin_by_Size_Type)))
     print *, 'Bin:',Bin_loop, RedshiftBinLimits(Bin_Loop,:)
     print *, '      Global Mean Size:', mean_size_global
     print *, '      Variance in Sizes:', size_variance(BCat%Cat(Bin_Loop),  trim(adjustl(Bin_by_Size_Type)))

     call Average_Size_in_RA_Dec_Grid(BCat%Cat(Bin_Loop), RAGrid, DecGrid, AvSize, OccupationGrid = nGrid, KappaEst = Kappa, Smoothed_OccupationGrid = Smoothed_nGrid, Size_Type = trim(adjustl(Bin_by_Size_Type)))

     call Average_Size_in_RA_Dec_Grid_Errors(BCat%Cat(Bin_Loop), mean_size_global, RAGrid ,DecGrid, Error_AvSize, Error_Kappa, by_Size_Type = trim(adjustl(Bin_by_Size_Type)))


     !-----Mass Estimation---!
     !-Get Mass Estimates-!                                                                                                 
     call Mass_Estimate_CircularAperture_byPixel(Kappa, Error_Kappa, RAGrid, DecGrid, Cluster_Pos, (/1.e0_double/60.e0_double/), Cluster_Masses, Cluster_Masses_Error, Source_Redshift = Mean_Discrete(BCat%Cat(Bin_Loop)%Redshift))
     
     print *, '---------------------------------------------------------------------------------------------------------------------'
     print *, 'Cluster Masses are (in 10^18 M_Sun/h):'
     print *, 'A901a:', Cluster_Masses(1), ' +- ', Cluster_Masses_Error(1)
     print *, 'A901b:', Cluster_Masses(2),' +- ',Cluster_Masses_Error(2)
     print *, 'A902:', Cluster_Masses(3),' +- ',Cluster_Masses_Error(3)
     print *, 'SW Group:', Cluster_Masses(4),' +- ',Cluster_Masses_Error(4)
     print *, '---------------------------------------------------------------------------------------------------------------------'

     !-- Output and Plotting Routine--!
     write(BinString, '(I1)') Bin_Loop
     Bin_Output_Directory = trim(adjustl(Output_Directory))//'Bin'//trim(adjustl(BinString))//'/'
     call system('mkdir '//trim(adjustl(Bin_Output_Directory)))

     Output_FilenameDir_AvSize = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_AvSize))
     Output_FilenameDir_KappaEst = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_KappaEst))
     Output_FilenameDir_KappaError = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_KappaError))
     Output_FilenameDir_SizeError = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_SizeError))
     Output_FilenameDir_OccGrid = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_OccGrid))
     Output_FilenameDir_SmoothedOccGrid = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_SmoothedOccGrid))

     call Output_toFile(RAGrid, DecGrid, AvSize, Output_FilenameDir_AvSize, Kappa, Output_FilenameDir_KappaEst, Error_Kappa, Output_FilenameDir_KappaError, Error_AvSize, Output_FilenameDir_SizeError, nGrid, Output_FilenameDir_OccGrid, Smoothed_nGrid, Output_FilenameDir_SmoothedOccGrid)
     call Run_Plotting_Routine(Output_FilenameDir_AvSize, Output_FilenameDir_KappaEst, Output_FilenameDir_KappaError, Output_FilenameDir_OccGrid)
     
     deallocate(AvSize, nGrid, Kappa, Smoothed_nGrid, Error_AvSize, Error_Kappa, Cluster_Masses, Cluster_Masses_Error)
     print *, 'Finished Bin:', Bin_Loop
     print *, 'Press Enter to Continue...'
     read(*,*)
     Cluster_Masses = 0.e0_Double; Cluster_Masses_Error = 0.e0_double
  end do

contains

  subroutine Run_Plotting_Routine(Filename1, Filename2, Error_Filename2, OccFilename)
    character(*)::Filename1, Filename2, Error_Filename2, OccFilename !-Filename1 is the average size, Filename2 is the Kappa Estimator, bit gridded in RA and Dec

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


  subroutine Output_toFile(RA, Dec, AverageSize, File_AverageSize, Kappa, File_Kappa,  KappaError, File_KappaError, AverageSizeError, File_AverageSizeError, OccGrid, File_OccGrid, SmoothedOccGrid, File_SmoothedOccGrid)
    real(double),intent(in)::RA(:), Dec(:)
    real(double),intent(in):: AverageSize(:,:), Kappa(:,:), SmoothedOccGrid(:,:), AverageSizeError(:,:), KappaError(:,:)
    integer,intent(in)::OccGrid(:,:)

    character(*), intent(in):: File_AverageSize, File_Kappa, File_KappaError, File_AverageSizeError, File_OccGrid, File_SmoothedOccGrid
    character(5)::fmtstring
    integer::i

    open(82,file = File_AverageSize)
    write(fmtstring,'(I5)') size(Dec)+1
    write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') 0.e0_double, Dec(:)
    write(fmtstring,'(I5)') size(AverageSize,2)+1
    do i = 2, size(AverageSize,1) + 1
       write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') RA(i-1), AverageSize(i-1,:)
    end do
    close(82)
    print *, 'Output to file: ', File_AverageSize
     
    open(56, file = File_Kappa)
    write(fmtstring,'(I5)') size(Dec)+1
    write(56, '('//trim(adjustl(fmtstring))//'(e14.7,x))') 0.e0_double, Dec(:)
    write(fmtstring,'(I5)') size(Kappa,2)+1
    do i = 2, size(Kappa,1) + 1
       write(56, '('//trim(adjustl(fmtstring))//'(e14.7,x))') RA(i-1), Kappa(i-1,:)
    end do
    close(56)
    print *, 'Output to file: ', File_Kappa
    
    open(87, file = File_KappaError)
    write(fmtstring,'(I5)') size(Dec)+1
    write(87, '('//trim(adjustl(fmtstring))//'(e14.7,x))') 0.e0_double, Dec(:)
    write(fmtstring,'(I5)') size(KappaError,2)+1
    do i = 2, size(KappaError,1) + 1
       write(87, '('//trim(adjustl(fmtstring))//'(e14.7,x))') RA(i-1), KappaError(i-1,:)
    end do
    close(87)
    print *, 'Output to file: ', File_KappaError
    
    open(87, file = File_AverageSizeError)
    write(fmtstring,'(I5)') size(Dec)+1
    write(87, '('//trim(adjustl(fmtstring))//'(e14.7,x))') 0.e0_double, Dec(:)
    write(fmtstring,'(I5)') size(AverageSizeError,2)+1
    do i = 2, size(AverageSizeError,1) + 1
       write(87, '('//trim(adjustl(fmtstring))//'(e14.7,x))') RA(i-1), AverageSizeError(i-1,:)
    end do
    close(87)
    print *, 'Output to file: ', File_AverageSizeError

    open(82,file = File_OccGrid)
    write(fmtstring,'(I5)') size(Dec)+1
    write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') 0.e0_double, Dec(:)
    write(fmtstring,'(I5)') size(OccGrid,2)+1
    do i = 2, size(OccGrid,1) + 1
       write(82, '('//trim(adjustl(fmtstring))//'(I5,x))') RA(i-1), OccGrid(i-1,:)
    end do
    close(82)
    print *, 'Output to file: ', File_OccGrid

    open(82,file = File_SmoothedOccGrid)
    write(fmtstring,'(I5)') size(Dec)+1
    write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') 0.e0_double, Dec(:)
    write(fmtstring,'(I5)') size(SmoothedOccGrid,2)+1
    do i = 2, size(SmoothedOccGrid,1) + 1
       write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') RA(i-1), SmoothedOccGrid(i-1,:)
    end do
    close(82)
    print *, 'Output to file: ', File_SmoothedOccGrid

  end subroutine Output_toFile


end program Mass_Mapping_RedshiftBins
