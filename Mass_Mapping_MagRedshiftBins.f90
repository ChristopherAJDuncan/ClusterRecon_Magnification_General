!---Code for the estiamtion of the convergence from Binning by ***ABSOLUTE** magnitude---!
program Mass_Mapping_MagBins
  use Catalogues; use Param_Types; use Convergence_Estimation; use Statistics, only: get_variance, mean_discrete; use MAss_Estimation; use Size_Histograms, only:Size_Histogram_Plotter, Size_Histogram_Catalogue, Size_Histogram_Circular_Aperture
  implicit none

  character(120)::Catalogue_Filename = 'Catalogues/STAGES_COMBO17_gsMag_size_matched.pzcat'
  integer,dimension(13)::Catalogue_Cols = (/-1,10,11,-1,-1,8,-1,-1,-1,14,16,17,5/) !-ntile, RA, Dec, xpos, ypos, Mag, MagErr, Flux, FluxErr, Size, g1, g2, redshift-!       

  character(120)::Output_Directory = 'Mass_Mapping_MagRedshiftBins_Output/'
  character(120)::Output_Filename_AvSize = 'Average_Size_in_RA_Dec_Grid.dat'
  character(120)::Output_Filename_KappaEst = 'Average_Size_KappaEst_in_RA_Dec_Grid.dat'
  character(120)::Output_Filename_KappaError = 'Average_Size_KappaError_in_RA_Dec_Grid.dat'
  character(120)::Output_Filename_SizeError = 'Average_Size_Error_in_RA_Dec_Grid.dat'
  character(120)::Output_Filename_OccGrid = 'nGrid.dat'
  character(120)::Output_Filename_SmoothedOccGrid = 'nGrid_Smoothed.dat' 
  character(120)::Output_FilenameDir_AvSize, Output_FilenameDir_KappaEst, Output_FilenameDir_KappaError, Output_FilenameDir_SizeError, Output_FilenameDir_OccGrid, Output_FilenameDir_SmoothedOccGrid
  character(10)::bin_by_sIze_Type = 'Physical'

  logical:: plot_Bin = .False.
  logical:: get_Size_Histograms = .true.
  character(120)::Bin_Output_Directory
  character(120)::Size_Histogram_Filename
  logical::here

  type(Catalogue)::Cat

  integer::nRA = 128, nDec = 128 !-Default 256-!

  !--Redshift Binning Decalrations--!                                                                                                                                                                                      
  real(double)::lmag, hmag
  type(Binned_Catalogue),allocatable:: BBCat(:)
  integer, parameter::nMagBin = 1, nRedBin = 1 !5, 4
  real(double),allocatable:: MagBinLimits(:,:), RedshiftLimits(:,:)
  character::BinString
  !----------------------------------!                                                                                                                                                                                                        

  !--Average Size Binning---!                  
  !-_Bin labels go as (Redshift, Magnitude, RA, Dec)-!
  real(double)::mean_size_global
  real(double),allocatable::RAGrid(:), DecGrid(:), AvSize(:,:), Kappa(:,:), AverageSize_Bin(:,:,:,:), Kappa_Bin(:,:,:,:)
  integer,allocatable::nGrid(:,:), nGrid_Bin(:,:,:,:)
  real(double),allocatable::Smoothed_nGrid(:,:), Smoothed_nGrid_Bin(:,:,:,:)
  real(double),allocatable::Error_AvSize(:,:), Error_Kappa(:,:), Error_AvSize_Bin(:,:,:,:), Error_Kappa_Bin(:,:,:,:)
  real(double)::beta(nMagBin, nRedBin)

  !--Cluster Information--!                                                                                                                                                                                                                   
  real(double),allocatable::Cluster_Pos(:,:) !-Defines Center of Apertures for Mass Estimation-!                                                                                                                                              
  real(double),allocatable::Cluster_Masses(:), Cluster_Masses_Error(:) !-Cluster-!       

  integer::Bin_Loop_Z, Bin_Loop_M, i
  integer::Occupation_limit = 100.0

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
  call Cut_By_PhotoMetricRedshift(Cat, 0.21e0_double) !--Cut out foreground--!     

  if(get_Size_Histograms) then
     call Size_Histogram_Catalogue(Cat%Physical_Sizes, trim(adjustl(Output_Directory))//'Size_Histogram_FullCatalogue.dat')
     call Size_Histogram_Plotter(trim(adjustl(Output_Directory))//'Size_Histogram_FullCatalogue.dat')
  end if

  
  print *, 'There are:', count(Cat%Redshift < 0.e0_Double), ' unassigned redshifts in the main catalogue'

  call Calculate_Bin_Limits_by_equalNumber(Cat%Redshift, nRedBin, RedshiftLimits)
  call Calculate_Bin_Limits_by_equalNumber(Cat%Mag, nMagBin, MagBinLimits)
  Verbose = .true.
  call bin_catalogue_by_Redshift_and_Magnitude(Cat, RedshiftLimits, MagBinLimits, BBCat)
  Verbose = .false.

  !--Set up beta measurement--!
  beta  =1.e4_double

  !--Set up _Bin arrays, which contain all the information--!
  allocate(AverageSize_Bin(nRedBin,nMagBin, nRA+1, nDec+1)); AverageSize_Bin = 0.e0_double
  allocate(Kappa_Bin(nRedBin,nMagBin, nRA+1, nDec+1)); Kappa_Bin = 0.e0_double
  allocate(nGrid_Bin(nRedBin,nMagBin, nRA+1, nDec+1)); nGrid_Bin = 0
  allocate(Smoothed_nGrid_Bin(nRedBin,nMagBin, nRA+1, nDec+1)); Smoothed_nGrid_Bin = 0.e0_double
  !--Set Errors to default to artificially large number - only changed if analysis is carried out - e.g. in occupation limit--!
  allocate(Error_AvSize_Bin(nRedBin,nMagBin, nRA+1, nDec+1)); Error_AvSize_Bin = 1.e30_double
  allocate(Error_Kappa_Bin(nRedBin,nMagBin, nRA+1, nDec+1)); Error_Kappa_Bin = 1.e30_double

  call construct_RA_Dec_Gridding(Cat, nRA, nDec, RAGrid, DecGrid)
  do Bin_Loop_Z = 1, nRedBin
     do Bin_Loop_M = 1, nMagBin

        print *, 'Doing Bin_Loop:', Bin_Loop_Z, Bin_Loop_M

        !--Check for unassigned redshifts: ERROR HANDLING--!
        if(any(BBCat(Bin_Loop_Z)%Cat(Bin_Loop_M)%Redshift < 0.e0_double)) then           
           print *, 'Of:', size(BBCat(Bin_Loop_Z)%Cat(Bin_Loop_M)%Redshift), ' galaxies in bin, ', count(BBCat(Bin_Loop_Z)%Cat(Bin_Loop_M)%Redshift < 0.e0_double) , 'are less than zero...'
           print *, BBCat(Bin_Loop_Z)%Cat(Bin_Loop_M)%Redshift
           STOP 'Binned catalogue has unassigned redshifts'
        end if

        !--Do not evaluate if there are too few galaxies in the sample--!
        if(BBCat(Bin_Loop_Z)%Occupation(Bin_Loop_M) <= Occupation_Limit) cycle
        
        mean_size_global = global_mean_size(BBCat(Bin_Loop_Z)%Cat(Bin_Loop_M), trim(adjustl(Bin_by_Size_Type)))
        print *, '      Global Mean Size:', mean_size_global
        print *, '      Variance in Sizes:', size_variance(BBCat(Bin_Loop_Z)%Cat(Bin_Loop_M),  trim(adjustl(Bin_by_Size_Type)), mean_size_global)
        
        call Average_Size_in_RA_Dec_Grid(BBCat(Bin_Loop_Z)%Cat(Bin_Loop_M), RAGrid, DecGrid, AvSize, OccupationGrid = nGrid, KappaEst = Kappa, Smoothed_OccupationGrid = Smoothed_nGrid, Size_Type = trim(adjustl(Bin_by_Size_Type)), beta_correction = beta(Bin_Loop_M, Bin_Loop_Z))
        
        call Average_Size_in_RA_Dec_Grid_Errors(BBCat(Bin_Loop_Z)%Cat(Bin_Loop_M), mean_size_global, RAGrid ,DecGrid, Error_AvSize, Error_Kappa, by_Size_Type = trim(adjustl(Bin_by_Size_Type)))


        !--Store for each bin--!
        AverageSize_Bin(Bin_Loop_Z,Bin_Loop_M,:,:) =  AvSize
        Error_AvSize_Bin(Bin_Loop_Z,Bin_Loop_M,:,:) = Error_AvSize
        Kappa_Bin(Bin_Loop_Z,Bin_Loop_M,:,:) = Kappa
        Error_Kappa_Bin(Bin_Loop_Z,Bin_Loop_M,:,:) = Error_Kappa
        nGrid_Bin(Bin_Loop_Z,Bin_Loop_M, :,:) = nGrid
        Smoothed_nGrid_Bin(Bin_Loop_Z,Bin_Loop_M,:,:) = Smoothed_nGrid
        
        deallocate(AvSize, Error_AvSize, Kappa, Error_Kappa, nGrid, Smoothed_nGrid)
        
        
        !-----Mass Estimation---!                                                                            
        call Mass_Estimate_Circular_Aperture_Catalogue(BBCat(Bin_Loop_Z)%Cat(Bin_Loop_M), 0, Cluster_Pos, (/1.e0_double/60.e0_double/))
        PRINT *, 'Reading:'
        READ(*,*)
        call Mass_Estimate_CircularAperture(Kappa_Bin(Bin_Loop_Z,Bin_Loop_M,:,:), Error_Kappa_Bin(Bin_Loop_Z,Bin_Loop_M,:,:), RAGrid, DecGrid, Cluster_Pos, (/1.e0_double/60.e0_double/), Cluster_Masses, Cluster_Masses_Error, Source_Redshift = Mean_Discrete(BBCat(Bin_Loop_Z)%Cat(Bin_Loop_M)%Redshift), Convergence_Map_Occupation = nGrid_Bin(Bin_Loop_Z,Bin_Loop_M, :,:))
        
        print *, '---------------------------------------------------------------------------------------------------------------------'
        print *, 'Cluster Masses are (in 10^18 M_Sun/h):'
        print *, 'A901a:', Cluster_Masses(1), ' +- ', Cluster_Masses_Error(1)
        print *, 'A901b:', Cluster_Masses(2),' +- ',Cluster_Masses_Error(2)
        print *, 'A902:', Cluster_Masses(3),' +- ',Cluster_Masses_Error(3)
        print *, 'SW Group:', Cluster_Masses(4),' +- ',Cluster_Masses_Error(4)
        print *, '---------------------------------------------------------------------------------------------------------------------'
        
        Cluster_Masses = 0.e0_double; Cluster_Masses_Error = 0.e0_double
        deallocate(Cluster_Masses, Cluster_Masses_Error)
        
        !--Output and Plotting--!
        write(BinString, '(I1)') Bin_Loop_Z
        Bin_Output_Directory = trim(adjustl(Output_Directory))//'RedshiftBin'//trim(adjustl(BinString))//'/'
        call system('mkdir '//trim(adjustl(Bin_Output_Directory)))
        write(BinString, '(I1)') Bin_Loop_M
        Bin_Output_Directory = trim(adjustl(Bin_Output_Directory))//trim(adjustl(BinString))//'/'
        call system('mkdir '//trim(adjustl(Bin_Output_Directory)))
        
        Output_FilenameDir_AvSize = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_AvSize))
        Output_FilenameDir_KappaEst = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_KappaEst))
        Output_FilenameDir_KappaError = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_KappaError))
        Output_FilenameDir_SizeError = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_SizeError))
        Output_FilenameDir_OccGrid = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_OccGrid))
        Output_FilenameDir_SmoothedOccGrid = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_SmoothedOccGrid))
        
        call Output_toFile(RAGrid, DecGrid, AverageSize_Bin(Bin_Loop_Z,Bin_Loop_M,:,:), Output_FilenameDir_AvSize, Kappa_Bin(Bin_Loop_Z,Bin_Loop_M,:,:), Output_FilenameDir_KappaEst, Error_Kappa_Bin(Bin_Loop_Z,Bin_Loop_M,:,:), Output_FilenameDir_KappaError, Error_AvSize_Bin(Bin_Loop_Z,Bin_Loop_M,:,:), Output_FilenameDir_SizeError, nGrid_Bin(Bin_Loop_Z,Bin_Loop_M,:,:), Output_FilenameDir_OccGrid, Smoothed_nGrid_Bin(Bin_Loop_Z,Bin_Loop_M,:,:), Output_FilenameDir_SmoothedOccGrid)
        if(plot_Bin) call Run_Plotting_Routine(Output_FilenameDir_AvSize, Output_FilenameDir_KappaEst, Output_FilenameDir_KappaError, Output_FilenameDir_OccGrid)

        if(get_Size_Histograms) then
           Size_Histogram_Filename = trim(adjustl(Bin_Output_Directory))//'Size_Histogram_Catalogue.dat'
           call Size_Histogram_Catalogue(BBCat(Bin_Loop_Z)%Cat(Bin_Loop_M)%Physical_Sizes, Size_Histogram_Filename)
           do i = 1, size(Cluster_Pos,1)
              write(Size_Histogram_Filename, '(I1)') i
              Size_Histogram_Filename = trim(adjustl(Bin_Output_Directory))//'Size_Histogram_Cluster_'//trim(adjustl(Size_Histogram_Filename))//'.dat'
              call Size_Histogram_Circular_Aperture(BBCat(Bin_Loop_Z)%Cat(Bin_Loop_M), Cluster_Pos(i,:), 1.e0_double/60.e0_double, Size_Histogram_Filename)
              call Size_Histogram_Plotter(trim(adjustl(Bin_Output_Directory))//'Size_Histogram_Catalogue.dat', Size_Histogram_Filename)
           end do
        end if
        
        print *, 'Finished Bin (Z,M):', Bin_Loop_Z, Bin_Loop_M
        print *, 'Press Enter to Continue...'
        !read(*,*)
        
     end do
  end do
  
  !--Combine Results--!
  print *, 'Combining Information'
  call Combine_Bin_Information(AverageSize_Bin, Error_AvSize_Bin, Kappa_Bin, Error_Kappa_Bin, nGrid_Bin, Smoothed_nGrid_Bin, AvSize, Error_AvSize, Kappa, Error_Kappa, nGrid, Smoothed_nGrid, Occupation_Limit)
  print *, 'Results combined'
  call Mass_Estimate_CircularAperture(Kappa, Error_Kappa, RAGrid, DecGrid, Cluster_Pos, (/1.e0_double/60.e0_double/), Cluster_Masses, Cluster_Masses_Error, Source_Redshift = Mean_Discrete(Cat%Redshift), Convergence_Map_Occupation = nGrid)

  print *, '---------------------------------------------------------------------------------------------------------------------'
  print *, 'Cluster Masses are (in 10^18 M_Sun/h):'
  print *, 'A901a:', Cluster_Masses(1), ' +- ', Cluster_Masses_Error(1)
  print *, 'A901b:', Cluster_Masses(2),' +- ',Cluster_Masses_Error(2)
  print *, 'A902:', Cluster_Masses(3),' +- ',Cluster_Masses_Error(3)
  print *, 'SW Group:', Cluster_Masses(4),' +- ',Cluster_Masses_Error(4)
  print *, '---------------------------------------------------------------------------------------------------------------------'
  
  Bin_Output_Directory = trim(adjustl(Output_Directory))//'Combined/'
  call system('mkdir '//trim(adjustl(Bin_Output_Directory)))
  
  Output_FilenameDir_AvSize = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_AvSize))
  Output_FilenameDir_KappaEst = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_KappaEst))
  Output_FilenameDir_KappaError = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_KappaError))
  Output_FilenameDir_SizeError = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_SizeError))
  Output_FilenameDir_OccGrid = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_OccGrid))
  Output_FilenameDir_SmoothedOccGrid = trim(adjustl(Bin_Output_Directory))//trim(adjustl(Output_Filename_SmoothedOccGrid))
  
  call Output_toFile(RAGrid, DecGrid, AvSize, Output_FilenameDir_AvSize, Kappa, Output_FilenameDir_KappaEst, Error_Kappa, Output_FilenameDir_KappaError, Error_AvSize, Output_FilenameDir_SizeError, nGrid, Output_FilenameDir_OccGrid, Smoothed_nGrid, Output_FilenameDir_SmoothedOccGrid)
  print *, 'Combined Result output to File, press Enter to plot:'
  read(*,*)
  call Run_Plotting_Routine(Output_FilenameDir_AvSize, Output_FilenameDir_KappaEst, Output_FilenameDir_KappaError, Output_FilenameDir_OccGrid)
  
  print *, 'Finished.'

contains

  subroutine Combine_Bin_Information(AverageSize_Bin, AverageSize_Error_Bin, Kappa_Bin, Kappa_Error_Bin, OccGrid_Bin, Smoothed_OccGrid_Bin, AverageSize, AverageSize_Error, Convergence, Convergence_Error, OccGrid, Smoothed_OccGrid, Occ_Limit)
    !--Combines the infomration from all bins to a global result--!

    real(double),intent(in), dimension(:,:,:,:):: AverageSize_Bin, AverageSize_Error_Bin, Kappa_Bin,Kappa_Error_Bin, Smoothed_OccGrid_Bin
    integer, intent(in):: OccGrid_Bin(:,:,:,:)
    integer,intent(in),optional::Occ_Limit !-If occupation falls below this limit, then not included, or error artificially set to !-inf-!

    real(double),intent(out), dimension(:,:),allocatable :: AverageSize, AverageSize_Error, Convergence, Convergence_Error,  Smoothed_OccGrid 
    integer, intent(out), dimension(:,:),allocatable :: OccGrid

    real(double),allocatable,dimension(:,:,:,:)::Weight_Size_Bin, Weight_Kappa_Bin
    real(double),allocatable, dimension(:,:)::SumWeight_Size, SumWeight_Kappa

    integer::i,j, nBinZ, nBinM

    integer::Method = 2

    allocate(AverageSize(size(AverageSize_Bin,3), size(AverageSize_bin,4))); AverageSize = 0.e0_double
    allocate(AverageSize_Error(size(AverageSize_Error_Bin,3), size(AverageSize_Error_bin,4))); AverageSize_Error = 0.e0_double
    allocate(Convergence(size(Kappa_Bin,3), size(Kappa_bin,4))); Convergence = 0.e0_double
    allocate(Convergence_Error(size(Kappa_Bin,3), size(Kappa_bin,4))); Convergence_Error = 0.e0_double
    allocate(Smoothed_OccGrid(size(Smoothed_OccGrid_Bin,3), size(Smoothed_OccGrid_bin,4))); Smoothed_OccGrid = 0.e0_double
    allocate(OccGrid(size(OccGrid_Bin,3), size(OccGrid_bin,4))); OccGrid = 0
    
    nBinZ =  size(AverageSize_Bin,1); nBinM = size(AverageSize_Bin,2)

    select case (Method)
    case(1) !-Mean across all bins-!
    do i = 1, nBinZ
       do j =1, nBinM
          if(present(Occ_limit) .and. sum(OccGrid_Bin(i,j,:,:)) <= Occ_Limit) cycle
          OccGrid = OccGrid + OccGrid_Bin(i,j,:,:)
          Smoothed_OccGrid = Smoothed_OccGrid + Smoothed_OccGrid_Bin(i,j,:,:)
          
          AverageSize = AverageSize + AverageSize_Bin(i,j,:,:)
          Convergence = Convergence + Kappa_Bin(i,j,:,:)
          
          !--Add Errors in Quadrature--!
          AverageSize_Error = AverageSize_Error + AverageSize_Error_Bin(i,j,:,:)*AverageSize_Error_Bin(i,j,:,:)
          Convergence_Error  = Convergence_Error + Kappa_Error_Bin(i,j,:,:)*Kappa_Error_Bin(i,j,:,:)
       end do
    end do
    !--Convert to Averages--!
    AverageSize = AverageSize/(nBinZ*nBinM)
    Convergence = Convergence/(nBinZ*nBinM)

    AverageSize_Error = dsqrt(AverageSize_Error)/(nBinZ*nBinM)
    Convergence_Error = dsqrt(Convergence_Error)/(nBinZ*nBinM)

    case(2) !-Inverse Variance Weighting--! THIS NEEDS DONE!!!!
       allocate(Weight_Size_Bin(size(AverageSize_Error_Bin,1),size(AverageSize_Error_Bin,2),size(AverageSize_Error_Bin,3), size(AverageSize_Error_Bin,4)))
       allocate(Weight_Kappa_Bin(size(Kappa_Error_Bin,1),size(Kappa_Error_Bin,2),size(Kappa_Error_Bin,3), size(Kappa_Error_Bin,4)))
       allocate(SumWeight_Size(size(Weight_Size_Bin,3), size(Weight_Size_Bin,4))); SumWeight_Size = 0.e0_double
       allocate(SumWeight_Kappa(size(Weight_Kappa_Bin,3),size(Weight_Kappa_Bin,4))); SumWeight_Kappa = 0.e0_double
       do i = 1, nBinZ
          do j = 1, nBinM
             if(present(Occ_limit) .and. sum(OccGrid_Bin(i,j,:,:)) <= Occ_Limit) cycle
             OccGrid = OccGrid + OccGrid_Bin(i,j,:,:)
             Smoothed_OccGrid = Smoothed_OccGrid + Smoothed_OccGrid_Bin(i,j,:,:)
             
             Weight_Size_Bin(i,j,:,:) = 1.e0_double/(AverageSize_Error_Bin(i,j,:,:)**2.)
             Weight_Kappa_Bin(i,j,:,:) = 1.e0_double/(Kappa_Error_Bin(i,j,:,:)**2.)
             
             AverageSize = AverageSize + AverageSize_Bin(i,j,:,:)*Weight_Size_Bin(i,j,:,:)
             Convergence = Convergence + Kappa_Bin(i,j,:,:)*Weight_Kappa_Bin(i,j,:,:)
             
             SumWeight_Size = SumWeight_Size + Weight_Size_Bin(i,j,:,:)
             SumWeight_Kappa = SumWeight_Kappa + Weight_Kappa_Bin(i,j,:,:)
          end do
       end do
       AverageSize_Error = dsqrt(1.e0_double/SumWeight_Size)
       Convergence_Error = dsqrt(1.e0_double/SumWeight_Kappa)

       AverageSize = AverageSize/SumWeight_Size
       Convergence = Convergence/SumWeight_Kappa
       
       deallocate(Weight_Size_Bin, Weight_Kappa_Bin, SumWeight_Size, SumWeight_Kappa)
       
    case default 
       print *, 'Combine_Bin_Information - Invalid method'
    end select


  end subroutine Combine_Bin_Information

  subroutine Run_Plotting_Routine(Filename1, Filename2, Error_Filename2, OccFilename)
    character(*)::Filename1, Filename2, Error_Filename2, OccFilename !-Filename1 is the average size, Filename2 is the Kappa Estimator, bit gridded in RA and Dec                                                                             

    character(50)::Plotter_Name = 'Catalogue_Average_Size_in_RA_Dec_Grid_Plotter.py' !!!                                                                                                                                                     \
                                                                                                                                                                                                                                              
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



end program Mass_Mapping_MagBins
