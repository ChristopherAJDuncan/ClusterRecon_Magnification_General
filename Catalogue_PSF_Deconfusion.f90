program Catalogue_PSF_Deconfusion
  use Catalogues; use Param_Types
  implicit none

  type(Catalogue):: Cat
  character(120):: Catalogue_Filename = 'Catalogues/STAGES_shear.cat'

  integer::Cut_Type = 1 !-0:None, 1:By Minimum Pixel Size, 2:Subtract By Quadrature-!

  call Catalogue_Readin(Cat, Catalogue_Filename, 'FR')
  call PSF_Correction(Cat, Cut_Type)

  !-Output-!
  call Catalogue_Output(Cat, 'Catalogues/STAGES_PSF_Corrected_byPixelSize.cat', .true.)

end program Catalogue_PSF_Deconfusion
