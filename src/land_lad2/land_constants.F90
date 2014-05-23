module land_constants_mod

use constants_mod, only : rdgas, rvgas, wtmair

implicit none
private

! ==== public interfaces =====================================================
integer, public, parameter :: &
     NBANDS   = 2, & ! number of spectral bands for short-wave radiation calculations
     BAND_VIS = 1, & ! visible radiation (wavelenght range?)
     BAND_NIR = 2    ! near infra-red radiation (wavelenght range?)

real, public, parameter :: d622 = rdgas/rvgas
real, public, parameter :: d378 = 1.0-d622
real, public, parameter :: d608 = d378/d622

real, public, parameter :: Rugas = 8.314472 ! universal gas constant, J K-1 mol-1

real, public, parameter :: seconds_per_year = 86400.0*365.0
real, public, parameter :: mol_C = 12.0e-3 ! molar mass of carbon, kg/mol
real, public, parameter :: mol_air = wtmair/1000.0 ! molar mass of air, kg/mol
real, public, parameter :: mol_CO2 = 44.00995e-3 ! molar mass of CO2,kg/mol
real, public, parameter :: mol_h2o = 18.0e-3 ! molar mass of water, kg/mol

real, public, parameter :: Pa_per_m = 10000 ! actually 9800
real, public, parameter :: MPa_per_m = 0.01 ! actually 0.0098
real, public, parameter :: mmol_per_kg_h2o = 1000/18*1000
real, public, parameter :: mmol_per_kg_air = 1000/29*1000


! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     version = '$Id: land_constants.F90,v 1.1 2010/08/15 14:41:58 pjp Exp $', &
     tagname = '$Name: no_fms_b_pjp $'

end module
