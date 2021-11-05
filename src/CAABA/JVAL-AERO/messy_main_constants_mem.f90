!****************************************************************************
!                Time-stamp: <2018-07-04 12:39:38 joec_pa>
!****************************************************************************

! Definitions of machine precision constants as Fortran PARAMETERs for MESSy
! Definitions of physical constants as Fortran PARAMETERs for MESSy

! Authors:
! Rolf Sander,     MPICH, 2004: original code
! Patrick Joeckel, MPICH, 2004: preprocessor-directives removed; the 
!                               BASEMODEL now may use the constants of the
!                               Modular Earth Submodel System ...

MODULE messy_main_constants_mem

  IMPLICIT NONE
  INTRINSIC :: SELECTED_INT_KIND, SELECTED_REAL_KIND, TINY
  PUBLIC !is already default

  CHARACTER(LEN=*), PARAMETER :: modstr = 'MESSy'
#ifdef _VCSREV_
  CHARACTER(LEN=*), PARAMETER :: modver = 'd2.53.0.24'//_VCSREV_
#else
  CHARACTER(LEN=*), PARAMETER :: modver = 'd2.53.0.24'
#endif

  ! DISTRIBUTION MODIFICATION (please do NOT modify the next line!)
  ! OPTIONS: #DAU-NGRADE#

  ! MACHINE PRECISION CONSTANTS
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)
  INTEGER, PARAMETER :: wp = dp
  REAL(DP), PARAMETER :: TINY_DP = TINY(0._dp) ! mz_rs_20060114
  REAL(DP), PARAMETER :: HUGE_DP = HUGE(0._dp) ! mz_rs_20100409
  ! mz_rs_20101124+
  ! practically, for many applications, HUGE_DP and BIG_DP are both like
  ! infinity. Sometimes, however, BIG_DP creates less numerical
  ! problems, e.g. when it is necessary to calculate BIG_DP**2:
  REAL(DP), PARAMETER :: BIG_DP = 1E40_dp
  ! mz_rs_20101124-
  !
  INTEGER, PARAMETER :: nerr = 6 ! mz_ab_20100503: for use in mo_mpi.f90
  ! op_bk_20160929+
  ! introduce an "abstract global patch", set to 0, this should be a saver way,
  ! than put "global" channels in patch 1. ICON patches can also be inactive ->
  ! this won't be the case for patch 0.
  INTEGER, PARAMETER :: global_patch = 0
  ! op_bk_20160929-

  ! FLAGS
  REAL(DP), PARAMETER :: FLAGGED_BAD = -1.0E+34_dp  ! FERRET

  ! mz_rs_20070904+
  ! STRINGS FOR UNIFORM OUTPUT:
  CHARACTER(LEN=*), PARAMETER :: HLINE1 = &
    '*************************************'// &
    '*************************************'
  CHARACTER(LEN=*), PARAMETER :: HLINE2 = &
    '-------------------------------------'// &
    '-------------------------------------'
  CHARACTER(LEN=*), PARAMETER :: HLINE3 = &
    '.....................................'// &
    '.....................................'
  ! mz_rs_20070904-

  ! STRING LENGTHs
  INTEGER, PARAMETER :: STRLEN_SHORT  = 8
  INTEGER, PARAMETER :: STRLEN_MEDIUM = 24
  INTEGER, PARAMETER :: STRLEN_LONG   = 64
  INTEGER, PARAMETER :: STRLEN_VLONG  = 80
  INTEGER, PARAMETER :: STRLEN_ULONG  = 256
  INTEGER, PARAMETER :: STRLEN_XLONG  = 512 ! mz_ab_20130227
  ! mz_rs_20100331+
  ! I'm not sure if 15 is really the upper limit for the length of
  ! KPP species. However, we currently don't have any species
  ! with more than 15 characters, and that works fine...
  INTEGER, PARAMETER :: STRLEN_KPPSPECIES =  15
  ! mz_rs_20100331-

  ! PHYSICAL CONSTANTS (CODATA Recommended Values, 2010, 
  ! http://physics.nist.gov/cuu/Constants/)
  REAL(dp), PARAMETER :: pi       = 3.14159265358979323846_dp
  REAL(dp), PARAMETER :: pi_icon  = 3.14159265358979323846264338327950288_wp
  REAL(dp), PARAMETER :: R_gas    = 8.3144621_dp      ! R [J/K/mol]
  REAL(dp), PARAMETER :: h_Planck = 6.62606957E-34_dp ! Planck constant [Js]
  REAL(dp), PARAMETER :: c_light  = 2.99792458E8_dp   ! speed of light [m/s]
  REAL(dp), PARAMETER :: stbo     = 5.670373E-8_dp    ! Stephan-Boltzmann constant [W/m2/K4]
  REAL(dp), PARAMETER :: N_A      = 6.02214129E23_dp ! Avogadro constant [1/mol]
  REAL(dp), PARAMETER :: N_A_kmol = 6.02214129E26_dp ! Avogadro constant [1/kmol] ! mz_ab_20130718
#ifndef MESSYIDTC
  REAL(dp), PARAMETER :: g        = 9.80665_dp   ! gravity acceleration [m/s2]
#else
  REAL(dp), PARAMETER :: g        = 9.80616_dp   ! gravity acceleration [m/s2]
#endif
  REAL(dp), PARAMETER :: T0       = 298.15_dp        ! standard temperature [K]
  REAL(dp), PARAMETER :: T0_INV   = 1._DP / T0       ! 1/T0 [1/K]
#ifndef MESSYIDTC
  REAL(dp), PARAMETER :: atm2Pa   = 101325._dp   ! conversion from [atm] to [Pa]
#else
  REAL(dp), PARAMETER :: atm2Pa   = 100000._dp   ! conversion from [atm] to [Pa]
#endif
  REAL(dp), PARAMETER :: cal2J    = 4.1868_dp        ! conversion from [cal] to [J]
  REAL(dp), PARAMETER :: k_B      = 1.3806488E-23_dp ! Boltzmann constant [J/K]
  REAL(dp), PARAMETER :: c_vKar   = 0.4_dp           !  Karman constant [?]

  ! standard atmosphere vertical gradient of the temperature in the troposphere 
  REAL(dp), PARAMETER :: stDTDZ = 0.0065_dp  ! [K/m] ! um_ak_20120103

  ! MXXX = molar mass of element XXX [g/mol]
  REAL(dp), PARAMETER :: MH  =   1.01_dp
! op_ff_20141016+
  REAL(dp), PARAMETER :: MD    =   2.01_dp
  REAL(dp), PARAMETER :: M13C  =  13.00_dp
  REAL(dp), PARAMETER :: M12C  =  12.00_dp
! op_ff_20141016-
  REAL(dp), PARAMETER :: MC  =  12.01_dp
  REAL(dp), PARAMETER :: MN  =  14.01_dp
  REAL(dp), PARAMETER :: M18O=  18.00_dp
  REAL(dp), PARAMETER :: MF  =  19.00_dp
  REAL(dp), PARAMETER :: MNa =  22.99_dp
  REAL(dp), PARAMETER :: MO  =  16.00_dp
  REAL(dp), PARAMETER :: MS  =  32.07_dp
  REAL(dp), PARAMETER :: MCl =  35.45_dp
  REAL(dp), PARAMETER :: MBr =  79.90_dp
  REAL(dp), PARAMETER :: MI  = 126.90_dp
  REAL(dp), PARAMETER :: MHg = 200.59_dp
! mz_ht_20160302+
  REAL(dp), PARAMETER :: MK  =  39.10_dp
  REAL(dp), PARAMETER :: MMg =  24.31_dp
  REAL(dp), PARAMETER :: MCa =  40.08_dp
! mz_ht_20160302-
  ! M_XXX = molar mass of compounds [g/mol]
  REAL(dp), PARAMETER :: M_O3  = MO*3._dp      ! molar mass of ozone [g/mol]
  REAL(dp), PARAMETER :: M_O2  = MO*2._dp      ! molar mass of oxygen [g/mol]
  REAL(dp), PARAMETER :: M_H2O = MH*2._dp + MO ! molar mass of H2O [g/mol]
  REAL(dp), PARAMETER :: M_HDO = MH + MD + MO  ! molar mass of HDO [g/mol]
  REAL(dp), PARAMETER :: M_HH18O=MH*2._dp+ M18O! molar mass of HH18O [g/mol]
  REAL(dp), PARAMETER :: M_N2  = MN*2._dp      ! molar mass of N2 [g/mol]

  ! DRY AIR AND WATER VAPOUR THERMODYNAMIC CONSTANTS
  REAL(dp), PARAMETER :: tmelt   = 273.15_dp    ! melting temp. of ice/snow [K]
  REAL(dp), PARAMETER :: ttrip   = 273.16_dp    ! triple point of water [K]
  REAL(dp), PARAMETER :: rho_H2O = 999.97_dp    ! density of H2O [kg/m3]
  ! fb_mk_20101021+
  REAL(dp), PARAMETER :: rho_sea = 1025._dp   ! density of sea water in [kg/m3]
  ! fb_mk_20101021-
  ! mz_rj_20131022+
  REAL(dp), PARAMETER :: rho_air = 1.225_dp   ! standard density of air [kg/m3]
  ! mz_rj_20131022-
  REAL(dp), PARAMETER :: M_air   = 28.970_dp    ! molar mass of dry air [g/mol]
#ifndef MESSYIDTC
  REAL(dp), PARAMETER :: cp_air  = 1005.46_dp   ! specific heat of dry air at
                                                ! constant pressure [J/K/kg]
#else
  REAL(dp), PARAMETER :: cp_air  = 1004.64_dp   ! specific heat of dry air at
                                                ! constant pressure [J/K/kg]
#endif
  ! mz_ap_20090519+
  REAL(dp), PARAMETER :: alv   = 2.5008e6_dp    ! latent heat for vaporisation 
  !                                             ! [J/kg]
  REAL(dp), PARAMETER :: als   = 2.8345e6_dp    ! latent heat for sublimation
  !                                             ! [J/kg]
  REAL(dp), PARAMETER :: alf   = als-alv        ! latent heat for fusion [J/kg]
  
  ! mz_ap_20090519-

  ! gas constant for dry air [J/K/kg]
#ifndef MESSYIDTC
  REAL(dp), PARAMETER :: rd      = 1000._dp * R_gas/M_air ! 287.05_dp
#else
  REAL(dp), PARAMETER :: rd      = 287.04_dp
#endif
  ! gas constant for water vapour
  REAL(dp), PARAMETER :: rv      = 1000._dp * R_gas/M_H2O ! 461.51_dp
  ! specific heat of water vapour at constant pressure [J/K/kg]
  REAL(dp), PARAMETER :: cpv     = 1869.46_dp
  ! dimensionless auxiliary constants
! op_re_20130718+
!!$  REAL(dp), PARAMETER :: vtmpc1  = rv/rd-1.0_dp
  REAL(dp), PARAMETER :: vtmpc1  = M_air / M_H2O - 1.0_dp
! op_re_20130718-
  REAL(dp), PARAMETER :: vtmpc2  = cpv/cp_air-1.0_dp
  REAL(dp), PARAMETER :: MM_eps  = M_H2O/M_air ! mz_hr_20070323

  ! cloud and radiation
  REAL(dp), SAVE     :: ceffmin = 10.0_dp    ! min eff.radius for ice cloud
  REAL(dp),PARAMETER :: ceffmax = 150.0_dp   ! max eff.radius for ice cloud
  REAL(dp), SAVE     :: ccwmin  = 1.0e-7_dp  ! cloud water limit for cover>0
  REAL(dp),PARAMETER :: cemiss  = 0.996_dp   ! LW emissivity 

  ! vdiff etc.
  REAL(dp), PARAMETER :: cvdifts = 1.5_dp  !  *factor for timestep weighting
  !                                        !   in *rhs* of *vdiff* and *scv*.

  ! PLANETARY PARAMETERS
#ifndef MESSYIDTC
  REAL(dp), PARAMETER :: radius_earth = 6371000.0_dp ! radius of the Earth [m]
#else
  REAL(dp), PARAMETER :: radius_earth = 6371299.0_dp ! radius of the Earth [m]
#endif
  REAL(dp), PARAMETER :: OneDay       = 86400.0_dp   ! one day [s]
  REAL(dp), PARAMETER :: OneSiderialDay = 86164.0_dp ! one day [s] ! mz_ab_20130129
#ifndef MESSYIDTC
  REAL(dp), PARAMETER :: omega = 0.7292E-4_dp   ! solid rotation velocity of the
                                                ! earth in 1/s ! mz_ab_20130129
#else
  REAL(dp), PARAMETER :: omega = 0.729212E-4_dp ! solid rotation velocity of the
                                                ! earth in 1/s ! mz_ab_20130129
#endif

  ! fu_kk_20061002+
  REAL(dp), PARAMETER :: solc  = 1365.0_dp           ! solar constant [W/m2]
  !REAL(dp), PARAMETER :: solc  = 1365.41_dp          ! solar constant [W/m2]
  ! fu_kk_20061002-
  ! *ratio: atmospheric height/radius of the earth.
  REAL(dp), PARAMETER :: crae = 0.1277e-02_dp

  ! fb_mk_20100212+
!!$  REAL(dp), PARAMETER :: alv   = 2.5008e6_dp ! latent heat for vaporisation 
!!$  !                                          ! in J/kg
!!$  REAL(dp), PARAMETER :: als   = 2.8345e6_dp ! latent heat for sublimation
!!$  !                                          ! in J/kg
!!$  REAL(dp), PARAMETER :: alf   = als-alv     ! latent heat for fusion in J/kg
  REAL(dp), PARAMETER :: clw   = 4186.84_dp  ! specific heat for liquid water
  !                                          ! J/K/kg
  REAL(dp), PARAMETER :: csw   = 3994._dp    ! specific heat for sea water
  !                                          ! J/K/kg
  REAL(dp), PARAMETER :: ctfreez = 271.38_dp ! temperature at which sea
                                             ! starts freezing/melting
  ! fb_mk_20100212-

  ! mz_ab_20090525+
  REAL(dp), PARAMETER:: AM = 1.673e-27     ! Atomic mass unit
  REAL(dp), PARAMETER:: ELCH =  1.602E-19  ! Electron charge

  REAL(dp), PARAMETER:: TWOPI = pi*2._dp      ! Pi*2.
  REAL(dp), PARAMETER:: PI_2  = pi*0.5_dp     ! Pi/2.
  REAL(dp), PARAMETER:: DTR   = pi/180._dp    ! Degrees to radians
  REAL(dp), PARAMETER:: RTD   = 180._dp/pi    ! Radians to degrees
  ! mz_ab_20090525-
  ! ub_ak_20170425+
  REAL(dp), PARAMETER:: RTD_icon = 180._dp/pi_icon    ! Radians to degrees
  REAL(dp), PARAMETER:: DTR_icon = pi_icon/180._dp    ! Radians to degrees
  ! ub_ak_20170425-

END MODULE messy_main_constants_mem

!*****************************************************************************
