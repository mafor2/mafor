!*****************************************************************************
!                Time-stamp: <2018-07-19 12:30:43 sander>
!*****************************************************************************

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program; if not, get it from:
! http://www.gnu.org/copyleft/gpl.html

!*****************************************************************************

!OPTION! -O nomove
MODULE messy_jval

  ! Photolysis routines to calculate J values
  ! Original code by Jochen Landgraf (MPICH, until 1998), see:
  !   Landgraf and Crutzen
  !   An efficient method for online calculations of photolysis and
  !   heating rates
  !   J. Atmos. Sci. 55, 863-878 (1998)

  ! submodel maintainer: Rolf Sander
  ! major contributions by Roland von Glasow, Patrick Joeckel,
  ! and Astrid Kerkweg
  ! see CHANGELOG for a detailed list of modifications
  
  ! -----------------------------------------------------
  ! MSK CHANGES FOR INTERFACING WITH THE MAFOR MODEL:
  ! 06.11.2020 LINE 1594
  !  !CALL jval_cal_uv
  !  IF (l_heating)  CALL jval_cal_uv
  ! -----------------------------------------------------
  !

  USE messy_main_constants_mem, ONLY: DP, k_B, HLINE2
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY, PTR_2D_ARRAY
  USE messy_cmn_photol_mem      ! IP_MAX, ip_*, jname

  IMPLICIT NONE
  PRIVATE
  ! NAME OF SUBMODEL
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'jval'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '14.2'

  PUBLIC :: aerosol_data
  PUBLIC :: combine_o3_fields
  PUBLIC :: jvalues
  PUBLIC :: jval_read_nml_ctrl
  PUBLIC :: jval_solar_time_control

  ! NAMELIST
  ! WHICH QUANTUM YIELD FOR CH3COCH3
  INTEGER, PUBLIC, SAVE :: qy_CH3COCH3 = 1
  ! SOLAR CYCLE DATA LINEARLY INTERPOLATED IN TIME
  ! AK: This needs to be dp, otherwise in ECHAM no p_bcast routine available
  REAL(DP), PUBLIC, SAVE  :: r_sol = 0.5_dp

  LOGICAL, PUBLIC, DIMENSION(IP_MAX), SAVE :: lp = .FALSE.

  ! pointer to photolysis rate coeff.
  TYPE(PTR_2D_ARRAY), PUBLIC, DIMENSION(:), POINTER, SAVE  :: jval_2d

  ! heating rates
  REAL(dp), PUBLIC, DIMENSION(:,:), POINTER :: rh_o2_2d, rh_o3_2d
  REAL(dp), PUBLIC, DIMENSION(:,:), POINTER :: fhuv_2d, fhuvdna_2d

  INTEGER, PARAMETER, PUBLIC :: MAXWAV =  7 ! for arrays for intervals 1-7
  INTEGER, PARAMETER :: dim55  = 55 ! dimension for ai_* and bi_*
  INTEGER, PARAMETER :: dim58  = 58 ! dimension for a0_*

  REAL :: &
    aext(MAXWAV,8,4),  & ! extinction coefficient [1/km]
    asca(MAXWAV,8,4),  & ! effective scattering cross section [cm^2/part.]
    aabs(MAXWAV,8,4),  & ! effective absorption cross section [cm^2/part.]
      ag(MAXWAV,8,4)     ! effective asymmetry factor
    !             ^
    !             1 = rural aerosol
    !             2 = maritime aerosol
    !             3 = urban aerosol
    !             4 = free troposphere aerosol

  ! reference: atlas2 spectrum and wmo86
  REAL, SAVE :: flux(MAXWAV) ! extraterrestic flux
  REAL, PARAMETER :: crray(MAXWAV) = &  ! rayleigh cross section [cm^2]
    (/ 3.1995E-25, 6.7791E-26, 5.4920E-26, 4.9733E-26, &
    4.2781E-26, 2.3131E-26, 3.6810E-27 /)
  ! integrated extraterrestric flux / extraterrestric flux (lambda(i))
  ! reference:  atlas2 spectrum
  REAL, SAVE :: f0(MAXWAV)
  ! integrated flux over the Schumann-Runge bands at toa
  REAL, SAVE :: SR_toa_flux
  REAL, SAVE :: phi_la

  ! **************************************************************************

  ! Some variables for subroutine jvalues are now module variables. This
  ! allows data interchange between jvalues and the jval_cal_* subroutines.

  INTEGER :: klev        ! number of levels
  INTEGER :: kproma_day  ! number of columns with daylight

  REAL, PARAMETER :: press_ref = 1.E4                        ! [Pa]
  REAL, PARAMETER :: dens_ref  = PRESS_REF/(k_B* 250.)*1.E-6 ! [molec./cm^3]

  ! coefficients for lyman-alpha parameterization from
  ! S. Chabrillat and G. Kockarts grl 24, p. 2659-2662 [2633]
  ! with the corrected values from grl 25 p.79  [2634]
  REAL, PARAMETER :: b_la(3) = (/ 6.84310E-01,  2.29821E-01,  8.65412E-02 /)
  REAL, PARAMETER :: c_la(3) = (/ 8.22114E-21,  1.77556E-20,  8.22112E-21 /)
  REAL, PARAMETER :: d_la(3) = (/ 6.00730E-21,  4.28569E-21,  1.28059E-20 /)
  REAL, PARAMETER :: e_la(3) = (/ 8.21666E-21,  1.63296E-20,  4.85121E-17 /)

  INTEGER, DIMENSION(:),   ALLOCATABLE  :: iu0 ! index for sorting kproma_day

  REAL,    DIMENSION(:,:),   ALLOCATABLE :: temp, press

  ! slant O2 column [m/cm^2]:
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: v2s_m

  ! slant O3 column [DU] (v3_du1 = intervals 0-2, v3_du2 = intervals 3-7):
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: v3_du1, v3_du2

  REAL,    DIMENSION(:,:),   ALLOCATABLE :: dlv2 ! log(v2)
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: dens ! density [molec/cm3]

  ! integrated actinic flux:
  REAL,    DIMENSION(:,:,:), ALLOCATABLE :: fint, finth

  ! correction for sza > sza_thr:
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: fj_corr

  REAL,    DIMENSION(:,:),   ALLOCATABLE :: tnorm_sr ! normalized T (T0 = 240 K)
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: tnorm    ! normalized T (T0 = 250 K)

  ! indices for lookup table:
  INTEGER, DIMENSION(:,:),   ALLOCATABLE :: i0, i1, i2, i3

  REAL,    DIMENSION(:,:),   ALLOCATABLE :: sig_top_o2, sig_top_o3
  REAL,    DIMENSION(:,:),   ALLOCATABLE :: r_m, r_o2

  REAL, DIMENSION(maxwav), PARAMETER :: fmax = &
       (/  2.5400E+12, 2.2200E+14, 7.4700E+13, 8.9700E+13, &
       1.2600E+14, 1.1100E+15, 2.6700E+15/)

  REAL, DIMENSION(maxwav), PARAMETER :: fmin = &
       (/  2.3400E+12, 2.2100E+14, 7.4300E+13, 8.9200E+13, &
       1.2600E+14, 1.1100E+15, 2.6700E+15/)

  REAL, DIMENSION(maxwav), PARAMETER :: f0max = &
       (/6.97E+01, 5.1333, 1.7116E+01, 8.930, 2.6699E+01, 18.62, 48.184/)
  REAL, DIMENSION(maxwav), PARAMETER :: f0min = &
       (/7.2679E+01,  5.0584E+00,  1.7139E+01,  8.9395E+00, &
       2.6714E+01,  1.8595E+01,  4.8184E+01/)

  !mz_ht_20151124+
  LOGICAL, PUBLIC, SAVE                        :: l_aero_inp=.FALSE.
  ! op_pj_20160825+
  LOGICAL, PUBLIC, SAVE                        :: l_export_to_smil = .FALSE.
  ! op_pj_20160825-
  REAL(dp), POINTER, DIMENSION(:,:,:), PUBLIC  :: zaer_sca => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:), PUBLIC  :: zaer_abs => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:), PUBLIC  :: zaer_ga  => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:), PUBLIC  :: zjv_asca => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:), PUBLIC  :: zjv_aabs => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:), PUBLIC  :: zjv_ga   => NULL()
  !mz_ht_20151124-

  ! **************************************************************************

CONTAINS

  ! **************************************************************************

  SUBROUTINE aerosol_data

    ! aerosol data from: Models for aerosols of the lower atmosphere
    ! and the effects of humidity variations on their optical properties
    ! E. P. Shettle and R. W. Fenn (1979)
    ! Environmental Research Paper No. 676
    ! different aerosol types:
    ! 1 = rural aerosol
    ! 2 = maritime aerosol
    ! 3 = urban aerosol
    ! 4 = free troposphere aerosol

    INTEGER :: i
    REAL, PARAMETER :: pn_ref(4) = (/ 15000., 4000., 20000., 5000./)

    ! rural model
    ! relative humidity = 0.0 %
    aext(:,1,1) = (/ 3.0331E-01, 2.6173E-01, 2.5400E-01, 2.5011E-01, 2.4395E-01, 2.1722E-01, 1.3731E-01 /)
    aabs(:,1,1) = (/ 9.2798E-02, 2.2387E-02, 1.6832E-02, 1.5121E-02, 1.3489E-02, 1.1261E-02, 8.4159E-03 /)
    ag(:,1,1)   = (/ 7.5290E-01, 6.8383E-01, 6.7781E-01, 6.7575E-01, 6.7346E-01, 6.6742E-01, 6.4481E-01 /)
    ! relative humidity = 50.0 %
    aext(:,2,1) = (/ 3.1458E-01, 2.7099E-01, 2.6296E-01, 2.5893E-01, 2.5255E-01, 2.2496E-01, 1.4230E-01 /)
    aabs(:,2,1) = (/ 9.4048E-02, 2.2410E-02, 1.6798E-02, 1.5083E-02, 1.3467E-02, 1.1366E-02, 8.4195E-03 /)
    ag(:,2,1)   = (/ 7.5485E-01, 6.8880E-01, 6.8292E-01, 6.8088E-01, 6.7856E-01, 6.7231E-01, 6.5029E-01 /)
    ! relative humidity = 70.0 %
    aext(:,3,1) = (/ 3.3914E-01, 2.9169E-01, 2.8296E-01, 2.7858E-01, 2.7167E-01, 2.4184E-01, 1.5346E-01 /)
    aabs(:,3,1) = (/ 9.7106E-02, 2.2834E-02, 1.7015E-02, 1.5237E-02, 1.3562E-02, 1.1400E-02, 8.5197E-03 /)
    ag(:,3,1)   = (/ 7.5859E-01, 6.9748E-01, 6.9220E-01, 6.9041E-01, 6.8844E-01, 6.8315E-01, 6.6063E-01 /)
    ! relative humidity = 80.0 %
    aext(:,4,1) = (/ 4.6086E-01, 3.9492E-01, 3.8325E-01, 3.7747E-01, 3.6840E-01, 3.2934E-01, 2.1149E-01 /)
    aabs(:,4,1) = (/ 1.0993E-01, 2.4472E-02, 1.7877E-02, 1.5895E-02, 1.4081E-02, 1.2031E-02, 8.8703E-03 /)
    ag(:,4,1)   = (/ 7.6932E-01, 7.2720E-01, 7.2359E-01, 7.2237E-01, 7.2103E-01, 7.1720E-01, 6.9661E-01 /)
    ! relative humidity = 90.0 %
    aext(:,5,1) = (/ 6.7903E-01, 6.1749E-01, 5.5327E-01, 5.1363E-01, 4.4529E-01, 2.0029E-01, 2.6598E-01 /)
    aabs(:,5,1) = (/ 1.2497E-01, 2.5991E-02, 1.8534E-02, 1.6356E-02, 1.4461E-02, 1.2863E-02, 9.1785E-03 /)
    ag(:,5,1)   = (/ 7.7428E-01, 7.5237E-01, 7.5024E-01, 7.4945E-01, 7.4847E-01, 7.4511E-01, 7.2836E-01 /)
    ! relative humidity = 95.0 %
    aext(:,6,1) = (/ 8.0780E-01, 7.0467E-01, 6.8690E-01, 6.7817E-01, 6.6449E-01, 6.0441E-01, 4.0403E-01 /)
    aabs(:,6,1) = (/ 1.3235E-01, 2.6829E-02, 1.8930E-02, 1.6640E-02, 1.4677E-02, 1.3195E-02, 9.3323E-03 /)
    ag(:,6,1)   = (/ 7.7486E-01, 7.6184E-01, 7.6042E-01, 7.5986E-01, 7.5910E-01, 7.5619E-01, 7.4127E-01 /)
    ! relative humidity = 98.0 %
    aext(:,7,1) = (/ 1.0511E+00, 9.3767E-01, 9.1720E-01, 9.0699E-01, 8.9083E-01, 8.1871E-01, 5.6955E-01 /)
    aabs(:,7,1) = (/ 1.4878E-01, 2.9403E-02, 2.0470E-02, 1.7882E-02, 1.5667E-02, 1.4020E-02, 9.8336E-03 /)
    ag(:,7,1)   = (/ 7.7704E-01, 7.7311E-01, 7.7252E-01, 7.7225E-01, 7.7183E-01, 7.6996E-01, 7.5918E-01 /)
    ! relative humidity = 99.0 %
    aext(:,8,1) = (/ 1.2932E+00, 1.1771E+00, 1.1548E+00, 1.1435E+00, 1.1253E+00, 1.0428E+00, 7.4704E-01 /)
    aabs(:,8,1) = (/ 1.6170E-01, 3.0955E-02, 2.1241E-02, 1.8451E-02, 1.6103E-02, 1.4596E-02, 1.0031E-02 /)
    ag(:,8,1)   = (/ 7.7794E-01, 7.7937E-01, 7.7928E-01, 7.7918E-01, 7.7898E-01, 7.7782E-01, 7.7005E-01 /)
    ! maritime model
    ! relative humidity = 0.0 %
    aext(:,1,2) = (/ 1.1714E-01, 1.0744E-01, 1.0552E-01, 1.0455E-01, 1.0299E-01, 9.6330E-02, 7.7600E-02 /)
    aabs(:,1,2) = (/ 2.3293E-02, 4.7427E-03, 3.3183E-03, 2.8928E-03, 2.5072E-03, 2.0902E-03, 1.3660E-03 /)
    ag(:,1,2)   = (/ 7.4802E-01, 7.0052E-01, 6.9643E-01, 6.9505E-01, 6.9353E-01, 6.8918E-01, 6.7501E-01 /)
    ! relative humidity = 50.0 %
    aext(:,2,2) = (/ 1.2518E-01, 1.1512E-01, 1.1335E-01, 1.1248E-01, 1.1111E-01, 1.0502E-01, 8.4545E-02 /)
    aabs(:,2,2) = (/ 2.3540E-02, 4.7178E-03, 3.2864E-03, 2.8636E-03, 2.4880E-03, 2.1214E-03, 1.3557E-03 /)
    ag(:,2,2)   = (/ 7.5592E-01, 7.1030E-01, 7.0599E-01, 7.0441E-01, 7.0253E-01, 6.9753E-01, 6.9120E-01 /)
    ! relative humidity = 70.0 %
    aext(:,3,2) = (/ 1.4988E-01, 1.3954E-01, 1.3763E-01, 1.3667E-01, 1.3514E-01, 1.2846E-01, 1.0760E-01 /)
    aabs(:,3,2) = (/ 2.4193E-02, 4.7563E-03, 3.2808E-03, 2.8459E-03, 2.4610E-03, 2.0977E-03, 1.3526E-03 /)
    ag(:,3,2)   = (/ 7.6810E-01, 7.3178E-01, 7.2840E-01, 7.2718E-01, 7.2576E-01, 7.2243E-01, 7.2140E-01 /)
    ! relative humidity = 80.0 %
    aext(:,4,2) = (/ 2.6744E-01, 2.5365E-01, 2.5127E-01, 2.5009E-01, 2.4827E-01, 2.4046E-01, 2.1678E-01 /)
    aabs(:,4,2) = (/ 2.6889E-02, 4.8768E-03, 3.2419E-03, 2.7725E-03, 2.3777E-03, 2.1241E-03, 1.3276E-03 /)
    ag(:,4,2)   = (/ 7.9436E-01, 7.7966E-01, 7.7799E-01, 7.7730E-01, 7.7636E-01, 7.7352E-01, 7.7178E-01 /)
    ! relative humidity = 90.0 %
    aext(:,5,2) = (/ 3.8330E-01, 3.6488E-01, 3.6140E-01, 3.5965E-01, 3.5686E-01, 3.4492E-01, 3.1084E-01 /)
    aabs(:,5,2) = (/ 2.9979E-02, 4.9768E-03, 3.1745E-03, 2.6765E-03, 2.2897E-03, 2.2327E-03, 1.3043E-03 /)
    ag(:,5,2)   = (/ 8.0012E-01, 7.9255E-01, 7.9193E-01, 7.9173E-01, 7.9151E-01, 7.9080E-01, 7.8588E-01 /)
    ! relative humidity = 95.0 %
    aext(:,6,2) = (/ 5.1615E-01, 4.9529E-01, 4.9145E-01, 4.8953E-01, 4.8648E-01, 4.7318E-01, 4.3302E-01 /)
    aabs(:,6,2) = (/ 3.1363E-02, 5.0292E-03, 3.1472E-03, 2.6332E-03, 2.2444E-03, 2.2567E-03, 1.3038E-03 /)
    ag(:,6,2)   = (/ 8.0759E-01, 8.0372E-01, 8.0348E-01, 8.0342E-01, 8.0340E-01, 8.0318E-01, 7.9781E-01 /)
    ! relative humidity = 98.0 %
    aext(:,7,2) = (/ 7.8803E-01, 7.6471E-01, 7.6086E-01, 7.5899E-01, 7.5608E-01, 7.4288E-01, 6.9431E-01 /)
    aabs(:,7,2) = (/ 3.2908E-02, 5.1181E-03, 3.1457E-03, 2.6121E-03, 2.2174E-03, 2.2908E-03, 1.3037E-03 /)
    ag(:,7,2)   = (/ 8.1493E-01, 8.1833E-01, 8.1815E-01, 8.1795E-01, 8.1751E-01, 8.1522E-01, 8.0982E-01 /)
    ! relative humidity = 99.0 %
    aext(:,8,2) = (/ 1.1277E+00, 1.1062E+00, 1.1025E+00, 1.1006E+00, 1.0976E+00, 1.0838E+00, 1.0284E+00 /)
    aabs(:,8,2) = (/ 3.4070E-02, 5.1842E-03, 3.1481E-03, 2.6024E-03, 2.2079E-03, 2.3432E-03, 1.3072E-03 /)
    ag(:,8,2)   = (/ 8.2128E-01, 8.2690E-01, 8.2699E-01, 8.2690E-01, 8.2662E-01, 8.2476E-01, 8.1909E-01 /)
    ! urban model
    ! relative humidity = 0.0 %
    aext(:,1,3) = (/ 3.3148E-01, 2.9477E-01, 2.8744E-01, 2.8369E-01, 2.7766E-01, 2.5096E-01, 1.6723E-01 /)
    aabs(:,1,3) = (/ 1.3662E-01, 1.0755E-01, 1.0369E-01, 1.0197E-01, 9.9459E-02, 8.9558E-02, 6.0876E-02 /)
    ag(:,1,3)   = (/ 7.7487E-01, 7.2345E-01, 7.1743E-01, 7.1491E-01, 7.1141E-01, 6.9828E-01, 6.5684E-01 /)
    ! relative humidity = 50.0 %
    aext(:,2,3) = (/ 3.4926E-01, 3.0928E-01, 3.0148E-01, 2.9751E-01, 2.9115E-01, 2.6305E-01, 1.7487E-01 /)
    aabs(:,2,3) = (/ 1.4041E-01, 1.0979E-01, 1.0577E-01, 1.0398E-01, 1.0139E-01, 9.1213E-02, 6.1876E-02 /)
    ag(:,2,3)   = (/ 7.7764E-01, 7.2861E-01, 7.2286E-01, 7.2045E-01, 7.1711E-01, 7.0448E-01, 6.6367E-01 /)
    ! relative humidity = 70.0 %
    aext(:,3,3) = (/ 4.6116E-01, 4.0134E-01, 3.9057E-01, 3.8521E-01, 3.7675E-01, 3.3989E-01, 2.2481E-01 /)
    aabs(:,3,3) = (/ 1.6134E-01, 1.2160E-01, 1.1665E-01, 1.1452E-01, 1.1148E-01, 9.9897E-02, 6.7307E-02 /)
    ag(:,3,3)   = (/ 7.8806E-01, 7.5152E-01, 7.4701E-01, 7.4507E-01, 7.4232E-01, 7.3156E-01, 6.9530E-01 /)
    ! relative humidity = 80.0 %
    aext(:,4,3) = (/ 7.0148E-01, 6.0491E-01, 5.8850E-01, 5.8046E-01, 5.6793E-01, 5.1364E-01, 3.4072E-01 /)
    aabs(:,4,3) = (/ 1.9374E-01, 1.3916E-01, 1.3278E-01, 1.3012E-01, 1.2643E-01, 1.1285E-01, 7.5244E-02 /)
    ag(:,4,3)   = (/ 7.9357E-01, 7.7371E-01, 7.7093E-01, 7.6965E-01, 7.6776E-01, 7.5984E-01, 7.3021E-01 /)
    ! relative humidity = 90.0 %
    aext(:,5,3) = (/ 1.0358E+00, 9.0438E-01, 8.8195E-01, 8.7096E-01, 8.5379E-01, 7.7850E-01, 5.2733E-01 /)
    aabs(:,5,3) = (/ 2.2461E-01, 1.5611E-01, 1.4843E-01, 1.4529E-01, 1.4105E-01, 1.2581E-01, 8.3559E-02 /)
    ag(:,5,3)   = (/ 7.9206E-01, 7.8587E-01, 7.8450E-01, 7.8377E-01, 7.8259E-01, 7.7704E-01, 7.5405E-01 /)
    ! relative humidity = 95.0 %
    aext(:,6,3) = (/ 1.4636E+00, 1.3097E+00, 1.2821E+00, 1.2683E+00, 1.2466E+00, 1.1491E+00, 8.0463E-01 /)
    aabs(:,6,3) = (/ 2.5345E-01, 1.7323E-01, 1.6439E-01, 1.6081E-01, 1.5601E-01, 1.3902E-01, 9.2399E-02 /)
    ag(:,6,3)   = (/ 7.8870E-01, 7.9251E-01, 7.9212E-01, 7.9176E-01, 7.9103E-01, 7.8710E-01, 7.7032E-01 /)
    ! relative humidity = 98.0 %
    aext(:,7,3) = (/ 2.2504E+00, 2.0973E+00, 2.0654E+00, 2.0488E+00, 2.0219E+00, 1.8952E+00, 1.4029E+00 /)
    aabs(:,7,3) = (/ 2.9543E-01, 1.9976E-01, 1.8929E-01, 1.8507E-01, 1.7944E-01, 1.5968E-01, 1.0601E-01 /)
    ag(:,7,3)   = (/ 7.8442E-01, 7.9688E-01, 7.9744E-01, 7.9746E-01, 7.9725E-01, 7.9526E-01, 7.8555E-01 /)
    ! relative humidity = 99.0 %
    aext(:,8,3) = (/ 2.9810E+00, 2.8574E+00, 2.8252E+00, 2.8076E+00, 2.7782E+00, 2.6340E+00, 2.0363E+00 /)
    aabs(:,8,3) = (/ 3.3013E-01, 2.2188E-01, 2.1001E-01, 2.0522E-01, 1.9882E-01, 1.7643E-01, 1.1639E-01 /)
    ag(:,8,3)   = (/ 7.8264E-01, 7.9839E-01, 7.9940E-01, 7.9963E-01, 7.9972E-01, 7.9897E-01, 7.9345E-01 /)
    ! free troposphere model
    ! relative humidity = 0.0 %
    aext(:,1,4) = (/ 9.6917E-02, 8.2981E-02, 8.0379E-02, 7.9069E-02, 7.6992E-02, 6.7999E-02, 4.1242E-02 /)
    aabs(:,1,4) = (/ 2.9018E-02, 5.9614E-03, 4.1862E-03, 3.6543E-03, 3.1696E-03, 2.6301E-03, 1.7271E-03 /)
    ag(:,1,4)   = (/ 7.4651E-01, 6.7635E-01, 6.7031E-01, 6.6827E-01, 6.6602E-01, 6.5979E-01, 6.3029E-01 /)
    ! relative humidity = 50.0 %
    aext(:,2,4) = (/ 1.0052E-01, 8.5881E-02, 8.3187E-02, 8.1836E-02, 7.9699E-02, 7.0451E-02, 4.2745E-02 /)
    aabs(:,2,4) = (/ 2.9357E-02, 5.9318E-03, 4.1461E-03, 3.6171E-03, 3.1448E-03, 2.6693E-03, 1.7141E-03 /)
    ag(:,2,4)   = (/ 7.4861E-01, 6.8186E-01, 6.7579E-01, 6.7365E-01, 6.7116E-01, 6.6398E-01, 6.3632E-01 /)
    ! relative humidity = 70.0 %
    aext(:,3,4) = (/ 1.0825E-01, 9.2398E-02, 8.9474E-02, 8.8007E-02, 8.5687E-02, 7.5677E-02, 4.6051E-02 /)
    aabs(:,3,4) = (/ 3.0195E-02, 5.9803E-03, 4.1383E-03, 3.5940E-03, 3.1103E-03, 2.6408E-03, 1.7100E-03 /)
    ag(:,3,4)   = (/ 7.5245E-01, 6.9062E-01, 6.8517E-01, 6.8330E-01, 6.8118E-01, 6.7508E-01, 6.4685E-01 /)
    ! relative humidity = 80.0 %
    aext(:,4,4) = (/ 1.4685E-01, 1.2477E-01, 1.2085E-01, 1.1890E-01, 1.1585E-01, 1.0271E-01, 6.3278E-02 /)
    aabs(:,4,4) = (/ 3.3598E-02, 6.1340E-03, 4.0904E-03, 3.5023E-03, 3.0054E-03, 2.6727E-03, 1.6797E-03 /)
    ag(:,4,4)   = (/ 7.6350E-01, 7.2093E-01, 7.1716E-01, 7.1587E-01, 7.1438E-01, 7.0974E-01, 6.8392E-01 /)
    ! relative humidity = 90.0 %
    aext(:,5,4) = (/ 2.1349E-01, 1.8287E-01, 1.7767E-01, 1.7512E-01, 1.7115E-01, 1.5382E-01, 9.7287E-02 /)
    aabs(:,5,4) = (/ 3.7513E-02, 6.2539E-03, 3.9994E-03, 3.3760E-03, 2.8912E-03, 2.8144E-03, 1.6493E-03 /)
    ag(:,5,4)   = (/ 7.6918E-01, 7.4714E-01, 7.4493E-01, 7.4409E-01, 7.4303E-01, 7.3917E-01, 7.1877E-01 /)
    ! relative humidity = 95.0 %
    aext(:,6,4) = (/ 2.5805E-01, 2.2351E-01, 2.1757E-01, 2.1465E-01, 2.1008E-01, 1.9000E-01, 1.2284E-01 /)
    aabs(:,6,4) = (/ 3.9242E-02, 6.3251E-03, 3.9701E-03, 3.3259E-03, 2.8370E-03, 2.8409E-03, 1.6490E-03 /)
    ag(:,6,4)   = (/ 7.6997E-01, 7.5695E-01, 7.5552E-01, 7.5494E-01, 7.5416E-01, 7.5099E-01, 7.3275E-01 /)
    ! relative humidity = 98.0 %
    aext(:,7,4) = (/ 3.2929E-01, 2.9127E-01, 2.8443E-01, 2.8102E-01, 2.7563E-01, 2.5151E-01, 1.6782E-01 /)
    aabs(:,7,4) = (/ 4.1205E-02, 6.4373E-03, 3.9673E-03, 3.2981E-03, 2.8015E-03, 2.8823E-03, 1.6499E-03 /)
    ag(:,7,4)   = (/ 7.7007E-01, 7.6643E-01, 7.6581E-01, 7.6552E-01, 7.6504E-01, 7.6274E-01, 7.4798E-01 /)
    ! relative humidity = 99.0 %
    aext(:,8,4) = (/ 4.0049E-01, 3.6140E-01, 3.5393E-01, 3.5014E-01, 3.4407E-01, 3.1647E-01, 2.1729E-01 /)
    aabs(:,8,4) = (/ 4.2656E-02, 6.5175E-03, 3.9683E-03, 3.2845E-03, 2.7889E-03, 2.9497E-03, 1.6547E-03 /)
    ag(:,8,4)   = (/ 7.6981E-01, 7.7196E-01, 7.7188E-01, 7.7176E-01, 7.7151E-01, 7.6990E-01, 7.5824E-01 /)

    DO i = 1,4
      asca(:,:,i) = aext(:,:,i) - aabs(:,:,i)   ! [1/km]
      ! scaling for particle density 1 part./cm^3 = effective cross sections
      asca(:,:,i) = asca(:,:,i)/pn_ref(i)*1.E-5
      aabs(:,:,i) = aabs(:,:,i)/pn_ref(i)*1.E-5
    ENDDO

  END SUBROUTINE aerosol_data

  ! **************************************************************************

  SUBROUTINE combine_o3_fields(input_O3, press_2d,  &
                               press_h, o3_h, v3_h, &
                               o3_2d, v3_2d)

    USE messy_main_constants_mem, ONLY: g, N_A, M_air

    INTRINSIC :: ABS, MAX, MIN, SIZE

    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: input_O3, press_2d
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: press_h, o3_h, v3_h
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: o3_2d
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: v3_2d

    ! for calculation of O3 column:
    INTEGER  :: kt0, kt1, j1_0, j2_0, j1_1, j2_1
    REAL(dp) :: av3_0, dpre_0, dpre_h_0, aro3_0
    REAL(dp) :: av3_1, dpre_1, dpre_h_1, aro3_1
    REAL(dp) :: dq1, dp1, dq2, dp2, press0
    INTEGER :: nlev_h, j, jp, jk
    REAL(dp), PARAMETER :: &
         sp = N_A * 1000./(M_air*g) * 1.E-4 ! [part./cm^2 * 1/Pa]

    nlev_h      = SIZE(press_h,2)
    o3_2d(:,3:) = input_O3(:,2:)

    ! calculate column density of O3

    ! O3 columns and slant columns
    ! use haloe O3 vertical columns above toa
    ! where is the press0 and press_2d(jp,1) in the haloe grid? ->
    ! kt0 next haloe level above press0
    ! kt1 next haloe level above press_2d(jp,1)
    jp_loop: DO jp = 1,SIZE(input_O3,1)
      press0 = 0.37 * press_2d(jp,1) ! 0.37 = 1/e
      dq1 = press_h(jp,1) - press0
      dp1 = press_h(jp,1) - press_2d(jp,1)
      kt0 = 1
      kt1 = 1
      DO j = 2,nlev_h
        dq2 = press_h(jp,j) - press0
        dp2 = press_h(jp,j) - press_2d(jp,1)
        IF (ABS(dq2) <= ABS(dq1)) THEN
          dq1 = dq2
          kt0 = j
        ENDIF
        IF (ABS(dp2) <= ABS(dp1)) THEN
          dp1 = dp2
          kt1 = j
        ENDIF
      ENDDO
      IF (dq1 <= 0.) kt0 = kt0 + 1
      IF (dp1 <= 0.) kt1 = kt1 + 1

      ! linear interpolation from haloe grid to press0 and press_2d(jp,1)
      ! level 0
      j1_0     = MAX(1,MIN(nlev_h-1,kt0))
      j2_0     = MAX(2,MIN(nlev_h,kt0+1))
      dpre_0   = press_h(jp,j2_0) - press0
      dpre_h_0 = press_h(jp,j2_0) - press_h(jp,j1_0)
      av3_0    = (v3_h(jp,j1_0) - v3_h(jp,j2_0)) / dpre_h_0
      aro3_0   = (o3_h(jp,j1_0) - o3_h(jp,j2_0)) / dpre_h_0
      ! level 1
      j1_1     = MIN(nlev_h-1,kt1)
      j2_1     = MIN(nlev_h,kt1+1)
      dpre_1   = press_h(jp,j2_1) - press_2d(jp,1)
      dpre_h_1 = press_h(jp,j2_1) - press_h(jp,j1_1)
      av3_1    = (v3_h(jp,j1_1) - v3_h(jp,j2_1)) / dpre_h_1
      aro3_1   = (o3_h(jp,j1_1) - o3_h(jp,j2_1)) / dpre_h_1
      ! o3_2d is a local variable with the same contents as relo3_2d
      ! (or the tracer O3) apart from 2 differences:
      ! - o3_2d also contains a dummy level 1 above top layer
      ! - level 2 of o3_2d is always from haloe data, not echam
      o3_2d(jp,1) = MAX(0._dp, aro3_0 * dpre_0 + o3_h(jp,j2_0))
      o3_2d(jp,2) = MAX(0._dp, aro3_1 * dpre_1 + o3_h(jp,j2_1))
      ! calculate v3_2d
      v3_2d(jp,1) = MAX(0._dp, av3_0 * dpre_0 + v3_h(jp,j2_0))
      v3_2d(jp,2) = MAX(0._dp, av3_1 * dpre_1 + v3_h(jp,j2_1))
      DO jk = 2,SIZE(press_2d,2)
        v3_2d(jp,jk+1) = v3_2d(jp,jk) + &
          sp * (press_2d(jp,jk)-press_2d(jp,jk-1)) &
          * 0.5 * (o3_2d(jp,jk+1) + o3_2d(jp,jk))
      ENDDO
      ! Note that v3_2d and o3_2d have nlev+1 layers with indices 1...nlev+1.
      ! Level 1 is above top layer; level 2 equals echam level 1 and so on

    ENDDO jp_loop

  END SUBROUTINE combine_o3_fields

  ! **************************************************************************

  SUBROUTINE jvalues(              &
    v3_2d,                         &
    cossza_1d, press_2d, relo3_2d, &
    rhum_2d, temp_2d, albedo_1d,   &
    aclc_2d, slf_1d, clp_2d,       &
    lmidatm, l_heating, pbllev)

    ! Note that although relo3_2d has the dimension 1:klev+1, the value
    ! relo3_2d(1) is not used at all here. Also, note that relo3_2d is
    ! _only_ used for the heating rates. For the calculation of the
    ! J-values, only v3_2d is used.

    ! interface for calculation of J-values and heating rates
    ! (a) the input data are sorted in j = 1,kproma_day, so that
    !     for j = 1,kproma_day u0(j) > u0lim (daytime).
    ! (b) corrections is done for u0 because of
    !     spherical geometry of the earth.
    ! (c) for some interpolation temperature and pressure
    !     on a virtual level "0" is necessary
    !     (see e.g. press(j,0) and temp(j,0) in aero_2d).
    ! (d) aerosols are taken in to account very rudimentary
    !     in the lowest pbllev levels a mixture of
    !     rural and maritime aerosol is asumed
    ! (e) initialzation of J-values rates
    ! (f) resort the output of the band model to j = 1,kproma

    ! correction of the air mass factor
    ! F. Kasten and T. Young, Revised optical air mass tabels and
    ! approximation formula (1989) Applied Optics Vol 28, no. 22 p. 4735
    ! and J. Lenoble, Atmospheric Radiative Transfer (1993), p. 236

    USE messy_main_constants_mem, ONLY: N_A, R_gas, pi, g, k_B, M_air

    INTRINSIC :: ACOS, AINT, COS, EXP, INT, LOG, MAX, MIN, SIZE, SQRT

    !-------------------------------------------------------------------------

    ! I/O

    REAL, DIMENSION(:,:), INTENT(IN) :: v3_2d         ! ozone column
    REAL, DIMENSION(:),   INTENT(IN) :: cossza_1d
    REAL, DIMENSION(:,:), INTENT(IN) :: press_2d
    REAL, DIMENSION(:,:), INTENT(IN) :: relo3_2d
    REAL, DIMENSION(:,:), INTENT(IN) :: rhum_2d
    REAL, DIMENSION(:,:), INTENT(IN) :: temp_2d
    REAL, DIMENSION(:),   INTENT(IN) :: albedo_1d
    REAL, DIMENSION(:,:), INTENT(IN) :: aclc_2d
    REAL, DIMENSION(:),   INTENT(IN) :: slf_1d
    REAL, DIMENSION(:,:), INTENT(IN) :: clp_2d
    LOGICAL,              INTENT(IN) :: lmidatm
    LOGICAL,              INTENT(IN) :: l_heating
    INTEGER,              INTENT(IN) :: pbllev   ! number of levels in pbl

    !-------------------------------------------------------------------------

    ! LOCAL PARAMETERs

    REAL, PARAMETER :: M_air_SI = M_air/1000.  ! molar mass of air [kg/mol]
    REAL, PARAMETER :: relo2 = 0.2095          ! O2 mixing ratio

    ! limiting values
    ! solar zenith angle threshold when correction becomes necessary
    REAL, PARAMETER :: sza_thr = 87.5
    ! if cos(sza) < u0lim => night
    ! changed to 94.5deg, important for polar lower stratosphere
    REAL, PARAMETER :: u0lim   = -0.078459095
    ! if cloud fraction < fraclim => clear sky
    REAL, PARAMETER :: fraclim = 0.01

    REAL, PARAMETER :: part_rural = 15000.
    REAL, PARAMETER :: part_sea   = 4000.

    ! O3 and O2 cross sections for T = 250 K
    REAL, PARAMETER :: crs_o3(7) = (/ &
      3.6448E-19,  1.8004E-18,  2.7700E-19, &
      1.0500E-19,  2.6000E-20,  0.0000E+00, 4.5500E-21/)
    REAL, PARAMETER :: crs_o2(7) = (/ &
      7.6300E-24,  0.0000E+00,  0.0000E+00, &
      0.0000E+00,  0.0000E+00,  0.0000E+00, 0.0000E+00/)

    ! polynomial coeff. to calculate tau_0 above 1 Pa
    REAL, PARAMETER :: ct(3) = (/ 1.3625E-01, 2.7173E+00, 1.7216E+01 /)

    ! change the units (part./cm^2 -> dobson units)
    ! Boltzmann constant k = 1.38E-23 [J/K], normal conditions T0 = 273 K,
    ! p0 = 1000 hPa, so constant = k*t0/p0 = 3.767E-20 cm^3
    REAL, PARAMETER :: constant = 3.767E-20

    ! look-up table
    INTEGER, PARAMETER :: ifil(115) = (/ &
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
      11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
      21,  22,  23,  24,  25,  26,  27,  28,  29,  30, &
      31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
      41,  41,  41,  41,  41,  42,  42,  42,  42,  42, &
      43,  43,  43,  43,  43,  44,  44,  44,  44,  44, &
      45,  45,  45,  45,  45,  46,  46,  46,  46,  46, &
      47,  47,  47,  47,  47,  48,  48,  48,  48,  48, &
      49,  49,  49,  49,  49,  50,  50,  50,  50,  50, &
      51,  51,  51,  51,  51,  52,  52,  52,  52,  52, &
      53,  53,  53,  53,  53,  54,  54,  54,  54,  54, &
      55,  55,  55,  55,  55 /)

    !---------------------------------------------------------------------------
    ! O2
    !---------------------------------------------------------------------------

    ! T-const parameters:
    REAL, PARAMETER :: a0_1_1_O2(dim58) = (/ &
      -8.2183E-22, -7.5708E-22, -6.9570E-22, -6.3768E-22, -5.8243E-22, -5.2913E-22, &
      -4.7808E-22, -4.2947E-22, -3.8362E-22, -3.4122E-22, -3.0251E-22, -2.6790E-22, &
      -2.3760E-22, -2.1142E-22, -1.8694E-22, -1.6567E-22, -1.4694E-22, -1.2963E-22, &
      -1.1419E-22, -1.0028E-22, -8.7287E-23, -7.5680E-23, -6.5376E-23, -5.5981E-23, &
      -4.7882E-23, -4.0766E-23, -3.4783E-23, -2.9760E-23, -2.5508E-23, -2.1894E-23, &
      -1.8788E-23, -1.6086E-23, -1.3715E-23, -1.1618E-23, -9.7779E-24, -8.1663E-24, &
      -6.7897E-24, -5.6263E-24, -4.6643E-24, -3.8835E-24, -3.2620E-24, -2.7741E-24, &
      -2.3911E-24, -2.0890E-24, -1.8410E-24, -1.6271E-24, -1.4312E-24, -1.2481E-24, &
      -1.0736E-24, -9.1030E-25, -7.6036E-25, -6.2472E-25, -5.0350E-25, -3.9726E-25, &
      -3.0975E-25, -2.5067E-25, -2.3447E-25, -2.5631E-25/)
    REAL, PARAMETER :: a0_1_2_O2(dim58) = (/ &
      3.7793E-20, 3.4931E-20, 3.2206E-20, 2.9618E-20, 2.7143E-20, 2.4745E-20, &
      2.2437E-20, 2.0230E-20, 1.8139E-20, 1.6198E-20, 1.4417E-20, 1.2818E-20, &
      1.1412E-20, 1.0192E-20, 9.0463E-21, 8.0468E-21, 7.1626E-21, 6.3423E-21, &
      5.6073E-21, 4.9425E-21, 4.3187E-21, 3.7592E-21, 3.2605E-21, 2.8039E-21, &
      2.4087E-21, 2.0600E-21, 1.7656E-21, 1.5175E-21, 1.3066E-21, 1.1266E-21, &
      9.7134E-22, 8.3568E-22, 7.1620E-22, 6.1008E-22, 5.1660E-22, 4.3441E-22, &
      3.6393E-22, 3.0413E-22, 2.5449E-22, 2.1405E-22, 1.8173E-22, 1.5626E-22, &
      1.3619E-22, 1.2030E-22, 1.0720E-22, 9.5872E-23, 8.5447E-23, 7.5670E-23, &
      6.6315E-23, 5.7532E-23, 4.9435E-23, 4.2083E-23, 3.5489E-23, 2.9688E-23, &
      2.4892E-23, 2.1643E-23, 2.0749E-23, 2.1959E-23/)
    REAL, PARAMETER :: a0_2_1_O2(dim58) = (/ &
      5.4893E-24, 4.9910E-24, 4.5286E-24, 4.1033E-24, 3.7097E-24, 3.3368E-24, &
      2.9851E-24, 2.6540E-24, 2.3434E-24, 2.0576E-24, 1.7978E-24, 1.5661E-24, &
      1.3778E-24, 1.2351E-24, 1.0531E-24, 9.0950E-25, 7.9314E-25, 6.7431E-25, &
      5.7891E-25, 4.9984E-25, 4.1944E-25, 3.5393E-25, 2.9938E-25, 2.4600E-25, &
      2.0324E-25, 1.6431E-25, 1.3426E-25, 1.1071E-25, 9.1312E-26, 7.5333E-26, &
      6.2129E-26, 5.1112E-26, 4.1850E-26, 3.4012E-26, 2.7413E-26, 2.1889E-26, &
      1.7368E-26, 1.3718E-26, 1.0836E-26, 8.5894E-27, 6.8820E-27, 5.5888E-27, &
      4.5734E-27, 3.7968E-27, 3.1593E-27, 2.6040E-27, 2.1105E-27, 1.6725E-27, &
      1.2896E-27, 9.7613E-28, 7.2546E-28, 5.4032E-28, 4.0244E-28, 2.9906E-28, &
      2.2089E-28, 1.6328E-28, 1.1931E-28, 1.0172E-28/)
    REAL, PARAMETER :: a0_2_2_O2(dim58) = (/ &
      -2.5135E-22, -2.2933E-22, -2.0879E-22, -1.8982E-22, -1.7220E-22, -1.5541E-22, &
      -1.3951E-22, -1.2448E-22, -1.1032E-22, -9.7233E-23, -8.5282E-23, -7.4577E-23, &
      -6.5840E-23, -5.9190E-23, -5.0670E-23, -4.3922E-23, -3.8430E-23, -3.2797E-23, &
      -2.8256E-23, -2.4477E-23, -2.0618E-23, -1.7460E-23, -1.4820E-23, -1.2225E-23, &
      -1.0139E-23, -8.2311E-24, -6.7526E-24, -5.5892E-24, -4.6273E-24, -3.8315E-24, &
      -3.1713E-24, -2.6182E-24, -2.1514E-24, -1.7549E-24, -1.4196E-24, -1.1379E-24, &
      -9.0640E-25, -7.1882E-25, -5.7007E-25, -4.5372E-25, -3.6493E-25, -2.9743E-25, &
      -2.4422E-25, -2.0337E-25, -1.6971E-25, -1.4028E-25, -1.1402E-25, -9.0638E-26, &
      -7.0114E-26, -5.3249E-26, -3.9713E-26, -2.9679E-26, -2.2178E-26, -1.6533E-26, &
      -1.2250E-26, -9.0811E-27, -6.6538E-27, -5.6796E-27/)
    REAL, PARAMETER :: a0_3_1_O2(dim58) = (/ &
      -1.1330E-26, -1.0224E-26, -9.2100E-27, -8.2949E-27, -7.4669E-27, -6.6897E-27, &
      -5.9624E-27, -5.2828E-27, -4.6454E-27, -4.0602E-27, -3.5298E-27, -3.0568E-27, &
      -2.7412E-27, -2.5965E-27, -2.1167E-27, -1.7930E-27, -1.5636E-27, -1.2740E-27, &
      -1.0735E-27, -9.2508E-28, -7.4518E-28, -6.1609E-28, -5.1750E-28, -4.0736E-28, &
      -3.2786E-28, -2.5062E-28, -1.9757E-28, -1.5981E-28, -1.2914E-28, -1.0429E-28, &
      -8.4200E-29, -6.7867E-29, -5.4443E-29, -4.3384E-29, -3.4274E-29, -2.6808E-29, &
      -2.0823E-29, -1.6120E-29, -1.2487E-29, -9.7168E-30, -7.6841E-30, -6.1691E-30, &
      -4.9535E-30, -4.0688E-30, -3.3467E-30, -2.7017E-30, -2.1353E-30, -1.6398E-30, &
      -1.2109E-30, -8.8107E-31, -6.1509E-31, -4.3973E-31, -3.1073E-31, -2.1638E-31, &
      -1.4626E-31, -1.0029E-31, -6.1246E-32, -5.2422E-32/)
    REAL, PARAMETER :: a0_3_2_O2(dim58) = (/ &
      5.1805E-25, 4.6918E-25, 4.2415E-25, 3.8334E-25, 3.4624E-25, 3.1127E-25, &
      2.7840E-25, 2.4754E-25, 2.1848E-25, 1.9168E-25, 1.6728E-25, 1.4543E-25, &
      1.3078E-25, 1.2404E-25, 1.0159E-25, 8.6370E-26, 7.5543E-26, 6.1814E-26, &
      5.2270E-26, 4.5176E-26, 3.6541E-26, 3.0319E-26, 2.5547E-26, 2.0194E-26, &
      1.6315E-26, 1.2530E-26, 9.9199E-27, 8.0548E-27, 6.5334E-27, 5.2960E-27, &
      4.2914E-27, 3.4715E-27, 2.7949E-27, 2.2353E-27, 1.7725E-27, 1.3918E-27, &
      1.0853E-27, 8.4357E-28, 6.5613E-28, 5.1263E-28, 4.0693E-28, 3.2785E-28, &
      2.6415E-28, 2.1761E-28, 1.7949E-28, 1.4530E-28, 1.1517E-28, 8.8709E-29, &
      6.5721E-29, 4.7975E-29, 3.3613E-29, 2.4108E-29, 1.7090E-29, 1.1939E-29, &
      8.0961E-30, 5.5681E-30, 3.4127E-30, 2.9238E-30/)
    REAL, PARAMETER :: c0_1_1_O2(dim58) = (/ &
      1.6313E-01, 1.6313E-01, 1.6313E-01, 1.5033E-01, 1.2271E-01, 8.0699E-02, &
      2.7953E-02, -3.1080E-02, -8.9210E-02, -1.3973E-01, -1.7625E-01, -1.9421E-01, &
      -1.8968E-01, -1.6239E-01, -1.1548E-01, -5.5354E-02, 1.1115E-02, 7.2937E-02, &
      1.2309E-01, 1.5084E-01, 1.5209E-01, 1.2658E-01, 7.7491E-02, 1.3151E-02, &
      -5.3841E-02, -1.0972E-01, -1.4160E-01, -1.4005E-01, -1.0220E-01, -3.1586E-02, &
      6.0093E-02, 1.5487E-01, 2.3362E-01, 2.7658E-01, 2.7465E-01, 2.2069E-01, &
      1.2617E-01, 2.9936E-03, -1.2934E-01, -2.5313E-01, -3.5165E-01, -4.1483E-01, &
      -4.4026E-01, -4.2883E-01, -3.8954E-01, -3.3760E-01, -2.8310E-01, -2.3620E-01, &
      -1.9962E-01, -1.7123E-01, -1.4599E-01, -1.1691E-01, -7.8799E-02, -3.1119E-02, &
      1.5819E-02, 4.0198E-02, 1.2886E-02, -5.1621E-02 /)
    REAL, PARAMETER :: c0_1_2_O2(dim58) = (/ &
      -6.6121E+00, -6.6121E+00, -6.6121E+00, -6.0413E+00, -4.8040E+00, -2.9134E+00, &
      -5.2923E-01, 2.1509E+00, 4.8016E+00, 7.1154E+00, 8.7953E+00, 9.6253E+00, &
      9.4146E+00, 8.1431E+00, 5.9477E+00, 3.1218E+00, -1.5497E-02, -2.9459E+00, &
      -5.3333E+00, -6.6594E+00, -6.7195E+00, -5.4899E+00, -3.1141E+00, 1.2847E-02, &
      3.2821E+00, 6.0202E+00, 7.5888E+00, 7.5119E+00, 5.6346E+00, 2.1181E+00, &
      -2.4658E+00, -7.2238E+00, -1.1192E+01, -1.3366E+01, -1.3268E+01, -1.0516E+01, &
      -5.6771E+00, 6.5426E-01, 7.4827E+00, 1.3895E+01, 1.9018E+01, 2.2316E+01, &
      2.3648E+01, 2.3047E+01, 2.0973E+01, 1.8220E+01, 1.5321E+01, 1.2816E+01, &
      1.0856E+01, 9.3281E+00, 7.9653E+00, 6.3888E+00, 4.3157E+00, 1.7124E+00, &
      -8.5978E-01, -2.2006E+00, -6.9303E-01, 2.8807E+00 /)
    REAL, PARAMETER :: c0_2_1_O2(dim58) = (/ &
      1.0067E-01, 1.0067E-01, 1.0067E-01, 8.8865E-02, 6.1607E-02, 1.8703E-02, &
      -2.9928E-02, -8.2111E-02, -1.2935E-01, -1.6186E-01, -1.8279E-01, -1.8464E-01, &
      -1.7619E-01, -1.4626E-01, -1.0857E-01, -6.3256E-02, -1.8559E-02, 2.6886E-02, &
      7.1887E-02, 1.1214E-01, 1.1752E-01, 1.1554E-01, 8.7678E-02, 5.0931E-02, &
      1.8629E-02, 4.8977E-03, 1.7167E-02, 4.9338E-02, 8.7722E-02, 1.1847E-01, &
      1.2938E-01, 1.2377E-01, 1.0500E-01, 8.3863E-02, 3.8984E-02, -2.2387E-02, &
      -9.2748E-02, -1.8991E-01, -2.5700E-01, -2.9953E-01, -2.7072E-01, -2.0761E-01, &
      -1.3340E-01, -5.4193E-02, -1.7096E-02, -1.2603E-02, -1.9185E-02, -1.6969E-02, &
      2.0296E-03, 3.2798E-02, 7.0451E-02, 9.8498E-02, 9.8113E-02, 6.3805E-02, &
      3.0444E-03, -5.1674E-02, -8.3758E-02, -1.0412E-01 /)
    REAL, PARAMETER :: c0_2_2_O2(dim58) = (/ &
      -4.2499E+00, -4.2499E+00, -4.2499E+00, -3.7233E+00, -2.5021E+00, -5.7147E-01, &
      1.6266E+00, 3.9958E+00, 6.1499E+00, 7.6388E+00, 8.6017E+00, 8.6869E+00, &
      8.2949E+00, 6.9003E+00, 5.1363E+00, 3.0066E+00, 8.9690E-01, -1.2572E+00, &
      -3.3992E+00, -5.3234E+00, -5.5815E+00, -5.4860E+00, -4.1376E+00, -2.3517E+00, &
      -7.7538E-01, -1.0253E-01, -7.0618E-01, -2.2954E+00, -4.1993E+00, -5.7305E+00, &
      -6.2762E+00, -5.9946E+00, -5.0484E+00, -3.9789E+00, -1.6990E+00, 1.4309E+00, &
      5.0334E+00, 1.0028E+01, 1.3489E+01, 1.5692E+01, 1.4195E+01, 1.0900E+01, &
      7.0114E+00, 2.8451E+00, 8.8633E-01, 6.4822E-01, 9.9838E-01, 8.8006E-01, &
      -1.3828E-01, -1.7936E+00, -3.8269E+00, -5.3470E+00, -5.3261E+00, -3.4529E+00, &
      -1.2319E-01, 2.8863E+00, 4.6573E+00, 5.7852E+00 /)

    !---------------------------------------------------------------------------
    ! O3
    !---------------------------------------------------------------------------

    ! interval 0
    REAL, PARAMETER :: a0_1_1_o3(dim58) = (/ &
      -5.0266E-21, -5.4393E-21, -5.8779E-21, -6.3301E-21, -6.7861E-21, -7.2489E-21, &
      -7.7158E-21, -8.1908E-21, -8.6591E-21, -9.1470E-21, -9.6073E-21, -1.0100E-20, &
      -1.1045E-20, -1.2851E-20, -1.2777E-20, -1.3153E-20, -1.3900E-20, -1.3768E-20, &
      -1.4009E-20, -1.4565E-20, -1.4284E-20, -1.4350E-20, -1.4658E-20, -1.4278E-20, &
      -1.4194E-20, -1.3586E-20, -1.3282E-20, -1.3232E-20, -1.3161E-20, -1.3038E-20, &
      -1.2881E-20, -1.2684E-20, -1.2426E-20, -1.2129E-20, -1.1816E-20, -1.1479E-20, &
      -1.1166E-20, -1.0886E-20, -1.0644E-20, -1.0432E-20, -1.0255E-20, -1.0076E-20, &
      -9.8596E-21, -9.5929E-21, -9.2271E-21, -8.7438E-21, -8.1461E-21, -7.4984E-21, &
      -6.7978E-21, -6.1194E-21, -5.4602E-21, -4.8508E-21, -4.2780E-21, -3.7223E-21, &
      -3.1841E-21, -2.6694E-21, -2.1642E-21, -1.7583E-21 /)
    REAL, PARAMETER :: a0_1_2_o3(dim58) = (/ &
      6.5642E-19, 6.7466E-19, 6.9413E-19, 7.1430E-19, 7.3473E-19, 7.5556E-19, &
      7.7666E-19, 7.9823E-19, 8.1958E-19, 8.4192E-19, 8.6310E-19, 8.8585E-19, &
      9.2970E-19, 1.0139E-18, 1.0104E-18, 1.0281E-18, 1.0633E-18, 1.0571E-18, &
      1.0686E-18, 1.0951E-18, 1.0817E-18, 1.0848E-18, 1.0997E-18, 1.0812E-18, &
      1.0772E-18, 1.0474E-18, 1.0324E-18, 1.0299E-18, 1.0264E-18, 1.0203E-18, &
      1.0124E-18, 1.0025E-18, 9.8956E-19, 9.7450E-19, 9.5861E-19, 9.4142E-19, &
      9.2542E-19, 9.1102E-19, 8.9855E-19, 8.8755E-19, 8.7836E-19, 8.6901E-19, &
      8.5767E-19, 8.4363E-19, 8.2432E-19, 7.9871E-19, 7.6691E-19, 7.3232E-19, &
      6.9477E-19, 6.5827E-19, 6.2267E-19, 5.8965E-19, 5.5848E-19, 5.2814E-19, &
      4.9865E-19, 4.7034E-19, 4.4246E-19, 4.1997E-19 /)
    REAL, PARAMETER :: a0_2_1_o3(dim58) = (/ &
      1.9659E-23, 2.1210E-23, 2.2803E-23, 2.4516E-23, 2.6124E-23, 2.7802E-23, &
      2.9458E-23, 3.1153E-23, 3.2770E-23, 3.4473E-23, 3.5961E-23, 3.7649E-23, &
      4.6093E-23, 6.6596E-23, 5.9393E-23, 5.8976E-23, 6.3617E-23, 5.7209E-23, &
      5.6189E-23, 5.9276E-23, 5.3213E-23, 5.1644E-23, 5.3096E-23, 4.7637E-23, &
      4.5648E-23, 3.8769E-23, 3.5279E-23, 3.4507E-23, 3.3714E-23, 3.2727E-23, &
      3.1749E-23, 3.0697E-23, 2.9494E-23, 2.8250E-23, 2.7014E-23, 2.5721E-23, &
      2.4544E-23, 2.3512E-23, 2.2609E-23, 2.1747E-23, 2.1044E-23, 2.0311E-23, &
      1.9467E-23, 1.8548E-23, 1.7366E-23, 1.5952E-23, 1.4283E-23, 1.2631E-23, &
      1.0890E-23, 9.2987E-24, 7.8249E-24, 6.5515E-24, 5.4049E-24, 4.3587E-24, &
      3.4177E-24, 2.6250E-24, 1.9674E-24, 1.5040E-24 /)
    REAL, PARAMETER :: a0_2_2_o3(dim58) = (/ &
      -1.2080E-21, -1.2765E-21, -1.3473E-21, -1.4236E-21, -1.4957E-21, -1.5712E-21, &
      -1.6461E-21, -1.7230E-21, -1.7967E-21, -1.8747E-21, -1.9432E-21, -2.0212E-21, &
      -2.4130E-21, -3.3684E-21, -3.0313E-21, -3.0117E-21, -3.2308E-21, -2.9270E-21, &
      -2.8785E-21, -3.0260E-21, -2.7350E-21, -2.6594E-21, -2.7297E-21, -2.4643E-21, &
      -2.3673E-21, -2.0302E-21, -1.8585E-21, -1.8204E-21, -1.7811E-21, -1.7319E-21, &
      -1.6830E-21, -1.6302E-21, -1.5696E-21, -1.5066E-21, -1.4438E-21, -1.3779E-21, &
      -1.3176E-21, -1.2646E-21, -1.2179E-21, -1.1733E-21, -1.1368E-21, -1.0985E-21, &
      -1.0543E-21, -1.0059E-21, -9.4355E-22, -8.6860E-22, -7.7982E-22, -6.9157E-22, &
      -5.9824E-22, -5.1265E-22, -4.3307E-22, -3.6405E-22, -3.0167E-22, -2.4455E-22, &
      -1.9298E-22, -1.4939E-22, -1.1308E-22, -8.7411E-23 /)

    REAL, PARAMETER :: c0_1_1_o3(dim58) = (/ &
      -7.0191E-03, -7.0191E-03, -7.0191E-03, -8.1775E-03, -9.0525E-03, -1.0315E-02, &
      -1.1296E-02, -1.1692E-02, -1.2369E-02, -1.2926E-02, -1.2667E-02, -1.3050E-02, &
      -1.2699E-02, -1.2539E-02, -1.2559E-02, -1.2045E-02, -1.2295E-02, -1.2351E-02, &
      -1.2501E-02, -1.2481E-02, -1.2809E-02, -1.2848E-02, -1.2163E-02, -1.1147E-02, &
      -1.0137E-02, -7.6925E-03, -5.2963E-03, -2.6730E-03, -2.2028E-04, 2.1178E-03, &
      3.7890E-03, 4.4592E-03, 3.6196E-03, 1.6453E-03, -1.3751E-03, -4.4430E-03, &
      -8.0093E-03, -1.0779E-02, -1.2009E-02, -1.1375E-02, -8.2903E-03, -3.7170E-03, &
      3.1828E-03, 9.6868E-03, 1.6113E-02, 2.0800E-02, 2.2935E-02, 2.2651E-02, &
      2.1574E-02, 2.0290E-02, 2.0017E-02, 2.1864E-02, 2.5309E-02, 2.9336E-02, &
      3.2886E-02, 3.4905E-02, 3.6604E-02, 4.0160E-02 /)
    REAL, PARAMETER :: c0_1_2_o3(dim58) = (/ &
      3.1254E-01, 3.1254E-01, 3.1254E-01, 3.6421E-01, 4.0341E-01, 4.6024E-01, &
      5.0458E-01, 5.2252E-01, 5.5342E-01, 5.7890E-01, 5.6701E-01, 5.8472E-01, &
      5.6840E-01, 5.6096E-01, 5.6190E-01, 5.3773E-01, 5.4953E-01, 5.5220E-01, &
      5.5935E-01, 5.5836E-01, 5.7414E-01, 5.7598E-01, 5.4284E-01, 4.9348E-01, &
      4.4419E-01, 3.2440E-01, 2.0651E-01, 7.6914E-02, -4.4740E-02, -1.6118E-01, &
      -2.4473E-01, -2.7838E-01, -2.3607E-01, -1.3616E-01, 1.7271E-02, 1.7373E-01, &
      3.5633E-01, 4.9868E-01, 5.6214E-01, 5.2933E-01, 3.6891E-01, 1.3019E-01, &
      -2.3136E-01, -5.7347E-01, -9.1277E-01, -1.1612E+00, -1.2748E+00, -1.2596E+00, &
      -1.2019E+00, -1.1328E+00, -1.1181E+00, -1.2182E+00, -1.4056E+00, -1.6255E+00, &
      -1.8200E+00, -1.9310E+00, -2.0248E+00, -2.2218E+00 /)
    REAL, PARAMETER :: c0_2_1_o3(dim58) = (/ &
      8.0025E-04, 8.0025E-04, 8.0025E-04, -4.2989E-03, -1.7255E-03, -1.4436E-03, &
      -2.4649E-03, -2.5963E-03, -2.9353E-03, -2.8235E-03, -5.4407E-05, -3.7774E-04, &
      7.6838E-05, 4.3952E-03, 3.2172E-03, 3.2044E-03, 5.0126E-03, 5.1538E-03, &
      5.0230E-03, 5.7908E-03, 7.5288E-03, 5.7935E-03, 6.8681E-03, 5.7070E-03, &
      5.8363E-03, 6.7479E-03, 3.9212E-03, 4.5167E-03, 1.3028E-03, 1.6390E-03, &
      6.7867E-04, -1.5781E-03, -1.9784E-03, -1.4371E-03, -1.6170E-03, -2.1123E-03, &
      -3.3934E-04, 1.9431E-03, 5.8769E-03, 1.2693E-02, 1.5833E-02, 2.1979E-02, &
      2.2427E-02, 1.9646E-02, 1.3546E-02, 1.1175E-02, 5.9491E-03, 2.3821E-03, &
      4.9171E-03, 8.4691E-03, 4.6355E-03, -6.7929E-03, -1.3726E-02, -3.1842E-02, &
      -4.3998E-02, -3.8830E-02, -2.5365E-02, -2.4250E-03 /)
    REAL, PARAMETER :: c0_2_2_o3(dim58) = (/ &
      -5.4338E-03, -5.4338E-03, -5.4338E-03, 2.2199E-01, 1.0670E-01, 9.4016E-02, &
      1.4018E-01, 1.4614E-01, 1.6160E-01, 1.5648E-01, 2.9102E-02, 4.4040E-02, &
      2.2948E-02, -1.7829E-01, -1.2316E-01, -1.2256E-01, -2.0790E-01, -2.1460E-01, &
      -2.0837E-01, -2.4507E-01, -3.2849E-01, -2.4485E-01, -2.9687E-01, -2.4043E-01, &
      -2.4675E-01, -2.9141E-01, -1.5234E-01, -1.8176E-01, -2.2346E-02, -3.9092E-02, &
      8.9268E-03, 1.2222E-01, 1.4239E-01, 1.1500E-01, 1.2414E-01, 1.4940E-01, &
      5.8625E-02, -5.8691E-02, -2.6168E-01, -6.1477E-01, -7.7800E-01, -1.0989E+00, &
      -1.1223E+00, -9.7603E-01, -6.5396E-01, -5.2830E-01, -2.5028E-01, -5.9807E-02, &
      -1.9569E-01, -3.8678E-01, -1.7977E-01, 4.3966E-01, 8.1682E-01, 1.8059E+00, &
      2.4721E+00, 2.1879E+00, 1.4446E+00, 1.7371E-01 /)

    !---------------------------------------------------------------------------
    ! heating rate O2
    !---------------------------------------------------------------------------

    ! interval 0
    REAL, PARAMETER :: a0_1_1_h_O2(dim58) = (/ &
      -1.8499E-17, -1.7023E-17, -1.5631E-17, -1.4320E-17, -1.3069E-17, -1.1865E-17, &
      -1.0715E-17, -9.6169E-18, -8.5858E-18, -7.6296E-18, -6.7583E-18, -5.9790E-18, &
      -5.2968E-18, -4.7079E-18, -4.1576E-18, -3.6800E-18, -3.2602E-18, -2.8722E-18, &
      -2.5272E-18, -2.2168E-18, -1.9271E-18, -1.6687E-18, -1.4396E-18, -1.2309E-18, &
      -1.0510E-18, -8.9307E-19, -7.6057E-19, -6.4948E-19, -5.5557E-19, -4.7594E-19, &
      -4.0769E-19, -3.4848E-19, -2.9663E-19, -2.5090E-19, -2.1085E-19, -1.7589E-19, &
      -1.4605E-19, -1.2092E-19, -1.0017E-19, -8.3372E-20, -6.9982E-20, -5.9493E-20, &
      -5.1295E-20, -4.4804E-20, -3.9449E-20, -3.4854E-20, -3.0635E-20, -2.6684E-20, &
      -2.2938E-20, -1.9461E-20, -1.6248E-20, -1.3355E-20, -1.0790E-20, -8.5388E-21, &
      -6.6786E-21, -5.4138E-21, -5.0405E-21, -5.4427E-21 /)
    REAL, PARAMETER :: a0_1_2_h_O2(dim58) = (/ &
      8.5044E-16, 7.8519E-16, 7.2342E-16, 6.6493E-16, 6.0888E-16, 5.5472E-16, &
      5.0272E-16, 4.5288E-16, 4.0586E-16, 3.6207E-16, 3.2199E-16, 2.8598E-16, &
      2.5433E-16, 2.2688E-16, 2.0113E-16, 1.7868E-16, 1.5887E-16, 1.4048E-16, &
      1.2405E-16, 1.0922E-16, 9.5312E-17, 8.2858E-17, 7.1768E-17, 6.1624E-17, &
      5.2848E-17, 4.5108E-17, 3.8590E-17, 3.3101E-17, 2.8444E-17, 2.4478E-17, &
      2.1066E-17, 1.8093E-17, 1.5480E-17, 1.3166E-17, 1.1132E-17, 9.3487E-18, &
      7.8209E-18, 6.5290E-18, 5.4583E-18, 4.5882E-18, 3.8920E-18, 3.3444E-18, &
      2.9148E-18, 2.5734E-18, 2.2906E-18, 2.0471E-18, 1.8227E-18, 1.6117E-18, &
      1.4109E-18, 1.2238E-18, 1.0503E-18, 8.9357E-19, 7.5402E-19, 6.3109E-19, &
      5.2915E-19, 4.5959E-19, 4.3898E-19, 4.6127E-19 /)
    REAL, PARAMETER :: a0_2_1_h_O2(dim58) = (/ &
      1.2394E-19, 1.1257E-19, 1.0210E-19, 9.2549E-20, 8.3622E-20, 7.5186E-20, &
      6.7265E-20, 5.9764E-20, 5.2789E-20, 4.6322E-20, 4.0467E-20, 3.5241E-20, &
      3.0990E-20, 2.7780E-20, 2.3672E-20, 2.0437E-20, 1.7819E-20, 1.5137E-20, &
      1.2990E-20, 1.1214E-20, 9.4013E-21, 7.9279E-21, 6.7009E-21, 5.4998E-21, &
      4.5377E-21, 3.6605E-21, 2.9863E-21, 2.4579E-21, 2.0230E-21, 1.6656E-21, &
      1.3712E-21, 1.1262E-21, 9.2054E-22, 7.4715E-22, 6.0127E-22, 4.7977E-22, &
      3.8022E-22, 3.0030E-22, 2.3703E-22, 1.8811E-22, 1.5045E-22, 1.2204E-22, &
      1.0039E-22, 8.3465E-23, 6.9243E-23, 5.7108E-23, 4.6401E-23, 3.6710E-23, &
      2.8284E-23, 2.1692E-23, 1.6032E-23, 1.1940E-23, 8.9192E-24, 6.7102E-24, &
      4.9096E-24, 3.6569E-24, 2.6802E-24, 2.2165E-24 /)
    REAL, PARAMETER :: a0_2_2_h_O2(dim58) = (/ &
      -5.6743E-18, -5.1721E-18, -4.7073E-18, -4.2811E-18, -3.8812E-18, -3.5016E-18, &
      -3.1436E-18, -2.8030E-18, -2.4849E-18, -2.1887E-18, -1.9195E-18, -1.6780E-18, &
      -1.4807E-18, -1.3312E-18, -1.1389E-18, -9.8684E-19, -8.6327E-19, -7.3616E-19, &
      -6.3397E-19, -5.4908E-19, -4.6206E-19, -3.9104E-19, -3.3166E-19, -2.7328E-19, &
      -2.2633E-19, -1.8335E-19, -1.5018E-19, -1.2407E-19, -1.0250E-19, -8.4706E-20, &
      -6.9986E-20, -5.7685E-20, -4.7322E-20, -3.8548E-20, -3.1137E-20, -2.4941E-20, &
      -1.9844E-20, -1.5736E-20, -1.2471E-20, -9.9369E-21, -7.9788E-21, -6.4958E-21, &
      -5.3613E-21, -4.4712E-21, -3.7202E-21, -3.0770E-21, -2.5074E-21, -1.9899E-21, &
      -1.5383E-21, -1.1837E-21, -8.7802E-22, -6.5621E-22, -4.9190E-22, -3.7129E-22, &
      -2.7262E-22, -2.0372E-22, -1.4980E-22, -1.2411E-22 /)
    REAL, PARAMETER :: a0_3_1_h_O2(dim58) = (/ &
      -2.5636E-22, -2.3107E-22, -2.0811E-22, -1.8770E-22, -1.6885E-22, -1.5120E-22, &
      -1.3484E-22, -1.1938E-22, -1.0509E-22, -9.1793E-23, -7.9820E-23, -6.9138E-23, &
      -6.1986E-23, -5.8746E-23, -4.7865E-23, -4.0550E-23, -3.5361E-23, -2.8790E-23, &
      -2.4253E-23, -2.0906E-23, -1.6821E-23, -1.3902E-23, -1.1669E-23, -9.1764E-24, &
      -7.3772E-24, -5.6257E-24, -4.4298E-24, -3.5769E-24, -2.8841E-24, -2.3240E-24, &
      -1.8733E-24, -1.5069E-24, -1.2068E-24, -9.6039E-25, -7.5676E-25, -5.9189E-25, &
      -4.5907E-25, -3.5550E-25, -2.7491E-25, -2.1442E-25, -1.6854E-25, -1.3466E-25, &
      -1.0970E-25, -9.0267E-26, -7.3631E-26, -5.9278E-26, -4.7190E-26, -3.6173E-26, &
      -2.6421E-26, -1.9961E-26, -1.3627E-26, -9.7720E-27, -6.7900E-27, -4.9931E-27, &
      -3.2057E-27, -2.3090E-27, -1.3269E-27, -1.1296E-27 /)
    REAL, PARAMETER :: a0_3_2_h_O2(dim58) = (/ &
      1.1722E-20, 1.0604E-20, 9.5845E-21, 8.6741E-21, 7.8296E-21, 7.0352E-21, &
      6.2960E-21, 5.5939E-21, 4.9425E-21, 4.3333E-21, 3.7826E-21, 3.2891E-21, &
      2.9572E-21, 2.8062E-21, 2.2970E-21, 1.9532E-21, 1.7082E-21, 1.3968E-21, &
      1.1808E-21, 1.0208E-21, 8.2475E-22, 6.8406E-22, 5.7602E-22, 4.5486E-22, &
      3.6705E-22, 2.8123E-22, 2.2239E-22, 1.8026E-22, 1.4589E-22, 1.1800E-22, &
      9.5465E-23, 7.7073E-23, 6.1949E-23, 4.9481E-23, 3.9136E-23, 3.0727E-23, &
      2.3927E-23, 1.8604E-23, 1.4445E-23, 1.1312E-23, 8.9260E-24, 7.1575E-24, &
      5.8498E-24, 4.8275E-24, 3.9492E-24, 3.1884E-24, 2.5454E-24, 1.9571E-24, &
      1.4343E-24, 1.0868E-24, 7.4475E-25, 5.3582E-25, 3.7361E-25, 2.7549E-25, &
      1.7754E-25, 1.2823E-25, 7.4012E-26, 6.3085E-26 /)

    REAL, PARAMETER :: a1_1_h_O2(dim55) = (/ &
      2.2019E-21, 2.0669E-21, 1.9421E-21, 1.7521E-21, 1.5821E-21, 1.4046E-21, &
      1.2432E-21, 1.1063E-21, 9.7442E-22, 8.7812E-22, 7.7660E-22, 7.1075E-22, &
      6.3191E-22, 5.8937E-22, 5.2695E-22, 4.9718E-22, 4.4878E-22, 4.2814E-22, &
      3.8663E-22, 3.7471E-22, 3.3779E-22, 3.3166E-22, 2.9899E-22, 2.9546E-22, &
      2.6873E-22, 2.6540E-22, 2.4191E-22, 2.4079E-22, 2.1877E-22, 2.2043E-22, &
      2.0100E-22, 2.0127E-22, 1.8433E-22, 1.8692E-22, 1.6986E-22, 1.7344E-22, &
      1.5765E-22, 1.6155E-22, 1.4597E-22, 1.5146E-22, 1.3272E-22, 1.1680E-22, &
      1.0094E-22, 9.1645E-23, 8.1263E-23, 7.4985E-23, 6.7353E-23, 6.2797E-23, &
      5.4034E-23, 4.9597E-23, 4.5948E-23, 4.2281E-23, 3.9107E-23, 3.6248E-23, &
      3.3680E-23 /)
    REAL, PARAMETER :: b1_1_h_O2(dim55) = (/ &
      6.4259E-20, 6.4664E-20, 6.5351E-20, 6.6871E-20, 6.8656E-20, 7.0963E-20, &
      7.3465E-20, 7.5929E-20, 7.8633E-20, 8.0848E-20, 8.3436E-20, 8.5280E-20, &
      8.7685E-20, 8.9089E-20, 9.1304E-20, 9.2436E-20, 9.4396E-20, 9.5283E-20, &
      9.7172E-20, 9.7744E-20, 9.9608E-20, 9.9934E-20, 1.0175E-19, 1.0195E-19, &
      1.0357E-19, 1.0378E-19, 1.0532E-19, 1.0539E-19, 1.0695E-19, 1.0682E-19, &
      1.0829E-19, 1.0827E-19, 1.0963E-19, 1.0942E-19, 1.1088E-19, 1.1056E-19, &
      1.1199E-19, 1.1163E-19, 1.1312E-19, 1.1258E-19, 1.1446E-19, 1.1626E-19, &
      1.1825E-19, 1.1953E-19, 1.2110E-19, 1.2212E-19, 1.2346E-19, 1.2432E-19, &
      1.2607E-19, 1.2702E-19, 1.2784E-19, 1.2871E-19, 1.2951E-19, 1.3026E-19, &
      1.3097E-19 /)
    REAL, PARAMETER :: a1_2_h_O2(dim55) = (/ &
      -2.8240E-24, -6.1030E-25, 9.4608E-26, 1.5331E-24, 2.0466E-24, 2.5907E-24, &
      2.8439E-24, 2.7861E-24, 2.8797E-24, 2.5420E-24, 2.5759E-24, 2.1666E-24, &
      2.2047E-24, 1.7922E-24, 1.8569E-24, 1.4870E-24, 1.5572E-24, 1.2428E-24, &
      1.3272E-24, 1.0399E-24, 1.1505E-24, 8.8128E-25, 9.9792E-25, 7.6595E-25, &
      8.5693E-25, 6.6778E-25, 7.6372E-25, 5.8658E-25, 6.8579E-25, 5.1267E-25, &
      6.0524E-25, 4.6585E-25, 5.3884E-25, 4.1184E-25, 4.9356E-25, 3.6429E-25, &
      4.4576E-25, 3.3411E-25, 4.1076E-25, 2.9950E-25, 3.0936E-25, 2.3719E-25, &
      2.1037E-25, 1.6574E-25, 1.5159E-25, 1.2310E-25, 1.1649E-25, 9.6794E-26, &
      1.1037E-25, 1.0128E-25, 9.2941E-26, 8.6866E-26, 8.0930E-26, 7.6328E-26, &
      7.1344E-26 /)
    REAL, PARAMETER :: b1_2_h_O2(dim55) = (/ &
      -1.4765E-22, -1.5429E-22, -1.5816E-22, -1.6967E-22, -1.7506E-22, -1.8214E-22, &
      -1.8606E-22, -1.8502E-22, -1.8694E-22, -1.7917E-22, -1.8004E-22, -1.6858E-22, &
      -1.6974E-22, -1.5613E-22, -1.5842E-22, -1.4437E-22, -1.4721E-22, -1.3369E-22, &
      -1.3753E-22, -1.2374E-22, -1.2933E-22, -1.1506E-22, -1.2153E-22, -1.0808E-22, &
      -1.1358E-22, -1.0166E-22, -1.0795E-22, -9.5902E-23, -1.0290E-22, -9.0259E-23, &
      -9.7248E-23, -8.6376E-23, -9.2251E-23, -8.1710E-23, -8.8697E-23, -7.7321E-23, &
      -8.4694E-23, -7.4311E-23, -8.1631E-23, -7.0727E-23, -7.1719E-23, -6.3563E-23, &
      -6.0197E-23, -5.4039E-23, -5.1909E-23, -4.7265E-23, -4.6104E-23, -4.2402E-23, &
      -4.5123E-23, -4.3188E-23, -4.1307E-23, -3.9862E-23, -3.8375E-23, -3.7164E-23, &
      -3.5791E-23 /)

    REAL, PARAMETER :: c0_1_1_h_O2(dim58) = (/ &
      1.6544E-01, 1.6544E-01, 1.6544E-01, 1.5192E-01, 1.2308E-01, 8.0673E-02, &
      2.6507E-02, -3.2853E-02, -9.1680E-02, -1.4242E-01, -1.7924E-01, -1.9693E-01, &
      -1.9252E-01, -1.6599E-01, -1.1973E-01, -5.9841E-02, 6.2923E-03, 6.7739E-02, &
      1.1680E-01, 1.4474E-01, 1.4599E-01, 1.2135E-01, 7.1299E-02, 8.6124E-03, &
      -5.9255E-02, -1.1442E-01, -1.4553E-01, -1.4349E-01, -1.0447E-01, -3.2614E-02, &
      6.0645E-02, 1.5808E-01, 2.3844E-01, 2.8349E-01, 2.8131E-01, 2.2803E-01, &
      1.3163E-01, 6.2099E-03, -1.2865E-01, -2.5422E-01, -3.5418E-01, -4.1854E-01, &
      -4.4227E-01, -4.2963E-01, -3.8962E-01, -3.3586E-01, -2.8065E-01, -2.3321E-01, &
      -1.9654E-01, -1.6878E-01, -1.4393E-01, -1.1473E-01, -7.5851E-02, -2.7431E-02, &
      2.0787E-02, 4.5287E-02, 1.9248E-02, -4.4049E-02 /)
    REAL, PARAMETER :: c0_1_2_h_O2(dim58) = (/ &
      -6.7124E+00, -6.7124E+00, -6.7124E+00, -6.1096E+00, -4.8175E+00, -2.9091E+00, &
      -4.6080E-01, 2.2341E+00, 4.9166E+00, 7.2408E+00, 8.9341E+00, 9.7514E+00, &
      9.5470E+00, 8.3107E+00, 6.1459E+00, 3.3309E+00, 2.0938E-01, -2.7032E+00, &
      -5.0387E+00, -6.3741E+00, -6.4340E+00, -5.2461E+00, -2.8239E+00, 2.2265E-01, &
      3.5346E+00, 6.2379E+00, 7.7682E+00, 7.6675E+00, 5.7322E+00, 2.1537E+00, &
      -2.5093E+00, -7.4007E+00, -1.1451E+01, -1.3730E+01, -1.3620E+01, -1.0902E+01, &
      -5.9665E+00, 4.8011E-01, 7.4387E+00, 1.3944E+01, 1.9141E+01, 2.2501E+01, &
      2.3744E+01, 2.3080E+01, 2.0967E+01, 1.8118E+01, 1.5180E+01, 1.2647E+01, &
      1.0682E+01, 9.1883E+00, 7.8462E+00, 6.2635E+00, 4.1487E+00, 1.5049E+00, &
      -1.1374E+00, -2.4849E+00, -1.0475E+00, 2.4591E+00 /)
    REAL, PARAMETER :: c0_2_1_h_O2(dim58) = (/ &
      1.0429E-01, 1.0429E-01, 1.0429E-01, 8.9880E-02, 5.9879E-02, 1.8582E-02, &
      -3.3983E-02, -8.4684E-02, -1.3260E-01, -1.6597E-01, -1.8676E-01, -1.8794E-01, &
      -1.7532E-01, -1.5047E-01, -1.1289E-01, -6.8522E-02, -2.1417E-02, 2.7059E-02, &
      6.9956E-02, 1.0657E-01, 1.1728E-01, 1.1740E-01, 8.1899E-02, 5.5565E-02, &
      2.2491E-02, 6.4790E-03, 2.1389E-02, 5.3682E-02, 9.1237E-02, 1.1849E-01, &
      1.3130E-01, 1.2411E-01, 1.0848E-01, 7.6345E-02, 3.6820E-02, -2.1606E-02, &
      -1.0270E-01, -1.8904E-01, -2.6165E-01, -2.9980E-01, -2.7215E-01, -2.0665E-01, &
      -1.1997E-01, -4.6731E-02, -8.1621E-03, -6.5224E-03, -1.3729E-02, -1.3854E-02, &
      3.2013E-04, 3.6318E-02, 7.3051E-02, 1.0182E-01, 1.0176E-01, 5.7931E-02, &
      7.1388E-04, -5.6235E-02, -8.9184E-02, -1.0821E-01 /)
    REAL, PARAMETER :: c0_2_2_h_O2(dim58) = (/ &
      -4.4111E+00, -4.4111E+00, -4.4111E+00, -3.7683E+00, -2.4243E+00, -5.6593E-01, &
      1.8100E+00, 4.1118E+00, 6.2966E+00, 7.8252E+00, 8.7815E+00, 8.8358E+00, &
      8.2505E+00, 7.0922E+00, 5.3337E+00, 3.2483E+00, 1.0250E+00, -1.2728E+00, &
      -3.3147E+00, -5.0647E+00, -5.5790E+00, -5.5848E+00, -3.8664E+00, -2.5866E+00, &
      -9.7255E-01, -1.8798E-01, -9.2156E-01, -2.5168E+00, -4.3795E+00, -5.7367E+00, &
      -6.3771E+00, -6.0161E+00, -5.2284E+00, -3.6025E+00, -1.5947E+00, 1.3851E+00, &
      5.5370E+00, 9.9751E+00, 1.3722E+01, 1.5698E+01, 1.4260E+01, 1.0841E+01, &
      6.2990E+00, 2.4465E+00, 4.1005E-01, 3.2315E-01, 7.0652E-01, 7.1321E-01, &
      -4.6522E-02, -1.9832E+00, -3.9668E+00, -5.5259E+00, -5.5226E+00, -3.1297E+00, &
      5.7547E-03, 3.1380E+00, 4.9567E+00, 6.0108E+00 /)

    !---------------------------------------------------------------------------
    ! heating rate O3
    !---------------------------------------------------------------------------

    ! interval 0
    REAL, PARAMETER :: a0_1_1_h_o3(dim58) = (/ &
      -1.2137E-22, -1.3140E-22, -1.4188E-22, -1.5273E-22, -1.6372E-22, -1.7487E-22, &
      -1.8603E-22, -1.9730E-22, -2.0862E-22, -2.1991E-22, -2.3132E-22, -2.4281E-22, &
      -2.6548E-22, -3.0836E-22, -3.0667E-22, -3.1524E-22, -3.3291E-22, -3.2967E-22, &
      -3.3507E-22, -3.4769E-22, -3.4126E-22, -3.4203E-22, -3.4898E-22, -3.3950E-22, &
      -3.3699E-22, -3.2199E-22, -3.1457E-22, -3.1286E-22, -3.1057E-22, -3.0750E-22, &
      -3.0359E-22, -2.9838E-22, -2.9225E-22, -2.8504E-22, -2.7745E-22, -2.6947E-22, &
      -2.6210E-22, -2.5539E-22, -2.4952E-22, -2.4461E-22, -2.4048E-22, -2.3635E-22, &
      -2.3147E-22, -2.2516E-22, -2.1654E-22, -2.0541E-22, -1.9180E-22, -1.7648E-22, &
      -1.6050E-22, -1.4439E-22, -1.2942E-22, -1.1511E-22, -1.0167E-22, -8.8738E-23, &
      -7.6149E-23, -6.4073E-23, -5.2003E-23, -4.2481E-23 /)
    REAL, PARAMETER :: a0_1_2_h_o3(dim58) = (/ &
      1.4659E-20, 1.5102E-20, 1.5568E-20, 1.6051E-20, 1.6544E-20, 1.7046E-20, &
      1.7550E-20, 1.8061E-20, 1.8578E-20, 1.9095E-20, 1.9619E-20, 2.0150E-20, &
      2.1202E-20, 2.3201E-20, 2.3121E-20, 2.3524E-20, 2.4358E-20, 2.4205E-20, &
      2.4462E-20, 2.5065E-20, 2.4757E-20, 2.4793E-20, 2.5130E-20, 2.4669E-20, &
      2.4547E-20, 2.3812E-20, 2.3447E-20, 2.3362E-20, 2.3249E-20, 2.3096E-20, &
      2.2900E-20, 2.2639E-20, 2.2330E-20, 2.1965E-20, 2.1579E-20, 2.1172E-20, &
      2.0795E-20, 2.0450E-20, 2.0147E-20, 1.9893E-20, 1.9678E-20, 1.9462E-20, &
      1.9207E-20, 1.8875E-20, 1.8420E-20, 1.7830E-20, 1.7106E-20, 1.6288E-20, &
      1.5431E-20, 1.4565E-20, 1.3756E-20, 1.2980E-20, 1.2249E-20, 1.1543E-20, &
      1.0853E-20, 1.0189E-20, 9.5229E-21, 8.9954E-21 /)
    REAL, PARAMETER :: a0_2_1_h_o3(dim58) = (/ &
      4.7696E-25, 5.1460E-25, 5.5300E-25, 5.9348E-25, 6.3343E-25, 6.7431E-25, &
      7.1414E-25, 7.5407E-25, 7.9369E-25, 8.3277E-25, 8.7149E-25, 9.0974E-25, &
      1.1137E-24, 1.6029E-24, 1.4309E-24, 1.4183E-24, 1.5282E-24, 1.3746E-24, &
      1.3491E-24, 1.4189E-24, 1.2764E-24, 1.2346E-24, 1.2684E-24, 1.1360E-24, &
      1.0876E-24, 9.2160E-25, 8.3909E-25, 8.1852E-25, 7.9720E-25, 7.7394E-25, &
      7.5008E-25, 7.2306E-25, 6.9453E-25, 6.6430E-25, 6.3409E-25, 6.0376E-25, &
      5.7610E-25, 5.5069E-25, 5.2866E-25, 5.0869E-25, 4.9187E-25, 4.7467E-25, &
      4.5521E-25, 4.3299E-25, 4.0531E-25, 3.7197E-25, 3.3457E-25, 2.9418E-25, &
      2.5496E-25, 2.1692E-25, 1.8378E-25, 1.5360E-25, 1.2673E-25, 1.0258E-25, &
      8.0461E-26, 6.2046E-26, 4.6123E-26, 3.5804E-26 /)
    REAL, PARAMETER :: a0_2_2_h_o3(dim58) = (/ &
      -2.9174E-23, -3.0837E-23, -3.2542E-23, -3.4347E-23, -3.6137E-23, -3.7977E-23, &
      -3.9777E-23, -4.1590E-23, -4.3397E-23, -4.5187E-23, -4.6968E-23, -4.8735E-23, &
      -5.8198E-23, -8.0994E-23, -7.2944E-23, -7.2355E-23, -7.7542E-23, -7.0259E-23, &
      -6.9046E-23, -7.2383E-23, -6.5541E-23, -6.3528E-23, -6.5162E-23, -5.8731E-23, &
      -5.6366E-23, -4.8234E-23, -4.4174E-23, -4.3158E-23, -4.2100E-23, -4.0942E-23, &
      -3.9749E-23, -3.8392E-23, -3.6955E-23, -3.5425E-23, -3.3890E-23, -3.2344E-23, &
      -3.0927E-23, -2.9621E-23, -2.8484E-23, -2.7450E-23, -2.6576E-23, -2.5678E-23, &
      -2.4658E-23, -2.3489E-23, -2.2028E-23, -2.0260E-23, -1.8271E-23, -1.6114E-23, &
      -1.4012E-23, -1.1965E-23, -1.0176E-23, -8.5399E-24, -7.0783E-24, -5.7600E-24, &
      -4.5477E-24, -3.5348E-24, -2.6559E-24, -2.0842E-24 /)

    REAL, PARAMETER :: a1_1_h_o3(dim55) = (/ &
      -2.6170E-21, -2.3714E-21, -2.1105E-21, -1.8362E-21, -1.5606E-21, -1.3370E-21, &
      -1.1074E-21, -9.5439E-22, -7.8330E-22, -6.8886E-22, -5.6589E-22, -5.1038E-22, &
      -4.2194E-22, -3.9035E-22, -3.2535E-22, -3.0790E-22, -2.5785E-22, -2.4914E-22, &
      -2.0979E-22, -2.0582E-22, -1.7405E-22, -1.7329E-22, -1.4647E-22, -1.4801E-22, &
      -1.2521E-22, -1.2786E-22, -1.0827E-22, -1.1177E-22, -9.4436E-23, -9.8422E-23, &
      -8.3150E-23, -8.7506E-23, -7.3654E-23, -7.8396E-23, -6.5949E-23, -7.0585E-23, &
      -5.9312E-23, -6.4041E-23, -5.3704E-23, -5.8376E-23, -4.7468E-23, -4.0009E-23, &
      -3.2526E-23, -2.8732E-23, -2.4054E-23, -2.1918E-23, -1.8685E-23, -1.7406E-23, &
      -1.3725E-23, -1.2469E-23, -1.1357E-23, -1.0366E-23, -9.4866E-24, -8.6851E-24, &
      -7.9774E-24 /)
    REAL, PARAMETER :: b1_1_h_o3(dim55) = (/ &
      6.5534E-20, 6.4797E-20, 6.3362E-20, 6.1168E-20, 5.8273E-20, 5.5367E-20, &
      5.1808E-20, 4.9054E-20, 4.5547E-20, 4.3374E-20, 4.0239E-20, 3.8684E-20, &
      3.5987E-20, 3.4945E-20, 3.2637E-20, 3.1974E-20, 2.9947E-20, 2.9572E-20, &
      2.7782E-20, 2.7592E-20, 2.5987E-20, 2.5947E-20, 2.4459E-20, 2.4548E-20, &
      2.3168E-20, 2.3335E-20, 2.2052E-20, 2.2290E-20, 2.1068E-20, 2.1359E-20, &
      2.0206E-20, 2.0546E-20, 1.9431E-20, 1.9824E-20, 1.8760E-20, 1.9168E-20, &
      1.8148E-20, 1.8588E-20, 1.7600E-20, 1.8058E-20, 1.6962E-20, 1.6119E-20, &
      1.5180E-20, 1.4656E-20, 1.3952E-20, 1.3604E-20, 1.3037E-20, 1.2796E-20, &
      1.2058E-20, 1.1791E-20, 1.1540E-20, 1.1304E-20, 1.1084E-20, 1.0873E-20, &
      1.0678E-20 /)
    REAL, PARAMETER :: a1_2_h_o3(dim55) = (/ &
      2.2344E-25, -1.9499E-24, -3.1949E-24, -4.0979E-24, -4.8478E-24, -4.6437E-24, &
      -4.9602E-24, -4.1908E-24, -4.2950E-24, -3.3922E-24, -3.4410E-24, -2.6204E-24, &
      -2.6772E-24, -1.9974E-24, -2.0731E-24, -1.5289E-24, -1.6252E-24, -1.1887E-24, &
      -1.2924E-24, -9.4153E-25, -1.0454E-24, -7.5764E-25, -8.6219E-25, -6.1992E-25, &
      -7.1943E-25, -5.1542E-25, -6.0862E-25, -4.3240E-25, -5.2051E-25, -3.6836E-25, &
      -4.4937E-25, -3.1544E-25, -3.9208E-25, -2.7231E-25, -3.4280E-25, -2.3761E-25, &
      -3.0238E-25, -2.0755E-25, -2.6845E-25, -1.8293E-25, -1.9049E-25, -1.3634E-25, &
      -1.1541E-25, -8.5676E-26, -7.6228E-26, -5.8217E-26, -5.3957E-26, -4.2072E-26, &
      -4.7639E-26, -4.1364E-26, -3.6349E-26, -3.2286E-26, -2.8904E-26, -2.6126E-26, &
      -2.3719E-26 /)
    REAL, PARAMETER :: b1_2_h_o3(dim55) = (/ &
      1.7664E-22, 1.8316E-22, 1.9001E-22, 1.9723E-22, 2.0511E-22, 2.0245E-22, &
      2.0736E-22, 1.9351E-22, 1.9565E-22, 1.7488E-22, 1.7613E-22, 1.5315E-22, &
      1.5488E-22, 1.3245E-22, 1.3514E-22, 1.1445E-22, 1.1836E-22, 9.9587E-23, &
      1.0431E-22, 8.7463E-23, 9.2708E-23, 7.7458E-23, 8.3260E-23, 6.9208E-23, &
      7.5229E-23, 6.2376E-23, 6.8481E-23, 5.6498E-23, 6.2710E-23, 5.1602E-23, &
      5.7719E-23, 4.7272E-23, 5.3441E-23, 4.3501E-23, 4.9528E-23, 4.0271E-23, &
      4.6132E-23, 3.7313E-23, 4.3130E-23, 3.4749E-23, 3.5508E-23, 2.9390E-23, &
      2.6762E-23, 2.2659E-23, 2.1237E-23, 1.8302E-23, 1.7554E-23, 1.5320E-23, &
      1.6436E-23, 1.5099E-23, 1.3968E-23, 1.3001E-23, 1.2154E-23, 1.1423E-23, &
      1.0760E-23 /)

    REAL, PARAMETER :: c0_1_1_h_o3(dim58) = (/ &
      -7.9855E-03, -7.9855E-03, -7.9855E-03, -9.3673E-03, -1.0438E-02, -1.1734E-02, &
      -1.2552E-02, -1.3566E-02, -1.4145E-02, -1.4611E-02, -1.4620E-02, -1.4621E-02, &
      -1.4551E-02, -1.4312E-02, -1.3990E-02, -1.3951E-02, -1.3725E-02, -1.3948E-02, &
      -1.4112E-02, -1.4251E-02, -1.4297E-02, -1.4233E-02, -1.3704E-02, -1.2551E-02, &
      -1.0726E-02, -8.6198E-03, -5.9638E-03, -2.8771E-03, 1.2781E-04, 2.6547E-03, &
      4.2461E-03, 4.8933E-03, 3.9358E-03, 1.7863E-03, -1.5202E-03, -5.4635E-03, &
      -9.3239E-03, -1.2417E-02, -1.3942E-02, -1.3267E-02, -1.0020E-02, -4.5719E-03, &
      2.5238E-03, 1.0490E-02, 1.7436E-02, 2.2612E-02, 2.5110E-02, 2.5261E-02, &
      2.4053E-02, 2.2477E-02, 2.2555E-02, 2.4513E-02, 2.8877E-02, 3.3852E-02, &
      3.8176E-02, 4.0812E-02, 4.3180E-02, 4.7724E-02 /)
    REAL, PARAMETER :: c0_1_2_h_o3(dim58) = (/ &
      3.5437E-01, 3.5437E-01, 3.5437E-01, 4.1600E-01, 4.6396E-01, 5.2227E-01, &
      5.5928E-01, 6.0531E-01, 6.3172E-01, 6.5306E-01, 6.5348E-01, 6.5350E-01, &
      6.5028E-01, 6.3913E-01, 6.2405E-01, 6.2220E-01, 6.1157E-01, 6.2210E-01, &
      6.2991E-01, 6.3658E-01, 6.3877E-01, 6.3570E-01, 6.1010E-01, 5.5403E-01, &
      4.6498E-01, 3.6179E-01, 2.3111E-01, 7.8634E-02, -7.0411E-02, -1.9625E-01, &
      -2.7582E-01, -3.0831E-01, -2.6005E-01, -1.5129E-01, 1.6683E-02, 2.1779E-01, &
      4.1544E-01, 5.7443E-01, 6.5314E-01, 6.1815E-01, 4.4930E-01, 1.6492E-01, &
      -2.0689E-01, -6.2594E-01, -9.9266E-01, -1.2670E+00, -1.3999E+00, -1.4080E+00, &
      -1.3432E+00, -1.2584E+00, -1.2626E+00, -1.3688E+00, -1.6062E+00, -1.8778E+00, &
      -2.1148E+00, -2.2597E+00, -2.3905E+00, -2.6422E+00 /)
    REAL, PARAMETER :: c0_2_1_h_o3(dim58) = (/ &
      -1.6165E-03, -1.6165E-03, -1.6165E-03, -2.8232E-03, -1.9496E-03, -2.9353E-03, &
      -3.6045E-03, -2.5226E-03, -3.2650E-03, -1.6102E-03, -1.0882E-03, -2.6250E-04, &
      1.5112E-03, 1.6614E-03, 3.8094E-03, 4.2192E-03, 5.7008E-03, 5.8037E-03, &
      6.4124E-03, 6.8792E-03, 7.4584E-03, 7.0290E-03, 6.6888E-03, 7.5922E-03, &
      5.7979E-03, 6.1172E-03, 5.7426E-03, 4.6160E-03, 2.5640E-03, 5.2821E-04, &
      -7.7150E-04, -1.7893E-03, -3.3838E-03, -2.3750E-03, -2.6544E-03, -2.6404E-03, &
      -3.7134E-04, 2.0935E-03, 6.6356E-03, 1.2796E-02, 1.8694E-02, 2.3208E-02, &
      2.4539E-02, 2.3023E-02, 1.6766E-02, 8.2739E-03, 6.8675E-03, 4.5913E-03, &
      7.0203E-03, 6.2734E-03, 5.8744E-03, -4.4369E-03, -1.8659E-02, -3.6573E-02, &
      -4.7990E-02, -4.6081E-02, -2.7760E-02, -5.5219E-03 /)
    REAL, PARAMETER :: c0_2_2_h_o3(dim58) = (/ &
      1.0145E-01, 1.0145E-01, 1.0145E-01, 1.5527E-01, 1.1613E-01, 1.6049E-01, &
      1.9073E-01, 1.4162E-01, 1.7547E-01, 9.9681E-02, 7.5667E-02, 3.7520E-02, &
      -4.4779E-02, -5.1778E-02, -1.5231E-01, -1.7156E-01, -2.4150E-01, -2.4637E-01, &
      -2.7535E-01, -2.9766E-01, -3.2546E-01, -3.0477E-01, -2.8830E-01, -3.3221E-01, &
      -2.4465E-01, -2.6029E-01, -2.4186E-01, -1.8620E-01, -8.4424E-02, 1.6956E-02, &
      8.1942E-02, 1.3304E-01, 2.1340E-01, 1.6235E-01, 1.7655E-01, 1.7583E-01, &
      5.9656E-02, -6.7038E-02, -3.0141E-01, -6.2051E-01, -9.2723E-01, -1.1629E+00, &
      -1.2326E+00, -1.1528E+00, -8.2247E-01, -3.7240E-01, -2.9758E-01, -1.7603E-01, &
      -3.0622E-01, -2.6604E-01, -2.4449E-01, 3.1438E-01, 1.0881E+00, 2.0662E+00, &
      2.6918E+00, 2.5868E+00, 1.5755E+00, 3.4351E-01 /)

    !---------------------------------------------------------------------------
    ! tau_i = effective optical depth
    !---------------------------------------------------------------------------

    ! lookup tables

    REAL, PARAMETER :: taua1_1(dim55) = (/ &
      9.0291E-02, 8.0836E-02, 7.2346E-02, 6.4955E-02, 5.8596E-02, 5.3170E-02, &
      4.8646E-02, 4.4776E-02, 4.1517E-02, 3.8734E-02, 3.6428E-02, 3.4335E-02, &
      3.2609E-02, 3.1047E-02, 2.9760E-02, 2.8467E-02, 2.7479E-02, 2.6459E-02, &
      2.5659E-02, 2.4815E-02, 2.4143E-02, 2.3457E-02, 2.2924E-02, 2.2285E-02, &
      2.1837E-02, 2.1334E-02, 2.0953E-02, 2.0441E-02, 2.0154E-02, 1.9720E-02, &
      1.9450E-02, 1.9079E-02, 1.8800E-02, 1.8512E-02, 1.8282E-02, 1.7964E-02, &
      1.7786E-02, 1.7494E-02, 1.7336E-02, 1.7085E-02, 1.6546E-02, 1.5734E-02, &
      1.5083E-02, 1.4513E-02, 1.4052E-02, 1.3613E-02, 1.3269E-02, 1.2923E-02, &
      1.2679E-02, 1.2418E-02, 1.2188E-02, 1.1981E-02, 1.1782E-02, 1.1607E-02, &
      1.1458E-02 /)
    REAL, PARAMETER :: taub1_1(dim55) = (/ &
      3.8714E-03, 3.2236E-02, 7.8932E-02, 1.3806E-01, 2.0482E-01, 2.7537E-01, &
      3.4548E-01, 4.1516E-01, 4.8196E-01, 5.4596E-01, 6.0477E-01, 6.6336E-01, &
      7.1601E-01, 7.6758E-01, 8.1327E-01, 8.6239E-01, 9.0242E-01, 9.4625E-01, &
      9.8267E-01, 1.0232E+00, 1.0571E+00, 1.0935E+00, 1.1230E+00, 1.1601E+00, &
      1.1872E+00, 1.2189E+00, 1.2439E+00, 1.2787E+00, 1.2989E+00, 1.3306E+00, &
      1.3510E+00, 1.3799E+00, 1.4024E+00, 1.4263E+00, 1.4460E+00, 1.4740E+00, &
      1.4900E+00, 1.5172E+00, 1.5323E+00, 1.5568E+00, 1.6111E+00, 1.7028E+00, &
      1.7845E+00, 1.8632E+00, 1.9325E+00, 2.0041E+00, 2.0645E+00, 2.1296E+00, &
      2.1785E+00, 2.2339E+00, 2.2859E+00, 2.3352E+00, 2.3851E+00, 2.4309E+00, &
      2.4721E+00 /)
    REAL, PARAMETER :: taua1_2(dim55) = (/ &
      2.9405E-04, 2.7624E-04, 2.5622E-04, 2.3168E-04, 2.0689E-04, 1.8370E-04, &
      1.6043E-04, 1.4216E-04, 1.2519E-04, 1.1208E-04, 9.7739E-05, 9.0017E-05, &
      7.9158E-05, 7.3840E-05, 6.4827E-05, 6.2338E-05, 5.4973E-05, 5.3183E-05, &
      4.6960E-05, 4.6053E-05, 4.1414E-05, 4.0288E-05, 3.5728E-05, 3.6471E-05, &
      3.2043E-05, 3.1950E-05, 2.8408E-05, 2.9726E-05, 2.5810E-05, 2.6529E-05, &
      2.3376E-05, 2.4114E-05, 2.2203E-05, 2.1965E-05, 1.9841E-05, 2.0660E-05, &
      1.8329E-05, 1.9302E-05, 1.7117E-05, 1.7776E-05, 1.5480E-05, 1.3571E-05, &
      1.1671E-05, 1.0559E-05, 9.1366E-06, 8.6354E-06, 7.5732E-06, 7.2476E-06, &
      5.9600E-06, 5.5182E-06, 5.0433E-06, 4.5682E-06, 4.3248E-06, 3.9847E-06, &
      3.5182E-06 /)
    REAL, PARAMETER :: taub1_2(dim55) = (/ &
      9.1337E-03, 9.1871E-03, 9.2972E-03, 9.4935E-03, 9.7538E-03, 1.0055E-02, &
      1.0416E-02, 1.0745E-02, 1.1093E-02, 1.1394E-02, 1.1760E-02, 1.1976E-02, &
      1.2307E-02, 1.2483E-02, 1.2803E-02, 1.2897E-02, 1.3196E-02, 1.3273E-02, &
      1.3556E-02, 1.3599E-02, 1.3834E-02, 1.3893E-02, 1.4146E-02, 1.4103E-02, &
      1.4371E-02, 1.4377E-02, 1.4609E-02, 1.4519E-02, 1.4796E-02, 1.4743E-02, &
      1.4981E-02, 1.4923E-02, 1.5077E-02, 1.5097E-02, 1.5279E-02, 1.5207E-02, &
      1.5418E-02, 1.5327E-02, 1.5536E-02, 1.5471E-02, 1.5702E-02, 1.5918E-02, &
      1.6156E-02, 1.6309E-02, 1.6524E-02, 1.6605E-02, 1.6792E-02, 1.6853E-02, &
      1.7111E-02, 1.7205E-02, 1.7312E-02, 1.7425E-02, 1.7486E-02, 1.7576E-02, &
      1.7704E-02 /)
    REAL, PARAMETER :: taua1_3(dim55) = (/ &
      -1.6566E-07, -3.5179E-08, 3.7882E-08, 1.1172E-07, 1.5734E-07, 1.8300E-07, &
      2.0711E-07, 1.9758E-07, 1.9658E-07, 1.7454E-07, 1.7817E-07, 1.4641E-07, &
      1.4890E-07, 1.1952E-07, 1.2673E-07, 9.7519E-08, 1.0423E-07, 8.0564E-08, &
      8.8805E-08, 6.7977E-08, 7.3188E-08, 5.8858E-08, 6.6132E-08, 4.7268E-08, &
      5.6996E-08, 4.3654E-08, 5.0776E-08, 3.5295E-08, 4.4550E-08, 3.2703E-08, &
      4.0261E-08, 2.9426E-08, 3.2833E-08, 2.7080E-08, 3.1572E-08, 2.3149E-08, &
      2.8543E-08, 2.0524E-08, 2.5594E-08, 1.9131E-08, 1.9420E-08, 1.4776E-08, &
      1.2911E-08, 1.0184E-08, 9.6719E-09, 7.3939E-09, 7.1189E-09, 5.6010E-09, &
      6.7829E-09, 6.0132E-09, 5.5800E-09, 5.3370E-09, 4.6488E-09, 4.4309E-09, &
      4.4105E-09 /)
    REAL, PARAMETER :: taub1_3(dim55) = (/ &
      -9.8469E-06, -1.0238E-05, -1.0640E-05, -1.1231E-05, -1.1710E-05, -1.2044E-05, &
      -1.2417E-05, -1.2246E-05, -1.2225E-05, -1.1718E-05, -1.1811E-05, -1.0922E-05, &
      -1.0997E-05, -1.0028E-05, -1.0284E-05, -9.1739E-06, -9.4458E-06, -8.4281E-06, &
      -8.8030E-06, -7.8033E-06, -8.0664E-06, -7.3070E-06, -7.7107E-06, -6.6165E-06, &
      -7.2051E-06, -6.3645E-06, -6.8311E-06, -5.7783E-06, -6.4308E-06, -5.5659E-06, &
      -6.1366E-06, -5.2914E-06, -5.5657E-06, -5.0882E-06, -5.4723E-06, -4.7311E-06, &
      -5.2192E-06, -4.4734E-06, -4.9576E-06, -4.3242E-06, -4.3533E-06, -3.8286E-06, &
      -3.5945E-06, -3.2182E-06, -3.1411E-06, -2.7698E-06, -2.7215E-06, -2.4361E-06, &
      -2.6731E-06, -2.5092E-06, -2.4115E-06, -2.3537E-06, -2.1813E-06, -2.1240E-06, &
      -2.1183E-06 /)

    REAL, PARAMETER ::  taua2(dim55) = (/ &
      1.4994E-01, 1.3078E-01, 1.1405E-01, 1.0032E-01, 8.9440E-02, 8.1120E-02, &
      7.4800E-02, 7.0080E-02, 6.6440E-02, 6.3680E-02, 6.1520E-02, 5.9840E-02, &
      5.8400E-02, 5.7320E-02, 5.6360E-02, 5.5560E-02, 5.4840E-02, 5.4280E-02, &
      5.3720E-02, 5.3280E-02, 5.2840E-02, 5.2480E-02, 5.2120E-02, 5.1840E-02, &
      5.1520E-02, 5.1280E-02, 5.1080E-02, 5.0800E-02, 5.0640E-02, 5.0440E-02, &
      5.0240E-02, 5.0080E-02, 4.9960E-02, 4.9800E-02, 4.9680E-02, 4.9560E-02, &
      4.9440E-02, 4.9360E-02, 4.9240E-02, 4.9120E-02, 4.8912E-02, 4.8592E-02, &
      4.8376E-02, 4.8208E-02, 4.8096E-02, 4.8032E-02, 4.7920E-02, 4.7920E-02, &
      4.7840E-02, 4.7840E-02, 4.7840E-02, 4.7840E-02, 4.7840E-02, 4.7840E-02, &
      4.7760E-02 /)
    REAL, PARAMETER :: taub2(dim55) = (/ &
      6.1088E-03, 6.3602E-02, 1.5558E-01, 2.6544E-01, 3.7968E-01, 4.8784E-01, &
      5.8580E-01, 6.7076E-01, 7.4538E-01, 8.0886E-01, 8.6394E-01, 9.1098E-01, &
      9.5490E-01, 9.9054E-01, 1.0246E+00, 1.0550E+00, 1.0842E+00, 1.1083E+00, &
      1.1337E+00, 1.1549E+00, 1.1771E+00, 1.1962E+00, 1.2161E+00, 1.2324E+00, &
      1.2517E+00, 1.2669E+00, 1.2800E+00, 1.2990E+00, 1.3103E+00, 1.3249E+00, &
      1.3400E+00, 1.3525E+00, 1.3621E+00, 1.3754E+00, 1.3856E+00, 1.3962E+00, &
      1.4071E+00, 1.4145E+00, 1.4260E+00, 1.4377E+00, 1.4586E+00, 1.4948E+00, &
      1.5219E+00, 1.5451E+00, 1.5620E+00, 1.5724E+00, 1.5920E+00, 1.5920E+00, &
      1.6081E+00, 1.6081E+00, 1.6081E+00, 1.6081E+00, 1.6081E+00, 1.6081E+00, &
      1.6301E+00 /)
    REAL, PARAMETER :: taua3(dim55) = (/ &
      1.4627E-02, 1.2811E-02, 1.1330E-02, 1.0177E-02, 9.2931E-03, 8.6222E-03, &
      8.1183E-03, 7.7196E-03, 7.4002E-03, 7.1526E-03, 6.9397E-03, 6.7612E-03, &
      6.6134E-03, 6.4776E-03, 6.3723E-03, 6.2660E-03, 6.1781E-03, 6.0942E-03, &
      6.0288E-03, 5.9625E-03, 5.8986E-03, 5.8466E-03, 5.8010E-03, 5.7508E-03, &
      5.7069E-03, 5.6709E-03, 5.6332E-03, 5.5991E-03, 5.5671E-03, 5.5351E-03, &
      5.5094E-03, 5.4792E-03, 5.4553E-03, 5.4313E-03, 5.4095E-03, 5.3874E-03, &
      5.3674E-03, 5.3474E-03, 5.3296E-03, 5.3115E-03, 5.2628E-03, 5.1941E-03, &
      5.1366E-03, 5.0886E-03, 5.0520E-03, 5.0160E-03, 4.9760E-03, 4.9601E-03, &
      4.9281E-03, 4.9121E-03, 4.8882E-03, 4.8802E-03, 4.8681E-03, 4.8482E-03, &
      4.8403E-03 /)
    REAL, PARAMETER :: taub3(dim55) = (/ &
      5.1550E-04, 4.6892E-02, 1.2180E-01, 2.0895E-01, 2.9795E-01, 3.8227E-01, &
      4.5824E-01, 5.2830E-01, 5.9245E-01, 6.4837E-01, 7.0178E-01, 7.5102E-01, &
      7.9549E-01, 8.3976E-01, 8.7671E-01, 9.1672E-01, 9.5196E-01, 9.8770E-01, &
      1.0172E+00, 1.0488E+00, 1.0808E+00, 1.1082E+00, 1.1333E+00, 1.1623E+00, &
      1.1887E+00, 1.2112E+00, 1.2358E+00, 1.2589E+00, 1.2813E+00, 1.3045E+00, &
      1.3239E+00, 1.3473E+00, 1.3665E+00, 1.3863E+00, 1.4049E+00, 1.4243E+00, &
      1.4423E+00, 1.4608E+00, 1.4778E+00, 1.4955E+00, 1.5443E+00, 1.6217E+00, &
      1.6937E+00, 1.7598E+00, 1.8149E+00, 1.8735E+00, 1.9435E+00, 1.9735E+00, &
      2.0375E+00, 2.0715E+00, 2.1255E+00, 2.1445E+00, 2.1748E+00, 2.2271E+00, &
      2.2490E+00 /)

    ! coefficients for polynomials
    REAL, PARAMETER :: &
      B_4_O3(4) = (/ 9.6253E-03, 2.5826E-03, -2.9300E-07, 3.6315E-11 /), &
      B_5_O3(3) = (/ 8.5935E-03, 3.1464E-04, -3.0217E-08             /), &
      B_7_O3(2) = (/ 1.2885E-03, 6.1319E-05                          /)

    !-------------------------------------------------------------------------

    ! LOCAL variables

    REAL, DIMENSION(SIZE(slf_1d),0:SIZE(temp_2d,2)) :: &
      v2,    & ! O2 column density             [part./cm^2]
      v2s,   & ! slant O2 column density       [part./cm^2]
      v3,    & ! O3 column density             [part./cm^2]
      v3s,   & ! slant O3 column density       [part./cm^2]
      relo3    ! O3 mixing ratio               [mol/mol]

    REAL, DIMENSION(SIZE(slf_1d),SIZE(temp_2d,2)) :: &
      v2s2,     & ! v2s (interval 1, ma-corr)
      dh_O2,    & ! O2 heating rate                   [K/s]
      dh_o3,    & ! O3 heating rate                   [K/s]
      dv2,      & ! diff. O2 column density           [part./cm^2]
      dv3,      & ! diff. O3 column density           [part./cm^2]
      aclc,     & ! cloud fraction                    [1]
      clp,      & ! cloud liquid water path per layer [g/m^2]
      rhum,     & ! relative humudity                 [%]
      part_col, & ! aerosol part. column per layer    [part/cm^2]
      tau_sca,  & ! scattering optical depth (slingo cloud parameter)
      tau_abs,  & ! absorption optical depth (slingo cloud parameter)
      gcld        ! asymmetry factor (slingo cloud parameter)

    REAL, DIMENSION(SIZE(slf_1d)) :: &
      slf, & ! fraction of sea surface
      u0     ! cosine of solar zenith angle

    INTEGER, DIMENSION(SIZE(temp_2d,2)) :: &
      iaer  ! aerosol type (1-2)
      ! 1 = mixture between rural and maritime aerosol with slf
      ! 2 = free troposphere aerosol

    REAL, DIMENSION(SIZE(slf_1d),0:SIZE(temp_2d,2),MAXWAV) :: &
      fact    ! actinic flux [photons/(cm^2 s)]
    REAL, DIMENSION(SIZE(slf_1d),0:SIZE(temp_2d,2),3:5) :: &
      facth   ! flux [photons/(cm^2 s)]

    REAL, DIMENSION(SIZE(slf_1d),SIZE(temp_2d,2),MAXWAV) :: &
      taus_clr,   & ! scat. opt. depth of layer (clear sky)
      taua_clr,   & ! abs. opt. depth of layer (clear sky)
      g_clr,      & ! asymmerty factor of layer (clear sky)
      taus_cld,   & ! scat. opt. depth of layer (clouds)
      taua_cld,   & ! abs. opt. depth of layer (clouds)
      g_cld,      & ! asymmerty factor of layer (clouds)
      taer_sca,   & ! scattering optical depth of aerosol (Shettle, Fenn)
      taer_abs,   & ! absorption optical depth of aerosol (Shettle, Fenn)
      gaer          ! asymmetry factor of aerosol (Shettle, Fenn)

    REAL, DIMENSION(0:MAXWAV) :: sig_h_o3, sig_h_O2 ! heating rates

    REAL, DIMENSION(SIZE(slf_1d),MAXWAV) :: &
      albedo, & ! surface albedo
      ftop      ! absorption correction above model
    REAL, DIMENSION(SIZE(slf_1d),SIZE(temp_2d,2),4) :: &
      bb   ! initialization coefficients of pifm

    INTEGER :: kproma
    INTEGER :: i, j, k, l, km, kp

    REAL :: zen, dz, aa1
    REAL :: coszet, cossza_thr, part_surf
    REAL :: ta_O2, ta_o3, temp_avg, domin, a1_sc, a2_sc, a4_sc
    REAL :: b1_sc, b2_sc, b4_sc, a1_ab, a2_ab, a4_ab, b1_ab, b2_ab, b4_ab
    REAL :: a1_g, a2_g, a4_g, b1_g, b2_g, b4_g, bsa, baa, ga
    REAL :: dv2s1 ! dv2s (interval 0)
    REAL :: dv3s1 ! dv3s (interval 0)
    REAL :: v3s1  ! v3s (intervals 1-2)
    REAL :: v3s2  ! v3s (intervals 3,4,5,7)
    REAL :: dlv2_h, v3s_du_h
    REAL :: v3_du, sig0_o2k, sig0_o3k, dtau_0

    ! effective optical depths for the different intervals (tau_6 is not used)
    REAL :: tau_0(SIZE(slf_1d)), tau_1, tau_2, tau_3, tau_4, tau_5, tau_7

    !-------------------------------------------------------------------------

    kproma = SIZE(slf_1d)
    klev   = SIZE(temp_2d,2)

    ALLOCATE(v2s_m(kproma,klev))
    ALLOCATE(v3_du1(kproma,klev))
    ALLOCATE(v3_du2(kproma,klev))
    ALLOCATE(tnorm_sr(kproma,klev))
    ALLOCATE(tnorm(kproma,klev))
    ALLOCATE(dlv2(kproma,klev))
    ALLOCATE(dens(kproma,klev))
    ALLOCATE(r_m(kproma,klev))
    ALLOCATE(r_o2(kproma,klev))

    ALLOCATE(fint(kproma,klev,0:MAXWAV))
    ALLOCATE(finth(kproma,klev,3:5))
    ALLOCATE(fj_corr(kproma,8))

    ALLOCATE(iu0(kproma))
    ALLOCATE(i0(kproma,klev))
    ALLOCATE(i1(kproma,klev))
    ALLOCATE(i2(kproma,klev))
    ALLOCATE(i3(kproma,klev))

    ALLOCATE(sig_top_o2(kproma,0:klev))
    ALLOCATE(sig_top_o3(kproma,0:klev))
    ALLOCATE(temp(kproma,0:klev))
    ALLOCATE(press(kproma,0:klev))

    ! calculate filter function for sza le 0.02
    kproma_day = 1
    DO j = 1,kproma
      IF (cossza_1d(j) >= u0lim) THEN
        iu0(kproma_day) = j
        kproma_day = kproma_day + 1
        ! get latitude indices
      ENDIF
    ENDDO
    kproma_day = kproma_day -1
    ! aerosol: boundary layer assumed for the lowest pbllev layers
    iaer(:) = 2
    iaer(klev-pbllev+1:klev) = 1

    cossza_thr = COS(sza_thr*pi/180.)

    DO j = 1,kproma_day
      v3(j,:)     = v3_2d(iu0(j),:)
      press(j,1:) = press_2d(iu0(j),:)
      relo3(j,:)  = relo3_2d(iu0(j),:)
      temp(j,1:)  = temp_2d(iu0(j),:)
      rhum(j,:)   = MAX(0.,MIN(100.,rhum_2d(iu0(j),:)))
      aclc(j,:)   = aclc_2d(iu0(j),:)
      clp(j,:)    = clp_2d(iu0(j),:)
      albedo(j,:) = albedo_1d(iu0(j))
      slf(j)      = slf_1d(iu0(j))       ! land sea fraction

      ! Rodgers, The radiative heat budget of the troposphere and lower
      ! stratosphere, Res. Rep. n0. a2, planatary circualtions project (1967)
      coszet = MAX(cossza_1d(iu0(j)),cossza_thr)
      u0(j)  = SQRT(coszet**2*1224.+1.)*1./35.

      ! pressure and temperature at virtual level 0
      ! above level 0 an only absorbing astmosphere is assumed
      press(j,0) = 0.37 * press(j,1) ! 0.37 = 1/e
      temp(j,0)  = (temp(j,2)-temp(j,1))/(press(j,2)-press(j,1)) * &
                   (-0.63) * press(j,1) + temp(j,1)
      ! assume constant mixing ratio above toa
      ! fix ozone mixing ratio in the first nt_lev layers by haloe climatology
      ! assume constant mixing ratio above toa (level 1)

      ! aerosol particle column per layer [part./cm^2]
      ! Christoph's assumption about aerosol particle density:
      ! n(z) = par_surf * (press(z)/press(surface))^3
      ! part_surf: aerosol particle  at the surface [part./cm^3]
      part_surf = part_rural * slf(j) + part_sea * (1.-slf(j))
      DO k = 1,klev
        aa1 = -R_gas/(M_air_SI*g) *LOG(press(j,k-1)/press(j,k))
        dz = aa1*0.5*(temp(j,k-1)+temp(j,k))
        part_col(j,k) = part_surf * 0.5 * &
          ((press(j,k)/press(j,klev))**3 + (press(j,k-1)/press(j,klev))**3) &
          * dz * 100.
      ENDDO
    ENDDO

    IF (l_heating)       rh_o3_2d(:,:)      = 0.0_dp
    IF (l_heating)       rh_o2_2d(:,:)      = 0.0_dp

    !-------------------------------------------------------------------------

    CALL column_cal

    !-------------------------------------------------------------------------

    CALL flux_cal

    ! ------------------------------------------------------------------------

    ! do first the calculation of the optical depth tau_0 above
    ! the model atmosphere

    ! initialization of tau_0
    k = 1
    DO j = 1,kproma_day

      dlv2_h   = MAX(44.5,MIN(56.,LOG(v2s(j,0))))
      v3s_du_h = MIN(300., v3s(j,0)* constant * 1.E+3)
      i0(j,k)       = MIN(58,INT(AINT((dlv2_h-44.5)/0.2) + 1.00001))
      tnorm_sr(j,k) = (temp(j,0)-240.) / 240.
      sig_top_O2(j,0) = p2(a0_1_1_O2(i0(j,k))*dlv2_h + a0_1_2_O2(i0(j,k)), &
        a0_2_1_O2(i0(j,k))*dlv2_h + a0_2_2_O2(i0(j,k)), &
        a0_3_1_O2(i0(j,k))*dlv2_h + a0_3_2_O2(i0(j,k)), &
        v3s_du_h) * &
        p2(1.,c0_1_1_O2(i0(j,k))*dlv2_h+c0_1_2_O2(i0(j,k)), &
        c0_2_1_O2(i0(j,k))*dlv2_h+c0_2_2_O2(i0(j,k)), &
        tnorm_sr(j,k))
      sig_top_o3(j,0) = p1(a0_1_1_o3(i0(j,k))*dlv2_h + a0_1_2_o3(i0(j,k)), &
        a0_2_1_o3(i0(j,k))*dlv2_h + a0_2_2_o3(i0(j,k)), &
        v3s_du_h) * &
        p2(1.,c0_1_1_o3(i0(j,k))*dlv2_h+c0_1_2_o3(i0(j,k)), &
        c0_2_1_o3(i0(j,k))*dlv2_h+c0_2_2_o3(i0(j,k)), &
        tnorm_sr(j,k))
      tau_0(j)     = p2(ct(1),ct(2),ct(3),(dlv2_h-47.)/47.)
    ENDDO

    ! basic species, minimal set, definition of arrays for tables and fits
    DO k = 1,klev ! altitude loop

      DO j  = 1,kproma_day ! longitude loop

        ! change of units [v2s_m] = meter, [v3_du] = DU
        v2s_m(j,k) = v2s(j,k) * constant * 1.E-2
        v3_du      = v3s(j,k) * constant * 1.E+3
        dlv2(j,k)  = LOG(v2s(j,k))

        ! allowed ranges for columns
        ! interval: 0
        IF (dlv2(j,k)>= 56.) THEN
          dlv2(j,k)  = 56.
          dv2s1 = 1.E+25
        ELSE
          dv2s1     = v2s(j,k)-v2s(j,k-1)
        ENDIF

        ! definition range of polynomial
        ! interval: 1
        IF (v2s_m(j,k)>= 500.) THEN
          v2s_m(j,k) = 500.
          v2s2(j,k)  = 1.3602E+24
        ELSE
          v2s2(j,k)   = v2s(j,k)
        ENDIF

        ! interval: 0 - 2
        IF (v3_du>= 300.) THEN
          dv3s1     = 1.E+19
          v3_du1(j,k) = 300.
          v3s1      = 8.161E+18
        ELSE
          dv3s1     = v3s(j,k)-v3s(j,k-1)
          v3_du1(j,k) = v3_du
          v3s1      = v3s(j,k)
        ENDIF

        ! interval: 3 - 7
        IF (v3_du>= 3000.) THEN
          v3_du2(j,k) = 3000.
          v3s2      = 8.161E+19
        ELSE
          v3_du2(j,k) = v3_du
          v3s2      = v3s(j,k)
        ENDIF

        ! scaling of temperature variable
        tnorm_sr(j,k) = (temp(j,k)-240.) / 240. ! for interval 0
        tnorm(j,k)    = (temp(j,k)-250.) / 250. ! for interval 1-7
        dens(j,k)  =  PRESS(J,K)/(k_B*TEMP(J,K))*1.E-6

        ! indices for lookup table
        ! mz_rs_20140314+
        ! At very high altitudes, dlv2 can be smaller than 44.5 and the
        ! resulting index i0 will be out of the valid range. This bug
        ! fix sets i0=1 in that case.
        i0(j,k)  = MAX(1,MIN(58,INT(AINT((dlv2(j,k)-44.5)/0.2) + 1.00001)))
        !i0(j,k)  = MIN(58,INT(AINT((dlv2(j,k)-44.5)/0.2) + 1.00001))
        ! mz_rs_20140314-
        i1(j,k)  = MIN(INT((v3_du-0.5)/2.5) +1,115)
        i1(j,k)  = ifil(i1(j,k))
        i2(j,k)  = MIN(INT((v3_du-0.5)/2.5) +1,115)
        i2(j,k)  = ifil(i2(j,k))
        i3(j,k)  = MIN(INT((v3_du-5.)/25.) +1,115)
        i3(j,k)  = ifil(i3(j,k))

        ! calculation of optical depths and integrated actinic flux
        ! effective optical depths and integrated actinic fluxes
        sig0_O2k = &
          p2(a0_1_1_O2(i0(j,k))*dlv2(j,k) + a0_1_2_O2(i0(j,k)), &
          a0_2_1_O2(i0(j,k))*dlv2(j,k) + a0_2_2_O2(i0(j,k)), &
          a0_3_1_O2(i0(j,k))*dlv2(j,k) + a0_3_2_O2(i0(j,k)), &
          v3_du1(j,k)) * &
          p2(1.,c0_1_1_O2(i0(j,k))*dlv2(j,k)+c0_1_2_O2(i0(j,k)), &
          c0_2_1_O2(i0(j,k))*dlv2(j,k)+c0_2_2_O2(i0(j,k)), &
          tnorm_sr(j,k))

        sig0_o3k = &
          p1(a0_1_1_o3(i0(j,k))*dlv2(j,k) + a0_1_2_o3(i0(j,k)), &
          a0_2_1_o3(i0(j,k))*dlv2(j,k) + a0_2_2_o3(i0(j,k)), &
          v3_du1(j,k)) * &
          p2(1.,c0_1_1_o3(i0(j,k))*dlv2(j,k)+c0_1_2_o3(i0(j,k)), &
          c0_2_1_o3(i0(j,k))*dlv2(j,k)+c0_2_2_o3(i0(j,k)), &
          tnorm_sr(j,k))
        dtau_0     = &
          0.5*(sig_top_O2(j,k-1) + sig0_O2k) * dv2s1 + &
          0.5*(sig_top_o3(j,k-1) + sig0_o3k) * dv3s1
        tau_0(j)   = tau_0(j) + dtau_0

        sig_top_O2(j,k) = sig0_O2k
        sig_top_o3(j,k) = sig0_o3k

        ! calculate the effective optical depths tau_i
        ! intervals i = 1..3 from lookup tables
        tau_1 = &
          p1(taub1_1(i1(j,k)),taua1_1(i1(j,k)),v3_du1(j,k)) + &
          p1(taub1_2(i1(j,k)),taua1_2(i1(j,k)),v3_du1(j,k)) * &
          v2s_m(j,k) + &
          p1(taub1_3(i1(j,k)),taua1_3(i1(j,k)),v3_du1(j,k)) * &
          v2s_m(j,k)**2
        tau_2 = p1(taub2(i2(j,k)),taua2(i2(j,k)),v3_du1(j,k))
        tau_3 = p1(taub3(i3(j,k)),taua3(i3(j,k)),v3_du2(j,k))
        ! for intervals i = 4..7 from the polynomials with B_i_O3(j)
        tau_4 = p2(B_4_O3(1),B_4_O3(2),B_4_O3(3),v3_du2(j,k))
        tau_5 = p1(B_5_O3(1),B_5_O3(2),v3_du2(j,k))
        ! tau_6 is not used
        tau_7 = p1(B_7_O3(1),B_7_O3(2),v3_du2(j,k))
        fint(j,k,0) = EXP(-tau_0(j))* SR_toa_flux
        fint(j,k,1) = EXP(-tau_1 + crs_o3(1)*v3s1 + &
          crs_O2(1)*v2s2(j,k))*f0(1)*fact(j,k,1)
        fint(j,k,2) = EXP(-tau_2 + crs_o3(2)*v3s1)*f0(2)*fact(j,k,2)
        IF (v3_du2(j,k) <= 1000.) THEN
          fint(j,k,3) = EXP(-tau_3 + crs_o3(3)*v3s2)*f0(3)*fact(j,k,3)
          finth(j,k,3) = EXP(-tau_3 + crs_o3(3)*v3s2)*f0(3)*facth(j,k,3)
        ELSE
          fint(j,k,3) = 0.
          finth(j,k,3) = 0.
        ENDIF
        IF (v3_du2(j,k) <= 2500.) THEN
          fint(j,k,4) = EXP(-tau_4 + crs_o3(4)*v3s2)*f0(4)*fact(j,k,4)
          finth(j,k,4) = EXP(-tau_4 + crs_o3(4)*v3s2)*f0(4)*facth(j,k,4)
        ELSE
          fint(j,k,4) = 0.
          finth(j,k,4) = 0.
        ENDIF
        fint(j,k,5) = EXP(-tau_5 + crs_o3(5)*v3s2)*f0(5)*fact(j,k,5)
        finth(j,k,5) = EXP(-tau_5 + crs_o3(5)*v3s2)*f0(5)*facth(j,k,5)
        fint(j,k,6) = f0(6)*fact(j,k,6)
        fint(j,k,7) = EXP(-tau_7 + crs_o3(7)*v3s2)*f0(7)*fact(j,k,7)

      ENDDO
    ENDDO

    !-------------------------------------------------------------------------

    ! reduction factors r_O2 und r_m for lyman-alpha line
    ! S. Chabrillat and G. Kockarts grl 24, p. 2659-2662 [2633]
    ! with the corrected values from grl 25 p.79 [2634]

    ! zero for non-ma
    r_m(:,:)  = 0. ! for CO2, CH4, SF6, H2SO4, and H2O
    r_O2(:,:) = 0. ! only for O2
    ! lyman-alpha photolysis for ma (only relevant in the mesosphere)
    IF (lmidatm) THEN
      DO i = 1,3
        DO k = 1,klev
          DO j = 1,kproma_day
            r_m(j,k)  = r_m(j,k)  + b_la(i) * EXP(-c_la(i)*v2s2(j,k))
            r_O2(j,k) = r_O2(j,k) + d_la(i) * EXP(-e_la(i)*v2s2(j,k))
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    !-------------------------------------------------------------------------

    j_loop: DO j = 1,kproma_day
      ! correction for zenith angles > sza_thr
      ! ref2642 from Lamago et al. (ACP 3, 1981-1990, 2003)
      fj_corr(j,:) = 1.
      zen = ACOS(cossza_1d(iu0(j)))*180./pi
      IF (zen>sza_thr) THEN
        fj_corr(j,1) = EXP(    1.2*((82.-zen)/5.5+1.))
        fj_corr(j,2) = EXP(1.3*1.2*((82.-zen)/5.5+1.))
        fj_corr(j,3) = EXP(1.4*1.2*((82.-zen)/5.5+1.))
        fj_corr(j,4) = EXP(1.5*1.2*((82.-zen)/5.5+1.))
        fj_corr(j,5) = EXP(2.0*1.2*((82.-zen)/5.5+1.))
        fj_corr(j,6) = EXP(3.0*1.2*((82.-zen)/5.5+1.))
        fj_corr(j,7) = EXP(4.0*1.2*((82.-zen)/5.5+1.))
        fj_corr(j,8) = EXP(5.0*1.2*((82.-zen)/5.5+1.))
      ENDIF
    ENDDO j_loop

    !-------------------------------------------------------------------------

    ! UV heating rates (only far UV)

    IF (l_heating) THEN
      DO k = 1,klev
        DO j = 1,kproma_day

          ! O2
          sig_h_O2(0) = p2(                                      &
            a0_1_1_h_O2(i0(j,k))*dlv2(j,k)+a0_1_2_h_O2(i0(j,k)), &
            a0_2_1_h_O2(i0(j,k))*dlv2(j,k)+a0_2_2_h_O2(i0(j,k)), &
            a0_3_1_h_O2(i0(j,k))*dlv2(j,k)+a0_3_2_h_O2(i0(j,k)), &
            v3_du1(j,k))                                         &
            * p2(1.,c0_1_1_h_O2(i0(j,k))*dlv2(j,k)               &
            + c0_1_2_h_O2(i0(j,k)),                              &
            c0_2_1_h_O2(i0(j,k))*dlv2(j,k)                       &
            + c0_2_2_h_O2(i0(j,k)),tnorm_sr(j,k) )
          sig_h_O2(1) = p1(                                      &
            b1_1_h_O2(i1(j,k)),a1_1_h_O2(i1(j,k)),v3_du1(j,k) )  &
            + p1( b1_2_h_O2(i1(j,k)),                            &
            a1_2_h_O2(i1(j,k)),v3_du1(j,k) )                     &
            * v2s_m(j,k)
          dh_O2(j,k) = (sig_h_O2(0) * fint(j,k,0) +              &
            sig_h_O2(1) * fint(j,k,1) ) * relo2
          rh_O2_2d(iu0(j),k) = REAL(MAX(0.0, dh_O2(j,k) *fj_corr(j,8))) ! 5
          ! O3
          sig_h_o3(0) = p1(a0_1_1_h_o3(i0(j,k))*dlv2(j,k)        &
            + a0_1_2_h_o3(i0(j,k)),                              &
            a0_2_1_h_o3(i0(j,k))*dlv2(j,k)                       &
            + a0_2_2_h_o3(i0(j,k)),                              &
            v3_du1(j,k)) *                                       &
            p2(1.,c0_1_1_h_o3(i0(j,k))*dlv2(j,k)+                &
            c0_1_2_h_o3(i0(j,k)),                                &
            c0_2_1_h_o3(i0(j,k))*dlv2(j,k)+                      &
            c0_2_2_h_o3(i0(j,k)),tnorm_sr(j,k))
          sig_h_o3(1) = p1(b1_1_h_o3(i1(j,k))                    &
            ,a1_1_h_o3(i1(j,k)),v3_du1(j,k))+                    &
            p1(b1_2_h_o3(i1(j,k)),a1_2_h_o3(i1(j,k))             &
            ,v3_du1(j,k))*                                       &
            v2s_m(j,k)
          dh_o3(j,k)   =   (sig_h_o3(0)    * fint(j,k,0) +       &
            sig_h_o3(1)    * fint(j,k,1) )*relo3(j,k)*1.e6
          rh_o3_2d(iu0(j),k) = REAL(MAX(0.0, dh_o3(j,k) *fj_corr(j,6))) ! 3
        ENDDO
      ENDDO
    ENDIF

    !-------------------------------------------------------------------------

    !MSK 06.11.2020
    !CALL jval_cal_uv
    IF (l_heating)  CALL jval_cal_uv
    !MSK END

    CALL jval_cal

    DEALLOCATE(v2s_m)
    DEALLOCATE(v3_du1)
    DEALLOCATE(v3_du2)
    DEALLOCATE(tnorm_sr)
    DEALLOCATE(tnorm)
    DEALLOCATE(dlv2)
    DEALLOCATE(dens)
    DEALLOCATE(r_m)
    DEALLOCATE(r_O2)

    DEALLOCATE(fint)
    DEALLOCATE(finth)
    DEALLOCATE(fj_corr)

    DEALLOCATE(iu0)
    DEALLOCATE(i0)
    DEALLOCATE(i1)
    DEALLOCATE(i2)
    DEALLOCATE(i3)

    DEALLOCATE(sig_top_O2)
    DEALLOCATE(sig_top_o3)
    DEALLOCATE(temp)
    DEALLOCATE(press)

    !-------------------------------------------------------------------------

  CONTAINS

    !-------------------------------------------------------------------------

    SUBROUTINE column_cal

      ! calculate O2 and O3 column densities

      INTEGER :: j, k
      REAL, PARAMETER :: sp = N_A/(M_air_SI*g)*1.E-4 ! [part./cm^2 * 1/Pa]

      DO j = 1,kproma_day
        v3s(j,:) = v3(j,:) / u0(j)
      ENDDO
      ! differential O3 columns and slant columns
      DO k = 1,klev
        DO j = 1,kproma_day
          dv3(j,k)  = v3(j,k) - v3(j,k-1)
        ENDDO
      ENDDO

      ! O2
      DO j = 1,kproma_day
        v2(j,:)  = sp * relo2 * press(j,:)
        v2s(j,:) = v2(j,:) / u0(j)
        ! differential columns and slant columns
        DO k = 1,klev
          dv2(j,k) = v2(j,k) - v2(j,k-1)
        ENDDO
      ENDDO

    END SUBROUTINE column_cal

    !-------------------------------------------------------------------------

    SUBROUTINE flux_cal

      ! run practical improved flux method (Zdunkwoski)
      ! to calculate actinic fluxes

      ! reference temperatures
      REAL, PARAMETER :: t1 = 226.
      REAL, PARAMETER :: t2 = 263.

      REAL, PARAMETER :: cr_O2(MAXWAV) = & ! O2 cross section
        (/ 7.6300E-24, 0.0000E+00, 0.0000E+00, 0.0000E+00, &
        0.0000E+00, 0.0000E+00, 0.0000E+00/)
      REAL, PARAMETER :: cr_o3(MAXWAV) = & ! O3 cross section (T = 298 K)
        (/ 3.6500E-19, 1.7900E-18, 2.7700E-19, 1.0500E-19, &
        2.6000E-20, 0.0000E+00, 4.5500E-21/)
      REAL, PARAMETER :: x1_o3(MAXWAV) = & ! T-dependence of a0_o3
        (/-2.7027E-23, 5.4054E-22, 0.0000E+00, 0.0000E+00, &
        0.0000E+00, 0.0000E+00, 0.0000E+00/)
      REAL, PARAMETER :: x2_o3(MAXWAV) = & ! T-dependence of a0_o3
        (/-4.1828E-25, 8.3655E-24, 0.0000E+00, 0.0000E+00, &
        0.0000E+00, 0.0000E+00, 0.0000E+00/)

      CALL slingo    ! cloud parameter (slingo parameterization)

      CALL aero_2d    ! aerosol parameter

      ! correction for absorption above model toa
      DO l = 1,MAXWAV
        j_loop1: DO j = 1,kproma_day
          ta_O2 = cr_O2(l) * v2(j,0)
          ta_o3 = (((temp(j,0)-t2)*x2_o3(l) +x1_o3(l))*(temp(j,0)-t1) &
            +cr_o3(l)) * v3(j,0)
          ftop(j,l) = EXP(-1./u0(j)*(ta_O2 + ta_o3))
        ENDDO j_loop1
      ENDDO

      ! optical depths of each layer
      ! mz_ht_20151124+
      ! SAVE PROPERTIES OF INTERNAL AEROSOL FOR OUTPUT AND COMPARISON ...
      IF (l_export_to_smil) THEN   ! op_pj_20160825
         DO l = 1,MAXWAV
            k_loop1b: DO k  = 1,klev
               j_loop2b: DO j = 1,kproma_day
                  zjv_asca(iu0(j),k,l) = taer_sca(j,k,l)
                  zjv_aabs(iu0(j),k,l) = taer_abs(j,k,l)
                  zjv_ga(iu0(j),k,l)   = gaer(j,k,l)
               END DO j_loop2b
            END DO k_loop1b
         END DO
      END IF ! op_pj_20160825
      ! ... THEN OVERWRITE THEM ...
      IF (l_aero_inp) THEN
         taer_sca = 0._dp
         taer_abs = 0._dp
         gaer     = 0._dp
         DO l = 1,MAXWAV
            k_loop1a: DO k  = 1,klev
               j_loop2a: DO j = 1,kproma_day
                  taer_sca(j,k,l) = zaer_sca(iu0(j),k,l)
                  taer_abs(j,k,l) = zaer_abs(iu0(j),k,l)
                  gaer(j,k,l)     = zaer_ga(iu0(j),k,l)
               END DO j_loop2a
            END DO k_loop1a
         END DO
      ENDIF
      ! mz_ht_20151124-

      DO l = 1,MAXWAV
        k_loop1: DO k  = 1,klev
          j_loop2: DO j = 1,kproma_day

            temp_avg  = 0.5 * (temp(j,k-1) + temp(j,k))

            ta_O2     = cr_O2(l) * dv2(j,k)
            ta_o3     = (((temp_avg-t2)*x2_o3(l) +x1_o3(l)) &
              *(temp_avg-t1)+cr_o3(l)) * dv3(j,k)
            taua_clr(j,k,l) = ta_o3 + ta_O2 + taer_abs(j,k,l)
            IF (k>1) THEN
              taus_clr(j,k,l) = crray(l)*1./0.21*dv2(j,k) + &
                taer_sca(j,k,l)
              g_clr(j,k,l)   = gaer(j,k,l)*taer_sca(j,k,l) / &
                taus_clr(j,k,l)
            ELSE
              taus_clr(j,k,l) = 0.
              g_clr(j,k,l)   = 0.
            ENDIF

            ! for calculation of g_clr see eg.: Radiation and cloud processes
            ! in the atmosphere, K.N. Liou, p.155 eq. 3.8.7

            taua_cld(j,k,l) = tau_abs(j,k) ! slingo option
            taus_cld(j,k,l) = tau_sca(j,k) ! slingo option
            g_cld(j,k,l)    = gcld(j,k)    ! slingo option

          ENDDO j_loop2
        ENDDO k_loop1
      ENDDO

      CALL pifm

      k_loop2: DO k = 0,klev
        j_loop3: DO j = 1,kproma_day
          fact(j,k,:) = fact(j,k,:) * ftop(j,:)
          facth(j,k,3:5) = facth(j,k,3:5) * ftop(j,3:5)
        ENDDO j_loop3
      ENDDO k_loop2

    END SUBROUTINE flux_cal

    !-------------------------------------------------------------------------

    SUBROUTINE slingo

      ! A. Slingo's data for cloud particle radiative properties (from 'a gcm
      ! parameterization for the shortwave properties of water clouds' jas
      ! vol. 46 may 1989 pp 1419-1427)
      ! here only for the spectral range 250nm - 690 nm

      REAL, PARAMETER :: &
        abar =  2.817E-02, & ! a coefficient for extinction optical depth
        bbar =  1.305,     & ! b coefficient for extinction optical depth
        cbar = -5.62E-08,  & ! c coefficient for single particle scat albedo
        dbar =  1.63E-07,  & ! d coefficient for single particle scat albedo
        ebar =  0.829,     & ! e coefficient for asymmetry parameter
        fbar =  2.482E-03, & ! f coefficient for asymmetry parameter
        cer  =  10.

      REAL :: &
        tau, & ! total optical depth
        wc,  & ! single scattering albedo
        cwc    ! co single scattering albedo 1-wc

      ! set cloud extinction optical depth, single scatter albedo,
      ! asymmetry parameter, and forward scattered fraction:
      ! do not let single scatter albedo be 1; delta-eddington solution
      ! for non-conservative case:
      wc   = MIN(1.-cbar-dbar*cer, 0.999999)
      cwc  = cbar + dbar * cer
      DO k = 1,klev
        DO j = 1,kproma_day
          tau           = clp(j,k) * (abar + bbar / cer)
          gcld(j,k)     = ebar + fbar * cer
          tau_abs(j,k)  = cwc * tau
          tau_sca(j,k)  = wc * tau
        ENDDO
      ENDDO

    END SUBROUTINE slingo

    !-------------------------------------------------------------------------

    SUBROUTINE aero_2d

      ! aerosol data from: Models for aerosols of the lower atmosphere
      ! and the effects of humidity variations on their optical properties
      ! E. P. Shettle and R. W. Fenn (1979)
      ! Environmental Research Paper No. 676

      ! relative humidity of lookup table [%]
      REAL, PARAMETER :: rh_ref(8) = (/ 0., 50., 70., 80., 90., 95., 98., 99./)

      INTEGER :: i_ref(kproma) ! interval index for relative humidity

      DO k = 1,klev

        DO j = 1,kproma_day
          ! DO i = 1,7
          ! IF (rh_ref(i)<= rhum(j,k)) i_ref(j) = i
          ! ENDDO
          IF (rh_ref(1)<= rhum(j,k)) i_ref(j) = 1
          IF (rh_ref(2)<= rhum(j,k)) i_ref(j) = 2
          IF (rh_ref(3)<= rhum(j,k)) i_ref(j) = 3
          IF (rh_ref(4)<= rhum(j,k)) i_ref(j) = 4
          IF (rh_ref(5)<= rhum(j,k)) i_ref(j) = 5
          IF (rh_ref(6)<= rhum(j,k)) i_ref(j) = 6
          IF (rh_ref(7)<= rhum(j,k)) i_ref(j) = 7
        ENDDO

        DO l = 1,MAXWAV
          DO j = 1,kproma_day

            domin = 1./(rh_ref(i_ref(j))-rh_ref(i_ref(j)+1))

            a1_sc = (asca(l,i_ref(j),1)-asca(l,i_ref(j)+1,1)) * domin
            a2_sc = (asca(l,i_ref(j),2)-asca(l,i_ref(j)+1,2)) * domin
            a4_sc = (asca(l,i_ref(j),4)-asca(l,i_ref(j)+1,4)) * domin

            b1_sc = asca(l,i_ref(j),1) - a1_sc*rh_ref(i_ref(j))
            b2_sc = asca(l,i_ref(j),2) - a2_sc*rh_ref(i_ref(j))
            b4_sc = asca(l,i_ref(j),4) - a4_sc*rh_ref(i_ref(j))

            a1_ab = (aabs(l,i_ref(j),1)-aabs(l,i_ref(j)+1,1)) * domin
            a2_ab = (aabs(l,i_ref(j),2)-aabs(l,i_ref(j)+1,2)) * domin
            a4_ab = (aabs(l,i_ref(j),4)-aabs(l,i_ref(j)+1,4)) * domin

            b1_ab = aabs(l,i_ref(j),1) - a1_ab*rh_ref(i_ref(j))
            b2_ab = aabs(l,i_ref(j),2) - a2_ab*rh_ref(i_ref(j))
            b4_ab = aabs(l,i_ref(j),4) - a4_ab*rh_ref(i_ref(j))

            a1_g = (ag(l,i_ref(j),1)-ag(l,i_ref(j)+1,1)) * domin
            a2_g = (ag(l,i_ref(j),2)-ag(l,i_ref(j)+1,2)) * domin
            a4_g = (ag(l,i_ref(j),4)-ag(l,i_ref(j)+1,4)) * domin

            b1_g = ag(l,i_ref(j),1) - a1_g*rh_ref(i_ref(j))
            b2_g = ag(l,i_ref(j),2) - a2_g*rh_ref(i_ref(j))
            b4_g = ag(l,i_ref(j),4) - a4_g*rh_ref(i_ref(j))

            bsa = 0.
            baa = 0.
            ga  = 0.

            IF (iaer(k) == 1) THEN
              bsa = (a1_sc*rhum(j,k) + b1_sc) * slf(j) + &
                (a2_sc*rhum(j,k) + b2_sc) * (1.-slf(j))
              baa = (a1_ab*rhum(j,k) + b1_ab) * slf(j) + &
                (a2_ab*rhum(j,k) + b2_ab) * (1.-slf(j))
              ga  = (a1_g *rhum(j,k) + b1_g)  * slf(j) + &
                (a2_g *rhum(j,k) + b2_g)  * (1.-slf(j))
            ENDIF
            IF (iaer(k) == 2) THEN
              bsa = (a4_sc*rhum(j,k) + b4_sc)
              baa = (a4_ab*rhum(j,k) + b4_ab)
              ga  = (a4_g *rhum(j,k) + b4_g)
            ENDIF

            taer_sca(j,k,l) = bsa * part_col(j,k)
            taer_abs(j,k,l) = baa * part_col(j,k)
            gaer(j,k,l) = ga
          ENDDO
        ENDDO

      ENDDO

    END SUBROUTINE aero_2d

    !-------------------------------------------------------------------------

    SUBROUTINE pifm

      ! practical improved flux method (pifm)
      ! to calculate actinic fluxes
      ! Zdunkowski,Welch,Korb: Beitr. Phys. Atmosph. vol. 53, p. 147 ff

      ! This version is not suitable for calculation for conserving
      ! scattering (w0 = 1). w0 is limited to w0 <= 1. - 1.E-15.
      ! for w0 = 1, al(4) and al(5) have to be calculated differently.

      INTRINSIC :: ABS

      ! arrays for matrix inversion
      REAL, DIMENSION(kproma,klev) :: &
        tu1, tu2, tu3, tu4, tu5, tu6, tu7, tu8, tu9

      REAL :: &
        al(kproma,2*klev,5),   & ! matrix coefficient
        rw(kproma,3*(klev+1)), & ! flux array for cloudy sky
        rf(kproma,3*(klev+1)), & ! flux array for clear sky
        sd,                    & ! parallel solar flux
        fd,                    & ! downward diffuse flux
        fu                       ! upward diffuse flux

      INTEGER :: nlev3p1,k3
      REAL :: f,tautot,w0,p_1,smooth1,smooth2
      REAL :: b0,bu0,alph1,alph2,alph3,alph4
      REAL :: eps,factor,ueps2,gam1,gam2,e,rm,tauscat,gg,ha,hb,hc,hd,gb,gc,gd
      REAL :: td1,td2,td3,td4,td5,td6,td7,tds1,tds2,tds3,tus1

      REAL, PARAMETER :: u = 2.      ! diffusivity factor
      REAL, PARAMETER :: delu0  = 1.E-3
      REAL, PARAMETER :: resonc = 1.E-6
      REAL, PARAMETER :: w0min  = 1.E-7

      !-----------------------------------------------------------------------

      nlev3p1 = 3 * (klev + 1)

      ! initiatization of pifm

      CALL pifmini

      ! option: maximum overlap

      ! calculation of the matrix coefficients a1,...,a5

      DO l  = 1,MAXWAV      ! wavel. loop

        ! first: clear sky

        DO k  = 1,klev      ! altitude loop
          DO j  = 1,kproma_day         ! longitude loop

            tautot = taus_clr(j,k,l)+taua_clr(j,k,l)

            ! single scattering albedo

            IF (tautot > 0.) THEN
              w0 = taus_clr(j,k,l) / tautot
            ELSE
              w0 = 0.
            ENDIF

            w0 = MIN(w0,1.-w0min)

            ! p_1: first expansion coefficient of the phase function
            p_1 = 3.*g_clr(j,k,l)

            ! f: fraction of radiation contained in diffraction peak
            f = g_clr(j,k,l)**2

            ! b0:  fractional mean backward scattering coefficient
            ! of diffuse light
            ! bu0: backward scattering coefficient of primary scattered
            ! parallel solar light
            ! for small p_1 smooth1,smooth2 manage the smooth change of
            ! b0 and bu0 to 0

            IF (p_1<= 0.1) THEN
              smooth1 = 1.33333333-p_1*3.3333333
              smooth2 = 10.*p_1
            ELSE
              smooth1 = 1.
              smooth2 = 1.
            ENDIF

            b0 = (3.-p_1)/8.  *smooth1
            bu0 = 0.5-u0(j)/4.*(p_1-3.*f)/(1.-f)  *smooth2

            ! alpha coefficient
            alph1 = u*(1.-(1.-b0)*w0)
            alph2 = u*b0*w0
            alph3 = w0*bu0*(1.-f)
            alph4 = w0*(1.-bu0)*(1.-f)

            ! epsilon and gamma coefficient
            eps = SQRT(alph1**2-alph2**2)
            factor = 1.-w0*f

            ! check for resonance condition in gam1 and gam2, if fulfil then
            ! chance u0(j) and calculate ueps2, bu0, alph3, alph4 again.
            ueps2 = (u0(j)*eps)**2
            IF (ABS(ueps2-factor**2)<resonc) THEN
              IF (ueps2<factor**2) THEN
                u0(j) = u0(j)-delu0
              ELSE
                u0(j) = u0(j)+delu0
              ENDIF
              ueps2 = (u0(j)*eps)**2
              bu0 = 0.5-u0(j)/4.*(p_1-3.*f)/(1.-f)  *smooth2
              alph3 = w0*bu0*(1.-f)
              alph4 = w0*(1.-bu0)*(1.-f)
            ENDIF

            gam1 = ( factor*alph3-u0(j)*(alph1*alph3+alph2*alph4) ) / &
              (factor**2-ueps2)
            gam2 = (-factor*alph4-u0(j)*(alph1*alph4+alph2*alph3) ) / &
              (factor**2-ueps2)

            e = EXP(-eps*tautot)
            rm = alph2/(alph1+eps)

            al(j,k,4) = e*(1.-rm**2)/(1.-e**2 * rm**2)
            al(j,k,5) = rm*(1.-e**2)/(1.-e**2 * rm**2)
            al(j,k,1) = EXP(-factor*tautot/u0(j))
            al(j,k,2) = -al(j,k,4)*gam2-al(j,k,5)*gam1*al(j,k,1) +gam2*al(j,k,1)
            al(j,k,3) = -al(j,k,5)*gam2-al(j,k,4)*gam1*al(j,k,1)+gam1

          ENDDO
        ENDDO

        ! second: cloudy sky

        DO k  = klev+1,2*klev ! altitude loop

          DO j = 1,kproma_day   ! longitude loop

            al(j,k,1) = al(j,k-klev,1)
            al(j,k,2) = al(j,k-klev,2)
            al(j,k,3) = al(j,k-klev,3)
            al(j,k,4) = al(j,k-klev,4)
            al(j,k,5) = al(j,k-klev,5)

          ENDDO

          DO j = 1,kproma_day    ! longitude loop

            IF (aclc(j,k-klev) <= fraclim) CYCLE

            tauscat = taus_clr(j,k-klev,l) + taus_cld(j,k-klev,l)
            tautot  = taua_clr(j,k-klev,l) + taua_cld(j,k-klev,l) &
              + tauscat
            gg      = g_cld(j,k-klev,l)*taus_cld(j,k-klev,l) / &
              tauscat

            IF (tautot > 0.) THEN
              w0 = tauscat/tautot
            ELSE
              w0 = 0.
            ENDIF

            w0 = MIN(w0,1.-w0min)

            p_1 = 3.*gg
            f = gg**2

            IF (p_1<= 0.1) THEN
              smooth1 = 1.33333333-p_1*3.3333333
              smooth2 = 10.*p_1
            ELSE
              smooth1 = 1.
              smooth2 = 1.
            ENDIF

            b0 = (3.-p_1)/8.  *smooth1
            bu0 = 0.5-u0(j)/4.*(p_1-3.*f)/(1.-f)  *smooth2

            alph1 = u*(1.-(1.-b0)*w0)
            alph2 = u*b0*w0
            alph3 = w0*bu0*(1.-f)
            alph4 = w0*(1.-bu0)*(1.-f)

            eps = SQRT(alph1**2-alph2**2)
            factor = 1.-w0*f

            ueps2 = (u0(j)*eps)**2
            IF (ABS(ueps2-factor**2)<resonc) THEN
              IF (ueps2<factor**2) THEN
                u0(j) = u0(j)-delu0
              ELSE
                u0(j) = u0(j)+delu0
              ENDIF
              ueps2 = (u0(j)*eps)**2
              bu0 = 0.5-u0(j)/4.*(p_1-3.*f)/(1.-f)  *smooth2
              alph3 = w0*bu0*(1.-f)
              alph4 = w0*(1.-bu0)*(1.-f)
            ENDIF

            gam1 = ( factor*alph3-u0(j)*(alph1*alph3+alph2*alph4)) / &
              (factor**2-ueps2)
            gam2 = (-factor*alph4-u0(j)*(alph1*alph4+alph2*alph3)) / &
              (factor**2-ueps2)

            e = EXP(-eps*tautot)
            rm = alph2/(alph1+eps)

            al(j,k,4) = e*(1.-rm**2)/(1.-e**2 * rm**2)
            al(j,k,5) = rm*(1.-e**2)/(1.-e**2 * rm**2)
            al(j,k,1) = EXP(-factor*tautot/u0(j))
            al(j,k,2) = -al(j,k,4)*gam2-al(j,k,5)*gam1*al(j,k,1) + gam2*al(j,k,1)
            al(j,k,3) = -al(j,k,5)*gam2-al(j,k,4)*gam1*al(j,k,1) + gam1

          ENDDO
        ENDDO

        ! matrix inversion

        DO j  = 1,kproma_day ! longitude loop

          ! direct solution of the first four equations

          rf(j,1) = u0(j)*flux(l)
          rw(j,1) = 0.
          rf(j,2) = 0.
          rw(j,2) = 0.

          ! 5th to 10th equation: bring matrix elements on the left of the main
          ! diagonal to the rhs:  save elements on the right of the main
          ! diagonal in array -tu(l,1)

          rf(j,3) = al(j,1,3) * bb(j,1,1) * rf(j,1)
          rf(j,4) = al(j,1,1) * bb(j,1,1) * rf(j,1)
          rf(j,5) = al(j,1,2) * bb(j,1,1) * rf(j,1)
          rw(j,3) = al(j,1+klev,3) * (1.-bb(j,1,1)) * rf(j,1)
          rw(j,4) = al(j,1+klev,1) * (1.-bb(j,1,1)) * rf(j,1)
          rw(j,5) = al(j,1+klev,2) * (1.-bb(j,1,1)) * rf(j,1)

          tu1(j,1) = 0.
          tu2(j,1) = al(j,1,4)      * bb(j,1,2)
          tu3(j,1) = al(j,1,4)      * (1.-bb(j,1,4))
          tu4(j,1) = al(j,1+klev,4) * (1.-bb(j,1,2))
          tu5(j,1) = al(j,1+klev,4) * bb(j,1,4)
          tu6(j,1) = al(j,1,5)      * bb(j,1,2)
          tu7(j,1) = al(j,1,5)      * (1.-bb(j,1,4))
          tu8(j,1) = al(j,1+klev,5) * (1.-bb(j,1,2))
          tu9(j,1) = al(j,1+klev,5) * bb(j,1,4)


        ENDDO

        ! blocks of 6 equations: eliminate left matrix elements, save right
        ! matrix elements in array -tu(l,i), calculate rhs.

        DO k = 2,klev
          k3 = 3*k
          DO j  = 1,kproma_day ! longitude loop

            ha  = bb(j,k,1)*tu6(j,k-1) + (1.-bb(j,k,3))*tu8(j,k-1)
            hb  = bb(j,k,1)*tu7(j,k-1) + (1.-bb(j,k,3))*tu9(j,k-1)
            hc  = (1.-bb(j,k,1))*tu6(j,k-1) + bb(j,k,3)*tu8(j,k-1)
            hd  = (1.-bb(j,k,1))*tu7(j,k-1) + bb(j,k,3)*tu9(j,k-1)
            ga  = bb(j,k,1)*rf(j,k3-2) + (1.-bb(j,k,3))*rw(j,k3-2)
            gb  = bb(j,k,1)*rf(j,k3-1) + (1.-bb(j,k,3))*rw(j,k3-1)
            gc  = (1.-bb(j,k,1))*rf(j,k3-2) + bb(j,k,3)*rw(j,k3-2)
            gd  = (1.-bb(j,k,1))*rf(j,k3-1) + bb(j,k,3)*rw(j,k3-1)

            td1        = 1./(1.-al(j,k,5)*ha)
            rf(j,k3)   = td1*(al(j,k,3)*ga + al(j,k,5)*gb)
            tu1(j,k)   = td1*al(j,k,5)*hb
            tu2(j,k)   = td1*al(j,k,4)*bb(j,k,2)
            tu3(j,k)   = td1*al(j,k,4)*(1.-bb(j,k,4))
            td2        = al(j,k+klev,5)*hc
            td3        = 1./(1.-al(j,k+klev,5)*hd-td2*tu1(j,k))
            rw(j,k3)   = td3*(al(j,k+klev,3)*gc + &
                         al(j,k+klev,5)*gd+td2*rf(j,k3))
            tu4(j,k)   = td3*(al(j,k+klev,4)*(1.-bb(j,k,2)) + &
                         td2*tu2(j,k))
            tu5(j,k)   = td3*(al(j,k+klev,4)*bb(j,k,4) + &
                         td2*tu3(j,k))
            rf(j,k3+1) = al(j,k,1)*ga
            rw(j,k3+1) = al(j,k+klev,1)*gc
            td4        = al(j,k,4)*ha
            td5        = al(j,k,4)*hb
            rf(j,k3+2) = al(j,k,2)*ga +al(j,k,4)*gb+td4* &
                         rf(j,k3) + td5*rw(j,k3)
            tu6(j,k)   = al(j,k,5)*bb(j,k,2) + td4*tu2(j,k) + &
                         td5* tu4(j,k)
            tu7(j,k)   = al(j,k,5)*(1.-bb(j,k,4))+td4*tu3(j,k)+ &
                         td5*tu5(j,k)
            td6        = al(j,k+klev,4)*hc
            td7        = al(j,k+klev,4)*hd
            tu8(j,k)   = al(j,k+klev,5)*(1.-bb(j,k,2)) + &
                         td6*tu2(j,k) + td7*tu4(j,k)
            tu9(j,k)   = al(j,k+klev,5)*bb(j,k,4) + td6*tu3(j,k)+ &
                         td7*tu5(j,k)
            rw(j,k3+2) = al(j,k+klev,2)*gc + al(j,k+klev,4)*gd + &
                         td6*rf(j,k3) + td7*rw(j,k3)

          ENDDO
        ENDDO

        ! last two equations: the same as before

        DO j = 1,kproma_day ! longitude loop
          tds1          = 1. / (1.-albedo(j,l)*tu6(j,klev))
          rf(j,nlev3p1) = tds1 * (albedo(j,l) * rf(j,nlev3p1-2)+ &
                          albedo(j,l) * rf(j,nlev3p1-1))
          tus1          = tds1*albedo(j,l)*tu7(j,klev)
          tds2          = albedo(j,l)*tu8(j,klev)
          tds3          = 1./(1.-albedo(j,l)*tu9(j,klev)-tds2*tus1)
          rw(j,nlev3p1) = tds3*(albedo(j,l)*rw(j,nlev3p1-2) + &
                          albedo(j,l)*rw(j,nlev3p1-1) + tds2*rf(j,nlev3p1))
          rf(j,nlev3p1) = rf(j,nlev3p1) + tus1*rw(j,nlev3p1)
        ENDDO

        ! now we have created an upper triangular matrix the elements of which
        ! are -tu(l,i), 0, or 1 (in the main diagonal). the 0 and 1 elements
        ! are not stored in an array. let us solve the system now and store the
        ! results in the arrays rf (fluxes clear sky) and rw (fluxes cloudy sky)

        DO k = klev,1,-1
          k3 = 3*k
          DO j  = 1,kproma_day ! longitude loop

            rw(j,k3+2) = rw(j,k3+2) + tu8(j,k)*rf(j,k3+3) + &
                         tu9(j,k)*rw(j,k3+3)
            rf(j,k3+2) = rf(j,k3+2) + tu6(j,k)*rf(j,k3+3) + &
                         tu7(j,k)*rw(j,k3+3)
            rw(j,k3)   = rw(j,k3) + tu4(j,k)*rf(j,k3+3)   + &
                         tu5(j,k)*rw(j,k3+3)
            rf(j,k3)   = rf(j,k3) + tu2(j,k)*rf(j,k3+3)   + &
                         tu3(j,k)*rw(j,k3+3) + tu1(j,k)*rw(j,k3)

            sd = rf(j,k3+1) + rw(j,k3+1)
            fd = rf(j,k3+2) + rw(j,k3+2)
            fu = rf(j,k3+3) + rw(j,k3+3)

            ! actinic flux shall not be caculated at toa
            fact(j,k,l) = MAX(0.e0 , sd/u0(j) + u * fd + u * fu)
            IF (l>2 .AND. l<6) facth(j,k,l) = MAX(0.e0 ,sd+fd)
          ENDDO
        ENDDO

        DO j  = 1,kproma_day ! longitude loop
          sd = rf(j,1) + rw(j,1)
          fd = rf(j,2) + rw(j,2)
          fu = rf(j,3) + rw(j,3)
          fact(j,0,l) = MAX(0.e0 , sd/u0(j) + u * fd + u * fu)
          IF (l>2 .AND. l<6) facth(j,0,l) = MAX(0.e0 ,sd+fd)
        ENDDO

      ENDDO

    END SUBROUTINE pifm

    !-------------------------------------------------------------------------

    SUBROUTINE pifmini

      ! initialization of practical improved flux method
      ! Zdunkowski et al 1982
      ! and
      ! J.F. Geleyn and A. Hollingsworth Contrib. to Atm. Phys. 52 no.1 p.1
      !-----------------------------------------------------------------------
      ! optical parameters of clouds and aerosol
      ! calculate continuity matrices for fractional cloudiness which find the
      ! flux going through a level i; levels are numbered from top to bottom.
      ! the following notation is used
      !
      ! s(cld,b) = bb(1)*s(clr,t)      + bb(3)*s(cld,t)
      ! s(clr,b) = (1-bb(1))*s(clr,t)  + (1-bb(3))*s(cld,t)
      ! fu(cld,t) = bb(2)*fu(clr,b)     + bb(4)*fu(cld,b)
      ! fu(fr,t) = (1-bb(2))*fu(clr,b) + (1-bb(4))*fu(cld,b)
      !
      ! where
      !
      ! cld = cloudy
      ! clr = clear
      ! b   = bottom of level k
      ! t   = top of level k
      ! s   = direct or diffuse downward flux
      ! fu  = diffuse upward flux

      ! first for highest level k = 1
      DO j  = 1,kproma_day ! longitude loop
        bb(j,1,1) = 1.-aclc(j,1)
        bb(j,1,2) = 1.
        bb(j,1,3) = 1.
        bb(j,1,4) = 1.
      ENDDO

      ! now for levels 2 to klev-1
      DO k = 2,klev-1
        km = k-1
        kp = k+1

        DO j  = 1,kproma_day ! longitude loop

          IF (aclc(j,km)<1.) THEN
            bb(j,k,1) = (1. - MAX(aclc(j,k),aclc(j,km))) / (1. - aclc(j,km))
          ELSE
            bb(j,k,1) = 1.
          ENDIF

          IF (aclc(j,km)>0.) THEN
            bb(j,k,3) = MIN(aclc(j,k),aclc(j,km)) / aclc(j,km)
          ELSE
            bb(j,k,3) = 1.
          ENDIF

          IF (aclc(j,kp)<1.) THEN
            bb(j,k,2) = (1. - MAX(aclc(j,k),aclc(j,kp))) / (1. - aclc(j,kp))
          ELSE
            bb(j,k,2) = 1.
          ENDIF

          IF (aclc(j,kp)>0.) THEN
            bb(j,k,4) =  MIN(aclc(j,k),aclc(j,kp)) / aclc(j,kp)
          ELSE
            bb(j,k,4) = 1.
          ENDIF

        ENDDO
      ENDDO

      ! finally for lowest level k = klev
      k  = klev
      km = klev-1

      DO j  = 1,kproma_day ! longitude loop

        IF (aclc(j,km)<1.) THEN
          bb(j,k,1) = (1. - MAX(aclc(j,k),aclc(j,km))) / (1. - aclc(j,km))
        ELSE
          bb(j,k,1) = 1.
        ENDIF

        IF (aclc(j,km)>0.) THEN
          bb(j,k,3) = MIN(aclc(j,k),aclc(j,km)) / aclc(j,km)
        ELSE
          bb(j,k,3) = 1.
        ENDIF

        bb(j,klev,2) = 1.
        bb(j,klev,4) = 1.

      ENDDO

    END SUBROUTINE pifmini

    !-------------------------------------------------------------------------

  END SUBROUTINE jvalues

  ! **************************************************************************

  SUBROUTINE jval_cal_uv

    INTEGER :: j, k
    REAL    :: dj,djb
    REAL, DIMENSION(3:5) :: coeff_uv, coeff_dna
    REAL, PARAMETER :: a3_uv(3) = (/6.5612E-19, -2.3893E-24,  2.6419E-28/)
    REAL, PARAMETER :: a4_uv(2) = (/6.3890E-19, -6.9346E-25/)
    REAL, PARAMETER :: a5_uv(4) = (/1.3218E-19,-6.9323E-23, &
         1.3036E-26,-8.3734E-31/)
    REAL, PARAMETER :: a3_dna(3) = (/8.3302E-20, -2.6503E-23, &
         3.0102E-27/)
    REAL, PARAMETER :: a4_dna(3) = (/7.8083E-21,-2.2999E-24, 2.3155E-28/)
    REAL, PARAMETER :: a5_dna(4) = (/1.3206E-22,-7.8328E-26, &
         1.6011E-29,-1.0860E-33/)

    fhuv_2d(:,:)    = 0.0_dp
    fhuvdna_2d(:,:) = 0.0_dp

    DO k = 1,klev
       DO j = 1,kproma_day
          coeff_uv(3)=p2(a3_uv(1),a3_uv(2),a3_uv(3),v3_du2(j,k))
          coeff_dna(3)=p2(a3_dna(1),a3_dna(2),a3_dna(3),v3_du2(j,k))
          coeff_uv(4)=p1(a4_uv(1),a4_uv(2),v3_du2(j,k))
          coeff_dna(4)=p2(a4_dna(1),a4_dna(2),a4_dna(3),v3_du2(j,k))
          coeff_uv(5)=p3(a5_uv(1),a5_uv(2),a5_uv(3),a5_uv(4),v3_du2(j,k))
          coeff_dna(5)=p3(a5_dna(1),a5_dna(2),a5_dna(3),a5_dna(4),v3_du2(j,k))
          dj=1.e4*(coeff_uv(3)*finth(j,k,3)+coeff_uv(4)*finth(j,k,4)+ &
               coeff_uv(5)*finth(j,k,5))
          djb=1.e4*(coeff_dna(3)*finth(j,k,3)+coeff_dna(4)*finth(j,k,4)+ &
               coeff_dna(5)*finth(j,k,5))
          fhuv_2d(iu0(j),k)=REAL(MAX(0.0,dj),dp)
          fhuvdna_2d(iu0(j),k)=REAL(MAX(0.0,djb),dp)
       ENDDO
    ENDDO

  END SUBROUTINE jval_cal_uv

  ! **************************************************************************

  ! include the subroutines jval_cal and jval_cal_* which are
  ! dynamically generated by jvpp:
#include "messy_jval_jvpp.inc"

  ! **************************************************************************

  FUNCTION p1(c0,c1,x)
    IMPLICIT NONE
    REAL :: p1,c0,c1,x
    p1 = c0 + c1*x
  END FUNCTION p1

  ! **************************************************************************

  FUNCTION p2(c0,c1,c2,x)
    IMPLICIT NONE
    REAL :: p2,c0,c1,c2,x
    p2 = c0 + (c1+c2*x)*x
  END FUNCTION p2

  ! **************************************************************************

  FUNCTION p3(c0,c1,c2,c3,x)
    IMPLICIT NONE
    REAL :: p3,c0,c1,c2,c3,x
    p3 = c0 + (c1+(c2+c3*x)*x)*x
  END FUNCTION p3

  ! **************************************************************************

  SUBROUTINE jval_read_nml_ctrl(status, iou)

    ! read CTRL namelist, check it, and initialize global switches/variables
    ! Author: Patrick Joeckel, (MPICH), Feb 2003

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    INTEGER, INTENT(OUT) :: status ! error status

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'jval_read_nml_ctrl'
    LOGICAL                      :: lex          ! file exists?
    INTEGER                      :: fstat        ! file status

    NAMELIST /CTRL/ r_sol, qy_CH3COCH3

    ! initialize
    status = 1 ! error

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML = CTRL, IOSTAT = fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! diagnose namelist and set global switches
    ! check namelist entries

    WRITE(*,*) HLINE2
    WRITE(*,*) 'SOLAR CYCLE PARAMETER r_sol = ',r_sol
    WRITE(*,*) 'NOTE: THIS IS POSSIBLY OBSOLETE DEPENDING ON CPL NAMELIST.'
    WRITE(*,*) HLINE2

    IF ((qy_CH3COCH3<1).OR.(qy_CH3COCH3>3)) THEN
      WRITE(*,*) HLINE2
      WRITE(*,*) 'ERROR: you must select: 1 <= qy_CH3COCH3 <= 3'
      WRITE(*,*) HLINE2
      RETURN
    ELSE
      WRITE(*,*) HLINE2
      WRITE(*,*) 'qy_CH3COCH3 = ', qy_CH3COCH3
      WRITE(*,*) HLINE2
    END IF
    ! end check namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! no error

  END SUBROUTINE jval_read_nml_ctrl

! ****************************************************************************

  SUBROUTINE jval_solar_time_control(status, cdisse, val)

    IMPLICIT NONE

    ! I/O
    INTEGER,                INTENT(OUT)          :: status
    REAL(DP),               INTENT(IN)           :: cdisse
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: val

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'jval_solar_time_control'

    status = 1 ! ERROR

    external_sol: IF (PRESENT(val)) THEN

      ! Fluxes must be adjusted to the Sun-Earth distance using the
      ! factor cdisse ("DIStance-Sun-Earth" in AU).
      ! f0 is a ratio of fluxes, cdisse cancels out.
      SELECT CASE(SIZE(val))
      CASE(1)
        ! The 10.7 cm solar radio flux (F10.7cm, in sfu = solar flux
        ! units) changes from about 70 sfu (1E-22 W m-2 Hz -1) to about
        ! 270 sfu during the solar cycle (see Fig. 1 in Tapping, Space
        ! Weather, 11, 394-406 (2013)). Thus, the current position in
        ! the solar cycle r_sol is:
        r_sol = (val(1) - 70.0) / 200.0
        flux(:) = (r_sol*fmax(:) + (1.-r_sol)*fmin(:))*cdisse
        SR_toa_flux = (r_sol * 9.6701E+12 + (1.-r_sol) * 8.7990E+12)*cdisse
        phi_la = (r_sol * 7.E+11     + (1.-r_sol) * 3.E+11)*cdisse
        f0(:) = r_sol*f0max + (1.-r_sol)*f0min
        !
      CASE(16)
        flux(:) = val(3:9)*cdisse
        SR_toa_flux =  val(2)*cdisse
        phi_la = val(1)*cdisse
        f0(:) = val(10:16)
        !
      CASE DEFAULT
        !
        WRITE(*,*) substr, ': ERROR IN CHOICE OF jval_solar (CPL)'//&
          &'! Number of parameters is ',SIZE(val),&
          ' but only 1 or 16 are currently implemented.'
        RETURN
        !
      END SELECT

    ELSE

      ! constant r_sol (CTRL namelist)
      flux(:) = (r_sol*fmax(:) + (1.-r_sol)*fmin(:))*cdisse
      SR_toa_flux = (r_sol * 9.6701E+12 + (1.-r_sol) * 8.7990E+12)*cdisse
      phi_la = (r_sol * 7.E+11     + (1.-r_sol) * 3.E+11)*cdisse
      ! ratio of fluxes; cdisse cancels out !
      f0(:) = r_sol*f0max + (1.-r_sol)*f0min

    ENDIF external_sol

    status = 0

  END SUBROUTINE jval_solar_time_control

! ****************************************************************************

END MODULE messy_jval

! ****************************************************************************
