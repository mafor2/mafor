! <gde_constants.f90 - A component of the Multicomponent
!                     Aerosol Dynamics Model MAFOR>
!*****************************************************************************! 
!* 
!*    Copyright (C) 2011-2021  Matthias Steffen Karl
!*
!*    Contact Information:
!*          Dr. Matthias Karl
!*          Sulzbrackring 13
!*          21037 Hamburg
!*          Germany
!*          email:  mattkar@googlemail.com
!*
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!*
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!*
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*
!*    The MAFOR code is intended for research and educational purposes. 
!*    Users preparing publications resulting from the usage of MAFOR are 
!*    requested to cite:
!*    1.  Karl, M., Gross, A., Pirjola, L., Leck, C., A new flexible
!*        multicomponent model for the study of aerosol dynamics
!*        in the marine boundary layer, Tellus B, 63(5),1001-1025,
!*        doi:10.1111/j.1600-0889.2011.00562.x, 2011.
!*    2.  Karl, M., Kukkonen, J., Keuken, M.P., Lutzenkirchen, S.,
!*        Pirjola, L., Hussein, T., Modelling and measurements of urban
!*        aerosol processes on the neighborhood scale in Rotterdam,
!*        Oslo and Helsinki, Atmos. Chem. Phys., 16,
!*        4817-4835, doi:10.5194/acp-16-4817-2016, 2016.
!*
!*****************************************************************************!
!*    All routines written by Matthias Karl
!* 
!*****************************************************************************!

    !----------------------------------------------------------------------!
    !                                                                      !
    !****                                                                  !
    ! Definitions of machine precision constants as                        !
    !                                  Fortran parameters for MAFOR        !
    ! Definitions of physical constants as Fortran parameters for MAFOR    !
    !                                                                      !
    !----------------------------------------------------------------------!


MODULE gde_constants

  IMPLICIT NONE
  INTRINSIC :: SELECTED_REAL_KIND
  ! PUBLIC is already default

! KPP DP - Double precision kind
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)

  ! Upper limit of name length of KPP species
  INTEGER, parameter :: STRLEN_KPPSPECIES =  15


  ! PHYSICAL CONSTANTS
  real(dp), parameter :: pi      = 3.14159265358979323846_dp
  real(dp), parameter :: R_gas   = 8.314409_dp  ! R [J/K/mol]
  ! Stephan-Boltzmann constant
  real(dp), parameter :: stbo    = 5.67E-8_dp   ! [W/m2/K4]
  real(dp), parameter :: N_A     = 6.022045E23_dp ! Avogadro constant [1/mol]
  real(dp), parameter :: g       = 9.80665_dp   ! gravity acceleration [m/s2]
  real(dp), parameter :: T0      = 298.15_dp    ! standard temperature [K]
  real(dp), parameter :: T0_INV  = 1._DP / T0   ! 1/T0 [1/K]
  real(dp), parameter :: atm2Pa  = 101325._dp   ! conversion from [atm] to [Pa]
  real(dp), parameter :: cal2J   = 4.1868_dp    ! conversion from [cal] to [J]
  real(dp), parameter :: k_B     = 1.380662E-23_dp ! Boltzmann constant [J/K]
  real(dp), parameter :: c_vKar  = 0.4_dp       !  Karman constant [?]
  
  ! MXXX = molar mass of element XXX [g/mol]
  real(dp), parameter :: MH  =   1.01_dp
  real(dp), parameter :: MC  =  12.01_dp
  real(dp), parameter :: MN  =  14.01_dp
  real(dp), parameter :: MF  =  19.00_dp
  real(dp), parameter :: MNa =  22.99_dp
  real(dp), parameter :: MO  =  16.00_dp
  real(dp), parameter :: MS  =  32.07_dp
  real(dp), parameter :: MCl =  35.45_dp
  real(dp), parameter :: MBr =  79.90_dp
  real(dp), parameter :: MI  = 126.90_dp
  real(dp), parameter :: MHg = 200.59_dp
  real(dp), parameter :: MCa =   40.0_dp
  real(dp), parameter :: MCO3=   60.0_dp
  real(dp), parameter :: MSi =   28.1_dp
  ! M_XXX = molar mass of compounds [g/mol]
  real(dp), parameter :: M_O3  = MO*3._dp      ! molar mass of ozone [g/mol]
  real(dp), parameter :: M_H2O = MH*2._dp + MO ! molar mass of H2O [g/mol]
  
  ! DRY AIR AND WATER VAPOUR THERMODYNAMIC CONSTANTS
  real(dp), parameter :: tmelt   = 273.15_dp    ! melting temp. of ice/snow [K]
  real(dp), parameter :: ttrip   = 273.16_dp    ! triple point of water [K]
  real(dp), parameter :: rho_H2O = 999.97_dp    ! density of H2O [kg/m3]
  real(dp), parameter :: M_air   = 28.970_dp    ! molar mass of dry air [g/mol]
  real(dp), parameter :: cp_air  = 1005.46_dp   ! specific heat of dry air at
                                                ! constant pressure [J/K/kg]
                                                
  ! mz_ap_20090519+
  real(dp), parameter :: alv   = 2.5008e6_dp    ! latent heat for vaporisation 
  !                                             ! [J/kg]
  real(dp), parameter :: als   = 2.8345e6_dp    ! latent heat for sublimation
  !                                             ! [J/kg]
  real(dp), parameter :: alf   = als-alv        ! latent heat for fusion [J/kg]
  
  ! mz_ap_20090519-

  ! gas constant for dry air [J/K/kg]
  real(dp), parameter :: rd      = 1000._dp * R_gas/M_air ! 287.05_dp
  ! gas constant for water vapour
  real(dp), parameter :: rv      = 1000._dp * R_gas/M_H2O ! 461.51_dp
  ! specific heat of water vapour at constant pressure [J/K/kg]
  real(dp), parameter :: cpv     = 1869.46_dp
  ! dimensionless auxiliary constants
  real(dp), parameter :: vtmpc1  = rv/rd-1.0_dp
  real(dp), parameter :: vtmpc2  = cpv/cp_air-1.0_dp
  real(dp), parameter :: MM_eps  = M_H2O/M_air ! mz_hr_20070323

  ! cloud and radiation
  real(dp), SAVE     :: ceffmin = 10.0_dp    ! min eff.radius for ice cloud
  real(dp),parameter :: ceffmax = 150.0_dp   ! max eff.radius for ice cloud
  real(dp), SAVE     :: ccwmin  = 1.0e-7_dp  ! cloud water limit for cover>0
  real(dp),parameter :: cemiss  = 0.996_dp   ! LW emissivity 

  ! PLANETARY parameterS
  real(dp), parameter :: radius_earth = 6371000.0_dp ! radius of the Earth [m]
  real(dp), parameter :: OneDay       = 86400.0_dp   ! one day [s]
  ! fu_kk_20061002+
  real(dp), parameter :: solc  = 1365.0_dp           ! solar constant [W/m2]
  !real(dp), parameter :: solc  = 1365.41_dp          ! solar constant [W/m2]
  ! fu_kk_20061002-
  ! *ratio: atmospheric height/radius of the earth.
  real(dp), parameter :: crae = 0.1277e-02_dp
  
  ! mz_ab_20090525+
  real(dp), parameter:: AM = 1.673e-27     ! Atomic mass unit
  real(dp), parameter:: ELCH =  1.602E-19  ! Electron charge

  real(dp), parameter:: TWOPI = pi*2._dp      ! Pi*2.
  real(dp), parameter:: PI_2  = pi*0.5_dp     ! Pi/2.
  real(dp), parameter:: DTR   = pi/180._dp    ! Degrees to radians
  real(dp), parameter:: RTD   = 180._dp/pi    ! Radians to degrees
  ! mz_ab_20090525-

  ! MAFOR PHYSICAL CONSTANTS
  real( dp), parameter :: Rv_h2o  = 0.4614       ! Rv (J/g*K)
  real( dp), parameter :: AVOng   = 1.e09_dp/N_A
  
  ! MAFOR CHEMISTRY CONSTANTS
  real( dp), parameter :: M_H2SO4 = 98.07             ! molar mass of H2SO4 [g/mol]
  real( dp), parameter :: MB      = 98.08*1.661e-27   ! sulphuric acid (kg/molec)
  real( dp), parameter :: MVOC    = 118.*1.661e-27    ! some organic (kg/molec)
  real( dp), parameter :: MAH     = 118.*1.661e-27    ! organic C4 dimer (kg/molec)
  real( dp), parameter :: MAN     = 63.*1.661e-27*2   ! nitric acid (kg/molec)
  real( dp), parameter :: MNH     = 18.*1.661e-27     ! ammonia (kg/molec)
  real( dp), parameter :: MAC     = 450.*1.661e-27    ! nucleating vapor (kg/molec)
  real( dp), parameter :: M_ABS   = 132.26            ! molar mass of ammonium bisulfate [g/mol] 
  real( dp), parameter :: M_msa   = 96.11
  real( dp), parameter :: M_nit   = 63.
  real( dp), parameter :: M_nh3   = 17.
  real( dp), parameter :: M_dma   = 45.08
  real( dp), parameter :: M_ca    = 450.00
  real( dp), parameter :: M_hcl   = 36.
  real( dp), parameter :: RHOH2O  = 999.9668          ! density of water [kg/m3]  
  
  ! MAFOR PLANETARY parameterS
  real( dp), parameter :: InitSpr = 80.               ! first day of spring [d]
  real( dp), parameter :: Cancer  = 23.45 * (pi/180.) ! angle of earth
  real( dp), parameter :: latitu  = 45. * (pi/180.)   ! latitude

END MODULE gde_constants

!*****************************************************************************
