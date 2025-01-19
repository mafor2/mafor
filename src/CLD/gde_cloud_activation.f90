! <gde_cloud_activation.f90 - A component of the Multicomponent
!                     Aerosol Dynamics Model MAFOR>
!*****************************************************************************! 
!* 
!*    Copyright (C) 2011-2024 Matthias Steffen Karl
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
!*    1.  Karl, M., Pirjola, L., Grönholm, T., Kurppa, M., Anand, S., 
!*        Zhang, X., Held, A., Sander, R., Dal Maso, M., Topping, D., 
!*        Jiang, S., Kangas, L., and Kukkonen, J., Description and 
!*        evaluation of the community aerosol dynamics model MAFOR v2.0,
!*        Geosci. Model Dev., 15, 
!*        3969-4026, doi:10.5194/gmd-15-3969-2022, 2022.
!*
!*!*****************************************************************************!
!*    All routines written by Matthias Karl
!* 
!*******************************************************************************!
!*
!*    References Cloud Module
!*    -----------------------
!*    H. Abdul-Razzak, S.J. Ghan, C. Rivera-Carpio,
!*        A parameterization of aerosol activation,
!*        1. Single aerosol type, J. Geophys. Res.,
!*        103, D6, 6123-6131, 1998.
!*    H. Abdul-Razzak and S.J. Ghan,
!*        A parameterization of aerosol activation,
!*        2. Multiple aerosol types, J. Geophys. Res.,
!*        105, D5, 6837-6844, 2000.
!*    H. Abdul-Razzak and S.J. Ghan,
!*        A parameterization of aerosol activation,
!*        3. Sectional representation, J. Geopphys. Res.,
!*        107, D3, 4026, doi:10.1029/2001JD000483, 2002.
!*    W.R. Leaitch, J.W. Strapp, G.A. Isaac,
!*        Cloud droplet nucleation and cloud scavenging
!*        of aerosol sulphate in polluted atmospheres,
!*        Tellus, 38B, 328-344, 1986.
!*    A. Kerkweg, S. Wurzler, T. Reisin, A. Bott,
!*        On the cloud processing of aerosol particles:
!*        An entraining air-parcel model with two-dimensional
!*        spectral cloud microphysics and a new formulation
!*        of the collection kernel,
!*        Q. J. R. Meteorolo. Soc., 129, 1-19, 2003
!*    J.H. Seinfeld and S.N. Pandis,
!*        Atmospheric Chemistry and Physics, From Air
!*        Pollution to Climate Change, 2nd Edition,
!*        John Wiley & Sons, Inc., Hoboken New Jersey, 2006.
!*    S.N Pandis, J.H. Seinfeld, C. Pilinis,
!*        Chemical composition differences in fog and
!*        cloud droplets of different sizes,
!*        Atmospheric Environment, 24A(7), 1957-1969, 1990a.
!*    S.N. Pandis, J.H. Seinfeld, C. Pilinis,
!*        The smog-fog-smog cycle and acid deposition
!*        J. Geophys. Res., 95, D11, 18,489-18,500, 1990b.
!*    W.A. Hoppel, J.W. Fitzgerald, G.M. Frick,
!*        R.E. Larson, Aerosol size distributions and
!*        optical properties found in the marine
!*        boundary layer over the Atlantic Ocean,
!*        J. Geophys. Res., 95, D4, 3659-3683, 1990.
!*    D. Topping, P. Connolly, G. McFiggans,
!*        Cloud droplet number enhanced by co-condensation
!*        of organic vapours, Supplementary Information,
!*        Nature Geoscience, 6, 443-446, 2013.
!*    A. Laaksonen, P. Korhonen, M. Kulmala, R.J. Charlson,
!*        Modification of the Koehler equation to include
!*        soluble trace gases and slightly soluble substances,
!*        J. Atmos. Sci., 55, 853-862, 1998.
!*
!*******************************************************************************!
module gde_cloud_activation

  use messy_mecca_kpp_Parameters
  use messy_mecca_kpp_global, only : APN
  use messy_mecca_kpp_global, only : ind_Clm_a, ind_Nap_a, ind_SO4mm_a
  use messy_mecca_kpp_global, only : ind_HSO4m_a, ind_NO3m_a, ind_CH3SO3m_a
  use messy_mecca_kpp_global, only : ind_DMAp_a, ind_NH4p_a, ind_DOC_a
  use messy_mecca_kpp_global, only : ind_HC2O4m_a, ind_C2O4mm_a
  use messy_mecca_kpp_global, only : ind_CH3COCOOm_a
  use messy_mecca_kpp_global, only : ind_C2H5C2O4m_a, ind_C2H4C2O4mm_a
  use messy_mecca_kpp_global, only : ind_so2_a, ind_h2o2_a, ind_hno3_a
  use messy_mecca_kpp_global, only : ind_IO3m_a
  use messy_mecca_kpp_global, only : ind_hcho_a, ind_SUCCAC_a, ind_ADIPAC_a
  use messy_mecca_kpp_global, only : ind_GLUTARAC_a, ind_OXALAC_a

  use gde_sensitiv,     only       : IDEPO
  use gde_sensitiv,     only       : IWETD
  use gde_sensitiv,     only       : ICOAG
  use gde_sensitiv,     only       : ICOND

  use gde_input_data,   only       : MMAX,AMAX
  use gde_input_data,   only       : NU,NA,AI,AS,CS
  use gde_input_data,   only       : aqmax
  use gde_input_data,   only       : NSOA
  
  use gde_input_data,   only       : A_OR1,A_OR2,A_OR3
  use gde_input_data,   only       : A_OR4,A_OR5,A_OR6
  use gde_input_data,   only       : A_OR7,A_OR8,A_OR9
  use gde_input_data,   only       : A_SAL,A_CHL
  use gde_input_data,   only       : A_EBC,A_DUS
  use gde_input_data,   only       : A_XXX,A_WAT
  
  use gde_input_data,   only       : SVI,MSA,NO3,DMA,NH4
  use gde_input_data,   only       : IOD,OXA,SUC,SSA

  use gde_input_data,   only       : DENAM,DENDU,DENXX,DENEC,DENALK
  use gde_input_data,   only       : lwcpart
  use gde_input_data,   only       : Lthjump,Lvpjump
  use gde_input_data,   only       : alphah2o,alphaCh2o
  use gde_input_data,   only       : nucomin,massmin
  use gde_input_data,   only       : CONVM
  use gde_input_data,   only       : surf_h2o_std,v2,v3
  use gde_toolbox,      only       : surf_succin

  use gde_toolbox,      only       : molec2ug

  use gde_deposition,   only       : settling
  use gde_condensation, only       : condensation_incloud
  use gde_coagulation,  only       : coagulation_coeff
  use gde_coagulation,  only       : coagulation_target
  use gde_coagulation,  only       : collection_kernel
  use gde_coagulation,  only       : coagulation


  private

  public  :: cloud_driver
  public  :: ccnactivation
  public  :: cloud_solver

  private :: insoluble_core
  private :: kinetic_factor
  private :: water_diffusivity
  private :: thermal_conductivity
  private :: curvat_effect
  private :: solute_effect


contains

subroutine cloud_driver(vupdra,temp,press,supersat,dwldtsum,hour_time, &
                        dir,dTdt,dSdt )
    !----------------------------------------------------------------------
    !
    !****  Driver for fog/cloud cycle
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate cloud droplet activation
    !      During each hour with activation conditions (incloud=2), 
    !      one cloud cycle occurs.
    !      The cloud cycle has three stages:
    !      1) Conditioning
    !      2) Growth
    !      3) Dissipation (descend)
    !      The onset of the growth phase is defined by RH>rhactive.
    !      The cloud activation starts at S>Scrit
    !      In this simple cloud scheme, the liquid water content of 
    !      the cloud increases linearly with adiabatic cooling 
    !      (updraft) for 20 min, then stays on a plateau for 20 min, 
    !      and thereafter decreases linearly with warming for 20 min.
    !      The updraft (and downdraft) velocity is constant.
    !      Similar as in Korhonen et al.(2005), Atmos. Chem. Phys.,
    !      5, 2561-2570, www.atmos-chem-phys.org/acp/5/2561/.
    !      Supersaturation at beginning of cloud cycle is zero.
    !
    !
    !      interface
    !      ---------
    !      The updraft velocity is given by the user in dispers.dat
    !      input is DPA and DPAW, the routine computes
    !      dSdt and dDPAW/dt
    !
    !        input:
    !           vupdra    updraft velocity               [m/s]
    !           temp      air temperature                [K]
    !           press     air pressure                   [Pa]
    !           supersat  supersaturation ratio          [-]
    !           dwldtsum  change lwc                     [m^3/m^3]
    !           hour_time time within hour               [s]
    !
    !        output:
    !           dir      sign of tendency               [-]
    !           dTdt     tendency temperature           [K/s]
    !           dSdt     tendency supersaturation       [1/s]
    !
    !
    !      method
    !      ------
    !      The basic equation for the simple cloud scheme
    !      is given by the time rate of change of supersaturation S
    !      for an air parcel rising adiabatically at uniform speed V
    !      from Leaitch et al. (1986):
    !          dS/dt = alpha*V - gamma*dMw/dt
    !      dMw/dt denotes the water condensation rate during the
    !      aerosol activation and subsequent growth process.
    !      dr/dt is the grwoth rate of droplet radius:
    !          dr/dt = (G/r)*( S - A/r + B*ap^3/r^3 )
    !      From this Smax and fraction of actitvated particles is
    !      evaluated.
    !
    !
    !      reference
    !      ---------
    !      W.R. Leaitch, J.W. Strapp, G.A. Isaac,
    !         Cloud droplet nucleation and cloud scavenging
    !         of aerosol sulphate in polluted atmospheres,
    !         Tellus, 38B, 328-344, 1986.
    !      J.H. Seinfeld and S.N. Pandis,
    !         Atmospheric Chemistry and Physics, From Air
    !         Pollution to Climate Change, 2nd Edition,
    !         John Wiley & Sons, Inc., Hoboken New Jersey, 2006.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    use gde_constants,    only : g,cp_air
    use gde_constants,    only : MM_eps,rd

    implicit none

    ! input
    real( dp), intent(in)                            :: vupdra      ! [m/s]
    real( dp), intent(in)                            :: temp        ! [K]
    real( dp), intent(in)                            :: press       ! [Pa]
    real( dp), intent(in)                            :: supersat    ! [-]
    real( dp), intent(in)                            :: dwldtsum    ! [m^3/m^3]
    real( dp), intent(in)                            :: hour_time   ! [s]

    ! output
    real( dp), intent(out)                           :: dir         ! [-]
    real( dp), intent(out)                           :: dTdt        ! [K/s]
    real( dp), intent(out)                           :: dSdt        ! [1/s]

    ! local
    real( dp)                                        :: lathe
    real( dp)                                        :: psat
    real( dp)                                        :: alphav
    real( dp)                                        :: gammaw


      !----------------------------------------------------------
      ! Cloud cycle (T-profile)
      ! The temperature profile is adopted from Case 1 of
      ! Pandis and Seinfeld, 1990a.
      !
      ! 20 min updraft
      if (hour_time.le.60.*20) then
        dir = 1._dp
      endif
      ! 20 min plateau
      ! this period is with dS/dt=0. The beginning of this period
      ! could also be defined as reaching the maximum supersaturation.
      if ((hour_time.gt.60.*20).and.(hour_time.le.60.*40)) then
        dir = 0._dp
      endif
      ! 20 min downdraft
      if (hour_time.gt.60.*40) then
        dir = -1.0_dp
      endif

      !----------------------------------------------------------
      !
      ! CONSTANTS AND UNITS CHECK
      !   g [m/s2], M_H2O [g/mol], M_air in [g/mol]
      !   R_gas   [J/K/mol] (=universal gas constant)
      !   rd      [J/K/kg]  (=gas constant of dry air)    
      !   cp_air  [J/K/kg]
      !   rho_H2O [kg/m3]
      !   lathe   [J/g]
      !   [J]=[Nm]=[m^2 kg s^-2], [Pa]=[Nm^-2]
      !
      !----------------------------------------------------------

      ! latent heat of condensation [J g^-1]
      lathe  = (2501000.-(2370.*(temp-273.15)))/1000.

      ! saturation partial pressure of gaseous H2O (in Pa)
      psat   = 6.982616467E+5_dp &
               - 1.888612677E+4_dp*temp    + 2.132971127E+2_dp*temp**2 &
               - 1.288405237E+0_dp*temp**3 + 4.393186046E-3_dp*temp**4 &
               - 8.023554873E-6_dp*temp**5 + 6.136820929E-9_dp*temp**6

      !----------------------------------------------------------
      ! size-invariant parameter alphav
      !   unit [1/m]
      ! Abdul-Razzak et al., Part 1, EQ(11)
      ! Leaitch et al. (1986) EQ(1)
      ! rd: gas constant of dry air [J/K/kg]
      ! MM_eps: ratio MW(H2O)/MW(air) = 0.622
      ! cp_air: specific heat capacity of dry air
      !         at constant pressure [J/K/kg]
      !----------------------------------------------------------

      alphav = ( (MM_eps*g*lathe*1.e3)                    /  &
                 (cp_air*rd*temp**2) )                    -  &
               ( g/(rd*temp) )

      !----------------------------------------------------------
      ! size-invariant parameter gammaw
      !   unit [-]
      ! Abdul-Razzak et al., Part 1, EQ(12)
      ! Leaitch et al. (1986) EQ(1)
      ! pressure should also change, but for now is constant
      !----------------------------------------------------------

      gammaw = ( press/(MM_eps*psat) )                    + &
               ( (MM_eps*(lathe*1.e3)**2*(1+supersat))    / &
                  (cp_air*rd*temp**2) )

      !----------------------------------------------------------
      ! Time rate of change of supersaturation S
      ! for an air parcel rising adiabatically at uniform speed V
      ! from Leaitch et al. (1986):
      !     dS/dt = alphav*V - gammaw*dwL/dt
      ! Seinfeld and Pandis (2006) EQ (17.79) p. 789
      ! Leaitch et al. (1986) EQ(1)
      ! *
      ! The equation is valid under the assumption that there
      ! is no entrainment into the rising/cooling air parcel.
      ! This means that in absence of water condensation, the 
      ! saturation varies linearly with the updraft velocity.
      ! Supersaturation in a cloud is the result of a balance
      ! between the cooling rate and the liquid water increase.
      ! The liquid water increase is limited by the mass
      ! transport of particles, which depends on particle size
      ! distribution and state of activation.
      ! The cooling rate obtained from this equation is used in
      ! Pandis et al. (1990a) and in Abdul-Razzak et al.
      ! NOTE: 
      ! Case 1 in Pandis et al. (1990a) shows that maximum
      ! supersaturation in a radiative cooling fog
      ! should be S = 0.11% = 0.0011 = 1.1e-3 after 6 minutes
      ! and then decrease
      !----------------------------------------------------------

      dSdt = dir*alphav*vupdra - gammaw*dwldtsum


      !ADIABATIC COOLING: CHANGE TEMPERATURE dT/dt
      ! dT/dz = - lapserate = - g/cp_wet 
      ! = -9.76 degC/km (dry to approx wet)
      ! dT/dt = V * dT/dz

      dTdt = (-1._dp)*(dir*g*vupdra)/cp_air


      !print *,'driver dSdt', dSdt, dir, dir*alphav*vupdra, gammaw*dwldtsum
      !print *,'hour_time',hour_time

      !!! TEST STOP
      !!! Stop after 5 min
      ! if (hour_time.gt. 5.*60) stop
      !!! Koehler test
      ! if (hour_time.gt. 10.*60) stop
      ! Stop after 1 min of activation
      !if (hour_time.gt.21.*60) stop
      ! Stop on plateau
      !if (hour_time.gt.30.*60) stop
      ! Stop in downdraft
      !if (hour_time.gt.50.*60) stop

  end subroutine cloud_driver



subroutine cloud_solver(caq,DTIME,IMAX,press,temp,mbh,DPcrit,    &
                        DPA,DPAW,VPT,ROOPW,M_oc,                 &
                        lwc,caqold,                              &
                        IAG,N,MASS )
    !----------------------------------------------------------------------
    !
    !****  Solver of cloud processes
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      solver of microphysical processes during cloud/fog
    !      activated when incloud=2 and RH>rhactive.
    !      Assuming a non-precipitating cloud.
    !
    !
    !      interface
    !      ---------
    !      The updraft velocity is given by the user in dispers.dat
    !      input is DPA and DPAW, the routine computes changes in
    !      MASS and N during cloud phases.
    !
    !        input:
    !           DTIME    time step                      [s]
    !           caq      aqueous phase concentration    [molec/cm^3]
    !           press    air pressure                   [Pa]
    !           temp     air temperature                [K]
    !           mbh      mixing layer height            [m]
    !           DPcrit   critical particle diameter     [m]
    !           DPA      dry particle diameter          [m]
    !           DPAW     wet particle diameter          [m]
    !           VPT      particle volume                [m^3]
    !           ROOPW    particle density               [kg/m^3]
    !           lwc      liquid water content           [vol(H2O)/vol(air)]
    !           caqold   aq. phase conc old time step   [molec/cm^3]
    !           M_oc     molecular weight SOA species   [g/mol]
    !
    !        in/out:
    !           IAG      coagulation target class, integer matrix
    !           N        number concentration in bin    [1/m^3]
    !           MASS     mass concentration in bin      [ng/m^3]
    !
    !
    !      method
    !      ------
    !      In-Cloud processes:
    !      - Settling of cloud droplets
    !      - Interstitial aerosol scavenging by droplets
    !      - Condensation of dissolved gases
    !      - Collision/Coalescence
    !
    !      Not considered: wet scavenging of particles
    !      by below cloud scavenging
    !
    !      Aqueous aerosol components
    !      in acid-base reactions:
    !      SVI = 1
    !      MSA = 2
    !      NO3 = 3
    !      DMA = 4
    !      NH4 = 5
    !      IOD = 6
    !      OXA = 7
    !      SUC = 8
    !      SAL = 9
    !
    !      reference
    !      ---------
    !      see header
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    use gde_constants,    only : M_H2SO4,M_msa,M_nit,M_dma
    use gde_constants,    only : M_nh3,MNa,MCl,M_hio3
    use gde_constants,    only : pi

    implicit none

    ! input
    integer, intent(in)                                :: IMAX
    real( dp), dimension(nspec), intent(in)            :: caq       ! [molec/cm^3]
    real( dp), intent(in)                              :: DTIME     ! [s]
    real( dp), intent(in)                              :: temp      ! [K]
    real( dp), intent(in)                              :: press     ! [Pa]
    real( dp), intent(in)                              :: mbh       ! [m]
    real( dp), intent(in)                              :: DPcrit    ! [m]
    real( dp), dimension(MMAX,0:(IMAX+1)), intent(in)  :: DPA       ! [m]
    real( dp), dimension(MMAX,0:(IMAX+1)), intent(in)  :: DPAW      ! [m]
    real( dp), dimension(MMAX,IMAX), intent(in)        :: VPT       ! [m^3]
    real( dp), dimension(MMAX,IMAX), intent(in)        :: ROOPW     ! [kg/m^3]
    real( dp), dimension(aqmax,APN), intent(in)        :: caqold    ! [molec/cm^3]
    real( dp), dimension(APN), intent(in)              :: lwc       ! [vol(H2O)/vol(air)]
    real( dp), dimension(NSOA), intent(in)             :: M_oc      ! [g/mol]


    ! output
    real( dp), dimension(MMAX,MMAX,IMAX,IMAX), intent(in out) :: IAG
    real( dp), dimension(MMAX,IMAX,AMAX),intent(inout) :: MASS      ! [ng/m^3]
    real( dp), dimension(MMAX,IMAX),intent(inout)      :: N         ! [1/m^3]

    ! local
    ! maximum growth in aq. phase mode for volatile gases
    real( dp),dimension(aqmax,APN)                     :: caqdiff
    real( dp),dimension(aqmax)                         :: mwaq
    real( dp), dimension(MMAX)                         :: sumdpcs
    real( dp)                                          :: wsurface
    real( dp)                                          :: coags
    real( dp)                                          :: coagscav

    integer                                            :: zkc
    integer                                            :: L
    integer                                            :: I,M,J,O
    integer                                            :: MM,K

    ! allocatable
    real( dp),allocatable,dimension(:,:,:)             :: fluxm
    real( dp),allocatable,dimension(:,:,:)             :: fluxcm
    real( dp),allocatable,dimension(:,:)               :: flux
    real( dp),allocatable,dimension(:,:)               :: fluxc
    real( dp),allocatable,dimension(:,:,:,:)           :: kcoag
    real( dp),allocatable,dimension(:,:,:,:)           :: kcoll
    real( dp),allocatable,dimension(:,:,:)             :: massold
    real( dp),allocatable,dimension(:,:)               :: depo
    real( dp),allocatable,dimension(:,:)               :: vterm

      !----------------------------------------------------------
      ! Allocate aerosol terms

      if (.not. allocated(fluxm))         ALLOCATE(fluxm(MMAX,IMAX,AMAX))
      if (.not. allocated(flux))          ALLOCATE(flux(MMAX,IMAX))
      if (.not. allocated(fluxcm))        ALLOCATE(fluxcm(MMAX,IMAX,aqmax))
      if (.not. allocated(fluxc))         ALLOCATE(fluxc(MMAX,IMAX))
      if (.not. allocated(kcoag))         ALLOCATE(kcoag(MMAX,MMAX,IMAX,IMAX))
      if (.not. allocated(kcoll))         ALLOCATE(kcoll(MMAX,MMAX,IMAX,IMAX))
      if (.not. allocated(massold))       ALLOCATE(massold(MMAX,IMAX,AMAX))
      if (.not. allocated(depo))          ALLOCATE(depo(MMAX,IMAX))
      if (.not. allocated(vterm))         ALLOCATE(vterm(MMAX,IMAX))


      ! Initialisation
      fluxm(:,:,:)    = 0._dp
      fluxcm(:,:,:)   = 0._dp
      flux(:,:)       = 0._dp
      fluxc(:,:)      = 0._dp
      caqdiff(:,:)    = 0._dp
      kcoag(:,:,:,:)  = 0._dp
      kcoll(:,:,:,:)  = 0._dp
      massold(:,:,:)  = 0._dp
      depo(:,:)       = 0._dp
      vterm(:,:)      = 0._dp
      caqdiff(:,:)    = 0._dp
      coags           = 0._dp 
      coagscav        = 0._dp


      ! ******************************************************************
      !     IN-CLOUD AEROSOL PROCESSES
      ! ******************************************************************


      ! CONDENSATIONAL GROWTH
      ! cloud solver is called first after salts are dissolved
      !
      ! THE CAQ CHANGE CAN BE TRANSLATED TO MASS CHANGE -> VPNEW
      ! TO DRIVE THE CONDENSATION FLUXES.
      !
      ! Calculate dCaq/dt [molec/cm^3]
      ! dCaq can be negative if species is depleted
      ! Partitioning of gases only if LWC is large enough
      ! Aqueous aerosol components in acid-base reactions:

      mwaq = (/ M_H2SO4,     &     ! SVI = 1
                M_msa,       &     ! MSA = 2
                M_nit,       &     ! NO3 = 3
                M_dma,       &     ! DMA = 4
                M_nh3+1.,    &     ! NH4 = 5
                M_hio3,      &     ! IOD = 6
                M_oc(1),     &     ! OXA = 7
                M_oc(2),     &     ! SUC = 8
                MNa+MCl /)         ! SAL = 9


      do zkc=1,APN
      !A_SUL (1)
           caqdiff(SVI,zkc) = caq(ind_SO4mm_a(zkc))             + & 
                              caq(ind_HSO4m_a(zkc))             - &
                              caqold(SVI,zkc)

      !A_MSA (2)
           caqdiff(MSA,zkc) = caq(ind_CH3SO3m_a(zkc))           - &
                              caqold(MSA,zkc)

      !A_NIT (3)
           caqdiff(NO3,zkc) = caq(ind_NO3m_a(zkc))              - &
                              caqold(NO3,zkc)

      !A_AMI (4)
           caqdiff(DMA,zkc) = caq(ind_DMAp_a(zkc))              - &
                              caqold(DMA,zkc)

      !A_NH4 (5)
           caqdiff(NH4,zkc) = caq(ind_NH4p_a(zkc))              - &
                              caqold(NH4,zkc)

      !A_IO3 (6)
           caqdiff(IOD,zkc) = caq(ind_IO3m_a(zkc))              - &
                              caqold(IOD,zkc)

! Organics: 
!   OXALAC, MGLYOAC, SUCCAC
      !A_OR1 (7)
           caqdiff(OXA,zkc) = caq(ind_HC2O4m_a(zkc))            + &
                              caq(ind_C2O4mm_a(zkc))            + &
                              caq(ind_CH3COCOOm_a(zkc))         - &
                              caqold(OXA,zkc)

      !A_OR2 (8)
           caqdiff(SUC,zkc) = caq(ind_C2H5C2O4m_a(zkc))         + &
                              caq(ind_C2H4C2O4mm_a(zkc))        + &
                              caq(ind_ADIPAC_a(zkc))            + &
                              caq(ind_GLUTARAC_a(zkc))          - &
                              caqold(SUC,zkc)

      !A_SAL (9)
           caqdiff(SSA,zkc) = caq(ind_Clm_a(zkc))               + &
                              caq(ind_Nap_a(zkc))               - &
                              caqold(SSA,zkc)

      ! Convert dCaq to mass change [ng/m^3]
           do L=1,aqmax

             caqdiff(L,zkc) = caqdiff(L,zkc)                    * &
                              molec2ug( mwaq(L) )*1.e3_dp

           enddo

      ! Growth Limit for NO3 (HNO3g) and SSA (HClg)
      ! ! NH4Cl(s) == NH3(g) + HCl(g)
      ! ! NH4NO3(s) == NH3(g) + HNO3(g)
      ! ! NaCl(s) + NO3m == NaNO3(s) + Clm

           caqdiff(SSA,zkc) = min( caqdiff(SSA,zkc), 0.3_dp )   !0.2
           caqdiff(NO3,zkc) = min( caqdiff(NO3,zkc), 0.3_dp )   !0.2

      ! Set growth of oxalic acid to zero, too volatile
           caqdiff(OXA,zkc) = 0.0

          ! if (zkc==3) then
          !   print*,'S(VI) diff [ng/m^3]',zkc,caqdiff(SVI,zkc)
          !   print*,'ORGC  diff [ng/m^3]',zkc,caqdiff(OXA,zkc)+caqdiff(SUC,zkc)
          !   print*,'SALT  diff [ng/m^3]',zkc,caqdiff(SSA,zkc)
          !   print*,'NITR  diff [ng/m^3]',zkc,caqdiff(NO3,zkc)
          !   print*,'AMMO  diff [ng/m^3]',zkc,caqdiff(NH4,zkc)
          ! endif

      enddo  ! aq. phase modes


      ! IN-CLOUD CONDENSATION
      ! Condensation flux must be called for aqueous aerosol
      ! species (resulting from acid-base reactions)
      if (ICOND.eq.1) then

        call condensation_incloud(IMAX,DTIME,VPT,ROOPW,DPA,    &
                                  N,MASS,  caqdiff, DPcrit,    &
                                  fluxc,fluxcm)
      endif



      ! DROPLET DEPOSITION (SETTLING)
      ! settling of large droplets according to
      ! smog-fog-smog cycle by Pandis et al. (1990b)
      ! vterm(M,I) [m/s]
      ! terminal velocity is also needed for
      ! collision/coalesence.
      call settling(DPAW,press,temp,ROOPW,IMAX, vterm)


      ! COLLISION / COALESENCE
      ! the collection kernel K describes the interaction
      ! of two colliding particles. The collision partners
      ! are assummed to fall with their terminal velocities
      ! Vinf. The smaller collision partner (s) will be 
      ! collected by the larger collision partner (b) as soon
      ! as it is inside the swept volume of the collector.
      if (ICOAG.eq.1) then

      ! coagulation coefficient for interstitial aerosol 
      ! kcoag(M,O,I,J) [m^3/s]
      ! removal of interstitial (non-activated) particles
      ! by larger droplets through Brownian diffusion
      ! according to Hoppel et al. (1990).
         call coagulation_coeff(temp,ROOPW,DPAW,kcoag,IMAX)

      ! collection kernel
      ! Kcoll = Ecoa*Ecol*|Vinf_b - Vinf_s|*pi*(rb+rs)**2
      !   Kcoll: collection kernel [m^3/s]
      ! includes kcoag(M,O,I,J) for interst. scavenging
         call collection_kernel(IMAX,temp,press,DPAW,vterm,       &
                                DPcrit,kcoag*1._dp,   kcoll)

      ! collection equation -> number and mass fluxes
         if (DPcrit.lt.2.e-6_dp) then
      ! droplet-particle collection
           call coagulation(temp,DTIME,ROOPW,DPAW,VPT,N,MASS,     &
                         IMAX,IAG,kcoll, coags,fluxm,flux)
         else
      ! particle coagulation
           call coagulation(temp,DTIME,ROOPW,DPAW,VPT,N,MASS,     &
                         IMAX,IAG,kcoag, coags,fluxm,flux)
         endif

      ! compute updated coagulation target classes
         call coagulation_target(VPT,IMAX, IAG)
      endif



      ! ******************************************************************
      !     CHANGE OF PARTICLE NUMBER AND MASS DUE IN-CLOUD PROCESSES
      ! ******************************************************************


      ! Condensational Growth
      ! 1 distribute dm(aq) over weighted particle surface area
      ! 2 condensation of aq. phase mass [ng/m^3]
      ! 3 condensation flux of mass between bins [ng/m^3]
      ! 4 lower bound of mass per bin [ng/m^3]

      if (ICOND.eq.1) then

       ! summation of droplet surface area
       ! ns(Dp)dDp surface size distribution
       sumdpcs(:)=0._dp
        do M=AI,CS
          do I=1,IMAX
             sumdpcs(M) = sumdpcs(M)                           +  &
                          pi*DPA(M,I)*DPA(M,I)*N(M,I)
          enddo
        enddo

        do M=AI,CS
          do I=1,IMAX

          ! weight of particle surface area
            wsurface = pi*DPA(M,I)*DPA(M,I)*N(M,I)/sumdpcs(M)
          ! uniform distribution for AI mode
            if (M==AI) wsurface = 1._dp/IMAX

          ! mass fluxes of the aqueous phase species (1 to aqmax-1)
          ! aqmax = MMAX-2
            do L=1, aqmax-1
              MASS(M,I,L)    = MASS(M,I,L) + caqdiff(L,M-2) * wsurface
              MASS(M,I,L)    = MASS(M,I,L) + fluxcm(M,I,L)
              MASS(M,I,L)    = max(MASS(M,I,L),massmin)
            enddo

          ! mass flux of seasalt (Chlorine)
            MASS(M,I,A_SAL)    = MASS(M,I,A_SAL)               +  &
                                 caqdiff(aqmax,M-2) * wsurface
            MASS(M,I,A_SAL)    = MASS(M,I,A_SAL) + fluxcm(M,I,aqmax)
            MASS(M,I,A_SAL)    = max(MASS(M,I,A_SAL),massmin)

         ! condensation flux N [1/m^3]
         ! only activated particles shall grow
         ! apply non-negative constraint
            if ( DPA(M,I).gt.DPcrit ) then
              N(M,I)=N(M,I)+fluxc(M,I)
              N(M,I)=max(N(M,I),0.0_dp)    
              ! write(6,*) 'Fcond',M,I, fluxc(M,I), N(M,I)
            endif

          enddo
        enddo

      endif



      ! COLLISION / COALESENCE
      ! Number and Mass loss of 
      ! -  non-activated particles with DPA < DPcrit
      ! -  smaller collision partner
      ! Mass gain (and number shift)
      ! -  by all activated particle sizes
      ! -  larger collision partner
      ! Collection flux in [1/m^3/s]
      if (ICOAG.eq.1) then

        do M=NU,CS
          do I=1,IMAX
             IF (N(M,I) .GT. nucomin) THEN
               do K=1,AMAX   
                 MASS(M,I,K)=MASS(M,I,K)+FLUXM(M,I,K)*DTIME  
                 MASS(M,I,K)=max(MASS(M,I,K),0.0_dp)             
               end do 
             ENDIF
             N(M,I)=N(M,I)+FLUX(M,I)*DTIME
             ! apply non-negative constraint because N(M,I) has already
             ! been changed by condensation flux
             N(M,I)=max(N(M,I),0.0_dp)
             ! write(6,*) 'Fcoll',M,I,N(M,I),FLUX(M,I)*DTIME
          enddo
        enddo

      endif


      ! Droplet deposition
      ! Only particles with diameter >100 nm are falling
      ! Deposition rate depo [1/s]
      ! mbh = cloud height   !  250.0 m 
      if (IDEPO.ge.1) then

        do M=AI,CS
          do I=1,IMAX

            depo(M,I) = vterm(M,I) / max(mbh, 250.0)

            do K=1,AMAX
                MASS(M,I,K)=MASS(M,I,K)                         - &
                MASS(M,I,K)*depo(M,I)*DTIME   
            enddo                                                       
            N(M,I)=N(M,I)-N(M,I)*depo(M,I)*DTIME
            !write(6,*) 'Ndepo',M,I,DPA(M,I),DPAW(M,I),vterm(M,I),depo(M,I)
          enddo
        enddo

      endif



      ! Deallocate aerosol terms
      deallocate(fluxm)
      deallocate(fluxcm)
      deallocate(flux)
      deallocate(fluxc)
      deallocate(kcoag)
      deallocate(kcoll)
      deallocate(massold)
      deallocate(depo)
      deallocate(vterm)


  end subroutine cloud_solver



subroutine ccnactivation(caq,IMAX,temp,press,      &
                         DPA,DPAW,lwc,ROOP,MASS,sforg,surfin,  &
                         M_oc,vupdra,supersat,dir,NTOT,N,      &
                         ddpwdt,dmwdt,dwldtsum,Nfa,Mfa,DPcrit )
    !----------------------------------------------------------------------
    !
    !****  Activation of cloud droplets
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate cloud droplet activation:
    !      All particles whose critical supersaturation is lower
    !      than maximum supersaturation Smax will activate.
    !      The parameterizations of Abdul-Razzad and Ghan are used
    !      the number of activated particles from Smax.
    !      All smaller particles that are not activated remain
    !      in the cloud interstitial. Nothing happens with them
    !      (has to be reconsidered later, according Kerminen, 2001).
    !      increase of water mass dMw/dt
    !      Supersaturation at beginning of cloud cycle is zero.
    !
    !
    !      interface
    !      ---------
    !      The updraft velocity is given by the user in dispers.dat
    !      input is DPA and DPAW, the routine computes
    !      dSdt and dDPAW/dt
    !
    !        input:
    !           caq      aqueous phase concentration    [molec/cm^3]
    !           temp     air temperature                [K]
    !           press    air pressure                   [Pa]
    !           DPA      dry particle diameter          [m]
    !           DPAW     wet particle diameter          [m]
    !           lwc      liquid water content           [vol(H2O)/vol(air)]
    !           ROOP     particle density               [kg/m^3]
    !           MASS     mass concentration in bin      [ng/m^3]
    !           NTOT     total particle number in mode  [1/m^3]
    !           N        number concentration in bin  ! [1/m^3]
    !           sforg    surface tension organic        [kg/s^2]
    !           M_oc     molecular weight SOA species   [g/mol]
    !           vupdra   updraft velocity               [m/s]
    !           supersat supersaturation ratio          [-]
    !           dir      sign of tendency               [-]
    !
    !        output:
    !           ddpwdt   tendency wet diameter          [m/s]
    !           dmwdt    tendency mass concentration    [ng/s]
    !           dwldtsum change lwc                     [m^3/m^3]
    !           Nfa      activated number fraction      [-]
    !           Mfa      activated mass fraction        [-]
    !           DPcrit   critical diameter              [m]
    !
    !      method
    !      ------
    !      The basic equation for the simple cloud scheme
    !      is given by the time rate of change of supersaturation S
    !      for an air parcel rising adiabatically at uniform speed V
    !      from Leaitch et al. (1986):
    !          dS/dt = alpha*V - gamma*dMw/dt
    !      dMw/dt denotes the water condensation rate during the
    !      aerosol activation and subsequent growth process.
    !      dr/dt is the grwoth rate of droplet radius:
    !          dr/dt = (G/r)*( S - A/r + B*ap^3/r^3 )
    !      From this Smax and fraction of actitvated particles is
    !      evaluated.
    !
    !
    !      reference
    !      ---------
    !      see header
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    use gde_constants,    only : g,M_air,cp_air
    use gde_constants,    only : rho_H2O,pi

    implicit none

    ! input
    integer, intent(in)                              :: IMAX
    integer, intent(in)                              :: surfin
    real( dp), dimension(nspec), intent(in)          :: caq         ! [molec/cm^3]
    real( dp), intent(in)                            :: temp        ! [K]
    real( dp), intent(in)                            :: press       ! [Pa]
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA         ! [m]
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPAW        ! [m]
    real( dp), dimension(MMAX,IMAX,AMAX),intent(in)  :: MASS        ! [ng/m^3]
    real( dp), dimension(MMAX,IMAX),intent(in)       :: N           ! [1/m^3]
    real( dp), dimension(APN),intent(in)             :: lwc         ! [vol(H2O)/vol(air)]
    real( dp), dimension(MMAX,IMAX),intent(in)       :: ROOP        ! [kg/m^3]
    real( dp), dimension(NU:CS), intent(in)          :: NTOT        ! [1/m^3]
    real( dp), intent(in)                            :: sforg       ! [kg/s^2]
    real( dp), dimension(NSOA), intent(in)           :: M_oc        ! [g/mol]
    real( dp), intent(in)                            :: vupdra      ! [m/s]
    real( dp), intent(in)                            :: supersat    ! [-]
    real( dp), intent(in)                            :: dir         ! [-]


    ! output
    real( dp),dimension(MMAX,IMAX), intent(out)      :: ddpwdt      ! [m/s]
    real( dp),dimension(MMAX,IMAX), intent(out)      :: dmwdt       ! [ng/s]
    real( dp), intent(out)                           :: dwldtsum    ! [m^3/m^3]
    real( dp), intent(out)                           :: Nfa         ! [-]
    real( dp), intent(out)                           :: Mfa         ! [-]
    real( dp), intent(out)                           :: DPcrit      ! [m]


    ! specific gas constant for dry air [J/(kgK)]
    real( dp), parameter                             :: R=287.058   ! [J/kg/K]
    real( dp), parameter                             :: dlim=5.5e-8 !6.0e-8 ! [m]

    ! local
    real( dp)                                        :: MWaero
    real( dp)                                        :: psat
    real( dp)                                        :: satenv
    real( dp)                                        :: scalef
    real( dp)                                        :: dplimit
    real( dp)                                        :: dwldt
    real( dp)                                        :: supermax
    real( dp)                                        :: A,B
    real( dp)                                        :: dpmaxi
    real( dp)                                        :: rho_air
    real( dp)                                        :: frac_nact
    real( dp)                                        :: frac_mact
    real( dp)                                        :: Nact
    real( dp)                                        :: Nsum
    real( dp)                                        :: Mact
    real( dp)                                        :: Msum
    logical  :: scritfound=.false.
    integer  :: I,M


    ! allocatable
    real( dp),allocatable,dimension(:,:)             :: DPAL
    real( dp),allocatable,dimension(:,:)             :: DPAU
    real( dp),allocatable,dimension(:,:)             :: curvat
    real( dp),allocatable,dimension(:,:)             :: hygrosc
    real( dp),allocatable,dimension(:,:)             :: gassol
    real( dp),allocatable,dimension(:,:)             :: kineticg
    real( dp),allocatable,dimension(:,:)             :: watdiff
    real( dp),allocatable,dimension(:,:)             :: vterm
    real( dp),allocatable,dimension(:,:)             :: DPU
    real( dp),allocatable,dimension(:,:)             :: masssol
    real( dp),allocatable,dimension(:,:)             :: masstot
    real( dp),allocatable,dimension(:,:)             :: densu
    real( dp),allocatable,dimension(:,:)             :: epsm
    real( dp),allocatable,dimension(:,:)             :: massh2o
    real( dp),allocatable,dimension(:,:)             :: xmolh2o
    real( dp),allocatable,dimension(:,:)             :: SCll
    real( dp),allocatable,dimension(:,:)             :: SCul
    real( dp),allocatable,dimension(:,:)             :: SCmm
    real( dp),allocatable,dimension(:,:)             :: sateq


      if (.not. allocated(DPAL))          ALLOCATE(DPAL(MMAX,0:(IMAX+1)))
      if (.not. allocated(DPAU))          ALLOCATE(DPAU(MMAX,0:(IMAX+1)))
      if (.not. allocated(DPU))           ALLOCATE( DPU(MMAX,0:(IMAX+1)))
      if (.not. allocated(curvat))        ALLOCATE(curvat(MMAX,IMAX))
      if (.not. allocated(hygrosc))       ALLOCATE(hygrosc(MMAX,IMAX))
      if (.not. allocated(gassol))        ALLOCATE(gassol(MMAX,IMAX))
      if (.not. allocated(kineticg))      ALLOCATE(kineticg(MMAX,IMAX))
      if (.not. allocated(watdiff))       ALLOCATE(watdiff(MMAX,IMAX))
      if (.not. allocated(vterm))         ALLOCATE(vterm(MMAX,IMAX))
      if (.not. allocated(masssol))       ALLOCATE(masssol(MMAX,IMAX))
      if (.not. allocated(masstot))       ALLOCATE(masstot(MMAX,IMAX))
      if (.not. allocated(densu))         ALLOCATE(densu(MMAX,IMAX))
      if (.not. allocated(epsm))          ALLOCATE(epsm(MMAX,IMAX))
      if (.not. allocated(massh2o))       ALLOCATE(massh2o(MMAX,IMAX))
      if (.not. allocated(xmolh2o))       ALLOCATE(xmolh2o(MMAX,IMAX))
      if (.not. allocated(SCll))          ALLOCATE(SCll(MMAX,IMAX))
      if (.not. allocated(SCul))          ALLOCATE(SCul(MMAX,IMAX))
      if (.not. allocated(SCmm))          ALLOCATE(SCmm(MMAX,IMAX))
      if (.not. allocated(sateq))         ALLOCATE(sateq(MMAX,IMAX))

      !----------------------------------------------------------
      !
      ! CONSTANTS AND UNITS CHECK
      !   g [m/s2], M_H2O [g/mol], M_air in [g/mol]
      !   R_gas   [J/K/mol] (=universal gas constant)
      !   rd      [J/K/kg]  (=gas constant of dry air)    
      !   cp_air  [J/K/kg]
      !   rho_H2O [kg/m3]
      !   lathe   [J/g]
      !   [J]=[Nm]=[m^2 kg s^-2], [Pa]=[Nm^-2]
      !
      !----------------------------------------------------------

      !----------------------------------------------------------
      ! droplet growth rate for water condensation/evaporation
      !   calculate diameter change rate dDp/dt
      !   this will be used to calculate:
      !   (a) liquid water concentration
      !   (b) change rate of liquid water concentration
      !   (c) new wet diameter DPAW (of droplets)
      !   (d) change rate of supersaturation
      ! dr/dt is the grwoth rate of droplet radius:
      !      dr/dt = (G/r)*( S - A/r + B*ap^3/r^3 )
      ! r: droplet radius
      ! ap: radius of dry particle (!)
      ! G: kinetic factor
      ! A: curvature effect
      ! B: solute effect (hygroscopicity)
      !----------------------------------------------------------

      ! Initialisation
      ddpwdt(:,:)   = 0._dp
      dmwdt(:,:)    = 0._dp

      dwldtsum=0.0_dp
      Nact=0._dp
      Nsum=0._dp
      Mact=0._dp
      Msum=0._dp


      ! saturation partial pressure of gaseous H2O (in Pa)
      psat   = 6.982616467E+5_dp &
               - 1.888612677E+4_dp*temp    + 2.132971127E+2_dp*temp**2 &
               - 1.288405237E+0_dp*temp**3 + 4.393186046E-3_dp*temp**4 &
               - 8.023554873E-6_dp*temp**5 + 6.136820929E-9_dp*temp**6

      ! Equilibrium saturation ratio per bin
      sateq(:,:) = 1.0_dp



      ! Calculate diameter of insoluble core DPU &
      ! Mass fraction of soluble material
      call insoluble_core(IMAX,MASS,ROOP,DPA,M_oc,masssol,densu,masstot, &
                          epsm,massh2o,xmolh2o,DPU)


      !----------------------------------------------------------
      ! Curvature effect A for water
      ! curvat is in unit [-]
      ! can be approximated by 0.66e-6/temp [m^-1]
      ! Abdul-Razzak et al., Part 1, EQ(5)
      !
      !       Input is aqueous phase concentration 
      !       from KPP in molec cm^-3(air)
      !
      !----------------------------------------------------------

      call curvat_effect(DPAW,temp,sforg,surfin,IMAX, caq,lwc,            &
                         massh2o,xmolh2o,curvat)


      !----------------------------------------------------------
      ! Solute effect B (hygroscopicity)
      !   unit [-]
      ! Abdul-Razzak et al., Part 1, EQ(6)
      ! SOLUTES:
      ! 1 sodium chloride     NaCl               (A_SAL)
      ! 2 sodium sulfate      (Na)2SO4           (A_SAL, A_SUL)
      ! 3 ammonium sulfate    (NH4)2SO4          (A_NH4, A_SUL)
      ! 4 ammonium bisulfate  NH4HSO4            (A_NH4, A_SUL)
      ! 5 ammonium nitrate    NH4NO3             (A_NH4, A_NO3)
      ! 6 organic acids       SuccAcid, OxalAcid (A_OR1, A_OR2)
      ! Dissolved gases:
      ! SO2, H2O2, HCHO, HNO3
      ! SUCCAC, ADIPAC
      !
      !       Input is aqueous phase concentration 
      !       from KPP in molec cm^-3(air)
      !
      !----------------------------------------------------------

      call solute_effect(DPA,DPU,DPAW,IMAX,lwc,massh2o,epsm,densu,N,      &
                         caq,hygrosc,gassol)


      !----------------------------------------------------------
      ! Critical supersaturation & critical diameter
      !
      ! The maximum supersaturation (Smax) reached inside a cloud
      ! is an important parameter.
      ! All particles whose critical supersaturation is lower
      ! than maximum supersaturation will activate.
      ! Smax can be evaluated by dS/dt = 0.
      ! Abdul-Razzak et al. (1998) present a mathematical solution
      ! EQ(31), to determine Smax for the aerosol size distribution
      ! of a single component.
      ! Here we want to evaluate the state of activation
      ! at the current environmental supersaturation.
      ! We set Smax to the current environmental supersaturation.
      !
      ! The iterative procedure to calculate the critical diameter
      ! (DPcrit) is described in
      ! Abdul-Razzak et al., Part 3, EQ(4) and EQs(12)-(15)
      !
      ! Calculate the critical supersaturation of each size
      ! section at the upper and lower limit of the dry diameter
      ! DPA(M,I) is the dry diameter in the middle of the section.
      ! DPcrit is the critical diameter (dry particle).
      !----------------------------------------------------------

      !----------------------------------------------------------
      ! A is used in [m] and has to multiplied by DPAW
      ! A [m] can be estimated as:
      ! A = (0.66/temp)*1.e-6
      ! B is used in [-] and is adjusted by liquid water volume
      ! Literature estimate of B in [-] 
      ! (Seinfeld and Pandis, 2006, p.770)
      ! B = (1.e-18*3.44e13*2.*2.0e-10/98.)/(DPAW(M,I)**3-DPU(M,I)**3)
      !----------------------------------------------------------

      supermax = supersat

      DPcrit = DPA(MMAX,IMAX+1)
      scritfound=.false.

      do M=1,MMAX
        do I=1,IMAX

          A = curvat(M,I)*DPAW(M,I)
          B = hygrosc(M,I)*(DPAW(M,I)**3-DPU(M,I)**3)/DPA(M,I)**3


          if (I.eq.1) then
            if (M.eq.1) then
              DPAL(M,I) = DPA(M,I)
              DPAU(M,I) = DPA(M,I)+0.5_dp*(DPA(M,I+1)-DPA(M,I))
            else
              DPAL(M,I) = DPA(M,I)-0.5_dp*(DPA(M,I)-DPA(M-1,IMAX))
              DPAU(M,I) = DPA(M,I)+0.5_dp*(DPA(M,I+1)-DPA(M,I))
            endif
          endif
          if ((I.gt.1).and.(I.lt.IMAX)) then
            DPAL(M,I) = DPA(M,I)-0.5_dp*(DPA(M,I)-DPA(M,I-1))
            DPAU(M,I) = DPA(M,I)+0.5_dp*(DPA(M,I+1)-DPA(M,I))
          endif
          if (I.eq.IMAX) then
           if (M.eq.MMAX) then
              DPAL(M,I) = DPA(M,I)-0.5_dp*(DPA(M,I)-DPA(M,I-1))
              DPAU(M,I) = DPA(M,I)
            else
              DPAL(M,I) = DPA(M,I)-0.5_dp*(DPA(M,I)-DPA(M,I-1))
              DPAU(M,I) = DPA(M,I)+0.5_dp*(DPA(M+1,1)-DPA(M,I))
            endif
          endif 

          ! Abdul-Razzak et al., Part 3, EQ(4)
          ! for the critical radius
          SCll(M,I) = (2./sqrt(B)) * (2.*A/(3.*DPAL(M,I)))**1.5
          SCul(M,I) = (2./sqrt(B)) * (2.*A/(3.*DPAU(M,I)))**1.5
          SCmm(M,I) = (2./sqrt(B)) * (2.*A/(3.* DPA(M,I)))**1.5


          if ((SCul(M,I).le.supermax).and.(.not.scritfound)) then
            DPcrit = DPA(M,I)
            !print *,'DPcrit (dry) found',DPcrit
            scritfound=.true.
          endif

          !print *,'Aest, A, Best, B ',M,I,(0.66/temp)*1.e-6, A,     &
          !  (1.e-18*3.44e13*2.*2.0e-10/98.)/(DPAW(M,I)**3-DPU(M,I)**3), B

          !print *,'Scrit',M,I,SCll(M,I),SCul(M,I),SCmm(M,I),supermax

        enddo
      enddo

      !!! print *,'Senv',supersat+1._dp

      !----------------------------------------------------------
      ! Equilibrium saturation ratio due to curvature and 
      ! solute effects
      !
      ! According to Laaksonen et al. (JAS 1998, EQ.23), soluble
      ! gases should be treated as third term in Koehler equation:
      ! Seq = 1 + A/Dp - B/Dp^3 - [ v * n_A/n_w ]
      !----------------------------------------------------------

      do M=1,MMAX
         do I=1,IMAX
           sateq(M,I) = 1._dp + curvat(M,I) - hygrosc(M,I) - gassol(M,I)
       !    print *,'sateq',M,I,DPAW(M,I),curvat(M,I),hygrosc(M,I),sateq(M,I)
      !!! koehler test start
          !if ((M==2).and.(I==17))print *,'0.05um',DPA(M,I),DPcrit,DPAW(M,I),sateq(M,I)
          !if ((M==3).and.(I== 6))print *,'0.10um',DPA(M,I),DPcrit,DPAW(M,I),sateq(M,I)
          !if ((M==3).and.(I==16))print *,'0.50um',DPA(M,I),DPcrit,DPAW(M,I),sateq(M,I)
      !!! koehler test end
         enddo
      enddo
      !!! koehler test
      !if(supersat+1._dp > 1.007) stop

      !----------------------------------------------------------
      ! Water mass change per bin
      ! Kinetic factor G for water vapour condensation
      !   unit  [m^2 s^-1]
      !      Abdul-Razzak et al., Part 1, EQ(16)
      ! Water diffusivity coefficient
      !   unit  [m^2 s^-1]
      !----------------------------------------------------------

      call kinetic_factor(DPAW,temp,press,IMAX,sateq,kineticg)

      call water_diffusivity(DPAW,temp,press,IMAX,watdiff)

      ! Environmental saturation ratio
      !   S_env = supersaturation + 1
      satenv  = supersat + 1.0_dp

      !----------------------------------------------------------
      ! TUNING FACTOR FOR dDp/dt
      !  Tuned to fit decay rate in Arctic sea fog,
      !  SF1 in: Zhao et al. (2022), Atm. Environ.,
      !  doi: 10.1016/j.atmosenv.2022.118943
      !  larger scalef causes smaller S and LWC
      !----------------------------------------------------------
      !scalef  = 0.125_dp  ! f=1/8
      !scalef  = 0.250_dp  ! f=1/4
      scalef  = 0.333_dp  ! f=1/3  SF1
      !scalef  = 1._dp

      ! Lower diameter limit for water condensation/evaporation
      dplimit = min( Dpcrit,dlim )


      !----------------------------------------------------------
      ! Change of wet diameter ddpwdt
      !   unit [m/s]
      ! Seinfeld and Pandis (2006) EQ (17.70) p. 785
      ! change is negative if the equilibrium saturation ratio is
      ! larger than the environmental saturation ratio
      ! dDp/dt is output of this routine and is used in mafor.f90
      ! to change DPAW(M,I)
      ! dMw/dt is output of this routine and is used in mafor.f90
      ! to change MASS(M,I,A_WAT) and LWC
      ! NOTE:
      ! The rate of growth of droplets is inversely proportional
      ! to their diameters so smaller droplets grow faster than
      ! larger ones.
      !----------------------------------------------------------

       do M=1,MMAX
        do I=1,IMAX

            ddpwdt(M,I) = satenv - sateq(M,I)

            ! Dp*dDp/dt
            ddpwdt(M,I) = ddpwdt(M,I) * kineticg(M,I)

            ! Check if ddpwdt is okay:
            ! Dp*dDp/dt != 3.5e-6 * (Senv-Seq) *1.e-4
            ! print *,'ddpwdt',M,I,DPA(M,I),DPAW(M,I),ddpwdt(M,I),  &
            !                    3.5e-6*(satenv -sateq(M,I))*1.e-4

            ! dDp/dt
            ddpwdt(M,I) = ddpwdt(M,I) / DPAW(M,I)
           
            ddpwdt(M,I) = ddpwdt(M,I)*scalef

            ! Change of water mass (ng) per bin per delta_time
            ! dmH2O/dt [Topping et al., SI, 2013]
            dmwdt(M,I)  = (satenv-sateq(M,I)) *ROOP(M,I) * watdiff(M,I)
            dmwdt(M,I)  = dmwdt(M,I)*2._dp*pi *psat * DPAW(M,I)
            dmwdt(M,I)  = dmwdt(M,I)*kineticg(M,I)*CONVM

            ! Change of water mass per particle (ng/m^3/s)
            dmwdt(M,I)  = dmwdt(M,I)*N(M,I)*3._dp

            ! Check if dmw/dt is consistent with dDp/dt
            ! Pandis et al. (1990), Eq.(7)
            ! print *,'dmwdt',M,I,DPAW(M,I),ddpwdt(M,I),            &
            !         dmwdt(M,I)*2./(max(N(M,I),nucomin)*rho_H2O*   &
            !         pi*DPAW(M,I)**2*CONVM)


          ! lower limit for critical diameter
            if ( DPA(M,I).lt.dplimit )   then
               ddpwdt(M,I) = 0._dp
               dmwdt(M,I)  = 0._dp        
            endif

          ! slow evaporation of AI, AS at constant T (05.11.2023)
            if (dir.eq.0._dp) then
               ddpwdt(AI,I) = ddpwdt(AI,I)*0.02_dp
               ddpwdt(AS,I) = ddpwdt(AS,I)*0.02_dp
               dmwdt(AI,I)  = dmwdt(AI,I) *0.02_dp
               dmwdt(AS,I)  = dmwdt(AS,I) *0.02_dp
            endif

          ! no evaporation in updraft/cooling
            if (dir.eq.1._dp) then
               ddpwdt(M,I) = max(ddpwdt(M,I),0._dp)
               dmwdt(M,I)  = max(dmwdt(M,I),0._dp)
            endif

          ! no condensation in downdraft/warming
            if (dir.eq.-1._dp) then
               ddpwdt(M,I) = min(ddpwdt(M,I),0._dp)
               dmwdt(M,I)  = min(dmwdt(M,I),0._dp)
            endif

        !  write(6,'(2I3,5ES12.4)') M,I,dir,DPA(M,I),DPAW(M,I),satenv -sateq(M,I),dmwdt(M,I)

        enddo
      enddo


      !----------------------------------------------------------
      ! Change of liquid water mixing ratio [vol(H2O)/vol(air)]
      ! in unit [m^3/m^3]
      !
      ! Abdul-Razzak et al., Part 1, EQ(13)
      ! dwL/dt = 4*pi*rho_H2O*SUM_i,imax(r**2*dr/dt*N_i(->act))
      !      r = Dp*0.5, dr/dt = 0.5*Dp/dt
      !
      ! Only activate particles count in the summation:
      !     N(M,I)--> DPA(M,I) > Dpcrit
      !----------------------------------------------------------

      do M=AI,CS
        do I=1,IMAX

          !!!if (DPA(M,I).gt.Dpcrit) then
          if (DPA(M,I).gt.dplimit) then

             if ((M.eq.CS).and.(N(M,I).lt.1._dp)) then

               dwldt = 4._dp * pi                             * &
                       (NTOT(CS)/IMAX)*(0.5_dp*DPAW(M,I))**2

             else

               dwldt = 4._dp * pi                             * &
                        N(M,I) * (0.5_dp*DPAW(M,I))**2

             endif

             dwldt = dwldt * 0.5_dp*ddpwdt(M,I)

          else

             dwldt = 0._dp

          endif


      ! summation of liquid water mixing ratio over size bins
      ! for particles above critical diameter
          if (DPA(M,I).ge.3.0e-8) then

            dwldtsum = dwldtsum + dwldt

          endif

        enddo
      enddo

      ! dwldtsum is the change of liquid water that is used
      ! to update the supersaturation dS/dt in cloud_driver 

      ! MULTIPLY WL SUM WITH RATIO RHO_H2O/RHO_AIR
      ! air density [kg/m^3] 1.2928
       rho_air   = press/R/temp
       dwldtsum  = dwldtsum * rho_H2O/rho_air

      ! print *,'wlsum',dwldtsum

      !----------------------------------------------------------
      ! Determine fraction of activated particles in size bin
      ! frac_nact: number fraction activated
      ! frac_mact: mass fraction activated
      ! NOTE: activated mass fraction is for the soluble aerosol
      ! Divide three cases when comparing Scrit to Smax
      !----------------------------------------------------------

      ! Critical diameter DPcrit
      dpmaxi = DPcrit
      Nact   = 0._dp
      Nsum   = 0._dp

      do M=1,MMAX
        do I=1,IMAX

          if (SCul(M,I).gt.supermax) then
            frac_nact = 0._dp
            frac_mact = 0._dp
          endif
          if ((SCll(M,I).ge.supermax).and.(SCul(M,I).le.supermax)) then
            frac_nact = log(supermax/SCll(M,I))/log(SCul(M,I)/SCll(M,I))
            dpmaxi = max(dpmaxi,DPAL(M,I))
            frac_mact = (dpmaxi**3-DPAL(M,I)**3)/(DPAU(M,I)**3-DPAL(M,I)**3)
          endif
          if (SCll(M,I).lt.supermax) then
            frac_nact = 1._dp
            frac_mact = 1._dp
          endif
        ! Fraction of activated particles and mass
        ! Only soluble mass gets activated
          Nact = Nact + N(M,I)*frac_nact
          Nsum = Nsum + N(M,I)
          Mact = Mact + masssol(M,I)*frac_mact
          Msum = Msum + masstot(M,I)

          !print *,'N act',M,I,N(M,I),frac_nact,Nact,Nsum

        enddo
      enddo
      Nfa = Nact/max(Nsum,nucomin)
      Mfa = Mact/max(Msum,massmin)



! Deallocate aerosol terms
       deallocate(DPU)
       deallocate(DPAL)
       deallocate(DPAU)
       deallocate(kineticg)
       deallocate(watdiff)
       deallocate(curvat)
       deallocate(hygrosc)
       deallocate(gassol)
       deallocate(vterm)
       deallocate(masssol)
       deallocate(masstot)
       deallocate(densu)
       deallocate(epsm)
       deallocate(massh2o)
       deallocate(xmolh2o)
       deallocate(SCll)
       deallocate(SCul)
       deallocate(SCmm)
       deallocate(sateq)


  end subroutine ccnactivation





subroutine insoluble_core(IMAX,MASS,ROOP,DPA,M_oc,masssol,densu,masstot, &
                          epsm,massh2o,xmolh2o,DPU)
    !----------------------------------------------------------------------
    !
    !****  Determine insoluble core of droplets
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Calculate diameter of insoluble core DPU &
    !      Mass fraction of soluble material
    !
    !      Insoluble aerosol components are:
    !      A_OR7
    !      A_OR8
    !      A_OR9
    !      A_EBC
    !      A_DUS
    !      A_XXX (may be treated as non-ionic electrolytes later)
    !
    !      interface
    !      ---------
    !
    !        input:
    !           MASS     mass concentration in bin      [ng/m^3]
    !           DPA      dry particle diameter          [m]
    !           ROOP     particle density               [kg/m^3]
    !           M_oc     molecular weight SOA species   [g/mol] 
    !
    !        output:
    !           masssol  soluble mass concentration     [ng/m^3]
    !           masstot  total mass per size bin        [ng/m^3]
    !           densu    density of insoluble mass      [kg/m^3]
    !           epsm     mass fraction soluble material [-]
    !           massh2o  mass of water per size bin     [ng/m^3]
    !           xmolh2o  mole fraction of water         [mol/mol]
    !           DPU      diameter of insolube core      [m^2 s^-1]
    !
    !      method
    !      ------
    !      Seinfeld and Pandis (2006) EQ (17.39) p. 775
    !
    !      reference
    !      ---------
    !      Seinfeld and Pandis (2006) EQ (17.39) p. 775
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    use gde_constants,    only : M_H2O

    implicit none

    ! input
    integer, intent(in)                              :: IMAX
    real( dp), dimension(MMAX,IMAX,AMAX),intent(in)  :: MASS     ! [ng/m^3]
    real( dp), dimension(MMAX,IMAX),intent(in)       :: ROOP     ! [kg/m^3]
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA      ! [m]
    real( dp), dimension(NSOA), intent(in)           :: M_oc     ! [g/mol]

    ! output
    real( dp),dimension(MMAX,IMAX),intent(out)       :: masssol
    real( dp),dimension(MMAX,IMAX),intent(out)       :: densu
    real( dp),dimension(MMAX,IMAX),intent(out)       :: epsm
    real( dp),dimension(MMAX,IMAX),intent(out)       :: massh2o
    real( dp),dimension(MMAX,IMAX),intent(out)       :: xmolh2o 
    real( dp),dimension(MMAX,0:(IMAX+1)),intent(out) :: DPU
    real( dp),dimension(MMAX,IMAX),intent(out)       :: masstot

    !local
    real( dp)                                        :: sumsolm
    real( dp)                                        :: sumunsolm
    real( dp)                                        :: moleh2o
    real( dp)                                        :: moleor1
    real( dp)                                        :: moleor2
    real( dp)                                        :: moleor3
    real( dp)                                        :: moleor4
    real( dp)                                        :: moleor5
    real( dp)                                        :: moleor6

    integer  :: I,M,Q

      do M=1,MMAX             
        do I=1,IMAX

          sumsolm=0._dp
          sumunsolm=0._dp
         ! Soluble species
          do Q=1,A_OR6
            sumsolm = sumsolm + MASS(M,I,Q)
          enddo
         ! Insoluble species
          do Q=A_OR7,A_XXX
            sumunsolm = sumunsolm + MASS(M,I,Q)
          enddo
          sumsolm   = sumsolm + MASS(M,I,A_SAL)
          sumunsolm = sumunsolm - MASS(M,I,A_SAL) - MASS(M,I,A_CHL)

         ! Total mass per size bin [ng/m^3]
          masstot(M,I) = sumsolm+sumunsolm
          masssol(M,I) = sumsolm

         ! Mass of water per size bin [ng/m^3]
          massh2o(M,I) = MASS(M,I,A_WAT)

         ! Mole fraction of water in mixture with organics [mol/mol]
          if (massh2o(M,I).eq.0._dp) then
            xmolh2o(M,I) = 1._dp
          else
            moleh2o     = massh2o(M,I)*1.e-9/M_H2O
            moleor1     = MASS(M,I,A_OR1)*1.e-9/M_oc(1)
            moleor2     = MASS(M,I,A_OR2)*1.e-9/M_oc(2)
            moleor3     = MASS(M,I,A_OR3)*1.e-9/M_oc(3)
            moleor4     = MASS(M,I,A_OR4)*1.e-9/M_oc(4)
            moleor5     = MASS(M,I,A_OR5)*1.e-9/M_oc(5)
            moleor6     = MASS(M,I,A_OR6)*1.e-9/M_oc(6)

            xmolh2o(M,I) = 1._dp-( (2.*moleor1+2.*moleor2   + &
                                    2.*moleor3+2.*moleor4   + &
                                    2.*moleor5+2.*moleor6)/moleh2o )
          endif

         ! Mass fraction of soluble material &
         ! Density of unsoluble material [kg/m^3]
          if (sumunsolm.eq.0._dp) then
            epsm(M,I)  = 1._dp
            densu(M,I) = 1200._dp
          else if ((sumsolm+sumunsolm).eq.0._dp) then
            epsm(M,I)  = 0._dp
            densu(M,I) = 1200._dp  
          else
            epsm(M,I)  = sumsolm/(sumsolm+sumunsolm)
            ! maximum 0.9999
            epsm(M,I)  = min(epsm(M,I),0.9999_dp)

            densu(M,I) =  (MASS(M,I,A_OR7)*DENALK                  + &
                           MASS(M,I,A_OR8)*DENALK                  + &
                           MASS(M,I,A_OR9)*DENALK                  + &
                           MASS(M,I,A_XXX)*DENXX                   + &
                           MASS(M,I,A_EBC)*DENEC                   + &
                           MASS(M,I,A_DUS)*DENDU)/sumunsolm
          endif

         ! Diameter of insoluble core
          DPU(M,I) = DPA(M,I)                                      / &
                     ( (densu(M,I)/ROOP(M,I))                      * & 
                       (epsm(M,I)/(1._dp-epsm(M,I))+1._dp ) )**(1/3)
          !print *,'dpu',M,I,xmolh2o(M,I),epsm(M,I),DPU(M,I)

        enddo
      enddo


  end subroutine insoluble_core


subroutine water_diffusivity(DPAW,temp,press,IMAX,watdiff)
    !----------------------------------------------------------------------
    !
    !****  Modified water diffusivity
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate modified water diffusivity coefficient
    !
    !      water accomodation of 0.96 is used (M.Z. Jacobson)
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature              [K]
    !           press    air pressure                 [Pa]
    !           DPAW     wet particle diameter        [m]
    !        output:
    !           watdiff  water diffusivity coeff.     [m^2 s^-1]
    !
    !      method
    !      ------
    !      Seinfeld and Pandis, p. 801, Eq. (15.66), 1998.
    !
    !      reference
    !      ---------
    !      Seinfeld and Pandis, p. 801, Eq. (15.66), 1998.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    use gde_constants,    only : pi,M_H2O,R_gas,N_A,M_air

    implicit none

    ! input
    INTEGER, intent(in)                              :: IMAX
    REAL( dp), intent(in)                            :: temp     ! [K]
    REAL( dp), intent(in)                            :: press    ! [Pa]
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPAW     ! [m]

    REAL( dp),dimension(MMAX,IMAX), intent(out)      :: watdiff

    REAL( dp), PARAMETER                   :: colldpH2O=3.11e-10 ! [m]
    REAL( dp), PARAMETER                   :: T0=288.15          ! [K]
    REAL( dp)                              :: diffH2O,rho_air
    real( dp)                              :: rp

    integer  :: I,M

      watdiff(:,:) = 0._dp
      
      !density of dry air [g m^-3]
      !  R_gas is in [J/K/mol], M_air in [g/mol]
      ! [J]=[Nm]=[m^2kgs^-2] and [Pa]=[Nm^-2]
      rho_air = press*M_air/(R_gas*T0)

      !molecular water diffusivity (in m^2 s^-1)
      diffH2O = (5./(16.*N_A*colldpH2O*colldpH2O*rho_air*1.e-3_dp))   * &
                SQRT(((R_gas*temp*M_air*1.e-3_dp)/(2*pi))             * &
                ((M_air+M_H2O)/M_H2O))

      
      !modified water diffusivity (in m^2 s^-1)
      !watdiff*1.e4 with alphaC=1 comparable to 
      !Seinfeld and Pandis, p. 784, Fig. 17.13, 2006.

      do M=1,MMAX                
        do I=1,IMAX

          rp=0.5*DPAW(M,I)
          watdiff(M,I) = diffH2O /                            &
                        (   ( rp/(rp+Lvpjump) )                      +  &
                           ( diffH2O/(rp*alphaCh2o) )                *  &
                        SQRT( (2*pi*M_H2O*1.e-3_dp)/(R_gas*temp) )  )

          !print *,'diffH',M,I,DPAW(M,I),diffH2O,watdiff(M,I)*1.e4_dp ! Dv in cm^2/s

        end do
      end do


  end subroutine water_diffusivity


subroutine thermal_conductivity(DPAW,temp,press,IMAX,thermcon)
    !----------------------------------------------------------------------
    !
    !****  Modified thermal conductivity
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate modified thermal conductivity coefficient
    !
    !      thermal water accomodation of 0.96 is used (M.Z. Jacobson)
    !
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature              [K]
    !           press    air pressure                 [Pa]
    !           DPAW     wet particle diameter        [m]
    !        output:
    !           thermcon  thermal conductivity coeff. [J m^-1 s^-1 K^-1]
    !
    !      method
    !      ------
    !      Seinfeld and Pandis, p. 804, Eq. (15.76), 1998.
    !
    !      reference
    !      ---------
    !      Seinfeld and Pandis, p. 804, Eq. (15.76), 1998.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    use gde_constants,    only : cp_air,pi,M_air,R_gas

    implicit none

    ! input
    INTEGER, intent(in)                              :: IMAX
    REAL( dp), intent(in)                            :: temp     ! [K]
    REAL( dp), intent(in)                            :: press    ! [Pa]
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPAW     ! [m]

    REAL( dp),dimension(MMAX,IMAX), intent(out)      :: thermcon


    REAL( dp), PARAMETER                   :: T0=288.15         ! [K]
    REAL( dp)                              :: Kdry,rho_air,cp_wet
    real( dp)                              :: rp

    integer  :: I,M

      thermcon(:,:) = 0._dp
      
      !density of dry air [g m^-3]
      !  R_gas is in [J/K/mol], M_air in [g/mol]
      ! [J]=[Nm]=[m^2kgs^-2] and [Pa]=[Nm^-2]
      rho_air = press*M_air/(R_gas*T0)
      
      !specific heat of wet air [J g^-1 K^-1]
      cp_wet = cp_air*1.e-3*(1.+0.856*0.03)

      !thermal conductivity of dry air [J m^-1 s^-1 K^-1]
      !Seinfeld and Pandis (2006) p. 786
      Kdry = 1.e-3*(4.39+0.071*temp)

      !modified thermal conductivity (in J m^-1 s^-1 K^-1)

      do M=1,MMAX                
        do I=1,IMAX

          rp=0.5*DPAW(M,I)

          thermcon(M,I) = Kdry /                                &
                        (   ( rp/(rp+Lthjump) )                        +  &
                           ( Kdry/(rp*alphah2o*cp_wet*rho_air) )       *  &
                        SQRT( (2*pi*M_air*1.e-3_dp)/(R_gas*temp) )  )

          !print *,'Kdry',M,I,DPAW(M,I),Kdry,thermcon(M,I)


        end do
      end do


  end subroutine thermal_conductivity


subroutine kinetic_factor(DPAW,temp,press,IMAX,sateq,kineticg)
    !----------------------------------------------------------------------
    !
    !****  Kinetic factor for water vapour condensation
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate kinetic factor G [m^2 s^-1]
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature              [K]
    !           press    air pressure                 [Pa]
    !           DPAW     wet particle diameter        [m]
    !           sateq    equil. saturation ratio      [-]
    !        output:
    !           kineticg kinetic factor               [m^2/s]
    !
    !      method
    !      ------
    !      Seinfeld and Pandis, p. 803, Eq. (15.74), 1998.
    !
    !      reference
    !      ---------
    !      Seinfeld and Pandis, p. 803, Eq. (15.74), 1998.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    use gde_constants,    only : Rv_h2o,rho_H2O

    implicit none

    ! input
    INTEGER, intent(in)                              :: IMAX
    REAL( dp), intent(in)                            :: temp     ! [K]
    REAL( dp), intent(in)                            :: press    ! [Pa]
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPAW     ! [m]
    REAL( dp),dimension(MMAX,IMAX), intent(in)       :: sateq    ! [-]

    REAL( dp),dimension(MMAX,IMAX), intent(out)      :: kineticg

    REAL( dp),dimension(MMAX,IMAX)         :: watdiff,thermcon
    REAL( dp)                              :: psat,lathe
    REAL( dp)                              :: termg1,termg2
    REAL( dp)                              :: termg3
    
    integer  :: I,M

      kineticg(:,:) = 0._dp
            
      !saturation partial pressure of gaseous H2O (in Pa)
      ! [Pa]=[N/m^2]=[J/m^3], [J]=[Nm]
      psat = 6.112*exp((17.67*(temp-273.15))&
           /((temp-273.15)+243.5))*100.
      
      !latent heat of condensation (in J g^-1)
      lathe = (2501000.-(2370.*(temp-273.15)))/1000.

      !modified water diffusivity  (in m^2 s^-1)
      CALL water_diffusivity(DPAW,temp,press,IMAX,watdiff)  

      !modified thermal conductivity (in J m^-1 s^-1 K^-1)
      CALL thermal_conductivity(DPAW,temp,press,IMAX,thermcon)



      !second term [dimless]
      termg2 = (lathe/(Rv_h2o*temp))-1._dp

      ! Note:   Rv_H2O [J/g/K]   is R_gas/M_H2O  
      !         rho_H2O [kg/m^3] is density of H2O
      !         
      do M=1,MMAX                
        do I=1,IMAX
          
          termg1 = (rho_H2O*Rv_h2o*1.e3_dp*temp)  /  &
                   (4.0_dp*psat*watdiff(M,I))                ! [m^-2*s]

          termg3 = (rho_H2O*lathe*1.e3_dp*sateq(M,I))   /  &
                   (4.0_dp*thermcon(M,I)*temp)               ! [m^-2*s]

          !kinetic factor G
          !  unit [m^2/s]   
          kineticg(M,I) = 1./(termg1+(termg2*termg3))

          ! Seinfeld and Pandis (2006) p. 786: at 283 K
          ! term1         = 1.21x10^9  m^-2*s   @ 10 um diameter
          ! term2 * term3 = 1.60*10^9  m^-2*s   @ 10 um diameter
          ! kineticg      = 3.6e-10    m^2/s
          ! computed kineticg = 2.94e-10 m^2/s at 277.5 K

         !print *,'Kinetic G',M,I,DPAW(M,I),termg1,termg2*termg3,kineticg(M,I)

        end do
      end do

      
  end subroutine kinetic_factor


subroutine curvat_effect(DPAW,temp,sforg,surfin,IMAX,caq,lwc,  &
                         massh2o,xmolh2o, curvat)
    !----------------------------------------------------------------------
    !
    !****  Curvature effect
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate curvature coefficient [dimless]
    !      "Kelvin effect"
    !      Seinfeld and Pandis (2006) EQ (17.28) p. 770
    !      curvature = A = approx 0.66e-6/T/Dp
    !
    !      Surface tension of mixture of dissolved organics
    !      and water
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature              [K]
    !           sforg    surface tension organic      [kg/s^2]
    !           DPAW     wet particle diameter        [m]
    !           lwc      liquid water content         [vol(H2O)/vol(air)]
    !           massh2o  mass concentration of water  [ng/m^3]
    !           xmolh2o  molar ratio water            [mol/mol]
    !           caq      aq. phase concentration      [molec/cm^3]
    !
    !        output:
    !           curvat  curvature coeff. A            [-]
    !
    !      method
    !      ------
    !      H. Abdul-Razzak, S.J. Ghan, C. Rivera-Carpio,
    !        A parameterization of aerosol activation,
    !        1. Single aerosol type, J. Geophys. Res.,
    !        103, D6, 6123-6131, 1998.
    !
    !      reference
    !      ---------
    !      H. Abdul-Razzak, S.J. Ghan, C. Rivera-Carpio,
    !        A parameterization of aerosol activation,
    !        1. Single aerosol type, J. Geophys. Res.,
    !        103, D6, 6123-6131, 1998.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    use gde_constants,    only : Rv_h2o,rho_H2O,M_H2O,N_A

    implicit none

    integer, intent(in)                              :: IMAX
    integer, intent(in)                              :: surfin
    real( dp), intent(in)                            :: temp     ! [K]
    real( dp), intent(in)                            :: sforg    ! [kg/s^2]
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPAW     ! [m]
    real( dp), dimension(APN),intent(in)             :: lwc      ! [vol(H2O)/vol(air)]
    real( dp), dimension(MMAX,IMAX),intent(in)       :: massh2o  ! [ng/m^3]
    real( dp), dimension(MMAX,IMAX),intent(in)       :: xmolh2o  ! [mol/mol]
    real( dp), dimension(nspec), intent(in)          :: caq      ! [molec/cm^3]

    real( dp),dimension(MMAX,IMAX), intent(out)      :: curvat   ! [-]

    real( dp), dimension(MMAX)             :: sumdp
    real( dp), dimension(MMAX)             :: maqDOC

    real( dp)                              :: surf_H2O
    real( dp)                              :: sigma_H2O
    real( dp)                              :: sigma_org
    real( dp)                              :: surftp

    real( dp)                              :: mpH2O
    real( dp)                              :: xmH2O
    real( dp)                              :: naqH2O
    real( dp)                              :: naqDOC


    integer  :: I,M,zkc


      curvat(:,:) = 1.e-12_dp

      ! summation droplet diameter
       do M=AI,CS              
         sumdp(M)=0._dp
         do I=1,IMAX
            sumdp(M) = sumdp(M) + DPAW(M,I)
         enddo
       enddo

      ! Dissolved Organic Matter: soluble organic vapors
      ! aqueous phase modes
      ! convert conc_aq to mola_aq
      !  molec cm^-3(air) to mol(DOC)/gH2O
      !  lwc(:) is lwc of modes AI-CS
      ! total molality in modes AI-CS
      ! INDEX zkc has to match aerosol mode M
      ! zkc   M
      ! 1     3
      ! 2     4
      ! 3     5
      ! caq: replace index for solute by zkc+1 if more than one aqueous mode
      maqDOC(:)=0._dp
      
      do zkc=1,APN
        maqDOC(zkc+2) = (caq(ind_OXALAC_a(zkc))+caq(ind_SUCCAC_a(zkc))+    &
                    caq(ind_GLUTARAC_a(zkc))+caq(ind_ADIPAC_a(zkc))+       &
                    caq(ind_HC2O4m_a(zkc))+caq(ind_C2O4mm_a(zkc))+         &
                    caq(ind_C2H5C2O4m_a(zkc))+caq(ind_C2H4C2O4mm_a(zkc)))* &
                    1.e3*v2/(N_A*lwc(1)*rho_H2O)
      enddo

      ! surface tension of water (in kg/s^2)
      surf_H2O = surf_h2o_std - 0.155*(temp - 273.15)
      sigma_H2O=surf_H2O*1.e-3_dp

      ! surface tension of organics [kg/s^2]
      if (surfin.eq.1) then
        sigma_org=surf_succin(temp)
      else
        sigma_org=sforg
      endif

      ! curvature effect in NU mode              
      ! surface tension of the mixture [kg s^-2]
      ! Seinfeld and Pandis (2006) EQ(17.A.34), p. 818
      
      do M=NU,NA
        do I=1,IMAX

          if (xmolh2o(M,I).eq.1._dp) then
            surftp = sigma_H2O
          else
            surftp = xmolh2o(M,I)*( sigma_H2O-sigma_org ) + &
                     sigma_org
          endif
          curvat(M,I) = 4._dp*surftp                      / &
                        (Rv_h2o*1.e3*temp*rho_H2O         * &
                        DPAW(M,I))

                                                 ! curvature in [m]    control in [m]
          !print*,'kelvin A',M,I,surftp,curvat(M,I),curvat(M,I)*DPAW(M,I),0.66e-6/temp

        enddo
      enddo

      ! curvature effect in AI-CS mode
      do M=AI,CS
        do I=1,IMAX

          mpH2O  = massh2o(M,I)*1.e-9_dp     ![gH2O in bin]

        ! Weighting with droplet mass distribution
        ! Largest droplet has largest mass
          naqDOC  = maqDOC(M) * DPAW(M,I)/sumdp(M)

        ! Solute's mol, n_s  [mol]
          naqDOC  = naqDOC * mpH2O

        ! Water's mol, n_w   [mol]
          naqH2O  = mpH2O / M_H2O  

        ! mole ratio X of water in bin [mol/mol]
          if (naqH2O.eq.0._dp) then
            xmH2O     = 1._dp
            sigma_org = 0._dp
          else
            xmH2O = naqH2O/(naqH2O+naqDOC)
          endif

        ! surface tension of the mixture [kg s^-2]
        ! Seinfeld and Pandis (2006) EQ(17.A.34), p. 818
          surftp = xmH2O*( sigma_H2O-sigma_org ) + sigma_org

        ! curvature effect of mixture (water+dissolved organics)
        ! dim. curvat: [kg/s2] / ([m2*kg/(1e-3kg*K*s2)]*K*[kg/m3]*m)
        ! curvat is dimensionless

          curvat(M,I) = 4._dp*surftp                   / &
                        (Rv_h2o*1.e3*temp*rho_H2O      * &
                         DPAW(M,I))

                                                 ! curvature in [m]    control in [m]
          !print*,'kelvin A',M,I,surftp,curvat(M,I),curvat(M,I)*DPAW(M,I),0.66e-6/temp

        enddo
      enddo


  end subroutine curvat_effect


subroutine solute_effect(DPA,DPU,DPAW,IMAX,lwc,massh2o,epsm,densu,  &
                         N,caq,    hygrosc,gassol)
    !----------------------------------------------------------------------
    !
    !****  Solute effect (Hygroscopicity) B
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate solute coefficient B [dimless]
    !      "Raoult effect"
    !      Seinfeld and Pandis (2006) EQ (17.29) p. 770
    !      hygrosc = B/(Dp^3-Dpu^3) = 
    !               approx 1e-18*(3.44e-13*2*mass[g]/M_ABS)/ &
    !                      (Dp^3-Dpu^3)
    !
    !      Solutes:
    !      NaCl       c(ind_Clm_a)
    !      (Na)2SO4   c(ind_Nap_a)   [Na+]>[Cl-]
    !      (NH4)2SO4  c(ind_SO4mm_a) [SO42-]-[Na+]
    !      NH4HSO4    c(ind_HSO4m_a)
    !      NH4NO3     c(ind_NO3m_a)
    !      NH4IO3     c(ind_IO3m_a)
    !      OXALATE    c(ind_C2O4mm_a)
    !      SUCCINATE  c(C2H4C2O4mm_a)
    !      Dissolved gases:
    !      SO2,H2O2,HCHO,HNO3,SUCCAC,ADIPAC
    !
    !      interface
    !      ---------
    !      aqueous phase concentrations from KPP
    !
    !        input:
    !           DPA      dry particle diameter        [m]
    !           DPU      insoluble particle diameter  [m]
    !           DPAW     wet particle diameter        [m]
    !           lwc      liquid water content         [vol(H2O)/vol(air)]
    !           massh2o  mass concentration of water  [ng/m^3]
    !           epsm     mass fraction of soluble     [-]
    !           densu    density of unsoluble         [kg/m^3]
    !           N        number conc. per bin         [1/m^3]
    !           caq      aq. phase concentration      [molec/cm^3]
    !
    !        output:
    !           hygrosc  solute coeff.                [-]
    !           gassol   soluble gas effect           [-]
    !
    !      method
    !      ------
    !      H. Abdul-Razzak, S.J. Ghan, C. Rivera-Carpio,
    !        A parameterization of aerosol activation,
    !        1. Single aerosol type, J. Geophys. Res.,
    !        103, D6, 6123-6131, 1998.
    !
    !      reference
    !      ---------
    !      H. Abdul-Razzak, S.J. Ghan, C. Rivera-Carpio,
    !        A parameterization of aerosol activation,
    !        1. Single aerosol type, J. Geophys. Res.,
    !        103, D6, 6123-6131, 1998.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    use gde_constants,    only : pi,M_H2O,rho_H2O,M_ABS
    use gde_constants,    only : N_A


    implicit none

    integer, intent(in)                              :: IMAX
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA      ! [m]
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPU      ! [m]
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPAW     ! [m]
    real( dp), dimension(APN),intent(in)             :: lwc      ! [vol(H2O)/vol(air)]
    real( dp), dimension(MMAX,IMAX),intent(in)       :: massh2o  ! [ng/m^3]
    real( dp), dimension(MMAX,IMAX),intent(in)       :: epsm     ! [-]
    real( dp), dimension(MMAX,IMAX),intent(in)       :: densu    ! [kg/m^3]
    real( dp), dimension(MMAX,IMAX),intent(in)       :: N        ! [1/m^3]
    real( dp), dimension(nspec), intent(in)          :: caq      ! [molec/cm^3]

    real( dp),dimension(MMAX,IMAX), intent(out)      :: hygrosc  ! [-]
    real( dp),dimension(MMAX,IMAX), intent(out)      :: gassol   ! [-]

    real( dp), parameter                   :: hygmin = 1.e-12_dp

    real( dp), dimension(MMAX)             :: sumdpcs
    real( dp), dimension(MMAX)             :: maqCL
    real( dp), dimension(MMAX)             :: maqNA
    real( dp), dimension(MMAX)             :: maqSO4
    real( dp), dimension(MMAX)             :: maqHSO
    real( dp), dimension(MMAX)             :: maqNH4
    real( dp), dimension(MMAX)             :: maqNO3
    real( dp), dimension(MMAX)             :: maqIO3
    real( dp), dimension(MMAX)             :: maqOXA
    real( dp), dimension(MMAX)             :: maqSUC
    real( dp), dimension(MMAX)             :: maqSO2
    real( dp), dimension(MMAX)             :: maqH2O2
    real( dp), dimension(MMAX)             :: maqHCHO
    real( dp), dimension(MMAX)             :: maqHNO3
    real( dp), dimension(MMAX)             :: maqSUCA
    real( dp), dimension(MMAX)             :: maqADPA

    real( dp)                              :: msolute
    real( dp)                              :: gsolute
    real( dp)                              :: nsolute
    real( dp)                              :: nsolbin
    real( dp)                              :: nsolgas
    real( dp)                              :: mpH2O
    real( dp)                              :: wsurface


    integer  :: I,M
    integer  :: zkc


      hygrosc(:,:) = hygmin
      gassol(:,:)  = 0._dp

      ! summation of droplet surface area
      ! ns(Dp)dDp surface size distribution
      sumdpcs(:)=0._dp
      do M=AI,CS
        do I=1,IMAX
          sumdpcs(M) = sumdpcs(M)                     +  &
                       pi*DPA(M,I)*DPA(M,I)*N(M,I)
        enddo
      enddo
 
      ! AQUEOUS PHASE MODES
      ! convert conc_aq to mola_aq
      ! Cl: molec cm^-3(air) to mol(NaCl)/gH2O
      !   lwc(3) is lwc of coarse mode
      ! total molality in coarse mode
      ! INDEX zkc has to match aerosol mode M
      ! zkc   M
      ! 1     3
      ! 2     4
      ! 3     5
      ! caq: replace index for solute by zkc+1 if more than one aqueous mode

      do zkc=1,APN 

      ! solute NaCl
        maqCL(zkc+2)  = (caq(ind_Clm_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O)
        maqNa(zkc+2)  = (caq(ind_Nap_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O)

      !!! solute Na2SO4
      !!    maqNA(zkc+2)  = ((caq(ind_Nap_a(zkc))-caq(ind_Clm_a(zkc)))*1.e3*v3)/ &
      !!                  (N_A*lwc(zkc)*rho_H2O)
      !!    maqNA(zkc+2)  = max(maqNA(zkc+1),0._dp)

      ! solute (NH4)2SO4
      !!  maqSO4(zkc+2) = ((caq(ind_SO4mm_a(zkc))-caq(ind_Nap_a(zkc))        + &
      !!                    caq(ind_Clm_a(zkc)))*1.e3*v3)/(N_A*lwc(zkc)*rho_H2O)
      !!  maqSO4(zkc+2) = max(maqSO4(zkc+1),0._dp)
        maqSO4(zkc+2) = (caq(ind_SO4mm_a(zkc))*1.e3*v3)/(N_A*lwc(zkc)*rho_H2O)
        maqNH4(zkc+2) = (caq(ind_NH4p_a(zkc) )*1.e3*(v2+v3)*0.5)/(N_A*lwc(zkc)*rho_H2O)

      ! solute NH4HSO4
        maqHSO(zkc+2) = (caq(ind_HSO4m_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O)

      ! solute NH4NO3
        maqNO3(zkc+2) = (caq(ind_NO3m_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O)

      ! solute NH4IO3
        maqIO3(zkc+2) = (caq(ind_IO3m_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O)

      ! solute OXALAC (from VOC + OH)
        maqOXA(zkc+2) = (caq(ind_HC2O4m_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O) +&
                        (caq(ind_C2O4mm_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O)

      ! solute SUCCAC (from BSOV)
        maqSUC(zkc+2) = (caq(ind_C2H5C2O4m_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O) +&
                        (caq(ind_C2H4C2O4mm_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O)

      ! solute ADIPAC (from BLOV): no acid-base reaction

      ! solute gases (SO2, H2O2, HCHO, HNO3, SUCCAC, ADIPAC)
        maqSO2 (zkc+2) = (caq(ind_so2_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O)
        maqH2O2(zkc+2) = (caq(ind_h2o2_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O)
        maqHCHO(zkc+2) = (caq(ind_hcho_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O)
        maqHNO3(zkc+2) = (caq(ind_hno3_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O)
        maqSUCA(zkc+2) = (caq(ind_SUCCAC_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O)
        maqADPA(zkc+2) = (caq(ind_ADIPAC_a(zkc))*1.e3*v2)/(N_A*lwc(zkc)*rho_H2O)

      enddo

      ! NUCLEATION mode
      ! consider only NH4HSO4 as mole solute n_s
      ! Seinfeld and Pandis (2006), EQ(17.40), page 775
      ! that take into account unsoluble material
      ! To calculate hygrosc use the equation 
      ! Seinfeld and Pandis (2006), EQ(17.29), page 770

      do M=NU,NA
        do I=1,IMAX

       ! mol solute per bin
          nsolute = ( epsm(M,I)*v2*pi*DPA(M,I)**3 )            /  &
                    ( 6._dp*M_ABS*1.e-3                        *  &
                      ( (epsm(M,I)/DENAM)                      +  &
                        ( (1._dp-epsm(M,I))/densu(M,I) ) ) )

       ! nsolute per bin
          nsolute = nsolute*0.1_dp

       ! DPU is diameter of unsoluble core
                   
          hygrosc(M,I) = ( 6._dp * M_H2O*1.e-3_dp * nsolute )  /  &
                         (pi * rho_H2O                         *  &
                         (DPAW(M,I)**3 - DPU(M,I)**3) )

        
          hygrosc(M,I) = max(hygrosc(M,I),hygmin)
          gassol(M,I)  = 0._dp

          !print *,'Raoult B',M,I,nsolute,(DPAW(M,I)**3-DPU(M,I)**3),hygrosc(M,I), &
          !        (1.e-18*3.44e13*2.*2.0e-10/M_ABS)/(DPAW(M,I)**3-DPU(M,I)**3)  !! estimate

        enddo
      enddo


      ! solute effect in AI-CS mode
      ! DPU is diameter of unsoluble core
      do M=AI,CS
        msolute = 0._dp
        gsolute = 0._dp
        ! dissolved salts in [mol/gH2O]
        msolute = maqCL(M)+maqNA(M) +maqSO4(M)+maqHSO(M)+maqNH4(M)+maqNO3(M)
        msolute = msolute+maqIO3(M)
        msolute = msolute+maqOXA(M)+maqSUC(M)
        ! dissolved gases in [mol/gH2O]
        gsolute = maqSO2(M)+maqH2O2(M)+maqHCHO(M)+maqHNO3(M)       
        gsolute = gsolute+maqSUCA(M)+maqADPA(M)

       ! print *,'msolute',M,lwc(M-1),msolute

        do I=1,IMAX

        ! Note that massh2o is weighted according to
        ! droplet number distribution in gde_addwater
          mpH2O  = massh2o(M,I)*1.e-9_dp     ![gH2O in bin]

        ! Weighting with droplet mass distribution
        ! Uniform distribution
        !!!  nsolbin = msolute/IMAX

        ! Particle with largest surface area have
        ! largest mass
        ! weight of particle surface area
          wsurface = pi*DPA(M,I)*DPA(M,I)*N(M,I)/sumdpcs(M)
          nsolbin = msolute * wsurface

        ! Solute's mol, n_s  [mol]
          nsolbin = nsolbin*mpH2O

        ! Solute mol per particle
        !!!nsolbin = nsolbin  / max((ntota(M)/IMAX),nucomin)
          nsolbin = nsolbin  / max( N(M,I), nucomin )

          !print *,'Raoult B',M,I,DPAW(M,I),nsolbin,ntota(M)/IMAX


        ! Dissolved gas term: gassol
        ! According to Laaksonen et al. (JAS 1998), soluble gases
        ! should be treated as third term in Koehler equation with
        ! Seq = 1 + A/Dp - B/Dp^3 - [ v * n_A/n_w ] (EQ.23)
        ! Solute (gas), n_g [mol/gH2O]
          nsolgas = gsolute * wsurface
        ! Solute (gas) mol per mole water
          gassol(M,I) = nsolgas * M_H2O


        ! Avoid diameter increase of AS over CS
          if (M.eq.AI) nsolbin = nsolbin * 0.6_dp
          if (M.eq.AS) nsolbin = nsolbin * 0.6_dp 

      !!! Koehler test start
          !nsolbin = 3.422E-18                    ! NaCl      bin1
          !nsolbin = 3.422E-18*8.                 ! NaCl      bin2
          !if (M.eq.AS) nsolbin = 3.422E-18*60.   ! NaCl      bin3
          !nsolbin = 2.268E-18                    ! (NH4)2SO4 bin1
          !nsolbin = 2.268E-18*8.                 ! (NH4)2SO4 bin2
          !if (M.eq.AS) nsolbin = 2.268E-18*60.   ! (NH4)2SO4 bin3
      !!! Koehler test end

        ! Solute term: B*(DP^3-DPu^3)

          hygrosc(M,I) = ( 6._dp * M_H2O*1.e-3_dp              *  & 
                           nsolbin )                           /  &
                         (pi * rho_H2O                         *  &
                         (DPAW(M,I)**3 - DPU(M,I)**3) )

          !!print *,'DP',M,I,DPA(M,I),DPU(M,I),DPAW(M,I)
          !print *,'Raoult B',M,I,                         &
          !         DPAW(M,I)**3, DPU(M,I)**3,             &
          !         nsolbin      ,                         &
          !         hygrosc(M,I)* (DPAW(M,I)**3-DPU(M,I)**3),           &
          !         (1.e-18*3.44e13*2.*2.0e-16/M_ABS)    !! simple estimate

        ! Minimum value if no solute
          hygrosc(M,I) = max(hygrosc(M,I),hygmin)

        enddo
      enddo


  end subroutine solute_effect


  end module gde_cloud_activation
