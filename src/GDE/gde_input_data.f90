! <gde_input_data.f90 - A component of the Multicomponent
!                     Aerosol Dynamics Model MAFOR>
!*****************************************************************************! 
!* 
!*    Copyright (C) 2011-2024  Matthias Steffen Karl
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
!*****************************************************************************!
!*    All routines written by Matthias Karl
!* 
!*****************************************************************************!
module gde_input_data

   use messy_mecca_kpp_Parameters, only  : nspec, dp

  implicit none

  !
  ! definition of the chemistry: #reactions, order of species, etc.
  ! parameters needed for chemistry and rate constants
  ! in messy_mecca_kpp_Parameters

!----------------------------------------------------------------
!   number of gas phase input tracers
     integer, parameter    :: gspec=81

!----------------------------------------------------------------
!   GDE Aerosol Parameters
!   Aerosol Modes, number of modes: 5
!
! IMPORTANT NOTE:
! If any of the aerosol indices below are changed
! do the corresponding changes in MESA/gde_thermo_data.f90
!----------------------------------------------------------------

! Name of aerosol mode
     integer, parameter       :: NU =1
     integer, parameter       :: NA =2  ! new nanoparticle mode
     integer, parameter       :: AI =3
     integer, parameter       :: AS =4
     integer, parameter       :: CS =5  
! Number of aerosol modes
     integer, parameter       :: MMAX=5
! Number of aerosol components 
     integer, parameter       :: AMAX=21
! Number of condensable aerosol components
     integer, parameter       :: QMAX=16
! Number of SOA components
     integer, parameter       :: NSOA=9
! Number of input aerosol components
     integer, parameter       :: iamax=10
! Number of aqueous aerosol species
     integer, parameter       :: aqmax=9

! INDEX declaration aerosol components
!  condensable components       
     integer, parameter       :: A_SUL = 1
     integer, parameter       :: A_MSA = 2
     integer, parameter       :: A_NIT = 3
     integer, parameter       :: A_AMI = 4
     integer, parameter       :: A_NH4 = 5
     integer, parameter       :: A_IO3 = 6
     integer, parameter       :: A_OR1 = 7
     integer, parameter       :: A_OR2 = 8
     integer, parameter       :: A_OR3 = 9
     integer, parameter       :: A_OR4 = 10
     integer, parameter       :: A_OR5 = 11
     integer, parameter       :: A_OR6 = 12
     integer, parameter       :: A_OR7 = 13
     integer, parameter       :: A_OR8 = 14
     integer, parameter       :: A_OR9 = 15
     integer, parameter       :: A_CHL = 16
!  non-volatile components
     integer, parameter       :: A_SAL = 17 ! 11
     integer, parameter       :: A_EBC = 18
     integer, parameter       :: A_DUS = 19
     integer, parameter       :: A_XXX = 20
     integer, parameter       :: A_WAT = 21

!   Aerosol Components from input
!        number of components: 10
     integer, parameter       :: SU = 1
     integer, parameter       :: OC = 2
     integer, parameter       :: AM = 3
     integer, parameter       :: NI = 4
     integer, parameter       :: MS = 5
     integer, parameter       :: SA = 6
     integer, parameter       :: XX = 7   !PBA
     integer, parameter       :: EC = 8
     integer, parameter       :: DU = 9
     integer, parameter       :: WA = 10

!   Aqueous aerosol species
     integer, parameter       :: SVI = 1
     integer, parameter       :: MSA = 2
     integer, parameter       :: NO3 = 3
     integer, parameter       :: DMA = 4
     integer, parameter       :: NH4 = 5
     integer, parameter       :: IOD = 6
     integer, parameter       :: OXA = 7
     integer, parameter       :: SUC = 8
     integer, parameter       :: SSA = 9

!----------------------------------------------------------------
! mafor static arrays
! Static arrays
     real( dp), dimension(nspec),save       :: emis
     real( dp), dimension(nspec),save       :: emis0
     real( dp), dimension(nspec),save       :: vdry
     real( dp),dimension(iamax), save       :: DEN
     real( dp),dimension(MMAX,iamax), save  :: BGMCTOT
     real( dp),dimension(MMAX,iamax), save  :: EMMCTOT

     real( dp),save           :: DPMIN
     real( dp),save           :: DPMAX
     real( dp),save           :: radiusm(NU:CS) 
     real( dp),save           :: NTOT(NU:CS),VTOT(NU:CS),ROOPWAV(NU:CS)
     real( dp),save           :: MSULFTOT(NU:CS),MORGCTOT(NU:CS),MSALTTOT(NU:CS)
     real( dp),save           :: MAMMOTOT(NU:CS),MNITRTOT(NU:CS),MMSAPTOT(NU:CS)
     real( dp),save           :: MIODATOT(NU:CS)
     real( dp),save           :: MXXXXTOT(NU:CS),MECBCTOT(NU:CS),MDUSTTOT(NU:CS)
     real( dp),save           :: MH2OTOT(NU:CS)
     real( dp),save           :: MORG1TOT(NU:CS),MORG2TOT(NU:CS),MORG3TOT(NU:CS)
     real( dp),save           :: MORG4TOT(NU:CS),MORG5TOT(NU:CS),MORG6TOT(NU:CS)
     real( dp),save           :: MORG7TOT(NU:CS),MORG8TOT(NU:CS),MORG9TOT(NU:CS)
     real( dp),save           :: ppmtot(NU:CS)
     real( dp),save           :: mtotin(NU:CS)
     real( dp),save           :: MTOT(NU:CS)
     real( dp),save           :: MTOTW(NU:CS)
     real( dp),save           :: BGMD(NU:CS),BGSIG(NU:CS),BGNTOT(NU:CS)
     real( dp),save           :: EGMD(NU:CS),ESIG(NU:CS),ENTOT(NU:CS)
     real( dp),save           :: DNVAPC
     real( dp),save           :: VVAP,NNUC,JNUC,N3,NEU1,NEU2,NEU3,NEU4,NEU5,NEU6,NEU7
     real( dp),save           :: CTNH4,CTNIT,CTSO4,CTCHL
     real( dp),save           :: KPEQ
     real( dp),save           :: GROC,GRSU,GRMS,GRTOT,DCORG,DCSU,CSSULT,CSORGT,COAGS
     real( dp),save           :: TENDE
     real( dp),save           :: DEB_NTOT,DEB_CTOTS1,DEB_CTOTS2,DEB_MTOT
     real( dp),save           :: DEB_CTOTA1,DEB_CTOTA2
     real( dp),save           :: DEB_CTOTO11,DEB_CTOTO12
     real( dp),save           :: DEB_CTOTO21,DEB_CTOTO22
     real( dp),save           :: model_time,dummy_time
     real( dp),save           :: hour_time,hour_timec,hour_timet,hour_timeo
     real( dp),save           :: hour_after_fog
     real( dp),save           :: hour_timemin,dilut_time,dilut2_time
     real( dp),save           :: cpu_time_start, cpu_time_end
! Parameterized particle number concentration
     real( dp),save           :: PNC1,PNC2,PNC3,vd1,vd2,vd3
     real( dp),save           :: kcoa1,kcoa2,kcoa3,zmblm
     real( dp),save           :: PNC4,PNC5,PNC6,vd4,vd5,vd6
     real( dp),save           :: kcoa4,kcoa5,kcoa6


!----------------------------------------------------------------
! mafor input variables
!
! real
! (general)
     real( dp),save           :: edms,eso2,eh2o2,cnh3,rain,cdms
     real( dp),save           :: camidoh,co3  
     real( dp),save           :: lat_deg,daynr,daytime,starttime
     real( dp),save           :: u10,lon_deg
     real( dp),save           :: RH,zmbl,owf
     real( dp),save           :: supersat
     real( dp),save           :: VEN,BSOV_FT
! (chamber)     
     real( dp),save           :: alphanit
     real( dp),save           :: cno,cno2,cco3,jno2m,cipn
! (plume)    
     real( dp),save           :: dila,dilcoef
     real( dp),save           :: dilrate
     real( dp),save           :: dilstore
     real( dp),save           :: emisratp 
     real( dp),save           :: zmbl_old,temp_old
     real( dp),save           :: wplm_old,wplm
     real( dp),save           :: TAIR
     real( dp),save           :: areapl
! (nucleation)
     real( dp),save           :: natot,fnuc     
! integer/char/logic
     integer,save             :: iday,imonth,runtime
     integer,save             :: mday,month,jday,gn,zkc
     integer,save             :: incloud

!----------------------------------------------------------------
! Particle phase parameters

! Numerical particle number conc. threshold [#/m^3]
     real( dp), parameter     :: nucomin=1.e-7_dp 
! Particle density [g cm^-3]
!  (H2SO4 particles)
     real( dp), parameter     :: pdens=1.8_dp
! Conversion kg into ng
     real( dp), parameter     :: CONVM=1.E12_dp
! Numerical component mass conc. threshold  [ng/m^3]     
     real( dp),PARAMETER      :: massmin=1.E-12_dp

! Aerosol Dynamics
!  mean free path of gas [m]
     real( dp), parameter     :: LAM=6.6E-8_dp
!  viscosity [kg/(ms)]
     real( dp), parameter     :: VIS=1.81E-5_dp
!  nucleation parameter
     real( dp), parameter     :: CNUC=1.E-20_dp

! Condensing Vapor
!   mass of a vapor molecule [kg]
     real( dp), parameter     :: MVAP=3.E-26_dp
!
!   vapor density [kg/m3]   H2SO4
! densities from Karl et al. 2007, JGR
     real( dp), parameter     :: DENV=1770._dp
! density organic aerosol [kg/m3]
!  from "organic.dat"
!     REAL( dp), parameter :: DENOC=1570._dp
!   density ammonium aerosol [kg/m3]
! prel. data from Barbara D'Anna    
     real( dp), parameter     :: DENAM=1300._dp
!   density nitrate aerosol [kg/m3]
! prel. data from Barbara D'Anna  
     real( dp), parameter     :: DENNI=1300._dp
!   density MSAp aerosol [kg/m3]
     real( dp), parameter     :: DENMS=1770._dp
!   density iodate aerosol [kg/m3]
     real( dp), parameter     :: DENIO=4629._dp
!   density sea salt aerosol [kg/m3]
     real( dp), parameter     :: DENSA=2240._dp
!   density biological organic aerosol [kg/m3]
     real( dp), parameter     :: DENXX=1150._dp
!   density soot aerosol [kg/m3] Lemmetty et al. AST 2008
     real( dp), parameter     :: DENEC=1200._dp
!   density mineral dust [kg/m3]
!      check literature!!!
     real( dp), parameter     :: DENDU=1400._dp     
!   density C16-C30 n-alkanes [kg/m3]
     real( dp), parameter     :: DENALK=900._dp
!   density generic organics [kg/m^3] 
     real( dp), parameter     :: DENOC=1570._dp
!
! diffusion coefficient of vapor [m2/s]
     real( dp), parameter     :: DIFV=2.5E-5_dp

! equilibrium concentration of vapor [1/m3]
! is now calculated as p0sat of H2SO4 
!     REAL( dp), parameter :: NSVAP=1.E+10_dp
! saturation vapour pressure of VOC [1/m3]
! Pirjola and Kulmala 2001: <= 1.0e12  
!     REAL( dp), parameter :: NSVAPORG=1.0e12_dp     
! is now calculated as p0sat of succinic acid

!----------------------------------------------------------------
! Cloud droplet activation
   ! threshold LWD for partitioning [m^3/m^3]
     real( dp), parameter     :: lwcpart=1.e-10_dp
   ! threshold RH for activation (99%) [-]
     real( dp), parameter     :: rhactiv=0.99_dp
     real( dp), parameter     :: rhstart=0.98_dp
   ! thermal jump length [m]
     real( dp), parameter     :: Lthjump=2.16E-7_dp
   ! vapor jump length [m]
     real( dp), parameter     :: Lvpjump=1.096E-7_dp
   ! thermal accomodation coefficient H2O
     real( dp), parameter     :: alphah2o=0.96_dp
   ! condensation coefficient H2O
     real( dp), parameter     :: alphaCh2o=1.00_dp
   ! surface tension pure H2O [g s^-2] = [dyn cm^-1]
     real( dp), parameter     :: surf_h2o_std = 76.1
   ! van't Hoff coefficients
     real( dp), parameter     :: v2=2._dp
     real( dp), parameter     :: v3=3._dp

!----------------------------------------------------------------
! Aqeuous phase parameters
  ! parameters for exchange gas phase - aqueous phase

    include 'amine.inc'

! Water Vapor
!   collision diameter of H2O [cm]
     real( dp), parameter     :: diamh2o=3.11E-8_dp
!   density of air [g/cm3]
     real( dp), parameter     :: rho_air=0.00123_dp



end module gde_input_data
