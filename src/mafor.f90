! <mafor.f90 - A component of the Multicomponent
!              Aerosol FORmation model MAFOR>
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
!  PROGRAM: MAFOR
!
!          Multicomponent Aerosol FORmation model
!
!  VERSION: 2.0.0
!
!  PURPOSE:    Lagrangian type sectional aerosol box model
!              which includes gas phase and aqueous phase 
!              chemistry in addition to aerosol dynamics. 
!
!
!  REFERENCES:
!    1.  Karl, M., Gross, A., Pirjola, L., Leck, C., A new flexible
!        multicomponent model for the study of aerosol dynamics
!        in the marine boundary layer, Tellus B, 63(5),1001-1025,
!        doi:10.1111/j.1600-0889.2011.00562.x, 2011.
!    2 . Karl, M., Leck, C., Gross, A., Pirjola, L., A study of
!        new particle formation in the marine boundary layer over
!        the central Arctic Ocean using a flexible multicomponent
!        aerosol dynamic model, Tellus B, 64, 17158,
!        doi:10.3402/tellusb.v64i0.17158, 2012a.
!    3.  Karl, M., Dye, C., Schmidbauer, N., Wisthaler, A., Mikoviny, T., 
!        D'Anna, B., M?ller, M., Clemente, E., Munoz, A., Porras, R., 
!        Rodenas, M., Vazquez, M. and T. Brauers, Study of OH-initiated 
!        degradation of 2-aminoethanol, Atmos. Chem. Phys., 12, 
!        1881-1901, 2012b.
!    4.  Karl, M., Kukkonen, J., Keuken, M.P., Lutzenkirchen, S.,
!        Pirjola, L., Hussein, T., Modelling and measurements of urban
!        aerosol processes on the neighborhood scale in Rotterdam,
!        Oslo and Helsinki, Atmos. Chem. Phys., 16,
!        4817-4835, doi:10.5194/acp-16-4817-2016, 2016.
!
!  EXTERNAL CODES AND LICENSES:
!        The MAFOR source code includes:
!        * Gas-phase chemistry and aqueous phase chemistry using the
!          CAABA/MECCA v4.0 / Module Efficiently Calculating the
!          Chemistry of the Atmosphere;
!        * Photolysis routines to calculate J values using the MESSY-JVAL
!            code of Max-Planck Institute for Chemistry, Mainz, Germany;
!        * A package with the Rosenbrock ROS3 solver, KPP-2.2.3
!            symbolic chemistry Kinetics PreProcessor 
!            (http://www.cs.vt.edu/~asandu/Software/KPP)
!            for dynamic creation of Fortran solver code;
!        * A package on thermodynamics equilibrium solver for 
!            ammonium nitrate, MESA.
!
!        The basic gas phase chemistry and aqueous phase chemistry is based 
!        on the community atmospheric chemistry box model CAABA/MECCA-4.0 
!        (Sander et al., 2019):
!        Sander, R., A. Baumgaertner, D. Cabrera-Perez, F. Frank,
!        S. Gromov, J.-U. Grooss, H. Harder, V. Huijnen, P. Joeckel,
!        V.A. Karydis, K. E. Niemeyer, A. Pozzer, H. Riede, M.G. Schultz,
!        D. Taraborelli, and S. Tauer, The community atmospheric
!        chemistry box model CAABA/MECCA-4.0,
!        Geosci. Model Dev., 12, 1365-1385, doi:10.5194/gmd-12-1365-2019, 2019.
!
!        MESSY-JVAL photolysis routines to calculate J values:
!        Original code by Jochen Landgraf (MPICH, until 1998), see:
!        Landgraf, J. and P.J. Crutzen (1998), An efficient method for online 
!        calculations of photolysis and heating rates, J. Atmos. Sci., 55, 
!        863-878.
!        JVAL and the JVal PreProcessor (JVPP) are community models
!        published under the GNU General Public License.
!        The JVAL-14 code was modularized according to the MESSy standard:
!        Sander, R., P. Joeckel, O. Kirner, A.T. Kunert, J. Landgraf,
!        and A. Pozzer, The photolysis module JVAL-14, compatible with the 
!        MESSy standard, and the JVal Preprocessor (JVPP),
!        Geosci. Model Dev., 7, 2653-2662, doi:10.5194/gmd-7-2653-2014, 2014.
!
!        KPP is distributed under the GPL License, (C) 1995-1997, V. Damian
!        and A. Sandu, CGRER, Univ. Ioawa, (C) 1997-2007, A. Sandu, Michigan
!        Tech, Virginia Tech with contributions from: R. Sander, Max-Planck
!        Institute for Chemistry, Mainz, Germany. Any changes to KPP must
!        be reported to the developers of KPP.
!        
!        MESA is part of the MOSAIC code distributed with the community
!        developed model WRF-Chem. Primary Author: Rahul A. Zaveri,
!        Co-Investigators: Richard C. Easter, William I. Gustafson Jr.. 
!        MOSAIC code was prepared by the Battelle Memorial Institute under
!        Contract No. DE-AC05-76RL0 1830 with the Department of Energy (DOE).
!        The following Terms of Usage apply: The MOSAIC code is intended for
!        research and educational purposes. Users are requested to contact
!        the Primary Author regarding the use of MOSAIC code for any 
!        commercial application. Users preparing publications resulting from 
!        the usage of MOSAIC are requested to cite:
!        Zaveri, R.A., R.C. Easter, J.D. Fast, and L.K. Peters (2008), Model
!        for Simulating Aerosol Interactions and Chemistry (MOSAIC), 
!        J. Geophys. Res., 113, D13204, doi:10.1029/2007JD008782.
!
!
!*****************************************************************************! 
!  DESCRIPTION:
!        This FORTRAN90-code solves the general dynamic equation
!        (GDE) for a size distribution of dry multi-component
!        aerosol.
!        The size distribution representation is fixed sectional
!        and the time integration Euler forward.
!        Accumulation and coarse mode of a marine aerosol
!
!        The aerosol is organized in four modes (MMAX=4):
!        NU: nucleation mode
!        AI: Aitken mode
!        AS: accumulation mode
!        CS: coarse mode
!        Each mode is split into a user-defined number of bins IMAX
!        Total number of bins is NB=IMAX*MMAX
!
!        The dynamical processes included are
!            - nucleation             INUC=1
!            - condensation           ICOND=1
!                   H2SO4 and MSA              (ICONS=1,2) 
!                   organic vapor              (ICONO=1)
!                   amine and nitrate          (ICONA=1)  
!                   ammonium nitrate           (ICONX=1)      
!            - coagulation            ICOAG>=1
!            - dry deposition         IDEPO>=1
!            - wet deposition         IWETD=1
! 
!
!        The size distribution representation is fixed sectional
!        and the time integration Euler forward.
!        Accumulation and coarse mode of a marine aerosol
!
!        Note: Introduction of the APC scheme of M.Z. Jacobson in
!              MAFOR V1.1 made the integration routine unconditionally
!              stable in terms of condensation/evaporation.
!        Code: original provided by
!               19.8.2006 Kari E. J. Lehtinen
!              and freely distributed
!        Start date V1.0: 31.10.2011
!
!        Computation of gas/aquoeus phase concentration in molec/cm3(air)
!        Treatment of the aqueous phase per mode
!        Chemistry Integration:
!        KPP-2.2.3 with Rosenbrock solver rosenbrock_posdef
!
!        The modal structure is
!        1) 1-10nm, 10-100nm, 100-1000nm, 1000-10000nm or
!        2) 1-10nm, 10-30nm, 30-100nm, 100-1000nm 
!   
!        If the model is run "in-cloud", the largest aerosol mode (CS) is
!        treated as water drops determined by the prescribed liquid water
!        content (lwc).
!
!        KPP aqueous phase: APN=1 and species *_a01 for aqueous phase 
!        chemistry in CS mode.
!
!        CCN ACTIVATION (preliminary)
!          Cloud droplet activation happens if relative humidity is above 99%
!          and the incloud flag is set to 1 by the user and if the computed
!          supersaturation is greater than equlibrium supersaturation.
!          During cloud droplet activation, the aerosol dynamics processes
!          are stopped.
!          During each hour with activation conditions one cloud cycle occurs.
!          In this simple cloud scheme, the liquid water content of the cloud
!          increases linearly (with height) until t_incloud=0.5 h and thereafter 
!          decreases linearly. The updraft (and downdraft) velocity is constant.
!          To activate, cloud flag and IAQP option must be 1.
!
!
!    ATTENTION: (!) Change of Henry's law constants has to be done in
!                   ChemProp in folder CAABA/mecca/tracer/chemprop
!               (!) xnom7dmso is set to zero
!
!        Sulphuric acid production SO2+OH:
!        SO2 + OH     = HSO3      : k_3rd(temp,cair,3.E-31,-3.3,1.5E-12,0.,0.6)
!        HSO3 + O2    = HO2 + SO3 : 1.3E-12*exp(-330/temp)
!        SO3 + H2O    = H2SO4     : 1.2E-15
!
!
!        Gas-phase chemistry included:
!            - v1.0:  inorganic photochemistry  (St-Tr-G-N)
!            - v1.0:  CH4 chemistry
!            - v1.0:  DMS chemistry             (Tr-G-S)
!            - v1.0:  MEA chemistry (=AMIN)     (Tr-G-A)
!            - v1.1:  C2,C3,C4-NMHC             (Tr-G-C)
!            - v1.1:  Isoprene chemistry        (Tr-G-C)
!            - v1.1:  MMA,DMA,TMA chemistry     (Tr-G-A)
!            - v1.2:  aromatic NMHC             (Tr-G-C)
!            - v1.3:  MIM2 isoprene chemistry   (Tr-G-C)
!            - v1.3:  RO2 lump species
!            - v1.4:  DEA,TEA chemistry         (Tr-G-A)
!            - v1.7:  AMP chemistry             (Tr-G-A)
!            - v1.8:  a-pinene chemistry        (Tr-G-CB)
!            - v1.9.6: Mainz Organic Mechanism [Sander et al., 2019]
!                      incl. iodine reactions   (Tr-(G||Aa)-I)
!            - v1.9.7: incl. chlorine reactions (Tr-(G||Aa)-I-Cl)
!
!         SOA is represented by 9 SOA precursors in the gas phase:
!              BSOV -  SVOC, secondary oxidized biogenic, "COV"
!              BLOV -  LVOC, secondary oxidized biogenic, in Nucl 8-10+12
!              BELV - ELVOC, secondary oxidized biogenic
!              ASOV -  SVOC, secondary oxidized aromatic
!              ALOV -  LVOC, secondary oxidized aromatic
!              AELV - ELVOC, secondary oxidized aromatic
!              PIOV -  IVOC, primary emitted n-alkane (default C21)
!              PSOV -  SVOC, primary emitted n-alkane (default C26)
!              PELV - ELVOC, primary emitted n-alkane (default C34)
!         and 9 SOA components in the particle phase:
!             A_OR1: secondary oxidized biogenic, succinic acid
!             A_OR2: secondary oxidized biogenic, carboxylic acid
!             A_OR3: secondary oxidized biogenic, extr.-low-volatile
!             A_OR4: secondary oxidized aromatic
!             A_OR5: secondary oxidized aromatic
!             A_OR6: secondary oxidized aromatic, extr.-low-volatile
!             A_OR7: primary emitted n-alkane
!             A_OR8: primary emitted n-alkane
!             A_OR9: primary emitted n-alkane, extr.-low-volatile
!         Additional organic carbon:
!             A_XXX: primary biological aerosol  
!
!*****************************************************************************! 
!
!  LIST OF FILE UNITS:
!
!     Input
!            8 ingeod.dat       meteo
!            9 sensitiv.dat     processes
!           20 inchem.dat       gas phase
!           21 inaero.dat       aerosol
!           22 organic.dat      organic vapours
!           23 incham.dat       chamber parameters
!           24 monitor.dat      chamber data 
!           25 emitpar.dat      particle emission
!           26 inbgair.dat      background aerosolplumedilr
!           27 inaqchem.dat     aqueous phase
!           28 dispers.dat      dispersion parameters
!
!     Output
!           11 soadis           soa distribution
!           12 debug.res        debug info output
!           13 total_n.res      total aerosol
!           14 size_dis.res     size distribution
!           15 size_dism.res    mass conc size distribution
!           16 concout.res      gas phase conc
!           17 aerconc.res      aerosol phase conc
!           18 wetdp.res        wet diameter
!           19 plume.res        plume dispersion output
!
!*****************************************************************************! 
!
!  HISTORY:
!         Author: Matthias Steffen Karl (MSK)
!
!
!  CHANGES COMPARED TO VERSION 2.0.0
!
!         dd.mm.yyyy  name   routine and one-line description of modification
!         -------------------------------------------------------------------
!
!
!        
!****************************************************************************

      program mafor

! interface to kpp
      use messy_mecca_kpp,            only : icntrl
      use messy_mecca_kpp_util,       only : initialize_indexarrays
! aqueous phase species _a
      use messy_mecca_kpp_global
! gas-phase species
      use messy_mecca_kpp_parameters, only : ind_N2,ind_O2,ind_CO2,ind_H2O
      use messy_mecca_kpp_parameters, only : ind_DMS,ind_H2O2,ind_NH3
      use messy_mecca_kpp_parameters, only : ind_O3,ind_NO,ind_NO2,ind_SO2
      use messy_mecca_kpp_parameters, only : ind_H2SO4,ind_CH3SO3H
      use messy_mecca_kpp_parameters, only : ind_IPN,ind_N2O5
      use messy_mecca_kpp_parameters, only : ind_ClNO3,ind_HNO3,ind_HCl
      use messy_mecca_kpp_parameters, only : ind_MEA,ind_MMA,ind_TMA
      use messy_mecca_kpp_parameters, only : ind_CH2NCH3,ind_AMP,ind_DMA
 ! condensable organics
      use messy_mecca_kpp_parameters, only : ind_BSOV,ind_BLOV,ind_BELV
      use messy_mecca_kpp_parameters, only : ind_ASOV,ind_ALOV,ind_AELV
      use messy_mecca_kpp_parameters, only : ind_PIOV,ind_PSOV,ind_PELV


      use gde_constants,      only    : R_gas,N_A,k_B,RHOH2O,M_H2O,M_msa,pi
      use gde_constants,      only    : M_nit,OneDay,MB,MAN,MAH,M_H2SO4,AVOng
      use gde_constants,      only    : M_nh3,MH,MC,MO,MCl,STRLEN_KPPSPECIES
      use gde_constants,      only    : rho_H2O,MNa,R_gas
      use gde_input_data
      use gde_toolbox,        only    : conch2o,acidps
      use gde_toolbox,        only    : satpress_alkane,molec2ug   
      use gde_transfer,       only    : transfer_coeff,aero_init_gasaq

      use gde_init_gas,       only    : readchem
      use gde_init_gas,       only    : readchamber
      use gde_init_gas,       only    : readmonit
      use gde_init_gas,       only    : readorganic
      use gde_init_gas,       only    : fco, fcoj
      use gde_init_gas,       only    : IAM,IOZ,F_HONO,KP_NIT,K_DIL,CAMI
      use gde_init_gas,       only    : L_MEA, L_NO2, L_HNO3, L_O3
      use gde_init_gas,       only    : V_CHAM,S_CHAM,CWIN
      use gde_init_gas,       only    : surf_org,surfin
      use gde_init_gas,       only    : gamma_oc1_m,gamma_oc2_m,gamma_oc3_m
      use gde_init_gas,       only    : gamma_oc4_m,gamma_oc5_m,gamma_oc6_m
      use gde_init_gas,       only    : gamma_oc7_m,gamma_oc8_m,gamma_oc9_m
      use gde_init_aqchem,    only    : readaqphase
      use gde_init_aero,      only    : readinaero
      use gde_init_aero,      only    : reademitpar
      use gde_init_aero,      only    : readbackgr
      use gde_init_aero,      only    : DPMMIN, GMD, SIG, NUM

      use gde_chem_gas,       only    : chemistry_solver
      use gde_chem_gas,       only    : emis_drydep,check_range
      use gde_chem_gas,       only    : chamber_dil,chamber_loss

      use gde_photo,          only    : photo_dealloc
      use gde_addwater,       only    : wetdiameter
      use gde_addwater,       only    : water_content
      use gde_aerosol_props,  only    : getDensity
      use gde_aerosol_props,  only    : getTotalMass
      use gde_aerosol_props,  only    : initSizeDistribution
      use gde_aerosol_props,  only    : initNumberMass
      use gde_aerosol_props,  only    : initBGNumberMass
      use gde_aerosol_props,  only    : initEMNumberMass

      use gde_plume,          only    : plumeDisp
      use gde_plume,          only    : readdispers
      use gde_plume,          only    : initplume
      use gde_plume,          only    : plumearea
      use gde_plume,          only    : hmix_st,dst_st,hsta_st,ta_st
      use gde_plume,          only    : BGH2O
      use gde_plume,          only    : vupdra,sst,sal

      use gde_seaflux,        only    : seasaltemis
      use gde_coagulation,    only    : coagulation_target
      use gde_nucleation,     only    : nucleationratio
      use gde_sensitiv,       only    : readsens,ICOAG,ICOND,INUC,INUCMEC
      use gde_sensitiv,       only    : IDEPO,IWETD,ICHAM,IKELV
      use gde_sensitiv,       only    : ICONS,ICONA,ICONO,IPEMI,ICONX,ICHEM
      use gde_sensitiv,       only    : IDMS,IOZO,IORG,IDIL,INANO
      use gde_sensitiv,       only    : ICONW,IAQC,IAQP,IDEB      
      use gde_aerosol_solver, only    : aerosol_solver
      use gde_aerosol_solver, only    : interface_mosaic

          
     implicit none


!   KPP Solver
     INTEGER, PARAMETER :: NBL = 1 ! N_block_length
     REAL(dp), DIMENSION(NBL,nspec) :: cbl     
     INTEGER            :: status          ! error status
     REAL( dp)          :: DTIME           ! time step [s]


!   Define Variables
     integer            :: IMAX
     integer            :: IAS1,IAS3,IAS4,IAS5,IAS6
     integer            :: IAS61,IAS62,IAS63,IAS7
     integer            :: IAS81,IAS82,IAS9
     integer            :: IAS10,IAS11,IAS12,IAS13,IAS14,IAS15,IAS16
     integer            :: IAS17,IAS18,IAS19,IAS20,IAS21,IAS22,IAS23
     integer            :: stat1
! Allocatable arrays
     REAL( dp),allocatable,dimension(:)       :: NVAP,NVAPO
     REAL( dp),allocatable,dimension(:,:)     :: ROOP,ROOPW
     REAL( dp),allocatable,dimension(:,:)     :: DLOGDP,DLINDP,N,LINDIS,LOGDIS,MLOGDIS
     REAL( dp),allocatable,dimension(:,:)     :: BGN,BGLOGDIS,DIFFDILN
     REAL( dp),allocatable,dimension(:,:)     :: EN
     REAL( dp),allocatable,dimension(:,:)     :: VPT,MPT,MPTW
     REAL( dp),allocatable,dimension(:,:)     :: DPA,DPAW
     REAL( dp),allocatable,dimension(:,:)     :: saltemis
     REAL( dp),allocatable,dimension(:,:,:)   :: MASS
     REAL( dp),allocatable,dimension(:,:,:)   :: CCOND
     REAL( dp),allocatable,dimension(:,:,:)   :: BGMASS
     REAL( dp),allocatable,dimension(:,:,:)   :: DIFFDILM
     REAL( dp),allocatable,dimension(:,:,:)   :: EMASS
     REAL( dp),allocatable,dimension(:,:,:,:) :: IAG
!MESA start
     real( dp),allocatable,dimension(:,:)     :: HYST
     real( dp),allocatable,dimension(:,:)     :: henry_eff_hno3
     real( dp),allocatable,dimension(:,:)     :: henry_eff_hcl
     real( dp),allocatable,dimension(:,:)     :: henry_eff_nh3
     real( dp),allocatable,dimension(:,:)     :: svmc_min_hno3
     real( dp),allocatable,dimension(:,:)     :: svmc_min_hcl
     real( dp),allocatable,dimension(:,:)     :: svmc_min_nh3
     real( dp),allocatable,dimension(:,:)     :: kh_nh3
     real( dp),allocatable,dimension(:,:)     :: cioncharge
     integer,allocatable,dimension(:,:)       :: flag_dissolution
     real( dp), save                          :: Keq_nh4no3_0
     real( dp), save                          :: Keq_nh4cl_0
!MESA end
!CCN start
     real( dp),allocatable,dimension(:,:)     :: ddpwdt
     real( dp),allocatable,dimension(:,:)     :: dmwdt
     real( dp),save                           :: dtedt
     real( dp),save                           :: dsupdt
     real( dp),save                           :: Nfa
     real( dp),save                           :: Mfa
     real( dp),save                           :: DPcrit
!CCN end


! Compounds
! c:  gas/aquoeus phase concentration in molec/cm3(air)

! SOA components
     real( dp), dimension(NSOA)        :: M_oc
     real( dp), dimension(NSOA)        :: nmo
     real( dp), dimension(NSOA)        :: foc
     real( dp), dimension(NSOA)        :: hvap
     real( dp), dimension(NSOA)        :: csat0
     real( dp), dimension(NSOA)        :: csate
     real( dp), dimension(NSOA)        :: csoagas
     real( dp), dimension(NSOA)        :: casoa
! End SOA components

! Aqueous phase of particles
     real( dp), dimension(aqmax,APN)   :: caqold
     REAL( dp), dimension(APN)         :: radius
     REAL( dp), dimension(0:nspec)     :: alpha_T0
     REAL( dp), dimension(0:nspec)     :: alpha_Tdep
     REAL( dp), dimension(0:nspec)     :: Henry_T0
     REAL( dp), dimension(0:nspec)     :: Henry_Tdep
     REAL( dp), dimension(0:nspec)     :: molar_mass
     CHARACTER(STRLEN_KPPSPECIES),dimension(0:nspec) :: chem_name
     LOGICAL, dimension(APN)   :: loghet
     real( dp) :: zhetT
     real( dp) :: lwcm,phm     ! from ingeod.dat    
     real( dp) :: dlwc
     integer   :: jb
     integer   :: snom7nno
     logical   :: wascloud=.false.
! End aqueous phase     
          
    
!*************************************************************

! Integer and Char Variables
     integer               :: I,M,K,Q,S
     integer               :: fx,jk,jl,jn,jm,jo,jq

!*************************************************************

     integer,dimension(12) :: imn=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
! fix for tab reading
     character(len=180)    :: inchemline
     character(len=1)      :: tab = ACHAR(9)
     integer               :: it
! end fix
     character(len=20)     :: fb2='(2X,e12.4)'
     character(len=10)     :: fb1
     character(len=40)     :: fb
     character(len=20)     :: fc2='(3X,e12.4)'
     character(len=10)     :: fc1
     character(len=40)     :: fc
     !character(len=10)     :: fc3='(3X,a10)'
     !character(len=40)     :: fa
     logical               :: ndis=.true.
     logical               :: firstloop=.true.

!*************************************************************
!*************************************************************
! OUTPUT SECTION

     OPEN(UNIT=12,FILE="debug.res",   STATUS="REPLACE")
     OPEN(UNIT=13,FILE="total_n.res", STATUS="REPLACE")
     OPEN(UNIT=14,FILE="size_dis.res",STATUS="REPLACE")
     OPEN(UNIT=15,FILE="size_dism.res",STATUS="REPLACE")
     OPEN(UNIT=16,FILE="concout.res", STATUS="REPLACE")
     OPEN(UNIT=17,FILE="aerconc.res", STATUS="REPLACE")
     OPEN(UNIT=18,FILE="wetdp.res",   STATUS="REPLACE")
     OPEN(UNIT=19,FILE="plume.res",   STATUS="REPLACE")
     OPEN(UNIT=11,FILE="soadis.res",  STATUS="REPLACE")

!*************************************************************
!*************************************************************
! INPUT SECTION

!   Default values and units
! METEO
    temp     = 298.15    ! temperature [K]
    press    = 101300.0  ! pressure [Pa]
    RH       = 0.50      ! relative humidty [rh frac]
    zmbl     = 1000.     ! mixing height [m]
    vupdra   = 0.2       ! updraft velocity [m/s]
    supersat = -0.02     ! supersaturation (RH=0.98)

! GENERAL
    incloud=0            ! in-cloud flag
    hour_time=0.         ! hour counter [s]
    hour_timec=0.        ! 1min counter [s]
    hour_timet=0.        ! 1min counter [s]
    hour_timeo=0.        ! output time counter [s]
    hour_timemin=0.      ! display time counter [s]
    
! CHAMBER
    IOZ      = 0         ! read O3 from monitor.dat
    F_HONO   = 0.        ! HONO chamber source
    alphanit = 1.        ! nitrate condensation (accom. coeff.)
    fco      = 1.        ! nitrate condensation
    fcoj     = 1.        ! nitrate condensation
    fkoh_mea = 1.        ! k(AMIN+OH) scaling
   
! ORGANIC VAPOURS
    surfin   = 1         ! surface tension organics (input/function)
    surf_org = 0.050     ! surface tension organics  [kg/s^2]
    
! PLUME DISPERSION
!  default and input values
!  Pohjola et al. Atm. Environ., 37, 339-351,2003
    dila     = 86.49     ! dispersion parameter A   [-]
    dilcoef  = 0.92332   ! dilution coefficient     [-]
    dilrate  = 0.0       ! dilution rate            [1/s]
    hmix_st  = 1.50      ! mixing height @ station  [m]
    dst_st   = 12.0      ! distance @ station       [m]
    hsta_st  = 10.0      ! stack height             [m]
    ta_st    = 400.0     ! temperature @ station    [K]
    BGH2O    = 1.E16     ! [molecules cm^-3]
    emisratp = 1.0       ! emission ratio line1/line2
!  plume height, width and temperature    
    temp_old = temp
    zmbl_old = zmbl
    wplm_old = 3.5
    wplm     = wplm_old
    dilut_time=0.        ! dilution time counter [s]
    dilstore = 0.

!   PN parameterization
    vd1      = 0.00528    ! [m/s]
    !vd2      = 0.00117    ! [m/s] PNC1:3
    !vd3      = 0.00024    ! [m/s] PNC1:3
    vd2      = 0.00181    ! [m/s] PNC1:6
    vd3      = 0.00068    ! [m/s] PNC1:6
    vd4      = 0.00039    ! [m/s] PNC1:6
    vd5      = 0.00023    ! [m/s] PNC1:6
    vd6      = 0.00024    ! [m/s] PNC1:6
    kcoa1    = 4.51E-15   ! [m^3/s]
    !kcoa2    = 3.10E-15   ! [m^3/s] PNC1:3
    !kcoa3    = 8.82E-16   ! [m^3/s] PNC1:3
    kcoa2    = 1.50E-14   ! [m^3/s] PNC1:6
    kcoa3    = 5.40E-15   ! [m^3/s] PNC1:6
    kcoa4    = 6.26E-15   ! [m^3/s] PNC1:6
    kcoa5    = 2.28E-15   ! [m^3/s] PNC1:6
    kcoa6    = 8.69E-16   ! [m^3/s] PNC1:6
                
! NUCLEATION
    natot    = 2.        ! number of molec in crit. cluster
    fnuc     = 1.        ! scaling factor of nucleation constant

!*************************************************************
!*************************************************************

! DRIVER
   !!! Read General input for simulation
   !!! incloud has to be 0 or 1 for the full run
      open(8,file='ingeod.dat',status='old', iostat=stat1)
! open error handling    
      if (stat1.ne.0) then
        write(6,*) 'File ingeod.dat cannot be opened !'
        stop
      endif
    
! fix for tab reading
      read(8,'(a)') inchemline
      do it=1,LEN(inchemline)
      ! replace each tab in the line with a space
        if( inchemline(it:it).EQ.tab) inchemline(it:it)=' '
      end do
      read(inchemline,*) runtime, iday, imonth, starttime, lat_deg,       &
          lon_deg, temp, press, RH, zmbl, incloud, u10, rain, edms, eso2, &
          eh2o2,cnh3,camidoh,cdms,co3,fnuc,lwcm,phm,dila,dilcoef
! end fix for tab reading
    
      if (RH .ge. rhactiv) then
        write(6,*) 'WARNING: RH >=0.99; cloud droplet activation'
      endif
      if (RH .gt. 1.10) then
        write(6,*) 'STOP: RH>1.1 in ingeod.dat'
        stop
      endif

      TAIR = temp

! TIMES
     ! parse date to find day number (julian day)
      mday = 0
      do  month = 1, imonth-1
         mday = mday + imn(month)
      end do
      jday = mday + iday
      daynr = float(jday - 1)
      write(6,*) 'DOY',daynr

! PROCESSES
   !!!  Read sensitivities file
   !!!  Unit 9 is reserved for this
   !!!  Specify, which aerosol processes are considered
      call readsens()

! SIMULATION TIME STEP
      DTIME = 5.0_dp
      if ((IDIL==1).or.(IDIL==2).or.(IDIL>4))  DTIME=0.1_dp
      if ((IDIL==3).or.(IDIL==4))   DTIME=0.01_dp
      model_time=starttime*3600._dp          !start time [s]
      TENDE  = model_time + runtime*3600._dp + DTIME !final time [s]
      if (IDIL==3)   TENDE=model_time + 60._dp + DTIME
      if (IDIL==4)   TENDE=model_time + 120._dp + DTIME

! GAS CHEMISTRY

     !!! Read Chemistry input for simulation
     !!! Initialize all tracers
     !!! emis in (cm^-2s^-1), vdry in (cms^-1)
      emis(:) = 0._dp
      vdry(:) = 0._dp
      c(:)    = 0._dp

     ! read kpp ctrl namelist (mafor.nml)
     !    CALL initialize_kpp_ctrl(status, iou, modstr)
     ! icntrl(3) = solver-specific method:
     ! default value 0
     !    icntrl(3) = 1 ! ros2: L-stable method, 2 stages, order 2
     !    icntrl(3) = 2 ! ros3: L-stable method, 3 stages, order 3 (recommended)
     !    icntrl(3) = 3 ! ros4: L-stable method, 4 stages, order 4
     !    icntrl(3) = 4 ! rodas3: stiffly-stable method, 4 stages, order 3
     !    icntrl(3) = 5 ! rodas4: stiffly-stable method, 6 stages, order 4
      icntrl(3) = 2

      CFACTOR = 1.000000e+00_dp
      DO i = 1, NVAR
        VAR(i) = (0.)*CFACTOR
      END DO

     !!! input concentrations, emission and deposition      
      IAM     = 1     ! by default use MEA

     !!! read inchem.dat
      call readchem(emis,vdry,c)
      emis0(:)=emis(:)


     !!! water vapour pressure and concentration
      ! cair: c(air) in [mcl/cc]
      cair    = (N_A/1.E6_dp) * press / (R_gas*temp)
      c(ind_H2O) = conch2o(temp,RH,press,cair)
      ! fixed concentrations: 1=O2, 2=N2, 3=CO2, 4=H2Oas
      ! in molec/cm^3
      FIX(1)   = 0.209460*2.46e19*298./temp
      FIX(2)   = 0.780840*2.46E19*298./temp
      FIX(3)   = 355.e-6*2.46E19*298./temp
      c(ind_O2) = FIX(1)
      c(ind_N2) = FIX(2)
      c(ind_CO2)= FIX(3)

! AQ. PHASE CHEMISTRY

     !!! Initialize indexarrays for aq. phase species
      call initialize_indexarrays


! SIMULATION INPUT

    ! INPUT CONDENSING ORGANIC VAPOURS
     !!! returns density of organic and soot particles
      call readorganic(DEN(OC),DEN(EC),M_oc,nmo,foc,hvap,csat0)

      if (ICHAM==1) then
    ! INPUT CHAMBER AND MONITOR DATA
        call readchamber()
        call readmonit(cco3,cno,cno2,cipn,temp,jno2m)
        ! convert ppb to molec cm-3
        IF (IOZ.EQ.1)  c(ind_O3)=cco3*(cair*1.E-9)
        c(ind_NO)  = cno* (cair*1.E-9)
        c(ind_NO2) = cno2*(cair*1.E-9)
        c(ind_IPN) = cipn*(cair*1.E-9)
      else    
    ! INPUT DISPERSION AND DEPOSITION DATA
        call readdispers()
!PLUME
        if (IDIL>0) then
        ! initial plume height, width and temperature
          call initplume(u10,temp_old,wplm,zmbl)
        ! call calculate plume area (for output in plume.res)
          call plumearea(wplm,zmbl,areapl) 
          wplm_old=wplm
          zmbl_old=zmbl
        endif
! PLUME
    ! INPUT PRESCRIBED CONCENTRATIONS
        ! from ingeod.dat
        if (IORG==1) then
          c(ind_BSOV)=camidoh   ! prescribe "COV"
        endif
        if (IDMS==1) then
          c(ind_DMS)=cdms
          c(ind_NH3)=cnh3
        endif
        if (IOZO .EQ. 1) c(ind_O3)=co3
      endif

  

! AEROSOL INPUT

   !!! Read Aerosol input
   !!! The modal structure is
   !!! 1) 1-10nm, 10-100nm, 100-1000nm, 1000-10000nm or
   !!! 2) 1-10nm, 10-50nm, 50-100nm, 100-1000nm 
   !!!  
      call readinaero(IMAX,DPMAX,ndis,MSULFTOT,MORGCTOT,         &
                      MAMMOTOT,MNITRTOT,MMSAPTOT,MSALTTOT,       &
                      MXXXXTOT,MECBCTOT,MDUSTTOT)



    !!! read Particle emission input
    !!! PBA emission rate (ng m^-2 s^-1)
      if (IPEMI.eq.1) then
        call reademitpar(EMMCTOT)
      endif


    !!! Read Background Aerosol input
      if (IDIL>0) then
        call readbackgr(BGMCTOT)
      endif

 

! ALLOCATE DYNAMIC ARRAYS

    allocate(ROOP(MMAX,IMAX),ROOPW(MMAX,IMAX), STAT=IAS1)
    allocate(CCOND(MMAX,IMAX,QMAX),  STAT=IAS61)
    allocate(NVAP(QMAX),NVAPO(QMAX),  STAT=IAS63)
    allocate(DLOGDP(MMAX,IMAX),DLINDP(MMAX,IMAX),N(MMAX,IMAX),LINDIS(MMAX,IMAX),LOGDIS(MMAX,IMAX), STAT=IAS7)
    allocate(VPT(MMAX,IMAX),MPT(MMAX,IMAX),MPTW(MMAX,IMAX),MLOGDIS(MMAX,IMAX), STAT=IAS81)
    allocate(MASS(MMAX,IMAX,AMAX), STAT=IAS9)
    allocate(DPA(MMAX,0:(IMAX+1)), DPAW(MMAX,0:(IMAX+1)), STAT=IAS10)
    allocate(IAG(MMAX,MMAX,IMAX,IMAX), STAT=IAS12)
    allocate(BGLOGDIS(MMAX,IMAX),BGN(MMAX,IMAX),DIFFDILN(MMAX,IMAX), STAT=IAS13)
    allocate(DIFFDILM(MMAX,IMAX,AMAX), STAT=IAS14)
    allocate(BGMASS(MMAX,IMAX,AMAX), STAT=IAS15)   
    allocate(EMASS(MMAX,IMAX,AMAX),EN(MMAX,IMAX), STAT=IAS16)
    allocate(saltemis(MMAX,IMAX), STAT=IAS6)
!MESA start
    allocate(henry_eff_hno3(mmax,imax),henry_eff_hcl(mmax,imax),henry_eff_nh3(mmax,imax), STAT=IAS17)
    allocate(svmc_min_hno3(mmax,imax),svmc_min_hcl(mmax,imax),svmc_min_nh3(mmax,imax),   STAT=IAS18)
    allocate(flag_dissolution(mmax,imax),                        STAT=IAS19)       
    allocate(hyst(mmax,imax),                                    STAT=IAS20)
    allocate(cioncharge(mmax,imax), kh_nh3(mmax,imax),           STAT=IAS21) 
!MESA end
!CCN start
    allocate(ddpwdt(mmax,imax),                                  STAT=IAS22)
    allocate(dmwdt(mmax,imax),                                   STAT=IAS23)
!CCN end


!*************************************************************
!*************************************************************
! INITIALIZE AEROSOL
!
! marine aerosol: accumulation and coarse mode [Jaenicke 1993]
! based on dry radius
!   ndis  : flag for number or mass conc distribution
!   GMD   : geometric mean diameter of mode [m]
!   DPMIN : diameter of smallest bin [m]
!   DPMAX : diameter of largest bin [m]
!   SIG   : geometric standard deviation of mode
!   NUM   : initial number concentration [part/m^3] !not used
!   MSULFTOT : H2SO4 particle conc [ng/m^3]
!   MORGCTOT : Organic particle conc [ng/m^3]
!   MAMMOTOT : NH4 particle conc [ng/m^3]
!   MNITRTOT : NO3 particle conc [ng/m^3]
!   MMSAPTOT : MSAp particle conc [ng/m^3]
!   MSALTTOT : Sea salt particle conc [ng/m^3]
!   MXXXXTOT : Biological particle conc [ng/m^3]
!   MECBCTOT : Elemental carbon particle conc [ng/m^3]
!   MDUSTTOT : Mineral dust particle conc [ng/m^3]


    ! total number particles [partic/m3]
    NTOT(:)   = 0._dp
    ENTOT(:)  = 0._dp
    BGNTOT(:) = 0._dp

    ! particle density [kg/m3]; later calculate for mixture
    DEN(SU)=DENV
    !DEN(OC) from readorganic
    DEN(AM)=DENAM
    DEN(NI)=DENNI
    DEN(MS)=DENMS
    DEN(SA)=DENSA
    DEN(XX)=DENXX   !PBA
    !DEN(EC) from readorganic
    DEN(DU)=DENDU
    DEN(WA)=RHOH2O

    MASS(:,:,:)       = 0._dp      ! mass conc. of compound in bin [ng/m^3]    
    N(NU:CS,1:IMAX)   = 0._dp      ! number conc. in bin     [partic/m^3]
    VPT(NU:CS,1:IMAX) = 0._dp      ! particle volume in bin  [m^3]
    MH2OTOT(NU:CS)    = 0._dp      ! H2O particle conc       [ng/m^3]
    DPMIN = 1.00E-9_dp             ! min. diam. of the size distribution [m]
    VVAP=MVAP/DENV                 ! vapor molecule volume [m^3]
    GRTOT=0._dp                    ! growth rate
    CSORGT=0._dp                   ! condensation sink organic [1/s]
    CSSULT=0._dp                   ! condensation sink sulfate [1/s]
    COAGS=0._dp                    ! coagulation sink [1/s]
    DNVAPC=0._dp                   ! amine loss to particles  
       
    ! condensing compound (vapor)
    NVAP(:) = 0._dp
    ! initial H2SO4 vapor conc.   [1/m^3]
    NVAP(A_SUL)=c(ind_H2SO4)*1.e6_dp
    ! initial organic vapor conc. [1/m^3]
    NVAP(A_OR1)=c(ind_BSOV) *1.e6_dp
    NVAP(A_OR2)=c(ind_BLOV) *1.e6_dp
    NVAP(A_OR3)=c(ind_BELV) *1.e6_dp
    NVAP(A_OR4)=c(ind_ASOV) *1.e6_dp
    NVAP(A_OR5)=c(ind_ALOV) *1.e6_dp
    NVAP(A_OR6)=c(ind_AELV) *1.e6_dp
    NVAP(A_OR7)=c(ind_PIOV) *1.e6_dp
    NVAP(A_OR8)=c(ind_PSOV) *1.e6_dp
    NVAP(A_OR9)=c(ind_PELV) *1.e6_dp
    ! for SOA output, SOA_gas     [ug/m^3]
    csoagas(1)=NVAP(A_OR1)*1.e-6*molec2ug(M_oc(1))
    csoagas(2)=NVAP(A_OR2)*1.e-6*molec2ug(M_oc(2))
    csoagas(3)=NVAP(A_OR3)*1.e-6*molec2ug(M_oc(3))
    csoagas(4)=NVAP(A_OR4)*1.e-6*molec2ug(M_oc(4))
    csoagas(5)=NVAP(A_OR5)*1.e-6*molec2ug(M_oc(5))
    csoagas(6)=NVAP(A_OR6)*1.e-6*molec2ug(M_oc(6))
    csoagas(7)=NVAP(A_OR7)*1.e-6*molec2ug(M_oc(7))
    csoagas(8)=NVAP(A_OR8)*1.e-6*molec2ug(M_oc(8))
    csoagas(9)=NVAP(A_OR9)*1.e-6*molec2ug(M_oc(9))
    ! for SOA output, Csat0(T)    [ug/m^3]
    do S=1,NSOA
      csate(S) = csat0(S)*EXP((hvap(S)/R_gas) *     &
                ((1._dp/298._dp)-(1._dp/temp)))
    enddo

    if (IDIL==3)  DPMIN = 1.50E-9_dp 

! total mass conc. in each mode      
    do M=NU,CS
     ppmtot(M) = MSULFTOT(M) + MORGCTOT(M) + MAMMOTOT(M) + MNITRTOT(M) + &
                 MMSAPTOT(M) + MSALTTOT(M) + MXXXXTOT(M) + MECBCTOT(M) + &
                 MDUSTTOT(M)
    enddo

! BIN SIZE DISTRIBUTION
   !!! initialise bin size distribution with one lognormal mode
   !!! DPA is dry diameter [m]
   !!! VPT is dry volume   [m^3]
   !!! DLOGDP logarithmic width, DLINDP linear width of size bin

    CALL initSizeDistribution(IMAX,DPMIN,DPMAX,DPA,VPT,DLOGDP,DLINDP)


! INITIAL NUMBER AND MASS IN SIZE BINS

   !!! initial number and mass in bin
   !!! M.J. Jacobson, p.457, eq. 13.24
   !!! this means: GMD is geometric-mean number diameter
   !!! we assume average particle density (DEN) to be
   !!! constant in all bins.
   !!! Calculate mass conc of compoound in each bin I
   !!! from log-normal mass distribution
   !!! thus GMD is geometric-mean mass diameter!!!

    CALL initNumberMass(IMAX,DPA,DLINDP,VPT,DEN,                    &
                            MSULFTOT, MMSAPTOT, MNITRTOT, MAMMOTOT, &
                            MSALTTOT, MECBCTOT, MDUSTTOT, MXXXXTOT, & 
                            MORGCTOT, MASS, N)


! aerosol mass concentration from input
    do M=NU,CS
     mtotin(M) = 0._dp
     do I=1,IMAX


        mtotin(M) = mtotin(M)       + MASS(M,I,A_SUL)                   + &  
                    MASS(M,I,A_OR1) + MASS(M,I,A_OR2) + MASS(M,I,A_OR3) + &
                    MASS(M,I,A_OR4) + MASS(M,I,A_OR5) + MASS(M,I,A_OR6) + &
                    MASS(M,I,A_OR7) + MASS(M,I,A_OR8) + MASS(M,I,A_OR9) + &
                    MASS(M,I,A_AMI) + MASS(M,I,A_NH4) + MASS(M,I,A_NIT) + &
                    MASS(M,I,A_MSA) + MASS(M,I,A_SAL) + MASS(M,I,A_CHL) + &
                    MASS(M,I,A_XXX) + MASS(M,I,A_EBC) + MASS(M,I,A_DUS)

      enddo
    enddo


! Background aerosol
    IF (IDIL>0) THEN

          call initBGNumberMass(IMAX,DPA,DLINDP,VPT,DEN,            &
                                BGMCTOT, BGMASS, BGN)

          do M=NU,CS
            do I=1,IMAX
              BGNTOT(M)=BGNTOT(M)+BGN(M,I)
              BGLOGDIS(M,I)=BGN(M,I)/DLOGDP(M,I)
            enddo
          enddo

    ENDIF


! Emitted aerosol
    IF (IPEMI .EQ. 1) THEN

          call initEMNumberMass(IMAX,DPA,DLINDP,VPT,DEN,            &
                                EMMCTOT,  EMASS,  EN)

    ENDIF


! CORRECTION OF MASS AND NUMBER TO MAINTAIN INITIAL MASS (each mode)
    if(ndis) then
      do M=NU,CS
        do I=1,IMAX
          do K=1,AMAX-1
            MASS(M,I,K)=MASS(M,I,K)*(ppmtot(M)/max(mtotin(M),massmin))
          end do
          N(M,I) = N(M,I)*(ppmtot(M)/max(mtotin(M),massmin))
        end do
      end do
    endif    

! some output diagnostics        
    do M=NU,CS
      do I=1,IMAX
        if (DPA(M,I).gt.3e-9) THEN        
         NTOT(M)=NTOT(M)+N(M,I)
        endif
        if (IPEMI .EQ. 1)  ENTOT(M)=ENTOT(M)+EN(M,I)
        LOGDIS(M,I)=N(M,I)/DLOGDP(M,I)     
      end do
      if(M.eq.NU) then
        write( 6,fmt='(a,i2,a,f13.5)' ) '    Mode', M,' init number conc. (>3nm)   [#/cm^3]',NTOT(M)*1.e-6
      else
        write( 6,fmt='(a,i2,a,f13.5)' ) '    Mode', M,' init number conc.          [#/cm^3]',NTOT(M)*1.e-6
      endif
      if (IDEB == 1) then
        if(M.eq.NU) then
          write(12,fmt='(a,i2,a,f13.5)' ) '    Mode', M,' init number conc. (>3nm)   [#/cm^3]',NTOT(M)*1.e-6
        else
          write(12,fmt='(a,i2,a,f13.5)' ) '    Mode', M,' init number conc.          [#/cm^3]',NTOT(M)*1.e-6
        endif
        if (IPEMI==1) write(12, '(A34, I2, E12.2)') 'number flux (# m^-2 s^-1) in mode',M,ENTOT(M)
        if (IDIL>0)   write(12,*) 'init bg number (# cm^-3) in mode',M,BGNTOT(M)*1.e-6
      endif
    end do


! WATER ACTIVITY
    ! Water activity from sulphate - ammonium - nitrate. MH2O in ng/m^3
    ! Munkelwitz and Tang 1994 / Seinfeld and Pandis 1997
    ! Sea salt hygroscopicity same  
    ! If in-cloud, MH2O and N of CS mode is calculated from LWC(CS)
    ! total water per mode in molec/cm^3
    ! convert MASS(M,I,A_WAT) into c(ind_H2O_a*) in molec/cm^3(air)
    ! and  liquid water content [m3/m3] lwc(as)

    MASS(NU,:,A_WAT)=0._dp
    DPAW(:,:)=DPA(:,:)
    lwc(:)=0._dp
    dmwdt(:,:)=0._dp


    CALL water_content(incloud,firstloop,IMAX,RH,lwcm,MASS,N,     &
             dmwdt,DPA,GMD,SIG,DLINDP,VPT,DPAW,MH2OTOT,NTOT(CS),  &
             lwc,c(ind_H2O_a(1)),xaer,wascloud  )



! TOTAL MASS CONCENTRATION
! Total dry mass and wet mass in kg/m^3
! Total mode masses in ng/m^3
    CALL getTotalMass(IMAX,MASS,MPT,MPTW,                         &
                      MTOT,MTOTW,MSULFTOT,MMSAPTOT,               &
                      MNITRTOT,MAMMOTOT,MSALTTOT,MECBCTOT,        &
                      MDUSTTOT,MXXXXTOT,MORGCTOT,MORG1TOT,        &
                      MORG2TOT,MORG3TOT,MORG4TOT,MORG5TOT,        &
                      MORG6TOT,MORG7TOT,MORG8TOT,MORG9TOT,        &
                      casoa   )

    do M=NU,CS
      do I=1,IMAX
        MLOGDIS(M,I)=MPTW(M,I)/DLOGDP(M,I)
      enddo
    enddo 


! DIAMETER OF WET PARTICLE
! only for output dpwet.res
! Wet diameter/volume is used to calculate aerosol processes
      if (ICONW .ge. 1) then
        call wetdiameter(IMAX,incloud,RH,MTOT,MTOTW,DPA,DPAW)
        do M=NU,CS
          do I=1,IMAX
            VPT(M,I) =(pi/6.)*DPAW(M,I)**3.  ! only firstloop              
          end do
        end do
      else
        do M=NU,CS
          do I=1,IMAX
            DPAW(M,I)=DPA(M,I)
          end do
        end do      
      endif

! AVERAGE DENSITY
! density in kg/m^3

      call getDensity(IMAX,temp,DEN(EC),DPA,MASS,MPT,MPTW,ROOP,ROOPW)

! COAGULATION COEFFICIENTS
      if (ICOAG .ge. 1) then
        call coagulation_target(VPT,IMAX,IAG)
      endif


! NUCLEATION
    IF (INUC .EQ. 1) CALL nucleationratio(DPA(NU,1),VVAP,NNUC)
    JNUC=0._dp
    N3=0._dp
    NEU1=NTOT(NU)
    ! NEU2: 25-100m particles
    NEU2=0._dp
    NEU3=0._dp
    NEU4=0._dp
    NEU5=0._dp
    NEU6=0._dp
    NEU7=0._dp    
    do I=1,IMAX
      if (DPA(NU,I).GT.3.E-9_dp) N3=N3+N(NU,I)
      if (DPA(AI,I).LT.25.E-9_dp) NEU1=NEU1+N(AI,I)
      if (DPA(AI,I).GE.25.E-9_dp) NEU2=NEU2+N(AI,I)
      if (DPA(AS,I).LT.1.E-7_dp)  NEU2=NEU2+N(AS,I)
      if ((DPA(AI,I).GE.24.E-9_dp).and.(DPA(AI,I).LT.50.E-9_dp)) NEU3=NEU3+N(AI,I)
      if ((DPA(AS,I).GE.45.E-9_dp).and.(DPA(AS,I).LT.50.E-9_dp)) NEU3=NEU3+N(AS,I)
      if ((DPA(AI,I).GE.50.E-9_dp).and.(DPA(AI,I).LT.75.E-9_dp)) NEU4=NEU4+N(AI,I)
      if ((DPA(AS,I).GE.50.E-9_dp).and.(DPA(AS,I).LT.75.E-9_dp)) NEU4=NEU4+N(AS,I)
      if (DPA(AI,I).GE.75.E-9_dp) NEU5=NEU5+N(AI,I)
      if (DPA(AS,I).GE.75.E-9_dp) NEU5=NEU5+N(AS,I)
      if (DPA(CS,I).LT.2.E-7_dp)  NEU6=NEU6+N(CS,I)
      if (DPA(CS,I).GE.2.E-7_dp)  NEU7=NEU7+N(CS,I)
    end do

! PN PARAMETERIZATION
    PNC1 = NEU1 - NTOT(NU)
    !PNC2 = NEU2             !PNC1:3
    !PNC3 = NTOT(CS)         !PNC1:3
    PNC2 = NEU3
    PNC3 = NEU4
    PNC4 = NEU5
    PNC5 = NEU6
    PNC6 = NEU7            
    zmblm=sqrt((hmix_st)**2+(dila*(1800.*u10*1.e-3)**dilcoef)**2)

!*************************************************************
!*************************************************************
! AQUEOUS PHASE CHEMISTRY

!Aqueous phase initialisation start
    !!! initialize gas phase to aerosol transfer
    CALL aero_init_gasaq(chem_name,Henry_T0,Henry_Tdep,alpha_T0,alpha_Tdep,molar_mass)
    !do jl=1,NSPEC
    !  write(6,*) 'henry',jl,chem_name(jl),Henry_T0(jl),Henry_Tdep(jl),alpha_T0(jl), &
    !                      alpha_Tdep(jl),molar_mass(jl)
    !end do   


    !!! initialize transfer reactions
    !!! forward + backward transfer gas/aq
    k_exf(:,:) = 0._dp
    k_exb(:,:) = 0._dp

    !!! default: no aqueous phase reactions
    xaer(:)=0.0
    
    !!! default: no partitioning gas/aq
    !!! switch 1/0 for H2SO4 and MSA gas/aq equilibrium
    !!! switch 1/0 for SO2, DMS and DMSO heterogeneous uptake
    !!! switch 1/0 for OH and HO2 and O2,O3 gas/aq equilibrium
    !!! switch 1/0 for NO2/NO3/N2O5/HONO/HNO3/HNO4/NH3 liquid uptake
    !!! at current: exchange of O3 and SO2 back to gas phase blocked
    xnom7sulf = 0
    xnom7msap = 0
    xnom7so2  = 0
    xnom7dmso = 0
    xnom7hox  = 0
    xnom7ox   = 0
    xnom7h2o2 = 0
    xnom7nox  = 0
    xnom7co2  = 0
    xnom7cl   = 0
    !!! CCM amine switches
    xnom7amin = 0
    xnom7nno  = 0
    xaqcamin  = 0
    xaqcnno   = 0
    xgcnno    = 1
    xgcamin   = 1
    xaeq(:)   = 0.0
    snom7nno  = 0

    !!! aqueous phase chemistry on
    IF (IAQC.EQ.1) THEN
      xaer(:)  = 1.0
      ya_soan1 = 0.0
      ya_soan2 = 0.0 
      ya_soan5 = 0.0
      xaqcamin = 1
      xaqcnno  = 1
    END IF

    !!! aqueous phase partitioning on
    IF (IAQP.EQ.1) THEN
      ya_soan1 = 0.0
      ya_soan2 = 0.0 
      ya_soan5 = 0.0
      xnom7so2  = 1
      xnom7hox  = 1
      xnom7ox   = 1
      xnom7h2o2 = 1
      xnom7sulf = 1
      xnom7msap = 1
      xnom7dmso = 0
      xnom7nox  = 1
      xnom7co2  = 1
      xnom7cl   = 1
      xnom7amin = 1
      xnom7nno  = 1
      snom7nno  = xnom7nno
      !!! do not allow heterogeneous chemistry at the moment
      do zkc=1,APN
        loghet(zkc)=.false.
      end do
    END IF

    !!! mean dry radius [m] needed when 
    !!! calculating mass transfer coeff.
    do M=NU,CS
      radiusm(M) = GMD(M)/2._dp
    enddo

    !!! aqueous modes: 1-3
    do jb=1,APN
      radius(jb) = radiusm(jb+1)
    end do

    if (IDEB == 1) then
      do jb=1,APN
         write(12,*) 'lwc',jb,lwc(jb),radius(jb)
      end do
    endif


    !!! cvfac: conversion factor dm^3(aq)/mol => cm^3(air)/molec
    do zkc=1,APN
      cvfac(zkc) = 1.E3 / ( N_A * lwc(zkc) )
    end do  

    !!! aerosol water is a fixed species in KPP
    FIX(4) = c(ind_H2O_a(1))
    !FIX(5) = c(ind_H2O_a(2))
    !FIX(6) = c(ind_H2O_a(3))

    !!! initialize c(:) for aqueous phase
    !!! read Initial aqueous phase concentrations
    ! aq. concentrations in coarse mode in molec/cm^3

    if ((IAQC.eq.1).OR.(IAQP.eq.1)) then

        call readaqphase(c)
      ! coarse mode pH
        if (incloud.eq.1) then
          c(ind_Hp_a(1))=10**((-1)*phm)   ![mol/dm^3]
          c(ind_Hp_a(1))=c(ind_Hp_a(1))*N_A*lwc(1)/1000.
        end if

    endif

!Aqueous phase initialisation end

!MESA start
! initialize thermodynamic solver
    NVAP(A_NH4)=c(ind_NH3) *1.e6_dp
    NVAP(A_NIT)=c(ind_HNO3)*1.e6_dp
    NVAP(A_CHL)=c(ind_HCl) *1.e6_dp
    call interface_mosaic( IMAX,temp,press,RH, alphanit,  &
                fcoj,mass,nvap,mpt,dpaw,n,hyst,              &    
                henry_eff_hno3,henry_eff_hcl,henry_eff_nh3,  &
                svmc_min_hno3,svmc_min_hcl,svmc_min_nh3,     &
                cioncharge,kh_nh3,flag_dissolution,          &
                Keq_nh4no3_0,Keq_nh4cl_0)
!MESA end

!CCN start
     dtedt  = 0._dp
     dsupdt = 0._dp
     Nfa    = 0._dp
     Mfa    = 0._dp
     DPcrit = DPA(MMAX,IMAX)
     caqold(:,:)=0._dp
     supersat = RH - 1._dp
!CCN end



!*************************************************************
!*************************************************************
! INITIALIZE OUTPUT

    daytime = daynr + (model_time/OneDay)

    ! formating size distribution output
    write(fb1,'(I6)') 4*IMAX+1
    fb="("//trim(adjustl(fb1))//fb2//")"
    ! formating concentration output
    write(fc1,'(I6)') NSPEC+1+APN+1
    fc="("//trim(adjustl(fc1))//fc2//")"
    !fa="("//trim(adjustl(fc1))//fc3//")"
        
    ! write first time step (t=0)
    write(13,'(29e17.6)') model_time,N3,NTOT(NU),NTOT(AI),NTOT(AS),NTOT(CS),    &
            NEU1,NEU2,DNVAPC,COAGS,CSORGT,GRTOT,JNUC,NEU3,NEU4,NEU5,NEU6,NEU7,  &
            PNC1,PNC2,PNC3,PNC4,PNC5,PNC6,                                      &
            supersat,Nfa,Mfa,DPcrit,temp
    !write(16,fa) 'time',(chem_name(jl),jl=1,nspec),'lwc1','lwc2','lwc3','lwc4','jr_no2'        
    write(16,fc) model_time,(c(jl),jl=1,nspec),(lwc(zkc),zkc=1,APN),0.00
    write(14,fb) model_time, (DPA(NU,I), I=1,IMAX),(DPA(AI,I), I=1,IMAX),       &
                 (DPA(AS,I), I=1,IMAX),(DPA(CS,I), I=1,IMAX)                          
    write(14,fb) model_time, (DLOGDP(NU,I), I=1,IMAX),(DLOGDP(AI,I), I=1,IMAX), &
                 (DLOGDP(AS,I), I=1,IMAX),(DLOGDP(CS,I), I=1,IMAX)
    if (IDIL>0) then
      write(14,fb) model_time,(BGLOGDIS(NU,I),I=1,IMAX),(BGLOGDIS(AI,I),I=1,IMAX), &
                 (BGLOGDIS(AS,I),I=1,IMAX),(BGLOGDIS(CS,I),I=1,IMAX)
    endif
    write(14,fb) model_time, (LOGDIS(NU,I), I=1,IMAX),(LOGDIS(AI,I), I=1,IMAX), &
                 (LOGDIS(AS,I), I=1,IMAX),(LOGDIS(CS,I), I=1,IMAX)
    write(15,fb) model_time,(MLOGDIS(NU,I),I=1,IMAX),(MLOGDIS(AI,I), I=1,IMAX), &
                 (MLOGDIS(AS,I), I=1,IMAX),(MLOGDIS(CS,I), I=1,IMAX)
    if ((IDIL==3).or.(IDIL==4)) then
        do K=1,AMAX        
            write(15,fb) model_time,((MASS(NU,I,K)/DLOGDP(NU,I)),I=1,IMAX),   &
                ((MASS(AI,I,K)/DLOGDP(AI,I)),I=1,IMAX),                       &
                ((MASS(AS,I,K)/DLOGDP(AS,I)),I=1,IMAX),                       &
                ((MASS(CS,I,K)/DLOGDP(CS,I)),I=1,IMAX) 
        end do              
    endif        
    write(18,fb) model_time, (DPAW(NU,I), I=1,IMAX),(DPAW(AI,I), I=1,IMAX),       &
                 (DPAW(AS,I), I=1,IMAX),(DPAW(CS,I), I=1,IMAX)
    write(17,'(41e12.4)') model_time,                                 &
               MSULFTOT(NU),MSULFTOT(AI),MSULFTOT(AS),MSULFTOT(CS),   &
               MMSAPTOT(NU),MMSAPTOT(AI),MMSAPTOT(AS),MMSAPTOT(CS),   &
               MXXXXTOT(NU),MXXXXTOT(AI),MXXXXTOT(AS),MXXXXTOT(CS),   &
               MORGCTOT(NU),MORGCTOT(AI),MORGCTOT(AS),MORGCTOT(CS),   &
               MAMMOTOT(NU),MAMMOTOT(AI),MAMMOTOT(AS),MAMMOTOT(CS),   &
               MNITRTOT(NU),MNITRTOT(AI),MNITRTOT(AS),MNITRTOT(CS),   &
               MECBCTOT(NU),MECBCTOT(AI),MECBCTOT(AS),MECBCTOT(CS),   &
               MDUSTTOT(NU),MDUSTTOT(AI),MDUSTTOT(AS),MDUSTTOT(CS),   &
               MSALTTOT(NU),MSALTTOT(AI),MSALTTOT(AS),MSALTTOT(CS),   &
               MH2OTOT(NU), MH2OTOT(AI), MH2OTOT(AS), MH2OTOT(CS)
! soa component output
    write(11,'(28e17.6)') model_time , (csate(S), S=1,NSOA),          &
               (csoagas(S), S=1,NSOA), (casoa(S), S=1,NSOA)
! plume dispersion output
    if (IDIL>0) then
      write(19,'(12e12.4)') model_time,dilut_time,temp_old,zmbl_old,  &
               wplm_old,areapl,dilrate,                               &
               NVAP(A_SUL)*1.e-6*molec2ug(M_H2SO4),                   &
               csoagas(7),csoagas(8),csoagas(9),                      &
               (ENTOT(NU)+ENTOT(AI)+ENTOT(AS)+ENTOT(CS))

               
    endif


    model_time=model_time+DTIME
    hour_time=hour_time+DTIME
    hour_timec=hour_timec+DTIME
    hour_timet=hour_timet+DTIME
    hour_timeo=hour_timeo+DTIME
    hour_timemin=hour_timemin+DTIME

!PLUME
    ! dilution time [s] between emission point and station
    ! dilution calculation starts 1s after emission point
    if (IDIL>0) then
      dilut_time=max((dst_st/u10),1.0_dp)
      if ((IDIL==3).or.(IDIL==4))  dilut_time=0.0+DTIME
    else
      dilut_time=0.0
    endif
!PLUME


!*************************************************************
!*************************************************************
! EVERY TIME STEP DO CHEMISTRY
!    Note: Aerosol exchange done after chemistry
! READ NEW METEO INPUT EVERY HOUR!!!

    call cpu_time(cpu_time_start)

    if (ICONW==2) write(6,*) 'MOSAIC coupling is activated'
    write(6,*) 'aerosol dynamics running ./. please wait...'  
     
    DO WHILE (model_time < TENDE)
 
       jx(:)=0._dp
! Initialise aerosol source terms
       do M=NU,CS
        do I=1,IMAX
          CCOND(M,I,:)  = 0._dp     
          saltemis(M,I) = 0._dp
        end do
       end do
       NVAPO(:) = 0._dp
       CSORGT   = 0._dp
       CSSULT   = 0._dp
       COAGS    = 0._dp
       JNUC     = 0._dp
       CTNH4    = 0._dp
       CTNIT    = 0._dp
       CTSO4    = 0._dp
       GRTOT    = 0._dp

!MESA start
    if (ICONW.eq.2) then
      do M=NU,CS
       do I=1,IMAX
         HYST(m,i)      = 0._dp
         do q=1,qmax
           CCOND(m,i,q) = 0._dp
         end do
       end do
      end do
    endif
!MESA end
       
! New meteo input every hour
! use: lat_deg,lon_deg,temp,press,RH,zmbl
      IF (hour_time .GT. 3600.1) THEN
        hour_time=0.
        read(8,*) runtime, iday, imonth, starttime, lat_deg,  &
                lon_deg,temp, press,RH,zmbl,incloud,u10,rain, &
                edms,eso2,eh2o2,cnh3,camidoh,cdms,co3,fnuc,   &
                lwcm,phm,dila,dilcoef
        emis(ind_DMS) = edms
        emis(ind_SO2) = eso2
        emis(ind_H2O2)= eh2o2
        TAIR = temp
        if (RH .ge. 0.99)   write(6,*) 'WARNING: RH >=0.99; program may stop'
        if (RH .gt. 1.10) then
          write(6,*) 'STOP: RH>1.1 in ingeod.dat'
          stop
        endif
      ENDIF

 ! PLUME DISPERSION
       if (IDIL>0) then

        ! time step for type 4 dilution is 0.5 s
        if (firstloop) dilstore=dilut_time
        if (dilstore+0.49.le.dilut_time) then
          dilstore = dilut_time
        endif

        call plumeDisp(DTIME,TAIR,temp_old,zmbl_old,wplm_old,     & 
                       u10,dilut_time,dilstore,dila,dilcoef,      &
                       zmbl,wplm,temp,dilrate,emisratp,  c    )

 ! Plume height not > BL height
        zmbl = min(zmbl,2000.0)

        ! call calculate plume area (for output in plume.res)
        call plumearea(wplm,zmbl,areapl) 

 ! Temperature correction, according to the Ideal Gas Law
 ! applied to aerosol particles in the plume
        do M=NU,CS
          do I=1,IMAX       
            do K=1,AMAX                  
              MASS(M,I,K)=MASS(M,I,K) *(temp_old/temp)
            end do
            N(M,I)=N(M,I) *(temp_old/temp)
          end do
        end do

        ! PN PARAMETERIZATION
        PNC1 = PNC1 *(temp_old/temp)
        PNC2 = PNC2 *(temp_old/temp)
        PNC3 = PNC3 *(temp_old/temp)
        PNC4 = PNC4 *(temp_old/temp)
        PNC5 = PNC5 *(temp_old/temp)
        PNC6 = PNC6 *(temp_old/temp)

        ! new temp is old temp of next time step
        temp_old=temp
        zmbl_old=zmbl
        wplm_old=wplm

      endif


        
! Prescribed gas concentrations
! For arctic nucleation runs
      IF (IDMS .EQ. 1) THEN
        c(ind_DMS)=cdms
        c(ind_NH3)=cnh3        
        if ((IDEB==1).and.(hour_timemin .EQ. 120.)) write(12,*) 'DMS',c(ind_DMS),'NH3',c(ind_NH3)
      ENDIF
      IF (IOZO .EQ. 1) THEN
        c(ind_O3)=co3
        if ((IDEB==1).and.(hour_timemin .EQ. 120.)) write(12,*) 'O3',c(ind_O3)
      ENDIF

! Particle emission from sea surface
! emissions rates [particles/(m2*s)]
      if (IPEMI .EQ. 2) then
        CALL seasaltemis(IMAX,DPA,u10,sst,sal,saltemis)
      endif

! Chamber dilution
      IF (ICHAM .EQ. 1) THEN
        CALL chamber_dil(DTIME,K_DIL,c)
        CALL chamber_loss(DTIME,temp,press,M_oc,csat0,hvap,     &
                    L_MEA,L_NO2,L_HNO3,L_O3,V_CHAM,S_CHAM,CWIN,c)
      ENDIF

! New monitor input every minute
      IF (ICHAM .EQ. 1) THEN
        IF (hour_timec .EQ. 60.) THEN
          hour_timec=0.
          call readmonit(cco3,cno,cno2,cipn,temp,jno2m)
        ENDIF
        ! convert ppb to molec cm-3
        IF (IOZ.EQ.1)  c(ind_O3)=cco3*(cair*1.E-9)
        c(ind_NO)  = cno* (cair*1.E-9)
        c(ind_NO2) = cno2*(cair*1.E-9)
        c(ind_IPN) = cipn*(cair*1.E-9)
      ELSE    
        !read organic vapor from ingeod.dat
        IF (IORG.EQ.1) THEN
          c(ind_BSOV)=camidoh
          if ((IDEB==1).and.(hour_timemin .EQ. 120.)) write(12,*) 'SOAN1',c(ind_BSOV)
        ENDIF  
      ENDIF

! Water pressure and concentration
! cair = c(air) in [mcl/cc]
       cair    = (N_A/1.E6_dp) * press / (R_gas*temp)
       ! H2O
       if (IDIL==3) then  
         c(ind_H2O) = c(ind_H2O) - ( dilrate*((c(ind_H2O)-BGH2O) *DTIME) )
       else
         c(ind_H2O) = conch2o(temp,RH,press,cair)
       endif


! Calculate new gas phase concentration after emission and dry deposition
       CALL emis_drydep(DTIME,emis,vdry,rain,zmbl,c)
!       do jl=1,nspec
!        write(6,*) 'emis_vd',jm,emis(jm),vdry(jm),c(jm)
!       end do

! Modulate gas emissions for plume dilution type 4
       if (IDIL==4) then
         emis(:)=emis0(:)*emisratp
       endif


! CLOUD MODULE
! Exchange between aerosol and gas phase
!   for all species (even those that don't partition)
! Calculate gas phase transfer coefficient to droplets [m s^-1]
! Calculate Henry coefficients [dimensionless]
! Forward/backward reaction rate coefficients
! LWC and pH are prescribed in "fog simulation"    
       ! Treatment of CS droplet mode vs particle mode
       if (incloud .eq. 1) then

          if (RH.ge.rhactiv) then
         ! LWC calculated elsewhere
          else
            supersat=RH-1._dp
            Nfa=0._dp
            Mfa=0._dp
            DPcrit=DPA(MMAX,IMAX)
         ! LWC in CS mode from ingeod.dat
         !   lwc(3)=lwcm
            lwc(1)=lwcm
          endif

!AQ partitioning start
          if (IAQP.eq.1) then

            if (RH.ge.rhactiv) then

              caqold(:,:)=0._dp
              if (hour_time.eq.tdiss) then
         ! Aerosol mass components to dissolve
         ! Do one time at beginning of simulation hour
         ! mXXXXtot(cs)   :: ng m^-3(air)
         ! ind_XXm_a(cs)  :: molec cm^-3(air)
               do zkc=1,APN
                 c(ind_SO4mm_a(zkc))=MSULFTOT(zkc+1)*1.e-15_dp*N_A/M_H2SO4
                 c(ind_HSO4m_a(zkc))=0._dp
                 c(ind_CH3SO3m_a(zkc))=MMSAPTOT(zkc+1)*1.e-15_dp*N_A/M_msa
                 c(ind_NO3m_a(zkc))=MNITRTOT(zkc+1)*1.e-15_dp*N_A/M_nit
                 c(ind_DMAp_a(zkc))=0._dp
                 c(ind_NH4p_a(zkc))=MAMMOTOT(zkc+1)*1.e-15_dp*N_A/M_nh3
                 c(ind_HC2O4m_a(zkc))=0._dp
                 c(ind_C2O4mm_a(zkc))=MORG1TOT(zkc+1)*1.e-15_dp*N_A/M_oc(1)
                 c(ind_C2H5C2O4m_a(zkc))=0._dp
                 c(ind_C2H4C2O4mm_a(zkc))=MORG2TOT(zkc+1)*1.e-15_dp*N_A/M_oc(2)
                 c(ind_Clm_a(zkc))=0.5*MSALTTOT(zkc+1)*1.e-15_dp*N_A/(MNa+MCl)
                 c(ind_Nap_a(zkc))=0.5*MSALTTOT(zkc+1)*1.e-15_dp*N_A/(MNa+MCl)
               enddo
               caqold(:,:)=0._dp
              endif
          ! Save old aqueous phase concentrations
              if (hour_time.gt.tdiss) then
                do zkc=1,APN
                  caqold(1,zkc)=c(ind_SO4mm_a(zkc))+c(ind_HSO4m_a(zkc))
                  caqold(2,zkc)=c(ind_CH3SO3m_a(zkc))
                  caqold(3,zkc)=c(ind_NO3m_a(zkc))
                  caqold(4,zkc)=c(ind_DMAp_a(zkc))
                  caqold(5,zkc)=c(ind_NH4p_a(zkc))
                  caqold(6,zkc)=c(ind_HC2O4m_a(zkc))+c(ind_C2O4mm_a(zkc))
                  caqold(7,zkc)=c(ind_C2H5C2O4m_a(zkc))+c(ind_C2H4C2O4mm_a(zkc))
                  caqold(8,zkc)=c(ind_Clm_a(zkc))+c(ind_Nap_a(zkc))
                enddo
              endif

            endif


         ! set water content and pH
            do zkc=1,APN
         ! prescribe pH from ingeod.dat
              c(ind_Hp_a(zkc))=10**((-1)*phm)   ![mol/dm^3]
              c(ind_Hp_a(zkc))=c(ind_Hp_a(zkc))*N_A*lwc(zkc)/1000.
              c(ind_OHm_a(zkc))=10**((-1)*(14-phm))   ![mol/dm^3]
              c(ind_OHm_a(zkc))=c(ind_OHm_a(zkc))*N_A*lwc(zkc)/1000.
         ! partitioning of gases only if LWC is large enough
              if (lwc(zkc).gt.1.e-10_dp) then
                 xaer(zkc)   = 1.0
              else
                 xaer(zkc)   = 0.0
              endif
            enddo


            xnom7nno  = snom7nno
         ! amine acid-base reactions in CS mode
            !xaeq(1)   = 0.0
            !xaeq(2)   = 0.0
            !xaeq(3)   = 1.0
            xaeq(1)   = 1.0

         ! Henry partitioning
            CALL transfer_coeff(radius,temp,press,xaer,lwc,Henry_T0,  &
                 Henry_Tdep,alpha_T0,alpha_Tdep,molar_mass,k_exf,k_exb)


            do zkc=1,APN

            ! heterogeneous chemistry (not included)
              if ((ind_H2O_a(zkc)/=0).AND.(loghet(zkc))) then
                zhetT = c(ind_H2O_a(zkc))
                ! Cl and Br not included in chemistry
                !if (ind_Clm_a(zkc)/=0) zhetT = zhetT + 5.E2 * c(ind_Clm_a(zkc))
                !if (ind_Brm_a(zkc)/=0) zhetT = zhetT + 3.E5 * c(ind_Brm_a(zkc))
                k_exf_N2O5(zkc)  = k_exf(zkc,IND_N2O5)  / zhetT
                !k_exf_ClNO3(zkc) = k_exf(zkc,IND_ClNO3) / zhetT
                !k_exf_BrNO3(zkc) = k_exf(zkc,IND_BrNO3) / zhetT
                k_exf_ClNO3(zkc) = 0._dp
                k_exf_BrNO3(zkc) = 0._dp
              else
                k_exf_N2O5(zkc)  = 0._dp
                k_exf_ClNO3(zkc) = 0._dp
                k_exf_BrNO3(zkc) = 0._dp
              endif


            ! update cvfac with new LWC
              cvfac(zkc) = 1.E3 / ( N_A * lwc(zkc) )

           end do  
           !do jl=1,nspec
           !    write(6,*) 'kext',jl,k_exf(CS,jl),k_exb(CS,jl)
           !end do
         endif

   !AQ partitioning end

       endif
! CLOUD MODULE END

! NOT INCLOUD CHEMISTRY
       if (incloud .eq. 0) then

         supersat=RH-1._dp
         Nfa=0._dp
         Mfa=0._dp
         DPcrit=DPA(MMAX,IMAX)


!MESA start
! Thermodynamics: new parameters for hno3 every 120 sec
! Stability and accuracy at long time step in accordance with
!    M.Z. Jacobson,  AEROSOL SCIENCE TECHNOLOGY, VOL. 39, 
!          NO. 2, 92-103, DOI: 10.1080/027868290904546, 2005.
!--------------------------------------------------------------
! CALL MESA INTERFACE EVERY 120 SEC
         if (hour_timet.eq.120.) then
         
             hour_timet=0.

             NVAP(A_NH4)=c(ind_NH3 )*1.e6_dp
             NVAP(A_NIT)=c(ind_HNO3)*1.e6_dp
             NVAP(A_CHL)=c(ind_HCl )*1.e6_dp
             call interface_mosaic( IMAX,temp,press,RH, alphanit,  &
                      fcoj,mass,nvap,mpt,dpaw,n,hyst,              &    
                      henry_eff_hno3,henry_eff_hcl,henry_eff_nh3,  &
                      svmc_min_hno3,svmc_min_hcl,svmc_min_nh3,     &
                      cioncharge,kh_nh3,flag_dissolution,          &
                      Keq_nh4no3_0,Keq_nh4cl_0)

         endif
!MESA end



! After cloud end
         if (wascloud) then
!AQ partitioning start 
! Drive gases out after evaporation of cloud
           if (IAQP.eq.1) then
             c(ind_FeOH2p_a(:)) = 0._dp
             xnom7nno    = 1                 !allow partitioning nitrosamines
             xaeq(:)     = 0._dp             !amine acid-base off
             xaer(:)     = 0._dp             !switch off aq. chemistry
             k_exf(:,:)  = 0._dp             !stop uptake (this was commented out)
             k_exb(:,:)  = k_exb(:,:)*1.e9   !drive gas out of droplets
             ! Henry partitioning
             CALL transfer_coeff(radius,temp,press,xaer,lwc,Henry_T0,  &
                 Henry_Tdep,alpha_T0,alpha_Tdep,molar_mass,k_exf,k_exb)
           endif
!AQ partitioning end
           caqold(:,:)=0._dp
          ! wascloud=.false.   ! set to false in water_content
         endif

       endif


       daytime = daynr + (model_time/OneDay)


! Gas phase chemistry integration

       if (ICHEM.EQ.1) then

         !  do jl=1,nspec
         !      write(6,*) 'c',jl,chem_name(jl),c(jl),k_exf(3,jl),k_exb(3,jl)
         !  end do

          call chemistry_solver(c,DTIME,firstloop,model_time,daytime,lat_deg, &
                 temp,press,cair,jno2m,RH,F_HONO,fco,lwc,cvfac,xaer,k_exf,    &
                 k_exb,k_exf_N2O5,k_exf_ClNO3,k_exf_BrNO3,   fcoj,jx) 


       endif



!!! Integration of aerosol processes !!!

! number concentration [1/m^3] threshold for numerical stability
       do M=NU,CS
        do I=1,IMAX
          IF (N(M,I) .LT. nucomin) THEN
            do K=1,AMAX-1
              MASS(M,I,K)=0._dp
            end do 
            N(M,I)=0._dp 
          ENDIF
        end do
       end do

! mass concentration [ng/m^3] threshold for numerical stability
       do M=NU,CS
        do I=1,IMAX
          do K=1,AMAX-1
            if ( MASS(M,I,K).lt.massmin ) MASS(M,I,K)=0._dp
          end do 
        end do
       end do

! Condensable vapours 
       NVAP(A_SUL)=c(ind_H2SO4)*1.e6_dp   ! [molec/m3]
       NVAP(A_MSA)=c(ind_CH3SO3H)*1.e6_dp
       NVAP(A_OR1)=c(ind_BSOV)*1.e6_dp 
       NVAP(A_OR2)=c(ind_BLOV)*1.e6_dp 
       NVAP(A_OR3)=c(ind_BELV)*1.e6_dp
       NVAP(A_OR4)=c(ind_ASOV)*1.e6_dp
       NVAP(A_OR5)=c(ind_ALOV)*1.e6_dp
       NVAP(A_OR6)=c(ind_AELV)*1.e6_dp
       NVAP(A_OR7)=c(ind_PIOV)*1.e6_dp
       NVAP(A_OR8)=c(ind_PSOV)*1.e6_dp
       NVAP(A_OR9)=c(ind_PELV)*1.e6_dp
       NVAP(A_NH4)=c(ind_NH3 )*1.e6_dp
       NVAP(A_NIT)=c(ind_HNO3)*1.e6_dp
       NVAP(A_CHL)=c(ind_HCl )*1.e6_dp
       ! Run for a single amine
       IF (IAM .EQ. 1)  NVAP(A_AMI)=c(ind_MEA)*1.e6_dp
       IF (IAM .EQ. 2)  NVAP(A_AMI)=c(ind_MMA)*1.e6_dp
       IF (IAM .EQ. 3)  NVAP(A_AMI)=c(ind_DMA)*1.e6_dp
       IF (IAM .EQ. 4)  NVAP(A_AMI)=c(ind_TMA)*1.e6_dp                     
       IF (IAM .EQ. 5)  NVAP(A_AMI)=c(ind_CH2NCH3)*1.e6_dp
       IF (IAM .EQ. 6)  NVAP(A_AMI)=c(ind_AMP)*1.e6_dp 


! Cloud droplet activation happens if relative humidity is above 99%
! and the incloud flag is set to 1 by the user and if the computed
! supersaturation is greater than equlibrium supersaturation.

      !  print*,'hour time(s) ',hour_time
      !  print *,'S ',supersat
      !  print *,'LWC(AI-CS) ',lwc(1),lwc(2),lwc(3)
      !   print *,'NH4+',c(ind_Hp_a02),c(ind_NH3_a02),c(ind_NH4p_a02)

       if ((incloud.EQ.1).and.(RH.ge.rhactiv)) then

           write(6,*) 'MAFOR stops because RH exceeds RHcrit !'
           stop


        ! ccnactivation returns change of droplet diameter
        ! and change of supersaturation
        !
        ! call ccnactivation(c,DTIME,IMAX,temp,press,RH,zmbl,rain,    &
        !                   DPA,DPAW,VPT,lwc,ROOP,MASS,surf_org,     &
        !                   surfin,M_oc(1),M_oc(2),hour_time,vupdra,supersat,  &
        !                   caqold,NTOT,N,IAG,COAGS,ddpwdt,dmwdt,dtedt,    &
        !                   dsupdt,Nfa,Mfa,DPcrit )
        !
        ! new supersaturation
        ! supersat = supersat + dsupdt*DTIME
        ! supersaturation should not go below RH=94%
        ! supersat = max(supersat,-0.06_dp)
        ! cloud maximum supersaturation: 0.11 (11%)
        ! supersat = min(supersat,0.11_dp)
        ! new temperature
        ! temp = temp + dtedt*DTIME

       else


! Compute aerosol processes
         call aerosol_solver(DTIME,IMAX,press,temp,ROOP,ROOPW,DPA,DPAW,VPT,N, &
                  MASS,IAG,NVAP,                                              &
                  zmbl,rain,hsta_st,u10,cair,RH,daytime,lat_deg,              &
                  hour_timemin,jx(ip_NO2),alphanit,fcoj,                      &
                  CAMI,KP_NIT,fnuc,INUCMEC,                                   &
                  DEN(OC),surfin,surf_org,M_oc,nmo,foc,hvap,csat0,            &
                  gamma_oc1_m,gamma_oc2_m,gamma_oc3_m,gamma_oc4_m,            &
                  gamma_oc5_m,gamma_oc6_m,gamma_oc7_m,gamma_oc8_m,            &
                  gamma_oc9_m,                                                &
                  MTOT,MORGCTOT,MECBCTOT,MDUSTTOT,                            &
                  MORG1TOT,MORG2TOT,MORG3TOT,MORG4TOT,MORG5TOT,               &
                  MORG6TOT,MORG7TOT,MORG8TOT,MORG9TOT,                        &
                  henry_eff_hno3,henry_eff_hcl,svmc_min_hno3,svmc_min_hcl,    & 
                  svmc_min_nh3, cioncharge,kh_nh3,flag_dissolution,           &
                  Keq_nh4no3_0,Keq_nh4cl_0,                                   &
                  COAGS,VVAP,JNUC,natot,GRTOT,CSSULT,CSORGT,NVAPO,csate  )

!debug
          !write(6,*) 'nvap biog org',NVAP(A_OR1),NVAP(A_OR2),NVAP(A_OR3)
          !write(6,*) 'nvap arom org',NVAP(A_OR4),NVAP(A_OR5),NVAP(A_OR6)
          !write(6,*) 'nvap prim org',NVAP(A_OR7),NVAP(A_OR8),NVAP(A_OR9)

       endif


! dilution with background particles [#/cm3]        
       IF (IDIL>0) THEN
         do M=NU,CS
          do I=1,IMAX
            do K=1,AMAX
              DIFFDILM(M,I,K)=MASS(M,I,K)-BGMASS(M,I,K)
            end do
            DIFFDILN(M,I)=N(M,I)-BGN(M,I)
          end do
         end do
       ENDIF


! -----------------------------------------------------  
! UPDATE MASS AND NUMBER CONCENTRATION BY ALL PROCESSES
! -----------------------------------------------------      
    
! DILUTION (with background particles)
       IF (IDIL>0) THEN
          do M=NU,CS
            do I=1,IMAX       
              do K=1,AMAX                  
                  MASS(M,I,K)=MASS(M,I,K)      -        &
                      dilrate*DIFFDILM(M,I,K)*DTIME
              end do
              N(M,I)=N(M,I)-dilrate*DIFFDILN(M,I)*DTIME
            end do
          end do                 
       ENDIF

! SEA-SALT EMISSION   
       IF (IPEMI.EQ.2) THEN
         do M=NU,CS
          do I=1,IMAX
            MASS(M,I,A_CHL)=MASS(M,I,A_CHL)        +     &
                 0.54*(saltemis(M,I)/zmbl) *CONVM  *     &
                 DEN(SA)*VPT(M,I)*DTIME
            MASS(M,I,A_SAL)=MASS(M,I,A_SAL)        +     &
                 0.46*(saltemis(M,I)/zmbl) *CONVM  *     &
                 DEN(SA)*VPT(M,I)*DTIME                                                      
            N(M,I)=N(M,I)+(saltemis(M,I)/zmbl)*DTIME
          end do
         end do
       ENDIF

! EMISSION (continuous EN and EMASS)   
       IF (IPEMI.EQ.1) THEN
         do M=NU,CS
          do I=1,IMAX
            do K=1,AMAX
              MASS(M,I,K)=MASS(M,I,K)              +     &
                    emisratp*(EMASS(M,I,K)/zmbl)*DTIME
            end do
            N(M,I)=N(M,I)+emisratp*(EN(M,I)/zmbl)*DTIME
          end do
         end do
       ENDIF


! LAST CHECK FOR NEGATIVE NUMBER CONCENTRATIONS
       do M=NU,CS
        do I=1,IMAX
          IF (N(M,I) .LT. 0._dp) THEN
            write(6,*) 'WARNING: negative N before Nucl.',M,I,N(M,I)
            write(6,*) 'program continues...'
            if (IDEB == 1) then
              write(12,fmt='(a,i2,i2,f15.2)') '    WARNING: negative N before Nucleation ',M,I,N(M,I)
            endif
            !stop
          ENDIF
          N(M,I) = max(N(M,I),0.0_dp)
        end do
       end do

! PN PARAMETERIZATION
       PNC1 = PNC1 - dilrate* (PNC1- BGNTOT(AI)) *DTIME
       PNC2 = PNC2 - dilrate* (PNC2- 0.6*BGNTOT(AS)) *DTIME
       PNC3 = PNC3 - dilrate* (PNC3- 0.2*BGNTOT(AS)) *DTIME
       PNC4 = PNC4 - dilrate* (PNC4- 0.2*BGNTOT(AS)) *DTIME
       PNC5 = PNC5 - dilrate* (PNC5- 0.5*BGNTOT(CS)) *DTIME
       PNC6 = PNC6 - dilrate* (PNC6- 0.5*BGNTOT(CS)) *DTIME
       PNC1 = PNC1 - PNC1* PNC1*kcoa1 *DTIME  - (vd1/zmblm)* PNC1 *DTIME
       PNC2 = PNC2 - PNC2* PNC2*kcoa2 *DTIME  - (vd2/zmblm)* PNC2 *DTIME
       PNC3 = PNC3 - PNC3* PNC3*kcoa3 *DTIME  - (vd3/zmblm)* PNC3 *DTIME
       PNC4 = PNC4 - PNC4* PNC4*kcoa4 *DTIME  - (vd4/zmblm)* PNC4 *DTIME
       PNC5 = PNC5 - PNC5* PNC5*kcoa5 *DTIME  - (vd5/zmblm)* PNC5 *DTIME
       PNC6 = PNC6 - PNC6* PNC6*kcoa6 *DTIME  - (vd6/zmblm)* PNC6 *DTIME

! END UPDATE NUMBER AND MASS CONCENTRATION

! Amine loss rate to particles
       IF (ICONA.EQ.1) THEN
         DNVAPC  = (NVAPO(A_AMI)-NVAP(A_AMI))/(NVAPO(A_AMI)*DTIME)
       ENDIF
! New gas phase concentration after condensation/nucleation
       c(ind_H2SO4)  = NVAP(A_SUL)*1.e-6_dp   ! [molec/cm3]
       c(ind_CH3SO3H)= NVAP(A_MSA)*1.e-6_dp 
       c(ind_BSOV)   = NVAP(A_OR1)*1.e-6_dp
       c(ind_BLOV)   = NVAP(A_OR2)*1.e-6_dp
       c(ind_BELV)   = NVAP(A_OR3)*1.e-6_dp
       c(ind_ASOV)   = NVAP(A_OR4)*1.e-6_dp
       c(ind_ALOV)   = NVAP(A_OR5)*1.e-6_dp
       c(ind_AELV)   = NVAP(A_OR6)*1.e-6_dp
       c(ind_PIOV)   = NVAP(A_OR7)*1.e-6_dp
       c(ind_PSOV)   = NVAP(A_OR8)*1.e-6_dp
       c(ind_PELV)   = NVAP(A_OR9)*1.e-6_dp
       c(ind_HNO3)   = NVAP(A_NIT)*1.e-6_dp
       c(ind_NH3 )   = NVAP(A_NH4)*1.e-6_dp
       c(ind_HCl )   = NVAP(A_CHL)*1.e-6_dp
       IF (IAM .EQ. 1)  c(ind_MEA)     = NVAP(A_AMI)*1.e-6_dp
       IF (IAM .EQ. 2)  c(ind_MMA)     = NVAP(A_AMI)*1.e-6_dp
       IF (IAM .EQ. 3)  c(ind_DMA)     = NVAP(A_AMI)*1.e-6_dp
       IF (IAM .EQ. 4)  c(ind_TMA)     = NVAP(A_AMI)*1.e-6_dp       
       IF (IAM .EQ. 5)  c(ind_CH2NCH3) = NVAP(A_AMI)*1.e-6_dp 
       IF (IAM .EQ. 6)  c(ind_AMP)     = NVAP(A_AMI)*1.e-6_dp
! for SOA output, SOA_gas     [ug/m^3]
       csoagas(1)=NVAP(A_OR1)*1.e-6*molec2ug(M_oc(1))
       csoagas(2)=NVAP(A_OR2)*1.e-6*molec2ug(M_oc(2))
       csoagas(3)=NVAP(A_OR3)*1.e-6*molec2ug(M_oc(3))
       csoagas(4)=NVAP(A_OR4)*1.e-6*molec2ug(M_oc(4))
       csoagas(5)=NVAP(A_OR5)*1.e-6*molec2ug(M_oc(5))
       csoagas(6)=NVAP(A_OR6)*1.e-6*molec2ug(M_oc(6))
       csoagas(7)=NVAP(A_OR7)*1.e-6*molec2ug(M_oc(7))
       csoagas(8)=NVAP(A_OR8)*1.e-6*molec2ug(M_oc(8))
       csoagas(9)=NVAP(A_OR9)*1.e-6*molec2ug(M_oc(9))


! Calculate new NTOT
       do M=NU,CS
         NTOT(M)=0._dp
         do I=1,IMAX
           NTOT(M)=NTOT(M)+N(M,I)
         enddo
       enddo


! Calculate new water content
      ! Water activity from sulphate - ammonium - nitrate. MH2O in ng/m^3
      ! Munkelwitz and Tang 1994 / Seinfeld and Pandis 1997
      ! Sea salt hygroscopicity same  
      ! If in-cloud, MH2O and N of CS mode is calculated from LWC(CS)

       if ((incloud.EQ.1).and.(RH.ge.rhactiv)) then
      ! cloud activation m(H2O) changed dynamically
      ! do we need a lower limit for mH2O?
         do M=NU,CS
          do I=1,IMAX       
            dmwdt(M,I) = dmwdt(M,I)*DTIME
          enddo
         enddo
       endif

       CALL water_content(incloud,firstloop,IMAX,RH,lwcm,MASS,N,      &
                dmwdt,DPA,GMD,SIG,DLINDP,VPT,DPAW,MH2OTOT,NTOT(CS),   &
                lwc,c(ind_H2O_a(1)),xaer,wascloud  )


! Calculate aerosol evolution parameters

! TOTAL MASS
! Total dry mass and wet mass in kg/m^3
! Total mode masses in ng/m^3
       CALL getTotalMass(IMAX,MASS,MPT,MPTW,                      &
                      MTOT,MTOTW,MSULFTOT,MMSAPTOT,               &
                      MNITRTOT,MAMMOTOT,MSALTTOT,MECBCTOT,        &
                      MDUSTTOT,MXXXXTOT,MORGCTOT,MORG1TOT,        &
                      MORG2TOT,MORG3TOT,MORG4TOT,MORG5TOT,        &
                      MORG6TOT,MORG7TOT,MORG8TOT,MORG9TOT, casoa   )


! AVERAGE DENSITY
! Density in kg/m^3

       call getDensity(IMAX,temp,DEN(EC),DPA,MASS,MPT,MPTW,ROOP,ROOPW)


! DIAMETER OF WET PARTICLE
! Wet diameter/volume is used to calculate aerosol processes
       IF (ICONW .ge. 1) THEN
         if ((incloud.EQ.1).and.(RH.ge.rhactiv)) then
          ! cloud activation DPAW changed dynamically
          do M=NU,CS
            do I=1,IMAX
              DPAW(M,I) = DPAW(M,I) + ddpwdt(M,I)*DTIME
              ! lower limit DPA*1.1
              DPAW(M,I) = max(DPAW(M,I),1.1*DPA(M,I))
            enddo
          enddo
         else
           CALL wetdiameter(IMAX,incloud,RH,MTOT,MTOTW,DPA,DPAW)
         endif
       ELSE
         do M=NU,CS
           do I=1,IMAX
             DPAW(M,I)=DPA(M,I)             
           end do
         end do          
       ENDIF

! UPDATE SIZE DISTRIBUTION
       do M=NU,CS
         do I=1,IMAX
          LOGDIS(M,I)=N(M,I)/DLOGDP(M,I)
          MLOGDIS(M,I)=MPT(M,I)/DLOGDP(M,I)
         end do
       end do
       do M=NU,CS
        VTOT(M)=0._dp
        do I=1,IMAX
          VTOT(M) = VTOT(M)+VPT(M,I)    
         end do
       end do

       N3=0._dp
       NEU1=NTOT(NU)
       ! NEU2: 25-100m particles
       NEU2=0._dp
       NEU3=0._dp
       NEU4=0._dp
       NEU5=0._dp
       NEU6=0._dp
       NEU7=0._dp    
       do I=1,IMAX
         if (DPA(NU,I).GT.3.E-9_dp) N3=N3+N(NU,I)
         if (DPA(AI,I).LT.25.E-9_dp) NEU1=NEU1+N(AI,I)
         if (DPA(AI,I).GE.25.E-9_dp) NEU2=NEU2+N(AI,I)
         if (DPA(AS,I).LT.1.E-7_dp)  NEU2=NEU2+N(AS,I)
         if ((DPA(AI,I).GE.24.E-9_dp).and.(DPA(AI,I).LT.50.E-9_dp)) NEU3=NEU3+N(AI,I)
         if ((DPA(AS,I).GE.45.E-9_dp).and.(DPA(AS,I).LT.50.E-9_dp)) NEU3=NEU3+N(AS,I)
         if ((DPA(AI,I).GE.50.E-9_dp).and.(DPA(AI,I).LT.75.E-9_dp)) NEU4=NEU4+N(AI,I)
         if ((DPA(AS,I).GE.50.E-9_dp).and.(DPA(AS,I).LT.75.E-9_dp)) NEU4=NEU4+N(AS,I)
         if (DPA(AI,I).GE.75.E-9_dp) NEU5=NEU5+N(AI,I)
         if (DPA(AS,I).GE.75.E-9_dp) NEU5=NEU5+N(AS,I)
         if (DPA(CS,I).LT.2.E-7_dp)  NEU6=NEU6+N(CS,I)
         if (DPA(CS,I).GE.2.E-7_dp)  NEU7=NEU7+N(CS,I)
      end do


! OUTPUT EVERY 0.1 SEC!
       if ((IDIL==3).or.(IDIL==4)) then
        if (hour_timeo.GT.0.0999) then
          write(16,fc) model_time,(c(jq),jq=1,nspec),(lwc(zkc),zkc=1,APN),jx(ip_NO2)
          write(14,fb) model_time,(LOGDIS(NU,I),I=1,IMAX),(LOGDIS(AI,I), I=1,IMAX), &
                 (LOGDIS(AS,I), I=1,IMAX),(LOGDIS(CS,I), I=1,IMAX)
          write(18,fb) model_time, (DPAW(NU,I), I=1,IMAX),(DPAW(AI,I), I=1,IMAX),   &
                 (DPAW(AS,I), I=1,IMAX),(DPAW(CS,I), I=1,IMAX)
          write(13,'(29e17.6)') model_time,N3,NTOT(NU),NTOT(AI),NTOT(AS),NTOT(CS),    &
                 NEU1,NEU2,DNVAPC,COAGS,CSORGT,GRTOT,JNUC,                            &
                 NEU3,NEU4,NEU5,NEU6,NEU7,PNC1,PNC2,PNC3,PNC4,PNC5,PNC6,     &
                 supersat,Nfa,Mfa,DPcrit,temp
          write(15,fb) model_time,(MLOGDIS(NU,I),I=1,IMAX),(MLOGDIS(AI,I),I=1,IMAX),&
                 (MLOGDIS(AS,I), I=1,IMAX),(MLOGDIS(CS,I), I=1,IMAX)
          do K=1,AMAX        
             write(15,fb) model_time,((MASS(NU,I,K)/DLOGDP(NU,I)),I=1,IMAX),   &
                 ((MASS(AI,I,K)/DLOGDP(AI,I)),I=1,IMAX),                       &
                 ((MASS(AS,I,K)/DLOGDP(AS,I)),I=1,IMAX),                       &
                 ((MASS(CS,I,K)/DLOGDP(CS,I)),I=1,IMAX) 
          end do
          write(17,'(41e17.6)') model_time,                             &
               MSULFTOT(NU),MSULFTOT(AI),MSULFTOT(AS),MSULFTOT(CS),   &
               MMSAPTOT(NU),MMSAPTOT(AI),MMSAPTOT(AS),MMSAPTOT(CS),   &
               MXXXXTOT(NU),MXXXXTOT(AI),MXXXXTOT(AS),MXXXXTOT(CS),   &
               MORGCTOT(NU),MORGCTOT(AI),MORGCTOT(AS),MORGCTOT(CS),   &
               MAMMOTOT(NU),MAMMOTOT(AI),MAMMOTOT(AS),MAMMOTOT(CS),   &
               MNITRTOT(NU),MNITRTOT(AI),MNITRTOT(AS),MNITRTOT(CS),   &
               MECBCTOT(NU),MECBCTOT(AI),MECBCTOT(AS),MECBCTOT(CS),   &               
               MDUSTTOT(NU),MDUSTTOT(AI),MDUSTTOT(AS),MDUSTTOT(CS),   &
               MSALTTOT(NU),MSALTTOT(AI),MSALTTOT(AS),MSALTTOT(CS),   &
               MH2OTOT(NU), MH2OTOT(AI), MH2OTOT(AS), MH2OTOT(CS)
          write(19,'(12e12.4)') model_time,dilut_time,temp_old,       & 
               zmbl_old,wplm_old,areapl,dilrate,                      &
               NVAP(A_SUL)*1.e-6*molec2ug(M_H2SO4),                   &
               csoagas(7),csoagas(8),csoagas(9),                      &
               (ENTOT(NU)+ENTOT(AI)+ENTOT(AS)+ENTOT(CS))*emisratp
          hour_timeo=0.
        endif
       endif

! OUTPUT EVERY 1 SEC!
       if ((IDIL==1).or.(IDIL==2).or.(IDIL>4)) then
        if ((hour_timeo.GT.0.999).and.(hsta_st>30.0) ) then
          write(16,fc) model_time,(c(jq),jq=1,nspec),(lwc(zkc),zkc=1,APN),jx(ip_NO2)
          write(14,fb) model_time,(LOGDIS(NU,I),I=1,IMAX),(LOGDIS(AI,I), I=1,IMAX), &
                  (LOGDIS(AS,I), I=1,IMAX),(LOGDIS(CS,I), I=1,IMAX)
          write(18,fb) model_time, (DPAW(NU,I), I=1,IMAX),(DPAW(AI,I), I=1,IMAX),   &
                  (DPAW(AS,I), I=1,IMAX),(DPAW(CS,I), I=1,IMAX)
          write(13,'(29e17.6)') model_time,N3,NTOT(NU),NTOT(AI),NTOT(AS),NTOT(CS),  &
                  NEU1,NEU2,DNVAPC,COAGS,CSORGT,GRTOT,JNUC,                         &
                  NEU3,NEU4,NEU5,NEU6,NEU7,PNC1,PNC2,PNC3,PNC4,PNC5,PNC6,           &
                  supersat,Nfa,Mfa,DPcrit,temp
          write(15,fb) model_time,(MLOGDIS(NU,I),I=1,IMAX),(MLOGDIS(AI,I),I=1,IMAX),&
                  (MLOGDIS(AS,I), I=1,IMAX),(MLOGDIS(CS,I), I=1,IMAX)
          do K=1,AMAX        
              write(15,fb) model_time,((MASS(NU,I,K)/DLOGDP(NU,I)),I=1,IMAX),   &
                  ((MASS(AI,I,K)/DLOGDP(AI,I)),I=1,IMAX),                       &
                  ((MASS(AS,I,K)/DLOGDP(AS,I)),I=1,IMAX),                       &
                  ((MASS(CS,I,K)/DLOGDP(CS,I)),I=1,IMAX) 
          end do
          write(19,'(12e12.4)') model_time,dilut_time,temp_old,       & 
               zmbl_old,wplm_old,areapl,dilrate,                      &
               NVAP(A_SUL)*1.e-6*molec2ug(M_H2SO4),                   &
               csoagas(7),csoagas(8),csoagas(9),                      &
               (ENTOT(NU)+ENTOT(AI)+ENTOT(AS)+ENTOT(CS))*emisratp
          hour_timeo=0.
        endif
       endif


! OUTPUT EVERY 10 SEC!
       if (hour_timeo .GT. 9.9999 ) then
         hour_timeo=0.

! write total number conc. and volume-mean-diameter current time step
         write(13,'(29e17.6)') model_time,N3,NTOT(NU),NTOT(AI),NTOT(AS),NTOT(CS),  &
               NEU1,NEU2,DNVAPC,COAGS,CSORGT,GRTOT,JNUC,                           &                 
               NEU3,NEU4,NEU5,NEU6,NEU7,PNC1,PNC2,PNC3,PNC4,PNC5,PNC6,             &
               supersat,Nfa,Mfa,DPcrit,temp
! write size distribution current time step
         write(14,fb) model_time,(LOGDIS(NU,I),I=1,IMAX),(LOGDIS(AI,I), I=1,IMAX), &
                 (LOGDIS(AS,I), I=1,IMAX),(LOGDIS(CS,I), I=1,IMAX)
! write mass conc. size distribution current time step (kg/m^3)
         write(15,fb) model_time,(MLOGDIS(NU,I),I=1,IMAX),(MLOGDIS(AI,I), I=1,IMAX),&
                 (MLOGDIS(AS,I), I=1,IMAX),(MLOGDIS(CS,I), I=1,IMAX)
! write wet diameter of particles                 
         write(18,fb) model_time, (DPAW(NU,I), I=1,IMAX),(DPAW(AI,I), I=1,IMAX),    &
                 (DPAW(AS,I), I=1,IMAX),(DPAW(CS,I), I=1,IMAX)                 
! write concentrations current time step
         write(16,fc) model_time,(c(jq),jq=1,nspec),(lwc(zkc),zkc=1,APN),jx(ip_NO2)
! write aerosol concentrations current time step
         write(17,'(41e12.4)') model_time,                              &
               MSULFTOT(NU),MSULFTOT(AI),MSULFTOT(AS),MSULFTOT(CS),   &
               MMSAPTOT(NU),MMSAPTOT(AI),MMSAPTOT(AS),MMSAPTOT(CS),   &
               MXXXXTOT(NU),MXXXXTOT(AI),MXXXXTOT(AS),MXXXXTOT(CS),   &
               MORGCTOT(NU),MORGCTOT(AI),MORGCTOT(AS),MORGCTOT(CS),   &
               MAMMOTOT(NU),MAMMOTOT(AI),MAMMOTOT(AS),MAMMOTOT(CS),   &
               MNITRTOT(NU),MNITRTOT(AI),MNITRTOT(AS),MNITRTOT(CS),   &
               MECBCTOT(NU),MECBCTOT(AI),MECBCTOT(AS),MECBCTOT(CS),   &               
               MDUSTTOT(NU),MDUSTTOT(AI),MDUSTTOT(AS),MDUSTTOT(CS),   &
               MSALTTOT(NU),MSALTTOT(AI),MSALTTOT(AS),MSALTTOT(CS),   &
               MH2OTOT(NU), MH2OTOT(AI), MH2OTOT(AS), MH2OTOT(CS)
! write plume dispersion parameters
          write(19,'(12e12.4)') model_time,dilut_time,temp_old,       & 
               zmbl_old,wplm_old,areapl,dilrate,                      &
               NVAP(A_SUL)*1.e-6*molec2ug(M_H2SO4),                   &
               csoagas(7),csoagas(8),csoagas(9),                      &
               (ENTOT(NU)+ENTOT(AI)+ENTOT(AS)+ENTOT(CS))*emisratp
               
! write soa component distribution
         write(11,'(28e17.6)') model_time , (csate(S), S=1,NSOA),     &
                       (csoagas(S), S=1,NSOA), (casoa(S), S=1,NSOA)
       endif



! Next time step
       firstloop=.false.
       model_time=model_time+DTIME
       hour_time=hour_time+DTIME
       hour_timec=hour_timec+DTIME
       hour_timet=hour_timet+DTIME
       hour_timeo=hour_timeo+DTIME
       dilut_time=dilut_time+DTIME

 !   print *,'diltime(t)',dilut_time,emisratp,dilrate

       if (IDEB == 1) then
         if (hour_timemin .EQ. 120.) then
          write(12,*) 'daytime   ',daytime
          write(12,*) 'model time',model_time
          write(12,*) 'JNUC [m3/s]         NO2 [molec/cm3]'
          write(12,'(4ES12.4)') JNUC,c(ind_NO2)   

          hour_timemin=0.
         endif
         hour_timemin=hour_timemin+DTIME
       endif  


    END DO

    if (IDEB == 1) then
      call cpu_time(cpu_time_end)
      write(12,*) 'CPU time',cpu_time_end-cpu_time_start
    endif  

! END INTEGRATION LOOP
!*************************************************************

! Close input files
    close(24)
! Close output files
    close(11)
    close(12)
    close(13)
    close(14)
    close(16)
    close(17)     
    close(18)
    close(19)


! Deallocate all arrays
    if (ICHEM.eq.1) CALL photo_dealloc()  

    deallocate(ROOP,ROOPW)
    deallocate(CCOND)
    deallocate(NVAP)
    deallocate(NVAPO)
    deallocate(DLOGDP,DLINDP,N,LINDIS,LOGDIS,MLOGDIS)
    deallocate(VPT,MPT,MPTW)
    deallocate(MASS) 
    deallocate(DPA,DPAW)
    deallocate(BGLOGDIS,DIFFDILM,DIFFDILN)
    deallocate(BGMASS,BGN)
    deallocate(EMASS,EN)     
    deallocate(IAG)
!MESA start
    deallocate(HYST)
    deallocate(henry_eff_hno3,henry_eff_hcl,henry_eff_nh3)
    deallocate(svmc_min_hno3,svmc_min_hcl,svmc_min_nh3)
    deallocate(cioncharge,kh_nh3)
    deallocate(flag_dissolution)
!MESA end
!CCN start
    deallocate(ddpwdt)
    deallocate(dmwdt)
!CCN end

   end program mafor
