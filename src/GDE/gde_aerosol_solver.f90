! <gde_aerosol_solver.f90 - A component of the Multicomponent
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
module gde_aerosol_solver

  use gde_sensitiv,   only      : IDEB
  use gde_sensitiv,   only      : IDEPO,ICHAM,IWETD,ICOAG
  use gde_sensitiv,   only      : INUC
  use gde_sensitiv,   only      : ICOND,IKELV,ICONS,ICONO,ICONA
  use gde_sensitiv,   only      : ICONX,INANO
  use gde_sensitiv,   only      : ICONW

  use gde_constants,  only      : M_H2SO4,M_nh3,M_nit,M_msa,M_hio3
  use gde_constants,  only      : M_ca,M_hcl
  use gde_constants,  only      : MB,MAH,MAN,MIO
  use gde_constants,  only      : pi

! for kelvin effect of n-alkanes
  use gde_init_gas, only        : gamma_oc7,gamma_oc8
! for chamber studies
  use gde_init_gas, only        : V_CHAM,S_SED,S_DIF
  use gde_init_gas, only        : DILPAR

  use gde_input_data, only      : MMAX,AMAX,QMAX
  use gde_input_data, only      : NSOA
  use gde_input_data, only      : A_SUL,A_NH4,A_AMI,A_NIT,A_XXX
  use gde_input_data, only      : A_OR1,A_OR2,A_OR3,A_OR4,A_OR5
  use gde_input_data, only      : A_OR6,A_OR7,A_OR8,A_OR9
  use gde_input_data, only      : A_MSA,A_CHL,A_IO3
  use gde_input_data, only      : A_WAT
  use gde_input_data, only      : NU,AI,AS,CS
  use gde_input_data, only      : massmin
  use gde_input_data, only      : nucomin
  use gde_input_data, only      : CTNH4,CTNIT,CTSO4,CTCHL
  use gde_input_data, only      : KPEQ
  use gde_input_data, only      : DCSU,DCORG
  use gde_input_data, only      : CONVM
  use gde_input_data, only      : DENV
  use gde_input_data, only      : DEB_CTOTS1,DEB_CTOTS2
  use gde_input_data, only      : DEB_CTOTA1,DEB_CTOTA2
  use gde_input_data, only      : DEB_NTOT,DEB_MTOT
  use gde_input_data, only      : DEB_CTOTO11,DEB_CTOTO12
  use gde_input_data, only      : DEB_CTOTO21,DEB_CTOTO22

  use gde_toolbox,    only      : molec2ug
  use gde_toolbox,    only      : newton
  use gde_toolbox,    only      : machineps

  use gde_deposition,     only  : depositpar
  use gde_deposition,     only  : depocanopy
  use gde_deposition,     only  : depozhang01  
  use gde_deposition,     only  : depofrough
  use gde_deposition,     only  : depositwall
  use gde_deposition,     only  : wetscavbulk
  use gde_deposition,     only  : wetscavsize
  use gde_coagulation,    only  : coagulation_coeff
  use gde_coagulation,    only  : coagulation
  use gde_coagulation,    only  : coagulation_target
  use gde_nucleation,     only  : nucleation
  use gde_nucleation,     only  : nucleationratio
  use gde_transfer,       only  : kelvin_nit
  use gde_transfer,       only  : kelvin_sulf
  use gde_transfer,       only  : kelvin_msap
  use gde_transfer,       only  : kelvin_org
  use gde_transfer,       only  : kelvin_alkane
  use gde_transfer,       only  : nano_koehler
  use gde_condensation,   only  : condensloss
  use gde_condensation,   only  : condensation
  use gde_condensation,   only  : apc_update_gasc
  use gde_condensation,   only  : apc_update_pmass


  private

  public :: aerosol_solver
  public :: interface_mosaic

! KPP DP - Double precision kind
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)


contains


  subroutine aerosol_solver(DTIME,IMAX,press,temp,DENSPAR,DENSPARW,DPA,DPAW,VPT,N,     &
                           MASS,IAG,NVAP,                                              &
                           mbl,rain,hsta_st,u10,cair,RH,daytime,lat_deg,               &
                           incloud,owf,hour_timemin,                                   &
                           jrno2,alphanit,fcoj,                                        &
                           CAMI,KP_NIT,fnuc,INUCMEC,                                   &
                           DENOC,surfin,surf_org,                                      &
                           M_oc,nmo,foc,hvap,csat0,                                    &
                           gamma_oc1_m,gamma_oc2_m,gamma_oc3_m,gamma_oc4_m,            &
                           gamma_oc5_m,gamma_oc6_m,gamma_oc7_m,gamma_oc8_m,            &
                           gamma_oc9_m,                                                &
                           MTOT,MORGCTOT,MECBCTOT,MDUSTTOT,                            &
                           MORG1TOT,MORG2TOT,MORG3TOT,MORG4TOT,MORG5TOT,               &
                           MORG6TOT,MORG7TOT,MORG8TOT,MORG9TOT,                        &
                           henry_eff_hno3,henry_eff_hcl,svmc_min_hno3,svmc_min_hcl,    &
                           svmc_min_nh3, cioncharge, kh_nh3,                           &
                           flag_dissolution, Keq_nh4no3_0,Keq_nh4cl_0,                 &
                           coags,vvap,jnuc,natot,                                      &
                           GRTOT,CSSULT,CSORGT,nvapo,csate     )
    !----------------------------------------------------------------------
    !     
    ! Main routine of MAFOR aerosol dynamics solver
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      call the aerosol dynamic processes
    !      compute new gas phase concentration after condensation
    !      compute new particle mass concentrations after condensation
    !      update particle mass and number nucl/coag/depo
    !
    !      interface
    !      ---------
    !      
    !
    !      method
    !      ------
    !      deposition
    !      coagulation
    !      Analytical Predictor for Condensation
    !      mass balance checks
    !      Euler differences to solve N and MASS
    !      growth rate
    !
    !      
    !      external
    !      --------
    !      none
    !      
    !      reference
    !      ---------
    !      none
    !
    !------------------------------------------------------------------

    implicit none

    ! input
    integer, intent(in)                        :: IMAX
    real( dp), intent(in)                      :: DTIME
    real( dp), intent(in)                      :: temp
    real( dp), intent(in)                      :: press
    real( dp), intent(in)                      :: mbl
    real( dp), intent(in)                      :: hsta_st
    real( dp), intent(in)                      :: u10
    real( dp), intent(in)                      :: rain
    real( dp), intent(in)                      :: cair
    real( dp), intent(in)                      :: RH
    real( dp), intent(in)                      :: jrno2
    real( dp), intent(in)                      :: daytime
    real( dp), intent(in)                      :: lat_deg
    real( dp), intent(in)                      :: owf
    integer, intent(in)                        :: incloud
    real( dp), intent(in)                      :: hour_timemin
    real( dp), intent(in)                      :: CAMI,KP_NIT,fnuc
    real( dp), intent(in)                      :: alphanit,fcoj
    real( dp), intent(in)                      :: DENOC
    real( dp), intent(in)                      :: surf_org
    real( dp), dimension(NSOA), intent(in)     :: M_oc
    real( dp), dimension(NSOA), intent(in)     :: nmo
    real( dp), dimension(NSOA), intent(in)     :: foc
    real( dp), dimension(NSOA), intent(in)     :: hvap
    real( dp), dimension(NSOA), intent(in)     :: csat0
    real( dp), intent(in)                      :: gamma_oc1_m(NU:CS)
    real( dp), intent(in)                      :: gamma_oc2_m(NU:CS)
    real( dp), intent(in)                      :: gamma_oc3_m(NU:CS)
    real( dp), intent(in)                      :: gamma_oc4_m(NU:CS)
    real( dp), intent(in)                      :: gamma_oc5_m(NU:CS)
    real( dp), intent(in)                      :: gamma_oc6_m(NU:CS)
    real( dp), intent(in)                      :: gamma_oc7_m(NU:CS)
    real( dp), intent(in)                      :: gamma_oc8_m(NU:CS)
    real( dp), intent(in)                      :: gamma_oc9_m(NU:CS)
    real( dp), intent(in)                      :: vvap

    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPAW
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPA
    real( dp), DIMENSION(MMAX,IMAX), intent(in)       :: DENSPARW
    real( dp), DIMENSION(MMAX,IMAX), intent(in)       :: DENSPAR
    real( dp), DIMENSION(MMAX,IMAX), intent(in)       :: VPT
    real( dp), dimension(MMAX), intent(in)            :: MTOT
    real( dp), dimension(MMAX), intent(in)            :: MORG1TOT,MORG2TOT
    real( dp), dimension(MMAX), intent(in)            :: MORG3TOT,MORG4TOT
    real( dp), dimension(MMAX), intent(in)            :: MORG5TOT,MORG6TOT
    real( dp), dimension(MMAX), intent(in)            :: MORG7TOT,MORG8TOT
    real( dp), dimension(MMAX), intent(in)            :: MORG9TOT
    real( dp), dimension(MMAX), intent(in)            :: MORGCTOT
    real( dp), dimension(MMAX), intent(in)            :: MECBCTOT,MDUSTTOT

! mesa start
    real( dp),dimension(MMAX,IMAX),intent(in)         :: henry_eff_hno3
    real( dp),dimension(MMAX,IMAX),intent(in)         :: henry_eff_hcl
    real( dp),dimension(MMAX,IMAX),intent(in)         :: svmc_min_hno3
    real( dp),dimension(MMAX,IMAX),intent(in)         :: svmc_min_hcl
    real( dp),dimension(MMAX,IMAX),intent(in)         :: svmc_min_nh3
    real( dp),dimension(MMAX,IMAX),intent(in)         :: cioncharge
    real( dp),dimension(MMAX,IMAX),intent(in)         :: kh_nh3
    integer,  dimension(MMAX,IMAX),intent(in)         :: flag_dissolution
    real( dp),intent(in)                              :: Keq_nh4no3_0
    real( dp),intent(in)                              :: Keq_nh4cl_0
! mesa end

    integer, intent(in)                        :: INUCMEC
    integer, intent(in)                        :: surfin

    ! in/out
    real( dp), dimension(MMAX,MMAX,IMAX,IMAX), intent(in out) :: IAG
    real( dp), dimension(MMAX,IMAX,AMAX),intent(in out)       :: MASS
    real( dp), DIMENSION(MMAX,IMAX), intent(in out)           :: N
    real( dp), dimension(QMAX), intent(in out)                :: NVAP


    ! output
    real( dp), intent(out)                            :: coags
    real( dp), intent(out)                            :: jnuc
    real( dp), intent(out)                            :: natot
    real( dp), dimension(QMAX),intent(out)            :: nvapo
    real( dp), intent(out)                            :: CSSULT
    real( dp), intent(out)                            :: CSORGT
    real( dp), intent(out)                            :: GRTOT
    real( dp), dimension(NSOA),intent(out)            :: csate

    ! local
    integer                     :: M,I,Q,K,J

    ! for kelvin effect of n-alkanes
    real( dp)                   :: sum_gamma(NU:CS)
    real( dp)                   :: gamma_oc7_tmp
    real( dp)                   :: gamma_oc8_tmp
    ! wet particle scavenging [1/s]
    real( dp)                   :: wetsc
    ! in-cloud interstitial scavenging [1/s]
    real( dp)                   :: coagscav
    ! condensation
    real( dp)                   :: DNVAPN
    real( dp)                   :: NVAPOLDN
    real( dp)                   :: nnuc
    real( dp)                   :: GRSU,GRMS,GROC    ! growth rate [nm/hr]
    real( dp)                   :: sumlwc
    ! iteration
    real( dp)                   :: mtol
    real( dp)                   :: xout
    real( dp)                   :: ci
    real( dp)                   :: sf
    real( dp)                   :: mnh4
    logical                     :: deriv_zero = .false.
    ! condensation
    real( dp), dimension(QMAX)  :: nsv
    real( dp), dimension(QMAX)  :: FVAP,LVAP
    real( dp), dimension(QMAX)  :: CAT,CAT2
    real( dp), dimension(QMAX)  :: DNVAP
    real( dp), dimension(QMAX)  :: knuc,kgro,ktn1 
    real( dp), dimension(QMAX)  :: DIFFMASS1
    real( dp), dimension(QMAX)  :: DIFFMASS2
    real( dp), dimension(QMAX)  :: mwarray

    real( dp),allocatable,dimension(:,:,:)   :: MMX
    real( dp),allocatable,dimension(:,:,:)   :: EXCESS
    real( dp),allocatable,dimension(:,:,:)   :: LOSS
    real( dp),allocatable,dimension(:,:,:)   :: TRANS
    real( dp),allocatable,dimension(:,:,:)   :: KOND
    real( dp),allocatable,dimension(:,:,:)   :: FLUXCM
    real( dp),allocatable,dimension(:,:)     :: FLUXCMC
    real( dp),allocatable,dimension(:,:)     :: FLUXC
    real( dp),allocatable,dimension(:,:,:,:) :: kcoag     
    real( dp),allocatable,dimension(:,:,:)   :: FLUXM
    real( dp),allocatable,dimension(:,:)     :: FLUX
    real( dp),allocatable,dimension(:,:,:)   :: cccond
    real( dp),allocatable,dimension(:,:)     :: depo
    real( dp),allocatable,dimension(:,:)     :: depowall
    real( dp),allocatable,dimension(:,:)     :: wetdep
    real( dp),allocatable,dimension(:,:,:)   :: keffect
    real( dp),allocatable,dimension(:,:)     :: keffectwat
    real( dp),allocatable,dimension(:,:)     :: keffectni,keffectsu
    real( dp),allocatable,dimension(:,:)     :: keffectms,keffectoc
    real( dp),allocatable,dimension(:,:)     :: keffectalk1,keffectalk2


! Allocate aerosol terms

       if (.not. allocated(MMX))          ALLOCATE(MMX(MMAX,IMAX,AMAX))
       if (.not. allocated(EXCESS))       ALLOCATE(EXCESS(MMAX,IMAX,QMAX))
       if (.not. allocated(LOSS))         ALLOCATE(LOSS(MMAX,IMAX,QMAX))
       if (.not. allocated(TRANS))        ALLOCATE(TRANS(MMAX,IMAX,QMAX))       
       if (.not. allocated(KOND))         ALLOCATE(KOND(MMAX,IMAX,QMAX))
       if (.not. allocated(CCCOND))       ALLOCATE(CCCOND(MMAX,IMAX,QMAX))       
       if (.not. allocated(FLUXCM))       ALLOCATE(FLUXCM(MMAX,IMAX,QMAX))
       if (.not. allocated(FLUXCMC))      ALLOCATE(FLUXCMC(MMAX,IMAX))
       if (.not. allocated(FLUXC))        ALLOCATE(FLUXC(MMAX,IMAX))
       if (.not. allocated(kcoag))        ALLOCATE(kcoag(MMAX,MMAX,IMAX,IMAX) )
       if (.not. allocated(FLUXM))        ALLOCATE(FLUXM(MMAX,IMAX,AMAX))
       if (.not. allocated(FLUX))         ALLOCATE(FLUX(MMAX,IMAX))
       if (.not. allocated(depo))         ALLOCATE(depo(MMAX,IMAX))
       if (.not. allocated(depowall))     ALLOCATE(depowall(MMAX,IMAX))       
       if (.not. allocated(wetdep))       ALLOCATE(wetdep(MMAX,IMAX))
       if (.not. allocated(keffect))      ALLOCATE(keffect(MMAX,IMAX,QMAX))       
       if (.not. allocated(keffectwat))   ALLOCATE(keffectwat(MMAX,IMAX))
       if (.not. allocated(keffectni))    ALLOCATE(keffectni(MMAX,IMAX))
       if (.not. allocated(keffectsu))    ALLOCATE(keffectsu(MMAX,IMAX))
       if (.not. allocated(keffectms))    ALLOCATE(keffectms(MMAX,IMAX))
       if (.not. allocated(keffectoc))    ALLOCATE(keffectoc(MMAX,IMAX))
       if (.not. allocated(keffectalk1))  ALLOCATE(keffectalk1(MMAX,IMAX))
       if (.not. allocated(keffectalk2))  ALLOCATE(keffectalk2(MMAX,IMAX))

! Initialize all aerosol terms
! that are output and that are local or locally allocated
       depo(:,:)      = 0._dp 
       depowall(:,:)  = 0._dp
       wetdep(:,:)    = 0._dp
       kcoag(:,:,:,:) = 0._dp
       flux(:,:)      = 0._dp
       fluxm(:,:,:)   = 0._dp
       coags          = 0._dp
       jnuc           = 0._dp
       natot          = 0._dp
       nnuc           = 0._dp
       keffect(:,:,:) = 1._dp
       LVAP(:)        = 0._dp
       FVAP(:)        = 0._dp
       DNVAP(:)       = 0._dp
       DNVAPN         = 0._dp
       NVAPOLDN       = 0._dp
       NVAPO(:)       = 0._dp
       nsv(:)         = 0._dp
       CAT(:)         = 0._dp
       CAT2(:)        = 0._dp
       DIFFMASS1(:)   = 0._dp
       DIFFMASS2(:)   = 0._dp
       MMX(:,:,:)     = 0._dp
       EXCESS(:,:,:)  = 0._dp
       LOSS(:,:,:)    = 0._dp
       KOND(:,:,:)    = 0._dp
       CCCOND(:,:,:)  = 0._dp
       FLUXCM(:,:,:)  = 0._dp
       FLUXCMC(:,:)   = 0._dp
       FLUXC(:,:)     = 0._dp
       CSSULT         = 0._dp
       CSORGT         = 0._dp
       GRSU           = 0._dp
       GRMS           = 0._dp
       GROC           = 0._dp
       GRTOT          = 0._dp


! PARTICLE DEPOSITION
! uses wet diameter DPAW
! particle dry deposition rate [1/s]
! particle deposition only when box is in contact with ground

       if ( (IDEPO.EQ.1).AND.(mbl.ge.2.0*hsta_st) ) then
         CALL depositpar(press,temp,DENSPAR,DPAW,mbl,owf,depo,IMAX)
       endif

       if ( (IDEPO.EQ.2).AND.(mbl.ge.2.0*hsta_st) ) then
         CALL depocanopy(press,temp,DENSPAR,DPAW,mbl,u10,depo,IMAX)
       endif

       if ( (IDEPO.EQ.3).AND.(mbl.ge.2.0*hsta_st) ) then
         CALL depofrough(press,temp,DENSPAR,DPAW,mbl,depo,IMAX)
       endif

       if ( (IDEPO.EQ.4).AND.(mbl.ge.2.0*hsta_st) ) then
         CALL depozhang01(temp,DENSPAR,DPAW,mbl,depo,IMAX)
       endif


! particle deposition rate to chamber walls [1/s]
       if ((ICHAM.EQ.1).AND.(IDEPO.GE.1)) then
         CALL depositwall(press,temp,DPA,           &
                                      DILPAR,V_CHAM,S_SED,S_DIF,   &
                                      depowall,IMAX)  
       endif

! particle wet scavenging rate [1/s]
       if (rain .gt. 0.0) then
         ! bulk aerosol
         if (IWETD.EQ.1) then
           wetsc = wetscavbulk(rain)
           do M=AI,CS
             do I=1,IMAX
               wetdep(M,I) = wetsc
             enddo
           enddo 
         endif
         ! size-resolved aerosol
         if (IWETD.EQ.2) CALL wetscavsize(IMAX,DPAW,press,temp,DENSPAR,rain,wetdep)
       endif


! PARTICLE COAGULATION
! coagulation of particles

       if (ICOAG .ge. 1) then
       
! coagulation coefficient kcoag(M,O,I,J)
! (also used for interstitial aerosol scavenging)
           CALL coagulation_coeff(temp,DENSPARW,DPAW,kcoag,IMAX)

! coagulation by Brownian diffusion
           CALL coagulation(temp,DTIME,DENSPARW,DPAW,VPT,N,MASS,  &
                          IMAX,IAG,kcoag,coags,fluxm,flux)

! compute updated coagulation target classes      
           CALL coagulation_target(VPT,IMAX,IAG)

       endif


! GAS-PARTICLE EQUILIBRIUM
! for amine-hno3 system
       if (ICHAM .EQ. 1) THEN
          KPEQ=KP_NIT
       else
          KPEQ=6.84E21_dp  ! NH4NO3 at 298 K
       endif

! NUCLEATION (first)
! Nucleation of particles, N(1)
!  JNUC: nucleation rate J in m^-3s^-1
       if (INUC == 1) then
         CALL nucleation(INUCMEC,cair,NVAP(A_SUL),NVAP(A_NH4),    &
                         NVAP(A_AMI),NVAP(A_NIT),                 &
                         NVAP(A_OR2),NVAP(A_IO3),                 &
                         temp,RH,jrno2,CAMI,                      &
                         KPEQ,fnuc,coags,daytime,lat_deg,natot,jnuc)
                         
         if ((IDEB == 1) .and. (hour_timemin == 120.)) then
            write(12,*) 'Jnuc', jnuc       
         endif
       endif


! CONDENSATION     
! Multi-component condensation
! condensation of vapours to particles:
! - sulphuric acid           NVAP(A_SUL)
! - methane sulphonic acid   NVAP(A_MSA)
! - iodic acid               NVAP(A_IO3)
! - ammonium nitrate         NVAP(A_NH4),NVAP(A_NIT)
! - amine nitrate            NVAP(A_AMI),NVAP(A_NIT)
! - ammonium chloride        NVAP(A_NH4),NVAP(A_CHL)
! - BSOA-1, BSOA-2           NVAP(A_OR1),NVAP(A_OR2)
! - BSOA-3                   NVAP(A_OR3)
! - ASOA-1, ASOA-2           NVAP(A_OR4),NVAP(A_OR5)
! - ASOA-3                   NVAP(A_OR6)
! - PIOA-1, PIOA-2           NVAP(A_OR7),NVAP(A_OR8)
! - PIOA-3                   NVAP(A_OR9)

! 1) Compute Kelvin Effect and Raoult Effect

       if (ICOND .EQ. 1) then

         ! Set Kelvin effect to 1 before any calculation

           keffectsu(NU:CS,1:IMAX)   = 1.0_dp
           keffectms(NU:CS,1:IMAX)   = 1.0_dp
           keffectoc(NU:CS,1:IMAX)   = 1.0_dp
           keffectalk1(NU:CS,1:IMAX) = 1.0_dp
           keffectalk2(NU:CS,1:IMAX) = 1.0_dp
           keffectni(NU:CS,1:IMAX)   = 1.0_dp

         ! Calculate Kelvin effect if IKELV==1 

         if (IKELV==1) then
         !   First for sulfuric acid, MSA and iodic acid
           if ((ICONS == 1).or.(ICONS == 2)) then
             CALL kelvin_sulf(DPA,temp,keffectsu,IMAX)
             CALL kelvin_msap(DPA,temp,keffectms,IMAX)
           endif
         !  Second for organic vapour
           if (ICONO .EQ. 1) then
             if (INANO.EQ.1) then
               CALL nano_koehler(DPA,temp,surfin,surf_org,M_oc(1),DENOC, &
                               MASS,keffectwat,keffectoc,IMAX)
             else   
               CALL kelvin_org(DPA,temp,surfin,surf_org,M_oc(1),DENOC,MASS,keffectoc,IMAX)
             endif   
         ! primary emitted alkanes with nC=21 and nC=26
             CALL kelvin_alkane(DPA,temp,21,gamma_oc7,keffectalk1,IMAX)             
             CALL kelvin_alkane(DPA,temp,26,gamma_oc8,keffectalk2,IMAX) 
         ! re-calculate molar fraction of primary emitted OC components
         ! (gamma_oc7 and gamma_oc8)
             do M=NU,CS
                 sum_gamma(M) =(gamma_oc1_m(M)*M_oc(1)+gamma_oc2_m(M)*M_oc(2)  &
                              + gamma_oc3_m(M)*M_oc(3)+gamma_oc4_m(M)*M_oc(4)  &
                              + gamma_oc5_m(M)*M_oc(5)+gamma_oc4_m(M)*M_oc(6)  &
                              + gamma_oc7_m(M)*M_oc(7)+gamma_oc8_m(M)*M_oc(8)  &
                              + gamma_oc9_m(M)*M_oc(9)    )
                 gamma_oc7_tmp = MORG7TOT(M)/(max(MORGCTOT(M),massmin))
                 gamma_oc8_tmp = MORG8TOT(M)/(max(MORGCTOT(M),massmin))
                 gamma_oc7(M)=(gamma_oc7_tmp/M_oc(7)) * sum_gamma(M)
                 gamma_oc8(M)=(gamma_oc8_tmp/M_oc(8)) * sum_gamma(M)

                 !print *,"aerosol_solver gamma_oc7 ",M,sum_gamma(M),gamma_oc7(M)
                 !print *,"aerosol_solver gamma_oc8 ",M,sum_gamma(M),gamma_oc8(M)                        
             end do
           endif   
         !   Third for amine-nitrate          
           if (ICONA .EQ. 1) then
               CALL kelvin_nit(DPA,temp,keffectni,IMAX)
           endif
         !   Fourth for ammonium nitrate & ammonium chloride
           if (ICONX .EQ. 1) then
               CALL kelvin_nit(DPA,temp,keffectni,IMAX)
           endif

         endif
       endif


! 2) Calculate saturation vapour concentrations, 
!    condensation coefficients, excess and loss

! Kelvin effect into one variable

       if (ICOND .EQ. 1) then

         do M=NU,CS
           do I=1,IMAX
         ! A_SUL = 1, A_MSA = 2, A_NIT = 3
         ! A_AMI = 4, A_NH4 = 5, A_IO3 = 6
             keffect(M,I,1)  = keffectsu(M,I)
             keffect(M,I,2)  = keffectms(M,I)
             keffect(M,I,3)  = keffectni(M,I)           
             keffect(M,I,4)  = keffectni(M,I)
             keffect(M,I,5)  = keffectni(M,I)
             keffect(M,I,6)  = keffectsu(M,I)
         ! soa components
             keffect(M,I,7)  = keffectoc(M,I)
             keffect(M,I,8)  = keffectoc(M,I)
             keffect(M,I,9)  = 1.0_dp
             keffect(M,I,10)  = keffectoc(M,I)
             keffect(M,I,11) = keffectoc(M,I)
             keffect(M,I,12) = 1.0_dp
             keffect(M,I,13) = keffectalk1(M,I)
             keffect(M,I,14) = keffectalk2(M,I)
             keffect(M,I,15) = 1.0_dp
         ! hcl
             keffect(M,I,16) = keffectni(M,I)
           end do
         end do

! Calculate totals for NH4, SO4, NO3, CHL
         CTSO4=0._dp
         CTNH4=0._dp
         CTNIT=0._dp
         CTCHL=0._dp
         do M=NU,CS
           do I=1,IMAX
             if (M.LT.CS) CTSO4=CTSO4+mass(M,I,A_SUL)
             CTNH4=CTNH4+mass(M,I,A_NH4)
             CTNIT=CTNIT+mass(M,I,A_NIT)
             CTCHL=CTCHL+mass(M,I,A_CHL)
           end do
         end do

         ! add gas-phase concentration in ng/m^3
         CTSO4=CTSO4+nvap(A_SUL)*1.e-3*molec2ug(M_H2SO4)
         CTNH4=CTNH4+nvap(A_NH4)*1.e-3*molec2ug(M_nh3)
         CTNIT=CTNIT+nvap(A_NIT)*1.e-3*molec2ug(M_nit)
         CTCHL=CTCHL+nvap(A_CHL)*1.e-3*molec2ug(M_hcl)

         ! sum up the water content for APD
         ! sum lwc (ng/m3)
         sumlwc=0._dp
         do M=NU,CS
           do I=1,IMAX
              sumlwc=sumlwc+mass(m,i,a_wat)
           end do
         end do


! coupled to MESA for APD (Analytical Predictor of Dissolution)
         CALL condensloss(temp,press,RH,DPAW,N,IMAX,MASS,alphanit,fcoj,KPEQ,               &
                          DENOC, M_oc,nmo,foc,hvap,csat0,                                  &
                          ctso4,ctnh4,ctnit,ctchl,NVAP,                                    &
                          MORG1TOT,MORG2TOT,MORG3TOT,MORG4TOT,MORG5TOT,                    &
                          MORG6TOT,MORG7TOT,MORG8TOT,MORG9TOT,                             &
                          MORGCTOT,MECBCTOT,MDUSTTOT,MTOT,                                 &
                          keffect,henry_eff_hno3,henry_eff_hcl,svmc_min_hno3,              &
                          svmc_min_hcl,svmc_min_nh3, flag_dissolution,                     &
                          Keq_nh4no3_0,Keq_nh4cl_0,sumlwc, DTIME,                          &
                          GRSU,GRMS,DCSU,GROC,DCORG,cccond,                                &
                          nsv,EXCESS,LOSS,TRANS,csate)

       endif

! 3) Sum up Condensation Sink

       if (ICOND .EQ. 1) then
         if ((ICONS == 1).or.(ICONS == 2)) then
            do M=NU,CS
             do I=1,IMAX
               CSSULT=CSSULT+(cccond(M,I,A_SUL)*N(M,I)/(DCSU*4.*pi))           
             end do
            end do
         endif
         if (ICONO .EQ. 1) then         
            do M=NU,CS
             do I=1,IMAX              
               CSORGT=CSORGT + ( (cccond(M,I,A_OR1)*N(M,I)/(DCORG*4.*pi)) + &
               (cccond(M,I,A_OR2)*N(M,I)/(DCORG*4.*pi))+(cccond(M,I,A_OR3)*N(M,I)/(DCORG*4.*pi)) +&
               (cccond(M,I,A_OR4)*N(M,I)/(DCORG*4.*pi))+(cccond(M,I,A_OR5)*N(M,I)/(DCORG*4.*pi)) +&
               (cccond(M,I,A_OR6)*N(M,I)/(DCORG*4.*pi))+(cccond(M,I,A_OR7)*N(M,I)/(DCORG*4.*pi)) +&
               (cccond(M,I,A_OR8)*N(M,I)/(DCORG*4.*pi))+(cccond(M,I,A_OR9)*N(M,I)/(DCORG*4.*pi)) )/9
             end do
            end do
         endif
         ! Calculate Condensation Sink of organic vapour
         if ( ICONO == 1) CSORGT=4*pi*DCORG*CSORGT
         if ((ICONS == 1).or.(ICONS == 2))  CSORGT=CSORGT + (4*pi*DCSU*CSSULT)
       endif

! 4) Compute condensation flux
  
       if (ICOND .EQ. 1) then
         CALL condensation(VPT,DTIME,cccond,excess,loss,trans,         & 
                           N,IMAX,MASS,LVAP,FVAP,FLUXCM,FLUXC,FLUXCMC ) 
       endif

! 5) JACOBSON APC SHEME

!   M.J. Jacobson, p.544, eq. 16.70
!   Compute Ca_sum for preservation of mass-balance
!   Ctot:CATxxxO  [ng/m^3] --> [ug m^-3]
!   1) calculate Ca_sum (MOLD)
      if (ICOND .EQ. 1) then
         do M=NU,CS
          do I=1,IMAX
            do Q=1,QMAX
              CAT(Q)=CAT(Q)+mass(M,I,Q)*1.e-3
            end do                                            
          end do
         end do
         do Q=1,QMAX
           nvapo(Q)=NVAP(Q)
         end do

     !   write(6,*) 'gdesolver nvap 0',nvapo(A_NH4),nvap(A_NH4)

! compute new gas phase concentration after condensation
         CALL apc_update_gasc(M_oc,nvapo,DTIME,CAT,FVAP,LVAP,DNVAP,NVAP)

     !   write(6,*) 'gdesolver nvap 1',nvapo(A_NH4),nvap(A_NH4)

! compute new condensation rate KOND [molec/m^3s]
         do M=NU,CS
          do I=1,IMAX
            do Q=1,QMAX 
              KOND(M,I,Q)=N(M,I)*cccond(M,I,Q)*(NVAP(Q)-nsv(Q)*keffect(M,I,Q))
            end do
          end do
        end do

! for ammonium nitrate condensation?
        if (ctnh4 > ctnit) then
           KOND(:,:,a_nh4)=KOND(:,:,a_nh4)*0.4
        endif

! If incloud=1 no condensation to CS mode
        if (incloud .eq. 1) then
           KOND(CS,:,:)   = 0._dp
           FLUXCM(CS,:,:) = 0._dp
        endif

! compute new particle mass concentrations after condensation
        CALL apc_update_pmass(DENOC,DENSPARW,DTIME,IMAX,M_oc,KOND,FLUXCM,MASS,MMX)


! 6) Mass balance correction 
! 2020-06-25 added by MSK
! 2021-04-04 revised by MSK
!      We place here the upper limit for the new mass concentration
!      to bound the final gas and aerosol concentraton between
!      0 and C_tot = CAT + NVAPO
! Correct gas phase concentration again for mass balance reasons
! for the secondary inorganic aerosols

! molecular weights
        mwarray(A_SUL) = M_H2SO4
        mwarray(A_MSA) = M_msa
        mwarray(A_NIT) = M_nit
        mwarray(A_AMI) = M_nit
        mwarray(A_NH4) = M_nh3
        mwarray(A_IO3) = M_hio3
        mwarray(A_CHL) = M_hcl
        mwarray(A_OR1) = M_oc(1)
        mwarray(A_OR2) = M_oc(2)
        mwarray(A_OR3) = M_oc(3)
        mwarray(A_OR4) = M_oc(4)
        mwarray(A_OR5) = M_oc(5)
        mwarray(A_OR6) = M_oc(6)
        mwarray(A_OR7) = M_oc(7)
        mwarray(A_OR8) = M_oc(8)
        mwarray(A_OR9) = M_oc(9)

! aerosol mass summation
! use DIFFMASS1 (old mass) and DIFFMASS2 (prelim. new mass)
! CAT2 in ug/m3;  DIFFMASS in ng/m3

        do M=NU,CS
          do I=1,IMAX
            do Q=1,QMAX
              CAT2(Q)     =CAT2(Q)     +mmx(M,I,Q)*1.e-3
              DIFFMASS1(Q)=DIFFMASS1(Q)+mass(M,I,Q)
              DIFFMASS2(Q)=DIFFMASS2(Q)+mmx(M,I,Q)
            enddo
          end do
        end do

      !      print *,'nvap1',nvap(a_nh4),nvapo(a_nh4),CAT(a_nh4),CAT2(a_nh4)

! new gas phase concentration
        do Q=1,QMAX

! special for semi-volatile secondary inorganic
! A_NH4 (NH4+)
          !if (Q==A_NH4)  then

      !      nvap(a_nh4) = nvapo(a_nh4) + (CAT(a_nh4)-CAT2(a_nh4))  * &
      !                    1.e6*(1./molec2ug(M_nh3))
      !               print *,'nvap2',nvap(a_nh4)
! A_NIT (NO3-)
          !else 
          if ( (Q==A_NIT).and.(ICONX==1) ) then
            nvap(a_nit) = nvapo(a_nit) + (CAT(a_nit)-CAT2(a_nit))  * &
                          1.e6*(1./molec2ug(M_nit))
! A_CHL (Cl-)
          else if (Q==A_CHL) then
            nvap(a_chl) = nvapo(a_chl) + (CAT(a_chl)-CAT2(a_chl))  * &
                          1.e6*(1./molec2ug(M_hcl))
          else
! all others only if CAT > CAT2
            if ( (CAT2(Q).lt.CAT(Q)).and.(CAT2(Q).gt.1.e-7) ) then
              nvap(Q) = nvapo(Q) + ( CAT(Q)-CAT2(Q) )*1.e6*(1./molec2ug( mwarray(Q) ))
            endif
          endif

          nvap(Q) = max(nvap(Q),0._dp)


! C_tot (ng/m3)
! gas (molec/m3) * 1.e-6 = gas (molec/cm3)
! gas (molec/cm3) * molecug = gas (ug/m3)
! gas (ug/m3) * 1.e3 = gas (ng/m3)
          DIFFMASS1(Q)=DIFFMASS1(Q)+nvapo(Q)*1.e-3*molec2ug(mwarray(Q))
          DIFFMASS2(Q)=DIFFMASS2(Q)+nvap(Q) *1.e-3*molec2ug(mwarray(Q))

! final aerosol mass concentration
      !   if (Q==5) then
      !     print *,'nh4 diffmass',DIFFMASS1(Q),DIFFMASS2(Q)
      !   else
           if (DIFFMASS2(Q).gt.DIFFMASS1(Q)) then
             do M=NU,CS
               do I=1,IMAX
                  mmx(m,i,q)=mmx(m,i,q)* (DIFFMASS1(q)/DIFFMASS2(q))
               end do
             end do
           endif
      !  endif

        enddo

! 2021-04-04 end final mass correction


!!----------------------------------------------------------
!! PNG solver: equilibration of NH3 (ICONW==2)
!!    M.Z. Jacobson,  AEROSOL SCIENCE TECHNOLOGY, VOL. 39, 
!!          NO. 2, 92-103, DOI: 10.1080/027868290904546, 2005.
!!   Newton-Raphson iteration
!!   x       : NVAPO(A_NH4)            [mlc/m^3]
!!   ctot    : ctnh4                   [ng/m^3]
!!   cioncharge = concentration of all 
!!           ions weighted by charge   [mlc/m^3]
!!   H'(NH3) : henry_eff_nh3           [-]
!!   K'(NH3) : kh_nh3                  [-]
!! NOTE: does not give stable results if there are
!!       emissions of ammonia
!!----------------------------------------------------------
!        if (ICONW .EQ. 2) then
!       ! compute the tolerance
!           mtol = 1.e5_dp * machineps(0.5_dp)
!           deriv_zero = .false.
!       ! iteration of new NH3 concentration
!           CALL newton( nvapo(A_NH4),henry_eff_nh3,kh_nh3,   &
!                        cioncharge, ctnh4, IMAX,mtol,mtol,   &
!                        .true.,deriv_zero, xout )
!           write(*,'(a33, ES12.4)') 'the solution of f(x) = 0 is x = ', xout
!       ! set new NH3 concentration
!           NVAP(A_NH4) = max(xout,0._dp)
!           sf=1.0_dp
!       ! insert NVAP(A_NH4) in Eq. (31) to overwrite MMX(A_NH4)
!           do M=NU,CS
!             do I=1,IMAX
!               ! net charge must not be positive
!               ci   = min(cioncharge(m,i), 0._dp)
!               mnh4 = -ci * NVAP(A_NH4)                 *  &
!                       henry_eff_nh3(m,i)*kh_nh3(m,i)*sf    /  &
!                       ( NVAP(A_NH4)*henry_eff_nh3(m,i)*kh_nh3(m,i)*sf + 1._dp )
!               ! convert molec/m^3 to ng/m^3
!               MMX(m,i,A_NH4) = max(mnh4*1.e-3_dp*molec2ug(M_nh3), 0._dp)
!               !print *,'mmx',m,i,mass(m,i,A_NH4),MMX(m,i,A_NH4)
!            end do
!           end do
!         endif  ! ICONW2
    
      endif   ! ICOND=1



! NUCLEATION
! Jacobson (2002), Equation (33)
!  mass transfer coefficient knuc [s^-1]
!  JNUC: nucleation rate J in m^-3s^-1
! not used at the moment
       kgro(:)=0.0
       knuc(:)=0.0
       ktn1(:)=0.0
       if (INUC .EQ. 1) then
         if (jnuc > 0) then
           do Q=1,QMAX
             if ( nvap(Q) > nsv(Q)*keffect(NU,1,Q) )  then
               knuc(Q) = jnuc*DENSPAR(NU,1)*vvap*nnuc / MB
               knuc(Q) = knuc(Q) * 1/(nvap(Q) - nsv(Q)*keffect(NU,1,Q) )
             else
               knuc(Q) = 0.0
            endif
           end do
         else
           knuc(:) = 0.0           
         endif
       endif   
       if (ICOND .EQ. 1) then
         do Q=1,QMAX 
           kgro(Q) = cccond(NU,1,Q)*n(NU,1)
         end do  
       endif
       do Q=1,QMAX
         ktn1(Q) = knuc(Q) + kgro(Q)
       end do

! -----------------------------------------------------
! UPDATE MASS AND NUMBER CONCENTRATION BY ALL PROCESSES
! -----------------------------------------------------

! CONDENSATION
       if (IDEB == 1) then
        DEB_NTOT =0.
        DEB_CTOTS1=0.
        DEB_CTOTA1=0.
        DEB_CTOTO11=0.
        DEB_CTOTO21=0.
        DEB_MTOT =0.
        do M=NU,CS
          do I=1,IMAX
            DEB_NTOT=DEB_NTOT+n(M,I)
            DEB_CTOTS1=DEB_CTOTS1+mass(M,I,A_SUL)
            DEB_CTOTA1=DEB_CTOTA1+mass(M,I,A_NH4)
            DEB_CTOTO11=DEB_CTOTO11+mass(M,I,A_OR1)
            DEB_CTOTO21=DEB_CTOTO21+mass(M,I,A_OR2)
            do K=1,AMAX 
              DEB_MTOT=DEB_MTOT+mass(M,I,K)
            end do
          end do
        end do
        DEB_NTOT = DEB_NTOT*1.e-6 
       endif    

       
       if (ICOND == 1) then 
         if (IDEB == 1) then
           write(12,*) 'Number and mass conservation for condensation:'
           write(12,fmt='(a,f17.5)') '    NTOT before condensation          [#/cm^3]', DEB_NTOT
           DEB_CTOTS1 = DEB_CTOTS1  + nvapo(A_SUL)*1.e-3*molec2ug(M_H2SO4)
           DEB_CTOTA1 = DEB_CTOTA1  + nvapo(A_NH4)*1.e-3*molec2ug(M_nh3)
           DEB_CTOTO11 = DEB_CTOTO11+ nvapo(A_OR1)*1.e-3*molec2ug(M_oc(1))
           DEB_CTOTO21 = DEB_CTOTO21+ nvapo(A_OR2)*1.e-3*molec2ug(M_oc(2))
         endif

         do M=NU,CS
          do I=1,IMAX
           do Q=1,QMAX
             MASS(M,I,Q)=MMX(M,I,Q)
           end do
           ! move core particles
           MASS(M,I,A_XXX) = mass(M,I,A_XXX) + fluxcmc(M,I)                             
           MASS(M,I,A_XXX) = max(mass(M,I,A_XXX),0.0_dp)
           ! numbers                          
           N(M,I)=N(M,I)+fluxc(M,I)           
          end do
         end do


         if (IDEB == 1) then
           ! test of number conversation for condensation
           DEB_NTOT=0.
           DEB_CTOTS2=0.
           DEB_CTOTA2=0.
           DEB_CTOTO12=0.
           DEB_CTOTO22=0.
           do M=NU,CS
             do I=1,IMAX
               DEB_NTOT=DEB_NTOT+n(M,I)
               DEB_CTOTS2=DEB_CTOTS2+mass(M,I,A_SUL)
               DEB_CTOTA2=DEB_CTOTA2+mass(M,I,A_NH4)
               DEB_CTOTO12=DEB_CTOTO12+mass(M,I,A_OR1) 
               DEB_CTOTO22=DEB_CTOTO22+mass(M,I,A_OR2) 
             end do
           end do

           DEB_CTOTS2  = DEB_CTOTS2+  nvap(A_SUL)*1.e-3*molec2ug(M_H2SO4)
           DEB_CTOTA2  = DEB_CTOTA2+  nvap(A_NH4)*1.e-3*molec2ug(M_nh3)
           DEB_CTOTO12 = DEB_CTOTO12+ nvap(A_OR1)*1.e-3*molec2ug(M_oc(1))
           DEB_CTOTO22 = DEB_CTOTO22+ nvap(A_OR2)*1.e-3*molec2ug(M_oc(2))
           write(12,fmt='(a,f17.5)') '    NTOT after condensation           [#/cm^3]', DEB_NTOT*1e-6
           write(12,fmt='(a,f13.5)') '    CTOT(sulfate) before condensation [ng/m^3]', DEB_CTOTS1
           write(12,fmt='(a,f13.5)') '    CTOT(sulfate) after condensation  [ng/m^3]', DEB_CTOTS2
           write(12,fmt='(a,f13.5)') '    CTOT(NH4  ) before condensation   [ng/m^3]', DEB_CTOTA1
           write(12,fmt='(a,f13.5)') '    CTOT(NH4  ) after condensation    [ng/m^3]', DEB_CTOTA2
           write(12,fmt='(a,f13.5)') '    CTOT(SOA-1) before condensation   [ng/m^3]', DEB_CTOTO11
           write(12,fmt='(a,f13.5)') '    CTOT(SOA-1) after condensation    [ng/m^3]', DEB_CTOTO12
           write(12,fmt='(a,f13.5)') '    CTOT(SOA-2) before condensation   [ng/m^3]', DEB_CTOTO21
           write(12,fmt='(a,f13.5)') '    CTOT(SOA-2) after condensation    [ng/m^3]', DEB_CTOTO22
              !     write(6,*) 'gdesolver NH4 balnce',DEB_CTOTA1,DEB_CTOTA2

         endif
       endif

         ! write(6,*) 'gedsolver nvap 2',nvapo(A_NH4),nvap(A_NH4)


! NUCLEATION
       if (INUC .EQ. 1) then 
         ! NNUC: volume ratio particle/vapor
         CALL nucleationratio(DPA(NU,1),vvap,nnuc) 
       ! Compute number and mass concentration after nucleation
         ! organic/sulfuric acid nucleation    
         if ( (INUCMEC.EQ.8).OR.(INUCMEC.EQ.9).OR.(INUCMEC.EQ.10).OR.(INUCMEC.EQ.12) ) then
            natot=1.0 
            NVAP(A_OR2)=nvap(A_OR2)-jnuc*nnuc*natot*DTIME
            MASS(NU,1,A_OR2) = MASS(NU,1,A_OR2)                      +       &
                       natot*MAH*jnuc*DTIME*CONVM
         endif
         ! amine/nitric acid nucleation 
         if (INUCMEC .EQ. 6) then
            NVAPOLDN=nvap(A_AMI)
            !do I=1,IMAX
              MASS(NU,1,A_AMI) = mass(NU,1,A_AMI)                   +       &
                       natot*MAN*jnuc*DTIME*CONVM
              MASS(NU,1,A_NIT) = mass(NU,1,A_NIT)                   +       &
                       natot*MAN*jnuc*DTIME*CONVM
              N(NU,1)=N(NU,1)+jnuc*DTIME
           
              NVAP(A_AMI)=nvap(A_AMI)-jnuc*nnuc*natot*DTIME
              NVAP(A_NIT)=nvap(A_NIT)-jnuc*nnuc*natot*DTIME  
            !enddo         
            DNVAPN=DNVAPN+(NVAPOLDN-nvap(A_AMI))
         ! iodic acid / sulfuric acid nucleation    
         else if (INUCMEC .EQ. 14) then
            NVAP(A_IO3)=nvap(A_IO3)-jnuc*nnuc*natot*DTIME
            MASS(NU,1,A_IO3) = MASS(NU,1,A_IO3)                     +       &
                       natot*MIO*jnuc*DTIME*CONVM
            ! for the sulfuric acid / iodic acid nucleation
            nvap(A_SUL) = nvap(A_SUL)-jnuc*nnuc*DTIME
            mass(NU,1,A_SUL) = mass(NU,1,A_SUL)                     +       &
                      natot*MB*jnuc*DTIME*CONVM
            N(NU,1)=N(NU,1)+jnuc*DTIME
         ! iodic acid nucleation    
         else if (INUCMEC .EQ. 15) then
            NVAP(A_IO3)=nvap(A_IO3)-jnuc*nnuc*natot*DTIME
            MASS(NU,1,A_IO3) = MASS(NU,1,A_IO3)                     +       &
                       natot*MIO*jnuc*DTIME*CONVM
            N(NU,1)=N(NU,1)+jnuc*DTIME
         else
       !   sulfuric acid nucleation
       !       new check on NVAP(A_SUL), 17.02.2013 
             if( jnuc*nnuc*DTIME > nvap(A_SUL) ) then
                jnuc = nvap(A_SUL)/(nnuc*DTIME)
                nvap(A_SUL) = 0.0
             else
                nvap(A_SUL) = nvap(A_SUL)-jnuc*nnuc*DTIME
             endif    
             mass(NU,1,A_SUL) = mass(NU,1,A_SUL)                    +       &
                       natot*MB*jnuc*DTIME*CONVM
       ! New Number concentration in NU,1    
           n(NU,1)=n(NU,1)+jnuc*DTIME
       ! Jacobson (2002), Equation (34)
       ! typically produces too little (NU,1) particles
       !    N(NU,1) = N(NU,1) + max( (                                       & 
       !       (1/molec2ug(M_H2SO4))*( MASS(NU,1,A_SUL)-MMX(NU,1,A_SUL) )*   &
       !       1.e3*(MB/(DENV*VVAP*NNUC)) *knuc(a_sul)/ktn1(a_sul)), 0.0)+&
       !       max(((1/molec2ug(M_oc2  ))*(MASS(NU,1,A_OR2)-MMX(NU,1,A_OR2))*&
       !      1.e3* (MAH/(DENOC*VVAP*NNUC)) * knuc(a_or2)/ktn1(a_or2) ), 0.0) 
         endif
       endif



! COAGULATION
       if (ICOAG .ge. 1) then
         if (IDEB == 1) then
           write(12,*) 'Mass conservation for coagulation:'
           write(12,fmt='(a,f13.5)') '    MTOT before coagulation           [ng/m^3]', DEB_MTOT
         endif       

         if (incloud .eq. 1) then
    ! Collection of interstatial particles in the non-precipitating cloud
    ! where the coarse mode droplets act as collectors and the particles
    ! in NU to AS mode are captured by Brownian diffusion and removed.
    ! Particles smaller than 10 nm will be scavenged in few minutes.
    ! Seinfeld & Pandis book 1998, ch. 15.7.4
    ! Limit this to not activated particles (Dp_dry < 100 nm)

           do M=NU,AS
             do I=1,IMAX
               if (DPA(M,I).lt.1.e-7) then
                 do J=1,IMAX
                   coagscav = kcoag(M,CS,I,J) * N(CS,J)
                   do K=1,AMAX
                     MASS(M,I,K)= MASS(M,I,K)             -      &
                              MASS(M,I,K)*coagscav*DTIME
                   end do
                   N(M,I)=N(M,I) - N(M,I)*coagscav*DTIME
               
               !  write(6,*) 'Kcoag',M,I,CS,J,DPAW(M,I),coagscav
                 enddo
               endif
             enddo
           enddo

         else

    ! Particle coagulation by Brownian diffusion

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

             !  write(6,*) 'Fcoag',M,I,N(M,I),FLUX(M,I)
            end do
           end do

         endif


    ! Test of mass conversation for coagulation
         if (IDEB == 1) then         
           DEB_MTOT=0.
           do M=NU,CS
             do I=1,IMAX
               do K=1,AMAX 
                 DEB_MTOT=DEB_MTOT+mass(M,I,K)
               end do
             end do
           end do
           write(12,fmt='(a,f13.5)') '    MTOT after coagulation            [ng/m^3]', DEB_MTOT
         endif

       endif


! DRY DEPOSITION
       if ((IDEPO.EQ.1).OR.(IDEPO.EQ.2).OR.(IDEPO.EQ.3)) then
         do M=NU,CS
          do I=1,IMAX
           do K=1,AMAX
             MASS(M,I,K)=MASS(M,I,K)                    -        &
                   MASS(M,I,K)*DEPO(M,I)*DTIME
           end do           
           N(M,I)=N(M,I)-N(M,I)*DEPO(M,I)*DTIME
          end do
         end do
       endif


! WET SCAVENGING
       if (IWETD.GE.1) then
         do M=NU,CS
          do I=1,IMAX
            do K=1,AMAX
             MASS(M,I,K)=MASS(M,I,K)                     -       &
                    MASS(M,I,K)*wetdep(M,I)*DTIME
            end do
            N(M,I)=N(M,I)-N(M,I)*wetdep(M,I)*DTIME
          end do
         end do
       endif


! DEPOSITION TO CHAMBER WALLS
       if ((ICHAM.EQ.1).AND.(IDEPO.GE.1)) then
         do M=NU,CS
          do I=1,IMAX
           do K=1,AMAX
             MASS(M,I,K)=MASS(M,I,K)                     -       &
                   MASS(M,I,K)*depowall(M,I)*DTIME   
           end do
           N(M,I)=N(M,I)-N(M,I)*depowall(M,I)*DTIME
          end do
         end do
       ENDIF



 ! Calculate Growth Rate of 1nm particle [nm/hr]
       ! According to Eq. 2-5 in Verheggen and Mozurkevich, ACP,6,2927-2942,2006
       ! GR=GR*(NVAP-NSVAP)*3600.*1.e9
       ! GR[m/s]=GR[m4/(s*molec)]*(NVAP-NSVAP)[molec/m3]  
       if ((ICONS == 1).or.(ICONS == 2)) then
         GRSU=GRSU*(NVAP(A_SUL)-nsv(A_SUL))*3600.*1.e9*1.
         GRMS=GRMS*(NVAP(A_MSA)-nsv(A_MSA))*3600.*1.e9*1.
       endif
       GROC=GROC*(NVAP(A_OR1)-nsv(A_OR1))*3600.*1.e9*1.    !mass fraction = 1
       GRTOT=GROC+GRSU+GRMS



! end update MASS and N


! Deallocate aerosol terms
       deallocate(mmx)
       deallocate(excess)
       deallocate(loss)
       deallocate(trans)
       deallocate(kond)
       deallocate(cccond)
       deallocate(fluxcm)
       deallocate(fluxcmc)
       deallocate(fluxc)
       deallocate(kcoag)
       deallocate(fluxm)
       deallocate(flux)
       deallocate(depo)
       deallocate(depowall)
       deallocate(wetdep)
       deallocate(keffect)
       deallocate(keffectwat)
       deallocate(keffectoc)
       deallocate(keffectni)
       deallocate(keffectsu)
       deallocate(keffectms)
       deallocate(keffectalk1)
       deallocate(keffectalk2)

  end subroutine aerosol_solver


  subroutine interface_mosaic(IMAX,temp,press,RH, alphanit,fcoj,    &
                      mass,nvap,mpt,dpaw,n,hyst,                    &
                      henry_eff_hno3,henry_eff_hcl,henry_eff_nh3,   &
                      svmc_min_hno3,svmc_min_hcl,svmc_min_nh3,      &
                      cioncharge,kh_nh3,flag_dissolution,           &
                      Keq_nh4no3_0,Keq_nh4cl_0 )
    !----------------------------------------------------------------------
    !     
    ! Interface for MOSAIC solver
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      interface to the MESA library
    !
    !      interface
    !      ---------
    !         output:
    !            henry_eff_hno3      per each bin
    !            henry_eff_nh3       per each bin
    !            svmc_min_hno3       per each bin
    !            svmc_min_nh3        per each bin
    !            cioncharge          per each bin
    !            kh_nh3              per each bin
    !            flag_dissolution    per each bin
    !            HYST ???
    ! 
    !      method
    !      ------
    !      dissolution and growth of semi-volatile inorganic compounds
    !
    !      
    !      external
    !      --------
    !      none
    !      
    !      reference
    !      ---------
    !      none
    !
    !------------------------------------------------------------------

! interface to mosaic
      use gde_mosaic_therm,           only : mosaic_control
      use gde_condensation,           only : nitcondens

    implicit none

    ! input
      integer, intent(in)                        :: IMAX
      real( dp), intent(in)                      :: temp
      real( dp), intent(in)                      :: press
      real( dp), intent(in)                      :: RH
      real( dp), intent(in)                      :: alphanit
      real( dp), intent(in)                      :: fcoj

      real( dp), dimension(QMAX),intent(in)            :: NVAP
      real( dp), dimension(MMAX,IMAX,AMAX),intent(in)  :: MASS
      real( dp), dimension(MMAX,IMAX),intent(in)       :: MPT    
      real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPAW
      real( dp), dimension(MMAX,IMAX),intent(in)       :: N

    ! inout
     real( dp),dimension(MMAX,IMAX),intent(inout)      :: HYST    !public

    ! output
     real( dp),dimension(MMAX,IMAX),intent(out)  :: henry_eff_hno3
     real( dp),dimension(MMAX,IMAX),intent(out)  :: henry_eff_hcl
     real( dp),dimension(MMAX,IMAX),intent(out)  :: henry_eff_nh3
     real( dp),dimension(MMAX,IMAX),intent(out)  :: svmc_min_hno3
     real( dp),dimension(MMAX,IMAX),intent(out)  :: svmc_min_hcl
     real( dp),dimension(MMAX,IMAX),intent(out)  :: svmc_min_nh3
     real( dp),dimension(MMAX,IMAX),intent(out)  :: kh_nh3
     real( dp),dimension(MMAX,IMAX),intent(out)  :: cioncharge
     integer,  dimension(MMAX,IMAX),intent(out)  :: flag_dissolution
     real( dp), intent(out)                      :: Keq_nh4no3_0
     real( dp), intent(out)                      :: Keq_nh4cl_0

    ! local
     integer                                    :: m,i,q
     real( dp), dimension(MMAX,IMAX,QMAX)       :: ccond
     real( dp), dimension(MMAX,IMAX)            :: CCONDNIT


         if (ICONW.eq.2) then

             ! calculate mass transfer coefficient ccond
             call nitcondens(imax,press,temp,dpaw,alphanit,fcoj,CCONDNIT)

             do m=1,mmax
               do i=1,imax
                 do q=1,QMAX
                   ccond(m,i,q) = 0._dp
                 end do
                 ccond(m,i,A_NIT) = ccondnit(m,i)
                 ccond(m,i,A_NH4) = ccondnit(m,i)
               end do
             end do

             !write(6,*) '1hno3 nh3',NVAP(A_NIT),NVAP(A_NH4)
             !write(6,*) '1no3  nh4',MASS(3,6,A_NH4),MASS(3,6,A_NIT)

             call mosaic_control(imax,RH,temp, mass,nvap,mpt,dpaw,n,  &
                      ccond,hyst,henry_eff_hno3,henry_eff_hcl,        &
                      henry_eff_nh3,svmc_min_hno3,svmc_min_hcl,       & 
                      svmc_min_nh3, cioncharge,kh_nh3,                &
                      flag_dissolution,Keq_nh4no3_0,Keq_nh4cl_0)


             !write(6,*) '2hno3 nh3',NVAP(A_NIT),NVAP(A_NH4)
             !write(6,*) '2no3  nh4',MASS(3,6,A_NH4),MASS(3,6,A_NIT)

         else

! default values if MOSAIC solver is not used
             do m=1,mmax
               do i=1,imax
                 henry_eff_hno3(m,i)   = 0._dp
                 henry_eff_hcl(m,i)    = 0._dp
                 henry_eff_nh3(m,i)    = 0._dp
                 svmc_min_hno3(m,i)    = 0._dp
                 svmc_min_hcl(m,i)     = 0._dp
                 svmc_min_nh3(m,i)     = 0._dp
                 cioncharge(m,i)       = 0._dp
                 flag_dissolution(m,i) = 0
               end do
             end do

         endif

  end subroutine interface_mosaic


  end module gde_aerosol_solver
