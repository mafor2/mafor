! <gde_cloud_microphys.f90 - A component of the Multicomponent
!                     Aerosol Dynamics Model MAFOR>
!*****************************************************************************! 
!* 
!*    Copyright (C) 2011-2026 Matthias Steffen Karl
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
!*    References Cloud Solver
!*    -----------------------
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
!*
!*******************************************************************************!
module gde_cloud_microphys

  use messy_mecca_kpp_Parameters
  use messy_mecca_kpp_global, only : APN
  use messy_mecca_kpp_global, only : ind_Clm_a, ind_Nap_a
  use messy_mecca_kpp_global, only : ind_HSO4m_a, ind_SO4mm_a
  use messy_mecca_kpp_global, only : ind_NO3m_a, ind_CH3SO3m_a
  use messy_mecca_kpp_global, only : ind_DMAp_a, ind_NH4p_a
  use messy_mecca_kpp_global, only : ind_IO3m_a
  use messy_mecca_kpp_global, only : ind_HC2O4m_a, ind_C2O4mm_a
  use messy_mecca_kpp_global, only : ind_CH3COCOOm_a
  use messy_mecca_kpp_global, only : ind_C2H5C2O4m_a
  use messy_mecca_kpp_global, only : ind_C2H4C2O4mm_a
  use messy_mecca_kpp_global, only : ind_ADIPAC_a
  use messy_mecca_kpp_global, only : ind_GLUTARAC_a

  use gde_constants,  only         : pi,N_A,k_B
  use gde_constants,  only         : LAM,VIS
  use gde_constants,  only         : CONVM
  use gde_constants,  only         : nucomin, massmin

  use gde_sensitiv,     only       : IDEPO
  use gde_sensitiv,     only       : IWETD
  use gde_sensitiv,     only       : ICOAG
  use gde_sensitiv,     only       : ICOND

  use gde_input_data,   only       : MMAX,AMAX
  use gde_input_data,   only       : NU,NA,AI,AS,CS
  use gde_input_data,   only       : aqmax
  use gde_input_data,   only       : NSOA
  use gde_input_data,   only       : SVI,MSA,NO3,DMA,NH4
  use gde_input_data,   only       : IOD,OXA,SUC,CHL,SSA
  use gde_input_data,   only       : A_SAL

  use gde_toolbox,      only       : molec2ug

  use gde_deposition,   only       : settling
  use gde_coagulation,  only       : coagulation_coeff
  use gde_coagulation,  only       : coagulation_target
  use gde_coagulation,  only       : coagulation


  private

  public  :: cloud_solver

  private :: condensation_incloud
  private :: collection_kernel


contains

subroutine cloud_solver(caq,DTIME,IMAX,press,temp,mbh,DPcrit,    &
                        DPA,DPAW,VPT,ROOPW,M_oc,rp0in,Dfracin,   &
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
    !           rp0in    primary spherule radius        [nm]
    !           Dfracin  fractal dimension              [-]
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
    !      CHL = 9
    !      SAL = 10
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

    implicit none

    ! input
    integer, intent(in)                                :: IMAX
    real( dp), dimension(nspec), intent(in)            :: caq       ! [molec/cm^3]
    real( dp), intent(in)                              :: DTIME     ! [s]
    real( dp), intent(in)                              :: temp      ! [K]
    real( dp), intent(in)                              :: press     ! [Pa]
    real( dp), intent(in)                              :: mbh       ! [m]
    real( dp), intent(in)                              :: DPcrit    ! [m]
    real( dp), intent(in)                              :: rp0in     ! [nm]
    real( dp), intent(in)                              :: Dfracin   ! [-]
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

    real( dp), parameter                               :: dlim=6.0e-8 ! [m]

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
                MCl,         &     ! CHL = 9         
                MNa /)             ! SOD = 10 = not condensate


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

      !A_CHL (9)
           caqdiff(CHL,zkc) = caq(ind_Clm_a(zkc))               - &
                              caqold(CHL,zkc)

      !A_SAL (10)
           caqdiff(SSA,zkc) = caq(ind_Nap_a(zkc))               - &
                              caqold(SSA,zkc)

      ! Convert dCaq to mass change [ng/m^3]
           do L=1,aqmax

             caqdiff(L,zkc) = caqdiff(L,zkc)                    * &
                              molec2ug( mwaq(L) )*1.e3_dp

           enddo

      ! Growth Limit for NO3 (HNO3g) and CHL (HClg)
      ! ! NH4Cl(s) == NH3(g) + HCl(g)
      ! ! NH4NO3(s) == NH3(g) + HNO3(g)
      ! ! NaCl(s) + NO3m == NaNO3(s) + Clm

         !  caqdiff(CHL,zkc) = min( caqdiff(CHL,zkc), 0.3_dp )
           caqdiff(NO3,zkc) = min( caqdiff(NO3,zkc), 0.6_dp )  ! 0.3

      ! Set growth of oxalic acid to zero, too volatile
           caqdiff(OXA,zkc) = 0.0
      ! TEST FOR MSA-BUDGET / without MSA aq. phase condensation
      !     caqdiff(MSA,zkc) = 0.0

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
         call coagulation_coeff(temp,rp0in,Dfracin,ROOPW,DPAW,kcoag,IMAX)

      ! collection kernel
      ! Kcoll = Ecoa*Ecol*|Vinf_b - Vinf_s|*pi*(rb+rs)**2
      !   Kcoll: collection kernel [m^3/s]
      ! includes kcoag(M,O,I,J) for interst. scavenging
         call collection_kernel(IMAX,temp,press,DPAW ,vterm,      &
                               DPcrit,kcoag,   kcoll)

      ! collection equation -> number and mass fluxes
         if (DPcrit.lt.4.e-6_dp) then
      ! droplet-particle collection
           call coagulation(temp,DTIME, DPAW,VPT,N,MASS,          &
                         IMAX,IAG,kcoll, coags,fluxm,flux)
         else
      ! particle coagulation
           call coagulation(temp,DTIME, DPAW,VPT,N,MASS,          &
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

      if ( (ICOND.eq.1) .and. (DPcrit.lt.4.e-6_dp) ) then

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

            if ( DPA(M,I).gt.dlim ) then

          ! mass fluxes of the aqueous phase species (1 to aqmax-1)
          ! APN = MMAX-2
              do L=1, aqmax-1
                MASS(M,I,L)    = MASS(M,I,L) + caqdiff(L,M-2) * wsurface
                MASS(M,I,L)    = MASS(M,I,L) + fluxcm(M,I,L)
                MASS(M,I,L)    = max(MASS(M,I,L),massmin)
              enddo

          ! mass flux of seasalt (Na+)
              MASS(M,I,A_SAL)    = MASS(M,I,A_SAL)               +  &
                                   caqdiff(aqmax,M-2) * wsurface
              MASS(M,I,A_SAL)    = MASS(M,I,A_SAL) + fluxcm(M,I,aqmax)
              MASS(M,I,A_SAL)    = max(MASS(M,I,A_SAL),massmin)

            endif

         ! condensation flux N [1/m^3]
         ! only activated particles shall grow
         ! apply non-negative constraint
            if ( DPA(M,I).gt.dlim ) then
           ! if ( DPA(M,I).gt.DPcrit ) then
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

      ! Reduce settling rate for large particles (SF1)
            if (DPAW(M,I).ge.1.e-5_dp) then
               depo(M,I) = depo(M,I)*0.5_dp
            endif

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


    !----------------------------------------------------------------------
    !----------------  INCLOUD CONDENSATION       -------------------------
    !----------------------------------------------------------------------

  subroutine condensation_incloud( IMAX,DTIME,VPT,ROOP,DPA,    &
                 N,MASS,caqdiff, DPCRIT,                       &
                 FLUXC,FLUXCM )
    !********************************************************************
    !
    !     C  O  N  D  E  N  S  A  T  I  O  N   I N   C L O U D
    !
    !********************************************************************
    !
    ! Condensation in cloud module
    !
    !   Computes condensation flux for aqueous phase modes
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !  INPUT
    !  -----
    !
    !  DPCRIT: critical diameter                  [m]
    !     VPT: particle volume in bin             [m^3]
    !     DPA: particle diameter                  [m^3]
    !    ROOP: wet particle density               [kg/m^3]
    !       N: particle number conc. in bin       [1/m^3]
    !    MASS: component mass conc. in bin        [ng/m^3]
    ! CAQDIFF: concentration change dCaq          [ng/m^3]
    ! 
    !
    !  OUTPUT
    !  ------
    !
    !  FLUXC:  condensation flux N in bin      [1/m^3] 
    ! FLUXCM:  condensation flux m in bin      [ng/m^3]
    !
    !  reference
    !  ---------
    !  S.N Pandis, J.H. Seinfeld, C. Pilinis,
    !      Chemical composition differences in fog and
    !      cloud droplets of different sizes,
    !      Atmospheric Environment, 24A(7), 1957-1969, 1990a
    !
    !********************************************************************   
    
     implicit none

    ! input
    integer, intent(in)                               :: IMAX
    real( dp), intent(in)                             :: DTIME
    real( dp), intent(in)                             :: DPCRIT

    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPA
    real( dp), dimension(MMAX,IMAX),intent(in)        :: VPT
    real( dp), dimension(MMAX,IMAX),intent(in)        :: N
    real( dp), dimension(MMAX,IMAX),intent(in)        :: ROOP
    real( dp),dimension(aqmax,APN),intent(in)         :: caqdiff

    real( dp), dimension(MMAX,IMAX,AMAX-1),intent(in) :: MASS
         
    ! output
    real( dp), dimension(MMAX,IMAX),intent(out)       :: FLUXC
    real( dp), dimension(MMAX,IMAX,aqmax),intent(out) :: FLUXCM

    ! local
    real( dp)                                         :: VPNEW
    real( dp)                                         :: JCOND
    real( dp),dimension(aqmax,APN)                    :: cdiff

    integer     :: I,M
    integer     :: L


       ! initialization
         FLUXC(:,:)=0._dp
         FLUXCM(:,:,:)=0._dp


       ! Calculate change of volume in a size bin
       ! Based on Pandis et al. (1990a) EQ (7)
       ! JCOND = Dp^2*dDp/dt = dV/dt
       ! convert kg to ng with CONVM
       ! AI mode is aqueous phase mode 1

         do M=AI,MMAX             
           do I=1,IMAX

             VPNEW = VPT(M,I)
             JCOND = 0._dp

          ! JCOND has unit [m^6/ng]
          ! Conversion of ROOP from kg/m^3 into ng/m^3
             if (N(M,I).gt.nucomin) then
               JCOND = 2._dp/(pi*ROOP(M,I)*CONVM*N(M,I))
             endif


          ! No growth below critical diameter
             if (DPA(M,I).le.DPCRIT) then
               VPNEW=VPT(M,I)
             else

          ! New volume after condensation
                do L=1,aqmax

          ! Weighting of mass change rate dCaq for each bin
                  !!cdiff(L,M-1) = caqdiff(L,M-1) * (DPAW(M,I)/sumdpcs(M))
                  cdiff(L,M-1) = caqdiff(L,M-1) / IMAX

          ! Condensation of aqueous species L
          ! cdiff: change of mass rate per bin [ng/m^3/s]
          ! zkc=M-1
                  VPNEW = VPNEW + JCOND * cdiff(L,M-1) *DTIME

               enddo
             endif

             ! Lower bound of volume growth
             VPNEW=max( VPT(M,I),VPNEW )
             ! Upper bound of volume growth
             if (I.eq.IMAX) then
               if (M.eq.CS) then
                  VPNEW=min( VPT(M,I)*1.10,VPNEW )
               else
                  VPNEW=min( VPT(M+1,1),VPNEW )
               endif
             else
              VPNEW=min( VPT(M,I+1),VPNEW )
             endif

          !   print *,'cond',M,I,N(M,I),JCOND,cdiff(1,M-1),VPT(M,I),VPNEW


       ! Subroutine condensation_incloud
       ! Compute condensation flux for
       ! aqueous phase modes
       ! Aerosol species
       !      SVI = 1   => A_SUL (1)
       !      MSA = 2   => A_MSA (2)
       !      NO3 = 3   => A_NIT (3)
       !      DMA = 4   => A_AMI (4)
       !      NH4 = 5   => A_NH4 (5)
       !      IOD = 6   => A_IO3 (6)
       !      OXA = 7   => A_OR1 (7)
       !      SUC = 8   => A_OR2 (8)
       !      CHL = 9   => A_CHL (16)
       !      SAL =10   => A_SAL (17)

             if (I.eq.1) then
          ! number fluxes
               FLUXC(M,I)=FLUXC(M,I)-N(M,I)
               FLUXC(M,I)=FLUXC(M,I)+N(M,I)*(VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I))
               FLUXC(M,I+1)=FLUXC(M,I+1)+N(M,I)*(VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I))
          ! mass fluxes of the aqueous phase species (1 to aqmax-1)  
               do L=1,aqmax-1               
                 FLUXCM(M,I,L)=FLUXCM(M,I,L)-MASS(M,I,L)
                 FLUXCM(M,I,L)=FLUXCM(M,I,L)+MASS(M,I,L)*(VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I))
                 FLUXCM(M,I+1,L)=FLUXCM(M,I+1,L)+MASS(M,I,L)*(VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I))
               enddo
          ! mass flux of seasalt
               FLUXCM(M,I,aqmax)=FLUXCM(M,I,aqmax)-MASS(M,I,A_SAL)
               FLUXCM(M,I,aqmax)=FLUXCM(M,I,aqmax)+MASS(M,I,A_SAL)*(VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I))
               FLUXCM(M,I+1,aqmax)=FLUXCM(M,I+1,aqmax)+MASS(M,I,A_SAL)*(VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I))
            endif

             if ((I.gt.1).and.(I.lt.IMAX)) then
          ! number fluxes
               FLUXC(M,I)=FLUXC(M,I)-N(M,I)
               FLUXC(M,I)=FLUXC(M,I)+N(M,I)*(VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I))
               FLUXC(M,I+1)=FLUXC(M,I+1)+N(M,I)*(VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I))
          ! mass fluxes of the aqueous phase species (1 to aqmax-1)  
               do L=1,aqmax-1             
                  FLUXCM(M,I,L)=FLUXCM(M,I,L)-MASS(M,I,L)
                  FLUXCM(M,I,L)=FLUXCM(M,I,L)+MASS(M,I,L)*(VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I))
                  FLUXCM(M,I+1,L)=FLUXCM(M,I+1,L)+MASS(M,I,L)*(VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I))
               enddo
          ! mass flux of seasalt
               FLUXCM(M,I,aqmax)=FLUXCM(M,I,aqmax)-MASS(M,I,A_SAL)
               FLUXCM(M,I,aqmax)=FLUXCM(M,I,aqmax)+MASS(M,I,A_SAL)*(VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I))
               FLUXCM(M,I+1,aqmax)=FLUXCM(M,I+1,aqmax)+MASS(M,I,A_SAL)*(VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I))
             endif


             if (I.eq.IMAX) then
          ! number fluxes
               FLUXC(M,I)=FLUXC(M,I)-N(M,I)
               if (M.lt.CS) then
                 FLUXC(M,I)=FLUXC(M,I)+N(M,I)*(VPT(M+1,1)-VPNEW)/(VPT(M+1,1)-VPT(M,I))
                 FLUXC(M+1,1)=FLUXC(M+1,1)+N(M,I)*(VPNEW-VPT(M,I))/(VPT(M+1,1)-VPT(M,I))
               endif
          ! mass fluxes of the aqueous phase species (1 to aqmax-1)  
               do L=1,aqmax-1
                  FLUXCM(M,I,L)=FLUXCM(M,I,L)-MASS(M,I,L)
                 if (M.lt.CS) then 
                   FLUXCM(M,I,L)=FLUXCM(M,I,L)+MASS(M,I,L)*(VPT(M+1,1)-VPNEW)/(VPT(M+1,1)-VPT(M,I))
                   FLUXCM(M+1,1,L)=FLUXCM(M+1,1,L)+MASS(M,I,L)*(VPNEW-VPT(M,I))/(VPT(M+1,1)-VPT(M,I))
                 endif
               enddo
          ! mass flux of seasalt
               FLUXCM(M,I,aqmax)=FLUXCM(M,I,aqmax)-MASS(M,I,A_SAL)
               if (M.lt.CS) then 
                 FLUXCM(M,I,aqmax)=FLUXCM(M,I,aqmax)+       &
                                   MASS(M,I,A_SAL)*(VPT(M+1,1)-VPNEW)/(VPT(M+1,1)-VPT(M,I))
                 FLUXCM(M+1,1,aqmax)=FLUXCM(M+1,1,aqmax)+   &
                                   MASS(M,I,A_SAL)*(VPNEW-VPT(M,I))/(VPT(M+1,1)-VPT(M,I))
               endif
             endif

        enddo
      enddo


  end subroutine condensation_incloud


  subroutine collection_kernel(IMAX,temp,press,DPA,vterm,      &
                               DPcrit,KCOA,  KCOL )
    !----------------------------------------------------------------------
    !
    !****  COLLECTION KERNEL
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      The collection kernel K describes the interaction
    !      of two colliding particles. The collision partners
    !      are assummed to fall with their terminal velocities
    !      Vinf. The smaller collision partner (s) will be 
    !      collected by the larger collision partner (b) as soon
    !      as it is inside the swept volume of the collector
    !      KCOL = Ecoa*Ecol*|Vinf_B - Vinf_s|*pi*(RB+rs)**2
    !      Ecoa: coalescence efficieny, Table 4, Kirkweg (2003)
    !      Ecol: collision efficiency, Table 2, Kirkweg (2003)
    !
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature                [K]
    !           press    air pressure                   [Pa]
    !           DPA      wet particle diameter          [m]
    !           VPT      particle volume                [m^3]
    !           ROOPW    particle density               [kg/m^3]
    !           vterm    terminal velocity              [m/s]
    !           DPcrit   critical particle diameter     [m]
    !           KCOA     coagulation coefficient        [m^3/s]
    !
    !
    !        output:
    !           KCOL     collection kernel              [m^3/s]
    !
    !
    !      method
    !      ------
    !      Collision/Coalescence as function of size
    !      of the collision partners
    !      For mode AI to CS (excluding NU mode)
    !
    !      reference
    !      ---------
    !      Ventilation coefficient:
    !      J.M. Straka,
    !      Cloud and Precipitation Microphysics,
    !      Principles and Paramterizations,
    !      Cambridge University Press, Cambridge, UK,
    !      p.116-118, 2009.
    !
    !      B.A. Tinsley and L. Zhou
    !      Parameterization of aerosol scavenging
    !      due to atmospheric ionization,
    !      J. Geophys. Res., 120, 8389-8410, 2015.
    !
    !      Collision efficiency:
    !      K. Young,
    !      Microphysical Cloud Processes,
    !      Oxford University Press, New York, USA, 1993.
    !
    !      Coalescence efficiency:
    !      K.V. Beard and H.T. Ochs
    !      Collection and Coalescence Efficiencies for Accretion,
    !      J. Geophys. Res., 89, D5, 7165-7169, 1984.
    !
    !      Collection kernel:
    !      A. Kerkweg, S. Wurzler, T. Reisin, A. Bott,
    !      On the cloud processing of aerosol particles:
    !      An entraining air-parcel model with two-dimensional
    !      spectral cloud microphysics and a new formulation
    !      of the collection kernel,
    !      Q.J.R. Meteorolo. Soc., 129, 1-19, 2003
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

    implicit none

    ! input
    integer, intent(in)                                :: IMAX
    real( dp), intent(in)                              :: temp      ! [K]
    real( dp), intent(in)                              :: press     ! [Pa]
    real( dp), dimension(MMAX,0:(IMAX+1)), intent(in)  :: DPA       ! [m]
    real( dp), dimension(MMAX,IMAX)                    :: vterm     ! [m/s]
    real( dp), intent(in)                              :: DPcrit    ! [m]
    real( dp), dimension(MMAX,MMAX,IMAX,IMAX), intent(in)  :: KCOA  ! [m^3/s]

    ! output
    real( dp), dimension(MMAX,MMAX,IMAX,IMAX), intent(out) :: KCOL  ! [m^3/s]


    ! local
    ! specific gas constant for dry air [J/(kgK)]
    real( dp), parameter                               :: R=287.058
    real( dp),dimension(MMAX,MMAX,IMAX,IMAX)           :: ECOL
    real( dp),dimension(MMAX,MMAX,IMAX,IMAX)           :: ECOA
    real( dp), dimension(MMAX,IMAX)                    :: DIF
    real( dp), dimension(MMAX,IMAX)                    :: VEN

    real( dp)                                          :: rho_air
    real( dp)                                          :: MYY
    real( dp)                                          :: NYY
    real( dp)                                          :: KNA
    real( dp)                                          :: CC
    real( dp)                                          :: DIFFCO
    real( dp)                                          :: NRE
    real( dp)                                          :: NSC
    real( dp)                                          :: NPE
    real( dp)                                          :: VEP

    real( dp)                                          :: rc
    real( dp)                                          :: rlim
    real( dp)                                          :: rs
    real( dp)                                          :: RB
    real( dp)                                          :: difs
    real( dp)                                          :: vterms
    real( dp)                                          :: acoa
    real( dp)                                          :: bcoa
    real( dp)                                          :: xcoa
    real( dp)                                          :: betacoa

    integer                                            :: I,M,J,O
    integer                                            :: MM,K


      ! Initialisation
      KCOL(:,:,:,:)   = 0._dp
      ECOL(:,:,:,:)   = 0._dp
      ECOA(:,:,:,:)   = 0._dp
      DIF(:,:)        = 0._dp
      VEN(:,:)        = 1._dp


      ! Critical radius
      rc      = DPcrit*0.5_dp
      !!! rlim    = 4.00e-7_dp      ! 400 nm
      !!! rlim    = 6.00e-7_dp      ! 600 nm
      rlim    = 8.00e-7_dp      ! 800 nm
      rlim    = MIN(rlim,rc)

      ! Air density
      ! [kg/m^3] 1.2928
      rho_air = press / R / temp      

      ! Dynamic viscosity of air
      ! [kg/m/s]
      MYY= (1.832e-5_dp*(temp**(1.5_dp))*406.4_dp)/            &
          (5093._dp*(temp+110.4_dp))

      ! Kinematic viscosity of air
      ! [m^2/s]
      NYY     = MYY / rho_air


      ! First calculate in all bins
      ! Particle diffusion coefficient
      ! for Brownian Diffusion DIF
      ! and Ventilation coefficient VEN

      do M=1,MMAX
        do I=1,IMAX

          ! Cunnningham slip-flow correction [-]
          ! CC = 1 + Kn_a(A + B*exp(-C/Kn_a))
          ! with A,B,C from Rogak and Flagan 
          ! [J. Coll. Interface Sci., 151, 1992]
          ! LAM and VIS in gde_input_data
          ! Kn_a is Knudsen number of air 
          KNA = 2._dp*LAM / DPA(M,I)
          CC  = 1._dp+KNA*(1.257_dp+0.4_dp*exp(-1.1_dp/KNA))

          ! Particle Diffusion Coefficient   [m^2/s]
          DIF(M,I) = ( k_B * temp * CC / (3.*pi*VIS*DPA(M,I)) )

          ! Vapor diffusivity [m^2/s]
          DIFFCO  = (k_B*temp*CC) /                            &
                    (6._dp*pi*MYY*0.5_dp*DPA(M,I))

          ! Reynold number [-]
          NRE = (vterm(M,I)*DPA(M,I)) / NYY

          ! Schmidt number [-]
          NSC = NYY / DIFFCO

          ! Ventilation coefficent [-]
          ! Straka (2009)
          ! Two conditions depending on NRE and NSC
          ! The transition gives a smooth VEN(Dp) curve
          if (NSC**(1./3.)*sqrt(NRE) .lt. 1.4) then
            VEN(M,I) = 1._dp + 0.108_dp*                       &
                       (NSC**(1./3.)*sqrt(NRE))**2.
          else
            VEN(M,I) = 0.78_dp + 0.308_dp*                     &
                       NSC**(1./3.)*sqrt(NRE)
          endif

       !print *,'DpA CASE DIF',M,I,DPA(M,I),NSC**(1./3.)*sqrt(NRE),VEN(M,I),DIF(M,I)

         enddo       ! I
      enddo         ! M


      ! Second calculate the collection kernel KCOL
      ! KCOL = Ecoa*Ecol*|Vinf_B - Vinf_s|*pi*(RB+rs)**2
      !  Ecoa: coalesence efficieny, Table 4, Kirkweg (2003)
      !  Ecol: collision efficiency, Table 2, Kirkweg (2003)
      ! Accretion is the growth of large partner by 
      ! sweeping the smaller partner.
      ! Large collision partner X(O,J)
      !  with radius RB
      ! Small collision partner X(M,I)
      !  with radius rs
      ! Outer loop: small partner
      ! Inner loop: large partner
      ! Coalescence:
      !   small partner rs: 1.0 - 30 um
      !   large partner RB: >50 um
      !   for small partner rs < 1um -> Ecoa=1
      !   for large partner RB < 50 um -> Ecoa=1

      do M=NU,CS
        do I=1,IMAX

          rs     = DPA(M,I)*0.5_dp
          difs   = DIF(M,I)
          vterms = vterm(M,I)

          do O=AI,CS
            do J=1,IMAX
!      O=CS
!      J=15

              ! Radius of larger collision partner
              RB      = DPA(O,J)*0.5_dp
              ! Initialize parameters of Ecoa
              betacoa = 0._dp
              bcoa    = 0._dp
              acoa    = 0._dp
              xcoa    = 0._dp

              ! Particle ventilation factor [-]
              ! From Tinsley and Zhou (2015), Eq.5 + Eq.6
              ! NPE: cube root of Peclet number
              NPE = ( (2._dp*vterm(O,J)*RB) / difs )**(1./3.)
              VEP = 1._dp + 0.530_dp*NPE * exp(-1.1_dp/NPE)


              IF (RB.GT.rs) THEN

                ! Collision efficiency [-]
                ! Ecol = 4*RB*DIF*VEN /
                !        (rs+RB)**2*|Vinf_B-Vinf_s|

                ECOL(M,O,I,J) = 4._dp*RB*difs*VEN(O,J)       /  &
                      ( (RB+rs)**2.*ABS(vterm(O,J)-vterms) )

                ! Coalescence efficiency [-]

                ! RB >= 50 um
                if (RB.ge.5.e-5) then
                  if ((rs.ge.1.e-6).and.(rs.le.3.e-5)) then

                    betacoa = LOG(rs*1.e6)                   +  &
                          0.44_dp*LOG(RB*1.e6/200._dp)
                    bcoa    = 0.0946_dp * betacoa - 0.319_dp
                    acoa    = SQRT(bcoa**2.0 + 0.00441_dp)
                    xcoa    = (acoa-bcoa)**(1._dp/3._dp)     -  &
                            (acoa+bcoa)**(1._dp/3._dp)
                    ECOA(M,O,I,J) = xcoa + 0.459_dp

                  else if (rs.gt.3.e-5) then  ! rs > 30 um
                    ECOA(M,O,I,J) = 0.6_dp
                  else                        ! rs < 1 um 
                    ECOA(M,O,I,J) = 1._dp
                  endif

                ! RB < 50 um
                else
                  if (rs.gt.3.e-5) then       ! rs > 30 um
                    ECOA(M,O,I,J) = 0.6_dp
                  else 
                    ECOA(M,O,I,J) = 1._dp
                  endif
                endif

              ! Collection kernel [m^3/s]

                KCOL(M,O,I,J) = ECOA(M,O,I,J) * ECOL(M,O,I,J)  *  &
                        ABS(vterm(O,J)-vterms) * pi * (RB+rs)**2.

              ! Collection kernel for interstitial scavenging
              ! if critcial radius is below 2 um
              ! upper radius limit rlim (= 800 nm)
                if ((rs.le.rlim).and.(rc.lt.2.e-6_dp)) then

              ! MSK 2025-11-07
              ! Pruppacher and Klett (1997): 
              ! Collision kernel for ventilated Brownian diffusion
              ! KCOL = 4pi * D_particle * VEN_particle * RB
              ! Particle diffusion coefficient of smaller particle
              ! This is the same as KCOL above when
              ! ECOL as above
              ! ECOA = 1 (stick always)

                  KCOL(M,O,I,J) = 4._dp*pi * difs * VEP * RB

              ! Size-dependent linear scaling function
                  if (rs.gt.2.0e-8_dp) then
                     KCOL(M,O,I,J) = KCOL(M,O,I,J)* (2.50e8*rs -4.00_dp)
                  endif


                endif

              ! No collection until critcial radius is below 1 um
                if (rc.ge.2.e-6_dp) then
                  KCOL(M,O,I,J) = 0._dp
                endif

              ! RB < rs
              ELSE
                ECOL(M,O,I,J) = 0._dp
                ECOA(M,O,I,J) = 0._dp
                KCOL(M,O,I,J) = 0._dp
              ENDIF

             ! if ((M==AI).and.(I==10).and.(rc.lt.2.e-6_dp)) then
             ! write(6,*) 'Kcol',O,J,rs/RB,ECOL(M,O,I,J),KCOA(M,O,I,J),KCOL(M,O,I,J)
             ! endif

           ! Print KCOL and KCOA for one collector size
           ! comment O-J loop to print KCOL
           ! comment rs-scaling
           ! and replace DPAW by DPA in the call to collection_kernel
           !  IF (RB.GT.rs) THEN
           !    write(6,'(I2,A1,I2,6ES12.4)') M,' ',I,RB,rs,KCOL(M,O,I,J),KCOA(M,O,I,J)
           !  ENDIF

            enddo   ! J
          enddo     ! O

        enddo       ! I
      enddo         ! M          


  end subroutine collection_kernel


  end module gde_cloud_microphys
