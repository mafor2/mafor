! <gde_aerosol_props.f90 - A component of the Multicomponent
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
module gde_aerosol_props

   use gde_sensitiv,   only  : ICONA,ICOAG

   use gde_input_data, only  : MMAX, AMAX, iamax
   use gde_input_data, only  : NSOA
   use gde_input_data, only  : NU,NA,AI,AS,CS
   use gde_input_data, only  : A_SUL,A_MSA,A_IO3
   use gde_input_data, only  : A_NIT,A_NH4,A_AMI
   use gde_input_data, only  : A_OR1,A_OR2,A_OR3
   use gde_input_data, only  : A_OR4,A_OR5,A_OR6
   use gde_input_data, only  : A_OR7,A_OR8,A_OR9
   use gde_input_data, only  : A_CHL
   use gde_input_data, only  : A_SAL,A_XXX,A_EBC,A_DUS
   use gde_input_data, only  : A_WAT

   use gde_input_data, only  : SU,OC,AM,NI,MS,SA,XX,EC,DU

   use gde_input_data, only  : CONVM
   use gde_input_data, only  : DENSA,DENXX,DENDU,DENEC
   use gde_input_data, only  : DENAM,DENNI
   use gde_input_data, only  : massmin

   use gde_constants,  only  : pi,RHOH2O

   use gde_init_aero,  only  : GMD,SIG
   use gde_init_aero,  only  : BGMD,BGSIG
   use gde_init_aero,  only  : EGMD,ESIG

   use gde_init_gas,   only  : rp0,Dfrac
   use gde_init_gas,   only  : gamma_oc1_m,gamma_oc2_m,gamma_oc3_m
   use gde_init_gas,   only  : gamma_oc4_m,gamma_oc5_m,gamma_oc6_m
   use gde_init_gas,   only  : gamma_oc7_m,gamma_oc8_m,gamma_oc9_m

   use gde_toolbox,    only  : roolm


implicit none

   INTRINSIC :: SELECTED_REAL_KIND

   public :: initsizedistribution
   public :: initNumberMass
   public :: initBGNumberMass
   public :: initEMNumberMass
   public :: getdensity
   public :: gettotalmass


! KPP DP - Double precision kind
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)



contains


  subroutine initSizeDistribution(IMAX,DPMIN,DPMAX,DPA,VPT,DLOGDP,DLINDP)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      initialize the lognormal size distribution
    !
    !      interface
    !      ---------
    !
    !        input:
    !           DPMIN   smallest diameter               [m]
    !           DPMAX   largest diameter                [m]
    !
    !        output:
    !           DPA     dry particle diameter           [m]
    !           VPT     dry particle volume             [m^3]
    !           DLOGDP  logarithmic width of bin        [-]
    !           DLINDP  linear width of bin             [m]
    !
    !      method
    !      ------
    !      initialize the lognormal size distribution
    !
    !      reference
    !      ---------
    !      none
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

  implicit none

     integer, intent(in)                                :: IMAX
     real( dp), intent(out)                             :: DPMIN
     real( dp), intent(out)                             :: DPMAX

     real( dp), dimension(MMAX,0:(IMAX+1)),intent(out)  :: DPA
     real( dp), DIMENSION(MMAX,IMAX), intent(out)       :: VPT
     real( dp), DIMENSION(MMAX,IMAX), intent(out)       :: DLOGDP
     real( dp), DIMENSION(MMAX,IMAX), intent(out)       :: DLINDP

! local
     real( dp), dimension(MMAX)        :: VRAT
     real( dp)                         :: QDPM


     integer                           :: cm
     integer                           :: M,I

!!! initialise bin size distribution with one lognormal mode
!   DPA is dry diameter [m]

        QDPM=EXP(LOG(DPMAX/DPMIN)/(DBLE(IMAX*MMAX)-1.)) 

! each mode to get IMAX bins
        cm=0
        do M=NU,CS
          do I=1, IMAX
! Bin diameter
            DPA(M,I)=DPMIN*(QDPM**(DBLE(cm*IMAX+I)-1._dp))
! Bin volume (dry)
            VPT(M,I)=(pi/6._dp)*DPA(M,I)**3._dp
          !  print *,'Mode,bin Dp Vp ',M,I,DPA(M,I),VPT(M,I)
         end do
         cm=cm+1
        end do

! first and last bin, volume ratio
        do M=NU,CS
          DPA(M,0)=DPA(M,1)/(DPA(M,2)/DPA(M,1))
          DPA(M,IMAX+1)=DPA(M,IMAX)*(DPA(M,IMAX)/DPA(M,IMAX-1))
          VRAT(M)=(DPA(M,IMAX)/DPA(M,1))**(3./(IMAX-1))
        end do

! all other bins
        do M=NU,CS
          do I=1,IMAX

! Logarithmic width of bin
            DLOGDP(M,I)=0.5_dp*LOG(DPA(M,I+1)/DPA(M,I-1))
! Linear width of bin
           !DLINDP(M,I)=0.5_dp*(DPA(M,I+1)-DPA(M,I-1))
            DLINDP(M,I)= DPA(M,I)*(2._dp**(1/3._dp))       *  &
                        ((VRAT(M))**(1/3._dp)-1._dp)       /  &
                        ((1._dp+VRAT(M))**(1/3._dp))                 

          end do
        end do

  end subroutine initSizeDistribution


!------------------------------------------------------------------

  subroutine initNumberMass(IMAX,DPA,DLINDP,VPT,DEN,                &
                            MSULFTOT, MMSAPTOT, MNITRTOT, MAMMOTOT, &
                            MSALTTOT, MECBCTOT, MDUSTTOT, MXXXXTOT, & 
                            MORGCTOT,   MASS, N)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      initialize number and mass of aerosol compounds in size bins
    !
    !      interface
    !      ---------
    !
    !        input:
    !           DPA     dry particle diameter           [m]
    !           MASS    mass concentration              [ng/m^3]
    !           DLINDP  linear width of bin             [m]
    !           VPT     dry particle volume             [m^3]
    !           DEN     density of compound             [kg/m^3]
    !
    !        output:
    !           MASS    mass concentration              [ng/m^3]
    !           N       number concentration            [part/m^3]
    !
    !      method
    !      ------
    !      Calculate mass conc of compoound in each bin I
    !      from log-normal mass distribution
    !      GMD is geometric-mean mass diameter!!!
    !      we assume average particle density (DEN) to be
    !      constant in all bins.
    !
    !      reference
    !      ---------
    !      none
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

  implicit none

     integer, intent(in)                                 :: IMAX
     real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)    :: DPA
     real( dp), DIMENSION(MMAX,IMAX), intent(in)         :: DLINDP
     real( dp), DIMENSION(MMAX,IMAX), intent(in)         :: VPT
     real( dp), dimension(iamax), intent(in)             :: DEN

     real( dp), dimension(MMAX),intent(in)               :: MSULFTOT
     real( dp), dimension(MMAX),intent(in)               :: MSALTTOT
     real( dp), dimension(MMAX),intent(in)               :: MAMMOTOT
     real( dp), dimension(MMAX),intent(in)               :: MNITRTOT
     real( dp), dimension(MMAX),intent(in)               :: MMSAPTOT
     real( dp), dimension(MMAX),intent(in)               :: MXXXXTOT
     real( dp), dimension(MMAX),intent(in)               :: MECBCTOT
     real( dp), dimension(MMAX),intent(in)               :: MDUSTTOT
     real( dp), dimension(MMAX),intent(in)               :: MORGCTOT


!in/out
     real( dp), dimension(MMAX,IMAX,AMAX),intent(in out) :: MASS
     real( dp), DIMENSION(MMAX,IMAX), intent(in out)     :: N

! local
     REAL( dp),allocatable,dimension(:,:,:)              :: VCONC
     REAL( dp),allocatable,dimension(:,:,:)              :: LINDISM

     real( dp)        :: deffsoot
     real( dp)        :: Dp0
     integer          :: M,I


! Allocate aerosol terms
        if (.not. allocated(VCONC))        ALLOCATE(VCONC(MMAX,IMAX,iamax))
        if (.not. allocated(LINDISM))      ALLOCATE(LINDISM(MMAX,IMAX,iamax))



        do M=1,MMAX
          do I=1,IMAX

            LINDISM(M,I,SU)=(MSULFTOT(M)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M)))) &
                           *EXP(-0.5_dp*(LOG(DPA(M,I)/GMD(M))/LOG(SIG(M)))**2._dp)
            LINDISM(M,I,OC)=(MORGCTOT(M)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M)))) &
                           *EXP(-0.5_dp*(LOG(DPA(M,I)/GMD(M))/LOG(SIG(M)))**2._dp)
            LINDISM(M,I,AM)=(MAMMOTOT(M)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M)))) &
                           *EXP(-0.5_dp*(LOG(DPA(M,I)/GMD(M))/LOG(SIG(M)))**2._dp)
            LINDISM(M,I,NI)=(MNITRTOT(M)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M)))) &
                           *EXP(-0.5_dp*(LOG(DPA(M,I)/GMD(M))/LOG(SIG(M)))**2._dp)
            LINDISM(M,I,MS)=(MMSAPTOT(M)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M)))) &
                           *EXP(-0.5_dp*(LOG(DPA(M,I)/GMD(M))/LOG(SIG(M)))**2._dp)
            LINDISM(M,I,SA)=(MSALTTOT(M)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M)))) &
                           *EXP(-0.5_dp*(LOG(DPA(M,I)/GMD(M))/LOG(SIG(M)))**2._dp)
            LINDISM(M,I,XX)=(MXXXXTOT(M)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M)))) &
                           *EXP(-0.5_dp*(LOG(DPA(M,I)/GMD(M))/LOG(SIG(M)))**2._dp)
            LINDISM(M,I,EC)=(MECBCTOT(M)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M)))) &
                           *EXP(-0.5_dp*(LOG(DPA(M,I)/GMD(M))/LOG(SIG(M)))**2._dp)
            LINDISM(M,I,DU)=(MDUSTTOT(M)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M)))) &
                           *EXP(-0.5_dp*(LOG(DPA(M,I)/GMD(M))/LOG(SIG(M)))**2._dp)


! Add contribution of next lower mode
            if ((M.gt.1).and.(M.lt.MMAX)) then           
              LINDISM(M,I,SU)=LINDISM(M,I,SU) +                                 &
                      (MSULFTOT(M-1)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M-1))))   &
                      *EXP(-0.5_dp*(LOG(DPA(M,I)/(0.3*GMD(M-1)))/LOG(SIG(M-1)))**2._dp)
              LINDISM(M,I,OC)=LINDISM(M,I,OC) +                                 &
                      (MORGCTOT(M-1)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M-1))))   &
                      *EXP(-0.5_dp*(LOG(DPA(M,I)/(0.3*GMD(M-1)))/LOG(SIG(M-1)))**2._dp)
              LINDISM(M,I,AM)=LINDISM(M,I,AM) +                                 &
                      (MAMMOTOT(M-1)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M-1))))   &
                      *EXP(-0.5_dp*(LOG(DPA(M,I)/(0.3*GMD(M-1)))/LOG(SIG(M-1)))**2._dp)
              LINDISM(M,I,NI)=LINDISM(M,I,NI)  +                                &
                      (MNITRTOT(M-1)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M-1))))   &
                      *EXP(-0.5_dp*(LOG(DPA(M,I)/(0.3*GMD(M-1)))/LOG(SIG(M-1)))**2._dp)
              LINDISM(M,I,MS)=LINDISM(M,I,MS)  +                                &
                      (MMSAPTOT(M-1)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M-1))))   &
                      *EXP(-0.5_dp*(LOG(DPA(M,I)/(0.3*GMD(M-1)))/LOG(SIG(M-1)))**2._dp)
              LINDISM(M,I,SA)=LINDISM(M,I,SA)  +                                &
                      (MSALTTOT(M-1)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M-1))))   &
                      *EXP(-0.5_dp*(LOG(DPA(M,I)/(0.3*GMD(M-1)))/LOG(SIG(M-1)))**2._dp)
              LINDISM(M,I,XX)=LINDISM(M,I,XX)  +                                &
                      (MXXXXTOT(M-1)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M-1))))   &
                      *EXP(-0.5_dp*(LOG(DPA(M,I)/(0.3*GMD(M-1)))/LOG(SIG(M-1)))**2._dp)
              LINDISM(M,I,EC)=LINDISM(M,I,EC)  +                                &
                      (MECBCTOT(M-1)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M-1))))   &
                      *EXP(-0.5_dp*(LOG(DPA(M,I)/(0.3*GMD(M-1)))/LOG(SIG(M-1)))**2._dp)
              LINDISM(M,I,DU)=LINDISM(M,I,DU)  +                                &
                      (MDUSTTOT(M-1)/(SQRT(2._dp*pi)*DPA(M,I)*LOG(SIG(M-1))))   &
                      *EXP(-0.5_dp*(LOG(DPA(M,I)/(0.3*GMD(M-1)))/LOG(SIG(M-1)))**2._dp)
            endif


! Mass concentrations in ng/m3
            MASS(M,I,A_SUL)=LINDISM(M,I,SU)*DLINDP(M,I)

! Divide ORGC between SOA components
            MASS(M,I,A_OR1)=gamma_oc1_m(M)*LINDISM(M,I,OC)*DLINDP(M,I)
            MASS(M,I,A_OR2)=gamma_oc2_m(M)*LINDISM(M,I,OC)*DLINDP(M,I)
            MASS(M,I,A_OR3)=gamma_oc3_m(M)*LINDISM(M,I,OC)*DLINDP(M,I)
            MASS(M,I,A_OR4)=gamma_oc4_m(M)*LINDISM(M,I,OC)*DLINDP(M,I)
            MASS(M,I,A_OR5)=gamma_oc5_m(M)*LINDISM(M,I,OC)*DLINDP(M,I)
            MASS(M,I,A_OR6)=gamma_oc6_m(M)*LINDISM(M,I,OC)*DLINDP(M,I)
            MASS(M,I,A_OR7)=gamma_oc7_m(M)*LINDISM(M,I,OC)*DLINDP(M,I)
            MASS(M,I,A_OR8)=gamma_oc8_m(M)*LINDISM(M,I,OC)*DLINDP(M,I)
            MASS(M,I,A_OR9)=gamma_oc9_m(M)*LINDISM(M,I,OC)*DLINDP(M,I)

! Divide between ammonium and amminium
            if(icona.eq.1) then
              MASS(M,I,A_AMI)=LINDISM(M,I,AM)*DLINDP(M,I)
              MASS(M,I,A_NH4)=0.0
            else
              MASS(M,I,A_AMI)=0.0
              MASS(M,I,A_NH4)=LINDISM(M,I,AM)*DLINDP(M,I)
            endif

            MASS(M,I,A_NIT)=LINDISM(M,I,NI)*DLINDP(M,I)
            MASS(M,I,A_MSA)=LINDISM(M,I,MS)*DLINDP(M,I)

! Divide SALT in chlorine and NaCl (=Na)
            MASS(M,I,A_CHL)=0.54*LINDISM(M,I,SA)*DLINDP(M,I)    ! chlorine
            MASS(M,I,A_SAL)=0.46*LINDISM(M,I,SA)*DLINDP(M,I)    ! sodium

            MASS(M,I,A_XXX)=LINDISM(M,I,XX)*DLINDP(M,I)
            MASS(M,I,A_EBC)=LINDISM(M,I,EC)*DLINDP(M,I) 
            MASS(M,I,A_DUS)=LINDISM(M,I,DU)*DLINDP(M,I)


! Volume concentration in m3/m3
            VCONC(M,I,SU)=MASS(M,I,A_SUL)/(CONVM*DEN(SU))
            VCONC(M,I,OC)=( MASS(M,I,A_OR1)+MASS(M,I,A_OR2)+MASS(M,I,A_OR3)  &
                          + MASS(M,I,A_OR4)+MASS(M,I,A_OR5)+MASS(M,I,A_OR6)  &
                          + MASS(M,I,A_OR7)+MASS(M,I,A_OR8)+MASS(M,I,A_OR9) ) &
                          / (CONVM*DEN(OC))
            if(icona.eq.1) then
              VCONC(M,I,AM)=MASS(M,I,A_AMI)/(CONVM*DEN(AM))
            else
              VCONC(M,I,AM)=MASS(M,I,A_NH4)/(CONVM*DEN(AM))
            endif
            VCONC(M,I,NI)=MASS(M,I,A_NIT)/(CONVM*DEN(NI))
            VCONC(M,I,MS)=MASS(M,I,A_MSA)/(CONVM*DEN(MS))
            VCONC(M,I,SA)=( MASS(M,I,A_SAL) + MASS(M,I,A_CHL) )  &
                          /(CONVM*DEN(SA))
            VCONC(M,I,XX)=MASS(M,I,A_XXX)/(CONVM*DEN(XX))
            VCONC(M,I,DU)=MASS(M,I,A_DUS)/(CONVM*DEN(DU))

! COAG=2 Correction for fractal geometry of soot
!           rp0     primary spherule radius         [nm]
!           Dfrac   fractal dimension               [-]
            Dp0=2._dp*rp0
            if ( (ICOAG.eq.2).or.(ICOAG.eq.5) ) then
              deffsoot = min( (DEN(EC)*(DPA(M,I)*1.e9/Dp0)**(Dfrac-2.7)), DEN(EC))    
              deffsoot = max(deffsoot,1200.0_dp)
            else
              deffsoot = DEN(EC)
            endif

            VCONC(M,I,EC)=MASS(M,I,A_EBC)/(CONVM*deffsoot)


! Calculated particle number concentration per bin
! consistent GMD (mass-based)
            N(M,I) = ( VCONC(M,I,SU) + VCONC(M,I,OC) + VCONC(M,I,AM)  +  &
                       VCONC(M,I,NI) + VCONC(M,I,MS) + VCONC(M,I,SA)  +  &
                       VCONC(M,I,XX) + VCONC(M,I,EC) + VCONC(M,I,DU)) / VPT(M,I) 

            !print *,'init ',M,I,VPT(M,I),N(M,I)

           end do
        end do

! Deallocate aerosol terms
        deallocate(VCONC)
        deallocate(LINDISM)


  end subroutine initNumberMass

!------------------------------------------------------------------

  subroutine initBGNumberMass(IMAX,DPA,DLINDP,VPT,DEN,                &
                              BGMCTOT, BGMASS, BGN)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      initialize number and mass of aerosol compounds in size bins
    !      for the background air aerosol
    !
    !      interface
    !      ---------
    !
    !        input:
    !           DPA     dry particle diameter           [m]
    !           BGMCTOT mass concentration per mode     [ng/m^3]
    !           DLINDP  linear width of bin             [m]
    !           VPT     dry particle volume             [m^3]
    !           DEN     density of compound             [kg/m^3]
    !
    !        output:
    !           BGMASS  mass concentration in BG        [ng/m^3]
    !           BGN     number concentration in BG      [part/m^3]
    !
    !      method
    !      ------
    !      Calculate mass conc of compoound in each bin I
    !      from log-normal mass distribution
    !      GMD is geometric-mean mass diameter!!!
    !      we assume average particle density (DEN) to be
    !      constant in all bins.
    !
    !      reference
    !      ---------
    !      none
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

  implicit none

     integer, intent(in)                                 :: IMAX
     real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)    :: DPA
     real( dp), DIMENSION(MMAX,IMAX), intent(in)         :: DLINDP
     real( dp), DIMENSION(MMAX,IMAX), intent(in)         :: VPT
     real( dp), dimension(iamax), intent(in)             :: DEN

     real( dp), dimension(MMAX,iamax), intent(in)        :: BGMCTOT


!in/out
     real( dp), dimension(MMAX,IMAX,AMAX),intent(in out) :: BGMASS
     real( dp), DIMENSION(MMAX,IMAX), intent(in out)     :: BGN

! local
     REAL( dp),allocatable,dimension(:,:,:)              :: BGVCONC
     REAL( dp),allocatable,dimension(:,:,:)              :: BGLINDISM

     integer          :: M,I,K


! Allocate aerosol terms
        if (.not. allocated(BGVCONC))      ALLOCATE(BGVCONC(MMAX,IMAX,iamax))
        if (.not. allocated(BGLINDISM))    ALLOCATE(BGLINDISM(MMAX,IMAX,iamax))



        do M=1,MMAX
          do I=1,IMAX

! calculate linear bin width for background aerosol (except for water)
            do K=1,iamax-1
              BGLINDISM(M,I,K)=( BGMCTOT(M,K) /(SQRT(2.*pi)*DPA(M,I)*LOG(BGSIG(M)))) &
                               *EXP(-0.5*(LOG(DPA(M,I)/BGMD(M))/LOG(BGSIG(M)))**2.)
            enddo

! Add contribution of next lower mode
            if ((M.gt.1).and.(M.lt.MMAX)) then
              do K=1,iamax-1
                BGLINDISM(M,I,K)=BGLINDISM(M,I,K)  +                           & 
                     ( BGMCTOT(M-1,K) /(SQRT(2.*pi)*DPA(M,I)*LOG(BGSIG(M-1)))) &
                     *EXP(-0.5*(LOG(DPA(M,I)/(0.3*BGMD(M-1)))/LOG(BGSIG(M-1)))**2.)
              enddo
            endif


! Mass concentrations in ng/m3
            BGMASS(M,I,A_SUL)=BGLINDISM(M,I,SU)*DLINDP(M,I)
            BGMASS(M,I,A_OR1)=gamma_oc1_m(M)*BGLINDISM(M,I,OC)*DLINDP(M,I)
            BGMASS(M,I,A_OR2)=gamma_oc2_m(M)*BGLINDISM(M,I,OC)*DLINDP(M,I)
            BGMASS(M,I,A_OR3)=gamma_oc3_m(M)*BGLINDISM(M,I,OC)*DLINDP(M,I)
            BGMASS(M,I,A_OR4)=gamma_oc4_m(M)*BGLINDISM(M,I,OC)*DLINDP(M,I)
            BGMASS(M,I,A_OR5)=gamma_oc5_m(M)*BGLINDISM(M,I,OC)*DLINDP(M,I)
            BGMASS(M,I,A_OR6)=gamma_oc6_m(M)*BGLINDISM(M,I,OC)*DLINDP(M,I)
            BGMASS(M,I,A_OR7)=gamma_oc7_m(M)*BGLINDISM(M,I,OC)*DLINDP(M,I)
            BGMASS(M,I,A_OR8)=gamma_oc8_m(M)*BGLINDISM(M,I,OC)*DLINDP(M,I)
            BGMASS(M,I,A_OR9)=gamma_oc9_m(M)*BGLINDISM(M,I,OC)*DLINDP(M,I)


! Divide between ammonium and amminium
            if(icona.eq.1) then
              BGMASS(M,I,A_AMI)=BGLINDISM(M,I,AM)*DLINDP(M,I)
              BGMASS(M,I,A_NH4)=0.0
            else
              BGMASS(M,I,A_AMI)=0.0
              BGMASS(M,I,A_NH4)=BGLINDISM(M,I,AM)*DLINDP(M,I)
            endif
            BGMASS(M,I,A_NIT)=BGLINDISM(M,I,NI)*DLINDP(M,I)
            BGMASS(M,I,A_MSA)=BGLINDISM(M,I,MS)*DLINDP(M,I)
! Divide SALT in chlorine and NaCl (=Na)
            BGMASS(M,I,A_CHL)=0.54*BGLINDISM(M,I,SA)*DLINDP(M,I)    ! chlorine
            BGMASS(M,I,A_SAL)=0.46*BGLINDISM(M,I,SA)*DLINDP(M,I)    ! sodium
            BGMASS(M,I,A_XXX)=BGLINDISM(M,I,XX)*DLINDP(M,I)
            BGMASS(M,I,A_EBC)=BGLINDISM(M,I,EC)*DLINDP(M,I) 
            BGMASS(M,I,A_DUS)=BGLINDISM(M,I,DU)*DLINDP(M,I)


! Volume concentration in m3/m3
            BGVCONC(M,I,SU)=BGMASS(M,I,A_SUL)/(CONVM*DEN(SU))
            BGVCONC(M,I,OC)=(BGMASS(M,I,A_OR1)+BGMASS(M,I,A_OR2)+  &
                             BGMASS(M,I,A_OR3)+BGMASS(M,I,A_OR4)+  &
                             BGMASS(M,I,A_OR5)+BGMASS(M,I,A_OR6)+  &
                             BGMASS(M,I,A_OR7)+BGMASS(M,I,A_OR8)+  &
                             BGMASS(M,I,A_OR9) ) &
                             /(CONVM*DEN(OC))

            if(icona.eq.1) then
              BGVCONC(M,I,AM)=BGMASS(M,I,A_AMI)/(CONVM*DEN(AM))
            else
              BGVCONC(M,I,AM)=BGMASS(M,I,A_NH4)/(CONVM*DEN(AM))
            endif
            BGVCONC(M,I,NI)=BGMASS(M,I,A_NIT)/(CONVM*DEN(NI))
            BGVCONC(M,I,MS)=BGMASS(M,I,A_MSA)/(CONVM*DEN(MS))
            BGVCONC(M,I,SA)=( BGMASS(M,I,A_SAL) + BGMASS(M,I,A_CHL) )  &
                          /(CONVM*DEN(SA))
            BGVCONC(M,I,XX)=BGMASS(M,I,A_XXX)/(CONVM*DEN(XX))
            BGVCONC(M,I,EC)=BGMASS(M,I,A_EBC)/(CONVM*DEN(EC))     
            BGVCONC(M,I,DU)=BGMASS(M,I,A_DUS)/(CONVM*DEN(DU)) 


! Calculated particle number concentration per bin
! consistent GMD (mass-based)
            BGN(M,I)=(BGVCONC(M,I,SU)+BGVCONC(M,I,OC)+BGVCONC(M,I,AM)  +  &
                      BGVCONC(M,I,NI)+BGVCONC(M,I,MS)+BGVCONC(M,I,SA)  +  &
                      BGVCONC(M,I,XX)+BGVCONC(M,I,EC)+BGVCONC(M,I,DU) )/VPT(M,I)                    

           end do
        end do


! Deallocate aerosol terms
        deallocate(BGVCONC)
        deallocate(BGLINDISM)


  end subroutine initBGNumberMass

!------------------------------------------------------------------

  subroutine initEMNumberMass(IMAX,DPA,DLINDP,VPT,DEN,                &
                              EMMCTOT, EMASS, EN)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      initialize number and mass of aerosol compounds in size bins
    !      for the emitted aerosol particles
    !
    !      interface
    !      ---------
    !
    !        input:
    !           DPA     dry particle diameter           [m]
    !           EMMCTOT mass concentration per mode     [ng/m^3]
    !           DLINDP  linear width of bin             [m]
    !           VPT     dry particle volume             [m^3]
    !           DEN     density of compound             [kg/m^3]
    !
    !        output:
    !           EMASS   mass concentration in BG        [ng/m^3]
    !           EN      number concentration in BG      [part/m^3]
    !
    !      method
    !      ------
    !      Calculate mass conc of compoound in each bin I
    !      from log-normal mass distribution
    !      GMD is geometric-mean mass diameter!!!
    !      we assume average particle density (DEN) to be
    !      constant in all bins.
    !
    !      reference
    !      ---------
    !      none
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

  implicit none

     integer, intent(in)                                 :: IMAX
     real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)    :: DPA
     real( dp), DIMENSION(MMAX,IMAX), intent(in)         :: DLINDP
     real( dp), DIMENSION(MMAX,IMAX), intent(in)         :: VPT
     real( dp), dimension(iamax), intent(in)             :: DEN

     real( dp), dimension(MMAX,iamax), intent(in)        :: EMMCTOT


!in/out
     real( dp), dimension(MMAX,IMAX,AMAX),intent(in out) :: EMASS
     real( dp), DIMENSION(MMAX,IMAX), intent(in out)     :: EN

! local
     REAL( dp),allocatable,dimension(:,:,:)              :: EVCONC
     REAL( dp),allocatable,dimension(:,:,:)              :: ELINDISM

     integer          :: M,I,K


! Allocate aerosol terms
        if (.not. allocated(EVCONC))       ALLOCATE(EVCONC(MMAX,IMAX,iamax))
        if (.not. allocated(ELINDISM))     ALLOCATE(ELINDISM(MMAX,IMAX,iamax))


! calculate linear bin width for emitted aerosol (except for water)

        do M=1,MMAX
          do I=1,IMAX
            do K=1,iamax-1
              ELINDISM(M,I,K)=( EMMCTOT(M,K) /(SQRT(2.*pi)*DPA(M,I)*LOG(ESIG(M)))) &
                               *EXP(-0.5*(LOG(DPA(M,I)/EGMD(M))/LOG(ESIG(M)))**2.)
            enddo


! Emitted mass in ng m^-2 s^-1
            EMASS(M,I,A_SUL)=ELINDISM(M,I,SU)*DLINDP(M,I)    
            EMASS(M,I,A_OR1)=gamma_oc1_m(M)*ELINDISM(M,I,OC)*DLINDP(M,I)
            EMASS(M,I,A_OR2)=gamma_oc2_m(M)*ELINDISM(M,I,OC)*DLINDP(M,I)
            EMASS(M,I,A_OR3)=gamma_oc3_m(M)*ELINDISM(M,I,OC)*DLINDP(M,I)
            EMASS(M,I,A_OR4)=gamma_oc4_m(M)*ELINDISM(M,I,OC)*DLINDP(M,I)
            EMASS(M,I,A_OR5)=gamma_oc5_m(M)*ELINDISM(M,I,OC)*DLINDP(M,I)
            EMASS(M,I,A_OR6)=gamma_oc6_m(M)*ELINDISM(M,I,OC)*DLINDP(M,I)
            EMASS(M,I,A_OR7)=gamma_oc7_m(M)*ELINDISM(M,I,OC)*DLINDP(M,I)
            EMASS(M,I,A_OR8)=gamma_oc8_m(M)*ELINDISM(M,I,OC)*DLINDP(M,I)
            EMASS(M,I,A_OR9)=gamma_oc9_m(M)*ELINDISM(M,I,OC)*DLINDP(M,I)

! Divide between ammonium and amminium
            if(icona.eq.1) then
              EMASS(M,I,A_AMI)=ELINDISM(M,I,AM)*DLINDP(M,I)
              EMASS(M,I,A_NH4)=0.0
            else
              EMASS(M,I,A_AMI)=0.0
              EMASS(M,I,A_NH4)=ELINDISM(M,I,AM)*DLINDP(M,I)
            endif
            EMASS(M,I,A_NIT)=ELINDISM(M,I,NI)*DLINDP(M,I)
            EMASS(M,I,A_MSA)=ELINDISM(M,I,MS)*DLINDP(M,I)     
! Divide SALT in chlorine and NaCl (=Na)
            EMASS(M,I,A_CHL)=0.54*ELINDISM(M,I,SA)*DLINDP(M,I)    ! chlorine
            EMASS(M,I,A_SAL)=0.46*ELINDISM(M,I,SA)*DLINDP(M,I)    ! sodium
            EMASS(M,I,A_XXX)=ELINDISM(M,I,XX)*DLINDP(M,I)    
            EMASS(M,I,A_EBC)=ELINDISM(M,I,EC)*DLINDP(M,I)
            EMASS(M,I,A_DUS)=ELINDISM(M,I,DU)*DLINDP(M,I) 


! Volume concentration in m3/m3
            EVCONC(M,I,SU)=EMASS(M,I,A_SUL)/(CONVM*DEN(SU))
            EVCONC(M,I,OC)=(EMASS(M,I,A_OR1)+EMASS(M,I,A_OR2)+    &
                            EMASS(M,I,A_OR3)+EMASS(M,I,A_OR4)+    &
                            EMASS(M,I,A_OR5)+EMASS(M,I,A_OR6)+    &
                            EMASS(M,I,A_OR7)+EMASS(M,I,A_OR8)+    &
                            EMASS(M,I,A_OR9)  ) &
                            /(CONVM*DEN(OC))

            if(icona.eq.1) then
              EVCONC(M,I,AM)=EMASS(M,I,A_AMI)/(CONVM*DEN(AM))
            else
              EVCONC(M,I,AM)=EMASS(M,I,A_NH4)/(CONVM*DEN(AM))
            endif     
            EVCONC(M,I,NI)=EMASS(M,I,A_NIT)/(CONVM*DEN(NI))
            EVCONC(M,I,MS)=EMASS(M,I,A_MSA)/(CONVM*DEN(MS))
            EVCONC(M,I,SA)=( EMASS(M,I,A_SAL) + EMASS(M,I,A_CHL) )  &
                          /(CONVM*DEN(SA))
            EVCONC(M,I,XX)=EMASS(M,I,A_XXX)/(CONVM*DEN(XX))
            EVCONC(M,I,EC)=EMASS(M,I,A_EBC)/(CONVM*DEN(EC))                
            EVCONC(M,I,DU)=EMASS(M,I,A_DUS)/(CONVM*DEN(DU))


! Calculated particle number concentration per bin
! consistent GMD (mass-based)
!   unit:  # m^-2 s^-1
            EN(M,I)=(EVCONC(M,I,SU)+EVCONC(M,I,OC)+EVCONC(M,I,AM)        +  &              
                     EVCONC(M,I,NI)+EVCONC(M,I,MS)+EVCONC(M,I,SA)        +  &
                     EVCONC(M,I,XX)+EVCONC(M,I,EC)+EVCONC(M,I,DU) )/VPT(M,I) 
           ! print *,'EN',M,I,ELINDISM(M,I,XX),EN(M,I)

           end do
        end do


! Deallocate aerosol terms
        deallocate(EVCONC)
        deallocate(ELINDISM)


  end subroutine initEMNumberMass


!------------------------------------------------------------------

  subroutine getDensity(IMAX,temp,denecin,DPA,MASS,MPT,MPTW,ROOP,ROOPW)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate the average particle density for all size bins
    !
    !      interface
    !      ---------
    !
    !        input
    !           MASS      mass conc. per bin              [ng/m^3]
    !           MPT       total dry mass conc. per bin    [kg/m^3]
    !           MPTW      total wet mass conc. per bin    [kg/m^3]
    !           DPA       diameter of dry particles       [m]
    !           temp      air temperature                 [K]
    !           denecin   particle density of soot from organic.dat  [kg/m^3]
    !        output:
    !           ROOP      particle density of dry aerosol [kg/m^3]
    !           ROOPW     particle density of dry aerosol [kg/m^3]
    !
    !      method
    !      ------
    !      ROOP=(Msa+Morg+Mwa)*roo(sa+org+wa)+Mdu*roo(du))/mass
    !      implicitly assumed that the densities of sa and org are the same
    !      mass fraction  MF = (MSULF+Morg)/(MSULF+Morg+MH2O)
    !      MF is dim.-less fraction
    !      density roo(sa+org+wa)=ROOLM(XM,T)
    !
    !      reference
    !      ---------
    !      none
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

  implicit none

! input
     integer, intent(in)                             :: IMAX
     REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA     ! [m]
     real( dp), dimension(MMAX,IMAX,AMAX),intent(in) :: MASS     ! [ng/m^3]
     real( dp), dimension(MMAX,IMAX),intent(in)      :: MPT      ! [kg/m^3]
     real( dp), dimension(MMAX,IMAX),intent(in)      :: MPTW     ! [kg/m^3]
     real( dp),intent(in)                            :: temp     ! [K]
     real( dp),intent(in)                            :: denecin  ! [kg/m^3]

! output
     real( dp), dimension(MMAX,IMAX),intent(out)     :: ROOP     ! [kg/m^3]
     real( dp), dimension(MMAX,IMAX),intent(out)     :: ROOPW    ! [kg/m^3]

! local
     real( dp), dimension(MMAX,IMAX)   :: massocs
     real( dp), dimension(MMAX,IMAX)   :: MF
     real( dp), dimension(MMAX,IMAX)   :: deffpsoot
     real( dp)                         :: Dp0
     real( dp)                         :: ROOP1,ROOPW1
     integer                           :: M,I


! initialization
        deffpsoot(:,:)=1200._dp
        Dp0=2._dp*rp0

        do M=1,MMAX
         do I=1,IMAX
         
! effective density of soot particles [kg/m^3]
! COAG=2 or 5: Correction for fractal geometry of soot
           if ( (ICOAG.eq.2).or.(ICOAG.eq.5) ) then
             deffpsoot(M,I) = min( (denecin*(DPA(M,I)*1.e9/Dp0)**(Dfrac-2.7_dp)), denecin)    
             deffpsoot(M,I) = max(deffpsoot(M,I),1200.0_dp)
           else
             deffpsoot(M,I) = denecin
           endif

           !write(6,*) 'deffsoot',denecin,deffpsoot(M,I)

! sum of SULF and OC components per size bin
           massocs(M,I) = MASS(M,I,A_SUL) + MASS(M,I,A_MSA)                      &
                         +  MASS(M,I,A_OR1) + MASS(M,I,A_OR2) + MASS(M,I,A_OR3)   &
                         +  MASS(M,I,A_OR4) + MASS(M,I,A_OR5) + MASS(M,I,A_OR6)   &
                         +  MASS(M,I,A_OR7) + MASS(M,I,A_OR8) + MASS(M,I,A_OR9) 
         end do
        end do



        do M=1,MMAX
         do I=1,IMAX
      !!! predefined density is that of elemental carbon
           if (MPT(M,I).eq. 0.0_dp) then

             ROOP(M,I) = deffpsoot(M,I)
             ROOPW(M,I)= deffpsoot(M,I)
         
           else if ( massocs(M,I).gt.massmin )  then

      !!! density of dry particles (sulf+oc)
              MF(M,I)  = 1.0_dp
              ROOP1    = roolm(MF(M,I),temp)

      !!! density of wet particles (sulf+oc+wat)
      ! roo(sa+org+wa)=ROOLM(XM,T) where XM=(Msa+Morg)/(Msa+Morg+Mwa)
      ! ROOP=(Msa+Morg+Mwa)*roo(sa+org+wa)+Mdu*roo(du))/mass 
      ! implicitly assumed that the densities of sa and org are the same

              MF(M,I)   =  massocs(M,I) / (massocs(M,I) + MASS(M,I,A_WAT))

              ROOPW1    = roolm(MF(M,I),temp)

              ROOP(M,I) = ((  massocs(M,I)*ROOP1                                + &
                              MASS(M,I,A_AMI)*DENAM                             + &
                              MASS(M,I,A_NH4)*DENNI                             + &
                              MASS(M,I,A_NIT)*DENNI                             + &
                              MASS(M,I,A_SAL)*DENSA                             + &
                              MASS(M,I,A_CHL)*DENSA                             + &
                              MASS(M,I,A_XXX)*DENXX                             + &
                              MASS(M,I,A_EBC)*deffpsoot(M,I)                    + &
                              MASS(M,I,A_DUS)*DENDU )    *1.e-12_dp )           / &
                              MPT(M,I)


              ROOPW(M,I)= (( (massocs(M,I) + MASS(M,I,A_WAT))*ROOPW1            + &
                              MASS(M,I,A_AMI)*DENAM                             + &
                              MASS(M,I,A_NH4)*DENNI                             + &
                              MASS(M,I,A_NIT)*DENNI                             + &
                              MASS(M,I,A_SAL)*DENSA                             + &
                              MASS(M,I,A_CHL)*DENSA                             + &
                              MASS(M,I,A_XXX)*DENXX                             + &
                              MASS(M,I,A_EBC)*deffpsoot(M,I)                    + &
                              MASS(M,I,A_DUS)*DENDU )    *1.e-12_dp )           / &
                              MPTW(M,I)


           else
           
             ! ROOP(M,I)  = 1000._dp
             ! ROOPW(M,I) = 1000._dp

              ROOP(M,I) = ((  MASS(M,I,A_AMI)*DENAM                             + &
                              MASS(M,I,A_NH4)*DENNI                             + &
                              MASS(M,I,A_NIT)*DENNI                             + &
                              MASS(M,I,A_SAL)*DENSA                             + &
                              MASS(M,I,A_CHL)*DENSA                             + &
                              MASS(M,I,A_XXX)*DENXX                             + &
                              MASS(M,I,A_EBC)*deffpsoot(M,I)                    + &
                              MASS(M,I,A_DUS)*DENDU )    *1.e-12_dp )           / &
                              MPT(M,I)


              ROOPW(M,I)= ((  MASS(M,I,A_WAT)*RHOH2O                            + &
                              MASS(M,I,A_AMI)*DENAM                             + &
                              MASS(M,I,A_NH4)*DENNI                             + &
                              MASS(M,I,A_NIT)*DENNI                             + &
                              MASS(M,I,A_SAL)*DENSA                             + &
                              MASS(M,I,A_CHL)*DENSA                             + &
                              MASS(M,I,A_XXX)*DENXX                             + &
                              MASS(M,I,A_EBC)*deffpsoot(M,I)                    + &
                              MASS(M,I,A_DUS)*DENDU )    *1.e-12_dp )           / &
                              MPTW(M,I)

           endif
          
           ROOP(M,I)  = max(ROOP(M,I),1000._dp)
           ROOPW(M,I) = max(ROOPW(M,I),1000._dp)
!debug
          !write(6,*) 'DENS',M,I,MPT(M,I),massocs(M,I),ROOP(M,I),DPA(M,I),ROOPW(M,I)


         end do
       end do

  end subroutine getDensity

!------------------------------------------------------------------

  subroutine getTotalMass(IMAX,MASS,MPT,MPTW, MTOT, MTOTW,        &
                          MSULFTOT, MMSAPTOT, MIODATOT,           &
                          MNITRTOT, MAMMOTOT,                     &
                          MSALTTOT, MECBCTOT, MDUSTTOT, MXXXXTOT, & 
                          MORGCTOT, MORG1TOT, MORG2TOT, MORG3TOT, &
                          MORG4TOT, MORG5TOT, MORG6TOT, MORG7TOT, &
                          MORG8TOT, MORG9TOT, casoa )
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate the total mass per aerosol mode
    !
    !      interface
    !      ---------
    !
    !        input:
    !           MASS    particle mass of component      [ng/m^3]
    !
    !        output:
    !           MPT     total dry mass                  [kg/m^3]
    !           MPTW    total wet mass                  [kg/m^3]
    !           MTOT    total dry mass                  [ng/m^3]
    !           MTOTW   total wet mass                  [ng/m^3]
    !           MxxxxTOT compound total mass            [ng/m^3]
    !           CASOA   total SOA-x mass                [ug/m^3]
    !
    !      method
    !      ------
    !      summation of the mass of aerosol components in each size bin
    !      to get the total dry mass and total wet mass per bin
    !      and the mass of component per aerosol mode
    !
    !      reference
    !      ---------
    !      none
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

  implicit none

     integer, intent(in)                             :: IMAX
     real( dp), dimension(MMAX,IMAX,AMAX),intent(in) :: MASS     ! [ng/m^3]

     real( dp), dimension(MMAX,IMAX),intent(out)     :: MPT      ! [kg/m^3]
     real( dp), dimension(MMAX,IMAX),intent(out)     :: MPTW     ! [kg/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MTOT     ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MTOTW    ! [ng/m^3]

     real( dp), dimension(MMAX),intent(out)          :: MSULFTOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MSALTTOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MAMMOTOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MNITRTOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MMSAPTOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MIODATOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MXXXXTOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MECBCTOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MDUSTTOT ! [ng/m^3]

     real( dp), dimension(MMAX),intent(out)          :: MORGCTOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MORG1TOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MORG2TOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MORG3TOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MORG4TOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MORG5TOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MORG6TOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MORG7TOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MORG8TOT ! [ng/m^3]
     real( dp), dimension(MMAX),intent(out)          :: MORG9TOT ! [ng/m^3]

     real( dp), dimension(NSOA),intent(out)          :: casoa    ! [ug/m^3]

! local

     integer                           :: M,I

! TOTAL MASS CONCENTRATION
! in kg/m^3


        do M=1,MMAX

! initialize variables in mode M
          MTOT(M)=0._dp
          MTOTW(M)=0._dp
          MSULFTOT(M)=0._dp
          MAMMOTOT(M)=0._dp
          MNITRTOT(M)=0._dp
          MMSAPTOT(M)=0._dp
          MIODATOT(M)=0._dp
          MSALTTOT(M)=0._dp
          MXXXXTOT(M)=0._dp
          MECBCTOT(M)=0._dp
          MDUSTTOT(M)=0._dp
          MORGCTOT(M)=0._dp
          MORG1TOT(M)=0._dp
          MORG2TOT(M)=0._dp
          MORG3TOT(M)=0._dp
          MORG4TOT(M)=0._dp
          MORG5TOT(M)=0._dp
          MORG6TOT(M)=0._dp
          MORG7TOT(M)=0._dp
          MORG8TOT(M)=0._dp
          MORG9TOT(M)=0._dp


          do I=1,IMAX

! Add up discretized mass concentrations [ng/m^3]
            MSULFTOT(M)=MSULFTOT(M)+MASS(M,I,A_SUL)
            MAMMOTOT(M)=MAMMOTOT(M)+MASS(M,I,A_AMI)+MASS(M,I,A_NH4)
            MNITRTOT(M)=MNITRTOT(M)+MASS(M,I,A_NIT)
            MMSAPTOT(M)=MMSAPTOT(M)+MASS(M,I,A_MSA)
            MIODATOT(M)=MIODATOT(M)+MASS(M,I,A_IO3)
            MSALTTOT(M)=MSALTTOT(M)+MASS(M,I,A_SAL)+MASS(M,I,A_CHL)
            MXXXXTOT(M)=MXXXXTOT(M)+MASS(M,I,A_XXX)
            MECBCTOT(M)=MECBCTOT(M)+MASS(M,I,A_EBC)
            MDUSTTOT(M)=MDUSTTOT(M)+MASS(M,I,A_DUS)
            MORG1TOT(M)=MORG1TOT(M)+MASS(M,I,A_OR1)
            MORG2TOT(M)=MORG2TOT(M)+MASS(M,I,A_OR2)
            MORG3TOT(M)=MORG3TOT(M)+MASS(M,I,A_OR3)
            MORG4TOT(M)=MORG4TOT(M)+MASS(M,I,A_OR4)
            MORG5TOT(M)=MORG5TOT(M)+MASS(M,I,A_OR5)
            MORG6TOT(M)=MORG6TOT(M)+MASS(M,I,A_OR6)
            MORG7TOT(M)=MORG7TOT(M)+MASS(M,I,A_OR7)
            MORG8TOT(M)=MORG8TOT(M)+MASS(M,I,A_OR8)
            MORG9TOT(M)=MORG9TOT(M)+MASS(M,I,A_OR9)

! Total mass per bin (dry and wet) [kg/m^3]

            MPT(M,I)=( MASS(M,I,A_SUL)+                   &
                       MASS(M,I,A_OR1)+MASS(M,I,A_OR2)+MASS(M,I,A_OR3) + &
                       MASS(M,I,A_OR4)+MASS(M,I,A_OR5)+MASS(M,I,A_OR6) + &
                       MASS(M,I,A_OR7)+MASS(M,I,A_OR8)+MASS(M,I,A_OR9) + &
                       MASS(M,I,A_AMI)+MASS(M,I,A_NH4)                 + &
                       MASS(M,I,A_NIT)+MASS(M,I,A_MSA)+MASS(M,I,A_IO3) + &
                       MASS(M,I,A_SAL)+MASS(M,I,A_CHL)+MASS(M,I,A_XXX) + &
                       MASS(M,I,A_EBC)+MASS(M,I,A_DUS) )*1.e-12_dp

            MPTW(M,I)=( MASS(M,I,A_SUL)+                  &
                       MASS(M,I,A_OR1)+MASS(M,I,A_OR2)+MASS(M,I,A_OR3) + &
                       MASS(M,I,A_OR4)+MASS(M,I,A_OR5)+MASS(M,I,A_OR6) + &
                       MASS(M,I,A_OR7)+MASS(M,I,A_OR8)+MASS(M,I,A_OR9) + &
                       MASS(M,I,A_AMI)+MASS(M,I,A_NH4)                 + &
                       MASS(M,I,A_NIT)+MASS(M,I,A_MSA)+MASS(M,I,A_IO3) + &
                       MASS(M,I,A_SAL)+MASS(M,I,A_CHL)+MASS(M,I,A_XXX) + &
                       MASS(M,I,A_EBC)+MASS(M,I,A_DUS)                 + &
                       MASS(M,I,A_WAT) )*1.e-12_dp


             MTOT(M)=MTOT(M)+MPT(M,I)*1.e12
             MTOTW(M)=MTOTW(M)+MPTW(M,I)*1.e12

          end do
! Sum organic carbon mass
          MORGCTOT(M) = MORG1TOT(M)+MORG2TOT(M)+MORG3TOT(M)      + &             
                        MORG4TOT(M)+MORG5TOT(M)+MORG6TOT(M)      + &
                        MORG7TOT(M)+MORG8TOT(M)+MORG9TOT(M)
! debug
! write(6,*) "totmass orgc M", M, MORGCTOT(M)
! write(6,*) "biognic orgc M", M, MORG1TOT(M),MORG2TOT(M),MORG3TOT(M)
! write(6,*) "aromatc orgc M", M, MORG4TOT(M),MORG5TOT(M),MORG6TOT(M)
! write(6,*) "primary orgc M", M, MORG7TOT(M),MORG8TOT(M),MORG9TOT(M)
        end do

! Finally calculate the total SOA(part) mass concentrations [ug/m^3]
        casoa(1) = 1.e-3_dp*(MORG1TOT(NU)+MORG1TOT(NA)+MORG1TOT(AI)+MORG1TOT(AS)+MORG1TOT(CS))
        casoa(2) = 1.e-3_dp*(MORG2TOT(NU)+MORG2TOT(NA)+MORG2TOT(AI)+MORG2TOT(AS)+MORG2TOT(CS))
        casoa(3) = 1.e-3_dp*(MORG3TOT(NU)+MORG3TOT(NA)+MORG3TOT(AI)+MORG3TOT(AS)+MORG3TOT(CS))
        casoa(4) = 1.e-3_dp*(MORG4TOT(NU)+MORG4TOT(NA)+MORG4TOT(AI)+MORG4TOT(AS)+MORG4TOT(CS))
        casoa(5) = 1.e-3_dp*(MORG5TOT(NU)+MORG5TOT(NA)+MORG5TOT(AI)+MORG5TOT(AS)+MORG5TOT(CS))
        casoa(6) = 1.e-3_dp*(MORG6TOT(NU)+MORG6TOT(NA)+MORG6TOT(AI)+MORG6TOT(AS)+MORG6TOT(CS))
        casoa(7) = 1.e-3_dp*(MORG7TOT(NU)+MORG7TOT(NA)+MORG7TOT(AI)+MORG7TOT(AS)+MORG7TOT(CS))
        casoa(8) = 1.e-3_dp*(MORG8TOT(NU)+MORG8TOT(NA)+MORG8TOT(AI)+MORG8TOT(AS)+MORG8TOT(CS))
        casoa(9) = 1.e-3_dp*(MORG9TOT(NU)+MORG9TOT(NA)+MORG9TOT(AI)+MORG9TOT(AS)+MORG9TOT(CS))


  end subroutine getTotalMass

!------------------------------------------------------------------

  subroutine restoreCoarseMode(IMAX,DPA,DLINDP,DEN,                     &
                               MSULFTOT,MMSAPTOT,MAMMOTOT,MORGCTOT,     &
                               MSALTTOT,MDUSTTOT,  MASS,N)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      restore number and mass of aerosol compounds in size bins
    !      of the coarse mode after fog/cloud evaporation
    !
    !      interface
    !      ---------
    !
    !        input:
    !           DPA     dry particle diameter           [m]
    !           MASS    mass concentration              [ng/m^3]
    !           DLINDP  linear width of bin             [m]
    !           DEN     density of compound             [kg/m^3]
    !
    !        output:
    !           MASS    mass concentration              [ng/m^3]
    !           N       number concentration            [part/m^3]
    !
    !      method
    !      ------
    !      Calculate mass conc of compoound in each bin I
    !      from log-normal mass distribution
    !      GMD is geometric-mean mass diameter!!!
    !      we assume average particle density (DEN) to be
    !      constant in all bins.
    !
    !      reference
    !      ---------
    !      none
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

  implicit none

     integer, intent(in)                                 :: IMAX
     real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)    :: DPA
     real( dp), DIMENSION(MMAX,IMAX), intent(in)         :: DLINDP
     real( dp), dimension(iamax), intent(in)             :: DEN

     real( dp), dimension(MMAX),intent(in)               :: MSULFTOT
     real( dp), dimension(MMAX),intent(in)               :: MMSAPTOT
     real( dp), dimension(MMAX),intent(in)               :: MAMMOTOT
     real( dp), dimension(MMAX),intent(in)               :: MORGCTOT
     real( dp), dimension(MMAX),intent(in)               :: MSALTTOT
     real( dp), dimension(MMAX),intent(in)               :: MDUSTTOT

!in/out
     real( dp), dimension(MMAX,IMAX,AMAX),intent(in out) :: MASS
     real( dp), DIMENSION(MMAX,IMAX), intent(in out)     :: N

! local
     real( dp)                                         :: VDRY
     REAL( dp),allocatable,dimension(:,:)              :: VCONCCS
     REAL( dp),allocatable,dimension(:,:)              :: LINDISMCS

     integer          :: I


! Allocate aerosol terms
        if (.not. allocated(VCONCCS))      ALLOCATE(VCONCCS(IMAX,iamax))
        if (.not. allocated(LINDISMCS))    ALLOCATE(LINDISMCS(IMAX,iamax))



        do I=1,IMAX

            LINDISMCS(I,DU)=(MDUSTTOT(CS)/(SQRT(2._dp*pi)*DPA(CS,I)*LOG(SIG(CS)))) &
                           *EXP(-0.5_dp*(LOG(DPA(CS,I)/GMD(CS))/LOG(SIG(CS)))**2._dp)
            LINDISMCS(I,SA)=(MSALTTOT(CS)/(SQRT(2._dp*pi)*DPA(CS,I)*LOG(SIG(CS)))) &
                           *EXP(-0.5_dp*(LOG(DPA(CS,I)/GMD(CS))/LOG(SIG(CS)))**2._dp)
            LINDISMCS(I,SU)=(MSULFTOT(CS)/(SQRT(2._dp*pi)*DPA(CS,I)*LOG(SIG(CS)))) &
                           *EXP(-0.5_dp*(LOG(DPA(CS,I)/GMD(CS))/LOG(SIG(CS)))**2._dp)
            LINDISMCS(I,MS)=(MMSAPTOT(CS)/(SQRT(2._dp*pi)*DPA(CS,I)*LOG(SIG(CS)))) &
                           *EXP(-0.5_dp*(LOG(DPA(CS,I)/GMD(CS))/LOG(SIG(CS)))**2._dp)
            LINDISMCS(I,AM)=(MAMMOTOT(CS)/(SQRT(2._dp*pi)*DPA(CS,I)*LOG(SIG(CS)))) &
                           *EXP(-0.5_dp*(LOG(DPA(CS,I)/GMD(CS))/LOG(SIG(CS)))**2._dp)
            LINDISMCS(I,OC)=(MORGCTOT(CS)/(SQRT(2._dp*pi)*DPA(CS,I)*LOG(SIG(CS)))) &
                           *EXP(-0.5_dp*(LOG(DPA(CS,I)/GMD(CS))/LOG(SIG(CS)))**2._dp)

! Mass concentrations in ng/m3
            MASS(CS,I,A_DUS)=LINDISMCS(I,DU)*DLINDP(CS,I)
            MASS(CS,I,A_SAL)=LINDISMCS(I,SA)*DLINDP(CS,I)
            MASS(CS,I,A_SUL)=LINDISMCS(I,SU)*DLINDP(CS,I)
            MASS(CS,I,A_MSA)=LINDISMCS(I,MS)*DLINDP(CS,I)
            MASS(CS,I,A_NH4)=LINDISMCS(I,AM)*DLINDP(CS,I)
            MASS(CS,I,A_OR2)=LINDISMCS(I,OC)*DLINDP(CS,I)

! Volume concentration in m3/m3
            VCONCCS(I,DU)=MASS(CS,I,A_DUS)/(CONVM*DEN(DU))
            VCONCCS(I,SA)=MASS(CS,I,A_SAL)/(CONVM*DEN(SA))
            VCONCCS(I,SU)=MASS(CS,I,A_SUL)/(CONVM*DEN(SU))
            VCONCCS(I,MS)=MASS(CS,I,A_MSA)/(CONVM*DEN(MS))
            VCONCCS(I,AM)=MASS(CS,I,A_NH4)/(CONVM*DEN(AM))
            VCONCCS(I,OC)=MASS(CS,I,A_OR2)/(CONVM*DEN(OC))

! Calculate dry volume in m3
! VDRY is better than VPT (wet diameter) to return to the
! original number size distribution
            VDRY=(pi/6._dp)*DPA(CS,I)**3._dp

! Calculated particle number concentration per bin
! consistent GMD (mass-based)
            N(CS,I) = ( VCONCCS(I,DU)+VCONCCS(I,SA)+      &
                        VCONCCS(I,SU)+VCONCCS(I,MS)+      &
                        VCONCCS(I,AM)+VCONCCS(I,OC) )   / VDRY

          !  print *,'restore ',I,VPT(CS,I),N(CS,I)

        end do


! Deallocate aerosol terms
        deallocate(VCONCCS)
        deallocate(LINDISMCS)


  end subroutine restoreCoarseMode

!------------------------------------------------------------------

end module gde_aerosol_props
