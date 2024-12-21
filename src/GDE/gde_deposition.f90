! <gde_deposition.f90 - A component of the Multicomponent
!                     Aerosol Dynamics Model MAFOR>
!*****************************************************************************! 
!* 
!*    Copyright (C) 2011-2023 Matthias Steffen Karl
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
!*    except:
!*    routine DEPOSITPAR written by Liisa Pirjola
!*    routine DIFFPAR written by Liisa Pirjola
!*    routine DEPOZHANG01 written by Mona Kurppa
!*    routine DEPOFROUGH written by Tareq Hussein, and belonging
!*    functions MIUAIR and LAMBDA written by Tareq Hussein
!*
!*    Default values for deposition parameters:
!*       ustar    = 1.17      ! friction velocity [m/s]
!*       znot     = 0.001     ! surface roughness [m]
!*       ADEP     = 1.7       ! dry dep. parameter A
!*       BDEP     = 51.8      ! dry dep. parameter B
!*       ZCAP     = 0.30      ! canopy height [m]
!*       DCOL     = 0.005     ! collection size [m]
!*       Fplus    = 0.0       ! roughness parameter [-]
!* 
!*****************************************************************************!
module gde_deposition

    use gde_constants,  only      : pi,k_B,g,c_vKar,M_air,R_gas
    use gde_input_data, only      : MMAX
    use gde_plume,      only      : znot,ustar
    use gde_plume,      only      : ADEP,BDEP,ZCAP,dcol,Fplus


    private
   
    public :: depositwall
    public :: depositpar
    public :: depocanopy
    public :: depozhang01
    public :: settling
    public :: depofrough
    public :: wetscavbulk
    public :: wetscavsize

! KPP DP - Double precision kind
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)

  contains


subroutine depositwall(press,temp,DPA,DIL,VC,AS,AD,depowall,IMAX)
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
    !      calculates particle diffusion and sedimentation to walls
    !      in chamber experiments
    !
    !      interface
    !      ---------
    !        input:
    !           press    [Pa]
    !           temp     [K]
    !           DPA      [m]
    !           DIL      [1/s]
    !           VC       [m3]
    !           AS       [m2]
    !           AD       [m2]
    !        output:
    !           depowall [1/s]
    !
    !      method
    !      ------
    !      COSIMA - a computer program simulating the dynamics of 
    !      fractal aerosols
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      Parameterization taken from:
    !      K.-H. Naumann, Aerosol Science, 34, 1371-1397, 2003
    !      COSIMA - a computer program simulating the dynamics of 
    !      fractal aerosols
    !
    !------------------------------------------------------------------

    implicit none

    INTEGER, intent(in)                               :: IMAX
    real( dp), intent(in)                             :: temp
    real( dp), intent(in)                             :: press
    real( dp), intent(in)                             :: DIL
    real( dp), intent(in)                             :: VC,AS,AD
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPA

    REAL( dp), DIMENSION(MMAX,IMAX), intent(out)      :: depowall

  ! local
    REAL( dp), DIMENSION(MMAX,0:(IMAX+1))             :: RP

    REAL( dp)                              :: MYY,DIFFCO,CC
    REAL( dp)                              :: delta,pres
    REAL( dp)                              :: adif,adil,ased
    REAL( dp)                              :: DENSPARCONST

    !  APA to tune delta in Eq. 40 of Naumann 2003. Use APA=1
    real( dp), parameter                   :: APA = 1.0_dp
    !  KD [m] adopted from Bunz and Dlugi, 1991
    real( dp), parameter                   :: KD=0.005_dp
    !   A [dimless] adopted from  Bunz and Dlugi, 1991 
    real( dp), parameter                   :: A=0.274_dp 

    INTEGER                                :: M,I

      pres=press/101325._dp
      !!!DENS= 1.2929*273.15/temp*pres/1.      !air density in kg/m^3, P in atm
      DENSPARCONST= 1000.  !kg/m^3  for all particles (1.assumption)


!!! geometry from input incham.dat
!!!      AD=130._dp       !m^2      EUPHORE, diffusion surface (diameter: 9.10m)
!!!      AS=65._dp        !m^2      EUPHORE, sedimentation surface
!!!      VC=177._dp       !m^3      EUPHORE, Volume chamber  
!!!      DIL=7.E-6_dp     !1/s      EUPHORE, dilution rate (fixed at the moment) 

!      P = 1.        !atm

      do M=1,MMAX
       do I=1,IMAX

         RP(M,I)=DPA(M,I)*0.5_dp
        ! get MYY in kg/m/s, DIFFCO in m^2/s
         CALL diffpar(RP(M,I),pres,temp,MYY,DIFFCO,CC)
        ! Calculate diffusive boundary layer thickness
        ! D_0 in Eq. 40 has the value of 1 (and the same unit as D) and its purpose is
        ! to make the term dimensionless [pers. commun. K.-H. Naummann, 08.06.2009].
         delta=(DIFFCO/1._dp)**(A*APA)
         delta=KD*delta
         adif=(DIFFCO*AD)/(delta*VC)      !1/s Naumann 2003
         !adif=CDIFF*SQRT(DIFFCO)         !1/s Verheggen 2006 CDIFF=3.6E-3 or higher
         ased=4._dp*pi*DENSPARCONST*(RP(M,I)**3)*g*DIFFCO*AS
         ased=ased/(3._dp*k_B*temp*VC)    !1/s
         adil=DIL                         !1/s                                  
         depowall(M,I)=adif+ased+adil     !1/s

      !    write(6,*) 'wall',M,I,delta,adif,ased,adil,depowall(M,I)

       end do
      end do

  end subroutine depositwall


subroutine depositpar(press,temp,DENSPAR,DPA,mbh,owf,depo,IMAX)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Dr. Liisa Pirjola
    !      Docent
    !      Department of Physics
    !      University of Helsinki
    !      P.O Box 64, FI-00014 Helsinki, Finland
    !      Department of Technology
    !      Metropolia University of Applied Sciences
    !      P.O. Box 4071, FI-01600 Vantaa, Finland
    !
    !      purpose
    !      -------
    !      calculates particle dry deposition rate
    !
    !      interface
    !      ---------
    !
    !        input:
    !           press    [Pa]
    !           temp     [K]
    !           mbh      [m]
    !           owf      [-]
    !           denspar  [kg/m3]
    !           DPA      [m]
    !         output:
    !           depo     [1/s]
    !
    !      method
    !      ------
    !      Parameterization for dry deposition velicity 
    !      from Schack et al. (1985).
    !      Taken from MONO32
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      Schack Jr., C. J., Pratsinis, S. E. and Friedlander, S. K. 1985. 
    !        A general correlation for deposition of suspended particles from
    !        turbulent gases to completely rough surfaces,
    !        Atmos. Environ. 19, 953-960.    
    !
    !      note: znot from input is used as if it was in "cm"
    !------------------------------------------------------------------

    implicit none

    integer, intent(in)                               :: IMAX
    real( dp), intent(in)                             :: temp,press,mbh
    real( dp), intent(in)                             :: owf
   ! REAL( dp), intent(in)                             :: ustar,znot
   ! REAL( dp), intent(in)                             :: ADEP,BDEP    
    real( dp), dimension(MMAX,IMAX), intent(in)       :: DENSPAR
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPA
    real( dp), dimension(MMAX,IMAX), intent(out)      :: depo
    real( dp), dimension(MMAX,0:(IMAX+1))             :: RP

    real( dp)                              :: DENS,pres
    real( dp)                              :: Sc,MYY,DIFFCO,X,VS,VG,CC
    real( dp)                              :: DENSPA
    real( dp)                              :: znotcm

    integer                                :: M,I

      pres=press/101325._dp
      DENS= 1.2929*273.15/temp*pres/1.      !air density in kg/m^3, P in atm
      DENSPA= 1400.  !kg/m^3  for all particles (1.assumption)
      !  water (Sehmel and Sutter, 1974)
      !USTAR=1.17    !m/s
      !ZNOT = 0.001   !m
      !ADEP=1.7
      !BDEP=51.8
!      P = 1.        !atm
      znotcm=znot
      ! set ZNOT over water/ice (water: 0.001 m/s)
      ! value from dispers.dat is taken
      ! if 0.<= owf <= 0.3
      if (owf.eq.1._dp) then
        znotcm=0.001
      else if ((owf.gt.0.3).and.(owf.lt.1._dp)) then
        znotcm=0.0005
      endif
      do M=1,MMAX
       do I=1,IMAX
         RP(M,I)=DPA(M,I)*0.5_dp
         ! get MYY in kg/m/s, DIFFCO in m^2/s
         CALL diffpar(RP(M,I),pres,temp,MYY,DIFFCO,CC)
         MYY=MYY/DENS     !m^2/s
         Sc = MYY/DIFFCO  !laaduton
         X = 2._dp*RP(M,I)*(ustar/znotcm/MYY)**(0.5_dp)*Sc**(1/3._dp)
         !VS = 4._dp/3._dp*pi*(RP(M,I))**3*DENSPA*g*DIFFCO/k_B/temp     !m/s
         VS = 4._dp/3._dp*pi*(RP(M,I))**3*DENSPAR(M,I)*g*DIFFCO/k_B/temp     !m/s
         VG = DIFFCO/(2*RP(M,I))*(ADEP*X + BDEP*X**3) + VS     !m/s
         !write(6,*) DPA(M,I),VG     
         depo(M,I) = VG/mbh       !1/s
       end do
      end do

      !stop

  end subroutine depositpar  

subroutine depozhang01(temp,DENSPAR,DPA,mbh,depo,IMAX)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Mona Kurppa
    !      Atmospheric Composition Research
    !      Finnish Meteorological Institute
    !      Helsinki, Finland
    !
    !      purpose
    !      -------
    !      calculates particle dry deposition rate to vegetation surface
    !      uses particle wet diameter (m)
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     [K]
    !           mbh      [m]
    !           DPA      [m]
    !        !   denspar [kg/m3]
    !
    !        !   ustar   [m/s]  - friction velocity
    !        !   znot    [m]    - surface roughness length
    !        output:
    !           depo  [1/s]
    !    
    !      method
    !      ------
    !      Parameterization for dry deposition velicity for smooth and 
    !      rough surfaces.
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      Zhang L., Gong, S., Padro, J., Barrie, L., 2001.
    !        A size-seggregated particle dry deposition scheme
    !        for an atmospheric aerosol module,
    !        Atmos. Environ., 35, 549-560.
    !
    !------------------------------------------------------------------

    implicit none

    INTEGER, intent(in)                               :: IMAX
    REAL( dp), intent(in)                             :: temp
    REAL( dp), intent(in)                             :: mbh

    REAL( dp), DIMENSION(MMAX,IMAX), intent(in)       :: DENSPAR
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPA

    REAL( dp), DIMENSION(MMAX,IMAX), intent(out)      :: depo

    real( dp), parameter :: adn = 1.3                ! air density
    real( dp), parameter :: alpha_d = 1.5            !
    real( dp), parameter :: am_airmol = 4.8096E-26   ! Average mass of an air molecule
    real( dp), parameter :: dcol = 0.002             ! collector size (m)
    real( dp), parameter :: gamma_d = 0.56           ! 

    ! SALSA configuration in street canyon case
    ! real( dp), parameter :: pdn = 1500.0         ! particle density (kg/m^3)
    ! real( dp), parameter :: ustar = 0.52        ! friction velocity (m/s)
    ! real( dp), parameter :: znot = 0.4          ! roughness length (m)  
    ! real( dp), parameter :: zC = 10.0           ! canopy height (m)

    real( dp)                              :: avis   ! molecular viscocity
    real( dp)                              :: kvis   ! kinematic viscocity
    real( dp)                              :: lambda !
    real( dp)                              :: mdiff  !
    real( dp)                              :: ra     ! aerodynamic resistance
    real( dp)                              :: rs     ! quasi-laminar resistance
    real( dp)                              :: Sc     ! Schmidt number
    real( dp)                              :: St     ! Stokes number
    real( dp)                              :: va     !
    real( dp)                              :: CC     ! slip correction
    real( dp)                              :: Kn     ! Knudsen number
    real( dp)                              :: vc     ! settling velocity
    real( dp)                              :: DENSPA ! particle density
    real( dp)                              :: zC     ! reference height
      
    INTEGER                                :: M,I

      DENSPA= 1400.  !kg/m^3  for all particles (1.assumption)
      
      ! Molecular viscosity of air (Eq. 4.54)
      avis = 1.8325E-5_dp * ( 416.16_dp / ( temp + 120.0_dp ) ) *  &
             ( temp / 296.16_dp )**1.5

      ! Kinematic viscocity
      kvis = avis / adn

      ! Thermal velocity of an air molecule (Eq. 15.32)
      va = SQRT( 8.0 * k_B * temp / ( pi * am_airmol ) )

      ! Mean free path (m) (Eq. 15.24)
      lambda = 2.0 * avis / ( adn * va )

      ! Aerodynamic resistance
      ! evaluated at canopy top = zC (here 10 m)
      zC = 10.0_dp
      if (mbh <= znot) then
         zC=znot
      else if (mbh <= zC) then
         zC=mbh
      else
        zC = 10.0_dp
      endif
      ra = LOG(zC / znot) / (0.4_dp * ustar)  !znot = z0

      do M=1,MMAX
       do I=1,IMAX

         ! Knudsen number (Eq. 15.23)
         Kn = MAX( 1.0E-2, lambda / ( DPA(M,I) * 0.5_dp ) ) ! To avoid underflow

         ! Cunningham slip-flow correction
         CC = 1.0 + Kn * ( 1.257_dp + 0.4_dp * EXP( -1.1_dp / Kn ) )

         ! Critical fall speed i.e. settling velocity  (Eq. 20.4)
         !vc = MIN(1.0_dp, (DPA(M,I))**2 * (DENSPA - adn) * g * CC / (18.0_dp * avis))
         vc = MIN(1.0_dp, (DPA(M,I))**2 * (DENSPAR(M,I) - adn) * g * CC / (18.0_dp * avis))

         ! Stokes number
         !!! St = vc * ustar**2 / kvis  ! for surfaces with bluff roughness elements
         ! for vegetated surfaces
         St =  vc * ustar / (g * dcol)

         ! Particle diffusivity coefficient (Eq. 15.29)
         mdiff = (k_B * temp * CC) / (3.0_dp * pi * avis * DPA(M,I) )

         ! Schmidt number
         Sc = kvis / mdiff

         ! The overall quasi-laminar resistance for particles (Zhang et al., Eq. 5)
         ! EPSILON(1.0) is the smallest epsilon for which 1.0 + epsilon > 1.0
         ! Rs  = 1 / 3.0 * ustar * R1 * (EB + EIM + EIN)
         ! R1  = EXP(-St**0.5)
         ! EB  = Sc**(-gamma_d)
         ! EIM = St / (alpha_d + St))**2
         ! EIN = 0.5*(DPA/dcol)**2
         rs = MAX(EPSILON(1.0), (3.0 * ustar * EXP(-St**0.5) *                &
                                  ( Sc**(-gamma_d)                         +  &  ! EB
                                    (St / (alpha_d + St))**2               +  &  ! EIM
                                    0.5 * (DPA(M,I) / dcol)**2                &  ! EIN
                                  )                                           &
                                 )  )

         rs = 1.0 / rs

         ! Total deposition velocity (m/s)
         depo(M,I) = vc + 1.0 / (ra + rs)

         !write(6,*) DPA(M,I), depo(M,I)

         depo(M,I) = depo(M,I)/mbh    !1/s
         
       end do
      end do

     ! stop

  end subroutine depozhang01

subroutine depocanopy(press,temp,DENSPAR,DPA,mbh,u10,depo,IMAX)
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
    !      calculates particle dry deposition rate to a rough surface
    !
    !      interface
    !      ---------
    !
    !        input:
    !           press [Pa]
    !           temp  [K]
    !           mbh   [m]
    !           denspar [kg/m3]
    !           utop  [m/s]
    !           dcol  [m]    - collector size
    !           ustar [m/s]  - friction velocity
    !           znot  [m]    - surface roughness length
    !           ZCAP  [m]    - canopy height
    !        output:
    !           depo  [1/s]
    !    
    !      method
    !      ------
    !      Parameterization for dry deposition velicity for smooth and 
    !      rough surfaces.
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      Kouznetsov, R. and Sovief, M. 2012.
    !      A methodology for evaluation of vertical dispersion and dry deposition
    !      of atmospheric aerosol, J. Geophys. Res., 117, D01202, 1-17,
    !      doi:10.1029/2011JD016366.
    !
    !------------------------------------------------------------------

    implicit none

    INTEGER, intent(in)                               :: IMAX
    REAL( dp), intent(in)                             :: temp,press,mbh,u10
  !  REAL( dp), intent(in)                             :: ustar,znot
  !  REAL( dp), intent(in)                             :: dcol,ZCAP
    REAL( dp), DIMENSION(MMAX,IMAX), intent(in)       :: DENSPAR
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPA
    REAL( dp), DIMENSION(MMAX,IMAX), intent(out)      :: depo
    REAL( dp), DIMENSION(MMAX,0:(IMAX+1))             :: RP

    real( dp)                              :: DENS,pres
    real( dp)                              :: Sc,MYY,DIFFCO,VS,CC
    real( dp)                              :: restar,vdif,vint,vimp
    real( dp)                              :: stoken,effimp,vdzo,ra
    real( dp)                              :: acap,ztop,taup
    real( dp)                              :: recol,stockeff
    real( dp)                              :: DENSPA

    INTEGER                                :: M,I

      pres=press/101325._dp
      DENS= 1.2929*273.15/temp*pres/1.      !air density in kg/m^3, P in atm
      DENSPA= 1400.  !kg/m^3  for all particles (1.assumption)
      !USTAR=1.17    !m/s
      !ZNOT = 0.001   !m
      acap=(ustar/u10)*dcol
      if (ZCAP <= znot) then
         ztop=znot
      else
         ztop=ZCAP
      endif      
!      P = 1.        !atm
      do M=1,MMAX
       do I=1,IMAX
         RP(M,I)=DPA(M,I)*0.5
         ! get MYY in kg/m/s, DIFFCO in m^2/s
         CALL diffpar(RP(M,I),pres,temp,MYY,DIFFCO,CC)
         ! Relaxation time
         !taup   = CC * DENSPA * DPA(M,I)**2. /(18.*MYY) 
         taup   = CC * DENSPAR(M,I) * DPA(M,I)**2. /(18.*MYY)
         MYY=MYY/DENS     !m^2/s kinematic viscosity
         Sc = MYY/DIFFCO
         ! Canopy Reynolds number Re*
         restar = ustar * acap/MYY 
         ! Collector Reynolds number
         recol  = restar * (u10/ustar)**2.
         ! Stoke number
         stoken = 2*taup * u10/dcol 
         ! Diffusion controlled deposition vdif
         vdif   = ustar * 2. * restar**(-0.5) * Sc**(-2./3.)
         ! Interception controlled deposition vint
         vint   = ustar * 80. * (DPA(M,I)/acap)**2. * restar**(0.5)
         ! effective Stoke number
         stockeff = stoken - recol**(-0.5)
         if (stockeff > 0.15) then
           effimp = exp( (-0.1/(stockeff-0.15)) - (1./sqrt(stockeff-0.15)) )
         else
           effimp = 0.0
         endif    
         ! Impaction controlled deposition vimp
         vimp   = ustar * (2*ustar/u10) * effimp  * (stoken - ((ustar*restar**(-0.5))/u10) )
         ! Settling velocity VS       
         !VS     = 4./3.*pi*(RP(M,I))**3*DENSPA *g*DIFFCO/k_B/temp     !m/s
         VS     = 4./3.*pi*(RP(M,I))**3*DENSPAR(M,I) *g*DIFFCO/k_B/temp     !m/s
         ! Deposition velocity within canopy layer (at roughness lenth z0)
         vdzo   = vdif + vint + vimp + VS
         ! Aerodynamic resistance ra = int( dz/K(z) ) from z0 to z1
         ! crude assumption of a logarithmic profile
         ra = (1/(c_vKar*ustar)) * log10(ztop/znot)
         ! Deposition velocity above canopy        
         depo(M,I) = (1./vdzo)*exp((-1.)*VS*ra) + (1/VS)*(1.-exp((-1.)*VS*ra))
         depo(M,I) = 1/depo(M,I)      !m/s
         !write(6,*) 'depocanopy',M,I,depo(M,I),vdif,vint,vimp,VS
         !write(6,*) DPA(M,I),depo(M,I)
         depo(M,I) = depo(M,I)/mbh    !1/s
       end do
      end do

      !stop

  end subroutine depocanopy

subroutine diffpar(RP,pres,temp,MYY,DIFFCO,CC)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Dr. Liisa Pirjola
    !      Docent
    !      Department of Physics
    !      University of Helsinki
    !      P.O Box 64, FI-00014 Helsinki, Finland
    !      Department of Technology
    !      Metropolia University of Applied Sciences
    !      P.O. Box 4071, FI-01600 Vantaa, Finland
    !
    !      purpose
    !      -------
    !      calculates Slip correction factor Cc and
    !      calculates Brownian diffusion coefficient
    !
    !      interface
    !      ---------
    !
    !        input:
    !           press  [Pa]
    !           temp   [K]
    !           RP     [m]    - radius
    !        output:
    !           MYY    [kg/m/s]
    !           DIFFCO [m^2/s]
    !           CC     [-] - slip correction factor
    !
    !      method
    !      ------
    !      taken from MONO32
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

    REAL( dp), intent(in)                  :: temp,pres,RP
    REAL( dp), intent(out)                 :: MYY,DIFFCO,CC

    REAL( dp)                              :: LAMB,KN

       LAMB= (6.73e-8_dp*temp*(1.+(110.4_dp/temp)))/ &
          (296._dp*pres*1.373_dp)
       MYY= (1.832e-5_dp*(temp**(1.5_dp))*406.4_dp)/ &
          (5093._dp*(temp+110.4_dp))

       KN= LAMB/RP
       CC = 1._dp + (KN*(1.142_dp+(0.558_dp*EXP((-.999_dp)/KN))))
       DIFFCO = (k_B*temp*CC)/(6._dp*pi*MYY*RP)

  end subroutine diffpar


subroutine depofrough(press,temp,DENSPAR,DPA,mbh,depo,IMAX)
    !--------------------------------------------------------------------------
    !      author
    !      -------
    !      Tareq Hussein
    !      Professor
    !      University of Jordan
    !      School of Science
    !      Department of Physics
    !      Amman, 11942 Jordan
    !      Mobile: +962 779 483608
    !      Tel:    +962 6 5355000, ext: 22060
    !      Fax:    +962 6 5300253
    !      e-mail: t.hussein@ju.edu.jo
    !      Currently:
    !      Visiting Professor
    !      University of Helsinki
    !      Institute for Atmospheric and Earth System Research (INAR)
    !      PL 64, FI-00014 UHEL
    !      Helsinki, Finland
    !      Mobile: +358 400 867890
    !
    ! This is a routine extracted from the main deposition routine to be used
    ! as a stand alone based on the Model developed by Tareq Hussein and
    ! Published in J Aerosol Sci Tech 2012 46 0044-0059 (Hussein et al.)
    !
    ! Inputs:
    ! Dp     = particle diameter, and it can be a row matrix [m].
    !     ==== DPA(M,I)
    ! rohp   = particle density [kg/m3] for each particle size bin. This can be
    !         a matrix of the same length as Dp or just a single value. 
    !     ==== DENSPAR(M,I)
    ! u      = Friction velocity near the surface [m/s]
    !     ====      ustar [m/s]  - friction velocity
    ! F_plus = dimensionless roughness height of the surface [--]   Fplus
    ! T      = Temperature [K]
    ! P      = Presure [Pa]
    ! g      = Acceleration of gravity [m/s2]  (constant)
    !
    ! Outputs are the deposition velocity [m/s] onto three different
    ! orientations of sufraces: facing up as GROUND, facing down as CEILING,
    ! and vertical as WALLS.
    !
    !   original matlab code by:
    !   tareq.hussein@helsinki.fi
    !
    !      Hussein, T., Smolik, J., Kerminen, V.-M., Kulmala, M., 2012.
    !      Modelling dry deposition of aerosol particles onto rough surfaces,
    !      Aerosol Science and Technology, 46, 44-59, 
    !      doi: 10.1080/02786826.2011.605814.
    !-------------------------------------------------------------------------

    implicit none

    INTEGER, intent(in)                    :: IMAX
    REAL( dp), intent(in)                  :: temp,press,mbh
  !  REAL( dp), intent(in)                  :: ustar
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPA
    REAL( dp), DIMENSION(MMAX,IMAX), intent(in)       :: DENSPAR
  !  REAL( dp), intent(in)                             :: Fplus

    REAL( dp), DIMENSION(MMAX,IMAX), intent(out)      :: depo  

    real( dp)                              :: DENSPA   
    real( dp)                              :: miu,lam,niu,roh
    real( dp)                              :: d1,d2,f1,f2
    real( dp)                              :: DRP,Ypluscbl,Yplusnull
    real( dp)                              :: CC,DIFFCO,tauP,tauPplus
    real( dp)                              :: CD,cdrex,rex
    real( dp)                              :: xcd,ypol    
    real( dp)                              :: VG,vgplus, Sc
!    real( dp)                              :: splus,Dplus
    real( dp)                              :: tauL
    real( dp)                              :: niutplus,vy2plus
    real( dp)                              :: pp,Gint,dyplus
    
    ! Number of points for polynomial fitting
    integer, parameter                     :: nre=25
    ! specific gas constant for dry air [J/(kg·K)]
    real( dp), parameter                   :: R = 287.058
    ! Ymax: boundary condition is sometimes 30 or 200 according to literature!
    real( dp), parameter                   :: Ymax = 100.0
    ! 600 - 1200 points for diamters < 5 um. 12000 for bigger particles!
    ! we use number of points: 450
    integer, parameter                     :: npoints = 450
    
    real( dp), dimension(nre)              :: RE,rlin,CDRE,xa,ya
    real( dp), dimension(npoints)          :: Yplus, ylin
    real( dp), dimension(npoints)          :: niutplusi,tauLi,dyplusi
    real( dp), dimension(npoints)          :: xi, dxi    
    real( dp), dimension(npoints)          :: vy2plusi,vpy2plusi !,Dplusi 
    real( dp), dimension(npoints)          :: Fdep,Gdep

    integer                                :: M,I,n,j,k

      !Fplus = 0.0  ! perfectly smooth surface, this is INPUT
      !Fplus = 0.5  ! water
      !Fplus = 0.55 ! gravel
      !Fplus = 1.60 ! artifical grass
      !Fplus = 0.2  ! ocean
      !DENSPA= 1500.  !kg/m^3  for all particles for TEST
      DENSPA= 1400.  !kg/m^3  for all particles for ocean TEST

      roh   = press/R/temp                ! air density [kg/m^3] 1.2928
      miu   = miuair(temp)                ! dynamic viscosity of air [kg/ms]
      niu   = miu / roh                   ! kinematic viscosity of air [m2/s]
      lam   = lambda(temp,press)          ! mean free path of air [m]
      
   
      !-- Creating the Y+ grid:
      d1   = LOG10(1.e-9 * ustar / niu)
      d2   = LOG10(Ymax)
      do n=2,npoints-1
          ylin(n) = d1 + (d2/npoints)*n - (d1/npoints)*n 
      end do
      ylin(1)       = d1
      ylin(npoints) = d2
      Yplus = 10**ylin       

      !-- Creating the Re-CDRe grid:
      f1   = -16.0
      f2   = 3.0
      do n=2,nre-1
          rlin(n)  = f1 + (f2/nre)*n - (f1/nre)*n 
      end do
      rlin(1)    = f1
      rlin(nre) = f2
      do j=1,nre
        RE(j)   = 10**rlin(j)     
        CDRE(j) = 0.0
      end do

      do j=1,nre
         if ( RE(j) .le. 0.1 ) then
           CD      = 24./RE(j)
           CDRE(j) = CD*RE(j)*RE(j)
         else if ( RE(j) .le. 2.0 ) then  
           CD      = (24./RE(j)) * ( 1. + (3./16.)*RE(j) + &
                     (9./160.)*RE(j)*RE(j) * LOG(2*RE(j)) )
           CDRE(j) = CD*RE(j)*RE(j)
         else if ( RE(j) .le. 500. ) then
           CD      = (24./RE(j)) * (1. + 0.15*(RE(j)**0.687))
           CDRE(j) = CD*RE(j)*RE(j)
         else if ( RE(j) .le. 2.e5) then
           CD      = 0.44
           CDRE(j) = CD*RE(j)*RE(j)
         endif                             
      end do

      ! nre order polynomial curve fit CDRE-RE
      do j=1,nre
          xa(j) = LOG(CDRE(j))
          ya(j) = LOG(RE(j))
          !write(6,*) 're',j,RE(j),CDRE(j),LOG(RE(j)),LOG(CDRE(j))  
      end do

      !% Boundary limits and surface orientation parameters
      !    i_surface  = [-1 0 1];           1ceiling 2walls 3ground        ![ceiling wall floor]
      ! BC at the top of the boundary
      Ypluscbl = Ymax
      ! Air wall normal fluctuating velocity intensity
      ! adopted from Guha (1997) after Zhao and Wu (2006a) 
      ! This is needed from the calculations of the turbophoresis
      ! process to the deposition velocity of indoor aerosol particles
      vy2plus  =  (0.005* Ypluscbl*Ypluscbl ) /   &
                  (1. + 0.002923 * ( Ypluscbl**(2.128)))
      vy2plus  =  vy2plus**2            
      ! dimensionless air turbulent viscosity niutplus
      ! i.e the air fluid turbulent viscosity to the kinematic 
      ! viscosity of air: niuT / niu   by Johansen (1991)
      niutplus = 0.4 * Ypluscbl
      ! Lagrangian time-scale of the fluid      
      tauL     = ( niutplus * niu )/( vy2plus * ustar*ustar )
      
      ! Particle Size Dependant Parameters
      do M=1,MMAX
       do I=1,IMAX
          DRP=DPA(M,I)
          ! BC at the rough surface
          Yplusnull = (DRP/2.) * (ustar/niu) + Fplus
          ! Cunningham slip correction factor
          ! for spherical particles according to Allen and Raabe (1982) for oil droplets
          ! and Allen and Raabe (1982) for solid particles. 2.1%
          CC       = 1. + ( (lam/DRP) * (2.34 + (1.05*EXP((-0.39)* (DRP/lam) ))) )    
          ! particle diffusion coefficient [m^2/s]
          DIFFCO   = (k_B * temp * CC) / (3. * miu * pi * DRP)
          Sc       = miu/DIFFCO
          ! Particle relaxation time 
          !tauP     = CC * DENSPA * DRP**2. /(18.*miu)  !TEST
          tauP     = CC * DENSPAR(M,I) * DRP**2. /(18.*miu)
          tauPplus = (tauP * ustar * ustar) / niu

          !-- Settling Velocity calculation
          ! if DRP < 10 nm we neglect settling
          if (DRP .ge. 1.e-8) then
            !cdrex = ( 4.*(DRP**3) * roh*(DENSPA - roh) *g*CC ) / (3.*miu*miu)   !TEST
            cdrex = ( 4.*(DRP**3) * roh*(DENSPAR(M,I) - roh) *g*CC ) / (3.*miu*miu)
            xcd=LOG(cdrex)
            ! polint returns value ypol of a polynomial of degree nre to the curve CDRE-RE 
            ! evaluated at xcd=log(cdrex)
            call polint(nre, xa, ya, xcd, ypol)
            ! the exponential of ypol is the desired Reynold value
            ! it is approx.: cdre = 24*re  =>  re = cdre/24  (see above)
            rex = EXP(ypol)
          else
            cdrex = 0.0
            xcd = 0.0
            rex = 0.0
          endif
          ! Gravitational settling velocity [m/s]
          VG = (miu / roh / DRP) * rex
          !VGREF     = 4./3.*pi*(0.5*DRP)**3*DENSPA *g*DIFFCO/k_B/temp     !m/s  !usual vg
          vgplus   = VG/ustar
          !write(6,*) 'vgplus',m,i,DRP,rex,VG,vgplus       

          ! Correction according to Liu and Agarwal (1974) 
          ! for the near wall air turbulences (not used)
          !splus    = (0.69 * ustar*ustar * tauP)/niu 
          !Yplusnull= Yplusnull+splus
           
          ! D+ = D/v + Ep/v ... Dimensionless Brownian and turbulent diffusivities
          !Dplus    =   DIFFCO/niu + (tauL/( (tauL+tauP)*niutplus ))
          !Dplus replaced by (1./Sc) + x * niutplus
                
          ! Calculate dyplusi, vy2plusi, niutplusi,tauLi on Yplus space    
          do k=1,npoints-1
             dyplusi(k)   = Yplus(k+1) - Yplus(k)
             vy2plusi(k)  = (0.005* Yplus(k)*Yplus(k) ) /   &
                            (1. + 0.002923 * ( Yplus(k)**(2.128) ))
             vy2plusi(k)  =  vy2plusi(k)**2
             ! niutplus: smooth from LaiNazaroff2000, rough from Johansen1991 
             niutplusi(k) = 1.e-32  !not zero
             if ( (Yplus(k).ge.0.0).and.(Yplus(k).lt.3.) ) then
                niutplusi(k) = (Yplus(k)/11.15)**(3.)
             else if ( (Yplus(k).ge.3.0).and.(Yplus(k).lt.52.108) ) then 
                niutplusi(k) = (Yplus(k)/11.4)**(2.) - 0.049774;
             else if ( (Yplus(k).ge.52.108).and.(Yplus(k).le.200.) ) then                  
                niutplusi(k) = 0.4 * Yplus(k)
             endif 
             tauLi(k)     =  (niutplusi(k)*niu) / (vy2plusi(k)*ustar*ustar)
          ! Turbophoresis term vtplus
          ! Eq. (3): Vt = -tauP * d(Vpy^2)/dy
          ! Vpy^2 = Particle wall normal fluctuating velocity intensity
          !        Vpy^2 = Vy^2 * (1 + tauP/tauL)^(-1)   valid for tauPplus<138
          ! To consider turbophoresis:
          !        Vpy^2 = Vy^2 * (1 + (tauP/tauL)^0.5 + (tauP/tauL)^1.5 )^(-1)
          !
             xi(k)        =  tauLi(k) / ( tauLi(k) + tauP ) 
             vpy2plusi(k) =  vy2plusi(k) * (  1./( 1. + (tauP/tauLi(k))**0.5   +  &
                               (tauP/tauLi(k))**1.5 )  )
             !Dplus replaced by (1./Sc) + x * niutplus
             !Dplusi(k)    =  DIFFCO/niu + (tauLi(k)/( (tauLi(k)+tauP) *           & 
             !                  max(niutplusi(k),1.e-32) ))
          end do
                         
          Fdep(:) = 0.0
          Gdep(:) = 0.0
          Gint    = 0.0
          dyplus  = 0.0
       
          ! Following is evaluation for the minor intergal F: pp is the integral
          ! For all Yplus greater equal Yplusnull
          do k=1,npoints-1 
            if (Yplus(k) .ge. Yplusnull) then 
              ! F_smooth  and F_rough
              !     p = ( i_surface * vg_plus  +  tauP_plus .* dX ./ dY_plus ) ./ ...
              !    ( 1/Sc + x(1:end-1) .* niuT_plus(1:end-1) );         
              !dxi(k)   = (xi(k+1)*vy2plusi(k+1)) - (xi(k)*vy2plusi(k))
              ! include turbophoresis term:
              dxi(k)   = vpy2plusi(k+1) - vpy2plusi(k)
              pp       = ( vgplus + tauPplus * (abs(dxi(k))/dyplusi(k)) )   / &
                   ( (1./Sc) + xi(k) * niutplusi(k) )
              ! Eq.(7a)
              Fdep(k)  = Fdep(k-1) + (pp*dyplusi(k))
              !write(6,*) 'Fdepk',k,vgplus,vtplus,dyplusi(k),pp,Fdep(k)
            endif             
          end do  ! k loop
            
          ! Following is evaluation for the main integral: G is the integrand
          ! until the second last element of Fdep.
          !matlab: G = exp( F - F(end) ) ./ D_plus;
          do k=1,npoints-2
            if (Yplus(k) .ge. Yplusnull) then
              ! G_smooth and G_rough
              !p = exp( F - F(end) ) ./ ( 1/Sc + tauL ./ ( tauL + tauP ) .* niuT_plus );
              Gdep(k) = EXP( Fdep(k) - Fdep(npoints-2) )  /  &
                       ( (1./Sc) + tauLi(k) / (tauLi(k) + tauP) * niutplusi(k) ) 
              Gint    = Gint + Gdep(k)  * dyplusi(k)
              !write(6,*) 'Gdepk',k,Gdep(k),Fdep(k),Gint
            endif
          end do
          depo(M,I) = 1./  Gint                  ! [-]
          depo(M,I) = depo(M,I)
          !write(6,*) tauPplus,depo(M,I)
          depo(M,I) = depo(M,I) * ustar         ! [m/s]
          !write(6,*) 'depo',m,i,DRP,depo(M,I)
          !write(6,*) DRP,depo(M,I)
          depo(M,I) = depo(M,I)/mbh             ! [1/s]

        end do
      end do !  end of M,I loop
      
     ! stop

    end subroutine depofrough


subroutine settling(DPAW,pres,temp,DENSPAR,IMAX,vterm  )
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
    !      calculates terminal velocity for settling
    !      of particles. Considers drag force at higher
    !      Reynold numbers, for the larger droplets.
    !
    !      interface
    !      ---------
    !
    !        input:
    !           press    [Pa]
    !           temp     [K]
    !           DPAW     [m]
    !           DENSPAR  [kg/m^3]
    !        output:
    !           vterm    [m/s]
    !
    !      method
    !      ------
    !      Creating the Re-CDRe grid is taken from the
    !      routine DEPOFROUGH
    !      For large particle the formulation of
    !      Beard and Pruppacher (1969) is used
    !      
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      J.H. Seinfeld and S.N. Pandis,
    !        Atmospheric Chemistry and Physics, From Air
    !        Pollution to Climate Change, 2nd Edition,
    !        John Wiley & Sons, Inc., Hoboken New Jersey, 2006,
    !        Pages 406-411.
    !      K.V. Beard and H.R. Pruppacher,
    !        A determination of the terminal velocity and
    !        drag of small water drops by means of a wind
    !        tunnel, J. Atmos. Sci., 26, 1066-1072, 1996.
    !
    !------------------------------------------------------------------

    implicit none

    ! input
    integer, intent(in)                              :: IMAX
    real( dp), intent(in)                            :: temp     ! [K]
    real( dp), intent(in)                            :: pres     ! [Pa]
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPAW     ! [m]
    real( dp), dimension(MMAX,IMAX), intent(in)      :: DENSPAR  ! [kg/m^3]

    ! output
    real( dp), dimension(MMAX,IMAX), intent(out)     :: vterm    ! [m/s]

    ! local
    ! specific gas constant for dry air [J/(kgK)]
    real( dp), parameter                             :: R=287.058   ! [J/kg/K]

    ! Number of points for polynomial fitting
    integer, parameter                               :: nre=25
    real( dp), dimension(nre)                        :: rlin
    real( dp), dimension(nre)                        :: RE
    real( dp), dimension(nre)                        :: CDRE
    real( dp), dimension(nre)                        :: xa,ya
    real( dp)                                        :: f1,f2
    real( dp)                                        :: CD
    real( dp)                                        :: cdrex,rex
    real( dp)                                        :: xcd,ypol 

    real( dp)                                        :: rho_air
    real( dp)                                        :: DRP
    real( dp)                                        :: LAMB
    real( dp)                                        :: MYY
    real( dp)                                        :: KN
    real( dp)                                        :: CC
    integer                                          :: I,M
    integer                                          :: j,n


      ! air density [kg/m^3] 1.2928
       rho_air = pres / R / temp
       
       ! free mean path [m]
       LAMB= (6.73e-8_dp*temp*(1.+(110.4_dp/temp)))/ &
          (296._dp*pres*1.373_dp)

       ! dynamic viscosity of air [kg/m/s]
       MYY= (1.832e-5_dp*(temp**(1.5_dp))*406.4_dp)/ &
          (5093._dp*(temp+110.4_dp))

       !-- Creating the Re-CDRe grid:
       f1   = -16.0
       f2   = 3.0
       do n=2,nre-1
           rlin(n)  = f1 + (f2/nre)*n - (f1/nre)*n
       end do
       rlin(1)   = f1
       rlin(nre) = f2
       do j=1,nre
         RE(j)   = 10**rlin(j)     
         CDRE(j) = 0.0
       end do

       do j=1,nre
          if ( RE(j) .le. 0.1 ) then
            CD      = 24./RE(j)
            CDRE(j) = CD*RE(j)*RE(j)
          else if ( RE(j) .le. 2.0 ) then  
            CD      = (24./RE(j)) * ( 1. + (3./16.)*RE(j) + &
                      (9./160.)*RE(j)*RE(j) * LOG(2*RE(j)) )
            CDRE(j) = CD*RE(j)*RE(j)
          else if ( RE(j) .le. 500. ) then
            CD      = (24./RE(j)) * (1. + 0.15*(RE(j)**0.687))
            CDRE(j) = CD*RE(j)*RE(j)
          else if ( RE(j) .le. 2.e5) then
            CD      = 0.44
            CDRE(j) = CD*RE(j)*RE(j)
          endif
       end do

       ! nre order polynomial curve fit CDRE-RE
       do j=1,nre
           xa(j) = LOG(CDRE(j))
           ya(j) = LOG(RE(j))
           !write(6,*) 're',j,RE(j),CDRE(j),LOG(RE(j)),LOG(CDRE(j))  
       end do


       ! Calculate settling velocity [m/s]
       do M=1,MMAX                
         do I=1,IMAX

         ! Droplet diameter
           DRP=DPAW(M,I)

         ! Knudsen number
           KN= 2._dp*LAMB/DRP

         ! Cunningham-Milliken slip correction factor
         ! as in Fitzgerald et al., 1998
           CC = 1._dp + KN * ( 1.37_dp + 0.4_dp * EXP( -1.1_dp / KN ) )

         !-- Settling Velocity calculation

           if (DRP .lt. 5.e-8) then
         ! if DRP < 50 nm we neglect settling

             cdrex = 0.0
             rex   = 0.0
             vterm(M,I) = 0.0

           elseif (DRP .lt. 1.e-5) then
         ! if 50 nm < DRP < 10 um
         ! Seinfeld and Pandis (2006), EQ 9.42

             cdrex = 0.0
             rex   = 0.0
             vterm(M,I) = (1._dp/18._dp)                      *  &
                          (DRP**2 *DENSPAR(M,I) * g*CC)
             vterm(M,I) = vterm(M,I) / MYY 

           else
         ! for large particles, the formulation of 
         ! Beard and Pruppacher (1969) is used.
         ! Reynold number based on CDRE-RE grid

             cdrex = ( 4.*(DRP**3) *rho_air *                    &
                     (DENSPAR(M,I) - rho_air) *g*CC )         /  &  
                     (3._dp*MYY*MYY)

             xcd=LOG(cdrex)

             ! polint returns value ypol of a polynomial of 
             ! degree nre to the curve CDRE-RE 
             ! evaluated at xcd=log(cdrex)
             call polint(nre, xa, ya, xcd, ypol)

             ! the exponential of ypol is the desired Reynold value
             ! it is approx.: cdre = 24*re  =>  re = cdre/24 
             ! (see above)
             rex = EXP(ypol)

             ! Beard and Pruppacher (1969)
             vterm(M,I) = (MYY * rex)/(rho_air*DRP)

           endif

      
           !write(6,*) 'vterm',M,I,DRP,cdrex,rex,vterm(M,I)

         enddo
       enddo


  end subroutine settling


  subroutine wetscavsize(IMAX,DPAW,pres,temp,DENSPAR,rain,wetscav )
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
    !      calculates size-dependent precipitation scavenging
    !      of particles based on raindrop-aerosol collection
    !      efficiency E. For a fix raindrop diameter of 1 mm.
    !      Based on the correlation for E that fits experiments
    !      presented by Slinn, 1983.
    !      Average MBL volume cloud fraction: 5-10%
    !
    !      interface
    !      ---------
    !
    !        input:
    !           press    [Pa]
    !           temp     [K]
    !           DPAW     [m]
    !           DENSPAR  [kg/m^3]
    !           rain     [mm/h]
    !        output:
    !           wetscav  [1/s]
    !
    !      method
    !      ------
    !      J.H. Seinfeld and S.N. Pandis,
    !        Atmospheric Chemistry and Physics, From Air
    !        Pollution to Climate Change, 2nd Edition,
    !        John Wiley & Sons, Inc., Hoboken New Jersey, 2006.
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      Slinn, W.G.N., 1983.
    !      Precipitation scavenging, In: Atmospheric Sciences and
    !      Power Production - 1979, Chap. 11 Divsion of Biomedical
    !      Environmental Research, U.S. Department of Energy, 
    !      Washington, DC.
    !
    !------------------------------------------------------------------

    implicit none

    ! input
    integer, intent(in)                              :: IMAX
    real( dp), intent(in)                            :: temp     ! [K]
    real( dp), intent(in)                            :: pres     ! [Pa]
    real( dp), intent(in)                            :: rain     ! [mm/h]
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPAW     ! [m]
    real( dp), dimension(MMAX,IMAX), intent(in)      :: DENSPAR  ! [kg/m^3]

    real( dp), dimension(MMAX,IMAX), intent(out)     :: wetscav  ! [1/s]

    real( dp), dimension(MMAX,IMAX)                  :: vterm    ! [m/s]
    real( dp)                                        :: rainfall
    real( dp)                                        :: rho_air
    real( dp)                                        :: miu_air
    real( dp)                                        :: niu_air
    real( dp)                                        :: miu_water
    real( dp)                                        :: LAMBDA
    real( dp)                                        :: CC
    real( dp)                                        :: KN
    real( dp)                                        :: DIFFCO
    real( dp)                                        :: Sc
    real( dp)                                        :: Re
    real( dp)                                        :: St
    real( dp)                                        :: taup
    real( dp)                                        :: omega
    real( dp)                                        :: sstar
    real( dp)                                        :: dpaero
    real( dp)                                        :: dpratio
    real( dp)                                        :: EC


    ! specific gas constant for dry air [J/(kgK)]
    real( dp), parameter                             :: R      = 287.058     ! [J/kg/K]

    ! velocity of rain drop [m/s] for drop diameter 0.5 mm
    ! as in Table 20.1 in Seinfeld and Pandis, 1998
    real( dp), parameter                             :: dprain = 0.5*1.e-3   ! [m]
    real( dp), parameter                             :: urain  = 1.59        ! [m/s] @ 0.5 mm
    ! velocity of rain drop [m/s] for drop diameter 0.1 mm
    ! scavenging coefficient is larger for smaller rain drops
   ! real( dp), parameter                             :: dprain = 0.1*1.e-3   ! [m]
   ! real( dp), parameter                             :: urain  = 0.26        ! [m/s] @ 0.1 mm

    ! average cloud volume for precipation
    real( dp), parameter                             :: FC     = 0.10

    integer  :: I,M

      ! rainfall intensity [m/s]
      rainfall=rain*1.e-3/3600.
      !rainfall= 1.0*1.e-3/3600. !test 1 mm/h

      ! free mean path [m]
      LAMBDA= (6.73e-8_dp*temp*(1.+(110.4_dp/temp)))/ &
              (296._dp*pres*1.373_dp)
      ! air density [kg/m^3] 1.2928
      rho_air   = pres/R/temp
      ! dynamic viscosity of air [kg/m/s]
      miu_air   = miuair(temp)
      ! kinetic viscosity of air [m^2/s]
      niu_air   = miu_air / rho_air
      ! dynamic viscocity of water [kg/m/s]
      ! Atkins Physical Chemistry p. 761
      miu_water = 0.891e-3_dp

      ! Reynold number of raindrop based on its radius [-]
      ! 22.386 @ 0.5 mm
      Re        = (dprain * urain * rho_air) / (2.*miu_air)

      ! S*
      sstar   = 1.2_dp + (1./12.)*log(1._dp + Re)
      sstar   = sstar / (1._dp + log(1._dp + Re))

      ! Viscocity ratio
      omega     = miu_water / miu_air

      ! get terminal velocity [m/s]
      CALL settling(DPAW,pres,temp,DENSPAR,IMAX,vterm)


      ! Particle Size Dependent Parameters
      ! Collection Efficiency Ec
      do M=1,MMAX
        do I=1,IMAX

          dpaero = DPAW(M,I)

         ! Ratio of particle-to-rain drop
          dpratio = dpaero / dprain

         ! Knudsen number [m/m]
          KN = 2._dp*LAMBDA / dpaero

         ! Slip correction factor
          CC = 1._dp + (KN*(1.142_dp+(0.558_dp*EXP((-.999_dp)/KN))))

         ! particle diffusion coefficient [m^2/s]
          DIFFCO = (k_B * temp * CC) / (3._dp * miu_air * pi * dpaero)
          
         ! Schmidt number of collected particle [-]
          Sc     = niu_air/DIFFCO

         ! Relaxation time of particles [s]
         taup    = CC * DENSPAR(M,I) * dpaero**2. / (18._dp * miu_air)

         ! Stoke number of collected particle [-]
         St      = (2*taup * (urain - vterm(M,I))) / dprain

         ! Collision Efficiency [-]
         EC      =    &
              ! Brownian Diffusion
                   (4._dp/(Re*Sc)) * ( 1._dp + 0.4_dp*Re**(1./2.)*Sc**(1./3.)   +  &
                   0.16_dp * Re**(1./2.) * Sc**(1./2.) )                        +  &
              ! Interception
                   4.*dpratio*(1._dp/omega + (1._dp + 2.*Re**(1./2.))*dpratio)        
         if (St>sstar) then
              ! Impaction
            EC   = EC + ( (St - sstar) / (St - sstar + (2./3.)) )**(3./2.)      *  &
                       ! scale by squared density ratio water/particle
                        (1000._dp/DENSPAR(M,I))**(1./2.)
         endif

         ! Scavenging coefficient [1/s]
         ! For monodisperse aerosols and raindrops that have same diameter dprain
         wetscav(M,I) = FC * (3./2.) * EC * rainfall / dprain

         ! if (rainfall>0.) print *,M,I,dpaero,EC,wetscav(M,I)*3600.

        enddo
      enddo


  end subroutine wetscavsize


  ELEMENTAL REAL(dp) function wetscavbulk(rain)
    !------------------------------------------------------------------
    !
    ! Wet scavenging of particles in s^-1
    !   
    ! Follows the bulk Precipitation Scavenging for particles
    ! as described in Fitzgerald et al., 1998
    ! We consider only in-cloud scavenging of aerosols in those
    ! sections that contain particles large enough to be activated
    ! during cloud formation.
    ! Thus nucleation mode aerosols are not scavenged.
    ! input: rainfall rate in mm hr^-1
    ! output: scavenging efficiency in s^-1
    !     average MBL volume cloud fraction: 5-10%
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      reference
    !      ---------
    !      Fitzgerald, J.W., Hoppel, W.A., Gelbard, F., 1998.
    !      A one-dimensional sectional model to simulate multicomponent
    !      aerosol dynamics in the marine boundary layer.
    !      1. Model description,
    !      Journal of Geophysical Research, 103, D13, 16085-16102.
    !
    !------------------------------------------------------------------
    !
      implicit none
        REAL( dp), intent(in)     :: rain
        REAL( dp), parameter      :: FC=0.10


        wetscavbulk=FC*3.49E-04*rain**(0.79)     


    END FUNCTION wetscavbulk


subroutine polint(PMAX, xa,ya,x,y)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      NUMERICAL RECIPES
    !
    !      purpose
    !      -------
    !      Polynomial Interpolation and Extrapolation
    !
    !      interface
    !      ---------
    !
    !
    !      method
    !      ------
    !      NUMERICAL RECIPES
    !      Given arrays xa and ya, each of length n, and given a value x, 
    !      this routine returns a value y, and an error estimate dy. 
    !      If P(x) is the polynomial of degree N-1 such that
    !      P(xai) = yai, i = 1, ... , n, then the returned value y = P(x).
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      Book Numerical Recipes, Ch. 3.1
    !
    !
    !------------------------------------------------------------------
    implicit none

    REAL( dp), intent(in)                     :: x 
    integer, intent(in)                       :: PMAX

    REAL( dp), dimension(PMAX), intent(in)    :: xa,ya
    REAL( dp), intent(out)                    :: y

    REAL( dp)                                 :: den,dif,dift,ho,hp,w,dy
    
    real( dp), dimension(PMAX)                :: c,d
    integer                                   :: i,m,ns

          ns=1
          dif=abs(x-xa(1))
          
          ! Here we find the index ns of the closest table entry
          do i=1,PMAX
            dift=abs(x-xa(i))
            if (dift .lt. dif) then
              ns=i
              dif=dift
            endif
            ! initialize the tableau of c's and d's.
            c(i)=ya(i) 
            d(i)=ya(i)
          end do
          
          ! This is the initial approximation to y.
          y=ya(ns)
          ns=ns-1
          
          ! For each column of the tableau
          ! we loop over the current c's and d's and update them.
          ! After each column in the tableau is completed, we decide
          ! which correction, c or d, we want to add to our accumulating
          ! value of y, i.e., which path to take through
          ! the tableau forking up or down. We do this in such a
          ! way as to take the most "straight line" route through the
          ! tableau to its apex, updating ns accordingly to keep track
          ! of where we are. This route keeps the partial approximations
          ! centered (insofar as possible) on the target x.
          ! The last dy added is thus the error indication.
          do m=1,PMAX-1 
            do i=1,PMAX-m 
              ho=xa(i)-x
              hp=xa(i+m)-x
              w=c(i+1)-d(i)
              den=ho-hp
              ! if(den.eq.0.)pause 'failure in polint'
              ! This error can occur only if two input xa's are (to within roundoff) identical.
              if (den .eq. 0.0) then
                den= den+1.e-32
              endif
              den=w/den
              d(i)=hp*den 
              c(i)=ho*den
            end do
            if (2*ns .lt. PMAX-m) then
              dy=c(ns+1)
            else
              dy=d(ns)
              ns=ns-1
            endif
            y=y+dy
          end do

  end subroutine polint


  ELEMENTAL REAL(dp) function miuair(temp)
    !------------------------------------------------------------------
    !
    ! This routine calculates the dynamic viscosity of air 
    ! (known as mu or eeta)
    !   
    ! Tareq Hussein
    ! October 12, 2012. Amman-Jordan
    !
    !      author
    !      -------
    !      Tareq Hussein
    !      Professor
    !      University of Jordan
    !      School of Science
    !      Department of Physics
    !      Amman, 11942 Jordan
    !      Mobile: +962 779 483608
    !      Tel:    +962 6 5355000, ext: 22060
    !      Fax:    +962 6 5300253
    !      e-mail: t.hussein@ju.edu.jo
    !      Currently:
    !      Visiting Professor
    !      University of Helsinki
    !      Institute for Atmospheric and Earth System Research (INAR)
    !      PL 64, FI-00014 UHEL
    !      Helsinki, Finland
    !      Mobile: +358 400 867890
    !
    ! input:  temp::   temperature [K]
    ! output: miu_air:: dynamic viscosity [Pa s]
    !      
    !------------------------------------------------------------------
    !  
      implicit none
        REAL( dp), intent(in)     :: temp
        ! reference viscosity in [Pa s] at T0
        REAL( dp), parameter      :: miu0 = 1.827e-5
        ! reference temperature [K]
        REAL( dp), parameter      :: T0  = 291.15
        ! Sutherland's constant [K]
        REAL( dp), parameter      :: C = 120
 
    
        miuair = miu0 * (T0 + C) 
        miuair = miuair / ( (temp + C) * (temp/T0)**(1.5) )    

    END FUNCTION miuair
    
  ELEMENTAL REAL(dp) function lambda(temp,press)
    !------------------------------------------------------------------
    !
    ! This routine calculates the mean free path of air molecules.
    ! (known as mu or eeta)
    !   
    ! Tareq Hussein
    ! October 12, 2012. Amman-Jordan
    !
    !      author
    !      -------
    !      Tareq Hussein
    !      Professor
    !      University of Jordan
    !      School of Science
    !      Department of Physics
    !      Amman, 11942 Jordan
    !      Mobile: +962 779 483608
    !      Tel:    +962 6 5355000, ext: 22060
    !      Fax:    +962 6 5300253
    !      e-mail: t.hussein@ju.edu.jo
    !      Currently:
    !      Visiting Professor
    !      University of Helsinki
    !      Institute for Atmospheric and Earth System Research (INAR)
    !      PL 64, FI-00014 UHEL
    !      Helsinki, Finland
    !      Mobile: +358 400 867890
    !
    ! input:  temp::   temperature [K]
    !         press::  air pressure [Pa]
    ! output: lambda == mean free path [m]
    !      
    !------------------------------------------------------------------
    !  
      implicit none
        REAL( dp), intent(in)     :: temp,press
        ! molecular diameter of an air molecule [m]
        REAL( dp), parameter      :: dm = 0.00037e-6
        !density of air molecules [1/m3]
        REAL( dp)                 :: nair  
        
        nair = press / k_B / temp 
        lambda = 1 / ( sqrt(2.0) * nair * pi * dm * dm );    

    END FUNCTION lambda


end module gde_deposition
