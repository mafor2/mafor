! <gde_toolbox.f90 - A component of the Multicomponent
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
!*    except:
!*    function  SUNDIS             written by Sasha Madronich
!*    functions ROOLQ and ROOLQD   written by Hanna Vehkamaeki
!*    functions SURTEN and WATERPS written by Hanna Vehkamaeki
!*    function  ACIDPS             written by Liisa Pirjola
!* 
!*****************************************************************************!
module gde_toolbox

 use gde_constants,   only : R_gas,N_A,pi,M_nh3

  implicit none

  private

  public :: k_3rd, k_arr, vmean_func, henry_func
  public :: zfarr
  public :: molec2ug
  public :: molecdiff
  public :: conch2o, k_SIV_H2O2,roolm
  public :: satpress_alkane
  public :: st_org, rho_org, surf_succin
  public :: surten, waterps, acidps
  public :: roolq,roolqd
  public :: satps_sulf
  public :: sulfhydrates
  public :: newton,machineps

! KPP DP - Double precision kind
  INTRINSIC :: SELECTED_REAL_KIND
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)

contains

!---------------------------------------------------------------------------
    ELEMENTAL REAL(dp) function conch2o(temp,RH,press,cair)

    !------------------------------------------------------------------
    implicit none
    REAL( dp) , INTENT(IN)     :: temp,RH,press
    REAL( dp) , INTENT(IN)     :: cair      ! air concentration [molecules/cm3]
    REAL( dp)                  :: psat

    ! saturation partial pressure of gaseous H2O [Pa]
    psat    = 6.982616467E+5_dp &
               - 1.888612677E+4_dp*temp    + 2.132971127E+2_dp*temp**2 &
               - 1.288405237E+0_dp*temp**3 + 4.393186046E-3_dp*temp**4 &
               - 8.023554873E-6_dp*temp**5 + 6.136820929E-9_dp*temp**6
    conch2o    = cair * RH * psat / press

  END FUNCTION conch2o

  ELEMENTAL REAL(dp) function k_SIV_H2O2(k_298,tdep,cHp,temp)
    ! special rate function for S(IV) + H2O2

    REAL( dp), INTENT(IN) :: k_298 ! k at T = 298.15K
    REAL( dp), INTENT(IN) :: tdep  ! temperature dependence
    REAL( dp), INTENT(IN) :: cHp   ! c(H+)
    REAL( dp), INTENT(IN) :: temp  ! temperature

    INTRINSIC EXP

    k_SIV_H2O2 = k_298 * EXP(tdep*(1._dp/temp-3.3540E-3_dp)) &
      * cHp / (cHp+0.1_dp)

  END FUNCTION k_SIV_H2O2

    ELEMENTAL REAL(dp) function k_3rd(temp,cair,k0_300K,n,kinf_300K,m,fc)

    !------------------------------------------------------------------
    !
    !****  Three body reaction rate constant
    !****
    !
    !------------------------------------------------------------------
    implicit none

    REAL( dp) , INTENT(IN)     :: temp      ! temperature [K]
    REAL( dp) , INTENT(IN)     :: cair      ! air concentration [molecules/cm3]
    REAL( dp) ,     INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL( dp) ,     INTENT(IN) :: n         ! exponent for low pressure limit
    REAL( dp) ,     INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL( dp) ,     INTENT(IN) :: m         ! exponent for high pressure limit
    REAL( dp) ,     INTENT(IN) :: fc        ! broadening factor (usually fc=0.6)
    REAL( dp)                  :: zt_help, k0_T, kinf_T, k_ratio

    zt_help = 300._dp/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    k_3rd   = k0_T/(1._dp+k_ratio)*fc**(1._dp/(1._dp+log10(k_ratio)**2))

  END FUNCTION k_3rd

  !---------------------------------------------------------------------------

  ELEMENTAL REAL(dp) function k_arr (k_298,tdep,temp)

    !------------------------------------------------------------------
    !
    !****  Arrhenius function
    !****  for equilibrium constants
    !
    !------------------------------------------------------------------
    implicit none
    REAL( dp) ,     INTENT(IN) :: k_298 ! k at T = 298.15K
    REAL( dp) ,     INTENT(IN) :: tdep  ! temperature dependence
    REAL( dp) ,     INTENT(IN) :: temp  ! temperature

    k_arr = k_298 * exp(tdep*(1./temp-3.3540E-3_dp))  ! 1/298.15=3.3540e-3

  END FUNCTION k_arr



  ELEMENTAL REAL(dp) function zfarr(rx1,er,ztrec)
    !------------------------------------------------------------------
    !
    !****  ZFARR calculation of Arrhenius expression for rate constants
    !
    !------------------------------------------------------------------
    !
    implicit none
    REAL( dp) ,intent(in) :: rx1,er,ztrec
    !
    zfarr=rx1*exp(er*ztrec)
    !
  end function zfarr

  ELEMENTAL REAL(dp) function henry_func(pa, pb,temp)
    !------------------------------------------------------------------
    !
    !****  calculation of T-dependent Henry coefficient
    !
    !------------------------------------------------------------------

      REAL( dp) ,intent(in)  :: pa, pb, temp

      ! 1/298.=3.3557e-3
      henry_func = pa*exp(pb*((1._dp/temp)-3.3557E-3_dp))

    END FUNCTION henry_func


  ELEMENTAL REAL(dp) function vmean_func(molweight,temp)

    !------------------------------------------------------------------
    !
    ! mean molecular speed from Maxwell-Boltzmann distribution:
    ! vmean=sqrt(8*R_gas*T/(M*pi))      (M in kg/mol)
    ! vmean in m/s
    ! sqrt(8*R_gas/pi)=4.60138
    !
    !------------------------------------------------------------------
    !
    implicit none
      REAL( dp), intent(in)     :: molweight, temp

      vmean_func = sqrt(temp/molweight)*4.60138_dp

    END FUNCTION vmean_func

  ELEMENTAL REAL(dp) function st_org(temp)
    !--------------------------------------------------surf_org(----------------
    !
    !  T-dependent saturation pressure of organic vapour
    !
    !------------------------------------------------------------------
    !
        implicit none
          REAL( dp), intent(in)     :: temp

          st_org = 0.05073_dp - temp*0.08381e-3_dp

    END FUNCTION st_org

  ELEMENTAL REAL(dp) function surf_succin(temp)
    !------------------------------------------------------------------
    !
    !  T-dependent surface tension of organic vapour
    !  represented by succinic acid  
    !  acc. to Hyvarinen et al., J. Chem. Eng. Data,51,255-260, 2006
    !  unit of surface tension: kg/s2 or N/m
    ! 
    !------------------------------------------------------------------
    !
        implicit none
          REAL( dp), intent(in)     :: temp

          surf_succin = 83.45_dp - temp*0.12_dp
          surf_succin = surf_succin*1.e-3_dp

    END FUNCTION surf_succin

  ELEMENTAL REAL(dp) function rho_org(temp)
    !------------------------------------------------------------------
    !
    !  Density of the liquid hexanol - water solution
    !      from "Data for phase transitions in aerosol systems"
    !              section 2.2; value for hexanol is used
    !
    !------------------------------------------------------------------
    !
        implicit none
          REAL( dp), intent(in)     :: temp

          rho_org = 1022.56_dp - temp*0.6947_dp
      

    END FUNCTION rho_org


  ELEMENTAL REAL(dp) function satpress_alkane(carbonn,temp)
    !------------------------------------------------------------------
    !
    !  Vapour pressure of n-alkanes as function of 
    !  number of carbon atoms
    !      from
    !      Lemmon, E.W. and Goodwin, A.R.H.								
    !      Critical properties and vapor pressure equation for 
    !      alkanes CnH2n+2: Normal alkanes with n<=36 and isomers 
    !      for n=4 through n=9.
    !      J. Phys. Chem. Ref. Data, 29(1), 2000								
    !
    !------------------------------------------------------------------
    !
        implicit none

          REAL( dp), intent(in)     :: temp
          INTEGER, intent(in)       :: carbonn
          
          REAL( dp), parameter      :: a1=1200.
          REAL( dp), parameter      :: a2=7.2353461
          REAL( dp), parameter      :: a3=-0.31819703
          REAL( dp), parameter      :: a4=0.43600696
          REAL( dp), parameter      :: a5=-0.26905663
          REAL( dp), parameter      :: b1=0.73499318
          REAL( dp), parameter      :: b2=2.1151684
          REAL( dp), parameter      :: b3=-0.36229342         
          REAL( dp), parameter      :: b4=0.69123121
          REAL( dp), parameter      :: b5=-0.22056059
          REAL( dp), parameter      :: b6=-2.8890416
          REAL( dp), parameter      :: c1=3.0
          REAL( dp), parameter      :: c2=1.1556787
          REAL( dp), parameter      :: c3=-0.053002297
          REAL( dp), parameter      :: c4=0.61691228
          REAL( dp), parameter      :: N1= -5.968127116539
          REAL( dp), parameter      :: N2=  1.244160804719
          REAL( dp), parameter      :: N3= -0.501340221969
          REAL( dp), parameter      :: N4= -1.316444840726
          REAL( dp), parameter      :: N5= -4.589426828138
          REAL( dp), parameter      :: N6=  1.365162911745
          REAL( dp), parameter      :: N7= -5.162331161063
          REAL( dp), parameter      :: N8= -5.678478208052
          REAL( dp), parameter      :: N9= -2.798018429501
          REAL( dp), parameter      :: N10= 3.010247883995
          REAL( dp), parameter      :: N11=-7.634352378867
          REAL( dp), parameter      :: N12=-0.178136926230

          REAL( dp)                 :: TC,PC,TR,tau,omega,PR0,PR1,PR2
          REAL( dp)                 :: lnsatp
          ! odd/even switch?
          INTEGER                   :: EVENODD,ODDEVEN
          
            ! check even/odd carbon number
            IF (mod(carbonn,2) .GT.0) THEN   !odd
              EVENODD=0
              ODDEVEN=1
            ELSE
              EVENODD=1
              ODDEVEN=0
            ENDIF  
                      
            ! Eq. 20 Calculation of critical temperature (K)
            TC=a1-EXP( a2+(a3*carbonn**a4)+ (EVENODD*a5/(carbonn**5)) )  
            ! Eq. 21 Calculation of critical pressure (MPa)
            PC=b1+EXP( b2+(b3*carbonn**b4)+(b5/carbonn)+   &
                 ((ODDEVEN*b6)/((carbonn+1)**4)) )    
            ! Eq.24 Calculation of omega		
            omega=c1-( EXP(c2+(c3*(carbonn**c4))) )    
            ! Eq. 2 Calculation of Tr and tau
            TR=temp/TC     
            tau=1-TR
            ! Eq. 13 Calculation of Parameter Pr0
            PR0=N1*tau+N2*(tau**1.5)+N3*(tau**2.5)+N4*(tau**5)
            ! Eq. 14 Calculation of Parameter Pr1
            PR1=N5*tau+N6*(tau**1.5)+N7*(tau**2.5)+N8*(tau**5)
            ! Eq. 15 Calculation of Parameter Pr2
            PR2=N9*tau+N10*(tau**1.5)+N11*(tau**2.5)+N12*(tau**5)
            
            ! Eq.12 Calculation of n-alkane vapour pressure (Pa)
            lnsatp=(1/TR)*(PR0+omega*PR1+(omega**2)*PR2)
            satpress_alkane=EXP(lnsatp)*PC*1.e6


    END FUNCTION satpress_alkane

  ELEMENTAL REAL(dp) function roolm(x2,temp)

    !------------------------------------------------------------------
    !
    ! Particle density in kg/m^3
    !
    !------------------------------------------------------------------
    !
    implicit none
        REAL( dp), intent(in)     :: x2, temp
        REAL( dp)                 :: ROO0,ROO60,ROO

          IF (x2.LT.0.6) THEN
           ROO0 = 0.99894 + x2*(0.74823 + x2*(-4.07622D-3 + x2*0.31788))
           ROO60 = 0.98299 +x2*(0.60819 + x2*(0.23326 + x2*0.15419))
          ELSE
           ROO0 = 0.47352 + x2*(4.90399 + x2*(-11.91650 + x2*(15.05760  -&
              x2*6.66837)))
          ROO60 =  0.25052 + x2*(5.73314 + x2*(-13.13814 + x2*(15.56578 -&
              x2*6.61870)))
          END IF
          ROO = ROO0 + (ROO60 - ROO0)*(temp - 273.15_dp)/60.0_dp
          roolm = 1000.0_dp*ROO

    END FUNCTION roolm

    ELEMENTAL REAL(dp) function molec2ug(mw)
    !------------------------------------------------------------------
    !
    !****  molec2ug converts molec/cm3 --> ug/m3. 
    !****  1/molec2ug converts ug/m3 --> molec/cm3
    ! input: MW in g/mol
    !
    !------------------------------------------------------------------
    !
    implicit none
    REAL( dp)  , parameter :: avog=6.022e23_dp
    REAL( dp)  , parameter :: AVOug=1.e12_dp/avog
    REAL( dp) , intent(in) :: mw

    molec2ug=mw*AVOug

    end function molec2ug

    ELEMENTAL REAL(dp) function molecdiff(mworg,sigma,temp,pres)
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
    !      calculate the molecular diffusion coefficient in [cm2/s]
    !
    !      interface
    !      ---------
    !
    !          input:
    !             mworg   molecular weight of vapor [g/mol]
    !             sigma   collision diameter        [Angstroem = 1e-10 m]
    !             temp    air temperature           [K]
    !             pres    air pressure              [atm]
    !
    !      method
    !      ------
    !      molecular diffusion coefficient in [cm2/s]
    !      Chapmans-Enskog equation with the first order 
    !      approximation of the collision parameter, Omega = 1
    !      following
    !      Wyslouzil et al.(1991), Equation 21:
    !
    !      D_AB = 0.00266*temp**(3/2) / p*sqrt(M_AB)*sigma_AB**2
    !      with
    !           M_AB = 2*(1/M_A + 1/M_B)**(-1)
    !           sigma_AB = (sigma_A + sigma_B)/2
    !      in the code below DC0 = D_AB
    !      A is index for trace gas, B is index for air
    !      sigma is the collision diameter, estimated from
    !      the molecular volume of the liquid, in Angstroem
    !
    !      for Air:     
    !           M_B = M_air = 28.965 g/mol
    !           sigma_B     = 51.96**(1/3) = 3.711 Angstroem 
    !      Ref: R.C. Reid, J.M. Prausnitz and B.E. Poling,
    !         The Properties of Gases and Liquids,
    !         McGraw-Hill, New York, 1987.
    !      For MSA:
    !           M_A         = 96.11 g/mol
    !           sigma_A     = 40.3**(1/3) = 3.428 Angstroem
    ! 
    !      Value of DC0 should be around 0.04 cm2/s to 0.05 cm2/s.
    !
    !      reference
    !      ---------
    !      Wyslouzil, B.E., Seinfeld, J.H., Flagan, R.C. 
    !      and Okuyama, K.
    !      Binary nucleation in acid-water systems. I. 
    !      Methanesulfonic acid-water, 
    !      J. Phys. Chem, 94(10), 6827-6850, 1991
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    !
    implicit none
 
    real( dp), intent(in)  :: mworg
    real( dp), intent(in)  :: sigma
    real( dp), intent(in)  :: temp
    real( dp), intent(in)  :: pres
    real( dp)              :: DC0


         DC0    = 1.e-3_dp*temp**1.75_dp*(1.0_dp/28.965_dp +        &
                  1.0_dp/mworg)**0.5_dp/pres
      molecdiff = DC0/(sigma + 51.96_dp**(1.0_dp/3.0_dp))**2.0_dp


    end function molecdiff

!------------------------------------------------------------------
!!! NUMERICAL METHOD FUNCTIONS !!!
!------------------------------------------------------------------

    subroutine newton(start, hk, kh, cion, ctot, imax, delta, eps,      &
                      print_out,deriv_zero,newtonx)
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
    !      Newton-Raphson iteration, root finding method
    !
    !      interface
    !      ---------
    !          input:
    !            start    initial x = NVAP(A_NH4)            [molec/m^3]
    !            hk       henry  H'(NH3)                     [-]
    !            kh       dissociation  K'(NH3)              [-]
    !            cion     ion charge concentration           [molec/m^3]
    !            ctot     total concentration                [ng/m^3]
    !
    !      method
    !      ------
    !      Newton-Raphson iteration
    !      Approximate the solution of f(x) = 0
    !      f       ! function
    !      df      ! derivate of function f
    !
    !      reference
    !      ---------
    !      https://www.tek-tips.com/viewthread.cfm?qid=1667584
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    !
    implicit none

    integer, intent(in)    ::  imax
    real( dp), intent(in)  ::  start
    real( dp), intent(in)  ::  ctot
    real( dp),dimension(4,imax),intent(in)  :: hk
    real( dp),dimension(4,imax),intent(in)  :: kh
    real( dp),dimension(4,imax),intent(in)  :: cion
    real( dp), intent(in)  ::  delta, eps
    logical, intent(in)    ::  print_out


    real( dp), intent(out) ::  newtonx
    logical, intent(out)   ::  deriv_zero 

    real( dp) :: fn
    real( dp) :: dfx
    real( dp) :: x
    real( dp) :: x_old 
    real( dp) :: ct
    integer   :: k

    k = 0
    x = start

! convert ng/m^3 to molec/m^3
    ct = (ctot*1.e-3)/molec2ug(M_nh3)
    ct = ct*1.e6

    dfx =  derivnh3(x,hk,kh,cion,imax)
    fn = fn_ammonia(x,hk,kh,cion,ct,imax)

    if (print_out) then
      ! print header and starting values
      write(*,'(A4, A12, A12, A12)') 'k', 'x', 'f(x)', 'df(x)'
      write(*,'(i4, ES12.4, ES12.4, ES12.4)') k, x, fn, dfx
    end if

    if (ABS(fn) < eps) then 
      ! solution at start point - terminate computation
      ! do nothing and return
    else

      do ! repeat-until loop
! derivate is defined in separate function
!        df = deriv(f, x)
        dfx = derivnh3(x,hk,kh,cion,imax)
        if (dfx .eq. 0) then 
          ! derivative is zero - terminate computation
          deriv_zero = .true.
          return
        endif
        ! compute 
        x_old = x
! f(x) is defined in separate function
        fn = fn_ammonia(x,hk,kh,cion,ct,imax)
        x = x - fn/dfx
        k = k + 1         ! iteration number 
        if (print_out) then
          ! print iteratzion result
          write(*,'(i4, ES12.4, ES12.4, ES12.4)') k, x, fn, dfx
        endif
        ! stop after max. 10 iterations
        if (k == 10) then
          exit
        endif
        if ((ABS(x - x_old) < delta) .or.   &
            (ABS(fn) < eps)) then
          ! exit repeat-until loop
          exit
        endif

      end do
    endif

    ! return value
    newtonx = x   
    return

  end subroutine newton


    ELEMENTAL REAL(dp) function machineps(start)
    !------------------------------------------------------------------
    !
    !****  calculate machine epsilon
    !
    !------------------------------------------------------------------
    !
    implicit none

    real( dp), intent(in) :: start ! starting point

    machineps = start
    do while (machineps + 1.0_dp > 1.0_dp)
      machineps = machineps / 2.0_dp
    end do   

  end function machineps 

    REAL(dp) function fn_ammonia(x,HK,KH,cion,ctot,imax)
    !------------------------------------------------------------------
    !
    !****  NH3 equilibration function
    !    M.Z. Jacobson,  AEROSOL SCIENCE TECHNOLOGY, VOL. 39, 
    !    NO. 2, 92-103, DOI: 10.1080/027868290904546, 2005.
    !    Equation (34)
    !    x = C_gas(NH3)
    !    HK = H'  = dim-less Henry-constant NH3
    !    KH = KH' = dim-less K_dissociation NH3_aq
    !
    !------------------------------------------------------------------
    !
    implicit none

    integer, intent(in)      :: imax
    real( dp), intent(in)    :: x
    real( dp), intent(in)    :: ctot
    real( dp),dimension(4,imax),intent(in)  :: HK
    real( dp),dimension(4,imax),intent(in)  :: KH
    real( dp),dimension(4,imax),intent(in)  :: cion

    real( dp), parameter :: sf = 1.0_dp
    real( dp)     :: ci,he,term1
    integer       :: i,m


      fn_ammonia = x
      do m=1,4
        do i=1,imax

        ! cion must not be positive
          ci = min(cion(m,i), 0._dp)

    !   write(*,'(i4, i4, ES12.4, ES12.4, ES12.4)') m, i, HK(m,i), KH(m,i),ci

          he    = HK(m,i) *sf
          term1 = x*he*KH(m,i)
          fn_ammonia = fn_ammonia +                   &
                       x*he       -                   &
                       ( ci*term1 / (term1 + 1._dp) )

        end do
      end do
      fn_ammonia = fn_ammonia - ctot

  end function fn_ammonia


    REAL(dp) function derivnh3(x,HK,KH,cion,imax)
    !------------------------------------------------------------------
    !
    !****  NH3 equilibration function derivative
    !    M.Z. Jacobson,  AEROSOL SCIENCE TECHNOLOGY, VOL. 39, 
    !    NO. 2, 92-103, DOI: 10.1080/027868290904546, 2005.
    !    Equation (35)
    !    x = C_gas(NH3)
    !    HK = H'  = dim-less Henry-constant NH3
    !    KH = KH' = dim-less K_dissociation NH3_aq
    !
    !------------------------------------------------------------------
    !
    implicit none

    integer, intent(in)    :: imax
    real( dp), intent(in)  :: x
    real( dp),dimension(4,imax),intent(in)  :: HK
    real( dp),dimension(4,imax),intent(in)  :: KH
    real( dp),dimension(4,imax),intent(in)  :: cion

    real( dp), parameter :: sf = 1.0_dp
    real( dp)     :: ci,he,term2
    integer       :: i,m


      derivnh3 = 1.00_dp
      do m=1,4
        do i=1,imax

        ! cion must not be positive
          ci = min(cion(m,i), 0._dp)
          he       = HK(m,i) *sf
          term2    = x*he*KH(m,i) + 1._dp
          derivnh3 = derivnh3 +  he             -    &
                     ( ci*he*KH(m,i) / term2 )  +    &
                     ( (ci*x*(he*KH(m,i))**2)   /    &
                       term2**2 )

        end do
      end do

  end function derivnh3


!------------------------------------------------------------------
!!! NUCLEATION RELATED FUNCTIONS !!!
!------------------------------------------------------------------

  subroutine sulfhydrates(t,pw,hydpw,kprod,rhoh,rhohm)
    !*******************************************************************
    !  Routine for evaluating hydrate coefficients and
    !   hydrates in gas phase
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      reference
    !      ---------
    !      Jaecker-Voirol et al JCP 87 4849-4852, 1987
    !      Wyslouzil et al., 1991
    !*******************************************************************
    implicit none

        REAL( dp), intent(in)                :: t,pw
        REAL( dp), intent(out)               :: hydpw
        REAL( dp), dimension(5), intent(out) :: kprod,rhoh,rhohm
        integer,parameter:: hydimax=5
        REAL( dp), dimension(hydimax)        :: kh,khm,kprodm
        REAL( dp)                            :: pstand
        REAL( dp)                            :: hydpwm
        integer                              :: i1,i2

       pstand=1.01325E5_dp

! Zeleznik equilibrium constant for reactions
       kh(1)= EXP(6159.1796_dp/t-14.3343_dp) ! water+h2so4 =1-hydrate
       kh(2)= EXP(5813.8776_dp/t-15.5140_dp) ! i-hydrate+water= i+1-hydrate
       kh(3)= EXP(5246.0_dp/t-14.90_dp)                ! i-hydrate+water= i+1-hydrate
       kh(4)= EXP(4934.0_dp/t-14.47_dp)                ! i-hydrate+water= i+1-hydrate
       kh(5)= EXP(4763.0_dp/t-14.21_dp)                ! i-hydrate+water= i+1-hydrate
! same for MSA-water, Wyslouzil et al., 1991
! t dependence missing at the moment.
       khm(1) = EXP(1478.72/t)  ! 142.9 at 298 K water+msa =1-hydrate
       khm(2) = EXP(1025.24/t)  ! 31.2 at 298 K
       khm(3) = EXP(819.0/t)    ! 15.6 at 298 K
       khm(4) = EXP(706.0/t)    ! 10.7 at 298 K
       khm(5) = EXP(648.0/t)    ! 8.80 at 298 K


       DO i1=1,hydimax
          !product of the equilibrium constants k1*k2*k3*..*ki1
                kprod(i1)=1.
                kprodm(i1)=1.
         DO i2 = 1,i1
                kprod(i1)= kh(i2)*kprod(i1)
                kprodm(i1)= khm(i2)*kprodm(i1)
         END DO
       END DO

        ! jaecker-voirol et al JCP 87 4849-4852, 1987:
        !  pa(total)/pa(free)=hydpw=1+k1(pw/p0)+..+k1*k2*..*ki1(pw/p0)^i1+..+
          hydpw = 1.0
          hydpwm=1.0
          DO i1 = 1,hydimax
            hydpw = hydpw + kprod(i1)*(pw/pstand)**i1
            hydpwm = hydpwm + kprodm(i1)*(pw/pstand)**i1
          END DO
        !concentration of the hydrates (divided by rhoa(tot))
        !rho(i1)/rhoa(tot)= k1*k2*..*ki1(pw/p0)^i1/hydpw
          DO i1 = 1,hydimax
            rhoh(i1) = kprod(i1)*(pw/pstand)**i1/hydpw
            rhohm(i1) = kprodm(i1)*(pw/pstand)**i1/hydpwm
          END DO

  end subroutine sulfhydrates


  ELEMENTAL REAL(dp) function surten(xmole,temp)
    !------------------------------------------------------------------
    !  Surface tension of the cluster
    !      xmole = mole fraction of h2so4
    !      temp = temperature [K]
    !      surten = surface tension [N/m]
    !
    !      about 230-323 K , x=0,...,1
    !      (valid down to the solid phase limit temp, depends on molefraction)
    !
    !         author
    !         -------
    !      Hanna Vehkamaeki
    !      Professor
    !      University of Helsinki
    !      Institute for Atmospheric and Earth System Research (INAR)
    !      PL 64, FI-00014 UHEL
    !      Helsinki, Finland
    !      Phone:  +358 50415747
    !      Email:  hann.vehkamaki@helsinki.fi 
    !
    !      reference
    !      ---------
    !      Noppel, M., Vehkamaki, H., Kulmala, M., An improved model
    !      for hydrate formation in sulfuric-acid - water nucleation
    !      J. Chem. Phys., Vol. 116, No. 1, 2002.
    !------------------------------------------------------------------
    !
     implicit none

          REAL( dp), intent(in)     :: xmole, temp
          REAL( dp)                 :: xmass, a,b
          REAL( dp)                 :: ma,mb

          ma=18.09*1.661E-27_dp
          mb=98.08*1.661E-27_dp

          xmass = xmole*mb/((1.0-xmole)*ma + xmole*mb) !mass fraction of h2so4

        ! low temperature surface tension
          a = 0.11864+xmass*(-0.11651+xmass*(0.76852+xmass *           &
          (-2.40909+xmass*(2.95434-xmass*1.25852))))
          b = -0.00015709+xmass*(0.00040102+xmass*(-0.00239950+xmass * &
          (0.007611235+xmass*(-0.00937386+xmass*0.00389722))))
          surten=a+temp*b

    END FUNCTION surten

  ELEMENTAL REAL(dp) function waterps(temp)
    !------------------------------------------------------------------
    !  Saturation vapour pressure of pure water in Pa
    !  Temperature temp in kelvins
    !
    !         author
    !         -------
    !      Hanna Vehkamaeki
    !      Professor
    !      University of Helsinki
    !      Institute for Atmospheric and Earth System Research (INAR)
    !      PL 64, FI-00014 UHEL
    !      Helsinki, Finland
    !      Phone:  +358 50415747
    !      Email:  hann.vehkamaki@helsinki.fi 
    !
    !      reference
    !      ---------
    !      Preining, O., Wagner, P.E., Pohl, F.G., Szymanski, W.
    !      Heterogeneous nucleation and droplet growth, report, 
    !      Inst. of Exp. Phys., Univ. of Vienna, 1981.
    !------------------------------------------------------------------
    !
     implicit none

        REAL( dp), intent(in)     :: temp

        waterps = 77.344913_dp - 7235.4247_dp/temp - 8.2*LOG(temp) + &
                0.0057113_dp*temp  ![Pa]
        waterps = EXP(waterps)

    END FUNCTION waterps

  ELEMENTAL REAL(dp) function acidps(temp)
    !------------------------------------------------------------------
    !  Saturation vapour pressure of pure sulphuric acid in Pa
    !  Temperature temp in kelvins
    !
    !         author
    !         -------
    !      Dr. Liisa Pirjola
    !      Docent
    !      Department of Physics
    !      University of Helsinki
    !      P.O Box 64, FI-00014 Helsinki, Finland
    !      Department of Technology
    !      Metropolia University of Applied Sciences
    !      P.O. Box 4071, FI-01600 Vantaa, Finland
    !
    !      reference
    !      ---------
    !      Kulmala, M., Laaksonen, A., Binary nucleation of water–sulfuric
    !      acid system: comparison of classical theories with different 
    !      H2SO4 saturation vapor pressures
    !      J. Chem. Phys., Vol. 93, No. 1, 696-701, 1990.
    !------------------------------------------------------------------
    !
     implicit none

        REAL( dp), intent(in)     :: temp
        REAL( dp)                 :: lpar,pstand

          pstand = 1.01325E5_dp
          lpar   = -11.695+LOG(pstand)  ! Zeleznik
          acidps = 1/360.15-1.0/temp
          acidps = acidps + 0.38/545.*(1.0+LOG(360.15/temp)-360.15/temp)
          acidps = 10156.0*acidps +lpar
          acidps = EXP(acidps)          ! Pa

    END FUNCTION acidps
    
    ELEMENTAL REAL(dp) function satps_sulf(temp,rhin,csul)
    !------------------------------------------------------------------
    !  Saturation vapour pressure of pure sulphuric acid in Pa
    !  Temperature temp xa  [K]
    !  Relative humidity rh [-]
    !  Concentration of h2so4 vapour csul [1/m^3]
    !  xa = mole fraction of h2so4
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      reference
    !      ---------
    !  P. Bolsaitis and J. F. Elliott			
    !    Thermodynamic Activities and equilibrium partial 			
    !    pressures for aqueous sulfuric acid solutions			
    !  J Phys Chem Eng Data 1990, 35, 69-85.
    !
    !------------------------------------------------------------------
    !
     implicit none

        REAL( dp), intent(in)     :: temp,rhin,csul
        
        REAL( dp)                 :: logps,xa,rhoa,t,rh

          rhoa=csul
          t=temp
          rh=rhin
          if (csul < 1.e10) then
            !print *,'Warning: sulfuric acid concentration < 1e10 1/m3, using 1e10 1/m3'
             rhoa=1.e10
          endif
          if (csul > 1.e17) then
            !print *,'Warning: sulfuric acid concentration > 1e17 1/m3, using 1e17 1/m3'
             rhoa=1.e17
          endif
          if (temp < 190.15) then
            !print *,'Warning: temperature < 190.15 K, using 190.15 K'
             t=190.15
          endif
          if (temp > 300.15) then
             !print *,'Warning: temperature > 300.15 K, using 300.15 K'
              t=300.15
          endif
          if (rhin < 0.0001) then
             !print *,'Warning: saturation ratio of water < 0.0001, using 0.0001'
              rh=0.0001
          endif
          if (rhin > 1.0) then
             !print *,'Warning: saturation ratioof water > 1 using 1'
              rh=1.0
          endif

          xa=  0.7409967177282139 - 0.002663785665140117*t +                   &
            0.002010478847383187*LOG(rh) - 0.0001832894131464668*t*LOG(rh)+    &
            0.001574072538464286*LOG(rh)**2 -                                  &
            0.00001790589121766952*t*LOG(rh)**2 +                              &
            0.0001844027436573778*LOG(rh)**3 -                                 &
            1.503452308794887e-6*t*LOG(rh)**3 -                                &
            0.003499978417957668*LOG(rhoa/1.e6) +                              &
            0.0000504021689382576*t*LOG(rhoa/1.e6)
          logps  = -21.4481_dp + 0.0487_dp*t + 7.3243_dp*log(xa)           &
                   - 0.0135_dp*t * log(xa)
          satps_sulf = 10**(logps)        !bar
          satps_sulf = satps_sulf*1.e5    !Pa

    END FUNCTION satps_sulf

  ELEMENTAL REAL(dp) function roolq(xmole,temp)
    !------------------------------------------------------------------
    !  Density of h2o-h2so4 solution
    !  x = mole fraction of h2so4
    !  temp = temperature [K]
    !  roolq = density [kg/m^3]
    !
    !  about 220-373 K , x=0,...,1
    !  (valid down to the solid phase limit temp, depends on molefraction)
    !
    !         author
    !         -------
    !      Hanna Vehkamaeki
    !      Professor
    !      University of Helsinki
    !      Institute for Atmospheric and Earth System Research (INAR)
    !      PL 64, FI-00014 UHEL
    !      Helsinki, Finland
    !      Phone:  +358 50415747
    !      Email:  hann.vehkamaki@helsinki.fi 
    !
    !      reference
    !      ---------
    !      Noppel, M., Vehkamaki, H., Kulmala, M., An improved model
    !      for hydrate formation in sulfuric-acid - water nucleation
    !      J. Chem. Phys., Vol. 116, No. 1, 2002.
    !------------------------------------------------------------------
    !
     implicit none

          REAL( dp), intent(in)     :: xmole, temp
        REAL( dp)                 :: xmass, a,b,c,ma,mb

        ma=18.09*1.661E-27_dp
          mb=98.08*1.661E-27_dp

          xmass = xmole*mb/((1.0-xmole)*ma + xmole*mb) !mass fraction of h2so4

          a = 0.7681724+xmass*(2.1847140+xmass*(7.1630022+xmass *           &
          (-44.31447+xmass*(88.75606+xmass*(-75.73729+xmass*23.43228)))))
          b = 1.808225e-3+xmass*(-9.294656e-3+xmass*(-0.03742148+xmass *    &
          (0.2565321+xmass*(-0.5362872+xmass *                            &
          (0.4857736-xmass*0.1629592)))))
          c = -3.478524e-6+xmass*(1.335867e-5+xmass*(5.195706e-5 +          &
          xmass*(-3.717636e-4+xmass*(7.990811e-4+xmass *                  &
          (-7.458060e-4+xmass*2.58139e-4)))))
          roolq=a+temp*(b+c*temp) ! g/cm^3
          roolq= roolq*1.0E3_dp   ! kg/m^3

    END FUNCTION roolq

  ELEMENTAL REAL(dp) function roolqd(xmole,temp)
    !------------------------------------------------------------------
    !  Calculates  the derivative of the density with respect to mole fraction
    !  roolqd = drho/dx [kg/m^4]
    !  xmole = mole fraction of h2so4, t in K
    !
    !         author
    !         -------
    !      Hanna Vehkamaeki
    !      Professor
    !      University of Helsinki
    !      Institute for Atmospheric and Earth System Research (INAR)
    !      PL 64, FI-00014 UHEL
    !      Helsinki, Finland
    !      Phone:  +358 50415747
    !      Email:  hann.vehkamaki@helsinki.fi 
    !
    !      reference
    !      ---------
    !      Noppel, M., Vehkamaki, H., Kulmala, M., An improved model
    !      for hydrate formation in sulfuric-acid - water nucleation
    !      J. Chem. Phys., Vol. 116, No. 1, 2002.
    !
    !------------------------------------------------------------------
    !
     implicit none

          REAL( dp), intent(in)     :: xmole, temp
        REAL( dp)                 :: roop,room, eps

        room=0.0
        roop=0.0
          eps=0.00001
          IF (xmole .LE. 1.0-eps .AND. xmole .GE. 0.0+eps) THEN
         roop=roolq(xmole+eps,temp)
         room=roolq(xmole-eps,temp)
          ELSE IF (xmole .GT. 1.0-eps ) THEN
         roop=roolq(xmole,temp)
         room=roolq(xmole-2*eps,temp)
          ELSE IF (xmole .LT. 0.0+eps ) THEN
         roop=roolq(xmole+2.*eps,temp)
         room=roolq(xmole,temp)
          END IF
          roolqd=(roop-room)/(2.*eps)

    END FUNCTION roolqd


!------------------------------------------------------------------

end module gde_toolbox
