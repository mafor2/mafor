! <gde_transfer.f90 - A component of the Multicomponent
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
module gde_transfer

  use messy_mecca_kpp_Global

  use gde_input_data,   only       : MMAX,AMAX
  use gde_input_data,   only       : A_OR1,A_OR2,A_OR3,A_OR4,A_OR5
  use gde_input_data,   only       : A_OR6,A_OR7,A_OR8,A_OR9
  use gde_input_data,   only       : A_SUL,A_NH4
  use gde_input_data,   only       : A_WAT
  use gde_input_data,   only       : xmAMIN,fcAMIN,vhAMIN,a_AMIN,b_AMIN
  use gde_input_data,   only       : surf_h2o_std
  use gde_input_data,   only       : pdens

  use gde_toolbox,      only       : surf_succin,surten
  use gde_toolbox,      only       : roolq

  private

  public :: transfer_coeff,aero_init_gasaq
  public :: kelvin_nit
  public :: kelvin_sulf,kelvin_msap
  public :: kelvin_org,kelvin_alkane
  public :: nano_koehler


contains

    !***************************************************************************

    !----------------------------------------------------------------------
    !----------------  TRANSFER TO AQUEOUS PHASE    -----------------------
    !----------------------------------------------------------------------

  subroutine aero_init_gasaq(chem_name,Henry_T0,Henry_Tdep,alpha_T0,alpha_Tdep,molar_mass)
    !----------------------------------------------------------------------
    !     
    !   Initialize the gas-phase / aqueous phase exchange
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      initialize the gas-phase / aqueous phase exchange
    !
    !      interface
    !      ---------
    !      interface to the CAABA/MECCA model
    !
    !      method
    !      ------
    !      subroutine mecca_aero_init_gasaq
    !      in messy_mecca_aero.f90
    !      ChemProp is generated in CAABA/mecca/tracer/chemprop
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

    use messy_main_tracer
    use messy_mecca_kpp_monitor,  only: SPC_NAMES
    use gde_constants,            only: STRLEN_KPPSPECIES
        
    implicit none
    
      REAL( dp),dimension(0:nspec),intent(out)  :: alpha_T0
      REAL( dp),dimension(0:nspec),intent(out)  :: alpha_Tdep
      REAL( dp),dimension(0:nspec),intent(out)  :: Henry_T0
      REAL( dp),dimension(0:nspec),intent(out)  :: Henry_Tdep
      REAL( dp),dimension(0:nspec),intent(out)  :: molar_mass                        
      CHARACTER(STRLEN_KPPSPECIES),dimension(0:nspec),intent(out) :: chem_name

      INTEGER :: icp ! Index ChemProp
      INTEGER :: jn


      do jn = 1,NSPEC
        icp = get_chemprop_index(TRIM(SPC_NAMES(jn)))
        chem_name(jn) = TRIM(SPC_NAMES(jn))
        if (icp /=0) then
          Henry_T0(jn)   = chemprop(icp)%cask_r(R_Henry_T0)
          Henry_Tdep(jn) = chemprop(icp)%cask_r(R_Henry_Tdep)
          alpha_T0(jn)   = chemprop(icp)%cask_r(R_alpha_T0)
          alpha_Tdep(jn) = chemprop(icp)%cask_r(R_alpha_Tdep)
          molar_mass(jn) = chemprop(icp)%cask_r(R_MOLARMASS) / 1000. ! [kg/mol]
        endif
      end do


  end subroutine aero_init_gasaq



  subroutine transfer_coeff(radius,temp,press,xaer,lwc,Henry_T0,Henry_Tdep,   &
                 alpha_T0,alpha_Tdep,molar_mass,k_exf,k_exb)
    !----------------------------------------------------------------------
    !     
    !****  Calculation of gas phase transfer coefficient to the aqueous phase
    !      in [m s^-1]
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculation of transfer coefficient according to Schwartz 1986
    !
    !      interface
    !      ---------
    !
    !        input:
    !          radius       droplet radius                     [m]      
    !          xaer         flag for aq. chem.
    !          lwc          liquid water content
    !          alpha_T0     T-dependence Henry coef.           [-]
    !          alpha_Tdep   T-dependence Henry coef.           [-]
    !          molar_mass   molar mass                         [kg/mol] 
    !         
    !        output:
    !          k_exf        transfer coef. forward reaction 
    !          k_exb        transfer coef. backward reaction
    !
    !      method
    !      ------
    ! Calculation of the transfer coefficient for the mean mode actual
    ! radius
    ! transfer coefficient after Schwarz, 1986 (see Sander & Crutzen '96)
    ! k_mt=1/(r**2/(3*D_gas)+4*r/(3*v_mean*alpha))
    ! if : D_gas=lambda*v_mean/3. then we can rewrite k_mt as:
    ! k_mt=vmean/(r**2/lambda+4*r/(3*alpha))
    !      
    !      reference
    !      ---------
    !      none
    !
    !      modifications
    !      -------------
    !      none
    !      LATER REPLACE RADIUS BY SIZE BIN RADIUS
    !
    !------------------------------------------------------------------
     use gde_constants, ONLY: R_gas,atm2Pa,T0,T0_INV
     
    implicit none

    ! input
      REAL( dp), intent(in)           :: press,temp
 
    ! mecca_mbl:
      REAL( dp), dimension(APN),intent(in)    :: radius    ! [m]      
      REAL( dp), dimension(APN),intent(in)    :: xaer
      REAL( dp), dimension(APN),intent(in)    :: lwc

      REAL( dp),dimension(0:nspec),intent(in) :: alpha_T0
      REAL( dp),dimension(0:nspec),intent(in) :: alpha_Tdep
      REAL( dp),dimension(0:nspec),intent(in) :: Henry_T0
      REAL( dp),dimension(0:nspec),intent(in) :: Henry_Tdep
      REAL( dp),dimension(0:nspec),intent(in) :: molar_mass 
    
    ! output
      REAL( dp),dimension(APN,nspec),intent(out) :: k_exf
      REAL( dp),dimension(APN,nspec),intent(out) :: k_exb

      integer   :: jn,jb
      
      REAL( dp) :: kmtrans
      REAL( dp) :: zlambda
      REAL( dp) :: vmean
      REAL( dp) :: alpha
      REAL( dp) :: henry      

      REAL, PARAMETER :: zrc_t = 1.E-9     ! radius threshold for kmt calc.
       

    ! multiplying factor to make the error of ignoring the need of
    ! integration over the mode minimal with little numerical effort
    ! deduced for box-model study
    !intfac(NU) = 1.0
    !intfac(AI) = 1.0   
    !intfac(AS) = 0.84
    !intfac(CS) = 0.55

    ! zlambda = mean free path
    ! 101325 Pa * 6.6E-8 m / 293.15 K = 2.28E-5 Pa*m/K is from
    ! equation (10-106) in Pruppacher and Klett = ref0064
    !                  zlambda = 2.28E-5 * temp / press
    zlambda = 2.28E-5 * temp / press

    ! alpha_Tdep,henry_To,henry_Tdep,molar_mass: initialize
    !      status = get_gasaq(TRIM(SPC_NAMES(jn)), &
    !    Henry_T0(jn), Henry_Tdep(jn), alpha_T0(jn), alpha_Tdep(jn), &
    !    molar_mass(jn))
    
    ! calculate exchange rate coefficients:
    k_exf(:,:) = 0._dp
    k_exb(:,:) = 0._dp
    ! loop over species
    DO jn=1,nspec
      ! mean molecular speed from Maxwell-Boltzmann distribution:
      ! vmean = SQRT(8*R_gas*T/(M*pi)) = sqrt(4.60138*T/M) [m/s]
      ! M is in kg/mol
      if (molar_mass(jn) > 0.001) then
        vmean = 4.60138 * sqrt(temp/molar_mass(jn))
      else
        vmean = 0.
      endif
      ! calculate accommodation coefficients alpha:
      alpha = 1._dp / (1._dp+(1._dp/alpha_T0(jn)-1._dp) * &
              exp(alpha_Tdep(jn)*((1._dp/T0)-(1._dp/temp))))
      ! calculate temperature dependent Henry's law constants:
      if (henry_T0(jn) > 0.) then
        henry = henry_T0(jn) * exp(henry_Tdep(jn)*((1._dp/temp)-T0_INV))
  ! debug Henry coefficient
  !     write(6,*) 'transfer',jn,vmean,alpha,henry
      else
        henry=0.
      endif
      if (alpha > 0.) then
        DO jb = 1, APN ! loop over modes
          if (radius(jb) >= zrc_t) then
            ! calculate mass transfer coefficients kmt after Schwarz, 1986
            ! k_mt=1/(r**2/(3*D_gas)+4*r/(3*v_mean*alpha))
            ! if D_gas=lambda*v_mean/3. then we can rewrite k_mt as:
            ! k_mt=vmean/(r**2/lambda+4*r/(3*alpha))
            ! 8.11.2011 MSK correct kmtrans (unit: 1/s)
             kmtrans = vmean / ( (radius(jb)*radius(jb)/  &
                        zlambda)+(4.*radius(jb)/(3.*alpha)) )     
            ! reaction rate coefficients:
            k_exf(jb,jn) = xaer(jb) * kmtrans * lwc(jb) ! forward
            if (henry > 0.) then
              k_exb(jb,jn) = xaer(jb) * kmtrans * atm2Pa &
                / (R_gas * 1.E3_dp * temp * henry) ! backward
            endif
            
           ! if (jb==3) write(6,*) 'kexf',jb,radius(jb),kmtrans,lwc(jb),k_exf(jb,jn)

          endif
        END DO
      endif
    END DO    

   ! write(6,*) 'radius', radius(AS), radius(CS)


    ! ------------------------------------------------------------------------
     
  end subroutine transfer_coeff


    !----------------------------------------------------------------------
    !----------------  TRANSFER TO SOLID PARTICLES  -----------------------
    !----------------------------------------------------------------------


subroutine kelvin_sulf(DPA,temp,kelffect,IMAX)
    !----------------------------------------------------------------------
    !
    !****  Kelvin effect for liquid sulfuric acid
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate Kelvin effect of sulfuric acid (dimensionless)
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature              [K]
    !           DPAI     particle diameter            [m]
    !
    !        output:
    !           kelffect Kelvin effect of H2SO4       [-]
    !
    !      method
    !      ------
    !      see Vehkamaki et al., J. Geophys. Res., 107(D22), doi:10.1029/
    !                            2002JD002184, 2002.          
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
    use gde_constants,  only : M_H2O,pi,R_gas,M_H2SO4
    use gde_toolbox,    only : surten,roolq,roolqd

    implicit none

    ! input
    INTEGER, intent(in)                              :: IMAX
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA
    REAL( dp), intent(in)                            :: temp

    ! output
    REAL( dp), dimension(MMAX,IMAX),intent(out)      :: kelffect

    REAL( dp), DIMENSION(MMAX,0:(IMAX+1))  :: RP
    REAL( dp)                              :: xmole,exponente 

    INTEGER                                :: M,I

    ! for Sulfuric Acid
    REAL( dp)                              :: st_h2so4        ! kg/s2
    REAL( dp)                              :: rho_h2so4       ! kg/m3

     xmole=1.0_dp
     ! calculate surface tension (kg/s2)
     st_h2so4=surten(xmole,temp)
     ! calculate liquid density  (kg/m3)
     rho_h2so4=roolq(xmole,temp)
     !write(6,*) 'sulf',st_h2so4, rho_h2so4
 
     do M=1,MMAX
      do I=1,IMAX
        RP(M,I)=DPA(M,I)*0.5_dp  ! [m]
        exponente=2._dp*st_h2so4*M_H2SO4*1.E-3_dp
        exponente=exponente/(R_gas*temp*rho_h2so4*RP(M,I))
        kelffect(M,I)=exp(exponente)
        !write(6,*) M,I,kelffect(M,I)
      end do
     end do

  end subroutine kelvin_sulf

subroutine kelvin_msap(DPA,temp,kelffect,IMAX)
    !----------------------------------------------------------------------
    !
    !****  Kelvin effect for methanesulfonic acid (MSA)
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate Kelvin effect of MSA (dimensionless)
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature              [K]
    !           DPAI     particle diameter            [m]
    !
    !        output:
    !           kelffect Kelvin effect of MSA         [-]
    !
    !      method
    !      ------
    !      see Kreidenweis and Seinfeld, Atm. Environ., 22, 2, pp. 283-296, 1988 
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
    use gde_constants,  only : M_H2O,pi,R_gas

    implicit none

    ! input
    INTEGER, intent(in)                              :: IMAX
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA
    REAL( dp), intent(in)                            :: temp
    ! output
    REAL( dp), dimension(MMAX,IMAX),intent(out)      :: kelffect

    REAL( dp), DIMENSION(MMAX,0:(IMAX+1))  :: RP
    REAL( dp)                              :: exponente 

    INTEGER                                :: M,I

    ! for MSA
    REAL( dp), PARAMETER                   :: M_msa   = 96.11 ! g/mol
    REAL( dp), PARAMETER                   :: st_msa  = 0.053 ! kg/s2
    REAL( dp), PARAMETER                   :: rho_msa = 1507  ! kg/m3 

     do M=1,MMAX
      do I=1,IMAX
        RP(M,I)=DPA(M,I)*0.5_dp  ! [m]
        exponente=2._dp*st_msa*M_msa*1.E-3_dp
        exponente=exponente/(R_gas*temp*rho_msa*RP(M,I))
        kelffect(M,I)=exp(exponente)
      end do
     end do

  end subroutine kelvin_msap

subroutine kelvin_nit(DPA,temp,kelffect,IMAX)
    !----------------------------------------------------------------------
    !
    !****  Kelvin effect for solid ammonium nitrate
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate Kelvin effect of ammonium nitrate (dimensionless)
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature              [K]
    !           DPAI     particle diameter            [m]
    !
    !        output:
    !           kelffect Kelvin effect of NH4NO3      [-]
    !
    !      method
    !      ------
    !      see Stelson and Seinfeld, Atm. Env., 41, S126-S135, 2007  
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
    use gde_constants,  only : M_H2O,pi,R_gas

    implicit none

    ! input
    INTEGER, intent(in)                              :: IMAX
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA
    REAL( dp), intent(in)                            :: temp
    ! output
    REAL( dp), dimension(MMAX,IMAX),intent(out)      :: kelffect

    REAL( dp), DIMENSION(MMAX,0:(IMAX+1))  :: RP
    REAL( dp)                              :: exponente 

    INTEGER                                :: M,I

    ! for Ammonium nitrate
    REAL( dp), PARAMETER                   :: M_nh4   = 80.043 ! g/mol
    REAL( dp), PARAMETER                   :: st_nh4  = 0.1084 ! kg/s2
    REAL( dp), PARAMETER                   :: rho_nh4 = 1725   ! kg/m3 

     do M=1,MMAX
      do I=1,IMAX
        RP(M,I)=DPA(M,I)*0.5_dp  ! [m]
        exponente=2._dp*st_nh4*M_nh4*1.E-3_dp
        exponente=exponente/(R_gas*temp*rho_nh4*RP(M,I))
        kelffect(M,I)=exp(exponente)
      end do
     end do

  end subroutine kelvin_nit
  
  subroutine kelvin_org(DPA,temp,surfin,sforg,M_org1,DENOC,masscl,kelffect,IMAX)
    !----------------------------------------------------------------------
    !
    !****  Kelvin effect for organic vapours
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate Kelvin effect of organics (dimensionless)
    !      organic vapour represented by succinic acid
    !      12.11.2011: System definition: OV - H2SO4 - H2O
    !                  Activity coeff. f_i assumed to be unity
    !                  Raoult effect considered once OV is "dissolved"
    !                  Surface tension: sigma = X_i * f_i * sigma_i
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature              [K]
    !           DPAI     particle diameter            [m]
    !           masscl   mass conc. per bin           [ng/m^3]
    !           sforg    surface tension organic      [kg/s^2]
    !           M_org1   molecular weight organic     [g/mol]
    !           DENOC    density organic              [kg/m^3]
    ! 
    !        output:
    !           kelffect Kelvin effect of organic     [-]
    !
    !      method
    !      ------
    !      none
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
    use gde_constants,  only : M_H2O,pi,R_gas,M_H2SO4

    implicit none

    ! input
    INTEGER, intent(in)                              :: IMAX,surfin
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA
    REAL( dp), dimension(MMAX,IMAX,AMAX),intent(in)  :: masscl
    REAL( dp), intent(in)                            :: temp
    REAL( dp), intent(in)                            :: sforg,M_org1,DENOC
    ! output
    REAL( dp), dimension(MMAX,IMAX),intent(out)      :: kelffect

    REAL( dp), PARAMETER                   :: massmin=1.e-6
    REAL( dp), PARAMETER                   :: ST0=273.15_dp
    REAL( dp), PARAMETER                   :: dens_wat=1000.
    REAL( dp)                              :: gamma_os=1._dp


    REAL( dp), DIMENSION(MMAX,0:(IMAX+1))  :: RP
    REAL( dp)                              :: mos,msul,mwat
    REAL( dp)                              :: sigma_sul,sigma_org,sigma_wat
    REAL( dp)                              :: surft,rhot,rho_h2so4
    REAL( dp)                              :: x_os,x_wat,x_sul
    REAL( dp)                              :: exponente 


    INTEGER                                :: M,I
    
     ! calculate surface tension (kg/s^2)
     IF (surfin.EQ.1) THEN
       sigma_org=surf_succin(temp)
     ELSE
       sigma_org=sforg
     ENDIF
     !write(6,*) 'org',sigma_org
     sigma_wat=surf_h2o_std-(0.155_dp*(temp-ST0))
     sigma_wat=sigma_wat*1e-3_dp
    
     do M=1,MMAX
       do I=1,IMAX
         ! calculate mass (ng/m^3)
         mos =masscl(M,I,A_OR1)+masscl(M,I,A_OR2)+masscl(M,I,A_OR3)+  &
              masscl(M,I,A_OR4)+masscl(M,I,A_OR5)+masscl(M,I,A_OR6)+  &
              masscl(M,I,A_OR7)+masscl(M,I,A_OR8)+masscl(M,I,A_OR9)
         msul=masscl(M,I,A_SUL)
         mwat=masscl(M,I,A_WAT)
         ! calculate molar fractions (-) and
         ! surface tension of particle (kg/s^2)
         if (mos.gt.massmin) then
           x_os=(mos/M_org1) /((mos/M_org1) + (msul/M_H2SO4) &
             + (mwat/M_H2O))
           x_wat=(mwat/M_H2O)/((mos/M_org1) + (msul/M_H2SO4) &
              + (mwat/M_H2O))
           x_sul=1.-x_os-x_wat
!MSK 21.11.2020 minimum x_os = 0.1_dp
           x_os=max(x_os,0.1_dp)
           sigma_sul=surten(x_sul,temp)
           surft=x_wat*sigma_wat + x_os*sigma_org + &
                 x_sul*sigma_sul
           surft=min(surft,0.065_dp)
        ! calculate liquid density  (kg/m3)
           rho_h2so4=roolq(x_sul,temp)
           rhot=x_sul*rho_h2so4+x_wat*dens_wat+x_os*DENOC
           rhot=max(rhot,1200._dp)         
        else
           x_os=1.0_dp      !to neglect Raoult effect
           x_sul=0.0_dp
           x_wat=0.0_dp
           gamma_os=1.0_dp
           surft=0.070_dp   !Zhang and Wexler (2002)
           rhot=1000._dp    !water droplet
        endif
        RP(M,I)=DPA(M,I)*0.5_dp  ! [m]
        exponente=2._dp*surft*M_org1*1.E-3_dp
        exponente=exponente/(R_gas*temp*rhot*RP(M,I))
!MSK 21.11.2020 set x_os=1 for larger particles
        if (RP(M,I)>3.E-8) x_os=1.0_dp
        kelffect(M,I)=x_os*gamma_os*exp(exponente)
!MSK 21.11.2020 Ke not < 1.1 for NU and AI mode
        if (M<3) kelffect(M,I)=max(kelffect(M,I),1.10_dp)

	  ! write(6,*) 'ke',M,I,x_os,rhot,exp(exponente),kelffect(M,I)
      end do
     end do

  end subroutine kelvin_org
  
  subroutine kelvin_alkane(DPA,temp,carbonn,gammao,kelffect,IMAX)
    !----------------------------------------------------------------------
    !
    !****  Kelvin effect for n-alkanes (C16-C30)
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate Kelvin effect of alkanes (dimensionless)
    !      input:
    !      carbon number of n-alkane
    !      gamma is mole fraction of n-alkanes (0,...,1)
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature              [K]
    !           DPA      particle diameter            [m]
    !           gammao   mole fraction of alkane      [-]
    ! 
    !        output:
    !           kelffect Kelvin effect of alkanes     [-]
    !
    !      method
    !      ------
    !      see Zhang and Wexler et al.,
    !      Evolution of particle number distribution near roadways.
    !      Part II: the 'Road-toAmbient' process
    !      Atmos. Environ., 38,6655-6665, 2004
    !   
    !      surface tension of n-alkanes (@20degC):
    !      Jasper, J.J.; Kerr, E.R., Gregorich, F., The orthobaric surface 
    !      tensions and thermodynamic properties of the liquid surfaces of 
    !      the n-alkanes, C6 to C28, J. Am. Chem. Soc., 1953, 75, 5252-5254
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
    use gde_constants,  only : M_H2O,pi,R_gas

    implicit none

    ! input
    INTEGER, intent(in)                              :: IMAX, carbonn
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA
    REAL( dp), dimension(MMAX),intent(in)            :: gammao
    REAL( dp), intent(in)                            :: temp

    ! output
    REAL( dp), dimension(MMAX,IMAX),intent(out)      :: kelffect

    REAL( dp), DIMENSION(MMAX,0:(IMAX+1))  :: RP
    REAL( dp)                              :: exponente,exponente1
    REAL( dp)                              :: M_alk,surf_alk,Rv,dens_alkc
    !!!REAL( dp), PARAMETER                   :: surf_alk=0.030   ! C24H50
    !!!REAL( dp), PARAMETER                   :: dens_alk=1500.   ! kg/m^3
    ! C16H34 through C32H62 liquid densities in kg/m^3
    !   Data from CPC Handbook of Chemistry and Physics (2010)
    INTEGER, PARAMETER                     :: nalk = 15
    REAL( dp),dimension(nalk),parameter    :: dens_alk= (/ &
      770.1,   778.0,   776.8,   785.5,   788.6,   &
      791.9,   794.4,   778.5,   799.1,   801.2,   &
      778.3,   779.6,   806.7,   808.3,   809.7       /)

    INTEGER                                :: M,I,C
 
     C=carbonn-15 
     M_alk=14*carbonn+2                   ! [g/mol]
     surf_alk=2.1758+0.3762*carbonn
     surf_alk=(1/(surf_alk/M_alk))*1.e-3  ! [kg/s^2]
     dens_alkc=dens_alk(C)*1.e-6        ! [kg/cm^3]
     ! Rv=R_gas*1.E7 is 8.31E7 gcm^2 s^-2 mol^-1K^-1
     Rv=R_gas*1.e7
     
     do M=1,MMAX
      do I=1,IMAX
        RP(M,I)=DPA(M,I)*0.5*100  ! [cm]
        exponente1=2.*surf_alk*M_alk
	    ! reduce Ke ("heterogenous reactions")
        exponente1=exponente1*0.4
        exponente=exponente1/(Rv*temp*dens_alkc*RP(M,I))
        kelffect(M,I)=gammao(M)*exp(exponente)
        kelffect(M,I)=max(kelffect(M,I),1.0_dp)
	    ! write(6,*) M,I,kelffect(M,I)
      end do
     end do

  end subroutine kelvin_alkane


  subroutine sat_kelvin_amine(temp,DPAI,xaqAMIN,srkelvini)
    !----------------------------------------------------------------------
    !
    !****  Saturation ratio: Kelvin effect
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate Kelvin part of saturation ratio (dimensionless)
    !         for each size bin
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature              [K]
    !           xaqAMIN  molality of AMIN             [mol/kg]
    !           DPAI     particle diameter            [m]
    !        output:
    !           srkelvini saturation ratio            [-]
    !
    !      method
    !      ------
    !      Particle surface tension obtained from fit to Szyszkowski-Langmuir
    !      see Facchini et al., Nature 1999
    !      xaq must be MOLALITY in mol/kgC here!!!
    !      uses average particle surface tension and density (i.e. water)
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
    use gde_constants,    only : M_H2O

    implicit none

    ! input
    REAL( dp), intent(in)       :: temp
    REAL( dp), intent(in)       :: xaqAMIN      ! [mol/kg]
    REAL( dp), intent(in)       :: DPAI         ! [m]
    REAL( dp), intent(out)      :: srkelvini

    ! surface tension in [g s^-1] which is [dyn cm^-1]
    REAL( dp), parameter   :: rgassi = 8.31e+7_dp    ! gcm^2 s^-2mol^-1K^-1
    REAL( dp)              :: surfpa, surf_h2o,xaqca

    ! surface tension of water
    surf_h2o = surf_h2o_std - 0.155_dp*(temp - 273.15_dp)

    ! molality of amine in mol/kgC
    xaqca = xaqAMIN*fcAMIN

    ! mean surface tension of particle
    surfpa = surf_h2o - a_AMIN*temp*log(1.+b_AMIN*xaqca)

    ! kelvin part of saturation vapor ratio
    ! diameter [m] -> radius [cm]
    srkelvini = 1._dp+((2._dp*surfpa*(M_H2O))/(DPAI*0.5_dp*100._dp*rgassi*temp*pdens))
    if (srkelvini.LT.0.0_dp) then
        write(6,*) 'skelvin',srkelvini
        stop
    endif

  end subroutine sat_kelvin_amine


subroutine sat_raoult_amine(DPAI,numberi,caqAMIN,srraoulti)
    !----------------------------------------------------------------------
    !
    !****  Saturation ratio: Raoult effect
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate Raoult part of saturation ratio (dimensionless)
    !
    !      interface
    !      ---------
    !
    !        input:
    !           caqAMIN  aq. phase conc AMIN          [molec cm^-3]
    !           DPAI     particle diameter            [m]
    !           numberi  number conc                  [1/m^3]
    !
    !        output:
    !           srraoulti saturation ratio            [-]
    !
    !      method
    !      ------
    !      input is aqueous phase conc of amine in molec cm^-3 (air)
    !      volume of dissolved compound is neglected with respect to large
    !      water volume
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
    use gde_constants,    only : M_H2O,pi,rho_H2O

    implicit none

    ! input
    REAL( dp), intent(in)                :: caqAMIN
    REAL( dp), intent(in)                :: DPAI          ! [m]
    REAL( dp), intent(in)                :: numberi       ! [1/m3]
    REAL( dp), intent(out)               :: srraoulti
    REAL( dp)                            :: masspAMIN,terma,termb

    ! convert amine conc [molec cm^-3] to mass per particle [g particle^-1]
    ! Raoult part of saturation vapor ratio
    ! number 1/m3  diameter m -> radius m  caqAMIN molec/cm3 -> ?

      masspAMIN  = (caqAMIN*xmAMIN)/numberi
      terma      = 3._dp*M_H2O*0.001_dp*vhAMIN*masspAMIN
      termb      = 4._dp*pi*((DPAI*0.5_dp)**3)*rho_H2O*xmAMIN
      srraoulti  = 1._dp-(terma/termb)
      if (srraoulti.LT.0.0_dp) then
        write(6,*) 'sraoult',srraoulti,terma,termb,masspAMIN,DPAI
        stop
      endif

  end subroutine sat_raoult_amine


    !----------------------------------------------------------------------
    !----------------  NANO-KOEHLER MODEL         -------------------------
    !----------------------------------------------------------------------


 subroutine nano_koehler(DPA,temp,surfin,sforg,M_os,DENOC,masscl,se_wat,se_os,IMAX)
    !----------------------------------------------------------------------
    !
    !****  Nano-Koehler model
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate Nano-Koehler effect (dimensionless)
    !      input:
    !      -
    !      - 
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature              [K]
    !           DPA      particle diameter            [m]
    !           masscl   mass conc. per bin           [ng/m^3]
    !           sforg    surface tension organic      [kg/s^2]
    !           M_os     molecular weight organic     [g/mol]
    !           DENOC    density organic              [kg/m^3]
    ! 
    !        output:
    !           se_wat   equil. saturation water      [-]
    !           se_os    equil. saturation organic    [-]
    !
    !      method
    !      ------
    !      calculate Nano-Koehler effect in os-abs-wat mixture
    !
    !      reference
    !      ---------
    !      Kulmala, M., Kerminen, V.-M., Anttila, T., Laaksonen, A.,
    !      and O'Dowd, C.D. (2004).
    !      Organic aerosol formation via sulphate cluster activation
    !      J. Geophys. Res., 109, D04205, doi:10.1029/2003JD003961.
    !   
    !      Surface tension of aqueous solution of organic acids:
    !      Alvarez, E., G. Vazquez, M. Sanchez-Vilas, B. Sanjurjo, and 
    !      J. M. Navaza (1997).
    !      Surface tension of organic acids + water binary mixtures from
    !      20 degC to 50 degC.
    !      J. Chem. Eng. Data, 42, 957-960.
    !
    !      Surface tension of aqueous solution of ammonium sulphate:
    !      Korhonen, P., Laaksonen, A., Batris, E., and Viisanen, Y. (1998).
    !      Thermodynamics for highly concentrated water - ammonium
    !      sulfate solutions.
    !      J. Aerosol Sci., 29, Suppl 1, S379-S380. 
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    use gde_constants,  only : pi,R_gas,M_H2O,M_ABS

    implicit none

    ! input
    INTEGER, intent(in)                              :: IMAX,surfin
    REAL( dp), dimension(MMAX,IMAX,AMAX),intent(in)  :: masscl
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA
    REAL( dp), intent(in)                            :: temp,sforg,M_os,DENOC

    ! output
    REAL( dp), dimension(MMAX,IMAX),intent(out)      :: se_wat,se_os

    REAL( dp), DIMENSION(MMAX,0:(IMAX+1))  :: RP
    REAL( dp), PARAMETER                   :: M_wat=18.0154
    REAL( dp), PARAMETER                   :: dens_wat=1000.
    REAL( dp), PARAMETER                   :: dens_os=1570.
    REAL( dp), PARAMETER                   :: massmin=1.e-6
    REAL( dp)                              :: x_os,x_pb,x_wat,w_abs
    REAL( dp)                              :: sigma_org,surf_pb,surf2,surft
    REAL( dp)                              :: gamma_os,gamma_pb,gamma_wat
    REAL( dp)                              :: exponent_wat,exponent_os
    REAL( dp)                              :: exponente
    REAL( dp)                              :: ke_wat,ke_os
    REAL( dp)                              :: mos,mabs,mwat
    REAL( dp), PARAMETER                   :: A0=0.071912134 
    REAL( dp), PARAMETER                   :: A1=0.022387172
    REAL( dp), PARAMETER                   :: A2=-0.079960647
    REAL( dp), PARAMETER                   :: A3=0.698516114
    REAL( dp), PARAMETER                   :: A4=-2.361475486
    REAL( dp), PARAMETER                   :: A5=4.291669494
    REAL( dp), PARAMETER                   :: A6=-3.664358167
    REAL( dp), PARAMETER                   :: A7=1.143836868
! Nano-Koehler model parameter (Kulmala et al., 2004)
    REAL( dp), PARAMETER                   :: AA=-190.
    REAL( dp), PARAMETER                   :: BB=100.
    REAL( dp), PARAMETER                   :: CC=0.3
    REAL( dp), PARAMETER                   :: DD=0.8

    INTEGER                                :: M,I

! surface tension of organic vapour (os)
     IF (surfin.EQ.1) THEN
       sigma_org=surf_succin(temp)
     ELSE
       sigma_org=sforg
     ENDIF
! mole fraction os and pb in a os-abs-wat mixture
! mole fraction of water (wat) in a pb cluster
! weight fraction ammonium bisulfate (abs) in inorganic cluster (pb)
! mass concentrations (from input MASS(M,I)) in ng/m^3
! fractions are calculated per each size bin
! nano_koehler applies only to nucleation mode
!
! ATTENTION:
! MASS(NU,I,A_WAT) IS STILL ZERO
! => CHANGE IN GDEBOX
     do M=1,MMAX
      do I=1,IMAX
       mos =masscl(M,I,A_OR1)+masscl(M,I,A_OR2)
       mabs=masscl(M,I,A_SUL)+masscl(M,I,A_NH4)
       mwat=masscl(M,I,A_WAT)
       x_os=(mos/M_os)/max(((mos/M_os) + 2*(mabs/M_ABS) &
              + (mwat/M_H2O)),massmin)
       x_pb=1.-x_os 
       x_wat=(mwat/M_H2O)/max( ( (mwat/M_H2O)  &
             + 2*(mabs/M_ABS) ),massmin)
       !lower limit: x_wat=0.25       
       x_wat=max(x_wat,0.25_dp)          
       w_abs=(mabs)/max((mabs + mwat),massmin)
       !upper limit: w_abs=0.9
       w_abs=min(w_abs,0.9_dp)      
!! surface tension aqueous ammonium bisulfate solution (=inorganic cluster pb)
       surf_pb=A0+A1*w_abs+A2*w_abs**2+A3*w_abs**3+A4*w_abs**4+ &
            A5*w_abs**5+A6*w_abs**6+A7*w_abs**7    ! [kg/s^2]
!! surface tension of pseudobinary solution (pb+os)
!! Error in equation (4) in Kulmala et al. (2004). 
!!   correct is: sigma_pb-(sigma_pb-sigma_os)*sur2*xos
       surf2=1.+( (CC*x_pb)/(1.-(DD*x_pb)) )
       surft=surf_pb-((surf2*x_os)*(surf_pb-sigma_org))
!! activity coefficients
       gamma_os=exp(  (AA*BB)/( R_gas*temp*(1.+(BB*x_os/x_pb))**2 )  )
       gamma_pb=exp(  AA/( R_gas*temp*(1.+(x_pb/(BB*x_os)))**2 )  )    
       gamma_wat=gamma_pb*(1.-0.305*w_abs-0.294*w_abs**2-0.443   &
                *w_abs**3)/x_wat
       !write(6,*) M,I,surft,gamma_pb,gamma_wat,w_abs
       RP(M,I)=DPA(M,I)*0.5      ! [m]
       exponent_wat=(2.*surft*M_wat*1.E-3)/(R_gas*temp*dens_wat*RP(M,I))
       exponent_os =(2.*surft*M_os *1.E-3)/(R_gas*temp*dens_os *RP(M,I))
       ke_wat=exp(exponent_wat)
       ke_os=exp(exponent_os)
       ! equilibrium saturation ratio 
       se_wat(M,I)=x_wat*gamma_wat*ke_wat
       se_os(M,I)=x_os*gamma_os*ke_os
       !write(6,*) 'Nano',M,I,se_wat(M,I),se_os(M,I)
      end do
     end do
     ! overwrite se_os for AI:CS
     do M=2,MMAX
      do I=1,IMAX
        RP(M,I)=DPA(M,I)*0.5    ! [m]
        exponente=2.*sigma_org*M_os*1.E-3
        exponente=exponente/(R_gas*temp*DENOC*RP(M,I))
        se_os(M,I)=exp(exponente)
      end do
     end do
    
  end subroutine nano_koehler



  end module gde_transfer
