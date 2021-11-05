! <gde_plume.f90 - A component of the Multicomponent
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
module gde_plume

  use messy_mecca_kpp_Parameters, only : dp, NSPEC

! gas-phase species
  use messy_mecca_kpp_parameters, only : ind_O3,ind_NO,ind_NO2,ind_SO2
  use messy_mecca_kpp_parameters, only : ind_H2SO4,ind_SO3
  use messy_mecca_kpp_parameters, only : ind_N2O5,ind_HNO3,ind_NH3
! condensable organics
  use messy_mecca_kpp_parameters, only : ind_BSOV,ind_BLOV,ind_BELV
  use messy_mecca_kpp_parameters, only : ind_ASOV,ind_ALOV,ind_AELV
  use messy_mecca_kpp_parameters, only : ind_PIOV,ind_PSOV,ind_PELV

  use gde_sensitiv,  only       : IDIL, ICHAM, IDEB

  use gde_constants, only       : M_H2SO4,pi

  use gde_init_aero, only       : BGNO,BGNO2,BGSO2,BGO3
  use gde_init_aero, only       : BGNH3,BGSULF,BGPIOV,BGPSOV

  use gde_toolbox,   only       : molecdiff
  use gde_toolbox,   only       : acidps


implicit none

   private

   public :: plumedisp
   public :: readdispers
   public :: initplume
   public :: plumearea

! deposition
   public :: ADEP,BDEP,ZCAP,dcol,Fplus
   public :: ustar,znot
! dispersion
   public :: hmix_st,dst_st,hsta_st,ta_st
   public :: DR_fin,T_fin,tau_c,tau_d,BGH2O
   public :: dil2_b,dil2_c,dil2_d,dil2_e,dil2_f
   public :: u0,sigw,tend1,tbeg2,tend2
! clouds
   public :: vupdra,sst,sal

   real( dp),save     :: ADEP,BDEP,ZCAP,DCOL,Fplus
   real( dp),save     :: ustar,znot
   real( dp),save     :: dil2_b,dil2_c,dil2_d,dil2_e,dil2_f
   real( dp),save     :: hmix_st,dst_st,hsta_st,ta_st
   real( dp),save     :: DR_fin,T_fin,tau_c,tau_d,BGH2O
   real( dp),save     :: u0,sigw,tend1,tbeg2,tend2
   real( dp),save     :: vupdra,sst,sal

contains

  subroutine plumedisp(DT,tair,told,zplum_old,wplum_old,u10,dilut_time,  &
                       dilstore,dila,dilcoef,zplum,wplum,tnew,dilrate,   &
                       emisratp,   cgas  )
    !----------------------------------------------------------------------
    !     
    !  Control of the plume dispersion module
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      call of subroutines of the plume module
    !      calculates plume height and T at current timestep
    !
    !      interface
    !      ---------
    !
    !        input:
    !          DT               model time step                 [s]
    !          firstloop        flag for initialisation
    !          tair             ambient air temperature         [K]
    !          told             plume temperature, old timestep [K]
    !          zplum_old        plume height, old timestep      [m]
    !          wplum_old        plume width, old timestep       [m]
    !          u10              wind speed                      [s]
    !          dilut_time       time passed in plume            [s]
    !          dilstore         time step of 0.5 s              [s]
    !          dila             dispersion parameter A          [-]
    !          dilcoef          dilution coefficient B          [-]
    !  
    !        output:
    !          zplum            plume height, current time      [m]
    !          wplum            plume width,  current time      [m]
    !          tnew             plume temperature, current time [K]
    !          dilrate          dilution rate of concentration  [1/s]
    !          emisratp         emission rate source1/source2   [-]
    !
    !        in/out:
    !           cgas            concentration of gases          [molec./cm3]
    !
    !      method
    !      ------
    !      IDIL = 1
    !           plume dispersion type 1
    !           Pohjola et al. Atm. Environ., 37, 339-351, 2003
    !           The dilution of the traffic-generated aerosol
    !           by background air is approximated by fitting a 
    !           power-law function y = a*x**(-b) 
    !           where x = u*t with u the wind speed and t 
    !           the time.
    !           a is dila and b is dilcoef in input file ingeod.dat
    !           dila     = 86.49     ! dispersion parameter A   [-]
    !           dilcoef  = 0.92332   ! dilution coefficient     [-]
    !
    !      IDIL = 2
    !           plume dispersion type 2
    !           dil2_a,dil2_b,dil2_c for concentration
    !           dil2_d,dile_e,dil2_f for temperature and zmbl
    !           Jari Walden, 2013
    !           dil2_b   = 0.202     ! Jari parameter b   [s]
    !           dil2_c   = -1.85     ! Jari parameter c   [-]
    !           dil2_d   = 9.2       ! Jari parameter d   [-]
    !           dil2_e   = 0.25      ! Jari parameter e   [s]
    !           dil2_f   = -1.95     ! Jari parameter f   [-]
    !
    !      IDIL = 3
    !           plume dispersion type 3
    !           Ronkko et al. 2013 diesel exhaust after-treatment
    !           DR_fin   = 12.0;     ! [-]
    !           T_fin    = 300.0;    ! [K]
    !           tau_c    = 0.03      ! [s]
    !           tau_d    = 0.12      ! [s]
    !           BGH2O    = 1.E16     ! [molecules cm^-3]
    !
    !      IDIL = 4
    !           plume dispersion type 4
    !           Urban Case, two stage dilution
    !           Stage 1 dilution rate according to
    !           Vignati et al., Sci. Total Environ., 235, 1999
    !           u0       = 0.23      ! [m/s]
    !           sigw     = 0.29      ! [m/s]
    !           tend1    = 7.0       ! [m]  width line 1
    !           tbeg2    = 13.0      ! [m]  distance to line 2
    !           tend2    = 22.5      ! [m]  width street
    !
    !      IDIL = 5
    !           plume dispersion type 5
    !           Gaussian plume with semi-elliptic cross section
    !           Jana Moldanova, MOCCA model, 2021
    !
    !      IDIL = 6
    !           plume dispersion type 6
    !           Modified Gaussian plume, Konopka 1995 formulation
    !           Jana Moldanova, MOCCA model, 2021
    !
    !      IDIL = 7
    !           plume dispersion type 7
    !           Ship plume in marine BL, open sea
    !           Chosson et al. (2008)
    !           Jana Moldanova, MOCCA model, 2021
    !
    !
    !      reference
    !      ---------
    !      Chosson, F., Paoli, R., and B. Cuenot,
    !        Ship plume dispersion rates in convective boundary layers
    !        for chemistry models, Atmos. Chem. Phys. 8, 4841-4853
    !        www.atmos-chem-phys.net/8/4841/2008/, 2008.
    !
    !      Konopka, P., Analytical Gaussian Solutions for Anisotropic
    !        Diffusion in a Linear Shear Flow, 
    !        J. Non-Equilib. Thermondyn., 20, 78-91, 1995.
    !
    !      Pohjola, M.A., Pirjola, L., Karppinen, A., Harkonen, J.,
    !        Korhonen, H., Hussein, T., Ketzel, M., and J. Kukkonen,
    !        Evaluation and modelling of the size fractionated aerosol
    !        particle number concentration measurements nearby a
    !        major road in Helsinki - Part I: Modelling results
    !        within the LIPIKA project, Atmos. Chem. Phys. 7, 4065-4080,
    !        www.atmos-chem-phys.net/7/4065/2007/, 2007.
    !
    !      Vignati, E., Berkowicz, R., Palmgren, F., Lyck, E., and
    !        Hummelshoj, P.,
    !        Transformation of size distributions of emitted particles 
    !        in streets.
    !        The Science of The Total Environment 235, 37-49, 1999.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

  implicit none

! input
     real( dp), intent(in)                      :: DT           ![s]
     real( dp), intent(in)                      :: tair         ![K]
     real( dp), intent(in)                      :: told         ![K]
     real( dp), intent(in)                      :: zplum_old    ![m]
     real( dp), intent(in)                      :: wplum_old    ![m]
     real( dp), intent(in)                      :: u10          ![m/s]
     real( dp), intent(in)                      :: dilut_time   ![s]
     real( dp), intent(in)                      :: dilstore     ![s]
     real( dp), intent(in)                      :: dila
     real( dp), intent(in)                      :: dilcoef


     real( dp),dimension(NSPEC), intent(in out) :: cgas         ![molec/cm3]

! output
     real( dp), intent(out)                     :: zplum        ![m]
     real( dp), intent(out)                     :: wplum        ![m]
     real( dp), intent(out)                     :: tnew         ![K]
     real( dp), intent(out)                     :: dilrate      ![1/s]
     real( dp), intent(out)                     :: emisratp     ![-]

! local
     real( dp)                                  :: wline1       ![m]
     real( dp)                                  :: wline2       ![m]
     real( dp)                                  :: ch2so4       ![molec/cm3]
     real( dp), parameter                       :: ratio_tc=0.744
!----------------------------------------------------------------

 ! Re-calculate mixing height of traffic air parcel
        call plumeheight(zplum_old,wplum_old,dilut_time,dila,dilcoef, &
                         u10,zplum,wplum)

 ! Re-calculate plume temperature (tnew)
        call plumetemp(DT,tair,dilut_time,told,zplum,zplum_old,tnew )

 ! Re-calculate dilution rate for concentration
        call plumedilr(dilut_time,dilstore,dilcoef,u10,DT,wplum_old,  &
                       wplum,dilrate)

 ! Wall-loss in diesel exhaust chamber
        if (IDIL==3) then
          if (dilut_time.ge.tau_d) then
          !write(6,*) 'before wall loss', c(ind_h2so4)
              ch2so4=cgas(ind_H2SO4)
              if (ICHAM .eq. 2)  call wall_loss_h2so4(DT,tnew,ch2so4)
              cgas(ind_H2SO4)=ch2so4
          !write(6,*) 'after wall loss', c(ind_h2so4)
          endif
        endif

 ! Debug: log the plume parameters --> into plume output !!!
        if (IDEB == 1) then
           write(12,fmt='(a,f8.2,f8.2,f8.2,f8.2)') 'z_PL Temp t_dil r_dil', &
                    zplum,tnew,dilut_time,dilrate
        endif


 ! Emission ratio source1/source2
 ! no emission accounted between tend1 and tbeg2 or after tend2
        emisratp = 1.0
        wline1   = tend1*u10
        wline2   = (tend2 - tbeg2)*u10
        if (IDIL==4) then
           if (dilut_time.gt.tend2*u10) then
              emisratp=0.0
           elseif ( (dilut_time.gt.tend1*u10).and. &
                    (dilut_time.le.tbeg2*u10) ) then
              emisratp=0.0
           elseif ( (dilut_time.gt.tbeg2*u10).and. &
                    (dilut_time.le.tend2*u10) ) then
           ! line source 2
              emisratp=(wline1/wline2)*ratio_tc
           else
              emisratp=1.0
           endif
        endif
                          
                          
 ! Temperature correction of gases, according to the Ideal Gas Law
 ! applied to the gases in the plume
        cgas(ind_NO)    = cgas(ind_NO)    * (told/tnew)
        cgas(ind_NO2)   = cgas(ind_NO2)   * (told/tnew)
        cgas(ind_SO2)   = cgas(ind_SO2)   * (told/tnew)
        cgas(ind_H2SO4) = cgas(ind_H2SO4) * (told/tnew)
        cgas(ind_SO3)   = cgas(ind_SO3)   * (told/tnew)
        cgas(ind_O3)    = cgas(ind_O3)    * (told/tnew)
        cgas(ind_N2O5)  = cgas(ind_N2O5)  * (told/tnew)
        cgas(ind_NH3)   = cgas(ind_NH3)   * (told/tnew)
        cgas(ind_HNO3)  = cgas(ind_HNO3)  * (told/tnew)

 ! applied to organic vapors
        cgas(ind_BSOV)  = cgas(ind_BSOV)  * (told/tnew)
        cgas(ind_BLOV)  = cgas(ind_BLOV)  * (told/tnew)
        cgas(ind_BELV)  = cgas(ind_BELV)  * (told/tnew)
        cgas(ind_ASOV)  = cgas(ind_ASOV)  * (told/tnew)
        cgas(ind_ALOV)  = cgas(ind_ALOV)  * (told/tnew)
        cgas(ind_AELV)  = cgas(ind_AELV)  * (told/tnew)
        cgas(ind_PIOV)  = cgas(ind_PIOV)  * (told/tnew)
        cgas(ind_PSOV)  = cgas(ind_PSOV)  * (told/tnew)
        cgas(ind_PELV)  = cgas(ind_PELV)  * (told/tnew)


 ! New gas-phase concentrations after dilution

 ! dilution with zero background concentration
        cgas(ind_BSOV)  = cgas(ind_BSOV) - ( dilrate * cgas(ind_BSOV)*DT)
        cgas(ind_BLOV)  = cgas(ind_BLOV) - ( dilrate * cgas(ind_BLOV)*DT)
        cgas(ind_BELV)  = cgas(ind_BELV) - ( dilrate * cgas(ind_BELV)*DT)
        cgas(ind_ASOV)  = cgas(ind_ASOV) - ( dilrate * cgas(ind_ASOV)*DT)
        cgas(ind_ALOV)  = cgas(ind_ALOV) - ( dilrate * cgas(ind_ALOV)*DT)
        cgas(ind_AELV)  = cgas(ind_AELV) - ( dilrate * cgas(ind_AELV)*DT)
        cgas(ind_PELV)  = cgas(ind_PELV) - ( dilrate * cgas(ind_PELV)*DT)

        cgas(ind_SO3)   = cgas(ind_SO3)  - ( dilrate * cgas(ind_SO3)*DT)

 ! dilution of NO, NO2, SO2, NH3, H2SO4, PIOV, PSOV
 ! with backround concentration from inbgair.dat
 ! (note: exhaust concentration higher than background)
        cgas(ind_NO)    = cgas(ind_NO)   - ( dilrate * ((cgas(ind_NO)-BGNO)    *DT) )
        cgas(ind_NO2)   = cgas(ind_NO2)  - ( dilrate * ((cgas(ind_NO2)-BGNO2)  *DT) )
        cgas(ind_SO2)   = cgas(ind_SO2)  - ( dilrate * ((cgas(ind_SO2)-BGSO2)  *DT) )
        cgas(ind_PIOV)  = cgas(ind_PIOV) - ( dilrate * ((cgas(ind_PIOV)-BGPIOV)*DT) )
        cgas(ind_PSOV)  = cgas(ind_PSOV) - ( dilrate * ((cgas(ind_PSOV)-BGPSOV)*DT) )
        cgas(ind_NH3)   = cgas(ind_NH3)  - ( dilrate * ((cgas(ind_NH3)-BGNH3)*DT)   )
        cgas(ind_H2SO4) = cgas(ind_H2SO4)- ( dilrate * ((cgas(ind_H2SO4)-BGSULF)*DT))

 ! entrainment of O3
        cgas(ind_O3)    = cgas(ind_O3)   + ( dilrate * (max((BGO3-cgas(ind_O3)),0.0_dp)*DT) )

 ! approximated dilution of H2SO4
 !       cgas(ind_H2SO4) = cgas(ind_H2SO4)- ( dilrate * 0.85*((cgas(ind_H2SO4)-1.0e05)*DT) )


  end subroutine plumedisp


!------------------------------------------------------------------

  subroutine readdispers()
    !----------------------------------------------------------------------
    !     
    !   Read input for plume dispersion and deposition
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Read input for plume dispersion and deposition
    !
    !      interface
    !      ---------
    !      dispers.dat
    !
    !      method
    !      ------
    !      read ascii file with space or tab separated entries
    !
    !      reference
    !      ---------
    !      none
    !
    !      modifications     real( dp), intent(out)                     :: cco3
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

  implicit none

! local
     integer               :: stat6


! read dispers.dat

       open(28,file='dispers.dat',status='old', iostat=stat6)
! open error handling    
       if (stat6.ne.0) then
          write(6,*) 'File dispers.dat cannot be opened !'
          stop
       end if 
       read(28,*) hmix_st,dst_st,hsta_st,ta_st
       read(28,*) dil2_b,dil2_c,dil2_d,dil2_e,dil2_f
       read(28,*) DR_fin,T_fin,tau_c,tau_d,BGH2O
       read(28,*) u0,sigw,tend1,tbeg2,tend2
       read(28,*) ustar,znot,ADEP,BDEP
       read(28,*) ZCAP,dcol,Fplus
       read(28,*) vupdra,sst,sal

       if (vupdra.gt.0.5_dp) then
         write(6,*) 'STOP: V_updraft must be < 0.5 m/s in dispers.dat'
         stop
       endif

       close(28)


  end subroutine readdispers


  subroutine plumeheight(zplum_old,wplum_old,dilut_time,dila,dilcoef,u10,&
                         zplum,wplum)
    !----------------------------------------------------------------------
    !     
    !  Calculate the plume height H_pl and plume width W_pl
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Calculate height of the expanding air parcel
    !      Dispersion parameters from HIWAY-2, US EPA for
    !      different stability conditions
    !      a' = dila is from dispers.dat
    !      b' = 0.91 for stable
    !      For circular plume cross section (type 1,2,4),
    !      the width of the plume is the same as its height.
    !      For dispersion type 3, plume width and height are 
    !      constant.
    !      For semi-elliptic plume cross section (type 5-7)
    !      the width is calculated here and the plume height
    !      is approximated with the formulation of 
    !      von Glasow et al. (2003) 
    !
    !      interface
    !      ---------
    !
    !        input:
    !          zplum_old        plume height, old timestep      [m]
    !          wplum_old        plume width, old timestep       [m]    
    !          dilut_time       time passed in plume            [s]
    !          dila             dispersion parameter A          [-]
    !          dilcoef          dilution coefficient B          [-]
    !          u10              wind speed                      [s]
    !
    !        output:
    !          zplum            plume height, new timestep      [m]
    !          wplum            plume width, new timestep       [m]
    !
    !      method
    !      ------
    !      NILU TR 12/2003, p.29: road-to-ambient type 1,2,4
    !      plume area constant for dispersion type 3
    !      MOCCA model and von Glasow et al. (2003) type 5-7
    !
    !      reference
    !      ---------
    !      Karl, M., Kukkonen, J., Keuken, M.P., Lutzenkirchen, S.,
    !        Pirjola, L., and T. Hussein, Modeling and measurements
    !        on the neighborhood scale in Rotterdam, Oslo and Helsinki,
    !        Atmos. Chem. Phys., 16, 4817-4835, 
    !        www.atmos-chem-phys.net/16/4817/2016/, 2016.
    !      von Glasow, R., Lawrence, M.G., Sander, R., and 
    !        P.J. Crutzen, Modeling the chemical effects of ship
    !        exhaust in the cloud-free marine boundary layer,
    !        Atmos. Chem. Phys., 3, 233-250,
    !        www.atmos-chem-phys.org/acp/3/233/, 2003.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

  implicit none

! input
     real( dp), intent(in)                      :: dilut_time    ![s]
     real( dp), intent(in)                      :: dila
     real( dp), intent(in)                      :: dilcoef
     real( dp), intent(in)                      :: u10           ![m/s]
     real( dp), intent(in)                      :: zplum_old     ![m]
     real( dp), intent(in)                      :: wplum_old     ![m]

! output
     real( dp), intent(out)                     :: zplum         ![m]
     real( dp), intent(out)                     :: wplum         ![m]

! local
     real( dp), parameter                       :: alpha   = 0.1 ![-]
     real( dp), parameter                       :: hkerb   = 3.5 ![m]
     real( dp), parameter                       :: APARpld = 0.2
     real( dp), parameter                       :: APARvgl = 0.75
     real( dp), parameter                       :: BPARpld = 0.6
     real( dp), parameter                       :: TILL    = 3.0 ![m]
     real( dp), parameter                       :: t0      = 1.0 ![s]
     real( dp), parameter                       :: plwith0 = 6.3 ![m]
     real( dp), parameter                       :: plhight0= 3.5 ![m] 
     real( dp)                                  :: s0,st
     real( dp)                                  :: dilt2

!----------------------------------------------------------------


        SELECT CASE (IDIL)
            CASE (1)

                 ! z_pl = sqrt( z_pl(0)**2 + (a'*(1e-3*u10*t)**b')**2 )
                 ! a' = dila
                 ! b' = dilcoef
                 !
                 zplum = sqrt((hmix_st)**2 +                           &
                        (dila*(dilut_time*u10*1.e-3)**dilcoef)**2)
                 wplum = zplum

            CASE (2)

                 ! z_pl = sqrt( z_pl(0)**2 + (a'*(1e-3*u10*t)**b')**2 )
                 ! a' = dila
                 ! b' = dilcoef
                 !
                 zplum = sqrt((hmix_st)**2 +                           &
                        (dila*(dilut_time*u10*1.e-3)**dilcoef)**2)
                 wplum = zplum

            CASE (3)

                 ! do not change plume height
                 !
                 zplum = zplum_old
                 wplum = wplum_old

            CASE (4)

                 if (dilut_time.le.tend1*u10) then

                 ! Stage 1:
                 ! S0 = pi*hmix_st**2
                 ! S(t) = (sqrt(S0) + t*sigw)**2 - (t*alpha*u0)**2
                 ! z_pl = sqrt( S(t)/pi )
                 !
                   dilt2 = 0.0
                   s0 = pi*hmix_st**2
                   st = ( sqrt(s0) + dilut_time*sigw )**2  -          &
                       ( dilut_time*alpha*u0 )**2
                   zplum = sqrt(st/pi)

                 else if ((dilut_time.gt.tend1*u10).and.   &
                          (dilut_time.le.tbeg2*u10)) then  
                 ! no change over "tram tracks"
                   dilt2 = 0.0
                   zplum = zplum_old

                 else if ((dilut_time.gt.tbeg2*u10).and.    &
                          (dilut_time.le.tend2*u10)) then

                   dilt2 = dilut_time - (tbeg2-tend1)*u10
                   s0 = pi*hmix_st**2
                   st = ( sqrt(s0) + dilt2*sigw )**2  -               &
                       ( dilt2*alpha*u0 )**2
                   zplum = sqrt(st/pi)

                 else if (dilut_time.gt.tend2*u10) then

                 ! Stage 2:
                 ! start with hmix_st --> should be zplum at end of stage 1
                 ! z_pl = sqrt( z_pl(0)**2 + (a'*(1e-3*u10*t)**b')**2 )
                 ! a' = dila
                 ! b' = dilcoef = 0.91 (stable)
                 ! here dilution time starts at 0 again
                 !
                   dilt2 = dilut_time - tend2*u10
                   zplum = sqrt((hkerb)**2 +                         &
                          ( dila*(dilt2*u10*1.e-3)**(0.91) )**2)

                 endif
                 ! plume width same as plume height
                 wplum = zplum

            CASE (5)

                 ! Gaussian dispersion, type 5
                 wplum = 1000.*APARpld*                              &
                         EXP( BPARpld*LOG( (dilut_time*u10+TILL)     &
                         /1000.) )
                 zplum = plhight0 * ( (dilut_time+t0)/t0 )**BPARpld

            CASE (6)

                 ! Modified Gaussian / Konopka, type 6
                 wplum = 1000.*APARpld*                              &
                         EXP( BPARpld*LOG( (dilut_time*u10+TILL)     &
                         /1000.) )
                 zplum = plhight0 * ( (dilut_time+t0)/t0 )**BPARpld

            CASE (7)

                 ! Chosson, open sea MBL, type 7
                 ! original von Glasow et al. (2003)
                 wplum = plwith0  * ( (dilut_time+t0)/t0 )**APARvgl
                 zplum = plhight0 * ( (dilut_time+t0)/t0 )**BPARpld

            CASE DEFAULT
 
              write(6,*) 'no valid plume dispersion option'
               STOP

         END SELECT

      !   write(6,*) 'type ',IDIL,' pl height ',zplum,' pl width ',wplum  

  end subroutine plumeheight


  subroutine plumedilr(dilut_time,dilstore,dilcoef,u10,DT,wplum_old,&
                       wplum,dilrate)
    !----------------------------------------------------------------------
    !     
    !  Calculate the dilution rate for concentration
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Calculate the dilution rate for concentration
    !      for the dilution of the plume with background air
    !
    !      interface
    !      ---------
    !
    !        input:
    !          dilut_time       time passed in plume            [s]
    !          dilstore         time step of 0.5 s              [s]
    !          dilcoef          dilution coefficient B          [-]
    !          u10              wind speed                      [s]
    !          DT               model time step                 [s]
    !          wplum_old        plume width, old timestep       [m]
    !          wplum            plume width, new timestep       [m]
    !
    !        output:
    !          dilrate          plume dilution rate             [1/s]
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

  implicit none

! input
     real( dp), intent(in)                      :: dilut_time    ![s]
     real( dp), intent(in)                      :: dilstore      ![s]
     real( dp), intent(in)                      :: dilcoef       ![-]
     real( dp), intent(in)                      :: u10           ![m/s]
     real( dp), intent(in)                      :: DT            ![s]
     real( dp), intent(in)                      :: wplum_old     ![m]
     real( dp), intent(in)                      :: wplum         ![m]

! output
     real( dp), intent(out)                     :: dilrate       ![1/s]

! local
     real( dp), parameter                       :: alpha   = 0.1 ![-]
     ! Konopka
     real( dp), parameter                       :: kdiffh  = 25. ![m^2/s]
     real( dp), parameter                       :: kdiffv  = 10. ![m^2/s]
     real( dp), parameter                       :: plwith0 = 6.3 ![m]
     real( dp), parameter                       :: plhight0= 3.5 ![m] 
     ! PSD2DIM relation between st. dev. and plum with/height
     real( dp), parameter                       :: PSD2DIM = 2.2
     ! Chosson
     real( dp), parameter                       :: BPARcho = 1.659
     real( dp), parameter                       :: CPARcho = 1.133
     
     real( dp)                                  :: dilut2_time
     real( dp)                                  :: ddrdt
     real( dp)                                  :: DR
     real( dp)                                  :: s0
     real( dp)                                  :: dilt2
     real( dp)                                  :: qplume
     real( dp)                                  :: plstdevh0
     real( dp)                                  :: plstdevv0
     real( dp)                                  :: dwidthdt

!----------------------------------------------------------------


        SELECT CASE (IDIL)
            CASE (1)
                 ! Plume dispersion type 1
                 !
                 dilrate = dilcoef/dilut_time

            CASE (2)
                 ! Plume dispersion type 2
                 !
                 if ((dil2_c > 0.).or.(dil2_c < -2.)) then
                    write(6,*) 'Type 2: dil2_c has to be between 0 and -2'
                    stop
                 endif
                 dilut2_time = max(dilut_time+dil2_b,1.0_dp)
                 dilrate     = min( ((-dil2_c)/dilut2_time),1.90_dp )     

            CASE (3)
                 ! Plume dispersion type 3
                 ! zero dilution after tau_d is reached
                 !
                 if ( tau_d .le. 0.0 ) then
                    write(6,*) 'Type 3: tau_d has to be greater than 0'
                    stop
                 endif
                 dilrate = log(DR_fin)/tau_d
                 if (dilut_time.ge.tau_d)  dilrate = 0.0

            CASE (4)
                 ! Plume dispersion type 4
                 !
                 if (dilut_time.le.tend1*u10) then

                 ! Stage 1: 
                 !   dDR/dt = ( -2a**2*u0**2*t +
                 !            2sigw*(sqrt(s0)+sigw*t) ) /s0
                 !   DR= 1 + dDR/dt *dilut_time
                 !   dilrate = dDR/dt / DR**2
                 !
                   dilt2   = 0.0
                   s0      = pi*hmix_st**2
                   ddrdt   = ( -2.0*alpha**2 *u0**2 *dilstore         + &
                             2.0*sigw* (sqrt(s0)+sigw*dilstore) ) / s0
                   DR      = 1.0 + ddrdt*dilstore
                   dilrate = ddrdt / DR**2

                 else if ((dilut_time.gt.tend1*u10).and.   &
                          (dilut_time.le.tbeg2*u10)) then

                 ! no change over "tram tracks"
                   dilt2 = tend1
                   s0    = pi*hmix_st**2
                   ddrdt = ( -2.0*alpha**2 *u0**2 *dilt2                + &
                             2.0*sigw* (sqrt(s0)+sigw*dilt2) ) / s0
                   DR    = 1.0 + ddrdt*dilt2
                   dilrate = ddrdt / DR**2

                 else if ((dilut_time.gt.tbeg2*u10).and.    &
                          (dilut_time.le.tend2*u10)) then

                   dilt2 = dilstore - (tbeg2-tend1)*u10
                   s0    = pi*hmix_st**2
                   ddrdt = ( -2.0*alpha**2 *u0**2 *dilt2                + &
                             2.0*sigw* (sqrt(s0)+sigw*dilt2) ) / s0
                   DR    = 1.0 + ddrdt*dilt2
                   dilrate = ddrdt / DR**2

                 else if (dilut_time.gt.tend2*u10) then

                   ddrdt   = 1.0
                   dilt2   = dilstore + 1.0 - tend2*u10 -0.3
                 ! Stage 2:
                 ! Plume dispersion type 1
                 !
                   dilrate = dilcoef/dilt2

                 endif

              !   print *,'dilrt',dilut_time,dilstore,dilt2,dilrate,ddrdt

            CASE (5)
                 ! Plume dispersion type 5
                 ! Gaussian dispersion
                 ! plume dispersion lambda(t)=1/A*dA/dt
                 !
                 dwidthdt = (wplum - wplum_old) / DT
                 dilrate  = dwidthdt/wplum

            CASE (6)
                 ! Plume dispersion type 6
                 ! Modified Gaussian / Konopka
                 ! plstdevh0: plume st.dev. horizontal
                 ! plstdevv0: plume st.dev. vertical
                 !
                 plstdevh0 = plwith0 /PSD2DIM
                 plstdevv0 = plhight0/PSD2DIM
                 qplume    = 4.*kdiffh*kdiffv*dilut_time +                &
                            2.*kdiffh*plstdevv0**2.      +                &
                            2.*kdiffv*plstdevh0**2
                 dilrate   = qplume / ( dilut_time*qplume +               &   
                            plstdevh0**2. * plstdevv0**2. )

            CASE (7)
                 ! Plume dispersion type 7
                 ! Chosson et al. (2008), open sea MBL
                 ! for buoyancy flux of 250 m^4/s^3
                 ! dilution rate in 1/min converted to 1/s
                 dilrate  = BPARcho/ (dilut_time/60.)**CPARcho
                 dilrate  = dilrate / 60. 

            CASE DEFAULT
 
              write(6,*) 'no valid plume dispersion option'
               STOP

         END SELECT

     !    write(6,*) 'type ',IDIL,' diltime ',dilut_time,' dilrate ',dilrate  


  end subroutine plumedilr


  subroutine plumetemp(DT,TAIR,dilut_time,t_old,zplum,zplum_old,t_new)
    !----------------------------------------------------------------------
    !     
    !  Calculate the in-plume temperature
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Calculate new temperature in plume for the current time step
    !      homogeneous temperature in the plume cross-section
    !      TAIR is the surrounding ambient air temperature
    !
    !      interface
    !      ---------
    !
    !        input:
    !          DT               model time step                 [s]
    !          TAIR             ambient air temperature         [K]
    !          dilut_time       time passed in plume            [s]
    !          t_old            plume temperature, old timestep [K]
    !          zplum            plume height, new timestep      [m]
    !          zplum_old        plume height, old timestep      [m]
    !
    !        output:
    !          t_new            plume temperature, new timestep [K]
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

  implicit none

! input
     real( dp), intent(in)                      :: DT          ![s]
     real( dp), intent(in)                      :: TAIR        ![K]
     real( dp), intent(in)                      :: dilut_time  ![s]
     real( dp), intent(in)                      :: t_old       ![K]
     real( dp), intent(in)                      :: zplum       ![m]
     real( dp), intent(in)                      :: zplum_old   ![m]

! output
     real( dp), intent(out)                     :: t_new       ![K]

! local
     real(dp)                                   :: DRPL
     real(dp)                                   :: dilut2_time

!----------------------------------------------------------------


        SELECT CASE (IDIL)
            CASE (1)
                 ! T_new = (T_old-T_air)/DR + T_air
                 !
                 DRPL=zplum/zplum_old
                 t_new=((t_old-TAIR)/DRPL)+TAIR
                 t_new=MAX(t_new,TAIR)

            CASE (2)
                 ! T_new = dil2_d*(dt + dil2_e)**dil2_f + T_air
                 !
                 dilut2_time=max(dilut_time+dil2_e,1.0_dp)
                 t_new = TAIR + dil2_d*(dilut2_time)**dil2_f
                 t_new=max(t_new,TAIR)

            CASE (3)
                 ! T_new = T_old-( (T_old-T_fin) * (1/tau_c) *DT )
                 !
                 t_new = t_old - ( (t_old-T_fin)*(1/tau_c) * DT )
                 t_new=max(t_new,T_fin)

            CASE (4)
                 ! same as type 1
                 DRPL=zplum/zplum_old
                 t_new=((t_old-TAIR)/DRPL)+TAIR
                 t_new=MAX(t_new,TAIR)

            CASE (5)
                 ! same as type 1
                 DRPL=zplum/zplum_old
                 t_new=((t_old-TAIR)/DRPL)+TAIR
                 t_new=MAX(t_new,TAIR)

            CASE (6)
                 ! same as type 1
                 DRPL=zplum/zplum_old
                 t_new=((t_old-TAIR)/DRPL)+TAIR
                 t_new=MAX(t_new,TAIR)

            CASE (7)
                 ! same as type 1
                 DRPL=zplum/zplum_old
                 t_new=((t_old-TAIR)/DRPL)+TAIR
                 t_new=MAX(t_new,TAIR)

            CASE DEFAULT
 
              write(6,*) 'no valid plume dispersion option'
               STOP

         END SELECT

  end subroutine plumetemp


  subroutine initplume( u10,t_old,wplum,zplum ) 
    !----------------------------------------------------------------------
    !     
    !  Calculate the initial plume dispersion
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate the cross section area of the plume
    !      
    !      Type 5: DTdisp is the start time of plume dispersion
    !              in MOCCA: 10s, but here: 1s
    !
    !      interface
    !      ---------
    !
    !        input:
    !          u10              wind speed                      [s]
    !
    !        output:
    !          t_old            plume temperature, old timestep [K]
    !          wplum            plume width, new timestep       [m]
    !          zplum            plume height, new timestep      [m]
    !
    !      method
    !      ------
    !      MOCCA plume dispersion, Jana Moldanova, IVL, 2021
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
     real( dp), intent(in)                      :: u10         ![m/s]

! output
     real( dp), intent(out)                     :: t_old       ![K]
     real( dp), intent(out)                     :: wplum       ![m]
     real( dp), intent(out)                     :: zplum       ![m]
     
! local
! Type 5 parameters
     real( dp), parameter                       :: APARpld = 0.2
     real( dp), parameter                       :: BPARpld = 0.6
     real( dp), parameter                       :: TILL    = 3.0
     real( dp), parameter                       :: DTdisp  = 1.0
! Type 6 parameters
     real( dp), parameter                       :: plwith0 = 6.3 !m
     real( dp), parameter                       :: plhight0= 3.5 !m 

     real( dp)                                  :: dilfct
     real( dp)                                  :: PLwidth


! Plume initial temperature from dispers.dat
        t_old=ta_st

        SELECT CASE (IDIL)
            CASE (1)
                 ! Plume dispersion type 1
                 zplum=hmix_st
                 wplum=zplum

            CASE (2)
                 ! Plume dispersion type 2
                 zplum=hmix_st
                 wplum=zplum

            CASE (3)
                 ! Plume dispersion type 3
                 zplum=hmix_st
                 wplum=zplum

            CASE (4)
                 ! Plume dispersion type 4
                 zplum=hmix_st
                 wplum=zplum

            CASE (5)
                 ! Plume dispersion type 5
                 zplum=hmix_st
                 wplum= 1000. * APARpld *                            &    
                       EXP(BPARpld*LOG(((DTdisp)*u10+TILL)/1000.))

            CASE (6)
                 ! Plume dispersion type 6
                 zplum=plhight0
                 wplum=plwith0

            CASE (7)
                 ! Plume dispersion type 7
                 zplum=hmix_st
                 wplum= 1000. * APARpld *                            &    
                       EXP(BPARpld*LOG(((DTdisp)*u10+TILL)/1000.))

            CASE DEFAULT
 
              write(6,*) 'no valid plume dispersion option'
               STOP

         END SELECT

  end subroutine initplume


  subroutine plumearea(wplum,zplum,PLA) 
    !----------------------------------------------------------------------
    !     
    !  Calculate the cross section area of plume
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      For circular plume cross section (type 1,2,3,4):
    !      Apl = 0.5*pi * zplum**2
    !
    !      For semi-elliptic plume cross section (type 5-7)
    !      Apl = (pi/8) * wplum * zplum
    !
    !      interface
    !      ---------
    !
    !        input:
    !          wplum            plume width, new timestep       [m]
    !          zplum            plume height, new timestep      [m]
    !
    !        output:
    !          PLA              plume area                      [m^2]
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

  implicit none

! input
     real( dp), intent(in)                      :: wplum       ![m]
     real( dp), intent(in)                      :: zplum       ![m]

! output
     real( dp), intent(out)                     :: PLA         ![m^2]

! local


        SELECT CASE (IDIL)
            CASE (1)
                 ! Plume dispersion type 1
                 PLA = 0.5*pi * zplum**2

            CASE (2)
                 ! Plume dispersion type 2
                 PLA = 0.5*pi * zplum**2

            CASE (3)
                 ! Plume dispersion type 3
                 PLA = 0.5*pi * zplum**2

            CASE (4)
                 ! Plume dispersion type 4
                 PLA = 0.5*pi * zplum**2

            CASE (5)
                 ! Plume dispersion type 5
                 PLA = (pi/8) * wplum * zplum

            CASE (6)
                 ! Plume dispersion type 6
                 PLA = (pi/8) * wplum * zplum

            CASE (7)
                 ! Plume dispersion type 7
                 PLA = (pi/8) * wplum * zplum

            CASE DEFAULT
 
              write(6,*) 'no valid plume dispersion option'
               STOP

         END SELECT

  end subroutine plumearea


  subroutine wall_loss_h2so4(time_step_len,temp,ch2so4)
    !----------------------------------------------------------------------
    !     
    !      Calculation of wall losses for h2so4 in ageing chamber
    !      for diesel exhaust
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculation of wall loss of Cgas
    !      this routine is called from within the plume dispersion
    !      module (not available for chamber experiments)
    !      input: temperature in K
    !      input/output the gas phase concentration
    !
    !      interface
    !      ---------
    !
    !        input:
    !          temp             temperature                   [K]
    !          time_step_len    time step length              [s]
    !
    !      method
    !      ------
    !      Vouitsis, E., Ntziachristos, L., and Z. Samaras,
    !         Modelling of diesel exhaust aerosol during 
    !         laboratory sampling, Atmos. Env. 39,1335, 
    !         doi:10.1016/j.atmosenv.2004.11.011, 2005
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
    REAL( dp), intent(in)     :: time_step_len        ! [s]
    REAL( dp), intent(in)     :: temp                 ! [K]
    !
    REAL( dp), intent(in out) :: ch2so4

    REAL( dp) :: p,rho_air,vf,viscosity,te,psa
    REAL( dp) :: sigma_A
    REAL( dp) :: D,DC,Re,Sc,Sh,WL
    REAL( dp) :: fct

        fct = time_step_len
	
        p = 1._dp                          ! atm
        rho_air = 1.2928_dp*273.15_dp/temp ! kg/m3
        ! fluid velocity in the ageing chamber, 55 lpm
        D = 5.5_dp                         ! cm tube diameter
        vf = 0.386_dp                      ! m/s
        viscosity = 1.8e-5_dp              ! Pas for air

    ! Saturation vapour pressure at the wall surface for H2SO4
    ! (Vehkamäki et al., Environmental Science and Technology 37, 3392, 2003)
        te=295._dp                         ! K wall temperature
        psa = acidps(te)                   ! Pa

    ! diffusion coefficient (cm^2/s) for H2SO4
    !   sigma is the collision diameter, estimated from
    !   the molecular volume of the liquid, in Angstroem
    !   for sulfphuric acid:
    !   sigma_A= 19.7**(1/3)   in Angstroem
        sigma_A = 19.7_dp**(1.0_dp/3.0_dp)
    ! Value of DC should be around 0.04 cm2/s to 0.05 cm2/s
        DC = molecdiff(M_H2SO4,sigma_A,temp,p)
       
    ! write(6,*) "Diffusion coefficient", DC
        Re = rho_air*vf*D*1.e-2_dp/viscosity       ! Reynolds number
        Sc = viscosity/rho_air/DC*1.e4_dp          ! Schmidt number
        Sh = 0.0096_dp*Re**0.913_dp*Sc**0.346_dp   ! Scherwood number
    ! Wall loss WL in cm^-3 s^-1
    ! Voutsis et al., 2005, eq (15)
    ! Cgas(h2so4) in mlc/cm^3
        WL = 4._dp*DC*( ch2so4    - &
             (psa*1.e-6_dp/1.38e-23_dp/te) )/D**2._dp*Sh

        if (Re.ge.2000._dp) write(6,*) 'turbulent flow'
        if (WL.lt.0._dp) WL=0._dp

        ch2so4 = ch2so4 - fct*WL


  end subroutine wall_loss_h2so4


end module gde_plume
