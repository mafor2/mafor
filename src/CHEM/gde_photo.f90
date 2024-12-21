! <gde_photo.f90 - A component of the Multicomponent
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
!*    routine photrates written by Rolf Sander and Hella Riede
!* 
!*****************************************************************************!
module gde_photo

  use messy_cmn_photol_mem
  use messy_mecca_kpp_parameters
  use messy_mecca_kpp_global
  use messy_main_tools

  private

  public :: photrates, photo_mbl, photo_dealloc

contains


  subroutine photrates(daynr,Lati,jr)
    !----------------------------------------------------------------------
    !     
    !****  photolysis rate calculation
    !
    !      author
    !      -------
    !      Rolf Sander,    MPICH, Mainz, 2003-2007
    !      Hella Riede,    MPICH, Mainz, 2007
    !
    !      Air Chemistry Department
    !      Max Planck Instiute of Chemistry (MPICH)
    !      P.O. Box 3060
    !      55020 Mainz, Germany
    !      Email: rolf.sander@mpic.de
    !
    !      purpose
    !      -------
    !      calculation of photolysis rates according to JVAL
    !
    !      interface
    !      ---------
    !      from messy_sappho.f90, part of CAABA/MECCA v4.0
    !      compatible with the MESSy standard, and the
    !      JVal PreProcessor (JVPP)
    !
    !      method
    !      ------
    !      Landgraf, J. and P.J. Crutzen
    !         An Efficient Method for Online Calculations of Photolysis
    !         and Heating Rates.
    !         Bull. Am. Met. Soc., 25, 863-878, 1998
    !
    !      J values from PAPER model by Landgraf et al.
    !      J value as a function of photon: J = A*exp(-B/(C+photon))
    !      photon=sin(PSI)=cos(THETA)
    !      THETA=zenith angle, PSI=90-THETA
    !      temp profile of atmos. is: data/prof.AFGL.midl.sum
    !      surface albedo :  0.00
    !      ozone column (Dobson) :300.00
    !      cloud cover: OFF
    !      paper model was started on: 00-01-18
    !
    !  
    !      reference
    !      ---------
    !      Landgraf, J. and P.J. Crutzen
    !      An Efficient Method for Online Calculations of Photolysis
    !      and Heating Rates.
    !      Bull. Am. Met. Soc., 25, 863-878, 1998
    !
    !      Sander, R., P. Joeckel, O. Kirner, A.T. Kunert, J. Landgraf,
    !      and A. Pozzer
    !      The photolysis module JVAL-14, compatible with the MESSy
    !      standard, and the JVal Preprocessor (JVPP).
    !      Geosci. Model Dev., 7, 2653–2662, 2014
    !      doi:10.5194/gmd-7-2653-2014
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    use gde_constants, only: pi,Cancer,InitSpr,OneDay

    implicit none

    ! local
    REAL( dp), intent(in)                   :: Lati, daynr
    REAL( dp),dimension(IP_MAX),intent(out) :: jr                             
    REAL( dp)                               :: Latitu
    REAL( dp)                               :: DayREAL
    REAL( dp)       :: SoDecli, SINPSI, PHOTON, FCT, DN  
    REAL( dp), PARAMETER :: DUSK = 0.0721347_dp   ! = 5.E-2/LOG(2.) = PHOTON at dusk

      ! day as a REAL( dp) value (at start of spring DayREAL( dp) = 0):
      ! InitSpr is JD 80., OneDay is 86400 s
      ! Day REAL( dp) is Julian Day - 80
      ! DayREAL( dp) = (model_time)/OneDay-InitSpr
      DayREAL = daynr - InitSpr
      !write(6,*) 'photo',DayREAL,Lati

      ! latitude as radiant
      Latitu=Lati * (pi/180._dp)   

      ! seasonal cycle:
      SoDecli = Cancer * sin(2._dp * pi * DayREAL / 365.25_dp)

      ! diurnal cycle of psi, the solar elevation angle
      SINPSI = sin(Latitu) * sin(SoDecli) &
           - cos(Latitu) * cos(SoDecli) * cos(2 * pi * DayREAL)

      ! PHOTON is approximately the positive part of SINPSI but
      !   - avoid sharp switching on and off
      !   - add some light at dawn before sunrise and at dusk after sunset
      PHOTON = DUSK*log(1._dp+exp(SINPSI/DUSK))
      ! write(6,*) 'photon', PHOTON

      ! DN is a day/night switch to set photolyses to about 0 at night
      FCT=50._dp
      DN = exp(FCT*SINPSI)/(exp(-FCT*SINPSI)+exp(FCT*SINPSI))

      jr(:)=0.0

      jr(ip_O1D   )= DN*4.4916E-04*exp(-3.5807E+00 /(3.2382E-01+PHOTON)) !J01
      jr(ip_O3P   )= DN*4.3917E-04*exp(-1.6665E-01 /(5.0375E-02+PHOTON)) !J02
      jr(ip_H2O2  )= DN*1.8182E-05*exp(-1.2565E+00 /(1.7904E-01+PHOTON)) !J03
      jr(ip_NO2   )= DN*1.2516E-02*exp(-4.5619E-01 /(6.9151E-02+PHOTON)) !J04
      jr(ip_NOO2  )= DN*2.2959E-02*exp(-1.0873E-01 /(2.1442E-02+PHOTON)) !J05
      jr(ip_NO2O  )= DN*1.9176E-01*exp(-1.2557E-01 /(1.9375E-02+PHOTON)) !J06
      jr(ip_N2O5  )= DN*1.1375E-04*exp(-1.1092E+00 /(1.7004E-01+PHOTON)) !J07
      jr(ip_HNO3  )= DN*3.4503E-06*exp(-2.1412E+00 /(2.6143E-01+PHOTON)) !J08
      jr(ip_CH3OOH)= DN*1.2858E-05*exp(-1.1739E+00 /(1.7044E-01+PHOTON)) !J09
      jr(ip_CHOH  )= DN*8.6933E-05*exp(-1.3293E+00 /(1.6405E-01+PHOTON)) !J10
      jr(ip_MGLYOX) = 2.*jr(ip_CHOH) ! based on code from Domenico      
      jr(ip_COH2  )= DN*8.9896E-05*exp(-8.4745E-01 /(1.2056E-01+PHOTON)) !J11
      jr(ip_HOCl  )= DN*4.9596E-04*exp(-8.4303E-01 /(1.3121E-01+PHOTON)) !J12
      jr(ip_Cl2O2 )= DN*2.8605E-03*exp(-8.8027E-01 /(1.4813E-01+PHOTON)) !J13
      jr(ip_ClNO2 )= DN*8.0989E-04*exp(-9.7319E-01 /(1.4656E-01+PHOTON)) !J14
      jr(ip_ClNO3 )= DN*7.4610E-05*exp(-7.2496E-01 /(1.2847E-01+PHOTON)) !J15
      jr(ip_Cl2   )= DN*3.9398E-03*exp(-6.3121E-01 /(1.0305E-01+PHOTON)) !J16
      jr(ip_HOBr  )= DN*3.2114E-03*exp(-4.6104E-01 /(8.9870E-02+PHOTON)) !J17
      jr(ip_BrNO2 )= DN*7.6992E-03*exp(-3.4144E-01 /(6.2243E-02+PHOTON)) !J18
      jr(ip_BrNO3 )= DN*1.9670E-03*exp(-5.2431E-01 /(9.9522E-02+PHOTON)) !J19
      jr(ip_Br2   )= DN*3.6899E-02*exp(-2.1097E-01 /(3.3855E-02+PHOTON)) !J20
      jr(ip_BrCl  )= DN*1.3289E-02*exp(-3.0734E-01 /(5.3134E-02+PHOTON)) !J21
      jr(ip_IO    )= DN*3.6977E-01*exp(-2.6523E-01 /(3.7541E-02+PHOTON)) !J22
      jr(ip_HOI   )= DN*1.2681E-02*exp(-3.7515E-01 /(6.1696E-02+PHOTON)) !J23
      jr(ip_INO3  )= DN*5.7985E-03*exp(-5.2214E-01 /(1.0423E-01+PHOTON)) !J24
      jr(ip_CH3I  )= DN*2.6536E-05*exp(-1.7937E+00 /(2.4260E-01+PHOTON)) !J25
      jr(ip_I2    )= DN*1.6544E-01*exp(-1.3475E-01 /(2.1859E-02+PHOTON)) !J26
      jr(ip_BrO   )= DN*6.9639E-02*exp(-7.4569E-01 /(1.0859E-01+PHOTON)) !J28
      jr(ip_ICl   )= DN*2.6074E-02*exp(-1.9215E-01 /(3.0177E-02+PHOTON)) !J29
      jr(ip_IBr   )= DN*7.4918E-02*exp(-1.5931E-01 /(2.3807E-02+PHOTON)) !J30
      jr(ip_INO2  )= DN*5.5440E-03*exp(-6.6060E-01 /(1.0034E-01+PHOTON)) !J31
      jr(ip_C3H7I )= DN*8.9493E-05*exp(-1.8899E+00 /(2.5853E-01+PHOTON)) !J32
      jr(ip_CH2ClI)= DN*4.9715E-04*exp(-1.4267E+00 /(2.0610E-01+PHOTON)) !J33
      jr(ip_CH2I2 )= DN*2.0184E-02*exp(-1.0349E+00 /(1.4961E-01+PHOTON)) !J34
      jr(ip_OClO  )= DN*1.0933E-01*exp(-4.4797E-01 /(7.2470E-02+PHOTON)) !J89
      jr(ip_HNO4  )= DN*2.1532E-05*exp(-1.9648E+00 /(2.1976E-01+PHOTON)) !J91
      jr(ip_HONO  )= DN*2.9165E-03*exp(-5.1317E-01 /(7.4940E-02+PHOTON)) !J92
      jr(ip_CH3Br )= DN*7.3959E-10*exp(-2.9326E+01 /(1.0132E-01+PHOTON)) !J99


   ! For use in Arctic and Antarctic scenarios 
   ! May be implemented later
   ! SUBROUTINE photo_ff
   !   ! quick and dirty conversion of J values by Landgraf et al.
   !   ! DAY8090 = ratio between first day of spring (JD80) and 1 April (JD90);
   !   DAY8090 = 0.674_DP;
   !   ! photon/COS(rad_lat) reaches about 1 on noon of first day of spring;
   !   ! PI/2. converts from 12h-mean to noon maximum;
   !
   !   jr(ip_O1D)    = 6.4697E-04 / EXP( 2.8200    /( 1.1388E-01+photon)) ! J01
   !   jr(ip_O3P)    = 7.2132E-04 / EXP( 1.7075    /( 1.9411E-01+photon)) ! J02
   !   jr(ip_H2O2)   = 5.3105E-05 / EXP( 1.3588    /( 1.8467E-01+photon)) ! J03
   !   jr(ip_NO2)    = 3.9861E-02 / EXP( 7.4565E-01/( 1.2747E-01+photon)) ! J04
   !   jr(ip_NOO2)   = 2.2000E-02 / EXP( 1.0000E-02/( 1.0000E-03+photon)) ! J05
   !   jr(ip_NO2O)   = 1.8000E-01 / EXP( 1.0000E-02/( 1.0000E-03+photon)) ! J06
   !   jr(ip_N2O5)   = 2.1621E-04 / EXP( 1.2521    /( 1.8418E-01+photon)) ! J07
   !   jr(ip_HNO3)   = 4.7367E-06 / EXP( 2.1397    /( 1.9219E-01+photon)) ! J08
   !   jr(ip_CH3OOH) = 4.7946E-05 / EXP( 1.3607    /( 1.8532E-01+photon)) ! J09
   !   jr(ip_CHOH)   = 2.1913E-04 / EXP( 1.4250    /( 1.7789E-01+photon)) ! J10
   !   jr(ip_COH2)   = 2.7881E-04 / EXP( 1.1487    /( 1.6675E-01+photon)) ! J11
   !   jr(ip_HOCl)   = 1.8428E-03 / EXP( 1.2575    /( 1.8110E-01+photon)) ! J12
   !   jr(ip_Cl2O2)  = 1.8302E-02 / EXP( 9.9560E-01/( 1.5859E-01+photon)) ! J13
   !   jr(ip_ClNO2)  = 2.4861E-03 / EXP( 1.2157    /( 1.7917E-01+photon)) ! J14
   !   jr(ip_ClNO3)  = 2.1778E-04 / EXP( 9.4755E-01/( 1.5296E-01+photon)) ! J15
   !   jr(ip_Cl2)    = 1.3150E-02 / EXP( 9.1210E-01/( 1.4929E-01+photon)) ! J16
   !   jr(ip_HOBr)   = 4.0575E-03 / EXP( 8.9913E-01/( 1.5034E-01+photon)) ! J17
   !   jr(ip_BrNO2)  = 5.0826E-02 / EXP( 6.5344E-01/( 1.1219E-01+photon)) ! J18
   !   jr(ip_BrNO3)  = 6.0737E-03 / EXP( 7.5901E-01/( 1.2785E-01+photon)) ! J19
   !   jr(ip_Br2)    = 6.9038E-02 / EXP( 4.9563E-01/( 8.3580E-02+photon)) ! J20
   !   jr(ip_BrCl)   = 3.4235E-02 / EXP( 6.3421E-01/( 1.0965E-01+photon)) ! J21
   !   jr(ip_IO)     = 1.2E-01 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J22
   !   jr(ip_HOI)    = 2.0E-03 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J23
   !   jr(ip_INO3)   = 1.4E-03 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J24
   !   jr(ip_CH3I)   = 9.1E-07 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J25
   !   jr(ip_I2)     = 6.1E-02 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J26
   !   jr(ip_BrO)    = 2.3319E-01 / EXP( 1.1094    /( 1.6736E-01+photon)) ! J28
   !   jr(ip_ICl)    = 9.3E-03 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J29
   !   jr(ip_IBr)    = 2.7E-02 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J30
   !   jr(ip_INO2)   = jx(ip_INO3) ! assumed by R.V.                      ! J31
   !   jr(ip_C3H7I)  = 2.8E-06 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J32
   !   jr(ip_CH2ClI) = 2.9E-05 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J33
   !   jr(ip_CH2I2)  = 2.1E-03 *PI/2. *DAY8090* photon/COS(rad_lat)       ! J34
   !   jr(ip_OClO)   = 3.2837E-01 / EXP( 6.9733E-01/( 1.1907E-01+photon)) ! J89
   !   jr(ip_HNO4)   = 6.0657E-05 / EXP( 1.5901    /( 1.8577E-01+photon)) ! J91
   !   jr(ip_HONO)   = 1.0218E-02 / EXP( 8.7261E-01/( 1.4850E-01+photon)) ! J92
   !   jr(ip_CH3Br)  = 0.                                                 ! J99
   !   jr(ip_CH3CHO) = jx(ip_CHOH)+jx(ip_COH2) ! mz_rs_20050911 assumed
   ! END SUBROUTINE photo_ff

  end subroutine photrates


  subroutine photo_mbl(firstloop,daynr,Lati,press,jr)
    !----------------------------------------------------------------------
    !     
    ! Calculate photolysis rates for all defined photolysis reactions (IP_MAX)
    !  in the marine boundary layer.
    !  Uses the routine JVAL from the PAPER model by Landgraf and Crutzen.
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculation of jr(IP_MAX)
    !      cossza=cos(THETA)=sin(PSI)
    !      input: number of day, latitude, air pressure
    !      output: jr(IP_MAX)
    !
    !      interface
    !      ---------
    !      from messy_jval_box.f90, part of CAABA/MECCA v4.0
    !      compatible with the MESSy standard, and the
    !      JVal PreProcessor (JVPP)
    !
    !      method
    !      ------
    !      interface to the MECCA model's JVAL module
    !      
    !      reference
    !      ---------
    !      Sander, R., P. Joeckel, O. Kirner, A.T. Kunert, J. Landgraf,
    !      and A. Pozzer
    !      The photolysis module JVAL-14, compatible with the MESSy
    !      standard, and the JVal Preprocessor (JVPP).
    !      Geosci. Model Dev., 7, 2653–2662, 2014
    !      doi:10.5194/gmd-7-2653-2014
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    use messy_jval,       only: jval_2d
    use messy_jval,       only: jvalues
    use messy_jval,       only: jval_solar_time_control
    !use messy_jval,       only: jval_read_nml_ctrl
    use messy_jval,       only: aerosol_data
    use messy_jval,       only: lp
    use messy_jval,       only: qy_CH3COCH3
    use messy_main_tools, only: nn_index
    use gde_constants,    only: pi,Cancer,InitSpr,OneDay
            
    implicit none

    ! input
    REAL( dp), intent(in)   :: Lati, daynr
    REAL( dp), intent(in)   :: press              ! in Pa  
    LOGICAL,   intent(in)   :: firstloop
    ! output
    REAL( dp),dimension(IP_MAX),intent(out) :: jr

    INTEGER, PARAMETER   :: nsza =  1
    INTEGER, PARAMETER   :: nlev = 19 ! number of pressure levels
    INTEGER, PARAMETER   :: jrow =  1 ! must be jrow=1 for box model
       
    REAL( dp),dimension(nsza,nlev+1)    :: v3
    REAL( dp),dimension(nsza,nlev+1)    :: relo3
    REAL( dp),dimension(nsza,nlev)      :: jpress     
    REAL( dp),dimension(nsza,nlev)      :: jrhum
    REAL( dp),dimension(nsza,nlev)      :: jtemp    
    REAL( dp),dimension(nsza,nlev)      :: aclc
    REAL( dp),dimension(nsza,nlev)      :: clp
    REAL( dp),dimension(nsza)           :: albedo
    REAL( dp),dimension(nsza)           :: slf    
    REAL( dp) :: cossza
    REAL( dp) :: Latitu
    REAL( dp) :: DayREAL
    REAL( dp) :: SoDecli

    LOGICAL   :: lmidatm
    LOGICAL   :: l_heating
    INTEGER   :: pbllev        ! number of levels in pbl
    INTEGER   :: jval_clev     ! current pressure level in jval
    INTEGER   :: status
    integer   :: ip
    integer   :: IAS1  
  
      albedo(:)  = 0.07_dp     ! DEFAULT_ALBEDO
      aclc(:,:)  = 0._dp       ! assume clear sky
      slf(:)     = 0._dp       ! 0 = sea
      ! clp = cloud liquid water path per layer [g/m^2]
      clp(:,:)   = 0._dp       ! assume clear sky
      lmidatm    = .FALSE.
      l_heating  = .FALSE.
      pbllev     = 5           ! number of levels in pbl

      ! GET COSSZA
      ! day as a REAL( dp) value (at start of spring DayREAL( dp) = 0):
      ! InitSpr is JD 80., OneDay is 86400 s
      ! Day REAL( dp) is Julian Day - 80
      ! DayREAL( dp) = (model_time)/OneDay-InitSpr
      DayREAL = daynr - InitSpr
      ! latitude as radiant
      Latitu=Lati * (pi/180._dp)   
      ! seasonal cycle:
      SoDecli = Cancer * sin(2._dp * pi * DayREAL / 365.25_dp)
      ! diurnal cycle of psi, the solar elevation angle
      cossza = sin(Latitu) * sin(SoDecli) &
           - cos(Latitu) * cos(SoDecli) * cos(2 * pi * DayREAL)


      ! array of photolysis rates with 3 dims
      IF (firstloop) THEN
        ALLOCATE(jval_2d(IP_MAX))
        DO ip = 1, IP_MAX
          ALLOCATE(jval_2d(ip)%ptr(nsza,nlev), STAT=IAS1)
          jval_2d(ip)%ptr(:,:) = 0.
        END DO
      ENDIF

        ! GLOBAL COLUMNS
        ! global average values are extracted with ferret from messy and
        ! jval_diag streams using e.g.: "list jrhum[i=@ave,j=@ave,l=1]"
        ! nlev = 19
        ! array values from DEFAULT_V3 etc. in caaba_mem.f90
        !
        ! vertical ozone column [mcl/cm2]
        v3(1,:)    = (/ &
          3.366E+17, 1.437E+18, 4.085E+18, 5.428E+18, 6.157E+18, 6.583E+18, &
          6.860E+18, 7.070E+18, 7.227E+18, 7.343E+18, 7.436E+18, 7.523E+18, &
          7.605E+18, 7.678E+18, 7.740E+18, 7.788E+18, 7.822E+18, 7.844E+18, &
          7.857E+18, 7.862E+18 /)
        ! relative ozone, i.e. ozone mixing ratio [mol/mol]
        relo3(1,:) = (/ &
          7.182E-06, 8.319E-06, 4.172E-06, 2.041E-06, 9.525E-07, 4.334E-07, &
          2.571E-07, 1.514E-07, 9.760E-08, 5.775E-08, 5.064E-08, 4.394E-08, &
          3.980E-08, 3.636E-08, 3.209E-08, 2.807E-08, 2.479E-08, 2.242E-08, &
          2.105E-08, 2.065E-08 /)
        ! pressure [Pa]
        jpress(1,:) = (/ &
          1000., 3000., 5040., 7339., 10248., 14053., 18935., 24966., 32107., &
          40212., 49027., 58204., 67317., 75897., 83472., 89631., 94099.,     &
          96838., 98169. /)
        ! relative humidity [%]
        jrhum(1,:)  = (/ &
          0.23, 1.57, 3.52, 11.73, 24.55, 25.31, 27.45, 36.46, 44.52, 46.27,  &
          46.48, 49.18, 51.73, 57.95, 72.82, 80.71, 81.66, 77.65, 76.18 /)
        ! temperature [K]
        jtemp(1,:)  = (/ &
          230.6, 218.2, 211.7, 207.0, 205.6, 210.9, 218.1, 225.8, 235.7, 246.6, &
          256.4, 264.2, 270.6, 275.4, 278.2, 280.9, 283.2, 284.9, 285.7 /)

      ! calculate pressure level in jpress according to current pressure
      CALL nn_index(jpress(1,:), press, jval_clev)     

      !! read jval ctrl namelist: jval.nml
      !! defines r_sol and qy_CH3COCH3
      ! CALL jval_read_nml_ctrl(status, 999)
      ! IF (status /= 0) STOP 1

      ! use r_sol [0,...1] in CTRL for solar cycle
      ! orbital parameter is set to 1.0 AU here (no orbital variation)
      ! no external solar cycle data provided here
      CALL jval_solar_time_control(status, 1.0_dp)
 
      ! QUANTUM YIELD FOR CH3COCH3:
      !    qy_ch3coch3 = 1 ! Gierzack & ECHAM5 (old IUPAC) (default)
      !    qy_ch3coch3 = 2 ! BLITZ 2004
      !    qy_ch3coch3 = 3 ! IUPAC
      qy_CH3COCH3 = 2


      ! intialize aerosol data
      CALL aerosol_data ! aerosol optical data (Shettle, Fenn)
      
      ! retrieve jval_gp for all photolysis reactions
      lp(:)=.TRUE.

      !write(6,*) 'gde_photo',status,jval_clev

      ! calculate jvalues
      ! messy_jval.f90 wants REAL, not REAL(DP)
      CALL jvalues(                                    &
        REAL(v3), REAL((/ cossza /)),                  &
        REAL(jpress), REAL(relo3),                     &
        REAL(jrhum), REAL(jtemp), REAL(albedo),        &
        REAL(aclc), REAL(slf), REAL(clp),              &
        lmidatm, l_heating, pbllev)

 
      do ip = 1, IP_MAX
        jr(ip) = jval_2d(ip)%ptr(1,jval_clev)
        !write(6,*) 'JR',ip,jr(ip)
      end do



  end subroutine photo_mbl

  subroutine photo_dealloc()
    !----------------------------------------------------------------------
    !     
    !   Deallocate the array of photolysis frequencies
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      deallocate the array of photolysis frequencies
    !
    !      interface
    !      ---------
    !      none
    !
    !      method
    !      ------
    !      free memory
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

    use messy_jval,       only: jval_2d
    implicit none
    
    integer   :: ip
    
  
     do ip = 1, IP_MAX
       deallocate(jval_2d(ip)%ptr)
       nullify(jval_2d(ip)%ptr)
     end do
     deallocate(jval_2d)
     nullify(jval_2d)

  end subroutine photo_dealloc

 end module gde_photo
