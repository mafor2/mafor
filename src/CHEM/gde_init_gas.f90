! <gde_init_gas.f90 - A component of the Multicomponent
!                     Aerosol Dynamics Model MAFOR>
!*****************************************************************************! 
!* 
!*    Copyright (C) 2011-2025  Matthias Steffen Karl
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
module gde_init_gas

   use messy_mecca_kpp_Parameters

   use messy_mecca_kpp_Global, only : ya_soan1, ya_soan2, fkoh_mea
   use messy_mecca_kpp_Global, only : wall
      
   use gde_input_data, only  : gspec
   use gde_input_data, only  : edms,eso2,eh2o2,cnh3

implicit none

   public :: readchem
   public :: readchamber
   public :: readmonit

! chamber
   public :: IAM, IOZ
   public :: F_HONO, KP_NIT, CAMI, K_DIL
   public :: L_MEA, L_NO2, L_HNO3, L_O3
   public :: fco, fcoj
   public :: DILPAR
   public :: V_CHAM,S_CHAM,S_SED,S_DIF
   public :: CWIN

   integer,save       :: IAM, IOZ
   real( dp),save     :: F_HONO, KP_NIT
   real( dp),save     :: K_DIL, CAMI
   real( dp),save     :: L_MEA,  L_NO2, L_HNO3, L_O3
   real( dp),save     :: fco,fcoj
   real( dp),save     :: DILPAR
   real( dp),save     :: V_CHAM,S_CHAM,S_SED,S_DIF
   real( dp),save     :: CWIN


contains

  subroutine readchem(emi,vdr,cg)
    !----------------------------------------------------------------------
    !     
    !   Read Chemistry input for simulation
    !   Initialize all tracers
    !   emis in (cm^-2s^-1), vdry in (cms^-1)
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Read Chemistry input for simulation
    !      Initialize all tracers
    !
    !      interface
    !      ---------
    !      inchem.dat
    !
    !      emis in (cm^-2s^-1), vdry in (cms^-1)
    !
    !      method
    !      ------
    !      read ascii file with space or tab separated entries
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

     real( dp),dimension(nspec), intent(out)    :: emi
     real( dp),dimension(nspec), intent(out)    :: vdr

     real( dp),dimension(nspec), intent(in out) :: cg


! local
     real( dp), dimension(gspec)     :: emisg,vdryg,cgasg

     integer               :: stat
     integer               :: gn
     integer               :: it,jk
! fix for tab reading
     character(len=180)    :: inchemline
     character(len=1)      :: tab = ACHAR(9)
     character(len=12)     :: dummyname

       open(20,file='inchem.dat',status='old', iostat=stat)
! open error handling    
       if (stat.ne.0) then
          write(6,*) 'File inchem.dat cannot be opened !'
          stop
       end if
       do jk=1,gspec
! fix for tab reading
         read(20,'(a)') inchemline
         do it=1,LEN(inchemline)
         ! replace each tab in the line with a space
           if( inchemline(it:it).EQ.tab) inchemline(it:it)=' '
         end do
         read(inchemline,*) dummyname, gn,cgasg(jk),emisg(jk),vdryg(jk)
! end fix for tab reading
       end do

       close(20)


!************     Initial Gas Phase (from inchem.dat)    ******************

! Emissions
	
       emi(ind_NO)      = emisg(1)
       emi(ind_NO2)     = emisg(2)
       emi(ind_HNO3)    = emisg(3)
       emi(ind_NH3)     = emisg(72)
       emi(ind_O3)      = emisg(6)
       emi(ind_H2SO4)   = emisg(5)
       emi(ind_CH4)     = emisg(7)
       emi(ind_HCOOH)   = emisg(8)
       emi(ind_HCHO)    = emisg(9)
       !!!emis(ind_H2O2)  = emisg(10)
       emi(ind_DMS)     = edms
       emi(ind_SO2)     = eso2
       emi(ind_H2O2)    = eh2o2
       emi(ind_CH3OOH)  = emisg(11)
       emi(ind_HONO)    = emisg(12)
       emi(ind_PAN)     = emisg(13)
       emi(ind_N2O5)    = emisg(14)
       emi(ind_HNO4)    = emisg(15)
       emi(ind_NO3)     = emisg(16)
       emi(ind_HCl)     = emisg(18)
       emi(ind_I2)      = emisg(31)
       emi(ind_CH3I)    = emisg(34)
       emi(ind_DMSO)    = emisg(38)
       emi(ind_CH3SO3H) = emisg(41)
       emi(ind_CO)      = emisg(42)
       emi(ind_C2H6)    = emisg(44)
       emi(ind_C2H4)    = emisg(45)
       emi(ind_C3H8)    = emisg(46)
       emi(ind_C3H6)    = emisg(47)
       emi(ind_NC4H10)  = emisg(48)
       emi(ind_MVK)     = emisg(49)
       emi(ind_MEK)     = emisg(50)
       emi(ind_C5H8)    = emisg(51)
       emi(ind_LTMB)    = emisg(52)
       emi(ind_TOLUENE) = emisg(53)
       emi(ind_LXYL)    = emisg(54)
       emi(ind_BPINENE) = emisg(55)
       emi(ind_CAMPHENE)= emisg(56)
       emi(ind_CARENE)  = emisg(57)
       emi(ind_SABINENE)= emisg(58)
       emi(ind_MMA)     = emisg(59)
       emi(ind_DMA)     = emisg(60)
       emi(ind_TMA)     = emisg(61)
       emi(ind_CH2NCH3) = emisg(62)
       emi(ind_MEA)     = emisg(63)
       emi(ind_DEA)     = emisg(64)
       emi(ind_TEA)     = emisg(65)
       emi(ind_AMP)     = emisg(66)
       emi(ind_TME)     = emisg(67)
       emi(ind_IPN)     = emisg(68)
       emi(ind_APINENE) = emisg(69)
       emi(ind_CHEX)    = emisg(70)
       emi(ind_SO3)     = emisg(71)
       ! (72) NH3 see above
       emi(ind_BSOV)    = emisg(73)
       emi(ind_BLOV)    = emisg(74)
       emi(ind_BELV)    = emisg(75)
       emi(ind_ASOV)    = emisg(76)
       emi(ind_ALOV)    = emisg(77)
       emi(ind_AELV)    = emisg(78)
       emi(ind_PIOV)    = emisg(79)
       emi(ind_PSOV)    = emisg(80)
       emi(ind_PELV)    = emisg(81)
    
    
! Dry Deposition

       vdr(ind_NO)      = vdryg(1)
       vdr(ind_NO2)     = vdryg(2)
       vdr(ind_HNO3)    = vdryg(3)
       vdr(ind_NH3)     = vdryg(72)
       vdr(ind_SO2)     = vdryg(4)
       vdr(ind_H2SO4)   = vdryg(5)
       vdr(ind_O3)      = vdryg(6)
       vdr(ind_CH4)     = vdryg(7)
       vdr(ind_HCOOH)   = vdryg(8)
       vdr(ind_HCHO)    = vdryg(9)
       vdr(ind_H2O2)    = vdryg(10)
       vdr(ind_CH3OOH)  = vdryg(11)
       vdr(ind_HONO)    = vdryg(12)
       vdr(ind_PAN)     = vdryg(13)
       vdr(ind_N2O5)    = vdryg(14)
       vdr(ind_HNO4)    = vdryg(15)
       vdr(ind_NO3)     = vdryg(16)
       vdr(ind_DMS)     = vdryg(17) 
       vdr(ind_HCl)     = vdryg(18)
       vdr(ind_I2)      = vdryg(31)
       vdr(ind_CH3I)    = vdryg(34)
       vdr(ind_DMSO)    = vdryg(38)
       vdr(ind_CH3SO3H) = vdryg(41)
       vdr(ind_CO)      = vdryg(42)
       vdr(ind_C2H6)    = vdryg(44)
       vdr(ind_C2H4)    = vdryg(45)
       vdr(ind_C3H8)    = vdryg(46)
       vdr(ind_C3H6)    = vdryg(47)
       vdr(ind_NC4H10)  = vdryg(48)
       vdr(ind_MVK)     = vdryg(49)
       vdr(ind_MEK)     = vdryg(50)
       vdr(ind_C5H8)    = vdryg(51)
       vdr(ind_LTMB)    = vdryg(52)
       vdr(ind_TOLUENE) = vdryg(53)
       vdr(ind_LXYL)    = vdryg(54)
       vdr(ind_BPINENE) = vdryg(55)
       vdr(ind_CAMPHENE)= vdryg(56)
       vdr(ind_CARENE)  = vdryg(57)
       vdr(ind_SABINENE)= vdryg(58)
       vdr(ind_MMA)     = vdryg(59)
       vdr(ind_DMA)     = vdryg(60)
       vdr(ind_TMA)     = vdryg(61)
       vdr(ind_CH2NCH3) = vdryg(62)
       vdr(ind_MEA)     = vdryg(63)
       vdr(ind_DEA)     = vdryg(64)
       vdr(ind_TEA)     = vdryg(65)
       vdr(ind_AMP)     = vdryg(66)
       vdr(ind_TME)     = vdryg(67)
       vdr(ind_IPN)     = vdryg(68)
       vdr(ind_APINENE) = vdryg(69)
       vdr(ind_CHEX)    = vdryg(70)
       vdr(ind_SO3)     = vdryg(71)
       ! (72) NH3 see above
       vdr(ind_BSOV)    = vdryg(73)
       vdr(ind_BLOV)    = vdryg(74)
       vdr(ind_BELV)    = vdryg(75)
       vdr(ind_ASOV)    = vdryg(76)
       vdr(ind_ALOV)    = vdryg(77)
       vdr(ind_AELV)    = vdryg(78)
       vdr(ind_PIOV)    = vdryg(79)
       vdr(ind_PSOV)    = vdryg(80)
       vdr(ind_PELV)    = vdryg(81)


! Initial Concentrations

       cg(ind_NO)       = cgasg(1)
       cg(ind_NO2)      = cgasg(2)
       cg(ind_HNO3)     = cgasg(3)
       ! NH3 from ingeod.dat
       cg(ind_NH3)      = cgasg(72)
       ! initial NH3 must not be zero
       cg(ind_NH3)      = max(1.0_dp,cg(ind_NH3))
       cg(ind_SO2)      = cgasg(4)
       cg(ind_H2SO4)    = cgasg(5)
       cg(ind_O3)       = cgasg(6)
       cg(ind_CH4)      = cgasg(7)
       cg(ind_HCOOH)    = cgasg(8)
       cg(ind_HCHO)     = cgasg(9)
       cg(ind_H2O2)     = cgasg(10)
       cg(ind_CH3OOH)   = cgasg(11)
       cg(ind_HONO)     = cgasg(12)
       cg(ind_PAN)      = cgasg(13)
       cg(ind_N2O5)     = cgasg(14)
       cg(ind_HNO4)     = cgasg(15)
       cg(ind_NO3)      = cgasg(16)
       cg(ind_DMS)      = cgasg(17)
       cg(ind_HCl)      = cgasg(18)
       cg(ind_I2)       = cgasg(31)
       cg(ind_CH3I)     = cgasg(34)
       cg(ind_DMSO)     = cgasg(38)
       cg(ind_CH3SO3H)  = cgasg(41)
       cg(ind_CO)       = cgasg(42)
       ! C2-C6 NMHC
       cg(ind_C2H6)     = cgasg(44)
       cg(ind_C2H4)     = cgasg(45)
       cg(ind_C3H8)     = cgasg(46)
       cg(ind_C3H6)     = cgasg(47)
       cg(ind_NC4H10)   = cgasg(48)
       cg(ind_MVK)      = cgasg(49)
       cg(ind_MEK)      = cgasg(50)
       ! aromatics & biogenics
       cg(ind_C5H8)     = cgasg(51)
       cg(ind_LTMB)     = cgasg(52)
       cg(ind_TOLUENE)  = cgasg(53)
       vdr(ind_LXYL)    = cgasg(54)
       cg(ind_BPINENE)  = cgasg(55)
       cg(ind_CAMPHENE) = cgasg(56)
       cg(ind_CARENE)   = cgasg(57)
       cg(ind_SABINENE) = cgasg(58)
       ! AMINES
       cg(ind_MMA)      = cgasg(59)
       cg(ind_DMA)      = cgasg(60)
       cg(ind_TMA)      = cgasg(61)
       cg(ind_CH2NCH3)  = cgasg(62)
       cg(ind_MEA)      = cgasg(63)
       cg(ind_DEA)      = cgasg(64)
       cg(ind_TEA)      = cgasg(65)
       cg(ind_AMP)      = cgasg(66)
	! NEW
       cg(ind_TME)      = cgasg(67)
       cg(ind_IPN)      = cgasg(68)
       cg(ind_APINENE)  = cgasg(69)
       cg(ind_CHEX)     = cgasg(70)
       cg(ind_SO3)      = cgasg(71)
       ! (72) NH3 see above
       cg(ind_BSOV)     = cgasg(73)
       cg(ind_BLOV)     = cgasg(74)
       cg(ind_BELV)     = cgasg(75)
       cg(ind_ASOV)     = cgasg(76)
       cg(ind_ALOV)     = cgasg(77)
       cg(ind_AELV)     = cgasg(78)
       cg(ind_PIOV)     = cgasg(79)
       cg(ind_PSOV)     = cgasg(80)
       cg(ind_PELV)     = cgasg(81)


       if ( cg(ind_SO3).gt.9.1e11 ) then
         write(6,*) 'STOP: initial SO3 too high (must be <= 9.0e11 cm^-3)'
         stop
       endif

  end subroutine readchem

!------------------------------------------------------------------

  subroutine readchamber()
    !----------------------------------------------------------------------
    !     
    !   Read chamber input for simulation
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Read chamber input for simulation
    !
    !      interface
    !      ---------
    !      incham.dat
    !
    !      CHAMBER DATA
    !        IOZ=1 read o3 from monitor.dat
    !        Amine-specific: IAM,KP_NIT,fco,fkoh_mea,CAMI
    !        K_DIL dilution rate for gases and particles
    !        IAM=1: MEA
    !        IAM=2: MMA
    !        IAM=3: DMA
    !        IAM=4: TMA
    !        IAM=5: MMI
    !        IAM=6: AMP
    !        Chamber-specific: V_CHAM,S_CHAM,S_SED,S_DIF,CWIN
    !
    !      method
    !      ------
    !      read ascii file with space or tab separated entries
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

! local
     integer               :: stat4
     integer               :: stat5


! read incham.dat

       open(23,file='incham.dat',status='old', iostat=stat4)
! open error handling    
       if (stat4.ne.0) then
          write(6,*) 'File incham.dat cannot be opened !'
          stop
       end if  
       read(23,*) IAM, IOZ,  KP_NIT, F_HONO, ya_soan1, ya_soan2, fco
       read(23,*) fkoh_mea, CAMI, K_DIL, L_MEA,  L_NO2, L_HNO3, L_O3
! chamber geometry
       read(23,*) DILPAR,V_CHAM,S_CHAM,S_SED,S_DIF,CWIN

       ! check molar yields
       if ((ya_soan1+ya_soan2).gt.1.0) then
         write(6,*) 'STOP: molar yield > 1.0 in incham.dat'
         stop      
       endif
       if (IAM.gt.5) then
         write(6,*) 'STOP: amine number must be <6 in incham.dat'
         stop
       endif
       if (V_CHAM.eq.0) then
         write(6,*) 'STOP: chamber volume  must be >0 m3 incham.dat'
         stop
       endif

       wall=1

       close(23)


! open monitor.dat
       open(24,file='monitor.dat',status='old', iostat=stat5)
! open error handling    
       if (stat5.ne.0) then
          write(6,*) 'File monitor.dat cannot be opened !'
          stop
       end if 


  end subroutine readchamber

!------------------------------------------------------------------

  subroutine readmonit(cco3,cno,cno2,cipn,temp,jno2m)
    !----------------------------------------------------------------------
    !     
    !   Read monitors input for simulation
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Read Chemistry input for simulation
    !      Initialize all tracers
    !
    !      interface
    !      ---------
    !      monitor.dat
    !
    !      CHAMBER DATA
    !      IOZ=1 read o3 from monitor.dat
    !      concentration (in ppb) times series
    !      from monitor data, called every 10 min
    !      scale calc J values to jno2 meas.
    !
    !      method
    !      ------
    !      read ascii file with space or tab separated entries
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

     real( dp), intent(out)                     :: cco3
     real( dp), intent(out)                     :: cno
     real( dp), intent(out)                     :: cno2
     real( dp), intent(out)                     :: cipn
     real( dp), intent(out)                     :: temp
     real( dp), intent(out)                     :: jno2m

      ! monitor.dat is opened in readcham()
      ! because read from monitor.dat happens continously

! read monitor.dat

       read(24,*) cco3, cno, cno2, jno2m, temp, cipn
       ! check for negative input values
       IF ( (cco3.LT.0.) .OR. (cno.LT.0.) .OR. (cno2.LT.0.)   &
            .OR. (jno2m.LT.0.) &
            .OR. (temp.LT.0.) .OR. (cipn.LT.0.) ) THEN
         write(6,*) 'STOP: negative input value in monitor.dat'
         stop
       ENDIF

       if (jno2m.LT.1.e-6) then
         fcoj=0.
       else
         fcoj=fco
       endif


  end subroutine readmonit

!------------------------------------------------------------------

end module gde_init_gas
