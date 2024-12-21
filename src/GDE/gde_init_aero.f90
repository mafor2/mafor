! <gde_init_aero.f90 - A component of the Multicomponent
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
module gde_init_aero

   use gde_input_data, only  : MMAX
   use gde_input_data, only  : NU,NA,AI,AS,CS
   use gde_input_data, only  : iamax
   use gde_input_data, only  : SU,OC,AM,NI,MS,SA,XX,EC,DU


implicit none

   INTRINSIC :: SELECTED_REAL_KIND

   public :: readinaero
   public :: reademitpar
   public :: readbackgr

   public :: DPMMIN, GMD, SIG, NUM
   public :: BGNO, BGNO2, BGSO2, BGO3
   public :: BGNH3, BGSULF
   public :: BGPIOV, BGPSOV
   public :: BGMD, BGSIG
   public :: EGMD, ESIG

! KPP DP - Double precision kind
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)

   real( dp),save     :: DPMMIN(NU:CS)
   real( dp),save     :: GMD(NU:CS)
   real( dp),save     :: SIG(NU:CS)
   real( dp),save     :: NUM(NU:CS)
   real( dp),save     :: BGMD(NU:CS)
   real( dp),save     :: BGSIG(NU:CS)
   real( dp),save     :: BGNO,BGNO2,BGSO2,BGO3
   real( dp),save     :: BGNH3,BGSULF
   real( dp),save     :: BGPIOV,BGPSOV
   real( dp),save     :: EGMD(NU:CS)
   real( dp),save     :: ESIG(NU:CS)


contains

  subroutine readinaero(IMAX,DPMAX,ndis,MSULFTOT,MORGCTOT,        &
                        MAMMOTOT,MNITRTOT,MMSAPTOT,MSALTTOT,      &
                        MXXXXTOT,MECBCTOT,MDUSTTOT )

!-------------------------
!---  inaero.dat
!-------------------------

   !!! Read Aerosol input
   !!! The modal structure is
   !!!    NU      NA       AI        AS          CS
   !!! 1) 1-10nm, 10-30nm, 30-100nm, 100-1000nm, 1000-10000nm or
   !!! 2) 1-5nm , 5-20nm,  20-50nm,  50-100nm,   100-1000nm

  implicit none

     integer, intent(out)                       :: IMAX
     logical, intent(out)                       :: ndis

     real( dp), intent(out)                     :: DPMAX
     real( dp), intent(out)                     :: MSULFTOT(NU:CS)
     real( dp), intent(out)                     :: MORGCTOT(NU:CS)
     real( dp), intent(out)                     :: MAMMOTOT(NU:CS)
     real( dp), intent(out)                     :: MNITRTOT(NU:CS)
     real( dp), intent(out)                     :: MMSAPTOT(NU:CS)
     real( dp), intent(out)                     :: MSALTTOT(NU:CS)
     real( dp), intent(out)                     :: MXXXXTOT(NU:CS)
     real( dp), intent(out)                     :: MECBCTOT(NU:CS)
     real( dp), intent(out)                     :: MDUSTTOT(NU:CS)

! local
     integer               :: stat
     integer               :: M

!! fix for tab reading
     character(len=180)    :: inchemline
     character(len=1)      :: tab = ACHAR(9)
     integer               :: it



       open(21,file='inaero.dat',status='old', iostat=stat)
! open error handling    
       if (stat.ne.0) then
          write(6,*) 'File inaero.dat cannot be opened !'
          stop
       end if

       read(21,*) DPMAX, IMAX

       write(6,*) 'IMAX',IMAX

! Minimum mode GMD [m]
        DPMMIN(NU)=1.5e-9_dp    
        IF (DPMAX.GE.1.0e-5_dp) THEN
         DPMMIN(NA)=9.0e-9_dp
         DPMMIN(AI)=3.0E-8_dp
         DPMMIN(AS)=1.0e-7_dp
         DPMMIN(CS)=1.5e-6_dp
       ELSE  !DPMAX=1.e-6
         DPMMIN(NA)=9.0e-9_dp
         DPMMIN(AI)=2.0E-8_dp
         DPMMIN(AS)=3.1e-8_dp
         DPMMIN(CS)=1.0e-7_dp
       ENDIF

       do M=NU,CS

! fix for tab reading
          read(21,'(a)') inchemline
          do it=1,LEN(inchemline)
          ! replace each tab in the line with a space
            if( inchemline(it:it).EQ.tab) inchemline(it:it)=' '
          end do
          read(inchemline,*) ndis,GMD(M),SIG(M),NUM(M),               &
                     MSULFTOT(M),MORGCTOT(M),MAMMOTOT(M),MNITRTOT(M), &
                     MMSAPTOT(M),MSALTTOT(M),MXXXXTOT(M),MECBCTOT(M), &
                     MDUSTTOT(M)
! end fix for tab reading

! check the lognormal distribution parameters
          if (SIG(M).lt.1.09) then
            write(6,*) 'STOP SIGMA too small in inaero.dat. Mode: ',M
            stop
          endif

          if (SIG(M).gt.2.20) then
            write(6,*) 'STOP SIGMA >2.2 not accepted in inaero.dat. Mode: ',M
            stop
          endif

          if (GMD(M).lt.DPMMIN(M)) then
            write(6,*) 'STOP GMD too small in inaero.dat. Mode: ',M
            stop
          endif

       end do


       close(21)


  end subroutine readinaero

!------------------------------------------------------------------

  subroutine reademitpar(EMMCTOT )

!-------------------------
!---  emitpar.dat
!-------------------------

    !!! read Particle emission input
    !!! Particle emission rate (ng m^-2 s^-1)
    !!! contiuous emission during the simulation

  implicit none


     real( dp), dimension(MMAX,iamax), intent(out) :: EMMCTOT



! local
     integer               :: stat
     integer               :: M

!! fix for tab reading
     character(len=180)    :: inchemline
     character(len=1)      :: tab = ACHAR(9)
     integer               :: it



       open(25,file='emitpar.dat',status='old', iostat=stat)
! open error handling    
       if (stat.ne.0) then
          write(6,*) 'File emitpar.dat cannot be opened !'
          stop
       end if

       do M=NU,CS
! fix for tab reading
         read(25,'(a)') inchemline
         do it=1,LEN(inchemline)
           ! replace each tab in the line with a space
           if( inchemline(it:it).EQ.tab) inchemline(it:it)=' '
         end do
         read(inchemline,*) EGMD(M), ESIG(M),                       & 
                EMMCTOT(M,SU),EMMCTOT(M,OC),EMMCTOT(M,AM),EMMCTOT(M,NI),  &
                EMMCTOT(M,MS),EMMCTOT(M,SA),EMMCTOT(M,XX),EMMCTOT(M,EC),  &
                EMMCTOT(M,DU)

! end fix for tab reading
       end do



       close(25)


  end subroutine reademitpar

!------------------------------------------------------------------

  subroutine readbackgr(BGMCTOT)

!-------------------------
!---  inbgair.dat
!-------------------------

    !!! Read Background Aerosol input +
    !!! Read background gas concentration (molec/cm^3)

  implicit none

     real( dp), dimension(MMAX,iamax), intent(out) :: BGMCTOT


! local
     integer               :: stat
     integer               :: M

!! fix for tab reading
     character(len=180)    :: inchemline
     character(len=1)      :: tab = ACHAR(9)
     integer               :: it



       open(26,file='inbgair.dat',status='old', iostat=stat)
! open error handling    
       if (stat.ne.0) then
          write(6,*) 'File inbgair.dat cannot be opened !'
          stop
       end if

       do M=NU,CS
! fix for tab reading
        read(26,'(a)') inchemline
        do it=1,LEN(inchemline)
          ! replace each tab in the line with a space
          if( inchemline(it:it).EQ.tab) inchemline(it:it)=' '
        end do
        read(inchemline,*) BGMD(M),BGSIG(M),                     &
                BGMCTOT(M,SU),BGMCTOT(M,OC),BGMCTOT(M,AM),BGMCTOT(M,NI),  &
                BGMCTOT(M,MS),BGMCTOT(M,SA),BGMCTOT(M,XX),BGMCTOT(M,EC),  &
                BGMCTOT(M,DU)

! end fix for tab reading
       end do
       read(26,'(a)') inchemline
        do it=1,LEN(inchemline)
          ! replace each tab in the line with a space
          if( inchemline(it:it).EQ.tab) inchemline(it:it)=' '
        end do
       read(inchemline,*) BGNO,BGNO2,BGSO2,BGO3,BGNH3,BGSULF,BGPIOV,BGPSOV


       close(26)


  end subroutine readbackgr

!------------------------------------------------------------------




!------------------------------------------------------------------

end module gde_init_aero
