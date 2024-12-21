! <module_fitdis_io.for - A component of FITAERO
!                 Fit utility for number size distributions>
!**********************************************************************!
!*
!**********************************************************************! 
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
!***********************************************************************! 
!*
!***********************************************************************
!***
!***      FITAERO
!***      Aerosol size distribution fit tool
!***      for MAFOR's inaero.dat
!***********************************************************************

      module module_fitaero_io

!***********************************************************************
!***  Module module_fitdis_io declares variables and parameters
!***  for input/output interfaces of FITAERO.
!***  Module contains all routines for reading and writing files
!***********************************************************************

      implicit none

!***********************************************************************

!     Declarations:

! file units input
      integer             :: funit_run
      integer             :: funit_log
      integer             :: funit_in_sizedis
      integer             :: funit_in_massfra

! file units output
      integer             :: funit_out_sizedis
      integer             :: funit_out_inaero

! file names

      character (len=256) :: fname_run
      character (len=256) :: fname_log
      character (len=256) :: fname_in_sizedis
      character (len=256) :: fname_in_massfra

      character (len=256) :: fname_outpath
      character (len=256) :: fname_out_sizedis
      character (len=256) :: fname_out_inaero

! file exist flags

      logical             :: fe_log
      logical             :: fe_in_sizedis
      logical             :: fe_in_massfra

! file formatted read

      character(len=1)    :: tab = ACHAR(9)  ! TAB delimiter


!     Routines and Functions:

! ************
      contains

!**********************************************************************

!      subroutine nxtdat(un,leof)

! The subroutine prepares for reading the next uncommented 
! line of data from file.

!***********************************************************************

!      implicit none

! Scalar arguments

!      integer             :: un
!      character(len=256)  :: TXTSTR
!      logical             :: leof

! UN     - Fileunit
! TXTSTR - Textstring
! LEOF   - If end of file then true else false 

! If fileunit is nonpositive then just return

!      IF (un .LE. 0) RETURN

! Fileunit is positive

!      leof = .FALSE.

! Read lines from file

!  100 CONTINUE
!      READ (un,1000,END=999) TXTSTR
!      IF (TXTSTR(1:1) .EQ. '*') GOTO 100
!      BACKSPACE(un)
!
!      RETURN
!
!  999 CONTINUE
!
!      leof = .TRUE.
!      return

! 1000 FORMAT (A256)

! End of subroutine nxtdat

!      end subroutine nxtdat

!***********************************************************************
! Functions

!***********************************************************************

      integer function get_file_unit ()

! The function returns a unit number that is not in use

!***********************************************************************

      implicit none

      integer   :: lu_max, lu, m
      integer   :: iostat
      character(len=256) :: fn
      logical   :: opened
!
!      m = lu_max
      if (m < 1) m = 97
! Start with unit 12 (10 is metadata, 11 is logfile)
!      m=12

      do lu = m,1,-1
       
         inquire (unit=lu, opened=opened, iostat=iostat, name=fn)
         if (iostat.ne.0) cycle
         if (.not.opened) exit

      end do
      print *,'next unit',lu

      get_file_unit = lu
      return

      end function get_file_unit

!***********************************************************************

! End of module module_fitaero_io

      end module module_fitaero_io
