! <get_user_input.for - A component of FITAERO
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

      subroutine get_user_input

!***********************************************************************
!***  Subroutine get_user_input reads the meta information
!***********************************************************************

!     Declarations of variables by using the MODULES feature:

      use module_fitaero_io
      use module_fitaero_exe


      implicit none

!***********************************************************************

!     Local declarations:

      integer :: i
      integer :: word_len
      
!***********************************************************************
!     Content of subroutine:


!     Definition of some simple input and output files:

      funit_run  = 10
      fname_run  = 'userfit.inp'

      open (unit=funit_run, file=fname_run,status='old')

!     Start reading parameters from the meta info file:

!     Files containing input-data on size distribution:
      read(funit_run,*) fname_in_sizedis
!
      print *,fname_in_sizedis

!     Files containing input-data on mass fractions:
      read(funit_run,*) fname_in_massfra
!
      print *,fname_in_massfra

!     Files containing output-data:
      read(funit_run,*) fname_outpath

!     Log file
      read(funit_run,*) fname_log
      !print *,'fname_log ',fname_log

!     Open Log file
      funit_log  = 11
      
      fe_log     = .false.
      if (fname_log(1:1) /= ' ') then
        open (unit=funit_log, file=fname_log,status='unknown')
        fe_log = .true.
      end if

!     Check if files exists
      INQUIRE(FILE=fname_in_sizedis,EXIST = fe_in_sizedis)

      if (.not.fe_in_sizedis) then
        if (fe_log) then
          write(funit_log,'(1X,A)') 'Missing input file sizedis'
        endif
        call stopit_fitaero('Missing input file sizedis')
      endif

      INQUIRE(FILE=fname_in_massfra,EXIST = fe_in_massfra)
      if (.not.fe_in_massfra) then
        if (fe_log) then
          write(funit_log,'(1X,A)') 'Missing input file massfrac'
        endif
        call stopit_fitaero('Missing input file massfrac')
      endif

!     The number of hours to compute:
      read(funit_run,*) n_obins

!     DPMAX in meter, as in inaero.dat:
      read(funit_run,*) dpmax

!     IMAX in meter, as in inaero.dat:
      read(funit_run,*) imax

!     Mass-correction option in aero.dat (T/F):
      read(funit_run,*) mcorr

!     RH in percent, convert to dimensionless
      read(funit_run,*) rhi
      rhi = rhi *0.01
      if (rhi.gt.0.99) then
        call stopit_fitaero('RH must not be >99%')      
      endif

      close(funit_run)

! **********************************************************************

!     Write the Log file:

      if (fe_log) then
        write(funit_log,'(1X,A)')         &
    '************************************************************'

        word_len = LEN_TRIM(fname_run)
        write(funit_log,'(1X,3A)') &
    'Input parameters read from the RUN_FILE: "', &
    fname_run(1:word_len),'"'

        write(funit_log,'(1X,A)') &
    '************************************************************'

        write(funit_log,*)
        write(funit_log,'(1X,A)') 'Files containing INPUT-data: '
        write(funit_log,'(1X,A)') '**************************** '
        write(funit_log,*)

        word_len = LEN_TRIM(fname_log)
        write(funit_log,'(1X,A58,2A)')  &
    'The name of the applied LOG_FILE = "', &
    fname_log(1:word_len),'"'


        write(funit_log,*)
        write(funit_log,'(1X,A)') 'User-defined steering parameters: '
        write(funit_log,'(1X,A)') '********************************* '
        write(funit_log,*)

        
        write(funit_log,'(1X,A)') &
    'END of input parameters read from the user-supplied META info:'
        write(funit_log,'(1X,A)') &
    '************************************************************'
      end if



      return
      end subroutine get_user_input
