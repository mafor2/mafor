! <gde_sensitiv.f90 - A component of the Multicomponent
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
module gde_sensitiv

implicit none

public :: readsens
public :: IDEPO,IWETD,ICOAG,ICOND,INUC,INUCMEC,ICHAM
public :: ICONS,ICONO,ICONA,ICONX,IPEMI,ICHEM,IKELV
public :: IDMS,IOZO,IORG,IDEB,IDIL,ISOA,INANO
public :: ICONW,IAQC,IAQP


integer,save               :: IDEPO,IWETD,ICOAG,ICOND,INUC,INUCMEC,ICHAM
integer,save               :: ICONS,ICONO,ICONA,ICONX,IPEMI,ICHEM,IKELV
integer,save               :: IDMS,IOZO,IORG,IDEB,IDIL,ISOA,INANO
integer,save               :: ICONW,IAQP,IAQC


contains

  subroutine readsens()

      open(9,file='sensitiv.dat',status='old')
      read(9,*) IDEPO,IWETD,ICOAG,ICOND,INUC,INUCMEC,ICHAM
      read(9,*) ICONS,ICONO,ICONA,ICONX,IPEMI,ICHEM,IKELV
      read(9,*) IDMS,IOZO,IORG,IDEB,IDIL,ISOA,INANO
      read(9,*) ICONW,IAQP,IAQC

      close(9)

  end subroutine readsens

!------------------------------------------------------------------

end module gde_sensitiv
