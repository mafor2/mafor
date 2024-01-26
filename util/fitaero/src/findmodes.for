! <findmodes.for - A component of FITAERO
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

      subroutine findmodes(ndata,nobs,mstart,binlo,binup)

!***********************************************************************
!***  Subroutine findmodes finds the lowest and highest bin
!***  of the four modes in the observation data
!***********************************************************************

!     Declarations of variables by using the MODULES feature:

      use module_fitaero_exe


      implicit none

!***********************************************************************

       integer, intent(in)                  :: nobs
       integer, intent(in)                  :: mstart

       real, dimension(nobs,nv), intent(in) :: ndata

       integer, dimension(nm), intent(out)  :: binlo
       integer, dimension(nm), intent(out)  :: binup

!     Local declarations:

       integer                         :: i,j,m
       integer                         :: cm
       real                            :: qdpm
       real                            :: obsdp
       real                            :: obsdn
       real                            :: limait
       real                            :: limacc
       real, dimension(nm)             :: dplo
       real, dimension(nm)             :: dpup
       real, dimension(nobs-1)         :: diffd
       real, dimension(nm,imax)        :: dpmod

       logical :: loli2=.false.
       logical :: upli1=.false.
       logical :: upli2=.false.
       logical :: upli3=.false.

!***********************************************************************
!     Content of subroutine:
!
! BIN SIZE DISTRIBUTION IN MAFOR
! IMAX,DPMAX from User input
! initialise bin size distribution with one lognormal mode
! DPA is dry diameter [m]
!      QDPM=EXP(LOG(DPMAX/DPMIN)/(DBLE(IMAX*MMAX)-1.)) 
!      ! each mode to get IMAX bins
!      cm=0
!      do M=NU,CS
!        do I=1, IMAX
!         ! Bin diameter
!           DPA(M,I)=DPMIN*(QDPM**(DBLE(cm*IMAX+I)-1._dp))
!         ! Bin volume (dry)
!           VPT(M,I)=(pi/6._dp)*DPA(M,I)**3._dp
!        end do
!        cm=cm+1
!      end do
! THEN FIND CLOSEST BINUP
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! *** First guess initialisation of diameter limits

          ! Nucleation mode
          dplo(1) = dpmin*1.e9
          dpup(1) = 12.0
          binlo(1)= 1
          ! Aitken mode
          dplo(2) =  9.0
          dpup(2) = 20.0   !min
          ! Accumulation mode (guess)
          dplo(3) = 19.0 
          dpup(3) = 70.0   !min
          ! Coarse mode
          dplo(4) = 150.0
          dpup(4) = ndata(nobs,1)
          binup(4)= nobs


! *** Initialise bin size distribution with one lognormal mode
! *** Note DP is dry diameter [m] everywhere

          qdpm = exp( log(dpmax/dpmin)/(real(imax*nm)-1.)  )

! *** Calculate the model bin diameters
! *** Each mode to get IMAX bins
          cm=0
          do m=1,nm
            do i=1,imax
              dpmod(m,i)=dpmin*(qdpm**(real(cm*imax+i)-1.))
              dpmod(m,i)=dpmod(m,i)*1.e9
            end do
            cm=cm+1
          end do

          print *,'dpmod',dpmod(1,imax),dpmod(2,imax),dpmod(3,imax), &
                          dpmod(4,imax)

          ! Nucleation mode
          dplo(1) = dpmod(1,1)
          dpup(1) = dpmod(1,imax)
          ! Aitken mode
          dplo(2) = dpmod(2,1)
          dpup(2) = dpmod(2,imax)
          ! Accumulation mode
          dplo(3) = dpmod(3,1)
          dpup(3) = dpmod(3,imax)
          ! Coarse mode
          dplo(4) = dpmod(4,1)
          


! *** Find corresponding upper diameter in observed size distribution

          do j=2,nobs

              obsdp = ndata(j,1)

              ! determine dpup nucleation mode
              if ((.not.upli1).and.(mstart.eq.1)) then
                if ( obsdp .gt. dpup(1) ) then
                   dpup(1)=ndata(j-1,1)
                   binup(1)=j-1
                   binlo(2)=j
                   loli2=.true.
                   upli1=.true.
                endif
              endif

              ! determine dplo Aitken mode
              if ((.not.loli2).and.(mstart.eq.2)) then
                if ( obsdp .gt. dpup(1) ) then
                   dplo(2)=ndata(j-1,1)
                   binlo(2)=j-1
                   loli2=.true.
                endif
              endif

              ! determine dpup Aitken mode
              ! and dplo Accum mode
              if (.not.upli2) then
                if ( obsdp .gt. dpup(2) ) then
                   dpup(2)=ndata(j-1,1)
                   dplo(3)=ndata(j,1)
                   binup(2)=j-1
                   binlo(3)=j
                   upli2=.true.
                endif
              endif

              ! determine dpup Accum mode
              ! and dplo Coarse mode
              if (.not.upli3) then
                if ( obsdp .gt. dpup(3) ) then
                   dpup(3)=ndata(j-1,1)
                   dplo(4)=ndata(j,1)
                   binup(3)=j-1
                   binlo(4)=j
                   upli3=.true.
                endif
              endif

         enddo

         print *,'dplim(1)',dplo(1),dpup(1)
         print *,'dplim(2)',dplo(2),dpup(2)
         print *,'dplim(3)',dplo(3),dpup(3)
         print *,'dplim(4)',dplo(4),dpup(4)



      end subroutine findmodes
