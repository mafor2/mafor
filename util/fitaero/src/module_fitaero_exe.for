! <module_fitaero_exe.for - A component of FITAERO
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

      module module_fitaero_exe

!***********************************************************************
!***  Module module_fitdis_exe declares variables and parameters
!***  for controlling the program flow of FITDIS
!***********************************************************************

      implicit none

!***********************************************************************

!     Declarations:

!     User Input

! number of observations
      integer              :: n_obins
! IMAX, number of bins per mode (mafor)
      integer              :: imax
! DPMAX of size distribution, in [m] (mafor)
      real                 :: dpmax
! mass correction wanted (T/F) in inaero.dat
      logical              :: mcorr=.false.
! RH in first hour (dim.less)
      real                 :: rhi

! ----------------------------------------------------------------------

! number of variables  (x+y)
      integer, parameter   :: nv = 2

! number of parameters (gmd, sigma, mass)
      integer, parameter   :: np = 3
      integer, parameter   :: np1= np+1

! *** Aerosol parameters
      integer, parameter   :: nm = 5           ! number of modes
      integer, parameter   :: nc = 10          ! number of components
!   Aerosol Components from input, number of components: 10
      integer, parameter   :: SU = 1
      integer, parameter   :: OC = 2
      integer, parameter   :: AM = 3
      integer, parameter   :: NI = 4
      integer, parameter   :: MS = 5
      integer, parameter   :: SA = 6
      integer, parameter   :: XX = 7   !PBA
      integer, parameter   :: EC = 8
      integer, parameter   :: DU = 9
      integer, parameter   :: WA = 10
!  Aerosol modes
      integer, parameter   :: NU = 1
      integer, parameter   :: NA = 2
      integer, parameter   :: AI = 3
      integer, parameter   :: AS = 4
      integer, parameter   :: CS = 5

!   Physical constants
      real, parameter      :: pi    = 3.14159265358979323846
      real, parameter      :: conv  = 1.e12    ! convert kg to ng
      real, parameter      :: dpmin = 1.5e-9   ! m
      real, parameter      :: DENV  = 1770.0   ! density H2SO4 [kg/m3]
      real, parameter      :: DENOC = 1570.0   ! density organic [kg/m3]
      real, parameter      :: DENAM = 1300.0   ! density NH4 [kg/m3]
      real, parameter      :: DENNI = 1300.0   ! density nitrate [kg/m3]
      real, parameter      :: DENMS = 1770.0   ! density MSAp [kg/m3]
      real, parameter      :: DENSA = 2240.0   ! density seasalt [kg/m3]
      real, parameter      :: DENXX = 1150.0   ! density bioaero [kg/m3]
      real, parameter      :: DENEC = 1200.0   ! density soot [kg/m3]
      real, parameter      :: DENDU = 1400.0   ! density dust [kg/m3]



! *** Simplex parameters

! number of function calls
      integer              :: kount

! *** Simplex initial information
      integer, parameter   :: maxcnt1 = 500
! stop at 20 error evaluations for AI mode
      integer, parameter   :: maxcnt2 = 20
      real, parameter      :: errmin = 1.0e-03
      !scaling <0 expands, >0 shrinks, 0: no scaling
      real, parameter      :: scalf  = -0.1


! *** dimensional arrays
      real, dimension(np1)        :: erf
      real, dimension(nm,np)      :: beta
      real, dimension(nm,np)      :: betamax
      real, dimension(nm,np)      :: betamin
      real, dimension(nm,np1,np)  :: p
      real, dimension(nm,np)      :: c
      real, dimension(nm,np)      :: r

! Allocatable arrays:
      real, allocatable    :: in_array_dp(:)
      real, allocatable    :: in_array_dn(:)
      real, allocatable    :: in_array_nuf(:)
      real, allocatable    :: in_array_naf(:)
      real, allocatable    :: in_array_aif(:)
      real, allocatable    :: in_array_acf(:)
      real, allocatable    :: in_array_cof(:)

      real, allocatable    :: dataa(:,:)
      real, allocatable    :: nest(:,:)
      real, allocatable    :: numout(:,:)
      real, allocatable    :: dense(:,:)



! *** Density estimates
      real, dimension(nm), parameter :: dens   = (/            &
                           1000.0, 1200.0, 1400.0, 1600.0, 1500.0 /) !kg m^-3

! *** Three sets of parameter estimates
! *** For documentation see:
! *** ./docu/FITAERO_params.xls

! *** IF 2.0E-6 <= dpmax < 1.0E-5  [m]  (FINE + COARSE) ******************!

! *** Parameter (beta) estimates
      real, dimension(nm), parameter :: es1_gmd   = (/             &
                                        8.5, 25.0, 120.0, 270.0, 2800.0 /) 
      real, dimension(nm), parameter :: es1_sigma = (/             &
                                        1.40, 1.60, 1.90, 1.55, 1.55 /)
      real, dimension(nm), parameter :: es1_mm    = (/             &
                                        0.55, 100.0, 2400.0, 6.e3, 8.e3 /)

! *** Parameter (beta) upper limits
      real, dimension(nm), parameter :: max1_gmd   = (/             &
                                        13.0, 120.0, 200.0, 480.0, 5000.0 /)
      real, dimension(nm), parameter :: max1_sigma = (/             &
                                        1.60, 1.80, 2.15, 2.00, 2.00 /)
      real, dimension(nm), parameter :: max1_mm    = (/             &
              12.*es1_mm(1),10.*es1_mm(2),10.*es1_mm(3),10.*es1_mm(4),10.*es1_mm(5) /) 

! *** Parameter (beta) lower limits
      real, dimension(nm), parameter :: min1_gmd   = (/             &
                                        5.5, 12.0, 85.0, 250.0, 800.0 /) 
      real, dimension(nm), parameter :: min1_sigma = (/             &
                                        1.31, 1.45, 1.75, 1.45, 1.45 /)
      real, dimension(nm), parameter :: min1_mm    = (/             &
                                        0.01, 10.0, 350.0, 1000.0, 1000.0 /)


! *** IF dpmax <  2.0E-6 [m]  (FINE)       *******************************!

! *** Parameter (beta) estimates
      real, dimension(nm), parameter :: es2_gmd   = (/             &
                                        9.8, 20.0, 60.0, 120.0, 340.0 /)
      real, dimension(nm), parameter :: es2_sigma = (/             &
                                        1.35, 1.60,1.80, 1.55, 1.55 /)
      real, dimension(nm), parameter :: es2_mm    = (/             &
                                        0.004, 16.0, 250.0, 5.e3, 7.e3 /)

! *** Parameter (beta) upper limits
      real, dimension(nm), parameter :: max2_gmd   = (/             &
                                        11.0, 80.0, 140.0, 230.0, 600.0 /)
      real, dimension(nm), parameter :: max2_sigma = (/             &
                                        1.45, 1.80, 2.00, 2.05, 2.20 /)
      real, dimension(nm), parameter :: max2_mm    = (/             &
              100.*es2_mm(1),30.*es2_mm(2),60.*es2_mm(3),60.*es2_mm(4),10.*es2_mm(5) /) 

! *** Parameter (beta) lower limits
      real, dimension(nm), parameter :: min2_gmd   = (/             &
                                        6.0, 12.0, 55.0, 75.0, 270.0 /) 
      real, dimension(nm), parameter :: min2_sigma = (/             &
                                        1.30, 1.35, 1.70, 1.45, 1.45 /)
      real, dimension(nm), parameter :: min2_mm    = (/             &
                                        0.001, 13.0, 62.0, 140.0, 140.0 /)


! *** IF dpmax = 1.0E-5  [m] (PM10)        *******************************!

! *** Parameter (beta) estimates
      real, dimension(nm), parameter :: es3_gmd   = (/             &
                                        7.0, 25.0, 40.0, 250.0, 2400.0 /) 
      real, dimension(nm), parameter :: es3_sigma = (/             &
                                        1.40, 1.60, 1.90, 1.55, 1.55 /)
      real, dimension(nm), parameter :: es3_mm    = (/             &
                                        0.25, 20.0, 300.0, 1600.0, 2500.0 /)

! *** Parameter (beta) upper limits
      real, dimension(nm), parameter :: max3_gmd   = (/             &
                                        10.0, 100.0, 160.0, 580.0, 5000.0 /)
      real, dimension(nm), parameter :: max3_sigma = (/             &
                                        1.56, 1.80, 2.15, 2.00, 2.00 /)
      real, dimension(nm), parameter :: max3_mm    = (/             &
              10.*es3_mm(1),20.*es3_mm(2),5.*es3_mm(3),5.*es3_mm(4),15.*es3_mm(5) /) 

! *** Parameter (beta) lower limits
      real, dimension(nm), parameter :: min3_gmd   = (/             &
                                        5.3, 12.0,  40.0, 150.0, 1700.0 /) 
      real, dimension(nm), parameter :: min3_sigma = (/             &
                                        1.30, 1.45, 1.75, 1.45, 1.45 /)
      real, dimension(nm), parameter :: min3_mm    = (/             &
                                        0.02, 1.0, 40.0, 800.0, 1700.0 /)



!     Routines and Functions:

! ************
      contains

!***********************************************************************
! Subroutines

      subroutine errfunc(xbeta,ndata,ndens,nv,binlo,binup,m,kount,error)

! The subroutine returns the error function for the data set
!***********************************************************************

        implicit none

        integer, intent(in)     :: nv,m
        integer, intent(in out) :: kount

        integer,dimension(nm), intent(in)        :: binlo,binup
        real, dimension(nm,nc-1), intent(in)     :: ndens
        real, dimension(n_obins,nv), intent(in)  :: ndata
        real, dimension(nm,np), intent(in)       :: xbeta

        real, intent(out)       :: error

  ! Local
        integer          :: i
        integer          :: mini
        integer          :: maxi
        real             :: yobs
        real             :: ycalc
        real             :: resi
        real             :: yhat
        real             :: massbin
        real             :: masswbin 
        real             :: massmode
        real             :: masswmode
        real             :: swfmode
        real             :: shrmode
        real             :: sigma_rel
        real             :: sqsigma
        real             :: sigma_dn
        real             :: dnmax

        !print *,'lo up',m,binlo(m),binup(m)
        !print *,'errf xbeta',m,xbeta(m,:)

! Allow relative uncertainty of 10 % for observed dN
        sigma_rel=0.1
        error=0.0

!***********************************************************************
!* First loop to calculate water content in aerosol mode;             *!
!* MAFOR aerosol size distribution function:                          *!
!***********************************************************************

        mini=binlo(m)
        maxi=binup(m)
        massmode=0.0
        masswmode=0.0
        shrmode=1.0

        do 10 i= mini,maxi

          yobs = ndata(i,nv)

          call initaero(xbeta,ndata,ndens,shrmode,binlo,binup,m,i,  &
                        yhat,massbin,masswbin)

          massmode=massmode+massbin
          masswmode=masswmode+masswbin

 10     continue

        swfmode=(masswmode/(max(massmode,1.e-20)))**(0.3333)
        swfmode=min(swfmode,3.0)
        swfmode=max(swfmode,1.0)
        shrmode=1./swfmode

!*** Function to reduce the mass in mode 2 and 3
!*** correction of the shrinking factor for Aitken mode and Acc mode
!*** Changing shrmode has actually an effect on the mode width and mass
!*** MAybe correction depends on mass composition (fraction so4).
!*** If mf(SO4) < 0.15 increase, if larger decrease
!*** The sum of non-hygroscopic mass (=EC+DU) may also have an effect

        if ((m==NU).and.(rhi > 0.90)) then
          shrmode=shrmode*1.4
        else
        ! stena
          shrmode=shrmode*0.2
        endif

        if (m==NA) then
            ! init1_test
            if (ndens(NA,EC)>0.10) then
              if (rhi > 0.90) then
                shrmode=shrmode*0.70
              else if (rhi > 0.75) then
                shrmode=shrmode*1.50         ! RH=80%
              else if (rhi > 0.65) then
                shrmode=shrmode*1.50         ! bkgr2_test
              else if (rhi > 0.45) then
                shrmode=shrmode*1.40         ! RH=50%
              else
                shrmode=shrmode*1.50         ! RH=10%
              endif
            else if (ndens(NA,SU)==0.15) then
              if (rhi > 0.65) then
              ! stena
                shrmode=shrmode*1.46
              else
              ! aces_test
                shrmode=shrmode*1.50
              endif
            endif
        endif

        if (m==AI) then

            if (ndens(AI,SU)==0.15) then
              if (rhi > 0.65) then
              ! stena
                shrmode=shrmode*3.0
              else
              ! aces_test
                shrmode=shrmode*1.50
              endif
            endif

            if (ndens(AI,SU)<0.15) then
              if (rhi > 0.65) then
                shrmode=shrmode*1.08
              else if (rhi > 0.45) then
              ! traff1_test, bktr1_test
                shrmode=shrmode*0.97
              else
                shrmode=shrmode*0.97
              endif
            else
             !use DU in AI mode to divide between exhaust and background
             !exhaust AI: DU<0.1, EC>0.15
             ! EXHAUST
             ! bkgr1_test, 
             ! xprs1, xprs2
              if (ndens(AI,EC)>0.15) then
                if (rhi > 0.90) then
                  shrmode=shrmode*1.40
                else if (rhi > 0.85) then
                  shrmode=shrmode*1.20
                else if (rhi > 0.60) then
                  shrmode=shrmode*0.50       ! bkgr2_test
                else if (rhi > 0.45) then
                  shrmode=shrmode*1.10
                else
                  shrmode=shrmode*0.90
                endif
              ! BACKGROUND
              else if (ndens(AI,EC)>0.10) then
              ! init1_test
                if (rhi > 0.90) then
                  shrmode=shrmode*0.70
                else if (rhi > 0.75) then
                  shrmode=shrmode*1.45       ! RH=80%
                else if (rhi > 0.65) then
                  shrmode=shrmode*1.40
                else if (rhi > 0.45) then
                  shrmode=shrmode*1.45       ! RH=50%
                else
                  shrmode=shrmode*1.45       ! RH=10%
                endif
              else if (ndens(AI,EC)>0.05) then
              ! aces_test
                  shrmode=shrmode*1.30
              else  ! EC < 0.05
                if (rhi > 0.85) then
              ! amar2_test, acont_test
                  shrmode=shrmode*1.30
                else if (rhi > 0.75) then
              ! arctic_test
                  shrmode=shrmode*1.05
                endif
              endif
            endif

        endif

        if (m==AS) then


            if (ndens(AS,SU)<=0.15) then
               if (ndens(AS,EC)<0.15) then
                 if (rhi > 0.65) then
                 !stena_test
                   shrmode=shrmode*1.60
                 else if (rhi > 0.45) then
                 !aces_test
                   shrmode=shrmode*1.00
                 else
                   shrmode=shrmode*0.99
                 endif
               else ! EC>0.15
               ! bktr1_test, traff1_test
                  shrmode=shrmode*1.20
               endif
            endif

          ! use EC in AS mode to divide between exhaust and background
          ! BACKGROUND (high SU, low EC)
          ! 26.11.2020 MSK use SALT>0.1 in AS mode for marine
          !                and SULF>0.2 in AS mode for coastal
            if ((ndens(AS,SU)>0.15).and.   &
                (ndens(AS,EC)<0.15)   ) then
              if (rhi > 0.92) then
                if (ndens(AS,SA)>0.1) then
                  shrmode=shrmode*0.72       !marine
                else
                  shrmode=shrmode*1.55
                endif
              else if (rhi > 0.89) then
                if (ndens(AS,SA)>0.1) then
                  !amar2_test
                  shrmode=shrmode*0.91       ! marine
                else
                  if (ndens(AS,SU)>0.2) then
                  !acont_test
                    shrmode=shrmode*1.30     ! coastal
                  else
                   shrmode=shrmode*1.2
                  endif
                endif
              ! init1_test
              ! arctic_test (RH=87%)
              else if (rhi > 0.79) then
                shrmode=shrmode*1.00         ! RH=50%
              else if (rhi > 0.65) then
                shrmode=shrmode*0.98
              else if (rhi > 0.55) then
                shrmode=shrmode*1.05
              else if (rhi > 0.45) then
                shrmode=shrmode*1.00         ! RH=50%
              else if (rhi > 0.20) then
              ! xprs1, xprs2
                shrmode=shrmode*0.90
              else
                shrmode=shrmode*1.50         ! RH=10%
              endif
            endif
 
           ! EXHAUST (high SU, high EC)
           ! bkgr1_test, 
            if ((ndens(AS,SU)>0.15).and.   &
                (ndens(AS,EC)>=0.15)   ) then
              if (rhi > 0.89) then
                shrmode=shrmode*0.90
              else if (rhi > 0.60) then
                shrmode=shrmode*1.25         ! bkgr2_test
              else if (rhi > 0.45) then
                shrmode=shrmode*0.80
              else
                shrmode=shrmode*0.77
              endif
            endif

        endif


        if (m==CS) then
            if (rhi > 0.75) then
               shrmode=shrmode*1.10
            else if (rhi > 0.60) then
               shrmode=shrmode*1.08
            else if (rhi > 0.40) then
               shrmode=shrmode*1.20
            else if (rhi > 0.15) then
               shrmode=shrmode*1.00
            else
               shrmode=shrmode*1.00
            endif
        endif


        dnmax = maxval( ndata(mini:maxi,2) )


        do 20 i= mini,maxi

          yobs = ndata(i,nv)

!***********************************************************************
!* Change the next statment to change the function being fit          *!
!* MAFOR aerosol size distribution function:                          *!
!***********************************************************************

          call initaero(xbeta,ndata,ndens,shrmode,binlo,binup,m,i,  &
                        yhat,massbin,masswbin)

          ycalc=yhat


          ! Residual as diff OBS-MOD
          !resi=yobs-ycalc

          ! Chi2 Residual
          if (i.ge.maxi) then
             sigma_rel=sigma_rel*0.55
          endif
          sqsigma=yobs*yobs 
          sigma_dn=sqrt(sigma_rel*sigma_rel*sqsigma) 

          resi=(yobs-ycalc)/sigma_dn

          ! Control the fitting result
          !print *,i,yobs,ycalc,resi

          error=error+resi*resi

 20     continue

        kount=kount+1

        !print *,'massmode2',m,massmode,masswmode,swfmode,shrmode


        return

      end subroutine errfunc

!***********************************************************************
! Functions

!***********************************************************************

        ELEMENTAL REAL FUNCTION  poly4(A1,A2,A3,A4,X)
          IMPLICIT NONE
        REAL, intent(in)   :: A1,A2,A3,A4
        REAL, intent(in)   :: X
        poly4 = A1 + X *( A2 + X*( A3 + X*A4))
        END FUNCTION poly4

        ELEMENTAL REAL FUNCTION  poly6(A1,A2,A3,A4,A5,A6,X)
          IMPLICIT NONE
        REAL, intent(in)   :: A1,A2,A3,A4,A5,A6
        REAL, intent(in)   :: X
        poly6 = A1 + X*( A2 + X*( A3 + X*( A4 +   &
               X * ( A5 + X * (A6  )))))
        end function poly6

!    real function xxxx ()

! The function returns xxxx

!***********************************************************************

!      implicit none


!      return


!      end function xxxx

!***********************************************************************

! End of module module_fitaero_exe

      end module module_fitaero_exe
