! <initaero.for - A component of FITAERO
!                 Fit utility for number size distributions>
!**********************************************************************!
!*
!**********************************************************************! 
!* 
!*    Copyright (C) 2011-2020  Matthias Steffen Karl
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

      subroutine initaero(xbeta,ndata,ndens,shrinkf,binlo,binup,m,   &
                          ibin,yhat,massbin,masswbin)

!***********************************************************************
!***  Subroutine initaero calculates number size distribution
!***  This is the user-supplied function for the minimization
!***  problem
!***********************************************************************

!     Declarations of variables by using the MODULES feature:

      use module_fitaero_exe


      implicit none

!***********************************************************************

      integer, intent(in)                     :: ibin
      integer, intent(in)                     :: m
      integer,dimension(nm), intent(in)       :: binlo,binup
      real, dimension(n_obins,nv), intent(in) :: ndata
      real, dimension(nm,np), intent(in)      :: xbeta
      real, dimension(nm,nc-1), intent(in)    :: ndens
      real, intent(in)                        :: shrinkf

      real, intent(out)                       :: yhat
      real, intent(out)                       :: massbin
      real, intent(out)                       :: masswbin

!     Local declarations:
      integer              :: i,o
      integer              :: nobs
      real, dimension(nc)  :: massc
      real                 :: logmod,vrat
      real                 :: massdry
      real                 :: dpmini,dpmaxi
      real                 :: betagmd
      real                 :: dpobs
      real                 :: densf
      real                 :: watfr
      real, parameter      :: exmax    = -70.0


!***********************************************************************
!     Content of subroutine:
!
!     0) include basic constants from gde_input_data (use aerosol_data)
!     1) setup the mafor diameter grid for a mode -> DLOGDP(M,I)
!     2) estimate mass in mode
!     3) calculate number per bin
!     5) calculate LOGDIS(M,I)=N(M,I)/DLOGDP(M,I)
!     6) find dN/dlogDp for each observed Dp
!
! The aerosol size distribution for one mode
! BETA(1)  : GMD_M (nm)
! BETA(2)  : SIGMA
! BETA(3)  : MASS  (ng/m^3)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! initialisation
          vrat=1. 
          yhat=1.

        ! Volume ratio per aerosol mode
        ! VRAT(M)=(DPA(M,IMAX)/DPA(M,1))**(3./(IMAX-1))
        ! dpmin = ndata(1,1)       min-bin of mode
        ! dpmax = ndata(nobs,1)    max-bin of mode

          dpmini = ndata(binlo(m),1)   ! nm lower bin of mode
          dpmaxi = ndata(binup(m),1)   ! nm upper bin of mode
          nobs   = binup(m)-binlo(m)+1

        ! if only one data point in nucleation mode
          if ((nobs==1).and.(m==1)) then
              dpmaxi = 10.0
              dpmini = dpmin*1.e9
              nobs = 2
          endif

          vrat=(dpmaxi/dpmini)**(3./(nobs-1));


        ! Weighted density per aerosol mode
        ! H2SO4 ORG NH4 NO3 MSAp SALT BPAB BC DUST Total
          densf = ndens(m,SU)*DENV  + ndens(m,OC)*DENOC  &
                + ndens(m,AM)*DENAM + ndens(m,NI)*DENNI  &
                + ndens(m,MS)*DENMS + ndens(m,SA)*DENSA  &
                + ndens(m,XX)*DENXX + ndens(m,EC)*DENEC  &
                + ndens(m,DU)*DENDU

       !print *,'aero: vrat,dpmin,dpmax,nobs',m,vrat,dpmini,dpmaxi,nobs
       !print *,'aero: beta',xbeta(m,1),xbeta(m,2),xbeta(m,3)


        ! First we convert beta(1), so GMD is in [m]
        ! and observed diameter also to [m]
          betagmd = xbeta(m,1)*1.e-9        ! nm->m
          dpobs   = ndata(ibin,1)*1.e-9     ! nm->m

        ! Calculate lindism [ng/m^3/m]
          yhat = xbeta(m,3)                                  &
                /( dpobs * sqrt(2.*pi) * log(xbeta(m,2)) )
          !print*,'yhat0',yhat
          logmod = ((-0.5)*(log( dpobs/betagmd )  &
                    /log(xbeta(m,2)) )**2.0)
                
          !print*,'logmod',logmod
          if (logmod .lt. exmax) then
            yhat = 0.0
          else
            yhat = yhat * exp(logmod)
          endif

          !print*,'yhat1',yhat

        ! calculate mass = lindism * dlindp                
        ! dlindp :: Linear width of bin [in m]
        ! mass will be in units [ng/m^3]
          yhat = yhat * dpobs*(2.**(1./3.))                   &
                *(vrat**(1./3.)-1.) / ((1.+vrat)**(1./3.))    

        ! assuming that wet diameter was measured (dpobs)
        ! then masswbin is mass including water (MTOTW)
          !print*,'mass in bin',yhat
          masswbin=yhat


        ! estimate the water content in mass-% from equation
        ! based on Bian et al.(Atmos. Chem. Phys., 14, 6417–6426, 2014)
        ! and calcALWC.xls: lw(%)=87.7705*RH**2.2677
        !  watfr = 0.01*87.7705*rhi**(2.2677)
        ! it shows that LWC is less sensitive to RH
          watfr = 0.01*87.77*rhi**(2.8)


          massdry=masswbin*(1.-watfr)

        ! derive component masses using input mass fractions
          do o=1,nc-1
            massc(o)= ndens(m,o)*massdry
          enddo

        ! calculate water content in this bin
        ! relative humidity from user input
          call awater(rhi,massc(SU),massc(AM),massc(NI),massc(SA), &
                      massc(OC),massc(WA) )

          if (rhi.lt.0.05) massc(WA) = 0.0           

        ! update of dry mass
          massdry = masswbin - massc(WA)

        ! calculate water content again
          do o=1,nc-1
            massc(o)= ndens(m,o)*massdry
          enddo

          call awater(rhi,massc(SU),massc(AM),massc(NI),massc(SA), &
                      massc(OC),massc(WA) )

          if (rhi.lt.0.05) massc(WA) = 0.0  

        ! return total dry mass in bin
          massbin = masswbin - massc(WA)

        ! calculate vconc = mass/dens*CONV
        ! based on density of dry aerosol
        ! dens in [kg/m^3]
        ! conv convert kg to ng --> [ng/m^3]
        ! dimensionless [-]
          yhat = yhat/( densf*conv )

          !print*,'yhat3',yhat

        ! adjust dpobs from wet to dry (shrinking factor)
          dpobs = dpobs * shrinkf


        ! calculate n = vconc / vpt
        ! vpt :: volume dry of particles
        !      = (pi/6.)*(dpobs[m])**3 [m^3]
        ! need number in particles/cm^3
        ! therefore: muliply by 1.e-6
          yhat = yhat/( (pi/6.)*(dpobs)**3.0 )
          yhat = yhat*1.e-6 
       !   print*,'yhat4',yhat



      end subroutine initaero
