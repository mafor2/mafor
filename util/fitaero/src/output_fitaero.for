! <output_fitaero.for -  A component of FITAERO
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

      subroutine output_fitaero(para,ndata,ndens,num,wamass,shrf,  &
                                binlo,binup,m,mstart,nout)

!***********************************************************************
!***  Subroutine output_fitdis writes output of the estimated size
!***  distribution and of inaero.dat for MAFOR
!***********************************************************************

!     Declarations of variables by using the MODULES feature:

      use module_fitaero_io
      use module_fitaero_exe


      implicit none

!***********************************************************************

      integer,dimension(nm), intent(in)       :: binlo,binup
      real, dimension(n_obins,nv), intent(in) :: ndata
      real, dimension(nm,np), intent(in)      :: para
      real, dimension(nm,nc-1), intent(in)    :: ndens
      real, dimension(nm), intent(in)         :: num
      real, dimension(nm), intent(in)         :: wamass
      real, dimension(nm), intent(in)         :: shrf
      real, dimension(nm,n_obins), intent(in) :: nout
      integer, intent(in)                     :: m
      integer, intent(in)                     :: mstart

!     Local declarations:

      integer                           :: i,n,j

      real, dimension(nm)               :: gmdout
      real, dimension(nm)               :: sigout
      real, dimension(nm)               :: msulf, morgc, mammo
      real, dimension(nm)               :: mnitr, mmsap, msalt
      real, dimension(nm)               :: mxxxx, mecbc, mdust

      character(len=1)                  :: tab = ACHAR(9)  ! TAB delimiter

!***********************************************************************
!     Content of subroutine:



       if (m.eq.mstart) then

! ***   Open output files:

         fname_out_sizedis  = trim(fname_outpath)//'/out_sizedis.dat'
         funit_out_sizedis  = get_file_unit()
         open (funit_out_sizedis, file = fname_out_sizedis,   &
              status = 'unknown', form  = 'formatted', action = 'write')


         fname_out_inaero   = trim(fname_outpath)//'/inaero.dat'
         funit_out_inaero   = get_file_unit()
         open (funit_out_inaero, file = fname_out_inaero,     &
              status = 'unknown', form  = 'formatted', action = 'write')

! ***   Write a header line
         write (funit_out_sizedis,'(A41)')    '*   Dp[nm]   dN_obs[cm-3]     dN_est[cm-3]'
       endif

! ***   Write data to file

       do i= binlo(m),binup(m)

         write (funit_out_sizedis, 1000) ndata(i,1),ndata(i,2),nout(m,i)

       enddo

! ***  Calculate mid bin of AI mode
!       j = INT(binup(2)-binlo(2))-1


       do n= 1,nm
         print *,'Dp in n',n,para(n,1)
       enddo

! ***   Write inaero.dat
! ***   ndis GMD SIGMA NTOT H2SO4 ORG NH4 NO3 MSAp SALT BPAB BC DUST
       if (m.eq.nm) then

         write (funit_out_inaero, 2000) dpmax,imax


         do n= 1,nm

         print *,'water',para(n,3),wamass(n),para(n,3)*(1.0-wamass(n))


         ! Not sure yet if the water shall be subtracted

           msulf(n) = ndens(n,SU)*para(n,3) !*(1.0-wamass(n))
           morgc(n) = ndens(n,OC)*para(n,3) !*(1.0-wamass(n))
           mammo(n) = ndens(n,AM)*para(n,3) !*(1.0-wamass(n))
           mnitr(n) = ndens(n,NI)*para(n,3) !*(1.0-wamass(n))
           mmsap(n) = ndens(n,MS)*para(n,3) !*(1.0-wamass(n))
           msalt(n) = ndens(n,SA)*para(n,3) !*(1.0-wamass(n))
           mxxxx(n) = ndens(n,XX)*para(n,3) !*(1.0-wamass(n))
           mecbc(n) = ndens(n,EC)*para(n,3) !*(1.0-wamass(n))
           mdust(n) = ndens(n,DU)*para(n,3) !*(1.0-wamass(n))


           sigout(n) = para(n,2)
           gmdout(n) = para(n,1)*1.e-9  !GMD in [m]
           
! ***   Function to reduce the wet diameter giving higher dN/dlogDp 
! ***   in MAFOR. This applies only to mode 2 and 3

           if (n==1) then
             if (rhi>0.89) then
               gmdout(n) = gmdout(n)*0.85
             endif
           endif
! AI MODE
           if (n==2) then
             if ( (ndens(2,SU)>0.15).and.(rhi>0.15) ) then
             !use DU in AI mode to divide between exhaust and background
           !exhaust AI: DU<0.1, EC>0.15
               if ((ndens(2,DU)<0.1).and.(ndens(2,EC)>0.15)) then
                 if (rhi>0.89) then
                   gmdout(n) = gmdout(n)*0.75
                 else if (rhi>0.15) then
                   gmdout(n) = gmdout(n)*0.84
                 endif
               else   
               
           !background
                 if (dpmax==1.0E-5) then
                   gmdout(n) = gmdout(n)*1.00             
                 else if (dpmax<2.0E-6) then
                   gmdout(n) = gmdout(n)*0.92               
                 else
                   if (rhi>0.85) then
                      gmdout(n) = gmdout(n)*1.1
                      if (ndens(3,SU)<0.3) then
                        sigout(n) = sigout(n)*0.92   ! coastal
                      else
                        sigout(n) = sigout(n)*0.91   ! marine
                        gmdout(n) = min(gmdout(n), 0.50E-7)
                      endif
                   else if (rhi>0.45) then
                     if (para(2,1).gt.95.0) then
                       gmdout(n) = gmdout(n)*1.05
                     else
                       gmdout(n) = gmdout(n)*1.4
                     endif
                   else if (rhi>0.15) then
                     gmdout(n) = gmdout(n)*0.89
                   endif
                 endif

               endif
             else
               gmdout(n) = gmdout(n)*max(shrf(n),0.96)
             endif
! AS MODE
           elseif (n==3) then

             if (dpmax==1.0E-5) then
                gmdout(n) = gmdout(n)*1.00
             else if (dpmax<2.0E-6) then
                gmdout(n) = gmdout(n)*max(shrf(n),0.92)
             else
       !MSK 26.11.2020 MARINE & COASTAL
             ! SALT > 0.1 [marine case]
               if (ndens(3,SA)>0.1) then
                   gmdout(n) = max(gmdout(n)*0.66, 1.80E-7)
                   sigout(n) = sigout(n)*0.85
             ! SULF > 0.2 [coastal case] 
               else if (ndens(3,SU)>0.2) then
                   gmdout(n) = gmdout(n)*0.70
                   sigout(n) = sigout(n)*1.08
               else
               !!! do not change !!!
                   gmdout(n) = gmdout(n)*max(shrf(n),0.92)
               endif
             endif
! CS MODE
           elseif (n==4) then
             gmdout(n) = gmdout(n)*max(shrf(n),0.90)
           endif

           print *,'out',n,gmdout(n)*1.e9


           write (funit_out_inaero, 2100) mcorr,tab,                   &
                gmdout(n),tab,sigout(n),tab,                          &
                num(n),tab,msulf(n),tab,morgc(n),tab,                 &
                mammo(n),tab,mnitr(n),tab,mmsap(n),tab,               &
                msalt(n),tab,mxxxx(n),tab,mecbc(n),tab,mdust(n)

         enddo

       endif

      return


 1000 format(F11.2, 4X, F11.2, 4X, F11.2)
 2000 format(E9.2, 2X, I3)
 2100 format(1L1,A,E8.2,A,F0.3,A,E8.2,A, 9(F0.3,A))

      end subroutine output_fitaero
