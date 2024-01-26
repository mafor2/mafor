! <awater.for - A component of FITAERO
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

      subroutine awater(aw,Msu,Mamsu,Mamnit,Msalt,Morgc,wh2o)


!***********************************************************************
!***  Subroutine awater calculates water content per bin
!***********************************************************************
!
!      Water activity below DRH
!
!      purpose
!      -------
!      calculate water activity with polynomals from Tang and
!      Munkelwitz, JGR 99: 18801-18808, 1994
!      22.08.2012: added NaCl (=SALT)
!
!      method
!      ------
!
!      input is aerosol and gas phase conc of sulfate, ammonium, nitrate
!      in molec/m3 or ng/m3
!      output is water concentration on aerosol in ng/m3
!      definitions:
!      mso4, mnh4, and mno3 are the number of micromoles/(cubic meter of air)
!      for sulfate, ammonium, and nitrate respectively
!      irhx is the relative humidity (%)
!      wh2o is the returned water amount in nanograms / cubic meter of air
!      x is the molar ratio of ammonium to sulfate
!      y0,y1,y1.5, y2 are the water contents in mass of water/mass of solute
!      for pure aqueous solutions with x equal 1, 1.5, and 2 respectively.
!      y3 is the value of the mass ratio of water to solute for
!      a pure ammonium nitrate  solution.
!
!      external
!      --------
!      none
!
!      reference
!      ---------
!      none
!
!------------------------------------------------------------------

!     Declarations of variables by using the MODULES feature:

      use module_fitaero_exe


      implicit none

!***********************************************************************


      REAL, intent(in)                :: Msu,Mamsu,Mamnit,Msalt
      REAL, intent(in)                :: Morgc,aw
      REAL, intent(out)               :: wh2o

!     Local declarations:

      REAL                            :: mso4,mnh4,mno3,msal,awc
      REAL                            :: morg,torg
      REAL                            :: tso4,tnh4,tno3,tsal,x,null
      REAL                            :: mfs0,mfs1,mfs15,mfs2
      REAL                            :: c0(4),c1(4),c15(4),c2(4)
      REAL                            :: y, y0,y1,y15,y2,y3,y40,y5,y6
      REAL                            :: kSO4(6),kNO3(6),mfsSO4
      REAL                            :: kNaCl(6),kSucc(6)
      REAL                            :: mfsNO3,mfsNaCl,mfsSucc
      REAL                            :: u,y140,y1540,yc

! *** molecular weights:
      REAL, parameter ::                             &
                mwh     = 1.0,                       &
                mwso4   = 96.0636,                   &
                mwnh4   = 18.0985,                   &
                mwno3   = 62.0649,                   &
                mwcl    = 35.45,                     &
                mwna    = 22.99,                     &
                mwsuc   = 118.0,                     &
                mw2     = mwso4 + 2.0 * mwnh4,       &
                mwano3  = mwno3 + mwnh4,             &
                mwnacl  = mwcl + mwna


!***********************************************************************
!     Content of subroutine:

!     The polynomials use data for aw as a function of mfs from Tang and
!     Munkelwitz, JGR 99: 18801-18808, 1994.
!     The polynomials were fit to Tang's values of water activity as a
!     function of mfs.

! *** coefficients of polynomials fit to Tang and Munkelwitz data
!     now give mfs as a function of water activity.

      data c1/0.9995178, -0.7952896, 0.99683673, -1.143874/
      data c15/1.697092,-4.045936, 5.833688, -3.463783/
      data c2/2.085067, -6.024139, 8.967967, -5.002934/

! *** the following coefficients are a fit to the data in Table 1 of
!     Nair & Vohra, J. Aerosol Sci., 6: 265-271, 1975
!      data c0/0.8258941, -1.899205, 3.296905, -2.214749 /
! *** New data fit to data from
!       Nair and Vohra J. Aerosol Sci., 6: 265-271, 1975
!       Giaque et al. J.Am. Chem. Soc., 82: 62-70, 1960
!       Zeleznik J. Phys. Chem. Ref. Data, 20: 157-1200
      data c0/ 0.798079, -1.574367, 2.536686, -1.735297 /

! *** polynomials for ammonium nitrate and ammonium sulfate are from:
!     Chan et al.1992, Atmospheric Environment (26A): 1661-1673.

      data kNO3/0.2906, 6.83665, -26.9093,   &
               46.6983, -38.803, 11.8837/
      data kSO4/ 2.27515, -11.147, 36.3369,  &
             -64.2134, 56.8341, -20.0953/

! *** polynomials for sodium chloride are from:
!     Tang, I.N., Tridico, A.C., Fung, K.H.
!     Thermodynamic and optical properties of sea salt aerosols
!     J. Geophys. Res., 102(D19), 23,269-23,275, 1997
!     polynomial fit according to Chan et al. 1992 method

      data kNaCl/1.2617, -4.3239, 10.658,    &
                -16.248,  13.324, -4.6718/

! *** polynomials for sodium succinate are from:
!     Peng, C. and C. K. Chan
!     The water cycles of water-soluble organic salts of atmospheric importance
!     Atmos. Environ., 35, 1183-1192, 2001

      data kSucc/0.38225, 5.4575, -21.965,   &
                34.34, -23.764, 5.5627/


! *** check range of per cent relative humidity
!       aw water activity = fractional relative humidity
! calculate masses of so4, nh4 and no3 from the masses
! of h2so4,(nh4)2so4 and nh4no3 in ng/m^3
       mso4 = Msu/(2.*mwh+mwso4)*mwso4 + Mamsu/mw2*mwso4
       mnh4 = Mamsu/mw2*mwnh4 + Mamnit/mwano3*mwnh4
       mno3 = Mamnit/mwano3*mwno3
       msal = Msalt/mwnacl*mwcl
       morg = Morgc/mwsuc

      null=0.0
      mfs0=0.0
      mfs1=0.0
      mfs15=0.0
      tso4 = max( mso4 , null )
      tnh4 = max( mnh4 , null )
      tno3 = max( mno3 , null )
      tsal = max( msal , null )
      torg = max( morg , null )
      x = 0.0
! *** if there is non-zero sulfate calculate the molar ratio
      if (tso4 .gt. 0.0 ) then
        x = tnh4 / tso4
      else
! *** otherwise check for non-zero nitrate and ammonium
        if ( (tno3 .gt. 0.0) .and. (tnh4 .gt. 0.0) ) x = 10.0
      end if
      y  = 0.0
      y2 = 0.0
      y3 = 0.0
      y5 = 0.0
      y6 = 0.0

! *** begin screen on x for calculating wh2o
      if ( x .lt. 1.0 ) then
!
          mfs0 = poly4(c0(1),c0(2),c0(3),c0(4),aw)
          mfs1 = poly4(c1(1),c1(2),c1(3),c1(4),aw)
          y0 = (1.0 - mfs0 ) / mfs0
          y1 = (1.0 - mfs1 ) / mfs1
          y = (1.0- x) * y0 + x * y1
!
       else if ( x .lt. 1.5) then
!
         if ( aw .ge. 0.40 ) then
            mfs1  = poly4(c1(1),c1(2),c1(3),c1(4),aw)
            mfs15 = poly4(c15(1),c15(2),c15(3),c15(4),aw)
            y1  = (1.0 - mfs1 ) / mfs1
            y15 = (1.0 - mfs15) / mfs15
            y = 2.0 * ( y1 * (1.5 - x) + y15 *( x - 1.0) )
         else
! *** set up for crystalization
! *** Crystallization is done as follows:
!      For 1.5 <= x, crystallization is assumed to occur at rh = 0.4
!      For x <= 1.0, crystallization is assumed to occur at an rh < 0.01,
!      and since the code does not allow ar rh < 0.01, crystallization
!      is assumed not to occur in this range.
!      For 1.0 <= x <= 1.5 the crystallization curve is a straignt line
!      from a value of y15 at rh = 0.4 to a value of zero at y1. From
!      point B to point A in the diagram.
!      The algorithm does a double interpolation to calculate the amount of
!      water.
!
!        y1(0.40)               y15(0.40)
!         +                     + Point B
!
!         +--------------------+
!       x=1                   x=1.5
!      Point A
!
           awc = 0.80 * (x - 1.0) ! rh along the crystallization curve.
           u = 0.4
           y = 0.0
           if ( aw .ge. awc ) then ! interpolate using crystalization curve
               mfs1  = poly4(c1(1),c1(2),c1(3),c1(4),u)
               mfs15 = poly4(c15(1),c15(2),c15(3),c15(4),u)
               y140  = (1.0 - mfs1 ) / mfs1
               y1540 = (1.0 - mfs15) / mfs15
               y40 = 2.0 * ( y140 * (1.5 - x) + y1540 *( x - 1.0) )
               yc = 2.0 * y1540 * (x -1.0) ! y along crystallization curve
               y = y40 - (y40 - yc) * (0.40-aw) / (0.40 - awc)
            end if ! end of checking for aw
          end if ! end of checking on irh

       else if( x .lt. 1.9999) then
!
           y= 0.0
           if( aw .ge. 0.40) then
             mfs15 = poly4(c15(1),c15(2),c15(3),c15(4),aw)
             mfs2  = poly4(c2(1),c2(2),c2(3),c2(4),aw)
             y15 = (1.0 - mfs15) / mfs15
             y2  = (1.0 - mfs2) / mfs2
             y = 2.0 * (y15 * (2.0 - x) + y2 * (x - 1.5) )
           end if ! end of check for crystallization
!
      else ! 1.9999 < x

! regime where ammonium sulfate and ammonium nitrate are in solution.
!
! *** following cf&s for both ammonium sulfate and ammonium nitrate
! *** check for crystallization here. their data indicate a 40% value
!     is appropriate.
            y2 = 0.0
            y3 = 0.0
            if ( aw .ge. 0.40) then
              mfsSO4 = poly6(kSO4(1),kSO4(2),kSO4(3),kSO4(4),kSO4(5),kSO4(6),aw)
              mfsNO3 = poly6(kNO3(1),kNO3(2),kNO3(3),kNO3(4),kNO3(5),kNO3(6),aw)
              y2 = (1.0 - mfsSO4) / mfsSO4
              y3 = (1.0 - mfsNO3) / mfsNO3

            end if
!
      end if ! end of checking on x

! for sea-salt, represented by NaCl
! CRH(NaCl): 0.47 (crystallisation below RH=47%)
      if ( aw .ge. 0.47) then
         mfsNaCl = poly6(kNaCl(1),kNaCl(2),kNaCl(3),kNaCl(4),kNaCl(5),kNaCl(6),aw)
         y5 = (1.0 - mfsNaCl) / mfsNaCl
      end if

! for organics, represented by Na-succinate
! CRH(NaSucc): 0.48 (crystallisation below RH=48%)
      if ( aw .ge. 0.48) then
         mfsSucc = poly6(kSucc(1),kSucc(2),kSucc(3),kSucc(4),kSucc(5),kSucc(6),aw)
         y6 = (1.0 - mfsSucc) / mfsSucc
      end if

! *** now set up output of wh2o

!      wh2o units are nanograms (liquid water) / cubic meter of air
!
      if ( x .lt. 1.9999) then
         wh2o = y*(tso4 + tnh4)
         wh2o = wh2o + y5*tsal
         wh2o = wh2o + y6*torg
      else

! *** this is the case that all the sulfate is ammonium sulfate
!     and the excess ammonium forms ammonum nitrate

        wh2o = y2*tso4 + y3 * tno3
        wh2o = wh2o + y5*tsal
        wh2o = wh2o + y6*torg
      end if
      !write(6,*) y,y2,y3,y5



      end subroutine awater
