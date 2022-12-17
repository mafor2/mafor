! <gde_seaflux.f90 - A component of the Multicomponent
!                     Aerosol Dynamics Model MAFOR>
!*****************************************************************************! 
!* 
!*    Copyright (C) 2011-2022  Matthias Steffen Karl
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
module gde_seaflux


   use gde_input_data, only  : MMAX
   use gde_input_data, only  : NU,AI,AS,CS


implicit none

   INTRINSIC :: SELECTED_REAL_KIND

   public :: seasaltemis

! KPP DP - Double precision kind
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)



contains

  subroutine seasaltemis(IMAX,DPA,u10,sst,sal,owf,saltemis)

    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Calculate seasalt particle flux dF/dr
    !      according to Spada et al. (2013)
    !
    !      interface
    !      ---------
    !
    !        input:
    !           u10:  wind speed                     [m/s]
    !           DPA:  dry particle diameter          [m]
    !           sst:  sea surface temperature        [K]
    !           sal:  salinity in seawater           [g/kg]
    !           owf:  open water fraction            [-]
    !        output:
    !           saltemis: particle number flux       [1/(m^3s)]
    !
    !      method
    !      ------
    !      Combination of the seasalt flux parameterizations of
    !      MO86, SM93 and MA03 to cover the entire diameter spectrum
    !
    !      reference
    !      ---------
    !      Spada, M., Jorba, O., Perez Garcia-Pando, C., Janjic, Z.,
    !      and Baldasano, J. M. (2013)
    !      Modeling and evaluation of the global sea-salt aerosol
    !      distribution: sensitivity to emission schemes and resolution
    !      effects at coastal/orographic sites.
    !      Atmos. Chem. Phys., 13, 11735-11755, doi:10.5194/acp-13-11735-2013
    !
    !      MO86
    !      Monahan, E. C., Spiel, D. E., and Davidson, K. L. (1986)
    !      A Model of Marine Aerosol Generation via Whitecaps and Wave Disruption. 
    !      167–174, Oceanographic Sciences Library, Springer, 
    !      Dordrecht, the Netherlands, doi:10.1007/978-94-009-4668-2_16
    !
    !      SM93
    !      Smith, M. H., Park, P. M., and Consterdine, I. E. (1993)
    !      Marine aerosol concentrations and estimated fluxes over the sea.
    !      Q. J. Roy. Meteor. Soc., 119, 809–824, doi:10.1002/qj.49711951211
    !
    !      MA03
    !      Martensson, E. M., Nilsson, E. D., de Leeuw, G., Cohen, L. H., 
    !      and Hansson, H.-C. (2003) 
    !      Laboratory simulations and parameterization of the primary marine 
    !      aerosol production.
    !      J. Geophys. Res.-Atmos., 108, 4297, doi:10.1029/2002JD002263
    !
    !      Petters, M. D. and Kreidenweis, S. M. (2007)
    !      A single parameter representation of hygroscopic growth and cloud
    !      condensation nuclei activity.
    !      Atmos. Chem. Phys., 7, 1961-1971, 
    !      http://www.atmos-chem-phys.net/7/1961/2007/
    !
    !      Zinke, J., Nilsson, E. D., Zieger, P., and Salter, M. E. (2022)
    !      The effect of seawater salinity and seawater temperature on sea salt
    !      aerosol production. 
    !      J. Geophys. Res.-Atmos., 127, e2021JD036005,
    !      https://doi.org/10.1029/2021JD036005
    ! 
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------


  implicit none

     integer, intent(in)                             :: IMAX
     real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPA     ! [m]
     real( dp), intent(in)                           :: u10       ! [m/s]
     real( dp), intent(in)                           :: sst       ! [K]
     real( dp), intent(in)                           :: sal       ! [g/kg]
     real( dp), intent(in)                           :: owf       ! [-]

     real( dp), dimension(MMAX,IMAX),intent(out)     :: saltemis  ! [1/(m^2*s)]

! local
     real( dp), dimension(MMAX,IMAX)   :: DFDR
     real( dp), dimension(MMAX,IMAX)   :: RPD
     real( dp), dimension(MMAX,IMAX)   :: RP80
     real( dp), dimension(MMAX,IMAX)   :: DLOGDP
     real( dp)                         :: W
     real( dp)                         :: DR
     real( dp)                         :: kappa
     real( dp)                         :: rhwet
     real( dp)                         :: third

! MO86
     real( dp)                         :: B
     real( dp)                         :: mo86
! SM93
     real( dp), dimension(2)           :: r0,A,f
     real( dp), dimension(2)           :: sm93
     real( dp)                         :: sumsm93
! MA03
     real( dp), dimension(3)           :: c0,c1,c2,c3,c4
     real( dp), dimension(3)           :: d0,d1,d2,d3,d4
     real( dp)                         :: AM, BM

     integer                           :: M,I,K,L

! Reference salinity (north sea) [-]
     real( dp), parameter              :: SALREF = 35._dp 

! MA03 polynomial coefficients, three size ranges

        data c0/ -2.881e6 , -6.743e6 ,  2.181e6  /
        data c1/ -3.003e13,  1.183e14, -4.165e12 /
        data c2/ -2.867e21, -8.148e20,  3.132e18 /
        data c3/  5.932e28,  2.404e27, -9.841e23 /
        data c4/ -2.576e35, -2.452e33,  1.085e29 /

        data d0/  7.609e8 ,  2.279e9 , -5.800e8  /
        data d1/  1.829e16, -3.787e16,  1.105e15 /
        data d2/  6.791e23,  2.528e23, -8.297e20 /
        data d3/ -1.616e31, -7.310e29,  2.601e26 /
        data d4/  7.188e37,  7.368e35, -2.859e31 /

        ! Whitecap cover W as fraction [-]
        W = 3.84_dp * 1.e-6 * u10**(3.41_dp)

        ! SM93 parameters
        r0(1) = 2.1   ! um
        r0(2) = 9.2   ! um
        A(1)  = 10._dp**(0.0676*u10+2.43)
        A(2)  = 10._dp**(0.959*u10**(0.5_dp)-1.476)
        f(1)  = 3.1
        f(2)  = 3.3

        ! Hygroscopic Growth factor for NaCl
        ! kappa_mean in Petters and Kreidenfels (2007)
        ! Equation (3) therein:
        ! Vwet = (rh/(1-rh))*kappa*Vdry
        ! Dpwet**3 = Vwet*(6/pi)
        ! Dpwet**3 = (rh/(1-rh))*kappa*Dpdry**3
        ! Dpwet    = ((rh/(1-rh))*kappa)**(1/3) * Dpdry
        kappa = 1.12_dp
        rhwet = 0.80_dp
        third = 1._dp/3._dp


        do M=NU,CS
         do I=1,IMAX

        ! dry particle radius in m 
             RPD(M,I) = 0.5*DPA(M,I)
        ! wet particle radius in m at RH=80%  
             RP80(M,I) = RPD(M,I) * ( (rhwet/(1._dp-rhwet))*kappa )**third


        ! Calculation of particle number flux dF/dr in [m^-3 s^-1]
        ! Smaller particles: dF/dr = 0

             if ( DPA(M,I).lt.2.e-8 ) DFDR(M,I) = 0.0


        ! MA03, dF/dlog(dp) [m^-2 s^-1] for 0.02 < Dp < 2.8 um

             if ( (DPA(M,I).ge.2.e-8).and.(DPA(M,I).lt.1.45e-7) ) then 
               AM = c4(1)*DPA(M,I)**4 + c3(1)*DPA(M,I)**3 + c2(1)*DPA(M,I)**2  +  & 
                    c1(1)*DPA(M,I) + c0(1)
               BM = d4(1)*DPA(M,I)**4 + d3(1)*DPA(M,I)**3 + d2(1)*DPA(M,I)**2  +  &
                    d1(1)*DPA(M,I) + d0(1)
               DFDR(M,I) = W * (AM*sst + BM)
             endif

             if ( (DPA(M,I).ge.1.45e-7).and.(DPA(M,I).lt.4.19e-7) ) then 
               AM = c4(2)*DPA(M,I)**4 + c3(2)*DPA(M,I)**3 + c2(2)*DPA(M,I)**2  +  & 
                    c1(2)*DPA(M,I) + c0(2)
               BM = d4(2)*DPA(M,I)**4 + d3(2)*DPA(M,I)**3 + d2(2)*DPA(M,I)**2  +  &
                    d1(2)*DPA(M,I) + d0(2)
               DFDR(M,I) = W * (AM*sst + BM)
             endif

             if ( (DPA(M,I).ge.4.19e-7).and.(DPA(M,I).le.2.8e-6) ) then 
               AM = c4(3)*DPA(M,I)**4 + c3(3)*DPA(M,I)**3 + c2(3)*DPA(M,I)**2  +  & 
                    c1(3)*DPA(M,I) + c0(3)
               BM = d4(3)*DPA(M,I)**4 + d3(3)*DPA(M,I)**3 + d2(3)*DPA(M,I)**2  +  &
                    d1(3)*DPA(M,I) + d0(3)
               DFDR(M,I) = W * (AM*sst + BM)
             endif


             if ( (DPA(M,I).ge.2.e-8).and.(DPA(M,I).le.2.8e-6) ) then 
        ! Salinity adjustment for MA03
               DFDR(M,I) = DFDR(M,I) * ( sal/SALREF )**third
        ! Particle number flux through sea surface is [m^-2 s^-1]
               saltemis(M,I) = DFDR(M,I)
        ! Zinke et al. (2022): MA03 is roughly one order of magnitude higher
        !                      at high salinity (sal = 35 g/kg)
               saltemis(M,I) = saltemis(M,I) * 0.1_dp
        ! Scaled with the open water fraction
               saltemis(M,I) = saltemis(M,I) * owf
             endif


        ! In Coarse mode (CS): MO86 or SM93
             if (DPA(M,I).gt.2.8e-6) then 

        ! MO86, dF/dr80 [m^-2 s^-1 um^-1] for Dp > 2.8 um

               B = (0.380_dp-log10(RP80(M,I)*1.e6)) / 0.650_dp
               mo86 = 1.373 * (u10)**(3.41_dp) * (RP80(M,I)*1.e6 )**(-3._dp)   *  &
                      (1._dp + 0.057_dp*(RP80(M,I)*1.e6)**1.05)                *  &
                      10._dp**(1.19*exp((-1._dp)*B**2))


        ! SM93, dF/dr80 [m^-2 s^-1 um^-1] for Dp > 2.8 um

               do K=1,2
                 sm93(K) = A(K) * exp( (-1._dp)*f(K)*(log(RP80(M,I)*1.e6/r0(K)))**2 )
               enddo
               sumsm93 = sm93(1) + sm93(2)


               if (u10.lt.9.0_dp) then

                 DFDR(M,I) = mo86

               else

                 DFDR(M,I) = max(mo86,sumsm93)

               endif

        ! Salinity adjustment
               DFDR(M,I) = DFDR(M,I) * ( sal/SALREF )**third

        ! Convert to [m^-3 s^-1]
               DFDR(M,I) = DFDR(M,I)*1.e6

        ! Particle number flux through sea surface is F [m^-2 s^-1]
        ! dF/dlogDp = dF/dDry * ln10 * Ddry
               saltemis(M,I) = DFDR(M,I) * 2.302585093 * DPA(M,I)

        ! Test with scaling by 0.1
               saltemis(M,I) = saltemis(M,I) * 0.1_dp

        ! Scaled with the open water fraction
               saltemis(M,I) = saltemis(M,I) * owf 

             endif

         ! Write the sea-salt particle emissions
     !        write(6,'(I2,A1,I2,4ES12.4)') M,' ',I,DPA(M,I), saltemis(M,I)

         end do
       end do


  end subroutine seasaltemis

!------------------------------------------------------------------


end module gde_seaflux
