! <gde_seaflux.f90 - A component of the Multicomponent
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
   use gde_input_data, only  : NU,NA,AI,AS,CS
   use gde_constants,  only  : pi


implicit none

   INTRINSIC :: SELECTED_REAL_KIND

   public :: seasaltemis

! KPP DP - Double precision kind
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)



contains

  subroutine seasaltemis(IMAX,DPA,DLOGDP,u10,sst,sal,owf,saltemis)

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
    !      Calculate seasalt particle flux dF/dlogDp
    !      according to Salter et al. (2015).
    !      Replaces earlier routine with same name
    !      that calculated seasalt particle flux
    !      according to Spada et al. (2013).
    !
    !      interface
    !      ---------
    !
    !        input:
    !           u10:  wind speed                     [m/s]
    !           DPA:  dry particle diameter          [m]
    !           DLOGDP:                              [-]
    !           sst:  sea surface temperature        [K]
    !           sal:  salinity in seawater           [g/kg]
    !           owf:  open water fraction            [-]
    !        output:
    !           saltemis: particle number flux       [1/(m^3s)]
    !
    !      method
    !      ------
    !      Empirical seasalt flux parameterization that covers
    !      the entire diameter spectrum.
    !      For sectional representation dF/dlogDp is multiplied by dlogDp.
    !
    !      reference
    !      ---------
    !      SAL15
    !      Salter, M. D., Zieger, P., Acosta Navarro, J. C., Grythe, H.,
    !      Kirkevag, A., Riipinen, I., Nilsson, E. D. (2015)
    !      An empirically derived inorganic sea spray source function
    !      incorporating sea surface temperature.
    !      Atmos. Chem. Phys., 15, 11047-11066,
    !      https://doi.org/10.5194/acp-15-11047-2015
    !
    !      Spada, M., Jorba, O., Perez Garcia-Pando, C., Janjic, Z.,
    !      and Baldasano, J. M. (2013)
    !      Modeling and evaluation of the global sea-salt aerosol
    !      distribution: sensitivity to emission schemes and resolution
    !      effects at coastal/orographic sites.
    !      Atmos. Chem. Phys., 13, 11735-11755, doi:10.5194/acp-13-11735-2013
    !
    !      Zinke, J., Nilsson, E. D., Zieger, P., and Salter, M. E. (2022)
    !      The effect of seawater salinity and seawater temperature on 
    !      sea salt aerosol production. 
    !      J. Geophys. Res.-Atmos., 127, e2021JD036005,
    !      https://doi.org/10.1029/2021JD036005
    !
    !      Ovadnevaite, J., Ceburnis, D., Canagaratna, M., Berresheim, H.,
    !      Bialek, J., Martucci, G., Worsnop, D. R. and O’Dowd, O. (2012)
    !      On the effect of wind speed on submicron sea salt mass concentrations
    !      and source fluxes, J. Geophys. Res., 117, D16201,
    !      doi:10.1029/2011JD017379.
    !
    !      modifications
    !      -------------
    !      Salinity adjustment by multiplication with 
    !      (sal/SALREF)**third
    !      Multiplication of dlogF/dlogDp by dlogDp
    !      If OWF < 70% use Ovadnevaite et al. (2012) exponent 2.7
    !      for the wind dependence
    !
    !------------------------------------------------------------------


  implicit none

     integer, intent(in)                              :: IMAX
     real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA       ! [m]
     real( dp), dimension(MMAX,IMAX), intent(in)      :: DLOGDP    ! [-]
     real( dp), intent(in)                            :: u10       ! [m/s]
     real( dp), intent(in)                            :: sst       ! [K]
     real( dp), intent(in)                            :: sal       ! [g/kg]
     real( dp), intent(in)                            :: owf       ! [-]

     real( dp), dimension(MMAX,IMAX),intent(out)      :: saltemis  ! [1/(m^2*s)]

! local
     real( dp), dimension(MMAX,IMAX)   :: DFDR

     real( dp)                         :: W
     real( dp)                         :: third
     real( dp)                         :: sstc
     real( dp)                         :: logdp

!SAL15
! Number of modes in SAL15
     integer, parameter                :: nmod = 3

     real( dp), dimension(nmod)        :: effprod
     real( dp), dimension(nmod)        :: NP
     real( dp), dimension(nmod)        :: dpmod,sigmo
     real( dp), dimension(nmod)        :: a0,b0,c0,d0


     integer                           :: M,I
     integer                           :: j

! Reference salinity (north sea) [-]
     real( dp), parameter              :: SALREF = 35._dp 


! SAL15 Modal parameters for the SST dependence given
! for three aerosol modes. Table 2 in Salter et al. (2015)
! Modes
        data dpmod/ 0.095e-6, 0.60e-6, 1.50e-6   /
        data sigmo/ 2.10    , 1.72   , 1.60      /
! Ai to Di
        data a0/ -5.2168e5 ,  0.00    , 0.00      /
        data b0/ 3.31725e7 ,  7.374e5 , 1.4210e4  /
        data c0/ -6.95275e8, -2.4803e7, 1.4662e7  /
        data d0/ 1.0684e10 , 7.7373e8 , 1.7075e8  /


!SAL15
! Entrainment flux of air into the ocean water column
! according to Eq.(7) in Salter et al. (2015)
! W :: Entrainment flux [m^3 m^-2 s^-1]


        if (owf.lt.0.70) then

!Ovadnevaite JGR 2012 !EAC2024 abstract
           if (u10>3.7_dp) then
             W = 0.1_dp * 1.e-8 * u10**(2.7_dp)
           else
             W = 0._dp
           endif

        else

           W = 2.0_dp * 1.e-8 * u10**(3.41_dp)

        endif 

        third = 1._dp/3._dp

! Sea surface temperature in degree Celsius
        sstc = sst - 273.0_dp

! First calculate number production flux
! in each of the three modes
! according to Eq.(9) in Salter et al. (2015)
! NP :: Number production flux [m^3 m^-2 s^-1]

       do j=1,nmod

           NP(j) = W * ( a0(j)*sstc**3 + b0(j)*sstc**2 + c0(j)*sstc + d0(j) )

       end do

! No particles in NU mode
       do I=1,IMAX
           DFDR(NU,I)     = 0.0_dp
           saltemis(NU,I) = 0.0_dp
       end do

       do M=NA,CS
         do I=1,IMAX

        ! initialize DFDR
            DFDR(M,I) = 0.0_dp
            logdp = LOG(DPA(M,I))

        ! Particle number flux through sea surface is F [m^-2 s^-1]
        ! dF/dlogDp = dF/dDry * ln10 * Ddry
        ! Calculate dF/dlog(dp) [m^-2 s^-1] for all size bins
        ! according to Eq.(8) in Salter et al. (2015)
            do j=1,nmod
              effprod(j) =  NP(j) / (SQRT(2._dp*pi)*LOG(sigmo(j)))
              effprod(j) = effprod(j)                                    *  &
                           EXP(-0.5_dp*( ((logdp-LOG(dpmod(j)))**2._dp)  /  &
                           ((LOG(sigmo(j)))**2._dp) ))

              DFDR(M,I) = DFDR(M,I) + effprod(j)
            end do


        ! No sea salt flux for particles with Ddry < 10 nm
            if (DPA(M,I).lt.1.0e-8) DFDR(M,I) = 0._dp

        ! Particle number flux through sea surface is [m^-2 s^-1]
        ! DFDR is dF/dlogDp
        ! dF = dF/dlogDp * dlogDp
            saltemis(M,I) = DFDR(M,I) * DLOGDP(M,I)

        ! Salinity adjustment
        ! The sea salt number, surface area and mass (or volume) 
        ! emissions are all scaled by the salinity SAL.
        ! NOTE: However, when assuming that the water droplet generation 
        ! process is not affected by variable salinity, then the 
        ! dry particle diameter and mass emission should be scaled by 
        ! the salinity correction, but not number emission.
        ! Zink et al. (2022), Figure 1 shows for SAL>10 g/kg that
        ! number, surface area and volume distribution increases
        ! with increasing salinity. A shift of diameter is only seen
        ! for SAL<10 g/kg.
            saltemis(M,I) = saltemis(M,I) * ( sal/SALREF )**third

        ! Linear scaling with the open water fraction
            saltemis(M,I) = saltemis(M,I) * owf


         end do
       end do

! Write the sea-salt particle emission dF/dlogDp
! for comparison with literature
! DFDR(M,I) := dF/dlogDp
!       do M=NU,CS
!         do I=1,IMAX
!             write(6,'(I2,A1,I2,4ES12.4)') M,' ',I,DPA(M,I),DFDR(M,I),saltemis(M,I)
!         end do
!       end do
!    stop


  end subroutine seasaltemis

!------------------------------------------------------------------


end module gde_seaflux
