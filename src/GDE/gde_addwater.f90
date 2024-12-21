! <gde_addwater.f90 - A component of the Multicomponent
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
!*    except:
!*    routine AWATER written by Francis S. Binkowski, and belonging
!*    functions POLY4 and POLY6 written by Francis S. Binkowski
!* 
!*****************************************************************************!
module gde_addwater

    use messy_mecca_kpp_Global, only    : APN, dp

    use gde_constants,  only            : M_H2O,N_A,RHOH2O,pi

    use gde_toolbox,    only            : molec2ug
    use gde_sensitiv,   only            : IDEB

    use gde_input_data, only            : AMAX
    use gde_input_data, only            : MMAX
    use gde_input_data, only            : NU,NA,AI,AS,CS
    use gde_input_data, only            : A_SUL,A_NH4,A_AMI,A_NIT,A_XXX
    use gde_input_data, only            : A_OR1,A_OR2,A_SAL,A_CHL
    use gde_input_data, only            : A_WAT
    use gde_input_data, only            : CONVM
    use gde_input_data, only            : rhactiv
    use gde_input_data, only            : rho_air
    use gde_input_data, only            : surf_h2o_std

    use gde_init_aero,  only            : GMD

    private

    public :: water_content
    public :: wetdiameter
    private :: partlwc
    private :: awater
 
  contains


    subroutine water_content(incloud,firstloop,IMAX,RH,lwcm,MASS,N,dmwdt, &
               DPA,SIG,DLINDP,VPT,DPAW,MH2OTOT,NTOTCS,                &
               lwc,cwa01,cwa02,cwa03,xaer,wascloud1,wascloud2  )
    !----------------------------------------------------------------------
    !
    !****  central routine to compute aerosol water content and
    !      liquid water content in and out of clouds   
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate c(ind_H2O_a**)
    !      calculate lwc [m3/m3]
    !      calculate MH2OTOT
    !
    !      interface
    !      ---------
    !
    !        input:
    !          RH,lwcm
    !          dmwdt  (mass conc of water per bin)
    !
    !        output:
    !          c(ind_H2O_a**)
    !          lwc
    !          MH2OTOT
    !          
    !      method
    !      ------
    !      Water activity from sulphate - ammonium - nitrate. MH2O in ng/m^3
    !      Munkelwitz and Tang 1994 / Seinfeld and Pandis 1997
    !      Sea salt hygroscopicity same  
    !      If in-cloud, MH2O and N of CS mode is calculated from LWC(CS)
    !      Three cases are divided:
    !      1) wet aerosol, no fog/cloud
    !      2) cloud droplets, fog/cloud, prescribed LWC
    !      3) cloud droplets, RH>RHactive, fog/cloud, dynamic LWC
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

    implicit none

    ! input
    INTEGER, intent(in)                               :: IMAX
    real( dp), intent(in)                             :: RH,lwcm
    real( dp), intent(in)                             :: SIG(NU:CS)
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPA
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPAW
    real( dp), dimension(MMAX,IMAX),intent(in)        :: DLINDP
    real( dp), dimension(MMAX,IMAX),intent(in)        :: VPT
    real( dp), DIMENSION(MMAX,IMAX), intent(in)       :: dmwdt

    integer, intent(in)                               :: incloud
    logical, intent(in)                               :: firstloop

    ! in/out
    real( dp), dimension(MMAX,IMAX,AMAX),intent(in out) :: MASS
    real( dp), DIMENSION(MMAX,IMAX), intent(in out)     :: N
    real( dp),dimension(APN), intent(in out)            :: lwc
    real( dp),dimension(APN), intent(in out)            :: xaer
    logical, intent(in out)                             :: wascloud1
    logical, intent(in out)                             :: wascloud2

    ! output
    real( dp), intent(out)                            :: MH2OTOT(NU:CS)
    real( dp), intent(out)                            :: cwa01
    real( dp), intent(out)                            :: cwa02
    real( dp), intent(out)                            :: cwa03
    real( dp), intent(out)                            :: NTOTCS

    ! local
    real( dp), DIMENSION(MMAX,IMAX)                   :: MH2OA
    real( dp)                   :: lindismwa
    real( dp)                   :: lindisn
    !!!real( dp)                   :: cwa01,cwa02
    integer                     :: M,I
    integer                     :: cmode


! Initialize all terms
! that are output and that are local or locally allocated
      cwa01      = 0._dp
      cwa02      = 0._dp
      cwa03      = 0._dp


      if (incloud.EQ.1) then
      
! For fixed cloud water droplets in coarse mode ("SINTEF fog simulation")

          NTOTCS     = 0._dp
          cmode      = 3     ! 1 = CS mode

! aerosol water content NU ... AS
          do M=NU,AS
            MH2OTOT(M) = 0._dp 
            do I=1,IMAX
              CALL awater(RH,MASS(M,I,A_SUL),MASS(M,I,A_NH4),MASS(M,I,A_NIT),            &
                       MASS(M,I,A_SAL)+MASS(M,I,A_CHL),MASS(M,I,A_OR1)+MASS(M,I,A_OR2),  &
                       MASS(M,I,A_WAT))
              if (RH.LT.0.05) MASS(M,I,A_WAT)=0._dp           
              MH2OTOT(M)=MH2OTOT(M)+MASS(M,I,A_WAT)    
            end do
          end do

! droplet water for aqueous phase chemistry (AI,AS,CS)
! get xaer (aq. phase partitioning)
          CALL partlwc(MH2OTOT,cwa01,cwa02,cwa03,xaer,lwc)


! intiialize lwc
          lwc(:) = 1.e-32_dp

! use LWC for CS mode from ingeod.dat
! recalculate water content and lwc in coarse mode
! droplet distribution  CS

! existing aerosol mass in CS mode is lost
!!!          MASS(CS,1:IMAX,1:AMAX)=0._dp

! prescribe droplet distribution in CS mode
          lwc( cmode )=lwcm

          cwa03=RHOH2O*lwc( cmode )*N_A*1.E-3_dp/M_H2O

          MH2OTOT(CS)=cwa03*1.e3*molec2ug(M_H2O)

          NTOTCS=0._dp

! number of droplets
! physical correct estimate of NTOTCS:
!   NTOTCS=lwc( cmode )/( (RHOH2O/rho_air)*(pi/6.)*GMD(CS)**3.)
! we assume that mean diameter of fog droplets is GMD(CS)
! and set rho_air to 1.2928 kg/m3

          NTOTCS = lwc( cmode )/ ( (RHOH2O/1.2928_dp)*(pi/6.)*GMD(CS)**3.)


! distribute over bins of CS mode
          do I=1,IMAX
            lindismwa=(MH2OTOT(CS)/(SQRT(2._dp*pi)*DPA(CS,I)*LOG(SIG(CS)))) &
                  *EXP(-0.5_dp*(LOG(DPA(CS,I)/GMD(CS))/LOG(SIG(CS)))**2._dp)
            lindisn  =(NTOTCS/(SQRT(2._dp*pi)*DPA(CS,I)*LOG(SIG(CS)))) &
                  *EXP(-0.5_dp*(LOG(DPA(CS,I)/GMD(CS))/LOG(SIG(CS)))**2._dp)
            MASS(CS,I,A_WAT)= lindismwa*DLINDP(CS,I)
            N(CS,I)         = lindisn  *DLINDP(CS,I)

          end do

          !print *,'cloud=1, lwc, N(CS)',lwc( cmode ), NTOTCS

          wascloud1=.true.

          if (firstloop) then
            write(6,fmt='(a,f13.5)') 'init drop number in mode CS',NTOTCS*1.e-6
            write(6,fmt='(a,3es12.4)') 'incloud lwc',lwc( cmode ),cwa03,MH2OTOT(CS)
            if (IDEB == 1) then
              write(12,fmt='(a,f13.5)' ) 'Init drop number in mode CS',NTOTCS*1.e-6      
              write(12,fmt='(a,e10.4,a,e10.4,a,e10.4)' ) 'incloud lwc ',lwc( cmode ), &
                   '  c(H2O)aq ',cwa03,'  m(H2O) ',MH2OTOT(CS) 
            endif
          endif


      else if ((incloud.EQ.2).and.(RH.ge.rhactiv)) then

! For cloud activation when RH > RH_activation
!
! Cloud droplet activation happens if relative humidity is above 99%
! and the incloud flag is set to 1 by the user and if the computed
! supersaturation is greater than equlibrium supersaturation.
! Kinetic growth by water condensation in CS mode      
! Calculate next MASS(CS,I,A_WAT) by kinetic growth (not equilibrium)

! The water mass increase is calculated in gde_cloud_activation
! NEW WATER MASS [ng/m^3]
! Seinfeld and Pandis (2006), EQ(17.65)

          ! aerosol water content AI ... CS
          do M=AI,CS
            do I=1,IMAX
              !print *,'MH2OA',M,I,MASS(M,I,A_WAT),dmwdt(M,I)
              MH2OA(M,I)  = MASS(M,I,A_WAT) + dmwdt(M,I)
            end do
          end do


          ! first calculate water mass using the ZRS ruleset
          ! aerosol water content NU ... CS
          do M=NU,CS
            MH2OTOT(M) = 0._dp 
            do I=1,IMAX
              CALL awater(RH,MASS(M,I,A_SUL),MASS(M,I,A_NH4),MASS(M,I,A_NIT),        &
                       MASS(M,I,A_SAL)+MASS(M,I,A_CHL),MASS(M,I,A_OR1)+MASS(M,I,A_OR2),   &
                       MASS(M,I,A_WAT))
              MH2OTOT(M)=MH2OTOT(M)+MASS(M,I,A_WAT)
            end do
          end do


          ! use the greater water mass of the two methods
          ! aerosol water content AI ... CS
          do M=AI,CS
            MH2OTOT(M) = 0._dp 
            do I=1,IMAX

              !print *,'massH2O',M,I,MASS(M,I,A_WAT),MH2OA(M,I)
              MASS(M,I,A_WAT) = max(MASS(M,I,A_WAT),MH2OA(M,I))
              MH2OTOT(M) = MH2OTOT(M) + MASS(M,I,A_WAT)

            enddo
          enddo


! droplet water for aqueous phase chemistry
          CALL partlwc(MH2OTOT,cwa01,cwa02,cwa03,xaer,lwc)

          ! No LWC in Aitken mode if not present (N(AI) < 20#/cm3)
          if ( sum(N(AI,:))/1.e6 .lt. 20._dp ) then
            lwc(1)  = 1.e-11_dp
            xaer(1) = 0._dp
          endif

          lwc(1) = min(lwc(1),1.e-5_dp)
          lwc(2) = min(lwc(2),1.e-5_dp)
          lwc(3) = min(lwc(3),1.e-5_dp)

          wascloud2=.true.

          !print *,'caw01',cwa01
          !print *,'caw02',cwa02
          !print *,'caw03',cwa03

          if (IDEB == 1) then
             write(12,fmt='(a,e10.4,a,e10.4,a,e10.4)' ) 'incloud lwc(CS) ',lwc(3), &
                  '  c(H2O)aq ',cwa03,'  m(H2O) ',MH2OTOT(CS) 
          endif


      else
!For (incloud.eq.0)  or  (incloud.eq.2).and.(RH.lt.rhactiv)


! aerosol water content NU ... CS
        do M=NU,CS
          MH2OTOT(M) = 0._dp 
          do I=1,IMAX
            CALL awater(RH,MASS(M,I,A_SUL),MASS(M,I,A_NH4),MASS(M,I,A_NIT),            &
                     MASS(M,I,A_SAL)+MASS(M,I,A_CHL),MASS(M,I,A_OR1)+MASS(M,I,A_OR2),  &
                     MASS(M,I,A_WAT))
            if (RH.LT.0.05) MASS(M,I,A_WAT)=0._dp
            MH2OTOT(M)=MH2OTOT(M)+MASS(M,I,A_WAT)
          end do
        end do

! evaporation of fog ("SINTEF fog simulation") or cloud ("cloud acivation")
        if (wascloud1) then
            MH2OTOT(:)=0._dp
            MASS(CS,1:IMAX,A_WAT)=0._dp
            wascloud1=.false.
        endif
        if (wascloud2) then
            MH2OTOT(:)=0._dp
            wascloud2=.false.
        endif

! droplet water for aqueous phase chemistry
        CALL partlwc(MH2OTOT,cwa01,cwa02,cwa03,xaer,lwc)

      endif



    end subroutine water_content


    subroutine partlwc(MH2OTOT,ch2oai,ch2oas,ch2ocs,xaera,lwca)
    !----------------------------------------------------------------------
    !
    !****  convert MASS(M,I,A_WAT) into c(ind_H2O_a**) in molec/cm^3(air)
    !      and liquid water content [m3/m3] for nucleation, Aitken and 
    !      accumulation mode     
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate c(ind_H2O_a**)
    !      calculate lwc [m3/m3]
    !
    !      interface
    !      ---------
    !
    !        input:
    !           MH2OTOT [ng/m^3]
    !    
    !      method
    !      ------
    !
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

    implicit none

     REAL( dp), dimension(MMAX),intent(in )       :: MH2OTOT

     REAL( dp),intent(out)                        :: ch2oai
     REAL( dp),intent(out)                        :: ch2oas
     REAL( dp),intent(out)                        :: ch2ocs
     REAL( dp),dimension(APN),intent(out)         :: lwca
     real( dp),dimension(APN),intent(out)         :: xaera

     INTEGER                                      :: jb

      ch2oai = (MH2OTOT(AI)/(M_H2O*1E9_dp))* N_A/1.E6_dp
      ch2oas = (MH2OTOT(AS)/(M_H2O*1E9_dp))* N_A/1.E6_dp
      ch2ocs = (MH2OTOT(CS)/(M_H2O*1E9_dp))* N_A/1.E6_dp


! several aq. phase modes
! jb=1: AI=3
! jb=2: AS=4
! jb=3: CS=5
      do jb=1,APN
        lwca(jb) = 1.E3_dp*(MH2OTOT(jb+2)/1.E9_dp)/1.E6_dp /RHOH2O
        !print *,'lwc',jb,MH2OTOT(jb+2),lwca(jb)
      end do
 

      if (ch2oai.le.0.0) then
        ch2oai   = 1.0E-32_dp
        lwca(1)  = 1.0E-32_dp
        xaera(1) = 0.0
      end if
      if (ch2oas.le.0.0) then
        ch2oas   = 1.0E-32_dp
        lwca(2)  = 1.0E-32_dp
        xaera(2) = 0.0         
      end if
      if (ch2ocs.le.0.0) then
        ch2ocs   = 1.0E-32_dp
        lwca(3)  = 1.0E-32_dp
        xaera(3) = 0.0
!  coarse aq. phase mode
!        lwca(1)=1.0E-32_dp
!        xaera(1)=0.0
      end if

    end subroutine partlwc


    subroutine wetdiameter(IMAX,RH,temp,MTOT,MTOTW,DPA,DPAW)
    !----------------------------------------------------------------------
    !
    !****  Calculate diameter of wet particle
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculate wet diameter
    !      Instantaneous adjustment to the ambient RH with no dynamical 
    !      growth.
    !
    !      interface
    !      ---------
    !
    !        input:
    !           rel. humidty   [-]
    !           temperature    [K]
    !           total mass dry [ng/m^3]
    !           total mass wet [ng/m^3]
    !           dry diameter   [m]
    ! 
    !    
    !      method
    !      ------
    !
    !      New (wet) diameter [m] after water added to particles
    !      Calculate Swelling Factor as proxy for diameter change
    !      SWF=DPA(WET)/DPA(DRY)=(1+MPTW/MPT)**(1/3)
    !      Fitzgerald et al. (1998)
    !      "A one-dimensional sectional model to simulate multicomponent
    !      aerosol dynamics in the marine boundary layer"
    !      J. Geophys. Res., 103, D13, 16085-16102, 1998.
    !      Equation (10)
    !      Modified with correction for particle curvature
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      Fitzgerald et al. (1998)
    !      "A one-dimensional sectional model to simulate multicomponent
    !      aerosol dynamics in the marine boundary layer"
    !      J. Geophys. Res., 103, D13, 16085-16102, 1998.
    !
    !------------------------------------------------------------------

    implicit none

    integer, intent(in)                               :: IMAX
    real( dp), intent(in)                             :: RH
    real( dp), intent(in)                             :: temp
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPA
    real( dp), dimension(MMAX),intent(in)             :: MTOT,MTOTW
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(out) :: DPAW
        
    real( dp), dimension(MMAX)             :: SWF
    real( dp)                              :: surf_H2O
    real( dp)                              :: sigma_H2O
    real( dp)                              :: A = 4.33e-6_dp
    real( dp)                              :: surfac
    real( dp)                              :: kelv
    real( dp)                              :: swfsum
    real( dp)                              :: swftest
    real( dp)                              :: TT
    integer                                :: M,I

      ! surface tension of water (in kg/s^2)
      surf_H2O= surf_h2o_std-0.155*(temp-273.15)
      sigma_H2O=surf_H2O*1.e-3_dp
      surfac=A*sigma_H2O/temp

      ! calculate swelling factor per mode
      swfsum=0.0_dp
      do M=1,MMAX
        SWF(M)=(MTOTW(M)/ (max(MTOT(M),1.e-32_dp)) )**(1./3.)
        SWF(M)=min(SWF(M),3.0_dp)
        SWF(M)=max(SWF(M),1.0_dp)
        swfsum=swfsum+SWF(M)
      end do 

!!!      SWF(NU)=max(SWF(NU), 0.5*SWF(NU) + 0.5*SWF(NA) )
!!!      SWF(NA)=max(SWF(NU), 0.5*SWF(NA) + 0.5*SWF(AI) )
!!!      SWF(AI)=max(SWF(NA), 0.5*SWF(AI) + 0.5*SWF(AS) )
!!!      SWF(AS)=max(SWF(AI),SWF(AS))
!!!      SWF(CS)=max(SWF(AS),SWF(CS))
!!!      ! avoid particle accumulation at NU-NA barrier
!!!      SWF(NU)=SWF(NA)

     ! Assume uniform SWF in all modes
     ! SWF is modified by curvature
      do M=1,MMAX
        do I=1,IMAX
          SWF(M)=swfsum/MMAX
          kelv=1._dp + ( surfac/(0.5*DPA(M,I)))   ! curvature (1+A/r)
          swftest=max( SWF(M)/kelv,1.0_dp)
          DPAW(M,I)=DPA(M,I) *swftest
          !  print *,'wetdp',M,I,DPA(M,I),SWF(M),swftest,DPAW(M,I)
        end do
      end do

    end subroutine wetdiameter
    
    
subroutine awater(aw,Msu,Mamsu,Mamnit,Msalt,Morgc,wh2o)
    !----------------------------------------------------------------------
    !
    !****  Water activity below DRH
    !
    !      author
    !      -------
    !      Dr. Francis S. Binkowski
    !
    !      contact
    !      -------
    !      Dr. Francis S. Binkowski
    !      Research Professor 
    !      Center for Environmental Modeling for Policy Development 
    !      Institute for the Environment 
    !      The University of North Carolina at Chapel Hill 
    !      Chapel Hill, NC 27599-1105 
    !      Email: frank_binkowski@unc.edu 
    !      Phone: (919) 966-2231 
    !      Fax: (919) 843-3113 
    !      http://ie.unc.edu/people/binkowski/ 
    !      Contact email for future reference:
    !      francisbinkowski@gmail.com
    !
    !
    !      purpose
    !      -------
    !      calculate water activity with polynomals from Tang and
    !      Munkelwitz, JGR 99: 18801-18808, 1994
    !      22.08.2012: added NaCl     (=SALT)
    !                  & added NaSucc (=ORG1)
    !
    !      interface
    !      ---------
    !
    !        input:
    !           aw     water activity         [-]
    !           Msu    mass conc. sulfate     [ng/m3]
    !           Mamsu  mass conc. ammonium    [ng/m3]
    !           Mamnit mass conc. nitrate     [ng/m3]
    !           Msalt  mass conc. sea-salt    [ng/m3]
    !           Morgc  mass conc. WSOC        [ng/m3]
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

    implicit none

      REAL( dp), intent(in)                :: Msu,Mamsu,Mamnit,Msalt
      REAL( dp), intent(in)                :: Morgc,aw
      REAL( dp), intent(out)               :: wh2o

      REAL( dp)                            :: mso4,mnh4,mno3,msal,awc
      REAL( dp)                            :: morg,torg
      REAL( dp)                            :: tso4,tnh4,tno3,tsal,x,null
      REAL( dp)                            :: mfs0,mfs1,mfs15,mfs2
      REAL( dp)                            :: c0(4),c1(4),c15(4),c2(4)
      REAL( dp)                            :: y, y0,y1,y15,y2,y3,y40,y5,y6
      REAL( dp)                            :: kSO4(6),kNO3(6),mfsSO4
      REAL( dp)                            :: kNaCl(6),kSucc(6)
      REAL( dp)                            :: mfsNO3,mfsNaCl,mfsSucc
      REAL( dp)                            :: u,y140,y1540,yc

! *** molecular weights:
      REAL( dp), parameter ::                           &
                mwh     = 1.0_dp,                       &
                mwso4   = 96.0636_dp,                   &
                mwnh4   = 18.0985_dp,                   &
                mwno3   = 62.0649_dp,                   &
                mwcl    = 35.45_dp,                     &
                mwna    = 22.99_dp,                     &
                mwsuc   = 118.0_dp,                     &
                mw2     = mwso4 + 2.0_dp * mwnh4,       &
                mwano3  = mwno3 + mwnh4,                &
                mwnacl  = mwcl + mwna

!     The polynomials use data for aw as a function of mfs from Tang and
!     Munkelwitz, JGR 99: 18801-18808, 1994.
!     The polynomials were fit to Tang's values of water activity as a
!     function of mfs.

! *** coefficients of polynomials fit to Tang and Munkelwitz data
!     now give mfs as a function of water activity.

      data c1/0.9995178_dp, -0.7952896_dp, 0.99683673_dp, -1.143874_dp/
      data c15/1.697092_dp,-4.045936_dp, 5.833688_dp, -3.463783_dp/
      data c2/2.085067_dp, -6.024139_dp, 8.967967_dp, -5.002934_dp/

! *** the following coefficients are a fit to the data in Table 1 of
!     Nair & Vohra, J. Aerosol Sci., 6: 265-271, 1975
!      data c0/0.8258941, -1.899205, 3.296905, -2.214749 /
! *** New data fit to data from
!       Nair and Vohra J. Aerosol Sci., 6: 265-271, 1975
!       Giaque et al. J.Am. Chem. Soc., 82: 62-70, 1960
!       Zeleznik J. Phys. Chem. Ref. Data, 20: 157-1200
      data c0/ 0.798079_dp, -1.574367_dp, 2.536686_dp, -1.735297_dp /

! *** polynomials for ammonium nitrate and ammonium sulfate are from:
!     Chan et al.1992, Atmospheric Environment (26A): 1661-1673.

      data kNO3/0.2906_dp, 6.83665_dp, -26.9093_dp,   &
               46.6983_dp, -38.803_dp, 11.8837_dp/
      data kSO4/ 2.27515_dp, -11.147_dp, 36.3369_dp,  &
             -64.2134_dp, 56.8341_dp, -20.0953_dp/

! *** polynomials for sodium chloride are from:
!     Tang, I.N., Tridico, A.C., Fung, K.H.
!     Thermodynamic and optical properties of sea salt aerosols
!     J. Geophys. Res., 102(D19), 23,269-23,275, 1997
!     polynomial fit according to Chan et al. 1992 method

      data kNaCl/1.2617_dp, -4.3239_dp, 10.658_dp,    &
                -16.248_dp,  13.324_dp, -4.6718_dp/

! *** polynomials for sodium succinate are from:
!     Peng, C. and C. K. Chan
!     The water cycles of water-soluble organic salts of atmospheric importance
!     Atmos. Environ., 35, 1183-1192, 2001

      data kSucc/0.38225_dp, 5.4575_dp, -21.965_dp,   &
                34.34_dp, -23.764_dp, 5.5627_dp/


! *** check range of per cent relative humidity
!       aw water activity = fractional relative humidity
! calculate masses of so4, nh4 and no3 from the masses
! of h2so4,(nh4)2so4 and nh4no3 in ng/m^3
       mso4 = Msu/(2.*mwh+mwso4)*mwso4 + Mamsu/mw2*mwso4
       mnh4 = Mamsu/mw2*mwnh4 + Mamnit/mwano3*mwnh4
       mno3 = Mamnit/mwano3*mwno3
       msal = Msalt/mwnacl*mwcl
       morg = Morgc/mwsuc

      null=0.0_dp
      mfs0=0.0_dp
      mfs1=0.0_dp
      mfs15=0.0_dp
      tso4 = max( mso4 , null )
      tnh4 = max( mnh4 , null )
      tno3 = max( mno3 , null )
      tsal = max( msal , null )
      torg = max( morg , null )
      x = 0.0_dp
! *** if there is non-zero sulfate calculate the molar ratio
      if (tso4 .gt. 0.0_dp ) then
        x = tnh4 / tso4
      else
! *** otherwise check for non-zero nitrate and ammonium
        if ( (tno3 .gt. 0.0_dp) .and. (tnh4 .gt. 0.0_dp) ) x = 10.0_dp
      end if
      y  = 0.0_dp
      y2 = 0.0_dp
      y3 = 0.0_dp
      y5 = 0.0_dp
      y6 = 0.0_dp

! *** begin screen on x for calculating wh2o
      if ( x .lt. 1.0_dp ) then
!
          mfs0 = poly4(c0(1),c0(2),c0(3),c0(4),aw)
          mfs1 = poly4(c1(1),c1(2),c1(3),c1(4),aw)
          y0 = (1.0_dp - mfs0 ) / mfs0
          y1 = (1.0_dp - mfs1 ) / mfs1
          y = (1.0_dp- x) * y0 + x * y1
!
       else if ( x .lt. 1.5_dp) then
!
         if ( aw .ge. 0.40_dp ) then
            mfs1  = poly4(c1(1),c1(2),c1(3),c1(4),aw)
            mfs15 = poly4(c15(1),c15(2),c15(3),c15(4),aw)
            y1  = (1.0_dp - mfs1 ) / mfs1
            y15 = (1.0_dp - mfs15) / mfs15
            y = 2.0_dp * ( y1 * (1.5_dp - x) + y15 *( x - 1.0_dp) )
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
           awc = 0.80_dp * (x - 1.0_dp) ! rh along the crystallization curve.
           u = 0.4_dp
           y = 0.0_dp
           if ( aw .ge. awc ) then ! interpolate using crystalization curve
               mfs1  = poly4(c1(1),c1(2),c1(3),c1(4),u)
               mfs15 = poly4(c15(1),c15(2),c15(3),c15(4),u)
               y140  = (1.0_dp - mfs1 ) / mfs1
               y1540 = (1.0_dp - mfs15) / mfs15
               y40 = 2.0_dp * ( y140 * (1.5_dp - x) + y1540 *( x - 1.0_dp) )
               yc = 2.0_dp * y1540 * (x -1.0_dp) ! y along crystallization curve
               y = y40 - (y40 - yc) * (0.40_dp-aw) / (0.40_dp - awc)
            end if ! end of checking for aw
          end if ! end of checking on irh

       else if( x .lt. 1.9999_dp) then
!
           y= 0.0_dp
           if( aw .ge. 0.40_dp) then
             mfs15 = poly4(c15(1),c15(2),c15(3),c15(4),aw)
             mfs2  = poly4(c2(1),c2(2),c2(3),c2(4),aw)
             y15 = (1.0_dp - mfs15) / mfs15
             y2  = (1.0_dp - mfs2) / mfs2
             y = 2.0_dp * (y15 * (2.0_dp - x) + y2 * (x - 1.5_dp) )
           end if ! end of check for crystallization
!
      else ! 1.9999_dp < x

! regime where ammonium sulfate and ammonium nitrate are in solution.
!
! *** following cf&s for both ammonium sulfate and ammonium nitrate
! *** check for crystallization here. their data indicate a 40% value
!     is appropriate.
            y2 = 0.0_dp
            y3 = 0.0_dp
            if ( aw .ge. 0.40_dp) then
              mfsSO4 = poly6(kSO4(1),kSO4(2),kSO4(3),kSO4(4),kSO4(5),kSO4(6),aw)
              mfsNO3 = poly6(kNO3(1),kNO3(2),kNO3(3),kNO3(4),kNO3(5),kNO3(6),aw)
              y2 = (1.0_dp - mfsSO4) / mfsSO4
              y3 = (1.0_dp - mfsNO3) / mfsNO3

            end if
!
      end if ! end of checking on x

! for sea-salt, represented by NaCl
! CRH(NaCl): 0.47 (crystallisation below RH=47%)
      if ( aw .ge. 0.47_dp) then
         mfsNaCl = poly6(kNaCl(1),kNaCl(2),kNaCl(3),kNaCl(4),kNaCl(5),kNaCl(6),aw)
         y5 = (1.0_dp - mfsNaCl) / mfsNaCl
      end if

! for organics, represented by Na-succinate
! CRH(NaSucc): 0.48 (crystallisation below RH=48%)
      if ( aw .ge. 0.48_dp) then
         mfsSucc = poly6(kSucc(1),kSucc(2),kSucc(3),kSucc(4),kSucc(5),kSucc(6),aw)
         y6 = (1.0_dp - mfsSucc) / mfsSucc
      end if

! *** now set up output of wh2o

!      wh2o units are nanograms (liquid water) / cubic meter of air
!
      if ( x .lt. 1.9999_dp) then
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

    !------------------------------------------------------------------
    !
    !  Polynomial functions POLY4 and POLY6
    !
    !      author
    !      -------
    !      Written by 
    !      Dr. Francis S. Binkowski 
    !      Research Professor 
    !      Center for Environmental Modeling for Policy Development 
    !      Institute for the Environment 
    !      The University of North Carolina at Chapel Hill 
    !      Chapel Hill, NC 27599-1105 
    !      Email: frank_binkowski@unc.edu 
    !      Phone: (919) 966-2231 
    !      Fax: (919) 843-3113 
    !      http://ie.unc.edu/people/binkowski/ 
    !      Contact email for future reference:
    !      francisbinkowski@gmail.com
    !
    !------------------------------------------------------------------

        ELEMENTAL REAL(dp) FUNCTION  poly4(A1,A2,A3,A4,X)
          IMPLICIT NONE
        REAL( dp), intent(in)   :: A1,A2,A3,A4
        REAL( dp), intent(in)   :: X
        poly4 = A1 + X *( A2 + X*( A3 + X*A4))
        END FUNCTION poly4

        ELEMENTAL REAL(dp) FUNCTION poly6(A1,A2,A3,A4,A5,A6,X)
          IMPLICIT NONE
        REAL( dp), intent(in)   :: A1,A2,A3,A4,A5,A6
        REAL( dp), intent(in)   :: X
        poly6 = A1 + X*( A2 + X*( A3 + X*( A4 +   &
               X * ( A5 + X * (A6  )))))
        END FUNCTION poly6
    !------------------------------------------------------------------

end module gde_addwater
