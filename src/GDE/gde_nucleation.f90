! <gde_nucleation.f90 - A component of the Multicomponent
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
!*
!*    All routines written by Matthias Karl
!*    except the following routines:
!*    routine YU_TIMN_NUCLEATION written by Fangqun Yu
!*    routine ACDC_NUCLEATION written by Shuai Jiang
!*    routine TERNARY_FIT written by Joonas Merikanto,
!*            with corrections by Matthias Karl
!*    routine NEWBINAPARA written by Anni Maattanen et al.
!*    routine BINHOMOGENEOUS written by Hanna Vehkamaeki, and belonging
!*    routines PARTIAL_MOLAR_VOL, ACTIV_BINARY, ZELEZNIK 
!*            written by Hanna Vehkamaeki, and belonging
!*    functions ACTIVITY_ACID and ACTIVITY_WATER written by Hanna Vehkamaeki.
!*
!*****************************************************************************!
module gde_nucleation

  IMPLICIT NONE
  INTRINSIC :: SELECTED_REAL_KIND

  private
 
  public :: nucleation, nucleationratio
  
  !public :: kinetic, activation
  !public :: kinetic_iodine, activation_iodine
  !public :: organcluster
  !public :: aminecluster
  !public :: ioncluster, ioninduced
  !public :: yu_timn_nucleation
  !public :: ACDC_nucleation
  !public :: ternary_fit

! KPP DP - Double precision kind
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)
  
contains

 subroutine nucleation(INUCMEC,cair,ch2so4,cnh3,camin,cnit,csoa2,ciod, &
                       temp,RH,jrno2, CA,KP,fnuc,coags,daynr,          &
                       Lati,natot,JNUC)
    !********************************************************************
    !
    !     N  U  C  L  E  A  T  I  O  N
    !
    !********************************************************************
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    ! Nucleation module
    !
    !          - nucleation             INUC=1
    !               kinetic               H2SO4         (INUCMEC=1) or
    !               homogeneous           H2SO4-H2O     (INUCMEC=2)
    !               homogeneous           H2SO4-H2O-NH3 (INUCMEC=3)
    !               ion-mediated          H2SO4-H2O-NH3 (INUCMEC=4)
    !               activation            H2SO4         (INUCMEC=5)
    !               kinetic               HNO3-AMIN     (INUCMEC=6)
    !               combination           H2SO4-H2O     (INUCMEC=7)
    !               "OS1" activation      H2SO4-ORG     (INUCMEC=8)
    !               "OS2" kinetic         H2SO4-ORG     (INUCMEC=9)
    !               "OS3" total           H2SO4-ORG     (INUCMEC=10)
    !               neutral & ion-induced H2SO4-H2O     (INUCMEC=11)
    !               diesel-exhaust        H2SO4-ORG     (INUCMEC=12)
    !               ACDC cluster          H2SO4-H2O-NH3 (INUCMEC=13)
    !               HIO3 kinetic/activ    HIO3          (INUCMEC=14)
    !
    !  INPUT
    !  -----
    !
    ! INUCMEC: nucleation option switch
    !    cair: mixing ratio of air (M)        [1/m^3]
    !  ch2so4: concentration of h2so4 vapour  [1/m^3]
    !    cnh3: concentration of nh3 vapour    [1/m^3]
    !   camin: concentration of amine         [1/m^3]
    !    cnit: concentration of nitric acid   [1/m^3]
    !   csoa2: concentration of SOA-2         [1/m^3]
    !    ciod: concentration of HIO3          [1/m^3]
    !    temp: air temperature                [K]
    !      RH: rel. humidity                  [-]
    !   jrno2: phot. freq. NO2                [s^-1]
    !      CA: amine nucleation param 
    !      KP: amine nucleation param
    !    fnuc: scal. factor of nucl. constant [-]
    !   coags: coagulation sink               [s^-1]
    !   daynr: day of year                    [-]
    !    Lati: latitude                       [-] 
    ! 
    !  OUTPUT
    !  ------
    !
    !   JNUC: nucleation rate J in m^-3s^-1 
    !  natot: number of H2SO4 molecs in crit cluster
    !
    !********************************************************************
    !* Lehtinen et al. (2007) to extrapolate to 1 nm clusters
    !* for nucleation option 11 newbinapara
    !*   Lehtinen, K. E., Dal Maso, M., Kulmala, M., and Kerminen, V. M., 
    !*   Journal of Aerosol Science, 38(9), 988-994, 2007
    !*   dc=0.7               ! Target size (nm)
    !*   (in geometric diameter = mobility diameter -0.3nm)
    !*
    !********************************************************************

     implicit none

        integer, intent(in)      :: INUCMEC
        REAL( dp), intent(in)    :: ch2so4,cnh3,cair
        REAL( dp), intent(in)    :: camin,cnit,csoa2,ciod
        REAL( dp), intent(in)    :: temp,RH,jrno2,CA,KP
        REAL( dp), intent(in)    :: fnuc,daynr,Lati,coags        
        REAL( dp), intent(out)   :: JNUC,natot

        !Local
        REAL( dp)                :: nwtot,nntot
        REAL( dp)                :: rc,natoti,natott
        REAL( dp)                :: JNUCI,JNUCT

        !MAATTANEN ET AL. (2018)
        REAL( dp)                :: rc_n     ! radius of the charged critical cluster [m]
        REAL( dp)                :: rc_i     ! radius of the charged critical cluster [m] 
        REAL( dp)                :: jnuc_n   ! Neutral nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
        REAL( dp)                :: jnuc_i   ! Charged nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
        REAL( dp)                :: na_n     ! sulfuric acid molecules in the neutral critical cluster
        REAL( dp)                :: na_i     ! sulfuric acid molecules in the charged critical cluster
        !LEHTINEN ET AL. (2007)
        REAL( dp)                :: cs       ! H2SO4 condensation sink (h-1)
        REAL( dp)                :: gr       ! Particle growth rate (nm/h) 
        REAL( dp)                :: gammax
        REAL( dp)                :: m,dx,d1
        REAL( dp)                :: LKK_n,LKK_i

!comment
        REAL( dp)                :: jnuc_n_17,jnuc_i_17       ! Nucleation rates at 1.7nm mobility diameter



         SELECT CASE (INUCMEC)
            CASE (1)

            ! Kinetic nucleation of sulfuric acid
               CALL kinetic(ch2so4,cnh3,cair,JNUC,natot)
               !write(6,*) NVAP,c(KPP_NH3),cair,JNUC,natot

            CASE (2)

            ! Binary homogeneous nucleation of sulfuric acid and water
            ! based on Vehkamaeki et al.(2002 & 2003)
               CALL binhomogeneous(temp,RH,ch2so4,nwtot,natot,rc,JNUC)

            CASE (3)

            ! Ternary homogeneous nucleation, fitted by Joonas Merikanto, 2006
            !             including stable ammonium bisulfate formation
            !             output is logarithm of nucleation rate (in cm^-3s^-1)
               CALL ternary_fit(temp,RH,ch2so4,cnh3,cair,JNUC,natot,nntot,rc)
               JNUC=DEXP(JNUC)
            ! MSK 06.09.2017 must be multiplied by 10^6 
            ! to get nucleation rate J in m^-3s^-1
               JNUC=JNUC*1.e6
               write(6,*) ' JNUC        natot'
               write(6,'(4ES12.4)') JNUC, natot

            CASE (4)

            ! Revised based on Yu et al. (2020) look-up tables
            ! for ternary ion-mediated nucleation of H2SO4-H2O-NH3
            !      X0 = [H2SO4] in #/cm3  (5E5-5E9)
            !      Y0 = RH in % (0.5-99.5)
            !      Z0 = T (in K) (190-304)
            !      U0 = Q = ionization rate (ion-pairs cm-3 s-1) (2-100)
            !      V0 = S = surface area (um2/cm3) (20-200)
            !      W0 = B = NH3 in #/cm3  (1E5-1E12)
            ! fixed U0=2.0 ion-pairs/cm3s and V0=20.0 um2/cm3
            !   IF(ch2so4.GT.5.0E+11_dp) THEN
            ! MSK 30.12.2020: just TIMN
                 CALL yu_timn_nucleation(ch2so4,RH,temp,cnh3,JNUC,natot)
                 write(6,*) 'J_TIMN',JNUC,natot

            !     IF(JNUC.LT.1000._dp) CALL ioncluster(ch2so4,JNUC,natot)
            !   ELSE
            ! Ion-mediated nucleation, Yu and Turco 2001
            !             based on kinetic collision approach, QSSA assumption
            !     CALL ioncluster(ch2so4,JNUC,natot)
            !   ENDIF

            CASE (5)

            ! Activation nucleation, Kulmala et al., 2006
            !             activation of clusters of one molecule of H2SO4
               CALL activation(fnuc,ch2so4,JNUC,natot)

            CASE (6)

            ! Nucleation of nitric acid/amine, photolysis-triggered
               IF (jrno2.GT.1e-5) THEN
                 CALL aminecluster(CA,KP,cnit,camin,JNUC,natot)
               ELSE
                 JNUC=0._dp
                 natot=1._dp
               ENDIF

            CASE (7)

            ! Combined nucleation mechanism
               CALL ioncluster(ch2so4,JNUCI,natoti)
               CALL activation(fnuc,ch2so4,JNUCT,natott)
               natot=natott+natoti
               JNUC=JNUCI+JNUCT
               !write(6,*) ' JNUC        natot'
               !write(6,'(4ES12.4)') JNUC, natot

            CASE (8)

            ! Nucleation of sulfuric acid/organic acid clusters "OS1"
               natot=1.
               JNUC=0.7E-7*ch2so4+0.7E-7*csoa2

            CASE (9)

            ! Nucleation of sulfuric acid/organic acid clusters "OS2"
               natot=2.
               JNUC=3.4E-20*(ch2so4**2.+ch2so4*csoa2)

            CASE (10)

            ! Nucleation of sulfuric acid/organic acid clusters "OS3"
               CALL organcluster(fnuc,ch2so4,csoa2,temp,JNUC,natot)

            CASE (11)

            ! New parameterization of neutral and ion-induced sulfuric acid-
            ! water nucleation by Maattanen et al. (2018)
               CALL newbinapara(temp,RH,ch2so4,cair, jnuc_n,jnuc_i, na_n,na_i, rc_n,rc_i)

               ! Nucleation rates at 1.0nm using Lehtinen et al. (2007)
               m  = -1.6
               ! Particle growth rate (nm/h); Nieminen et al., 2010
               gr = (ch2so4*1.E-6_dp)/(661.1*(RH*100)**2-1.129E5*(RH*100)+1.549E7)
               ! Typical CS value in atmosphere in 1/h 
               cs = 10
               ! Target size (in geometric diameter = mobility diameter -0.3nm)
               dx= 0.7       ! to get 1 nm nucleated particles
               !!! Neutral case
               ! diameter of critical cluster
               ! MSK 25.02.2022: rc_n is already in nm
               d1=2.*rc_n
               ! gamma-factor in Lehtinen et al.
               gammax=max(0.0_dp,1.0_dp/(m+1)*( (dx/(d1))**(m+1)-1))
               ! final scaling factor in Lehtinen et al.
               LKK_n=exp(-gammax*d1*cs/gr)
               !!! Charged case
               ! MSK 25.02.2022: rc_i is already in nm
               d1=2.*rc_i
               gammax=max(0.0_dp,1.0_dp/(m+1)*( (dx/(d1))**(m+1)-1)) 
               LKK_i=exp(-gammax*d1*cs/gr)

            ! Nucleation rate of neutral and charged case should be added
            ! according to pers. commun. Maattanen (2020)
            ! and multiplied by 10^6 to get nucleation rate J in m^-3s^-1
               JNUC=jnuc_n*LKK_n + jnuc_i*LKK_i
               ! MSK 25.02.2022: na_n only added if jnuc_n>0
               natot=na_i
               if (jnuc_n > 0) natot=natot+na_n
               JNUC=JNUC*1.e6

             !  write(6,*) ' JNUC        natot'
             !  write(6,'(4ES12.4)') JNUC, natot
             
            CASE (12)

            ! Nucleation of sulfuric acid/organic acid clusters "OS3"
               CALL diesel(fnuc,ch2so4,csoa2,JNUC,natot)

            CASE (13)

            ! ACDC implementation in a look-up table way for ambient observation simulation
            ! by Shuai Jiang (2017)
!               write(6,*) 'ch2so4,cnh3,temp,RH,coags', ch2so4,cnh3,temp,RH,coags
               CALL ACDC_nucleation(ch2so4,cnh3,temp,RH,coags,JNUC,natot)
               write(6,*) ' JNUC        natot'
               write(6,'(4ES12.4)') JNUC, natot

            ! CASE (1x)
            ! ACDC implementation in a direct way for laminar flow tube simulation
            ! by Shuai Jiang (2017)
            !   CALL ACDC_direct(ch2so4,cnh3,temp,RH,coags,JNUC,natot)

            CASE (14)

            ! Kinetic nucleation of iodic acid
            !   CALL kinetic_iodine(ciod,JNUC,natot)
            ! Activation of clusters of one molecule HIO3
            !   CALL activation_iodine(ciod,JNUC,natot)
            ! Nucleation of sulfuric acid/iodic acid clusters
            ! Vuollekoski et al. (2009), Equation (3)
            ! They assume that diameter of nucleated particles is 1.5 nm
               natot=1.
               JNUC=1.e-18 * ch2so4 * ciod
            !   write(6,'(6ES12.4)') ciod, JNUC, natot

            CASE DEFAULT
 
              write(6,*) 'no valid nucleation option'
               STOP   

         END SELECT 

  end subroutine nucleation
  


  ! *********************************************************************
  ! *********************************************************************
  ! *********************************************************************


     subroutine nucleationratio(DNUC,VVAP,nucl)
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
    !      calculate volume ratio particle/vapor
    !
    !      interface
    !      ---------
    !
    !        input:
    !           DNUC   diameter nucl particle [m]
    !           VVAP   vapor molec volume     [m^3] 
    !
    !      method
    !      ------
    !      calculate volume ratio particle/vapor
    !
    !      reference
    !      ---------
    !      none
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

    use gde_constants,    only      : pi

    implicit none

    ! input
     REAL( dp), intent(in)    :: DNUC        ! diameter nucl particle [m]
     REAL( dp), intent(in)    :: VVAP        ! vapor molec volume [m^3]
     REAL( dp), intent(out)   :: nucl
     REAL( dp)                :: VNUC        ! nucl particle volume [m^3]

     VNUC=(pi/6._dp)*DNUC**3._dp
     nucl=VNUC/VVAP                          ! volume ratio particle/vapor [-]

   end subroutine nucleationratio
  

     subroutine kinetic(c2,c3,cair,jnuc,natot)
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
    !      kinetic nucleation: new parameterisation
    !
    !      interface
    !      ---------
    !
    !        input:
    !           c2: concentration of h2so4 vapour [1/m^3]
    !           c3: concentration of nh3 vapour [1/cm^3]
    !           natot: total number of h2so4 molecules 
    !                  in the critical cluster
    !
    !      method
    !      ------
    !      kinetic nucleation: empirical parameterisation
    !      based on Kulmala, ACP 2006.
    !      Original paramters from Lehtinen:
    !      JNUC=CNUC*NVAP**2._dp, with CNUC=1.E-20_dp
    !
    !      reference
    !      ---------
    !      Nucleation mechanism proposed in M. Kulmala ACPD,5,11277-11293,2005
    !      and Kulmala et al.,ACP,6,787-793,2006
    !      Parameter C_2 from Riipinen et al.,ACP,7,1899-1914,2007
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

     implicit none

        REAL( dp), intent(in)    :: c2,c3,cair
        REAL( dp), intent(out)   :: jnuc,natot
        REAL( dp)                :: C_2

        IF (c3 .GT. cair/1.e10_dp) THEN
          C_2=2.3e-20_dp         !m^3s^-1
        ELSE
          C_2=3.2e-20_dp         !m^3s^-1 Median value BACCI/QUEST IV (max=1.8E-19)
        ENDIF
        natot=2._dp
        jnuc=C_2*(c2**2._dp)       !m^-3s^-1

   end subroutine kinetic


     subroutine kinetic_iodine(c2,jnuc,natot)
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
    !      kinetic nucleation of iodic acid
    !
    !      interface
    !      ---------
    !
    !        input:
    !           c2: concentration of hio3 vapour [1/m^3]
    !           natot: total number of hio3 molecules 
    !                  in the critical cluster
    !
    !      method
    !      ------
    !      kinetic nucleation: empirical parameterisation
    !      OIO nucleation parameterization of Vuollekoski, JGR, 2009
    !
    !      reference
    !      ---------
    !      Nucleation mechanism proposed in M. Kulmala ACPD,5,11277-11293,2005
    !      and Kulmala et al.,ACP,6,787-793,2006
    !      Parameter K from:
    !      Vuollekoski, H., Kerminen, V.-M., Anttila, T., Sihto, S.-L.,
    !      Vana, M., Ehn, M., Korhonen, H., McFiggans, G., O'Dowd, C. D.
    !      Kulmala, M.: Iodine dioxide nucleation simulations in coastal
    !      and remote marine environments, J. Geophys. Res., 114, D02206, 
    !      doi:10.1029/2008JD010713, 2009.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

     implicit none

        REAL( dp), intent(in)    :: c2
        REAL( dp), intent(out)   :: jnuc,natot
        REAL( dp)                :: K

        ! kinetic nucleation of OIO, Vuollekoski (2009):
        ! best choice is K = 1.e-12 cm^3s^-1
        ! maximum K value from kinetic gas theory 
        ! K = 1.e-10 cm^3s^-1
        !K=1.0e-21_dp             !m^3s^-1
        K=1.0e-19_dp             !m^3s^-1

        natot=2._dp
        jnuc=K*(c2**2._dp)       !m^-3s^-1

   end subroutine kinetic_iodine


     subroutine activation(fnuc,c2,jnuc,natot)
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
    !      cluster activation: new parameterisation
    !
    !      interface
    !      ---------
    !
    !        input:
    !           c2: concentration of h2so4 vapour [1/m^3]
    !           natot: total number of h2so4 molecules 
    !                  in the critical cluster
    !
    !
    !      method
    !      ------
    !      cluster activation: empirical parameterisation
    !      based on Kulmala, ACP 2006.
    !
    !      reference
    !      ---------
    !      Nucleation mechanism proposed in M. Kulmala ACPD,5,11277-11293,2005
    !      and Kulmala et al.,ACP,6,787-793,2006
    !      Activation of clusters of one molecule of H2SO4
    !      Parameter C_1 from Riipinen et al.,ACP,7,1899-1914,2007
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

     implicit none

        REAL( dp), intent(in)    :: fnuc,c2
        REAL( dp), intent(out)   :: jnuc,natot
        REAL( dp)                :: C_1

        C_1=2.4e-7_dp        !s^-1  Median value BACCI/QUEST IV (max=2.0E-6)
        natot=1._dp
        jnuc=fnuc*C_1*c2          !m^-3s^-1     

   end subroutine activation


     subroutine activation_iodine(c2,jnuc,natot)
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
    !      cluster activation of iodic acid
    !
    !      interface
    !      ---------
    !
    !        input:
    !           c2: concentration of hio3 vapour [1/m^3]
    !           natot: total number of hio3 molecules 
    !                  in the critical cluster
    !
    !
    !      method
    !      ------
    !      cluster activation: empirical parameterisation
    !      OIO nucleation parameterization of Vuollekoski, JGR, 2009
    !
    !      reference
    !      ---------
    !      Nucleation mechanism proposed in M. Kulmala ACPD,5,11277-11293,2005
    !      and Kulmala et al.,ACP,6,787-793,2006
    !      Parameter A from:
    !      Vuollekoski, H., Kerminen, V.-M., Anttila, T., Sihto, S.-L.,
    !      Vana, M., Ehn, M., Korhonen, H., McFiggans, G., O'Dowd, C. D.
    !      Kulmala, M.: Iodine dioxide nucleation simulations in coastal
    !      and remote marine environments, J. Geophys. Res., 114, D02206, 
    !      doi:10.1029/2008JD010713, 2009.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

     implicit none

        REAL( dp), intent(in)    :: c2
        REAL( dp), intent(out)   :: jnuc,natot
        REAL( dp)                :: A

        ! cluster activation of OIO, Vuollekoski (2009):
        ! best choice is A = 3.e-6 s^-1
        !A=3.0e-6_dp        !s^-1
        ! better suited for arctic conditions
        A=0.5e-6_dp        !s^-1
        natot=1._dp
        jnuc=A*c2          !m^-3s^-1     

   end subroutine activation_iodine


     subroutine organcluster(fnuc,c2,c3,temp,jnuc,natot)
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
    !     sulfphuric acid - organic acid cluster nucleation
    !
    !      interface
    !      ---------
    !
    !        input:
    !           c2: concentration of h2so4 vapour    [1/m^3]
    !           c3: concentration of organic vapour  [1/m^3]
    !           natot: total number of molecules 
    !                  in the critical cluster
    !
    !      method
    !      ------
    !      sulfphuric acid - organic acid cluster nucleation
    !      Original paramters from Lehtinen:
    !      JNUC=CNUC*NVAP**2._dp, with CNUC=1.E-20_dp
    !
    !      reference
    !      ---------
    !      Zhao et al., J Phys Chem A, 113, 680-689, 2009
    !      Binding energy of 17-18 kcal/mol for h2so4 - benzoic acid and
    !                                           h2so4 - CPA clusters   
    ! 
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

     implicit none

        REAL( dp), intent(in)    :: fnuc,c2,c3,temp
        REAL( dp), intent(out)   :: jnuc,natot
        REAL( dp)                :: DENSPAR,MONOVOL,beta11
        REAL( dp)                :: clusterc
 
        ! calculate collision rate of two gas-phase h2so4-org (1:1) clusters
        !   beta11 = 4.53E-16 m^3s^-1 at 280K
        DENSPAR= 1000._dp
        MONOVOL=98.08_dp*1.661e-27/1.83e3
        beta11= 4*2.**0.5*(3./4./3.1416)**(1./6.)*(6.*1.38d-23*temp/   &
            denspar)**0.5*(MONOVOL)**(1./6.)    !m^3s^-1
        ! typical "field" collision rate
        beta11=1.e-20_dp              !m^3s^-1
        if (c3.gt.c2) then
          clusterc=c2+c3              !m^-3
          natot=2._dp
          beta11=0.5*beta11
        else
          ! kinetic h2so4 dimer nucleation
          clusterc=c2
          natot=2._dp
        endif
        jnuc=fnuc*beta11*(clusterc**2._dp)       !m^-3s^-1

   end subroutine organcluster


     subroutine diesel(fnuc,c2,c3,jnuc,natot)
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
    !      Empirical parameterization for nucleation in diesel exhaust
    !
    !      interface
    !      ---------
    !
    !        input:
    !           c2: concentration of h2so4 vapour    [1/m^3]
    !           c3: concentration of organic vapour  [1/m^3]
    !           natot: total number of molecules 
    !                  in the critical cluster
    !
    !      method
    !      ------
    !      Heteromolecular nucleation. K1 and K2 from L. Pirjola, 
    !      pers. commun., 2014
    !
    !      reference
    !      ---------
    !      Ronkko et al., Effect of gaseous sulphuric acid on diesel exhaust
    !      nanoparticle formation and characteristics, EST 47, 11882-11889, 
    !      2013
    ! 
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

     implicit none

        REAL( dp), intent(in)    :: fnuc,c2,c3
        REAL( dp), intent(out)   :: jnuc,natot
        real( dp)                :: K1,K2

        ! K1, K2 from least square fit to experimental diesel exhaust data
        K1 = 3.8e-23_dp   !m^3s^-1
        K2 = 5.6e-23_dp   !m^3s^-1
        
        jnuc  = (K1*c2*c2) + (K2*c2*c3)  !m^-3s^-1
        jnuc=fnuc*jnuc
        natot = 2._dp 


   end subroutine diesel


     subroutine aminecluster(CNUC,kpeq,c2,c3,jnuc,natot)
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
    !      kinetic nucleation of amine clusters
    !
    !      interface
    !      ---------
    !
    !        input:
    !           c2: concentration of hno3 vapour  [1/m^3]
    !           c3: concentration of amine vapour [1/m^3]
    !           kpeq: equlibrium partition coeff. for 
    !                 amine+HNO3 reaction [(1/cm^3)^2)
    !           natot: total number of amine molecules 
    !                  in the critical cluster
    !
    !
    !      method
    !      ------
    !      empirical parameterisation for amine clusters
    ! 
    !      Original paramters from Lehtinen:
    !      JNUC=CNUC*NVAP**2._dp, with CNUC=1.E-20_dp
    !
    !      reference
    !      ---------
    !      Nucleation mechanism proposed in M. Kulmala ACPD,5,11277-11293,2005
    !      and Kulmala et al.,ACP,6,787-793,2006
    !      Parameter C_2 from Riipinen et al.,ACP,7,1899-1914,2007
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

     implicit none

        REAL( dp), intent(in)    :: c2,c3,kpeq,CNUC
        REAL( dp), intent(out)   :: jnuc,natot
        REAL( dp)                :: kpeqm
 
        kpeqm=kpeq*1.e12_dp
        !C_2=1.E-28_dp     !1.E-25
        natot=1._dp

        ! Threshold: c(amin)*c(hno3) > Kp*10

        IF ((c2*c3) .gt. kpeqm) THEN
        !IF ((c2*c3) .gt. kpeqm*10.) THEN
          !jnuc=C_2*(c3**2._dp)       !m^-3s^-1
          jnuc=CNUC*(c2*c3)
        ELSE
          jnuc=0._dp
        ENDIF
        !write(6,*) 'ami',c2,c3,kpeqm


   end subroutine aminecluster
  

     subroutine ioncluster(c2,J_ion,natot)
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
    !      ion mediated nucleation mechanism
    !
    !      interface
    !      ---------
    !
    !        input:
    !           c2: concentration of h2so4 vapour [1/m^3]
    !           natot: total number of h2so4 molecules 
    !                  in the critical cluster
    !
    !
    !      method
    !      ------
    !      Based on kinetic collision approach, QSSA assumption
    !      Yu and Turco, JGR, 106(D5), 4797-4814, 2001
    !      Equation 13
    !
    !      Assumptions:
    !      - QSSA for cluster collision
    !      - 3 H2SO4 molecules in a critical cluster (n=3) 
    !      - GCR induced constant ionization rate Q=2 ion pairs cm^-3s^-1
    !      - assume same kf and kr for all cluster sizes
    !      - neglecting cluster dissociation (only ion-ion recombination)
    !
    !      wet diameter of hydrated H2SO4 molecule: 0.65 nm
    !      critical diameter of neutral cluster:    1.8 nm
    !
    !      reference
    !      ---------
    !      Yu and Turco, JGR, 106(D5), 4797-4814, 2001
    !      Karl, M., Gross, A., Pirjola, L., and C. Leck, A new flexible 
    !      multicomponent model for the study of aerosol dynamics in the 
    !      marine boundary layer, Tellus B, 63(5), 1001-1025, 
    !      doi: 10.1111/j.1600-0889.2011.00562.x, 2011.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

     implicit none

        REAL( dp), intent(in)    :: c2
        REAL( dp), intent(out)   :: J_ion,natot
        REAL( dp)                :: cacid,qgcr,kc_forw,kc_reco

        natot=3._dp
        qgcr=2.2_dp                    ! ion pairs cm^-3s^-1
        kc_forw=6.E-10_dp              ! cm^-3s^-1
        kc_reco=1.E-6_dp               ! cm^-3s^-1 
        cacid=c2*1.e-6_dp
       
        ! units: cm^-3s^-1
        J_ion = qgcr*(1./(1.+ (SQRT(qgcr)*SQRT(kc_reco)) /       &
            (kc_forw*cacid) ))**4

        J_ion = J_ion * 1.e6_dp      !m^-3s^-1

   end subroutine ioncluster

 

      subroutine yu_timn_nucleation(X0,Y0,Z0,W0,XJTIM,natot)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Fangqun Yu
    !
    !      contact
    !      -------
    !      Dr. Fangqun Yu
    !      Professor
    !      Atmospheric Sciences Research Center
    !      State University of New York at Albany
    !      251 Fuller Road Albany, New York 12203-3649
    !      Tel: 518-437-8767; Fax: 518-437-8758
    !      E-mail: fyu@albany.edu
    !
    !      purpose
    !      -------
    !      This subroutine is to calculate ternary ion-mediated nucleation (TIMN) 
    !      from lookup tables using multiple-variable interpolation scheme. 
    !
    !      The present lookup table should cover almost all the possible conditions in 
    !      the troposphere relevant to atmospheric nucleation. The range and resolution 
    !      in each parameter space can be extended in the future if needed.
    !
    !      interface
    !      ---------
    !
    !      INPUT (valid value range):
    !      X0 = [H2SO4] in #/cm3  (5E5-5E9)
    !      Y0 = RH in % (0.5-99.5)
    !      Z0 = T (in K) (190-304)
    !      U0 = Q = ionization rate (ion-pairs cm-3 s-1) (2-100)
    !      V0 = S = surface area (um2/cm3) (20-200)
    !      W0 = B = NH3 in #/cm3  (1E5-1E12)
    !
    !      OUTPUT:
    !      XJTIM: Ternary Ion-Mediated Nucleation rate (# cm-3 s-1)
    !      The nucleation rates were calculated at clusters containing ~10
    !      H2SO4 molecules (and associated NH3/H2O molecules, mass diameter ~1.5 nm)
    !
    !
    !      method
    !      ------
    !      This subroutine is to calculate ternary ion-mediated nucleation (TIMN) 
    !      from lookup tables using multiple-variable interpolation scheme. 
    !
    !      reference
    !      ---------
    !      1. Yu, F., Nadykto, A. B., Herb, J., Luo, G., Nazarenko, K. M.,
    !         Uvarova, L. A., H2SO4-H2O-NH3 ternary ion-mediated nucleation (TIMN):
    !         Kinetic-based model and comparison with CLOUD measurements, Atmos.
    !         Chem. Phys.,18, 17451-17474,
    !         https://doi.org/10.5194/acp-18-17451-2018, 2018. 
    !      2. Yu, F., Nadykto, A. B., Luo, G., and Herb, J., 
    !         H2SO4-H2O binary and H2SO4-H2O-NH3 ternary homogeneous and 
    !         ion-mediated nucleation: look-up table version 1.0 for 3-D modeling 
    !         application,
    !         Geosci. Model Dev., 13, 2663-2670, doi:10.5194/gmd-13-2663-2020, 2020. 
    !
    !
    !      modifications
    !      -------------
    !      MSK CHANGED DURING IMPLEMENTATION 22.08.2020
    !      Include file of XJTIMN needs much memory
    !      XJTIMN 6-DIM      : XJTIMN(MC,MRH,MT,MQ,MS,MB)
    !      XJTIMN = reshape((/ values... /), shape(XJTIMN))
    !      We reduce size by allowing only 1 Q and 1 S value
    !      Q=2.000E+00 => Q(1); S=2.000E+01 => S(1)
    !      => TIMN(MC,MRH,MT,1,1,MB)
    !      SIZE OF XJTIMN( MC,MRH,MT,1,1,MB) = 
    !              39*26*33 x 32 = 33462 x 32
    !      and XJTIMN 4-DIM : XJTIMN(MC,MRH,MT,MB)   
    !      The include file can be split in more subunits if needed
    !      MSK CHANGED 20.12.2020
    !      first read include file into XJTIMN(MC,MRH*MT*MB)
    !      then translate to XJTIMN2(MC,MRH,MT,MB)
    !
    ! NOTE:
    !     BIMN includes BHN, THN includes BHN, and TIMN includes both BIMN and THN
    !
    !     Nucleation rates are calculated at 1.7 nm mobility
    !     diameter (corresponding to mass diameter of 1.5 nm;
    !     Yu et al., 2018).
    !
    !------------------------------------------------------------------

    implicit none

    ! input
    real( dp),intent(in)       :: X0,Y0,Z0,W0

    ! output
    real( dp),intent(out)      :: XJTIM
    real( dp),intent(out)      :: natot

        ! Parameters
        integer, parameter :: MC  = 32
        integer, parameter :: MT  = 39
        integer, parameter :: MRH = 26
        integer, parameter :: MS  = 2
        integer, parameter :: MB  = 33
        integer, parameter :: MQ  = 8

        ! Lookup data
        ! check to leave Q and S as constant values
        real( dp),dimension(MC)              :: C
        real( dp),dimension(MRH)             :: RH 
        real( dp),dimension(MT)              :: T
        real( dp),dimension(MQ)              :: Q
        real( dp),dimension(MS)              :: S
        real( dp),dimension(MB)              :: B
!MSK 30.12.2020: first read into array with size MC x MRH*MT*MB
        !!!real( dp),dimension(MC,MRH,MT,MB)    :: XJTIMN
        real( dp),dimension(MC,MRH*MT*MB)    :: XJTIMN
        real( dp),dimension(MC,MRH,MT,MB)    :: XJTIMN2

        ! dummy array for xjtimn data 
        real( dp),dimension(MC*MRH*MT*MB)    :: dumxjtimn

        ! local variables

        real( dp) :: X,Y,Z,U,V,W
        real( dp) :: C11,B11,Q11,S11
        real( dp) :: VOL,FRACT
        real( dp) :: dx1,dx2,dy1,dy2,dz1,dz2,du1,du2,dv1,dv2,dw1,dw2
        real( dp) :: dx,dy,dz,du,dv,dw
        real( dp) :: XDH,XDT,X1,X2,Y1,Y2,XJ1,XJ2,YJ

        integer   :: IC1,IC2,JRH1,JRH2,KT1,KT2,IQ1,IQ2,IS1,IS2,IB1,IB2
        integer   :: IC, JRH, KT, IQ, IS,IB
        integer   :: IX

       ! The present lookup table should cover almost all the possible conditions in 
       ! the troposphere relevant to atmospheric nucleation. The range and resolution 
       ! in each parameter space can be extended in the future if needed.

       ! Avoid the input values to be changed due to out of the range reset
       ! MSK NOTE: convert RH and c(H2SO4) and c(NH3)
       ! MSK NOTE: check default values of ion. rate and surface area
       X = X0*1.E-6_dp    ! [H2SO4] in #/cm3
       Y = Y0*100._dp     ! RH in %
       Z = Z0             ! T in K
       U = 2.00E+00_dp    ! ionization rate (ion-pairs cm-3 s-1) (2-100)
       V = 2.00E+01_dp    ! surface area (um2/cm3) (20-200)
       W = W0*1.E-6       ! [NH3] in #/cm3


       !H2SO4 [molecules/cm^3]: C(MC)
       C(:)  = (/  &
         5.000E+05, 6.290E+05, 7.920E+05, 9.980E+05, 1.260E+06, 1.580E+06, 1.990E+06,  &
         2.510E+06, 3.150E+06, 3.970E+06, 5.000E+06, 6.290E+06, 7.920E+06, 9.980E+06,  &
         1.260E+07, 1.580E+07, 1.990E+07, 2.510E+07, 3.150E+07, 3.970E+07, 5.000E+07,  &
         6.290E+07, 7.920E+07, 9.980E+07, 1.260E+08, 1.580E+08, 1.990E+08, 2.510E+08,  &
         3.150E+08, 3.970E+08, 5.000E+08, 5.000E+09  /)        
       !RH  [%]: RH(MRH)gamma,m,dx,d1,
       RH(:) = (/  &
         5.000E-01, 4.000E+00, 8.000E+00, 1.200E+01, 1.600E+01, 2.000E+01, 2.400E+01,  &
         2.800E+01, 3.200E+01, 3.600E+01, 4.000E+01, 4.400E+01, 4.800E+01, 5.200E+01,  &
         5.600E+01, 6.000E+01, 6.400E+01, 6.800E+01, 7.200E+01, 7.600E+01, 8.000E+01,  &
         8.400E+01, 8.800E+01, 9.200E+01, 9.600E+01, 9.950E+01  /)
       !T   [K]: T(MT)
       T(:)  = (/  &
         1.900E+02, 1.930E+02, 1.960E+02, 1.990E+02, 2.020E+02, 2.050E+02, 2.080E+02,  &
         2.110E+02, 2.140E+02, 2.170E+02, 2.200E+02, 2.230E+02, 2.260E+02, 2.290E+02,  &
         2.320E+02, 2.350E+02, 2.380E+02, 2.410E+02, 2.440E+02, 2.470E+02, 2.500E+02,  &
         2.530E+02, 2.560E+02, 2.590E+02, 2.620E+02, 2.650E+02, 2.680E+02, 2.710E+02,  &
         2.740E+02, 2.770E+02, 2.800E+02, 2.830E+02, 2.860E+02, 2.890E+02, 2.920E+02,  &
         2.950E+02, 2.980E+02, 3.010E+02, 3.040E+02  /)
       !Q   [ion pairs/cm^3/s]: Q(MQ)
       Q(:)  = (/  &
         2.000E+00, 3.000E+00, 4.500E+00, 6.750E+00, 1.010E+01, 1.520E+01, 2.280E+01,  &
         1.000E+02  /)
       !S    [um/m^2]: S(MS)
       S(:)   = (/  &
         2.000E+01, 2.000E+02  /)
       !NH3 [molecules/cm^3]: B(MB)
       B(:)   = (/  &
         1.000E+05, 1.000E+08, 1.260E+08, 1.580E+08, 2.000E+08, 2.510E+08, 3.160E+08,  &
         3.980E+08, 5.010E+08, 6.310E+08, 7.940E+08, 1.000E+09, 1.260E+09, 1.580E+09,  &
         2.000E+09, 2.510E+09, 3.160E+09, 3.980E+09, 5.010E+09, 6.310E+09, 7.940E+09,  &
         1.000E+10, 1.260E+10, 1.580E+10, 2.000E+10, 2.510E+10, 3.160E+10, 3.980E+10,  &
         5.010E+10, 6.310E+10, 7.940E+10, 1.000E+11, 1.000E+12  /)


       !MSK CHANGE 22.08.2020
       !  include file needs much memory
       !XJTIMN 6-DIM      : XJTIMN(MC,MRH,MT,MQ,MS,MB)
       !XJTIMN = reshape((/ values... /), shape(XJTIMN))
       !  We reduce size by allowing only 1 Q and 1 S value
       !  Q=2.000E+00 => Q(1); S=2.000E+01 => S(1)
       !  => TIMN(MC,MRH,MT,1,1,MB)
       !  and XJTIMN 4-DIM : XJTIMN(MC,MRH,MT,MB)   
       !  The include file can be split in more subunits if needed
       ! SIZE OF XJTIMN( MC,MRH,MT,1,1,MB) = 39*26*33 x 32 = 33462 x 32
       ! reshape dummy array dumxjtimn in include file

       include 'yu_timn_j6d.inc'

       XJTIMN = reshape(dumxjtimn, shape(XJTIMN) ) 

!CHECK INPUT MATRIX
!        do IC=1,MC
!          do IX=1,10   !MRH*MT*MB
!           write(6,*) 'XJTIMN',IC,IX,XJTIMN(IC,IX)
!          enddo
!        enddo
!        stop

!MSK 30.12.2020: write input matrix into XJTIMN2
        do IC = 1,MC
          IX=1
          do KT = 1,MT
           do JRH = 1,MRH
             do IB =1, MB

               XJTIMN2(IC,JRH,KT,IB) = XJTIMN(IC,IX)
               !write(6,*) 'IX',IX,IC,JRH,KT,IB
               IX=IX+1

             enddo
           enddo
          enddo
        enddo


! check input values
       !write(6,*) 'ch2so4',X,'rh',Y,'T',Z,'cnh3',W

       ! Use the formula to calculate C and Q to get values with more digits, otherwise
       ! may cause problem when input C and Q are very clsoe to C(IC),Q(IQ) as
       ! IC and IQ are decided with formula

       ! MSK NOTE: the program will only stop if the look-up tables
       !           have been modified

       IF(C(1).NE.5.0E5.or.C(MC).NE.5.0E9) THEN
          WRITE(6,*)"STOP2: need to check JTIMN look-up table inputs"                  
          STOP
       ENDIF
       DO IC = 2, MC-1
          C11 = C(IC)                                                          
          C(IC) = C(IC-1)*10.**(0.1)
          IF(abs(1.-C11/C(IC)).GT.0.02) THEN                                  
           WRITE(6,*)"STOP3: need to check JTIMN look-up table inputs"                  
           STOP
          ENDIF                                                                
       ENDDO

       IF(B(1).NE.1.0E5.or.B(2).NE.1.0E8.or.B(MB).NE.1.0E12) THEN
          WRITE(6,*)"STOP4: need to check JTIMN look-up table inputs"                  
          STOP
       ENDIF
       DO IB = 3, MB-1
          B11 = B(IB)
          B(IB) = B(IB-1)*10.**(0.1)
          IF(abs(1.-B11/B(IB)).GT.0.02) THEN
           WRITE(6,*)"STOP5: need to check JTIMN look-up table inputs"                  
           STOP
          ENDIF
       ENDDO

       IF(Q(1).NE.2.0.or.Q(MQ).NE.100.0) THEN
          WRITE(6,*)"STOP6: need to check JTIMN look-up table inputs"                  
          STOP
       ENDIF
       DO IQ = 2, MQ-1
          Q11 = Q(IQ)                                                          
          Q(IQ) = Q(IQ-1)*1.5
          IF(abs(1.-Q11/Q(IQ)).GT.0.02) THEN
           WRITE(6,*)"STOP7: need to check JTIMN look-up table inputs"                  
           STOP
          ENDIF
       ENDDO

       ! If the inputed values are out of the lookup table valid ranges, set them to 
       ! boundary values for now for all variables except surface area. Care should be 
       ! taken if your inputted values are frequently out of the specified ranges.
 
       IF(X.LT.C(1)) THEN
           WRITE(6,10) X, C(1), C(1)
           X = C(1)
       ELSEIF(X.GT.C(MC)) THEN
           WRITE(6,11) X, C(MC), C(MC)
           X =C(MC)
       ENDIF

       IF(Y.LT.RH(1)) THEN
           WRITE(6,12) Y, RH(1), RH(1)
           Y =RH(1) 
       ELSEIF(Y.GT.RH(MRH)) THEN
           WRITE(6,13) Y, RH(MRH), RH(MRH)
           Y =RH(MRH)
       ENDIF

       IF(Z.LT.T(1)) THEN
           WRITE(6,14) Z, T(1), T(1)
           Z =T(1)
       ELSEIF(Z.GT.T(MT)) THEN
           WRITE(6,15) Z, T(MT), T(MT)
           Z =T(MT)
       ENDIF

       !MSK 22.08.2020: constant Q
       !IF(U.LT.Q(1)) THEN
       !    WRITE(6,16) U, Q(1), Q(1)
       !    U =Q(1)
       !ELSEIF(U.GT.Q(MQ)) THEN
       !    WRITE(6,17) U, Q(MQ), Q(MQ)
       !    U =Q(MQ)
       !ENDIF
!
! Allow extropolation for S
!        IF(V.LT.S(1)) THEN
!!           WRITE(6,18) V, S(1), S(1)
!           V =S(1)
!        ELSEIF(V.GT.S(MS)) THEN
!!           WRITE(6,19) V, S(MS), S(MS)
!           V =S(MS)
!        ENDIF

       IF(W.LT.B(1)) THEN
           WRITE(6,20) W, B(1), B(1)
           W =B(1)
       ELSEIF(W.GT.B(MB)) THEN
           WRITE(6,21) W, B(MB), B(MB)
           W =B(MB)
       ENDIF

       ! end of range check of input values

       IC1 =MAX0(INT(1.+10.*LOG10(X/5.E5)),1)
       IC2 = MIN0(IC1 + 1,MC)
       IF(IC2.EQ.MC) IC1=MC-1
        
       XDH = 4.
       IF(Y.LT.RH(2)) THEN
         JRH1 = 1.
       ELSE
         JRH1 = MAX0(INT((Y-RH(2))/XDH+2.),2)
       ENDIF
       JRH2 = MIN0(JRH1 + 1,MRH)
       IF(JRH2.EQ.MRH) JRH1=MRH-1

       XDT = 3.0
       KT1 = MAX0(INT((Z-190.0)/XDT)+1,1)
       KT2 = MIN0(KT1 + 1,MT)
       IF(KT2.EQ.MT) KT1=MT-1

       IF(U.LT.Q(1)) THEN
         IQ1 =1.
       ELSE
         IQ1 = MAX0(INT(1.+LOG10(U/Q(1))/LOG10(1.5)),1)
       ENDIF
       IQ2 = MIN0(IQ1 + 1,MQ)
       IF(IQ2.EQ.MQ) IQ1=MQ-1


       IF(W.LT.B(2)) THEN
         IB1 =1.
       ELSE
         IB1 = MAX0(INT(2.+10.*LOG10(W/B(2))),2)
       ENDIF
       IB2 = MIN0(IB1 + 1,MB)
       IF(IB2.EQ.MB) IB1=MB-1

       ! logJ log[H2SO4] interpolation
       dx1 = LOG10(X/C(IC1))
       dx2 = LOG10(C(IC2)/X)

       dy1 = LOG10(Y/RH(JRH1))
       dy2 = LOG10(RH(JRH2)/Y)

       dz1 = Z-T(KT1)
       dz2 = T(KT2)-Z

       !MSK 22.08.2020: constant Q, du1 and du2 not used
       !du1 = U - Q(IQ1)
       !du2 = Q(IQ2) - U

       ! logJ log[H2SO4] interpolation
       dw1 = LOG10(W/B(IB1))
       dw2 = LOG10(B(IB2)/W)


     ! write(6,*) 'lims',IC1,IC2,JRH1,JRH2,KT1,KT2,IQ1,IQ2


       !JTIMN
        XJ1 = 0.  
        XJ2 = 0.

       !MSK 22.08.2020: constant Q, du1 and du2 not used
        !VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)*(dw1+dw2)
        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(dw1+dw2)

     ! write(6,*) 'VOL',dx1,dx2,dy1,dy2,dz1,dz2,dw1,dw2,VOL
     ! write (6,*) 'XJTIMN',XJTIMN(1:MC,1,1,1)


        DO KT = KT1,KT2
          IF(KT.EQ.KT1) THEN
            dz = dz2
	  ELSE
            dz = dz1
          ENDIF
      	  DO JRH = JRH1,JRH2
            IF(JRH.EQ.JRH1) THEN
              dy = dy2
	    ELSE
              dy = dy1
            ENDIF
            DO IC = IC1,IC2
              IF(IC.EQ.IC1) THEN
                dx = dx2
	      ELSE
                dx = dx1
              ENDIF
 	      DO IB =IB1, IB2
                IF(IB.EQ.IB1) THEN
                  dw = dw2
	        ELSE
                  dw = dw1
                ENDIF

                !MSK 22.08.2020: constant Q
	        !DO IQ =IQ1, IQ2

                  !IF(IQ.EQ.IQ1) THEN
                  !  du = du2
	          !ELSE
                  !  du = du1
                  !ENDIF

                  !FRACT = dx*dy*dz*du*dw/VOL
                  FRACT = dx*dy*dz*dw/VOL

                  ! nucleation rate
                  !MSK 22.08.2020: constant S [=S(1)]
                  !XJ1 = XJ1 + FRACT*XJTIMN(IC,JRH,KT,IQ,1,IB)
                  !XJ2 = XJ2 + FRACT*XJTIMN(IC,JRH,KT,IQ,2,IB)
                  XJ1 = XJ1 + FRACT*XJTIMN2(IC,JRH,KT,IB)
	        !ENDDO

	      ENDDO
            ENDDO
	  ENDDO
	ENDDO

     !  write(6,*) 'XJ1 log10J 1/(cm3*s)', XJ1

       ! Log10J -->J
       XJ1 = 10.**XJ1
       !XJ2 = 10.**XJ2

       !MSK 22.08.2020: constant S
       !MSK: the interpolation according to inputed S
       !     can be commented if using only S(1)
       ! Interpolate to get J at inputed S
       !X1 = S(1)
       !Y1 = XJ1
       !X2 = S(2)
       !Y2 = XJ2
       !IF(Y1.GT.Y2) THEN   ! Extrapolation for S
       !  YJ = Y1*(Y2/Y1)**((X1-V)/(X1-X2))
       !ELSE  
       !  YJ=Y2
       !ENDIF
       !XJTIM = YJ

       XJTIM = XJ1

       ! Nucleation rate in m^-3s^-1
       XJTIM = XJTIM * 1.e6_dp

    !  write(6,*) 'J (1/(m3*s)', XJTIM

       !     Nucleation rates are calculated at 1.7 nm mobility
       !     diameter (corresponding to mass diameter of 1.5 nm;
       !     Yu et al., 2018).

       ! assume 3 molecules in critical cluster
       natot=3._dp


     ! Print warnings for input data that is out of the range
     ! of the parameterization 
 10     FORMAT("yu_timn_nucleation WARNING: INPUTED [H2SO4]=",ES9.2,"<",ES9.2, &
               " set it to ",ES9.2)
 11     FORMAT("yu_timn_nucleation WARNING: INPUTED [H2SO4]=",ES9.2,">",ES9.2, &
               " set it to ",ES9.2)
 12     FORMAT("yu_timn_nucleation WARNING: INPUTED RH =",F5.1,"% <",F5.1,     &
               "% set it to ",F5.1,"%")
 13     FORMAT("yu_timn_nucleation WARNING: INPUTED RH =",F5.1,"% >",F5.1,     &
               "% set it to ",F5.1,"%")
 14     FORMAT("yu_timn_nucleation WARNING: INPUTED T =",F6.1,"K <",F6.1,      &
               "K set it to ",F6.1,"K")
 15     FORMAT("yu_timn_nucleation WARNING: INPUTED T =",F6.1,"K >",F6.1,      &
               "K set it to ",F6.1,"K")
 16     FORMAT("yu_timn_nucleation WARNING: INPUTED Q =",F6.1," <",F6.1,       &
               " ion-pair/cm3s set it to ",F6.1)
 17     FORMAT("yu_timn_nucleation WARNING: INPUTED Q =",F6.1," >",F6.1,       &
               " ion-pair/cm3s set it to ",F6.1)
 18     FORMAT("yu_timn_nucleation WARNING: INPUTED S =",F6.1," <",F6.1,       &
               " um2/cm3 set it to ",F6.1)
 19     FORMAT("yu_timn_nucleation WARNING: INPUTED S =",F6.1," >",F6.1,       &
               " um2/cm3 set it to ",F6.1)
 20     FORMAT("yu_timn_nucleation WARNING: INPUTED [NH3]=",ES9.2,"<",ES9.2,   &
               " set it to ",ES9.2)
 21     FORMAT("yu_timn_nucleation WARNING: INPUTED [NH3]=",ES9.2,">",ES9.2,   &
               " set it to ",ES9.2)


   end subroutine yu_timn_nucleation




     subroutine ACDC_nucleation(X0,Y0,Z0,U0,V0,XJ0,XI0)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Shuai Jiang
    !
    !      contact
    !      -------
    !      Dr. Shuai Jiang
    !      Associate Professor
    !      School of Information Science and Technology
    !      University of Science and Technology of China (USTC)
    !      Hefei, Anhui, 230026
    !      P.R. China
    !      Email: shuaijiang@ustc.edu.cn
    !
    !      purpose
    !      -------
    !      This subroutine is to calculate rates from lookup tables
    !      using multiple-variable interpolation scheme based on
    !      ACDC derived rates
    !
    !      interface
    !      ---------
    !
    ! INPUT PARAMETERS AND VALID RANGES:
    ! X0 = Sulfuric acid vapor concentration (H2SO4]) in #/cm3  (1E4-3.16E9)
    ! Y0 = Ammonia concentration (NH3]) in #/cm3  (1E6-1E11)
    ! Z0 = Temperature (T) in K (180-320)
    ! U0 = Relative humidity (RH) in % (0-100)
    ! V0 = Coagulation Sink in #/s (1E-5-1E-1)
    !
    !  Parameters
    !  (1 ) MA  : NUMBER OF POINTS IN H2SO4 CONCENTRATION DIMENSION
    !  (2 ) MB  : NUMBER OF POINTS IN NH3 CONCENTRATION DIMENSION
    !  (3 ) MT   : NUMBER OF POINTS IN TEMPERATURE DIMENSION
    !  (4 ) MRH  : NUMBER OF POINTS IN RELATIVE HUMIDITY DIMENSION
    !  (5 ) MCS   : NUMBER OF POINTS IN COAGULATION SINK DIMENSION
    !  Arrays
    !  (6 ) A   : VALUES AT POINTS IN H2SO4 CONCENTRATION DIMENSION
    !  (7 ) B   : VALUES AT POINTS IN NH3 CONCENTRATION DIMENSION  
    !  (8 ) T   : VALUES AT POINTS IN TEMPERATURE DIMENSION
    !  (9 ) RH  : VALUES AT POINTS IN RELATIVE HUMIDITY DIMENSION
    !  (10) CS   : VALUES AT POINTS IN SURFACE AREA DIMENSION
    !
    ! OUTPUT:
    ! XJ0: Nucleation rate in  #/cm3s
    !
    !      method
    !      ------
    !      The multiple-variable interpolation scheme is originally 
    !      written by Prof. Fangqun Yu in State University of New York at Albany
    !      and here revised by Shuai Jiang (shuaijiang@ustc.edu.cn) 
    !      in University of Science and Technology of China (USTC)
    !      The ACDC derived rates now is obtained from Hanna Vehkamaki Group 
    !      for H2SO4-NH3-H2O system
    !
    !      reference
    !      ---------
    !      Henschel, H., Kurten, T., and Vehkamaki, H.,
    !      Computational Study on the effect of hydration on new 
    !      particle formation in the sulfuric acid/ammonia and 
    !      sulfuric acid/dimethylamine systems.
    !      J. Phys. Chem. A, 120, 1886-1896, 
    !      doi:10.1021/acs.jpca.5b11366, 2016.
    !
    !      Baranizadeh, E., Murphy, B.N., Julin, J., Falahat, S.
    !      Reddington, C.L., Arola, A., Ahlm, L., Mikkonen, S.,
    !      Fountoukis, C., Patoulias, D., Minikin, A., Hamburger, T.,
    !      Laaksonen, A., Pandis, S.N., Vehkamaki, H., Lehtinen, K.E.J., 
    !      Riipinen, I.,
    !      Implementation of state-of-the-art ternary new-particle 
    !      formation scheme to the regional chemical transport model 
    !      PMCAMx-UF in Europe.
    !      Geosci. Model Dev., 9, 2741–2754,
    !      doi:10.5194/gmd-9-2741-2016, 2016.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
   
    implicit none

    integer, parameter :: MA=12, MB=11, MT=15, MRH=21, MCS=9

    ! input
    REAL( dp),intent(in)       :: X0,Y0,Z0,U0,V0
    
    ! output
    REAL( dp),intent(out)      :: XJ0,XI0

    ! lookup data
    REAL( dp),dimension(MA)              :: A
    REAL( dp),dimension(MB)              :: B 
    REAL( dp),dimension(MT)              :: T
    REAL( dp),dimension(MRH)             :: RH
    REAL( dp),dimension(MCS)             :: CS
    REAL( dp),dimension(MB,MA,MCS,MRH,MT) :: XJACDC
!    REAL( dp),dimension(MA,MB,MT,MRH,MCS) :: XJACDC
!    REAL( dp),dimension(MB,MA,MCS,MRH,MT) :: XJACDC
    REAL( dp),dimension(MB*MA*MCS*MRH*MT) :: dumxjacdc

        ! local variables
        real( dp) :: X,Y,Z,U,V 
        real( dp) :: dx1,dx2,dy1,dy2,dz1,dz2 ,du1,du2,dv1,dv2
        real( dp) :: dx,dy,dz ,du,dv
        real( dp) :: RATIO, VOL, VOL3, FRACT
        real( dp) :: A11,B11,RH11,CS11,T11
        integer   :: IA1,IA2,JB1,JB2,KT1,KT2,IRH1,IRH2,ICS1,ICS2
        integer   :: IA,IB,IT,JRH,KCS
         
       X = X0*1.E-6_dp   ! [H2SO4] in #/cm3
       Y = Y0*1.E-6_dp    ! [NH3] in #/cm3
       Z = Z0
       U = U0
       V = V0
       
       !H2SO4 [molecules/cm^3]: C(MA)
       A(:)  = (/  &
         1.00e+04,  3.16e+04,  1.00e+05,  &
         3.16e+05,  1.00e+06,  3.16e+06,  & 
         1.00e+07,  3.16e+07,  1.00e+08,  &
         3.16e+08,  1.00e+09,  3.16e+09  /)        
       !NH3 [molecules/cm^3]: C(MB)
       B(:) = (/  &
         1.00e+06,  3.16e+06,  1.00e+07,  &
         3.16e+07,  1.00e+08,  3.16e+08,  &
         1.00e+09,  3.16e+09,  1.00e+10,  &
         3.16e+10,  1.00e+11  /)
       !T   [K]: T(MT)
       T(:)  = (/  &
         1.80e+02,  1.90e+02,  2.00e+02,  & 
         2.10e+02,  2.20e+02,  2.30e+02,  &
         2.40e+02,  2.50e+02,  2.60e+02,  &
         2.70e+02,  2.80e+02,  2.90e+02,  &
         3.00e+02,  3.10e+02,  3.20e+02  /)
       !Q   [ion pairs/cm^3/s]: Q(MQ)
       RH(:)  = (/  &
         0.00e+00,  5.00e-02,  1.00e-01,  &
         1.50e-01,  2.00e-01,  2.50e-01,  &
         3.00e-01,  3.50e-01,  4.00e-01,  &
         4.50e-01,  5.00e-01,  5.50e-01,  &
         6.00e-01,  6.50e-01,  7.00e-01,  &
         7.50e-01,  8.00e-01,  8.50e-01,  &
         9.00e-01,  9.50e-01,  1.00e+00  /)
       !S    [um/m^2]: S(MS)
       CS(:)   = (/  &
         1.00e-05,  3.16e-05,  1.00e-04,  &
         3.16e-04,  1.00e-03,  3.16e-03,  &
         1.00e-02,  3.16e-02,  1.00e-01  /)


       !XJACDC 5-DIM      : XJACDC(MA,MB,MT,MRH,MCS)
       !XJACDC = reshape((/ values... /), shape(XJIMN))
  
       include 'Lookup_table_J_5D-mk.inc'
       XJACDC = reshape(dumxjacdc, shape(XJACDC) ) 
       
!       write(6,*) "XJACDC", XJACDC(1,2,1,1,1)
       
       XJACDC=LOG10(XJACDC)
       


!      DO IT =1, MT
!         DO JRH = 1,MRH
!          DO KCS = 1,MCS
!           DO IA =1, MA
!            READ(31,201)(XJACDC(IA,IB,IT,JRH,KCS),IB = 1,MB)
!            DO IB=1, MB
!             XJACDC(IA,IB,IT,JRH,KCS)=LOG10(XJACDC(IA,IB,IT,JRH,KCS))
!            ENDDO
!           ENDDO
!          ENDDO
!         ENDDO
!      ENDDO

! Use the formula to get values with more 
! digits, otherwise may cause problem in interpretation
        DO IA = 2, MA
           A11 = A(IA)                                                          
           RATIO = 10.**(0.5)
           A(IA) = A(IA-1)*RATIO
! Double check to make sure that you get the right table.
           IF(abs(1.-A11/A(IA)).GT.0.02) THEN                                  
              write(6,*)"need check JACDC look-up table inputs"                  
              stop                                                              
           ENDIF                                                                
        ENDDO
        
        DO IB = 2, MB
           B11 = B(IB)                                                          
           RATIO = 10.**(0.5)
           B(IB) = B(IB-1)*RATIO
! Double check to make sure that you get the right table.
           IF(abs(1.-B11/B(IB)).GT.0.02) THEN                                  
              write(6,*)"need check JACDC look-up table inputs"                  
              stop                                                              
           ENDIF                                                                
        ENDDO
        
        DO KCS = 2, MCS
           CS11 = CS(KCS)                                                          
           RATIO = 10.**(0.5)
           CS(KCS) = CS(KCS-1)*RATIO
! Double check to make sure that you get the right table.
           IF(abs(1.-CS11/CS(KCS)).GT.0.02) THEN                                  
              write(6,*)"need check JACDC look-up table inputs"                  
              stop                                                              
           ENDIF                                                                
        ENDDO

        DO IT = 2, MT
           T11 = T(IT)                                                          
           RATIO = 10.
           T(IT) = T(IT-1)+RATIO
! Double check to make sure that you get the right table.
           IF(abs(1.-T11/T(IT)).GT.0.02) THEN                                  
              write(6,*)"need check JACDC look-up table inputs"                  
              stop                                                              
           ENDIF                                                                
        ENDDO       

        DO JRH = 2, MRH
           RH11 = RH(JRH)                                                          
           RATIO = 5.00e-02
           RH(JRH) = RH(JRH-1)+RATIO
! Double check to make sure that you get the right table.
           IF(abs(1.-RH11/RH(JRH)).GT.0.02) THEN                                  
              write(6,*)"need check JACDC look-up table inputs"                  
              stop                                                              
           ENDIF                                                                
        ENDDO
        
! Check of the input values 
        IF(X.LT.A(1)) THEN
           WRITE(6,10) X, A(1)
           XJ0 = 1.E-20
        ELSEIF(X.GT.A(MA)) THEN
           WRITE(6,11) X, A(MA), A(MA)
           X = A(MA)
        ENDIF

        IF(Y.LT.B(1)) THEN
           WRITE(6,12) Y, B(1)
           XJ0 = 1.E-20
           RETURN
        ELSEIF(Y.GT.B(MB)) THEN
           WRITE(6,13) Y, B(MB), B(MB)
           Y = B(MB)
        ENDIF

        IF(Z.LT.T(1)) THEN
           WRITE(6,14) Z, T(1), T(1)
           Z = T(1)
        ELSEIF(Z.GT.T(MT)) THEN
           WRITE(6,15) Z, T(MT)
           XJ0 = 1.E-20
           RETURN
        ENDIF

        IF(U.LT.RH(1)) THEN
           WRITE(6,16) U, RH(1)
           XJ0 = 1.E-20
           RETURN
        ELSEIF(U.GT.RH(MRH)) THEN
           WRITE(6,17) U, RH(MRH), RH(MRH)
           U =RH(MRH)
        ENDIF

        IF(V.LT.CS(1)) THEN
           WRITE(6,18) V, CS(1), CS(1)
           V = CS(1)
        ELSEIF(V.GT.CS(MCS)) THEN
           WRITE(6,19) V, CS(MCS)
           XJ0 = 1.E-20
           RETURN
        ENDIF        

 10     FORMAT("ACDC WARNING: INPUTED [A]=",ES9.2,"<",ES9.2, & 
          ", set JACDC to 1.E-20 cm-3s-1")
 11     FORMAT("ACDC WARNING: INPUTED [A]=",ES9.2,">",ES9.2, & 
          " set it to ",ES9.2)
 12     FORMAT("ACDC WARNING: INPUTED [B] =",F5.1,"<",F5.1, & 
          "set JACDC to 1.E-20 cm-3s-1")
 13     FORMAT("ACDC WARNING: INPUTED [B] =",F5.1,">",F5.1, &
          "set it to ",F5.1)
 14     FORMAT("ACDC WARNING: INPUTED T =",F6.1,"K <",F6.1, &
          "K set it to ",F6.1,"K")
 15     FORMAT("ACDC WARNING: INPUTED T =",F6.1,"K >",F6.1, &
          "K, set JACDC to 1.E-20 cm-3s-1")
 16     FORMAT("ACDC WARNING: INPUTED RH =",F6.1," <",F6.1, &
          "set JACDC to 1.E-20 cm-3s-1")
 17     FORMAT("ACDC WARNING: INPUTED RH =",F6.1," >",F6.1, &
          "set it to ",F6.1)
 18     FORMAT("ACDC WARNING: INPUTED CS =",F6.1,"s-1 <",F6.1, &
          "s-1 set it to ",F6.1)
 19     FORMAT("ACDC WARNING: INPUTED CS =",F6.1,"s-1 >",F6.1, &
          "s-1, set JACDC to 1.E-20 cm-3s-1")        
        
        IA1 = MAX0(INT(1.+2.*LOG10(X/A(1))),1)
        IA2 = MIN0(IA1 + 1,MA)
        IF(IA2.EQ.MA) IA1=MA-1
        
        JB1 = MAX0(INT(1.+2.*LOG10(Y/B(1))),1)
        JB2 = MIN0(JB1 + 1,MB)
        IF(JB2.EQ.MB) JB1=MB-1
        
        KT1 = MAX0(INT(Z/10.-17.),1)
        KT2 = MIN0(KT1 + 1,MT)
        IF(KT2.EQ.MT) KT1=MT-1
        
        IRH1 = MAX0(INT(1.+20.*U),1)
        IRH2 = MIN0(IRH1 + 1,MRH)
        IF(IRH2.EQ.MRH) IRH1=MRH-1  
        
        ICS1 = MAX0(INT(1.+2.*LOG10(V/CS(1))),1)
        ICS2 = MIN0(ICS1 + 1,MCS)
        IF(ICS2.EQ.MCS) ICS1=MCS-1      

       ! logJ log[H2SO4] interpolation
          dx1 = LOG10(X/A(IA1))   ! logJ log[H2SO4] interpolation
          dx2 = LOG10(A(IA2)/X)
          dy1 = LOG10(Y/B(JB1))
          dy2 = LOG10(B(JB2)/Y)
          dz1 = Z-T(KT1)
          dz2 = T(KT2)-Z

          du1 = U - RH(IRH1)
          du2 = RH(IRH2) - U
          dv1 = LOG10(V/CS(ICS1))
          dv2 = LOG10(CS(ICS2)/V)

       ! initialize output
        XJ0 = 0.

        VOL = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)*(du1+du2)*(dv1+dv2)
       ! VOL3 = (dx1+dx2)*(dy1+dy2)*(dz1+dz2)

                  
        DO IT = KT1,KT2
          IF(IT.EQ.KT1) THEN
            dz = dz2
                ELSE
            dz = dz1
          ENDIF
                DO JRH =IRH1, IRH2
            IF(JRH.EQ.IRH1) THEN
              du = du2
                ELSE
              du = du1
            ENDIF
                DO KCS =ICS1, ICS2
              IF(KCS.EQ.ICS1) THEN
                dv = dv2
                 ELSE
                dv = dv1
              ENDIF
              DO IA = IA1,IA2
                IF(IA.EQ.IA1) THEN
                  dx = dx2
                    ELSE
                  dx = dx1
                ENDIF                  
                DO IB = JB1,JB2
                  IF(IB.EQ.JB1) THEN
                    dy = dy2
                      ELSE
                    dy = dy1
                  ENDIF                        
                  FRACT = dx*dy*dz*du*dv/VOL 
                  XJ0 = XJ0 + FRACT*XJACDC(IB,IA,KCS,JRH,IT)
!                  write(6,*) "XJ0,FRACT,XJACDC,10.**XJ0", IA,IB,IT,JRH,KCS,XJ0,FRACT,10.**XJACDC(IB,IA,KCS,JRH,IT),10.**XJ0
!                  stop
                ENDDO
              ENDDO
            ENDDO
                ENDDO
              ENDDO

        ! Log10J -->J
        XJ0 = 10.**XJ0
        IF(X.LT.A(1)) THEN
          XJ0 = 0.0
        ENDIF

        XJ0 = XJ0 * 1.e6_dp      !m^-3s^-1
        
        XI0=2._dp
        
!     write(21,*) X0,Y0,Z0,U0,V0,XJ0
!     write(6,*) X0,Y0,Z0,U0,V0,XJ0
201 FORMAT(100(1PE9.2))
     write(6,*) "[A] in ACDC is (cm-3)", X0/1.00e6_dp
     write(6,*) "[B] in ACDC is (cm-3)", Y0/1.00e6_dp
     write(6,*) "T in ACDC is (K)", Z0
     write(6,*) "RH in ACDC is (%)", U0
     write(6,*) "CS in ACDC is (s-1) ", V0
     write(6,*) "jnuc in ACDC is (cm-3 s-1) ", XJ0/1.00e6_dp

    

   end subroutine ACDC_nucleation



 
     subroutine ternary_fit(t,rh,nas,nbs,cair,J_log,nacid,namm,r)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Joonas Merikanto
    !
    !      contact
    !      -------
    !      Dr. Joonas Merikanto 
    !      Climate System Research Unit
    !      Finnish Meteorological Institute
    !      Erik Palmenin Aukio 1
    !      FI-00560
    !      Helsinki, Finland
    !      Phone: +358405668512
    !      Email: joonas.merikanto@fmi.fi
    !
    !
    !      purpose
    !      -------
    !      Fortran 90 subroutine that calculates the parameterized composition
    !      and nucleation rate of critical clusters in h2o-h2so4-nh3 vapor
    !
    !      WARNING: The fit should not be used outside its limits of validity
    !      (limits indicated below)
    !
    !      interface
    !      ---------
    !
    ! IN:
    ! T:     temperature (K), limits 235-295 K
    ! rh:    relative humidity as fraction (eg. 0.5=50%) limits 0.05-0.95
    ! c2:    sulfuric acid concentration (molecules/cm3) limits 5x104 - 109 molecules/cm3
    ! c3:    ammonia mixing ratio (ppt) limits 0.1 - 1000 ppt
    !
    ! OUT:
    ! J_log: logarithm of nucleation rate (1/(s cm3))
    ! ntot:  total number of molecules in the critical cluster
    ! nacid: number of sulfuric acid molecules in the critical cluster
    ! namm:  number of ammonia molecules in the critical cluster
    ! r:     radius of the critical cluster (nm)
    !
    !
    !      method
    !      ------
    !      calculates the parameterized composition
    !      and nucleation rate of critical clusters in h2o-h2so4-nh3 vapor
    !
    !      reference
    !      ---------
    !      Merikanto, J., I. Napari, H. Vehkamäki, T. Anttila, and M. Kulmala, 2007.
    !      New parameterization of sulfuric acid-ammonia-water ternary nucleation 
    !      rates at tropospheric conditions, J. Geophys. Res., 112, D15207, 
    !      doi:10.1029/2006JD007977.
    !
    !      Merikanto, J., Napari, I., Vehkamäki, H., Anttila, T., Kulmala, 2009.
    !      Correction to "New parameterization of sulfuric acid-ammonia-water 
    !      ternary nucleation rates at tropospheric conditions", 
    !      J. Geophys. Res., 114, D09206, doi:10.1029/2009JD012136.
    !
    !      modifications
    !      -------------
    !      Remarks from Matthias Karl (MSK) on 2020-07-09:
    !
    !      The Fortran subroutine was distributed together with the article of
    !      Merikanto et al. (2009) under the name "ternary_parameterization"
    !      and is only slightly modified for the use in MAFOR.
    !
    !      The Correction (Merikanto et al., 2009) states that all coefficients
    !      should be given with 16 significant digits. This was already the case
    !      in the original code. The term i15 should be multiplied with RH,
    !      which was already the case in the original code, but not in Merikanto
    !      et al. (2007).
    ! 
    !      Changed by Matthias Karl (MSK):
    !      2017-04-30 set lower limit of j_log is -70. to avoid numerical problems
    !      2017-09-06 corrected unit conversion for NH3
    !      2020-07-09 (t_onset.gt.t) is the correct conditional
    !
    !------------------------------------------------------------------

    implicit none

       REAL( dp),intent(in)          :: t,rh,nas,nbs,cair
       REAL( dp),intent(out)         :: J_log,nacid,namm,r
       REAL( dp)                     :: J,t_onset,c2,c3

       ! convert sulfuric acid concentration from molec/m3 to molec/cm3
       c2 = nas*1.e-6_dp

       ! convert ammonia conc from molec/m3 to ppt
       c3 = nbs*1.e-6_dp
!MSK 06.09.2017
! c3 is in molec/cm^3; cair is in molec/cm^3
!       c3 = c3/(cair*1.e-10)
!MSK 06.09.2017 correct conversion to ppt needs 10^-12
       c3 = c3/(cair*1.e-12)
      !write(6,*) 'NH3 ppt', c3, nbs

       t_onset=143.6002929064716 + 1.0178856665693992*rh +                &
        10.196398812974294*LOG(c2) -                                      &
        0.1849879416839113*LOG(c2)**2 - 17.161783213150173*LOG(c3) +      &
        (109.92469248546053*LOG(c3))/LOG(c2) +                            &
        0.7734119613144357*LOG(c2)*LOG(c3)-0.15576469879527022*LOG(c3)**2
      !write(6,*) 't_onset',t_onset,t,c2,c3

       !-----------------------------------------------------------------
       ! Comment from MSK 2020-07-09:
       ! (t_onset.gt.t) is the correct conditional for allowing ternary
       ! nucleation.
       ! In Merikanto et al. (2007) page 6 it is stated:
       ! "Thus, within the validity regions of the parameterization, 
       ! significant ternary nucleation does not
       ! occur at temperatures higher than 295 K.", i.e. only when the 
       ! sulfuric acid concentration is > 10^9 cm^-3.
       ! In Merikanto et al. (2007) page 4 it is stated:
       ! "[21] If t_onset exceeds the temperature of interest,
       ! the nucleation rate is less than 5x10^(-6) cm^(-3)s^(-1)
       ! and should be set to zero.", i.e. if t_onset is below the air
       ! temperature then there will be no ternary nucleation.
       ! Example:
       !    NH3 790.9 ppt, H2SO4 1.e7 cm-3 --> t_onset = 268.014091 K
       !    NH3 790.9 ppt, H2SO4 1.e8 cm-3 --> t_onset = 282.974765 K
       !    NH3 790.9 ppt, H2SO4 1.e9 cm-3 --> t_onset = 297.238735 K
       ! => only with H2SO4=1.e9 cm-3 ternary nucleation occurs at 290 K
       !----------------------------------------------------------------- 

       if(t_onset.gt.t) then


       J_log=-12.861848898625231 +                                        &
!i=1 and i=6
        4.905527742256349*c3 -                 &
        358.2337705052991*rh -                                            &
        0.05463019231872484*c3*t + 4.8630382337426985*rh*t +              &
        0.00020258394697064567*c3*t**2 - 0.02175548069741675*rh*t**2 -    &
        2.502406532869512e-7*c3*t**3 + 0.00003212869941055865*rh*t**3 -   &
!i=5
        4.39129415725234e6/LOG(c2)**2 +                                   &
        (56383.93843154586*t)/LOG(c2)**2 -                                &
        (239.835990963361*t**2)/LOG(c2)**2 +                              &
        (0.33765136625580167*t**3)/LOG(c2)**2 -                           &

!i=15 MSK: must be multiplied by RH (see Correction, Merikanto et al., 2009)
        (629.7882041830943*rh)/(c3**3*LOG(c2)) +                          &
        (7.772806552631709*rh*t)/(c3**3*LOG(c2)) -                        &
        (0.031974053936299256*rh*t**2)/(c3**3*LOG(c2)) +                  &
        (0.00004383764128775082*rh*t**3)/(c3**3*LOG(c2)) +                &

!i=3
        1200.472096232311*LOG(c2) - 17.37107890065621*t*LOG(c2) +         &
        0.08170681335921742*t**2*LOG(c2) -                                &
        0.00012534476159729881*t**3*LOG(c2) -                             &
!i=4
        14.833042158178936*LOG(c2)**2+0.2932631303555295*t*LOG(c2)**2 -   &
        0.0016497524241142845*t**2*LOG(c2)**2 +                           &
        2.844074805239367e-6*t**3*LOG(c2)**2 -                            &
!i=7 and i=10
        231375.56676032578*LOG(c3) -                                      &
        100.21645273730675*rh*LOG(c3) + 2919.2852552424706*t*LOG(c3) +    &
        0.977886555834732*rh*t*LOG(c3)-12.286497122264588*t**2*LOG(c3) -  &
        0.0030511783284506377*rh*t**2*LOG(c3) +                           &
        0.017249301826661612*t**3*LOG(c3) +                               &
        2.967320346100855e-6*rh*t**3*LOG(c3) +                            &
!i=12
        (2.360931724951942e6*LOG(c3))/LOG(c2) -                           &
        (29752.130254319443*t*LOG(c3))/LOG(c2) +                          &
        (125.04965118142027*t**2*LOG(c3))/LOG(c2) -                       &
        (0.1752996881934318*t**3*LOG(c3))/LOG(c2) +                       &
!i=11
        5599.912337254629*LOG(c2)*LOG(c3) -                               &
        70.70896612937771*t*LOG(c2)*LOG(c3) +                             &
        0.2978801613269466*t**2*LOG(c2)*LOG(c3) -                         &
        0.00041866525019504*t**3*LOG(c2)*LOG(c3) +                        &
!i=8
        75061.15281456841*LOG(c3)**2 -                                    &
        931.8802278173565*t*LOG(c3)**2 +                                  &
        3.863266220840964*t**2*LOG(c3)**2 -                               &
        0.005349472062284983*t**3*LOG(c3)**2 -                            &
!i=16
        (732006.8180571689*LOG(c3)**2)/LOG(c2) +                          &
        (9100.06398573816*t*LOG(c3)**2)/LOG(c2) -                         &
        (37.771091915932004*t**2*LOG(c3)**2)/LOG(c2) +                    &
        (0.05235455395566905*t**3*LOG(c3)**2)/LOG(c2) -                   &
!i=18
        1911.0303773001353*LOG(c2)*LOG(c3)**2 +                           &
        23.6903969622286*t*LOG(c2)*LOG(c3)**2 -                           &
        0.09807872005428583*t**2*LOG(c2)*LOG(c3)**2 +                     &
        0.00013564560238552576*t**3*LOG(c2)*LOG(c3)**2 -                  &
!i=9
        3180.5610833308*LOG(c3)**3 + 39.08268568672095*t*LOG(c3)**3 -     &
        0.16048521066690752*t**2*LOG(c3)**3 +                             &
        0.00022031380023793877*t**3*LOG(c3)**3 +                          &
!i=17
        (40751.075322248245*LOG(c3)**3)/LOG(c2) -                         &
        (501.66977622013934*t*LOG(c3)**3)/LOG(c2) +                       &
        (2.063469732254135*t**2*LOG(c3)**3)/LOG(c2) -                     &
        (0.002836873785758324*t**3*LOG(c3)**3)/LOG(c2) +                  &
!i=19
        2.792313345723013*LOG(c2)**2*LOG(c3)**3 -                         &
        0.03422552111802899*t*LOG(c2)**2*LOG(c3)**3 +                     &
        0.00014019195277521142*t**2*LOG(c2)**2*LOG(c3)**3 -               &
        1.9201227328396297e-7*t**3*LOG(c2)**2*LOG(c3)**3 -                &
!i=2
        980.923146020468*LOG(rh) + 10.054155220444462*t*LOG(rh) -         &
        0.03306644502023841*t**2*LOG(rh) +                                &
        0.000034274041225891804*t**3*LOG(rh) +                            &
!i=13
        (16597.75554295064*LOG(rh))/LOG(c2) -                             &
        (175.2365504237746*t*LOG(rh))/LOG(c2) +                           &
        (0.6033215603167458*t**2*LOG(rh))/LOG(c2) -                       &
        (0.0006731787599587544*t**3*LOG(rh))/LOG(c2) -                    &
!i=14
        89.38961120336789*LOG(c3)*LOG(rh) +                               &
        1.153344219304926*t*LOG(c3)*LOG(rh) -                             &
        0.004954549700267233*t**2*LOG(c3)*LOG(rh) +                       &
        7.096309866238719e-6*t**3*LOG(c3)*LOG(rh) +                       &
!i=20
        3.1712136610383244*LOG(c3)**3*LOG(rh) -                           &
        0.037822330602328806*t*LOG(c3)**3*LOG(rh) +                       &
        0.0001500555743561457*t**2*LOG(c3)**3*LOG(rh) -                   &
        1.9828365865570703e-7*t**3*LOG(c3)**3*LOG(rh)


        J=EXP(J_log)

       ! write(6,*) 'ternary0', J, J_log,c2,c3,t,rh

       ! Comment from MSK:
       ! ntot:  total number of molecules in the critical cluster
       !        is not used further and therefore commented here.
       ! ntot=57.40091052369212 - 0.2996341884645408*t +                   &
       !  0.0007395477768531926*t**2 -                                     &
       !  5.090604835032423*LOG(c2) + 0.011016634044531128*t*LOG(c2) +     &
       !  0.06750032251225707*LOG(c2)**2 - 0.8102831333223962*LOG(c3) +    &
       !  0.015905081275952426*t*LOG(c3) -                                 &
       !  0.2044174683159531*LOG(c2)*LOG(c3) +                             &
       !  0.08918159167625832*LOG(c3)**2 -                                 &
       !  0.0004969033586666147*t*LOG(c3)**2 +                             &
       !  0.005704394549007816*LOG(c3)**3 + 3.4098703903474368*LOG(J) -    &
       !  0.014916956508210809*t*LOG(J) +                                  &
       !  0.08459090011666293*LOG(c3)*LOG(J) -                             &
       !  0.00014800625143907616*t*LOG(c3)*LOG(J) +                        &
       !  0.00503804694656905*LOG(J)**2

        r=3.2888553966535506e-10 - 3.374171768439839e-12*t +              &
         1.8347359507774313e-14*t**2 + 2.5419844298881856e-12*LOG(c2) -   &
         9.498107643050827e-14*t*LOG(c2) +                                &
         7.446266520834559e-13*LOG(c2)**2 +                               &
         2.4303397746137294e-11*LOG(c3) +                                 &
         1.589324325956633e-14*t*LOG(c3) -                                &
         2.034596219775266e-12*LOG(c2)*LOG(c3) -                          &
         5.59303954457172e-13*LOG(c3)**2 -                                &
         4.889507104645867e-16*t*LOG(c3)**2 +                             &
         1.3847024107506764e-13*LOG(c3)**3 +                              &
         4.141077193427042e-15*LOG(J) - 2.6813110884009767e-14*t*LOG(J)+  &
         1.2879071621313094e-12*LOG(c3)*LOG(J) -                          &
         3.80352446061867e-15*t*LOG(c3)*LOG(J) -                          &
         1.8790172502456827e-14*LOG(J)**2

        nacid=-4.7154180661803595 + 0.13436423483953885*t -               &
         0.00047184686478816176*t**2 -                                    &
         2.564010713640308*LOG(c2) + 0.011353312899114723*t*LOG(c2) +     &
         0.0010801941974317014*LOG(c2)**2+0.5171368624197119*LOG(c3)-     &
         0.0027882479896204665*t*LOG(c3) +                                &
         0.8066971907026886*LOG(c3)**2 -                                  &
         0.0031849094214409335*t*LOG(c3)**2 -                             &
         0.09951184152927882*LOG(c3)**3 +                                 &
         0.00040072788891745513*t*LOG(c3)**3+1.3276469271073974*LOG(J)-   &
         0.006167654171986281*t*LOG(J) -                                  &
         0.11061390967822708*LOG(c3)*LOG(J) +                             &
         0.0004367575329273496*t*LOG(c3)*LOG(J) +                         &
         0.000916366357266258*LOG(J)**2

        namm=71.20073903979772 - 0.8409600103431923*t +                   &
         0.0024803006590334922*t**2 +                                     &
         2.7798606841602607*LOG(c2) - 0.01475023348171676*t*LOG(c2) +     &
         0.012264508212031405*LOG(c2)**2 - 2.009926050440182*LOG(c3) +    &
         0.008689123511431527*t*LOG(c3) -                                 &
         0.009141180198955415*LOG(c2)*LOG(c3) +                           &
         0.1374122553905617*LOG(c3)**2 -                                  &
         0.0006253227821679215*t*LOG(c3)**2 +                             &
         0.00009377332742098946*LOG(c3)**3 + 0.5202974341687757*LOG(J) -  &
         0.002419872323052805*t*LOG(J) +                                  &
         0.07916392322884074*LOG(c3)*LOG(J) -                             &
         0.0003021586030317366*t*LOG(c3)*LOG(J) +                         &
         0.0046977006608603395*LOG(J)**2


!MSK 30.04.2017: lower limit of j_log is -70. to avoid numerical problems
         j_log=max(-70._dp,j_log)

       else
       ! Nucleation rate less that 5E-6, setting j_log arbitrary small
!MSK 30.04.2017: but not smaller than 1.e-32
         j_log=-70._dp
       endif
      
! output is logarithm of nucleation rate (in cm^-3s^-1)

   end subroutine ternary_fit




     subroutine newbinapara(t,satrat,ch2so4,airn, jnuc_n,jnuc_i, na_n,na_i, rc_n,rc_i)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      A. Määttänen, J. Merikanto,H. Henschel, J. Duplissy, R. Makkonen, 
    !      I. K. Ortega and H. Vehkamäki
    !
    !      Copyright (C)2018 Määttänen et al. 2018
    !
    !      Please cite both the JGR paper and the Zenodo repository in 
    !      publications that use the code
    !
    !      contact
    !      -------
    !      Dr. Anni Määttänen
    !      LATMOS (Laboratoire ATMospheres et Observations Spatiales)
    !      Tour 45-46, 4eme etage
    !      Boite 102 / 4 Place Jussieu
    !      FR-75005
    !      Paris, France
    !      Phone: +33(0)144274970
    !      Email: anni.maattanen@latmos.ipsl.fr
    !      Email: joonas.merikanto@fmi.fi
    !      Email: hanna.vehkamaki@helsinki.fi
    !
    !      purpose
    !      -------
    !      neutral and ion-induced sulfuric acid-water particle formation rate
    !      of critical clusters
    !
    !      interface
    !      ---------
    !
    !        input:
    !           t          ! temperature in K 
    !           satrat     ! saturatio ratio of water (between zero and 1)
    !           rhoa       ! sulfuric acid concentration in 1/cm3
    !           csi        ! Ion condensation sink (s-1)
    !           airn       ! Air molecule concentration in (cm-3)
    !           ipr        ! Ion pair production rate (cm-3 s-1)
    !
    !        output:
    !           jnuc_n     ! Neutral nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
    !           jnuc_i     ! Charged nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
    !           na_n       ! sulfuric acid molecules in the neutral critical cluster
    !           na_i       ! sulfuric acid molecules in the charged critical cluster
    !           rc_n       ! radius of the neutral critical cluster in nm 
    !           rc_i       ! radius of the charged critical cluster in nm 
    !
    !      method
    !      ------
    !      Calculates parametrized values for neutral and ion-induced sulfuric 
    !      acid-water particle formation rate of critical clusters, 
    !      number of particles in the critical clusters, the radii of the 
    !      critical clusters in H2O-H2SO4-ion system if temperature, saturation 
    !      ratio of water, sulfuric acid concentration, and, optionally, 
    !      either condensation sink due to pre-existing particles and 
    !      ion pair production rate,
    !      or atmospheric concentration of negative ions are given. 
    !
    !      The code calculates also the kinetic limit and the particle 
    !      formation rate above this limit (in which case we set ntot=1 and na=1)
    !
    !      reference
    !      ---------
    !      A. Määttänen, J. Merikanto, H. Henschel, J. Duplissy, R. Makkonen, 
    !      I. K. Ortega and H. Vehkamäki (2018), New parameterizations for 
    !      neutral and ion-induced sulfuric acid-water particle formation in 
    !      nucleation and kinetic regimes, J. Geophys. Res. Atmos., 
    !      122, doi:10.1002/2017JD027429.
    !
    !      A. Määttänen, J. Merikanto, H. Henschel, J. Duplissy, R. Makkonen, 
    !      I. K. Ortega, and H. Vehkamäki, (2018, April 13). 
    !      Revised release of a Fortran code including the particle formation 
    !      parameterizations published in Määttänen et al., JGR D, 2018 
    !      (Version v1.0). Journal of Geophysical Research Atmospheres. 
    !      Zenodo. http://doi.org/10.5281/zenodo.1217782
    !
    !      Brasseur, G., and A.  Chatel (1983),  paper  presented  at  the  
    !      9th  Annual  Meeting  of  the  European Geophysical Society, 
    !      Leeds, Great Britain, August 1982.
    !  
    !      Dunne, Eimear M., et al.(2016), Global atmospheric particle formation 
    !      from CERN CLOUD measurements,
    !      Science 354.6316, 1119-1124.   
    !
    !
    !      modifications
    !      -------------
    !      MSK CHANGES: set variables from double precision to real(dp)
    !      Use the same ion pair production rate as in yu_timn_nucleation
    !
    !------------------------------------------------------------------

     implicit none

        real( dp),intent(in) :: t          ! temperature in K 
        real( dp),intent(in) :: satrat     ! saturatio ratio of water (between zero and 1) = RH
        real( dp),intent(in) :: ch2so4     ! sulfuric acid concentration in 1/m3
        real( dp),intent(in) :: airn       ! Air molecule concentration in (cm-3)

        real( dp),intent(out) :: jnuc_n   ! Neutral nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
        real( dp),intent(out) :: jnuc_i   ! Charged nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
        real( dp),intent(out) :: na_n     ! sulfuric acid molecules in the neutral critical cluster
        real( dp),intent(out) :: na_i     ! sulfuric acid molecules in the charged critical cluster
        real( dp),intent(out) :: rc_n     ! radius of the charged critical cluster in nm 
        real( dp),intent(out) :: rc_i     ! radius of the charged critical cluster in nm 

        ! OPtional output
        real( dp) :: ntot_n               ! total number of molecules in the neutral critical cluster
        real( dp) :: ntot_i               ! total number of molecules in the charged critical cluster
        real( dp) :: x_n                  ! mole fraction of H2SO4 in the neutral critical cluster 
        real( dp) :: x_i                  ! mole fraction of H2SO4 in the charged critical cluster 
                                          ! (note that x_n=x_i in nucleation regime)
        real( dp) :: n_i                  ! number of ion pairs in air (cm-3) 
        real( dp) :: rhoatres             ! treshold concentration of h2so4 (1/cm^3) 
                                          ! for neutral kinetic nucleation
        logical   :: kinetic_n            ! true if kinetic neutral nucleation
        logical   :: kinetic_i            ! true if kinetic ion-induced nucleation

        ! Local
        real( dp) :: rhoa          ! sulfuric acid concentration in 1/cm3
        real( dp) :: csi           ! Ion condensation sink (s-1)
        real( dp) :: ipr           ! Ion pair production rate (cm-3 s-1)
        real( dp) :: x             ! mole fraction of H2SO4 in the critical cluster 
        real( dp) :: satratln      ! bounded water saturation ratio for neutral case (between 5.E-6 - 1.0)
        real( dp) :: satratli      ! bounded water saturation ratio for ion-induced case (between 1.E-7 - 0.95)
        real( dp) :: rhoaln        ! bounded concentration of h2so4 for neutral case (between 10^10 - 10^19 m-3)
        real( dp) :: rhoali        ! bounded concentration of h2so4 for ion-induced case (between 10^10 - 10^22 m-3)
        real( dp) :: tln           ! bounded temperature for neutral case (between 165-400 K)
        real( dp) :: tli           ! bounded temperature for ion-induced case (195-400 K)
        real( dp) :: kinrhotresn   ! threshold sulfuric acid for neutral kinetic nucleation   
        real( dp) :: kinrhotresi   ! threshold sulfuric acid for ion-induced kinetic nucleation 
        real( dp) :: jnuc_i1       ! Ion-induced rate for n_i=1 cm-3 
        real( dp) :: xloss         ! Ion loss rate 
        real( dp) :: recomb        ! Ion-ion recombination rate 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MSK 10.09.2020: Unit conversion of inputs
!                 and other initialisation
!                 Use the same ion pair production rate as in yu_timn_nucleation

        rhoa = ch2so4*1.E-6_dp      ! [H2SO4] in #/cm3
        csi  = 1.0/480.             ! Inverse lifetime of ions
        ipr  = 2.00E+00_dp          ! ionization rate (ion-pairs cm-3 s-1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !--- 0) Initializations:

        kinetic_n=.false.
        kinetic_i=.false.
        jnuc_n=0.0
        jnuc_i=0.0
        ntot_n=0.0
        ntot_i=0.0
        na_n=0.0
        na_i=0.0
        rc_n=0.0
        rc_i=0.0
        x=0.0
        x_n=0.0
        x_i=0.0
        satratln=satrat
        satratli=satrat
        rhoaln=rhoa
        rhoali=rhoa
        tln=t
        tli=t
        n_i=0.0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !   write(6,*) 'Input:'
  !   write(6,*) ' T(K)    RH(%)      SA(cm-3)     csi (s-1)    airn(cm-3)   ipr(s-1cm-3)'
  !   write(6,'(F7.1,F11.5,3ES13.4,F8.2)') t,satrat*100,rhoa,csi,airn, ipr
  !   write(6,*)


  !Temperature bounds
         if(t.le.165.) then
           write(6,*) 'Warning (newbinapara): temperature < 165.0 K, using 165.0 K in neutral nucleation calculation'
           tln=165.0
         end if
         if(t.le.195.) then
           write(6,*) 'Warning (newbinapara): temperature < 195.0 K, using 195.0 K in ion-induced nucleation calculation'
           tli=195.0
         end if
         if(t.ge.400.) then
           write(6,*) 'Warning (newbinapara): temperature > 400. K, using 400. K in nucleation calculations'
           tln=400.
           tli=400.
         end if

  ! Saturation ratio bounds
        if(satrat.lt.1.E-7) then
           write(6,*) 'Warning: (newbinapara) saturation ratio of'
           write(6,*) 'water < 1.e-7, using 1.e-7 in ion-induced nucleation calculation'
           satratli=1.E-7
        end if
        if(satrat.lt.1.E-5) then
           write(6,*) 'Warning (newbinapara): saturation ratio of'
           write(6,*) 'water < 1.e-5, using 1.e-5 in neutral nucleation calculation'
           satratln=1.E-5
        end if
        if(satrat.gt.0.95) then
           write(6,*) 'Warning (newbinapara): saturation ratio of'
           write(6,*) 'water > 0.95, using 0.95 in ion-induced nucleation calculation'
           satratli=0.95
        end if
        if(satrat.gt.1.0) then
           write(6,*) 'Warning (newbinapara): saturation ratio of'
           write(6,*) 'water > 1 using 1 in neutral nucleation calculation'
           satratln=1.0
        end if

  ! Sulfuric acid concentration bounds
        if(rhoa.le.1.e4) then
           write(6,*) 'Warning (newbinapara): sulfuric acid < 1e4 1/cm3, '
           write(6,*) 'using 1e4 1/cm3 in nucleation calculation'
           rhoaln=1.e4
           rhoali=1.e4
        end if
        if(rhoa.gt.1.e13) then
           write(6,*) 'Warning (newbinapara): sulfuric acid > 1e13 1/cm3, '
           write(6,*) 'using 1e13 1/cm3 in neutral nucleation calculation'
           rhoaln=1.e13
        end if
        if(rhoa.gt.1.e16) then
           write(6,*) 'Warning (newbinapara): sulfuric acid concentration '
           write(6,*) ' > 1e16 1/cm3, using 1e16 1/cm3 in ion-induced nucleation calculation'
           rhoali=1.e16
        end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Critical cluster composition (valid for both cases, bounds not used here) 
        x_n=  7.9036365428891719e-1 - 2.8414059650092153e-3*tln +  &
              1.4976802556584141e-2*LOG(satratln)               -  &
              2.4511581740839115e-4*tln*LOG(satratln)           +  &
              3.4319869471066424e-3 *LOG(satratln)**2           -  &
              2.8799393617748428e-5*tln*LOG(satratln)**2        +  &
              3.0174314126331765e-4*LOG(satratln)**3            -  &
              2.2673492408841294e-6*tln*LOG(satratln)**3        -  &
              4.3948464567032377e-3*LOG(rhoaln)                 +  &
              5.3305314722492146e-5*tln*LOG(rhoaln)

        x_i=  7.9036365428891719e-1 - 2.8414059650092153e-3*tli +  &
              1.4976802556584141e-2*LOG(satratli)               -  &
              2.4511581740839115e-4*tli*LOG(satratli)           +  &
              3.4319869471066424e-3 *LOG(satratli)**2           -  &
              2.8799393617748428e-5*tli*LOG(satratli)**2        +  &
              3.0174314126331765e-4*LOG(satratli)**3            -  &
              2.2673492408841294e-6*tli*LOG(satratli)**3        -  &
              4.3948464567032377e-3*LOG(rhoali)                 +  &
              5.3305314722492146e-5*tli*LOG(rhoali)
       
        x_n=MIN(MAX(x_n,1.e-30),1.) 
        x_i=MIN(MAX(x_i,1.e-30),1.) 



  !Neutral nucleation
  
  !Kinetic limit check
        if (satratln .ge. 1.e-2 .and. satratln .le. 1.) then
          kinrhotresn=exp(7.8920778706888086e+1 +                               & 
           7.3665492897447082*satratln - 1.2420166571163805e+4/tln            + &
           (-6.1831234251470971e+2*satratln)/tln - 2.4501159970109945e-2*tln  - &
           1.3463066443605762e-2*satratln*tln + 8.3736373989909194e-06*tln**2 - &
           1.4673887785408892*Log(satratln)                                   + & 
           (-3.2141890006517094e+1*Log(satratln))/tln                         + &
           2.7137429081917556e-3*tln*Log(satratln)) !1/cm3     
          if(kinrhotresn.lt.rhoaln) kinetic_n=.true.
        endif
        if (satratln .ge. 1.e-4  .and. satratln .lt. 1.e-2) then     
          kinrhotresn=exp(7.9074383049843647e+1 -                               &
           2.8746005462158347e+1*satratln - 1.2070272068458380e+4/tln         + &
           (-5.9205040320056632e+3*satratln)/tln - 2.4800372593452726e-2*tln  - &
           4.3983007681295948e-2*satratln*tln + 2.5943854791342071e-5*tln**2  - &
           2.3141363245211317*Log(satratln)                                   + &
           (9.9186787997857735e+1*Log(satratln))/tln                          + &
           5.6819382556144681e-3*tln*Log(satratln)) !1/cm3
          if(kinrhotresn.lt.rhoaln) kinetic_n=.true.
        endif
        if (satratln .ge. 5.e-6  .and. satratln .lt. 1.e-4) then
          kinrhotresn=exp(8.5599712000361677e+1 +                               & 
           2.7335119660796581e+3*satratln - 1.1842350246291651e+4/tln         + &
           (-1.2439843468881438e+6*satratln)/tln - 5.4536964974944230e-2*tln  + &
           5.0886987425326087*satratln*tln + 7.1964722655507067e-5*tln**2     - &
           2.4472627526306372*Log(satratln)                                   + &
           (1.7561478001423779e+2*Log(satratln))/tln                          + &
           6.2640132818141811e-3*tln*Log(satratln)) !1/cm3
          if(kinrhotresn.lt.rhoaln) kinetic_n=.true. 
        endif



        if(kinetic_n) then    
           ! Dimer formation rate
           jnuc_n=1.E6*(2.*0.3E-9)**2.*sqrt(8.*3.141593*1.38E-23*               &
                  (1./(1.661e-27*98.07)+1./(1.661e-27*98.07)))/2.*sqrt(t)*rhoa**2.
           ntot_n=1. !set to 1 
! The critical cluster contains one molecule, but the produced cluster contains 2 molecules
           na_n=1.   
           x_n=na_n/ntot_n  ! so also set this to 1
           rc_n=0.3E-9
        else
           jnuc_n= 2.1361182605986115e-1 + 3.3827029855551838 *tln              -  &
                   3.2423555796175563e-2*tln**2                                 +  &
                   7.0120069477221989e-5*tln**3 +8.0286874752695141/x_n         +  &
                   (-2.6939840579762231e-1)*LOG(satratln)                       +  &
                   1.6079879299099518*tln*LOG(satratln)                         +  &
                   (-1.9667486968141933e-2)*tln**2*LOG(satratln)                +  &
                   5.5244755979770844e-5*tln**3*LOG(satratln)                   +  &
                   (7.8884704837892468*LOG(satratln))/x_n                       +  &
                   4.6374659198909596*LOG(satratln)**2                          -  &
                   8.2002809894792153e-2*tln*LOG(satratln)**2                   +  &
                   8.5077424451172196e-4*tln**2*LOG(satratln)**2                +  &
                   (-2.6518510168987462e-6)*tln**3*LOG(satratln)**2             +  &
                   (-1.4625482500575278*LOG(satratln)**2)/x_n                   -  &
                   5.2413002989192037e-1*LOG(satratln)**3                       +  &
                   5.2755117653715865e-3*tln*LOG(satratln)**3                   +  &
                   (-2.9491061332113830e-6)*tln**2*LOG(satratln)**3             +  &
                   (-2.4815454194486752e-8)*tln**3*LOG(satratln)**3             +  &
                   (-5.2663760117394626e-2*LOG(satratln)**3)/x_n                +  &
                   1.6496664658266762*LOG(rhoaln)                               +  &
                   (-8.0809397859218401e-1)*tln*LOG(rhoaln)                     +  &
                   8.9302927091946642e-3*tln**2*LOG(rhoaln)                     +  &
                   (-1.9583649496497497e-5)*tln**3*LOG(rhoaln)                  +  &
                   (-8.9505572676891685*LOG(rhoaln))/x_n                        +  &
                   (-3.0025283601622881e+1)*LOG(satratln)*LOG(rhoaln)           +  &
                   3.0783365644763633e-1*tln*LOG(satratln)*LOG(rhoaln)          +  &
                   (-7.4521756337984706e-4)*tln**2*LOG(satratln)*LOG(rhoaln)    +  &
                   (-5.7651433870681853e-7)*tln**3*LOG(satratln)*LOG(rhoaln)    +  &
                   (1.2872868529673207*LOG(satratln)*LOG(rhoaln))/x_n           +  &
                   (-6.1739867501526535e-1)*LOG(satratln)**2*LOG(rhoaln)        +  &
                   7.2347385705333975e-3*tln*LOG(satratln)**2*LOG(rhoaln)       +  &
                   (-3.0640494530822439e-5)*tln**2*LOG(satratln)**2*LOG(rhoaln) + &
                   6.5944609194346214e-8*tln**3*LOG(satratln)**2*LOG(rhoaln)    +  &
                   (-2.8681650332461055e-2*LOG(satratln)**2*LOG(rhoaln))/x_n    +  &
                   6.5213802375160306*LOG(rhoaln)**2                            +  &
                   (-4.7907162004793016e-2)*tln*LOG(rhoaln)**2                  +  &
                   (-1.0727890114215117e-4)*tln**2*LOG(rhoaln)**2               +  &
                   5.6401818280534507e-7*tln**3*LOG(rhoaln)**2                  +  &
                   (5.4113070888923009e-1*LOG(rhoaln)**2)/x_n                   +  &
                   5.2062808476476330e-1*LOG(satratln)*LOG(rhoaln)**2           +  &
                   (-6.0696882500824584e-3)*tln*LOG(satratln)*LOG(rhoaln)**2    +  &
                   2.3851383302608477e-5*tln**2*LOG(satratln)*LOG(rhoaln)**2    +  &
                   (-1.5243837103067096e-8)*tln**3*LOG(satratln)*LOG(rhoaln)**2 +  &
                   (-5.6543192378015687e-2*LOG(satratln)*LOG(rhoaln)**2)/x_n    +  &
                   (-1.1630806410696815e-1)*LOG(rhoaln)**3                      +  &
                   1.3806404273119610e-3*tln*LOG(rhoaln)**3                     +  &
                   (-2.0199865087650833e-6)*tln**2*LOG(rhoaln)**3               +  &
                   (-3.0200284885763192e-9)*tln**3*LOG(rhoaln)**3               +  &
                   (-6.9425267104126316e-3*LOG(rhoaln)**3)/x_n
           jnuc_n=EXP(jnuc_n) 



           ntot_n =-3.5863435141979573e-3 - 1.0098670235841110e-1 *tln             +  &
                   8.9741268319259721e-4 *tln**2 - 1.4855098605195757e-6*tln**3    -  & 
                   1.2080330016937095e-1/x_n + 1.1902674923928015e-3*LOG(satratln) -  & 
                   1.9211358507172177e-2*tln*LOG(satratln)                         +  &
                   2.4648094311204255e-4*tln**2*LOG(satratln)                      -  &
                   7.5641448594711666e-7*tln**3*LOG(satratln)                      +  &
                   (-2.0668639384228818e-02*LOG(satratln))/x_n                     -  &
                   3.7593072011595188e-2*LOG(satratln)**2                          +  &
                   9.0993182774415718e-4 *tln*LOG(satratln)**2                     +  &
                   (-9.5698412164297149e-6)*tln**2*LOG(satratln)**2                +  &
                   3.7163166416110421e-8*tln**3*LOG(satratln)**2                   +  &
                   (1.1026579525210847e-2*LOG(satratln)**2)/x_n                    +  &
                   1.1530844115561925e-2 *LOG(satratln)**3                         -  &
                   1.8083253906466668e-4 *tln*LOG(satratln)**3                     +  &
                   8.0213604053330654e-7*tln**2*LOG(satratln)**3                   +  &
                   (-8.5797885383051337e-10)*tln**3*LOG(satratln)**3               +  &
                   (1.0243693899717402e-3*LOG(satratln)**3)/x_n                    +  &
                   (-1.7248695296299649e-2)*LOG(rhoaln)                            +  &
                   1.1294004162437157e-2*tln*LOG(rhoaln)                           +  &
                   (-1.2283640163189278e-4)*tln**2*LOG(rhoaln)                     +  &
                   2.7391732258259009e-7*tln**3*LOG(rhoaln)                        +  &
                   (6.8505583974029602e-2*LOG(rhoaln))/x_n                         +  &
                   2.9750968179523635e-1*LOG(satratln)*LOG(rhoaln)                 +  &
                   (-3.6681154503992296e-3)*tln*LOG(satratln)*LOG(rhoaln)          +  &
                   1.0636473034653114e-5*tln**2*LOG(satratln)*LOG(rhoaln)          +  &
                   5.8687098466515866e-9*tln**3*LOG(satratln)*LOG(rhoaln)          +  &
                   (-5.2028866094191509e-3*LOG(satratln)*LOG(rhoaln))/x_n          +  &
                   7.6971988880587231e-4*LOG(satratln)**2*LOG(rhoaln)              -  &
                   2.4605575820433763e-5*tln*LOG(satratln)**2*LOG(rhoaln)          +  &
                   2.3818484400893008e-7*tln**2*LOG(satratln)**2*LOG(rhoaln)       +  &
                   (-8.8474102392445200e-10)*tln**3*LOG(satratln)**2*LOG(rhoaln)   +  &
                   (-1.6640566678168968e-4*LOG(satratln)**2*LOG(rhoaln))/x_n       -  &
                   7.7390093776705471e-2*LOG(rhoaln)**2                            +  &
                   5.8220163188828482e-4*tln*LOG(rhoaln)**2                        +  &
                   1.2291679321523287e-6*tln**2*LOG(rhoaln)**2                     +  &
                   (-7.4690997508075749e-9)*tln**3*LOG(rhoaln)**2                  +  &
                   (-5.6357941220497648e-3*LOG(rhoaln)**2)/x_n                     +  &
                   (-4.7170109625089768e-3)*LOG(satratln)*LOG(rhoaln)**2           +  &
                   6.9828868534370193e-5*tln*LOG(satratln)*LOG(rhoaln)**2          +  &
                   (-3.1738912157036403e-7)*tln**2*LOG(satratln)*LOG(rhoaln)**2    +  &
                   2.3975538706787416e-10*tln**3*LOG(satratln)*LOG(rhoaln)**2      +  &
                   (4.2304213386288567e-4*LOG(satratln)*LOG(rhoaln)**2)/x_n        +  &
                   1.3696520973423231e-3*LOG(rhoaln)**3                            +  &
                   (-1.6863387574788199e-5)*tln*LOG(rhoaln)**3                     +  &
                   2.7959499278844516e-8*tln**2*LOG(rhoaln)**3                     +  &
                   3.9423927013227455e-11*tln**3*LOG(rhoaln)**3                    +  &
                   (8.6136359966337272e-5*LOG(rhoaln)**3)/x_n
             ntot_n=EXP(ntot_n)
     
             rc_n=EXP(-22.378268374023630+0.44462953606125100*x_n                  +  &
                      0.33499495707849131*LOG(ntot_n)) !in meters

             na_n=x_n*ntot_n
             if (na_n .lt. 1.) then
               write(6,*) 'Warning (newbinapara): number of acid molecules'
               write(6,*) '< 1 in nucleation regime, setting na_n=1'
               na_n=1.0
             endif
        endif
  

  ! Set the neutral nucleation rate to 0.0 if less than 1.0e-7      
        if(jnuc_n.lt.1.e-7) then
           jnuc_n=0.0
        endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Threshold neutral nucleation rate (j > 1/cm3s) parameterization (can be commented out if not needed)
        if (tln .ge. 310.) then
           rhoatres=exp(-2.8220714121794250 + 1.1492362322651116e+1*satratln           - &
                    3.3034839106184218e+3/tln + (-7.1828571490168133e+2*satratln)/tln  + & 
                    1.4649510835204091e-1*tln - 3.0442736551916524e-2*satratln*tln     - &
                    9.3258567137451497e-5*tln**2 - 1.1583992506895649e+1*Log(satratln) + & 
                    (1.5184848765906165e+3*Log(satratln))/tln                          + &
                    1.8144983916747057e-2*tln*Log(satratln)) !1/cm3
        endif
        if (tln .gt. 190. .and. tln .lt. 310.) then
           rhoatres=exp(-3.1820396091231999e+2 + 7.2451289153199676*satratln           + &
                    2.6729355170089486e+4/tln + (-7.1492506076423069e+2*satratln)/tln  + &
                    1.2617291148391978*tln - 1.6438112080468487e-2*satratln*tln        - &
                    1.4185518234553220e-3*tln**2 - 9.2864597847386694*Log(satratln)    + &
                    (1.2607421852455602e+3*Log(satratln))/tln                          + &
                    1.3324434472218746e-2*tln*Log(satratln)) !1/cm3
        endif
        if (tln .lt. 185. .and. tln .gt. 155.) then
           rhoatres=1.1788859232398459e+5 - 1.0244255702550814e+4*satratln             + &
                    4.6815029684321962e+3*satratln**2 -1.6755952338499657e+2*tln
        endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  ! Ion-induced nucleation parameterization 
 
        if(ipr.gt.0.0) then ! if the ion production rate is above zero
     
     ! Calculate the ion induced nucleation rate wrt. concentration of 1 ion/cm3
     
           kinrhotresi = 5.3742280876674478e1                                      -  &
                         6.6837931590012266e-3 *log(satratli)**(-2)                -  &
                         1.0142598385422842e-01 * log(satratli)**(-1)              -  &
                         6.4170597272606873e+00 * log(satratli)                    -  &
                         6.4315798914824518e-01 * log(satratli)**2                 -  &
                         2.4428391714772721e-02 * log(satratli)**3                 -  &
                         3.5356658734539019e-04 * log(satratli)**4                 +  &
                         2.5400015099140506e-05 * tli * log(satratli)**(-2)        -  &
                         2.7928900816637790e-04 * tli * log(satratli)**(-1)        +  &
                         4.4108573484923690e-02 * tli * log(satratli)              +  &
                         6.3943789012475532e-03 * tli * log(satratli)**(2)         +  &
                         2.3164296174966580e-04 * tli * log(satratli)**(3)         +  &
                         3.0372070669934950e-06 * tli * log(satratli)**4           +  &
                         3.8255873977423475e-06 * tli**2 * log(satratli)**(-1)     -  &
                         1.2344793083561629e-04 * tli**2 * log(satratli)           -  &
                         1.7959048869810192e-05 * tli**2 * log(satratli)**(2)      -  &
                         3.2165622558722767e-07 * tli**2 * log(satratli)**3        -  &
                         4.7136923780988659e-09 * tli**3 * log(satratli)**(-1)     +  &
                         1.1873317184482216e-07 * tli**3 * log(satratli)           +  &
                         1.5685860354866621e-08 * tli**3 * log(satratli)**2        -  &
                         1.4329645891059557e+04 * tli**(-1)                        +  &
                         1.3842599842575321e-01 * tli                              -  &
                         4.1376265912842938e-04 * tli**(2) + 3.9147639775826004e-07 * tli**3
     
           kinrhotresi=exp(kinrhotresi) !1/cm3


          ! write(6,*) 'ion kinetic limit',kinrhotresi,rhoali

           if(kinrhotresi.lt.rhoali) kinetic_i=.true.


           if(kinetic_i) then    
               jnuc_i1=1.0E6*(0.3E-9 + 0.487E-9)**2.*sqrt(8.*3.141593*1.38E-23*(1./(1.661e-27*98.07)+ &
                       1./(1.661e-27*98.07))) * sqrt(tli)*rhoali   !1/cm3s  
               ntot_i=1. !set to 1 
               na_i=1.
               x_i=na_i/ntot_i  ! so also set this to 1
               rc_i=0.487E-9

              ! write(6,*) 'ion kinetic 1: ',jnuc_i1,na_i,rc_i

           else
               jnuc_i1 = 3.0108954259038608e+01+tli*6.1176722090512577e+01+(tli**2)                     * & 
                         8.7240333618891663e-01+(tli**3)* (-4.6191788649375719e-03)                     + &
                         (tli**(-1))*8.3537059107024481e-01 + (1.5028549216690628e+01                   + &
                         tli*(-1.9310989753720623e-01)+(tli**2)*8.0155514634860480e-04+(tli**3)         * &
                         (-1.0832730707799128e-06)+(tli**(-1))*1.7577660457989019)                      * &
                         (LOG(satratli)**(-2)) + (-2.0487870170216488e-01 +  tli                        * &
                         1.3263949252910405e-03 +  (tli**2) * (-8.4195688402450274e-06)                 + &
                         (tli**3) * 1.6154895940993287e-08 + (tli**(-1)) * 3.8734212545203874e+01)      * &
                         (LOG(satratli)**(-2) * LOG(rhoali)) + (1.4955918863858371 +  tli               * &
                         9.2290004245522454e+01 +  (tli**2) * (-8.9006965195392618e-01)                 + &
                         (tli**3) * 2.2319123411013099e-03 + (tli**(-1)) * 4.0180079996840852e-03)      * &
                         (LOG(satratli)**(-1) * LOG(rhoali)**(-1)) + (7.9018031228561085 +  tli         * &
                         (-1.1649433968658949e+01) +  (tli**2) * 1.1400827854910951e-01                 + &
                         (tli**3) * (-3.1941526492127755e-04) + (tli**(-1)) *(-3.7662115740271446e-01)) * & 
                         (LOG(satratli)**(-1)) + (1.5725237111225979e+02 +  tli * (-1.0051649979836277) + &
                         (tli**2) * 1.1866484014507624e-03 + (tli**3) * 7.3557614998540389e-06          + &
                         (tli**(-1)) * 2.6270197023115189) * (LOG(satratli)**(-1) * LOG(rhoali))        + &
                         (-1.6973840122470968e+01 +  tli * 1.1258423691432135e-01 +  (tli**2)           * &
                         (-2.9850139351463793e-04) + (tli**3) * 1.4301286324827064e-07 + (tli**(-1))    * &
                         1.3163389235253725e+01) * (LOG(satratli)**(-1) * LOG(rhoali)**2)               + &
                         (-1.0399591631839757 +  tli * 2.7022055588257691e-03 +  (tli**2)               * &
                         (-2.1507467231330936e-06) + (tli**3) * 3.8059489037584171e-10 + (tli**(-1))    * &
                         1.5000492788553410e+02) * (LOG(satratli)**(-1) * LOG(rhoali)**3)               + &
                         (1.2250990965305315 +  tli * 3.0495946490079444e+01 +  (tli**2)                * &
                         2.1051563135187106e+01 + (tli**3) * (-8.2200682916580878e-02) + (tli**(-1))    * &
                         2.9965871386685029e-02) * (LOG(rhoali)**(-2)) + (4.8281605955680433 +  tli     * &
                         1.7346551710836445e+02 +  (tli**2) * (-1.0113602140796010e+01) + (tli**3)      * &
                         3.7482518458685089e-02 + (tli**(-1)) * (-1.4449998158558205e-01))              * &
                         (LOG(rhoali)**(-1)) + (2.3399230964451237e+02 +  tli *(-2.3099267235261948e+01)+ &
                         (tli**2) * 8.0122962140916354e-02 + (tli**3) * 6.1542576994557088e-05          + &
                         (tli**(-1)) * 5.3718413254843007) * (LOG(rhoali)) + (1.0299715519499360e+02    + &
                         tli * (-6.4663357203364136e-02) +  (tli**2) * (-2.0487150565050316e-03)        + &
                         (tli**3) * 8.7935289055530897e-07 + (tli**(-1)) * 3.6013204601215229e+01)      * &
                         (LOG(rhoali)**2) + (-3.5452115439584042 +  tli * 1.7083445731159330e-02        + &
                         (tli**2) * (-1.2552625290862626e-05) + (tli**3) * 1.2968447449182847e-09       + &
                         (tli**(-1)) * 1.5748687512056560e+02) * (LOG(rhoali)**3)                       + &
                         (2.2338490119517975 +  tli * 1.0229410216045540e+02 +  (tli**2)                * &
                         (-3.2103611955174052) + (tli**3) * 1.3397152304977591e-02 + (tli**(-1))        * &
                         (-2.4155187776460030e-02)) * (LOG(satratli)* LOG(rhoali)**(-2))                + &
                         (3.7592282990713963 +  tli * (-1.5257988769009816e+02) +  (tli**2)             * &
                         2.6113805420558802 + (tli**3) * (-9.0380721653694363e-03) + (tli**(-1))        * &
                         (-1.3974197138171082e-01)) * (LOG(satratli)* LOG(rhoali)**(-1))                + &
                         (1.8293600730573988e+01 +  tli * 1.8344728606002992e+01 +  (tli**2)            * &
                         (-4.0063363221106751e-01) + (tli**3) * 1.4842749371258522e-03 + (tli**(-1))    * &
                         1.1848846003282287) * (LOG(satratli)) + (-1.7634531623032314e+02 +  tli        * &
                         4.9011762441271278 +  (tli**2) * (-1.3195821562746339e-02) + (tli**3)          * &
                         (-2.8668619526430859e-05) + (tli**(-1)) * (-2.9823396976393551e-01))           * &
                         (LOG(satratli)* LOG(rhoali)) + (-3.2944043694275727e+01 +  tli                 * &
                         1.2517571921051887e-01 +  (tli**2) * 8.3239769771186714e-05 + (tli**3)         * &
                         2.8191859341519507e-07 + (tli**(-1)) * (-2.7352880736682319e+01))              * &
                         (LOG(satratli)* LOG(rhoali)**2) + (-1.1451811137553243 +  tli                  * &
                         2.0625997485732494e-03 +  (tli**2) * (-3.4225389469233624e-06) + (tli**3)      * &
                         4.4437613496984567e-10 + (tli**(-1)) * 1.8666644332606754e+02)                 * &
                         (LOG(satratli)* LOG(rhoali)**3) + (3.2270897099493567e+01 +  tli               * &
                         7.7898447327513687e-01 +  (tli**2) * (-6.5662738484679626e-03) + (tli**3)      * &
                         3.7899330796456790e-06 + (tli**(-1)) * 7.1106427501756542e-01)                 * &
                         (LOG(satratli)**2 * LOG(rhoali)**(-1)) + (-2.8901906781697811e+01 +  tli       * &
                         (-1.5356398793054860) +  (tli**2) * 1.9267271774384788e-02 + (tli**3)          * &
                         (-5.3886270475516162e-05) + (tli**(-1)) * 5.0490415975693426e-01)              * &
                         (LOG(satratli)**2) + (3.3365683645733924e+01 +  tli *(-3.6114561564894537e-01) + &
                         (tli**2) * 9.2977354471929262e-04 + (tli**3) * 1.9549769069511355e-07          + &
                         (tli**(-1)) * (-8.8865930095112855)) * (LOG(satratli)**2 * LOG(rhoali))        + &
                         (2.4592563042806375 +  tli * (-8.3227071743101084e-03) +  (tli**2)             * &
                         8.2563338043447783e-06 + (tli**3) * (-8.4374976698593496e-09) + (tli**(-1))    * &
                         (-2.0938173949893473e+02)) * (LOG(satratli)**2 * LOG(rhoali)**2)               + &
                         (4.4099823444352317e+01 +  tli * 2.5915665826835252 +  (tli**2)                * &
                         (-1.6449091819482634e-02) + (tli**3) * 2.6797249816144721e-05 + (tli**(-1))    * &
                         5.5045672663909995e-01)* satratli
               jnuc_i1=EXP(jnuc_i1)


               ntot_i = abs((-4.8324296064013375e+04 +  tli * 5.0469120697428906e+02 +  (tli**2)        * &
                       (-1.1528940488496042e+00) + (tli**(-1)) * (-8.6892744676239192e+02) + (tli**(3)) * &
                       4.0030302028120469e-04) + (-6.7259105232039847e+03 +  tli                        * & 
                       1.9197488157452008e+02 +  (tli**2) * (-1.3602976930126354e+00)                   + &
                       (tli**(-1)) * (-1.1212637938360332e+02) + (tli**(3)) * 2.8515597265933207e-03)   * &
                       LOG(satratli)**(-2) * LOG(rhoali)**(-2)                                          + &
                       (2.6216455217763342e+02 +  tli * (-2.3687553252750821e+00) +  (tli**2)           * &
                       7.4074554767517521e-03 + (tli**(-1)) *(-1.9213956820114927e+03) + (tli**(3))     * &
                       (-9.3839114856129453e-06)) * LOG(satratli)**(-2) + (3.9652478944137344e+00 +  tli* &
                       1.2469375098256536e-02 +  (tli**2) * (-9.9837754694045633e-05) + (tli**(-1))     * &
                       (-5.1919499210175138e+02) + (tli**(3)) * 1.6489001324583862e-07)                 * &
                       LOG(satratli)**(-2) * LOG(rhoali) + (2.4975714429096206e+02 +  tli               * &
                       1.7107594562445172e+02 +  (tli**2) * (-7.8988711365135289e-01) + (tli**(-1))     * &
                       (-2.2243599782483177e+01) + (tli**(3)) * (-1.6291523004095427e-04))              * &
                       LOG(satratli)**(-1) * LOG(rhoali)**(-2) + (-8.9270715592533611e+02 +  tli        * &
                       1.2053538883338946e+02 +  (tli**2) * (-1.5490408828541018e+00) + (tli**(-1))     * &
                       (-1.1243275579419826e+01) + (tli**(3)) * 4.8053105606904655e-03)                 * &
                       LOG(satratli)**(-1) * LOG(rhoali)**(-1) + (7.6426441642091631e+03 +  tli         * &
                       (-7.1785462414656578e+01) +  (tli**2) * 2.3851864923199523e-01 + (tli**(-1))     * &
                       8.5591775688708395e+01 + (tli**(3)) * (-3.7000473243342858e-04))                 * &
                       LOG(satratli)**(-1) + (-5.1516826398607911e+01 +  tli * 9.1385720811460558e-01   + &
                       (tli**2) * (-3.5477100262158974e-03) + (tli**(-1)) * 2.7545544507625586e+03      + &
                       (tli**(3)) * 5.4708262093640928e-06) * LOG(satratli)**(-1) * LOG(rhoali)         + &
                       (-3.0386767129196176e+02 +  tli * (-1.1033438883583569e+04) +  (tli**2)          * &
                       8.1296859732896067e+01 + (tli**(-1)) * 1.2625883141097162e+01 + (tli**(3))       * &
                       (-1.2728497822219101e-01)) * LOG(rhoali)**(-2) + (-3.3763494256461472e+03 +  tli * &
                       3.1916579136391006e+03 +  (tli**2) * (-2.7234339474441143e+01) + (tli**(-1))     * &
                       (-2.1897653262707397e+01) + (tli**(3)) * 5.1788505812259071e-02)                 * &
                       LOG(rhoali)**(-1) + (-1.8817843873687068e+03 +  tli * 4.3038072285882070e+00     + &
                       (tli**2) * 6.6244087689671860e-03 + (tli**(-1)) * (-2.7133073605696295e+03)      + &
                       (tli**(3))*(-1.7951557394285043e-05)) * LOG(rhoali) + (-1.7668827539244447e+02   + &
                       tli * 4.8160932330629913e-01 +  (tli**2)*(-6.3133007671100293e-04) + (tli**(-1)) * &
                       2.5631774669873157e+04 + (tli**(3)) * 4.1534484127873519e-07) * LOG(rhoali)**(2) + &
                       (-1.6661835889222382e+03 +  tli * 1.3708900504682877e+03 +  (tli**2)             * &
                       (-1.7919060052198969e+01) + (tli**(-1)) * (-3.5145029804436405e+01) + (tli**(3)) * &
                       5.1047240947371224e-02) * LOG(satratli)* LOG(rhoali)**(-2)                       + &
                       (1.0843549363030939e+04 +  tli * (-7.3557073636139577e+01) +  (tli**2)           * &
                       1.2054625131778862e+00 + (tli**(-1)) * 1.9358737917864391e+02 + (tli**(3))       * &
                       (-4.2871620775911338e-03)) * LOG(satratli)* LOG(rhoali)**(-1)                    + &
                       (-2.4269802549752835e+03 +  tli * 1.1348265061941714e+01 +  (tli**2)             * &
                       (-5.0430423939495157e-02) + (tli**(-1)) * 2.3709874548950634e+03 + (tli**(3))    * &
                       1.4091851828620244e-04) * LOG(satratli) + (5.2745372575251588e+02 +  tli         * &
                       (-2.6080675912627314e+00) +  (tli**2) * 5.6902218056670145e-03 + (tli**(-1))     * &
                       (-3.2149319482897838e+04) + (tli**(3))*(-5.4121996056745853e-06)) * LOG(satratli)* &
                       LOG(rhoali) + (-1.6401959518360403e+01 +  tli * 2.4322962162439640e-01           + &
                       (tli**2) * 1.1744366627725344e-03 + (tli**(-1)) * (-8.2694427518413195e+03)      + &
                       (tli**(3)) *(-5.0028379203873102e-06))* LOG(satratli)**(2)                       + &
                       (-2.7556572017167782e+03 +  tli * 4.9293344495058264e+01 +  (tli**2)             * &
                       (-2.6503456520676050e-01) + (tli**(-1)) * 1.2130698030982167e+03 + (tli**(3))    * &
                       4.3530610668042957e-04)* LOG(satratli)**2 * LOG(rhoali)**(-1)                    + &
                       (-6.3419182228959192e+00 +  tli * 4.0636212834605827e-02 +  (tli**2)             * &
                       (-1.0450112687842742e-04) + (tli**(-1)) * 3.1035882189759656e+02 + (tli**(3))    * &
                       9.4328418657873500e-08)* LOG(satratli)**(-3) + (3.0189213304689042e+03 +  tli    * &
                       (-2.3804654203861684e+01) +  (tli**2) * 6.8113013411972942e-02 + (tli**(-1))     * & 
                       6.3112071081188913e+02 + (tli**(3)) * (-9.4460854261685723e-05))* (satratli)     * &
                       LOG(rhoali) + (1.1924791930673702e+04 +  tli * (-1.1973824959206000e+02)         + &
                       (tli**2) * 1.6888713097971020e-01 + (tli**(-1)) * 1.8735938211539585e+02         + &
                       (tli**(3)) * 5.0974564680442852e-04)* (satratli) + (3.6409071302482083e+01       + &
                       tli * 1.7919859306449623e-01 +  (tli**2)*(-1.0020116255895206e-03) + (tli**(-1)) * &
                       (-8.3521083354432303e+03) + (tli**(3)) * 1.5879900546795635e-06)* satratli       * &
                       LOG(rhoali)**(2))

               rc_i = (-3.6318550637865524e-08 +  tli * 2.1740704135789128e-09   +  (tli**2)            * &
                      (-8.5521429066506161e-12) + (tli**3) * (-9.3538647454573390e-15))                 + &
                      (2.1366936839394922e-08 +  tli * (-2.4087168827395623e-10) +  (tli**2)            * &
                      8.7969869277074319e-13 + (tli**3) *(-1.0294466881303291e-15))* LOG(satratli)**(-2)* &
                      LOG(rhoali)**(-1) + (-7.7804007761164303e-10 +  tli * 1.0327058173517932e-11      + &
                      (tli**2) * (-4.2557697639692428e-14) + (tli**3) * 5.4082507061618662e-17)         * &
                      LOG(satratli)**(-2) + (3.2628927397420860e-12 +  tli * (-7.6475692919751066e-14)  + &
                      (tli**2) * 4.1985816845259788e-16 + (tli**3) * (-6.2281395889592719e-19))         * &
                      LOG(satratli)**(-2) * LOG(rhoali) + (2.0442205540818555e-09 +  tli                * &
                      4.0441858911249830e-08 +  (tli**2) * (-3.3423487629482825e-10)                    + &
                      (tli**3) * 6.8000404742985678e-13)* LOG(satratli)**(-1) * LOG(rhoali)**(-2)       + &
                      (1.8381489183824627e-08 +  tli * (-8.9853322951518919e-09) +  (tli**2)            * &
                      7.5888799566036185e-11 + (tli**3) *(-1.5823457864755549e-13))* LOG(satratli)**(-1)* &
                      LOG(rhoali)**(-1) + (1.1795760639695057e-07 +  tli * (-8.1046722896375875e-10)    + &
                      (tli**2) * 9.1868604369041857e-14 + (tli**3) * 4.7882428237444610e-15)            * &
                      LOG(satratli)**(-1) + (-4.4028846582545952e-09 +  tli * 4.6541269232626618e-11    + &
                      (tli**2) *(-1.1939929984285194e-13) + (tli**3) * 2.3602037016614437e-17)          * &
                      LOG(satratli)**(-1) * LOG(rhoali) + (2.7885056884209128e-11 +  tli                * &
                      (-4.5167129624119121e-13) +  (tli**2) * 1.6558404997394422e-15                    + &
                      (tli**3) * (-1.2037336621218054e-18))* LOG(satratli)**(-1) * LOG(rhoali)**2       + &
                      (-2.3719627171699983e-09 +  tli * (-1.5260127909292053e-07) +  (tli**2)           * &
                      1.7177017944754134e-09 + (tli**3) *(-4.7031737537526395e-12))* LOG(rhoali)**(-2)  + &
                      (-5.6946433724699646e-09 +  tli * 8.4629788237081735e-09 +  (tli**2)              * &
                      (-1.7674135187061521e-10) + (tli**3) * 6.6236547903091862e-13)* LOG(rhoali)**(-1) + &
                      (-2.2808617930606012e-08 +  tli * 1.4773376696847775e-10 +  (tli**2)              * &
                      (-1.3076953119957355e-13) + (tli**3) * 2.3625301497914000e-16)* LOG(rhoali)       + &
                      (1.4014269939947841e-10 +  tli * (-2.3675117757377632e-12) +  (tli**2)            * &
                      5.1514033966707879e-15 + (tli**3) *(-4.8864233454747856e-18))* LOG(rhoali)**2     + &
                      (6.5464943868885886e-11 +  tli * 1.6494354816942769e-08 +  (tli**2)               * &
                      (-1.7480097393483653e-10) + (tli**3) * 4.7460075628523984e-13)* LOG(satratli)     * &
                      LOG(rhoali)**(-2) + (8.4737893183927871e-09 +  tli * (-6.0243327445597118e-09)    + &
                      (tli**2) * 5.8766070529814883e-11 + (tli**3) * (-1.4926748560042018e-13))         * &
                      LOG(satratli)* LOG(rhoali)**(-1) + (1.0761964135701397e-07 +  tli                 * &
                      (-1.0142496009071148e-09) +  (tli**2) * 2.1337312466519190e-12                    + &
                      (tli**3) * 1.6376014957685404e-15)* LOG(satratli) + (-3.5621571395968670e-09      + &
                      tli * 4.1175339587760905e-11 +  (tli**2) * (-1.3535372357998504e-13)              + &
                      (tli**3) * 8.9334219536920720e-17)* LOG(satratli)* LOG(rhoali)                    + &
                      (2.0700482083136289e-11 +  tli * (-3.9238944562717421e-13) +  (tli**2)            * &
                      1.5850961422040196e-15 + (tli**3) *(-1.5336775610911665e-18))* LOG(satratli)      * &
                      LOG(rhoali)**2 + (1.8524255464416206e-09 +  tli * (-2.1959816152743264e-11)       + &
                      (tli**2) * (-6.4478119501677012e-14) + (tli**3) * 5.5135243833766056e-16)         * &
                      LOG(satratli)**2 * LOG(rhoali)**(-1) + (1.9349488650922679e-09 +  tli             * &
                      (-2.2647295919976428e-11) +  (tli**2) * 9.2917479748268751e-14                    + &
                      (tli**3) * (-1.2741959892173170e-16))* LOG(satratli)**2                           + &
                      (2.1484978031650972e-11 +  tli * (-9.3976642475838013e-14) +  (tli**2)            * &
                      (-4.8892738002751923e-16) + (tli**3) * 1.4676120441783832e-18)* LOG(satratli)**2  * &
                      LOG(rhoali) + (6.7565715216420310e-13 +  tli * (-3.5421162549480807e-15)          + &
                      (tli**2) * (-3.4201196868693569e-18) + (tli**3) * 2.2260187650412392e-20)         * &
                      LOG(satratli)**3 * LOG(rhoali)

               na_i=x_i*ntot_i
               if (na_i .lt. 1.) then
                  write(6,*) 'Warning (newbinapara): number of acid molecules < 1'
                  write(6,*) 'in nucleation regime, setting na_n=1'
                  na_n=1.0
              endif
           endif


           jnuc_i=jnuc_i1 
     ! Ion loss rate (1/s)
           xloss=csi+jnuc_i
     
     ! Recombination (here following Brasseur and Chatel, 1983)   
           recomb=6.0e-8*sqrt(300./tli)+6.0e-26*airn*(300./tli)**4
     
     ! Small ion concentration in air (1/cm3) (following Dunne et al., 2016)
     ! max function is to avoid n_i to go practically zero at very high J_ion 
           n_i=max(0.01,(sqrt(xloss**2.0+4.0*recomb*ipr)-xloss)/(2.0*recomb))


     ! Ion-induced nucleation rate
     ! Min function is to ensure that max function above does not cause J_ion to overshoot 
           jnuc_i=min(ipr,n_i*jnuc_i1)

     ! Set the ion-induced nucleation rate to 0.0 if less than 1.0e-7      
           if(jnuc_i.lt.1.e-7) then
              jnuc_i=0.0
           endif

        endif


  !  write(6,*) 'Output, neutral nucleation:'
  !  write(6,*) ' J_n (cm-3s-1)  Ntot_n     x         Na_n   rc_n(nm)       Kinetic  rhoatres(cm-3)'
  !  write(6,'(ES12.4,3F10.2,ES13.4,L6,ES18.4)') jnuc_n,ntot_n,x_n,na_n,rc_n*1.E9,kinetic_n,rhoatres
  !  write(6,*)
  !  write(6,*) 'Output, ion-induced nucleation:'
  !  write(6,*) ' J_i (cm-3s-1)  Ntot_i     x         Na_i   rc_i(nm)       Kinetic  Ions(cm-3) '
  !  write(6,'(ES12.4,3F10.2,ES13.4,L6,ES18.4)') jnuc_i,ntot_i,x_i,na_i,rc_i*1.E9,kinetic_i, n_i


   end subroutine newbinapara



     subroutine binhomogeneous(t,rh,rhoa,nwtot,natot,rc,jnuc)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Hanna Vehkamaeki
    !
    !      contact
    !      -------
    !      Dr. Hanna Vehkamaeki
    !      Professor
    !      University of Helsinki
    !      Institute for Atmospheric and Earth System Research (INAR)
    !      PL 64, FI-00014 UHEL
    !      Helsinki, Finland
    !      Phone:  +358 50415747
    !      Email:  hann.vehkamaki@helsinki.fi 
    !
    !   *
    !        (c) Hanna Vehkamaeki, 2002
    !   *
    !
    !      purpose
    !      -------
    !      Calculates binary nucleation rate using revised theory,
    !      stauffer+binder&stauffer kinetics and noppel hydrate correction
    !      for homogeneous nucleation
    !      (binheterogeneous for heterogeneous nucleation)
    !
    !      interface
    !      ---------
    !
    !        t     temperature [K]
    !        rh    relative humidity %/100
    !        rhoa  concentration of h2so4 vapour [1/m^3]
    !        x     mole fraction in teh core of the critical cluster
    !        nwtot total number of water molecules in the critical cluster
    !        natot total number of h2so4 molecules in the critical cluster
    !        rc    radius of the critical cluster core [m]
    !        jnuc  nucleation rate [1/m^3s]
    !
    !        rhow  concentration of water vapour [1/m^3]
    !        rha   relative acidity %/100
    !
    !
    !      method
    !      ------
    !
    !
    !      reference
    !      ---------
    !      Vehkamaki, H., Kulmala, M., Napari, I., Lehtinen, K. E. J., 
    !      Timmreck, C., Noppel, M., and A. Laaksonen, An improved 
    !      parameterization for sulphuric acid-water nucleation rates for 
    !      tropospheric and stratospheric conditions, 
    !      J. Geophys. Res., 107, D22, 4622, doi:10.1029/2002JD002184, 2002.
    !
    !      modifications
    !      -------------
    !      2019-11-01 added by Matthias S. Karl (MSK):
    !      Extension for higher temperatures, 300K to 400K
    !      Vehkamaki, H., Kulmala, M., Lehtinen, K.E.J.
    !      Modelling Bindary Homogeneous Nucleation of Water-Sulfuric
    !      Acid Vapours: Parameterisation for High Temperature
    !      Emissions.
    !      Environ. Sci. Technol. 2003, 37, 3392-3398
    !      2020-06-24 MSK renamed function names
    !      vola       --> partial_molar_vol
    !      activi     --> activ_binary
    !      hydraatit2 --> sulfhydrates
    !      activityw  --> activity_water
    !      activitya  --> activity_acid
    !
    !------------------------------------------------------------------

    use gde_constants,    only      : pi,k_B
    use gde_toolbox,      only      : waterps,acidps,surten,roolq
    use gde_toolbox,      only      : sulfhydrates

    implicit none

        REAL( dp), intent(in)    :: t,rh,rhoa
        REAL( dp), intent(out)   :: rc,jnuc,nwtot,natot

        REAL( dp)                :: zero, one, pstand, ma, mb,rhow,rha
        REAL( dp), DIMENSION(2)  :: act, molvol
        REAL( dp), dimension(5)  :: kprod,rhoh,rhohm
        REAL( dp)                :: g, actw,acta, vw, va, vc, x
        REAL( dp)                :: hyd,hydpw0x,sute,hydpw
        REAL( dp)                :: Raa,Rww,Raw,sumww,sumaa,sumaw
        REAL( dp)                :: mclu, mb2, ma2,mscale,sqrtmsca
        REAL( dp)                :: racid, rwater
        REAL( dp)                :: na, nw, xguess
        REAL( dp)                :: term1, term2,term3, term4,r
        REAL( dp)                :: zeld, theta, rav, detr, deltag
        REAL( dp)                :: atx,btx,ctx,dtx,etx,ftx,gtx
        REAL( dp)                :: htx,itx,jtx
        REAL( dp)                :: antx,bntx,cntx,dntx,entx,fntx,gntx
        REAL( dp)                :: hntx,intx,jntx
        REAL( dp)                :: nctot

          ma=18.09*1.661E-27_dp
          mb=98.08*1.661E-27_dp
          pstand=1.01325E5_dp  ! Pa
          zero=0.
          one=1.

          rhow = (rh*waterps(t))/k_B/t
          rha  =  rhoa*k_B*t/acidps(t)

!MSK 07.04.2018 no BNS for high T and low H2SO4
          if ((t.gt.305).and.((rhoa/1.e6).lt.2.e09)) then
              write(6,*) "No binary nucl, T:",t," > 305 K, H2SO4 too low",rhoa/1.e6
            return
          endif
!MSK 01.11.2019 BNS according to Vehkamaeki et al., 2003 (EST)
          if ((t.gt.305).and.((rhoa/1.e6).ge.2.e09)) then

            ! critical cluster mole fraction of sulfuric acid 
            ! [Eq. 10 in Vehkamaeki et al., 2003]
             x = 0.847012 - 0.0029656*t - 0.00662266*LOG(rhoa/1.e6) +        &
                0.0000587835*t*LOG(rhoa/1e6) + 0.0592653*LOG(rh)    -      &
                0.000363192*t*LOG(rh) + 0.0230074*LOG(rh)**2        -      &        
                0.0000851374*t*LOG(rh)**2 + 0.00217417*LOG(rh)**3   -      &
                7.923d-6*t*LOG(rh)**3

            ! coefficients a(T,x) ... i(T,x) are functions of
            ! temperature and critical cluster mole fraction x
             atx = -0.00156975 - 0.134245*t + 0.100507*t**2         -        &
                  0.000460103*t**3 + 0.187416/x**2 + 0.0104122/x
             btx = 0.00195077 + 0.168038*t - 0.0225755*t**2         +        &
                  0.0000827149*t**3 + 0.0025029/x**2 + 0.0155215/x
             ctx = 0.000154084 - 0.0280301*t                        +        &
                  0.00154587*t**2 - 4.52701d-6*t**3 +                      &
                  0.0915323/x**2 + 0.0711652/x
             dtx = -0.00509267 - 0.00796846*t                       +        &
                  0.0000446828*t**2 - 8.79425d-8*t**3               +      &
                  0.133991/x**2 + 0.831112/x
             etx = -0.0227223 - 1.56512*t + 0.00380717*t**2         +        &
                  0.0000164109*t**3 + 1.29499/x**2 + 0.0474821/x
             ftx = 0.00310646 + 0.304518*t - 0.000564012*t**2       -        &
                  2.03267d-6*t**3 - 0.351584/x**2 + 0.103749/x
             gtx = 0.077543 - 0.00196315*t                          -        &
                  0.0000130412*t**2 + 6.62369d-8*t**3               +      &
                  0.011347/x**2 + 0.0972804/x
             htx = -0.153143 + 0.0575392*t - 0.000306511*t**2       -        &
                  2.96097d-8*t**3 - 0.0982514/x**2 + 0.336286/x
             itx = -0.552173 - 0.00207043*t                         +        &
                  0.0000144032*t**2 + 8.83d-9*t**3                  +      &
                  0.0119833/x**2 - 0.0700025/x
             jtx = 0.126544 - 0.00136029*t + 5.90598d-6*t**2        -        &
                  4.1715d-9*t**3 + 0.00170807/x**2 - 0.0064323/x

            !print *,'BSN logRH',LOG(rh)
            !print *,'BSN logRO',LOG(rhoa/1.e6)
            !print *,'BSN a',atx
            !print *,'BSN b',btx*LOG(rh)
            !print *,'BSN c',ctx*LOG(rh)**2
            !print *,'BSN d',dtx*LOG(rh)**3
            !print *,'BSN e',etx*LOG(rhoa/1.e6) 
            !print *,'BSN f',ftx*LOG(rh)*LOG(rhoa/1.e6) 
            !print *,'BSN g',gtx*LOG(rh)**2*LOG(rhoa/1.e6)  
            !print *,'BSN h',htx*LOG(rhoa/1.e6)**2
            !print *,'BSN i',itx*LOG(rh)*LOG(rhoa/1.e6)**2
            !print *,'BSN j',jtx*LOG(rhoa/1.e6)**3
            !print *,'BSN exp',atx+btx*LOG(rh)+ctx*LOG(rh)**2+dtx*LOG(rh)**3+&
            !                  etx*LOG(rhoa/1.e6)+ftx*LOG(rh)*LOG(rhoa/1.e6)+&
            !                  gtx*LOG(rh)**2*LOG(rhoa/1.e6)+htx*LOG(rhoa/1.e6)**2+&
            !                  itx*LOG(rh)*LOG(rhoa/1.e6)**2+jtx*LOG(rhoa/1.e6)**3

             ! nucleation rate [1/cm^3s] [Eq. 11 in Vehkamaeki et al., 2003]
             jnuc = EXP( atx + btx*LOG(rh)                          +        &
                   ctx*LOG(rh)**2                                   +      &
                   dtx*LOG(rh)**3                                   +      & 
                   etx*LOG(rhoa/1.e6)                               +      &
                   ftx*LOG(rh)*LOG(rhoa/1.e6)                       +      & 
                   gtx*LOG(rh)**2*LOG(rhoa/1.e6)                    +      &
                   htx*LOG(rhoa/1.e6)**2                            +      &
                   itx*LOG(rh)*LOG(rhoa/1.e6)**2                    +      &
                   jtx*LOG(rhoa/1.e6)**3 )
             !print *,'BSN jnuc',t,rh,rhoa/1.e6,x,jnuc

             ! nucleation rate [1/m^3s] 
             jnuc = jnuc * 1.e6_dp


            ! very low nucleation rate doesn't increase the particle conc.
            ! of nucleation mode but may cause funny behaviour in radius!
            ! thus the limit 10^-5 cm^-3 s^-1
             if (jnuc.lt.10.) then

                 jnuc  = 0.0
                 natot = 0.0
                 nwtot = 0.0
                 rc    = 0.5d-9

             else

             ! coefficients A(T,x) ... I(T,x) are functions of
             ! temperature and critical cluster mole fraction x
               antx = 7.51024d-6 + 0.000502054*t                    -        &
                     0.0000368602*t**2 + 1.08256d-6*t**3            -      &
                     0.000270282/x
               bntx = -4.30048d-6 - 0.000730133*t                   +        &
                     0.000252062*t**2 - 1.01648d-6*t**3             -      &
                     0.00114283/x
               cntx = -4.42156d-6 - 0.0023486*t                     +        &
                     3.0065d-7*t**2 + 2.44797d-8*t**3               -      &
                     0.00250226/x
               dntx = -0.000167057 + 0.000207504*t                  -        &
                     1.13013d-6*t**2 + 1.80268d-9*t**3              -      &
                     0.0168245/x
               entx = 0.0000985954 + 0.00451285*t                   -        &
                     0.0000512557*t**2 + 4.60749d-8*t**3            -      &
                     0.00214318/x
               fntx = 0.0000636528 - 0.00288529*t                   +        &
                     6.51706d-6*t**2 + 2.32601d-8*t**3              -      &
                     0.0110319/x
               gntx = 0.000449239 + 0.0000689416*t                  -        &
                     3.50302d-7*t**2 + 1.07451d-10*t**3             +      &
                     0.00169646/x
               hntx = 0.000831844 - 5.35108d-6*t                    +        &
                     1.66432d-6*t**2 - 3.05108d-9*t**3              -      &
                     0.000306251/x
               intx = 0.00355374 + 0.0000306009*t                   -        &
                     2.11004d-7*t**2 - 2.11436d-11*t**3             +      &
                     0.00074989/x
               jntx = -0.00143534 + 7.856d-6*t                      -        &
                     3.45128d-8*t**2 + 5.21547d-11*t**3             -      &
                     0.000021423/x

               ! total number of molecules in the critical cluster
               nctot = EXP( antx + bntx*LOG(rh)                     +        &
                      cntx*LOG(rh)**2                               +      &
                      dntx*LOG(rh)**3                               +      &
                      entx*LOG(rhoa/1.e6)                           +      &
                      fntx*LOG(rh)*LOG(rhoa/1.e6)                   +      &
                      gntx*LOG(rh)**2*LOG(rhoa/1.e6)                +      &
                      hntx*LOG(rhoa/1.e6)**2                        +      &
                      intx*LOG(rh)*LOG(rhoa/1.e6)**2                +      &
                      jntx*LOG(rhoa/1.e6)**3 )
               natot = nctot*x
               nwtot = nctot*(1.-x)
               ! radius of critical cluster [nm]
               rc = EXP( -1.6525507 + 0.45852848*x                  +         &
                    0.33483673*LOG(nctot) )
               ! radius of critical cluster [m]
               rc = rc * 1.e-9
             endif

             return
          endif

!MSK: for T<305K
          ! critical cluster mole fraction
          x=  0.7409967177282139 - 0.002663785665140117*t +                   &
            0.002010478847383187*LOG(rh) - 0.0001832894131464668*t*LOG(rh)+ &
            0.001574072538464286*LOG(rh)**2 -                               &
            0.00001790589121766952*t*LOG(rh)**2 +                           &
            0.0001844027436573778*LOG(rh)**3 -                              &
            1.503452308794887e-6*t*LOG(rh)**3 -                             &
            0.003499978417957668*LOG(rhoa/1.e6) +                           &
            0.0000504021689382576*t*LOG(rhoa/1.e6)
        ! partial molar volumes [m^3] in crit. cluster
          CALL partial_molar_vol(x,t,molvol)
          vw=molvol(1)             ! water [m^3]
          va=molvol(2)             ! h2so4 [m^3]
          CALL activ_binary(x,t,act)     ! activities in the critical cluster

          actw = act(1) !water
          acta = act(2) !h2so4
          sute = surten(x,t)       !surface tension of the critical cluster [N/m]

        CALL sulfhydrates(t,rh*waterps(t),hydpw,kprod,rhoh,rhohm)
        ! pw0   [Pa] saturation vapour pressure of pure water
          hydpw0x = 1.0        ! no hydrate correction to equilibrium vapour pressure
          hyd = hydpw0x/hydpw  ! hydrate correction to actual vapour pressure

        ! radius of the critical cluster [m]
          rc = 2.0*sute*(x*va + (1.0-x)*vw) /  &
             ((x*LOG(rha/acta)+(1.0-x)*LOG(rh/actw)+x*LOG(hyd))*t*k_B)

      ! critical radius must be greater than zero
          IF(rc .GT. 0.0) THEN
                g = 4.0*pi*rc*rc*sute/(3.0*t*k_B)        !deltaG/kT= free energy of
                                                              !the critical cluster
                vc=4.0*pi/3.0*rc*rc*rc                   !volume of the cluster [m^3]
                natot=vc*x/(x*va+(1.-x)*vw)              !number of acids in the cluster core
                nwtot=vc*(1.-x)/(x*va+(1.-x)*vw)         !number of waters in the cluster core
!            write(6,*) 'natot,nwtot',natot,nwtot

                mclu= natot*mb+nwtot*ma             !mass of the cluster [kg]

        ! radii of water and h2so4 molecules [m]
                rwater=(3.*ma/(4.*pi*roolq(zero,t)))**(1./3.)
                racid=(3.*mb/(4.*pi*roolq(one,t)))**(1./3.)

                mscale=1.e27_dp
                sqrtmsca=3.162277660e13_dp
                mclu=mclu*mscale
                mb2=mb*mscale
                ma2=ma*mscale

      ! calculating the growth tensor R
      ! binder&stauffer  advances in physics 25 343-396 1976
                sumww = (rc+rwater)**2.0*(1./ma2+1./mclu)**0.5*rhow
                sumww =sumww*sqrtmsca
                sumaw=0.0
                sumaa =(rc+racid)**2.0*(1./mclu + 1./mb2)**(0.5)*rhoa
                sumaa=sumaa*sqrtmsca
                Rww=  sumww*(8.0*pi*k_B*t)**0.5      !impingement rate of water
                Raw = sumaw*(8.0*pi*k_B*t)**0.5      !"mixed" impingement rate
                Raa = sumaa*(8.0*pi*k_B*t)**0.5      !impingement rate of acids

      !*****************************************************************
      ! Kulmala&viisanen '91 approx
      !    Zeldovich factor
                zeld  = SQRT(sute/(t*k_B))*((1.0-x)*vw +  &
                x*va)/(2.0*pi*rc*rc)
                theta = ATAN(x/(1.0-x))
                detr  = Rww*Raa - Raw*Raw
                rav   = Rww*SIN(theta)*SIN(theta) +       &
                Raa*COS(theta)*COS(theta) -         &
                2.0*Raw*SIN(theta)*COS(theta)
                rav   = detr/rav             !average growth rate
                term1 = rav*zeld             !kinetic pre-factor
      !*****************************************************************
        ! energy surface to the right level,
        ! distribution set to go through "experimental" w(1,2) (noppel)
                na=1.      ! mole fraction of a 2-hydrate
                nw=2.
                xguess= na/(na+nw)
      ! partial molar volumes [m^3] in the critical cluster
                CALL partial_molar_vol(xguess,t,molvol)
                 vw=molvol(1) !water [m^3]
                va=molvol(2) !h2so4 [m^3]
                r=(3./(4.*pi)*(na*va+nw*vw))**(1./3.)
                CALL activ_binary(xguess,t,act)  !activities
                deltag= 4.*pi*r*r*surten(xguess,t)/(k_B*t)
                deltag=deltag- na*LOG(rha/(act(2)*hydpw))
                deltag=deltag- nw*LOG(rh/act(1))
                term2=rhoa/hydpw*(rhow*k_B*t/pstand)**2.*kprod(2)
                term3=EXP(deltag)
                term4=EXP(-g)
       !*****************************

                jnuc=term1*term2*term3*term4        !nucleation rate [1/m^3s]

            ! very low nucleation rate doesn't increase the particle conc.
            ! of nucleation mode but may cause funny behaviour in radius!
            ! thus the limit 10^-5 cm^-3 s^-1
                if ( (jnuc.lt.10.) .OR. (natot .lt. 0.0) .OR. (nwtot .lt. 0.0) ) then
                 jnuc  = 0.0
                 natot = 0.0
                 nwtot = 0.0
                 x     = 0.0
                 rc    = 0.5d-9
                endif

          ELSE    ! critical cluster radius = zero
              jnuc  = 0.0
              natot = 0.0
              nwtot = 0.0
              x     = 0.0
              rc    = 0.5d-9
          ENDIF

   end subroutine binhomogeneous

  subroutine partial_molar_vol(xmole,t,molvol)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Hanna Vehkamaeki
    !
    !      contact
    !      -------
    !      Dr. Hanna Vehkamaeki
    !      Professor
    !      University of Helsinki
    !      Institute for Atmospheric and Earth System Research (INAR)
    !      PL 64, FI-00014 UHEL
    !      Helsinki, Finland
    !      Phone:  +358 50415747
    !      Email:  hann.vehkamaki@helsinki.fi 
    !
    !   *
    !        (c) Hanna Vehkamaeki, 2002
    !   *
    !
    !      purpose
    !      -------
    !      Calculates partial molar volumes for water&sulphuric acid
    !
    !      interface
    !      ---------
    !
    !        input:
    !           molvol(i)  [m^3] i=1 water, i=2 h2so4
    !           xmole =mole fraction of h2so4, t in K,
    !
    !
    !      method
    !      ------
    !      Calculates partial molar volumes for water&sulphuric acid
    !
    !      reference
    !      ---------
    !      Vehkamaki, H., Kulmala, M., Napari, I., Lehtinen, K. E. J., 
    !      Timmreck, C., Noppel, M., and A. Laaksonen, An improved 
    !      parameterization for sulphuric acid-water nucleation rates for 
    !      tropospheric and stratospheric conditions, 
    !      J. Geophys. Res., 107, D22, 4622, doi:10.1029/2002JD002184, 2002.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

    use gde_toolbox,      only      : roolq,roolqd

    implicit none

        REAL( dp), intent(in)                :: xmole,t
        REAL( dp), intent(out), DIMENSION(2) :: molvol
        REAL( dp)                            :: rho, drho,ma,mb

        ma=18.09*1.661E-27_dp
          mb=98.08*1.661E-27_dp
          rho = roolq(xmole,t)     ! density [kg/m^3]

   ! Derivative of the density [kg/m^3], with respect to mole fraction

          drho = roolqd(xmole,t)
          molvol(1) = ma/rho + xmole*(ma*(1.-xmole)+mb*xmole) * &
                  drho/(rho*rho)
          molvol(2) = mb/rho -                                  &
                  (1-xmole)*(ma*(1.-xmole)+mb*xmole)*drho/(rho*rho)
          IF (molvol(1).LE.0.0 .OR. molvol(2).LE.0.0) THEN
                  molvol(1) = ma/rho
                molvol(2) = mb/rho
          ENDIF

  end subroutine partial_molar_vol


  subroutine activ_binary (x,t,act)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Hanna Vehkamaeki
    !
    !      contact
    !      -------
    !      Dr. Hanna Vehkamaeki
    !      Professor
    !      University of Helsinki
    !      Institute for Atmospheric and Earth System Research (INAR)
    !      PL 64, FI-00014 UHEL
    !      Helsinki, Finland
    !      Phone:  +358 50415747
    !      Email:  hann.vehkamaki@helsinki.fi 
    !
    !   *
    !        (c) Hanna Vehkamaeki, 2002
    !   *
    !
    !      purpose
    !      -------
    !      Chemical activities of water and sulphuric acid in liquid phase
    !      calculated using Zeleznik thermodynamics
    !
    !      interface
    !      ---------
    !
    !        input:
    !           x = mole fraction of h2so4
    !           t = temperature in K
    !           act(1), act(2) = activities of water and h2so4, respectively
    !
    !      method
    !      ------
    !      Chemical activities of water and sulphuric acid in liquid phase
    !      calculated using Zeleznik thermodynamics
    !
    !      reference
    !      ---------
    !      Vehkamaki, H., Kulmala, M., Napari, I., Lehtinen, K. E. J., 
    !      Timmreck, C., Noppel, M., and A. Laaksonen, An improved 
    !      parameterization for sulphuric acid-water nucleation rates for 
    !      tropospheric and stratospheric conditions, 
    !      J. Geophys. Res., 107, D22, 4622, doi:10.1029/2002JD002184, 2002.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

    implicit none

          REAL( dp), intent(in)                :: x,t
          REAL( dp), intent(out), DIMENSION(2) :: act

! write(6,*) x,t,'in activi'
          IF (x .GT. 1.0 .OR. x .LT. 0.0)  then
        print *, x,t
        STOP 'x>1 or x<0  in activi'
          end if
          CALL  zeleznik(x,t,act)

  end subroutine activ_binary

  subroutine zeleznik(x,T,act)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Hanna Vehkamaeki
    !
    !      contact
    !      -------
    !      Dr. Hanna Vehkamaeki
    !      Professor
    !      University of Helsinki
    !      Institute for Atmospheric and Earth System Research (INAR)
    !      PL 64, FI-00014 UHEL
    !      Helsinki, Finland
    !      Phone:  +358 50415747
    !      Email:  hann.vehkamaki@helsinki.fi 
    !
    !   *
    !        (c) Hanna Vehkamaeki, 2002
    !   *
    !
    !      purpose
    !      -------
    !      Water and sulfuric acid activities in liquid
    !      aqueous solutions.
    !
    !      interface
    !      ---------
    !
    !        input:
    !           x = mole fraction of h2so4
    !           t = temperature in K
    !
    !      method
    !      ------
    !
    !
    !      reference
    !      ---------
    !      Frank J. Zeleznik, Thermodynnamic properties
    !      of the aqueous sulfuric acid system to 220K-350K,
    !      mole fraction 0,...,1
    !      J. Phys. Chem. Ref. Data, Vol. 20, No. 6,pp.1157, 1991
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

    implicit none

        REAL( dp), intent(in)                :: x,t
        REAL( dp), intent(out), DIMENSION(2) :: act

          act(2) = activity_acid(x,t)
          act(1) = activity_water(x,t)

  end subroutine zeleznik


    !
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Hanna Vehkamaeki
    !
    !      contact
    !      -------
    !      Dr. Hanna Vehkamaeki
    !      Professor
    !      University of Helsinki
    !      Institute for Atmospheric and Earth System Research (INAR)
    !      PL 64, FI-00014 UHEL
    !      Helsinki, Finland
    !      Phone:  +358 50415747
    !      Email:  hann.vehkamaki@helsinki.fi 
    !
    !   *
    !        (c) Hanna Vehkamaeki, 2002
    !   *
    !
    !      purpose
    !      -------
    !      Functions activity_acid and activity_water
    !      and all related functions:
    !
    !      interface
    !      ---------
    !
    !        input:
    !
    !
    !      method
    !      ------
    !
    !
    !      reference
    !      ---------
    !      none

    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

  ELEMENTAL REAL(dp) function activity_acid(xal,T)
    !
     implicit none

        REAL( dp), intent(in)     :: xal, T
        REAL( dp)                 :: aaa

        aaa=1.E0_dp - 1.E-7_Dp
        activity_acid = EXP(lnAa(xal,T)-lnAa(aaa,T))

    END FUNCTION activity_acid

!------------------------------------------------------------------

  ELEMENTAL REAL(dp) function activity_water(xal,T)
    !
     implicit none

          REAL( dp), intent(in)     :: xal, T
        REAL( dp)                 :: toka,eka,x1

      x1=xal
        eka=-(                                                          &
       (2*m121(T)+e121(T)*log(x1)+e211(T)*(log(1-x1)+1))*x1           &
       +(2*m221(T)+e221(T)*(2*log(1-x1)+1))*(1-x1)                    &
       -(m111(T)+e111(T)*(log(x1)+1))*x1*x1                           &
      -(2*m121(T)+e121(T)*(log(x1)+1)                                 &
                 +e211(T)*(log(1-x1)+1))*x1*(1-x1)                    &
             -(m221(T)+e221(T)*(log(1-x1)+1))*(1-x1)**2               &
        +x1*(2*m122(T)+e122(T)*log(x1)+e212(T)*log(1-x1))*x1*(1-x1)   &
       +x1*(1-x1)*((2*m122(T)+e122(T)*log(x1)                         &
                             +e212(T)*(log(1-x1)+1))*x1               &
                    -(6*m122(T)+e122(T)*(3*log(x1)+1)                 &
                             +e212(T)*(3*log(1-x1)+1))*(1-x1)*x1)  )

      x1=1.E-12
        toka=-(                                                         &
       (2*m121(T)+e121(T)*log(x1)+e211(T)*(log(1-x1)+1))*x1           &
       +(2*m221(T)+e221(T)*(2*log(1-x1)+1))*(1-x1)                    &
       -(m111(T)+e111(T)*(log(x1)+1))*x1*x1                           &
      -(2*m121(T)+e121(T)*(log(x1)+1)                                 &
                 +e211(T)*(log(1-x1)+1))*x1*(1-x1)                    &
             -(m221(T)+e221(T)*(log(1-x1)+1))*(1-x1)**2               &
        +x1*(2*m122(T)+e122(T)*log(x1)+e212(T)*log(1-x1))*x1*(1-x1)   &
       +x1*(1-x1)*((2*m122(T)+e122(T)*log(x1)                         &
                             +e212(T)*(log(1-x1)+1))*x1               &
                    -(6*m122(T)+e122(T)*(3*log(x1)+1)                 &
                             +e212(T)*(3*log(1-x1)+1))*(1-x1)*x1)   )

          activity_water = EXP(eka-toka)

    END FUNCTION  activity_water

!start of functions related to zeleznik activities

          ELEMENTAL REAL(dp) FUNCTION lnAa(x1,T)
          IMPLICIT NONE
         REAL( dp), intent(in)   :: x1,T

          lnAa=-(                                                        &
           (2*m111(T)+e111(T)*(2*log(x1)+1))*x1                        &
         +(2*m121(T)+e211(T)*log(1-x1)+e121(T)*(log(x1)+1))*(1-x1)     &
         -(m111(T)+e111(T)*(log(x1)+1))*x1*x1                          &
         -(2*m121(T)+e121(T)*(log(x1)+1)+e211(T)*(log(1-x1)+1)         &
         -(2*m122(T)+e122(T)*log(x1)                                   &
                +e212(T)*log(1-x1))*(1-x1))*x1*(1-x1)                  &
         -(m221(T)+e221(T)*(log(1-x1)+1))*(1-x1)**2                    &
         -x1*(1-x1)*(                                                  &
                       (6*m122(T)+e122(T)*(3*log(x1)+1)                &
                               +e212(T)*(3*log(1-x1)+1)                &
                        )*x1*(1-x1)                                    &
                     -(2*m122(T)+e122(T)*(log(x1)+1)                   &
                                        +e212(T)*log(1-x1)             &
                         )*(1-x1))          )
          END FUNCTION lnAa


          ELEMENTAL REAL(dp) FUNCTION lnAw(x1,T)
          IMPLICIT NONE
         REAL( dp), intent(in)   :: x1,T

          lnAw=-(                                                        &
       (2*m121(T)+e121(T)*log(x1)+e211(T)*(log(1-x1)+1))*x1            &
       +(2*m221(T)+e221(T)*(2*log(1-x1)+1))*(1-x1)                     &
       -(m111(T)+e111(T)*(log(x1)+1))*x1*x1                            &
      -(2*m121(T)+e121(T)*(log(x1)+1)                                  &
                 +e211(T)*(log(1-x1)+1))*x1*(1-x1)                     &
             -(m221(T)+e221(T)*(log(1-x1)+1))*(1-x1)**2                &
        +x1*(2*m122(T)+e122(T)*log(x1)+e212(T)*log(1-x1))*x1*(1-x1)    &
       +x1*(1-x1)*((2*m122(T)+e122(T)*log(x1)                          &
                             +e212(T)*(log(1-x1)+1))*x1                &
                    -(6*m122(T)+e122(T)*(3*log(x1)+1)                  &
                             +e212(T)*(3*log(1-x1)+1))*(1-x1)*x1)  )
          END FUNCTION lnAw

          ELEMENTAL REAL(dp) FUNCTION m111(T)
          IMPLICIT NONE
         REAL( dp), intent(in)   :: T
          m111=-23.524503387D0+0.0406889449841D0*T-0.151369362907D-4     &
                    *T**2+2961.44445015D0/T+0.492476973663D0*log(T)
          END FUNCTION m111

          ELEMENTAL REAL(dp) FUNCTION m121(T)
          IMPLICIT NONE
         REAL( dp), intent(in)   :: T
          m121=1114.58541077D0-1.1833078936D0*T-0.00209946114412D0*T**2  &
                     -246749.842271D0/T+34.1234558134D0*log(T)
          END FUNCTION m121

          ELEMENTAL REAL(dp) FUNCTION m221(T)
          IMPLICIT NONE
         REAL( dp), intent(in)   :: T
          m221=-80.1488100747D0-0.0116246143257D0*T                      &
           +0.606767928954D-5*T**2                                     &
                 +3092.72150882D0/T+12.7601667471D0*log(T)
          END FUNCTION m221

          ELEMENTAL REAL(dp) FUNCTION m122(T)
          IMPLICIT NONE
         REAL( dp), intent(in)   :: T
          m122=888.711613784D0-2.50531359687D0*T                         &
           +0.000605638824061D0*T**2                                   &
                     -196985.296431D0/T+74.550064338D0*log(T)
          END FUNCTION m122

          ELEMENTAL REAL(dp) FUNCTION e111(T)
          IMPLICIT NONE
         REAL( dp), intent(in)   :: T
          e111=2887.31663295D0-3.32602457749D0*T                         &
           -0.2820472833D-2*T**2                                       &
                     -528216.112353D0/T+0.68699743564D0*log(T)
          END FUNCTION e111

          ELEMENTAL REAL(dp) FUNCTION e121(T)
          IMPLICIT NONE
         REAL( dp), intent(in)   :: T
          e121=-370.944593249D0-0.690310834523D0*T                      &
           +0.56345508422D-3*T**2                                     &
                     -3822.52997064D0/T+94.2682037574D0*log(T)
          END FUNCTION e121

          ELEMENTAL REAL(dp) FUNCTION e211(T)
          IMPLICIT NONE
         REAL( dp), intent(in)   :: T
          e211=38.3025318809D0-0.0295997878789D0*T                      &
           +0.120999746782D-4*T**2                                    &
                     -3246.97498999D0/T-3.83566039532D0*log(T)
          END FUNCTION e211

          ELEMENTAL REAL(dp) FUNCTION e221(T)
          IMPLICIT NONE
         REAL( dp), intent(in)   :: T
          e221=2324.76399402D0-0.141626921317D0*T                       &
           -0.00626760562881D0*T**2                                   &
                     -450590.687961D0/T-61.2339472744D0*log(T)
          END FUNCTION e221

          ELEMENTAL REAL(dp) FUNCTION e122(T)
          IMPLICIT NONE
         REAL( dp), intent(in)   :: T
          e122=-1633.85547832D0-3.35344369968D0*T                       &
           +0.00710978119903D0*T**2                                   &
                     +198200.003569D0/T+246.693619189D0*log(T)
          END FUNCTION e122

          ELEMENTAL REAL(dp) FUNCTION e212(T)
          IMPLICIT NONE
         REAL( dp), intent(in)   :: T
          e212=1273.75159848D0+1.03333898148D0*T                        &
           +0.00341400487633D0*T**2                                   &
                     +195290.667051D0/T-431.737442782D0*log(T)
          END FUNCTION e212

!------------------------------------------------------------------


  end module gde_nucleation
