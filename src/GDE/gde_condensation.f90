! <gde_condensation.f90 - A component of the Multicomponent
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
!*    routine CONDENSATION_COEFF2 written by Liisa Pirjola and Matthias Karl
!* 
!*****************************************************************************!
module gde_condensation

    use messy_mecca_kpp_global,    only : APN,dp

    use gde_constants,  only      : MB,MIO,MAN,MAH,MVOC,MAC,MNH
    use gde_constants,  only      : M_H2SO4,M_msa,M_hio3
    use gde_constants,  only      : M_nit,M_ca,M_nh3,M_hcl
    use gde_constants,  only      : M_H2O,M_air,MC
    use gde_constants,  only      : pi,k_B,N_A,R_gas,T0,RHOH2O

    use gde_input_data, only      : MMAX,QMAX,AMAX
    use gde_input_data, only      : NU,AI,AS,CS
    use gde_input_data, only      : aqmax
    use gde_input_data, only      : NSOA
    use gde_input_data, only      : DENV,DENMS,DENIO
    use gde_input_data, only      : DENXX,DENNI,DENAM
    use gde_input_data, only      : DENEC   
    use gde_input_data, only      : A_SUL,A_MSA,A_IO3
    use gde_input_data, only      : A_NIT,A_AMI,A_NH4
    use gde_input_data, only      : A_OR1,A_OR2,A_OR3,A_OR4,A_OR5
    use gde_input_data, only      : A_OR6,A_OR7,A_OR8,A_OR9
    use gde_input_data, only      : A_XXX,A_SAL,A_CHL
    use gde_input_data, only      : CONVM,massmin,MVAP,DIFV
    use gde_input_data, only      : nucomin
    use gde_sensitiv,   only      : ICONS,ICONA,ICONO,ICONX,ISOA 
    use gde_sensitiv,   only      : ICONW
    use gde_toolbox,    only      : molec2ug
    use gde_toolbox,    only      : molecdiff
    use gde_toolbox,    only      : waterps
    !use gde_toolbox,    only      : acidps ! not used
    use gde_toolbox,    only      : satps_sulf
    use gde_toolbox,    only      : sulfhydrates

    private
   
    public :: condensloss
    public :: condensation
    public :: apc_update_pmass
    public :: apc_update_gasc
    public :: condensation_incloud
    public :: condensation_coeff2
    public :: nitcondens

  contains
  
    subroutine condensloss(temp,press,RH,DPA,N,IMAX,MASS,alphanit,fcoj,KPEQ,                    &
                           DENOC,M_oc,nmo,foc,hvap,csat0,                                       &
                           ctso4,ctnh4,ctnit,ctchl,NVAP,                                        &
                           MORG1TOT,MORG2TOT,MORG3TOT,MORG4TOT,MORG5TOT,                        &
                           MORG6TOT,MORG7TOT,MORG8TOT,MORG9TOT,                                 &
                           MORGCTOT,MECBCTOT,MDUSTTOT,MTOT,  keffect,                           &
                           henry_eff_hno3,henry_eff_hcl,svmc_min_hno3,svmc_min_hcl,             & 
                           svmc_min_nh3,flag_dissolution,Keq_nh4no3_0,Keq_nh4cl_0,              &
                           sumlwc,DTIME,GRSU,GRMS,DCSU,GROC,DCORG,CCOND,                        &
                           NSV,EXCESS,LOSS,TRANS,csate)

                 
    !----------------------------------------------------------------------
    !
    !****  controls calculation of multicomponent condensation
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !
    !
    !      interface
    !      ---------
    !
    !        input:
    !           press [Pa]
    !           temp  [K]
    !           mbh   [m]
    !           denspar [kg/m3]
    !    
    !      method
    !      ------
    !        1) calculate saturation vapour density
    !        2) calculate SOA partitioning
    !        3) calculate excess and loss
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
    
    INTEGER, intent(in)                               :: IMAX
    real( dp), intent(in)                             :: temp,press,RH
    real( dp), dimension(MMAX), intent(in)            :: MORG1TOT,MORG2TOT,MORG3TOT
    real( dp), dimension(MMAX), intent(in)            :: MORG4TOT,MORG5TOT,MORG6TOT
    real( dp), dimension(MMAX), intent(in)            :: MORG7TOT,MORG8TOT,MORG9TOT
    real( dp), dimension(MMAX), intent(in)            :: MORGCTOT,MECBCTOT,MDUSTTOT
    real( dp), dimension(MMAX), intent(in)            :: MTOT
    real( dp), DIMENSION(QMAX), intent(in)            :: NVAP
    
    real( dp), dimension(MMAX,IMAX,QMAX),intent(in)   :: MASS, keffect
    real( dp), dimension(MMAX,IMAX), intent(in)       :: N
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPA
    ! (chamber)     
    real( dp),intent(in)                              :: alphanit,fcoj,KPEQ
    ! (organics)
    real( dp),intent(in)                              :: DENOC
    real( dp), dimension(NSOA), intent(in)            :: M_oc
    real( dp), dimension(NSOA), intent(in)            :: nmo
    real( dp), dimension(NSOA), intent(in)            :: foc
    real( dp), dimension(NSOA), intent(in)            :: hvap
    real( dp), dimension(NSOA), intent(in)            :: csat0
    ! (SIA)
    real( dp),intent(in)                              :: ctso4,ctnh4,ctnit,ctchl

! mesa start
    real( dp),dimension(MMAX,IMAX),intent(in)         :: henry_eff_hno3
    real( dp),dimension(MMAX,IMAX),intent(in)         :: henry_eff_hcl
    real( dp),dimension(MMAX,IMAX),intent(in)         :: svmc_min_hno3
    real( dp),dimension(MMAX,IMAX),intent(in)         :: svmc_min_hcl
    real( dp),dimension(MMAX,IMAX),intent(in)         :: svmc_min_nh3
    integer, dimension(MMAX,IMAX),intent(in)          :: flag_dissolution
    real( dp),intent(in)                              :: Keq_nh4no3_0
    real( dp),intent(in)                              :: Keq_nh4cl_0
    real( dp),intent(in)                              :: sumlwc
    real( dp),intent(in)                              :: DTIME
! mesa end

    !output
    real( dp), DIMENSION(MMAX,IMAX,QMAX), intent(out) :: CCOND,EXCESS,LOSS
    real( dp), DIMENSION(MMAX,IMAX,QMAX), intent(out) :: TRANS
    real( dp), DIMENSION(QMAX), intent(out)           :: NSV
    real( dp), intent(out)                            :: GRSU,GRMS,DCSU,GROC,DCORG
    real( dp), dimension(NSOA),intent(out)            :: csate

! local
    real( dp), dimension(MMAX,IMAX)     :: CCONDSUL,CCONDMSA,CCONDIOD
    real( dp), dimension(MMAX,IMAX)     :: CCONDNIT,CCONDCHL
    real( dp), dimension(MMAX,IMAX)     :: CCONDORG1,CCONDORG2,CCONDORG3
    real( dp), dimension(MMAX,IMAX)     :: CCONDORG4,CCONDORG5,CCONDORG6
    real( dp), dimension(MMAX,IMAX)     :: CCONDORG7,CCONDORG8,CCONDORG9

! soa parameters    
    real( dp), dimension(NSOA)          :: mwsoa
    real( dp), dimension(NSOA)          :: msize
    real( dp), dimension(NSOA)          :: ocfrac
    real( dp), dimension(NSOA)          :: ps
    real( dp)                           :: cmnon1,cmnon2,cmnon3
    real( dp)                           :: pmtot
    real( dp)                           :: pptot
    real( dp)                           :: fom
    real( dp)                           :: soa1tot,soa2tot
    real( dp)                           :: soa4tot,soa5tot
    real( dp)                           :: soa7tot,soa8tot
    real( dp)                           :: cs1,cs2
    real( dp)                           :: cs3,cs4
    real( dp)                           :: cs5,cs6

    real( dp)                           :: dcorg1,dcorg2,dcorg3,dcorg4,dcorg5
    real( dp)                           :: dcorg6,dcorg7,dcorg8,dcorg9
    real( dp)                           :: groc1,groc2,groc3,groc4,groc5 
    real( dp)                           :: groc6,groc7,groc8,groc9

! nh4-hno3-hcl
    real( dp)                           :: KPEQNH4
    real( dp)                           :: KPEQNCL
    real( dp)                           :: KPNIT
    real( dp)                           :: KPNCL
    real( dp)                           :: psnh4,psno3,pscl
    real( dp)                           :: nh3free,hno3tot
    real( dp)                           :: Ttrans
    real( dp), parameter                :: psmin_nh3  = 2.46e09_dp

    INTEGER                             :: M,I,Q,S



! Calculate saturation vapour concentrations in [molec/m^3]
! Compute condensation coefficients 

         ! 1) ammonium nitrate concentration. Use Kp in (molec/cm3)^2
         !KPEQNH4=6.84E21_dp  ! NH4NO3 at 298 K
         call eqnh4nitrate(NVAP(A_NIT)*1.e-6,NVAP(A_NH4)*1.e-6,    &
                           NVAP(A_CHL)*1.e-6, temp, KPEQNH4,       &
                           KPEQNCL, psno3,pscl,psnh4 )

         call nitcondens(IMAX,press,temp,DPA,alphanit,fcoj,CCONDNIT)

         if (iconw.eq.2) then
           KPNIT = Keq_nh4no3_0   ! from MESA
           KPNCL = Keq_nh4cl_0    ! from MESA
         else
           KPNIT = KPEQNH4
           KPNCL = KPEQNCL
         endif


         ! 2) alkylammonium nitrate concentration. Use Kp in (molec/cm3)^2
         if (ICONA == 1) then
             call eqnitrate(NVAP(A_NIT)*1.e-6,NVAP(A_AMI)*1.e-6,KPEQ,nsv(A_NIT),nsv(A_AMI))
             call nitcondens(IMAX,press,temp,DPA,alphanit,fcoj,CCONDNIT)
             KPNIT = KPEQ
         endif

         ! 3) sulphuric acid, methane sulphuric acid, iodic acid
         ! MSK 29.09.2013 New approximation of nsv(A_SUL) for aqueous h2s4o solution
         !    P. Bolsaitis and J. F. Elliott			
         !    Thermodynamic Activities and equilibrium partial 			
         !    pressures for aqueous sulfuric acid solutions			
         !    J Phys Chem Eng Data 1990, 35, 69-85.
         ! MSK 03.10.2022 Transition temperature below which MSA behaves like a ELVOC
         !    T_trans = a - b*RH + c*RH^2 - d*RH^3 + e*RH^4
         !    Hodshire, A. L., Campuzano-Jost, P., Kodros, J. K., 
         !    Croft, B., Nault, B. A., Schroder, J. C., 
         !    Jimenez, J. L., and Pierce, J. R.
         !    The potential role of methanesulfonic acid (MSA) 
         !    in aerosol formation and growth and the associated 
         !    radiative forcings
         !    Atmos Chem Phys 2019, 19, 3137-3160, 
         !    doi:10.5194/acp-19-3137-2019.

         ! transition temperature (in K)
         Ttrans = 252.0 - 6.19e-1_dp*RH*100._dp + 3.49e-2_dp*(RH*100._dp)**2 - &
                  5.6e-4_dp*(RH*100._dp)**3 + 3.32e-6_dp*(RH*100._dp)**4

         !SUL
         !old nsv(A_SUL)=acidps(temp) !Pa
         !old nsv(A_SUL)=nsv(A_SUL)/(k_B*1.e6*temp)        !molec/m^3    (factor 1E6 scaled)
         nsv(A_SUL) = satps_sulf(temp,RH,NVAP(A_SUL))  !Pa
         nsv(A_SUL)=nsv(A_SUL)/(k_B*temp)              !molec/m^3

         !MSA
         if (temp >= Ttrans) then
         ! MSA behaves as SVOC
           nsv(A_MSA)=(-8.00648e3_dp/temp)+2.14237_dp*LOG(temp)+7.45208_dp     !mmHg
           nsv(A_MSA)=EXP(nsv(A_MSA))                                          !mmHg
           nsv(A_MSA)=nsv(A_MSA)*(1.e5_dp/750.06_dp)                           !Pa
           nsv(A_MSA)=nsv(A_MSA)*1.e-6_dp              !low volatile vapor, scaled by 1.e-6    
         else
         ! MSA behaves as ELVOC
           nsv(A_MSA)=1.0e-11_dp                       !Pa
         endif
         nsv(A_MSA)=nsv(A_MSA)/(k_B*temp)              !molec/m^3

         !HIO3
         nsv(A_IO3)=2.2e-9_dp                          !Pa, EPISuite MPBPWIN v1.42
         nsv(A_IO3)=nsv(A_IO3)*1.e-2_dp                !EL vapor, scaled by 0.01
         nsv(A_IO3)=nsv(A_IO3)/(k_B*temp)              !molec/m^3

         call condensation_coeff2(temp,RH,press,DPA,GRSU,GRMS,DCSU,CCONDSUL,CCONDMSA,CCONDIOD,IMAX)


         ! 4) ammonium sulfate NH4HSO4, (NH4)2SO4 or ammonium nitrate
         !     seinfeld and pandis (1998):                               
         !     atmospheric chemistry and physics (0-471-17816-0)
         !     ctnh4,ctso4,ctnit are in ng/m^3
         !     ctnh4 < 2*ctso4 --> ammonium-poor
         !     ctnh4 > 2*ctso4 --> ammonium-rich
         !     if c(NH3_free)*c(HNO3) > Kp --> NH4NO3 condensation
         !
         !nh3free = ctnh4 - (ctso4*2.*M_nh3/M_H2SO4)
         !nh3free = max(0.0,nh3free)
         !write(6,*) 'ctnh4 ctso4 ctnit',ctnh4,ctso4,ctnit
         !nh3free = nh3free*1e-3/molec2ug(M_nh3)   !molec/cm^3
         !hno3tot = ctnit*1e-3/molec2ug(M_nit)     !molec/cm^3


         !NH3
         if ( ctnh4 < 2.0*ctso4 ) then
         ! 2 NH3(g) + H2SO4 --> (NH4)2SO4 is one-way
         ! NH4SO4, p0(NH3) not below psmin_nh3
           nsv(a_nh4) = max(nsv(a_sul),psmin_nh3)
           !print *,'pNH3 1',nsv(a_nh4),ctnh4,ctso4,ctnh4/ctso4
         else
           if ( NVAP(A_NIT)*NVAP(A_NH4)*1.e-12 > KPNIT ) then
             nsv(a_nh4) = psnh4
             !print *,'pNH3 2',nsv(a_nh4),ctnh4,ctso4,ctnh4/ctso4
           else if ( NVAP(A_CHL)*NVAP(A_NH4)*1.e-12 > KPNCL ) then
             !nsv(a_nh4) = psnh4
         ! TEST MSK 24.05.2023 condensation of NH4Cl
             nsv(a_nh4) = psnh4*0.995_dp
             !print *,'pNH3 3',nsv(a_nh4),ctnh4,ctso4,ctnh4/ctso4
           else
             nsv(a_nh4) = max(nsv(a_sul),psmin_nh3)
             !print *,'pNH3 4',nsv(a_nh4),ctnh4,ctso4,ctnh4/ctso4
           endif
         endif
      ! write(6,'(6ES12.4)') KPNIT, Keq_nh4no3_0,NVAP(A_NIT), NVAP(A_NH4), NVAP(A_NIT)*NVAP(A_NH4)*1.e-12

         !5) hydrochloric acid 
         call chlcondens(IMAX,press,temp,DPA,CCONDCHL)


         !6) Sat. pressure SOA-1, ..., SOA-9
         !   sat. pressure in Pa at given temperature T
         !   psat0(298K) = (C0(298K)*R_gas*298)/(1.e6*M_oc)
         !   psat0(T) = psat0(298K)*EXP((Hvap/R_gas)*((1/298K)-(1/temp)))
         ! csate is the effective saturation concentration C* (here C0)
         ! do for all SOA components (A_OR1, ..., A_OR9)
         ! S+A_OR1-1 is the index of SOA components
         do S=1,NSOA
           nsv(S+A_OR1-1)=(csat0(S)*R_gas*298._dp)/(1.e6_dp*M_oc(S))
           nsv(S+A_OR1-1)=nsv(S+A_OR1-1)*EXP((hvap(S)/R_gas)*  &
                          ((1._dp/298._dp)-(1._dp/temp)))
           !effective saturation concentration C*(T), assume gamma=1
           csate(S)=csat0(S)*EXP((hvap(S)/R_gas)*((1._dp/298._dp)-(1._dp/temp)))
         enddo
         !print *,"nsv(T) Pa ",nsv(A_OR1),nsv(A_OR2),nsv(A_OR3),nsv(A_OR4),nsv(A_OR5)
         !print *,"nsv(T) Pa ",nsv(A_OR6),nsv(A_OR7),nsv(A_OR8),nsv(A_OR9)


         !!! BELOW nsv is converted to sat. concentration in molec/m^3  !!!

         !   SOA PARTITIONING         
         if ((ISOA == 1).or.(ISOA==2)) then
             ! non-volatile SOC = 0.0 (BELV,AELV,PELV)
               cmnon1 = MORG3TOT(NU)+MORG3TOT(AI)+MORG3TOT(AS)+MORG3TOT(CS)
               cmnon2 = MORG6TOT(NU)+MORG6TOT(AI)+MORG6TOT(AS)+MORG6TOT(CS)
               cmnon3 = MORG9TOT(NU)+MORG9TOT(AI)+MORG9TOT(AS)+MORG9TOT(CS)
             ! fraction of absorptive material (organics) in total PM
               pmtot= max( (MTOT(NU)+MTOT(AI)+MTOT(AS)+MTOT(CS)),massmin )
               fom  = (MORGCTOT(NU)+MORGCTOT(AI)+MORGCTOT(AS)+MORGCTOT(CS)) / pmtot
               fom  = max(fom,0.3_dp)  ! at least 30% OM
               fom  = min(fom,1.0_dp)
             ! total SOA aerosol concentration
             ! S+A_OR1-1 is the index of SOA components
               do S=1,NSOA
                 mwsoa(S)  = M_oc(S)    ! g/mol
                 msize(S)  = nmo(S)     ! (nC+nO)
                 ocfrac(S) = foc(S)     ! (O:C ratio, carbon fraction)
                 ps(S)     = nsv(S+A_OR1-1)   ! Pa
               enddo
             ! biogenic SOA
               soa1tot = MORG1TOT(NU)+MORG1TOT(AI)+MORG1TOT(AS)+MORG1TOT(CS)
               soa2tot = MORG2TOT(NU)+MORG2TOT(AI)+MORG2TOT(AS)+MORG2TOT(CS)
             ! aromatic SOA
               soa4tot = MORG4TOT(NU)+MORG4TOT(AI)+MORG4TOT(AS)+MORG4TOT(CS)
               soa5tot = MORG5TOT(NU)+MORG5TOT(AI)+MORG5TOT(AS)+MORG5TOT(CS)
             ! nalkene SOA
               soa7tot = MORG7TOT(NU)+MORG7TOT(AI)+MORG7TOT(AS)+MORG7TOT(CS)
               soa8tot = MORG8TOT(NU)+MORG8TOT(AI)+MORG8TOT(AS)+MORG8TOT(CS)
             ! pre-existing adsorptive material
               pptot   = MECBCTOT(NU)+MECBCTOT(AI)+MECBCTOT(AS)+MECBCTOT(CS) &
                        +MDUSTTOT(NU)+MDUSTTOT(AI)+MDUSTTOT(AS)+MDUSTTOT(CS)

               call soapartition(temp,soa1tot,soa2tot,soa4tot,soa5tot,  &
                                 soa7tot,soa8tot                     ,  & 
                                 cmnon1,cmnon2,cmnon3                ,  &
                                 fom,ps,mwsoa,msize,ocfrac           ,  &
                                 pptot,cs1,cs2,cs3,cs4,cs5,cs6)

             ! saturation concentration in molec/m^3
             ! biogenic SOA
               NSV(A_OR1)=cs1
               NSV(A_OR2)=cs2
             ! aromatic SOA
               NSV(A_OR4)=cs3
               NSV(A_OR5)=cs4
             ! nalkene SOA
               NSV(A_OR7)=cs5
               NSV(A_OR8)=cs6
             ! non-volatile SOC (BELV,AELV,PELV)
             ! csat(T) in ug m^-3 --> nsv(T) in molec m^-3
             ! BELV (SOA-3)
               NSV(A_OR3)=csat0(3)*EXP((hvap(3)/R_gas)*((1._dp/298._dp)-(1._dp/temp)))
               NSV(A_OR3)=NSV(A_OR3)*1.e6/molec2ug(M_oc(3))
             ! AELV (SOA-6)
               NSV(A_OR6)=csat0(6)*EXP((hvap(6)/R_gas)*((1._dp/298._dp)-(1._dp/temp)))
               NSV(A_OR6)=NSV(A_OR6)*1.e6/molec2ug(M_oc(6))
             ! PELV (SOA-9)     
               NSV(A_OR9)=csat0(9)*EXP((hvap(9)/R_gas)*((1._dp/298._dp)-(1._dp/temp)))
               NSV(A_OR9)=NSV(A_OR9)*1.e6/molec2ug(M_oc(9))
             ! effective saturation concentration C*(T) [Donahue et al., 2006]
             ! (for elvoc we use C0) in ug/m^3
               csate(1)=cs1*molec2ug(M_oc(1))*1.e-6 
               csate(2)=cs2*molec2ug(M_oc(2))*1.e-6
               csate(4)=cs3*molec2ug(M_oc(4))*1.e-6
               csate(5)=cs4*molec2ug(M_oc(5))*1.e-6
               csate(7)=cs5*molec2ug(M_oc(7))*1.e-6
               csate(8)=cs6*molec2ug(M_oc(8))*1.e-6

         else
             ! saturation concentration in molec/m^3
               NSV(A_OR1)=NSV(A_OR1)/(k_B*temp) 
               NSV(A_OR2)=NSV(A_OR2)/(k_B*temp)
               NSV(A_OR3)=NSV(A_OR3)/(k_B*temp) 
               NSV(A_OR4)=NSV(A_OR4)/(k_B*temp)
               NSV(A_OR5)=NSV(A_OR5)/(k_B*temp)
               NSV(A_OR6)=NSV(A_OR6)/(k_B*temp)
               NSV(A_OR7)=NSV(A_OR7)/(k_B*temp)
               NSV(A_OR8)=NSV(A_OR8)/(k_B*temp)
               NSV(A_OR9)=NSV(A_OR9)/(k_B*temp)
         endif

        ! print *,"nsv(T) molec/m^3 ",nsv(A_OR1),nsv(A_OR2),nsv(A_OR3)
        ! print *,"nsv(T) molec/m^3 ",nsv(A_OR4),nsv(A_OR5),nsv(A_OR6)
        ! print *,"nsv(T) molec/m^3 ",nsv(A_OR7),nsv(A_OR8),nsv(A_OR9)


! condensation coefficient of organic vapors
         call orgcondens(press,temp,DPA,M_oc(1),DENOC,groc1,dcorg1,CCONDORG1,IMAX)
         call orgcondens(press,temp,DPA,M_oc(2),DENOC,groc2,dcorg2,CCONDORG2,IMAX)
         call orgcondens(press,temp,DPA,M_oc(3),DENOC,groc3,dcorg3,CCONDORG3,IMAX)
         call orgcondens(press,temp,DPA,M_oc(4),DENOC,groc4,dcorg4,CCONDORG4,IMAX)   
         call orgcondens(press,temp,DPA,M_oc(5),DENOC,groc5,dcorg5,CCONDORG5,IMAX)
         call orgcondens(press,temp,DPA,M_oc(6),DENOC,groc6,dcorg6,CCONDORG6,IMAX)
         call orgcondens(press,temp,DPA,M_oc(7),DENOC,groc7,dcorg7,CCONDORG7,IMAX)
         call orgcondens(press,temp,DPA,M_oc(8),DENOC,groc8,dcorg8,CCONDORG8,IMAX)
         call orgcondens(press,temp,DPA,M_oc(9),DENOC,groc9,dcorg9,CCONDORG9,IMAX)

! growth rate by organic vapors
         GROC=max(groc1,groc2)
         DCORG=(dcorg1+dcorg2+dcorg3+dcorg4+dcorg5+dcorg6+dcorg7+dcorg8+dcorg9)/9


! Initialize component-specific condensation rate
         do M=NU,CS
          do I=1,IMAX
            ccond(M,I,A_SUL)=CCONDSUL(M,I)
            ccond(M,I,A_MSA)=CCONDMSA(M,I)
            ccond(M,I,A_IO3)=CCONDIOD(M,I)
            ccond(M,I,A_AMI)=CCONDNIT(M,I)
!A_CHL
            if (ctchl > ctnit) then
              ccond(M,I,a_chl)=CCONDCHL(M,I)           !low-NO3
            else
              if (ctnh4 > 0.8*ctso4) then
             ! aerosol neutralized by ammonia
                  ccond(M,I,A_CHL)=CCONDCHL(M,I)*5.0   !high-N  (case 3)
              else
                if ( mass(m,i,a_chl).lt.mass(m,i,a_nit) ) then
                  ccond(M,I,A_CHL)=CCONDCHL(M,I)*40.0  !low-NH4 (case 3)
                else
                  ccond(M,I,A_CHL)=CCONDCHL(M,I)
                endif
              endif
            endif
!A_NO3
            if (( NVAP(A_NIT)*NVAP(A_NH4)*1.e-12 > KPNIT ) .or. &
                ( NVAP(A_NIT)*NVAP(A_AMI)*1.e-12 > KPNIT ) ) then
            !NH4NO3 condensation (case 7)
               ccond(m,i,a_nit)=ccondnit(m,i)
               !if((m==3).and.(i==2)) print *,'no3co1',KPNIT,  &
               !NVAP(A_NIT)*NVAP(A_NH4)*1.e-12,ccond(m,i,a_nit)
            else
              if (mass(m,i,a_chl).gt.10._dp) then
            !CL REPLACEMENT
                if (ctnh4 > 0.8*ctso4) then
            !high-N  (case 3)
                  ccond(m,i,a_nit)=ccondnit(m,i)*3.e-5
               !if((m==3).and.(i==2)) print *,'no3co2',KPNIT,  &
               !NVAP(A_NIT)*NVAP(A_NH4)*1.e-12,ccond(m,i,a_nit)
                else
            !low-NH4 (case 3)
                  ccond(m,i,a_nit)=ccondnit(m,i)*2.e-3
               !if((m==3).and.(i==2)) print *,'no3co3',KPNIT,   &
               !NVAP(A_NIT)*NVAP(A_NH4)*1.e-12,ccond(m,i,a_nit)
                endif
              else
            !no condensation
               ccond(m,i,a_nit)=0.0
               !if((m==3).and.(i==2)) print *,'no3co4',KPNIT,   &
               !NVAP(A_NIT)*NVAP(A_NH4)*1.e-12,ccond(m,i,a_nit)
              endif
            endif
!A_NH4
            if ( ctnh4 < 2.0*ctso4 ) then
              if ( ctnh4 > 1.5*ctso4 ) then
                ccond(m,i,a_nh4)=2*ccondsul(m,i)  !(NH4)2SO4
              else
                ccond(m,i,a_nh4)=ccondsul(m,i)    !NH4HSO4
              endif
            else
              if ( NVAP(A_NIT)*NVAP(A_NH4)*1.e-12 > KPNIT ) then
                ccond(m,i,a_nh4)=ccondnit(m,i)    !NH4NO3
              else if ( NVAP(A_CHL)*NVAP(A_NH4)*1.e-12 > KPNCL ) then
                ccond(m,i,a_nh4)=ccondnit(m,i)    !NH4Cl         
              else
                ccond(m,i,a_nh4)=0.0              !no condensation
              endif
            endif
            ccond(M,I,A_OR1)=CCONDORG1(M,I)
            ccond(M,I,A_OR2)=CCONDORG2(M,I)
            ccond(M,I,A_OR3)=CCONDORG3(M,I)
            ccond(M,I,A_OR4)=CCONDORG4(M,I)
            ccond(M,I,A_OR5)=CCONDORG5(M,I)
            ccond(M,I,A_OR6)=CCONDORG6(M,I)
            ccond(M,I,A_OR7)=CCONDORG7(M,I)
            ccond(M,I,A_OR8)=CCONDORG8(M,I)
            ccond(M,I,A_OR9)=CCONDORG9(M,I)
          end do
         end do

        ! if (iconw.eq.2) then
        !   write(6,*) 'KPnit',NVAP(A_NIT)*NVAP(A_NH4)*1.e-12,KPNIT,ccond(3,3,a_nit)
        !   write(6,*) 'sum LWC (ng/m3)',sumlwc
        ! endif


!--------------------------------------------------------------------- 
! compute gas-phase/particle condensation loss
!
! MESA-coupling: use SVMC from MESA for a_nh4 and a_nit
!    according to in Jacobson (2002); Eq. (16) 
!    M.Z. Jacobson,  JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 107, 
!          NO. D19, 4366, doi:10.1029/2001JD002044, 2002.
!
!    M.Z. Jacobson,  AEROSOL SCIENCE TECHNOLOGY, VOL. 39, 
!          NO. 2, 92-103, DOI: 10.1080/027868290904546, 2005.
!
! Solution for dissolutional growth for HNO3
! Analytical Predictor for Dissolution (APD)
! update with
! Predictor of Nonequilibrium Growth (PNG)
!
! EXCESS:  Equation (9) and (10) in Jacobson (2005)
!          and also Equation (11) in Jacobson (2002).
! LOSS:    Numerator of Eq. (12) in Jacobson (2005)
!          and also the numerator of Eq. (16) in Jacobson (2002).
!          Note the negative sign in the EXP term
! TRANS    Denominator of Eq. (12) in Jacobson (2005)
!          and also denominator of Eq. (16) in Jacobson (2002).
!          Note the negative sign in the EXP term
!
! conversion factor 1.e3 for mass: ng/m^3 --> mlc/m^3
!
! Apply ICONW2 only if LWC above low-LWC limit
! LWC limit in PNG solver (Jacobson, 2005): 0.01 ug/m^3
! LWC is sum of aerosol water in all bins = sumlwc
! sumlwc must be > 10 ng/m^3
!
!---------------------------------------------------------------------

         do m=1,mmax
           do i=1,imax
             if (iconw.eq.2) then
               if ( ctnh4 < 1.0*ctso4 ) then
              !!! if ( ctnh4 < 2.0*ctso4 ) then
                 nsv(a_nh4) = nsv(a_sul)
               else
                 nsv(a_nh4) = svmc_min_nh3(m,i)       !from MESA
                 nsv(a_nh4) = nsv(a_nh4) *0.05_dp     !tuning 
               endif
               nsv(a_nit) = svmc_min_hno3(m,i)        !from MESA
               nsv(a_chl) = svmc_min_hcl(m,i)         !from MESA
               nsv(a_chl) = nsv(a_chl)*10.0_dp        !tuning
             else
              ! nsv(a_nh4) is already determined above 
               nsv(a_nit) = psno3
               nsv(a_chl) = pscl
             endif
             if (nsv(a_nh4).eq.0.)then
               nsv(a_nh4) = nsv(a_sul)                  ! very low ps 
             endif

   ! debug
   !         write(6,'(2I3,4ES12.4)') m,i,nsv(a_nh4),nsv(a_nit),nsv(a_chl)
   !         write(6,'(2I3,4ES12.4)') m,i,ccond(m,i,a_nh4),ccond(m,i,a_nit),ccond(m,i,a_chl)
   !         write(6,'(2I3,4ES12.4)') m,i,svmc_min_nh3(m,i),nsv(a_nh4), ctnh4, ctso4
   ! debug

   !! APD Condensation/Evaporation
             do q=1,qmax

               !EXCESS
               excess(m,i,q) = nvap(q) - nsv(q) * keffect(m,i,q)

               !LOSS
               if ((excess(m,i,q).gt.0.) .or. (mass(m,i,q).gt.massmin)) then
                 loss(m,i,q)=n(m,i)*ccond(m,i,q)*nsv(q) * keffect(m,i,q)
               else
                 loss(m,i,q)=0.
               end if

               !TRANS
               trans(m,i,q) = n(m,i)*ccond(m,i,q)

             end do

    !! PNG Dissolution
             if ( (iconw.eq.2).and.(flag_dissolution(m,i)==1)   &
                   .and.(sumlwc.gt.10._dp) ) then

               !EXCESS
               ! HNO3 dissolution
               excess(m,i,a_nit) = nvap(a_nit)  -                      &
                   ( mass(m,i,a_nit)*1.e3*(1./molec2ug(M_nit)) *       &
                     keffect(m,i,a_nit) / henry_eff_hno3(m,i) )

               ! HCl dissolution
               excess(m,i,a_chl) = nvap(a_chl)  -                      &
                   ( mass(m,i,a_chl)*1.e3*(1./molec2ug(M_hcl)) *       &
                    keffect(m,i,a_chl) / henry_eff_hcl(m,i) )


               !LOSS                     
               ! HNO3 dissolution
               loss(m,i,a_nit) = mass(m,i,a_nit)*1.e3*(1./molec2ug(M_nit)) * &
                   ( 1. - EXP( (-1.)*((DTIME* n(m,i) * ccond(m,i,a_nit) *     &
                    keffect(m,i,a_nit) ) / henry_eff_hno3(m,i) ))  )

               ! HCl dissolution
               loss(m,i,a_chl) = mass(m,i,a_chl)*1.e3*(1./molec2ug(M_hcl)) * &
                   ( 1. - EXP( (-1.)*((DTIME* n(m,i) * ccond(m,i,a_chl) *     &
                    keffect(m,i,a_chl) ) / henry_eff_hcl(m,i) ))  )


               !TRANS
               ! HNO3 dissolution
               trans(m,i,a_nit) = (henry_eff_hno3(m,i) /keffect(m,i,a_nit)) * &
                   ( 1. - EXP( (-1.)*((DTIME* n(m,i) * ccond(m,i,a_nit) *    &
                    keffect(m,i,a_nit) ) / henry_eff_hno3(m,i) ))  )

               ! HCl dissolution
               trans(m,i,a_chl) = (henry_eff_hcl(m,i) /keffect(m,i,a_chl)) * &
                   ( 1. - EXP( (-1.)*((DTIME* n(m,i) * ccond(m,i,a_chl) *    &
                    keffect(m,i,a_chl) ) / henry_eff_hcl(m,i) ))  )
             endif
             
           !  write(6,'(2I3,4ES12.4)') m,i,excess(m,i,a_chl),loss(m,i,a_chl),trans(m,i,a_chl)

           end do
         end do


! switch off condensation/evaporation for compounds with zero ICON option 

         if (ICONS == 0) then
           excess(:,:,a_sul) = 0.0
           excess(:,:,a_msa) = 0.0
           excess(:,:,a_io3) = 0.0
           loss(:,:,a_sul  ) = 0.0
           loss(:,:,a_msa  ) = 0.0
           loss(:,:,a_io3  ) = 0.0
           trans(:,:,a_sul ) = 0.0
           trans(:,:,a_msa ) = 0.0
           trans(:,:,a_io3 ) = 0.0
         endif
         if (ICONO == 0) then
           excess(:,:,a_or1) = 0.0; excess(:,:,a_or2) = 0.0; excess(:,:,a_or3) = 0.0
           excess(:,:,a_or4) = 0.0; excess(:,:,a_or5) = 0.0; excess(:,:,a_or6) = 0.0
           excess(:,:,a_or7) = 0.0; excess(:,:,a_or8) = 0.0; excess(:,:,a_or9) = 0.0
           loss(:,:,a_or1)   = 0.0; loss(:,:,a_or2)   = 0.0; loss(:,:,a_or3)   = 0.0
           loss(:,:,a_or4)   = 0.0; loss(:,:,a_or5)   = 0.0; loss(:,:,a_or6)   = 0.0
           loss(:,:,a_or7)   = 0.0; loss(:,:,a_or8)   = 0.0; loss(:,:,a_or9)   = 0.0
           trans(:,:,a_or1)  = 0.0; trans(:,:,a_or2)  = 0.0; trans(:,:,a_or3)  = 0.0
           trans(:,:,a_or4)  = 0.0; trans(:,:,a_or5)  = 0.0; trans(:,:,a_or6)  = 0.0
           trans(:,:,a_or7)  = 0.0; trans(:,:,a_or8)  = 0.0; trans(:,:,a_or9)  = 0.0
         endif
         if (ICONA == 0) then
           excess(:,:,a_ami) = 0.0
           loss(:,:,a_ami  ) = 0.0
           trans(:,:,a_ami ) = 0.0
         endif
         if (ICONX == 0) then
           excess(:,:,a_nh4) = 0.0
           loss(:,:,a_nh4  ) = 0.0
           trans(:,:,a_nh4 ) = 0.0
           excess(:,:,a_chl) = 0.0
           loss(:,:,a_chl  ) = 0.0
           trans(:,:,a_chl ) = 0.0
         endif
         if ((ICONA == 0).and.(ICONX == 0)) then
           excess(:,:,a_nit) = 0.0
           loss(:,:,a_nit  ) = 0.0
           trans(:,:,a_nit ) = 0.0
         endif

 
    end subroutine condensloss


  subroutine condensation(VPT,DTIME,CCONDIN,EXCESS,LOSS,TRANS,N,IMAX,MASS, &
                          LVAP,FVAP,FLUXCM,FLUXC,FLUXCMC)
    
    !********************************************************************
    !
    !     C  O  N  D  E  N  S  A  T  I  O  N
    !
    !********************************************************************
    !
    ! Condensation module
    !
    !  1)   Initialize terms
    !  2)   Compute gas-phase/particle condensation loss 
    !  3)   Compute condensation flux for modes NU to AS
    !  4)   Compute condensation flux for mode  CS
    !
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !  INPUT
    !  -----
    !
    !     VPT: particle volume in bin         [m^3]
    !   DTIME: solver time step               [s]
    !   CCOND: condensation coefficient       [m^3/s]
    !  EXCESS: vapour conc. difference        [1/m^3]
    !    LOSS: particle condensation losss    [molec/(m^3s)]
    !       N: particle number conc. in bin   [1/m^3]
    !    MASS: component mass conc. in bin    [ng/m^3]
    ! 
    !
    !  OUTPUT
    !  ------
    !
    !  FLUXC:  condensation flux N in bin      [1/m^3] 
    ! FLUXCM:  condensation flux m in bin      [ng/m^3]
    ! FLUXCMC: flux of core particles in bin   [ng/m^3]
    !   LVAP:  condensation gas-phase loss     [1/s]
    !   FVAP:  condensation gas-phase prod.    [molec/(m^3s)]
    !
    !      method
    !      ------
    !      Analytical Predictor of Condensation
    !
    !      reference
    !      ---------
    !      Jacobson, M. Z., Fundamentals of Atmospheric Modeling, 
    !      Second Edition, Cambridge University Press, 2005
    !
    !      modifications
    !      -------------
    !      none
    !
    !********************************************************************   
    
     implicit none

    ! input
    integer, intent(in)                              :: IMAX
    real( dp), intent(in)                            :: DTIME
    real( dp), dimension(MMAX,IMAX),intent(in)       :: VPT,N
    real( dp), dimension(MMAX,IMAX,QMAX),intent(in)  :: EXCESS,LOSS
    real( dp), dimension(MMAX,IMAX,QMAX),intent(in)  :: CCONDIN
    real( dp), dimension(MMAX,IMAX,QMAX),intent(in)  :: TRANS
    real( dp), dimension(MMAX,IMAX,AMAX-1),intent(in):: MASS
         
    ! output
    real( dp), dimension(QMAX),intent(out)           :: FVAP,LVAP
 
    real( dp), dimension(MMAX,IMAX),intent(out)      :: FLUXC
    real( dp), dimension(MMAX,IMAX,QMAX),intent(out) :: FLUXCM
    real( dp), dimension(MMAX,IMAX),intent(out)      :: FLUXCMC

    real( dp), dimension(QMAX)                       :: JCOND,VVAPC    

    real( dp)   :: VPNEW
    integer     :: I,M,Q

! Initialize vapor molecule volume [m^3]
       VVAPC(A_SUL)=MVAP/DENV 
       VVAPC(A_MSA)=MVAP/DENV
       VVAPC(A_IO3)=MVAP/DENIO
       VVAPC(A_NIT)=MVAP/DENNI   ! VVAPN
       VVAPC(A_NH4)=MVAP/DENNI   ! VVAPN
       VVAPC(A_CHL)=MVAP/DENNI   ! VVAPN
       VVAPC(A_AMI)=MVAP/DENAM   ! VVAPA 
       VVAPC(A_OR1)=MVAP/DENXX   ! VVAPX
       VVAPC(A_OR2)=MVAP/DENXX   ! VVAPX
       VVAPC(A_OR3)=MVAP/DENXX   ! VVAPX
       VVAPC(A_OR4)=MVAP/DENXX   ! VVAPX       
       VVAPC(A_OR5)=MVAP/DENXX   ! VVAPX
       VVAPC(A_OR6)=MVAP/DENXX   ! VVAPX
       VVAPC(A_OR7)=MVAP/DENXX   ! VVAPX
       VVAPC(A_OR8)=MVAP/DENXX   ! VVAPX
       VVAPC(A_OR9)=MVAP/DENXX   ! VVAPX

! Initialize terms !!!!
! Initialize all prod+loss terms for gas phase 
       do Q=1,QMAX
         JCOND(Q)=0.
         FVAP(Q)=0.   
         LVAP(Q)=0.            
       end do

       ! initialization
       FLUXC(:,:)=0._dp
       FLUXCM(:,:,:)=0._dp
       FLUXCMC(:,:)=0._dp

! Compute condensation flux
         do M=NU,CS
          do I=1,IMAX
            VPNEW=0.0_dp
            do Q=1,QMAX
              JCOND(Q)=CCONDIN(M,I,Q)*EXCESS(M,I,Q)
            end do
            VPNEW=VPT(M,I)
            do Q=1,QMAX
              IF ((MASS(M,I,Q).GT.massmin)) THEN
                VPNEW=VPNEW + JCOND(Q)*VVAPC(Q)*DTIME
                ! if ((M.eq.2).and.(I.eq.IMAX)) print *,'VPNEW',M,I,Q,VPT(M,I),VPNEW
              ENDIF
            end do
            IF (VPNEW .LE. 0.0) THEN
              IF (ICONO .EQ. 1) THEN  !not including ELVOCs
                IF ((JCOND(A_OR1).LT.0.0).AND.(MASS(M,I,A_OR1).GT.massmin)) THEN
                   VPNEW=VPNEW-JCOND(A_OR1)*VVAPC(A_OR1)*DTIME
                ENDIF   
                IF ((JCOND(A_OR2).LT.0.0).AND.(MASS(M,I,A_OR2).GT.massmin)) THEN
                   VPNEW=VPNEW-JCOND(A_OR2)*VVAPC(A_OR2)*DTIME
                ENDIF   
                IF ((JCOND(A_OR4).LT.0.0).AND.(MASS(M,I,A_OR4).GT.massmin)) THEN
                   VPNEW=VPNEW-JCOND(A_OR4)*VVAPC(A_OR4)*DTIME
                ENDIF   
                IF ((JCOND(A_OR5).LT.0.0).AND.(MASS(M,I,A_OR5).GT.massmin)) THEN
                   VPNEW=VPNEW-JCOND(A_OR5)*VVAPC(A_OR5)*DTIME
                ENDIF
                IF ((JCOND(A_OR7).LT.0.0).AND.(MASS(M,I,A_OR7).GT.massmin)) THEN
                   VPNEW=VPNEW-JCOND(A_OR7)*VVAPC(A_OR7)*DTIME
                ENDIF   
                IF ((JCOND(A_OR8).LT.0.0).AND.(MASS(M,I,A_OR8).GT.massmin)) THEN
                   VPNEW=VPNEW-JCOND(A_OR8)*VVAPC(A_OR8)*DTIME
                ENDIF 
              ENDIF
              IF (ICONA .EQ. 1) THEN 
                IF (JCOND(A_AMI).LT.0.0) VPNEW=VPNEW-JCOND(A_AMI)*VVAPC(A_AMI)*DTIME
                IF (JCOND(A_NIT).LT.0.0) VPNEW=VPNEW-JCOND(A_NIT)*VVAPC(A_NIT)*DTIME
              ENDIF
            ENDIF
            !
            ! CHECK: VPT(M,I+1) must be greater than VPT(M,I)            
            ! VPNEW is bound between VPT(M,I-1) and VPT(M,I+1)
            !
            IF (I.EQ.1) THEN
              IF (VPT(M,I).GE.VPT(M,I+1)) THEN
                write(6,*) 'EMERGENCY STOP: i=1 vpt(i) >= vpt(i+1)',M,I
                stop 
              ENDIF
! MSK 18.11.2024 start new limits for i==1
              IF (M.NE.NU) THEN
                vpnew=max(vpnew,vpt(m,1)*0.50)
                vpnew=min(vpnew,vpt(m,i+1)*0.90)
              ELSE
                VPNEW=max(VPNEW,VPT(M-1,IMAX)*1.10)
                VPNEW=min(VPNEW,VPT(M,I+1)*0.90)
              ENDIF
! MSK 18.11.2024 end
            ENDIF
            IF ((I.EQ.IMAX).AND.(M.NE.CS)) THEN
              IF (VPT(M,I).GE.VPT(M+1,1)) THEN
                write(6,*) 'EMERGENCY STOP: i=imax vpt(i) >= vpt(i+1)',M,I
                stop 
              ENDIF             
              VPNEW=max(VPNEW,VPT(M,I-1)*1.10)
              VPNEW=min(VPNEW,VPT(M+1,1)*0.90)
            ENDIF
            IF ((I.GT.1).AND.(I.LT.IMAX)) THEN
              IF (VPT(M,I).GE.VPT(M,I+1)) THEN
                write(6,*) 'EMERGENCY STOP: 1<i<imax vpt(i) >= vpt(i+1)',M,I
                stop 
              ENDIF            
              VPNEW=max(VPNEW,VPT(M,I-1)*1.10)
              VPNEW=min(VPNEW,VPT(M,I+1)*0.90)
            ENDIF

            ! Gas phase terms
            do Q=1,QMAX
              FVAP(Q) = FVAP(Q) + LOSS(M,I,Q)
!MESA              LVAP(Q) = LVAP(Q) + N(M,I)*CCONDIN(M,I,Q)
              LVAP(Q) = LVAP(Q) + TRANS(M,I,Q)
            end do


            ! Compute condensation flux for modes NU to CS
            IF (I.EQ.1) THEN
               IF (M.EQ.NU) THEN    !special: [NU,1]
                 IF (VPNEW.GE.VPT(M,I)) THEN
                 ! number fluxes
                   FLUXC(M,I)=FLUXC(M,I)-N(M,I)
                   FLUXC(M,I)=FLUXC(M,I)+N(M,I)*( (VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I)) )
                   FLUXC(M,I+1)=FLUXC(M,I+1)+N(M,I)*( (VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I)) )
                 ! mass fluxes                
                 ! MSK 13.05.2023 Fixed brackets around volume ratio term
                   do Q=1,QMAX                              
                     FLUXCM(M,I,Q)=FLUXCM(M,I,Q)-MASS(M,I,Q)
                     FLUXCM(M,I,Q)=FLUXCM(M,I,Q)+MASS(M,I,Q)*( (VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I)) )
                     FLUXCM(M,I+1,Q)=FLUXCM(M,I+1,Q)+MASS(M,I,Q)*( (VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I)) )
                   end do
                   FLUXCMC(M,I)=FLUXCMC(M,I)-MASS(M,I,A_XXX)
                   FLUXCMC(M,I)=FLUXCMC(M,I)+MASS(M,I,A_XXX)*( (VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I)) )
                   FLUXCMC(M,I+1)=FLUXCMC(M,I+1)+MASS(M,I,A_XXX)*( (VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I)) )
                 ELSE  ! VPNEW LT VPT(M,I)                 
                 ! number fluxes
                    FLUXC(M,I)=FLUXC(M,I)-N(M,I)
                    FLUXC(M,I)=FLUXC(M,I)+N(M,I)*(VPNEW/VPT(M,I))
                 ! mass fluxes                
                    do Q=1,QMAX 
                      FLUXCM(M,I,Q)=FLUXCM(M,I,Q)-MASS(M,I,Q)
                      FLUXCM(M,I,Q)=FLUXCM(M,I,Q)+MASS(M,I,Q)*(VPNEW/VPT(M,I))
                    end do
                    FLUXCMC(M,I)=FLUXCMC(M,I)-MASS(M,I,A_XXX)
                    FLUXCMC(M,I)=FLUXCMC(M,I)+MASS(M,I,A_XXX)*(VPNEW/VPT(M,I))    
                 ENDIF           
               ELSE   ! AI,AS,CS modes
                 IF (VPNEW.GE.VPT(M,I)) THEN
                 ! number fluxes
                    FLUXC(M,I)=FLUXC(M,I)-N(M,I)
                    FLUXC(M,I)=FLUXC(M,I)+N(M,I)*( (VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I)) )
                    FLUXC(M,I+1)=FLUXC(M,I+1)+N(M,I)*( (VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I)) )
                 ! mass fluxes                
                    do Q=1,QMAX                              
                       FLUXCM(M,I,Q)=FLUXCM(M,I,Q)-MASS(M,I,Q)
                       FLUXCM(M,I,Q)=FLUXCM(M,I,Q)+MASS(M,I,Q)*( (VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I)) )
                       FLUXCM(M,I+1,Q)=FLUXCM(M,I+1,Q)+MASS(M,I,Q)*( (VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I)) )
                    end do
                    FLUXCMC(M,I)=FLUXCMC(M,I)-MASS(M,I,A_XXX)
                    FLUXCMC(M,I)=FLUXCMC(M,I)+MASS(M,I,A_XXX)*( (VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I)) )
                    FLUXCMC(M,I+1)=FLUXCMC(M,I+1)+MASS(M,I,A_XXX)*( (VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I)) )
                 ELSE  ! VPNEW LT VPT(M,I)
                 ! number fluxes
                    FLUXC(M,I)=FLUXC(M,I)-N(M,I)
                    FLUXC(M,I)=FLUXC(M,I)+N(M,I)*( (VPNEW-VPT(M-1,IMAX))/(VPT(M,I)-VPT(M-1,IMAX)) )
                    FLUXC(M-1,IMAX)=FLUXC(M-1,IMAX)+N(M,I)*( (VPNEW-VPT(M,I))/(VPT(M-1,IMAX)-VPT(M,I)) )
                 ! mass fluxes                
                    do Q=1,QMAX 
                      FLUXCM(M,I,Q)=FLUXCM(M,I,Q)-MASS(M,I,Q)
                      !BUG FIX 16.05.2013
                      fluxcm(m,i,q)=fluxcm(m,i,q)+mass(m,i,q)*( (vpnew-vpt(m-1,imax))/(vpt(m,i)-vpt(m-1,imax)) )
                      ! END BUG FIX                     
                      FLUXCM(M-1,IMAX,Q)=FLUXCM(M-1,IMAX,Q)+MASS(M,I,Q)*( (VPNEW-VPT(M,I)) &
                                        /(VPT(M-1,IMAX)-VPT(M,I)) )
                    end do
                    FLUXCMC(M,I)=FLUXCMC(M,I)-MASS(M,I,A_XXX)
                    FLUXCMC(m,i)=FLUXCMC(m,i)+mass(m,i,a_xxx)*( (vpnew-vpt(m-1,imax))/(vpt(m,i)-vpt(m-1,imax)) )
                    FLUXCMC(M-1,IMAX)=FLUXCMC(M-1,IMAX)+MASS(M,I,A_XXX)*( (VPNEW-VPT(M,I)) &
                                        /(VPT(M-1,IMAX)-VPT(M,I)) )         
                 ENDIF            
               ENDIF
            ENDIF

            IF ((I.GT.1).AND.(I.LT.IMAX)) THEN
               IF (VPNEW.GE.VPT(M,I)) THEN
               ! number fluxes
                 FLUXC(M,I)=FLUXC(M,I)-N(M,I)
                 FLUXC(M,I)=FLUXC(M,I)+N(M,I)*( (VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I)) )
                 FLUXC(M,I+1)=FLUXC(M,I+1)+N(M,I)*( (VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I)) )
               ! mass fluxes                
                 do Q=1,QMAX                              
                   FLUXCM(M,I,Q)=FLUXCM(M,I,Q)-MASS(M,I,Q)
                   FLUXCM(M,I,Q)=FLUXCM(M,I,Q)+MASS(M,I,Q)*( (VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I)) )
                   FLUXCM(M,I+1,Q)=FLUXCM(M,I+1,Q)+MASS(M,I,Q)*( (VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I)) )
                 end do
                 FLUXCMC(M,I)=FLUXCMC(M,I)-MASS(M,I,A_XXX)
                 FLUXCMC(M,I)=FLUXCMC(M,I)+MASS(M,I,A_XXX)*( (VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I)) )
                 FLUXCMC(M,I+1)=FLUXCMC(M,I+1)+MASS(M,I,A_XXX)*( (VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I)) )           
               ELSE  ! VPNEW LT VPT(M,I)
               ! number fluxes
                 FLUXC(M,I)=FLUXC(M,I)-N(M,I)
                 FLUXC(M,I)=FLUXC(M,I)+N(M,I)*( (VPNEW-VPT(M,I-1))/(VPT(M,I)-VPT(M,I-1)) )
                 FLUXC(M,I-1)=FLUXC(M,I-1)+N(M,I)*( (VPNEW-VPT(M,I))/(VPT(M,I-1)-VPT(M,I)) )
               ! mass fluxes                
                 do Q=1,QMAX 
                   FLUXCM(M,I,Q)=FLUXCM(M,I,Q)-MASS(M,I,Q)
                   FLUXCM(M,I,Q)=FLUXCM(M,I,Q)+MASS(M,I,Q)*( (VPNEW-VPT(M,I-1))/(VPT(M,I)-VPT(M,I-1)) )
                   FLUXCM(M,I-1,Q)=FLUXCM(M,I-1,Q)+MASS(M,I,Q)*( (VPNEW-VPT(M,I))/(VPT(M,I-1)-VPT(M,I)) )
                 end do
                 FLUXCMC(M,I)=FLUXCMC(M,I)-MASS(M,I,A_XXX)
                 FLUXCMC(M,I)=FLUXCMC(M,I)+MASS(M,I,A_XXX)*( (VPNEW-VPT(M,I-1))/(VPT(M,I)-VPT(M,I-1)) )
                 FLUXCMC(M,I-1)=FLUXCMC(M,I-1)+MASS(M,I,A_XXX)*( (VPNEW-VPT(M,I))/(VPT(M,I-1)-VPT(M,I)) )             
               ENDIF
            ENDIF
            
            IF (I.EQ.IMAX) THEN
               IF (M.EQ.CS) THEN    !special: [CS,IMAX]   
                  IF (VPNEW.GE.VPT(M,I)) THEN
                  ! number fluxes
                     FLUXC(M,I)=FLUXC(M,I)-N(M,I)
                     FLUXC(M,I)=FLUXC(M,I)+N(M,I)*VPNEW/VPT(M,I)
                  ! mass fluxes
                     do Q=1,QMAX
                       FLUXCM(M,I,Q)=FLUXCM(M,I,Q)-MASS(M,I,Q)
                       FLUXCM(M,I,Q)=FLUXCM(M,I,Q)+MASS(M,I,Q)*VPNEW/VPT(M,I)
                     end do
                     FLUXCMC(M,I)=FLUXCMC(M,I)-MASS(M,I,A_XXX)
                     FLUXCMC(M,I)=FLUXCMC(M,I)+MASS(M,I,A_XXX)*VPNEW/VPT(M,I)
                  ELSE  ! VPNEW LT VPT(M,I)
                  ! number fluxes
                     FLUXC(M,I)=FLUXC(M,I)-N(M,I)
                     FLUXC(M,I)=FLUXC(M,I)+N(M,I)*( (VPNEW-VPT(M,I-1))/(VPT(M,I)-VPT(M,I-1)) )
                     FLUXC(M,I-1)=FLUXC(M,I-1)+N(M,I)*( (VPNEW-VPT(M,I))/(VPT(M,I-1)-VPT(M,I)) )
                  ! mass fluxes                
                     do Q=1,QMAX 
                       FLUXCM(M,I,Q)=FLUXCM(M,I,Q)-MASS(M,I,Q)
                       FLUXCM(M,I,Q)=FLUXCM(M,I,Q)+MASS(M,I,Q)*( (VPNEW-VPT(M,I-1))/(VPT(M,I)-VPT(M,I-1)) )
                       FLUXCM(M,I-1,Q)=FLUXCM(M,I-1,Q)+MASS(M,I,Q)*( (VPNEW-VPT(M,I))/(VPT(M,I-1)-VPT(M,I)) )
                     end do
                     FLUXCMC(M,I)=FLUXCMC(M,I)-MASS(M,I,A_XXX)
                     FLUXCMC(M,I)=FLUXCMC(M,I)+MASS(M,I,A_XXX)*( (VPNEW-VPT(M,I-1))/(VPT(M,I)-VPT(M,I-1)) )
                     FLUXCMC(M,I-1)=FLUXCMC(M,I-1)+MASS(M,I,A_XXX)*( (VPNEW-VPT(M,I))/(VPT(M,I-1)-VPT(M,I)) )
                  ENDIF         
               ELSE   ! NU,AI,AS modes
                  IF (VPNEW.GE.VPT(M,I)) THEN
                  ! number fluxes
                     FLUXC(M,I)=FLUXC(M,I)-N(M,I)
                     FLUXC(M,I)=FLUXC(M,I)+N(M,I)*( (VPT(M+1,1)-VPNEW)/(VPT(M+1,1)-VPT(M,I)) )
                     FLUXC(M+1,1)=FLUXC(M+1,1)+N(M,I)*( (VPNEW-VPT(M,I))/(VPT(M+1,1)-VPT(M,I)) )
                  ! mass fluxes
                     do Q=1,QMAX
                      FLUXCM(M,I,Q)=FLUXCM(M,I,Q)-MASS(M,I,Q)
                      FLUXCM(M,I,Q)=FLUXCM(M,I,Q)+MASS(M,I,Q)*( (VPT(M+1,1)-VPNEW)/(VPT(M+1,1)-VPT(M,I)) )
                      FLUXCM(M+1,1,Q)=FLUXCM(M+1,1,Q)+MASS(M,I,Q)*( (VPNEW-VPT(M,I))/(VPT(M+1,1)-VPT(M,I)) )
                     end do
                     FLUXCMC(M,I)=FLUXCMC(M,I)-MASS(M,I,A_XXX)
                     FLUXCMC(M,I)=FLUXCMC(M,I)+MASS(M,I,A_XXX)*( (VPT(M+1,1)-VPNEW)/(VPT(M+1,1)-VPT(M,I)) )
                     FLUXCMC(M+1,1)=FLUXCMC(M+1,1)+MASS(M,I,A_XXX)*( (VPNEW-VPT(M,I))/(VPT(M+1,1)-VPT(M,I)) )  
                  ELSE  ! VPNEW LT VPT(M,I)
                  ! number fluxes
                      FLUXC(M,I)=FLUXC(M,I)-N(M,I)
                      FLUXC(M,I)=FLUXC(M,I)+N(M,I)*( (VPNEW-VPT(M,I-1))/(VPT(M,I)-VPT(M,I-1)) )
                      FLUXC(M,I-1)=FLUXC(M,I-1)+N(M,I)*( (VPNEW-VPT(M,I))/(VPT(M,I-1)-VPT(M,I)) )
                  ! mass fluxes                
                      do Q=1,QMAX 
                       FLUXCM(M,I,Q)=FLUXCM(M,I,Q)-MASS(M,I,Q)
                       FLUXCM(M,I,Q)=FLUXCM(M,I,Q)+MASS(M,I,Q)*( (VPNEW-VPT(M,I-1))/(VPT(M,I)-VPT(M,I-1)) )
                       FLUXCM(M,I-1,Q)=FLUXCM(M,I-1,Q)+MASS(M,I,Q)*( (VPNEW-VPT(M,I))/(VPT(M,I-1)-VPT(M,I)) )
                      end do
                      FLUXCMC(M,I)=FLUXCMC(M,I)-MASS(M,I,A_XXX)
                      FLUXCMC(M,I)=FLUXCMC(M,I)+MASS(M,I,A_XXX)*( (VPNEW-VPT(M,I-1))/(VPT(M,I)-VPT(M,I-1)) )
                      FLUXCMC(M,I-1)=FLUXCMC(M,I-1)+MASS(M,I,A_XXX)*( (VPNEW-VPT(M,I))/(VPT(M,I-1)-VPT(M,I)) )
                  ENDIF
               ENDIF
            ENDIF
             
          end do     ! end I-loop
         end do      ! end M-loop

  end subroutine condensation
 
  subroutine apc_update_gasc(M_oc,NVAPO,DTIME,CATO,FVAP,LVAP,DNVAP,NVAP)

    !********************************************************************
    !
    !     A  P  C  :  U  P  D  A  T  E    G  A  S    P  H  A  S  E 
    !
    !********************************************************************
    !
    ! JACOBSON APC SHEME
    !
    ! APC module to update particle mass concentrations after new gas phase
    !     concentrations are found.
    !     M.Z. Jacobson, p.544 ("Fundamentals of Atmospheric Modelling, 2005")
    !
    !    Compute new gas phase concentration after condensation
    !    constrained by the mass-balance equation
    !    M.J. Jacobson, p.544, eq. 16.71
    !    Upper limit: Cg_new=min(Cg_new,Ctot)
    !              Ca_sum in [ug m^-3]
    !              NVAP      [molec m^-3]
    !              molec2ug(MW) converts [molec cm^-3] --> [ug m^-3]
    !              Ctot      [molec m^-3]
    !    1) calculate Ctot=Ca_sum+Cg_old
    !    2) use upper limit to calculate Cg_new
    !    3) calculate delta_Cg=Cg_old-Cg_new
    !    4) for each bin calculate KOND_new
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !  INPUT
    !  -----
    !
    !  NVAPO:  old gas phase conc.             [molec/m^3]
    !   LVAP:  condensation gas-phase loss     [1/s]
    !   FVAP:  condensation gas-phase prod.    [molec/(m^3s)]
    !   CATO:  total component mass conc.      [ug/m^3]
    !   DTIME: solver time step                [s]    
    !
    !
    !  OUTPUT
    !  ------
    !
    !    NVAP: new gas phase conc.             [molec/m^3]
    !   DNVAP: gas phase conc. diff.           [molec/m^3]
    !
    !      method
    !      ------
    !      Analytical Predictor of Condensation
    !
    !      reference
    !      ---------
    !      Jacobson, M. Z., Fundamentals of Atmospheric Modeling, 
    !      Second Edition, Cambridge University Press, 2005
    !
    !      modifications
    !      -------------
    !      none
    !
    !********************************************************************   
     
     implicit none

    ! input
    REAL( dp), intent(in)                            :: DTIME
    real( dp), dimension(NSOA), intent(in)           :: M_oc    
    REAL( dp), dimension(QMAX),intent(in)            :: NVAPO,CATO,FVAP,LVAP    
             
    ! output
    REAL( dp), dimension(QMAX),intent(out)           :: NVAP
    REAL( dp), dimension(QMAX),intent(out)           :: DNVAP
    
    REAL( dp), dimension(QMAX)                       :: CAT 

! Initialize terms !!!!
        DNVAP(:)=0.
! Bug corrected 08.07.2014 (Eq. 16.71)
! NVAP(A_SUL) =(NVAP(A_SUL)+(FVAP(A_SUL)*DTIME))/(1._dp+LVAP(A_SUL)*DTIME)
! has to be:
! NVAP(A_SUL) =(NVAPO(A_SUL)+(FVAP(A_SUL)*DTIME))/(1._dp+LVAP(A_SUL)*DTIME)

        if ((ICONS == 1).or.(ICONS == 2)) then
          ! sulphuric acid
          CAT(A_SUL)  =(CATO(A_SUL)*1.e6*(1./molec2ug(M_H2SO4)))+NVAPO(A_SUL)
          NVAP(A_SUL) =(NVAPO(A_SUL)+(FVAP(A_SUL)*DTIME))/(1._dp+LVAP(A_SUL)*DTIME)
          NVAP(A_SUL) =min(NVAP(A_SUL),CAT(A_SUL))
          DNVAP(A_SUL)=NVAPO(A_SUL)-NVAP(A_SUL)
          ! MSA
          CAT(A_MSA)  =(CATO(A_MSA)*1.e6*(1./molec2ug(M_msa)))+NVAPO(A_MSA)
          NVAP(A_MSA) =(NVAPO(A_MSA)+(FVAP(A_MSA)*DTIME))/(1._dp+LVAP(A_MSA)*DTIME)
          NVAP(A_MSA) =min(NVAP(A_MSA),CAT(A_MSA))
          DNVAP(A_MSA)=NVAPO(A_MSA)-NVAP(A_MSA)
          ! iodic acid
          CAT(A_IO3)  =(CATO(A_IO3)*1.e6*(1./molec2ug(M_hio3)))+NVAPO(A_IO3)
          NVAP(A_IO3) =(NVAPO(A_IO3)+(FVAP(A_IO3)*DTIME))/(1._dp+LVAP(A_IO3)*DTIME)
          NVAP(A_IO3) =min(NVAP(A_IO3),CAT(A_IO3))
          DNVAP(A_IO3)=NVAPO(A_IO3)-NVAP(A_IO3)
        endif
        IF (ICONO .EQ. 1) THEN
          ! SOA-1
          CAT(A_OR1) =(CATO(A_OR1)*1.e6*(1./molec2ug(M_oc(1))))+NVAPO(A_OR1)
          NVAP(A_OR1) =(NVAPO(A_OR1)+(FVAP(A_OR1)*DTIME))/(1._dp+(LVAP(A_OR1)*DTIME))
          NVAP(A_OR1) =min(NVAP(A_OR1),CAT(A_OR1))
          DNVAP(A_OR1)=NVAPO(A_OR1)-NVAP(A_OR1)
          ! SOA-2
          CAT(A_OR2)  =(CATO(A_OR2)*1.e6*(1./molec2ug(M_oc(2))))+NVAPO(A_OR2)
          NVAP(A_OR2) =(NVAPO(A_OR2)+(FVAP(A_OR2)*DTIME))/(1._dp+(LVAP(A_OR2)*DTIME))
          NVAP(A_OR2) =min(NVAP(A_OR2),CAT(A_OR2))
          DNVAP(A_OR2)=NVAPO(A_OR2)-NVAP(A_OR2)
          ! SOA-3
          CAT(A_OR3)  =(CATO(A_OR3)*1.e6*(1./molec2ug(M_oc(3))))+NVAPO(A_OR3)
          NVAP(A_OR3) =(NVAPO(A_OR3)+(FVAP(A_OR3)*DTIME))/(1._dp+(LVAP(A_OR3)*DTIME))
          NVAP(A_OR3) =min(NVAP(A_OR3),CAT(A_OR3))        
          DNVAP(A_OR3)=NVAPO(A_OR3)-NVAP(A_OR3)
          ! SOA-4
          CAT(A_OR4)  =(CATO(A_OR4)*1.e6*(1./molec2ug(M_oc(4))))+NVAPO(A_OR4)
          NVAP(A_OR4) =(NVAPO(A_OR4)+(FVAP(A_OR4)*DTIME))/(1._dp+(LVAP(A_OR4)*DTIME))
          NVAP(A_OR4) =min(NVAP(A_OR4),CAT(A_OR4))
          DNVAP(A_OR4)=NVAPO(A_OR4)-NVAP(A_OR4)
          ! SOA-5
          CAT(A_OR5) =(CATO(A_OR5)*1.e6*(1./molec2ug(M_oc(5))))+NVAPO(A_OR5)
          NVAP(A_OR5) =(NVAPO(A_OR5)+(FVAP(A_OR5)*DTIME))/(1._dp+(LVAP(A_OR5)*DTIME))
          NVAP(A_OR5) =min(NVAP(A_OR5),CAT(A_OR5))
          DNVAP(A_OR5)=NVAPO(A_OR5)-NVAP(A_OR5)
          ! SOA-6
          CAT(A_OR6) =(CATO(A_OR6)*1.e6*(1./molec2ug(M_oc(6))))+NVAPO(A_OR6)
          NVAP(A_OR6) =(NVAPO(A_OR6)+(FVAP(A_OR6)*DTIME))/(1._dp+(LVAP(A_OR6)*DTIME))
          NVAP(A_OR6) =min(NVAP(A_OR6),CAT(A_OR6))
          DNVAP(A_OR6)=NVAPO(A_OR6)-NVAP(A_OR6)              
          ! SOA-7
          CAT(A_OR7)  =(CATO(A_OR7)*1.e6*(1./molec2ug(M_oc(7))))+NVAPO(A_OR7)
          NVAP(A_OR7) =(NVAPO(A_OR7)+(FVAP(A_OR7)*DTIME))/(1._dp+(LVAP(A_OR7)*DTIME))
          NVAP(A_OR7) =min(NVAP(A_OR7),CAT(A_OR7))
          DNVAP(A_OR7)=NVAPO(A_OR7)-NVAP(A_OR7)
          ! SOA-8
          CAT(A_OR8)  =(CATO(A_OR8)*1.e6*(1./molec2ug(M_oc(8))))+NVAPO(A_OR8)
          NVAP(A_OR8) =(NVAPO(A_OR8)+(FVAP(A_OR8)*DTIME))/(1._dp+(LVAP(A_OR8)*DTIME))
          NVAP(A_OR8) =min(NVAP(A_OR8),CAT(A_OR8))
          DNVAP(A_OR8)=NVAPO(A_OR8)-NVAP(A_OR8)
          ! SOA-9
          CAT(A_OR9)  =(CATO(A_OR9)*1.e6*(1./molec2ug(M_oc(9))))+NVAPO(A_OR9)
          NVAP(A_OR9) =(NVAPO(A_OR9)+(FVAP(A_OR9)*DTIME))/(1._dp+(LVAP(A_OR9)*DTIME))
          NVAP(A_OR9) =min(NVAP(A_OR9),CAT(A_OR9))
          DNVAP(A_OR9)=NVAPO(A_OR9)-NVAP(A_OR9)

        ENDIF
        IF (ICONA .EQ. 1) THEN
          ! amine
          CAT(A_AMI)  =(CATO(A_AMI)*1.e6*(1./molec2ug(M_nit)))+NVAPO(A_AMI)
          NVAP(A_AMI) =(NVAPO(A_AMI)+(FVAP(A_AMI)*DTIME))/(1._dp+LVAP(A_AMI)*DTIME)
          NVAP(A_AMI) =min(NVAP(A_AMI),CAT(A_AMI))
          DNVAP(A_AMI)=NVAPO(A_AMI)-NVAP(A_AMI)
          ! nitrate         
          CAT(A_NIT)  =(CATO(A_NIT)*1.e6*(1./molec2ug(M_nit)))+NVAPO(A_NIT)
          NVAP(A_NIT) =(NVAPO(A_NIT)+(FVAP(A_NIT)*DTIME))/(1._dp+LVAP(A_NIT)*DTIME)     
          NVAP(A_NIT) =min(NVAP(A_NIT),CAT(A_NIT))
          DNVAP(A_NIT)=NVAPO(A_NIT)-NVAP(A_NIT)                
        ENDIF
        IF (ICONX .EQ. 1) THEN
          ! ammonium (=NH4)
          ! 2020-06-25 corrected MW of NH4 (replaced M_nit by M_nh3)
          CAT(A_NH4)  =(CATO(A_NH4)*1.e6*(1./molec2ug(M_nh3)))+NVAPO(A_NH4)
          NVAP(A_NH4) =(NVAPO(A_NH4)+(FVAP(A_NH4)*DTIME))/(1._dp+LVAP(A_NH4)*DTIME)
          NVAP(A_NH4) =min(NVAP(A_NH4),CAT(A_NH4))
          DNVAP(A_NH4)=NVAPO(A_NH4)-NVAP(A_NH4)      
          ! nitrate         
          CAT(A_NIT)  =(CATO(A_NIT)*1.e6*(1./molec2ug(M_nit)))+NVAPO(A_NIT)
          NVAP(A_NIT) =(NVAPO(A_NIT)+(FVAP(A_NIT)*DTIME))/(1._dp+LVAP(A_NIT)*DTIME)
          NVAP(A_NIT) =min(NVAP(A_NIT),CAT(A_NIT))
          DNVAP(A_NIT)=NVAPO(A_NIT)-NVAP(A_NIT)
          ! chloride         
          CAT(A_CHL)  =(CATO(A_CHL)*1.e6*(1./molec2ug(M_hcl)))+NVAPO(A_CHL)
          NVAP(A_CHL) =(NVAPO(A_CHL)+(FVAP(A_CHL)*DTIME))/(1._dp+LVAP(A_CHL)*DTIME)  
          NVAP(A_CHL) =min(NVAP(A_CHL),CAT(A_CHL))
          DNVAP(A_CHL)=NVAPO(A_CHL)-NVAP(A_CHL)
        ENDIF

  end subroutine apc_update_gasc
    
  subroutine apc_update_pmass(DENOC,ROOPW,DTIME,IMAX,M_oc,KOND,FLUXCM,MASS,MMXO)

    !********************************************************************
    !
    !     A  P  C  :  U  P  D  A  T  E    M  A  S  S 
    !
    !********************************************************************
    !
    ! JACOBSON APC SHEME
    !
    ! APC module to update particle mass concentrations after new gas phase
    !     concentrations are found.
    !     M.J. Jacobson, p.544 ("Fundamentals of Atmospheric Modelling, 2005")
    !
    !    Compute number and mass concentration after condensation
    !    1) Calculate Caer_new
    !    2) Caer_new=max(Caer_new,0)
    !    New particle phase concentration. M.J. Jacobson, p.545, eq. 16.72
    !    convert NVAP [molec/m^-3] --> [ng/m^-3] 
    !    lower limit: Caer_new=max(Caer_new,0)
    !    upper limit: SUM Caer_new=min(SUM Caer_new,Ctot-Cg_new)
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !  INPUT
    !  -----
    !
    !   ROOPW: particle density in bin         [kg/m^3]
    !   DTIME: solver time step                [s]
    !    M_OC: molecular weight                [g/mol]
    !    KOND: condensation rate               [molec/m^3s]
    !  FLUXCM: condensation flux m in bin      [ng/m^3] 
    !    MASS: component mass conc. in bin     [ng/m^3]
    !    MOLD: component old mass conc. in bin [ng/m^3]
    !   DNVAP: gas phase conc. diff.           [molec/m^3]
    !
    !
    !  OUTPUT
    !  ------
    !
    !    MMXO: component new conc. in bin      [ng/m^3] 
    !
    !      method
    !      ------
    !      Analytical Predictor of Condensation
    !
    !      reference
    !      ---------
    !      Jacobson, M. Z., Fundamentals of Atmospheric Modeling, 
    !      Second Edition, Cambridge University Press, 2005
    !
    !      modifications
    !      -------------
    !      none
    !
    !********************************************************************   
    
     implicit none

    ! input
    integer, intent(in)                              :: IMAX
    real( dp), intent(in)                            :: DTIME
    real( dp), intent(in)                            :: DENOC
    real( dp), dimension(NSOA), intent(in)           :: M_oc
    real( dp), dimension(MMAX,IMAX),intent(in)       :: ROOPW
    real( dp), dimension(MMAX,IMAX,AMAX),intent(in)  :: MASS
    real( dp), dimension(MMAX,IMAX,QMAX),intent(in)  :: FLUXCM,KOND
             
    ! output
    real( dp), dimension(MMAX,IMAX,AMAX),intent(out) :: MMXO

    ! local
    real( dp), dimension(MMAX,IMAX,AMAX)             :: MMX
    real( dp), dimension(NSOA)                       :: MOC

    integer     :: M,I,Q


! Initialize terms !!!!
        do M=NU,CS
         do I=1,IMAX
          do Q=1,QMAX
           !output mass = input mass 
           MMXO(M,I,Q)=MASS(M,I,Q)           
          end do
         end do
        end do    

! Molecular weight MAC of organic components [kg/molec]

        MOC(:) = M_oc(:) *1.661e-27


! Compute mass concentration change by condensation
!
! 1) Calculate Caer_new
! 2) Caer_new=max(Caer_new,0)
! 3) Second check of Caer_new (only for ICONO and ICONA)
!
! MSK 14.11.2014: minimum mmxo set to massmin (not zero)
!
        do M=NU,CS
         do I=1,IMAX
           if ((ICONS == 1).or.(ICONS == 2)) then
             ! sulphuric acid and MSA            
             MMX(M,I,A_SUL)  = MASS(M,I,A_SUL) +   &
                             KOND(M,I,A_SUL)*MB*DTIME*CONVM*ROOPW(M,I)/DENV
             MMX(M,I,A_SUL)  = MMX(M,I,A_SUL) + FLUXCM(M,I,A_SUL)
             MMXO(M,I,A_SUL) = max(MMX(M,I,A_SUL), massmin)

             MMX(M,I,A_MSA)  = MASS(M,I,A_MSA) +   &
                             KOND(M,I,A_MSA)*MB*DTIME*CONVM*ROOPW(M,I)/DENMS
             MMX(M,I,A_MSA)  = MMX(M,I,A_MSA) + FLUXCM(M,I,A_MSA)           
             MMXO(M,I,A_MSA) = max(MMX(M,I,A_MSA), massmin)

             MMX(M,I,A_IO3)  = MASS(M,I,A_IO3) +   &
                             KOND(M,I,A_IO3)*MIO*DTIME*CONVM*ROOPW(M,I)/DENIO
             MMX(M,I,A_IO3)  = MMX(M,I,A_IO3) + FLUXCM(M,I,A_IO3)          
             MMXO(M,I,A_IO3) = max(MMX(M,I,A_IO3), massmin)
           endif
           IF (ICONO .EQ. 1) THEN
             ! biogenic secondary oxygenated
             MMX(M,I,A_OR1)  = MASS(M,I,A_OR1) +   &
                             KOND(M,I,A_OR1)*MOC(1)*DTIME*CONVM*ROOPW(M,I)/DENOC
             MMX(M,I,A_OR1)  = MMX(M,I,A_OR1) + FLUXCM(M,I,A_OR1)
             MMXO(M,I,A_OR1)  = max(MMX(M,I,A_OR1), massmin)

             MMX(M,I,A_OR2)  = MASS(M,I,A_OR2) +    &
                             KOND(M,I,A_OR2)*MOC(2)*DTIME*CONVM*ROOPW(M,I)/DENOC
             MMX(M,I,A_OR2)  = MMX(M,I,A_OR2) + FLUXCM(M,I,A_OR2)
             MMXO(M,I,A_OR2)  = max(MMX(M,I,A_OR2), massmin)

             MMX(M,I,A_OR3)  = MASS(M,I,A_OR3) +    &
                             KOND(M,I,A_OR3)*MOC(3)*DTIME*CONVM*ROOPW(M,I)/DENOC
             MMX(M,I,A_OR3)  = MMX(M,I,A_OR3) + FLUXCM(M,I,A_OR3)
!MSK 2024-11-05 bug fix: massmin limit causes bumps in size distribution                              
             !!!MMXO(M,I,A_OR3)  = max(MMX(M,I,A_OR3), massmin)
             MMXO(M,I,A_OR3)  = MMX(M,I,A_OR3)

             ! aromatic secondary oxygenated
             MMX(M,I,A_OR4)  = MASS(M,I,A_OR4) +   &
                             KOND(M,I,A_OR4)*MOC(4)*DTIME*CONVM*ROOPW(M,I)/DENOC
             MMX(M,I,A_OR4)  = MMX(M,I,A_OR4) + FLUXCM(M,I,A_OR4)
             MMXO(M,I,A_OR4)  = max(MMX(M,I,A_OR4), massmin)
             MMX(M,I,A_OR5)  = MASS(M,I,A_OR5) +    &
                             KOND(M,I,A_OR5)*MOC(5)*DTIME*CONVM*ROOPW(M,I)/DENOC
             MMX(M,I,A_OR5)  = MMX(M,I,A_OR5) + FLUXCM(M,I,A_OR5)
             MMXO(M,I,A_OR5)  = max(MMX(M,I,A_OR5), massmin)
             MMX(M,I,A_OR6)  = MASS(M,I,A_OR6) +    &
                             KOND(M,I,A_OR6)*MOC(6)*DTIME*CONVM*ROOPW(M,I)/DENOC
             MMX(M,I,A_OR6)  = MMX(M,I,A_OR6) + FLUXCM(M,I,A_OR6)       
!MSK 2024-11-05 bug fix: massmin limit causes bumps in size distribution
             !!!MMXO(M,I,A_OR6)  = max(MMX(M,I,A_OR6), massmin)
             MMXO(M,I,A_OR6)  = MMX(M,I,A_OR6)

             ! primary emitted (n-alkanes)
             MMX(M,I,A_OR7)  = MASS(M,I,A_OR7) +   &
                             KOND(M,I,A_OR7)*MOC(7)*DTIME*CONVM*ROOPW(M,I)/DENOC
             MMX(M,I,A_OR7)  = MMX(M,I,A_OR7) + FLUXCM(M,I,A_OR7)
             MMXO(M,I,A_OR7)  = max(MMX(M,I,A_OR7), massmin)
             MMX(M,I,A_OR8)  = MASS(M,I,A_OR8) +    &
                             KOND(M,I,A_OR8)*MOC(8)*DTIME*CONVM*ROOPW(M,I)/DENOC
             MMX(M,I,A_OR8)  = MMX(M,I,A_OR8) + FLUXCM(M,I,A_OR8)
             MMXO(M,I,A_OR8)  = max(MMX(M,I,A_OR8), massmin)
             MMX(M,I,A_OR9)  = MASS(M,I,A_OR9) +    &
                             KOND(M,I,A_OR9)*MOC(9)*DTIME*CONVM*ROOPW(M,I)/DENOC
             MMX(M,I,A_OR9)  = MMX(M,I,A_OR9) + FLUXCM(M,I,A_OR9)
!MSK 2024-11-05 bug fix: massmin limit causes bumps in size distribution
             !!! MMXO(M,I,A_OR9)  = max(MMX(M,I,A_OR9), massmin)
             MMXO(M,I,A_OR9)  = MMX(M,I,A_OR9)
           ENDIF
           IF (ICONA .EQ. 1) THEN        
             ! amine condensation
             MMX(M,I,A_AMI)  = MASS(M,I,A_AMI) +   &
                             KOND(M,I,A_AMI)*MAN*DTIME*CONVM*ROOPW(M,I)/DENAM 
             MMX(M,I,A_AMI)  = MMX(M,I,A_AMI) + FLUXCM(M,I,A_AMI)                              
             MMXO(M,I,A_AMI)  = max(MMX(M,I,A_AMI), massmin)

             ! amine-nitrate condensation    
             MMX(M,I,A_NIT)  = MASS(M,I,A_NIT) +   &
                             KOND(M,I,A_NIT)*MAN*DTIME*CONVM*ROOPW(M,I)/DENNI 
             MMX(M,I,A_NIT)  = MMX(M,I,A_NIT) + FLUXCM(M,I,A_NIT)                              
             MMXO(M,I,A_NIT)  = max(MMX(M,I,A_NIT), massmin)           
           ENDIF
           IF (ICONX .EQ. 1) THEN          
             ! ammonium nitrate condensation
             ! 2021-05-23 corrected MW of NH4 
             MMX(M,I,A_NH4)  = MASS(M,I,A_NH4) +   &
                             KOND(M,I,A_NH4)*MNH*DTIME*CONVM*ROOPW(M,I)/DENAM
             MMX(M,I,A_NH4)  = MMX(M,I,A_NH4) + FLUXCM(M,I,A_NH4)                             
             MMXO(M,I,A_NH4) = max(MMX(M,I,A_NH4), massmin)

             MMX(M,I,A_NIT)  = MASS(M,I,A_NIT) +   &
                             KOND(M,I,A_NIT)*MAN*DTIME*CONVM*ROOPW(M,I)/DENNI 
             MMX(M,I,A_NIT)  = MMX(M,I,A_NIT) + FLUXCM(M,I,A_NIT)                             
             MMXO(M,I,A_NIT) = max(MMX(M,I,A_NIT), massmin)

             MMX(M,I,A_CHL)  = MASS(M,I,A_CHL) +   &
                             KOND(M,I,A_CHL)*MAN*DTIME*CONVM*ROOPW(M,I)/DENNI 
             MMX(M,I,A_CHL)  = MMX(M,I,A_CHL) + FLUXCM(M,I,A_CHL)                             
             MMXO(M,I,A_CHL) = max(MMX(M,I,A_CHL), massmin)
           ENDIF
          end do
         end do


  end subroutine apc_update_pmass



    !----------------------------------------------------------------------
    !----------------  INCLOUD CONDENSATION       -------------------------
    !----------------------------------------------------------------------

  subroutine condensation_incloud( IMAX,DTIME,VPT,ROOP,DPA,    &
                 N,MASS,caqdiff, DPCRIT,                       &
                 FLUXC,FLUXCM )
    !********************************************************************
    !
    !     C  O  N  D  E  N  S  A  T  I  O  N   I N   C L O U D
    !
    !********************************************************************
    !
    ! Condensation in cloud module
    !
    !   Computes condensation flux for aqueous phase modes
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !  INPUT
    !  -----
    !
    !  DPCRIT: critical diameter                  [m]
    !     VPT: particle volume in bin             [m^3]
    !     DPA: particle diameter                  [m^3]
    !    ROOP: wet particle density               [kg/m^3]
    !       N: particle number conc. in bin       [1/m^3]
    !    MASS: component mass conc. in bin        [ng/m^3]
    ! CAQDIFF: concentration change dCaq          [ng/m^3]
    ! 
    !
    !  OUTPUT
    !  ------
    !
    !  FLUXC:  condensation flux N in bin      [1/m^3] 
    ! FLUXCM:  condensation flux m in bin      [ng/m^3]
    !
    !  reference
    !  ---------
    !  S.N Pandis, J.H. Seinfeld, C. Pilinis,
    !      Chemical composition differences in fog and
    !      cloud droplets of different sizes,
    !      Atmospheric Environment, 24A(7), 1957-1969, 1990a
    !
    !********************************************************************   
    
     implicit none

    ! input
    integer, intent(in)                               :: IMAX
    real( dp), intent(in)                             :: DTIME
    real( dp), intent(in)                             :: DPCRIT

    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in)  :: DPA
    real( dp), dimension(MMAX,IMAX),intent(in)        :: VPT
    real( dp), dimension(MMAX,IMAX),intent(in)        :: N
    real( dp), dimension(MMAX,IMAX),intent(in)        :: ROOP
    real( dp),dimension(aqmax,APN),intent(in)         :: caqdiff

    real( dp), dimension(MMAX,IMAX,AMAX-1),intent(in) :: MASS
         
    ! output
    real( dp), dimension(MMAX,IMAX),intent(out)       :: FLUXC
    real( dp), dimension(MMAX,IMAX,aqmax),intent(out) :: FLUXCM

    ! local
    real( dp)                                         :: VPNEW
    real( dp)                                         :: JCOND
    real( dp),dimension(aqmax,APN)                    :: cdiff

    integer     :: I,M
    integer     :: L


       ! initialization
         FLUXC(:,:)=0._dp
         FLUXCM(:,:,:)=0._dp


       ! Calculate change of volume in a size bin
       ! Based on Pandis et al. (1990a) EQ (7)
       ! JCOND = Dp^2*dDp/dt = dV/dt
       ! convert kg to ng with CONVM
       ! AI mode is aqueous phase mode 1

         do M=AI,MMAX             
           do I=1,IMAX

             VPNEW = VPT(M,I)
             JCOND = 0._dp

          ! JCOND has unit [m^6/ng]
          ! Conversion of ROOP from kg/m^3 into ng/m^3
             if (N(M,I).gt.nucomin) then
               JCOND = 2._dp/(pi*ROOP(M,I)*CONVM*N(M,I))
             endif


          ! No growth below critical diameter
             if (DPA(M,I).le.DPCRIT) then
               VPNEW=VPT(M,I)
             else

          ! New volume after condensation
                do L=1,aqmax

          ! Weighting of mass change rate dCaq for each bin
                  !!cdiff(L,M-1) = caqdiff(L,M-1) * (DPAW(M,I)/sumdpcs(M))
                  cdiff(L,M-1) = caqdiff(L,M-1) / IMAX

          ! Condensation of aqueous species L
          ! cdiff: change of mass rate per bin [ng/m^3/s]
          ! zkc=M-1
                  VPNEW = VPNEW + JCOND * cdiff(L,M-1) *DTIME

               enddo
             endif

             ! Lower bound of volume growth
             VPNEW=max( VPT(M,I),VPNEW )
             ! Upper bound of volume growth
             if (I.eq.IMAX) then
               if (M.eq.CS) then
                  VPNEW=min( VPT(M,I)*1.10,VPNEW )
               else
                  VPNEW=min( VPT(M+1,1),VPNEW )
               endif
             else
              VPNEW=min( VPT(M,I+1),VPNEW )
             endif

          !   print *,'cond',M,I,N(M,I),JCOND,cdiff(1,M-1),VPT(M,I),VPNEW


       ! Subroutine condensation_incloud
       ! Compute condensation flux for
       ! aqueous phase modes
       ! Aerosol species
       !      SVI = 1   => A_SUL (1)
       !      MSA = 2   => A_MSA (2)
       !      NO3 = 3   => A_NIT (3)
       !      DMA = 4   => A_AMI (4)
       !      NH4 = 5   => A_NH4 (5)
       !      IOD = 6   => A_IO3 (6)
       !      OXA = 7   => A_OR1 (7)
       !      SUC = 8   => A_OR2 (8)
       !      SAL = 9   => A_SAL (17)

             if (I.eq.1) then
          ! number fluxes
               FLUXC(M,I)=FLUXC(M,I)-N(M,I)
               FLUXC(M,I)=FLUXC(M,I)+N(M,I)*(VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I))
               FLUXC(M,I+1)=FLUXC(M,I+1)+N(M,I)*(VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I))
          ! mass fluxes of the aqueous phase species (1 to aqmax-1)  
               do L=1,aqmax-1               
                 FLUXCM(M,I,L)=FLUXCM(M,I,L)-MASS(M,I,L)
                 FLUXCM(M,I,L)=FLUXCM(M,I,L)+MASS(M,I,L)*(VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I))
                 FLUXCM(M,I+1,L)=FLUXCM(M,I+1,L)+MASS(M,I,L)*(VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I))
               enddo
          ! mass flux of seasalt
               FLUXCM(M,I,aqmax)=FLUXCM(M,I,aqmax)-MASS(M,I,A_SAL)
               FLUXCM(M,I,aqmax)=FLUXCM(M,I,aqmax)+MASS(M,I,A_SAL)*(VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I))
               FLUXCM(M,I+1,aqmax)=FLUXCM(M,I+1,aqmax)+MASS(M,I,A_SAL)*(VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I))
            endif

             if ((I.gt.1).and.(I.lt.IMAX)) then
          ! number fluxes
               FLUXC(M,I)=FLUXC(M,I)-N(M,I)
               FLUXC(M,I)=FLUXC(M,I)+N(M,I)*(VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I))
               FLUXC(M,I+1)=FLUXC(M,I+1)+N(M,I)*(VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I))
          ! mass fluxes of the aqueous phase species (1 to aqmax-1)  
               do L=1,aqmax-1             
                  FLUXCM(M,I,L)=FLUXCM(M,I,L)-MASS(M,I,L)
                  FLUXCM(M,I,L)=FLUXCM(M,I,L)+MASS(M,I,L)*(VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I))
                  FLUXCM(M,I+1,L)=FLUXCM(M,I+1,L)+MASS(M,I,L)*(VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I))
               enddo
          ! mass flux of seasalt
               FLUXCM(M,I,aqmax)=FLUXCM(M,I,aqmax)-MASS(M,I,A_SAL)
               FLUXCM(M,I,aqmax)=FLUXCM(M,I,aqmax)+MASS(M,I,A_SAL)*(VPT(M,I+1)-VPNEW)/(VPT(M,I+1)-VPT(M,I))
               FLUXCM(M,I+1,aqmax)=FLUXCM(M,I+1,aqmax)+MASS(M,I,A_SAL)*(VPNEW-VPT(M,I))/(VPT(M,I+1)-VPT(M,I))
             endif


             if (I.eq.IMAX) then
          ! number fluxes
               FLUXC(M,I)=FLUXC(M,I)-N(M,I)
               if (M.lt.CS) then
                 FLUXC(M,I)=FLUXC(M,I)+N(M,I)*(VPT(M+1,1)-VPNEW)/(VPT(M+1,1)-VPT(M,I))
                 FLUXC(M+1,1)=FLUXC(M+1,1)+N(M,I)*(VPNEW-VPT(M,I))/(VPT(M+1,1)-VPT(M,I))
               endif
          ! mass fluxes of the aqueous phase species (1 to aqmax-1)  
               do L=1,aqmax-1
                  FLUXCM(M,I,L)=FLUXCM(M,I,L)-MASS(M,I,L)
                 if (M.lt.CS) then 
                   FLUXCM(M,I,L)=FLUXCM(M,I,L)+MASS(M,I,L)*(VPT(M+1,1)-VPNEW)/(VPT(M+1,1)-VPT(M,I))
                   FLUXCM(M+1,1,L)=FLUXCM(M+1,1,L)+MASS(M,I,L)*(VPNEW-VPT(M,I))/(VPT(M+1,1)-VPT(M,I))
                 endif
               enddo
          ! mass flux of seasalt
               FLUXCM(M,I,aqmax)=FLUXCM(M,I,aqmax)-MASS(M,I,A_SAL)
               if (M.lt.CS) then 
                 FLUXCM(M,I,aqmax)=FLUXCM(M,I,aqmax)+       &
                                   MASS(M,I,A_SAL)*(VPT(M+1,1)-VPNEW)/(VPT(M+1,1)-VPT(M,I))
                 FLUXCM(M+1,1,aqmax)=FLUXCM(M+1,1,aqmax)+   &
                                   MASS(M,I,A_SAL)*(VPNEW-VPT(M,I))/(VPT(M+1,1)-VPT(M,I))
               endif
             endif

        enddo
      enddo


  end subroutine condensation_incloud



    !----------------------------------------------------------------------
    !----------------  CONDENSATION COEFFICIENTS  -------------------------
    !----------------------------------------------------------------------


  subroutine condensation_coeff2(temp,rh,press,DPA,UPT_SU,UPT_MS,DC0,   &
                                 DNB,DNBMSA,DNBIOD,IMAX)
    !----------------------------------------------------------------------
    !
    !****  Calculates condensation coefficients
    !
    !      author
    !      -------
    !      Liisa Pirjola and Matthias Karl
    !
    !      Dr. Liisa Pirjola
    !      Docent
    !      Department of Physics
    !      University of Helsinki
    !      P.O Box 64, FI-00014 Helsinki, Finland
    !      Department of Technology
    !      Metropolia University of Applied Sciences
    !      P.O. Box 4071, FI-01600 Vantaa, Finland
    !
    !      purpose
    !      -------
    !      calculation of condensation coefficient
    !      for variable accomodation coefficient alpha
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature                [K]
    !           press    air pressure                   [Pa]
    !           rh       relative humidity              [-]
    !           DPA      particle diameter              [m]
    !
    !        output:
    !           DNB      condensation coefficient H2SO4 [m^3 s^-1]
    !           DNBMSA   condensation coefficient MSA   [m^3 s^-1]
    !           DNBIOD   condensation coefficient HIO3  [m^3 s^-1]
    !           DC0      diffusion coefficient          [cm^2 s^-1]
    !           UPT_SU   uptake rate H2SO4              [m4/(s*molec)]
    !           UPT_MS   uptake rate MSA                [m4/(s*molec)]
    !
    !
    !      method
    !      ------
    !      condensational flux of sulphuric acid
    !      see eg. Pirjola and Kulmala, Atm. Res. 46, pp 321-347
    !
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      Pirjola, L., Kulmala, M. Modelling the formation of H2SO4-H2O
    !      particles in rural, urban and marine conditions,
    !      Atmospheric Research, 46, 321-347, 1998.
    !      Pirjola, L., Korhonen, H., Kulmala, M. Condensation/evaporation
    !      of insoluble organic vapor as functions of source rate and 
    !      saturation vapor pressure, J. Geophys. Res., 107, D11, 4108, 
    !      doi:10.1029/2011JD001228, 2002.
    !
    !------------------------------------------------------------------

    implicit none

    ! io
    INTEGER, intent(in)                              :: IMAX
    real( dp), intent(in)                            :: temp,press,rh
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA
    real( dp), dimension(MMAX,IMAX),intent(out)      :: DNB,DNBMSA,DNBIOD
    !REAL( dp), intent(out)                           :: GR
    real( dp), intent(out)                           :: UPT_SU,DC0
    real( dp), intent(out)                           :: UPT_MS

    INTEGER, PARAMETER                     :: hydimax=5
    real( dp)                              :: rhoht,pres,hydpw
    real( dp)                              :: N1,DNBI,wabs
    real( dp), DIMENSION(hydimax)          :: rhoh,kprod
    real( dp), DIMENSION(hydimax)          :: rhohm
    real( dp)                              :: DNBMSAI,DNBIODI
    real( dp)                              :: rhohmt
    real( dp), DIMENSION(hydimax+1)        :: dc,lamda,ca,knm,bm
    real( dp), DIMENSION(MMAX,0:(IMAX+1))  :: RP
    real( dp)                              :: CBAR

    INTEGER I,J,L,M

! new parameter list for H2SO4 and MSA (17.02.2013)
    real,  dimension(3) :: &
           alphs = (/ 0.5,   0.13,  0.5    /), &
           m_s   = (/ 98.08, 96.11, 175.91 /), &
           sigma = (/ 19.7,  40.3,  19.7   /), &
           as    = (/ 13.1,  15.7,  13.1   /) 
! overwrite with unity accommodation coefficients
     if (ICONS == 2) then
       alphs(1) = 1.0
       alphs(2) = 1.0
       alphs(3) = 1.0
     endif
! get hydrate coefficients
     wabs=rh*waterps(temp)
     CALL sulfhydrates(temp,wabs,hydpw,kprod,rhoh,rhohm)

! use pressure in [atm]
     pres=press/101325._dp

! diffusion coefficient DC(H2SO4) for output  (cm^2 s^-1)
     N1=0._dp
     DC0 = 1.e-3_dp*temp**1.75_dp*SQRT(1.0/M_air + 1.0/(M_S(1) + &
             N1*M_h2o))
     DC0 = DC0/(pres*(SIGMA(1)**(1.0_dp/3.0_dp) + (51.96_dp +     &
             N1*AS(1))**(1.0_dp/3.0_dp))**2.0_dp)

! total hydrate fraction (5 H2O molecules)
     rhoht=0.
     rhohmt=0.
     DNB=0.
     DNBMSA=0.
     DO L = 1,hydimax
         rhoht = rhoht + rhoh(L)
         rhohmt = rhohmt + rhohm(L)
     END DO

     RP(:,:)=0._dp
     do M=1,MMAX
      do I=1,IMAX
       RP(M,I)=DPA(M,I)*0.5_dp  ! [m]
       DNBI=0._dp
       DNBMSAI=0._dp
       DNBIODI=0._dp
! Loop for condensation of H2SO4, MSA and HIO3 as hydrates with up to five H2O molecules
! all molar weights given in g/mol
       DO J=1,3
       ! Corrected 30.01.2013 loop l=1,hydimax+1  
          do L=1, hydimax+1
! number of h2o molecules
       !     n1 = float(l)
            n1 = float(L-1)
! diffusion coefficient DC  (cm^2 s^-1)
            DC(L) = 1.e-3_dp*temp**1.75_dp*SQRT(1.0/M_air + 1.0/(M_S(J) + &
                  N1*M_h2o))
            DC(L) = DC(L)/(pres*(sigma(J)**(1.0_dp/3.0_dp) + (51.96_dp +     &
                  N1*AS(J))**(1.0_dp/3.0_dp))**2.0_dp)
! mean velocity in air (cm s^-1)
            CA(L) = 1e2_dp*(8.0_dp*8314.7_dp*temp / &
                  (pi*(m_s(J)+N1*M_H2O)))**(0.5_dp)
! mean free path (m) Seinfeld and Pandis (p. 457)
            LAMDA(L) = 32._dp*1.e-2_dp*DC(L)/(3._dp*pi*(1+((m_s(J)+       &
                     N1*M_H2O)/M_air))*CA(L))
            KNM(L) = LAMDA(L)/RP(M,I)
! transitional correction factor according to Fuchs and Sutugin
            BM(L) = (KNM(L)+1.0_dp)/(0.377_dp*KNM(L)+1.0_dp+4.0_dp*KNM(L)*KNM(L) / &
                 (3.0_dp*alphs(J)) + 4.0_dp*KNM(L)/(3.0_dp*alphs(J)))
          !write(6,*) 'J L DC CA LAMDA',J,L,DC(L),CA(L),LAMDA(L),BM(L),rhoh(L),rhoht
! DNBI and DNBMSAI in [cm^2/s]
            if (L==1) then  ! free acid
              if (J.eq.1) dnbi    = dnbi   +bm(L)*dc(L)*(1-rhoht)
              if (J.eq.2) dnbmsai = dnbmsai+bm(L)*dc(L)*(1-rhohmt)
              if (J.eq.3) dnbiodi = dnbiodi+bm(L)*dc(L)*(1-rhoht)
            else            ! hydrates
              if (J.eq.1) dnbi    = dnbi   +bm(L)*dc(L)*rhoh(L-1)
              if (J.eq.2) dnbmsai = dnbmsai+bm(L)*dc(L)*rhohm(L-1)
              if (J.eq.3) dnbiodi = dnbiodi+bm(L)*dc(L)*rhoh(L-1)
            end if

          END DO
       END DO
! DNB and DNBMSA in [m^3/s]
       DNB(M,I)    = 4.0e-4_dp*pi*RP(M,I)*DNBI
       DNBMSA(M,I) = 4.0e-4_dp*pi*RP(M,I)*DNBMSAI
       DNBIOD(M,I) = 4.0e-4_dp*pi*RP(M,I)*DNBIODI
      ! write(6,*) M,I,DNB(M,I),DNBIOD(M,I),'rad:',RP(M,I)
      end do
     end do


     ! Uptake rate of 1nm diameter particle (radius=0.5e-9m)
     ! According to Eq. 2-5 in Verheggen and Mozurkevich, ACP,6,2927-2942,2006
     !UPTR1=(1._dp/ALPHS(1))+(0.75_dp*RP(1,3)/LAMDA(1))-(0.47_dp*RP(1,3)/(RP(1,3)+LAMDA(1)))
     !UPTR2=(1._dp/ALPHS(2))+(0.75_dp*RP(1,3)/LAMDA(1))-(0.47_dp*RP(1,3)/(RP(1,3)+LAMDA(1)))
     !CBAR=SQRT(8._dp*k_B*temp/(pi*MVAP))    ! mean molecular speed
     !UPTAKE_SU=CBAR*MB/(4._dp*UPTR1*DENV)
     !UPTAKE_MS=CBAR*MB/(4._dp*UPTR2*DENMS)
     !   Growth Rate (not used)
     CBAR=SQRT(8._dp*k_B*temp/(pi*MVAP))
     UPT_SU=0.5_dp*CBAR*MVAP/DENV
     UPT_MS=0.5_dp*CBAR*MVAP/DENMS

  end subroutine condensation_coeff2

  subroutine orgcondens(press,temp,DPA,M_org1,DENOC,UPTAKE_ORG,DC0,ORGFLUX,IMAX)
    !----------------------------------------------------------------------
    !
    !****  Calculates organic condensation coefficients
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !
    !      calculates the flux of organic vapour using the diffusion
    !      coefficient and mean free path of sulphuric acid
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature                [K]
    !           press    air pressure                   [Pa]
    !           DPA      particle diameter              [m]
    !           M_org1   MW organic 1                   [g/mol]
    !           DENOC    particle density               [kg/m^3]
    !
    !        output:
    !           ORGFLUX  condensation coefficient       [m^3 s^-1]
    !           DC0      diffusion coefficient          [cm^2 s^-1]
    !           UPTAKE_ORG uptake rate organic          [m4/(s*molec)]
    !
    !      method
    !      ------
    !      Correction of Fuchs and Sutugin equation to consider molecule size
    !      and particle diffusion (DCP)
    !
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      Lehtinen, K.E.J. and M. Kulmala (2003). A model for particle formation
    !      and growth in the atmosphere with molecular resolution in size,
    !      Atmos. Chem. Phys., 3, 251-257.
    !
    !------------------------------------------------------------------

    implicit none

    ! input
    INTEGER, intent(in)                              :: IMAX
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA
    REAL( dp), intent(in)                            :: temp,press
    REAL( dp), intent(in)                            :: M_org1,DENOC
    ! output
    REAL( dp), dimension(MMAX,IMAX),intent(out)      :: ORGFLUX
    REAL( dp), intent(out)                           :: UPTAKE_ORG,DC0

    ! local
    REAL( dp), parameter                   :: ALPHA_ORG=1.0
    real( dp), parameter                   :: sigma_A = 5.85    ! Angstroem
    REAL( dp), DIMENSION(MMAX,0:(IMAX+1))  :: RP
    REAL( dp), dimension(MMAX,IMAX)        :: KNM0,BM0
    REAL( dp)                              :: LAMDA0,pres
    REAL( dp)                              :: UPTR,CBAR
    REAL( dp)                              :: RMOL0,VMOL0,DCP,M_part

    integer :: I,M

! use pressure in [atm]
     pres=press/101325._dp

! DC0 molecular diffusion coefficient in [cm2/s]
! Chapmans-Enskog equation with the first order approximation
! of the collision parameter, Omega = 1
! sigma is the collision diameter, estimated from
!           sigma_AB = (sigma_A + sigma_B)/2
!   A is index for trace gas, B is index for air
!         the molecular volume of the liquid, in Angstroem
! for Toluene: sigma_A = 5.85 Angstroem
!    Ref: J.R. Li, R.J. Kuppler and H.C. Zhou, 
!         Chem. Soc. Rev.,38, 1477-1504, 2009.
! Value of DC should be around 0.04 cm2/s to 0.05 cm2/s

     DC0    = molecdiff(M_org1,sigma_A,temp,pres)

     LAMDA0 = 3.0E-4_dp*DC0/(8.0_dp*8314.7_dp*temp/(pi*M_org1))**(0.5_dp)
     ! radius of molecule [m]
     VMOL0  = M_org1*1.E-3_dp/(N_A*DENOC)
     RMOL0  = 0.5_dp*(VMOL0*6.0_dp/pi)**(1.0_dp/3.0_dp)

     do M=1,MMAX
      do I=1,IMAX
        IF(M.EQ.1) THEN   !particle diffusion coefficient for NU
          M_part=RHOH2O*N_A*1.E3_dp*(1.0_dp/6.0_dp)*pi*DPA(M,I)**3.0_dp
          DCP=0.001_dp*temp**1.75*(1.0_dp/28.965_dp + 1.0_dp/M_part)**0.5_dp/pres
          DCP=DCP/(19.7_dp**(1.0_dp/3.0_dp) + 51.96_dp**(1.0_dp/3.0_dp))**2.0_dp
          DC0=DC0+DCP
          LAMDA0 = 3.0E-4_dp*DC0/( 8.0_dp*8314.7_dp*temp/(pi*M_org1)                &
                  + 8.0_dp*8314.7_dp*temp/(pi*M_part) )**(0.5_dp)
        ENDIF
        RP(M,I)=DPA(M,I)*0.5_dp  ! [m]
        KNM0(M,I)=LAMDA0/(RP(M,I)+RMOL0)
        BM0(M,I) = (KNM0(M,I) +1.0_dp)/(0.377_dp*KNM0(M,I)+1.0_dp+4.0_dp*KNM0(M,I)  &
                  * KNM0(M,I)/(3.0_dp*ALPHA_ORG)                                    &
                  + 4.0_dp*KNM0(M,I)/(3.0_dp*ALPHA_ORG))
        ORGFLUX(M,I) = BM0(M,I)*DC0
        ORGFLUX(M,I) = 4.0E-4_dp*pi*(RP(M,I)+RMOL0)*ORGFLUX(M,I)     
        !write(6,*) M,I,KNM0(M,I),BM0(M,I),ORGFLUX(M,I)
      end do
     end do

     ! Uptake rate of 1nm diameter particle (radius=0.5e-9m)
     ! According to Eq. 2-5 in Verheggen and Mozurkevich, ACP,6,2927-2942,2006
     ! LAMDA0 [m], CBAR [m/s], DENOC [kg/m3], MAH [kg/molec], UPTAKE_ORG [m4/(s*molec)]
     UPTR=(1._dp/ALPHA_ORG)+(0.75_dp*RP(1,1)/LAMDA0)-(0.47_dp*RP(1,1)/(RP(1,1)+LAMDA0))
     CBAR=SQRT(8._dp*k_B*temp/(pi*MVOC))    ! mean molecular speed
     UPTAKE_ORG=0.5*CBAR*MVOC/(4._dp*UPTR*DENOC)

  end subroutine orgcondens

  subroutine nitcondens(IMAX,press,temp,DPA,ALPHA_NIT,FC_NIT,NITFLUX)
    !----------------------------------------------------------------------
    !
    !****  Calculates ammoniumnitrate condensation coefficients
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !
    !   calculates the flux of ammonium/nitrate using the diffusion
    !   coefficient and mean free path of sulphuric acid
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature                [K]
    !           press    air pressure                   [Pa]
    !           DPA      particle diameter              [m]
    !           ALPHA_NIT accommodation coefficient
    !           FC_NIT   scaling factor
    !
    !        output:
    !           NITFLUX  condensation coefficient       [m^3 s^-1]
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

    ! input
    integer, intent(in)                              :: IMAX
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA
    real( dp), intent(in)                            :: temp,press,ALPHA_NIT,FC_NIT
    ! output
    real( dp), dimension(MMAX,IMAX),intent(out)      :: NITFLUX

    ! local
    REAL( dp), DIMENSION(MMAX,0:(IMAX+1))  :: RP
    REAL( dp), dimension(MMAX,IMAX)        :: KNM0,BM0
    real( dp), parameter                   :: Mwnit = 102._dp
    real( dp)                              :: sigma_nit
    real( dp)                              :: DC0,LAMDA0,pres
    !real( dp)                              :: UPTR,CBAR,UPTAKE_NIT

    integer :: I,M

! use pressure in [atm]
     pres=press/101325._dp

! DC0 molecular diffusion coefficient in [cm2/s]
! Chapmans-Enskog equation with the first order approximation
! of the collision parameter, Omega = 1
! sigma is the collision diameter, estimated from
!           sigma_AB = (sigma_A + sigma_B)/2
!   A is index for trace gas, B is index for air
!         the molecular volume of the liquid, in Angstroem
!   use sigma of sulfphuric acid:
!   sigma_A= 19.7**(1/3)   in Angstroem
! Value of DC should be around 0.04 cm2/s to 0.05 cm2/s

     sigma_nit = 19.7_dp**(1.0_dp/3.0_dp) 
     DC0    = molecdiff(Mwnit,sigma_nit,temp,pres)

     LAMDA0 = 3.0E-4_dp*DC0/(8.0_dp*8314.7_dp*temp/(pi*102._dp))**(0.5_dp)

     do M=1,MMAX
      do I=1,IMAX
        RP(M,I)=DPA(M,I)*0.5_dp  ! [m]
        KNM0(M,I)=LAMDA0/RP(M,I)
        BM0(M,I) = (KNM0(M,I) +1.0_dp)/(0.377_dp*KNM0(M,I)+1.0_dp+4.0_dp*KNM0(M,I)* &
           KNM0(M,I)/(3.0_dp*ALPHA_NIT)  + &
           4.0_dp*KNM0(M,I)/(3.0_dp*ALPHA_NIT))
        NITFLUX(M,I) = BM0(M,I)*DC0
        NITFLUX(M,I) = 4.0E-4_dp*pi*RP(M,I)*NITFLUX(M,I)*FC_NIT    
        !write(6,*) M,I,KNM0(M,I),BM0(M,I),NITFLUX(M,I)
      end do
     end do

     ! Uptake rate of 1nm diameter particle (radius=0.5e-9m)
     ! According to Eq. 2-5 in Verheggen and Mozurkevich, ACP,6,2927-2942,2006
     !UPTR=(1._dp/ALPHA_NIT)+(0.75_dp*RP(1,3)/LAMDA0)-(0.47_dp*RP(1,3)/(RP(1,3)+LAMDA0))
     !QCBAR=SQRT(8._dp*k_B*temp/(pi*MVAP))    ! mean molecular speed
     !UPTAKE_NIT=CBAR*MAN/(4._dp*UPTR*DENNI)

  end subroutine nitcondens

  subroutine chlcondens(IMAX,press,temp,DPA,CHLFLUX)
    !----------------------------------------------------------------------
    !
    !****  Calculates chlorine condensation coefficients
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !
    !   calculates the flux of hydrochloric acid using the diffusion
    !   coefficient and mean free path of sulphuric acid
    !   Expression is the same as for routine nitcondens
    !   except that the molecular weight of chlorine and the
    !   accommodation coefficient of 0.15 is used
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature                [K]
    !           press    air pressure                   [Pa]
    !           DPA      particle diameter              [m]
    !
    !        output:
    !           CHLFLUX  condensation coefficient       [m^3 s^-1]
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

    ! input
    integer, intent(in)                              :: IMAX
    real( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA
    real( dp), intent(in)                            :: temp,press
    ! output
    real( dp), dimension(MMAX,IMAX),intent(out)      :: CHLFLUX

    ! local
    REAL( dp), DIMENSION(MMAX,0:(IMAX+1))  :: RP
    REAL( dp), dimension(MMAX,IMAX)        :: KNM0,BM0
    real( dp)                              :: sigma_nit
    real( dp)                              :: DC0,LAMDA0,pres
! molecular weight
    real( dp), parameter                   :: Mw_hcl=36._dp
! accommodation coefficient
    real( dp), parameter                   :: alpha_hcl=0.15_dp


    integer :: I,M

! use pressure in [atm]
     pres=press/101325._dp


! DC0 molecular diffusion coefficient in [cm2/s]
! Chapmans-Enskog equation with the first order approximation
! of the collision parameter, Omega = 1
! sigma is the collision diameter, estimated from
!           sigma_AB = (sigma_A + sigma_B)/2
!   A is index for trace gas, B is index for air
!         the molecular volume of the liquid, in Angstroem
!   use sigma of sulfphuric acid:
!   sigma_A= 19.7**(1/3)   in Angstroem
! Value of DC should be around 0.04 cm2/s to 0.05 cm2/s

     sigma_nit = 19.7_dp**(1.0_dp/3.0_dp) 
     DC0    = molecdiff(Mw_hcl,sigma_nit,temp,pres)

     LAMDA0 = 3.0E-4_dp*DC0/(8.0_dp*8314.7_dp*temp/(pi*102._dp))**(0.5_dp)

     do M=1,MMAX
      do I=1,IMAX
        RP(M,I)=DPA(M,I)*0.5_dp  ! [m]
        KNM0(M,I)=LAMDA0/RP(M,I)
        BM0(M,I) = (KNM0(M,I) +1.0_dp)/(0.377_dp*KNM0(M,I)+1.0_dp+4.0_dp*KNM0(M,I)* &
           KNM0(M,I)/(3.0_dp*alpha_hcl)  + &
           4.0_dp*KNM0(M,I)/(3.0_dp*alpha_hcl))
        CHLFLUX(M,I) = BM0(M,I)*DC0
        CHLFLUX(M,I) = 4.0E-4_dp*pi*RP(M,I)*CHLFLUX(M,I) 
      end do
     end do


  end subroutine chlcondens

  
  subroutine eqnitrate(cnac,camin,kpmol,EQNIT,EQNH3)
    !----------------------------------------------------------------------
    !
    !****  Calculates molecular saturation concentration
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !
    !   calculates the molecular saturation concentration
    !   for condensation of alkylammonium nitrate salts
    !   formed from amine + nitric acid reaction
    !
    !      interface
    !      ---------
    !
    !      input is:
    !      temp:  temperature in K
    !      kpmol: dissociation constant in (molec/cm3)2
    !      cnac:  gas phase nitric acid concentration in molec/cm3
    !      camin: gas phase amin concentration in molec/cm3
    !
    !      output is:
    !      New equil. concentration (gas) on particle surface in molec./m3
    !      EQNIT, EQNH3
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

    ! input
    REAL( dp), intent(in)                  :: kpmol
    REAL( dp), intent(in)                  :: cnac,camin   ! molec/cm^3
    ! output
    REAL( dp), intent(out)                 :: EQNIT,EQNH3

    ! dissociation enthalpy, entropy, and heat capacity from kp_nitrate.xls
    ! calculated for methylammonium nitrate
    ! Murphy et al., ACP, 7, 2313-2337,2007
    REAL( dp),PARAMETER                    :: DH0=-179.69  ! kJ/mol
    REAL( dp),PARAMETER                    :: DS0=-313.31  ! J/molK
    REAL( dp),PARAMETER                    :: DCP=66.69    ! J/molK
    REAL( dp),PARAMETER                    :: R_mol=82.06  ! atm cm^3/molK

   !  REAL( dp)                              :: lnkp,kpppb,kpatm

    !! Calculation method with thermodyn. parameters
    !lnkp=((DS0-DCP)/R_gas) - ((DH0*1e3_dp-T0*DCP)/(R_gas*temp)) +  &
    !      ((DCP/R_gas)*LOG(temp/T0))
    !kpatm=1/(EXP(lnkp))
    !kpmol=(kpatm/((R_mol*temp)**2))*N_A**2

    !! use a fixed Kp for MEA:
    !! Exp. value for TEA-nitrate at 293 K: 1.85E-7 Pa^2
    !! Murphy et al., ACP, 7, 2313-2337,2007
    ! kpmol=1.13E22_dp

    !! use Kp(T) of NH4NO3 for MEA-nitrate
    !! Stelson et al., 1979
    !! kpmol=6.84E21_dp  ! NH4NO3 at 298 K
    !lnkp=84.6_dp-(24220._dp/temp)-6.1_dp*LOG(temp/298._dp)
    !kpppb=EXP(lnkp)
    !kpmol=kpppb*2.463E20*(298._dp/temp)

    ! Eq. conc. of ammonium in molec/cm^3
    EQNH3=((camin-cnac)/2) + SQRT( (((camin-cnac)*(camin-cnac))/4) + kpmol)
    ! Eq. conc. of nitric acid in molec/cm^3
    EQNIT=kpmol/EQNH3

    ! New equil. concentration (gas) on particle surface in molec./m3
    EQNIT=EQNIT*1.e6_dp
    EQNH3=EQNH3*1.e6_dp

  end subroutine eqnitrate


  subroutine eqnh4nitrate(cnac,camin,cchl,temp,kpmol_no3,kpmol_cl,  &
                          EQNIT,EQCHL,EQNH3)
    !----------------------------------------------------------------------
    !
    !****  Calculates molecular saturation concentration
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !
    !   calculates the molecular saturation concentration
    !   for condensation of ammonium nitrate salt
    !   formed from ammonia + nitric acid reaction
    !   extended for ammonium chloride salt as described in
    !   Jacobson (AST, 2005)
    !
    !      interface
    !      ---------
    !
    !      input is:
    !      temp:  temperature in K
    !      kpmol: dissociation constant in (molec/cm3)2
    !      cnac:  gas phase nitric acid concentration in molec/cm3
    !      camin: gas phase amin concentration in molec/cm3
    !      cchl:  gas phase hydrochloric acid in molec/cm3
    !
    !      output is:
    !      New equil. concentration (gas) on particle surface in molec./m3
    !      EQNIT, EQCHL, EQNH3
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
    !
    !    M.Z. Jacobson,  AEROSOL SCIENCE TECHNOLOGY, VOL. 39, 
    !          NO. 2, 92-103, DOI: 10.1080/027868290904546, 2005.
    !
    !------------------------------------------------------------------

    implicit none

    ! input
    real( dp), intent(in)                  :: cnac     ! molec/cm^3
    real( dp), intent(in)                  :: cchl     ! molec/cm^3
    real( dp), intent(in)                  :: camin    ! molec/cm^3
    real( dp), intent(in)                  :: temp     ! K

    ! output
    real( dp), intent(out)                 :: EQNIT
    real( dp), intent(out)                 :: EQCHL
    real( dp), intent(out)                 :: EQNH3
    real( dp), intent(out)                 :: kpmol_no3
    real( dp), intent(out)                 :: kpmol_cl

    ! local
    real( dp)                              :: tt
    real( dp)                              :: rt
    real( dp)                              :: lnkp
    real( dp)                              :: kpppb
    real( dp)                              :: kpatm
    real( dp)                              :: cnull
    real( dp)                              :: kpsum


    !! Gas-solid equilibrium Kp(NH4NO3)
    !! kpmol: dissociation constant in (molec/cm^3)^2
    !! Stelson et al., 1979
    !! kpmol=6.84E21_dp  ! NH4NO3 at 298 K
    !! Stelson Kp is about 1/4 smaller than in MESA
    lnkp=84.6_dp-(24220._dp/temp)-6.1_dp*LOG(temp/298._dp)
    kpppb=EXP(lnkp)
    kpmol_no3=kpppb*2.463E20*(298._dp/temp)

    !! Gas-solid equilibrium Kp(NH4Cl)
    !! kpmol=4.93E22_dp  ! NH4Cl at 298 K
    tt=298.15_dp/temp
    ! kpatm in [mol/kg/atm]
    kpatm=8.43E-17_dp*EXP(-71.0_dp*(tt-1._dp)+2.4_dp*(1._dp+LOG(tt)-tt))
    rt=82.056*temp                 ! [cm^3 atm/mol]
    kpmol_cl=kpatm/rt**2           ! [(mol/cm^3)^2]
    ! kpmol: dissociation constant in (molec/cm^3)^2
    !kpmol_cl=kpmol_cl *N_A**2
    kpmol_cl=kpmol_cl *N_A
    kpmol_cl=kpmol_cl *N_A
    ! reduced kp due to interaction with NO3
    kpmol_cl=kpmol_cl*0.01_dp

    !! Equil. conc. of ammonium in molec/cm^3
    cnull=camin-cnac-cchl
    kpsum=kpmol_no3+kpmol_cl
    if (cchl.eq.0._dp) kpsum=kpmol_no3
    if (cnac.eq.0._dp) kpsum=kpmol_cl

    EQNH3=( cnull/2 ) + 0.5_dp* SQRT( cnull*cnull + 4._dp*kpsum )

    !! Equil. conc. of nitric acid in molec/cm^3
    EQNIT=kpmol_no3/EQNH3
    EQCHL=kpmol_cl /EQNH3

    !! New equil. concentration (gas) 
    !! on particle surface in molec./m^3
    EQNIT=EQNIT*1.e6_dp
    EQCHL=EQCHL*1.e6_dp
    EQNH3=EQNH3*1.e6_dp

  end subroutine eqnh4nitrate




  subroutine soapartition(temp,           &
                        ca1_nul,ca2_nul,ca3_nul,ca4_nul,ca5_nul,ca6_nul,   &
                        cm1_in,cm2_in,cm3_in,fom_in,                       &
                        psatn,mwsoa,msize,ocfrc,pptot,                     &
                        csat1,csat2,csat3,csat4,csat5,csat6)

    !----------------------------------------------------------------------
    !
    !****  Calculates saturation concentration of SOA compounds
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !
    !      Solution from Seinfeld and Pandis (1998) book (13.5.2): 
    !      binary pseudo-ideal solution with pre-existing aerosol (OCp)
    !      2-D VBS set - prediction model for the activity coefficients
    !      from Donahue et al. (2011). 
    !
    !      interface
    !      ---------
    !
    !      input is:
    !!!      yamidoh: molar yield of SOA compound
    !!!      coh    : OH concentration, current time step         (molec/m3)
    !!!      cvoc   : parent VOC, current time step               (molec/m3)
    !!!      cg_null:	SOA gas phase concentration,  old time step (molec/m3)
    !      ca_null: SOA total aerosol concentration, old time step (ng/m3)
    !      cm_ini : Non-volatile OC concentration,   old time step (ng/m3)
    !      fom    : fraction of absorptive (organic) matter in total PM
    !      msize  : size of solvent (=nC + nO)
    !      ocfrc  : carbon fraction (O:C ratio)
    !      temp   : temperature         (K)
    !
    !      output is:
    !      cg_new : SOA gas phase concentration  (molec/m3)
    !      csat   : SOA saturation concentration (molec/m3)
    !
    !      method
    !      ------
    !
    !      call after chemistry solver
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !
    !   Donahue et al. (2011):
    !     Donahue, N. M., Epstein, S. A., Pandis, S. N., and
    !     Robinson, A. L., A two-dimensional volatility basis
    !     set: 1. organic-aerosol mixing thermodynamics,
    !     Atmos. Chem. Phys., 11, 3303-3318, 2011.
    !
    !   Seinfeld and Pandis (1998):                               
    !     Atmospheric Chemistry and Physics (0-471-17816-0)  
    !     chapter 13.5.2 Formation of binary ideal solution 
    !
    !   Odum et al. (1996):                                   
    !     Gas/particle partitioning and secondary organic    
    !     aerosol yields,  Environ. Sci. Technol. 30, 2580-2585. 
    !
    !   Bilde et al. (2003):
    !     Even-odd alternation of evaporation rates and vapor
    !     pressures of C3-C9 dicarboxylic acid aerosols, 
    !     Environ. Sci. Technol. 37, 1371-1378.
    !
    !
    !   References for the physical adsorption model: 
    !
    !   Pankow (1994):                                                  
    !     An absorption model of the gas/aerosol                      
    !     partitioning involved in the formation of                   
    !     secondary organic aerosol, Atmos. Environ. 28(2), 189-193.  
    !
    !   Ellison et al. (1999):
    !     Atmospheric processing of organic aerosols,
    !     J. Geophys. Res., 104(D9), 11,633-11,641.
    !
    !------------------------------------------------------------------

    implicit none

    ! input
    real( dp), intent(in)                  :: ca1_nul,ca2_nul,ca3_nul,ca4_nul
    real( dp), intent(in)                  :: ca5_nul,ca6_nul         ! ng m^-3
    real( dp), intent(in)                  :: cm1_in,cm2_in,cm3_in    ! ng m^-3
    real( dp),dimension(NSOA), intent(in)  :: psatn                   ! Pa
    real( dp),dimension(NSOA), intent(in)  :: mwsoa                   ! g/mol
    real( dp),dimension(NSOA), intent(in)  :: msize                   ! (nC + nO)
    real( dp),dimension(NSOA), intent(in)  :: ocfrc                   ! (O:C ratio)
    real( dp), intent(in)                  :: pptot                   ! ng m^-3
    real( dp), intent(in)                  :: fom_in                  ! ---
    real( dp), intent(in)                  :: temp                    ! K

    ! output
    real( dp), intent(out)                 :: csat1,csat2,csat3,csat4,csat5,csat6      

    integer                                :: i
    integer, parameter                     :: nsoal=6
    ! Pankow 1994: parameters for adsorption and absorption
    ! cbet = e^(Q1 - Qv)/RT determines shape of BET curve
    !    cbet = 100 polar surface, cbet = 0.1 unpolar surface
    !    Pankow: (Q1 - Qv) = 1.5 kcal/mol = 6276 J/mol
    !            (Q1 - Qv)/RT = 2.56, cbet = e^(Q1 - Qv)/RT = 12.9
    !    Pankow: C = 10^(-8.6) for urban particles
    !            C = nsites*atsp*temp*cbet/1600
    !            nsites = 7.e-11 cm^-2
    real( dp), parameter                   :: atsp   = 19.0     ! m^2 g^-1
    real( dp), parameter                   :: dpp    = 2.0e-7   ! m  (Dp=200 nm)  
    real( dp), parameter                   :: surfmol= 2.00e-19 ! m^2 = 20.0 A^2
    real( dp), parameter                   :: cbet   = 10.0
    real( dp), parameter                   :: gamma1 = 1.0      ! activity coeff.
    real( dp), parameter                   :: conini = 1.e-3    ! ug m^-3
    real( dp), parameter                   :: xinf   = 1.e-6
    !real( dp), parameter                   :: fom    = 0.3      ! fraction in OM
    real( dp), parameter                   :: bCO    = -0.3     ! carbon-oxygen non-ideality
    
    real( dp)                              :: npp,asurf,nsites,kadsr
    real( dp)                              :: cm_ini
    real( dp)                              :: xsoa1,xsoa2,xsoa3,xsoa4,xsoa5,xsoa6
    real( dp)                              :: xini1,xini2,xini3
    real( dp)                              :: xsolv,fsolv
    real( dp), dimension(nsoal)            :: Morg,Msiz,fcar
    real( dp), dimension(nsoal)            :: canull,psat,csat,xom
    real( dp), dimension(nsoal)            :: loggam,gammao

!--------------------------------------------------------------------------
    ! Initialization

    ! convert to ug m^-3
    canull(1)=ca1_nul*1.e-3
    canull(2)=ca2_nul*1.e-3
    canull(3)=ca3_nul*1.e-3    
    canull(4)=ca4_nul*1.e-3
    canull(5)=ca5_nul*1.e-3
    canull(6)=ca6_nul*1.e-3

    cm_ini=(cm1_in+cm2_in+cm3_in)*1.e-3 
    
    ! initialize
    ! molecular weight
    Morg(1)=mwsoa(1)     ! BSOV
    Morg(2)=mwsoa(2)     ! BLOV
    Morg(3)=mwsoa(4)     ! ASOV
    Morg(4)=mwsoa(5)     ! ALOV
    Morg(5)=mwsoa(7)     ! PIOV
    Morg(6)=mwsoa(8)     ! PSOV

    ! solute size
    Msiz(1)=msize(1)     ! BSOV
    Msiz(2)=msize(2)     ! BLOV
    Msiz(3)=msize(4)     ! ASOV
    Msiz(4)=msize(5)     ! ALOV
    Msiz(5)=msize(7)     ! PIOV
    Msiz(6)=msize(8)     ! PSOV

    ! carbon fraction 
    fcar(1)=ocfrc(1)     ! BSOV
    fcar(2)=ocfrc(2)     ! BLOV
    fcar(3)=ocfrc(4)     ! ASOV
    fcar(4)=ocfrc(5)     ! ALOV
    fcar(5)=ocfrc(7)     ! PIOV
    fcar(6)=ocfrc(8)     ! PSOV

    ! saturation pressure
    psat(1)=psatn(1)     ! BSOV
    psat(2)=psatn(2)     ! BLOV
    psat(3)=psatn(4)     ! ASOV
    psat(4)=psatn(5)     ! ALOV
    psat(5)=psatn(7)     ! PIOV
    psat(6)=psatn(8)     ! PSOV

    ! molar ratio (ug/m3*mol/g)
    xsoa1=canull(1)/Morg(1)
    xsoa2=canull(2)/Morg(2)
    xsoa3=canull(3)/Morg(3)
    xsoa4=canull(4)/Morg(4)
    xsoa5=canull(5)/Morg(5)
    xsoa6=canull(6)/Morg(6)

    xini1=(cm1_in*1.e-3)/mwsoa(3)    !BELV
    xini2=(cm2_in*1.e-3)/mwsoa(6)    !AELV
    xini3=(cm3_in*1.e-3)/mwsoa(9)    !PELV

    cm_ini=conini+cm_ini

!--------------------------------------------------------------------------


    ! Physical adsorption [Pankow 1994]
    !   pptot: mass conc primary particles (ng/m^3)
    !!!! aetot: mass conc total particles   (ng/m^3)
    ! Molecular composition of aerosol [Ellison 1999]
    !   Particle with Dp=200 nm 
    !   surface area: 1.3e-9 cm^2
    !   6.3e5 molecules at surface
    !   is struck by 15 OH radicals per second

   if (ISOA == 2) then
      npp    = pptot*1.e-12/(1./6.*pi*DENEC*dpp**3)  ! # m^-3_air
      !npp    = 1.                                ! # m^-3_air
      npp    = min(npp, 1.e5_dp)                     
      asurf  = pi*dpp**2                         ! m^2
      nsites = asurf/surfmol                     ! molecules / #
      ! only every 1/100000 site is free (=ca. 6 sites)
      nsites = nsites/1.e5                       ! sites / #
      nsites = nsites*npp                        ! sites m^-3_air
      nsites = nsites*dpp                        ! sites m^-2   
      nsites = nsites*1.e-4                      ! sites cm^-2
      nsites = min(nsites,1.e-4_dp)
      kadsr  = nsites*atsp*cbet*temp/1600.
      !write(6,*) 'N,SURF',npp,asurf
      !write(6,*) 'nsites',nsites
      !write(6,*) 'Kp_ads',kadsr/psat(2)
      !write(6,*) 'Kp_sor',(R_gas*temp*fom_in*760.)/(Morg(2)*gamma1*1.e6*psat(2))
    else
      kadsr  = 0.0
    endif
    
!--------------------------------------------------------------------------

    ! calculate activity coefficient
    ! based on Donahue et al. (2011), Equation (13)
    ! s = solvent (organic mixture)
    ! i = solute (invidiual SOA component)
    ! log10(gammao_i) = -2*bCO*Msiz_i*((fc_i)^2+(fc_s)^2-2*fc_i*fc_s)
    ! gammao should be 1 for bulk O:C of 0.25:1
    !     and >1 ... 2.5 for bulk O:C of 0.75:1
    ! ------------------------------------------
    ! total mass of solvent (organic mixture)
    ! add infinitesimal mass to avoid division by zero
    xsolv = xsoa1+xsoa2+xsoa3+xsoa4+xsoa5+xsoa6          +    & 
            xini1+xini2+xini3                            +    &
            xinf
    ! Also add the O:C ratio of the non-volatile OCp !
    ! carbon fraction of the solvent
    fsolv = (xsoa1/xsolv)*fcar(1)   + &
            (xsoa2/xsolv)*fcar(2)   + &
            (xsoa3/xsolv)*fcar(3)   + &
            (xsoa4/xsolv)*fcar(4)   + &
            (xsoa5/xsolv)*fcar(5)   + &
            (xsoa6/xsolv)*fcar(6)   + &
    ! non-volatile components
            (xini1/xsolv)*ocfrc(3)  + &
            (xini2/xsolv)*ocfrc(6)  + &
            (xini3/xsolv)*ocfrc(9)

    ! Limit fsolv to 0.3 to prevent very high gamma for PIOV and PSOV
    fsolv = min(fsolv,0.3_dp)

    do i=1,nsoal

       loggam(i) = (-2._dp)*bCO*Msiz(i)*( fcar(i)*fcar(i) +   &
                   fsolv*fsolv - (2._dp*fcar(i)*fsolv)  )
       gammao(i)  = 10._dp**(loggam(i))
!debug
!       print *,"log(gamma)",i,Msiz(i),fcar(i),loggam(i),gammao(i)
    end do
!debug
!    print *,'solvent (mol/ug/m3)',xsolv,fsolv


    ! calculate SOA saturation concentration (ug m^-3)
    ! input saturation vapor pressure (Pa) and MW (g/mol)
    do i=1,nsoal
      ! sat. concentration of (pure) non-interacting SOA
      ! csat_i (ug m^-3) is calculated from psat(T)_i (Pa)
      ! and the activity coefficient in the organic solution
      csat(i)=psat(i)*Morg(i)*gammao(i)*1.e6/(R_gas*temp*fom_in)

      if (ISOA == 2) then  ! physical adsorption
        csat(i)=1./( ( R_gas*temp*fom_in/(psat(i)*Morg(i)*gammao(i)*1.e6) )  &
                   + ( kadsr/psat(i)/760. ) )
      endif    

      ! sat. concentration of (n-)solution of SOA and OCp
      ! calculate mole fraction xom1 using aerosol mass
      !   aerosol mass concentration in ug m^-3
      !   _ini : pre-exisiting OCp
      ! xom1=(ca[1]/MW[1])/(SUM(ca[i]/M[i])+(cm_ini/MW_ini))
      ! csat1=csat1*xom1
      xom(i)=0._dp ! initialize
      IF (canull(i) .GT. 10.*conini) THEN    ! arbitrary threshold
        xom(i)=canull(i)/Morg(i)
        xom(i)=xom(i)/((xsoa1+xsoa2+xsoa3+xsoa4+xsoa5+xsoa6) +  &
               (cm_ini/MC))
        csat(i)=csat(i)*xom(i)
!debug
!          write(6,*) 'SOA part i',i,temp,xom(i),psat(i),csat(i)
      ELSE
        csat(i)=csat(i)*1.
!debug
!          write(6,*) 'SOA part i',i,temp,xom(i),psat(i),csat(i)
      END IF

    end do

!--------------------------------------------------------------------------

    ! return saturation concentration to _main
    csat1=csat(1)*1.e6/molec2ug(Morg(1))              ! molec m^-3   
    csat2=csat(2)*1.e6/molec2ug(Morg(2))              ! molec m^-3
    csat3=csat(3)*1.e6/molec2ug(Morg(3))              ! molec m^-3
    csat4=csat(4)*1.e6/molec2ug(Morg(4))              ! molec m^-3
    csat5=csat(5)*1.e6/molec2ug(Morg(5))              ! molec m^-3
    csat6=csat(6)*1.e6/molec2ug(Morg(6))              ! molec m^-3


  end subroutine soapartition
 
    
end module gde_condensation
