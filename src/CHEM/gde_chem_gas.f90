! <gde_chem_gas.f90 - A component of the Multicomponent
!                     Aerosol Dynamics Model MAFOR>
!*****************************************************************************! 
!* 
!*    Copyright (C) 2011-2021  Matthias Steffen Karl
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
!*****************************************************************************!
!*    All routines written by Matthias Karl
!* 
!*****************************************************************************! 
module gde_chem_gas

  use messy_cmn_photol_mem
  use messy_mecca_kpp_parameters
  use messy_mecca_kpp_global
  use messy_mecca_kpp

  use gde_constants, only       : k_B,pi,R_gas
  use gde_input_data, only      : NSOA
  use gde_toolbox,   only       : molecdiff
  use gde_photo,     only       : photo_mbl
  use gde_sensitiv,  only       : ICHAM

  private

  public :: emis_drydep, chamber_src, chamber_dil, chamber_loss
  public :: chemistry_solver
  public :: check_range


contains


  subroutine chemistry_solver(cgas,DTIME,firstloop,model_time,daytime,lat_deg, &
           temp,press,cair,jno2m,RH,F_HONO,fco,lwc,cvfac,xaer,k_exf,k_exb,    &
           k_exf_N2O5,k_exf_ClNO3,k_exf_BrNO3,   fcoj,jrnew)
    !----------------------------------------------------------------------
    !     
    ! Main routine to call the calculation of photolysis rates and
    !  do integration of chemistry
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      from input:
    !      input/output new gas phase concentration
    !
    !      interface
    !      ---------
    !
    !        input:
    !          k_exf        transfer coef. forward reaction 
    !          k_exb        transfer coef. backward reaction
    !          xaer         flag for aq. chem.
    !          cvfac        unit conversion factor
    !          lwc          liquid water content
    !          k_exf_N2O5   uptake coef. N2O5
    !          k_exf_ClNO3  uptake coef. ClNO3
    !          k_exf_BrNO3  uptake coef. BrNO3   
    !          DTIME        time step
    !          model_time   simulated time
    !          lat_deg      latitude degrees
    !          daytime      time of day
    !          temp         air temperature                    [K]
    !          press        air pressure                       [Pa]
    !          cair         c(air) (wet)                       [mcl/cm^3]
    !          jno2m        scaled JNO2                        [1/s]
    !          RH           rel. humidity                      [-]
    !          F_HONO       scaling of HONO source             [-]
    !          fco          light on/off                       [-]
    !
    !      method
    !      ------
    !      interface to the MECCA chemistry solver
    !      jx: photolysis rates in CAABA/MECCA
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

    implicit none

    ! input
    real( dp),dimension(nspec), intent(in out) :: cgas
    real( dp),dimension(APN,nspec),intent(in)  :: k_exf,k_exb  
    real( dp),dimension(APN),intent(in)        :: xaer,cvfac,lwc
    real( dp),dimension(APN),intent(in)        :: k_exf_N2O5
    real( dp),dimension(APN),intent(in)        :: k_exf_ClNO3,k_exf_BrNO3    
    real( dp), intent(in)                      :: DTIME
    real( dp), intent(in)                      :: model_time
    real( dp), intent(in)                      :: lat_deg,daytime
    real( dp), intent(in)                      :: temp,press,cair
    real( dp), intent(in)                      :: jno2m,RH,F_HONO,fco


    logical, intent(in)     :: firstloop

    ! output
    real( dp),dimension(IP_MAX),intent(out) :: jrnew
    real( dp), intent(out)  :: fcoj

    ! local
    integer             :: jo,jp
    real( dp)           :: RF
    real( dp)           :: dt

! KPP SOLVER
    integer, parameter  :: NBL = 1    ! N_block_length
    integer             :: status     ! error status
    real(dp), dimension(NBL,nspec) :: cbl    
   ! real( dp)           :: rtols = 1E-4_dp ! accuracy

!!! initialize kpp variables
      ! RTOL(:) = rtols   ! must be defined for ros3
      ! ATOL(:) = 1E-6_dp ! must be defined for ros3
!!! new rtol/atol for MOM chemistry
!!! in messy_mecca_kpp_initialize.f90:
       RTOL(:) = 1E-2_dp ! relative tolerance
       ATOL(:) = 1E1_dp  ! absolute tolerance

! Calculate photolysis rates [s^-1]

     !  write(6,*) 'cg o3 oh before',cgas(ind_o3),cgas(ind_oh)

       call photo_mbl(firstloop,daytime,lat_deg,press,jx)

       fcoj = fco

! CHAMBER STUFF
       IF (ICHAM .EQ. 1) THEN
         !adjust to measured jno2
         RF=jno2m/(MAX(jx(ip_NO2),1.E-14_dp))
         jx(ip_O1D)=RF*jx(ip_O1D)
         jx(ip_O3P)=RF*jx(ip_O3P)
         jx(ip_NO2)=RF*jx(ip_NO2)
         jx(ip_HNO3)=RF*jx(ip_HNO3)
         jx(ip_CHOH)=RF*jx(ip_CHOH)
         jx(ip_COH2)=RF*jx(ip_COH2)
         jx(ip_HONO)=RF*jx(ip_HONO)
         jx(ip_NOO2)=RF*jx(ip_NOO2)
         jx(ip_NO2O)=RF*jx(ip_NO2O) 
         jx(ip_N2O5)=RF*jx(ip_N2O5)
         jx(ip_H2O2)=RF*jx(ip_H2O2)
         jx(ip_CH3CO3H)=RF*jx(ip_CH3CO3H)
         jx(ip_CH3CHO)=RF*jx(ip_CH3CHO)
         jx(ip_CH3COCH3)=RF*jx(ip_CH3COCH3)
         jx(ip_MGLYOX)=RF*jx(ip_MGLYOX)
         jx(ip_MACR)=RF*jx(ip_MACR)
         jx(ip_MVK)=RF*jx(ip_MVK)
         jx(ip_GLYOX)=RF*jx(ip_GLYOX)
         jx(ip_HOCH2CHO)=RF*jx(ip_HOCH2CHO)
         !chamber source of HONO
         CALL chamber_src(DTIME,jx(ip_NO2 ),RH,temp,F_HONO,cgas )
         if (jno2m.LT.1.e-6) then
           fcoj=0.
         else
           fcoj=fco
         endif  
       ENDIF

!!! Integration of gas phase chemistry !!!

! Check rate constants before integration
!       do jp=1,nreact
!         write(6,*) 'r',jp,rconst(jp)
!       end do

! Main kpp call (see messy_mecca_box.f90, kpp_integrate)

       dt = DTIME
       ! call kpp fill routines for data transfer
       CALL fill_temp(status, SPREAD(temp,1,NBL))
       CALL fill_cair(status, SPREAD(cair,1,NBL))
       CALL fill_press(status, SPREAD(press,1,NBL))
       !IF (REQ_MCFCT) THEN
       !  CALL fill_mcexp(status, SPREAD(mcexp,1,NBL))
       !ENDIF
       ! heterogeneous reactions:
       !dummy_khet_St(:) = 0.
       !dummy_khet_Tr(:) = 0.
       !CALL fill_khet_Tr(status, SPREAD(dummy_khet_Tr,1,NBL))
       !CALL fill_khet_St(status, SPREAD(dummy_khet_St,1,NBL))
       CALL fill_jx(status, SPREAD(jx,1,NBL))
       CALL fill_lwc(status, SPREAD(lwc,1,NBL))
       CALL fill_cvfac(status, SPREAD(cvfac,1,NBL))
       CALL fill_xaer(status, SPREAD(xaer,1,NBL))
       CALL fill_k_exf(status, SPREAD(k_exf,1,NBL))
       CALL fill_k_exb(status, SPREAD(k_exb,1,NBL))
       CALL fill_k_exf_N2O5(status, SPREAD(k_exf_N2O5,1,NBL))
       CALL fill_k_exf_ClNO3(status, SPREAD(k_exf_ClNO3,1,NBL))
       CALL fill_k_exf_BrNO3(status, SPREAD(k_exf_BrNO3,1,NBL))


       !IF (l_tag) CALL mecca_tag_preprocess

! call kpp_integrate
       ! update_rconst is called inside kpp_integrate
       CALL check_range('before kpp:',cgas(:),model_time,cair)
       cgas(:) = MAX(cgas(:),0._dp) ! set negative values to zero
       cbl = SPREAD(cgas,1,NBL) ! add one dummy dimension
       CALL kpp_integrate(dt,cbl)  ! main kpp call
       cgas = cbl(1,:)          ! remove the dummy dimension
       CALL check_range('after kpp: ',cgas(:),model_time,cair)

       !IF (l_tag) CALL mecca_tag_postprocess
    
! End chemistry integration

       do jo = 1,IP_MAX
       ! write(6,*) 'jrate',jo,jx(jo)
        jrnew(jo) = jx(jo)
       end do

      ! write(6,*) 'cg o3 oh after',cgas(ind_o3),cgas(ind_oh)

! Avoid infinitesimal organic vapour concentrations 
       c(ind_BSOV) = max(c(ind_BSOV),1.e-32_dp)
       c(ind_BLOV) = max(c(ind_BLOV),1.e-32_dp)
       c(ind_BELV) = max(c(ind_BELV),1.e-32_dp)
       c(ind_ASOV) = max(c(ind_ASOV),1.e-32_dp)
       c(ind_ALOV) = max(c(ind_ALOV),1.e-32_dp)
       c(ind_AELV) = max(c(ind_AELV),1.e-32_dp)
       c(ind_PIOV) = max(c(ind_PIOV),1.e-32_dp)
       c(ind_PSOV) = max(c(ind_PSOV),1.e-32_dp)
       c(ind_PELV) = max(c(ind_PELV),1.e-32_dp)


  end subroutine chemistry_solver


  subroutine emis_drydep(time_step_len,emis,vdry,rain,zmbl,cgas)
    !----------------------------------------------------------------------
    !     
    ! Emission and deposition are currently calculated with Euler forward.
    !  This should be changed if numerical problems occur.
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculation of emission and dry deposition of Cgas
    !      from input: emission rate, v_dry, BL height
    !      input/output the gas phase concentration
    !
    !      interface
    !      ---------
    !
    !        input:
    !          emis         emission rate                     [molec(g)/(cm2*s)]
    !          vdry         dry deposition rate               [cm/s]
    !          time_step_len    time step length              [s]
    !          zmbl         mixing layer height               [m]    
    !          rain         rainfall rate                     [mm/h]
    !
    !      method
    !      ------
    !      interface to the MECCA chemistry solver
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

    implicit none

    ! input
    REAL( dp), intent(in)   :: emis(nspec),vdry(nspec)
    REAL( dp), intent(in)   :: time_step_len      ! in seconds
    REAL( dp), intent(in)   :: zmbl               ! in meter
    REAL( dp), intent(in)   :: rain               ! in mm/h
    !
    REAL( dp),dimension(nspec), intent(in out) :: cgas

    REAL( dp) :: fct

    !zmbl    = 1000._dp           ! height of boundary layer [m]
    ! emissions rates [molec(g)/(cm2*s)] conv. to [molec(g)/cm3]
    fct = time_step_len / (100._dp * zmbl)
    cgas(ind_O3)       = cgas(ind_O3)      +emis(ind_O3     )* fct
    cgas(ind_NO)       = cgas(ind_NO)      +emis(ind_NO     )* fct
    cgas(ind_NO2)      = cgas(ind_NO2)     +emis(ind_NO2    )* fct
    cgas(ind_CO)       = cgas(ind_CO)      +emis(ind_CO     )* fct
    cgas(ind_NH3)      = cgas(ind_NH3)     +emis(ind_NH3    )* fct
    cgas(ind_DMS)      = cgas(ind_DMS)     +emis(ind_DMS    )* fct
    cgas(ind_SO2)      = cgas(ind_SO2)     +emis(ind_SO2    )* fct
    cgas(ind_HCHO)     = cgas(ind_HCHO)    +emis(ind_HCHO   )* fct
    cgas(ind_H2O2)     = cgas(ind_H2O2)    +emis(ind_H2O2   )* fct
    cgas(ind_HNO3)     = cgas(ind_HNO3)    +emis(ind_HNO3   )* fct    
    cgas(ind_HONO)     = cgas(ind_HONO)    +emis(ind_HONO   )* fct 
    cgas(ind_H2SO4)    = cgas(ind_H2SO4)   +emis(ind_H2SO4  )* fct
    cgas(ind_CH4)      = cgas(ind_CH4)     +emis(ind_CH4    )* fct
    cgas(ind_HCOOH)    = cgas(ind_HCOOH)   +emis(ind_HCOOH  )* fct
    cgas(ind_CH3OOH)   = cgas(ind_CH3OOH)  +emis(ind_CH3OOH )* fct
    cgas(ind_PAN)      = cgas(ind_PAN)     +emis(ind_PAN    )* fct
    cgas(ind_N2O5)     = cgas(ind_N2O5)    +emis(ind_N2O5   )* fct
    cgas(ind_HNO4)     = cgas(ind_HNO4)    +emis(ind_HNO4   )* fct
    cgas(ind_NO3)      = cgas(ind_NO3)     +emis(ind_NO3    )* fct
    cgas(ind_HCl)      = cgas(ind_HCl)     +emis(ind_HCl    )* fct
    cgas(ind_CH3I)     = cgas(ind_CH3I)    +emis(ind_CH3I   )* fct
    cgas(ind_CH3SO3H)  = cgas(ind_CH3SO3H) +emis(ind_CH3SO3H)* fct
    cgas(ind_SO3)      = cgas(ind_SO3)     +emis(ind_SO3    )* fct
    
    ! VOC
    cgas(ind_C2H6)     = cgas(ind_C2H6)    +emis(ind_C2H6   )* fct
    cgas(ind_C3H8)     = cgas(ind_C3H8)    +emis(ind_C3H8   )* fct
    cgas(ind_C5H8)     = cgas(ind_C5H8)    +emis(ind_C5H8   )* fct
    cgas(ind_TOLUENE)  = cgas(ind_TOLUENE) +emis(ind_TOLUENE)* fct
    cgas(ind_C2H4)     = cgas(ind_C2H4)    +emis(ind_C2H4   )* fct
    cgas(ind_C3H6)     = cgas(ind_C3H6)    +emis(ind_C3H6   )* fct
    cgas(ind_NC4H10)   = cgas(ind_NC4H10)  +emis(ind_NC4H10 )* fct
    cgas(ind_MVK)      = cgas(ind_MVK)     +emis(ind_MVK    )* fct
    cgas(ind_MEK)      = cgas(ind_MEK)     +emis(ind_MEK    )* fct
    cgas(ind_LTMB)     = cgas(ind_LTMB)    +emis(ind_LTMB   )* fct
    cgas(ind_TME)      = cgas(ind_TME)     +emis(ind_TME    )* fct
    cgas(ind_IPN)      = cgas(ind_IPN)     +emis(ind_IPN    )* fct
    cgas(ind_APINENE)  = cgas(ind_APINENE) +emis(ind_APINENE)* fct
    cgas(ind_CHEX)     = cgas(ind_CHEX)    +emis(ind_CHEX   )* fct
    cgas(ind_BUT1ENE)  = cgas(ind_BUT1ENE) +emis(ind_BUT1ENE)* fct
    cgas(ind_LXYL)     = cgas(ind_LXYL)    +emis(ind_LXYL   )* fct
    cgas(ind_BPINENE)  = cgas(ind_BPINENE) +emis(ind_BPINENE)* fct
    cgas(ind_CAMPHENE) = cgas(ind_CAMPHENE)+emis(ind_CAMPHENE)* fct
    cgas(ind_CARENE)   = cgas(ind_CARENE)  +emis(ind_CARENE )* fct
    cgas(ind_SABINENE) = cgas(ind_SABINENE)+emis(ind_SABINENE)* fct

    ! AMINES
    cgas(ind_MEA)      = cgas(ind_MEA)     +emis(ind_MEA    )* fct  
    cgas(ind_DEA)      = cgas(ind_DEA)     +emis(ind_DEA    )* fct  
    cgas(ind_TEA)      = cgas(ind_TEA)     +emis(ind_TEA    )* fct  
    cgas(ind_MMA)      = cgas(ind_MMA)     +emis(ind_MMA    )* fct  
    cgas(ind_DMA)      = cgas(ind_DMA)     +emis(ind_DMA    )* fct  
    cgas(ind_TMA)      = cgas(ind_TMA)     +emis(ind_TMA    )* fct  
    cgas(ind_CH2NCH3)  = cgas(ind_CH2NCH3) +emis(ind_CH2NCH3)* fct
    cgas(ind_AMP)      = cgas(ind_AMP)     +emis(ind_AMP    )* fct
    
    ! SOA        
    cgas(ind_BSOV)     = cgas(ind_BSOV)    +emis(ind_BSOV   )* fct
    cgas(ind_BLOV)     = cgas(ind_BLOV)    +emis(ind_BLOV   )* fct   
    cgas(ind_BELV)     = cgas(ind_BELV)    +emis(ind_BELV   )* fct   
    cgas(ind_ASOV)     = cgas(ind_ASOV)    +emis(ind_BSOV   )* fct
    cgas(ind_ALOV)     = cgas(ind_ALOV)    +emis(ind_BLOV   )* fct   
    cgas(ind_AELV)     = cgas(ind_AELV)    +emis(ind_BELV   )* fct
    cgas(ind_PIOV)     = cgas(ind_PIOV)    +emis(ind_PIOV   )* fct
    cgas(ind_PSOV)     = cgas(ind_PSOV)    +emis(ind_PSOV   )* fct   
    cgas(ind_PELV)     = cgas(ind_PELV)    +emis(ind_PELV   )* fct


    ! deposition velocities [cm/s]
    fct = time_step_len / (zmbl * 100.)
    cgas(ind_O3)      = (1._dp-fct*vdry(ind_O3))      *cgas(ind_O3)
    cgas(ind_H2O2)    = (1._dp-fct*vdry(ind_H2O2))    *cgas(ind_H2O2)
    cgas(ind_NH3)     = (1._dp-fct*vdry(ind_NH3))     *cgas(ind_NH3)
    cgas(ind_NO2)     = (1._dp-fct*vdry(ind_NO2))     *cgas(ind_NO2)
    cgas(ind_N2O5)    = (1._dp-fct*vdry(ind_N2O5))    *cgas(ind_N2O5)
    cgas(ind_HNO3)    = (1._dp-fct*vdry(ind_HNO3))    *cgas(ind_HNO3)
    cgas(ind_CH3OOH)  = (1._dp-fct*vdry(ind_CH3OOH))  *cgas(ind_CH3OOH)
    cgas(ind_HCHO)    = (1._dp-fct*vdry(ind_HCHO))    *cgas(ind_HCHO)
    cgas(ind_HCOOH)   = (1._dp-fct*vdry(ind_HCOOH))   *cgas(ind_HCOOH)
    cgas(ind_HONO)    = (1._dp-fct*vdry(ind_HONO))    *cgas(ind_HONO)
    cgas(ind_HCl)     = (1._dp-fct*vdry(ind_HCl))     *cgas(ind_HCl)
    cgas(ind_CH3I)    = (1._dp-fct*vdry(ind_CH3I))    *cgas(ind_CH3I)
    cgas(ind_SO2)     = (1._dp-fct*vdry(ind_SO2))     *cgas(ind_SO2)
    cgas(ind_H2SO4)   = (1._dp-fct*vdry(ind_H2SO4))   *cgas(ind_H2SO4)
    cgas(ind_CH3SO3H) = (1._dp-fct*vdry(ind_CH3SO3H)) *cgas(ind_CH3SO3H)
    !cgas(ind_DMSO)    = (1._dp-fct*vdry(ind_DMSO)) * cgas(ind_DMSO)
    cgas(ind_DMS)     = (1._dp-fct*vdry(ind_DMS))     *cgas(ind_DMS)
    cgas(ind_SO3)     = (1._dp-fct*vdry(ind_SO3))     *cgas(ind_SO3)
    cgas(ind_PAN)     = (1._dp-fct*vdry(ind_PAN))     *cgas(ind_PAN)
    cgas(ind_HNO4)    = (1._dp-fct*vdry(ind_HNO4))    *cgas(ind_HNO4)
    cgas(ind_NO3)     = (1._dp-fct*vdry(ind_NO3))     *cgas(ind_NO3)
    cgas(ind_CH4)     = (1._dp-fct*vdry(ind_CH4))     *cgas(ind_CH4)
    cgas(ind_CO)      = (1._dp-fct*vdry(ind_CO))      *cgas(ind_CO)
 
    cgas(ind_IPN)     = (1._dp-fct*vdry(ind_IPN))     *cgas(ind_IPN)    
    cgas(ind_MMA)     = (1._dp-fct*vdry(ind_MMA))     *cgas(ind_MMA)
    cgas(ind_DMA)     = (1._dp-fct*vdry(ind_DMA))     *cgas(ind_DMA)
    cgas(ind_TMA)     = (1._dp-fct*vdry(ind_TMA))     *cgas(ind_TMA)
    cgas(ind_MEA)     = (1._dp-fct*vdry(ind_MEA))     *cgas(ind_MEA)
    cgas(ind_DEA)     = (1._dp-fct*vdry(ind_DEA))     *cgas(ind_DEA)
    cgas(ind_TEA)     = (1._dp-fct*vdry(ind_TEA))     *cgas(ind_TEA)
    cgas(ind_CH2NCH3) = (1._dp-fct*vdry(ind_CH2NCH3)) *cgas(ind_CH2NCH3)
    cgas(ind_AMP)     = (1._dp-fct*vdry(ind_AMP))     *cgas(ind_AMP)

    cgas(ind_C2H6)    = (1._dp-fct*vdry(ind_C2H6))    *cgas(ind_C2H6)   
    cgas(ind_C2H4)    = (1._dp-fct*vdry(ind_C2H4))    *cgas(ind_C2H4)   
    cgas(ind_C3H8)    = (1._dp-fct*vdry(ind_C3H8))    *cgas(ind_C3H8)
    cgas(ind_C3H6)    = (1._dp-fct*vdry(ind_C3H6))    *cgas(ind_C3H6)
    cgas(ind_NC4H10)  = (1._dp-fct*vdry(ind_NC4H10))  *cgas(ind_NC4H10)
    cgas(ind_MVK)     = (1._dp-fct*vdry(ind_MVK))     *cgas(ind_MVK)    
    cgas(ind_MEK)     = (1._dp-fct*vdry(ind_MEK))     *cgas(ind_MEK)
    cgas(ind_C5H8)    = (1._dp-fct*vdry(ind_C5H8))    *cgas(ind_C5H8)
    cgas(ind_LTMB)    = (1._dp-fct*vdry(ind_LTMB))    *cgas(ind_LTMB)
    cgas(ind_TOLUENE) = (1._dp-fct*vdry(ind_TOLUENE)) *cgas(ind_TOLUENE)
    cgas(ind_TME)     = (1._dp-fct*vdry(ind_TME))     *cgas(ind_TME)  
    cgas(ind_APINENE) = (1._dp-fct*vdry(ind_APINENE)) *cgas(ind_APINENE)
    cgas(ind_CHEX)    = (1._dp-fct*vdry(ind_CHEX))    *cgas(ind_CHEX)
    cgas(ind_BUT1ENE) = (1._dp-fct*vdry(ind_BUT1ENE)) *cgas(ind_BUT1ENE)
    cgas(ind_LXYL)    = (1._dp-fct*vdry(ind_LXYL))    *cgas(ind_LXYL)
    cgas(ind_BPINENE) = (1._dp-fct*vdry(ind_BPINENE)) *cgas(ind_BPINENE)
    cgas(ind_CAMPHENE)= (1._dp-fct*vdry(ind_CAMPHENE))*cgas(ind_CAMPHENE)
    cgas(ind_CARENE)  = (1._dp-fct*vdry(ind_CARENE))  *cgas(ind_CARENE)
    cgas(ind_SABINENE)= (1._dp-fct*vdry(ind_SABINENE))*cgas(ind_SABINENE)

    cgas(ind_BSOV)    = (1._dp-fct*vdry(ind_BSOV))    *cgas(ind_BSOV) 
    cgas(ind_BLOV)    = (1._dp-fct*vdry(ind_BLOV))    *cgas(ind_BLOV)    
    cgas(ind_BELV)    = (1._dp-fct*vdry(ind_BELV))    *cgas(ind_BELV)   
    cgas(ind_ASOV)    = (1._dp-fct*vdry(ind_ASOV))    *cgas(ind_ASOV) 
    cgas(ind_ALOV)    = (1._dp-fct*vdry(ind_ALOV))    *cgas(ind_ALOV)    
    cgas(ind_AELV)    = (1._dp-fct*vdry(ind_AELV))    *cgas(ind_AELV) 
    cgas(ind_PIOV)    = (1._dp-fct*vdry(ind_PIOV))    *cgas(ind_PIOV) 
    cgas(ind_PSOV)    = (1._dp-fct*vdry(ind_PSOV))    *cgas(ind_PSOV)  
    cgas(ind_PELV)    = (1._dp-fct*vdry(ind_PELV))    *cgas(ind_PELV) 


! Calculate wet deposition of SO2, cloud thickness 200m
    cgas(ind_SO2)     = (1.-(time_step_len*6.9E-03*rain/200.))*cgas(ind_SO2)

  end subroutine emis_drydep


  subroutine chamber_dil(time_step_len,kdil,cgas)
    !----------------------------------------------------------------------
    !     
    !      Calculation of dilution in a chamber
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculation of dilution of Cgas
    !      from input: first order dilution rate constant
    !      input/output the gas phase concentration
    !
    !      interface
    !      ---------
    !
    !        input:
    !          kdil         chamber dilution rate             [1/s]
    !          time_step_len    time step length              [s]
    !
    !      method
    !      ------
    !      interface to the MECCA chemistry solver
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

    implicit none

    ! input
    REAL( dp), intent(in)   :: kdil               ! in 1/s
    REAL( dp), intent(in)   :: time_step_len      ! in seconds

    !
    REAL( dp),dimension(nspec), intent(in out) :: cgas

    REAL( dp) :: fct

    fct = time_step_len

    ! Amines
    cgas(ind_MEA)       = (1._dp-fct*kdil) * cgas(ind_MEA)
    cgas(ind_H2NCHO)    = (1._dp-fct*kdil) * cgas(ind_H2NCHO)
    cgas(ind_MMA)       = (1._dp-fct*kdil) * cgas(ind_MMA)
    cgas(ind_DMA)       = (1._dp-fct*kdil) * cgas(ind_DMA)
    cgas(ind_TMA)       = (1._dp-fct*kdil) * cgas(ind_TMA)
    cgas(ind_CH3NHCHO)  = (1._dp-fct*kdil) * cgas(ind_CH3NHCHO)
    cgas(ind_TMADF)     = (1._dp-fct*kdil) * cgas(ind_TMADF)
    cgas(ind_AMP)       = (1._dp-fct*kdil) * cgas(ind_AMP)
    cgas(ind_AMPNNO2)   = (1._dp-fct*kdil) * cgas(ind_AMPNNO2)
    cgas(ind_NAMP)      = (1._dp-fct*kdil) * cgas(ind_NAMP)
    cgas(ind_H2NCOCH3)  = (1._dp-fct*kdil) * cgas(ind_H2NCOCH3)
    cgas(ind_H2NCOCH2OH)= (1._dp-fct*kdil) * cgas(ind_H2NCOCH2OH)
    cgas(ind_CH2NCH3)   = (1._dp-fct*kdil) * cgas(ind_CH2NCH3)
    cgas(ind_DMNNO2)    = (1._dp-fct*kdil) * cgas(ind_DMNNO2)
    cgas(ind_CH2NCH3)   = (1._dp-fct*kdil) * cgas(ind_CH2NCH3)

    ! VOC
    cgas(ind_LTMB)      = (1._dp-fct*kdil) * cgas(ind_LTMB)
    cgas(ind_TME)       = (1._dp-fct*kdil) * cgas(ind_TME)
    cgas(ind_IPN)       = (1._dp-fct*kdil) * cgas(ind_IPN)    
    cgas(ind_APINENE)   = (1._dp-fct*kdil) * cgas(ind_APINENE)
    cgas(ind_C5H8)      = (1._dp-fct*kdil) * cgas(ind_C5H8)    
    cgas(ind_CHEX)      = (1._dp-fct*kdil) * cgas(ind_CHEX)  
    cgas(ind_BUT1ENE)   = (1._dp-fct*kdil) * cgas(ind_BUT1ENE)
    cgas(ind_LXYL)      = (1._dp-fct*kdil) * cgas(ind_LXYL)
    cgas(ind_BPINENE)   = (1._dp-fct*kdil) * cgas(ind_BPINENE)
    cgas(ind_CAMPHENE)  = (1._dp-fct*kdil) * cgas(ind_CAMPHENE)
    cgas(ind_CARENE)    = (1._dp-fct*kdil) * cgas(ind_CARENE)
    cgas(ind_SABINENE)  = (1._dp-fct*kdil) * cgas(ind_SABINENE)

    ! Major products
    cgas(ind_MGLYOX)    = (1._dp-fct*kdil) * cgas(ind_MGLYOX) 
    cgas(ind_ACETOL)    = (1._dp-fct*kdil) * cgas(ind_ACETOL) 
    cgas(ind_GLYOX)     = (1._dp-fct*kdil) * cgas(ind_GLYOX)         
    cgas(ind_HCHO)      = (1._dp-fct*kdil) * cgas(ind_HCHO)
    cgas(ind_MVK)       = (1._dp-fct*kdil) * cgas(ind_MVK)
    cgas(ind_MACR)      = (1._dp-fct*kdil) * cgas(ind_MACR)
           
    cgas(ind_NO)        = (1._dp-fct*kdil) * cgas(ind_NO)
    cgas(ind_NO2)       = (1._dp-fct*kdil) * cgas(ind_NO2)
    cgas(ind_O3)        = (1._dp-fct*kdil) * cgas(ind_O3)

    ! SOA
    cgas(ind_BSOV)      = (1._dp-fct*kdil) * cgas(ind_BSOV)
    cgas(ind_BLOV)      = (1._dp-fct*kdil) * cgas(ind_BLOV)   
    cgas(ind_BELV)      = (1._dp-fct*kdil) * cgas(ind_BELV)  
    cgas(ind_ASOV)      = (1._dp-fct*kdil) * cgas(ind_ASOV) 
    cgas(ind_ALOV)      = (1._dp-fct*kdil) * cgas(ind_ALOV)    
    cgas(ind_AELV)      = (1._dp-fct*kdil) * cgas(ind_AELV) 
    cgas(ind_PIOV)      = (1._dp-fct*kdil) * cgas(ind_PIOV)
    cgas(ind_PSOV)      = (1._dp-fct*kdil) * cgas(ind_PSOV)  
    cgas(ind_PELV)      = (1._dp-fct*kdil) * cgas(ind_PELV)


  end subroutine chamber_dil


  subroutine chamber_loss(time_step_len,temp,press,Morg,Csat0,dhvap,    &
                          lmea,lno2,lhno3,lo3,vcham,scham,cwin,cgas)
    !----------------------------------------------------------------------
    !     
    !      Calculation of wall loss in a chamber
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculation of wall loss of Cgas
    !      as first order wall loss rate constant
    !      wall loss of SOA(gas) according to Zhang et al. (2014)
    !      for irreversible sticking of organic vapors
    !      from input: first order wall loss rate constant   [1/s]
    !                  vcham - chamber volume                [m3]
    !                  scham - total surface area of chamber [m2]
    !      input/output the gas phase concentration
    !
    !      interface
    !      ---------
    !
    !        input:
    !          lmea             wall loss rate                  [1/s]
    !          lno2             wall loss rate                  [1/s]
    !          lhno3            wall loss rate                  [1/s]
    !          lo3              wall loss rate                  [1/s]
    !          time_step_len    time step length                [s]
    !          temp             temperature                     [K]
    !          press            pressure                        [Pa]
    !          Morg             molecular weight SOA            [g/mol]
    !          Csat0            saturation concentration        [ug/m3]
    !          dhvap            enthalpy of vaporization        [J/mol]
    !          vcham            chamber volume                  [m3]
    !          scham            total surface area of chamber   [m2]
    !          cwin             equivalent wall organic aerosol [mg/m3]
    !
    !      method
    !      ------
    !      calculate Cgas change
    !
    !      reference
    !      ---------
    !      Wall loss of SOA(gas) compounds:
    !
    !      Zhang, X, Schwantes, R.H., McVay, R.C., Lignell, H.,
    !      Coggon, M.M., Flagan, R.C. and J.H. Seinfeld,
    !        Vapor wall deposition in Teflon chambers,
    !        Atmos. Chem. Phys., 15, 4197-4214,
    !        DOI:10.5194/acp-15-4197-2015, 2015.
    !
    !      Zhang, X., Cappa, C.D., Jathar, S.H., McVay, R.C.,
    !      Ensberg, J.J., Kleeman, M.J. and J.H. Seinfeld,
    !        Influence of vapor wall loss in laboratory chambers
    !        on yields of secondary organic aerosol,
    !        P. Natl. Acad. Sci. USA, 111, 5802-5807,
    !        DOI:10.1073/pnas.1404727111, 2014.
    !
    !      Crump, J.G. and J.H. Seinfeld, 
    !        Turbulent deposition and gravitational sedimentation
    !        of an aerosol in a vessel of arbitrary shape,
    !        J. Aerosol Sci., 12(5), 405-415, 1981.
    !
    !      Naumann, K.-H.,
    !        COSIMA - a computer program simulating the dynamics of 
    !        fractal aerosols
    !        Aerosol Science, 34, 1371-1397, 2003.
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

    implicit none

    ! input
    real( dp), intent(in)   :: temp,press
    real( dp), intent(in)   :: lmea,lno2,lhno3,lo3  ! in 1/s
    real( dp), intent(in)   :: vcham                ! in m3
    real( dp), intent(in)   :: scham                ! in m2
    real( dp), intent(in)   :: cwin                 ! in mg/m3
    real( dp), intent(in)   :: time_step_len        ! in seconds
    real( dp), dimension(NSOA), intent(in)     :: Morg
    real( dp), dimension(NSOA), intent(in)     :: Csat0
    real( dp), dimension(NSOA), intent(in)     :: dhvap

    ! in and out
    real( dp),dimension(nspec), intent(in out) :: cgas


    ! local
    real( dp), dimension(NSOA) :: cst
    real( dp), dimension(NSOA) :: aw
    real( dp), dimension(NSOA) :: lsoa
    real( dp), dimension(NSOA) :: wsoa
    real( dp), parameter       :: keddy = 40._dp         ! 1/s
    real( dp), parameter       :: sigma_A = 5.85         ! Angstroem  
    real( dp) :: fct
    real( dp) :: svrat
    real( dp) :: pres
    real( dp) :: DC0
    real( dp) :: CBAR
    real( dp) :: MWOC
    real( dp) :: logaw
    real( dp) :: dnom
    real( dp) :: cwinu
    integer   :: i


    fct = time_step_len

    !------------------------------------------------------------------
    !
    ! Approximation for trace gases' wall loss
    !
    !------------------------------------------------------------------

    cgas(ind_NO2)       = (1._dp-fct*lno2) * cgas(ind_NO2)
    cgas(ind_HNO3)      = (1._dp-fct*lhno3)* cgas(ind_HNO3)
    cgas(ind_O3)        = (1._dp-fct*lo3)  * cgas(ind_O3)
    ! Amines wall loss
    cgas(ind_MEA)       = (1._dp-fct*lmea) * cgas(ind_MEA)
    cgas(ind_MMA)       = (1._dp-fct*lmea) * cgas(ind_MMA)
    cgas(ind_DMA)       = (1._dp-fct*lmea) * cgas(ind_DMA)
    cgas(ind_TMA)       = (1._dp-fct*lmea) * cgas(ind_TMA)
    cgas(ind_CH2NCH3)   = (1._dp-fct*lmea) * cgas(ind_CH2NCH3)
    cgas(ind_AMP)       = (1._dp-fct*lmea) * cgas(ind_AMP)
    ! VOC wall loss same as for O3  
    cgas(ind_APINENE)   = (1._dp-fct*lo3)  * cgas(ind_APINENE) 
    cgas(ind_C5H8)      = (1._dp-fct*lo3)  * cgas(ind_C5H8) 
    cgas(ind_CHEX)      = (1._dp-fct*lo3)  * cgas(ind_CHEX) 
    cgas(ind_BUT1ENE)   = (1._dp-fct*lo3)  * cgas(ind_BUT1ENE)
    cgas(ind_LTMB)      = (1._dp-fct*lo3)  * cgas(ind_LTMB)
    cgas(ind_LXYL)      = (1._dp-fct*lo3)  * cgas(ind_LXYL)
    cgas(ind_BPINENE)   = (1._dp-fct*lo3)  * cgas(ind_BPINENE)
    cgas(ind_CAMPHENE)  = (1._dp-fct*lo3)  * cgas(ind_CAMPHENE)
    cgas(ind_CARENE)    = (1._dp-fct*lo3)  * cgas(ind_CARENE)
    cgas(ind_SABINENE)  = (1._dp-fct*lo3)  * cgas(ind_SABINENE)
    ! Major gas-phase products
    cgas(ind_MGLYOX)    = (1._dp-fct*lo3)  * cgas(ind_MGLYOX) 
    cgas(ind_ACETOL)    = (1._dp-fct*lo3)  * cgas(ind_ACETOL) 
    cgas(ind_GLYOX)     = (1._dp-fct*lo3)  * cgas(ind_GLYOX)         
    cgas(ind_HCHO)      = (1._dp-fct*lo3)  * cgas(ind_HCHO)
    cgas(ind_MVK)       = (1._dp-fct*lo3)  * cgas(ind_MVK)
    cgas(ind_MACR)      = (1._dp-fct*lo3)  * cgas(ind_MACR)

    !------------------------------------------------------------------
    !
    ! Approximation for SOA(gas) wall loss
    ! assuming deposition of the organic vapour to the
    ! chamber walls and treating it as first order loss rate,
    ! and evaporation from the wall using gas-to-wall partitioning
    ! Based on Zhang et al. (2014).
    !
    !------------------------------------------------------------------
    ! keddy is the coefficient of eddy diffusion in [1/s]
    ! introduced in Crump and Seinfeld (1981), which may be
    ! evaluated from the turbulent kinetic dissipation rate
    ! Here we evaluate keddy from the diffusion boundary layer
    ! thickness (delta) at the chamber wall.
    ! Equation 20 in Crump and Seinfeld (1981), they choose:
    !   delta = (3pi/4)*SQRT(DIFFCO/keddy)
    ! Equation 40 in Naumann (2003) expresses delta as:
    !   delta = k_D*(DIFFCO/D0)**a
    !   with k_D=0.005 m and a=0.274 and D0=1 m2/s
    ! Solving for keddy:
    !   keddy = (3/4pi)**2*DIFFCO/(k_D**2*DIFFCO**2a)
    ! DIFFCO is calculated in subroutine diffpar (gde_deposition)
    ! For now we use a constant value ke = 40 s
    ! that indicates homogeneous mixing of chamber air
    ! Variation of ke by factor 10 changes lsoa by only few %
    !------------------------------------------------------------------
    ! aw_i is the mass accommodation coefficient of 
    !    organic vapor onto the wall in [-]
    ! Zhang et al. (2015) give an empiric relation with C*
    ! in their Equation (23): 
    !      log10(aw_i)=-0.1919*log10C*_i - 6.32 
    ! Here we simply use C0(T) instead
    ! using C0 provided in organic.dat
    !------------------------------------------------------------------

    ! surface on volume ratio (S/V) in [1/m]
    svrat = scham/vcham

    ! use pressure in [atm]
    pres=press/101325._dp

    ! compute the wall loss rate k_w [1/s]
    do i=1,NSOA
! For Toluene: sigma_A = 5.85 Angstroem
! Ref: J.R. Li, R.J. Kuppler and H.C. Zhou, 
!      Chem. Soc. Rev.,38, 1477-1504, 2009.
! Value of D0 should be around 0.04 cm2/s to 0.05 cm2/s
       DC0    = molecdiff(Morg(i),sigma_A,temp,pres)
! Convert to m2/s
       DC0    = DC0 * 1.e-4_dp  
! Mean thermal speed [m/s]
! Molecular weight has to be in units kg/molec
       MWOC   = Morg(i)*1.661e-27
       CBAR   = SQRT(8._dp*k_B*temp/(pi*MWOC))
! Mass accommodation coefficient of organic vapor onto the wall
! Csat0 is from organic.dat [ug/m3]
       cst(i) = Csat0(i)*EXP((dhvap(i)/R_gas) *                 &
                ((1._dp/298._dp)-(1._dp/temp)))
       logaw  = -0.1919_dp*LOG10(cst(i)) - 6.32_dp
       aw(i)  = 10._dp**logaw
! Denominator of Eq. 22b in Zhang et al. (2015)
       dnom   = 1._dp + (pi/2._dp)*(  (aw(i)*CBAR)/(4._dp*SQRT(keddy*DC0)) )
       lsoa(i) = (svrat*0.25_sp*aw(i)*CBAR) / dnom
     ! print *,'lsoa in : ',i,aw(i),DC0,dnom,lsoa(i)  
    end do

    cgas(ind_BSOV)      = (1._dp-fct*lsoa(1)) * cgas(ind_BSOV)
    cgas(ind_BLOV)      = (1._dp-fct*lsoa(2)) * cgas(ind_BLOV)
    cgas(ind_BELV)      = (1._dp-fct*lsoa(3)) * cgas(ind_BELV)  
    cgas(ind_ASOV)      = (1._dp-fct*lsoa(4)) * cgas(ind_ASOV) 
    cgas(ind_ALOV)      = (1._dp-fct*lsoa(5)) * cgas(ind_ALOV)    
    cgas(ind_AELV)      = (1._dp-fct*lsoa(6)) * cgas(ind_AELV) 
    cgas(ind_PIOV)      = (1._dp-fct*lsoa(7)) * cgas(ind_PIOV)
    cgas(ind_PSOV)      = (1._dp-fct*lsoa(8)) * cgas(ind_PSOV)  
    cgas(ind_PELV)      = (1._dp-fct*lsoa(9)) * cgas(ind_PELV)


    ! compute the evaporation from the wall
    ! According to Zhang et al. (2014)
    ! using the equivalent wall organic matter (C_w)
    ! from incham.dat (typical value: 10 mg/m3)
    !   wsoa = lsoa/(K_w*C_w)        [1/s]
    ! K_w is the gas-wall partition coefficient,
    ! which is chosen equal to the gas-particle
    ! partition coefficient, with
    !    K_w = 1/C* (Donahue et al., 2006)
    ! here C* is approximated by C0(T)

    cwinu=cwin*1000._dp   ! convert mg/m3 to ug/m3
    do i=1,NSOA
       wsoa(i) = lsoa(i) / ( cwinu/(max(1.e-10,cst(i))) )
      ! print *,'wsoa in : ',i,cst(i),wsoa(i)  
    end do

    !print *,'BSOV after loss',cgas(ind_BSOV),cgas(ind_BLOV)
    cgas(ind_BSOV)      = cgas(ind_BSOV)   +wsoa(1)*fct *cgas(ind_BSOV)
    cgas(ind_BLOV)      = cgas(ind_BLOV)   +wsoa(2)*fct *cgas(ind_BLOV)
    cgas(ind_BELV)      = cgas(ind_BELV)   +wsoa(3)*fct *cgas(ind_BELV)
    cgas(ind_ASOV)      = cgas(ind_ASOV)   +wsoa(4)*fct *cgas(ind_ASOV)
    cgas(ind_ALOV)      = cgas(ind_ALOV)   +wsoa(5)*fct *cgas(ind_ALOV)
    cgas(ind_AELV)      = cgas(ind_AELV)   +wsoa(6)*fct *cgas(ind_AELV)
    cgas(ind_PIOV)      = cgas(ind_PIOV)   +wsoa(7)*fct *cgas(ind_PIOV)
    cgas(ind_PSOV)      = cgas(ind_PSOV)   +wsoa(8)*fct *cgas(ind_PSOV)
    cgas(ind_PELV)      = cgas(ind_PELV)   +wsoa(9)*fct *cgas(ind_PELV)
    !print *,'BSOV after evap',cgas(ind_BSOV),cgas(ind_BLOV)


  end subroutine chamber_loss


  subroutine chamber_src(time_step_len,jno2,RH,temp,fhono,cgas)
    !----------------------------------------------------------------------
    !     
    !      Calculation of wall source in a chamber
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Calculate HONO and HCHO wall source
    !      EUPHORE from Zador et al. 2006, JAC
    !      input: jno2, RH
    !      input/output the gas phase concentration
    !
    !      interface
    !      ---------
    !
    !        input:
    !          jno2             photoly. frequency NO2        [1/s]
    !          RH               rel. humidity                 [-]
    !          temp             temperature                   [K]
    !          fhono            scaling HONO source           [-]
    !          time_step_len    time step length              [s]
    !
    !      method
    !      ------
    !      calculate Cgas change
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

    implicit none

    ! input
    REAL( dp), intent(in)   :: jno2,RH,temp,fhono
    REAL( dp), intent(in)   :: time_step_len      ! in seconds
    !
    REAL( dp),dimension(nspec), intent(in out) :: cgas

    REAL( dp) :: fct,whonodry,whonowet,whcho

    ! whono,whcho in molec cm-3 s-1
    whonodry = jno2 *7.3E+21_dp*EXP(-8945._dp/temp)
    ! test: no temp dependence
    !whonodry = jno2 *7.3E+21*EXP(-8945./300.15)
    whonowet = whonodry
    IF (RH .gt. 0.05) whonowet = whonodry + (jno2*5.8E+08_dp*(RH*100._dp)**(0.36_dp))
    whcho    = jno2*3.1E+17_dp*EXP(-5686._dp/temp)    
    !write(6,*) jno2,whonodry,whonowet,whcho,cgas(ind_HONO)

    fct = time_step_len

    cgas(ind_HONO)     = cgas(ind_HONO)    +whonowet*fhono*fct
    cgas(ind_HCHO)     = cgas(ind_HCHO)    +whcho*fct

  end subroutine chamber_src



  subroutine check_range(infostring,conc,model_time,cair)

    !----------------------------------------------------------------------
    !     
    !      Check concentration range of solution
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Check concentration range of solution
    !
    !      interface
    !      ---------
    !
    !        input:
    !          conc             tracer concentration          [molec/cm^3]
    !          model_time       model time                    [s]
    !          cair             c(air) (wet)                  [mcl/cm^3]
    !
    !      method
    !      ------
    !      print a warning if a concentration is not in the correct range
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

    implicit none

      character(LEN=*), intent(in) :: infostring
      real(dp),         intent(in) :: conc(:) ! tracer concentration
      real(dp),         intent(in) :: model_time ! model time
      real(dp),         intent(in) :: cair
      integer :: jt


      INTRINSIC SIZE

      tracer_loop: DO jt=1,SIZE(conc)
        wrong_conc: IF ((conc(jt)<0._dp).OR.(conc(jt)>cair)) THEN
          WRITE(*,*) & 
            infostring , ', time =', model_time, &
            ', c =', conc(jt), ' mcl/cm3 for ', jt
          STOP
        ENDIF wrong_conc
      END DO tracer_loop

  end subroutine check_range
  

  end module gde_chem_gas
