!*****************************************************************************
! Time-stamp: <2017-06-08 14:15:15 joec_pa>
!*****************************************************************************

! MECCA-AERO: calculate gas/aerosol transfer rates for MECCA

! Authors:
! Astrid Kerkweg, MPICH, Mainz, 2003-...
! Rolf Sander,    MPICH, Mainz, 2003-...

!*****************************************************************************

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program; if not, get it from:
! http://www.gnu.org/copyleft/gpl.html

!*****************************************************************************

MODULE messy_mecca_aero

  USE messy_main_constants_mem,  ONLY: atm2Pa, R_gas, T0, T0_INV
  USE messy_mecca_kpp ! ind_*, update_rconst, kpp_integrate, APN
                      ! initialize_indexarrays, SPC_NAMES, EQN_NAMES,
                      ! EQN_TAGS, NSPEC, NREACT

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PUBLIC :: submodstr='mecca_aero'

  REAL(DP) :: alpha_T0(0:NSPEC), alpha_Tdep(0:NSPEC)
  REAL(DP) :: Henry_T0(0:NSPEC), Henry_Tdep(0:NSPEC)
  REAL(DP) :: molar_mass(0:NSPEC)

  PUBLIC :: mecca_aero_init_gasaq
  PUBLIC :: mecca_aero_calc_k_ex
  PUBLIC :: mecca_aero_calc_k_ex_ocean
  PUBLIC :: mecca_aero_diag

CONTAINS

  !***************************************************************************

  SUBROUTINE mecca_aero_init_gasaq(l_print)

    USE messy_main_tracer
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, HLINE2
    IMPLICIT NONE
    INTEGER :: icp ! Index ChemProp
    INTEGER :: jn
    LOGICAL, INTENT(IN), OPTIONAL :: l_print

    IF (PRESENT(l_print)) THEN
      IF (l_print) THEN
      PRINT *, HLINE2
      PRINT *, "         Henry's law coefficients "// &
        "and accommodation coefficients"
      PRINT *, HLINE2
      PRINT *, 'species           Henry_T0 Henry_Tdep'// &
        '   alpha_T0 alpha_Tdep         M'
      PRINT *, '                   [M/atm]        [K]'// &
        '        [1]        [K]   [g/mol]'
      PRINT *, HLINE2
      ENDIF
    ENDIF
    DO jn = 1,NSPEC
      icp = get_chemprop_index(TRIM(SPC_NAMES(jn)))
      IF (icp /=0) THEN
        Henry_T0(jn)   = chemprop(icp)%cask_r(R_Henry_T0)
        Henry_Tdep(jn) = chemprop(icp)%cask_r(R_Henry_Tdep)
        alpha_T0(jn)   = chemprop(icp)%cask_r(R_alpha_T0)
        alpha_Tdep(jn) = chemprop(icp)%cask_r(R_alpha_Tdep)
        molar_mass(jn) = chemprop(icp)%cask_r(R_MOLARMASS) / 1000. ! [kg/mol]
        IF (PRESENT(l_print)) THEN
          IF (l_print) THEN
            WRITE(*,'(1X,A15)', ADVANCE='NO') SPC_NAMES(jn)
            IF (Henry_T0(jn)>=0.) THEN
              WRITE(*,'(ES11.2)', ADVANCE='NO') Henry_T0(jn)
            ELSE
              WRITE(*,'(A)', ADVANCE='NO') '   --------'
            ENDIF
            IF (Henry_Tdep(jn)>=0.) THEN
              WRITE(*,'(F11.0)', ADVANCE='NO') Henry_Tdep(jn)
            ELSE
              WRITE(*,'(A)', ADVANCE='NO') '     ------'
            ENDIF
            IF (alpha_T0(jn)>=0.) THEN
              WRITE(*,'(ES11.2)', ADVANCE='NO') alpha_T0(jn)
            ELSE
              WRITE(*,'(A)', ADVANCE='NO') '   --------'
            ENDIF
            IF (alpha_Tdep(jn)>=0.) THEN
              WRITE(*,'(F11.0)', ADVANCE='NO') alpha_Tdep(jn)
            ELSE
              WRITE(*,'(A)', ADVANCE='NO') '     ------'
            ENDIF
            WRITE(*,'(F10.2)') 1000.*molar_mass(jn) ! [g/mol]
          ENDIF
        ENDIF
      ENDIF
    ENDDO

  END SUBROUTINE mecca_aero_init_gasaq

  !***************************************************************************

  SUBROUTINE mecca_aero_calc_k_ex( &
    radius, temp, press, loghet, xaer, lwc, C, &
    k_exf, k_exb, k_exf_N2O5, k_exf_ClNO3, k_exf_BrNO3)

    ! calculation of gas-aqueous-phase exchange coefficients k_ex

    REAL(DP), INTENT(IN)  :: radius(:)
    REAL(DP), INTENT(IN)  :: temp
    REAL(DP), INTENT(IN)  :: press
    LOGICAL,  INTENT(IN)  :: loghet(:)
    REAL(DP), INTENT(IN)  :: xaer(:)
    REAL(DP), INTENT(IN)  :: lwc(:)
    REAL(DP), INTENT(IN)  :: C(:)
    REAL(DP), INTENT(OUT) :: k_exf(:,:)
    REAL(DP), INTENT(OUT) :: k_exb(:,:)
    REAL(DP), INTENT(OUT) :: k_exf_N2O5(:)
    REAL(DP), INTENT(OUT) :: k_exf_ClNO3(:)
    REAL(DP), INTENT(OUT) :: k_exf_BrNO3(:)

    REAL(DP), PARAMETER :: zrc_t = 1.E-9_DP ! radius threshold for kmt calc.

    INTEGER  :: jn, jb
    REAL(DP) :: zhetT, kmt, lambda, vmean, alpha, henry

    ! ------------------------------------------------------------------------

    ! lambda = mean free path
    ! 101325 Pa * 6.6E-8 m / 293.15 K = 2.28E-5 Pa*m/K is from
    ! equation (10-106) in Pruppacher and Klett = ref0064
    lambda = 2.28E-5_DP * temp / press

    ! calculate exchange rate coefficients:
    k_exf(:,:) = 0._DP
    k_exb(:,:) = 0._DP
    DO jn = 1,NSPEC ! loop over species
      IF (molar_mass(jn)<=0.) CYCLE
      ! mean molecular speed from Maxwell-Boltzmann distribution:
      ! vmean = SQRT(8*R_gas*T/(M*pi)) = sqrt(4.60138*T/M) [m/s]
      vmean = 4.60138_DP * SQRT(temp/molar_mass(jn))
      ! calculate temperature dependent Henry's law constants:
      henry = henry_T0(jn) * EXP(henry_Tdep(jn)*((1._DP/temp)-T0_INV))
      IF (alpha_T0(jn) > 0._DP) THEN
        ! calculate accommodation coefficients alpha:
        alpha = 1._DP / (1._DP+(1._DP/alpha_T0(jn)-1._DP) * &
          EXP(alpha_Tdep(jn)*((1._DP/T0)-(1._DP/temp))))
        DO jb = 1, APN ! loop over modes
          IF (radius(jb) >= zrc_t) THEN
            ! calculate mass transfer coefficients kmt after Schwarz, 1986
            ! k_mt=1/(r**2/(3*D_gas)+4*r/(3*v_mean*alpha))
            ! if D_gas=lambda*v_mean/3. then we can rewrite k_mt as:
            ! k_mt=vmean/(r**2/lambda+4*r/(3*alpha))
            kmt = vmean / (radius(jb)*(radius(jb)/lambda+4./(3._DP*alpha)))
            !kmt = kmt * 0.84_DP ! for M7 accumulation mode
            !kmt = kmt * 0.55_DP ! for M7 coarse mode
            ! reaction rate coefficients:
            k_exf(jb,jn) = xaer(jb) * kmt * lwc(jb) ! forward
            IF (henry > 0._DP) THEN
              k_exb(jb,jn) = xaer(jb) * kmt * atm2Pa &
                / (R_gas * 1.E3_DP * temp * henry) ! backward
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    ! ------------------------------------------------------------------------

    DO jb = 1, APN
      IF ((ind_H2O_a(jb)/=0).AND.(loghet(jb))) THEN
        zhetT = C(ind_H2O_a(jb))
        IF (ind_Clm_a(jb)/=0) zhetT = zhetT + 5.E2_DP * C(ind_Clm_a(jb))
        IF (ind_Brm_a(jb)/=0) zhetT = zhetT + 3.E5_DP * C(ind_Brm_a(jb))
        k_exf_N2O5(jb)  = k_exf(jb,IND_N2O5)  / zhetT
        k_exf_ClNO3(jb) = k_exf(jb,IND_ClNO3) / zhetT
        k_exf_BrNO3(jb) = k_exf(jb,IND_BrNO3) / zhetT
      ELSE
        k_exf_N2O5(jb)  = 0._DP
        k_exf_ClNO3(jb) = 0._DP
        k_exf_BrNO3(jb) = 0._DP
      ENDIF
    ENDDO

  END SUBROUTINE mecca_aero_calc_k_ex

  !***************************************************************************

  SUBROUTINE mecca_aero_calc_k_ex_ocean(xaer, radius, temp, zmix, &
    k_exf, k_exb)

    ! calculation of gas-aqueous-phase exchange coefficients for
    ! exchange with the ocean surface

    REAL(dp), INTENT(IN)  :: xaer(:)
    REAL(dp), INTENT(IN)  :: radius(:)
    REAL(DP), INTENT(IN)  :: temp
    REAL(dp), INTENT(IN)  :: zmix
    REAL(dp), INTENT(INOUT) :: k_exf(:,:)
    REAL(dp), INTENT(INOUT) :: k_exb(:,:)

    REAL(dp) :: K_tot(0:NSPEC) ! [1/s]
    REAL(DP) :: henry, yhenry
    INTEGER  :: jn, zkc

    CALL mecca_aero_K_tot(K_tot) ! define K_tot for all species

    DO zkc = 1, APN ! loop over modes/bins
      IF (radius(zkc)<0.) THEN ! negative dummy value for ocean mixed layer
        DO jn = 1, NSPEC   ! loop over species
          ! calculate temperature dependent Henry's law constants:
          henry = henry_T0(jn) * EXP(henry_Tdep(jn)*((1._DP/temp)-T0_INV))
          yhenry = atm2Pa / (R_gas * 1.E3_DP * temp * henry)
          ! reaction rate coefficients
          k_exf(zkc,jn) = xaer(zkc) * K_tot(jn)          / zmix ! forward
          k_exb(zkc,jn) = xaer(zkc) * K_tot(jn) * yhenry / zmix ! backward
        ENDDO
      ENDIF
    ENDDO

  CONTAINS

    ! ------------------------------------------------------------------------

    SUBROUTINE mecca_aero_K_tot(K_tot)

      REAL(dp), INTENT(OUT) :: K_tot(0:NSPEC)

      ! maybe this subroutine can be replaced later by obtaining K_tot
      ! from the airsea submodel...

      K_tot(:) = 1E-6_dp ! [1/s] default value

      ! special values for some species:
      K_tot(ind_O3) = 1E-7_dp

    END SUBROUTINE mecca_aero_K_tot

    ! ------------------------------------------------------------------------

  END SUBROUTINE mecca_aero_calc_k_ex_ocean

  ! --------------------------------------------------------------------------

  SUBROUTINE mecca_aero_diag(zmr_Brsalt, zmr_Brorg, zmr_BrSScap, &
    zmrac_Brsalt, zmrac_Brorg, zmrac_BrSScap)

    IMPLICIT NONE

    REAL(DP), INTENT(IN)    :: zmr_BrSScap,  zmrac_BrSScap
    REAL(DP), INTENT(IN)    :: zmr_Brsalt,   zmr_Brorg
    REAL(DP), INTENT(INOUT) :: zmrac_Brsalt, zmrac_Brorg

    REAL(DP) :: loss
    REAL(DP) :: Brsum, zmrcapdiff

    ! Define ratio between Brorg and Brsalt
    ! loss1 = Brsalt /(Brorg + Brsalt)
    Brsum =  zmr_Brsalt + zmr_Brorg
    IF (Brsum < TINY(Brsum) .OR. Brsum < 10.e-25) THEN
       loss = 1._DP
    ELSE
       loss =  zmr_Brsalt / Brsum
    ENDIF
    ! calculate sum of captured bromine 
    ! as the tracer BrSScap is accumlated the difference is the
    ! sum of captured Bromine in one timestep.
    zmrcapdiff   = zmrac_BrSScap - zmr_BrSScap
    ! depending on the ratio of the origines of Br atomes they are
    ! captured and loose there prior identitiy
    zmrac_Brsalt = zmrac_Brsalt - loss * zmrcapdiff
    zmrac_Brorg  = zmrac_Brorg  - (1._DP-loss) * zmrcapdiff

  END SUBROUTINE mecca_aero_diag

  !***************************************************************************

END MODULE messy_mecca_aero

!*****************************************************************************
