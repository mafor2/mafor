! If you want to change this file, edit xpolymecca and
! template_messy_mecca_poly_si.f90, not smil/messy_mecca_poly_si.f90 !!!

MODULE messy_mecca_poly_si

  USE messy_main_constants_mem, ONLY: DP
  USE messy_mecca_kpp, ONLY: &
    REQ_HET, REQ_PHOTRAT, NSPEC_ori => NSPEC

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_poly, initialize, mr2c, c2mr, &
    kpp_integrate, mecca_new_tracer_gp, &
    fill_temp, fill_cair, fill_press, fill_jx, fill_khet_Tr, fill_khet_St
  PUBLIC :: REQ_HET, REQ_PHOTRAT
  INTEGER, PUBLIC, PARAMETER :: NMAXMECCA = 1
  INTEGER, PUBLIC, SAVE, DIMENSION(NMAXMECCA) :: NSPEC
  LOGICAL, PUBLIC, SAVE, DIMENSION(NMAXMECCA) :: l_fixed_step

  ! op_pj_20170904+
  ! for some reason (maybe the renaming of external routines at USE)
  ! this is required for the gfortran compiler
  INTERFACE kpp_integrate
     MODULE PROCEDURE kpp_integrate
  END INTERFACE kpp_integrate

  INTERFACE initialize
     MODULE PROCEDURE initialize
  END INTERFACE initialize

  INTERFACE fill_temp
     MODULE PROCEDURE fill_temp
  END INTERFACE fill_temp

  INTERFACE fill_cair
     MODULE PROCEDURE fill_cair
  END INTERFACE fill_cair

  INTERFACE fill_press
     MODULE PROCEDURE fill_press
  END INTERFACE fill_press

  INTERFACE fill_jx
     MODULE PROCEDURE fill_jx
  END INTERFACE fill_jx

  INTERFACE fill_khet_Tr
     MODULE PROCEDURE fill_khet_Tr
  END INTERFACE fill_khet_Tr

  INTERFACE fill_khet_St
     MODULE PROCEDURE fill_khet_St
  END INTERFACE fill_khet_St
  ! op_pj_20170904-

CONTAINS

  !***************************************************************************

  SUBROUTINE initialize_poly

    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,      ONLY: find_next_free_unit
    USE messy_mecca,           ONLY: modstr
    USE messy_mecca_kpp,       ONLY: &
      initialize_kpp_ctrl, l_fixed_step_ori => l_fixed_step, &
      nfsteps, icntrl, rcntrl, t_steps
    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'initialize_poly'
    INTEGER :: iou    ! I/O unit
    INTEGER :: status ! error status

    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      CALL initialize_kpp_ctrl(status, iou, modstr)
      IF (status /= 0) CALL error_bi('error in intitalize_kpp_ctrl',substr)
    ENDIF
    CALL p_bcast(l_fixed_step, p_io)
    CALL p_bcast(nfsteps,      p_io)
    CALL p_bcast(icntrl,       p_io)
    CALL p_bcast(rcntrl,       p_io)
    CALL p_bcast(t_steps,      p_io)
    NSPEC(1) = NSPEC_ori
    l_fixed_step(1) = l_fixed_step_ori

  END SUBROUTINE initialize_poly

  !***************************************************************************

  SUBROUTINE initialize
    USE messy_mecca_kpp, ONLY: initialize_ori => initialize
    CALL initialize_ori
  END SUBROUTINE initialize

  !***************************************************************************

  SUBROUTINE mr2c(jpm, zmrbc, c_air, conc, pind_H2O)
    USE messy_mecca_mem_si ! idt_*
    USE messy_mecca_kpp ! without ONLY to get all ind_*
    ! convert mixing ratio [mol/mol] to concentration [mcl/cc]
    INTEGER,  INTENT(IN)  :: jpm
    REAL(DP), INTENT(IN)  :: zmrbc(:,:), c_air(:)
    REAL(DP), INTENT(OUT) :: conc(:,:)
    INTEGER,  INTENT(OUT) :: pind_H2O

    pind_H2O = ind_H2O
    INCLUDE 'messy_mecca_mr2c_si.inc'

  END SUBROUTINE mr2c

  !***************************************************************************

  SUBROUTINE c2mr(jpm, zmrac, conc, c_air)
    USE messy_mecca_mem_si ! idt_*
    USE messy_mecca_kpp ! without ONLY to get all ind_*
    ! convert concentration [mcl/cc] to mixing ratio [mol/mol]
    INTEGER,  INTENT(IN)  :: jpm
    REAL(DP), INTENT(OUT) :: zmrac(:,:)
    REAL(DP), INTENT(IN)  :: conc(:,:), c_air(:)
    REAL(DP), DIMENSION(SIZE(c_air)) :: riac ! 1/c(air)

    riac(:) = 1._dp/c_air(:)
    INCLUDE 'messy_mecca_c2mr_si.inc'

  END SUBROUTINE c2mr

  !***************************************************************************

  SUBROUTINE kpp_integrate(jpm, time_step_len, conc, IERRF, &
    xNacc, xNrej, l_debug, PE)

    USE messy_mecca_kpp, ONLY: kpp_integrate_ori => kpp_integrate
    IMPLICIT NONE
    INTEGER, INTENT(IN)                             :: jpm
    REAL(DP),INTENT(IN)                             :: time_step_len
    REAL(DP),INTENT(INOUT), DIMENSION(:,:)          :: conc
    INTEGER, INTENT(OUT),   DIMENSION(:),  OPTIONAL :: IERRF
    INTEGER, INTENT(OUT),   DIMENSION(:),  OPTIONAL :: xNacc
    INTEGER, INTENT(OUT),   DIMENSION(:),  OPTIONAL :: xNrej
    LOGICAL, INTENT(IN),                   OPTIONAL :: l_debug
    INTEGER, INTENT(IN),                   OPTIONAL :: PE

    CALL kpp_integrate_ori(time_step_len, conc, IERRF=IERRF, &
      xNacc=xNacc, xNrej=xNrej, l_debug=l_debug, PE=PE)

  END SUBROUTINE kpp_integrate

  !***************************************************************************

  SUBROUTINE fill_temp(jpm, status, array)

    USE messy_mecca_kpp, ONLY: fill_temp_ori => fill_temp
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: jpm
    INTEGER, INTENT(OUT) :: status
    REAL(DP), INTENT(IN), DIMENSION(:) :: array

    CALL fill_temp_ori(status,array)

  END SUBROUTINE fill_temp

  !***************************************************************************

  SUBROUTINE fill_cair(jpm, status, array)

    USE messy_mecca_kpp, ONLY: fill_cair_ori => fill_cair
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: jpm
    INTEGER, INTENT(OUT) :: status
    REAL(DP), INTENT(IN), DIMENSION(:) :: array

    CALL fill_cair_ori(status,array)

  END SUBROUTINE fill_cair

  !***************************************************************************

  SUBROUTINE fill_press(jpm, status, array)

    USE messy_mecca_kpp, ONLY: fill_press_ori => fill_press
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: jpm
    INTEGER, INTENT(OUT) :: status
    REAL(DP), INTENT(IN), DIMENSION(:) :: array

    CALL fill_press_ori(status,array)

  END SUBROUTINE fill_press

  !***************************************************************************

  SUBROUTINE fill_jx(jpm, status, array)

    USE messy_mecca_kpp, ONLY: fill_jx_ori => fill_jx
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: jpm
    INTEGER, INTENT(OUT) :: status
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: array

    CALL fill_jx_ori(status,array)

  END SUBROUTINE fill_jx

  !***************************************************************************

  SUBROUTINE fill_khet_Tr(jpm, status, array)

    USE messy_mecca_kpp, ONLY: fill_khet_Tr_ori => fill_khet_Tr
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: jpm
    INTEGER, INTENT(OUT) :: status
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: array

    CALL fill_khet_Tr_ori(status,array)

  END SUBROUTINE fill_khet_Tr

  !***************************************************************************

  SUBROUTINE fill_khet_St(jpm, status, array)

    USE messy_mecca_kpp, ONLY: fill_khet_St_ori => fill_khet_St
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: jpm
    INTEGER, INTENT(OUT) :: status
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: array

    CALL fill_khet_St_ori(status,array)

  END SUBROUTINE fill_khet_St

  !***************************************************************************

  SUBROUTINE mecca_new_tracer_gp(c_pa_asm, i_pa_amode)

    ! ECHAM5/MESSy
    USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, LGTRSTR
    USE messy_main_tracer_bi,  ONLY: tracer_halt
    ! MESSy
    USE messy_mecca, ONLY: modstr
    USE messy_mecca_mem_si ! idt_*
    USE messy_main_tracer, ONLY: new_tracer, AIR, ON, OFF, &
      set_tracer, AEROSOL, AMOUNTFRACTION, MODAL, BIN, &
      I_ADVECT, I_CONVECT, I_VDIFF, I_WETDEP, I_DRYDEP, I_SEDI, &
      I_SCAV, I_MIX, I_FORCE_COL, I_INTEGRATE, I_TIMEFILTER, I_FORCE_INIT, &
      I_AEROSOL_METHOD, I_AEROSOL_MODE, I_AEROSOL_SOL, S_AEROSOL_MODEL, &
      R_MOLARMASS, R_PSS  , R_DRYREAC_SF, R_VINI, R_AEROSOL_DENSITY
    USE messy_main_constants_mem, ONLY: MH, MC, MN, MNa, MO, MS, MCl, MBr, &
      MI, MF, MHg
    ! When using a small (e.g. gas-phase only) chemistry mechanism,
    ! forcheck will say that AEROSOL, AMOUNTFRACTION, LGTRSTR, ON, and
    ! OFF are not used ("named constant not used"). Ignore this info and
    ! leave these variables in the ONLY lists!

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: c_pa_asm
    INTEGER,          INTENT(IN) :: i_pa_amode
    INTEGER :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_new_tracer_gp'
    CHARACTER(LEN=*), PARAMETER :: setname = GPTRSTR

    CALL start_message_bi(modstr, 'GRID POINT TRACER REQUEST', substr)
    ! The following INCLUDE statement contains all "CALL new_tracer" commands
    ! for all species that are used in the kpp chemistry scheme.
    ! Note that MECCA will never create a tracer for H2O.
    INCLUDE 'messy_mecca_trac_si.inc'
    CALL end_message_bi(modstr, 'GRID POINT TRACER REQUEST', substr)

  END SUBROUTINE mecca_new_tracer_gp

  !***************************************************************************

END MODULE messy_mecca_poly_si

!*****************************************************************************
