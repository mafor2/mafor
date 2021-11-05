! If you want to change this file, edit mecca/template_messy_mecca_kpp.f90 
! (and maybe kp4/templates/initialize_kpp_ctrl_template.f90 or
! kp4/src/create_kpp_module.C), but not smcl/messy_mecca_kpp.f90 !!!

! To survive the transport of variables from mecca to kpp via kp4, these
! rules must be followed (using mcexp as an example):
! - gas.eqn:
!   - declare variable inside a KPPPP_DIRECTIVE block:
!     !KPPPP_DIRECTIVE vector variable definition start
!       REAL(dp) :: mcexp(NREACT) ! Monte-Carlo factor
!     !KPPPP_DIRECTIVE vector variable definition end
! - messy_mecca_box.f90:
!   - declare the variable here:
!     REAL(DP) :: mcexp(NREACT)
! - template_messy_mecca_kpp.f90:
!   - add mcexp to the list of PRIVATE variables
!   - create a new fill subroutine:
!     SUBROUTINE fill_mcexp(status,array)
!     ...
!     END SUBROUTINE fill_mcexp
! - messy_mecca.f90:
!   - the new variable must not be used here! It can only be accessed
!     here as a SUBROUTINE (or FUNCTION) parameter from the SMIL level,
!     e.g.:
!     SUBROUTINE define_mcexp(mcexp)
!       REAL(DP), DIMENSION(:), INTENT(OUT) :: mcexp
!       ...
!     END SUBROUTINE define_mcexp

MODULE messy_mecca_kpp
  USE messy_mecca_kpp_function,   ONLY: A, CalcStoichNum
  USE messy_mecca_kpp_global      ! c, rconst, temp, press, dt, JX, khet_*,
                                  ! nsubsteps, logsteps, substep, method,
                                  ! rtol, atol, mcexp, ...
  USE messy_mecca_kpp_initialize, ONLY: initialize
  USE messy_mecca_kpp_integrator, ONLY: integrate, IERR_NAMES
  USE messy_mecca_kpp_parameters  ! ind_*
  USE messy_mecca_kpp_rates,      ONLY: update_rconst
  USE messy_mecca_kpp_monitor,    ONLY: SPC_NAMES, EQN_NAMES, EQN_TAGS
  USE messy_mecca_kpp_util,       ONLY: initialize_indexarrays

  IMPLICIT NONE
  SAVE

  ! the variables for which fill routines exist cannot be public:
  PRIVATE :: jx, cair, press, temp, C, cvfac, lwc, xaer, k_exf, k_exb, &
    k_exf_N2O5, k_exf_ClNO3, k_exf_BrNO3, mcexp
  !mz_hr_20130409+: fix variable leak of 'C' to mecca_init
  PRIVATE :: FIX, VAR
  !mz_hr_20130409- 

  ! start modified insertion from initialize_kpp_ctrl_template.f90

  ! NOTES:
  ! - ICNTRL RCNTRL are automatically defined by kpp
  ! - "USE messy_main_tools" is in Module_header of messy_mecca_kpp.f90

  ! FOR FIXED TIME STEP CONTROL
  ! ... max. number of fixed time steps (sum must be 1)
  INTEGER, PARAMETER         :: NMAXFIXSTEPS = 50
  ! ... switch for fixed time stepping
  LOGICAL, PUBLIC            :: l_fixed_step = .FALSE.
  INTEGER, PUBLIC            :: nfsteps = 1
  ! ... number of kpp control parameters
  INTEGER, PARAMETER, PUBLIC :: NKPPCTRL = 20
  !
  INTEGER,  DIMENSION(NKPPCTRL), PUBLIC     :: icntrl = 0
  REAL(DP), DIMENSION(NKPPCTRL), PUBLIC     :: rcntrl = 0.0_dp
  REAL(DP), DIMENSION(NMAXFIXSTEPS), PUBLIC :: t_steps = 0.0_dp

  ! END HEADER MODULE initialize_kpp_ctrl_template

CONTAINS

SUBROUTINE initialize_kpp_ctrl(status, iou, modstr)

  IMPLICIT NONE

  ! I/O
  INTEGER,          INTENT(OUT) :: status
  INTEGER,          INTENT(IN)  :: iou     ! logical I/O unit
  CHARACTER(LEN=*), INTENT(IN)  :: modstr  ! read <modstr>.nml

  ! LOCAL
  REAL(DP) :: tsum
  INTEGER  :: i

  CALL kpp_read_nml_ctrl(status, iou)
  IF (status /= 0) RETURN

  ! check fixed time steps
  tsum = 0.0_dp
  DO i=1, NMAXFIXSTEPS
     IF (t_steps(i) < TINY(0.0_DP)) EXIT
     tsum = tsum + t_steps(i)
  END DO

  nfsteps = i-1

  l_fixed_step = (nfsteps > 0) .AND. ( (tsum -1.0) < TINY(0.0_DP) )

  status = 0

CONTAINS

  SUBROUTINE kpp_read_nml_ctrl(status, iou)

    ! READ MECCA NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Astrid Kerkweg,  MPICH, June 2007
    !         Patrick Joeckel, MPICH, June 2007

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! LOCAL
    LOGICAL :: lex   ! file exists?
    INTEGER :: fstat ! file status
    CHARACTER(LEN=*), PARAMETER :: substr = 'kpp_read_nml_ctrl'

    NAMELIST /CTRL_KPP/ icntrl, rcntrl, t_steps

    ! INITIALIZE
    status = 1 ! error

    ! INPUT NAMELIST
    CALL read_nml_open(lex, substr, iou, 'CTRL_KPP', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_KPP, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_KPP', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    WRITE(*,*) 'solver-specific method:      icntrl(3) = ', icntrl(3)
    WRITE(*,*) 'max. number of kpp-substeps: icntrl(4) = ', icntrl(4)

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! no error

  END SUBROUTINE kpp_read_nml_ctrl

END SUBROUTINE initialize_kpp_ctrl

SUBROUTINE error_output(C,ierr,PE)

  USE messy_main_tools,      ONLY: str

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ierr
  INTEGER, INTENT(IN) :: PE
  REAL(dp), DIMENSION(:),INTENT(IN) :: C

  INTEGER,SAVE :: NUM =0
  INTEGER :: iou
  INTEGER :: i

  CHARACTER(LEN=250)  :: filename
  CHARACTER(LEN=1000) :: strnum
  CHARACTER(LEN=1000) :: strPE
  CHARACTER(256)      :: info

  LOGICAL             :: opened
  IF (ierr >= 0) RETURN

  NUM = NUM +1

  strnum=str(NUM)
  strPE=str(PE)

  WRITE(filename,*) 'mecca_PE'//TRIM(STRPE)//'_'//TRIM(STRNUM)//'.txt'

  iou = 0
  DO i=100,200
     INQUIRE(unit=i,opened=opened)
     IF (.NOT.opened) THEN
        iou = i
        EXIT
     END IF
  END DO

  IF (iou==0) THEN
     WRITE(info,'(a,i2.2,a,i2.2,a)') &
          'No unit in range < 100 : 200 > free.'
     RETURN
  ENDIF

  OPEN(IOU,FILE =TRIM(ADJUSTL(filename))&
       ,STATUS='NEW',ACTION= 'WRITE')

  ! ERROR STATUS
  WRITE(IOU,*) 'KPP ERRORSTATUS IERR:'
  WRITE(IOU,*) IERR

  ! SPECIES NAME
  WRITE(IOU,*) 'SELECTED MECHANISM'
  WRITE(IOU,*) WANTED
  WRITE(IOU,*)

  WRITE(IOU,*) 'NUMBER OF SPECIES:'
  WRITE(IOU,*) NSPEC
  WRITE(IOU,*)
  WRITE(IOU,*)  'NAMES OF SPECIES:'
  DO i=1,NSPEC
  WRITE(IOU,*)  SPC_NAMES(i)
  ENDDO
  WRITE(IOU,*)

  ! CONCENTRATIONS
  WRITE(IOU,*) 'concentrations (molec/cm^3(air) before KPP'
  DO i=1,NSPEC
     WRITE(IOU,*) SPC_NAMES(i), C(i)
  ENDDO
  WRITE(IOU,*)

  ! rates
  WRITE(IOU,*) 'rate contants'
  DO i=1,NREACT
     WRITE(IOU,*) RCONST(i)
  ENDDO
  WRITE(IOU,*)
  CLOSE(IOU)

END SUBROUTINE error_output

! end modified insertion from initialize_kpp_ctrl_template.f90

! start modified insertion from kpp_integrate from kp4/src/create_kpp_module.C

SUBROUTINE kpp_integrate (time_step,Conc,ierrf,xNacc,xNrej,istatus,&
     l_debug,PE)

  IMPLICIT NONE

  REAL(dp), INTENT(IN)                    :: time_step
  REAL(dp), INTENT(INOUT), DIMENSION(:,:) :: Conc
  INTEGER,  INTENT(OUT),   OPTIONAL       :: ierrf(:)
  INTEGER,  INTENT(OUT),   OPTIONAL       :: xNacc(:)
  INTEGER,  INTENT(OUT),   OPTIONAL       :: xNrej(:)
  INTEGER,  INTENT(INOUT), OPTIONAL       :: istatus(:)
  INTEGER,  INTENT(IN),    OPTIONAL       :: PE         ! mz_ak_20071206
  LOGICAL,  INTENT(IN),    OPTIONAL       :: l_debug    ! mz_ak_20071206

  ! this subroutine only treats the first (is=1) element:
  INTEGER, PARAMETER :: is = 1

  REAL(dp)               :: dt
  INTEGER, DIMENSION(20) :: istatus_u
  INTEGER                :: ierr_u

  IF (PRESENT (istatus)) istatus = 0

    C(:) = Conc(is,:)

    ! activate next line ONLY to check if all rate coefficients
    ! contain an uncertainty factor for Monte-Carlo simulations:
    ! CALL montecarlo_check

    CALL update_rconst

    dt = time_step

    ! integrate from t=0 to t=dt
    CALL integrate(0._dp,dt,icntrl,rcntrl,istatus_u=istatus_u,ierr_u=ierr_u)

    ! mz_ak_20071206+
    IF (PRESENT(l_debug) .AND. PRESENT(PE)) THEN
       IF (l_debug) CALL error_output(Conc(is,:),ierr_u,PE)
    ENDIF
    ! mz_ak_20071206-

    Conc(is,:) = C(:)

    ! Return Diagnostic Information

    IF(PRESENT(ierrf)) ierrf(is) = IERR_U
    IF(PRESENT(xNacc)) xNacc(is) = istatus_u(4)
    IF(PRESENT(xNrej)) xNrej(is) = istatus_u(5)

    IF (PRESENT (istatus)) THEN
      istatus(1:8) = istatus(1:8) + istatus_u(1:8)
    END IF

  CONTAINS

    !-------------------------------------------------------------------------

    SUBROUTINE montecarlo_check

      USE messy_main_constants_mem,   ONLY: HLINE2, TINY_DP

      INTEGER :: i
      REAL(DP) :: RCONST0(NREACT)

      ! jx(:), khet* and C(:) must be >0 to check reaction rate coefficients:
      jx(:)      = 1e-5
      khet_St(:) = 1e-5
      khet_Tr(:) = 1e-5
      C(:)       = 1e-12

      ! calculate rate coefficients for Monte-Carlo factors -1 and 1, then
      ! verify that the results are different:
      mcexp(:) = -1.
      CALL Update_RCONST()
      RCONST0(:) = RCONST(:)
      mcexp(:) = 1.
      CALL Update_RCONST()
      PRINT *, HLINE2
      PRINT *, 'Monte-Carlo Check (search for "WRONG"):'
      PRINT *, HLINE2
      DO i=1, NREACT
        WRITE(*,'(I4,A15)', ADVANCE='NO') i, TRIM(EQN_TAGS(i))
        IF (ABS(RCONST0(i)-RCONST(i)) < TINY_DP) THEN
          WRITE(*,'(1PE14.5,A)') RCONST(i), '  ***** WRONG *****'
        ELSE
          WRITE(*,'(A,2(1PE14.5))') '  OKAY       ', RCONST0(i), RCONST(i)
        ENDIF
      ENDDO
      PRINT *, HLINE2
      PRINT *, 'Finished montecarlo_check which is currently activated in'
      PRINT *, 'mecca/template_messy_mecca_kpp.f90'
      STOP 1

    END SUBROUTINE montecarlo_check

    !-------------------------------------------------------------------------

END SUBROUTINE kpp_integrate

! end modified insertion from kpp_integrate from kp4/src/create_kpp_module.C

! start modified insertion from kp4-generated fill routines

  SUBROUTINE fill_TEMP(status,array)
    INTEGER,INTENT(OUT)               :: status
    REAL (dp),INTENT(IN),DIMENSION(:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    TEMP = array(is)
  END SUBROUTINE fill_TEMP

  SUBROUTINE fill_cair(status,array)
    INTEGER,INTENT(OUT)               :: status
    REAL (dp),INTENT(IN),DIMENSION(:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    cair = array(is)
  END SUBROUTINE fill_cair

  SUBROUTINE fill_press(status,array)
    INTEGER,INTENT(OUT)               :: status
    REAL (dp),INTENT(IN),DIMENSION(:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    press = array(is)
  END SUBROUTINE fill_press

  SUBROUTINE fill_mcexp(status,array)
    INTEGER,INTENT(OUT)                 :: status
    REAL (dp),INTENT(IN),DIMENSION(:,:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    mcexp = array(is,:)
  END SUBROUTINE fill_mcexp

  SUBROUTINE fill_xaer(status,array)
    INTEGER,INTENT(OUT)                 :: status
    REAL (dp),INTENT(IN),DIMENSION(:,:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    xaer = array(is,:)
  END SUBROUTINE fill_xaer

  SUBROUTINE fill_cvfac(status,array)
    INTEGER,INTENT(OUT)                 :: status
    REAL (dp),INTENT(IN),DIMENSION(:,:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    cvfac = array(is,:)
  END SUBROUTINE fill_cvfac

  SUBROUTINE fill_lwc(status,array)
    INTEGER,INTENT(OUT)                 :: status
    REAL (dp),INTENT(IN),DIMENSION(:,:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    lwc = array(is,:)
  END SUBROUTINE fill_lwc

  SUBROUTINE fill_k_exf(status,array)
    INTEGER,INTENT(OUT)                   :: status
    REAL (dp),INTENT(IN),DIMENSION(:,:,:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    k_exf = array(is,:,:)
  END SUBROUTINE fill_k_exf

  SUBROUTINE fill_k_exb(status,array)
    INTEGER,INTENT(OUT)                   :: status
    REAL (dp),INTENT(IN),DIMENSION(:,:,:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    k_exb = array(is,:,:)
  END SUBROUTINE fill_k_exb

  SUBROUTINE fill_k_exf_N2O5(status,array)
    INTEGER,INTENT(OUT)                 :: status
    REAL (dp),INTENT(IN),DIMENSION(:,:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    k_exf_N2O5 = array(is,:)
  END SUBROUTINE fill_k_exf_N2O5

  SUBROUTINE fill_k_exf_ClNO3(status,array)
    INTEGER,INTENT(OUT)                 :: status
    REAL (dp),INTENT(IN),DIMENSION(:,:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    k_exf_ClNO3 = array(is,:)
  END SUBROUTINE fill_k_exf_ClNO3

  SUBROUTINE fill_k_exf_BrNO3(status,array)
    INTEGER,INTENT(OUT)                 :: status
    REAL (dp),INTENT(IN),DIMENSION(:,:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    k_exf_BrNO3 = array(is,:)
  END SUBROUTINE fill_k_exf_BrNO3

  SUBROUTINE fill_JX(status,array)
    INTEGER,INTENT(OUT)                 :: status
    REAL (dp),INTENT(IN),DIMENSION(:,:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    JX = array(is,:)
  END SUBROUTINE fill_JX

  SUBROUTINE fill_khet_Tr(status,array)
    INTEGER,INTENT(OUT)                 :: status
    REAL (dp),INTENT(IN),DIMENSION(:,:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    khet_Tr = array(is,:)
  END SUBROUTINE fill_khet_Tr

  SUBROUTINE fill_khet_St(status,array)
    INTEGER,INTENT(OUT)                 :: status
    REAL (dp),INTENT(IN),DIMENSION(:,:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    khet_St = array(is,:)
  END SUBROUTINE fill_khet_St

  ! mz_ab_20101119+
  SUBROUTINE fill_temp_ion(status,array)
    INTEGER,INTENT(OUT)                 :: status
    REAL (dp),INTENT(IN),DIMENSION(:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    temp_ion = array(is)
  END SUBROUTINE fill_temp_ion

  SUBROUTINE fill_temp_elec(status,array)
    INTEGER,INTENT(OUT)                 :: status
    REAL (dp),INTENT(IN),DIMENSION(:) :: array
    ! this subroutine only treats the first (is=1) element:
    INTEGER, PARAMETER :: is = 1
    status = 0
    temp_elec = array(is)
  END SUBROUTINE fill_temp_elec 
  ! mz_ab_20101119-

! end modified insertion from kp4-generated fill routines

END MODULE messy_mecca_kpp
