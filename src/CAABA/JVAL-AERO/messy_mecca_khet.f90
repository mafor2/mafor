!*****************************************************************************
!                Time-stamp: <2010-11-24 17:31:07 sander>
!*****************************************************************************

! AUTHORS
!   Rolf Sander,     MPICH, 2006

MODULE messy_mecca_khet

  USE messy_main_constants_mem, ONLY: &
    dp,      & ! kind parameter for real
    pi,      & ! pi
    R_gas      ! universal gas constant / (J/(mol*K))
  USE messy_mecca,              ONLY: modstr

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: modstr

  ! NAME OF SUBSUBMODEL
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodstr = modstr//'_khet'

  ! SUBROUTINES
  PUBLIC :: mecca_khet_read_nml_ctrl
  ! FUNCTIONS
  PUBLIC :: k_mt          ! rate coefficient for mass transfer

  ! CTRL_KHET NAMELIST
  LOGICAL, PUBLIC, SAVE :: &
    l_troposphere  = .FALSE., &
    l_stratosphere = .FALSE.

  REAL(dp), PUBLIC, PARAMETER :: &
    M_HNO3  = 63.012_dp, & ! molar mass / (g/mol)
    M_N2O5  = 108.00_dp    ! molar mass / (g/mol)

  REAL(dp), PUBLIC, PARAMETER :: & ! uptake coefficients
    gamma_N2O5 = 0.02_dp, &        ! Evans and Jacob, 2005
    gamma_HNO3 = 0.1_dp, &
    gamma_Hg   = 0.0_dp, &
    gamma_RGM  = 0.1_dp

CONTAINS

  !***************************************************************************

  SUBROUTINE mecca_khet_read_nml_ctrl(status, iou)

    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
      

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    !--- Local variables:

    NAMELIST /CTRL_KHET/ &
      l_troposphere, l_stratosphere

    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_khet_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE

    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    !--- 1) Read namelist:

    CALL read_nml_open(lex, substr, iou, 'CTRL_KHET', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_KHET, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_KHET', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    ! CHECK TRACER INITIALIZATION

    WRITE(*,*) '    Initialize module KHET'
    IF (l_troposphere) THEN
      WRITE(*,*) '  Heterogenous reaction rates for all tropospheric aerosol'
    ELSE
      WRITE(*,*) '  Heterogenous reaction rates not considered'//&
        ' for tropospheric aerosol, except N2O5'
    END IF
    IF (l_stratosphere) THEN
      WRITE(*,*) '  Heterogenous reaction rates for stratospheric aerosol'
    ELSE
      WRITE(*,*) '  Heterogenous reaction rates not considered'//&
        ' for stratospheric aerosol'
    END IF

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE mecca_khet_read_nml_ctrl

  !***************************************************************************

  ELEMENTAL REAL(dp) FUNCTION k_mt(M, T, A, gamma)

    ! calculate mass transfer coefficient k_mt

    IMPLICIT NONE
    INTRINSIC :: SQRT

    REAL(dp), INTENT (IN) :: M     ! molar mass [g/mol]
    REAL(dp), INTENT (IN) :: T     ! temperature [K]
    REAL(dp), INTENT (IN) :: A     ! aerosol surface [m2/m3]
    REAL(dp), INTENT (IN) :: gamma ! uptake coefficient [1]

    REAL(dp) :: v ! mean molecular speed (Maxwell-Boltzmann distribution) [m/s]

    v    = SQRT(8._dp*R_gas*T/(1E-3_dp*M*pi))    ! v = sqrt(8RT/(M*PI))   
    k_mt = 0.25_dp * v * gamma * A               ! [1/s]

  END FUNCTION k_mt

  !***************************************************************************

END MODULE messy_mecca_khet
