!*****************************************************************************
!                Time-stamp: <2019-03-13 14:55:59 sander>
!*****************************************************************************

! submodel MECCA
! Calculates chemistry
! written by:
!   Astrid Kerkweg, MPICH, Mainz, June 2003/Jan 2004
!   Rolf Sander,    MPICH, Mainz, 2003-2005

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

MODULE messy_mecca

  USE messy_main_constants_mem, ONLY: DP, HLINE2, STRLEN_SHORT

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mecca_read_nml_ctrl           ! read CTRL namelist and initialize
  PUBLIC :: steady_state_reached, steady_state_reached_old
  PUBLIC :: define_mcexp
  PUBLIC :: kpromanlev_to_nbl
  PUBLIC :: nbl_to_nblx
  INTERFACE nblx_to_nbl
    MODULE PROCEDURE nblx_to_nbl_int
    MODULE PROCEDURE nblx_to_nbl_real
  END INTERFACE
  PUBLIC :: nblx_to_nbl

  CHARACTER(LEN=*), PUBLIC, PARAMETER :: modstr = 'mecca' !> name of module
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: modver = '4.0'   !> module version
  LOGICAL, PUBLIC, SAVE :: l_aero     !> switch for aero chemistry

  ! GLOBAL CTRL-NAMELIST
  CHARACTER(LEN=STRLEN_SHORT), PUBLIC, SAVE :: mecca_aero = 'AUTO' !> for l_aero
  LOGICAL, PUBLIC, SAVE :: l_force_khet = .FALSE. !> switch for khet
  LOGICAL, PUBLIC, SAVE :: l_kpp_debug  = .FALSE. !> switch for kpp debugging
  LOGICAL, PUBLIC, SAVE :: l_tag        = .FALSE. !> switch for tagging
  INTEGER, PUBLIC, SAVE :: mcexp_seed   = 0       !> Monte-Carlo factor seed

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE mecca_read_nml_ctrl(status, iou)

    ! READ MECCA NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Astrid Kerkweg, MPICH, June 2003
    !         Rolf Sander, 2003, 2008

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_mecca_kpp,  ONLY: REQ_AEROSOL, REQ_MCFCT

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status !> error status
    INTEGER, INTENT(IN)  :: iou    !> logical I/O unit

    ! LOCAL
    LOGICAL :: lex   !> file exists?
    INTEGER :: fstat !> file status
    CHARACTER(LEN=*), PARAMETER :: substr = 'mecca_read_nml_ctrl'

    NAMELIST /CTRL/ mecca_aero, l_force_khet, l_kpp_debug, l_tag, &
      mcexp_seed

    ! INITIALIZE
    status = 1 !> error

    ! INPUT NAMELIST
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.NOT.lex) RETURN    !> <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  !> error while reading namelist

    SELECT CASE (TRIM(mecca_aero))
    CASE ('ON')
      WRITE(*,*) 'setting l_aero = T because mecca_aero = ON'
      l_aero = .TRUE.
    ! mz_rs_20150223+
    ! "OFF" is not allowed anymore
    ! CASE ('OFF')
    !   WRITE(*,*) 'setting l_aero = F because mecca_aero = OFF'
    ! l_aero = .FALSE.
    ! mz_rs_20150223-
    CASE ('AUTO')
      WRITE(*,*) 'setting l_aero = REQ_AEROSOL'
      l_aero = REQ_AEROSOL
    CASE DEFAULT
      WRITE(*,*) 'mecca_aero = ', TRIM(mecca_aero)
      WRITE(*,*) 'mecca_aero must be [ON/AUTO]'
      RETURN ! return with status=1
    END SELECT

    IF (REQ_AEROSOL.NEQV.l_aero) &
      WRITE(*,*) 'WARNING: REQ_AEROSOL and l_aero are different!'
    WRITE(*,*) 'mecca_aero   = ', TRIM(mecca_aero)
    WRITE(*,*) 'REQ_AEROSOL  = ', REQ_AEROSOL
    WRITE(*,*) 'l_aero       = ', l_aero
    WRITE(*,*) 'l_force_khet = ', l_force_khet
    WRITE(*,*) 'l_kpp_debug  = ', l_kpp_debug
    WRITE(*,*) 'l_tag        = ', l_tag
    IF (REQ_MCFCT) &
      WRITE(*,*) 'mcexp_seed   = ', mcexp_seed

    CALL read_nml_close(substr, iou, modstr)
    status = 0 !> no error

  END SUBROUTINE mecca_read_nml_ctrl

  ! --------------------------------------------------------------------------

  LOGICAL FUNCTION steady_state_reached_old(c,timesteplen)

    USE messy_mecca_kpp, ONLY: ind_OH, ind_HO2
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: c
    REAL(DP),               INTENT(IN) :: timesteplen
    REAL(DP), SAVE :: old_oh  = 0.
    REAL(DP), SAVE :: old_ho2 = 0.
    REAL(DP) :: change_oh, change_ho2
     
    steady_state_reached_old = .FALSE.

    IF (ind_OH >0) change_oh  = ABS((old_oh-c(ind_oh))/c(ind_oh))/timesteplen
    IF (ind_HO2>0) change_ho2 = ABS((old_ho2-c(ind_ho2))/c(ind_ho2))/timesteplen
    ! Steady state is defined here as less than a relative 1e-6 change
    ! per second. Note that this definition is probably only useful if
    ! the day/night cycle is switched off.
    IF ( (change_oh<1e-6).AND.(change_ho2<1e-6) ) THEN
      steady_state_reached_old = .TRUE.
    ENDIF
    IF (ind_OH >0) old_oh  = c(ind_oh)
    IF (ind_HO2>0) old_ho2 = c(ind_ho2)

  END FUNCTION steady_state_reached_old

  ! --------------------------------------------------------------------------

  !mz_hr_20160425-
  LOGICAL FUNCTION steady_state_reached(c,timesteplen)

    ! Steady state is defined as relative concentration change of less
    ! than 'ssthresh' per second for the concentrations of the species
    ! defined by 'ssind*', typically OH and HO2. For OH and HO2, this
    ! definition is probably only useful if the day/night cycle is
    ! switched off; see also caaba.f90, SUBROUTINE calc_sza, boolean
    ! l_steady_state_stop. Introduced new variable names were
    ! introducedmodified to allow easier change to other variables for
    ! steady state definition by just changing USE statement and the
    ! ssind* definitions. Note: Prior to the changes mentioned above,
    ! ind_* == 0 would have resulted in change_oh not to be initialized
    ! and not to be defined

    USE messy_mecca_kpp, ONLY: ind_OH, ind_HO2, SPC_NAMES

    IMPLICIT NONE

    INTRINSIC :: TRIM

    REAL(DP), DIMENSION(:), INTENT(IN) :: c
    REAL(DP),               INTENT(IN) :: timesteplen

    ! LOCAL
    INTEGER(DP) :: ssind1, ssind2
    REAL(DP) :: ssval1 = 0.
    REAL(DP) :: ssval2 = 0.
    REAL(DP), SAVE :: ssval1_old = 0.
    REAL(DP), SAVE :: ssval2_old = 0.
    REAL(DP) :: change_ssval1, change_ssval2
    REAL(DP), PARAMETER :: ssthresh = 1.e-5
    LOGICAL, SAVE :: l_init = .TRUE.

    steady_state_reached = .FALSE.

    ssind1 =  ind_OH
    ssind2 =  ind_HO2

    IF (ssind1 > 0) ssval1 = C(ssind1)
    IF (ssind2 > 0) ssval2 = C(ssind2)
    
    IF (ssval1 <= 0.0_dp .AND. ssval1_old <= 0.0_dp) THEN
      ! both values (close to) zero => relative change undefined, set to zero
      change_ssval1 = 0.0_dp
    ELSE
      ! use average of old and new in denominator to catch cases where one of
      ! the two is zero
      change_ssval1 = ABS((ssval1_old-ssval1)/((ssval1_old+ssval1)/2))/timesteplen
    ENDIF
    IF (ssval2 <= 0.0_dp .AND. ssval2_old <= 0.0_dp) THEN
      ! both values (close to zero), relative change meaningless, set to zero
      change_ssval2 = 0.0_dp
    ELSE
      ! use average of old and new in denominator to catch cases where one of
      ! the two is zero
      change_ssval2 = ABS((ssval2_old-ssval2)/((ssval2_old+ssval2)/2))/timesteplen
    ENDIF

    IF (l_init) THEN
      ! no difference to be determined; make sure program will not stop
      change_ssval1 = ssthresh*2._dp
      change_ssval2 = ssthresh*2._dp
      l_init = .FALSE.
    ELSE
      WRITE(*,'(ES10.2,A,ES10.2,A,A,A)') change_ssval1, " / ", ssthresh,&
           " (s-1 change in ", TRIM(SPC_NAMES(ssind1)), " / steady-state threshold)"
      WRITE(*,'(ES10.2,A,ES10.2,A,A,A)') change_ssval2, " / ", ssthresh,&
           " (s-1 change in ", TRIM(SPC_NAMES(ssind2)), " / steady-state threshold)"
    ENDIF

    IF ( (change_ssval1 <= ssthresh).AND.(change_ssval2 <= ssthresh) ) THEN
      steady_state_reached = .TRUE.
    ENDIF

    ssval1_old = ssval1
    ssval2_old = ssval2

  END FUNCTION steady_state_reached
  !mz_hr_20160425- 

  ! --------------------------------------------------------------------------

  SUBROUTINE define_mcexp(status, mcexp)

    USE messy_main_rnd,   ONLY: RND_MTW_GAUSS, rnd_init, rnd_number, rnd_finish
    IMPLICIT NONE
    INTEGER,                INTENT(OUT) :: status
    REAL(DP), DIMENSION(:), INTENT(OUT) :: mcexp
    INTEGER :: id_rnd

    ! assign a set of normally distributed random numbers to mcexp:
    CALL rnd_init(status, id_rnd, RND_MTW_GAUSS, mcexp_seed)
    IF (status/=0) RETURN
    CALL rnd_number(id_rnd, mcexp(:))
    CALL rnd_finish(id_rnd)

  END SUBROUTINE define_mcexp

  !---------------------------------------------------------------------------
  
  ! simliar to the f90 function PACK
  FUNCTION nbl_to_nblx(array_nbl, l_meccanum_nbl, nblx)

    ! I/O:
    REAL(DP), DIMENSION(:), INTENT(IN) :: array_nbl
    LOGICAL,  DIMENSION(:), INTENT(IN) :: l_meccanum_nbl
    INTEGER,                INTENT(IN) :: nblx
    REAL(DP), DIMENSION(nblx) :: nbl_to_nblx
    ! local:
    INTEGER :: jb, jbx

    jbx = 0
    DO jb = 1, SIZE(array_nbl)
      IF (l_meccanum_nbl(jb)) THEN
        jbx = jbx + 1
        nbl_to_nblx(jbx) = array_nbl(jb)
      ENDIF
    ENDDO

  END FUNCTION nbl_to_nblx

  !---------------------------------------------------------------------------
  
  ! simliar to the f90 function UNPACK
  SUBROUTINE nblx_to_nbl_real(array_nblx, array_nbl, l_meccanum_nbl)

    ! I/O:
    REAL(DP), DIMENSION(:), INTENT(IN)    :: array_nblx
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: array_nbl
    LOGICAL,  DIMENSION(:), INTENT(IN)    :: l_meccanum_nbl
    ! local:
    INTEGER :: jb, jbx

    jbx = 0
    DO jb = 1, SIZE(array_nbl)
      IF (l_meccanum_nbl(jb)) THEN
        jbx = jbx + 1
        array_nbl(jb) = array_nblx(jbx)
      ENDIF
    ENDDO

  END SUBROUTINE nblx_to_nbl_real

  !---------------------------------------------------------------------------
  
  ! simliar to the f90 function UNPACK
  SUBROUTINE nblx_to_nbl_int(array_nblx, array_nbl, l_meccanum_nbl)

    ! I/O:
    INTEGER, DIMENSION(:), INTENT(IN)    :: array_nblx
    INTEGER, DIMENSION(:), INTENT(INOUT) :: array_nbl
    LOGICAL, DIMENSION(:), INTENT(IN)    :: l_meccanum_nbl
    ! local:
    INTEGER :: jb, jbx

    jbx = 0
    DO jb = 1, SIZE(array_nbl)
      IF (l_meccanum_nbl(jb)) THEN
        jbx = jbx + 1
        array_nbl(jb) = array_nblx(jbx)
      ENDIF
    ENDDO

  END SUBROUTINE nblx_to_nbl_int

  !---------------------------------------------------------------------------
  
#define ALLATONCE
  ! simliar to the f90 function RESHAPE
  FUNCTION kpromanlev_to_nbl(array_2d)

    ! I/O:
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: array_2d
    ! kpromanlev_to_nbl has the same overall size as array_2d but is 1d:
    REAL(DP), DIMENSION(SIZE(array_2d)) :: kpromanlev_to_nbl
    ! local:
    INTEGER :: jb, jk, jp, kproma, nlev

    kproma = SIZE(array_2d,1)
    nlev = SIZE(array_2d,2)
#ifdef ALLATONCE
    DO jk = 1, nlev 
      kpromanlev_to_nbl(kproma*(jk-1)+1:kproma*jk) = array_2d(:,jk)
    ENDDO
#else
    jb = 0
    DO jk = 1, nlev 
      DO jp = 1, kproma
        jb = jb + 1
        kpromanlev_to_nbl(jb) = array_2d(jp,jk)
      ENDDO
    ENDDO
#endif

  END FUNCTION kpromanlev_to_nbl

!*****************************************************************************

END MODULE messy_mecca

!*****************************************************************************
