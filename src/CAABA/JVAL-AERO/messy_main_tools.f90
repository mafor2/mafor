! ************************************************************************
MODULE messy_main_tools
! ************************************************************************

  ! MESSy-SMCL tools
  !
  ! Authors: Patrick Joeckel,     MPICH, Mainz,            2003-...
  !          Hella Riede,         MPICH, Mainz,            2009-...
  !          Rolf Sander,         MPICH, Mainz,            2006-...
  !          Sabine Brinkop,      DLR,   Oberpfaffenhofen  2012-...
  !          Christopher Kaiser,  DLR,   Oberpfaffenhofen  2013-...

  USE messy_main_constants_mem, ONLY: SP, DP, STRLEN_LONG, STRLEN_VLONG

  IMPLICIT NONE
  PRIVATE

  ! mz_rs_20081222+
  ! 1D array also for integer:
  TYPE PTR_1D_ARRAY_INT
     INTEGER, DIMENSION(:), POINTER :: PTR => NULL()
  END TYPE PTR_1D_ARRAY_INT
  ! mz_rs_20081222-

  ! op_pj_20120912+
  TYPE PTR_0D_ARRAY
     REAL(DP), POINTER :: PTR => NULL()
  END TYPE PTR_0D_ARRAY
  ! op_pj_20120912-

  TYPE PTR_1D_ARRAY
     REAL(DP), DIMENSION(:), POINTER :: PTR => NULL()
  END TYPE PTR_1D_ARRAY

  TYPE PTR_2D_ARRAY
     REAL(DP), DIMENSION(:,:), POINTER :: PTR => NULL()
  END TYPE PTR_2D_ARRAY

  TYPE PTR_3D_ARRAY
     REAL(DP), DIMENSION(:,:,:), POINTER :: PTR => NULL()
  END TYPE PTR_3D_ARRAY

  TYPE PTR_4D_ARRAY
     REAL(DP), DIMENSION(:,:,:,:), POINTER :: PTR => NULL()
  END TYPE PTR_4D_ARRAY

  TYPE PTR_5D_ARRAY
     REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: PTR => NULL()
  END TYPE PTR_5D_ARRAY

  PUBLIC :: PTR_0D_ARRAY, PTR_1D_ARRAY, PTR_2D_ARRAY, PTR_3D_ARRAY &
          , PTR_4D_ARRAY, PTR_5D_ARRAY, PTR_1D_ARRAY_INT

  ! fb_mk_20120116+
  TYPE t_reset_par
     LOGICAL  :: l = .FALSE.
     REAL(dp) :: v = 0.0_dp
  END TYPE t_reset_par
  PUBLIC :: t_reset_par
  !fb_mk_20120116-

  ! mz_ht_20040414+
  ! FOR LOOKUP TABLE INITIALIZATION (CONVECTION at al.)
  ! lower and upper bound in K * 1000
  INTEGER,  PARAMETER, PUBLIC :: jptlucu1 =  50000  ! lookup table lower bound
  INTEGER,  PARAMETER, PUBLIC :: jptlucu2 = 400000  ! lookup table upper bound
  ! table - e_s*Rd/Rv
  REAL(dp), TARGET, SAVE,  PUBLIC :: tlucua(jptlucu1:jptlucu2)
  ! table - for derivative calculation
  REAL(dp), TARGET, SAVE,  PUBLIC :: tlucub(jptlucu1:jptlucu2)
  ! table - l/cp
  REAL(dp), TARGET, SAVE,  PUBLIC :: tlucuc(jptlucu1:jptlucu2)
  ! table
  REAL(dp), TARGET, SAVE,  PUBLIC :: tlucuaw(jptlucu1:jptlucu2)
  ! mz_ht_20040414-

  ! mz_rs_20060110+
  INTERFACE str
     MODULE PROCEDURE str_logical
     MODULE PROCEDURE str_integer
     MODULE PROCEDURE str_real_sp
     MODULE PROCEDURE str_real_dp
  END INTERFACE
  PUBLIC :: str
  ! mz_rs_20060110-
  ! ub_ak_20170406+
  PUBLIC :: lcase
  ! ub_ak_20170406-

  ! mz_hr_20160508+
  PUBLIC :: psat_mk  ! H2O psat [Pa], Murphy-Koop 2005, replaces psat, psatf
  PUBLIC :: spec2relhum      ! convert specific to relative humidity
  PUBLIC :: rel2spechum      ! convert relative to specific humidity
  PUBLIC :: spec2relhum_q    ! same as spec2relhum with RH = q / q_s
  PUBLIC :: cair_q           ! calculate cair from temp, press, spechum
  PUBLIC :: cair_c           ! calculate cair from temp, press, c(H2O)
  PUBLIC :: spechum2mr       ! calculate H2O mixing ratio from spechum
  PUBLIC :: mr2spechum       ! calculate spechum from H2O mixing ratio
  ! obsolete, use spechum2mr + rel2spechum PUBLIC :: rh2mr ! rel. hum. -> mix. ratio (calls relhum*2mr)
  ! obsolete, use spechum2mr + rel2spechum PUBLIC :: relhum2mr ! relative humidity -> mol(H2O)/mol(dryair)
  ! obsolete, use spechum2mr + rel2spechum PUBLIC :: relhumwmo2mr ! WMO relative humidity -> mol(H2O)/mol(dryair)
  ! obsolete, use cair_q + rel2spechum     PUBLIC :: cair_trad
  ! obsolete, use cair_q + rel2spechum     PUBLIC :: cair_wmo
  ! obsolete, use spec2relhum              PUBLIC :: spec2relhumwmo
  ! obsolete, use rel2spechum              PUBLIC :: rel2spechumwmo
  ! mz_hr_20160508- 
  ! op_bk_20170111+
  PUBLIC :: get_name_domain
  ! op_bk_20170111-

  INTERFACE iso2ind
     MODULE PROCEDURE iso2ind_1d
     MODULE PROCEDURE iso2ind_2d
  END INTERFACE

  INTERFACE ind2val
     MODULE PROCEDURE ind2val_1d
     MODULE PROCEDURE ind2val_2d
  END INTERFACE

! mz_ab_20100610+
  INTERFACE remap_bounds
     MODULE PROCEDURE remap_bounds1
     MODULE PROCEDURE remap_bounds2
     MODULE PROCEDURE remap_bounds3
     MODULE PROCEDURE remap_bounds4
  END INTERFACE
! mz_ab_20100610-

  ! mz_bk_20110707+
  INTERFACE str2num
     MODULE PROCEDURE str2num_real_dp
     MODULE PROCEDURE str2num_real_sp
     MODULE PROCEDURE str2num_integer
  END INTERFACE
  ! mz_bk_20110707-

  ! SUBROUTINES
  PUBLIC :: read_nml_open       ! Utilities ...
  PUBLIC :: read_nml_check      ! ... to simplify ...
  PUBLIC :: read_nml_close      ! ... namelist input

  PUBLIC :: iso2ind             ! find index
  PUBLIC :: ind2val             ! find value at index level
  PUBLIC :: int2str             ! convert integer to string
  PUBLIC :: strcrack            ! cracking strings into parts
  PUBLIC :: nn_index            ! look for nearest neighbour(s) in ordered list
                                ! (minimum distance)

  ! op_sb_20120724+ 
  PUBLIC :: nl_index            ! find index above/below in monotonous list
  ! op_sb_20120724-

  ! mz_ab_20100624+
  PUBLIC :: fliparray           ! reverses the order of 1-D array
  ! mz_ab_20100624-

  ! mz_ap_20070913+
  PUBLIC :: ns_index            ! look for sourrounding neighbour(s) in list
  ! mz_ap_20070913-

  ! mz_ht_20040414+
  PUBLIC :: init_convect_tables ! lookup table for convection et al. 
  ! mz_ht_20040414-

  PUBLIC :: match_wild          ! compare strings with wildcards
  PUBLIC :: str2chob            ! convert string to channel/object list
  PUBLIC :: bilin_weight        ! weights for bilinear interpolation

  PUBLIC :: find_next_free_unit

  ! mz_ab_20100610+
  PUBLIC :: remap_bounds  ! for obtaining pointers(-1:) etc, used for bounds
  ! mz_ab_20100610-

  ! mz_rs_20101118+
! op_ck_20131218+
!!$  PUBLIC :: density
  PUBLIC :: mass_density
! op_ck_20131218-
  PUBLIC :: layerthickness
  ! mz_rs_20101118-

  ! mz_ab_20100831+
  PUBLIC :: full2half           ! converts variable on full level pressures
                                ! to half level pressures
  PUBLIC :: spline1D            ! spline interpolation
  PUBLIC :: splint1D
  ! mz_ab_20100831-
  ! mz_ab_20111125+
  PUBLIC :: CalcTimeAngles
  ! mz_ab_20111125+

  ! mz_bk_20110707+
  PUBLIC :: str2num             ! convert string to numerical value
  ! mz_bk_20110707-
  ! mz_hr_20100704+ 
  PUBLIC :: ucase               ! turn a string into all uppercase
  ! mz_hr_20100704- 

  ! cy_ma_20120625+
  PUBLIC :: ERRFUNC             ! error function: evaluates erf(x), erfc(x), 
                                ! and exp(x*x)*erfc(x) for a real argument  x
  ! cy_ma_20120625-

  ! um_ak_20131111+
  PUBLIC :: calc_hybrid_coeff
  ! um_ak_20131111-

  ! op_mm_20140527+
  PUBLIC :: interpol_stag2mid
  PUBLIC :: interpol_mid2stag
  ! op_mm_20140527-

  ! um_ak_20140814+
  PUBLIC :: tautsp2D
  ! um_ak_20140814-
  ! op_pj_20161215+
  PUBLIC :: quick_sort
  ! op_pj_20161215-

  ! mz_ht_20140925+
  PUBLIC :: IDXSORT_DP          ! sorts the indices of a double precision array
                                !       beginning with the lowest value
  ! mz_ht_20140925-

  ! op_pj_20180523+
#ifdef LF
  PUBLIC :: isnan               ! workaround for Lahey/Fujitsu Compiler 8.1b
#endif
  ! op_pj_20180523-

CONTAINS

! -----------------------------------------------------------------------
  SUBROUTINE read_nml_open(lex, substr, iou, nmlstr, modstr, l_print)

    USE messy_main_blather, ONLY: start_message
    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    LOGICAL, INTENT(OUT)                   :: lex       ! file exists ?
    CHARACTER(LEN=*), INTENT(IN)           :: substr    ! calling routine
    INTEGER,          INTENT(IN)           :: iou       ! unit
    CHARACTER(LEN=*), INTENT(IN)           :: nmlstr    ! namelist
    CHARACTER(LEN=*), INTENT(IN)           :: modstr    ! module name
    LOGICAL,          INTENT(IN), OPTIONAL :: l_print

    ! LOCAL:
    LOGICAL :: zl_print

    IF (PRESENT(l_print)) THEN
       zl_print = l_print
    ELSE
       zl_print = .TRUE.  ! DEFAULT
    END IF

    CALL start_message(TRIM(modstr), 'INITIALISATION', substr, zl_print)

    ! CHECK IF FILE EXISTS
    INQUIRE(file=TRIM(modstr)//'.nml', exist=lex)
    IF (.NOT.lex) THEN
       WRITE(*,*) '*** WARNING: FILE '''//TRIM(modstr)//'.nml'&
            &//'''  NOT FOUND !'
       RETURN
    END IF

    ! OPEN FILE
    OPEN(iou,file=TRIM(modstr)//'.nml')
    IF (zl_print) &
         WRITE(*,*) 'Reading namelist '''//TRIM(nmlstr)//''''//&
         &' from '''//TRIM(modstr)//'.nml',''' (unit ',iou,') ...'

  END SUBROUTINE read_nml_open
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE read_nml_check(fstat, substr, iou, nmlstr, modstr, l_print)

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(IN)                    :: fstat     ! file status
    CHARACTER(LEN=*), INTENT(IN)           :: substr    ! calling routine
    INTEGER,          INTENT(IN)           :: iou       ! unit
    CHARACTER(LEN=*), INTENT(IN)           :: nmlstr    ! namelist
    CHARACTER(LEN=*), INTENT(IN)           :: modstr    ! module name
    LOGICAL,          INTENT(IN), OPTIONAL :: l_print

    ! LOCAL
    LOGICAL :: zl_print

    IF (PRESENT(l_print)) THEN
       zl_print = l_print
    ELSE
       zl_print = .TRUE.  ! DEFAULT
    END IF

    IF (fstat /= 0) THEN
       WRITE(*,*) '*** ERROR: READ ERROR in NAMELIST '''//TRIM(nmlstr)//''''&
            &//' in FILE '''//TRIM(modstr)//'.nml'//''' !'
       CLOSE(iou)
    ELSE
       IF (zl_print) WRITE(*,*) ' ... OK !'
    END IF

  END SUBROUTINE read_nml_check
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE read_nml_close(substr, iou, modstr, l_print)

    USE messy_main_blather, ONLY: end_message
    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)           :: substr    ! calling routine
    INTEGER,          INTENT(IN)           :: iou       ! unit
    CHARACTER(LEN=*), INTENT(IN)           :: modstr    ! module name
    LOGICAL,          INTENT(IN), OPTIONAL :: l_print

    ! LOCAL
    LOGICAL :: zl_print

    IF (PRESENT(l_print)) THEN
       zl_print = l_print
    ELSE
       zl_print = .TRUE.  ! DEFAULT
    END IF

    CLOSE(iou)

    CALL end_message(TRIM(modstr), 'INITIALISATION', substr, zl_print)

  END SUBROUTINE read_nml_close
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE iso2ind_1d(field, iso, k, f, lrev)
    
    ! ISOSURFACE TO INDEX
    ! OUTPUT IS THE LEVEL INDEX FOR A GIVEN ISOSURFACE
    ! NOTE:
    !   THIS ROUTINE IS WORKING PROPERLY ONLY IF THE 'FIELD'
    !   IS MONOTONIC
    !   (e.g., POTENTIAL TEMPERATURE, PRESSURE, POTENTIAL VORTICITY, ...)
    ! METHOD:
    !  lrev = .false. (default) -> SEARCH FROM LOW INDEX TO HIGH INDEX
    !  lrev = .true.            -> SEARCH FROM HIGH INDEX TO LOW INDEX
    !  ADJUST INDEX
    !  OPTIONAL: FRACTION OF BOX 'BELOW' ISO-SURFACE
    
    IMPLICIT NONE
    INTRINSIC :: NINT, ABS, MAX, MIN, PRESENT, REAL, SIZE, TINY

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN)            :: field  ! input field
    REAL(DP),               INTENT(IN)            :: iso    ! isosurface value
    INTEGER ,               INTENT(OUT)           :: k      ! isosurface level
    REAL(DP),               INTENT(OUT), OPTIONAL :: f      ! fraction in layer
    LOGICAL,                INTENT(IN),  OPTIONAL :: lrev   ! reverse order ?

    ! LOCAL
    INTEGER :: nk, jk, dk
    ! mz_rs_20080224: dk changed to integer
    REAL(DP) :: zf
    LOGICAL  :: llrev
    INTEGER  :: jstart, jstop, jstep

    k = 0
    nk = SIZE(field)    

    IF (PRESENT(lrev)) THEN
       llrev = lrev
    ELSE
       llrev = .FALSE. ! default
    END IF

    IF (llrev) THEN
       jstart = nk
       jstop  = 2
       jstep = -1
    ELSE
       jstart = 2
       jstop  = nk
       jstep = 1
    END IF

    DO jk = jstart, jstop, jstep
       IF ( ( (iso >= field(jk-1)) .AND. (iso <= field(jk)) ) .OR. &
            ( (iso <= field(jk-1)) .AND. (iso >= field(jk)) ) ) THEN
          k=jk
          EXIT
       END IF
    END DO

    IF ( k == 0 ) THEN   ! NOT FOUND
       IF (llrev) THEN
!!$       k = 2
          k = 1
          IF (PRESENT(f)) f = 1.0_dp
       ELSE
          k = nk
          IF (PRESENT(f)) f = 0.0_dp
       END IF
       RETURN
    END IF

    ! ADJUST INDEX
    ! CALCULATE FRACTION OF BOX 'BELOW' ISO-SURFACE
    !
    ! METHOD: LINEAR INTERPOLATION
    !
    ! THE FOLLOWING CONDITION MUST ALWAYS BE .TRUE.,
    ! SINCE THE FIRST LEVEL WITH 
    !    FIELD(k-1) <= ISO <= FLIELD(k)
    ! OR
    !    FIELD(k-1) >= ISO >= FLIELD(k)
    ! IS SEARCHED
    !
    IF ( ABS( field(k) - field(k-1) ) > TINY(0.0) ) THEN
       zf = ABS( (iso-field(k-1)) / (field(k)-field(k-1)) )    ! e [0,1)
    ELSE
       zf = 0.5_dp  ! SHOULD BE NEVER REACHED !!!
    END IF
    
    zf = MIN(1._dp,zf)
    zf = MAX(0._dp,zf)
    
    ! dk = INT(zf+0.5)
    dk = NINT(zf)
    ! zf e [0,0.5] -> dk = 0 -> ONE LEVEL ABOVE (-1)
    ! zf e (0.5,1) -> dk = 1 -> KEEP LEVEL

    k = k - 1 + dk
    
    ! CALCULATE FRACTION OF BOX 'BELOW' ISO-SURFACE
    ! ONE LEVEL ABOVE (dk = 0) -> zf e [0, 0.5]
    !                          -> f  e [0.5, 0]
    !       EXAMPLE: zf  = 0   -> ISO AT BOX MID         -> FRACT. = 0.5
    !                zf  = 0.5 -> ISO AT LOWER INTERFACE -> FRACT. = 0.0
    ! KEEP LEVEL      (dk = 1) -> zf e (0.5, 1)
    !                          -> f  e (1, 0.5)
    !       EXAMPLE: zf  = 0.5 -> ISO AT UPPER INTERFACE -> FRACT. = 1.0
    !                zf  = 1   -> ISO AT BOX MID         -> FRACT. = 0.5
    
    !                            FOR dk=1            FOR dk=0
    IF (PRESENT(f)) f = (1.5_dp-zf)*REAL(dk,DP) + (0.5_dp-zf)*REAL(1-dk,DP)
     
  END SUBROUTINE iso2ind_1d
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE iso2ind_2d(kproma, field, iso, k, f, lrev)
    
    ! as iso2ind_1d, however for vectors ...
    
    IMPLICIT NONE
    INTRINSIC :: NINT, ABS, MAX, MIN, PRESENT, REAL, SIZE, TINY

    ! I/O
    INTEGER,                  INTENT(IN)            :: kproma
    REAL(DP), DIMENSION(:,:), INTENT(IN)            :: field ! input field
    REAL(DP), DIMENSION(:),   INTENT(IN)            :: iso   ! isosurface value
    INTEGER,  DIMENSION(:),   INTENT(OUT)           :: k     ! isosurface level
    REAL(DP), DIMENSION(:),   INTENT(OUT), OPTIONAL :: f     ! fraction in layer
    LOGICAL,                  INTENT(IN),  OPTIONAL :: lrev  ! reverse order ?

    ! LOCAL
    INTEGER  :: nk, jk, jl, dk
    ! mz_rs_20080224: dk changed to integer
    REAL(DP) :: zf
    LOGICAL  :: llrev
    INTEGER  :: jstart, jstop, jstep
    INTEGER, DIMENSION(kproma) :: ilfound

    k = 0
    ilfound(:) = 0
    nk = SIZE(field,2)    

    IF (PRESENT(lrev)) THEN
       llrev = lrev
    ELSE
       llrev = .FALSE. ! default
    END IF

    IF (llrev) THEN
       jstart = nk
       jstop  = 2
       jstep = -1
    ELSE
       jstart = 2
       jstop  = nk
       jstep = 1
    END IF

    DO jk = jstart, jstop, jstep
       DO jl = 1, kproma
          IF ( ilfound(jl) == 1) CYCLE
          IF ( ( (iso(jl) >= field(jl,jk-1)) .AND. &
               (iso(jl) <= field(jl,jk)) ) .OR. &
               ( (iso(jl) <= field(jl,jk-1)) .AND. &
               (iso(jl) >= field(jl,jk)) ) ) THEN
             k(jl) = jk
             ilfound(jl) = 1 
          END IF
       END DO
    END DO

    DO jl = 1, kproma

       IF ( k(jl) == 0 ) THEN   ! NOT FOUND
          IF (llrev) THEN
!!$          k(jl) = 2
             k(jl) = 1
             IF (PRESENT(f)) f(jl) = 1.0_dp
          ELSE
             k(jl) = nk
             IF (PRESENT(f)) f(jl) = 0.0_dp
          END IF

       ELSE
          
          ! ADJUST INDEX
          ! CALCULATE FRACTION OF BOX 'BELOW' ISO-SURFACE
          !
          ! METHOD: LINEAR INTERPOLATION
          !
          ! THE FOLLOWING CONDITION MUST ALWAYS BE .TRUE.,
          ! SINCE THE FIRST LEVEL WITH 
          !    FIELD(k-1) <= ISO <= FLIELD(k)
          ! OR
          !    FIELD(k-1) >= ISO >= FLIELD(k)
          ! IS SEARCHED
          !
          IF ( ABS( field(jl,k(jl)) - field(jl,k(jl)-1) ) > TINY(0.0_dp) ) THEN
             zf = ABS( (iso(jl)-field(jl,k(jl)-1)) / &
                  (field(jl,k(jl))-field(jl,k(jl)-1)) )    ! e [0,1)
          ELSE
             zf = 0.5_dp  ! SHOULD BE NEVER REACHED !!!
          END IF
    
          zf = MIN(1._dp,zf)
          zf = MAX(0._dp,zf)
    
          ! dk = INT(zf+0.5)
          dk = NINT(zf)
          ! zf e [0,0.5] -> dk = 0 -> ONE LEVEL ABOVE (-1)
          ! zf e (0.5,1) -> dk = 1 -> KEEP LEVEL
    
          k(jl) = k(jl) - 1 + dk
    
          ! CALCULATE FRACTION OF BOX 'BELOW' ISO-SURFACE
          ! ONE LEVEL ABOVE (dk = 0) -> zf e [0, 0.5]
          !                          -> f  e [0.5, 0]
          !       EXAMPLE: zf  = 0   -> ISO AT BOX MID         -> FRACT. = 0.5
          !                zf  = 0.5 -> ISO AT LOWER INTERFACE -> FRACT. = 0.0
          ! KEEP LEVEL      (dk = 1) -> zf e (0.5, 1)
          !                          -> f  e (1, 0.5)
          !       EXAMPLE: zf  = 0.5 -> ISO AT UPPER INTERFACE -> FRACT. = 1.0
          !                zf  = 1   -> ISO AT BOX MID         -> FRACT. = 0.5

          !                            FOR dk=1            FOR dk=0
          IF (PRESENT(f)) f(jl) = (1.5_dp-zf)*REAL(dk,dp) + &
               (0.5_dp-zf)*REAL(1-dk,dp)
          
       END IF

    END DO
     
  END SUBROUTINE iso2ind_2d
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE ind2val_1d(val, field, k, f)

    ! CONVERT INDEX (AND FRACTION 'BELOW') IN A MONOTONOUS FIELD
    ! TO THE VALUE
    ! METHOD:
    !   - PICK OUT INDEX
    !   - LINEAR INTERPOLATION BETWEEN NEIGHBOURS, IF f IS PRESENT
    
    USE messy_main_constants_mem, ONLY: FLAGGED_BAD

    IMPLICIT NONE
    INTRINSIC :: PRESENT, SIZE, LBOUND, UBOUND

    ! I/O
    REAL(DP),               INTENT(OUT)          :: val   ! value
    REAL(DP), DIMENSION(:), INTENT(IN)           :: field ! field
    INTEGER,                INTENT(IN)           :: k     ! level
    REAL(DP),               INTENT(IN), OPTIONAL :: f     ! fraction
    
    ! LOCAL
    INTEGER  :: nk
    REAL(DP) :: ri, gf
    
    nk = SIZE(field)
    
    IF (PRESENT(f)) THEN
       ri  = 0.5_dp - f           ! e (-0.5,0.5) -> (top, bottom) of box
       IF (ri >= 0.0_dp) THEN
          IF (k == nk) THEN
             gf  = (field(k)-field(k-1))
          ELSE
             gf  = (field(k+1)-field(k))
          END IF
       ELSE
          IF (k == 1) THEN
             gf  = (field(k+1)-field(k))
          ELSE
             gf  = (field(k)-field(k-1))
          END IF
       END IF
       !
       val = field(k) + ri * gf
    ELSE
       ! op_pj_20180204+
       IF ((k < LBOUND(field,1)) .OR. (k > UBOUND(field,1))) THEN
          val = FLAGGED_BAD
       ELSE                              ! op_pj_20180204-
          val = field(k)
       ENDIF                             ! op_pj_20180204
    END IF
    
  END SUBROUTINE ind2val_1d
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE ind2val_2d(kproma, val, field, k, f)

    ! as ind2val_1d, however for vectors ...
    
    IMPLICIT NONE
    INTRINSIC :: PRESENT, SIZE

    ! I/O
    INTEGER,                  INTENT(IN)           :: kproma
    REAL(DP), DIMENSION(:),   INTENT(OUT)          :: val   ! value
    REAL(DP), DIMENSION(:,:), INTENT(IN)           :: field ! field
    INTEGER,  DIMENSION(:),   INTENT(IN)           :: k     ! level
    REAL(DP), DIMENSION(:),   INTENT(IN), OPTIONAL :: f     ! fraction
    
    ! LOCAL
    INTEGER  :: jl, nk
    REAL(DP) :: ri, gf
    
    nk = SIZE(field,2)
    
    IF (PRESENT(f)) THEN
       DO jl = 1, kproma
          ri  = 0.5_dp - f(jl)       ! e (-0.5,0.5) -> (top, bottom) of box
          IF (ri >= 0.0_dp) THEN
             IF (k(jl) == nk) THEN
                gf  = (field(jl,k(jl))-field(jl,k(jl)-1))
             ELSE
                gf  = (field(jl,k(jl)+1)-field(jl,k(jl)))
             END IF
          ELSE
             IF (k(jl) == 1) THEN
                gf  = (field(jl,k(jl)+1)-field(jl,k(jl)))
             ELSE
                gf  = (field(jl,k(jl))-field(jl,k(jl)-1))
             END IF
          END IF
          !
          val(jl) = field(jl,k(jl)) + ri * gf
       END DO
    ELSE
       DO jl = 1, kproma
          val(jl) = field(jl,k(jl))
       END DO
    END IF
    
  END SUBROUTINE ind2val_2d
! -----------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE int2str(str, ii, cpad, cerr)
    
    IMPLICIT NONE
    INTRINSIC :: MOD, LEN, PRESENT
    
    ! I/O
    CHARACTER(LEN=*), INTENT(OUT)          :: str   ! STRING
    INTEGER,          INTENT(IN)           :: ii    ! INTEGER
    CHARACTER,        INTENT(IN), OPTIONAL :: cpad  ! CHAR FOR PADDING
    CHARACTER,        INTENT(IN), OPTIONAL :: cerr  ! CHAR FOR ERROR
    
    ! LOCAL
    INTEGER              :: n, zi, zn, k
    INTEGER              :: rest
    CHARACTER            :: zcpad
    
    IF (PRESENT(cpad)) THEN
       zcpad = cpad
    ELSE
       zcpad = '0'      ! DEFAULT PADDING
    END IF
    
    n  = LEN(str)
    zi = ii
    zn = n
    
    DO
       rest = MOD(zi, 10)
       zi   = zi/10
       WRITE(str(zn:zn),'(i1)') rest
       zn = zn - 1
       IF (zi == 0) EXIT
       IF (zn == 0) EXIT
    END DO
    
    IF (PRESENT(cerr)) THEN
       IF ( (zn == 0) .AND. (zi /= 0) ) THEN
          DO k = 1, n
             WRITE(str(k:k),'(a1)') cerr
          END DO
       END IF
    END IF
    
    DO k = zn, 1, -1
       WRITE(str(k:k),'(a1)') zcpad
    END DO
    
  END SUBROUTINE int2str
! ---------------------------------------------------------------------

! mz_rs_20060110+
! ---------------------------------------------------------------------
  ! various types of the function str:

  CHARACTER(LEN=5) FUNCTION str_logical(zlogical, fmt)
    ! create string from logical
    IMPLICIT NONE
    INTRINSIC :: PRESENT
    LOGICAL, INTENT(in) :: zlogical
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: fmt
    IF (PRESENT(fmt)) THEN
      WRITE(str_logical,fmt) zlogical
    ELSE
      IF (zlogical) THEN
        str_logical = 'TRUE '
      ELSE
        str_logical = 'FALSE'
      ENDIF
    ENDIF
  END FUNCTION str_logical

  CHARACTER(LEN=STRLEN_LONG) FUNCTION str_integer(zinteger, fmt)
    ! create string from integer
    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, PRESENT
    INTEGER, INTENT(in) :: zinteger
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: fmt
    IF (PRESENT(fmt)) THEN
      WRITE(str_integer,fmt) zinteger
    ELSE
      WRITE(str_integer,*) zinteger
      str_integer = ADJUSTL(str_integer) ! remove leading spaces
    ENDIF
  END FUNCTION str_integer

  CHARACTER(LEN=STRLEN_LONG) FUNCTION str_real_sp(zreal_sp, fmt)
    ! create string from real_sp
    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, PRESENT
    REAL(sp), INTENT(in) :: zreal_sp
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: fmt
    IF (PRESENT(fmt)) THEN
      WRITE(str_real_sp,fmt) zreal_sp
    ELSE
      WRITE(str_real_sp,*) zreal_sp
      str_real_sp = ADJUSTL(str_real_sp) ! remove leading spaces
    ENDIF
  END FUNCTION str_real_sp

  CHARACTER(LEN=STRLEN_LONG) FUNCTION str_real_dp(zreal_dp, fmt)
    ! create string from real_dp
    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, PRESENT
    REAL(dp), INTENT(in) :: zreal_dp
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: fmt
    IF (PRESENT(fmt)) THEN
      WRITE(str_real_dp,fmt) zreal_dp
    ELSE
      WRITE(str_real_dp,*) zreal_dp
      str_real_dp = ADJUSTL(str_real_dp) ! remove leading spaces
    ENDIF
  END FUNCTION str_real_dp
! ---------------------------------------------------------------------
! mz_rs_20060110-

! ---------------------------------------------------------------------
  SUBROUTINE strcrack(str, ch, el, n)

    ! strcrack = string crack

    ! Split the string <str> into small pieces which are separated by
    ! the character <ch>. Delete trailing spaces from the resulting <n>
    ! pieces, then put them into the array <el>.

    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, ASSOCIATED, INDEX, LEN_TRIM, TRIM

    ! I/O
    CHARACTER(LEN=*),               INTENT(IN)  :: str ! string 2b cracked
    CHARACTER,                      INTENT(IN)  :: ch  ! separating char
    CHARACTER(LEN=*), DIMENSION(:), POINTER     :: el  ! field of substrings
    INTEGER,                        INTENT(OUT) :: n   ! # of substrings
    
    ! LOCAL
    INTEGER :: idx1, idx2, i
    
    ! INIT
    IF (ASSOCIATED(el)) DEALLOCATE(el)
    NULLIFY(el)
    n = 0
    
    ! EMPTY STRING
    IF ( (TRIM(str) == '') .OR. (TRIM(str) == ch) ) RETURN
    
    idx1 = 0
    idx2 = 0
    DO 
       idx1 = idx2 + 1
       IF (idx1 > LEN_TRIM(str(:))) EXIT
       IF (INDEX(TRIM(str(idx1:)), ch) == 0) THEN
          idx2 = LEN_TRIM(str(:)) + 1
       ELSE
          idx2 = idx2 + INDEX(TRIM(str(idx1:)), ch)
       END IF
       IF (idx1 == idx2) CYCLE
       
       n = n + 1
       
    END DO
    
    ! ALLOCATE SPACE
    ALLOCATE(el(n))
    DO i=1, n
       el(i) = ''
    END DO
    
    n = 0
    idx1 = 0
    idx2 = 0
    DO     
       idx1 = idx2 + 1
       IF (idx1 > LEN_TRIM(str(:))) EXIT
       IF (INDEX(TRIM(str(idx1:)), ch) == 0) THEN
          idx2 = LEN_TRIM(str(:)) + 1
       ELSE
          idx2 = idx2 + INDEX(TRIM(str(idx1:)), ch)
       END IF
       IF (idx1 == idx2) CYCLE
       
       n = n + 1
       
!       el(n) = str(idx1:idx2-1)
       el(n) = ADJUSTL(str(idx1:idx2-1))
       
    END DO

  END SUBROUTINE strcrack
! ---------------------------------------------------------------------  

! mz_ap_20070913+
! ---------------------------------------------------------------------
  SUBROUTINE ns_index(list, value, idx, idx2, l_periodic)

    IMPLICIT NONE
    INTRINSIC :: PRESENT, SIZE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN)            :: list
    REAL(DP),               INTENT(IN)            :: value
    INTEGER,                INTENT(OUT)           :: idx
    INTEGER,                INTENT(OUT), OPTIONAL :: idx2
    LOGICAL,                INTENT(IN),  OPTIONAL :: l_periodic 
    ! list(i) < value < list(i+1)
    ! l_periodic : periodic boundaries

    ! LOCAL
    INTEGER     :: n, i

    ! INIT
    n = SIZE(list)
    idx = 0

    ! we suppouse that both lists are ordered
    IF (value<list(1)) THEN
      DO i=1, n
         IF (value <= list(i)) THEN
            idx = i
         END IF
      END DO
    ELSE
      DO i=1, n
         IF (value >= list(i)) THEN
            idx = i
         END IF
      END DO
    END IF
    IF (PRESENT(idx2)) THEN
       IF (idx == n) THEN
          IF (PRESENT(l_periodic)) THEN
            idx2 = 1
          ELSE
            idx2 = n
          ENDIF
       ELSE
          idx2 = idx+1
       ENDIF
    END IF
    IF (idx == 0 ) THEN
      IF (PRESENT(l_periodic)) THEN 
        idx=n
        IF (PRESENT(idx2)) idx2=1
      ELSE 
        idx=1
        IF (PRESENT(idx2)) idx2=2
      ENDIF
    ENDIF

  END SUBROUTINE ns_index
! mz_ap_20070913-
! --------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE nn_index(list, value, idx, idx2, fac, status)

    IMPLICIT NONE
    INTRINSIC :: ABS, PRESENT, SIZE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN)            :: list
    REAL(DP),               INTENT(IN)            :: value
    INTEGER,                INTENT(OUT)           :: idx
    INTEGER,                INTENT(OUT), OPTIONAL :: idx2
    REAL(DP),               INTENT(IN),  OPTIONAL :: fac
    INTEGER,                INTENT(OUT), OPTIONAL :: status

    ! LOCAL
    REAL(DP)    :: zfac
    INTEGER     :: n, i
    REAL(DP)    :: dmin

    ! INIT
    n = SIZE(list)
    idx = 0
    IF (PRESENT(idx2)) idx2 = 0
    IF (PRESENT(fac)) THEN
       zfac = fac
    ELSE
       zfac = 2.0_dp     ! default
    END IF
    IF (PRESENT(status)) status = 0

    ! LIST MUST CONTAIN AT LEAST TWO ELEMENTS
    IF (n < 2) THEN
       idx = 1
       IF (PRESENT(idx2))   idx2 = 1
       IF (PRESENT(status)) status = 1
       RETURN
    END IF

    ! ALLOW zfac TIMES THE FIRST / LAST INTERVAL FOR EXTRAPOLATION ...
    IF (list(1) <= list(n)) THEN
       ! increasing ...
       IF (value < (list(1) - (list(2)-list(1))*zfac) ) THEN
          idx = 1
          IF (PRESENT(idx2))   idx2 = 1
          IF (PRESENT(status)) status = 2
          RETURN
       END IF
       IF (value > (list(n) + (list(n)-list(n-1))*zfac) ) THEN
          idx = n
          IF (PRESENT(idx2))   idx2 = n
          IF (PRESENT(status)) status = 3
          RETURN
       END IF
    ELSE
       ! decreasing ...
       IF (value > (list(1) - (list(2)-list(1))*zfac) ) THEN
          idx = 1
          IF (PRESENT(idx2))   idx2 = 1
          IF (PRESENT(status)) status = 4
          RETURN
       END IF
       IF (value < (list(n) + (list(n)-list(n-1))*zfac) ) THEN
          idx = n
          IF (PRESENT(idx2))   idx2 = n
          IF (PRESENT(status)) status = 5
          RETURN
       ENDIF
    END IF

    ! SEARCH MINIMUM
    dmin = ABS(list(1) - list(n))
    DO i=1, n
       IF (ABS(list(i)-value) <= dmin) THEN
          dmin = ABS(list(i)-value)
          idx = i
       END IF
    END DO

    ! NOTE: FOR IDX2 IT IS ASSUMED THAT THE LIST IS ORDERED
    IF (PRESENT(idx2)) THEN
       IF (idx == 1) THEN
          idx2 = 2
          RETURN
       END IF
       IF (idx == n) THEN
          idx2 = n-1
          RETURN
       END IF
       IF ( ABS(list(idx+1)-value) <= ABS(list(idx-1)-value) ) THEN
          idx2 = idx+1
       ELSE
          idx2 = idx-1
       END IF
    END IF

  END SUBROUTINE nn_index
! --------------------------------------------------------------------
! op_sb_20120724+
! ---------------------------------------------------------------------
  SUBROUTINE nl_index(list, value, idx, status)

    IMPLICIT NONE
    INTRINSIC :: PRESENT, SIZE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN)            :: list
    REAL(DP),               INTENT(IN)            :: value
    INTEGER,                INTENT(OUT)           :: idx
    INTEGER,                INTENT(OUT), OPTIONAL :: status

    ! LOCAL
    INTEGER     :: n, i

    ! INIT
    n = SIZE(list)
    idx = 0
  
    IF (PRESENT(status)) status = 0

    ! LIST MUST CONTAIN AT LEAST TWO ELEMENTS
    IF (n < 2) THEN
       idx = 1
       IF (PRESENT(status)) status = 1
       RETURN
    END IF

    IF (list(1) <= list(n)) THEN
       ! increasing ...
       IF (value < list(1) ) THEN
          idx = 1
          IF (PRESENT(status)) status = 2
          RETURN
       END IF
       IF (value > list(n) ) THEN
          idx = n
          IF (PRESENT(status)) status = 3
          RETURN
       END IF
       DO i = 1, n
         IF (value >= list(i)) then
          idx = i 
         ELSE
          RETURN
         ENDIF
       END DO    
    ELSE
       ! decreasing ...
       IF (value > list(1) ) THEN
          idx = 1
          IF (PRESENT(status)) status = 4
          RETURN
       END IF
       IF (value < list(n) ) THEN
          idx = n
          IF (PRESENT(status)) status = 5
          RETURN
       ENDIF
       DO i = 1, n
         IF (value <= list(i)) then
          idx = i 
       ELSE
          RETURN
         ENDIF
       END DO    
    END IF

  END SUBROUTINE nl_index
! --------------------------------------------------------------------
! op_sb_20120724-

! mz_ab_20100624+
! --------------------------------------------------------------------
  FUNCTION fliparray(list) RESULT(buffer)
    
    ! flips (reverses) array:
    ! (1,2,3,4) --> (4,3,2,1)

    IMPLICIT NONE
    INTRINSIC :: SIZE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN) :: list
    ! LOCAL
    REAL(DP), DIMENSION(:), ALLOCATABLE :: buffer
    ! mz_rs_20140310+
    ! FORCHECK SAYS:
    !    1042  REAL(DP), DIMENSION(:), ALLOCATABLE :: buffer
    !   BUFFER
    ! (file: messy_main_tools.f90, line:    1042)
    ! **[216 E] illegally specified automatic, static or allocatable
    !           (The object is a function result)
    ! mz_rs_20140310-

    INTEGER :: i, n

    n = SIZE(list,1)

    ALLOCATE(buffer(n))

    DO i=1,n
       buffer(i) = list(n-i+1)
    END DO

  END FUNCTION fliparray
! --------------------------------------------------------------------
! mz_ab_20100624-

! mz_ht_20042510+
! ---------------------------------------------------------------------  
  SUBROUTINE init_convect_tables

    ! Lookup tables for convective adjustment code
    !
    ! D. Salmond, CRAY (UK), August 1991, original code

    USE messy_main_constants_mem, ONLY: rd, rv, tmelt, cpd => cp_air &
                                      , alv, als ! op_pj_20100209

    IMPLICIT NONE
    INTRINSIC :: EXP, LOG, REAL
                       
    REAL(dp), PARAMETER :: zavl1 = -6096.9385_dp
    REAL(dp), PARAMETER :: zavl2 =    21.2409642_dp
    REAL(dp), PARAMETER :: zavl3 =    -2.711193_dp
    REAL(dp), PARAMETER :: zavl4 =     1.673952_dp
    REAL(dp), PARAMETER :: zavl5 =     2.433502_dp 

    REAL(dp), PARAMETER :: zavi1 = -6024.5282_dp
    REAL(dp), PARAMETER :: zavi2 =    29.32707_dp
    REAL(dp), PARAMETER :: zavi3 =     1.0613868_dp
    REAL(dp), PARAMETER :: zavi4 =    -1.3198825_dp
    REAL(dp), PARAMETER :: zavi5 =    -0.49382577_dp   
    
! op_pj_20100209+
!!$    REAL(dp), PARAMETER :: alv   = 2.5008e6_dp ! latent heat for
!!$    !                                          ! vaporisation in J/kg
!!$    REAL(dp), PARAMETER :: als   = 2.8345e6_dp ! latent heat for
!!$    !                                          ! sublimation in J/kg
! op_pj_20100209-
    ! Constants used for computation of saturation mixing ratio
    ! over liquid water (*c_les*) or ice(*c_ies*)
    REAL(dp),PARAMETER :: c3les = 17.269_dp           ! 
    REAL(dp),PARAMETER :: c3ies = 21.875_dp           ! 
    REAL(dp),PARAMETER :: c4les = 35.86_dp            ! 
    REAL(dp),PARAMETER :: c4ies = 7.66_dp             ! 
    REAL(dp),PARAMETER :: c5les = c3les*(tmelt-c4les) ! 
    REAL(dp),PARAMETER :: c5ies = c3ies*(tmelt-c4ies) ! 
   
    REAL(dp) :: z5alvcp, z5alscp, zalvdcp, zalsdcp
    REAL(dp) :: ztt, zldcp
!!$    REAL(dp) :: zcvm3, zcvm4, zcvm5
    REAL(dp) :: zcvm4, zcvm5
    REAL(dp) :: zavm1, zavm2, zavm3, zavm4, zavm5

    INTEGER :: it

    z5alvcp = c5les*alv/cpd
    z5alscp = c5ies*als/cpd

    zalvdcp = alv/cpd
    zalsdcp = als/cpd

    DO it = jptlucu1, jptlucu2
      ztt = 0.001_dp*REAL(it,dp)
      IF ((ztt-tmelt) > 0.0_dp) THEN
!!$        zcvm3 = c3les
        zcvm4 = c4les
        zcvm5 = z5alvcp
        zldcp = zalvdcp
        zavm1 = zavl1
        zavm2 = zavl2
        zavm3 = zavl3
        zavm4 = zavl4
        zavm5 = zavl5
      ELSE
!!$        zcvm3 = c3ies
        zcvm4 = c4ies
        zcvm5 = z5alscp
        zldcp = zalsdcp
        zavm1 = zavi1
        zavm2 = zavi2
        zavm3 = zavi3
        zavm4 = zavi4
        zavm5 = zavi5
      END IF
      tlucuc(it)  = zldcp
      tlucua(it)  = EXP((zavm1/ztt+zavm2+zavm3*0.01_dp* &
           ztt+zavm4*ztt*ztt*1.e-5_dp+zavm5*LOG(ztt)))*rd/rv
      tlucub(it)  = zcvm5*(1.0_dp/(ztt-zcvm4))**2
      tlucuaw(it) = EXP((zavl1/ztt+zavl2+zavl3*0.01_dp* &
           ztt+zavl4*ztt*ztt*1.e-5_dp+zavl5*LOG(ztt)))*rd/rv
    END DO
    
  END SUBROUTINE init_convect_tables
! ---------------------------------------------------------------------  
! mz_ht_20042510-

! ---------------------------------------------------------------------  
LOGICAL FUNCTION match_wild(pattern, string)

   ! compare given string for match to pattern which may
   ! contain wildcard characters:
   ! "?" matching any one character, and
   ! "*" matching any zero or more characters.
   ! Both strings may have trailing spaces which are ignored.
   ! Authors: Clive Page, userid: cgp  domain: le.ac.uk, 2003 (original code)
   !          Rolf Sander, 2005 (bug fixes and pattern preprocessing)
   ! Minor bug fixed by Clive Page, 2005 Nov 29.

   ! This program is free software; you can redistribute it and/or modify
   ! it under the terms of the GNU General Public License as published by
   ! the Free Software Foundation; either version 2 of the License, or
   ! (at your option) any later version.
   !
   ! This program is distributed in the hope that it will be useful,
   ! but WITHOUT ANY WARRANTY; without even the implied warranty of
   ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   ! GNU General Public License for more details.
   !
   ! You should have received a copy of the GNU General Public License
   ! along with this program; if not, write to the Free Software
   ! Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 
   ! 02110-1301  USA

   IMPLICIT NONE
   INTRINSIC :: INDEX, LEN, LEN_TRIM, REPEAT

   CHARACTER(LEN=*), INTENT(IN) :: pattern ! pattern may contain * and ?
   CHARACTER(LEN=*), INTENT(IN) :: string  ! string to be compared
   INTEGER :: lenp, lenp2, lens, n, p2, p, s
   INTEGER :: n_question, n_asterisk

   CHARACTER(LEN=LEN(pattern)) :: pattern2

   lens = LEN_TRIM(string)
   lenp = LEN_TRIM(pattern)

   ! If the pattern is empty, always return true
   IF (lenp == 0) THEN
     match_wild = .TRUE.
     RETURN
   ENDIF

   ! The pattern must be preprocessed. All consecutive occurences of
   ! one or more question marks ('?') and asterisks ('*') are sorted and
   ! compressed. The result is stored in pattern2.

   pattern2(:)=''
   p  = 1 ! current position in pattern
   p2 = 1 ! current position in pattern2
   DO
     IF ((pattern(p:p) == '?').OR.(pattern(p:p) == '*')) THEN
       ! a special character was found in the pattern
       n_question = 0
       n_asterisk = 0
       DO WHILE (p <= lenp)
         ! count the consecutive question marks and asterisks
         IF ((pattern(p:p) /= '?').AND.(pattern(p:p) /= '*')) EXIT
         IF (pattern(p:p) == '?') n_question = n_question + 1
         IF (pattern(p:p) == '*') n_asterisk = n_asterisk + 1
         p = p + 1
       ENDDO
       IF (n_question>0) THEN ! first, all the question marks
         pattern2(p2:p2+n_question-1) = REPEAT('?',n_question)
         p2 = p2 + n_question
       ENDIF
       IF (n_asterisk>0) THEN ! next, the asterisk (only one!)
         pattern2(p2:p2) = '*'
         p2 = p2 + 1
       ENDIF
     ELSE
       ! just a normal character
       pattern2(p2:p2) = pattern(p:p)
       p2 = p2 + 1
       p = p + 1
     ENDIF
     IF (p > lenp) EXIT
   ENDDO
   lenp2 = LEN_TRIM(pattern2)

   ! The modified wildcard in pattern2 is compared to the string:

   p2 = 1
   s = 1
   match_wild = .FALSE.
   DO
     IF (pattern2(p2:p2) == '?') THEN
       ! accept any char in string
       p2 = p2 + 1
       s = s + 1
     ELSEIF (pattern2(p2:p2) == "*") THEN
       p2 = p2 + 1
       IF (p2 > lenp2) THEN
         ! anything goes in rest of string
         match_wild = .TRUE.
         EXIT ! .TRUE.
       ELSE
         ! search string for char at p2
         n = INDEX(string(s:), pattern2(p2:p2))
         IF (n == 0) EXIT  ! .FALSE.
         s = n + s - 1
       ENDIF
     ELSEIF (pattern2(p2:p2) == string(s:s)) THEN
       ! single char match
       p2 = p2 + 1
       s = s + 1
     ELSE
       ! non-match
       EXIT ! .FALSE.
     ENDIF
     IF (p2 > lenp2 .AND. s > lens) THEN
       ! end of both pattern2 and string
       match_wild = .TRUE.
       EXIT ! .TRUE.
     ENDIF

     IF (s > lens .AND. p2 == lenp) THEN
       IF(pattern2(p2:p2) == "*") THEN
         ! "*" at end of pattern2 represents an empty string
         match_wild = .TRUE.
         EXIT ! .TRUE.
       END IF
     ENDIF
     IF (p2 > lenp2 .OR. s > lens) THEN
       ! end of either pattern2 or string
       EXIT ! .FALSE.
     ENDIF
   ENDDO

END FUNCTION match_wild
! ---------------------------------------------------------------------  

! ----------------------------------------------------------------------
  SUBROUTINE str2chob(status, str, n, c, o)

! op_pj_20100803+
!!$    USE messy_main_constants_mem, ONLY: STRLEN_ULONG
! op_pj_20100803-

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, TRIM, LEN

    ! I/O
    INTEGER,           INTENT(OUT)               :: status
    CHARACTER(LEN=*),  INTENT(IN)                :: str
    INTEGER,           INTENT(OUT)               :: n
    CHARACTER(LEN=*),  DIMENSION(:), POINTER     :: c
    CHARACTER(LEN=*),  DIMENSION(:), POINTER     :: o

    ! LOCAL
    INTEGER                                            :: nc, no, two, j, i
! op_pj_20100803+
!!$    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: tmp1 => NULL()
!!$    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: tmp2 => NULL()
!!$    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: tmp3 => NULL()
    CHARACTER(LEN=LEN(str)), DIMENSION(:), POINTER :: tmp1 !=> NULL()
    CHARACTER(LEN=LEN(str)), DIMENSION(:), POINTER :: tmp2 !=> NULL()
    CHARACTER(LEN=LEN(str)), DIMENSION(:), POINTER :: tmp3 !=> NULL()
! op_pj_20100803-

    ! INIT
    IF (ASSOCIATED(c)) DEALLOCATE(c)
    IF (ASSOCIATED(o)) DEALLOCATE(o)
    NULLIFY(c)
    NULLIFY(o)
    n = 0
    ! op_pj_20100811+
    NULLIFY(tmp1)
    NULLIFY(tmp2)
    NULLIFY(tmp3)
    ! op_pj_20100811-

    status = 1 ! ERROR

    ! COUNT OBJECTS
    CALL strcrack(TRIM(str), ';', tmp1, nc)           ! CHANNEL BLOCKS
    DO i=1, nc
       CALL strcrack(TRIM(tmp1(i)), ':', tmp2, two)   ! ONE CHANNEL PER BLOCK
       IF (two /= 2) RETURN      ! ERROR: 0 or more than 1 ':' in string
       CALL strcrack(TRIM(tmp2(2)), ',', tmp3, no)    ! OBJECTS PER CHANNEL
       n = n + no
    END DO

    ! ALLOCATE SPACE
    ALLOCATE(c(n))
    ALLOCATE(o(n))
    ! INIT
    DO i=1,n
       c(i) = ''
       o(i) = ''
    END DO

    ! PARSE STRING
    n = 0
    CALL strcrack(TRIM(str), ';', tmp1, nc)           ! CHANNEL BLOCKS
    DO i=1, nc
       CALL strcrack(TRIM(tmp1(i)), ':', tmp2, two)   ! ONE CHANNEL PER BLOCK
       CALL strcrack(TRIM(tmp2(2)), ',', tmp3, no)    ! OBJECTS PER CHANNEL
       DO j=1, no
          n = n + 1
          c(n) = TRIM(tmp2(1))
          o(n) = TRIM(tmp3(j))
       END DO
    END DO   

    IF (ASSOCIATED(tmp1)) THEN
      DEALLOCATE(tmp1) ; NULLIFY(tmp1)
    ENDIF
    IF (ASSOCIATED(tmp2)) THEN
      DEALLOCATE(tmp2) ; NULLIFY(tmp2)
    ENDIF
    IF (ASSOCIATED(tmp3)) THEN
      DEALLOCATE(tmp3) ; NULLIFY(tmp3)
    ENDIF

    status  = 0  ! NO ERROR

  END SUBROUTINE str2chob
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE bilin_weight(vn, v, w)
    
    IMPLICIT NONE

    ! I/O
    REAL(dp), DIMENSION(4, 2), INTENT(IN)  :: vn ! neighbours
    REAL(dp), DIMENSION(2),    INTENT(IN)  :: v  ! position
    REAL(dp), DIMENSION(4),    INTENT(OUT) :: w  ! weight

    ! LOCAL
    REAL(DP) :: t, u

    t = (v(1) - vn(1,1)) / (vn(2,1) - vn(1,1))
    u = (v(2) - vn(1,2)) / (vn(3,2) - vn(1,2))

    w(1) = (1._dp-t)*(1._dp-u)
    w(2) = t*(1._dp-u)
    w(3) = t*u
    w(4) = (1._dp-t)*u

  END SUBROUTINE bilin_weight
! ----------------------------------------------------------------------

! mz_hr_20100704+ 
! ----------------------------------------------------------------------
  REAL(dp) FUNCTION psat_mk(p_temp, p_l_psat_liquid)

    ! psat_mk replaces psat and psatf

    !> psat: saturation water vapor pressure [Pa] over liquid or ice as
    !> recommended by Murphy and Koop, Q. J. R. Meteorol. Soc. (2005),
    !> 131, pp. 1539-1565
    !> converging to zero, always > 0, defined at relatively low
    !>   temperatures as occurring in UTLS
    !> temperature in Kelvin [K]
    !> switch from ice to liquid above 273.15 K (= tmelt) in analogy
    !>   to convect_init_tables
    !> valid temperature range liquid formula: 123 K < T < 332 K
    !> valid temperature range ice formula:    110 K < T < 273.16 K

    !> always use l_psat_liquid=T for WMO relative humidity and for
    !> supercooling effects and supersaturation

    USE messy_main_constants_mem, ONLY: tmelt

    IMPLICIT NONE

    INTRINSIC LOG, EXP, TANH, PRESENT

    ! I/O
    REAL(dp), INTENT(IN)            :: p_temp     ! temperature [K]
    LOGICAL,  INTENT(IN), OPTIONAL  :: p_l_psat_liquid ! always over liquid surface?

    ! LOCAL
    REAL(dp), PARAMETER      :: lowlim_ice = 110._dp ! lower temp limit [K]
    REAL(dp), PARAMETER      :: lowlim_liq = 123._dp ! lower temp limit [K]
    REAL(dp), PARAMETER      :: uplim      = 332._dp ! upper temp limit [K]
    REAL(dp) :: lowlim
    LOGICAL  :: l_psat_liquid


    IF (PRESENT(p_l_psat_liquid)) THEN
      l_psat_liquid = p_l_psat_liquid
    ELSE
      l_psat_liquid = .FALSE.
    ENDIF

    IF (l_psat_liquid) THEN
      lowlim = lowlim_liq
    ELSE
      lowlim = lowlim_ice
    ENDIF

    ! define psat according to temperature and l_psat_liquid
    IF ((p_temp > tmelt) .OR. l_psat_liquid) THEN
      ! psat over liquid
      psat_mk = EXP( 54.842763_dp - 6763.22_dp/p_temp &
        - 4.21_dp*LOG(p_temp) + 0.000367_dp*p_temp    &
        + TANH(0.0415_dp*(p_temp-218.8_dp))           &
        * (53.878_dp - 1331.22_dp/p_temp              &
        - 9.44523_dp*LOG(p_temp) + 0.014025_dp*p_temp) )
    ELSE
      ! psat over ice
      psat_mk = EXP( 9.550426_dp - 5723.265_dp/p_temp &
        + 3.53068_dp*LOG(p_temp) - 0.00728332_dp*p_temp )
    ENDIF

    ! issue warning outside valid temperature range
    IF (p_temp > uplim) THEN
      WRITE(*,'(A,F6.1,A,F6.1,A,F6.1,A)') '  WARNING psat_mk:&
        & temperature (', p_temp, ' K) above valid temperature range [',&
        lowlim, ',', uplim, ']'
    ENDIF

    IF (p_temp < lowlim) THEN
      IF ((p_temp > tmelt) .OR. l_psat_liquid) THEN
        WRITE(*,'(A,F6.1,A,F6.1,A,F6.1,A)') '  WARNING psat_mk:&
          & temperature (', p_temp, ' K) below valid temperature range over&
          & liquid [', lowlim, ',', uplim, ']'
      ELSE
        WRITE(*,'(A,F6.1,A,F6.1,A,F6.1,A)') '  WARNING psat_mk:&
          & temperature (', p_temp, ' K) below valid temperature range over&
          & ice [', lowlim, ',', uplim, ']'
      ENDIF
    ENDIF

  END FUNCTION psat_mk
! ----------------------------------------------------------------------------
! mz_hr_20100704-

! mz_hr_20080221+
! ----------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE ucase(string)

  !> convert letters to uppercase
  !> from book Stephen Chapman "F90/95"

    IMPLICIT NONE

    INTRINSIC LGE, LLE, ACHAR, IACHAR, LEN

    ! I/O
    CHARACTER(LEN=*),        INTENT(INOUT)  :: string

    ! LOCAL
    INTEGER :: i 
    INTEGER :: length

    ! get len of str
    length = LEN(string)

    ! shift lower case to upper case
    DO i=1, length
      IF (LGE(string(i:i),'a') .AND. LLE(string(i:i),'z')) THEN
        string(i:i) = ACHAR(IACHAR(string(i:i)) - 32)
      ENDIF
    ENDDO
    
  END SUBROUTINE ucase
! ----------------------------------------------------------------------------
! mz_hr_20080221-

! mz_hr_20100607+
! ----------------------------------------------------------------------------
  ! merged spec2relhum and spec2relhumwmo here
  REAL(dp) FUNCTION spec2relhum(status, spechum, temp, press, &
                                l_hum_emac, l_psat_liquid, l_relhum_wmo)

  !> calculate relative humidity from specific humidity
  !> specific humidity q = (mass water vapor)/(mass humid air)
  !> l_relhum_wmo = FALSE
  !>   relhum = (partial pressure H2O)/(saturation partial pressure H2O)
  !>          = p(H2O)/psat(H2O)
  !> l_relhum_wmo = TRUE
  !>   calculate relative humidity as defined by WMO from specific humidity
  !>   NOTE: WMO relative humidity always defined over liquid
  !>   relative humidity   = omega_v/omega_vs
  !>     omega_v: mass mixing ratio of water vapor 
  !>              = (mass water vapor)/(mass dry air)
  !>     omega_vs: saturation mass mix. ratio of H2O (g) over liquid surface
  !>   MM_eps = (molar mass H2O)/(molar mass dry air) = M_H2O/M_air
  !>   [Jacobson, Fundamentals of Atmospheric Modeling, Cambridge University
  !>    Press, 1999, p. 23, 36]
  !> status = 1: WARNING
  !> status = 2: ERROR

    USE messy_main_constants_mem, ONLY: MM_eps, FLAGGED_BAD

    IMPLICIT NONE

    INTRINSIC :: ABS, INT, PRESENT

    !I/O
    INTEGER, INTENT(OUT) :: status     ! 1 = warning, 2 = error
    REAL(dp), INTENT(IN) :: spechum    ! kg/kg
    REAL(dp), INTENT(IN) :: temp       ! K
    REAL(dp), INTENT(IN) :: press      ! Pa
    !LOGICAL, INTENT(IN), OPTIONAL :: l_hum_emac
    LOGICAL, INTENT(IN), OPTIONAL :: l_hum_emac, l_psat_liquid, l_relhum_wmo

    !LOCAL
    REAL(dp) :: p_v, psat, omega_v, omega_vs
    LOGICAL  :: z_l_hum_emac, z_l_psat_liquid, z_l_relhum_wmo

    status = 0

    ! better readability by defining local logicals
    IF (PRESENT(l_hum_emac)) THEN
      z_l_hum_emac = l_hum_emac
    ELSE
      z_l_hum_emac = .FALSE.
    ENDIF

    IF (PRESENT(l_psat_liquid)) THEN
      z_l_psat_liquid = l_psat_liquid
    ELSE
      z_l_psat_liquid = .FALSE.
    ENDIF

    IF (PRESENT(l_relhum_wmo)) THEN
      z_l_relhum_wmo = l_relhum_wmo
    ELSE
      z_l_relhum_wmo = .FALSE.
    ENDIF

    IF (z_l_relhum_wmo) THEN
      z_l_psat_liquid = .TRUE. ! WMO: always above liquid
    ENDIF

    ! WARNING or ERROR status if spechum value out of bounds
    IF (spechum >= 1._dp) THEN
      ! WRITE(*,'(A, F10.3)') 'ERROR spec2relhum: spechum >= 1 : ', spechum
      WRITE(*,'(A,ES10.3,A)') 'ERROR spec2relhum: spechum >= 1 (',&
        spechum, '), relhum set to undefined'
      spec2relhum = FLAGGED_BAD
      status = 2
      RETURN
    ELSEIF (spechum < 0._dp) THEN
      WRITE(*,'(A,ES10.3,A)') 'ERROR spec2relhum: spechum < 0 (',&
        spechum, '), relhum set to undefined'
      spec2relhum = FLAGGED_BAD
      status = 2
      RETURN
    ENDIF
    
    ! determine saturation partial pressure of water vapor
    IF (z_l_hum_emac) THEN
      IF (z_l_psat_liquid) THEN
        psat = tlucuaw(INT(temp*1000._dp))/MM_eps
      ELSE
        !psat = tlucua(INT(z_temp*1000._dp))*Rv/Rd
        psat = tlucua(INT(temp*1000._dp))/MM_eps
      ENDIF
    ELSE
      psat = psat_mk(temp,z_l_psat_liquid)
    ENDIF

    ! water vapor mass mixing ratio in dry air
    omega_v = spechum/(1._dp-spechum)

    IF (z_l_relhum_wmo) THEN
      !> water vapor saturation mass mixing ratio
      omega_vs = MM_eps * psat / (press - psat)

      !> calculate WMO relative humidity
      spec2relhum = omega_v / omega_vs
    ELSE
      !> calculate water vapor partial pressure, taking into account that
      !>   pressure in the model is for humid air, i.e., this would be WRONG
      !>   since pressure NOT for dry air: p_v = spechum/MM_eps * press
      p_v = press * omega_v / (MM_eps + omega_v)

      !> calculate traditional relative humidity
      spec2relhum = p_v / psat
    ENDIF

    ! Warning for relative humidity > 1 (supersaturation)
    ! Error if negative relative humidity
    IF (spec2relhum > 1._dp) THEN
      WRITE(*,'(A,ES10.3,A)') '   WARNING spec2relhum: relhum > 1 (', &
        spec2relhum, ')'
      status = 1
    ELSEIF (spec2relhum < 0._dp) THEN
      WRITE(*,'(A,ES10.3,A)') '   ERROR spec2relhum: relhum < 0 (', &
        spec2relhum, '), relhum set to undefined'
      spec2relhum = FLAGGED_BAD
      status = 2
      RETURN
    ENDIF

  END FUNCTION spec2relhum
! --------------------------------------------------------------------- 
! mz_hr_20100607-

! op_ck_20131218+
  SUBROUTINE spec2relhum_q(status, relhum, spechum, temp, press, liq_only)

    ! Calculates RH from specific humidity q; RH defined as q / q_s, where q_s
    ! is specific humidity at saturation when keeping total pressure constant.

    USE messy_main_constants_mem, ONLY: vtmpc1

    IMPLICIT NONE
    INTRINSIC :: PRESENT, INT, MIN, MAX, EPSILON

    ! I/O
    INTEGER,  INTENT(OUT)          :: status
    REAL(dp), INTENT(OUT)          :: relhum
    REAL(dp), INTENT(IN)           :: spechum, temp, press
    LOGICAL,  INTENT(IN), OPTIONAL :: liq_only

    ! LOCAL
    LOGICAL  :: lonly
    INTEGER  :: it
    REAL(dp) :: qsat, psat

    ! Initialization
    status = 1 ! unspecified ERROR
    relhum = -1._dp
    it = 0
    qsat = -1._dp
    psat = -1._dp
    IF (PRESENT(liq_only)) THEN
       lonly = liq_only
    ELSE
       lonly = .FALSE.
    END IF
       
    ! Look up saturation vapor pressure (psat = tlucua? * M_air / M_H2O)
    ! and calculate saturation specific humidity (qsat)
    it = INT(temp * 1000._dp)
    IF (it < jptlucu1 .OR. it > jptlucu2) THEN
       status = 2 ! ERROR: Lookup overflow
       RETURN
    ELSE
       ! Obtain (and limit) saturation water vapor pressure from lookup table
       ! (adapted from messy_main_data_bi::main_data_global_start in MESSy2.50)
       IF (lonly) THEN
          psat = MIN(tlucuaw(it), 0.5_dp*press)
       ELSE
          psat = MIN(tlucua(it), 0.5_dp*press)
       END IF
    END IF

    ! Calculate (and limit) saturation specific humidity (adapted from
    ! messy_main_data_bi::main_data_global_start in MESSy2.50)
    qsat = MAX(2._dp*EPSILON(1._dp), psat/(press-vtmpc1*psat))

    status = 0 ! No ERROR

    relhum = spechum / qsat

  END SUBROUTINE spec2relhum_q
! op_ck_20131218-

! mz_hr_20100608+
! ---------------------------------------------------------------------
  REAL(dp) FUNCTION rel2spechum(status, relhum, temp, press, &
                                l_hum_emac, l_psat_liquid, l_relhum_wmo)

    !> calculates specific humidity from relative humidity
    !> for WMO definition, see, e.g., Jacobson, Fundamentals of Atmospheric
    !>   Modeling, Cambridge University Press, 1999
    !> see also inverse function spec2relhum

    USE messy_main_constants_mem, ONLY: MM_eps, FLAGGED_BAD

    IMPLICIT NONE

    INTRINSIC :: INT, PRESENT

    !I/O
    INTEGER, INTENT(OUT) :: status
    REAL(dp), INTENT(IN) :: relhum    ! [0-1]
    REAL(dp), INTENT(IN) :: temp      ! K
    REAL(dp), INTENT(IN) :: press     ! Pa
    LOGICAL, INTENT(IN), OPTIONAL :: l_hum_emac, l_psat_liquid, l_relhum_wmo

    !LOCAL
    !REAL(dp) :: omega_v, psat
    REAL(dp) :: psat, omega_v, omega_vs
    LOGICAL  :: z_l_hum_emac, z_l_psat_liquid, z_l_relhum_wmo

    status = 0

    IF (PRESENT(l_hum_emac)) THEN
      z_l_hum_emac = l_hum_emac
    ELSE
      z_l_hum_emac = .FALSE.
    ENDIF

    IF (PRESENT(l_psat_liquid)) THEN
      z_l_psat_liquid = l_psat_liquid
    ELSE
      z_l_psat_liquid = .FALSE.
    ENDIF

    IF (PRESENT(l_relhum_wmo)) THEN
      z_l_relhum_wmo = l_relhum_wmo
    ELSE
      z_l_relhum_wmo = .FALSE.
    ENDIF

    IF (z_l_relhum_wmo) THEN
      z_l_psat_liquid = .TRUE. ! WMO: always above liquid
    ENDIF

    ! WARNING or ERROR status if relhum out of bounds
    IF (relhum > 1._dp) THEN
      WRITE(*,'(A,ES10.3,A)') '   WARNING rel2spechum: relhum > 1 (',&
        relhum, ')'
      status = 1
    ELSEIF (relhum < 0._dp) THEN
      WRITE(*,'(A,ES10.3,A)') '   ERROR rel2spechum: relhum < 0 (',&
        relhum, '), spechum set to undefined'
      rel2spechum = FLAGGED_BAD
      status = 2
      RETURN
    ENDIF
    
    ! determine saturation partial pressure of water vapor
    IF (z_l_hum_emac) THEN
      IF (z_l_psat_liquid) THEN
        psat = tlucuaw(INT(temp*1000._dp))/MM_eps
      ELSE
        psat = tlucua(INT(temp*1000._dp))/MM_eps
      ENDIF
    ELSE
      psat = psat_mk(temp,z_l_psat_liquid)
    ENDIF

    IF (z_l_relhum_wmo) THEN
      !> water vapor saturation mass mixing ratio
      omega_vs = MM_eps * psat / (press - psat)

      !> calculate spechum with WMO definition of relhum
      rel2spechum = relhum * omega_vs / (1._dp + relhum * omega_vs)
    ELSE
      !> calculate spechum from relhum, taking into account that pressure
      !>   is defined for the moist atmosphere
      !   thus wrong: rel2spechum = MM_eps * relhum * psat / press
      omega_v = MM_eps * relhum * psat / (press - relhum * psat)
      rel2spechum = omega_v / (1._dp + omega_v)
    ENDIF

    ! Warning for relative humidity > 1 (supersaturation)
    ! Error if negative relative humidity
    IF (rel2spechum > 1._dp) THEN
      WRITE(*,'(A,ES10.3,A)') '   ERROR rel2spechum: spechum >= 1 (', &
        rel2spechum, '), spechum set to undefined'
      rel2spechum = FLAGGED_BAD
      status = 2
    ELSEIF (rel2spechum < 0._dp) THEN
      WRITE(*,'(A,ES10.3,A)') '   ERROR rel2spechum: spechum < 0 (', &
        rel2spechum, '), spechum set to undefined'
      rel2spechum = FLAGGED_BAD
      status = 2
      RETURN
    ENDIF

  END FUNCTION rel2spechum
! ----------------------------------------------------------------------------
! mz_hr_20100608-

!mz_hr_20160508+
! ----------------------------------------------------------------------------
  REAL(dp) FUNCTION spechum2mr(status, spechum)

    !> Calculates water mixing ratio (mol H2O)/(mol dryair) as function
    !>   of specific humidity.
    !> To calculate the water mixing ratio from relative humidity, first
    !>   convert relhum to spechum with the function rel2spechum

    USE messy_main_constants_mem, ONLY: FLAGGED_BAD, MM_eps

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! 1: warning, 2: error
    REAL(dp), INTENT(IN) :: spechum    ! [kg/kg]

    status = 0

    spechum2mr = 1._DP/MM_eps * spechum/(1._DP - spechum)

    IF (spechum2mr >= 1._dp) THEN
      WRITE(*,'(A,ES10.3,A)') 'ERROR spechum2mr: mixing ratio >= 1.0 (', &
        spechum2mr, ')'
      spechum2mr = FLAGGED_BAD
      status = 2
      RETURN
    ELSEIF (spechum2mr < 0._dp) THEN
      WRITE(*,'(A,ES10.3,A)') 'ERROR spechum2mr: mixing ratio < 0 (', &
        spechum2mr, ')'
      spechum2mr = FLAGGED_BAD
      status = 2
      RETURN
    ENDIF

  END FUNCTION spechum2mr
! ---------------------------------------------------------------------------- 
!mz_hr_20160508-

!mz_hr_20160508+
! ----------------------------------------------------------------------------
  REAL(dp) FUNCTION mr2spechum(status, mrH2O)

    !> Calculates specific humidity as function of water mixing ratio
    !> (mol H2O)/(mol dryair).
    !> To calculate the relative humidity from the water mixing ratio, use
    !>   spec2relhum afterwards

    USE messy_main_constants_mem, ONLY: FLAGGED_BAD, MM_eps

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! 1: warning, 2: error
    REAL(dp), INTENT(IN) :: mrH2O      ! [mol/mol]

    status = 0

    mr2spechum = mrH2O / (1._DP/MM_eps + mrH2O)

    IF (mr2spechum >= 1._dp) THEN
      WRITE(*,'(A,ES10.3,A)') 'ERROR mr2spechum: mixing ratio >= 1.0 (', &
        mr2spechum, ')'
      mr2spechum = FLAGGED_BAD
      status = 2
      RETURN
    ELSEIF (mr2spechum < 0._dp) THEN
      WRITE(*,'(A,ES10.3,A)') 'ERROR mr2spechum: mixing ratio < 0 (', &
        mr2spechum, ')'
      mr2spechum = FLAGGED_BAD
      status = 2
      RETURN
    ENDIF

  END FUNCTION mr2spechum
! ---------------------------------------------------------------------------- 
!mz_hr_20160508-

!mz_hr_20160503+
! ----------------------------------------------------------------------------
  REAL(DP) FUNCTION cair_q(spechum, temp, press, l_hum_emac)

  !> calculate cair from specific humidity

    USE messy_main_constants_mem, ONLY: R_gas, N_A, MM_eps

    IMPLICIT NONE

    INTRINSIC :: PRESENT

    !I/O
    LOGICAL, OPTIONAL, INTENT(IN) :: l_hum_emac ! treat humidity as in EMAC
    REAL(dp), INTENT(IN) :: spechum    ! kg/kg
    REAL(dp), INTENT(IN) :: temp       ! K
    REAL(dp), INTENT(IN) :: press      ! Pa

    !LOCAL
    LOGICAL :: z_l_hum_emac


    IF (PRESENT(l_hum_emac)) THEN
      z_l_hum_emac = l_hum_emac
    ELSE
      z_l_hum_emac = .FALSE.
    ENDIF

    IF (z_l_hum_emac) THEN
      !> calculate cair from specific humidity as in messy_mecca_e5
      cair_q = N_A/1.E6_DP * press/(R_gas*temp) &
                  * 1._DP/(1._DP + spechum*(1/MM_eps - 1._DP))
    ELSE
      cair_q = N_A/1.E6_DP * press/(R_gas*temp) &
        * (1._DP-spechum)/(1._DP + spechum*(1/MM_eps - 1._DP))
    ENDIF

  END FUNCTION cair_q
! ----------------------------------------------------------------------------
!mz_hr_20160503-

!mz_hr_20160513+
! ----------------------------------------------------------------------------
  REAL(DP) FUNCTION cair_c(c_H2O, temp, press, l_hum_emac)

  !> calculate cair from water vapor concentration

    USE messy_main_constants_mem, ONLY: R_gas, N_A, MM_eps

    IMPLICIT NONE

    INTRINSIC :: PRESENT

    !I/O
    LOGICAL, OPTIONAL, INTENT(IN) :: l_hum_emac ! treat humidity as in EMAC
    REAL(dp), INTENT(IN) :: c_H2O      ! kg/kg ! op_pj_20161006: [cm^-3] ???
    REAL(dp), INTENT(IN) :: temp       ! K
    REAL(dp), INTENT(IN) :: press      ! Pa

    !LOCAL
    LOGICAL :: z_l_hum_emac
    REAL(dp) :: spechum      ! tmp variable for l_hum_emac
    REAL(dp) :: mr_H2O_moist ! tmp variable for l_hum_emac


    IF (PRESENT(l_hum_emac)) THEN
      z_l_hum_emac = l_hum_emac
    ELSE
      z_l_hum_emac = .FALSE.
    ENDIF

    IF (z_l_hum_emac) THEN
      ! calculate cair from specific humidity as in messy_mecca_e5, derived
      !   from original formula as in cair_q
      mr_H2O_moist = c_H2O * 1.E6_DP/N_A * R_gas*temp/press
      spechum = (1._dp/(1._dp-mr_H2O_moist) - 1._dp)/(1._dp/MM_eps + &
        1._dp/(1._dp-mr_H2O_moist) - 1._dp)
      cair_c = N_A/1.E6_DP * press/(R_gas*temp) * 1._dp/(1._dp + &
        (1._dp/MM_eps-1._dp)*spechum)
      !print *, "EMAC cair_c = ", cair_c
    ELSE
      cair_c = N_A/1.E6_DP * press/(R_gas*temp) - c_H2O
      !print *, "norm cair_c = ", cair_c
    ENDIF

  END FUNCTION cair_c
! ----------------------------------------------------------------------------
!mz_hr_20160513-

! mz_ab_20100517+
! ----------------------------------------------------------------------------
  FUNCTION remap_bounds1(lb1,array) RESULT(ptr)
    INTEGER,                   INTENT(IN)          :: lb1
    REAL(dp), DIMENSION(lb1:), INTENT(IN), TARGET  :: array
    REAL(dp), DIMENSION(:),                POINTER :: ptr
    ptr => array
  END FUNCTION remap_bounds1
  
  FUNCTION remap_bounds2(lb1,lb2,array) RESULT(ptr)
    INTEGER,                        INTENT(IN)          :: lb1,lb2
    REAL(dp), DIMENSION(lb1:,lb2:), INTENT(IN), TARGET  :: array
    REAL(dp), DIMENSION(:,:),                   POINTER :: ptr
    ptr => array
  END FUNCTION remap_bounds2
  
  FUNCTION remap_bounds3(lb1,lb2,lb3,array) RESULT(ptr)
    INTEGER,                             INTENT(IN)          :: lb1,lb2,lb3
    REAL(dp), DIMENSION(lb1:,lb2:,lb3:), INTENT(IN), TARGET  :: array
    REAL(dp), DIMENSION(:,:,:),                      POINTER :: ptr
    ptr => array
  END FUNCTION remap_bounds3
  
  FUNCTION remap_bounds4(lb1,lb2,lb3,lb4,array) RESULT(ptr)
    INTEGER,                                  INTENT(IN)  :: lb1,lb2,lb3,lb4
    REAL(dp), DIMENSION(lb1:,lb2:,lb3:,lb4:), INTENT(IN), TARGET  :: array
    REAL(dp), DIMENSION(:,:,:,:),                         POINTER :: ptr
    ptr => array
  END FUNCTION remap_bounds4
! --------------------------------------------------------------------- 
! mz_ab_20100517-

! mz_ab_20100831+
! --------------------------------------------------------------------- 
  SUBROUTINE full2half(full,half,press,pressi)

    ! This subroutine returns a variable on half level pressures
    ! given a 3-D variable (but single jrow) on full level pressures.
    ! 
    ! half(:,k) lies between full(:,k) and full(:,k+1), so
    ! it is the lower boundary. Note that pressi is defined for upper
    ! boundary, therefore we use pressi(k+1) below.

    IMPLICIT NONE
    INTRINSIC :: size

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: full &
                                            , press, pressi
    REAL(dp), DIMENSION(:,:), INTENT(OUT) :: half 

    INTEGER :: nlev, k

    nlev = SIZE(press,2)

    DO k=1, nlev-1
       half(:,k)=(pressi(:,k+1)-press (:,k  ))/&
                   (press (:,k+1)-press (:,k  )) * full(:,k)   + &
                   (press (:,k+1)-pressi(:,k+1))/&
                   (press (:,k+1)-press (:,k  )) * full(:,k+1)
    ENDDO

    ! It is not possible to calculate this for lowest level, so
    ! this is set to the lowest-but-1 level
    half(:,nlev)=half(:,nlev-1)

  END SUBROUTINE full2half
! --------------------------------------------------------------------- 

! mz_ab_20111125+
! --------------------------------------------------------------------- 
    SUBROUTINE CalcTimeAngles(utsec, sda, ssa, rlt, csza, &
                              sza, daylen_sec, rlat, phi)

      USE messy_main_constants_mem, ONLY: DTR

      IMPLICIT NONE
      INTRINSIC :: ACOS, COS, SIN, SIZE

      REAL(dp), INTENT(IN)  :: utsec, daylen_sec
      REAL(dp), INTENT(IN)  :: sda
      REAL(dp), INTENT(IN)  :: rlat(:), phi(:)

      REAL(dp), INTENT(OUT) :: rlt(:), ssa(:)       &
                             , csza(:,:), sza(:,:)

      INTEGER :: lat_dim, lon_dim, l, m

      lat_dim = SIZE(rlat,1)
      lon_dim = SIZE(phi,1)

      LONLOOP: DO l = 1, lon_dim

         ! Solar siderial angle

         ssa(l)    = phi(l) + 360._dp*(utsec/daylen_sec) + 180._dp
         IF ( ssa(l) >= 360.0_dp ) ssa(l) = ssa(l) - 360.0_dp

         ! Local time

         rlt(l) = 180.0_dp + ssa(l)
         IF ( rlt(l) > 360.0_dp ) rlt(l) = rlt(l) - 360.0_dp
         rlt(l) = rlt(l)*DTR

         LATLOOP: DO m = 1 , lat_dim

            ! Calculate solar zenith angle

            csza(m,l)  = -COS(rlat(m))*COS(sda)*COS(rlt(l))+ &
                 SIN(rlat(m))*SIN(sda)
            sza(m,l)   =  ACOS(csza(m,l))

         END DO LATLOOP

      END DO LONLOOP


    END SUBROUTINE CalcTimeAngles
! --------------------------------------------------------------------- 
! mz_ab_20111125-

! --------------------------------------------------------------------- 
  SUBROUTINE Spline1D(X,Y,N,YP1,YPN,Y2,natspline)
    ! Numerical recipies Spline1d routine, see also partner split1d
    !
    ! GIVEN ARRAYS X AND Y OF LENGTH N CONTAINING A TABULATED FUNCTION
    ! IE YI=F(XI), WITH X1 < X2 < ... < XN, AND GIVEN VALUES YP1 AND
    ! YPN FOR THE FIRST DERIVATIVE OF THE INTERPOLATIONG FUNCTIONS AT
    ! POINTS 1 AND N RESPECTIVELY, THIS ROUTINE RETURNS AN ARRAY Y2 OF
    ! LENGTH N WHICH CONTAINS THE SECOND DERIVATIVES OF THE INTERPOLATING
    ! FUNCTION AT THE TABULATED POINTS XI.  IF YP2 AND/OR YPN ARE EQUAL
    ! TO 1E30 OR LARGER, THE ROUTINE IS SIGNALLED TO SET THE CORRESPONDING
    ! BOUNDARY CONDITIONS FOR A NATURAL SPLINE, WITH ZERO SECOND DERIVATIVES
    ! ON THAT BOUNDARY.

    IMPLICIT NONE
    INTRINSIC :: PRESENT

    REAL(dp), INTENT(IN) :: X(:), Y(:)
    INTEGER,  INTENT(IN) :: n
    REAL(dp), INTENT(IN) :: yp1,ypn
    REAL(dp), INTENT(OUT) :: Y2(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: natspline

    ! Local
    INTEGER, PARAMETER :: NMAX=400
    REAL(dp) :: U(NMAX)
!!$    REAL(dp) :: xa, y2a, ya
    REAL(dp) :: sig,p,qn,un
    INTEGER  :: i,k
!!$    INTEGER  :: klo,khi
!!$    REAL(dp) :: h,a,b
    LOGICAL  :: lnatspline = .FALSE.

    IF (PRESENT(natspline)) THEN
       lnatspline=natspline
    ELSE
       lnatspline=.FALSE.
    ENDIF

    IF (lnatspline) THEN
       Y2(1)=0._dp
       U(1)=0._dp
    ELSE
       Y2(1)=-0.5_dp
       U(1)=(3.0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
    ENDIF

    DO I=2,N-1
       SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
       P=SIG*Y2(I-1)+2._dp
       Y2(I)=(SIG-1.0_dp)/P
       U(I)=(6._dp*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))   &
            /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
    ENDDO

    IF (lnatspline) THEN
       QN=0._dp
       UN=0._dp
    ELSE
       QN=0.5_dp
       UN=(3._dp/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
    ENDIF

    Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1._dp)
    DO K=N-1,1,-1
       Y2(K)=Y2(K)*Y2(K+1)+U(K)
    ENDDO

    RETURN

  END SUBROUTINE Spline1D
! --------------------------------------------------------------------- 

! --------------------------------------------------------------------- 
  SUBROUTINE Splint1D(XA,YA,Y2A,N,X,Y,status)

    ! Numerical recipies spline routine. See also Spline1d
      
      ! GIVEN THE ARRAYS XA AND YA OF LENGTH N, WHICH TABULATE A FUNCTION
      ! (WITH THE XAI'S IN ORDER), AND GIVEN THE ARRAY Y2A, WHICH IS THE
      ! OUTPUT FROM SPLINE ABOVE, AND GIVEN A VALUE OF X, THIS ROUTINE
      ! RETURNS A CUBIC-SPLINE INTERPOLATED VALUE Y.
            
    REAL(dp), INTENT(IN)  :: xa(:), ya(:), y2a(:)
    INTEGER,  INTENT(IN)  :: n
    REAL(dp), INTENT(IN)  :: X
    REAL(dp), INTENT(OUT) :: Y
    INTEGER,  INTENT(OUT) :: status

!!$    INTEGER  :: nmax
!!$    REAL(dp) :: Y2,U
!!$    REAL(dp) :: yp1,ypn,sig,p,qn,un
!!$    INTEGER  :: i,k
    INTEGER  :: k
    INTEGER  :: klo,khi
    REAL(dp) :: h,a,b

    status = 1 ! ERROR
      
    KLO=1
    KHI=N
    DO WHILE (KHI-KLO > 1) 
       K=(KHI+KLO)/2
       IF (XA(K) > X) THEN
          KHI=K
       ELSE
          KLO=K
       ENDIF
    ENDDO
    
    H=XA(KHI)-XA(KLO)
    IF (H == 0 ) RETURN
    A=(XA(KHI)-X)/H
    B=(X-XA(KLO))/H
    Y=A*YA(KLO)+B*YA(KHI) +                                  &
         ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.0
    
    status = 0 ! NO ERROR
      
    RETURN
      
  END SUBROUTINE Splint1D
!---------------------------------------------------------------------
! mz_ab_20100831-

!---------------------------------------------------------------------
  FUNCTION find_next_free_unit(istart,istop) RESULT(unit)

    ! I/O
    INTEGER, INTENT(IN)  :: istart, istop
    INTEGER :: unit

    ! LOCAL
    LOGICAL        :: found, opened
    INTEGER        :: i
    CHARACTER(256) :: info

    found = .FALSE.
    DO i=istart,istop
       INQUIRE(unit=i,opened=opened)
       IF (.NOT.opened) THEN
          unit = i
          found = .TRUE.
          EXIT
       END IF
    END DO

    IF (.NOT. found) THEN
       WRITE(info,'(a,i2.2,a,i2.2,a)') &
         'No unit in range <',istart,':',istop,'> free.'
       !CALL error_bi(info,'find_next_free_unit')
    END IF

  END FUNCTION find_next_free_unit
! --------------------------------------------------------------------------

! mz_rs_20101118+
! op_ck_20131218+
!!$! --------------------------------------------------------------------------
!!$  ELEMENTAL REAL(dp) FUNCTION density(press, temp, sphum)
!!$
!!$    ! calculate air density in kg/m3
!!$
!!$    USE messy_main_constants_mem, ONLY: R_gas, M_air, M_H2O
!!$
!!$    REAL(dp), INTENT(in) :: press
!!$    REAL(dp), INTENT(in) :: temp
!!$    REAL(dp), INTENT(in) :: sphum
!!$
!!$    density = press /( R_gas * &
!!$      temp * (1._dp +(M_air / M_H2O -1._dp) *sphum))
!!$
!!$  END FUNCTION density
!!$!---------------------------------------------------------------------------
! --------------------------------------------------------------------------
  ELEMENTAL REAL(dp) FUNCTION mass_density(press, temp, sphum)

    ! calculate (wet) air density in kg m-3

    USE messy_main_constants_mem, ONLY: rd, vtmpc1

    IMPLICIT NONE

    REAL(dp), INTENT(in) :: press
    REAL(dp), INTENT(in) :: temp
    REAL(dp), INTENT(in) :: sphum

    mass_density = press / (rd * temp * (1._dp + vtmpc1 * sphum))

  END FUNCTION mass_density
!---------------------------------------------------------------------------
! op_ck_20131218-

! --------------------------------------------------------------------------
  ELEMENTAL REAL(dp) FUNCTION layerthickness(geopot_u, geopot_l)

    USE messy_main_constants_mem,  ONLY: g

    REAL(dp), INTENT(in) :: geopot_u
    REAL(dp), INTENT(in) :: geopot_l
    !LOCAL

    layerthickness  = (geopot_u-geopot_l)/ g
  END FUNCTION layerthickness
! --------------------------------------------------------------------------
! mz_rs_20101118-

! mz_bk_20110707+
  ! --------------------------------------------------------------------------
  ! str2num : converts string to number of kinds real(dp), real(sp) or integer
  ! procedures : str2num(character(len=*) str, real(dp) out, integer err)
  !              str2num(character(len=*) str, real(sp) out, integer err)
  !              str2num(character(len=*) str, integer out, integer err)

  SUBROUTINE str2num_real_dp(str, out, err)
    IMPLICIT NONE
    INTRINSIC :: PRESENT
    CHARACTER(LEN=*), INTENT(in)   :: str
    REAL(dp), INTENT(out)          :: out
    INTEGER, INTENT(out), OPTIONAL :: err

    !LOCAL
    INTEGER :: status

    IF (PRESENT(err)) err = 1

    READ(unit=str,fmt=*,IOSTAT=status) out
    IF (status /= 0) RETURN
    
    IF (PRESENT(err)) err=0
  END SUBROUTINE str2num_real_dp

  SUBROUTINE str2num_real_sp(str, out, err)
    IMPLICIT NONE
    INTRINSIC :: PRESENT
    CHARACTER(LEN=*), INTENT(in)   :: str
    REAL(sp), INTENT(out)          :: out
    INTEGER, INTENT(out), OPTIONAL :: err

    !LOCAL
    INTEGER :: status

    IF (PRESENT(err)) err = 1

    READ(unit=str,fmt=*,IOSTAT=status) out
    IF (status /= 0) RETURN
    
    IF (PRESENT(err)) err=0
  END SUBROUTINE str2num_real_sp

  SUBROUTINE str2num_integer(str, out, err)
    IMPLICIT NONE
    INTRINSIC :: PRESENT
    CHARACTER(LEN=*), INTENT(in)   :: str
    INTEGER, INTENT(out)           :: out
    INTEGER, INTENT(out), OPTIONAL :: err

    !LOCAL
    INTEGER :: status

    IF (PRESENT(err)) err = 1

    READ(unit=str,fmt=*,IOSTAT=status) out
    IF (status /= 0) RETURN
    
    IF (PRESENT(err)) err=0
  END SUBROUTINE str2num_integer
  ! --------------------------------------------------------------------------
! mz_bk_20110707-

! cy_ma_20120625+
!************************************************************************
! Marina Astitha - August 2009
SUBROUTINE ERRFUNC(X,ERESULT,JINT)
  ! This packet evaluates erf(x), erfc(x), and exp(x*x)*erfc(x)
  ! for a real argument  x.  It contains three function type
  ! subprograms: erf, erfc, and erfcx (or derf, derfc, and derfcx),
  ! and one subroutine type subprogram, calerf.  The calling
  ! statements for the primary entries are:
  
  ! y=erf(x)     (or   y=derf(x)),
  ! y=erfc(x)    (or   y=derfc(x)),
  ! and
  ! y=erfcx(x)   (or   y=derfcx(x)).
  
  ! Function                     Parameters for calerf
  ! Call              Arg                  Result          Jint
  ! 
  ! erf(arg)      any real argument         erf(arg)          0
  ! erfc(arg)     abs(arg)  <  xbig        erfc(arg)          1
  ! erfcx(arg)    xneg  <  arg  <  xmax   erfcx(arg)          2
  ! The main computation evaluates near-minimax approximations:
  ! from "Rational Chebyshev Approximations for the Error Function"
  ! by W. J. Cody, Math. Comp., 1969, pp. 631-638.  This
  ! transportable program uses rational functions that theoretically
  ! approximate  erf(x)  and  erfc(x)  to at least 18 significant
  ! decimal digits.  The accuracy achieved depends on the arithmetic
  ! system, the compiler, the intrinsic functions, and proper
  ! selection of the machine-dependent constants.
  
  ! Explanation of machine-dependent constants:
  ! xmin   = The smallest positive floating-point number.
  ! xinf   = The largest positive finite floating-point number.
  ! xneg   = The largest negative argument acceptable to erfcx;
  ! the negative of the solution to the equation
  ! 2*exp(x*x) = xinf.
  ! xsmall = Argument below which erf(x) may be represented by
  ! 2*x/sqrt(pi)  and above which  x*x  will not underflow.
  ! A conservative value is the largest machine number x
  ! such that   1.0 + x = 1.0   to machine precision.
  ! xbig   = Largest argument acceptable to erfc;  solution to
  ! the equation:  w(x)* (1-0.5/x**2) = xmin,  where
  ! w(x) = exp(-x*x)/[x*sqrt(pi)].
  ! xhuge  = Argument above which  1.0 - 1/(2*x*x) = 1.0  to
  ! machine precision.  a conservative value is
  ! 1/[2*sqrt(xsmall)]
  ! xmax   = Largest acceptable argument to erfcx; the minimum
  ! of xinf and 1/[sqrt(pi)*xmin].
  ! Approximate values for some important machines are:
  !                       xmin       xinf        xneg     xsmall
  ! CDC 7600      (s.p.)  3.13e-294   1.26e+322   -27.220  7.11e-15
  ! Cray-1        (s.p.)  4.58e-2467  5.45e+2465  -75.345  7.11e-15
  ! IEEE (IBM/XT,
  ! Sun, etc.)  (s.p.)  1.18e-38    3.40e+38     -9.382  5.96e-8
  ! IEEE (IBM/XT,
  ! Sun, etc.)  (d.p.)  2.23d-308   1.79d+308   -26.628  1.11d-16
  ! IBM 195       (d.p.)  5.40d-79    7.23e+75    -13.190  1.39d-17
  ! Univac 1108   (d.p.)  2.78d-309   8.98d+307   -26.615  1.73d-18
  ! Vax d-format  (d.p.)  2.94d-39    1.70d+38     -9.345  1.39d-17
  ! Vax g-format  (d.p.)  5.56d-309   8.98d+307   -26.615  1.11d-16
  
  !                        xbig       xhuge       xmax
  ! CDC 7600      (s.p.)  25.922      8.39e+6     1.80x+293
  ! Cray-1        (s.p.)  75.326      8.39e+6     5.45e+2465
  ! IEEE (IBM/XT,
  ! Sun, etc.)  (s.p.)   9.194      2.90e+3     4.79e+37
  ! IEEE (IBM/XT,
  ! Sun, etc.)  (d.p.)  26.543      6.71d+7     2.53d+307
  ! IBM 195       (d.p.)  13.306      1.90d+8     7.23e+75
  ! Univac 1108   (d.p.)  26.582      5.37d+8     8.98d+307
  ! Vax d-format  (d.p.)   9.269      1.90d+8     1.70d+38
  ! Vax g-format  (d.p.)  26.569      6.71d+7     8.98d+307
  
  ! Error returns:
  ! The program returns  erfc = 0      for  arg  >=  xbig;
  ! erfcx = xinf  for  arg  <  xneg;
  ! and
  ! erfcx = 0     for  arg  >=  xmax.
  ! Intrinsic functions required are:
  ! abs, aint, exp
  
  ! Author: W. J. Cody
  ! Mathematics And Computer Science Division
  ! Argonne National Laboratory
  ! Argonne, IL 60439
  ! Latest modification: March 19, 1990
  !
  ! Imported from the DEAD Model code (Charlie Zender)
  IMPLICIT NONE
  INTRINSIC :: ABS, AINT, EXP
  ! Parameters
  ! Mathematical constants
  REAL,PARAMETER::four=4.0e0
  REAL,PARAMETER::half=0.5e0
  REAL,PARAMETER::one=1.0e0
  REAL,PARAMETER::sixten=16.0e0
  REAL,PARAMETER::sqrpi=5.6418958354775628695e-1
  REAL,PARAMETER::thresh=0.46875e0
  REAL,PARAMETER::two=2.0e0
  REAL,PARAMETER::zero=0.0e0
  ! Machine-dependent constants
  REAL,PARAMETER::XBIG=9.194e0
  REAL,PARAMETER::XHUGE=2.90e3
  REAL,PARAMETER::XINF=3.40e+38
  REAL,PARAMETER::XMAX=4.79e37
  REAL,PARAMETER::XNEG=9.382e0
  REAL,PARAMETER::XSMALL=5.96e-8
  ! Coefficients for approximation to  erf  in first interval
  REAL,PARAMETER::a(5)=(/     &
       3.16112374387056560e+00,1.13864154151050156e+02, &
       3.77485237685302021e+02,3.20937758913846947e+03, &
       1.85777706184603153e-01 /)
  REAL,PARAMETER::b(4)=(/    &
       2.36012909523441209e+01,2.44024637934444173e+02,  &
       1.28261652607737228e+03,2.84423683343917062e+03 /)
  ! Coefficients for approximation to  erfc  in second interval
  REAL,PARAMETER::c(9)=(/   &
       5.64188496988670089e-01,8.88314979438837594e+00, &
       6.61191906371416295e+01,2.98635138197400131e+02, &
       8.81952221241769090e+02,1.71204761263407058e+03, &
       2.05107837782607147e+03,1.23033935479799725e+03, &
       2.15311535474403846e-08 /)
  REAL,PARAMETER::d(8)=(/   &
       1.57449261107098347e+01,1.17693950891312499e+02, &
       5.37181101862009858e+02,1.62138957456669019e+03, &
       3.29079923573345963e+03,4.36261909014324716e+03, &
       3.43936767414372164e+03,1.23033935480374942e+03 /)
  ! Coefficients for approximation to  erfc  in third interval
  REAL,PARAMETER::p(6)=(/  &
       3.05326634961232344e-01,3.60344899949804439e-01, &
       1.25781726111229246e-01,1.60837851487422766e-02, &
       6.58749161529837803e-04,1.63153871373020978e-02 /)
  REAL,PARAMETER::q(5)=(/   &
       2.56852019228982242e+00,1.87295284992346047e+00,&
       5.27905102951428412e-01,6.05183413124413191e-02,&
       2.33520497626869185e-03 /)
  ! Local
  INTEGER i,jint
  REAL arg,del,eresult,x,xden,xnum,y,ysq
  ! Main Code
  arg=x
  y=ABS(x)
  IF (y <= thresh) THEN
     ! Evaluate  erf  for  |x| <= 0.46875
     ysq=zero
     IF (y > xsmall) ysq=y*y
     xnum=a(5)*ysq
     xden=ysq
     DO i=1,3
        xnum=(xnum+a(i))*ysq
        xden=(xden+b(i))*ysq
     END DO
     eresult=x*(xnum+a(4))/(xden+b(4))
     IF (jint /= 0) eresult=one-eresult
     IF (jint == 2) eresult=EXP(ysq)*eresult
     GO TO 800
     ! Evaluate  erfc  for 0.46875 <= |x| <= 4.0
  ELSE IF (y <= four) THEN
     xnum=c(9)*y
     xden=y
     DO i=1,7
        xnum=(xnum+c(i))*y
        xden=(xden+d(i))*y
     END DO
     eresult=(xnum+c(8))/(xden+d(8))
     IF (jint /= 2) THEN
        ysq=AINT(y*sixten)/sixten
        del=(y-ysq)*(y+ysq)
        eresult=EXP(-ysq*ysq)*EXP(-del)*eresult
     END IF
     ! Evaluate  erfc  for |x| > 4.0
  ELSE
     eresult=zero
     IF (y >= xbig) THEN
        IF ((jint /= 2).or.(y >= xmax)) go to 300
        IF (y >= xhuge) THEN
           eresult=sqrpi/y
           go to 300
        END IF
     END IF
     ysq=one/(y*y)
     xnum=p(6)*ysq
     xden=ysq
     do i=1,4
        xnum=(xnum+p(i))*ysq
        xden=(xden+q(i))*ysq
     END do
     eresult=ysq*(xnum+p(5))/(xden+q(5))
     eresult=(sqrpi-eresult)/y
     IF (jint /= 2) THEN
        ysq=aint(y*sixten)/sixten
        del=(y-ysq)*(y+ysq)
        eresult=exp(-ysq*ysq)*exp(-del)*eresult
     END IF
  END IF
  ! Fix up for negative argument, erf, etc.
300 IF (jint == 0) THEN
     eresult=(half-eresult)+half
     IF (x < zero) eresult=-eresult
  ELSE IF (jint == 1) THEN
     IF (x < zero) eresult=two-eresult
  ELSE
     IF (x < zero) THEN
        IF (x < xneg) THEN
           eresult=xinf
        ELSE
           ysq=aint(x*sixten)/sixten
           del=(x-ysq)*(x+ysq)
           y=exp(ysq*ysq)*exp(del)
           eresult=(y+y)-eresult
        END IF
     END IF
  END IF
800 RETURN
END SUBROUTINE ERRFUNC
  ! --------------------------------------------------------------------------
! cy_ma_20120625-

  ! um_ak_20131111+ ! um_gg_2013111-
  ! --------------------------------------------------------------------------
  SUBROUTINE calc_hybrid_coeff(ke, vcoord, vcflat, irefatm, ivctype &
       , hyam, hyai, hybm, hybi)

    ! MESSy/SMCL
    USE messy_main_constants_mem,  ONLY: g, rd     

    IMPLICIT NONE
    INTRINSIC :: EXP, LOG, PRESENT, SQRT

    INTEGER,  INTENT(IN) :: ke        ! number of vertical levels
    INTEGER,  INTENT(IN) :: irefatm   ! indicator for reference atmosphere
    INTEGER,  INTENT(IN) :: ivctype   ! type of vertical coordinate
    REAL(dp), INTENT(IN) :: vcflat    ! coordinate where levels become flat
    ! incoming vertical coordinate ( [Pa] for ivctype==1, [m] else )
    ! For ivctype = 2 or 3, vcoord and vcflat are given in meters 
    ! defined on interface levels
    REAL(dp), INTENT(IN), DIMENSION(ke+1) :: vcoord    

    ! hybrid coefficients (at mid layer and interface)
    REAL(dp), INTENT(OUT), DIMENSION(ke),   OPTIONAL :: hyam, hybm
    REAL(dp), INTENT(OUT), DIMENSION(ke+1), OPTIONAL :: hyai, hybi
    
    ! LOCAL
    REAL(dp), DIMENSION(ke+1) :: vcoord_c, p0hl_c, hhl_c
    REAL(dp)                  :: vcflat_c, zgdrt, ztdbe, zbetf, t00
    REAL(dp), PARAMETER       :: t0sl=288.15_dp,dt0lp=42._dp,h_scal=1.e4_dp,delta_t=75._dp
    REAL(dp), PARAMETER       :: p0sl = 1.e5_dp
!!$    REAL(dp), DIMENSION(ke)   :: phyam, phybm
    REAL(dp), DIMENSION(ke+1) :: phyai, phybi

    INTEGER                   :: k

    IF (ivctype == 1) THEN
       DO k = 1, ke+1
          vcoord_c(k) = vcoord(k)
       ENDDO
       vcflat_c = vcflat

    ELSEIF ((ivctype == 2) .OR. (ivctype == 3)) THEN
       ! assume hsurf=hsurfs=0 -> ak=vcoord -> hhl=vcoord
       hhl_c(ke+1) = 0._dp
       DO  k = 1, ke
          hhl_c(k) = vcoord(k)
       ENDDO

       ! Compute the reference pressure at half levels
       IF (irefatm == 1) THEN
          zgdrt = g/rd/t0sl
          ztdbe = t0sl/dt0lp
          zbetf = 2.*dt0lp*zgdrt/t0sl
          vcflat_c  = EXP ( - ztdbe*(1._dp                              &
               - SQRT(1._dp - zbetf*vcflat)) )
          DO  k = 1, ke+1
             p0hl_c(k) = p0sl * EXP ( - ztdbe*(1._dp                    &
                  - SQRT(1._dp - zbetf*hhl_c(k))) )
          ENDDO
       ELSEIF (irefatm == 2) THEN
          t00   = t0sl-delta_t
          vcflat_c  = EXP ( - g/rd*h_scal/t00 * LOG(                    &
               (EXP(vcflat/h_scal)*t00 + delta_t)/(t00 + delta_t)) )
          DO  k = 1, ke+1
             p0hl_c(k) = p0sl * EXP ( - g/rd*h_scal/t00 * LOG( &
                  (EXP(hhl_c(k)/h_scal)*t00 + delta_t)/(t00 + delta_t)) )
          ENDDO
       ENDIF !irefatm
       vcoord_c(:) = p0hl_c(:)/p0sl

    ENDIF !ivctype

    DO k = 1, ke+1
       IF( vcoord_c(k) <= vcflat_c ) THEN
          phyai(k) = vcoord_c(k)*p0sl
          phybi(k) = 0._dp
       ELSE
          phyai(k) = vcflat_c*p0sl*(1.0_dp - vcoord_c(k))/(1._dp - vcflat_c)
          phybi(k) = (vcoord_c(k) - vcflat_c)/(1._dp - vcflat_c)
       ENDIF
    ENDDO
    IF (PRESENT(hyai)) hyai = phyai
    IF (PRESENT(hybi)) hybi = phybi

    IF (PRESENT(hyam)) THEN
       DO k = 1, ke
          hyam(k) = (phyai(k)+phyai(k+1))/2._dp
       ENDDO
    END IF
    IF (PRESENT(hybm)) THEN
       DO k = 1, ke
          hybm(k) = (phybi(k)+phybi(k+1))/2._dp
       END DO
    END IF

  END SUBROUTINE calc_hybrid_coeff
  ! --------------------------------------------------------------------------
  ! um_ak_20131111- ! um_gg_2013111-

  ! op_mm_20140610+
  !---------------------------------------------------------------------------
  ! This subroutine converts the field from the staggered COSMO Grid onto 
  ! the midpoints
  ! Select intdim = 1 for u-velocity, 2 for v-velocity, 3 for w-velocity
  ! Boundaries are treated as in cosmo src_output: Each loop goes from 
  ! 2:ie/2:je. Row 1 is set to Row 2 

  SUBROUTINE interpol_stag2mid(stag_field, mid_field, ie,je,ke, intdim)

    implicit none
  
    INTEGER, INTENT(IN)       ::ie,je,ke
    INTEGER, INTENT(IN)       ::intdim
    REAL(DP), INTENT(IN)      ::stag_field(:,:,:)
    REAL(DP), INTENT(OUT)     ::mid_field(:,:,:)

    !local
    INTEGER                   ::i,j,k
    
    mid_field=0.0_dp

    select case (intdim)

    case (1)
       do j= 1,je
          do k=1,ke
             do i=2, ie 
                mid_field(i,j,k)=0.5_dp*(stag_field(i,j,k)+stag_field(i-1,j,k))
             end do
          end do
       end do
       mid_field(1,:,:) = mid_field(2,:,:)
    case (2)
       do j=2,je
          do k=1,ke
             do i=1,ie
                mid_field(i,j,k)=0.5_dp*(stag_field(i,j,k)+stag_field(i,j-1,k))
             end do
          end do
       end do
       mid_field(:,1,:)=mid_field(:,2,:)
    case (3)
       do j=1,je
          do k=ke+1,2,-1
             do i=1,ie
! op_mm_20141027+
!!$             mid_field(i,j,k)=0.5_dp*(stag_field(i,j,k) &
!!$                  +stag_field(i,j,k-1))
                mid_field(i,j,k-1)=0.5_dp*(stag_field(i,j,k) &
                     +stag_field(i,j,k-1)) 
! op_mm_20141027-
             end do
          end do
       end do

    CASE DEFAULT

       print*, "Error in interpol_stag2mid. Unknown interpolation dimension "

    END SELECT

  END SUBROUTINE interpol_stag2mid
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! This subroutine converts fields from a midpoint grid to the staggered 
  ! COSMO Grid
  ! Select intdim = 1 for u-velocity, 2 for v-velocity
  
  SUBROUTINE interpol_mid2stag(mid_field, stag_field,ie,je,ke, intdim)

    implicit none

    INTEGER, INTENT(IN)       ::ie,je,ke
    INTEGER, INTENT(IN)       ::intdim
    REAL(DP), INTENT(OUT)     ::stag_field(:,:,:)
    REAL(DP), INTENT(IN)      ::mid_field(:,:,:)

    !local
    INTEGER                   ::i,j,k

    stag_field(:,:,:) = 0.0_dp ! op_pj_20170222

    select case (intdim)

    case (1)
       do j=1, je
          do k=1,ke
             do i=2,ie
                stag_field(i,j,k)=&
                     0.5_dp*(mid_field(i-1,j,k)+stag_field(i,j,k))
             end do
          end do
       end do
       stag_field(1,:,:) = stag_field(2,:,:)
    case (2)
       do j=2, je
          do k=1,ke
             do i=1,ie
                stag_field(i,j,k)= &
                     0.5_dp*(mid_field(i,j,k)+mid_field(i,j-1,k))
             end do
          end do
       end do
       stag_field(:,1,:) = stag_field(:,2,:)

    CASE DEFAULT

       print*, "Error in interpol_mid2stag. Unknown interpolation dimension "

    END SELECT

  END SUBROUTINE interpol_mid2stag
  !----------------------------------------------------------------------------
  ! op_mm_20140610-

! um_ak_20140814+
!==============================================================================
! adopted from COSMO model
!==============================================================================

!------------------------------------------------------------------------------

SUBROUTINE tautsp2D (TAU, GTAU, NTAU, NI, IMIN, IMAX, NTAUMAX, GAMMA,       &
                     S, BREAK, COEF, L_VEC, IFLAG)

!------------------------------------------------------------------------------
!
! Description:
!   Computes the tension spline for abscissae TAU(i) and ordinates GTAU(i)
!   (i = 1,...,NTAU), but vectorized over many "columns"
!
! From: Carl DeBoor: A Practical Guide to Splines
!
! Parameters:
!   TAU       Sequence of data points; must be strictly increasing (abscissae)
!   GTAU      Corresponding sequence of function values (ordinates)
!   NTAU      Number of data points; must be at least 4
!
!   NI        Number of columns
!   IMIN,IMAX Start- and end-index of columns
!   NTAUMAX   Maximum over the NTAUs in all columns
!
!   GAMMA     Indicates whether additional flexibility is desired:
!               0   :  no additional knots
!             0<  <3:  Under certain conditions on the given data at points
!                      i-1,i,i+1,i+2, a knot is added in the i-th interval
!                      (i=2,...,NTAU-2).
!             3<  <6:  Same, except that knots might also be added in intervals
!                      in which an inflection point would be permitted
!             The interpolant gets rounded with increasing GAMMA
!
! Output:
!   BREAK,COEF,L,K give the polynomial representation of the interpolant.
!   For BREAK(i) <= X <= BREAK(i+1) the interpolant has the following form:
!
!       F(X)=COEF(1,i)+DX(COEF(2,i)+DX/2(COEF(3,i)+DX/3(COEF(4,i)))
!                        with DX = X-BREAK(i) and i=1,...,L
! Workspace
!   S(NTAU,6)  WORK-ARRAY
!
! Return Codes:
!  IFLAG = 0  No Error
!  IFLAG = 2  Wrong Input
!
! Method:
!
!==============================================================================

INTEGER,  INTENT(IN)   ::  NI, IMIN, IMAX, NTAUMAX

INTEGER,  INTENT(IN)   ::  NTAU(NI)

REAL(DP), INTENT(IN)   :: &
    GTAU(NI,NTAUMAX),                    &
    TAU (NI,NTAUMAX),                    &
    GAMMA

REAL(DP),  INTENT(OUT)  :: &
    BREAK(NI,*),                         &
    COEF (NI,4,*),                       &
    S    (NI,NTAUMAX,6)

INTEGER, INTENT(OUT)  :: &
    L_VEC(NI),                           &
    IFLAG

! Local Variables
INTEGER  I,J,K,L,METHOD,NTAUM1, mb_err_indx_i,mb_err_indx_k

REAL(DP) C,D,DEL,DENOM,DIVDIF,ENTRY,FACTR2,GAM,      &
                  ONEMG3,SIXTH,TEMP,X,ALPH

REAL(DP)                     :: &
    RATIO_VEC   (NI),                    &
    Z_VEC       (NI),                    &
    ZETA_VEC    (NI),                    &
    ZT2_VEC     (NI),                    &
    ALPHA_VEC   (NI),                    &
    FACTOR_VEC  (NI),                    &
    ONEMZT_VEC  (NI),                    &
    ENTRY3      (NI)

!==============================================================================
!
   ALPH(X)= MIN (1.0_dp,ONEMG3/X)
!
!   Test of input parameters
!
      mb_err_indx_i = -1
      DO i = IMIN, IMAX
         IF (NTAU(i) .LT. 4) mb_err_indx_i = i
      ENDDO

      IF ( mb_err_indx_i /= -1 ) THEN
        WRITE (*,'(A,2I4,A)') '  NTAU =', NTAU(mb_err_indx_i), mb_err_indx_i, &
                              ':  MUST BE GREATER THAN 4'
        GO TO 999
      ENDIF

!
!   Computation of Delta(TAU) and of 1. and 2. derivation of data
!
      mb_err_indx_i = -1
      mb_err_indx_k = -1

      DO k = 1, NTAUMAX
         DO i = IMIN, IMAX
            IF (k <= NTAU(i)-1) THEN
               S(i,k,1)=TAU(i,k+1)-TAU(i,k)
               IF (S(i,k,1) .LE. 0.0_dp) THEN
                  mb_err_indx_i = i
                  mb_err_indx_k = k
               ELSE
                  S(i,k+1,4) = (GTAU(i,k+1) - GTAU(i,k))/S(i,k,1)
               ENDIF
            ENDIF
         ENDDO

         IF ( mb_err_indx_i /= -1 .OR. mb_err_indx_k /= -1) THEN
            WRITE (*,'(A,2I3,A,2E15.6,A)')                                    &
               ' Point ',mb_err_indx_i, mb_err_indx_k, ' and the following ', &
               TAU(mb_err_indx_i,mb_err_indx_k), TAU(mb_err_indx_i,mb_err_indx_k+1), &
               ' are in the wrong order! '
            GO TO 999
         ENDIF
      ENDDO

      DO k = 1, NTAUMAX
         DO i = IMIN, IMAX
            IF (k >= 2 .AND. k <= NTAU(i)-1) THEN
               S(i,k,4) = S(i,k+1,4) - S(i,k,4)
            ENDIF
         ENDDO
      ENDDO
!
!   Construct system of equations for 2. derivatives at TAU. At each interior
!   data point there is one continuity equation; at the first and the last
!   interior data point there is an additional one for a total of NTAU
!   equations in NTAU unknowns.
!
      DO i = IMIN, IMAX
         S(i,2,2) = S(i,1,1)/3.0_dp
      ENDDO

      SIXTH = 1.0_dp/6.0_dp
      METHOD = 2
      GAM = GAMMA
      IF(GAM .LE. 0.0_dp) METHOD = 1
      IF(GAM .GT. 3.0_dp) THEN
        METHOD = 3
        GAM = GAM - 3.0_dp
      ENDIF
      ONEMG3=1.0_dp - GAM/3.0_dp
!
!   Loops over k (interpolation points) and i (columns)
!
      DO k = 2, NTAUMAX
         DO i = IMIN, IMAX
            IF ( k <= NTAU(i)-2 ) THEN

               Z_VEC(i)=0.5_dp
               IF (METHOD /= 1) THEN
                  IF ( ((METHOD == 2) .AND.                                   &
                    (S(i,k,4)*S(i,k+1,4) >= 0.0_dp)) .OR. (METHOD == 3) ) THEN
                     TEMP = ABS(S(i,k+1,4))
                     DENOM = ABS(S(i,k,4)) + TEMP
                     IF (DENOM /= 0.0_dp) THEN
                        Z_VEC(i) = TEMP/DENOM
                        IF(ABS(Z_VEC(i)-0.5_dp).LE.SIXTH) Z_VEC(i)=0.5_dp
                     ENDIF
                  ENDIF
               ENDIF
               S(i,k,5) = Z_VEC(i)
!
!   Set up part of the k-th equation which depends on the k-th interval
!
               IF (Z_VEC(i)-0.5 .LT. 0._dp) THEN
                  ZETA_VEC(i) = GAM*Z_VEC(i)
                  ONEMZT_VEC(i) = 1.0_dp - ZETA_VEC(i)
                  ZT2_VEC(i) = ZETA_VEC(i)**2
                  ALPHA_VEC(i) = ALPH(ONEMZT_VEC(i))
                  FACTOR_VEC(i) = ZETA_VEC(i) /                       &
                                   (ALPHA_VEC(i)*(ZT2_VEC(i) - 1.0_dp) + 1.0_dp)
                  S(i,k,6) = ZETA_VEC(i)*FACTOR_VEC(i)/6.0_dp
                  S(i,k,2) = S(i,k,2) + S(i,k,1) *                  &
                                   ((1.0_dp-ALPHA_VEC(i)*ONEMZT_VEC(i))  &
                                     *FACTOR_VEC(i)*0.5_dp-S(i,k,6))
                  IF(S(i,k,2).LE.0.0) S(i,k,2) = 1.0_dp
                  S(i,k,3) = S(i,k,1)/6.0_dp
!
               ELSE IF (Z_VEC(i)-0.5 .EQ. 0._dp) THEN
!
                  S(i,k,2) = S(i,k,2) + S(i,k,1)/3.0_dp
                  S(i,k,3) = S(i,k,1)/6.0_dp
!
               ELSE
!
                  ONEMZT_VEC(i) = GAM*(1.0_dp - Z_VEC(i))
                  ZETA_VEC(i) = 1.0_dp - ONEMZT_VEC(i)
                  ALPHA_VEC(i) = ALPH(ZETA_VEC(i))
                  FACTOR_VEC(i) = ONEMZT_VEC(i) /                    &
                   (1.0_dp - ALPHA_VEC(i)*ZETA_VEC(i)*(1.0_dp + ONEMZT_VEC(i)))
                  S(i,k,6) = ONEMZT_VEC(i)*FACTOR_VEC(i)/6.0_dp
                  S(i,k,2) = S(i,k,2) + S(i,k,1)/3.0_dp
                  S(i,k,3) = S(i,k,6) * S(i,k,1)
               ENDIF
!
               IF (k == 2) THEN
                  S(i,1,5) = 0.5_dp
!
!   The first two equations enforce continuity of the first and of the third
!   derivative across TAU(2)
!
                  S(i,1,2) = S(i,1,1)/6.0_dp
                  S(i,1,3) = S(i,2,2)
                  ENTRY3(i) = S(i,2,3)

                  IF (Z_VEC(i)-0.5_dp .LT. 0._dp) THEN
                     FACTR2 = ZETA_VEC(i)*(ALPHA_VEC(i)                 &
                              *(ZT2_VEC(i)-1.0_dp)+1.0_dp)/             &
                        (ALPHA_VEC(i)*(ZETA_VEC(i)*ZT2_VEC(i)-1.0_dp) + 1.0_dp)
                     RATIO_VEC(i) = FACTR2*S(i,2,1)/S(i,1,2)
                     S(i,2,2) = FACTR2*S(i,2,1) + S(i,1,1)
                     S(i,2,3) = -FACTR2 * S(i,1,1)
!
                  ELSE IF (Z_VEC(i)-0.5_dp .EQ. 0._dp) THEN
!
                     RATIO_VEC(i) = S(i,2,1)/S(i,1,2)
                     S(i,2,2) = S(i,2,1) + S(i,1,1)
                     S(i,2,3) = -S(i,1,1)
!
                  ELSE
!
                     RATIO_VEC(i) = S(i,2,1)/S(i,1,2)
                     S(i,2,2) = S(i,2,1) + S(i,1,1)
                     S(i,2,3) = -S(i,1,1)*6.0_dp*ALPHA_VEC(i)*S(i,2,6)
                  ENDIF
!
!   At this point the first two equations read
!            DIAG(1)*X1 +    U(1)*X2 + ENTRY3*X3 = R(2)
!     -RATIO*DIAG(1)*X1 + DIAG(2)*X2 +   U(2)*X3 = 0
!   Eliminate first unknown from second equation
!
                  S(i,2,2) = RATIO_VEC(i)*S(i,1,3) + S(i,2,2)
                  S(i,2,3) = RATIO_VEC(i)*ENTRY3(i) + S(i,2,3)
                  S(i,1,4) = S(i,2,4)
                  S(i,2,4) = RATIO_VEC(i)*S(i,1,4)
!

               ELSE ! k > 2
!
                  S(i,k,2) = RATIO_VEC(i)*S(i,k-1,3) + S(i,k,2)
                  S(i,k,4) = RATIO_VEC(i)*S(i,k-1,4) + S(i,k,4)
               ENDIF  ! k == 2
!
!  Set up the part of the next equation which depends on the k-th interval
!
               IF (Z_VEC(i)-0.5_dp .LT. 0._dp) THEN
                  RATIO_VEC(i) = -S(i,k,6)*S(i,k,1)/S(i,k,2)
                  S(i,k+1,2) = S(i,k,1)/3.0_dp
!
               ELSE IF (Z_VEC(i)-0.5_dp .EQ. 0._dp) THEN
!
                  RATIO_VEC(i) = -S(i,k,1)/(6.0_dp*S(i,k,2))
                  S(i,k+1,2) = S(i,k,1)/3.0_dp
!
               ELSE
!
                  RATIO_VEC(i) = -(S(i,k,1)/6.0_dp)/S(i,k,2)
                  S(i,k+1,2)   = S(i,k,1) *                            &
                               ((1.0 - ZETA_VEC(i)*ALPHA_VEC(i))*      &
                               0.5_dp*FACTOR_VEC(i)-S(i,k,6))
               ENDIF
!
!   End of Loop over k (interpolation points) and i (columns)
!
            ENDIF ! k <= NTAU(i)-2
         ENDDO ! i = IMIN, IMAX
      ENDDO ! k = 2, NTAUMAX

      DO i = IMIN, IMAX
         k = NTAU(i)-1
         S(i,k,5) = 0.5_dp

!
!   The last two equations enforce continuity of third derivative and of
!   first derivative across TAU (NTAU-1)
!
         ENTRY = RATIO_VEC(i)*S(i,k-1,3) + S(i,k,2) + S(i,k,1)/3.0_dp
         S(i,k+1,2) = S(i,k,1)/6.0_dp
         S(i,k+1,4) = RATIO_VEC(i)*S(i,k-1,4) + S(i,k,4)
         IF (Z_VEC(i)-0.5_dp .LT. 0._dp) THEN
            RATIO_VEC(i) = S(i,k,1) * 6.0_dp * S(i,k-1,6) *              &
                             ALPHA_VEC(i)/S(i,k-1,2)
            S(i,k,2) = RATIO_VEC(i)*S(i,k-1,3) +S(i,k,1) + S(i,k-1,1)
            S(i,k,3) = -S(i,k-1,1)
!
         ELSE IF (Z_VEC(i)-0.5_dp .EQ. 0._dp) THEN
!
            RATIO_VEC(i) = S(i,k,1)/S(i,k-1,2)
            S(i,k,2) = RATIO_VEC(i)*S(i,k-1,3) + S(i,k,1) + S(i,k-1,1)
            S(i,k,3) = -S(i,k-1,1)
!
         ELSE
!
            FACTR2 = ONEMZT_VEC(i) * (ALPHA_VEC(i) *                    &
                            (ONEMZT_VEC(i)**2-1.0_dp)+1.0_dp)  /        &
                            (ALPHA_VEC(i)*(ONEMZT_VEC(i)**3-1.0_dp)+1.0_dp)
            RATIO_VEC(i) = FACTR2*S(i,k,1)/S(i,k-1,2)
            S(i,k,2) = RATIO_VEC(i)*S(i,k-1,3) + FACTR2*S(i,k-1,1)  &
                         + S(i,k,1)
            S(i,k,3) = -FACTR2*S(i,k-1,1)
         ENDIF
!
!  At this point the last two equations read
!           DIAG(k)*Xk +      U(k)*Xk+1 = R(k)
!    -RATIO*DIAG(k)*Xk + DIAG(k+1)*Xk+1 = R(k+1)
!  Eliminate Xk from last equation
!
         S(i,k,4) = RATIO_VEC(i)*S(i,k-1,4)
         RATIO_VEC(i) = -ENTRY/S(i,k,2)
         S(i,k+1,2) = RATIO_VEC(i)*S(i,k,3) + S(i,k+1,2)
         S(i,k+1,4) = RATIO_VEC(i)*S(i,k,4) + S(i,k+1,4)
      ENDDO ! i = IMIN, IMAX

!
!   Back Substitution
!
      DO i = IMIN, IMAX
         S(i,NTAU(i),4) = S(i,NTAU(i),4)/S(i,NTAU(i),2)
      ENDDO

      DO k = NTAUMAX,2,-1
         DO i = IMIN, IMAX
            IF (k <= NTAU(i)-1) THEN
               S(i,k,4) = (S(i,k,4) - S(i,k,3)*S(i,k+1,4))/S(i,k,2)
            ENDIF
         ENDDO
      ENDDO

      DO i = IMIN, IMAX
         S(i,1,4) = (S(i,1,4) - S(i,1,3)*S(i,2,4) -              &
                       ENTRY3(i)*S(i,3,4))/S(i,1,2)
      ENDDO
!
!   Construct polynomial pieces
!
      DO i = IMIN, IMAX
         BREAK(i,1) = TAU(i,1)
         L_VEC(i) = 1
      ENDDO

      DO k = 1, NTAUMAX
         DO i = IMIN, IMAX
            IF ( k <= NTAU(i)-1) THEN
               L = L_VEC(i)
               COEF(i,1,L) = GTAU(i,k)
               COEF(i,3,L) = S(i,k,4)
               DIVDIF = (GTAU(i,k+1) - GTAU(i,k))/S(i,k,1)
               Z_VEC(i) = S(i,k,5)
!
               IF (Z_VEC(i)-0.5_dp .LT. 0._dp) THEN
                  ! US avoid division by 0, if Z_VEC is veeeery small
                  ! by treating it as 0
                  IF(ABS(Z_VEC(i)) < 1E-50_dp) THEN
                     COEF(i,2,L) = DIVDIF
                     COEF(i,3,L) = 0.0_dp
                     COEF(i,4,L) = 0.0_dp
                  ELSE
                     ZETA_VEC(i) = GAM*Z_VEC(i)
                     ONEMZT_VEC(i) = 1.0_dp -ZETA_VEC(i)
                     C = S(i,k+1,4)/6.0_dp
                     D = S(i,k,4)*S(i,k,6)
                     L = L + 1
                     DEL = ZETA_VEC(i)*S(i,k,1)
                     BREAK(i,L) = TAU(i,k) + DEL
                     ZT2_VEC(i) = ZETA_VEC(i)**2
!                    ALPHA_VEC(i) = ALPH(ONEMZT_VEC(i))
                     ALPHA_VEC(i) = MIN(1.0_dp,ONEMG3/ONEMZT_VEC(i))
                     FACTOR_VEC(i) = ONEMZT_VEC(i)**2*ALPHA_VEC(i)
                     COEF(i,1,L) = GTAU(i,k) + DIVDIF*DEL +              &
                          S(i,k,1)**2*(D*ONEMZT_VEC(i)*(FACTOR_VEC(i)- &
                             1.0_dp) + C*ZETA_VEC(i)*(ZT2_VEC(i) - 1.0_dp))
                     COEF(i,2,L) = DIVDIF + S(i,k,1) *                   &
                         (D*(1.0_dp-3.0_dp*FACTOR_VEC(i)) + C*(3.0*ZT2_VEC(i)- &
                                    1.0))
                     COEF(i,3,L) = 6.0_dp*(D*ALPHA_VEC(i)*ONEMZT_VEC(i)   &
                                      + C*ZETA_VEC(i))
                     COEF(i,4,L) = 6.0_dp*(C - D*ALPHA_VEC(i))/S(i,k,1)
                     COEF(i,4,L-1) = COEF(i,4,L) -6.0_dp*D*                 &
                                       (1.0_dp-ALPHA_VEC(i))/(DEL*ZT2_VEC(i))
                     COEF(i,2,L-1) = COEF(i,2,L) -                       &
                                   DEL*(COEF(i,3,L)-DEL*0.5_dp*COEF(i,4,L-1))
                  ENDIF
!
               ELSE IF (Z_VEC(i)-0.5_dp .EQ. 0._dp) THEN
!
                  COEF(i,2,L) = DIVDIF - S(i,k,1) *                      &
                                  (2.0*S(i,k,4) + S(i,k+1,4))/6.0_dp
                  COEF(i,4,L) = (S(i,k+1,4) - S(i,k,4))/S(i,k,1)
!
               ELSE
!
                  ONEMZT_VEC(i) = GAM*(1.0_dp - Z_VEC(i))
                  IF(ONEMZT_VEC(i).EQ.0.0_dp) THEN
                     COEF(i,2,L) = DIVDIF
                     COEF(i,3,L) = 0.0_dp
                     COEF(i,4,L) = 0.0_dp
                  ELSE
                     ZETA_VEC(i) = 1.0_dp - ONEMZT_VEC(i)
!                    ALPHA_VEC(i) = ALPH(ZETA_VEC(i))
                     ALPHA_VEC(i) = MIN(1.0_dp,ONEMG3/ZETA_VEC(i))
                     C = S(i,k+1,4)*S(i,k,6)
                     D = S(i,k,4)/6.0_dp
                     DEL = ZETA_VEC(i)*S(i,k,1)
                     BREAK(i,L+1) = TAU(i,k) + DEL
                     COEF(i,2,L) = DIVDIF -S(i,k,1)*(2.0_dp*D + C)
                     COEF(i,4,L) = 6.0_dp*(C*ALPHA_VEC(i) - D)/S(i,k,1)
                     L = L + 1
                     COEF(i,4,L) = COEF(i,4,L-1) + 6.0_dp *                 &
                         (1.0_dp-ALPHA_VEC(i))*C/(S(i,k,1)*ONEMZT_VEC(i)**3)
                     COEF(i,3,L) = COEF(i,3,L-1) + DEL* COEF(i,4,L-1)
                     COEF(i,2,L) = COEF(i,2,L-1) + DEL*(COEF(i,3,L-1)  &
                                     +DEL*0.5_dp*COEF(i,4,L-1))
                     COEF(i,1,L) = COEF(i,1,L-1) + DEL*(COEF(i,2,L-1)  &
                                     +DEL*0.5_dp*                      &
                           (COEF(i,3,L-1) + (DEL/3.0_dp)*COEF(i,4,L-1)))
                  ENDIF
               ENDIF
!
               L = L + 1
               BREAK(i,L) = TAU(i,k+1)
               L_VEC(i) = L
            ENDIF
         ENDDO ! i = IMIN, IMAX
      ENDDO ! k = 1, NTAUMAX

      IFLAG = 0

      RETURN
!
 999  IFLAG = 2

END SUBROUTINE tautsp2D

!==============================================================================
!==============================================================================
! um_ak_20140814-

!---------------------------------------------------------------------------
RECURSIVE SUBROUTINE quick_sort(list, order)

  ! sort routine to arrange array elements from smallest to largest
  !
  ! grabbed from A millers web site http://users.bigpond.net.au/amiller/
  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified by Alan Miller to include an associated integer array which gives
  ! the positions of the elements in the original order.
  ! mvr modified integer array to intent inout - may now be any integer 
  !     array that gets sorted along with associated real array

  implicit none

  REAL(DP), DIMENSION (:), INTENT(INOUT) :: list
  INTEGER,  DIMENSION (:), INTENT(INOUT) :: order

  ! Local variable
  INTEGER :: i

  CALL quick_sort_1(1, SIZE(list))

CONTAINS

  RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

    implicit none
    INTEGER, INTENT(IN) :: left_end, right_end

    ! Local variables
    INTEGER             :: i, j, itemp
    REAL(DP)            :: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 6

    IF (right_end < left_end + max_simple_sort_size) THEN
       ! Use interchange sort for small lists
       CALL interchange_sort(left_end, right_end)

    ELSE
       ! Use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1

       DO
          ! Scan list from left end until element >= reference is found
          DO
             i = i + 1
             IF (list(i) >= reference) EXIT
          END DO
          ! Scan list from right end until element <= reference is found
          DO
             j = j - 1
             IF (list(j) <= reference) EXIT
          END DO

          IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          ELSE IF (i == j) THEN
             i = i + 1
             EXIT
          ELSE
             EXIT
          END IF
       END DO

       IF (left_end < j) CALL quick_sort_1(left_end, j)
       IF (i < right_end) CALL quick_sort_1(i, right_end)
    END IF

  END SUBROUTINE quick_sort_1

  SUBROUTINE interchange_sort(left_end, right_end)

    implicit none
    INTEGER, INTENT(IN) :: left_end, right_end

    ! Local variables
    INTEGER             :: i, j, itemp
    REAL(DP)            :: temp

    DO i = left_end, right_end - 1
       DO j = i+1, right_end
          IF (list(i) > list(j)) THEN
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          END IF
       END DO
    END DO

  END SUBROUTINE interchange_sort

END SUBROUTINE quick_sort
!---------------------------------------------------------------------------

! mz_ht_20140923+
!---------------------------------------------------------------------------
Subroutine IDXSORT_DP (XVALT, IRNGT)
!   Ranks array XVALT into index array IRNGT, using merge-sort
! __________________________________________________________
!   This version is not optimized for performance, and is thus
!   not as difficult to read as some other ones.
!   Michel Olagnon - April 2000
! __________________________________________________________
! __________________________________________________________
  Real (dp), Dimension (:), Intent (In) :: XVALT
  Integer, Dimension (:), Intent (Out) :: IRNGT
  ! __________________________________________________________
  !
  Integer, Dimension (:), Allocatable :: JWRKT
  Integer :: LMTNA, LMTNC
  Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
  !
  NVAL = Min (SIZE(XVALT), SIZE(IRNGT))
  If (NVAL <= 0) Then
    Return
  End If
  !
  !  Fill-in the index array, creating ordered couples
  !
  Do IIND = 2, NVAL, 2
    If (XVALT(IIND-1) < XVALT(IIND)) Then
      IRNGT (IIND-1) = IIND - 1
      IRNGT (IIND) = IIND
    Else
      IRNGT (IIND-1) = IIND
      IRNGT (IIND) = IIND - 1
    End If
  End Do
  If (Modulo (NVAL, 2) /= 0) Then
    IRNGT (NVAL) = NVAL
  End If
  !
  !  We will now have ordered subsets A - B - A - B - ...
  !  and merge A and B couples into     C   -   C   - ...
  !
  Allocate (JWRKT(1:NVAL))
  LMTNC = 2
  LMTNA = 2
  !
  !  Iteration. Each time, the length of the ordered subsets
  !  is doubled.
  !
  Do
    If (LMTNA >= NVAL) Exit
    IWRKF = 0
    LMTNC = 2 * LMTNC
    IWRK = 0
    !
    !   Loop on merges of A and B into C
    !
    Do
      IINDA = IWRKF
      IWRKD = IWRKF + 1
      IWRKF = IINDA + LMTNC
      JINDA = IINDA + LMTNA
      If (IWRKF >= NVAL) Then
        If (JINDA >= NVAL) Exit
        IWRKF = NVAL
      End If
      IINDB = JINDA
      !
      !   Shortcut for the case when the max of A is smaller
      !   than the min of B (no need to do anything)
      !
      If (XVALT(IRNGT(JINDA)) <= XVALT(IRNGT(JINDA+1))) Then
        IWRK = IWRKF
        Cycle
      End If
      !
      !  One steps in the C subset, that we create in the final rank array
      !
      Do
        If (IWRK >= IWRKF) Then
          !
          !  Make a copy of the rank array for next iteration
          !
          IRNGT (IWRKD:IWRKF) = JWRKT (IWRKD:IWRKF)
          Exit
        End If
        !
        IWRK = IWRK + 1
        !
        !  We still have unprocessed values in both A and B
        !
        If (IINDA < JINDA) Then
          If (IINDB < IWRKF) Then
            If (XVALT(IRNGT(IINDA+1)) > XVALT(IRNGT(IINDB+1))) &
              & Then
              IINDB = IINDB + 1
              JWRKT (IWRK) = IRNGT (IINDB)
            Else
              IINDA = IINDA + 1
              JWRKT (IWRK) = IRNGT (IINDA)
            End If
          Else
            !
            !  Only A still with unprocessed values
            !
            IINDA = IINDA + 1
            JWRKT (IWRK) = IRNGT (IINDA)
          End If
        Else
          !
          !  Only B still with unprocessed values
          !
          IRNGT (IWRKD:IINDB) = JWRKT (IWRKD:IINDB)
          IWRK = IWRKF
          Exit
        End If
        !
      End Do
    End Do
    !
    !  The Cs become As and Bs
    !
    LMTNA = 2 * LMTNA
  End Do
  !
  !  Clean up
  !
  Deallocate (JWRKT)
  Return
  !
End Subroutine IDXSORT_DP
!---------------------------------------------------------------------------
! mz_ht_20140923-

! op_bk_20170111+
!---------------------------------------------------------------------------
SUBROUTINE get_name_domain(status, string, cname, domain)
! split a string into a channel name (cname) and a domain number (domain),
! when string is of kind: <channel name>_D<domain number (2digit)>
! e.g.  string -> tr_O3_D02
!       cname  -> tr_O3
!       domain -> 2        
!       status -> 0 (domain number present)

  USE messy_main_constants_mem, ONLY: global_patch ! op_pj_20180126

  IMPLICIT NONE

  INTRINSIC :: INDEX, LEN_TRIM, TRIM

  INTEGER, INTENT(OUT)          :: status
  CHARACTER(LEN=*), INTENT(IN)  :: string
  CHARACTER(LEN=*), INTENT(OUT) :: cname
  INTEGER, INTENT(OUT)          :: domain

  INTEGER :: idx, len
  CHARACTER(LEN=2) :: str_num

  ! init status -> error: no domain number in string
  cname = ''
! op_pj_20180126+
!!$  domain = -1
  domain = global_patch
! op_pj_20180126-
  idx = INDEX(string, '_D')
  ! test if domain number is specified
  IF (idx == 0) THEN
     cname = TRIM(string)
     status=1
     RETURN
  END IF
  len = LEN_TRIM(string)
  ! test two digit domain number
  IF (len-(idx+2)+1 /= 2) THEN
     cname = TRIM(string)
     status=2
     RETURN
  END IF
  str_num = TRIM(string(idx+2:len))
  ! status=1 if it is not possible to convert to integer, otherwise status=0
  CALL str2num(str_num,domain,status)
  IF (status /= 0) THEN
! op_pj_20180126+
!!$  domain = -1
     domain = global_patch
! op_pj_20180126-
     cname = TRIM(string)
     status=3
  ELSE
     cname = TRIM(string(1:idx-1))
     status = 0
  END IF

END SUBROUTINE get_name_domain
!---------------------------------------------------------------------------
! op_bk_20170111-

! ub_ak_20170406+
!---------------------------------------------------------------------------
elemental subroutine lcase(word)
  ! convert a word to lower case
  character (len=*) , intent(in out) :: word
  integer                            :: i,ic,nlen
  nlen = len(word)
  do i=1,nlen
     ic = iachar(word(i:i))
     if (ic >= 65 .and. ic <= 90) word(i:i) = achar(ic+32)
  end do
end subroutine lcase
!---------------------------------------------------------------------------
! ub_ak_20170406-

! op_pj_20180523+
!---------------------------------------------------------------------------
! workaround for Lahey/Fujitsu Compiler 8.1b
#ifdef LF
LOGICAL FUNCTION ISNAN(x)
  IMPLICIT NONE
  REAL(dp), INTENT(in) :: x
  ISNAN = (x/=x)
END FUNCTION ISNAN
#endif
!---------------------------------------------------------------------------
! op_pj_20180523-

! ************************************************************************
END MODULE messy_main_tools
! ************************************************************************
