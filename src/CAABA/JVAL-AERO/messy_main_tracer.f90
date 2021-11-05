! **********************************************************************
MODULE messy_main_tracer
! **********************************************************************

  ! THIS MODULE IS THE MAIN SMCL MODULE OF THE GENERIC MESSy-SUBMODEL
  ! 'TRACER'

  USE messy_main_constants_mem, ONLY: DP, SP &
       , STRLEN_MEDIUM, STRLEN_LONG, STRLEN_VLONG, STRLEN_XLONG &
       , STRLEN_ULONG & ! op_pj_20170508
       , BIG_DP, MH, MC, MN, MF, MNa, MO, MS, MCl, MBr, MI, MHg &
       , MD, M13C, M12C, M18O, MK, MMg, MCa

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE

  CHARACTER(len=*), PARAMETER, PUBLIC :: modstr = 'tracer'
  CHARACTER(len=*), PARAMETER, PUBLIC :: modver = '2.5'

  PUBLIC :: DP, SP ! mz_ab_20151129 added SP

  ! RANK OF TRACER INDEX
  INTEGER, SAVE, PUBLIC :: TRRANK = 3 ! DEFAULT, CAN BE CHANGED BY BM

  ! GLOBAL NAMELIST SWITCHES (CTRL)
  LOGICAL, SAVE, PUBLIC :: l_family     = .FALSE.
  LOGICAL, SAVE, PUBLIC :: l_pdef       = .FALSE.

  ! SPECIAL STRING LENGTHs
  INTEGER, PARAMETER, PUBLIC  :: STRLEN_TRSET  = 15
  INTEGER, PARAMETER, PUBLIC  :: STRLEN_FNAME  = 2*STRLEN_MEDIUM + 1         

  ! GENERAL SWITCHES
  INTEGER, PARAMETER, PUBLIC :: UNDEFINED = -1
  INTEGER, PARAMETER, PUBLIC :: OFF       =  0
  INTEGER, PARAMETER, PUBLIC :: ON        =  1
  ! ... AEROSOL MODEL METHOD
  INTEGER, PARAMETER, PUBLIC :: MODAL     = 2
  INTEGER, PARAMETER, PUBLIC :: BIN       = 3

  ! TYPE
  INTEGER, PARAMETER, PUBLIC :: SINGLE  = 0
  INTEGER, PARAMETER, PUBLIC :: FAMILY  = 1
  INTEGER, PARAMETER, PUBLIC :: ISOTOPE = 2
  INTEGER, PARAMETER         :: MAX_TYPE = 2

  ! 'HOSTING' MEDIUM
  INTEGER, PARAMETER, PUBLIC :: AIR        = 1
  INTEGER, PARAMETER, PUBLIC :: AEROSOL    = 2
  INTEGER, PARAMETER, PUBLIC :: CLOUD      = 3
  INTEGER, PARAMETER, PUBLIC :: OCEAN      = 4
  INTEGER, PARAMETER, PUBLIC :: LAKE       = 5
  INTEGER, PARAMETER, PUBLIC :: RIVER      = 6
  INTEGER, PARAMETER, PUBLIC :: LANDICE    = 7
  INTEGER, PARAMETER, PUBLIC :: SEAICE     = 8
  INTEGER, PARAMETER, PUBLIC :: VEGETATION = 9
  INTEGER, PARAMETER         :: MAX_MEDIUM = 9

  ! QUANTITY
  INTEGER, PARAMETER, PUBLIC :: AMOUNTFRACTION = 1 ! = MOLAR MIXING RATIO
  INTEGER, PARAMETER, PUBLIC :: NUMBERDENSITY  = 2
  INTEGER, PARAMETER, PUBLIC :: CONCENTRATION  = 3
  INTEGER, PARAMETER         :: MAX_QUANTITY   = 3

  ! um_ak_20120716+
  ! INITIAL TYPES
  INTEGER, PARAMETER, PUBLIC :: T_INI_ZERO = 0 ! initialize with zero or const
  INTEGER, PARAMETER, PUBLIC :: T_INI_FILE = 1 ! initialize from file/MMD
  INTEGER, PARAMETER, PUBLIC :: T_INI_USER = 2 ! to be initialized by user

  ! LATERAL BOUNDARY CONDITIONS
  INTEGER, PARAMETER, PUBLIC :: T_LBC_ZERO = 0 ! initialize lateral BC to zero
  INTEGER, PARAMETER, PUBLIC :: T_LBC_FILE = 1 ! initialize tracer from 
                                               ! boundary condition file (lbf)
                                               ! or by MMD
  INTEGER, PARAMETER, PUBLIC :: T_LBC_CST  = 2 ! constant lateral bd. conditions
  INTEGER, PARAMETER, PUBLIC :: T_LBC_ZEROGRAD = 3 ! zero-gradient 
                                                   ! lateral boundary conditions
  INTEGER, PARAMETER, PUBLIC :: T_LBC_USER = 4 ! lateral boundary conditions
                                               ! defined by user  
  INTEGER, PARAMETER         :: T_LBC_MAX  = 4 ! maximum parameter value

  ! BOTTOM BOUNDARY CONDITIONS
  INTEGER, PARAMETER, PUBLIC :: T_BBC_ZEROFLUX = 0 ! zero flux at the bottom
  INTEGER, PARAMETER, PUBLIC :: T_BBC_ZEROVAL  = 1 ! zero value at the bottom
  INTEGER, PARAMETER, PUBLIC :: T_BBC_SURF_VAL = 2 ! bottom value provided as
                                                   ! surface field

  ! RELAXATION TYPES
  INTEGER, PARAMETER, PUBLIC :: T_RELAX_OFF    = 0 ! no relaxation
  INTEGER, PARAMETER, PUBLIC :: T_RELAX_FULL   = 1 ! full relaxation 
                                                   ! (at the 4 boundaries)
  INTEGER, PARAMETER, PUBLIC :: T_RELAX_INFLOW = 2 ! relaxation at inflow 
                                                   ! boundaries only
  INTEGER, PARAMETER         :: T_RELAX_MAX    = 2 ! maximum parameter value
 
  ! ADVECTION TYPES
  INTEGER, PARAMETER, PUBLIC :: T_ADV_OFF      = 0 ! no advection
  INTEGER, PARAMETER, PUBLIC :: T_ADV_ON       = 1 ! Default advection scheme 
                                                   ! (at the 4 boundaries)
  ! special for COSMO  (moisture variables)
  INTEGER, PARAMETER, PUBLIC :: T_ADV_2_LF     = 2 ! 
  INTEGER, PARAMETER, PUBLIC :: T_ADV_3_LF     = 3 ! semi-lagrange 
                                                   ! boundaries only
  INTEGER, PARAMETER         :: T_ADV_MAX      = 3 ! maximum parameter value

  ! um_ak_20120716-
  ! um_ak_20130514+
  ! TURBULENCE (VDIFF) PARAMETERISATION (COSMO)
  INTEGER, PARAMETER, PUBLIC :: T_TURB_OFF   = 0 ! no turbulent mixing
  INTEGER, PARAMETER, PUBLIC :: T_TURB_1D    = 1 ! 1-d (vert.) turbulent mixing
  INTEGER, PARAMETER, PUBLIC :: T_TURB_3D    = 2 ! 3-d (vert.) turbulent mixing
  INTEGER, PARAMETER         :: T_TURB_MAX   = 2 ! maximum parameter value
  ! DAMPING  (COSMO)
  INTEGER, PARAMETER, PUBLIC :: T_DAMP_OFF     = 0
  INTEGER, PARAMETER, PUBLIC :: T_DAMP_ON      = 1
  INTEGER, PARAMETER, PUBLIC :: T_DAMP_FORCED  = 2
  INTEGER, PARAMETER, PUBLIC :: T_DAMP_MAX     = 2
  ! um_ak_20130514-
  ! ub_ak_20170703+
  ! SEDIMENTATION (COSMO)
  INTEGER, PARAMETER, PUBLIC :: T_SEDIM_OFF = 0 ! no sedimentation
  INTEGER, PARAMETER, PUBLIC :: T_SEDIM_ON  = 1 ! apply sedimentation
  INTEGER, PARAMETER, PUBLIC :: T_SEDIM_MAX = 1 ! max. allowable parameter value
  ! PERTURBED PHYSICS
  INTEGER, PARAMETER, PUBLIC :: T_SPPTPERT_OFF  = 0 ! no sppt perturbation
  INTEGER, PARAMETER, PUBLIC :: T_SPPTPERT_ON   = 1 ! sppt perturbation
  INTEGER, PARAMETER, PUBLIC :: T_SPPTPERT_MAX  = 1 ! maximum allowable parameter value
  INTEGER, PARAMETER, PUBLIC :: T_BLOCKMEM_OFF = 0 ! no extra block memory required for tracer
  INTEGER, PARAMETER, PUBLIC :: T_BLOCKMEM_ON  = 1 ! block memory required
  INTEGER, PARAMETER, PUBLIC :: T_BLOCKMEM_MAX = 1 ! maximum allowable parameter value
  
  ! ub_ak_20170703-

  ! op_pj_20160824+
  INTEGER, PARAMETER :: NAMES_CASK_STRLEN = 15
  ! op_pj_20160824-

  ! =================================================================
  ! mz_rs_20160520+
#include "messy_main_tracer_chemprop_mafor.inc"
  ! mz_rs_20160520-
  ! =================================================================

  ! INTEGER CONTAINERS
  INTEGER, PARAMETER, PUBLIC :: I_ADVECT     =  MAX_CASK_I_CHEMPROP +  1  ! ADVECTION
  INTEGER, PARAMETER, PUBLIC :: I_CONVECT    =  MAX_CASK_I_CHEMPROP +  2  ! CONVECTION
  INTEGER, PARAMETER, PUBLIC :: I_VDIFF      =  MAX_CASK_I_CHEMPROP +  3  ! VERTICAL DIFFUSION
  INTEGER, PARAMETER, PUBLIC :: I_WETDEP     =  MAX_CASK_I_CHEMPROP +  4  ! WET DEPOSITION
  INTEGER, PARAMETER, PUBLIC :: I_DRYDEP     =  MAX_CASK_I_CHEMPROP +  5  ! DRY DEPOSITION
  INTEGER, PARAMETER, PUBLIC :: I_SEDI       =  MAX_CASK_I_CHEMPROP +  6  ! SEDIMENTATION
  INTEGER, PARAMETER, PUBLIC :: I_SCAV       =  MAX_CASK_I_CHEMPROP +  7  ! SCAVENGING
  INTEGER, PARAMETER, PUBLIC :: I_MIX        =  MAX_CASK_I_CHEMPROP +  8  ! TURBULENT MIXING
  INTEGER, PARAMETER, PUBLIC :: I_FORCE_COL  =  MAX_CASK_I_CHEMPROP +  9  ! FORCING IN COLUMN MODE
  INTEGER, PARAMETER, PUBLIC :: I_INTEGRATE  =  MAX_CASK_I_CHEMPROP + 10  ! TIME INTEGRATION
  INTEGER, PARAMETER, PUBLIC :: I_TIMEFILTER =  MAX_CASK_I_CHEMPROP + 11  ! TIME FILTER
  INTEGER, PARAMETER, PUBLIC :: I_FORCE_INIT =  MAX_CASK_I_CHEMPROP + 12  ! FORCE INIT AFTER RESTART
  INTEGER, PARAMETER, PUBLIC :: I_AEROSOL_METHOD = MAX_CASK_I_CHEMPROP + 13 ! MODAL OR BIN
  INTEGER, PARAMETER, PUBLIC :: I_AEROSOL_MODE   = MAX_CASK_I_CHEMPROP + 14 ! MODE OR BIN NUMBER
  INTEGER, PARAMETER, PUBLIC :: I_AEROSOL_SOL    = MAX_CASK_I_CHEMPROP + 15 ! SOLUBLE ON/OFF
  INTEGER, PARAMETER, PUBLIC :: I_AEROSOL_HETICE = MAX_CASK_I_CHEMPROP + 16 ! HIGHER ICE SCAV. FRAC.
  INTEGER, PARAMETER, PUBLIC :: I_HDIFF          = MAX_CASK_I_CHEMPROP + 17 ! HORIZONTAL DIFFUSION
  INTEGER, PARAMETER, PUBLIC :: I_RELAX          = MAX_CASK_I_CHEMPROP + 18 ! BOUNDARY DATA AVAILABLE
                                                                            ! i.e. relaxation possible
  INTEGER, PARAMETER, PUBLIC :: I_MMD_INIT       = MAX_CASK_I_CHEMPROP + 19
! op_pj_20100319+
  INTEGER, PARAMETER, PUBLIC :: I_TAG_REG_IDT    = MAX_CASK_I_CHEMPROP + 20 ! id of associated regular
                        !                                                   ! species
  INTEGER, PARAMETER, PUBLIC :: I_TAG_SPECIFIC   = MAX_CASK_I_CHEMPROP + 21 ! flag for special
  !                                                                         ! treatment
! op_pj_20100319-
  ! um_ak_20120716+
  INTEGER, PARAMETER, PUBLIC :: I_INITIAL        = MAX_CASK_I_CHEMPROP + 22 ! initialisation of trac
                                                                            ! treatment 
  INTEGER, PARAMETER, PUBLIC :: I_LBC            = MAX_CASK_I_CHEMPROP + 23 ! lateral boundary cond.
                                                                            ! treatment
  INTEGER, PARAMETER, PUBLIC :: I_DAMP           = MAX_CASK_I_CHEMPROP + 24 ! require damping?
  INTEGER, PARAMETER, PUBLIC :: I_GRIBTAB        = MAX_CASK_I_CHEMPROP + 25 ! grib table number
  INTEGER, PARAMETER, PUBLIC :: I_GRIBPARAM      = MAX_CASK_I_CHEMPROP + 26 ! grib parameter number
  ! um_ak_20120716-
  !
  ! mz_bs_20150702+
  INTEGER, PARAMETER, PUBLIC :: I_MTSKIP         = MAX_CASK_I_CHEMPROP + 27 ! handled by MTSKIP?
  ! mz_bs_20150702-
  ! ub_ak_20170215+
  INTEGER, PARAMETER, PUBLIC :: I_BBC            = MAX_CASK_I_CHEMPROP + 28 ! bottom boundary condition
  INTEGER, PARAMETER, PUBLIC :: I_SPPTPERT       = MAX_CASK_I_CHEMPROP + 29 ! PERTURBED PHYSICS
  INTEGER, PARAMETER, PUBLIC :: I_BLOCKMEM       = MAX_CASK_I_CHEMPROP + 30 ! BLOCKED MEMORY REQUIRED? 
  ! ub_ak_20170215-

  INTEGER, PARAMETER, PUBLIC :: MAX_CASK_I_PROCESS = 30
  !
  CHARACTER(LEN=NAMES_CASK_STRLEN), DIMENSION(MAX_CASK_I_PROCESS) &
       , PARAMETER, PUBLIC :: &
       NAMES_CASK_I_PROCESS = (/ &
       'advect         ', 'convect        ', 'vdiff          ' , &
       'wetdep         ', 'drydep         ', 'sedi           ' , &
       'scav           ', 'mix            ', 'force_col      ' , &
       'integrate      ', 'timefilter     ', 'force_init     ' , &
       'aerosol_method ', 'aerosol_mode   ', 'aerosol_sol    ' , &
       'aerosol_hetice ', 'hori_diff      ', 'relaxation     ' , &
       'mmd_init       ' , &
       'tag_reg_idt    ', 'tag_specific   '                    , &
       'initial type   ', 'lateral_bounds ', 'damping        ' , &
       'grib table num ', 'grib param num ',                     &
       'mtskip         ', 'bot_bound_cond ', 'SPPT           ' , &
       'blockphy       '/) ! ub_ak_20170703
  ! NOTES: - CASK_I(I_TAG_REG_IDT) will default to TRACER-ID
  !          (see subroutine new_tracer)
  INTEGER, DIMENSION(MAX_CASK_I_PROCESS), PARAMETER, PUBLIC :: &
       DEFAULT_CASK_I_PROCESS = &
       (/ ON, ON, ON, OFF, OFF, OFF, OFF, ON, OFF, ON, ON, OFF &
       , MODAL, 0, ON, OFF, OFF, T_RELAX_FULL, OFF, 0, 0 &
!       , T_LBC_FILE, T_INI_FILE, ON , 0 , 0, 0 /)
       , T_INI_FILE, T_LBC_FILE, ON , -1 , -1 &  ! um_ak_20130625
       , OFF, T_BBC_ZEROFLUX, T_SPPTPERT_OFF, T_BLOCKMEM_OFF /) ! ub_ak_20170705

  ! =================================================================

  ! STRING CONTAINERS
  INTEGER, PARAMETER, PUBLIC :: S_AEROSOL_MODEL    = MAX_CASK_S_CHEMPROP + 1
  INTEGER, PARAMETER, PUBLIC :: MAX_CASK_S_PROCESS = 1
  !
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(MAX_CASK_S_PROCESS), &
       PARAMETER, PUBLIC :: DEFAULT_CASK_S_PROCESS = (/ '' /)
  !
  CHARACTER(LEN=NAMES_CASK_STRLEN), DIMENSION(MAX_CASK_S_PROCESS) &
       , PARAMETER, PUBLIC :: &
       NAMES_CASK_S_PROCESS = (/ &
       'aerosol_model  ' &
       /)  

  ! =================================================================

  ! REAL CONTAINERS
  INTEGER, PARAMETER, PUBLIC :: R_VINI             = MAX_CASK_R_CHEMPROP + 1
  INTEGER, PARAMETER, PUBLIC :: MAX_CASK_R_PROCESS = 1
  !
  REAL(DP), DIMENSION(MAX_CASK_R_PROCESS), PARAMETER, PUBLIC :: &
       DEFAULT_CASK_R_PROCESS = (/ 0.0_dp /)
  !
  CHARACTER(LEN=NAMES_CASK_STRLEN), DIMENSION(MAX_CASK_R_PROCESS) &
       , PARAMETER, PUBLIC :: &
       NAMES_CASK_R_PROCESS = (/ &
       'vini           ' /)  

  ! mz_rs_20160520+
  ! combine PROCESS and CHEMPROP:
  ! 1) dimension:
  INTEGER, PUBLIC, PARAMETER :: &
       MAX_CASK_I = MAX_CASK_I_CHEMPROP + MAX_CASK_I_PROCESS
  INTEGER, PUBLIC, PARAMETER :: &
       MAX_CASK_R = MAX_CASK_R_CHEMPROP + MAX_CASK_R_PROCESS
  INTEGER, PUBLIC, PARAMETER :: &
       MAX_CASK_S = MAX_CASK_S_CHEMPROP + MAX_CASK_S_PROCESS

  ! 2) default values:
  INTEGER, DIMENSION(MAX_CASK_I_CHEMPROP), &
    PARAMETER, PUBLIC :: DEFAULT_CASK_I_CHEMPROP = -999
  REAL(DP), DIMENSION(MAX_CASK_R_CHEMPROP), &
    PARAMETER, PUBLIC :: DEFAULT_CASK_R_CHEMPROP = -999.999_DP
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(MAX_CASK_S_CHEMPROP), &
    PARAMETER, PUBLIC :: DEFAULT_CASK_S_CHEMPROP = '---'

  INTEGER, DIMENSION(MAX_CASK_I), PARAMETER, PUBLIC :: &
    DEFAULT_CASK_I = (/ DEFAULT_CASK_I_CHEMPROP, DEFAULT_CASK_I_PROCESS /)
  REAL(DP), DIMENSION(MAX_CASK_R), PARAMETER, PUBLIC :: &
    DEFAULT_CASK_R = (/ DEFAULT_CASK_R_CHEMPROP, DEFAULT_CASK_R_PROCESS /)
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(MAX_CASK_S), PARAMETER, PUBLIC :: &
    DEFAULT_CASK_S = (/ DEFAULT_CASK_S_CHEMPROP, DEFAULT_CASK_S_PROCESS /)

  ! 3) names:
  CHARACTER(LEN=NAMES_CASK_STRLEN), DIMENSION(MAX_CASK_I), &
    PARAMETER, PUBLIC :: NAMES_CASK_I = &
    (/ NAMES_CASK_I_CHEMPROP, NAMES_CASK_I_PROCESS /)
  CHARACTER(LEN=NAMES_CASK_STRLEN), DIMENSION(MAX_CASK_R), &
    PARAMETER, PUBLIC :: NAMES_CASK_R = &
    (/ NAMES_CASK_R_CHEMPROP, NAMES_CASK_R_PROCESS /)
  CHARACTER(LEN=NAMES_CASK_STRLEN), DIMENSION(MAX_CASK_S), &
    PARAMETER, PUBLIC :: NAMES_CASK_S = &
    (/ NAMES_CASK_S_CHEMPROP, NAMES_CASK_S_PROCESS /)
  ! mz_rs_20160520-

  ! =================================================================

  ! op_pj_20150811+
  ! overwrite tracer properties (cask contents) via CTRL-namelist
  TYPE T_TRACPROP_IO
     CHARACTER(LEN=10*STRLEN_TRSET+10) :: trset = ''
     CHARACTER(LEN=10*STRLEN_XLONG+10) :: trlist = ''
     ! LEN must be maximum of NAMES_CASK_[I,S,R]_STRLEN
     CHARACTER(LEN=15)              :: caskname = ''
     CHARACTER(LEN=STRLEN_MEDIUM)   :: cont = ''
  END TYPE T_TRACPROP_IO
  INTEGER, PARAMETER :: NMAXTRACPROP = 100
  TYPE(T_TRACPROP_IO), DIMENSION(NMAXTRACPROP), SAVE :: TPROP
  PUBLIC :: T_TRACPROP_IO, NMAXTRACPROP, TPROP ! only for MPI broadcast on BMIL
  ! op_pj_20150811-

  ! =================================================================

  ! SPECIAL ERROR NUMBERS
  INTEGER, PARAMETER, PUBLIC :: TR_EXIST  = 202  ! TRACER EXISTS
  INTEGER, PARAMETER, PUBLIC :: TR_NEXIST = 205  ! TRACER DOES NOT EXIST

  ! TRACER BLOCK SIZE required for COSMO halo exchange
  INTEGER, PUBLIC  :: n_trcr_block = 0  ! um_ak_20150316

  ! PUBLIC STRUCTURES
  PUBLIC :: t_ident         ! TRACER IDENTIFICATION
  PUBLIC :: t_meta          ! ADDITIONAL META INFORMATION
  !
  PUBLIC :: t_trinfo        ! META-STRUCT WITH ALL INFORMATION
  PUBLIC :: t_trinfo_tp     ! TRACER PROPERTIES
  PUBLIC :: t_trinfo_list   ! LIST OF META-STRUCT WITH ALL INFORMATION

  ! PUBLIC SUBROUTINES
  PUBLIC :: new_tracer_set   ! BML, BMIL     ! DEFINE NEW TRACER SET
  PUBLIC :: copy_tracer_set  ! BML, BMIL     ! COPY COMPLETE TRACER SET
  PUBLIC :: setup_tracer_set ! BML, BMIL     ! ALLOCATE MEMORY FOR TRACER SET
  PUBLIC :: get_tracer_set   ! BML, BMIL     ! SET REFERENCES TO TRACER SETS
  PUBLIC :: print_tracer_set ! BML, BMIL     ! PRINT TRACER SET SUMMARY
  PUBLIC :: print_tracer_set_val ! BML, BMIL ! PRINT TRACER VALUE RANGE
  PUBLIC :: clean_tracer_set ! BML, BMIL     ! REMOVE TRACER SET FROM MEMORY
  ! um_ak_20130527+
  PUBLIC :: get_tracer_set_info
  ! um_ak_20130527+
  !
  PUBLIC :: new_tracer       ! SMIL          ! DEFINE NEW TRACER IN SET
  PUBLIC :: set_tracer       ! SMIL          ! DEFINE TRACER PROPERTIES
! mz_pj_20070817+ will become obsolete
  PUBLIC :: new_tracer_old   ! SMIL          ! DEFINE NEW TRACER IN SET
! mz_pj_20070817-
! mz_pj_20090507+
  INTERFACE get_tracer
     MODULE PROCEDURE get_tracer_by_name
     MODULE PROCEDURE get_tracer_by_id
     ! op_pj_20100325+
     MODULE PROCEDURE get_tracer_i
     MODULE PROCEDURE get_tracer_r
     MODULE PROCEDURE get_tracer_s
     ! op_pj_20100325-
  END INTERFACE
! mz_pj_20090507-
  PUBLIC :: get_tracer       ! SMIL          ! GET TRACER INFORMATION FROM SET
  PUBLIC :: get_tracer_list  ! SMIL          ! GET TRACERS WITH SAME BASENAME
  PUBLIC :: tracer_iniflag                   ! GET/SET INITIALISATION FLAG
  !
  PUBLIC :: tracer_error_str                 ! RETURN STATUS INFORMATION
  PUBLIC :: param2string                     ! PARAMETER TO STRING CONVERSION
  PUBLIC :: full2base_sub                    ! fullname -> basename + subname
  PUBLIC :: main_tracer_read_nml_ctrl
  PUBLIC :: set_tracer_properties            ! op_pj_20150811
  PUBLIC :: get_tracer_set_id                ! get tracer set id

  INTERFACE get_tracer_set
     MODULE PROCEDURE get_tracer_set_by_name
     MODULE PROCEDURE get_tracer_set_by_id
  END INTERFACE

  INTERFACE set_tracer
     MODULE PROCEDURE set_tracer_i
     MODULE PROCEDURE set_tracer_r
     MODULE PROCEDURE set_tracer_d
     MODULE PROCEDURE set_tracer_s
     MODULE PROCEDURE set_tracer_refspec
  END INTERFACE

  ! mz_rs_20160521+
  PUBLIC :: get_chemprop, get_chemprop_index
  INTERFACE get_chemprop
     MODULE PROCEDURE get_chemprop_real
     MODULE PROCEDURE get_chemprop_int
     MODULE PROCEDURE get_chemprop_char
  END INTERFACE
  ! mz_rs_20160521-

  ! PRIVATE SUBROUTINES

  ! STRUCTURE DEFINITIONS
  TYPE t_ident    ! IDENTIFICATION
     CHARACTER(LEN=STRLEN_MEDIUM)     :: basename    = '' ! name of tracer
     CHARACTER(LEN=STRLEN_MEDIUM)     :: subname     = '' ! OPTIONAL subname
     CHARACTER(LEN=STRLEN_FNAME)      :: fullname    = '' ! name_subname
!!$  CHARACTER(LEN=STRLEN_VLONG)      :: longname    = '' ! ! op_pj_20170508
     CHARACTER(LEN=STRLEN_ULONG)      :: longname    = '' ! ! op_pj_20170508
     CHARACTER(LEN=STRLEN_MEDIUM)     :: unit        = '' !
     CHARACTER(LEN=STRLEN_MEDIUM)     :: submodel    = '' ! requesting submodel
     INTEGER                          :: idx              ! tracer index in set
     INTEGER                          :: idxblck     = 0  ! tracer index in blck set ! ub_ak_20170727
     INTEGER                          :: medium   = AIR   ! hosting medium
     INTEGER                          :: quantity = AMOUNTFRACTION
     INTEGER                          :: type     = SINGLE ! 
  END TYPE t_ident

  TYPE t_meta
     INTEGER,                      DIMENSION(MAX_CASK_I)  :: cask_i &
          = DEFAULT_CASK_I
     REAL(dp),                     DIMENSION(MAX_CASK_R)  :: cask_r &
          = DEFAULT_CASK_R
     CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(MAX_CASK_S)  :: cask_s &
          = DEFAULT_CASK_S
     LOGICAL :: linit = .FALSE.
  END TYPE t_meta

  TYPE t_trinfo
     TYPE(t_ident) :: ident       ! IDENTIFICATION
     TYPE(t_meta)  :: meta        ! ADDITIONAL META-INFORMATION
  END TYPE t_trinfo

  TYPE t_trinfo_list
     TYPE(t_trinfo)               :: info
     TYPE(t_trinfo_list), POINTER :: next
  END TYPE t_trinfo_list

  TYPE t_trinfo_tp
     TYPE(t_trinfo), POINTER :: tp
  END TYPE t_trinfo_tp

  ! ################## INTERNAL MEMORY AND POINTER MANAGEMENT ############
  TYPE T_TRACERSET
     ! NAME OF THE TRACER-SET FOR IDENTIFICATION
     CHARACTER(LEN=STRLEN_TRSET)              :: name
     !
     ! SWITCH TO ENABLE/DISABLE THIS TRACER SET, E.G., IF ITS
     ! EXISTENCE IS DEPENDENT ON A SPECIFIC SUBMODEL
     LOGICAL                                  :: l_enable = .TRUE.
     !
     ! SWITCH FOR TRIGGERING THIS SET TO BE INITIALIZED BY THE
     ! main_tracer_init_tracer (BMIL) SUBROUTINE
     LOGICAL                                  :: l_init = .TRUE.
     !
     ! NUMBER OF TRACERS IN SET
     INTEGER                                  :: ntrac  = 0
     ! NUMBER OF BLOCKED TRACERS IN SET (required for COSMO) ! ub_ak_20170727
     INTEGER                                  :: nblck  = 0  ! ub_ak_20170727
     !
     ! CONCATENATED LIST OF TRACER META INFORAMTION STRUCTURES
     TYPE(t_trinfo_list),             POINTER :: tilist => NULL()
     !
     ! ENUMERATED LIST OF TRACER META INFORAMTION STRUCTURES
     ! NOTE: THIS IS ONLY AVAILABLE AFTER CALL setup_tracer_set
     TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti     => NULL()
     !
     ! MEMORY FOR TRACER DATA
     REAL(DP), DIMENSION(:,:,:,:,:,:),  POINTER :: mem    => NULL()
     ! NUMBER OF (TIME-) LEVELS
     INTEGER                                    :: nt     = 0
     ! 3 STANDARD (TIME-) LEVELS OF MEMORY FOR TRACER DATA
     LOGICAL                                  :: l_tfstd = .FALSE.
     REAL(DP), DIMENSION(:,:,:,:,:),  POINTER :: xt      => NULL()
     REAL(DP), DIMENSION(:,:,:,:,:),  POINTER :: xtte    => NULL()
     REAL(DP), DIMENSION(:,:,:,:,:),  POINTER :: xtm1    => NULL()
     ! 'EXTENDED' MEMORY
     INTEGER                                  :: nex     = 0
     REAL(DP), DIMENSION(:,:,:,:,:),  POINTER :: xmem    => NULL()
     ! um_ak_20130521+
     ! 3d (4D)-TRACER FIELD INCLUDING TIME LEVELS
     REAL(DP), DIMENSION(:,:,:,:,:),  POINTER :: xt_tl   => NULL()
     ! um_ak_20130521-
     ! ub_ak_20170705+
     REAL(DP), DIMENSION(:,:,:,:),    POINTER :: xtblck   => NULL()
     REAL(DP), DIMENSION(:,:,:),      POINTER :: xtteblck => NULL()
     ! ub_ak_20170705-

  END TYPE T_TRACERSET
  PUBLIC :: T_TRACERSET
  !
  INTEGER, PARAMETER, PUBLIC :: NMAXSETID = 10  ! MAX. NUMBER OF TRACER SETS
  INTEGER,      SAVE, PUBLIC :: NSETID    = 0   ! ACT. NUMBER OF TRACER SETS
  TYPE(T_TRACERSET), DIMENSION(NMAXSETID), PUBLIC, SAVE :: TRSET
  ! ################## INTERNAL MEMORY AND POINTER MANAGEMENT ############

CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE new_tracer_set(status, setname, l_enable)

    IMPLICIT NONE
    INTRINSIC :: LEN, TRIM

    ! I/O
    INTEGER,               INTENT(OUT) :: status
    CHARACTER(LEN=*),      INTENT(IN)  :: setname
    LOGICAL,               INTENT(IN)  :: l_enable

    ! LOCAL
    INTEGER :: zid

    status = 100 ! ERROR: setname too long
    IF (LEN(setname) > STRLEN_TRSET) RETURN

    status = 101 ! ERROR: setname not unique
    DO zid = 1, NSETID
       IF (TRIM(trset(zid)%name) == TRIM(setname)) RETURN
    END DO

    status = 102 ! ERROR: no free tracer set available
    NSETID = NSETID + 1
    IF (NSETID > NMAXSETID) RETURN

    ! RESULT
    trset(NSETID)%name     = TRIM(setname)
    trset(NSETID)%l_enable = l_enable

    status = 0

  END SUBROUTINE new_tracer_set
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE copy_tracer_set(status, oldset, newset)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, LEN, TRIM

    ! I/O
    INTEGER,             INTENT(OUT) :: status
    CHARACTER(LEN=*),    INTENT(IN)  :: oldset
    CHARACTER(LEN=*),    INTENT(IN)  :: newset

    ! LOCAL
    INTEGER                      :: id_old, id_new
    TYPE(t_trinfo_list), POINTER :: ti_old => NULL()
    TYPE(t_trinfo_list), POINTER :: ti_new => NULL()
    TYPE(t_trinfo_list), POINTER :: te_new => NULL()

    status = 100 ! ERROR: setname too long
    IF (LEN(oldset) > STRLEN_TRSET) RETURN
    IF (LEN(newset) > STRLEN_TRSET) RETURN

    status = 101 ! ERROR: setname not unique
    DO id_new = 1, NSETID
       IF (TRIM(trset(id_new)%name) == TRIM(newset)) RETURN
    END DO

    status = 102 ! ERROR: no free tracer set available
    NSETID = NSETID + 1
    IF (NSETID > NMAXSETID) RETURN
    id_new = NSETID

    ! CHECK, IF OLD SET IS AVAILABLE; AND SET ID
    CALL get_tracer_set_id(status, oldset, id_old)
    IF (status /= 0) RETURN

    ! RESULT: NEW SET
    trset(id_new)%name = TRIM(newset)
    trset(id_new)%l_enable = trset(id_old)%l_enable

    ! COPY ALL TRACERS; LOOP OVER TRACERS IN OLD SET
    ti_old => trset(id_old)%tilist
    tracer_loop: DO 
       IF (.NOT. ASSOCIATED(ti_old)) EXIT

       ! ADD NEW TRACER TO NEW LIST
       ALLOCATE(ti_new)
       NULLIFY(ti_new%next)
       IF (trset(id_new)%ntrac == 0) THEN
          trset(id_new)%tilist => ti_new   ! SET POINTER TO FIRST ELEMENT
       ELSE
          te_new%next => ti_new            ! SET NEXT POINTER OF LAST ELEMENT
          !                                ! TO NEW ELEMENT
       END IF
       te_new => ti_new

       ! ADJUST TRACER INDEX
       trset(id_new)%ntrac = trset(id_new)%ntrac + 1

       ! SET TRACER INFORMATION
       ti_new%info = ti_old%info

       ti_old => ti_old%next
    END DO tracer_loop

    status = 0

  END SUBROUTINE copy_tracer_set
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE new_tracer(status, setname, basename, submodel             &
       , idx, subname, longname, unit, medium, quantity, type           &
       , cask_i, cask_r, cask_s, refspec) ! op_pj_20170314: refspec

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, INDEX, LEN, PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    !
    CHARACTER(LEN=*), INTENT(IN)            :: basename
    CHARACTER(LEN=*), INTENT(IN)            :: submodel
    !
    INTEGER,          INTENT(OUT), OPTIONAL :: idx
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: subname
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: longname
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: unit
    INTEGER,          INTENT(IN),  OPTIONAL :: medium
    INTEGER,          INTENT(IN),  OPTIONAL :: quantity
    INTEGER,          INTENT(IN),  OPTIONAL :: type
    !
    INTEGER,  DIMENSION(MAX_CASK_I),   INTENT(IN),  OPTIONAL :: cask_i
    REAL(DP), DIMENSION(MAX_CASK_R),   INTENT(IN),  OPTIONAL :: cask_r
    CHARACTER(LEN=STRLEN_MEDIUM), &
         DIMENSION(MAX_CASK_S), INTENT(IN), OPTIONAL :: cask_s
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL  :: refspec ! op_pj_20170314

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER   :: substr = 'new_tracer'
    INTEGER                        :: zid     = 0
    TYPE(t_trinfo_list), POINTER   :: ti      => NULL()
    TYPE(t_trinfo_list), POINTER   :: te      => NULL()
    CHARACTER(LEN=STRLEN_FNAME)    :: fullname
    INTEGER                        :: i
    CHARACTER(LEN=STRLEN_FNAME)    :: specname ! op_pj_20170314

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! TRACER
    status = 200 ! ERROR: tracer basename too long
    IF (LEN(basename) > STRLEN_MEDIUM) RETURN

    status = 208 ! ERROR: basename MUST NOT contain '_'
    IF (INDEX(basename, '_') /= 0) RETURN

    IF (PRESENT(subname)) THEN
       status = 201 ! ERROR: tracer subname too long
       IF (LEN(subname) > STRLEN_MEDIUM) RETURN
       IF (TRIM(subname) /= '') THEN
          fullname = TRIM(basename)//'_'//TRIM(subname)
       ELSE
          fullname = TRIM(basename)
       END IF
    ELSE
       fullname = TRIM(basename)
    END IF

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       IF (PRESENT(idx)) idx = 0
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (TRIM(ti%info%ident%fullname) == TRIM(fullname)) THEN
          IF (PRESENT(idx)) idx = ti%info%ident%idx ! op_pj_20171017
          status = TR_EXIST  ! ERROR: tracer exists already
          RETURN
       END IF
       te => ti
       ti => ti%next
    END DO

    ! ADD NEW TRACER TO LIST
    ALLOCATE(ti)
    NULLIFY(ti%next)
    IF (trset(zid)%ntrac == 0) THEN
       trset(zid)%tilist => ti         ! SET POINTER TO FIRST ELEMENT
    ELSE
       te%next => ti                   ! SET NEXT POINTER OF LAST ELEMENT
       !                               ! TO NEW ELEMENT
    END IF

    ! ADJUST TRACER INDEX
    trset(zid)%ntrac = trset(zid)%ntrac + 1
    IF (PRESENT(idx)) idx = trset(zid)%ntrac

    ! SET TRACER INFORMATION
    ! - IDENTIFICATION
    ! -- MANDATORY
    ti%info%ident%basename                   = TRIM(basename)    ! BASENAME
    !
    IF (LEN(submodel) > STRLEN_MEDIUM) THEN
       status = 204 ! ERROR: unit too long
       RETURN
    END IF
    ti%info%ident%submodel                   = TRIM(submodel)    ! SUBMODEL
    !
    ti%info%ident%fullname                   = TRIM(fullname)    ! FULLNAME
    ti%info%ident%idx                        = trset(zid)%ntrac  ! SET INDEX

    ! -- OPTIONAL
    !
    IF (PRESENT(subname)) THEN           ! SUBNAME
       IF (LEN(subname) > STRLEN_MEDIUM) THEN
          status = 201 ! ERROR: subname too long
          RETURN
       END IF
       ti%info%ident%subname    = TRIM(subname)
    END IF
    !
    !                                    ! LONGNAME
    IF (PRESENT(longname)) THEN
       IF (LEN(longname) > STRLEN_VLONG) THEN
          status = 207 ! ERROR: longname too long
          RETURN
       END IF
       ti%info%ident%longname  = TRIM(longname)
    END IF
    !
    IF (PRESENT(unit)) THEN              ! UNIT
       IF (LEN(unit) > STRLEN_MEDIUM) THEN
          status = 203 ! ERROR: unit too long
          RETURN
       END IF
       ti%info%ident%unit          = TRIM(unit)
    END IF
    !
    IF (PRESENT(medium)) THEN            ! MEDIUM
       IF ((medium < 0) .OR. (medium > MAX_MEDIUM) ) THEN
          status = 300  ! ERROR: unknown medium
          RETURN
       END IF
       ti%info%ident%medium   = medium
    END IF
    !
    IF (PRESENT(quantity)) THEN          ! QUANTITY
       IF ((quantity < 0) .OR. (quantity > MAX_QUANTITY) ) THEN
          status = 301  ! ERROR: unknown quantity
          RETURN
       END IF
       ti%info%ident%quantity = quantity
    END IF
    ! 
    IF (PRESENT(type)) THEN              ! TYPE
       IF ((type < 0) .OR. (type > MAX_TYPE  ) ) THEN
          status = 302  ! ERROR: unknown type
          RETURN
       END IF
       ti%info%ident%type = type
    END IF

    ! ASSIGN OPTIONAL PARAMETERS
    ! DEFAULT
    ti%info%meta%cask_i(:) = default_cask_i(:)
    ti%info%meta%cask_r(:) = default_cask_r(:)
    ti%info%meta%cask_s(:) = default_cask_s(:)
    ! mz_rs_20160521+
    ! op_pj_20170314+
    IF (PRESENT(refspec)) THEN
       specname = TRIM(ADJUSTL(refspec))
    ELSE
       specname = TRIM(ADJUSTL(basename))
    END IF
    ! op_pj_20170314-
    ! check if basename or refspec can be found in chemprop, then add values:
    DO i = 1, N_CHEMPROP
!!$    IF (TRIM(chemprop(i)%kppname) == TRIM(basename)) THEN ! op_pj_20170314
       IF (TRIM(chemprop(i)%kppname) == TRIM(specname)) THEN ! op_pj_20170314
          WRITE(*,*) substr,' INFO: using chemical properties of '//&
               &TRIM(chemprop(i)%kppname)//' for tracer '//TRIM(fullname)
          ti%info%meta%cask_i(1:MAX_CASK_I_CHEMPROP) = chemprop(i)%cask_i(:)
          ti%info%meta%cask_r(1:MAX_CASK_R_CHEMPROP) = chemprop(i)%cask_r(:)
          ti%info%meta%cask_s(1:MAX_CASK_S_CHEMPROP) = chemprop(i)%cask_s(:)
          EXIT
       ENDIF
    ENDDO
    ! mz_rs_20160521-
    ! SPECIFIC
    IF (PRESENT(cask_i)) ti%info%meta%cask_i(:) = cask_i(:)
    IF (PRESENT(cask_r)) ti%info%meta%cask_r(:) = cask_r(:)
    IF (PRESENT(cask_s)) ti%info%meta%cask_s(:) = cask_s(:)
    !
    ! special "dynamic default":
    IF (ti%info%meta%cask_i(I_TAG_REG_IDT) == 0) &
         ti%info%meta%cask_i(I_TAG_REG_IDT) = ti%info%ident%idx

    ! ub_ak_20170705+
    ! block memory requested?
    IF (ti%info%meta%cask_i(I_BLOCKMEM) == T_BLOCKMEM_ON) THEN
       trset(zid)%nblck = trset(zid)%nblck + 1
       ti%info%ident%idxblck = trset(zid)%nblck  ! SET BLOCK INDEX
    END IF
    ! ub_ak_20170705-

    status = 0

  END SUBROUTINE new_tracer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_tracer_s(status, setname, idx, flag, s)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    INTEGER,          INTENT(IN)            :: flag
    CHARACTER(LEN=*), INTENT(IN)            :: s

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    IF ( (flag < 1) .OR. (flag > MAX_CASK_S) ) THEN
       status = 702 ! STRING FLAG INDEX OUT OF RANGE
       RETURN
    END IF
    ti%info%meta%cask_s(flag) = s

    status = 0

  END SUBROUTINE set_tracer_s
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_tracer_i(status, setname, idx, flag, i)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    INTEGER,          INTENT(IN)            :: flag
    INTEGER,          INTENT(IN)            :: i

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    IF ( (flag < 1) .OR. (flag > MAX_CASK_I) ) THEN
       status = 700 ! INTEGER FLAG INDEX OUT OF RANGE
       RETURN
    END IF

    ! ub_ak_20170710+
    IF (flag == I_BLOCKMEM .AND. i == T_BLOCKMEM_ON) THEN
       ! add tracer to block memory list
       IF (ti%info%ident%idxblck <= 0) THEN
          trset(zid)%nblck = trset(zid)%nblck + 1 
          ti%info%ident%idxblck = trset(zid)%nblck
       END IF
    END IF
    ! ub_ak_20170710-

    ti%info%meta%cask_i(flag) = i

    status = 0

  END SUBROUTINE set_tracer_i
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_tracer_d(status, setname, idx, flag, r)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    INTEGER,          INTENT(IN)            :: flag
    REAL(DP),         INTENT(IN)            :: r

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    IF ( (flag < 1) .OR. (flag > MAX_CASK_R) ) THEN
       status = 701 ! REAL FLAG INDEX OUT OF RANGE
       RETURN
    END IF
    ti%info%meta%cask_r(flag) = r

    status = 0

  END SUBROUTINE set_tracer_d
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_tracer_r(status, setname, idx, flag, r)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    INTEGER,          INTENT(IN)            :: flag
    REAL    ,         INTENT(IN)            :: r

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    IF ( (flag < 1) .OR. (flag > MAX_CASK_R) ) THEN
       status = 701 ! REAL FLAG INDEX OUT OF RANGE
       RETURN
    END IF
    ti%info%meta%cask_r(flag) = REAL(r,dp)

    status = 0

  END SUBROUTINE set_tracer_r
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE set_tracer_refspec(status, setname, idx, refspec)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    CHARACTER(LEN=*), INTENT(IN)            :: refspec 

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid
    INTEGER                       :: spec_id

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    spec_id = get_chemprop_index(TRIM(refspec))
    IF (spec_id > 0) THEN
       ti%info%meta%cask_i(1:MAX_CASK_I_CHEMPROP) = &
            chemprop(spec_id)%cask_i(:)
       ti%info%meta%cask_r(1:MAX_CASK_R_CHEMPROP) = &
            chemprop(spec_id)%cask_r(:)
       ti%info%meta%cask_s(1:MAX_CASK_S_CHEMPROP) = & 
            chemprop(spec_id)%cask_s(:)
    ELSE
       status = 3001
       RETURN
    ENDIF

    status = 0

  END SUBROUTINE set_tracer_refspec
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE new_tracer_old(status, setname, basename, submodel         &
       , idx, subname, longname, unit, medium, quantity, type           &
       , nadvect, nconvect, nvdiff, nwetdep, ndrydep, nsedi, nscav      &
       , nmix                                                           &
       , nintegrate, ntimefilter                                        &
       , vini, lforce_init                                              &
       , nforce_col                                                     &
       , molarmass, henry, dryreac_sf                                   &
       , m_aerosol_method, m_aerosol_modelname, m_aerosol_density       &
       , m_aerosol_mode, m_aerosol_nsol                                 &
       )

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, INDEX, LEN, PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    !
    CHARACTER(LEN=*), INTENT(IN)            :: basename
    CHARACTER(LEN=*), INTENT(IN)            :: submodel
    !
    INTEGER,          INTENT(OUT), OPTIONAL :: idx
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: subname
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: longname
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: unit
    INTEGER,          INTENT(IN),  OPTIONAL :: medium
    INTEGER,          INTENT(IN),  OPTIONAL :: quantity
    INTEGER,          INTENT(IN),  OPTIONAL :: type
    !
    INTEGER,          INTENT(IN),  OPTIONAL :: nadvect 
    INTEGER,          INTENT(IN),  OPTIONAL :: nconvect
    INTEGER,          INTENT(IN),  OPTIONAL :: nvdiff
    INTEGER,          INTENT(IN),  OPTIONAL :: nwetdep
    INTEGER,          INTENT(IN),  OPTIONAL :: ndrydep
    INTEGER,          INTENT(IN),  OPTIONAL :: nsedi
    INTEGER,          INTENT(IN),  OPTIONAL :: nscav
    INTEGER,          INTENT(IN),  OPTIONAL :: nmix
    ! SPECIAL FOR COLUMN MODE
    INTEGER,          INTENT(IN),  OPTIONAL :: nforce_col
    ! NUMERICAL
    INTEGER,          INTENT(IN),  OPTIONAL :: nintegrate
    INTEGER,          INTENT(IN),  OPTIONAL :: ntimefilter
    ! INFRASTRUCTURE
    REAL(DP),         INTENT(IN),  OPTIONAL :: vini
    LOGICAL,          INTENT(IN),  OPTIONAL :: lforce_init
    ! SINGLE
    REAL(DP),         INTENT(IN),  OPTIONAL :: molarmass
    REAL(DP),         INTENT(IN),  OPTIONAL :: henry
    REAL(DP),         INTENT(IN),  OPTIONAL :: dryreac_sf
    ! MEDIUM AEROSOL
    INTEGER,          INTENT(IN),  OPTIONAL :: m_aerosol_method
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: m_aerosol_modelname
    REAL(DP),         INTENT(IN),  OPTIONAL :: m_aerosol_density
    INTEGER,          INTENT(IN),  OPTIONAL :: m_aerosol_mode
    INTEGER,          INTENT(IN),  OPTIONAL :: m_aerosol_nsol

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'new_tracer_old'
    INTEGER                        :: zid     = 0
    TYPE(t_trinfo_list), POINTER   :: ti      => NULL()
    TYPE(t_trinfo_list), POINTER   :: te      => NULL()
    CHARACTER(LEN=STRLEN_FNAME)    :: fullname
    INTEGER                        :: i

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! TRACER
    status = 200 ! ERROR: tracer basename too long
    IF (LEN(basename) > STRLEN_MEDIUM) RETURN

    status = 208 ! ERROR: basename MUST NOT contain '_'
    IF (INDEX(basename, '_') /= 0) RETURN

    IF (PRESENT(subname)) THEN
       status = 201 ! ERROR: tracer subname too long
       IF (LEN(subname) > STRLEN_MEDIUM) RETURN
       IF (TRIM(subname) /= '') THEN
          fullname = TRIM(basename)//'_'//TRIM(subname)
       ELSE
          fullname = TRIM(basename)
       END IF
    ELSE
       fullname = TRIM(basename)
    END IF

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       IF (PRESENT(idx)) idx = 0
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (TRIM(ti%info%ident%fullname) == TRIM(fullname)) THEN
          status = TR_EXIST  ! ERROR: tracer exists already
          RETURN
       END IF
       te => ti
       ti => ti%next
    END DO

    ! ADD NEW TRACER TO LIST
    ALLOCATE(ti)
    NULLIFY(ti%next)
    IF (trset(zid)%ntrac == 0) THEN
       trset(zid)%tilist => ti         ! SET POINTER TO FIRST ELEMENT
    ELSE
       te%next => ti                   ! SET NEXT POINTER OF LAST ELEMENT
       !                               ! TO NEW ELEMENT
    END IF

    ! ADJUST TRACER INDEX
    trset(zid)%ntrac = trset(zid)%ntrac + 1
    IF (PRESENT(idx)) idx = trset(zid)%ntrac

    ! SET TRACER INFORMATION
    ! - IDENTIFICATION
    ! -- MANDATORY
    ti%info%ident%basename                   = TRIM(basename)    ! BASENAME
    !
    IF (LEN(submodel) > STRLEN_MEDIUM) THEN
       status = 204 ! ERROR: unit too long
       RETURN
    END IF
    ti%info%ident%submodel                   = TRIM(submodel)    ! SUBMODEL
    !
    ti%info%ident%fullname                   = TRIM(fullname)    ! FULLNAME
    ti%info%ident%idx                        = trset(zid)%ntrac  ! SET INDEX

    ! -- OPTIONAL
    !
    IF (PRESENT(subname)) THEN           ! SUBNAME
       IF (LEN(subname) > STRLEN_MEDIUM) THEN
          status = 201 ! ERROR: subname too long
          RETURN
       END IF
       ti%info%ident%subname    = TRIM(subname)
    END IF
    !
    !                                    ! LONGNAME
    IF (PRESENT(longname)) THEN
       IF (LEN(longname) > STRLEN_VLONG) THEN
          status = 207 ! ERROR: longname too long
          RETURN
       END IF
       ti%info%ident%longname  = TRIM(longname)
    END IF
    !
    IF (PRESENT(unit)) THEN              ! UNIT
       IF (LEN(unit) > STRLEN_MEDIUM) THEN
          status = 203 ! ERROR: unit too long
          RETURN
       END IF
       ti%info%ident%unit          = TRIM(unit)
    END IF
    !
    IF (PRESENT(medium)) THEN            ! MEDIUM
       IF ((medium < 0) .OR. (medium > MAX_MEDIUM) ) THEN
          status = 300  ! ERROR: unknown medium
          RETURN
       END IF
       ti%info%ident%medium   = medium
    END IF
    !
    IF (PRESENT(quantity)) THEN          ! QUANTITY
       IF ((quantity < 0) .OR. (quantity > MAX_QUANTITY) ) THEN
          status = 301  ! ERROR: unknown quantity
          RETURN
       END IF
       ti%info%ident%quantity = quantity
    END IF
    ! 
    IF (PRESENT(type)) THEN              ! TYPE
       IF ((type < 0) .OR. (type > MAX_TYPE  ) ) THEN
          status = 302  ! ERROR: unknown type
          RETURN
       END IF
       ti%info%ident%type = type
    END IF

    ! ASSIGN OPTIONAL PARAMETERS
    ! DEFAULT
    ti%info%meta%cask_i(:) = default_cask_i(:)
    ti%info%meta%cask_r(:) = default_cask_r(:)
    ti%info%meta%cask_s(:) = default_cask_s(:)
    ! op_pj_20160825+
    ! check if basename can be found in chemprop, then add values:
    DO i = 1, N_CHEMPROP
       IF (TRIM(chemprop(i)%kppname) == TRIM(basename)) THEN
          WRITE(*,*) substr,' INFO: using chemical properties of '//&
               &TRIM(chemprop(i)%kppname)//' for tracer '//TRIM(fullname)
          ti%info%meta%cask_i(1:MAX_CASK_I_CHEMPROP) = chemprop(i)%cask_i(:)
          ti%info%meta%cask_r(1:MAX_CASK_R_CHEMPROP) = chemprop(i)%cask_r(:)
          ti%info%meta%cask_s(1:MAX_CASK_S_CHEMPROP) = chemprop(i)%cask_s(:)
          EXIT
       ENDIF
    ENDDO
    ! op_pj_20160825-

    ! PHYSICAL PROCESSES
    IF (PRESENT(nadvect))     ti%info%meta%cask_i(I_advect)     = nadvect 
    IF (PRESENT(nconvect))    ti%info%meta%cask_i(I_convect)    = nconvect
    IF (PRESENT(nvdiff))      ti%info%meta%cask_i(I_vdiff)      = nvdiff
    IF (PRESENT(nwetdep))     ti%info%meta%cask_i(I_wetdep)     = nwetdep
    IF (PRESENT(ndrydep))     ti%info%meta%cask_i(I_drydep)     = ndrydep
    IF (PRESENT(nsedi))       ti%info%meta%cask_i(I_sedi)       = nsedi
    IF (PRESENT(nscav))       ti%info%meta%cask_i(I_scav)       = nscav
    IF (PRESENT(nmix))        ti%info%meta%cask_i(I_mix)        = nmix
    ! SPECIAL FOR COLUMN MODE
    IF (PRESENT(nforce_col))  ti%info%meta%cask_i(I_force_col)  = nforce_col
    ! NUMERICAL
    IF (PRESENT(nintegrate))  ti%info%meta%cask_i(I_integrate)  = nintegrate
    IF (PRESENT(ntimefilter)) ti%info%meta%cask_i(I_timefilter) = ntimefilter
    ! INFRASTRUCTURE
    IF (PRESENT(vini))        ti%info%meta%cask_r(R_vini)        = vini
    IF (PRESENT(lforce_init)) THEN
       IF (lforce_init) THEN
          ti%info%meta%cask_i(I_force_init) = ON
       ELSE
          ti%info%meta%cask_i(I_force_init) = OFF
       END IF
    END IF
    ! SINGLE
    IF (PRESENT(molarmass))   ti%info%meta%cask_r(R_molarmass)  = molarmass
    IF (PRESENT(henry))       ti%info%meta%cask_r(R_pss  )      = henry
    IF (PRESENT(dryreac_sf))  ti%info%meta%cask_r(R_dryreac_sf) = dryreac_sf
    ! MEDIUM AEROSOL
    IF (PRESENT(m_aerosol_method))    ti%info%meta%cask_I(I_aerosol_method) = &
         m_aerosol_method
    IF (PRESENT(m_aerosol_modelname)) ti%info%meta%cask_S(S_aerosol_model) = &
         m_aerosol_modelname
    IF (PRESENT(m_aerosol_density))  ti%info%meta%cask_R(R_aerosol_density) = &
         m_aerosol_density
    IF (PRESENT(m_aerosol_mode))      ti%info%meta%cask_I(I_aerosol_mode) = &
         m_aerosol_mode
    IF (PRESENT(m_aerosol_nsol))      ti%info%meta%cask_I(I_aerosol_sol) = &
         m_aerosol_nsol

    ! special "dynamic default":
    IF (ti%info%meta%cask_i(I_TAG_REG_IDT) == 0) &
         ti%info%meta%cask_i(I_TAG_REG_IDT) = ti%info%ident%idx

     ! ub_ak_20170705+
    ! block memory requested?
    IF (ti%info%meta%cask_i(I_BLOCKMEM) == T_BLOCKMEM_ON) THEN
       trset(zid)%nblck = trset(zid)%nblck + 1
       ti%info%ident%idxblck = trset(zid)%nblck  ! SET BLOCK INDEX
    END IF
    ! ub_ak_20170705-
   !
    status = 0

  END SUBROUTINE new_tracer_old
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_set_by_name(status, setname, trlist, ti, ntrac &
!    , xt, xtte, xtm1, xmem, l_tfstd, l_init, l_enable)        ! um_ak_20130521
     , xt, xtte, xtm1, xmem, xt_tl, l_tfstd, l_init, l_enable &! um_ak_20130521
     , ntracblck,xtblck, xtteblck) ! ub_ak_20170705
    IMPLICIT NONE
    INTRINSIC :: LEN, PRESENT, TRIM

    ! I/O
    INTEGER,                         INTENT(OUT)           :: status
    CHARACTER(LEN=*),                INTENT(IN)            :: setname
    TYPE(t_trinfo_list),             POINTER,     OPTIONAL :: trlist
    TYPE(t_trinfo_tp), DIMENSION(:), POINTER,     OPTIONAL :: ti
    INTEGER,                         INTENT(OUT), OPTIONAL :: ntrac
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xt
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xtte
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xtm1
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xmem
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xt_tl !um_ak_20130521
    ! ub_ak_20170705+
    INTEGER,                         INTENT(OUT), OPTIONAL :: ntracblck
    REAL(DP), DIMENSION(:,:,:,:),    POINTER,     OPTIONAL :: xtblck 
    REAL(DP), DIMENSION(:,:,:),      POINTER,     OPTIONAL :: xtteblck 
    ! ub_ak_20170705-

    LOGICAL,                         INTENT(OUT), OPTIONAL :: l_tfstd
    LOGICAL,                         INTENT(OUT), OPTIONAL :: l_init
    LOGICAL,                         INTENT(OUT), OPTIONAL :: l_enable

    ! LOCAL
    INTEGER :: zid = 0

    status = 100 ! ERROR: setname too long
    IF (LEN(setname) > STRLEN_TRSET) RETURN

    status = 103 ! ERROR: set not available
    DO zid = 1, NSETID
       IF (TRIM(trset(zid)%name) == TRIM(setname)) THEN
          status = 0
          EXIT
       END IF
    END DO
    IF (status /= 0) RETURN

    IF (PRESENT(ntrac))  ntrac  =  trset(zid)%ntrac
    IF (PRESENT(trlist)) trlist => trset(zid)%tilist
    IF (PRESENT(ti))     ti     => trset(zid)%ti

    IF (PRESENT(xt))     xt     => trset(zid)%xt
    IF (PRESENT(xtte))   xtte   => trset(zid)%xtte
    IF (PRESENT(xtm1))   xtm1   => trset(zid)%xtm1
    IF (PRESENT(xmem))   xmem   => trset(zid)%xmem

    IF (PRESENT(xt_tl))  xt_tl  => trset(zid)%xt_tl ! um_ak_20130530
    ! ub_ak_20170705+
    IF (PRESENT(ntracblck)) ntracblck =  trset(zid)%nblck
    IF (PRESENT(xtblck))    xtblck    => trset(zid)%xtblck   ! ub_ak_20170705
    IF (PRESENT(xtteblck))  xtteblck  => trset(zid)%xtteblck ! ub_ak_20170705
    ! ub_ak_20170705-

    IF (PRESENT(l_tfstd))  l_tfstd  = trset(zid)%l_tfstd
    IF (PRESENT(l_init))   l_init   = trset(zid)%l_init
    IF (PRESENT(l_enable)) l_enable = trset(zid)%l_enable

    status = 0

  END SUBROUTINE get_tracer_set_by_name
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_set_by_id(status, setid, setname, trlist, ti, ntrac &
!    , xt, xtte, xtm1, xmem, l_tfstd, l_init, l_enable)        !um_ak_20130521
     , xt, xtte, xtm1, xmem, xt_tl, l_tfstd, l_init, l_enable &!um_ak_20130521
     , ntracblck, xtblck, xtteblck) ! ub_ak_20170705

    IMPLICIT NONE
    INTRINSIC :: PRESENT

    ! I/O
    INTEGER,                         INTENT(OUT)           :: status
    INTEGER,                         INTENT(IN)            :: setid
    CHARACTER(LEN=STRLEN_TRSET),     INTENT(OUT), OPTIONAL :: setname
    TYPE(t_trinfo_list),             POINTER,     OPTIONAL :: trlist
    TYPE(t_trinfo_tp), DIMENSION(:), POINTER,     OPTIONAL :: ti
    INTEGER,                         INTENT(OUT), OPTIONAL :: ntrac
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xt
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xtte
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xtm1
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xmem
    REAL(DP), DIMENSION(:,:,:,:,:),  POINTER,     OPTIONAL :: xt_tl!um_ak_20130521
    ! ub_ak_20170705+
    INTEGER,                         INTENT(OUT), OPTIONAL :: ntracblck
    REAL(DP), DIMENSION(:,:,:,:),    POINTER,     OPTIONAL :: xtblck
    REAL(DP), DIMENSION(:,:,:),      POINTER,     OPTIONAL :: xtteblck
    ! ub_ak_20170705-
    LOGICAL,                         INTENT(OUT), OPTIONAL :: l_tfstd
    LOGICAL,                         INTENT(OUT), OPTIONAL :: l_init
    LOGICAL,                         INTENT(OUT), OPTIONAL :: l_enable

    IF (setid > NSETID) THEN
       status = 104 ! TRACER-SET ID NOT EXISTENT
       RETURN
    ENDIF

    IF (PRESENT(setname)) setname   =  trset(setid)%name
    IF (PRESENT(trlist))  trlist    => trset(setid)%tilist
    IF (PRESENT(ti))      ti        => trset(setid)%ti
    IF (PRESENT(ntrac))   ntrac     =  trset(setid)%ntrac

    IF (PRESENT(xt))      xt        => trset(setid)%xt
    IF (PRESENT(xtte))    xtte      => trset(setid)%xtte
    IF (PRESENT(xtm1))    xtm1      => trset(setid)%xtm1
    IF (PRESENT(xmem))    xmem      => trset(setid)%xmem

    IF (PRESENT(xt_tl))   xt_tl     => trset(setid)%xt_tl ! um_ak_20130521
    ! ub_ak_20170705+
    IF (PRESENT(ntracblck)) ntracblck =  trset(setid)%nblck
    IF (PRESENT(xtblck))    xtblck    => trset(setid)%xtblck   ! ub_ak_20170705
    IF (PRESENT(xtteblck))  xtteblck  => trset(setid)%xtteblck ! ub_ak_20170705
    ! ub_ak_20170705-

    IF (PRESENT(l_tfstd))  l_tfstd  =  trset(setid)%l_tfstd
    IF (PRESENT(l_init))   l_init   =  trset(setid)%l_init
    IF (PRESENT(l_enable)) l_enable =  trset(setid)%l_enable

    status = 0

  END SUBROUTINE get_tracer_set_by_id
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_set_id(status, setname, id)

    IMPLICIT NONE
    INTRINSIC :: LEN, TRIM

    ! I/O
    INTEGER,                        INTENT(OUT)  :: status
    CHARACTER(LEN=*),               INTENT(IN)   :: setname
    INTEGER,                        INTENT(OUT)  :: id

    status = 100 ! ERROR: setname too long
    IF (LEN(setname) > STRLEN_TRSET) RETURN

    status = 103 ! ERROR: set not available
    DO id = 1, NSETID
       IF (TRIM(trset(id)%name) == TRIM(setname)) THEN
          status = 0
          EXIT
       END IF
    END DO
    IF (status /= 0) RETURN

    status = 0

  END SUBROUTINE get_tracer_set_id
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE setup_tracer_set(status, setname, dim, nt, l_tfstd, l_init &
       , nbounds, l2tls, dimblck, ntblck) ! ub_ak_20170705

    IMPLICIT NONE
    INTRINSIC :: NULL

    ! I/O
    INTEGER,               INTENT(OUT) :: status
    CHARACTER(LEN=*),      INTENT(IN)  :: setname ! name of tracer set
    INTEGER, DIMENSION(3), INTENT(IN)  :: dim     ! 3 free dimensions
    INTEGER,               INTENT(IN)  :: nt      ! number of (time-) levels
    LOGICAL,               INTENT(IN)  :: l_tfstd ! level 2 and 3 specific
    LOGICAL,               INTENT(IN)  :: l_init  ! allow init
    ! mz_ab_20100509+
    INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL :: nbounds    
    LOGICAL,               INTENT(IN), OPTIONAL :: l2tls ! ub_ak_20160617
    INTEGER,               INTENT(IN), OPTIONAL :: dimblck ! blocklength
    INTEGER,               INTENT(IN), OPTIONAL :: ntblck  ! number timelevel block
    ! mz_ab_20100509-

    ! LOCAL
    INTEGER :: zid, zstat
    TYPE(t_trinfo_list),   POINTER     :: til => NULL()
    INTEGER                            :: jt
    INTEGER, DIMENSION(3)              :: zdim ! mz_ab_20100509
    INTEGER                            :: ntb

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! SETUP ENUMERATED LIST OF TRACER META INFORMATION STRUCTURES
    ALLOCATE(trset(zid)%ti(trset(zid)%ntrac))
    til => trset(zid)%tilist
    DO jt=1, trset(zid)%ntrac
       IF (jt /= til%info%ident%idx) THEN
          status = 999 
          RETURN
       END IF
       trset(zid)%ti(jt)%tp => til%info
       til => til%next
    END DO

    ! CHECK LEVELS
    IF (nt < 1) THEN
       status = 105 ! WRONG NUMBER OF (TIME-) LEVELS
       RETURN
    END IF
    trset(zid)%nt = nt

    status = 1000 ! ERROR: MEMORY ALLOCATION FAILED
    ! mz_ab_20100528+
    IF (PRESENT(nbounds)) THEN
       zdim(:)=dim(:)+2*nbounds(:) 
    ELSE 
       zdim(:) = dim(:)
    END IF

    SELECT CASE (TRRANK)
    CASE(3)
       ALLOCATE(trset(zid)%mem( zdim(1), zdim(2) &
            , trset(zid)%ntrac, zdim(3), 1, nt)  &
            , STAT=zstat)
       ! ub_ak_20170705+
       ! ALLOCATE ADDITIONAL BLOCK MEMORY
       IF (trset(zid)%nblck  > 0) THEN
          IF (.NOT. PRESENT(dimblck)) THEN
             status = 1002 ! TRACER BLOCKING REQUESTED / dimblck undefined
             RETURN
          ELSE
             IF (dimblck <= 0) THEN
                status = 1003 ! TRACER BLOCKING REQUESTED / dimblck <= 0
             END IF
          END IF
          IF (PRESENT(ntblck)) THEN
             ntb = ntblck
          ELSE
             ntb = 2
          END IF
          ALLOCATE(trset(zid)%xtblck(dimblck,trset(zid)%nblck,zdim(3),ntb)&
               , STAT=zstat)
          IF (zstat /= 0) RETURN
          ALLOCATE(trset(zid)%xtteblck(dimblck,trset(zid)%nblck,zdim(3))&
               , STAT=zstat)
          IF (zstat /= 0) RETURN
          trset(zid)%xtblck(:,:,:,:)  = 0.0_DP
          trset(zid)%xtteblck(:,:,:)  = 0.0_DP
       ENDIF
       ! ub_ak_20170705-
    CASE(4)
       ALLOCATE(trset(zid)%mem(zdim(1),zdim(2)   &
            , zdim(3), trset(zid)%ntrac, 1, nt)  &
            , STAT=zstat)
       ! ub_ak_20170705+
       ! ALLOCATE ADDITIONAL BLOCK MEMORY
       IF (trset(zid)%nblck  > 0) THEN
          IF (.NOT. PRESENT(dimblck)) THEN
             status = 1002 ! TRACER BLOCKING REQUESTED / dimblck undefined
             RETURN
          ELSE
             IF (dimblck <= 0) THEN
                status = 1003 ! TRACER BLOCKING REQUESTED / dimblck <= 0
             END IF
          END IF
          IF (PRESENT(ntblck)) THEN
             ntb = ntblck
          ELSE
             ntb = 2
          END IF
          ALLOCATE(trset(zid)%xtblck(dimblck,zdim(3),trset(zid)%nblck,ntb)&
               , STAT=zstat)
          IF (zstat /= 0) RETURN
          ALLOCATE(trset(zid)%xtteblck(dimblck,zdim(3),trset(zid)%nblck)&
               , STAT=zstat)
          IF (zstat /= 0) RETURN
          trset(zid)%xtblck(:,:,:,:)  = 0.0_DP
          trset(zid)%xtteblck(:,:,:)  = 0.0_DP
       ENDIF
       ! ub_ak_20170705-
    END SELECT

    IF (zstat /= 0) RETURN
    trset(zid)%mem(:,:,:,:,:,:) = 0.0_DP

    trset(zid)%xt   => trset(zid)%mem(:,:,:,:,:,1)

    ! um_ak_20150309+
    IF  ((l_tfstd) .AND. (nt >= 5)) THEN
       NULLIFY(trset(zid)%xt)
       trset(zid)%l_tfstd = .TRUE.
       trset(zid)%xtte => trset(zid)%mem(:,:,:,:,:,1)
       trset(zid)%xtm1 => trset(zid)%mem(:,:,:,:,:,2)
       trset(zid)%xt   => trset(zid)%mem(:,:,:,:,:,3)
       trset(zid)%xmem => trset(zid)%mem(:,:,:,:,1,4:)
       trset(zid)%nex = 4
       ! ub_ak_20170215+
       IF (PRESENT(l2tls)) THEN
          IF (l2tls) THEN
             trset(zid)%xt_tl=> trset(zid)%mem(:,:,:,:,1,2:3)
          ELSE
             trset(zid)%xt_tl=> trset(zid)%mem(:,:,:,:,1,2:4)
          ENDIF
       ELSE
       ! ub_ak_20170215-
          trset(zid)%xt_tl=> trset(zid)%mem(:,:,:,:,1,2:4)
       ENDIF ! ub_ak_20170215
    ! um_ak_20150309-
    ELSE IF ((l_tfstd) .AND. (nt >= 3)) THEN
    ! um_ak_20130611-
!    IF ((l_tfstd) .AND. (nt >= 3)) THEN
       trset(zid)%l_tfstd = .TRUE.
       trset(zid)%xtte => trset(zid)%mem(:,:,:,:,:,2)
       trset(zid)%xtm1 => trset(zid)%mem(:,:,:,:,:,3)
       IF (nt > 3) THEN
          trset(zid)%xmem => trset(zid)%mem(:,:,:,:,1,4:)
          trset(zid)%nex = 4
       ELSE
          trset(zid)%xmem => NULL() 
          trset(zid)%nex = 0
       END IF
       NULLIFY(trset(zid)%xt_tl) ! um_ak_20150409
    ELSE
       trset(zid)%l_tfstd = .FALSE.
       trset(zid)%xtte => NULL()
       trset(zid)%xtm1 => NULL()
       IF (nt > 1) THEN
          trset(zid)%xmem => trset(zid)%mem(:,:,:,:,1,2:)
          trset(zid)%nex = 2
       ELSE
          trset(zid)%xmem => NULL() 
          trset(zid)%nex = 0
       END IF
       NULLIFY(trset(zid)%xt_tl) ! um_ak_20150409
    END IF

    trset(zid)%l_init = l_init

    status = 0

  END SUBROUTINE setup_tracer_set
  ! -------------------------------------------------------------------

  ! um_ak_20130527+
  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_set_info(status, setname, ntrac, nt, l_enabled &
      , nblck )! ub_ak_20170705
    ! ub_ak_20160914: l_enabled added

    IMPLICIT NONE
 
    ! I/O
    INTEGER,            INTENT(OUT) :: status
    CHARACTER(LEN=*),   INTENT(IN)  :: setname 
    INTEGER, OPTIONAL,  INTENT(OUT) :: ntrac ! current number of tracers
    INTEGER, OPTIONAL,  INTENT(OUT) :: nt    ! number of time-levels
    LOGICAL, OPTIONAL,  INTENT(OUT) :: l_enabled  ! is the tracer set enabled?
    INTEGER, OPTIONAL,  INTENT(OUT) :: nblck ! current number of tracers
    
    ! LOCAL
    INTEGER                         :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    IF (PRESENT(ntrac)) THEN
       ntrac = trset(zid)%ntrac
    END IF

    IF (PRESENT(nt)) THEN
       nt = trset(zid)%nt
    ENDIF

    IF (PRESENT(l_enabled)) THEN
      l_enabled = trset(zid)%l_enable
    ENDIF

    ! ub_ak_20170705+
    IF (PRESENT(nblck)) THEN
       nblck = trset(zid)%nblck
    END IF
    ! ub_ak_20170705-

    status = 0

  END SUBROUTINE get_tracer_set_info
  ! -------------------------------------------------------------------
  ! um_ak_20130527-

  ! -------------------------------------------------------------------
  SUBROUTINE clean_tracer_set(status, setname)
    
    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,               INTENT(OUT) :: status
    CHARACTER(LEN=*),      INTENT(IN)  :: setname
    
    ! LOCAL
    INTEGER                      :: zid, zstat
    TYPE(t_trinfo_list), POINTER :: ti => NULL()
    TYPE(t_trinfo_list), POINTER :: te => NULL()
    INTEGER                      :: jt
    
    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    status = 1001 ! ERROR: MEMORY DEALLOCATION FAILED
    IF (ASSOCIATED(trset(zid)%mem)) THEN
      DEALLOCATE(trset(zid)%mem,   STAT=zstat)
      IF (zstat /= 0) RETURN
    END IF
    
    ! ub_ak_20170705+
    IF (ASSOCIATED(trset(zid)%xtblck)) THEN
      DEALLOCATE(trset(zid)%xtblck,   STAT=zstat)
      IF (zstat /= 0) RETURN
    END IF
    IF (ASSOCIATED(trset(zid)%xtteblck)) THEN
      DEALLOCATE(trset(zid)%xtteblck,   STAT=zstat)
      IF (zstat /= 0) RETURN
    END IF
    ! ub_ak_20170705-

    ! DELETE ENUMERATED LIST OF TRACER META INFORMATION
    DO jt=1, trset(zid)%ntrac
       NULLIFY(trset(zid)%ti(jt)%tp)
    END DO
    IF (ASSOCIATED(trset(zid)%ti)) DEALLOCATE(trset(zid)%ti)

    ! DELETE TRACER INFO STRUCT
    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       te => ti%next
       DEALLOCATE(ti, STAT=zstat)
       IF (zstat /= 0) THEN
          status = 1001 ! ERRRO: memory deallocation failed
          RETURN
       END IF
       ti => te
    END DO
    
    status = 0

  END SUBROUTINE clean_tracer_set
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE print_tracer_set

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, TRIM

    ! LOCAL
    INTEGER                        :: zid, ntrac, i, l
    TYPE(t_trinfo_list), POINTER   :: ti => NULL()
    TYPE(t_trinfo_list), POINTER   :: te => NULL()
    CHARACTER(LEN=200)                 :: form = ''
    CHARACTER(LEN=NAMES_CASK_STRLEN)   :: tmpstr = ''
    CHARACTER(LEN=4*MAX_CASK_I+1)      :: oline = ''

    WRITE(form,*) '(1x,i3,1x,a15,1x,a8,1x,3(i1,1x),',MAX_CASK_I,'(1x,i3))'

    set_loop: DO zid = 1, NSETID

    ! NUMBER OF TRACERS
    ntrac = trset(zid)%ntrac

    WRITE(*,*) '============================================================='
    WRITE(*,*) 'TRACER-SET        : ', TRIM(trset(zid)%name)
    IF (trset(zid)%l_enable) THEN
       WRITE(*,*) 'STATUS            : enabled'
    ELSE
       WRITE(*,*) 'STATUS            : disabled'
    END IF
    WRITE(*,*) 'NUMBER OF TRACERS : ', ntrac
    WRITE(*,*)
    ! ub_ak_20170705+
    IF (trset(zid)%nblck > 0) THEN
       WRITE(*,*) 'NUMBER OF BLOCKED TRACERS : ',  trset(zid)%nblck
       WRITE(*,*)
    END IF
    ! ub_ak_20170705-

    IF (ntrac > 0) THEN
    WRITE(*,*) &
         '----IDENT-------------------------  ----META--------------------- '

    DO l=1, NAMES_CASK_STRLEN
       DO i=1, MAX_CASK_I
          tmpstr = NAMES_CASK_I(i)
!!$          WRITE(oline(2*i-1:2*i),'(a1,a3)') ' ',tmpstr(l:l)
          WRITE(oline((i-1)*4+1:4*i),'(a1,a3)') ' ',tmpstr(l:l)
       END DO
       WRITE(*,*) '                                   ', oline
    END DO
    WRITE(*,*) '                             M Q T                           '
    WRITE(*,*) 'idt fullname(1:15)  submodel                                 '
    WRITE(*,*) 
    END IF

    ! LOOP OVER TRACERS IN SET
    ti => trset(zid)%tilist
    tracer_loop: DO 
       IF (.NOT. ASSOCIATED(ti)) EXIT

       WRITE(*,form) &
          ti%info%ident%idx,             &
          ti%info%ident%fullname(1:15),  &
          ti%info%ident%submodel(1:8),   &
          ti%info%ident%medium,          &
          ti%info%ident%quantity,        &
          ti%info%ident%type,            &
          ti%info%meta%cask_i(:)

        ti => ti%next
    END DO tracer_loop

    WRITE(*,*) &
         '=================================================================='

    END DO set_loop

  END SUBROUTINE print_tracer_set
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE print_tracer_set_val

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, MINVAL, MAXVAL, TRIM, SIZE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: fstr1 = '(3x,a4,2x,e12.4,a5,e12.4,2x,a)'
    CHARACTER(LEN=*), PARAMETER  :: fstr2 = '(3x,i4.4,2x,e12.4,a5,e12.4,2x,a)'
    ! ub_ak_20170705+
    CHARACTER(LEN=*), PARAMETER  :: fstr3 = &
                                        '(3x,a8,i4.4,2x,e12.4,a5,e12.4,2x,a)'
    ! ub_ak_20170705-
    INTEGER                      :: zid, ntrac, jt, jtb
    TYPE(t_trinfo_list), POINTER :: ti => NULL()
    TYPE(t_trinfo_list), POINTER :: te => NULL()
    INTEGER                      :: n

    set_loop: DO zid=1, NSETID

    ! NUMBER OF TRACERS
    ntrac = trset(zid)%ntrac

    WRITE(*,*) '============================================================='
    WRITE(*,*) 'TRACER-SET        : ', TRIM(trset(zid)%name)
    IF (trset(zid)%l_enable) THEN
       WRITE(*,*) 'STATUS            : enabled'
    ELSE
       WRITE(*,*) 'STATUS            : disabled'
    END IF
    WRITE(*,*) 'NUMBER OF TRACERS : ', ntrac
    WRITE(*,*)
   ! ub_ak_20170705+
    IF (trset(zid)%nblck > 0) THEN
       WRITE(*,*) 'NUMBER OF BLOCKED TRACERS : ',  trset(zid)%nblck
       WRITE(*,*)
    END IF
    ! ub_ak_20170705-

    ! LOOP OVER TRACERS IN SET
    ti => trset(zid)%tilist
    tracer_loop: DO 
       IF (.NOT. ASSOCIATED(ti)) EXIT

       jt = ti%info%ident%idx
       WRITE(*,*) jt &
            , ti%info%ident%basename,' - ', ti%info%ident%subname, ' - '  &
            , ti%info%ident%submodel

       SELECT CASE (TRRANK)

       CASE(3) !---------------------------------------------------------

          WRITE(*,fstr1) 'XT  ' &
               , MINVAL(trset(zid)%xt(:,:,jt,:,:)), ' ... '  &
               , MAXVAL(trset(zid)%xt(:,:,jt,:,:))           &
               , ti%info%ident%unit

          IF (ASSOCIATED(trset(zid)%xtm1)) THEN
             WRITE(*,fstr1) 'XTM1' &
                  , MINVAL(trset(zid)%xtm1(:,:,jt,:,:)), ' ... ' &
                  , MAXVAL(trset(zid)%xtm1(:,:,jt,:,:))          &
                  , ti%info%ident%unit
          END IF

          IF (ASSOCIATED(trset(zid)%xtte)) THEN
             WRITE(*,fstr1) 'XTTE' &
                  , MINVAL(trset(zid)%xtte(:,:,jt,:,:)), ' ... ' &
                  , MAXVAL(trset(zid)%xtte(:,:,jt,:,:))          &
                  , ti%info%ident%unit

          END IF

          IF (ASSOCIATED(trset(zid)%xmem)) THEN
             DO n=1, SIZE(trset(zid)%xmem, 5)
                WRITE(*,fstr2) n &
                     , MINVAL(trset(zid)%xmem(:,:,jt,:,n)), ' ... ' &
                     , MAXVAL(trset(zid)%xmem(:,:,jt,:,n))          &
                     , ti%info%ident%unit

             END DO
          END IF
          ! ub_ak_20170705+
          IF (ASSOCIATED(trset(zid)%xtblck)) THEN
             jtb = ti%info%ident%idxblck
             IF (jtb > 0) THEN
                DO n=1, SIZE(trset(zid)%xtblck, 4)
                   WRITE(*,fstr3) 'xtblck  ',n &
                     , MINVAL(trset(zid)%xtblck(:,jtb,:,n)), ' ... ' &
                     , MAXVAL(trset(zid)%xtblck(:,jtb,:,n))          &
                     , ti%info%ident%unit
                   
                END DO
                WRITE(*,fstr3) 'xtteb1ck',n &
                     , MINVAL(trset(zid)%xtteblck(:,jtb,:)), ' ... ' &
                     , MAXVAL(trset(zid)%xtteblck(:,jtb,:))          &
                     , ti%info%ident%unit
             END IF                   
          END IF
          ! ub_ak_20170705-

       CASE(4) !---------------------------------------------------------

          WRITE(*,fstr1) 'XT  ' &
               , MINVAL(trset(zid)%xt(:,:,:,jt,:)), ' ... '  &
               , MAXVAL(trset(zid)%xt(:,:,:,jt,:))           &
               , ti%info%ident%unit

          IF (ASSOCIATED(trset(zid)%xtm1)) THEN
             WRITE(*,fstr1) 'XTM1' &
                  , MINVAL(trset(zid)%xtm1(:,:,:,jt,:)), ' ... ' &
                  , MAXVAL(trset(zid)%xtm1(:,:,:,jt,:))          &
                  , ti%info%ident%unit
          END IF

          IF (ASSOCIATED(trset(zid)%xtte)) THEN
             WRITE(*,fstr1) 'XTTE' &
                  , MINVAL(trset(zid)%xtte(:,:,:,jt,:)), ' ... ' &
                  , MAXVAL(trset(zid)%xtte(:,:,:,jt,:))          &
                  , ti%info%ident%unit

          END IF

          IF (ASSOCIATED(trset(zid)%xmem)) THEN
             DO n=1, SIZE(trset(zid)%xmem, 5)
                WRITE(*,fstr2) n &
                     , MINVAL(trset(zid)%xmem(:,:,:,jt,n)), ' ... ' &
                     , MAXVAL(trset(zid)%xmem(:,:,:,jt,n))          &
                     , ti%info%ident%unit

             END DO
          END IF

          ! ub_ak_20170705+
          IF (ASSOCIATED(trset(zid)%xtblck)) THEN
             jtb = ti%info%ident%idxblck
             IF (jtb > 0) THEN
                DO n=1, SIZE(trset(zid)%xtblck, 4)
                   WRITE(*,fstr3) 'xtblck  ',n &
                     , MINVAL(trset(zid)%xtblck(:,:,jtb,n)), ' ... ' &
                     , MAXVAL(trset(zid)%xtblck(:,:,jtb,n))          &
                     , ti%info%ident%unit
                   
                END DO
                WRITE(*,fstr3) 'xtteb1ck',n &
                     , MINVAL(trset(zid)%xtteblck(:,:,jtb)), ' ... ' &
                     , MAXVAL(trset(zid)%xtteblck(:,:,jtb))          &
                     , ti%info%ident%unit
             END IF
          END IF
          ! ub_ak_20170705-
       END SELECT

       ti => ti%next
    END DO tracer_loop

    END DO set_loop

    WRITE(*,*) '============================================================='

  END SUBROUTINE print_tracer_set_val
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_by_name(status, setname, basename &
       , subname, idx, fullname, longname, unit, submodel &
       , medium, quantity, type                           &
       , trinfo, pxt, pxtm1, pxtte, pxmem                 &
       , idxblck, pxtblck, pxtteblck) ! ub_ak_20170705

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, INDEX, LEN, PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    !
    CHARACTER(LEN=*), INTENT(IN)            :: basename
    !
    CHARACTER(LEN=*), INTENT(IN),  OPTIONAL :: subname
    INTEGER,          INTENT(OUT), OPTIONAL :: idx
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: fullname
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: longname
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: unit
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: submodel
    INTEGER,          INTENT(OUT), OPTIONAL :: medium
    INTEGER,          INTENT(OUT), OPTIONAL :: quantity
    INTEGER,          INTENT(OUT), OPTIONAL :: type
    !
    TYPE(t_trinfo),   INTENT(OUT), OPTIONAL :: trinfo
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxt
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxtm1
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxtte
    REAL(DP), DIMENSION(:,:,:,:), POINTER, OPTIONAL :: pxmem

    ! ub_ak_20170705+
    INTEGER,  INTENT(OUT)                , OPTIONAL :: idxblck 
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxtblck
    REAL(DP), DIMENSION(:,:),     POINTER, OPTIONAL :: pxtteblck
    ! ub_ak_20170705-

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    CHARACTER(LEN=STRLEN_FNAME)   :: zfullname
    INTEGER                       :: zid
    INTEGER                       :: zidx
    INTEGER                       :: zidxb ! ub_ak_20170705

    ! INIT
    IF (PRESENT(idx))        idx       = 0
    IF (PRESENT(fullname))   fullname  = ''
    IF (PRESENT(longname))   longname  = ''
    IF (PRESENT(unit))       unit      = ''
    IF (PRESENT(submodel))   submodel  = ''
    IF (PRESENT(medium))     medium    = 0
    IF (PRESENT(quantity))   quantity  = 0
    IF (PRESENT(type))       type      = 0
    IF (PRESENT(idxblck))    idxblck   = 0 ! ub_ak_20170705
    !
    ! trinfo
    !
    IF (PRESENT(pxt))   NULLIFY(pxt)
    IF (PRESENT(pxtm1)) NULLIFY(pxtm1)
    IF (PRESENT(pxtte)) NULLIFY(pxtte)
    IF (PRESENT(pxmem)) NULLIFY(pxmem)
    ! ub_ak_20170705+
    IF (PRESENT(pxtblck))   NULLIFY(pxtblck)
    IF (PRESENT(pxtteblck)) NULLIFY(pxtteblck)
    ! ub_ak_20170705-

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    status = 208 ! ERROR: basename MUST NOT contain '_'
    IF (INDEX(basename, '_') /= 0) RETURN

    status = 200 ! ERROR: tracer basename too long
    IF (LEN(basename) > STRLEN_MEDIUM) RETURN
    
    IF (PRESENT(subname)) THEN
       status = 201 ! ERROR: tracer subname too long
       IF (LEN(subname) > STRLEN_MEDIUM) RETURN
       IF (TRIM(subname) /= '') THEN
          zfullname = TRIM(basename)//'_'//TRIM(subname)
       ELSE
          zfullname = TRIM(basename)
       END IF
    ELSE
       zfullname = TRIM(basename)
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (TRIM(ti%info%ident%fullname) == TRIM(zfullname)) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    ! RETURN TRACER INFORMATION
    ! - IDENTIFICATION
    zidx = ti%info%ident%idx
    IF (PRESENT(idx)) idx = zidx
    ! ub_ak_20170705+
    zidxb = ti%info%ident%idxblck
    IF (PRESENT(idxblck)) idxblck = zidxb
    ! ub_ak_20170705-
    !
    IF (PRESENT(pxt)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xt)) THEN
          status = 600  ! ERROR: xt pointer not associated
          RETURN
       ELSE
          SELECT CASE (TRRANK)
          CASE(3)
             pxt => trset(zid)%xt(:,:,zidx,:,1)
          CASE(4)
             pxt => trset(zid)%xt(:,:,:,zidx,1)
          END SELECT
       END IF
    END IF
    !
    IF (PRESENT(pxtm1)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xtm1)) THEN
          status = 601  ! ERROR: xtm1 pointer not associated
          RETURN
       ELSE
          SELECT CASE(TRRANK)
          CASE(3)
             pxtm1 => trset(zid)%xtm1(:,:,zidx,:,1)
          CASE(4)
             pxtm1 => trset(zid)%xtm1(:,:,:,zidx,1)
          END SELECT
       END IF
    END IF
    !
    IF (PRESENT(pxtte)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xtte)) THEN
          status = 602  ! ERROR: xtte pointer not associated
          RETURN
       ELSE
          SELECT CASE(TRRANK)
          CASE(3)
             pxtte => trset(zid)%xtte(:,:,zidx,:,1)
          CASE(4)
             pxtte => trset(zid)%xtte(:,:,:,zidx,1)
          END SELECT
       END IF
    END IF
    !
    IF (PRESENT(pxmem)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xmem)) THEN
          status = 603  ! ERROR: xmem pointer not associated
          RETURN
       ELSE
          SELECT CASE(TRRANK)
          CASE(3)
             pxmem => trset(zid)%xmem(:,:,zidx,:,:)
          CASE(4)
             pxmem => trset(zid)%xmem(:,:,:,zidx,:)
          END SELECT
       END IF
    END IF

    ! ub_ak_20170705+
    IF (PRESENT(pxtblck)) THEN
       IF (zidxb > 0) THEN
          IF (.NOT. ASSOCIATED(trset(zid)%xtblck)) THEN
             status = 604  ! ERROR: xtblck pointer not associated
             RETURN
          ELSE
             SELECT CASE(TRRANK)
             CASE(3)
                pxtblck => trset(zid)%xtblck(:,zidxb,:,:)
             CASE(4)
                pxtblck => trset(zid)%xtblck(:,:,zidxb,:)
             END SELECT
          END IF
       ELSE
          NULLIFY(pxtblck)
       END IF
    END IF
    IF (PRESENT(pxtteblck)) THEN
       IF (zidxb > 0) THEN
          IF (.NOT. ASSOCIATED(trset(zid)%xtteblck)) THEN
             status = 604  ! ERROR: xtteblck pointer not associated
             RETURN
          ELSE
             SELECT CASE(TRRANK)
             CASE(3)
                pxtteblck => trset(zid)%xtteblck(:,zidxb,:)
             CASE(4)
                pxtteblck => trset(zid)%xtteblck(:,:,zidxb)
             END SELECT
          END IF
       ELSE
          NULLIFY(pxtteblck)
       END IF
    END IF
    ! ub_ak_20170705-

    ! - OPTIONAL INFORMATION
    IF (PRESENT(fullname)) THEN
       IF (LEN(fullname) > STRLEN_FNAME) THEN
          status = 206 ! ERROR: fullname too long
          RETURN
       END IF
       fullname = TRIM(ti%info%ident%fullname)
    END IF
    !
    !                                    ! LONGNAME
    IF (PRESENT(longname)) THEN
       IF (LEN(longname) > STRLEN_VLONG) THEN
          status = 207 ! ERROR: longname too long
          RETURN
       END IF
       longname = TRIM(ti%info%ident%longname)
    END IF
    !
    IF (PRESENT(unit)) THEN              ! UNIT
       IF (LEN(unit) > STRLEN_MEDIUM) THEN
          status = 203 ! ERROR: unit too long
          RETURN
       END IF
       unit = TRIM(ti%info%ident%unit)
    END IF
    !
    IF (PRESENT(submodel)) THEN          ! SUBMODEL
       IF (LEN(submodel) > STRLEN_MEDIUM) THEN
          status = 204 ! ERROR: submodel too long
          RETURN
       END IF
       submodel = TRIM(ti%info%ident%submodel)
    END IF
    !
    IF (PRESENT(medium)) THEN            ! MEDIUM
       IF ((ti%info%ident%medium < 0) .OR. &
            (ti%info%ident%medium > MAX_MEDIUM) ) THEN
          status = 300  ! ERROR: unknown medium
          RETURN
       END IF
       medium = ti%info%ident%medium
    END IF
    !
    IF (PRESENT(quantity)) THEN          ! QUANTITY
       IF ((ti%info%ident%quantity < 0) .OR. &
            (ti%info%ident%quantity > MAX_QUANTITY) ) &
            THEN
          status = 301  ! ERROR: unknown quantity
          RETURN
       END IF
       quantity = ti%info%ident%quantity
    END IF
    ! 
    IF (PRESENT(type)) THEN              ! TYPE
       IF ((ti%info%ident%type < 0) .OR. &
            (ti%info%ident%type > MAX_MEDIUM) ) THEN
          status = 302  ! ERROR: unknown type
          RETURN
       END IF
       type = ti%info%ident%type
    END IF

    IF (PRESENT(trinfo)) trinfo = ti%info

    status = 0

  END SUBROUTINE get_tracer_by_name
  ! -------------------------------------------------------------------

! mz_pj_20090507+
  ! -------------------------------------------------------------------
  ! ub_ak_20170707+
!  SUBROUTINE get_tracer_by_id(status, setname, id, basename & 
  SUBROUTINE get_tracer_by_id(status, setname, idtrac, basename & 
  ! ub_ak_20170707-
       , subname, fullname, longname, unit, submodel         &
       , medium, quantity, type                              &
       , trinfo, pxt, pxtm1, pxtte, pxmem                    &!)
       , idx, idxblck, pxtblck, pxtteblck) ! ub_ak_20170705

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, LEN, PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    ! ub_ak_20170705+
    !INTEGER,          INTENT(IN)            :: id
    INTEGER,          INTENT(IN)            :: idtrac
    ! ub_ak_20170705-
    !
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: basename
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: subname
    !
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: fullname
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: longname
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: unit
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: submodel
    INTEGER,          INTENT(OUT), OPTIONAL :: medium
    INTEGER,          INTENT(OUT), OPTIONAL :: quantity
    INTEGER,          INTENT(OUT), OPTIONAL :: type
    !
    TYPE(t_trinfo),   INTENT(OUT), OPTIONAL :: trinfo
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxt
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxtm1
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxtte
    REAL(DP), DIMENSION(:,:,:,:), POINTER, OPTIONAL :: pxmem

    ! ub_ak_20170705+
    ! tracer index if id <= 0; 
    ! NOTE: for the case, that the tracer is identified via its block index,
    ! we require an extra variable here, as idtrac (prior id)
    ! can not be INTENT(INOUT) as the subroutine get_tracer is
    ! often called in a loop over all tracers (which makes idtrac
    ! unchangable).
    INTEGER,  INTENT(OUT),                 OPTIONAL :: idx
    INTEGER,  INTENT(INOUT),               OPTIONAL :: idxblck
    REAL(DP), DIMENSION(:,:,:),   POINTER, OPTIONAL :: pxtblck
    REAL(DP), DIMENSION(:,:),     POINTER, OPTIONAL :: pxtteblck
    ! ub_ak_20170705-

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid
    INTEGER                       :: id

    ! INIT
    IF (PRESENT(basename))   basename  = ''
    IF (PRESENT(subname))    subname   = ''

    IF (PRESENT(fullname))   fullname  = ''
    IF (PRESENT(longname))   longname  = ''
    IF (PRESENT(unit))       unit      = ''
    IF (PRESENT(submodel))   submodel  = ''
    IF (PRESENT(medium))     medium    = 0
    IF (PRESENT(quantity))   quantity  = 0
    IF (PRESENT(type))       type      = 0
    !
    ! trinfo
    !
    IF (PRESENT(pxt))   NULLIFY(pxt)
    IF (PRESENT(pxtm1)) NULLIFY(pxtm1)
    IF (PRESENT(pxtte)) NULLIFY(pxtte)
    IF (PRESENT(pxmem)) NULLIFY(pxmem)
    ! ub_ak_20170705+
    IF (PRESENT(pxtblck))   NULLIFY(pxtblck)
    IF (PRESENT(pxtteblck)) NULLIFY(pxtteblck)
    ! ub_ak_20170705-

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist

    IF (idtrac <= 0 .AND. .NOT. PRESENT(idxblck)) THEN
       status = 400 ! ERROR: tracer id  or blck id must be present and > 0!
       RETURN
    END IF
    IF (idtrac > 0) THEN ! ub_ak_20170705
       DO
          IF (.NOT. ASSOCIATED(ti)) EXIT
          ! ub_ak_20170707+
          !IF (ti%info%ident%idx == id) EXIT   ! FOUND
          IF (ti%info%ident%idx == idtrac) EXIT   ! FOUND
          ! ub_ak_20170707-
          ti => ti%next
       END DO
       IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
          status = TR_NEXIST
          RETURN
       END IF
       id = idtrac
       ! ub_ak_20170705+
       IF (PRESENT(idxblck)) THEN
          IF (idxblck > 0) THEN
             IF (ti%info%ident%idxblck /= idxblck) THEN
                status = 401 ! ERROR: TRACER ID and BLOCK ID DO NOT MATCH 
                RETURN
             END IF
          ELSE
             idxblck =  ti%info%ident%idxblck
          END IF
       END IF
    ELSE
       IF (PRESENT(idxblck)) THEN
          DO
             IF (.NOT. ASSOCIATED(ti)) EXIT
             IF (ti%info%ident%idxblck == idxblck) EXIT   ! FOUND
             ti => ti%next
          END DO
          IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
             status = TR_NEXIST
             RETURN
          END IF
          id   = ti%info%ident%idx
       ELSE
          status = 400 ! ERROR: TRACER ID <= 0  / IDXBLCK not present
      END IF
   END IF
   ! ub_ak_20170705-


    ! RETURN TRACER INFORMATION
    ! - IDENTIFICATION
    IF (PRESENT(basename)) THEN
       IF (LEN(basename) > STRLEN_MEDIUM) THEN
          status = 200 ! ERROR: basename too long
          RETURN
       END IF
       basename = TRIM(ti%info%ident%basename)
    END IF

    IF (PRESENT(subname)) THEN
       IF (LEN(subname) > STRLEN_MEDIUM) THEN
          status = 201 ! ERROR: subname too long
          RETURN
       END IF
       subname = TRIM(ti%info%ident%subname)
    END IF

    IF (PRESENT(pxt)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xt)) THEN
          status = 600  ! ERROR: xt pointer not associated
          RETURN
       ELSE
          SELECT CASE(TRRANK)
          CASE(3)
             pxt => trset(zid)%xt(:,:,id,:,1)
          CASE(4)
             pxt => trset(zid)%xt(:,:,:,id,1)
          END SELECT
       END IF
    END IF
    !
    IF (PRESENT(pxtm1)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xtm1)) THEN
          status = 601  ! ERROR: xtm1 pointer not associated
          RETURN
       ELSE
          SELECT CASE(TRRANK)
          CASE(3)
             pxtm1 => trset(zid)%xtm1(:,:,id,:,1)
          CASE(4)
             pxtm1 => trset(zid)%xtm1(:,:,:,id,1)
          END SELECT
       END IF
    END IF
    !
    IF (PRESENT(pxtte)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xtte)) THEN
          status = 602  ! ERROR: xtte pointer not associated
          RETURN
       ELSE
          SELECT CASE(TRRANK)
          CASE(3)
             pxtte => trset(zid)%xtte(:,:,id,:,1)
          CASE(4)
             pxtte => trset(zid)%xtte(:,:,:,id,1)
          END SELECT
       END IF
    END IF
    !
    IF (PRESENT(pxmem)) THEN
       IF (.NOT. ASSOCIATED(trset(zid)%xmem)) THEN
          status = 603  ! ERROR: xmem pointer not associated
          RETURN
       ELSE
          SELECT CASE(TRRANK)
          CASE(3)
             pxmem => trset(zid)%xmem(:,:,id,:,:)
          CASE(4)
             pxmem => trset(zid)%xmem(:,:,:,id,:)
          END SELECT
       END IF
    END IF

    ! ub_ak_20170705+
    IF (PRESENT(pxtblck)) THEN
       IF (idxblck > 0) THEN
          IF (.NOT. ASSOCIATED(trset(zid)%xtblck)) THEN
             status = 604  ! ERROR: xtblck pointer not associated
             RETURN
          ELSE
             SELECT CASE(TRRANK)
             CASE(3)
                pxtblck => trset(zid)%xtblck(:,idxblck,:,:)
             CASE(4)
                pxtblck => trset(zid)%xtblck(:,:,idxblck,:)
             END SELECT
          END IF
       ELSE
          NULLIFY(pxtblck)
       END IF
    END IF
    IF (PRESENT(pxtteblck)) THEN
       IF (idxblck > 0) THEN
          IF (.NOT. ASSOCIATED(trset(zid)%xtteblck)) THEN
             status = 604  ! ERROR: xtteblck pointer not associated
             RETURN
          ELSE
             SELECT CASE(TRRANK)
             CASE(3)
                pxtteblck => trset(zid)%xtteblck(:,idxblck,:)
             CASE(4)
                pxtteblck => trset(zid)%xtteblck(:,:,idxblck)
             END SELECT
          END IF
       ELSE
          NULLIFY(pxtteblck)
       END IF
    END IF
    ! ub_ak_20170705-

    ! - OPTIONAL INFORMATION
    ! ub_ak_20170707+
    IF (PRESENT(idx)) idx = id
    ! ub_ak_20170707-

    IF (PRESENT(fullname)) THEN
       IF (LEN(fullname) > STRLEN_FNAME) THEN
          status = 206 ! ERROR: fullname too long
          RETURN
       END IF
       fullname = TRIM(ti%info%ident%fullname)
    END IF
    !
    !                                    ! LONGNAME
    IF (PRESENT(longname)) THEN
       IF (LEN(longname) > STRLEN_VLONG) THEN
          status = 207 ! ERROR: longname too long
          RETURN
       END IF
       longname = TRIM(ti%info%ident%longname)
    END IF
    !
    IF (PRESENT(unit)) THEN              ! UNIT
       IF (LEN(unit) > STRLEN_MEDIUM) THEN
          status = 203 ! ERROR: unit too long
          RETURN
       END IF
       unit = TRIM(ti%info%ident%unit)
    END IF
    !
    IF (PRESENT(submodel)) THEN          ! SUBMODEL
       IF (LEN(submodel) > STRLEN_MEDIUM) THEN
          status = 204 ! ERROR: submodel too long
          RETURN
       END IF
       submodel = TRIM(ti%info%ident%submodel)
    END IF
    !
    IF (PRESENT(medium)) THEN            ! MEDIUM
       IF ((ti%info%ident%medium < 0) .OR. &
            (ti%info%ident%medium > MAX_MEDIUM) ) THEN
          status = 300  ! ERROR: unknown medium
          RETURN
       END IF
       medium = ti%info%ident%medium
    END IF
    !
    IF (PRESENT(quantity)) THEN          ! QUANTITY
       IF ((ti%info%ident%quantity < 0) .OR. &
            (ti%info%ident%quantity > MAX_QUANTITY) ) &
            THEN
          status = 301  ! ERROR: unknown quantity
          RETURN
       END IF
       quantity = ti%info%ident%quantity
    END IF
    ! 
    IF (PRESENT(type)) THEN              ! TYPE
       IF ((ti%info%ident%type < 0) .OR. &
            (ti%info%ident%type > MAX_MEDIUM) ) THEN
          status = 302  ! ERROR: unknown type
          RETURN
       END IF
       type = ti%info%ident%type
    END IF

    IF (PRESENT(trinfo)) trinfo = ti%info

    status = 0

  END SUBROUTINE get_tracer_by_id
  ! -------------------------------------------------------------------
! mz_pj_20090507-

! op_pj_20100325+
  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_s(status, setname, idx, flag, s)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    INTEGER,          INTENT(IN)            :: flag
    CHARACTER(LEN=*), INTENT(OUT)           :: s

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF
    
    IF ( (flag < 1) .OR. (flag > MAX_CASK_S) ) THEN
       status = 702 ! STRING FLAG INDEX OUT OF RANGE
       RETURN
    END IF
    s = ti%info%meta%cask_s(flag)

    status = 0
    
  END SUBROUTINE get_tracer_s
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_i(status, setname, idx, flag, i)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    INTEGER,          INTENT(IN)            :: flag
    INTEGER,          INTENT(OUT)           :: i

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    IF ( (flag < 1) .OR. (flag > MAX_CASK_I) ) THEN
       status = 700 ! INTEGER FLAG INDEX OUT OF RANGE
       RETURN
    END IF
    i = ti%info%meta%cask_i(flag)

    status = 0

  END SUBROUTINE get_tracer_i
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_r(status, setname, idx, flag, r)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    INTEGER,          INTENT(IN)            :: idx
    INTEGER,          INTENT(IN)            :: flag
    REAL(DP),         INTENT(OUT)           :: r

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN
    
    ! CHECK IF TRACER SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    ! SEARCH FOR EXISITNG TRACER IN LIST
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == idx) EXIT   ! FOUND
       ti => ti%next
    END DO

    IF (.NOT. ASSOCIATED(ti)) THEN ! END OF LIST REACHED
       status = TR_NEXIST
       RETURN
    END IF

    IF ( (flag < 1) .OR. (flag > MAX_CASK_R) ) THEN
       status = 701 ! REAL FLAG INDEX OUT OF RANGE
       RETURN
    END IF
    r = ti%info%meta%cask_r(flag)

    status = 0

  END SUBROUTINE get_tracer_r
  ! -------------------------------------------------------------------
! op_pj_20100325-

  ! -------------------------------------------------------------------
  SUBROUTINE get_tracer_list(status, setname, basename, idxs, subnames)


    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, INDEX, LEN, PRESENT, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    CHARACTER(LEN=*), INTENT(IN)            :: basename
    INTEGER, DIMENSION(:), POINTER          :: idxs       ! INTENT(OUT)
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER, OPTIONAL :: subnames

    ! LOCAL
    TYPE(t_trinfo_list), POINTER  :: ti => NULL()
    INTEGER                       :: zid
    INTEGER                       :: n

    ! INIT
    IF (ASSOCIATED(idxs)) THEN
       DEALLOCATE(idxs)
       NULLIFY(idxs)
    END IF

    IF (PRESENT(subnames)) THEN
       IF (ASSOCIATED(subnames)) THEN
          DEALLOCATE(subnames)
          NULLIFY(subnames)
       END IF
    END IF

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    status = 208 ! ERROR: basename MUST NOT contain '_'
    IF (INDEX(basename, '_') /= 0) RETURN

    status = 200 ! ERROR: tracer basename too long
    IF (LEN(basename) > STRLEN_MEDIUM) RETURN
    
    ! SEARCH FOR EXISITNG TRACERS IN LIST
    n = 0
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (TRIM(ti%info%ident%basename) == TRIM(basename)) n=n+1 ! FOUND
       ti => ti%next
    END DO

    ! ALLOCATE OUTPUT
    ALLOCATE(idxs(n))
    idxs(:) = 0
    IF (PRESENT(subnames)) THEN
       ALLOCATE(subnames(n))
       subnames(:) = ''
    END IF

    ! SEARCH FOR EXISITNG TRACERS IN LIST
    n = 0
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (TRIM(ti%info%ident%basename) == TRIM(basename)) THEN
          n=n+1
          idxs(n) = ti%info%ident%idx
          IF (PRESENT(subnames)) subnames(n) = TRIM(ti%info%ident%subname)
       END IF
       ti => ti%next
    END DO

    status = 0

  END SUBROUTINE get_tracer_list
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
!!$  SUBROUTINE tracer_iniflag(status, setname, id, linit)
  SUBROUTINE tracer_iniflag(status, setname, id, lset, lget)

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: setname
    INTEGER,          INTENT(IN)  :: id
!!$    LOGICAL,          INTENT(IN)  :: linit
    LOGICAL,          INTENT(IN),  OPTIONAL  :: lset
    LOGICAL,          INTENT(OUT), OPTIONAL  :: lget

    ! LOCAL
    INTEGER                      :: zid
    TYPE(t_trinfo_list), POINTER :: ti => NULL()
    LOGICAL                      :: lfound

    IF (PRESENT(lget)) lget = .FALSE.

    ! TRACER SET
    CALL get_tracer_set_id(status, setname, zid)
    IF (status /= 0) RETURN

    ! CHECK, IF SET IS ACTIVE
    IF (.NOT. trset(zid)%l_enable) THEN
       status = 0
       RETURN
    END IF

    lfound = .FALSE.

    ! SEARCH FOR EXISITNG TRACER IN LIST WITH THIS ID
    ti => trset(zid)%tilist
    DO
       IF (.NOT. ASSOCIATED(ti)) EXIT
       IF (ti%info%ident%idx == id) THEN
          lfound = .TRUE.
          EXIT   ! FOUND
       END IF
       ti => ti%next
    END DO

    IF (.NOT. lfound) THEN
       status = TR_NEXIST ! TRACER DOES NOT EXIST
       RETURN
    END IF

    ! SET / GET FLAG
!!$    ti%info%meta%linit = linit
    IF (PRESENT(lset)) ti%info%meta%linit = lset
    IF (PRESENT(lget)) lget = ti%info%meta%linit

  END SUBROUTINE tracer_iniflag
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  FUNCTION tracer_error_str(status)

    USE messy_main_tools, ONLY: int2str

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=STRLEN_VLONG)           :: tracer_error_str
    INTEGER,                   INTENT(IN) :: status

    ! LOCAL
    CHARACTER(LEN=4) :: echar 

    CALL int2str(echar, status, cpad='0', cerr='*')

    SELECT CASE(echar)
       ! NO ERROR
    CASE('0000')
       tracer_error_str = 'E'//echar//': NO ERROR'

       ! TRACER-SET ERRORS
    CASE('0100')
       tracer_error_str = 'E'//echar//': TRACER-SET NAME TOO LONG'
    CASE('0101')
       tracer_error_str = 'E'//echar//': TRACER-SET NAME EXISTS ALREADY'
    CASE('0102')
       tracer_error_str = 'E'//echar//': TRACER-SET LIST FULL'
    CASE('0103')
       tracer_error_str = 'E'//echar//': TRACER-SET NAME NOT EXISTENT'
    CASE('0104')
       tracer_error_str = 'E'//echar//': TRACER-SET ID NOT EXISTENT'
    CASE('0105')
       tracer_error_str = 'E'//echar//': WRONG NUMBER OF (TIME-) LEVELS'

       ! TRACER ERRORS
    CASE('0200')
       tracer_error_str = 'E'//echar//': TRACER BASENAME TOO LONG'
    CASE('0201')
       tracer_error_str = 'E'//echar//': TRACER SUBNAME TOO LONG'
    CASE('0202')
       tracer_error_str = 'E'//echar//': TRACER EXISTS ALREADY'
    CASE('0203')
       tracer_error_str = 'E'//echar//': TRACER UNIT TOO LONG'
    CASE('0204')
       tracer_error_str = 'E'//echar//': TRACER SUBMODEL NAME TOO LONG'
    CASE('0205')
       tracer_error_str = 'E'//echar//': TRACER DOES NOT EXIST'
    CASE('0206')
       tracer_error_str = 'E'//echar//': TRACER FULLNAME TOO LONG'
    CASE('0207')
       tracer_error_str = 'E'//echar//': TRACER LONGNAME TOO LONG'
    CASE('0208')
     tracer_error_str = 'E'//echar//': TRACER BASENAME MUST NOT CONTAIN ''_'''

       ! TYPE ERRORS
    CASE('0300')
       tracer_error_str = 'E'//echar//': UNKNOWN TRACER MEDIUM'
    CASE('0301')
       tracer_error_str = 'E'//echar//': UNKNOWN TRACER QUANTITY'
    CASE('0302')
       tracer_error_str = 'E'//echar//': UNKNOWN TRACER TYPE'
       
       ! ub_ak_20170705+
       ! TRACER INDEX ERRORS 
     CASE('0400')
       tracer_error_str = 'E'//echar//': TRACER ID <=0  / IDXBLCK not present'
     CASE('0401')
       tracer_error_str = 'E'//echar//': TRACER ID and BLOCK ID DO NOT MATCH '
       ! ub_ak_20170705-

       ! POINTER ERRORS
    CASE('0600')
       tracer_error_str = 'E'//echar//': POINTER XT NOT ASSOCIATED'
    CASE('0601')
       tracer_error_str = 'E'//echar//': POINTER XTM1 NOT ASSOCIATED'
    CASE('0602')
       tracer_error_str = 'E'//echar//': POINTER XTTE NOT ASSOCIATED'
    CASE('0603')
       tracer_error_str = 'E'//echar//': POINTER XMEM NOT ASSOCIATED'
       ! ub_ak_20170705+
    CASE('0604')
       tracer_error_str = 'E'//echar//': POINTER XTBLCK NOT ASSOCIATED'
    CASE('0605')
       tracer_error_str = 'E'//echar//': POINTER XTTEBLCK NOT ASSOCIATED'
    CASE('0606')
       tracer_error_str = 'E'//echar//': TIMELEVEL FOR PTR XTBLCK NOT PRESENT'
       ! ub_ak_20170705-


       ! META INFORMATION ERRORS
    CASE('0700')
       tracer_error_str = 'E'//echar//': INTEGER FLAG INDEX OUT OF RANGE'
    CASE('0701')
       tracer_error_str = 'E'//echar//': REAL FLAG INDEX OUT OF RANGE'
    CASE('0702')
       tracer_error_str = 'E'//echar//': STRING FLAG INDEX OUT OF RANGE'
    CASE('0703')
       tracer_error_str = 'E'//echar//': CONTAINER UNKOWN'
       
       ! MEMORY ERRORS
    CASE('1000')
       tracer_error_str = 'E'//echar//': MEMORY ALLOCATION FAILED'
    CASE('1001')
       tracer_error_str = 'E'//echar//': MEMORY DEALLOCATION FAILED'
       ! ub_ak_20170705-
   CASE('1002')
       tracer_error_str = 'E'//echar//': TRACER BLOCKING REQUESTED / dimblck undefined'
   CASE('1003')
       tracer_error_str = 'E'//echar//': TRACER BLOCKING REQUESTED / dimblck <= 0'
       ! ub_ak_20170705+

       ! TRACER FAMILY ERROR
    CASE('2000')
       tracer_error_str = 'E'//echar//': TRACER MODE IS ALREADY ACTIVE'
    CASE('2001')
       tracer_error_str = 'E'//echar//': FAMILY MODE IS ALREADY ACTIVE'
    CASE('2010')
       tracer_error_str = 'E'//echar//': UNKOWN CONVERSION FLAG'

! op_pj_20170811+
       ! CHEMICAL PROPERTIES ERROR
    CASE('3001')
       tracer_error_str = 'E'//echar//': REFERENCE CHEMICAL SPECIES UNKNOWN'
! op_pj_20170811-

    CASE DEFAULT
       tracer_error_str = 'E'//echar//': UNKONW ERROR STATUS'
    END SELECT

  END FUNCTION tracer_error_str
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  FUNCTION PARAM2STRING(i, mode)

    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, TRIM

    ! I/O
    CHARACTER(LEN=STRLEN_MEDIUM)               :: PARAM2STRING
    INTEGER,                     INTENT(IN)    :: i
    CHARACTER(LEN=*),            INTENT(IN)    :: mode

    PARAM2STRING = ''

    SELECT CASE(TRIM(ADJUSTL(mode)))
       !
       CASE('switch')
          SELECT CASE(i)
          CASE(ON)
             PARAM2STRING = 'ON' 
          CASE(OFF)
             PARAM2STRING = 'OFF'
          CASE(MODAL)
             PARAM2STRING = 'MODAL'
          CASE(BIN)
             PARAM2STRING = 'BIN'
          CASE DEFAULT
             PARAM2STRING = 'UNDEFINED'
          END SELECT

       CASE('type')
          SELECT CASE(i)
          CASE(SINGLE)
             PARAM2STRING = 'SINGLE'
          CASE(FAMILY)
             PARAM2STRING = 'FAMILY'
          CASE(ISOTOPE)
             PARAM2STRING = 'ISOTOPE'
          CASE DEFAULT
             PARAM2STRING = 'UNDEFINED'
          END SELECT

       CASE('medium')
          SELECT CASE(i)
          CASE(AIR)
             PARAM2STRING = 'AIR'
          CASE(AEROSOL)
             PARAM2STRING = 'AEROSOL'
          CASE(CLOUD)
             PARAM2STRING = 'CLOUD'
          CASE(OCEAN)
             PARAM2STRING = 'OCEAN'
          CASE(LAKE)
             PARAM2STRING = 'LAKE'
          CASE(RIVER)
             PARAM2STRING = 'RIVER'
          CASE(LANDICE)
             PARAM2STRING = 'LANDICE'
          CASE(SEAICE)
             PARAM2STRING = 'SEAICE'
          CASE(VEGETATION)
             PARAM2STRING = 'VEGETATION'
          CASE DEFAULT
             PARAM2STRING = 'UNDEFINED'
          END SELECT

       CASE('quantity')
          SELECT CASE(i)
          CASE(AMOUNTFRACTION)
             PARAM2STRING = 'AMOUNTFRACTION'
          CASE(NUMBERDENSITY)
             PARAM2STRING = 'NUMBERDENSITY'
          CASE(CONCENTRATION)
             PARAM2STRING = 'CONCENTRATION'
          CASE DEFAULT
             PARAM2STRING = 'UNDEFINED'
          END SELECT

       CASE DEFAULT
          PARAM2STRING = 'UNDEFINED'
    END SELECT

  END FUNCTION PARAM2STRING
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE full2base_sub(status, fullname, basename, subname)

    IMPLICIT NONE

    INTRINSIC :: LEN, INDEX, TRIM

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    CHARACTER(LEN=*), INTENT(IN)  :: fullname
    CHARACTER(LEN=*), INTENT(OUT) :: basename
    CHARACTER(LEN=*), INTENT(OUT) :: subname

    ! LOCAL
    INTEGER :: si

    status = 206 ! ERROR: fullname too long
    IF (LEN(fullname) > STRLEN_FNAME) RETURN

    status = 200 ! ERROR: tracer basename too long
    IF (LEN(basename) > STRLEN_MEDIUM) RETURN

    status = 201 ! ERROR: tracer subname too long
    IF (LEN(basename) > STRLEN_MEDIUM) RETURN

    ! INIT
    basename = ''
    subname = ''

    si = INDEX(fullname, '_')
    IF (si == 0) THEN
       basename = TRIM(fullname)
       subname  = ''
    ELSE
       basename = TRIM(fullname(:si-1))
       subname  = TRIM(fullname(si+1:))
    END IF

    status = 0

  END SUBROUTINE full2base_sub
  ! -------------------------------------------------------------------

! op_pj_20150811+
  ! -------------------------------------------------------------------
  SUBROUTINE set_tracer_properties(extstatus, QTPROP, lprint, patch_id)

    USE messy_main_tools, ONLY: strcrack, ucase, get_name_domain, int2str

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT)                                    :: extstatus
    TYPE(T_TRACPROP_IO), DIMENSION(:), INTENT(IN), OPTIONAL :: QTPROP
    LOGICAL,                           INTENT(IN), OPTIONAL :: lprint
    INTEGER,                           INTENT(IN), OPTIONAL :: patch_id

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'set_tracer_properties'
    TYPE(T_TRACPROP_IO), DIMENSION(:), ALLOCATABLE :: ZTPROP
    INTEGER  :: nmax
    INTEGER  :: status
    INTEGER  :: i
    CHARACTER(LEN=2*STRLEN_MEDIUM+1), DIMENSION(:), POINTER :: trn => NULL()
    CHARACTER(LEN=2*STRLEN_TRSET+1),  DIMENSION(:), POINTER :: trs => NULL()
    INTEGER  :: n, j, m, k
    CHARACTER(LEN=STRLEN_MEDIUM)     :: basename    = '' ! name of tracer
    CHARACTER(LEN=STRLEN_MEDIUM)     :: subname     = '' ! OPTIONAL subname
    INTEGER  :: idx
    INTEGER  :: jc
    INTEGER  :: icont
    CHARACTER(LEN=STRLEN_MEDIUM)     :: scont
    REAL(DP) :: rcont
    LOGICAL  :: lfound
    LOGICAL  :: llprint
    ! op_pj_20170809+
    CHARACTER(LEN=STRLEN_TRSET) :: zts
    INTEGER                     :: ireqpatch
    CHARACTER(LEN=2)            :: pstr = '  '
    ! op_pj_20170809-
    INTEGER  :: ids ! species index ! op_pj_20170811

    extstatus = 0

    IF (PRESENT(QTPROP)) THEN
       nmax = SIZE(QTPROP)
       ALLOCATE(ZTPROP(nmax))
       ZTPROP(:) = QTPROP(:)
    ELSE
       nmax = SIZE(TPROP)
       ALLOCATE(ZTPROP(nmax))
       ZTPROP(:) = TPROP(:)
    END IF

    IF (PRESENT(lprint)) THEN
       llprint = lprint
    ELSE
       llprint = .TRUE.
    ENDIF

    nmlprop: DO i=1, nmax

       IF (TRIM(ZTPROP(i)%trset)    == '') CYCLE
       IF (TRIM(ZTPROP(i)%trlist)   == '') CYCLE
       IF (TRIM(ZTPROP(i)%caskname) == '') CYCLE
       IF (TRIM(ZTPROP(i)%cont)     == '') CYCLE

       ! op_pj_20160823+
       ! split tracer sets into individual names
       CALL strcrack(TRIM(ZTPROP(i)%trset), ';', trs, m)
       ! op_pj_20160823-

       ! split list of tracers into individual names
       CALL strcrack(TRIM(ZTPROP(i)%trlist), ';', trn, n)

       ntrac: DO j=1, n
          ! convert fullname into basename and subname
          CALL full2base_sub(status, TRIM(trn(j)), basename, subname)
          IF (status /= 0 .AND. llprint) &
               WRITE(*,*) substr//' (1): '//tracer_error_str(status)

          nset: DO k=1, m

             ! op_pj_20170809+
             ! SPECIAL HOOK TO TAKE INTO ACCOUNT THE PATCH_ID, IF ANY
             ! 1) check, if tracer_set is already specific with patch_id
             CALL get_name_domain(status, TRIM(trs(k)), zts, ireqpatch)
             IF (status == 0) THEN
                ! tracer set name is specific for patch ireqpatch
                IF (PRESENT(patch_id)) THEN
                   ! requested patch and actual patch mismatch
                   IF (ireqpatch /= patch_id) CYCLE
                ELSE
                   ! patch requested, but no patch available
                   CYCLE
                ENDIF
             ELSE
                ! tracer set name is NOT specific for patch
                IF (PRESENT(patch_id)) THEN
                   ! expand name for patch (implicit wildcard)
                   CALL int2str(pstr,patch_id)
                   zts = TRIM(zts)//'_D'//pstr
                !ELSE
                ! ! nothing to do: set name is not specific, no patch requested
                END IF
             ENDIF
             ! op_pj_20170809-

             ! get tracer ID
             CALL get_tracer(status, TRIM(zts), basename &
                  , subname=subname, idx=idx)
             IF (status /= 0) THEN
                IF (llprint) WRITE(*,*) substr//' (2): '//&
                     &tracer_error_str(status)//&
                     &' [ set/tracer = '//TRIM(zts)//'/'//TRIM(trn(j))//' ]'
                CYCLE
             END IF

             ! search for correct container ... and set contents ...
             lfound = .FALSE.

             ! ... integer
             DO jc = 1, MAX_CASK_I
                IF (TRIM(NAMES_CASK_I(jc)) == TRIM(ZTPROP(i)%caskname)) THEN

                   lfound = .TRUE.

                   ! SPECIAL CASES: ON, OFF, ...
                   scont = TRIM(ZTPROP(i)%cont)
                   CALL ucase(scont)
                   SELECT CASE(TRIM(scont))
                   CASE('ON')
                      icont = ON
                   CASE('OFF')
                      icont = OFF

                   CASE('MODAL')
                      icont = MODAL
                   CASE('BIN')
                      icont = BIN

                   CASE('T_INI_ZERO')
                      icont = T_INI_ZERO
                   CASE('T_INI_FILE')
                      icont = T_INI_FILE
                   CASE('T_INI_USER')
                      icont = T_INI_USER
                   CASE('T_LBC_ZERO')
                      icont = T_LBC_ZERO
                   CASE('T_LBC_FILE')
                      icont = T_LBC_FILE
                   CASE('T_LBC_CST')
                      icont = T_LBC_CST
                   CASE('T_LBC_ZEROGRAD')
                      icont = T_LBC_ZEROGRAD
                   CASE('T_LBC_USER')
                      icont = T_LBC_USER
                   CASE('T_BBC_ZEROFLUX')
                      icont = T_BBC_ZEROFLUX
                   CASE('T_BBC_ZEROVAL')
                      icont = T_BBC_ZEROVAL
                   CASE('T_BBC_SURF_VAL')
                      icont = T_BBC_SURF_VAL
                   CASE('T_RELAX_OFF')
                      icont = T_RELAX_OFF
                   CASE('T_RELAX_FULL')
                      icont = T_RELAX_FULL
                   CASE('T_RELAX_INFLOW')
                      icont = T_RELAX_INFLOW
                   CASE('T_ADV_OFF')
                      icont = T_ADV_OFF
                   CASE('T_ADV_ON')
                      icont = T_ADV_ON
                   CASE('T_ADV_2_LF')
                      icont = T_ADV_2_LF
                   CASE('T_ADV_3_LF')
                      icont = T_ADV_3_LF
                   CASE('T_TURB_OFF')
                      icont = T_TURB_OFF
                   CASE('T_TURB_1D')
                      icont = T_TURB_1D
                   CASE('T_TURB_3D')
                      icont = T_TURB_3D
                   CASE('T_DAMP_OFF')
                      icont = T_DAMP_OFF
                   CASE('T_DAMP_ON')
                      icont = T_DAMP_ON
                   CASE('T_DAMP_FORCED')
                      icont = T_DAMP_FORCED

                   CASE DEFAULT
                      READ(ZTPROP(i)%cont,*) icont
                   END SELECT
                   CALL set_tracer(status, TRIM(zts), idx, jc &
                        , icont)
                   IF (status /= 0 .AND. llprint) WRITE(*,*) substr//' (3): '//&
                        &tracer_error_str(status)//&
                        &' [ set/tracer = '//TRIM(zts)//'/',idx,' ]'
                END IF
             END DO

             ! ... string
             DO jc = 1, MAX_CASK_S
                IF (TRIM(NAMES_CASK_S(jc)) == TRIM(ZTPROP(i)%caskname)) THEN
                   lfound = .TRUE.
                   CALL set_tracer(status, TRIM(zts), idx, jc &
                        , TRIM(ZTPROP(i)%cont))
                   IF (status /= 0 .AND. llprint) WRITE(*,*) substr//' (4): '//&
                        &tracer_error_str(status)//&
                        &' [ set/tracer = '//TRIM(zts)//'/',idx,' ]'
                END IF
             END DO

             ! ... real
             DO jc = 1, MAX_CASK_R
                IF (TRIM(NAMES_CASK_R(jc)) == TRIM(ZTPROP(i)%caskname)) THEN
                   lfound = .TRUE.
                   READ(ZTPROP(i)%cont,*) rcont
                   CALL set_tracer(status, TRIM(zts), idx, jc &
                        , rcont)
                   IF (status /= 0 .AND. llprint) WRITE(*,*) substr//' (5): '//&
                        &tracer_error_str(status)//&
                        &' [ set/tracer = '//TRIM(zts)//'/',idx,' ]'
                END IF
             END DO

             IF (.NOT. lfound) THEN
! op_pj_20170811+
!!$             extstatus = 703
!!$             RETURN
                ! add special additional case(s)
                SELECT CASE(TRIM(ZTPROP(i)%caskname))
                CASE('refspec')
                   ids = get_chemprop_index(TRIM(ZTPROP(i)%cont))
                   IF (ids > 0) THEN
                      ! ... integer
                      DO jc = 1, MAX_CASK_I_CHEMPROP
                         CALL set_tracer(status, TRIM(zts), idx, jc, &
                              chemprop(ids)%cask_i(jc) )
                         IF (status /= 0 .AND. llprint) &
                              WRITE(*,*) substr//' (6): '//&
                              &tracer_error_str(status)//&
                              &' [ set/tracer = '//TRIM(zts)//'/',idx,' ]'
                      END DO
                      ! ... string
                      DO jc = 1, MAX_CASK_S_CHEMPROP
                         CALL set_tracer(status, TRIM(zts), idx, jc, &
                              chemprop(ids)%cask_s(jc) )
                         IF (status /= 0 .AND. llprint) &
                              WRITE(*,*) substr//' (7): '//&
                              &tracer_error_str(status)//&
                              &' [ set/tracer = '//TRIM(zts)//'/',idx,' ]'
                      END DO
                      ! ... real
                      DO jc = 1, MAX_CASK_R_CHEMPROP
                         CALL set_tracer(status, TRIM(zts), idx, jc, &
                              chemprop(ids)%cask_r(jc) )
                         IF (status /= 0 .AND. llprint) &
                              WRITE(*,*) substr//' (8): '//&
                              &tracer_error_str(status)//&
                              &' [ set/tracer = '//TRIM(zts)//'/',idx,' ]'
                      END DO
                   ELSE
                      ! REFERENCE CHEMICAL SPECIES UNKNOWN
                      status = 3001
                      RETURN
                   END IF
                CASE DEFAULT
                   extstatus = 703
                   RETURN
                END SELECT
! op_pj_20170811-
             END IF

          END DO nset

       END DO ntrac

       IF (ASSOCIATED(trn)) THEN
          DEALLOCATE(trn)
          NULLIFY(trn)
       END IF

       IF (ASSOCIATED(trs)) THEN
          DEALLOCATE(trs)
          NULLIFY(trs)
       END IF

    END DO nmlprop

    IF (ALLOCATED(ZTPROP)) DEALLOCATE(ZTPROP)

  END SUBROUTINE set_tracer_properties
  ! -------------------------------------------------------------------
! op_pj_20150811-

  ! -------------------------------------------------------------------
  SUBROUTINE main_tracer_read_nml_ctrl(status, iou)

    ! TRACER MODULE ROUTINE (CORE)
    !
    ! READ TRACER NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Jul 2003

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ l_family, l_pdef, TPROP ! op_pj_20150811 added TPROP

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    IF (l_family) THEN
       WRITE(*,*) '  TRACER FAMILIES              : ON'
    ELSE
       WRITE(*,*) '  TRACER FAMILIES              : OFF'
    END IF

    IF (l_pdef) THEN
       WRITE(*,*) '  TRACER PDEF                  : ON'
    ELSE
       WRITE(*,*) '  TRACER PDEF                  : OFF'
    END IF

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE main_tracer_read_nml_ctrl
  ! -------------------------------------------------------------------
  
  ! mz_rs_20170209+
  INTEGER FUNCTION get_chemprop_index(name)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: name
    INTEGER :: i
    get_chemprop_index = 0 ! dummy value for non-existing index
    DO i = 1, N_CHEMPROP
      IF (TRIM(chemprop(i)%kppname)==name) THEN
        get_chemprop_index = i
      ENDIF
    ENDDO
  END FUNCTION get_chemprop_index
  ! mz_rs_20170209-

  ! mz_rs_20160521+
  INTEGER FUNCTION get_chemprop_real(name, container, value)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: name
    INTEGER,          INTENT(IN)  :: container
    REAL(DP),         INTENT(OUT) :: value
    INTEGER :: i
    get_chemprop_real = 1 ! set status to error until species "name" is found
    DO i = 1, N_CHEMPROP
      IF (TRIM(chemprop(i)%kppname)==name) THEN
         value = chemprop(i)%cask_r(container)
         get_chemprop_real = 0 ! status = okay
      ENDIF
    ENDDO
  END FUNCTION get_chemprop_real

  INTEGER FUNCTION get_chemprop_int(name, container, value)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: name
    INTEGER,          INTENT(IN)  :: container
    INTEGER,          INTENT(OUT) :: value
    INTEGER :: i
    get_chemprop_int = 1 ! set status to error until species "name" is found
    DO i = 1, N_CHEMPROP
      IF (TRIM(chemprop(i)%kppname)==name) THEN
         value = chemprop(i)%cask_i(container)
         get_chemprop_int = 0 ! status = okay
      ENDIF
    ENDDO
  END FUNCTION get_chemprop_int

  INTEGER FUNCTION get_chemprop_char(name, container, value)
    IMPLICIT NONE
    CHARACTER(LEN=*),  INTENT(IN)  :: name
    INTEGER,           INTENT(IN)  :: container
    CHARACTER(LEN=24), INTENT(OUT) :: value
    INTEGER :: i
    get_chemprop_char = 1 ! set status to error until species "name" is found
    DO i = 1, N_CHEMPROP
      IF (TRIM(chemprop(i)%kppname)==name) THEN
         value = chemprop(i)%cask_s(container)
         get_chemprop_char = 0 ! status = okay
      ENDIF
    ENDDO
  END FUNCTION get_chemprop_char
  ! mz_rs_20160521-

! **********************************************************************
END MODULE messy_main_tracer
! **********************************************************************
