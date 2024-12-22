! -*- f90 -*- Created automatically by messy_main_tracer_chemprop.awk, DO NOT EDIT!

INTEGER, PUBLIC, PARAMETER :: S_casrn = 1
INTEGER, PUBLIC, PARAMETER :: R_molarmass = 1
INTEGER, PUBLIC, PARAMETER :: R_Henry_T0 = 2
INTEGER, PUBLIC, PARAMETER :: R_Henry_Tdep = 3
INTEGER, PUBLIC, PARAMETER :: R_Sechenov = 4
INTEGER, PUBLIC, PARAMETER :: R_alpha_T0 = 5
INTEGER, PUBLIC, PARAMETER :: R_alpha_Tdep = 6
INTEGER, PUBLIC, PARAMETER :: R_K_acid = 7
INTEGER, PUBLIC, PARAMETER :: R_K_acid2 = 8
INTEGER, PUBLIC, PARAMETER :: R_Vmol_bp = 9
INTEGER, PUBLIC, PARAMETER :: R_psat = 10
INTEGER, PUBLIC, PARAMETER :: R_psat_Tdep = 11
INTEGER, PUBLIC, PARAMETER :: R_pss = 12
INTEGER, PUBLIC, PARAMETER :: R_dryreac_sf = 13
INTEGER, PUBLIC, PARAMETER :: I_charge = 1
INTEGER, PUBLIC, PARAMETER :: R_aerosol_density = 14
INTEGER, PUBLIC, PARAMETER :: MAX_CASK_I_CHEMPROP = 1
INTEGER, PUBLIC, PARAMETER :: MAX_CASK_R_CHEMPROP = 14
INTEGER, PUBLIC, PARAMETER :: MAX_CASK_S_CHEMPROP = 1
CHARACTER(LEN=NAMES_CASK_STRLEN), DIMENSION(MAX_CASK_I_CHEMPROP), &
  PARAMETER, PUBLIC :: NAMES_CASK_I_CHEMPROP = (/ &
  'charge         ' /)
CHARACTER(LEN=NAMES_CASK_STRLEN), DIMENSION(MAX_CASK_R_CHEMPROP), &
  PARAMETER, PUBLIC :: NAMES_CASK_R_CHEMPROP = (/ &
  'molarmass      ', 'Henry_T0       ', 'Henry_Tdep     ', &
  'Sechenov       ', 'alpha_T0       ', 'alpha_Tdep     ', &
  'K_acid         ', 'K_acid2        ', 'Vmol_bp        ', &
  'psat           ', 'psat_Tdep      ', 'pss            ', &
  'dryreac_sf     ', 'aerosol_density' /)
CHARACTER(LEN=NAMES_CASK_STRLEN), DIMENSION(MAX_CASK_S_CHEMPROP), &
  PARAMETER, PUBLIC :: NAMES_CASK_S_CHEMPROP = (/ &
  'casrn          ' /)
TYPE, PUBLIC :: t_meta_chemprop
  CHARACTER(LEN=20) :: kppname
  INTEGER, DIMENSION(MAX_CASK_I_CHEMPROP) :: cask_i
  REAL(DP), DIMENSION(MAX_CASK_R_CHEMPROP) :: cask_r
  CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(MAX_CASK_S_CHEMPROP) :: cask_s
END TYPE t_meta_chemprop

TYPE(t_meta_chemprop), PARAMETER :: spc_001 = t_meta_chemprop( &
  'O1D', &
  (/ 0 /), & ! integer
  (/ MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_002 = t_meta_chemprop( &
  'O3P', &
  (/ 0 /), & ! integer
  (/ MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_003 = t_meta_chemprop( &
  'O2', &
  (/ 0 /), & ! integer
  (/ MO*2_DP, 1.3E-3_DP, 1500._DP, 9.99999_DP, 0.01_DP, 2000._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.0_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '7782-44-7               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_004 = t_meta_chemprop( &
  'O3', &
  (/ 0 /), & ! integer
  (/ MO*3_DP, 1.03E-2_DP, 2830._DP, 9.99999_DP, 0.002_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.01_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '10028-15-6              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_005 = t_meta_chemprop( &
  'H', &
  (/ 0 /), & ! integer
  (/ MH, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '12385-13-6              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_006 = t_meta_chemprop( &
  'H2', &
  (/ 0 /), & ! integer
  (/ MH*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '1333-74-0               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_007 = t_meta_chemprop( &
  'OH', &
  (/ 0 /), & ! integer
  (/ MO+MH, 3.0E1_DP, 4300._DP, 9.99999_DP, 0.01_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 25._DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '3352-57-6               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_008 = t_meta_chemprop( &
  'HO2', &
  (/ 0 /), & ! integer
  (/ MH+MO*2_DP, 3.9E3_DP, 5900._DP, 9.99999_DP, 0.5_DP, 0._DP, 3.5E-5_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.37E6_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '3170-83-0               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_009 = t_meta_chemprop( &
  'H2O', &
  (/ 0 /), & ! integer
  (/ MH*2_DP+MO, BIG_DP, 0._DP, 9.99999_DP, 0.0_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.0_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '7732-18-5               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_010 = t_meta_chemprop( &
  'H2O2', &
  (/ 0 /), & ! integer
  (/ MH*2_DP+MO*2_DP, 1.E5_DP, 6338._DP, 9.99999_DP, 0.077_DP, 3127._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 7.45E4_DP, &
  1.0_DP, 9.99999_DP /), & ! real
  (/ '7722-84-1               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_011 = t_meta_chemprop( &
  'H2OH2O', &
  (/ 0 /), & ! integer
  (/ MH*2_DP+MO+MH*2_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '25655-83-8              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_012 = t_meta_chemprop( &
  'N', &
  (/ 0 /), & ! integer
  (/ MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_013 = t_meta_chemprop( &
  'N2D', &
  (/ 0 /), & ! integer
  (/ MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_014 = t_meta_chemprop( &
  'N2', &
  (/ 0 /), & ! integer
  (/ MN*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_015 = t_meta_chemprop( &
  'NH3', &
  (/ 0 /), & ! integer
  (/ MN+MH*3_DP, 60.2_DP, 4160._DP, 9.99999_DP, 0.06_DP, 0._DP, -1.75E-5_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.02E4_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '7664-41-7               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_016 = t_meta_chemprop( &
  'N2O', &
  (/ 0 /), & ! integer
  (/ MN*2_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_017 = t_meta_chemprop( &
  'NO', &
  (/ 0 /), & ! integer
  (/ MN+MO, 1.9E-3_DP, 1480._DP, 9.99999_DP, 5.0E-5_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.E-3_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '10102-43-9              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_018 = t_meta_chemprop( &
  'NO2', &
  (/ 0 /), & ! integer
  (/ MN+MO*2_DP, 1.2E-2_DP, 2360._DP, 9.99999_DP, 0.0015_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E-2_DP, &
  1.0_DP, 9.99999_DP /), & ! real
  (/ '10102-44-0              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_019 = t_meta_chemprop( &
  'NO3', &
  (/ 0 /), & ! integer
  (/ MN+MO*3_DP, 3.8E-2_DP, 2000._DP, 9.99999_DP, 0.04_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.8_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '12033-49-7              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_020 = t_meta_chemprop( &
  'N2O5', &
  (/ 0 /), & ! integer
  (/ MN*2_DP+MO*5_DP, BIG_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.E4_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '10102-03-1              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_021 = t_meta_chemprop( &
  'HONO', &
  (/ 0 /), & ! integer
  (/ MH+MO+MN+MO, 4.9E1_DP, 4780._DP, 9.99999_DP, 0.04_DP, 0._DP, 5.1E-4_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 5.05E4_DP, 0.1_DP, &
  9.99999_DP /), & ! real
  (/ '7782-77-6               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_022 = t_meta_chemprop( &
  'HNO3', &
  (/ 0 /), & ! integer
  (/ MH+MN+MO*3_DP, 2.45E6_DP/1.5E1_DP, 8694._DP, 9.99999_DP, 0.5_DP, 0._DP, &
  15.4_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.E4_DP, &
  1.0_DP, 9.99999_DP /), & ! real
  (/ '7697-37-2               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_023 = t_meta_chemprop( &
  'HNO4', &
  (/ 0 /), & ! integer
  (/ MH+MN+MO*4_DP, 1.26E4_DP, 6900._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.E4_DP, &
  1.0_DP, 9.99999_DP /), & ! real
  (/ '26404-66-0              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_024 = t_meta_chemprop( &
  'NH2', &
  (/ 0 /), & ! integer
  (/ MH*2_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, 0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_025 = t_meta_chemprop( &
  'HNO', &
  (/ 0 /), & ! integer
  (/ MH+MN+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, 0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_026 = t_meta_chemprop( &
  'NHOH', &
  (/ 0 /), & ! integer
  (/ MH*2_DP+MN+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, 0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_027 = t_meta_chemprop( &
  'NH2O', &
  (/ 0 /), & ! integer
  (/ MH*2_DP+MN+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, 0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_028 = t_meta_chemprop( &
  'NH2OH', &
  (/ 0 /), & ! integer
  (/ MH*3_DP+MN+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, 0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_029 = t_meta_chemprop( &
  'LNITROGEN', &
  (/ 0 /), & ! integer
  (/ MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_030 = t_meta_chemprop( &
  'CH2OO', &
  (/ 0 /), & ! integer
  (/ MC+MH*2_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_031 = t_meta_chemprop( &
  'CH2OOA', &
  (/ 0 /), & ! integer
  (/ MC+MH*2_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_032 = t_meta_chemprop( &
  'CH3', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '2229-07-4               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_033 = t_meta_chemprop( &
  'CH3O', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '2143-68-2               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_034 = t_meta_chemprop( &
  'CH3O2', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MO*2_DP, 6._DP, 5600._DP, 9.99999_DP, 0.01_DP, 2000._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 6._DP, &
  1.0_DP, 9.99999_DP /), & ! real
  (/ '2143-58-0               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_035 = t_meta_chemprop( &
  'CH3OH', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MO+MH, 2.20E2_DP, 5200._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.2E2_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '67-56-1                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_036 = t_meta_chemprop( &
  'CH3OOH', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MO+MO+MH, 3.0E2_DP, 5322._DP, 9.99999_DP, 0.0046_DP, 3273._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 3.E2_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '3031-73-0               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_037 = t_meta_chemprop( &
  'CH4', &
  (/ 0 /), & ! integer
  (/ MC+MH*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '74-92-8                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_038 = t_meta_chemprop( &
  'CO', &
  (/ 0 /), & ! integer
  (/ MC+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '630-08-0                ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_039 = t_meta_chemprop( &
  'CO2', &
  (/ 0 /), & ! integer
  (/ MC+MO*2_DP, 3.4E-2_DP, 2400._DP, 9.99999_DP, 0.01_DP, 2000._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.0_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '124-38-9                ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_040 = t_meta_chemprop( &
  'HCHO', &
  (/ 0 /), & ! integer
  (/ MH+MC+MH+MO, 3.23E3_DP, 7100._DP, 9.99999_DP, 0.04_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 3.2E3_DP, 0.1_DP, &
  9.99999_DP /), & ! real
  (/ '50-00-0                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_041 = t_meta_chemprop( &
  'HCOOH', &
  (/ 0 /), & ! integer
  (/ MH+MC+MO+MO+MH, 8.9E3_DP, 6100._DP, 9.99999_DP, 0.014_DP, 3978._DP, &
  1.78E-4_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 6.7E6_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '64-18-6                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_042 = t_meta_chemprop( &
  'HOCH2O2', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 4.6E4_DP, &
  1.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_043 = t_meta_chemprop( &
  'HOCH2OH', &
  (/ 0 /), & ! integer
  (/ MC+MH*4_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.7E6_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_044 = t_meta_chemprop( &
  'HOCH2OOH', &
  (/ 0 /), & ! integer
  (/ MC+MH*4_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.7E6_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_045 = t_meta_chemprop( &
  'CH3NO3', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_046 = t_meta_chemprop( &
  'CH3O2NO2', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_047 = t_meta_chemprop( &
  'CH3ONO', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MO*2_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_048 = t_meta_chemprop( &
  'CN', &
  (/ 0 /), & ! integer
  (/ MC+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_049 = t_meta_chemprop( &
  'HCN', &
  (/ 0 /), & ! integer
  (/ MC+MH+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 7.6_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_050 = t_meta_chemprop( &
  'HOCH2O2NO2', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_051 = t_meta_chemprop( &
  'NCO', &
  (/ 0 /), & ! integer
  (/ MC+MN+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_052 = t_meta_chemprop( &
  'LCARBON', &
  (/ 0 /), & ! integer
  (/ MC, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_053 = t_meta_chemprop( &
  'C2H2', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_054 = t_meta_chemprop( &
  'C2H4', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_055 = t_meta_chemprop( &
  'C2H5O2', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*5_DP+MO*2_DP, 6._DP, 5600._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.0_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '3170-61-4               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_056 = t_meta_chemprop( &
  'C2H5OH', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*6_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.0E2_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_057 = t_meta_chemprop( &
  'C2H5OOH', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*6_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.4E2_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_058 = t_meta_chemprop( &
  'C2H6', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_059 = t_meta_chemprop( &
  'CH2CHOH', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*4_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.0E2_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_060 = t_meta_chemprop( &
  'CH2CO', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*2_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E6_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_061 = t_meta_chemprop( &
  'CH3CHO', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MC+MH+MO, 1.29E1_DP, 5890._DP, 9.99999_DP, 3.0E-2_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.3E1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '75-07-0                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_062 = t_meta_chemprop( &
  'CH3CHOHO2', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*5_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_063 = t_meta_chemprop( &
  'CH3CHOHOOH', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_064 = t_meta_chemprop( &
  'CH3CO', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_065 = t_meta_chemprop( &
  'CH3CO2H', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MC+MO*2_DP+MH, 4.1E3_DP, 6200._DP, 9.99999_DP, 2.0E-2_DP, &
  4079._DP, 1.75E-5_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  7.4E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '64-19-7                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_066 = t_meta_chemprop( &
  'CH3CO3', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_067 = t_meta_chemprop( &
  'CH3CO3H', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*4_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  8.4E2_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_068 = t_meta_chemprop( &
  'ETHGLY', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*6_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_069 = t_meta_chemprop( &
  'GLYOX', &
  (/ 0 /), & ! integer
  (/ MC+MH+MO+MC+MH+MO, 4.19E5_DP, 7480._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 4.2E5_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '107-22-2                ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_070 = t_meta_chemprop( &
  'HCOCH2O2', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_071 = t_meta_chemprop( &
  'HCOCO', &
  (/ 0 /), & ! integer
  (/ MH+MC+MO+MC+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_072 = t_meta_chemprop( &
  'HCOCO2H', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*2_DP+MO*3_DP, 1.09E4_DP, 4810._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  6.61E-4_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 6.6E7_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '298-12-4                ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_073 = t_meta_chemprop( &
  'HCOCO3', &
  (/ 0 /), & ! integer
  (/ MH+MC+MO+MC+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_074 = t_meta_chemprop( &
  'HCOCO3H', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*2_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.7E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_075 = t_meta_chemprop( &
  'HOCH2CH2O', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*5_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_076 = t_meta_chemprop( &
  'HOCH2CH2O2', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*5_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_077 = t_meta_chemprop( &
  'HOCH2CHO', &
  (/ 0 /), & ! integer
  (/ MH+MO+MC+MH*2_DP+MC+MH+MO, 4.E4_DP, 4630._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 4.10E4_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '141-46-8                ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_078 = t_meta_chemprop( &
  'HOCH2CO', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_079 = t_meta_chemprop( &
  'HOCH2CO2H', &
  (/ 0 /), & ! integer
  (/ MH+MO+MC+MH*2_DP+MC+MO*2_DP+MH, 2.83E4_DP, 4030._DP, 9.99999_DP, 0.1_DP, &
  0._DP, 1.48E-4_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.2E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_080 = t_meta_chemprop( &
  'HOCH2CO3', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_081 = t_meta_chemprop( &
  'HOCH2CO3H', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.2E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_082 = t_meta_chemprop( &
  'HOCHCHO', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_083 = t_meta_chemprop( &
  'HOOCH2CHO', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*4_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_084 = t_meta_chemprop( &
  'HOOCH2CO2H', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.2E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_085 = t_meta_chemprop( &
  'HOOCH2CO3', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_086 = t_meta_chemprop( &
  'HOOCH2CO3H', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*4_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_087 = t_meta_chemprop( &
  'HYETHO2H', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_088 = t_meta_chemprop( &
  'C2H5NO3', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*5_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_089 = t_meta_chemprop( &
  'C2H5O2NO2', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*5_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_090 = t_meta_chemprop( &
  'CH3CN', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 5.3_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_091 = t_meta_chemprop( &
  'ETHOHNO3', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*5_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.9E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_092 = t_meta_chemprop( &
  'NCCH2O2', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*2_DP+MO*2_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_093 = t_meta_chemprop( &
  'NO3CH2CHO', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_094 = t_meta_chemprop( &
  'NO3CH2CO3', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*2_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_095 = t_meta_chemprop( &
  'NO3CH2PAN', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MO*8_DP+MN*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 2.8_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_096 = t_meta_chemprop( &
  'PAN', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MO*5_DP+MN, 2.8_DP, 5730._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.8_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '2278-22-0               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_097 = t_meta_chemprop( &
  'PHAN', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_098 = t_meta_chemprop( &
  'ACETOL', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*6_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.7E2_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_099 = t_meta_chemprop( &
  'ALCOCH2OOH', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_100 = t_meta_chemprop( &
  'C2H5CHO', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*6_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.3E1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_101 = t_meta_chemprop( &
  'PROPACID', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*6_DP+MO*2_DP, 5.7E3_DP, 6800._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  1.35E-5_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 3.E5_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '79-09-4                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_102 = t_meta_chemprop( &
  'C2H5CO3', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*5_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_103 = t_meta_chemprop( &
  'PERPROACID', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.5E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_104 = t_meta_chemprop( &
  'C33CO', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*2_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_105 = t_meta_chemprop( &
  'C3H6', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_106 = t_meta_chemprop( &
  'C3H8', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*8_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_107 = t_meta_chemprop( &
  'CH3CHCO', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*4_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E6_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_108 = t_meta_chemprop( &
  'CH3COCH2O2', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*5_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_109 = t_meta_chemprop( &
  'CH3COCH3', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MC+MO+MC+MH*3_DP, 27.8_DP, 5530._DP, 9.99999_DP, 3.72E-3_DP, &
  6395._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.E1_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '67-64-1                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_110 = t_meta_chemprop( &
  'CH3COCO2H', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*4_DP+MO*3_DP, 3.11E5_DP, 5090._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  4.07E-3_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 4.3E8_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '127-17-3                ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_111 = t_meta_chemprop( &
  'CH3COCO3', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*3_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_112 = t_meta_chemprop( &
  'CH3COCO3H', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_113 = t_meta_chemprop( &
  'CHOCOCH2O2', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*3_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_114 = t_meta_chemprop( &
  'HCOCH2CHO', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*4_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_115 = t_meta_chemprop( &
  'HCOCH2CO2H', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.6E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_116 = t_meta_chemprop( &
  'HCOCH2CO3', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*3_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_117 = t_meta_chemprop( &
  'HCOCH2CO3H', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*4_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_118 = t_meta_chemprop( &
  'HCOCOCH2OOH', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_119 = t_meta_chemprop( &
  'HOC2H4CO2H', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.2E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_120 = t_meta_chemprop( &
  'HOC2H4CO3', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*5_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_121 = t_meta_chemprop( &
  'HOC2H4CO3H', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_122 = t_meta_chemprop( &
  'HOCH2COCH2O2', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*5_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_123 = t_meta_chemprop( &
  'HOCH2COCH2OOH', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_124 = t_meta_chemprop( &
  'HOCH2COCHO', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*4_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.1E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_125 = t_meta_chemprop( &
  'HYPERACET', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.7E2_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_126 = t_meta_chemprop( &
  'HYPROPO2', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*7_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_127 = t_meta_chemprop( &
  'HYPROPO2H', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.2E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_128 = t_meta_chemprop( &
  'IC3H7O2', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*7_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_129 = t_meta_chemprop( &
  'IC3H7OOH', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*8_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.3E2_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_130 = t_meta_chemprop( &
  'IPROPOL', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*8_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.3E2_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_131 = t_meta_chemprop( &
  'MGLYOX', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MC+MO+MC+MH+MO, 3.5E3_DP, 7500._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 3.7E3_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '78-98-8                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_132 = t_meta_chemprop( &
  'NC3H7O2', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*7_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_133 = t_meta_chemprop( &
  'NC3H7OOH', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*8_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.3E2_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_134 = t_meta_chemprop( &
  'NPROPOL', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*8_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.3E2_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_135 = t_meta_chemprop( &
  'PROPENOL', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*6_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.3E1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_136 = t_meta_chemprop( &
  'C32OH13CO', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*4_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_137 = t_meta_chemprop( &
  'C3DIALO2', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*3_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_138 = t_meta_chemprop( &
  'C3DIALOOH', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_139 = t_meta_chemprop( &
  'HCOCOHCO3', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*3_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_140 = t_meta_chemprop( &
  'HCOCOHCO3H', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*4_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_141 = t_meta_chemprop( &
  'METACETHO', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*4_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.7E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_142 = t_meta_chemprop( &
  'C3PAN1', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*5_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_143 = t_meta_chemprop( &
  'C3PAN2', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*3_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_144 = t_meta_chemprop( &
  'CH3COCH2O2NO2', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*5_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_145 = t_meta_chemprop( &
  'IC3H7NO3', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*7_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  0.83E0_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_146 = t_meta_chemprop( &
  'NC3H7NO3', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*7_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  0.83E0_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_147 = t_meta_chemprop( &
  'NOA', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*5_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_148 = t_meta_chemprop( &
  'PPN', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*5_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.8_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_149 = t_meta_chemprop( &
  'PR2O2HNO3', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*7_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.1E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_150 = t_meta_chemprop( &
  'PRONO3BO2', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*6_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_151 = t_meta_chemprop( &
  'PROPOLNO3', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*7_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.5E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_152 = t_meta_chemprop( &
  'HCOCOHPAN', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*3_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.0E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_153 = t_meta_chemprop( &
  'BIACET', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  7.4E1_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_154 = t_meta_chemprop( &
  'BIACETO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_155 = t_meta_chemprop( &
  'BIACETOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.3E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_156 = t_meta_chemprop( &
  'BIACETOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_157 = t_meta_chemprop( &
  'BUT1ENE', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_158 = t_meta_chemprop( &
  'BUT2OLO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_159 = t_meta_chemprop( &
  'BUT2OLO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*9_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_160 = t_meta_chemprop( &
  'BUT2OLOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*10_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_161 = t_meta_chemprop( &
  'BUTENOL', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 8.9_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_162 = t_meta_chemprop( &
  'C312COCO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*3_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_163 = t_meta_chemprop( &
  'C312COCO3H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_164 = t_meta_chemprop( &
  'C3H7CHO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 8.9_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_165 = t_meta_chemprop( &
  'C413COOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_166 = t_meta_chemprop( &
  'C44O2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_167 = t_meta_chemprop( &
  'C44OOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_168 = t_meta_chemprop( &
  'C4CODIAL', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_169 = t_meta_chemprop( &
  'CBUT2ENE', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_170 = t_meta_chemprop( &
  'CH3COCHCO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_171 = t_meta_chemprop( &
  'CH3COCHO2CHO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_172 = t_meta_chemprop( &
  'CH3COCOCO2H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.3E8_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_173 = t_meta_chemprop( &
  'CH3COOHCHCHO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_174 = t_meta_chemprop( &
  'CHOC3COO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_175 = t_meta_chemprop( &
  'CO23C3CHO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.6E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_176 = t_meta_chemprop( &
  'CO2C3CHO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.7E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_177 = t_meta_chemprop( &
  'CO2H3CHO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.1E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_178 = t_meta_chemprop( &
  'CO2H3CO2H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_179 = t_meta_chemprop( &
  'CO2H3CO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_180 = t_meta_chemprop( &
  'CO2H3CO3H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_181 = t_meta_chemprop( &
  'EZCH3CO2CHCHO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_182 = t_meta_chemprop( &
  'EZCHOCCH3CHO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_183 = t_meta_chemprop( &
  'HCOCCH3CHOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_184 = t_meta_chemprop( &
  'HCOCCH3CO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_185 = t_meta_chemprop( &
  'HCOCO2CH3CHO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_186 = t_meta_chemprop( &
  'HMAC', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.7E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_187 = t_meta_chemprop( &
  'HO12CO3C4', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  5.E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_188 = t_meta_chemprop( &
  'HVMK', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.7E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_189 = t_meta_chemprop( &
  'IBUTALOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_190 = t_meta_chemprop( &
  'IBUTDIAL', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.7E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_191 = t_meta_chemprop( &
  'IBUTOLBO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*9_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_192 = t_meta_chemprop( &
  'IBUTOLBOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*10_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_193 = t_meta_chemprop( &
  'IC4H10', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*10_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_194 = t_meta_chemprop( &
  'IC4H9O2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*9_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_195 = t_meta_chemprop( &
  'IC4H9OOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.1E2_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_196 = t_meta_chemprop( &
  'IPRCHO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 5.9E-1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_197 = t_meta_chemprop( &
  'IPRCO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_198 = t_meta_chemprop( &
  'IPRHOCO2H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.2E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_199 = t_meta_chemprop( &
  'IPRHOCO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_200 = t_meta_chemprop( &
  'IPRHOCO3H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_201 = t_meta_chemprop( &
  'MACO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_202 = t_meta_chemprop( &
  'MACO2H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*6_DP, 2.58E3_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  2.24E-5_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 5.8E5_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '79-41-4                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_203 = t_meta_chemprop( &
  'MACO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_204 = t_meta_chemprop( &
  'MACO3H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.4E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_205 = t_meta_chemprop( &
  'MACR', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 5.0E1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_206 = t_meta_chemprop( &
  'MACRO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_207 = t_meta_chemprop( &
  'MACRO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_208 = t_meta_chemprop( &
  'MACROH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  5.E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_209 = t_meta_chemprop( &
  'MACROOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  5.E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_210 = t_meta_chemprop( &
  'MBOOO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_211 = t_meta_chemprop( &
  'MEK', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.E1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_212 = t_meta_chemprop( &
  'MEPROPENE', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_213 = t_meta_chemprop( &
  'MPROPENOL', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 5.9E-1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_214 = t_meta_chemprop( &
  'MVK', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.3E1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_215 = t_meta_chemprop( &
  'NC4H10', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*10_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_216 = t_meta_chemprop( &
  'PERIBUACID', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.7E2_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_217 = t_meta_chemprop( &
  'TBUT2ENE', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_218 = t_meta_chemprop( &
  'TC4H9O2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*9_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_219 = t_meta_chemprop( &
  'TC4H9OOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  7.0E1_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_220 = t_meta_chemprop( &
  'BZFUCO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_221 = t_meta_chemprop( &
  'BZFUO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_222 = t_meta_chemprop( &
  'BZFUONE', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  36._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_223 = t_meta_chemprop( &
  'BZFUOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_224 = t_meta_chemprop( &
  'CO14O3CHO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.6E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_225 = t_meta_chemprop( &
  'CO14O3CO2H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.6E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_226 = t_meta_chemprop( &
  'CO2C4DIAL', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*2_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_227 = t_meta_chemprop( &
  'EPXC4DIAL', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.6E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_228 = t_meta_chemprop( &
  'EPXDLCO2H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.6E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_229 = t_meta_chemprop( &
  'EPXDLCO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*3_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_230 = t_meta_chemprop( &
  'EPXDLCO3H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_231 = t_meta_chemprop( &
  'HOCOC4DIAL', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.1E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_232 = t_meta_chemprop( &
  'MALANHY', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*2_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  36._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_233 = t_meta_chemprop( &
  'MALANHYO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*3_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_234 = t_meta_chemprop( &
  'MALANHYOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_235 = t_meta_chemprop( &
  'MALDALCO2H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.6E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_236 = t_meta_chemprop( &
  'MALDALCO3H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_237 = t_meta_chemprop( &
  'MALDIAL', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.6E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_238 = t_meta_chemprop( &
  'MALDIALCO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*3_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_239 = t_meta_chemprop( &
  'MALDIALO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_240 = t_meta_chemprop( &
  'MALDIALOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_241 = t_meta_chemprop( &
  'MALNHYOHCO', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*2_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_242 = t_meta_chemprop( &
  'MECOACEOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.1E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_243 = t_meta_chemprop( &
  'MECOACETO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_244 = t_meta_chemprop( &
  'BUT2OLNO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*9_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  8.8E1_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_245 = t_meta_chemprop( &
  'C312COPAN', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*3_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_246 = t_meta_chemprop( &
  'C4PAN5', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_247 = t_meta_chemprop( &
  'IBUTOLBNO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*9_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  8.8E1_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_248 = t_meta_chemprop( &
  'IC4H9NO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*9_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.4E-3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_249 = t_meta_chemprop( &
  'MACRNO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_250 = t_meta_chemprop( &
  'MPAN', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_251 = t_meta_chemprop( &
  'MVKNO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_252 = t_meta_chemprop( &
  'PIPN', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.9_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_253 = t_meta_chemprop( &
  'TC4H9NO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*9_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.4E-3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_254 = t_meta_chemprop( &
  'EPXDLPAN', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*3_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.8_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_255 = t_meta_chemprop( &
  'MALDIALPAN', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*3_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.8_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_256 = t_meta_chemprop( &
  'NBZFUO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*4_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_257 = t_meta_chemprop( &
  'NBZFUONE', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*3_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  20._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_258 = t_meta_chemprop( &
  'NBZFUOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*5_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.4E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_259 = t_meta_chemprop( &
  'NC4DCO2H', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*3_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.6E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_260 = t_meta_chemprop( &
  'LBUT1ENO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*9_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_261 = t_meta_chemprop( &
  'LBUT1ENOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*10_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_262 = t_meta_chemprop( &
  'LC4H9O2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*9_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_263 = t_meta_chemprop( &
  'LC4H9OOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.1E2_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_264 = t_meta_chemprop( &
  'LHMVKABO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_265 = t_meta_chemprop( &
  'LHMVKABOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  5.E6_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_266 = t_meta_chemprop( &
  'LMEKO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_267 = t_meta_chemprop( &
  'LMEKOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_268 = t_meta_chemprop( &
  'LBUT1ENNO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*9_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  8.8E1_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_269 = t_meta_chemprop( &
  'LC4H9NO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*9_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.4E-3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_270 = t_meta_chemprop( &
  'LMEKNO3', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  0.86E0_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_271 = t_meta_chemprop( &
  'C1ODC2O2C4OD', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*7_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_272 = t_meta_chemprop( &
  'C1ODC2O2C4OOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_273 = t_meta_chemprop( &
  'C1ODC2OOHC4OD', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_274 = t_meta_chemprop( &
  'C1ODC3O2C4OOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_275 = t_meta_chemprop( &
  'C1OOHC2O2C4OD', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_276 = t_meta_chemprop( &
  'C1OOHC2OOHC4OD', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_277 = t_meta_chemprop( &
  'C1OOHC3O2C4OD', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_278 = t_meta_chemprop( &
  'C4MDIAL', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_279 = t_meta_chemprop( &
  'C511O2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*7_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_280 = t_meta_chemprop( &
  'C511OOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_281 = t_meta_chemprop( &
  'C512O2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*7_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_282 = t_meta_chemprop( &
  'C512OOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_283 = t_meta_chemprop( &
  'C513CO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_284 = t_meta_chemprop( &
  'C513O2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*7_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_285 = t_meta_chemprop( &
  'C513OOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_286 = t_meta_chemprop( &
  'C514O2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*7_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_287 = t_meta_chemprop( &
  'C514OOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_288 = t_meta_chemprop( &
  'C59O2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_289 = t_meta_chemprop( &
  'C59OOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.E11_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_290 = t_meta_chemprop( &
  'C5H8', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.01_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_291 = t_meta_chemprop( &
  'CHOC3COCO3', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*5_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_292 = t_meta_chemprop( &
  'CHOC3COOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_293 = t_meta_chemprop( &
  'CO13C4CHO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_294 = t_meta_chemprop( &
  'CO23C4CHO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_295 = t_meta_chemprop( &
  'CO23C4CO3', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*5_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_296 = t_meta_chemprop( &
  'CO23C4CO3H', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_297 = t_meta_chemprop( &
  'DB1O', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_298 = t_meta_chemprop( &
  'DB1O2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_299 = t_meta_chemprop( &
  'DB1OOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_300 = t_meta_chemprop( &
  'DB2O2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_301 = t_meta_chemprop( &
  'DB2OOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_302 = t_meta_chemprop( &
  'HCOC5', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.7E2_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_303 = t_meta_chemprop( &
  'ISOPAOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_304 = t_meta_chemprop( &
  'ISOPBO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_305 = t_meta_chemprop( &
  'ISOPBOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_306 = t_meta_chemprop( &
  'ISOPBOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_307 = t_meta_chemprop( &
  'ISOPDO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_308 = t_meta_chemprop( &
  'ISOPDOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_309 = t_meta_chemprop( &
  'ISOPDOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_310 = t_meta_chemprop( &
  'MBO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E6_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_311 = t_meta_chemprop( &
  'MBOACO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_312 = t_meta_chemprop( &
  'MBOCOCO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_313 = t_meta_chemprop( &
  'ME3FURAN', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E6_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_314 = t_meta_chemprop( &
  'ACCOMECHO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.7E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_315 = t_meta_chemprop( &
  'ACCOMECO3', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*5_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_316 = t_meta_chemprop( &
  'ACCOMECO3H', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.1E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_317 = t_meta_chemprop( &
  'C24O3CCO2H', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  7.4E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_318 = t_meta_chemprop( &
  'C4CO2DBCO3', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*3_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_319 = t_meta_chemprop( &
  'C4CO2DCO3H', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*4_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_320 = t_meta_chemprop( &
  'C5134CO2OH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.1E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_321 = t_meta_chemprop( &
  'C54CO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.6E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_322 = t_meta_chemprop( &
  'C5CO14O2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*5_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_323 = t_meta_chemprop( &
  'C5CO14OH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.2E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_324 = t_meta_chemprop( &
  'C5CO14OOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.1E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_325 = t_meta_chemprop( &
  'C5DIALCO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*4_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_326 = t_meta_chemprop( &
  'C5DIALO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*5_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_327 = t_meta_chemprop( &
  'C5DIALOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.6E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_328 = t_meta_chemprop( &
  'C5DICARB', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.7E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_329 = t_meta_chemprop( &
  'C5DICARBO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*7_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_330 = t_meta_chemprop( &
  'C5DICAROOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_331 = t_meta_chemprop( &
  'MC3ODBCO2H', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.6E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_332 = t_meta_chemprop( &
  'MMALANHY', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*4_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  20._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_333 = t_meta_chemprop( &
  'MMALANHYO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*5_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_334 = t_meta_chemprop( &
  'MMALNHYOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_335 = t_meta_chemprop( &
  'TLFUO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*7_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_336 = t_meta_chemprop( &
  'TLFUONE', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.6E1_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_337 = t_meta_chemprop( &
  'TLFUOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_338 = t_meta_chemprop( &
  'C4MCONO3OH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_339 = t_meta_chemprop( &
  'C514NO3', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*7_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_340 = t_meta_chemprop( &
  'C5PAN9', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*5_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_341 = t_meta_chemprop( &
  'CHOC3COPAN', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*5_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_342 = t_meta_chemprop( &
  'DB1NO3', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_343 = t_meta_chemprop( &
  'ISOPBDNO3O2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_344 = t_meta_chemprop( &
  'ISOPBNO3', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  8.9E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_345 = t_meta_chemprop( &
  'ISOPDNO3', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  8.9E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_346 = t_meta_chemprop( &
  'NC4CHO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*7_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_347 = t_meta_chemprop( &
  'NC4OHCO3', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_348 = t_meta_chemprop( &
  'NC4OHCO3H', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_349 = t_meta_chemprop( &
  'NC4OHCPAN', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*8_DP+MN*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_350 = t_meta_chemprop( &
  'NISOPO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_351 = t_meta_chemprop( &
  'NISOPOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_352 = t_meta_chemprop( &
  'NMBOBCO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_353 = t_meta_chemprop( &
  'ACCOMEPAN', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*5_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.8_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_354 = t_meta_chemprop( &
  'C4CO2DBPAN', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*3_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.0E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_355 = t_meta_chemprop( &
  'C5COO2NO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*5_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.8_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_356 = t_meta_chemprop( &
  'NC4MDCO2H', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*5_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.6E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_357 = t_meta_chemprop( &
  'NTLFUO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_358 = t_meta_chemprop( &
  'NTLFUOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*7_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_359 = t_meta_chemprop( &
  'LC578O2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_360 = t_meta_chemprop( &
  'LC578OOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.E11_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_361 = t_meta_chemprop( &
  'LDISOPACO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_362 = t_meta_chemprop( &
  'LDISOPACO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_363 = t_meta_chemprop( &
  'LHC4ACCHO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_364 = t_meta_chemprop( &
  'LHC4ACCO2H', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.6E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_365 = t_meta_chemprop( &
  'LHC4ACCO3', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*7_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_366 = t_meta_chemprop( &
  'LHC4ACCO3H', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.2E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_367 = t_meta_chemprop( &
  'LIEPOX', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_368 = t_meta_chemprop( &
  'LISOPAB', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*1_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_369 = t_meta_chemprop( &
  'LISOPACO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_370 = t_meta_chemprop( &
  'LISOPACO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_371 = t_meta_chemprop( &
  'LISOPACOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_372 = t_meta_chemprop( &
  'LISOPCD', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*1_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_373 = t_meta_chemprop( &
  'LISOPEFO', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_374 = t_meta_chemprop( &
  'LISOPEFO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_375 = t_meta_chemprop( &
  'LMBOABO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*11_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_376 = t_meta_chemprop( &
  'LMBOABOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*12_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_377 = t_meta_chemprop( &
  'LME3FURANO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*7_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_378 = t_meta_chemprop( &
  'LZCO3C23DBCOD', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*5_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_379 = t_meta_chemprop( &
  'LZCO3HC23DBCOD', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_380 = t_meta_chemprop( &
  'LZCODC23DBCOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_381 = t_meta_chemprop( &
  'LC5PAN1719', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*7_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_382 = t_meta_chemprop( &
  'LISOPACNO3', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_383 = t_meta_chemprop( &
  'LISOPACNO3O2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_384 = t_meta_chemprop( &
  'LMBOABNO3', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*11_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_385 = t_meta_chemprop( &
  'LNISO3', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_386 = t_meta_chemprop( &
  'LNISOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 6.E3_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_387 = t_meta_chemprop( &
  'LNMBOABO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_388 = t_meta_chemprop( &
  'LNMBOABOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_389 = t_meta_chemprop( &
  'LZCPANC23DBCOD', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*5_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.E1_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_390 = t_meta_chemprop( &
  'C614CO', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_391 = t_meta_chemprop( &
  'C614O2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_392 = t_meta_chemprop( &
  'C614OOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*10_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_393 = t_meta_chemprop( &
  'CO235C5CHO', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_394 = t_meta_chemprop( &
  'CO235C6O2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_395 = t_meta_chemprop( &
  'CO235C6OOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*8_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_396 = t_meta_chemprop( &
  'BENZENE', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.8E-1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_397 = t_meta_chemprop( &
  'BZBIPERO2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*7_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_398 = t_meta_chemprop( &
  'BZBIPEROOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*8_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_399 = t_meta_chemprop( &
  'BZEMUCCO', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_400 = t_meta_chemprop( &
  'BZEMUCCO2H', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.6E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_401 = t_meta_chemprop( &
  'BZEMUCCO3', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_402 = t_meta_chemprop( &
  'BZEMUCCO3H', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_403 = t_meta_chemprop( &
  'BZEMUCO2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*7_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_404 = t_meta_chemprop( &
  'BZEMUCOOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*8_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_405 = t_meta_chemprop( &
  'BZEPOXMUC', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.6E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_406 = t_meta_chemprop( &
  'BZOBIPEROH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.0E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_407 = t_meta_chemprop( &
  'C5CO2DBCO3', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_408 = t_meta_chemprop( &
  'C5CO2DCO3H', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_409 = t_meta_chemprop( &
  'C5CO2OHCO3', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_410 = t_meta_chemprop( &
  'C5COOHCO3H', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_411 = t_meta_chemprop( &
  'C6125CO', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.7E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_412 = t_meta_chemprop( &
  'C615CO2O2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*7_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_413 = t_meta_chemprop( &
  'C615CO2OOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.1E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_414 = t_meta_chemprop( &
  'C6CO4DB', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_415 = t_meta_chemprop( &
  'C6H5O', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.9E3_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_416 = t_meta_chemprop( &
  'C6H5O2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_417 = t_meta_chemprop( &
  'C6H5OOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.9E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_418 = t_meta_chemprop( &
  'CATEC1O', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.6E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_419 = t_meta_chemprop( &
  'CATEC1O2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_420 = t_meta_chemprop( &
  'CATEC1OOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.6E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_421 = t_meta_chemprop( &
  'CATECHOL', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*4_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.6E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_422 = t_meta_chemprop( &
  'CPDKETENE', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*4_DP+MO*1_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_423 = t_meta_chemprop( &
  'PBZQCO', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*4_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.6E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_424 = t_meta_chemprop( &
  'PBZQO2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_425 = t_meta_chemprop( &
  'PBZQONE', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*4_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  20._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_426 = t_meta_chemprop( &
  'PBZQOOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_427 = t_meta_chemprop( &
  'PHENO2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*7_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_428 = t_meta_chemprop( &
  'PHENOL', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.9E3_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_429 = t_meta_chemprop( &
  'PHENOOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*8_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_430 = t_meta_chemprop( &
  'C614NO3', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*9_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_431 = t_meta_chemprop( &
  'BZBIPERNO3', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*7_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.9E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_432 = t_meta_chemprop( &
  'BZEMUCNO3', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*7_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.0E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_433 = t_meta_chemprop( &
  'BZEMUCPAN', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.8_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_434 = t_meta_chemprop( &
  'C5CO2DBPAN', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.7E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_435 = t_meta_chemprop( &
  'C5CO2OHPAN', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*8_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.0E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_436 = t_meta_chemprop( &
  'DNPHEN', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*4_DP+MO*5_DP+MN*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 2.3E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_437 = t_meta_chemprop( &
  'DNPHENO2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*10_DP+MN*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_438 = t_meta_chemprop( &
  'DNPHENOOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*10_DP+MN*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 2.3E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_439 = t_meta_chemprop( &
  'HOC6H4NO2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  70._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_440 = t_meta_chemprop( &
  'NBZQO2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*4_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_441 = t_meta_chemprop( &
  'NBZQOOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.4E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_442 = t_meta_chemprop( &
  'NCATECHOL', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.6E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_443 = t_meta_chemprop( &
  'NCATECO2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*9_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_444 = t_meta_chemprop( &
  'NCATECOOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*7_DP+MO*9_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_445 = t_meta_chemprop( &
  'NCPDKETENE', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*3_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_446 = t_meta_chemprop( &
  'NDNPHENO2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*4_DP+MO*12_DP+MN*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_447 = t_meta_chemprop( &
  'NDNPHENOOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*12_DP+MN*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 2.3E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_448 = t_meta_chemprop( &
  'NNCATECO2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*11_DP+MN*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_449 = t_meta_chemprop( &
  'NNCATECOOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*11_DP+MN*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 2.3E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_450 = t_meta_chemprop( &
  'NPHEN1O', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*4_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  70._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_451 = t_meta_chemprop( &
  'NPHEN1O2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*4_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_452 = t_meta_chemprop( &
  'NPHEN1OOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*5_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  70._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_453 = t_meta_chemprop( &
  'NPHENO2', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*6_DP+MO*8_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_454 = t_meta_chemprop( &
  'NPHENOOH', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*7_DP+MO*8_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.6E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_455 = t_meta_chemprop( &
  'C235C6CO3H', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_456 = t_meta_chemprop( &
  'C716O2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_457 = t_meta_chemprop( &
  'C716OOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*10_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_458 = t_meta_chemprop( &
  'C721O2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_459 = t_meta_chemprop( &
  'C721OOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*12_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_460 = t_meta_chemprop( &
  'C722O2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_461 = t_meta_chemprop( &
  'C722OOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*12_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_462 = t_meta_chemprop( &
  'CO235C6CHO', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_463 = t_meta_chemprop( &
  'CO235C6CO3', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_464 = t_meta_chemprop( &
  'MCPDKETENE', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*6_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_465 = t_meta_chemprop( &
  'ROO6R3O', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_466 = t_meta_chemprop( &
  'ROO6R3O2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_467 = t_meta_chemprop( &
  'ROO6R5O2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_468 = t_meta_chemprop( &
  'BENZAL', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*6_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 36._DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_469 = t_meta_chemprop( &
  'C6CO2OHCO3', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_470 = t_meta_chemprop( &
  'C6COOHCO3H', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_471 = t_meta_chemprop( &
  'C6H5CH2O2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_472 = t_meta_chemprop( &
  'C6H5CH2OOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.9E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_473 = t_meta_chemprop( &
  'C6H5CO3', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*5_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_474 = t_meta_chemprop( &
  'C6H5CO3H', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*6_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.4E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_475 = t_meta_chemprop( &
  'C7CO4DB', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.7E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_476 = t_meta_chemprop( &
  'CRESO2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*9_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_477 = t_meta_chemprop( &
  'CRESOL', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.9E3_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_478 = t_meta_chemprop( &
  'CRESOOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*10_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_479 = t_meta_chemprop( &
  'MCATEC1O', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_480 = t_meta_chemprop( &
  'MCATEC1O2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_481 = t_meta_chemprop( &
  'MCATEC1OOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.6E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_482 = t_meta_chemprop( &
  'MCATECHOL', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.6E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_483 = t_meta_chemprop( &
  'OXYL1O2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_484 = t_meta_chemprop( &
  'OXYL1OOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.9E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_485 = t_meta_chemprop( &
  'PHCOOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*6_DP+MO*2_DP, 1.4E4_DP, 6500._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  6.25E-5_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.5E7_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '65-85-0                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_486 = t_meta_chemprop( &
  'PTLQCO', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*6_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.4E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_487 = t_meta_chemprop( &
  'PTLQO2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_488 = t_meta_chemprop( &
  'PTLQONE', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*6_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  41._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_489 = t_meta_chemprop( &
  'PTLQOOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_490 = t_meta_chemprop( &
  'TLBIPERO2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*9_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_491 = t_meta_chemprop( &
  'TLBIPEROOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*10_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_492 = t_meta_chemprop( &
  'TLEMUCCO', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  3.1E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_493 = t_meta_chemprop( &
  'TLEMUCCO2H', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  7.4E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_494 = t_meta_chemprop( &
  'TLEMUCCO3', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_495 = t_meta_chemprop( &
  'TLEMUCCO3H', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.2E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_496 = t_meta_chemprop( &
  'TLEMUCO2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*9_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_497 = t_meta_chemprop( &
  'TLEMUCOOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*10_DP+MO*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_498 = t_meta_chemprop( &
  'TLEPOXMUC', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  41._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_499 = t_meta_chemprop( &
  'TLOBIPEROH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.0E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_500 = t_meta_chemprop( &
  'TOL1O', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.9E3_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_501 = t_meta_chemprop( &
  'TOLUENE', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.5E-1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_502 = t_meta_chemprop( &
  'C7PAN3', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*8_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_503 = t_meta_chemprop( &
  'C6CO2OHPAN', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*8_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.0E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_504 = t_meta_chemprop( &
  'C6H5CH2NO3', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  36._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_505 = t_meta_chemprop( &
  'DNCRES', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*6_DP+MO*5_DP+MN*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 2.3E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_506 = t_meta_chemprop( &
  'DNCRESO2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*10_DP+MN*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_507 = t_meta_chemprop( &
  'DNCRESOOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*10_DP+MN*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 2.3E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_508 = t_meta_chemprop( &
  'MNCATECH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.6E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_509 = t_meta_chemprop( &
  'MNCATECO2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*9_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_510 = t_meta_chemprop( &
  'MNCATECOOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*9_DP+MO*9_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_511 = t_meta_chemprop( &
  'MNCPDKETENE', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*5_DP+MO*2_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_512 = t_meta_chemprop( &
  'MNNCATCOOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*11_DP+MN*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 2.3E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_513 = t_meta_chemprop( &
  'MNNCATECO2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*11_DP+MN*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_514 = t_meta_chemprop( &
  'NCRES1O', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*6_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  70._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_515 = t_meta_chemprop( &
  'NCRES1O2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*6_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_516 = t_meta_chemprop( &
  'NCRES1OOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  70._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_517 = t_meta_chemprop( &
  'NCRESO2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*8_DP+MO*8_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_518 = t_meta_chemprop( &
  'NCRESOOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*9_DP+MO*8_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.6E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_519 = t_meta_chemprop( &
  'NDNCRESO2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*6_DP+MO*10_DP+MN*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_520 = t_meta_chemprop( &
  'NDNCRESOOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*12_DP+MN*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 2.3E3_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_521 = t_meta_chemprop( &
  'NPTLQO2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*6_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_522 = t_meta_chemprop( &
  'NPTLQOOH', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.4E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_523 = t_meta_chemprop( &
  'PBZN', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*5_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  70._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_524 = t_meta_chemprop( &
  'TLBIPERNO3', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*9_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  70._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_525 = t_meta_chemprop( &
  'TLEMUCNO3', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*9_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4.0E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_526 = t_meta_chemprop( &
  'TLEMUCPAN', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  2.8_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_527 = t_meta_chemprop( &
  'TOL1OHNO2', &
  (/ 0 /), & ! integer
  (/ MC*7_DP+MH*7_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  7.0E1_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_528 = t_meta_chemprop( &
  'C721CHO', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*12_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_529 = t_meta_chemprop( &
  'C721CO3', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_530 = t_meta_chemprop( &
  'C721CO3H', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*12_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_531 = t_meta_chemprop( &
  'C810O2', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_532 = t_meta_chemprop( &
  'C810OOH', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*14_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_533 = t_meta_chemprop( &
  'C811O2', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_534 = t_meta_chemprop( &
  'C812O2', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_535 = t_meta_chemprop( &
  'C812OOH', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*14_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_536 = t_meta_chemprop( &
  'C813O2', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_537 = t_meta_chemprop( &
  'C813OOH', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*14_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_538 = t_meta_chemprop( &
  'C85O2', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_539 = t_meta_chemprop( &
  'C85OOH', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*14_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_540 = t_meta_chemprop( &
  'C86O2', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_541 = t_meta_chemprop( &
  'C86OOH', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*14_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_542 = t_meta_chemprop( &
  'C89O2', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_543 = t_meta_chemprop( &
  'C89OOH', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*14_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_544 = t_meta_chemprop( &
  'C8BC', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*14_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E6_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_545 = t_meta_chemprop( &
  'C8BCCO', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*12_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E6_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_546 = t_meta_chemprop( &
  'C8BCO2', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_547 = t_meta_chemprop( &
  'C8BCOOH', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*12_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_548 = t_meta_chemprop( &
  'NORPINIC', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*12_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 6.166E-05_DP, 2.291E-06_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4E13_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_549 = t_meta_chemprop( &
  'EBENZ', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*10_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.2E-1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_550 = t_meta_chemprop( &
  'STYRENE', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*8_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 3.7E-1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_551 = t_meta_chemprop( &
  'STYRENO2', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*9_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_552 = t_meta_chemprop( &
  'STYRENOOH', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*10_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  70._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_553 = t_meta_chemprop( &
  'C721PAN', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*11_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_554 = t_meta_chemprop( &
  'C810NO3', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*14_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_555 = t_meta_chemprop( &
  'C89NO3', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*13_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_556 = t_meta_chemprop( &
  'C8BCNO3', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*11_DP+MO*3_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_557 = t_meta_chemprop( &
  'NSTYRENO2', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*8_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_558 = t_meta_chemprop( &
  'NSTYRENOOH', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*9_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  70._DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_559 = t_meta_chemprop( &
  'LXYL', &
  (/ 0 /), & ! integer
  (/ MC*8_DP+MH*10_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.7E-1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_560 = t_meta_chemprop( &
  'C811CO3', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_561 = t_meta_chemprop( &
  'C811CO3H', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*14_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_562 = t_meta_chemprop( &
  'C85CO3', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_563 = t_meta_chemprop( &
  'C85CO3H', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*12_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_564 = t_meta_chemprop( &
  'C89CO2H', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*14_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  6.6E7_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_565 = t_meta_chemprop( &
  'C89CO3', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_566 = t_meta_chemprop( &
  'C89CO3H', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*14_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_567 = t_meta_chemprop( &
  'C96O2', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_568 = t_meta_chemprop( &
  'C96OOH', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*16_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_569 = t_meta_chemprop( &
  'C97O2', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_570 = t_meta_chemprop( &
  'C97OOH', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*16_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_571 = t_meta_chemprop( &
  'C98O2', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_572 = t_meta_chemprop( &
  'C98OOH', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*16_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_573 = t_meta_chemprop( &
  'NOPINDCO', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*12_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_574 = t_meta_chemprop( &
  'NOPINDO2', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_575 = t_meta_chemprop( &
  'NOPINDOOH', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*14_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_576 = t_meta_chemprop( &
  'NOPINONE', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*14_DP+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E6_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_577 = t_meta_chemprop( &
  'NOPINOO', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*14_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_578 = t_meta_chemprop( &
  'NORPINAL', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*14_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_579 = t_meta_chemprop( &
  'NORPINENOL', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*14_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_580 = t_meta_chemprop( &
  'PINIC', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*14_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 6.166E-05_DP, 2.291E-06_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  4E13_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_581 = t_meta_chemprop( &
  'C811PAN', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*13_DP+MO*7_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_582 = t_meta_chemprop( &
  'C89PAN', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*13_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_583 = t_meta_chemprop( &
  'C96NO3', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*15_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_584 = t_meta_chemprop( &
  'C9PAN2', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*13_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_585 = t_meta_chemprop( &
  'LTMB', &
  (/ 0 /), & ! integer
  (/ MC*9_DP+MH*12_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.2E-1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_586 = t_meta_chemprop( &
  'APINAOO', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_587 = t_meta_chemprop( &
  'APINBOO', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_588 = t_meta_chemprop( &
  'APINENE', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*16_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_589 = t_meta_chemprop( &
  'BPINAO2', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_590 = t_meta_chemprop( &
  'BPINAOOH', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*18_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_591 = t_meta_chemprop( &
  'BPINENE', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*16_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_592 = t_meta_chemprop( &
  'C106O2', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_593 = t_meta_chemprop( &
  'C106OOH', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*16_DP+MO*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_594 = t_meta_chemprop( &
  'C109CO', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_595 = t_meta_chemprop( &
  'C109O2', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_596 = t_meta_chemprop( &
  'C109OOH', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*16_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_597 = t_meta_chemprop( &
  'C96CO3', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_598 = t_meta_chemprop( &
  'CAMPHENE', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*16_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_599 = t_meta_chemprop( &
  'CARENE', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*16_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_600 = t_meta_chemprop( &
  'MENTHEN6ONE', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*16_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_601 = t_meta_chemprop( &
  'OH2MENTHEN6ONE', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*17_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_602 = t_meta_chemprop( &
  'OHMENTHEN6ONEO2', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_603 = t_meta_chemprop( &
  'PERPINONIC', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*16_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  7.4E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_604 = t_meta_chemprop( &
  'PINAL', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*16_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_605 = t_meta_chemprop( &
  'PINALO2', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_606 = t_meta_chemprop( &
  'PINALOOH', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*14_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_607 = t_meta_chemprop( &
  'PINENOL', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*16_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_608 = t_meta_chemprop( &
  'PINONIC', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*16_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  7.4E5_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_609 = t_meta_chemprop( &
  'RO6R1O2', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_610 = t_meta_chemprop( &
  'RO6R3O2', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_611 = t_meta_chemprop( &
  'ROO6R1O2', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_612 = t_meta_chemprop( &
  'SABINENE', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*16_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_613 = t_meta_chemprop( &
  'BPINANO3', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*17_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_614 = t_meta_chemprop( &
  'C106NO3', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*15_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_615 = t_meta_chemprop( &
  'C10PAN2', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*15_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_616 = t_meta_chemprop( &
  'PINALNO3', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*13_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_617 = t_meta_chemprop( &
  'RO6R1NO3', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*17_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_618 = t_meta_chemprop( &
  'ROO6R1NO3', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*17_DP+MO*6_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_619 = t_meta_chemprop( &
  'LAPINABNO3', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*17_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_620 = t_meta_chemprop( &
  'LAPINABO2', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_621 = t_meta_chemprop( &
  'LAPINABOOH', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*18_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_622 = t_meta_chemprop( &
  'LNAPINABO2', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_623 = t_meta_chemprop( &
  'LNAPINABOOH', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*17_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_624 = t_meta_chemprop( &
  'LNBPINABO2', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*10_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_625 = t_meta_chemprop( &
  'LNBPINABOOH', &
  (/ 0 /), & ! integer
  (/ MC*10_DP+MH*17_DP+MO*5_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, 1.0E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_626 = t_meta_chemprop( &
  'LHAROM', &
  (/ 0 /), & ! integer
  (/ MC*11_DP+MH*14_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.2E-1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_627 = t_meta_chemprop( &
  'CHF3', &
  (/ 0 /), & ! integer
  (/ MC+MH+MF*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_628 = t_meta_chemprop( &
  'CHF2CF3', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH+MF*5_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_629 = t_meta_chemprop( &
  'CH3CF3', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MF*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_630 = t_meta_chemprop( &
  'CH2F2', &
  (/ 0 /), & ! integer
  (/ MC+MH*2_DP+MF*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_631 = t_meta_chemprop( &
  'CH3CHF2', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*4_DP+MF*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_632 = t_meta_chemprop( &
  'CCl4', &
  (/ 0 /), & ! integer
  (/ MC+MCl*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_633 = t_meta_chemprop( &
  'CF2Cl2', &
  (/ 0 /), & ! integer
  (/ MC+MF*2_DP+MCl*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_634 = t_meta_chemprop( &
  'CF2ClCF2Cl', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MF*4_DP+MCl*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_635 = t_meta_chemprop( &
  'CF2ClCFCl2', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MF*3_DP+MCl*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_636 = t_meta_chemprop( &
  'CF3CF2Cl', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MF*5_DP+MCl, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_637 = t_meta_chemprop( &
  'CFCl3', &
  (/ 0 /), & ! integer
  (/ MC+MF+MCl*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_638 = t_meta_chemprop( &
  'CH2Cl2', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*2_DP+MCl*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_639 = t_meta_chemprop( &
  'CH2FCF3', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*2_DP+MF*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_640 = t_meta_chemprop( &
  'CH3CCl3', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MCl*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_641 = t_meta_chemprop( &
  'CH3CFCl2', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MF+MCl*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, &
  0.1_DP, 0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  9.99999_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_642 = t_meta_chemprop( &
  'CH3Cl', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MCl, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.4E-2_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_643 = t_meta_chemprop( &
  'CHCl3', &
  (/ 0 /), & ! integer
  (/ MC+MH+MCl*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_644 = t_meta_chemprop( &
  'CHF2Cl', &
  (/ 0 /), & ! integer
  (/ MC+MH+MF*2_DP+MCl, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_645 = t_meta_chemprop( &
  'Cl', &
  (/ 0 /), & ! integer
  (/ MCl, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_646 = t_meta_chemprop( &
  'Cl2', &
  (/ 0 /), & ! integer
  (/ MCl*2_DP, 9.3E-2_DP, 2000._DP, 9.99999_DP, 0.038_DP, 6546._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 7.E-2_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '7782-50-5               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_647 = t_meta_chemprop( &
  'Cl2O2', &
  (/ 0 /), & ! integer
  (/ MCl*2_DP+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_648 = t_meta_chemprop( &
  'ClNO2', &
  (/ 0 /), & ! integer
  (/ MCl+MN+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_649 = t_meta_chemprop( &
  'ClNO3', &
  (/ 0 /), & ! integer
  (/ MCl+MN+MO*3_DP, BIG_DP, 0._DP, 9.99999_DP, 0.108_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.E30_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '14545-72-3              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_650 = t_meta_chemprop( &
  'ClO', &
  (/ 0 /), & ! integer
  (/ MCl+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_651 = t_meta_chemprop( &
  'HCl', &
  (/ 0 /), & ! integer
  (/ MH+MCl, 2._DP/1.7_DP, 9000._DP, 9.99999_DP, 0.074_DP, 3072._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1E14_DP, &
  1.0_DP, 9.99999_DP /), & ! real
  (/ '7647-01-0               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_652 = t_meta_chemprop( &
  'HOCl', &
  (/ 0 /), & ! integer
  (/ MH+MO+MCl, 6.6E2_DP, 5880._DP, 9.99999_DP, 0.5_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 6.7E2_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '7790-92-3               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_653 = t_meta_chemprop( &
  'OClO', &
  (/ 0 /), & ! integer
  (/ MCl+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_654 = t_meta_chemprop( &
  'LCHLORINE', &
  (/ 0 /), & ! integer
  (/ MCl, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_655 = t_meta_chemprop( &
  'Br', &
  (/ 0 /), & ! integer
  (/ MBr, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_656 = t_meta_chemprop( &
  'Br2', &
  (/ 0 /), & ! integer
  (/ MBr*2_DP, 7.25E-1_DP, 4390._DP, 9.99999_DP, 0.038_DP, 6546._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.7_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '7726-95-6               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_657 = t_meta_chemprop( &
  'BrCl', &
  (/ 0 /), & ! integer
  (/ MBr+MCl, 9.4E-1_DP, 5600._DP, 9.99999_DP, 0.038_DP, 6546._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0_DP, 0.1_DP, &
  9.99999_DP /), & ! real
  (/ '13863-41-7              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_658 = t_meta_chemprop( &
  'BrNO2', &
  (/ 0 /), & ! integer
  (/ MBr+MN+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_659 = t_meta_chemprop( &
  'BrNO3', &
  (/ 0 /), & ! integer
  (/ MBr+MN+MO*3_DP, BIG_DP, 0._DP, 9.99999_DP, 0.063_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E30_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '40423-14-1              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_660 = t_meta_chemprop( &
  'BrO', &
  (/ 0 /), & ! integer
  (/ MBr+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_661 = t_meta_chemprop( &
  'CF2ClBr', &
  (/ 0 /), & ! integer
  (/ MC+MF*2_DP+MCl+MBr, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_662 = t_meta_chemprop( &
  'CF3Br', &
  (/ 0 /), & ! integer
  (/ MC+MF*3_DP+MBr, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_663 = t_meta_chemprop( &
  'CH2Br2', &
  (/ 0 /), & ! integer
  (/ MC+MH*2_DP+MBr*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.1_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_664 = t_meta_chemprop( &
  'CH2ClBr', &
  (/ 0 /), & ! integer
  (/ MC+MH*2_DP+MCl+MBr, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.1_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_665 = t_meta_chemprop( &
  'CH3Br', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MBr, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.6E-1_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_666 = t_meta_chemprop( &
  'CHBr3', &
  (/ 0 /), & ! integer
  (/ MC+MH+MBr*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.7_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_667 = t_meta_chemprop( &
  'CHCl2Br', &
  (/ 0 /), & ! integer
  (/ MC+MH+MCl*2_DP+MBr, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 4.0E-1_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_668 = t_meta_chemprop( &
  'CHClBr2', &
  (/ 0 /), & ! integer
  (/ MC+MH+MCl+MBr*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 8.7E-1_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_669 = t_meta_chemprop( &
  'HBr', &
  (/ 0 /), & ! integer
  (/ MH+MBr, 1.3_DP, 10000._DP, 9.99999_DP, 0.032_DP, 3940._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.3E17_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '10035-10-6              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_670 = t_meta_chemprop( &
  'HOBr', &
  (/ 0 /), & ! integer
  (/ MH+MO+MBr, 1.3E3_DP, 5862._DP, 9.99999_DP, 0.5_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.1E1_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '13517-11-8              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_671 = t_meta_chemprop( &
  'LBROMINE', &
  (/ 0 /), & ! integer
  (/ MBr, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_672 = t_meta_chemprop( &
  'C3H7I', &
  (/ 0 /), & ! integer
  (/ MI+MH*7_DP+MC*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_673 = t_meta_chemprop( &
  'CH2ClI', &
  (/ 0 /), & ! integer
  (/ MC+MH*2_DP+MCl+MI, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.4E1_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_674 = t_meta_chemprop( &
  'CH2I2', &
  (/ 0 /), & ! integer
  (/ MI*2_DP+MC+MH*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_675 = t_meta_chemprop( &
  'CH3I', &
  (/ 0 /), & ! integer
  (/ MI+MC+MH*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_676 = t_meta_chemprop( &
  'HI', &
  (/ 0 /), & ! integer
  (/ MH+MI, BIG_DP, 0._DP, 9.99999_DP, 0.036_DP, 4130._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.E30_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '10034-85-2              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_677 = t_meta_chemprop( &
  'HIO3', &
  (/ 0 /), & ! integer
  (/ MH+MI+MO*3_DP, BIG_DP, 0._DP, 9.99999_DP, 0.01_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.E30_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_678 = t_meta_chemprop( &
  'HOI', &
  (/ 0 /), & ! integer
  (/ MH+MO+MI, 4.5E2_DP, 5862._DP, 9.99999_DP, 0.5_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.E30_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '14332-21-9              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_679 = t_meta_chemprop( &
  'I', &
  (/ 0 /), & ! integer
  (/ MI, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_680 = t_meta_chemprop( &
  'I2', &
  (/ 0 /), & ! integer
  (/ MI*2_DP, 3._DP, 4431._DP, 9.99999_DP, 0.01_DP, 2000._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 3.0_DP, 0.1_DP, &
  9.99999_DP /), & ! real
  (/ '7553-56-2               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_681 = t_meta_chemprop( &
  'I2O2', &
  (/ 0 /), & ! integer
  (/ MI*2_DP+MO*2_DP, BIG_DP, 0._DP, 9.99999_DP, 0.1_DP, 2000._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.0_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_682 = t_meta_chemprop( &
  'IBr', &
  (/ 0 /), & ! integer
  (/ MI+MBr, 2.4E1_DP, 5600._DP, 9.99999_DP, 0.018_DP, 2000._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.0_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '7789-33-5               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_683 = t_meta_chemprop( &
  'ICl', &
  (/ 0 /), & ! integer
  (/ MI+MCl, 1.1E2_DP, 5600._DP, 9.99999_DP, 0.018_DP, 2000._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.1E2_DP, 0.1_DP, &
  9.99999_DP /), & ! real
  (/ '7790-99-0               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_684 = t_meta_chemprop( &
  'INO2', &
  (/ 0 /), & ! integer
  (/ MI+MN+MO*2_DP, BIG_DP, 0._DP, 9.99999_DP, 0.1_DP, 2000._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 4.5_DP, 0.1_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_685 = t_meta_chemprop( &
  'INO3', &
  (/ 0 /), & ! integer
  (/ MI+MN+MO*3_DP, BIG_DP, 0._DP, 9.99999_DP, 0.1_DP, 2000._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.E30_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_686 = t_meta_chemprop( &
  'IO', &
  (/ 0 /), & ! integer
  (/ MI+MO, 4.5E2_DP, 5862._DP, 9.99999_DP, 0.5_DP, 2000._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.0_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_687 = t_meta_chemprop( &
  'IPART', &
  (/ 0 /), & ! integer
  (/ MI*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_688 = t_meta_chemprop( &
  'OIO', &
  (/ 0 /), & ! integer
  (/ MO+MI+MO, BIG_DP, 0._DP, 9.99999_DP, 0.01_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.0_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_689 = t_meta_chemprop( &
  'CH3SO2', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MS+MO*2_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_690 = t_meta_chemprop( &
  'CH3SO3', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MS+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '75-75-2                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_691 = t_meta_chemprop( &
  'CH3SO3H', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MS+MO*3_DP+MH, BIG_DP, 0._DP, 9.99999_DP, 0.076_DP, 1762._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.E30_DP, &
  1.0_DP, 9.99999_DP /), & ! real
  (/ '75-75-2                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_692 = t_meta_chemprop( &
  'DMS', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MS+MC+MH*3_DP, 5.4E-1_DP, 3460._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.0_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '75-18-3                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_693 = t_meta_chemprop( &
  'DMSO', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MS+MO+MC+MH*3_DP, 9.5E4_DP, 1300._DP, 9.99999_DP, 0.048_DP, &
  2578._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  5.E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '67-68-5                 ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_694 = t_meta_chemprop( &
  'CH3SOOH', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MS+MO*2_DP+MH, 1.0E8_DP, 1760._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.0_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '17696-73-0              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_695 = t_meta_chemprop( &
  'HPMTF', &
  (/ 0 /), & ! integer
  (/ MH+MO+MO+MC+MH*2_DP+MS+MC+MH+MO, 1.33E4_DP, 1300._DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  0.0_DP, 0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_696 = t_meta_chemprop( &
  'HOOCH2SCO', &
  (/ 0 /), & ! integer
  (/ MH+MO+MO+MC+MH*2_DP+MS+MC+MO, 1.33E4_DP, 1300._DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  0.0_DP, 0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_697 = t_meta_chemprop( &
  'H2SO4', &
  (/ 0 /), & ! integer
  (/ MH*2_DP+MS+MO*4_DP, 1.E11_DP, 0._DP, 9.99999_DP, 0.65_DP, 0._DP, 1E3_DP, &
  1.023E-2_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.3E15_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '7664-93-9               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_698 = t_meta_chemprop( &
  'OCS', &
  (/ 0 /), & ! integer
  (/ MC+MS+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_699 = t_meta_chemprop( &
  'S', &
  (/ 0 /), & ! integer
  (/ MS, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_700 = t_meta_chemprop( &
  'SF6', &
  (/ 0 /), & ! integer
  (/ MS+MF*6_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_701 = t_meta_chemprop( &
  'SH', &
  (/ 0 /), & ! integer
  (/ MS+MH, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_702 = t_meta_chemprop( &
  'SO', &
  (/ 0 /), & ! integer
  (/ MS+MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_703 = t_meta_chemprop( &
  'SO2', &
  (/ 0 /), & ! integer
  (/ MS+MO*2_DP, 1.3_DP, 2900._DP, 9.99999_DP, 0.11_DP, 0._DP, 1.23E-2_DP, &
  6.61E-8_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.45E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '7446-09-5               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_704 = t_meta_chemprop( &
  'SO3', &
  (/ 0 /), & ! integer
  (/ MS+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_705 = t_meta_chemprop( &
  'LSULFUR', &
  (/ 0 /), & ! integer
  (/ MS, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_706 = t_meta_chemprop( &
  'Hg', &
  (/ 0 /), & ! integer
  (/ MHg, 0.13_DP, 2700._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 0.13_DP, 0.1_DP, &
  9.99999_DP /), & ! real
  (/ '7439-97-6               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_707 = t_meta_chemprop( &
  'HgO', &
  (/ 0 /), & ! integer
  (/ MHg+MO, 3.2E6_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.4E7_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '21908-53-2              ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_708 = t_meta_chemprop( &
  'HgCl', &
  (/ 0 /), & ! integer
  (/ MHg+MCl, 2.4E7_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.4E7_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '7546-30-7               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_709 = t_meta_chemprop( &
  'HgCl2', &
  (/ 0 /), & ! integer
  (/ MHg+MCl*2_DP, 2.4E7_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.4E7_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '7487-94-7               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_710 = t_meta_chemprop( &
  'HgBr', &
  (/ 0 /), & ! integer
  (/ MHg+MBr, 2.4E7_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.4E7_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_711 = t_meta_chemprop( &
  'HgBr2', &
  (/ 0 /), & ! integer
  (/ MHg+MBr*2_DP, 2.4E7_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.4E7_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '7789-47-1               ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_712 = t_meta_chemprop( &
  'ClHgBr', &
  (/ 0 /), & ! integer
  (/ MCl+MHg+MBr, 2.4E7_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.4E7_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_713 = t_meta_chemprop( &
  'BrHgOBr', &
  (/ 0 /), & ! integer
  (/ MBr+MHg+MO+MBr, 2.4E7_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.4E7_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_714 = t_meta_chemprop( &
  'ClHgOBr', &
  (/ 0 /), & ! integer
  (/ MCl+MHg+MO+MBr, 2.4E7_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.4E7_DP, 1.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_715 = t_meta_chemprop( &
  'LTERP', &
  (/ 0 /), & ! integer
  (/ 136.24_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 4.9E-2_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_716 = t_meta_chemprop( &
  'LALK4', &
  (/ 0 /), & ! integer
  (/ 73.23_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 7.7E-4_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_717 = t_meta_chemprop( &
  'LALK5', &
  (/ 0 /), & ! integer
  (/ 106.97_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 2.5E-4_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_718 = t_meta_chemprop( &
  'LARO1', &
  (/ 0 /), & ! integer
  (/ 100.47_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.4E-1_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_719 = t_meta_chemprop( &
  'LARO2', &
  (/ 0 /), & ! integer
  (/ 113.93_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.4E-1_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_720 = t_meta_chemprop( &
  'LOLE1', &
  (/ 0 /), & ! integer
  (/ 61.68_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 5.0E-3_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_721 = t_meta_chemprop( &
  'LOLE2', &
  (/ 0 /), & ! integer
  (/ 79.05_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 5.0E-3_DP, &
  0.0_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_722 = t_meta_chemprop( &
  'LfPOG02', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_723 = t_meta_chemprop( &
  'LfPOG03', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_724 = t_meta_chemprop( &
  'LfPOG04', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_725 = t_meta_chemprop( &
  'LfPOG05', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_726 = t_meta_chemprop( &
  'LbbPOG02', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_727 = t_meta_chemprop( &
  'LbbPOG03', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_728 = t_meta_chemprop( &
  'LbbPOG04', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_729 = t_meta_chemprop( &
  'LfSOGsv01', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_730 = t_meta_chemprop( &
  'LfSOGsv02', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_731 = t_meta_chemprop( &
  'LbbSOGsv01', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_732 = t_meta_chemprop( &
  'LbbSOGsv02', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_733 = t_meta_chemprop( &
  'LfSOGiv01', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_734 = t_meta_chemprop( &
  'LfSOGiv02', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_735 = t_meta_chemprop( &
  'LfSOGiv03', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_736 = t_meta_chemprop( &
  'LfSOGiv04', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_737 = t_meta_chemprop( &
  'LbbSOGiv01', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_738 = t_meta_chemprop( &
  'LbbSOGiv02', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_739 = t_meta_chemprop( &
  'LbbSOGiv03', &
  (/ 0 /), & ! integer
  (/ 250._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_740 = t_meta_chemprop( &
  'LbSOGv01', &
  (/ 0 /), & ! integer
  (/ 180._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_741 = t_meta_chemprop( &
  'LbSOGv02', &
  (/ 0 /), & ! integer
  (/ 180._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_742 = t_meta_chemprop( &
  'LbSOGv03', &
  (/ 0 /), & ! integer
  (/ 180._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_743 = t_meta_chemprop( &
  'LbSOGv04', &
  (/ 0 /), & ! integer
  (/ 180._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_744 = t_meta_chemprop( &
  'LbOSOGv01', &
  (/ 0 /), & ! integer
  (/ 180._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_745 = t_meta_chemprop( &
  'LbOSOGv02', &
  (/ 0 /), & ! integer
  (/ 180._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_746 = t_meta_chemprop( &
  'LbOSOGv03', &
  (/ 0 /), & ! integer
  (/ 180._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_747 = t_meta_chemprop( &
  'LaSOGv01', &
  (/ 0 /), & ! integer
  (/ 150._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_748 = t_meta_chemprop( &
  'LaSOGv02', &
  (/ 0 /), & ! integer
  (/ 150._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_749 = t_meta_chemprop( &
  'LaSOGv03', &
  (/ 0 /), & ! integer
  (/ 150._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_750 = t_meta_chemprop( &
  'LaSOGv04', &
  (/ 0 /), & ! integer
  (/ 150._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_751 = t_meta_chemprop( &
  'LaOSOGv01', &
  (/ 0 /), & ! integer
  (/ 150._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_752 = t_meta_chemprop( &
  'LaOSOGv02', &
  (/ 0 /), & ! integer
  (/ 150._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_753 = t_meta_chemprop( &
  'LaOSOGv03', &
  (/ 0 /), & ! integer
  (/ 150._DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 1.0E5_DP, 0.0_DP, &
  9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_754 = t_meta_chemprop( &
  'ISO2', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_755 = t_meta_chemprop( &
  'ISON', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*9_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.7E4_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_756 = t_meta_chemprop( &
  'ISOOH', &
  (/ 0 /), & ! integer
  (/ MC*5_DP+MH*10_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  1.7E6_DP, 0.1_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_757 = t_meta_chemprop( &
  'LHOC3H6O2', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*7_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_758 = t_meta_chemprop( &
  'LHOC3H6OOH', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*8_DP+MO*3_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_759 = t_meta_chemprop( &
  'MVKO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*7_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_760 = t_meta_chemprop( &
  'MVKOOH', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*8_DP+MO*4_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_761 = t_meta_chemprop( &
  'NACA', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*3_DP+MO*4_DP+MN, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_762 = t_meta_chemprop( &
  'N2O3', &
  (/ 0 /), & ! integer
  (/ MN*2_DP+MO*3_DP, 6.0E-1_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_763 = t_meta_chemprop( &
  'N2O4', &
  (/ 0 /), & ! integer
  (/ MN*2_DP+MO*4_DP, 1.4_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_764 = t_meta_chemprop( &
  'MGLYOXAC', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MC+MO+MC+MO+MO+MH, 3.11E5_DP, 0._DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_765 = t_meta_chemprop( &
  'BSOV', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*6_DP+MO*4_DP, 1.E8_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_766 = t_meta_chemprop( &
  'BLOV', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*10_DP+MO*4_DP, 1.E8_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_767 = t_meta_chemprop( &
  'MMA', &
  (/ 0 /), & ! integer
  (/ MC+MH*3_DP+MN+MH*2_DP, 36.5_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_768 = t_meta_chemprop( &
  'DMA', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*7_DP+MN, 31.4_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_769 = t_meta_chemprop( &
  'TMA', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*9_DP+MN, 9.6_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_770 = t_meta_chemprop( &
  'MEA', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*7_DP+MO+MN, 6.18E6_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_771 = t_meta_chemprop( &
  'DEA', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*11_DP+MO*2_DP+MN, 2.58E7_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_772 = t_meta_chemprop( &
  'TEA', &
  (/ 0 /), & ! integer
  (/ MC*6_DP+MH*15_DP+MO*3_DP+MN, 1.42E9_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_773 = t_meta_chemprop( &
  'HNCO', &
  (/ 0 /), & ! integer
  (/ MH+MN+MC+MO, 2.50E1_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_774 = t_meta_chemprop( &
  'H2NCHO', &
  (/ 0 /), & ! integer
  (/ MH*2_DP+MN+MC+MH+MO, 7.19E5_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_775 = t_meta_chemprop( &
  'H2NCOCH2OH', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*5_DP+MO*2_DP+MN, 5.83E5_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_776 = t_meta_chemprop( &
  'CH3NHCHO', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*5_DP+MO+MN, 2.98E4_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_777 = t_meta_chemprop( &
  'DMNCHO', &
  (/ 0 /), & ! integer
  (/ MC*3_DP+MH*7_DP+MO+MN, 1.36E4_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_778 = t_meta_chemprop( &
  'MMNNO2', &
  (/ 0 /), & ! integer
  (/ MC+MH*4_DP+MO*2_DP+MN*2_DP, 6.85E2_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_779 = t_meta_chemprop( &
  'MEANNO2', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*6_DP+MO*3_DP+MN*2_DP, 4.45E6_DP, 0._DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_780 = t_meta_chemprop( &
  'MEANNO', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*6_DP+MO*2_DP+MN*2_DP, 1.00E7_DP, 0._DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_781 = t_meta_chemprop( &
  'NDMA', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*6_DP+MO+MN*2_DP, 2.75E2_DP, 0._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_782 = t_meta_chemprop( &
  'DMNNO2', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*6_DP+MO*2_DP+MN*2_DP, 3.13E2_DP, 0._DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_783 = t_meta_chemprop( &
  'NDELA', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*10_DP+MO*3_DP+MN*2_DP, 2.06E8_DP, 0._DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_784 = t_meta_chemprop( &
  'DEANNO2', &
  (/ 0 /), & ! integer
  (/ MC*4_DP+MH*10_DP+MO*4_DP+MN*2_DP, 2.37E9_DP, 0._DP, 9.99999_DP, 0.1_DP, &
  0._DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, &
  -9.99999_DP, -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_785 = t_meta_chemprop( &
  'ONE', &
  (/ 0 /), & ! integer
  (/ 1.0_DP, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_786 = t_meta_chemprop( &
  'O', &
  (/ 0 /), & ! integer
  (/ MO, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_787 = t_meta_chemprop( &
  'C', &
  (/ 0 /), & ! integer
  (/ MC, -9.99999_DP, -9.99999_DP, 9.99999_DP, 0.1_DP, 0._DP, 9.99999_DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, -9.99999_DP, &
  -9.99999_DP, 9.99999_DP /), & ! real
  (/ '?                       ' /) ) ! string

TYPE(t_meta_chemprop), PARAMETER :: spc_788 = t_meta_chemprop( &
  'OXL', &
  (/ 0 /), & ! integer
  (/ MC*2_DP+MH*2_DP+MO*4_DP, 3.2E6_DP, 7285._DP, 9.99999_DP, 0.1_DP, 0._DP, &
  9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 9.99999_DP, 3.26E6_DP, &
  0.1_DP, 9.99999_DP /), & ! real
  (/ '144-62-7                ' /) ) ! string

! combine all species into one array:
INTEGER, PUBLIC, PARAMETER :: N_CHEMPROP = 788
TYPE(t_meta_chemprop), PUBLIC, PARAMETER, DIMENSION(N_CHEMPROP) :: chemprop = &
  (/ spc_001, spc_002, spc_003, spc_004, spc_005, spc_006, spc_007, spc_008, &
  spc_009, spc_010, spc_011, spc_012, spc_013, spc_014, spc_015, spc_016, &
  spc_017, spc_018, spc_019, spc_020, spc_021, spc_022, spc_023, spc_024, &
  spc_025, spc_026, spc_027, spc_028, spc_029, spc_030, spc_031, spc_032, &
  spc_033, spc_034, spc_035, spc_036, spc_037, spc_038, spc_039, spc_040, &
  spc_041, spc_042, spc_043, spc_044, spc_045, spc_046, spc_047, spc_048, &
  spc_049, spc_050, spc_051, spc_052, spc_053, spc_054, spc_055, spc_056, &
  spc_057, spc_058, spc_059, spc_060, spc_061, spc_062, spc_063, spc_064, &
  spc_065, spc_066, spc_067, spc_068, spc_069, spc_070, spc_071, spc_072, &
  spc_073, spc_074, spc_075, spc_076, spc_077, spc_078, spc_079, spc_080, &
  spc_081, spc_082, spc_083, spc_084, spc_085, spc_086, spc_087, spc_088, &
  spc_089, spc_090, spc_091, spc_092, spc_093, spc_094, spc_095, spc_096, &
  spc_097, spc_098, spc_099, spc_100, spc_101, spc_102, spc_103, spc_104, &
  spc_105, spc_106, spc_107, spc_108, spc_109, spc_110, spc_111, spc_112, &
  spc_113, spc_114, spc_115, spc_116, spc_117, spc_118, spc_119, spc_120, &
  spc_121, spc_122, spc_123, spc_124, spc_125, spc_126, spc_127, spc_128, &
  spc_129, spc_130, spc_131, spc_132, spc_133, spc_134, spc_135, spc_136, &
  spc_137, spc_138, spc_139, spc_140, spc_141, spc_142, spc_143, spc_144, &
  spc_145, spc_146, spc_147, spc_148, spc_149, spc_150, spc_151, spc_152, &
  spc_153, spc_154, spc_155, spc_156, spc_157, spc_158, spc_159, spc_160, &
  spc_161, spc_162, spc_163, spc_164, spc_165, spc_166, spc_167, spc_168, &
  spc_169, spc_170, spc_171, spc_172, spc_173, spc_174, spc_175, spc_176, &
  spc_177, spc_178, spc_179, spc_180, spc_181, spc_182, spc_183, spc_184, &
  spc_185, spc_186, spc_187, spc_188, spc_189, spc_190, spc_191, spc_192, &
  spc_193, spc_194, spc_195, spc_196, spc_197, spc_198, spc_199, spc_200, &
  spc_201, spc_202, spc_203, spc_204, spc_205, spc_206, spc_207, spc_208, &
  spc_209, spc_210, spc_211, spc_212, spc_213, spc_214, spc_215, spc_216, &
  spc_217, spc_218, spc_219, spc_220, spc_221, spc_222, spc_223, spc_224, &
  spc_225, spc_226, spc_227, spc_228, spc_229, spc_230, spc_231, spc_232, &
  spc_233, spc_234, spc_235, spc_236, spc_237, spc_238, spc_239, spc_240, &
  spc_241, spc_242, spc_243, spc_244, spc_245, spc_246, spc_247, spc_248, &
  spc_249, spc_250, spc_251, spc_252, spc_253, spc_254, spc_255, spc_256, &
  spc_257, spc_258, spc_259, spc_260, spc_261, spc_262, spc_263, spc_264, &
  spc_265, spc_266, spc_267, spc_268, spc_269, spc_270, spc_271, spc_272, &
  spc_273, spc_274, spc_275, spc_276, spc_277, spc_278, spc_279, spc_280, &
  spc_281, spc_282, spc_283, spc_284, spc_285, spc_286, spc_287, spc_288, &
  spc_289, spc_290, spc_291, spc_292, spc_293, spc_294, spc_295, spc_296, &
  spc_297, spc_298, spc_299, spc_300, spc_301, spc_302, spc_303, spc_304, &
  spc_305, spc_306, spc_307, spc_308, spc_309, spc_310, spc_311, spc_312, &
  spc_313, spc_314, spc_315, spc_316, spc_317, spc_318, spc_319, spc_320, &
  spc_321, spc_322, spc_323, spc_324, spc_325, spc_326, spc_327, spc_328, &
  spc_329, spc_330, spc_331, spc_332, spc_333, spc_334, spc_335, spc_336, &
  spc_337, spc_338, spc_339, spc_340, spc_341, spc_342, spc_343, spc_344, &
  spc_345, spc_346, spc_347, spc_348, spc_349, spc_350, spc_351, spc_352, &
  spc_353, spc_354, spc_355, spc_356, spc_357, spc_358, spc_359, spc_360, &
  spc_361, spc_362, spc_363, spc_364, spc_365, spc_366, spc_367, spc_368, &
  spc_369, spc_370, spc_371, spc_372, spc_373, spc_374, spc_375, spc_376, &
  spc_377, spc_378, spc_379, spc_380, spc_381, spc_382, spc_383, spc_384, &
  spc_385, spc_386, spc_387, spc_388, spc_389, spc_390, spc_391, spc_392, &
  spc_393, spc_394, spc_395, spc_396, spc_397, spc_398, spc_399, spc_400, &
  spc_401, spc_402, spc_403, spc_404, spc_405, spc_406, spc_407, spc_408, &
  spc_409, spc_410, spc_411, spc_412, spc_413, spc_414, spc_415, spc_416, &
  spc_417, spc_418, spc_419, spc_420, spc_421, spc_422, spc_423, spc_424, &
  spc_425, spc_426, spc_427, spc_428, spc_429, spc_430, spc_431, spc_432, &
  spc_433, spc_434, spc_435, spc_436, spc_437, spc_438, spc_439, spc_440, &
  spc_441, spc_442, spc_443, spc_444, spc_445, spc_446, spc_447, spc_448, &
  spc_449, spc_450, spc_451, spc_452, spc_453, spc_454, spc_455, spc_456, &
  spc_457, spc_458, spc_459, spc_460, spc_461, spc_462, spc_463, spc_464, &
  spc_465, spc_466, spc_467, spc_468, spc_469, spc_470, spc_471, spc_472, &
  spc_473, spc_474, spc_475, spc_476, spc_477, spc_478, spc_479, spc_480, &
  spc_481, spc_482, spc_483, spc_484, spc_485, spc_486, spc_487, spc_488, &
  spc_489, spc_490, spc_491, spc_492, spc_493, spc_494, spc_495, spc_496, &
  spc_497, spc_498, spc_499, spc_500, spc_501, spc_502, spc_503, spc_504, &
  spc_505, spc_506, spc_507, spc_508, spc_509, spc_510, spc_511, spc_512, &
  spc_513, spc_514, spc_515, spc_516, spc_517, spc_518, spc_519, spc_520, &
  spc_521, spc_522, spc_523, spc_524, spc_525, spc_526, spc_527, spc_528, &
  spc_529, spc_530, spc_531, spc_532, spc_533, spc_534, spc_535, spc_536, &
  spc_537, spc_538, spc_539, spc_540, spc_541, spc_542, spc_543, spc_544, &
  spc_545, spc_546, spc_547, spc_548, spc_549, spc_550, spc_551, spc_552, &
  spc_553, spc_554, spc_555, spc_556, spc_557, spc_558, spc_559, spc_560, &
  spc_561, spc_562, spc_563, spc_564, spc_565, spc_566, spc_567, spc_568, &
  spc_569, spc_570, spc_571, spc_572, spc_573, spc_574, spc_575, spc_576, &
  spc_577, spc_578, spc_579, spc_580, spc_581, spc_582, spc_583, spc_584, &
  spc_585, spc_586, spc_587, spc_588, spc_589, spc_590, spc_591, spc_592, &
  spc_593, spc_594, spc_595, spc_596, spc_597, spc_598, spc_599, spc_600, &
  spc_601, spc_602, spc_603, spc_604, spc_605, spc_606, spc_607, spc_608, &
  spc_609, spc_610, spc_611, spc_612, spc_613, spc_614, spc_615, spc_616, &
  spc_617, spc_618, spc_619, spc_620, spc_621, spc_622, spc_623, spc_624, &
  spc_625, spc_626, spc_627, spc_628, spc_629, spc_630, spc_631, spc_632, &
  spc_633, spc_634, spc_635, spc_636, spc_637, spc_638, spc_639, spc_640, &
  spc_641, spc_642, spc_643, spc_644, spc_645, spc_646, spc_647, spc_648, &
  spc_649, spc_650, spc_651, spc_652, spc_653, spc_654, spc_655, spc_656, &
  spc_657, spc_658, spc_659, spc_660, spc_661, spc_662, spc_663, spc_664, &
  spc_665, spc_666, spc_667, spc_668, spc_669, spc_670, spc_671, spc_672, &
  spc_673, spc_674, spc_675, spc_676, spc_677, spc_678, spc_679, spc_680, &
  spc_681, spc_682, spc_683, spc_684, spc_685, spc_686, spc_687, spc_688, &
  spc_689, spc_690, spc_691, spc_692, spc_693, spc_694, spc_695, spc_696, &
  spc_697, spc_698, spc_699, spc_700, spc_701, spc_702, spc_703, spc_704, &
  spc_705, spc_706, spc_707, spc_708, spc_709, spc_710, spc_711, spc_712, &
  spc_713, spc_714, spc_715, spc_716, spc_717, spc_718, spc_719, spc_720, &
  spc_721, spc_722, spc_723, spc_724, spc_725, spc_726, spc_727, spc_728, &
  spc_729, spc_730, spc_731, spc_732, spc_733, spc_734, spc_735, spc_736, &
  spc_737, spc_738, spc_739, spc_740, spc_741, spc_742, spc_743, spc_744, &
  spc_745, spc_746, spc_747, spc_748, spc_749, spc_750, spc_751, spc_752, &
  spc_753, spc_754, spc_755, spc_756, spc_757, spc_758, spc_759, spc_760, &
  spc_761, spc_762, spc_763, spc_764, spc_765, spc_766, spc_767, spc_768, &
  spc_769, spc_770, spc_771, spc_772, spc_773, spc_774, spc_775, spc_776, &
  spc_777, spc_778, spc_779, spc_780, spc_781, spc_782, spc_783, spc_784, &
  spc_785, spc_786, spc_787, spc_788 /)
