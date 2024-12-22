!****************************************************************************
!                Time-stamp: <2023-12-11 21:10:48 sander>
!****************************************************************************

! Definitions that all photolysis submodels have in common

! Author:
! Rolf Sander, MPICH, 2008-2020

! NOTES:
! - The index ip_BrNO3 refers to the total photolysis of BrNO3. The
!   index ip_BrONO2 is only for the branch BrNO3 + hv -> BrO + NO2

MODULE messy_cmn_photol_mem

  IMPLICIT NONE

  ! ip_* = index of photolysis
  INTEGER, PUBLIC, PARAMETER :: &
    ip_O2           =   1, ip_O3P          =   2, ip_O1D          =   3, &
    ip_H2O2         =   4, ip_NO2          =   5, ip_NO2O         =   6, &
    ip_NOO2         =   7, ip_N2O5         =   8, ip_HNO3         =   9, &
    ip_HNO4         =  10, ip_PAN          =  11, ip_HONO         =  12, &
    ip_CH3OOH       =  13, ip_COH2         =  14, ip_CHOH         =  15, &
    ip_CH3CO3H      =  16, ip_CH3CHO       =  17, ip_CH3COCH3     =  18, &
    ip_MGLYOX       =  19, ip_HOCl         =  20, ip_OClO         =  21, &
    ip_Cl2O2        =  22, ip_ClNO3        =  23, ip_ClNO2        =  24, &
    ip_Cl2          =  25, ip_BrO          =  26, ip_HOBr         =  27, &
    ip_BrCl         =  28, ip_BrNO3        =  29, ip_BrNO2        =  30, &
    ip_Br2          =  31, ip_CCl4         =  32, ip_CH3Cl        =  33, &
    ip_CH3CCl3      =  34, ip_CFCl3        =  35, ip_CF2Cl2       =  36, &
    ip_CH3Br        =  37, ip_CF2ClBr      =  38, ip_CF3Br        =  39, &
    ip_CH3I         =  40, ip_C3H7I        =  41, ip_CH2ClI       =  42, &
    ip_CH2I2        =  43, ip_IO           =  44, ip_HOI          =  45, &
    ip_I2           =  46, ip_ICl          =  47, ip_IBr          =  48, &
    ip_INO2         =  49, ip_INO3         =  50, ip_SO2          =  51, &
    ip_SO3          =  52, ip_OCS          =  53, ip_CS2          =  54, &
    ip_H2O          =  55, ip_N2O          =  56, ip_NO           =  57, &
    ip_CO2          =  58, ip_HCl          =  59, ip_CHCl2Br      =  60, &
    ip_CHClBr2      =  61, ip_CH2ClBr      =  62, ip_CH2Br2       =  63, &
    ip_CHBr3        =  64, ip_SF6          =  65, ip_NO3NOO       =  66, &
    ip_ClONO2       =  67, ip_MACR         =  68, ip_MVK          =  69, &
    ip_GLYOX        =  70, ip_HOCH2CHO     =  71, ip_CH4          =  72, &
    ip_O2_b1b2      =  73, ip_O2_b1        =  74, ip_O2_b2        =  75, &
    ip_O3PO1D       =  76, ip_O3Pp         =  77, ip_H2O1D        =  78, &
    ip_N2           =  79, ip_N2_b1        =  80, ip_N2_b2        =  81, &
    ip_N2_b3        =  82, ip_NN2D         =  83, ip_NOp          =  84, &
    ip_Op_em        =  85, ip_O2p_em       =  86, ip_Op_O_em      =  87, &
    ip_N2p_em       =  88, ip_Np_N_em      =  89, ip_Np_N2D_em    =  90, &
    ip_N_N2D_em     =  91, ip_Op_em_b      =  92, ip_se_O2_b1     =  93, &
    ip_se_O2_b2     =  94, ip_se_N2_b1     =  95, ip_se_N2_b2     =  96, &
    ip_se_N2_b3     =  97, ip_se_N2_b4     =  98, ip_se_Op_em     =  99, &
    ip_O2_aurq      = 100, ip_N2_aurq      = 101, ip_H2SO4        = 102, &
    ip_C3O2         = 103, ip_CH3NO3       = 104, ip_CH3O2NO2     = 105, &
    ip_CH3ONO       = 106, ip_CH3O2        = 107, ip_HCOOH        = 108, &
    ip_HO2NO2       = 109, ip_OHNO3        = 110, ip_BrONO2       = 111, &
    ip_CH3OCl       = 112, ip_MEO2NO2dummy = 113, ip_CHF2Cl       = 114, &
    ip_F113         = 115, ip_C2H5NO3      = 116, ip_NOA          = 117, &
    ip_MEKNO3       = 118, ip_BENZAL       = 119, ip_HOPh3Me2NO2  = 120, &
    ip_HOC6H4NO2    = 121, ip_CH3CHO2VINY  = 122, ip_CH3COCO2H    = 123, &
    ip_IPRCHO2HCO   = 124, ip_C2H5CHO2HCO  = 125, ip_C2H5CHO2ENOL = 126, &
    ip_C3H7CHO2HCO  = 127, ip_C3H7CHO2VINY = 128, ip_PeDIONE24    = 129, &
    ip_PINAL2HCO    = 130, ip_PINAL2ENOL   = 131, ip_CF2ClCFCl2   = 132, &
    ip_CH3CFCl2     = 133, ip_CF3CF2Cl     = 134, ip_CF2ClCF2Cl   = 135, &
    ip_CHCl3        = 136, ip_CH2Cl2       = 137, ip_HO2          = 138, &
    ip_ClO          = 139, ip_HOOCCOOH     = 140, ip_CBrF2CBrF2   = 141, &
    ip_CH3CF2Cl     = 142, ip_MEK          = 143, ip_ACETOL       = 144, &
    ip_IC3H7NO3     = 145, ip_HI           = 146, ip_OIO          = 147, &
    ip_I2O2         = 148, ip_I2O3         = 149, ip_I2O4         = 150, &
    ip_INO          = 151, ip_Cl2O         = 152, ip_Cl2O3        = 153, &
    ip_ClNO         = 154, ip_ClONO        = 155

  ! The next time a new photolysis is added, it can be inserted here.

  ! IP_MAX must be set to the highest ip_* value from the definitions above:
  INTEGER, PUBLIC, PARAMETER :: IP_MAX = 155

  CHARACTER(LEN=12), PUBLIC, PARAMETER, DIMENSION(IP_MAX) :: jname = (/ &
    'O2          ', 'O3P         ', 'O1D         ', &
    'H2O2        ', 'NO2         ', 'NO2O        ', &
    'NOO2        ', 'N2O5        ', 'HNO3        ', &
    'HNO4        ', 'PAN         ', 'HONO        ', &
    'CH3OOH      ', 'COH2        ', 'CHOH        ', &
    'CH3CO3H     ', 'CH3CHO      ', 'CH3COCH3    ', &
    'MGLYOX      ', 'HOCl        ', 'OClO        ', &
    'Cl2O2       ', 'ClNO3       ', 'ClNO2       ', &
    'Cl2         ', 'BrO         ', 'HOBr        ', &
    'BrCl        ', 'BrNO3       ', 'BrNO2       ', &
    'Br2         ', 'CCl4        ', 'CH3Cl       ', &
    'CH3CCl3     ', 'CFCl3       ', 'CF2Cl2      ', &
    'CH3Br       ', 'CF2ClBr     ', 'CF3Br       ', &
    'CH3I        ', 'C3H7I       ', 'CH2ClI      ', &
    'CH2I2       ', 'IO          ', 'HOI         ', &
    'I2          ', 'ICl         ', 'IBr         ', &
    'INO2        ', 'INO3        ', 'SO2         ', &
    'SO3         ', 'OCS         ', 'CS2         ', &
    'H2O         ', 'N2O         ', 'NO          ', &
    'CO2         ', 'HCl         ', 'CHCl2Br     ', &
    'CHClBr2     ', 'CH2ClBr     ', 'CH2Br2      ', &
    'CHBr3       ', 'SF6         ', 'NO3NOO      ', &
    'ClONO2      ', 'MACR        ', 'MVK         ', &
    'GLYOX       ', 'HOCH2CHO    ', 'CH4         ', &
    'O2_b1b2     ', 'O2_b1       ', 'O2_b2       ', &
    'O3PO1D      ', 'O3Pp        ', 'H2O1D       ', &
    'N2          ', 'N2_b1       ', 'N2_b2       ', &
    'N2_b3       ', 'NN2D        ', 'NOp         ', &
    'Op_em       ', 'O2p_em      ', 'Op_O_em     ', &
    'N2p_em      ', 'Np_N_em     ', 'Np_N2D_em   ', &
    'N_N2D_em    ', 'Op_em_b     ', 'se_O2_b1    ', &
    'se_O2_b2    ', 'se_N2_b1    ', 'se_N2_b2    ', &
    'se_N2_b3    ', 'se_N2_b4    ', 'se_Op_em    ', &
    'O2_aurq     ', 'N2_aurq     ', 'H2SO4       ', &
    'C3O2        ', 'CH3NO3      ', 'CH3O2NO2    ', &
    'CH3ONO      ', 'CH3O2       ', 'HCOOH       ', &
    'HO2NO2      ', 'OHNO3       ', 'BrONO2      ', &
    'CH3OCl      ', 'MEO2NO2dummy', 'CHF2Cl      ', &
    'F113        ', 'C2H5NO3     ', 'NOA         ', &
    'MEKNO3      ', 'BENZAL      ', 'HOPh3Me2NO2 ', &
    'HOC6H4NO2   ', 'CH3CHO2VINY ', 'CH3COCO2H   ', &
    'IPRCHO2HCO  ', 'C2H5CHO2HCO ', 'C2H5CHO2ENOL', &
    'C3H7CHO2HCO ', 'C3H7CHO2VINY', 'PeDIONE24   ', &
    'PINAL2HCO   ', 'PINAL2ENOL  ', 'CF2ClCFCl2  ', &
    'CH3CFCl2    ', 'CF3CF2Cl    ', 'CF2ClCF2Cl  ', &
    'CHCl3       ', 'CH2Cl2      ', 'HO2         ', &
    'ClO         ', 'HOOCCOOH    ', 'CBrF2CBrF2  ', &
    'CH3CF2Cl    ', 'MEK         ', 'ACETOL      ', &
    'IC3H7NO3    ', 'HI          ', 'OIO         ', &
    'I2O2        ', 'I2O3        ', 'I2O4        ', &
    'INO         ', 'Cl2O        ', 'Cl2O3       ', &
    'ClNO        ', 'ClONO       ' /)

END MODULE messy_cmn_photol_mem

!*****************************************************************************
