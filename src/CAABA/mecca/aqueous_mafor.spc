{-----------------------------------------------------------------------------}
{------------------------------ aerosol mode: ## -----------------------------}
{-----------------------------------------------------------------------------}

{------------------------------- neutral species -----------------------------}

{------------------------------------- O -------------------------------------}

{------------------------------------- H -------------------------------------}

{------------------------------------- N -------------------------------------}
N2O3_a##       =              3O  + 2N  ; {@N_2O_3}          {dinitrogen trioxide}
N2O4_a##       =              4O  + 2N  ; {@N_2O_4}          {dinitrogen tetraoxide}
 
{------------------------------------- C -------------------------------------}

{1C}
CH2O2H2_a##    =   C +  4H +  2O        ; {@CH_2(OH)_2}      {}
MMA_a##        =   C +  5H        +  N  ; {@MMA}             {methylamine}
NH2CH2_a##     =   C +  4H        +  N  ; {@CH_2NH_2}        {methylamine radical}
HNCO_a##       =   C +   H +   O  +  N  ; {@HNCO}            {ioscyanic acid} 
H2NCHO_a##     =   C +  3H +   O  +  N  ; {@H2NCHO}          {formamide}
MMNNO2_a##     =   C +  2H +  2O  + 2N  ; {@MMNNO2}          {methylnitramine}
MSIA_a##       =   C +  4H +  S + 2O    ; {@MSIA}            {methyl sulfinic acid}

{2C}
OXALAC_a##     =  IGNORE                ; {@OXALAC}          {oxalic acid, 2C +  2H  +  4O}
HCOCO2H_a##    =  2C +  2H  +  3O       ; {@HCOCO_2H}        {oxoethanoic acid}
HOCH2CHO_a##   =  2C +  4H  +  2O       ; {@HOCH_2CHO}       {glycolaldehyde}
HOCH2CO2H_a##  =  2C +  4H  +  3O       ; {@HOCH_2CO_2H}     {hydroxyethanoic acid}
CH3CO3_a##     =  2C +  3H  +  3O       ; {@CH_3COO_2}       {peroxyacetyl radical}
GLYOX_a##      =  2C +  2H  +  2O       ; {@GLYOX}           {CHOCHO = glyoxal} 
DMA_a##        =  2C +  7H         +  N ; {@DMA}             {dimethylamine}
MEA_a##        =  2C +  7H  +   O  +  N ; {@MEA}             {ethanolamine}
MEANNO_a##     =  2C +  6H  +  2O  + 2N ; {@MEANNO}          {N-nitroso ethanolamine} 
MEANNO2_a##    =  2C +  6H  +  3O  + 2N ; {@MEANNO2}         {N-nitro ethanolamine}
NDMA_a##       =  2C +  6H  +   O  + 2N ; {@NDMA}            {N-nitroso dimethylamine}
DMNNO2_a##     =  2C +  6H  +  2O  + 2N ; {@DMNNO2}          {dimethylnitramine}
CH3NHCH2_a##   =  2C +  6H  +         N ; {@CH_3NHCH_2}      {methylamine methyl radical} 
CH3NHNHCH3_a## =  2C +  8H  +        2N ; {@CH_3NHNHCH_3}    {dimethylhydrazine} 
NH2C2H4NH2_a## =  2C +  8H  +        2N ; {@NH_2CH_2CH_2NH_2}  {ethylenediamine} 
NH2CH2CHOH_a## =  2C +  6H  +   O  +  N ; {@NH_2CH_2CHOH}    {ethanolamine radical} 
H2NCOCH2OH_a## =  2C +  5H  +  2O  +  N ; {@H2NCOCH2OH}      {2-hydroxy acetamide} 
CH3NHCHO_a##   =  2C +  5H  +   O  +  N ; {@CH_3NHCHO}       {N-methyl formamide} 
CH3NCO_a##     =  2C +  3H  +   O  +  N ; {@CH_3NCO}         {methyl isocyanic acid}
HPMTF_a##      =  2C +  4H  +  3O  +  S ; {@HPMTF}           {hydroperoxyl methyl thioformate}
HOOCH2SCO_a##  =  2C +  3H  +  3O  +  S ; {@HOOCH2SCO}       {}
 
{3C}
MGLYOX_a##     =  3C +  4H  +  2O       ; {@MGLYOX}          {methylglyoxal}
MGLYOAC_a##    =  3C +  4H  +  3O       ; {@MGLYOAC}         {methylglyoxylic acid}
DOC_a##        =  IGNORE                ; {@DOC}             {dissolved organic carbon DOC}
DOCO_a##       =  IGNORE                ; {@DOCO}            {oxidized DOC}
TMA_a##        =  3C +  9H         +  N ; {@TMA}             {trimethylamine}
DMNCH2_a##     =  3C +  8H         +  N ; {@(CH_3)_2NCH_2}   {dimethylamine methyl radical} 
DMNCHO_a##     =  3C +  7H  +   O  +  N ; {@DMNCHO}          {N,N-dimethyl formamide} 
MALONAC_a##    =  IGNORE                ; {@MALONAC}         {malonic acid, 3C +  4H  +  4O}
 
{4C}
DEA_a##        =  4C + 11H  +  2O  +  N ; {@DEA}             {diethanolamine} 
NDELA_a##      =  4C + 10H  +  3O  + 2N ; {@NDELA}           {N-nitroso diethanolamine}
DEANNO2_a##    =  4C + 10H  +  4O  + 2N ; {@DEANNO2}         {N-nitro diethanolamine}
DEAN_a##       =  4C + 10H  +  2O  +  N ; {@HOETNHCH_2CHOH}  {diethanolamine radical} 
SUCCAC_a##     =  IGNORE                ; {@SUCCAC}          {succinic acid, 4C +  6H  +  4O}

{5C}
GLUTARAC_a##   =  IGNORE                ; {@GLUTARAC}        {glutaric acid, 5C +  8H  +  4O}

{6C}
TEA_a##        =  6C + 15H  +  3O  +  N ; {@TEA}             {triethanolamine} 
DENCH2CHOH_a## =  6C + 14H  +  3O  +  N ; {@DENCH_2CHOH}     {triethanolamine radical}  
ADIPAC_a##     =  IGNORE                ; {@ADIPAC}          {adipic acid, 6C + 10H  +  4O}


{----------------------------------- ions ------------------------------------}

{------------------------------------- O -------------------------------------}

{------------------------------------- H -------------------------------------}

{------------------------------------- N -------------------------------------}

{------------------------------------- C -------------------------------------}

{1C}
MMAp_a##       =   C +  6H         +  N  + Pls  ; {@MMA^+}           {methylaminium}
MMNp_a##       =   C +  5H         +  N  + Pls  ; {@CH_3NH_2^+}      {methylamine N-radical cation} 
NH2CH2p_a##    =   C +  4H         +  N  + Pls  ; {@CH_2NH_2^+}      {iminium}
NH3CH2p_a##    =   C +  5H         +  N  + Pls  ; {@CH_2NH_3^+}      {methylaminium radical} 
NCOm_a##       =   C        +   O  +  N  + Min  ; {@NCO^-}           {isocyanate}


{2C}
HC2O4m_a##     =  IGNORE                 + Min  ; {@HC_2O_4^-}       {hydrogen oxalate,2C +   H  +  4O}
C2O4mm_a##     =  IGNORE                 + 2Min ; {@C_2O_4^<2->}     {oxalate, 2C        +  4O}
HCOCOOm_a##    =  2C +   H  +  3O        + Min  ; {@HCOCOO^-}        {}
MEAp_a##       =  2C +  8H  +   O  +  N  + Pls  ; {@MEA^+}           {ethanolaminium} 
DMAp_a##       =  2C +  8H         +  N  + Pls  ; {@DMA^+}           {dimethylaminium}
DMNp_a##       =  2C +  7H         +  N  + Pls  ; {@(CH_3)_2NH^+}    {dimethylamine N-radical cation} 
CH3NHCH2p_a##  =  2C +  6H         +  N  + Pls  ; {@CH_3NH^+CH_2}    {methyl iminium}
CH3NH2CH2p_a## =  2C +  7H         +  N  + Pls  ; {@CH_3NH_2^+CH_2}    {dimethylaminium radical} 
MENp_a##       =  2C +  7H  +   O  +  N  + Pls  ; {@HOCH_2CH_2NH_2^+}  {ethanolamine N-radical cation} 
NH3CH2CHOHp_a## = 2C +  7H  +   O  +  N  + Pls  ; {@HOCHCH_2NH_3^+}  {ethanolaminium radical} 
 
{3C}
CH3COCOOm_a##  =  3C +  3H   + 3O        + Min  ; {@CH_3COCOO^-}     {methylglyoxalate}
TMAp_a##       =  3C + 10H         +  N  + Pls  ; {@TMA^+}           {trimethylaminium}
TMNp_a##       =  3C +  9H         +  N  + Pls  ; {@(CH_3)_3N^+}     {trimethylamine N-radical cation} 
DMNCH2p_a##    =  3C +  8H         +  N  + Pls  ; {@(CH_3)_2N^+CH_2}  {dimethyl iminium}
DMNHCH2p_a##   =  3C +  9H         +  N  + Pls  ; {@(CH_3)_2NH^+CH_2} {trimethylaminium radical} 
 
{4C}
DEAp_a##       =  4C + 12H  +  2O  +  N  + Pls  ; {@DEA^+}                {diethanolaminium}
DENp_a##       =  4C + 13H  +  2O  +  N  + Pls  ; {@(HOET)_2NH^+}         {diethanolamine N-radical cation} 
DENHp_a##      =  4C + 12H  +  2O  +  N  + Pls  ; {@HOETNH_2CH_2CHOH^+}   {diethanolaminium radical} 
C2H5C2O4m_a##  =  IGNORE                 + Min  ; {@CH_2CH_2HC_2O_4^-}    {hydrogen succinate, 4C +  5H  +  4O}
C2H4C2O4mm_a## =  IGNORE                 + 2Min ; {@CH_2CH_2C_2O_4^<2->}  {succinate, 4C +  4H  +  4O}

{6C} 
TEAp_a##       =  6C + 16H  +  3O  +  N  + Pls ; {@TEA^+}                 {triethanolaminium}
TENp_a##       =  6C + 15H  +  3O  +  N  + Pls ; {@(HOET)_3N^+}           {triethanolamine N-radical cation} 
DENIMp_a##     =  6C + 15H  +  3O  +  N  + Pls ; {@(HOET)_2N^+CH_2CH_2OH} {diethanol iminium}
TENHp_a##      =  6C + 15H  +  3O  +  N  + Pls ; {@(HOET)_2NH^+CH_2CHOH}  {triethanolaminium radical} 
 

{------------------------------------M----------------------------------------}
