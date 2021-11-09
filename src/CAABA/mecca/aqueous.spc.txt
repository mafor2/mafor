{-----------------------------------------------------------------------------}
{------------------------------ aerosol mode: ## -----------------------------}
{-----------------------------------------------------------------------------}

{------------------------------- neutral species -----------------------------}

{------------------------------------- O -------------------------------------}

O2_a##         = 2O                   ; {@O_2}          {oxygen}
O3_a##         = 3O                   ; {@O_3}          {ozone}

{------------------------------------- H -------------------------------------}

OH_a##         =  H +  O              ; {@OH}           {hydroxyl radical}
HO2_a##        =  H + 2O              ; {@HO_2}         {perhydroxyl radical}
H2O_a##        = 2H +  O              ; {@H_2O}         {water}
H2O2_a##       = 2H + 2O              ; {@H_2O_2}       {hydrogen peroxide}

{------------------------------------- N -------------------------------------}

NH3_a##        = 3H      +  N         ; {@NH_3}         {ammonia}
NO_a##         =       O +  N         ; {@NO}           {nitric oxide}
NO2_a##        =      2O +  N         ; {@NO_2}         {nitrogen dioxide}
NO3_a##        =      3O +  N         ; {@NO_3}         {nitrogen trioxide}
HONO_a##       =  H + 2O +  N         ; {@HONO}         {nitrous acid}
HNO3_a##       =  H + 3O +  N         ; {@HNO_3}        {nitric acid}
HNO4_a##       =  H + 4O +  N         ; {@HNO_4}        {pernitric acid}

{------------------------------------- C -------------------------------------}

{1C}
CH3OH_a##      =   C +  4H +   O      ; {@CH_3OH}       {methanol}
HCOOH_a##      =   C +  2H +  2O      ; {@HCOOH}        {formic acid}
HCHO_a##       =   C +  2H +   O      ; {@HCHO}         {methanal (formaldehyde)}
CH3O2_a##      =   C +  3H +  2O      ; {@CH_3OO}       {methylperoxy radical}
CH3OOH_a##     =   C +  4H +  2O      ; {@CH_3OOH}      {}
CO2_a##        =   C       +  2O      ; {@CO_2}         {carbon dioxide}

{2C}
CH3CO2H_a##    =  2C +  4H +  2O      ; {@CH_3COOH}     {acetic acid}
PAN_a##        =  2C +  3H +  5O +  N ; {@PAN}          {peroxyacetylnitrate}
CH3CHO_a##     =  2C +  4H +   O      ; {@CH_3CHO}      {acetaldehyde}

{3C}
CH3COCH3_a##   =  3C +  6H +   O      ; {@CH_3COCH_3}   {acetone}

{------------------------------------- Cl ------------------------------------}

Cl_a##         = Cl                   ; {@Cl}           {chlorine atom}
Cl2_a##        = 2Cl                  ; {@Cl_2}         {molecular chlorine}
HCl_a##        = H + Cl               ; {@HCl}          {hydrogen chloride}
HOCl_a##       = H + O + Cl           ; {@HOCl}         {hypochlorous acid}

{------------------------------------- Br ------------------------------------}

Br_a##         = Br                   ; {@Br}           {bromine atom}
Br2_a##        = 2Br                  ; {@Br_2}         {molecular bromine}
HBr_a##        = H + Br               ; {@HBr}          {hydrogen bromide}
HOBr_a##       = H + O + Br           ; {@HOBr}         {hypobromous acid}
BrCl_a##       = Br + Cl              ; {@BrCl}         {bromine chloride}

{------------------------------------- I -------------------------------------}

I2_a##         = 2I                   ; {@I_2}          {molecular iodine}
IO_a##         = I + O                ; {@IO}           {iodine monoxide radical}
HOI_a##        = H + O + I            ; {@HOI}          {hypoiodous acid}
ICl_a##        = I + Cl               ; {@ICl}          {iodine chloride}
IBr_a##        = I + Br               ; {@IBr}          {iodine bromide}

{------------------------------------- S -------------------------------------}

SO2_a##        = S + 2O               ; {@SO_2}         {sulfur dioxide}
H2SO4_a##      = 2H + S + 4O          ; {@H_2SO_4}      {sulfuric acid}
DMS_a##        = 2C + 6H + S          ; {@DMS}          {dimethyl sulfide: CH3SCH3}
DMSO_a##       = 2C + 6H + S + O      ; {@DMSO}         {dimethyl sulfoxide: CH3SOCH3}

{------------------------------------- Hg ------------------------------------}

Hg_a##         = Hg                   ; {@Hg}           {mercury}
HgO_a##        = Hg + O               ; {@HgO}          {} 
HgOHOH_a##     = Hg + 2O + 2H         ; {@Hg(OH)_2}     {}
HgOHCl_a##     = Hg + O + H + Cl      ; {@Hg(OH)Cl}     {}
HgCl2_a##      = Hg + 2Cl             ; {@HgCl_2}       {}
HgBr2_a##      = Hg + 2Br             ; {@HgBr_2}       {}
HgSO3_a##      = Hg + S + 3O          ; {@HgSO_3}       {}
ClHgBr_a##     = Hg + Cl + Br         ; {@ClHgBr}       {}
BrHgOBr_a##    = Hg + O + 2Br         ; {@BrHgOBr}      {}
ClHgOBr_a##    = Hg + O + Cl + Br     ; {@ClHgOBr}      {}

{------------------------------------Fe---------------------------------------}

FeOH3_a##      = Fe + 3O + 3H         ; {@FeOH3}        {}
FeCl3_a##      = Fe + 3Cl             ; {@FeCl3}        {}
FeF3_a##       = Fe + 3F              ; {@FeF3}         {}

{----------------------------------- ions ------------------------------------}

{------------------------------------- O -------------------------------------}

O2m_a##        = 2O            + Min  ; {@O_2^-}        {}
OHm_a##        = H +  O        + Min  ; {@OH^-}         {}
HO2m_a##       = H + 2O        + Min  ; {@HO2^-}        {}
O2mm_a##       = 2O            + 2Min ; {@O2^<2->}      {}

{------------------------------------- H -------------------------------------}

Hp_a##         =  H             + Pls ; {@H^+}          {}

{------------------------------------- N -------------------------------------}

NH4p_a##       = N + 4H         + Pls ; {@NH_4^+}       {ammonium}
NO2m_a##       =      2O +  N   + Min ; {@NO_2^-}       {nitrite}
NO3m_a##       =      3O +  N   + Min ; {@NO_3^-}       {nitrate}
NO4m_a##       =      4O +  N   + Min ; {@NO_4^-}       {peroxy nitrate}

{------------------------------------- C -------------------------------------}

{1C}
CO3m_a##       = C + 3O         + Min ; {@CO_3^-}       {}
HCOOm_a##      = H + C + 2O     + Min ; {@HCOO^-}       {formate}
HCO3m_a##      = H + C + 3O     + Min ; {@HCO_3^-}      {hydrogen carbonate}

{2C}
CH3COOm_a##    = 2C + 3H + 2O   + Min ; {@CH_3COO^-}    {acetate}

{------------------------------------- Cl ------------------------------------}

Clm_a##        = Cl             + Min ; {@Cl^-}         {chloride}
Cl2m_a##       = 2Cl            + Min ; {@Cl_2^-}       {}
ClOm_a##       = Cl + O         + Min ; {@ClO^-}        {}
ClOHm_a##      = H + O + Cl     + Min ; {@ClOH^-}       {}

{------------------------------------- Br ------------------------------------}

Brm_a##        = Br             + Min ; {@Br^-}         {bromide}
Br2m_a##       = 2Br            + Min ; {@Br_2^-}       {}
BrOm_a##       = Br + O         + Min ; {@BrO^-}        {}
BrOHm_a##      = H + O + Br     + Min ; {@BrOH^-}       {}
BrCl2m_a##     = Br + 2Cl       + Min ; {@BrCl_2^-}     {}
Br2Clm_a##     = 2Br + Cl       + Min ; {@Br_2Cl^-}     {}

{------------------------------------- I -------------------------------------}

Im_a##         = I              + Min ; {@I^-}          {iodide}
IO2m_a##       = I + 2O         + Min ; {@IO_2^-}       {}
IO3m_a##       = I + 3O         + Min ; {@IO_3^-}       {iodate}
ICl2m_a##      = I + 2Cl        + Min ; {@ICl_2^-}      {}
IBr2m_a##      = I + 2Br        + Min ; {@IBr_2^-}      {}

{------------------------------------- S -------------------------------------}

SO3m_a##       = S + 3O          + Min ; {@SO_3^-}       {}
SO3mm_a##      = S + 3O         + 2Min ; {@SO_3^<2->}    {sulfite}
SO4m_a##       = S + 4O          + Min ; {@SO_4^-}       {}
SO4mm_a##      = S + 4O         + 2Min ; {@SO_4^<2->}    {sulfate}
SO5m_a##       = S + 5O          + Min ; {@SO_5^-}       {}
HSO3m_a##      = H + S + 3O      + Min ; {@HSO_3^-}      {hydrogen sulfite}
HSO4m_a##      = H + S + 4O      + Min ; {@HSO_4^-}      {hydrogen sulfate}
HSO5m_a##      = H + S + 5O      + Min ; {@HSO_5^-}      {}
CH3SO3m_a##    = C + 3H + S + 3O + Min ; {@CH_3SO_3^-}   {MSA anion}
CH2OHSO3m_a##  = C + 3H + S + 4O + Min ; {@CH_2OHSO_3^-} {}

{------------------------------------Hg---------------------------------------}

Hgp_a##        = Hg                +  Pls ; {@Hg^+}              {}
Hgpp_a##       = Hg                + 2Pls ; {@Hg^<2+>}           {}
HgOHp_a##      = Hg + O + H        +  Pls ; {@HgOH^+}            {}
HgClp_a##      = Hg + Cl           +  Pls ; {@HgCl^+}            {}
HgBrp_a##      = Hg + Br           +  Pls ; {@HgBr^+}            {}
HgSO32mm_a##   = Hg + 2S + 6O      + 2Min ; {@Hg(SO_3)_2^<2->}   {}

{------------------------------------Fe---------------------------------------}

Fepp_a##        = Fe             + 2Pls ; {@Fe^<2+>}         {Fe(II)}
FeOpp_a##       = Fe + O         + 2Pls ; {@FeO^<2+>}        {Fe(II)}
FeOHp_a##       = Fe + O + H     + Pls  ; {@FeOH^+}          {Fe(II)}
FeOH2p_a##      = Fe + 2O + 2H   + Pls  ; {@Fe(OH)_2^+}      {Fe(II)}
FeClp_a##       = Fe + Cl        + Pls  ; {@FeCl^+}          {Fe(II)}
Feppp_a##       = Fe             + 3Pls ; {@Fe^<3+>}         {Fe(III)}
FeHOpp_a##      = Fe + O + H     + 2Pls ; {@FeHO^<2+>}       {Fe(III)}
FeHO2pp_a##     = Fe + 2O + H    + 2Pls ; {@FeHO_2^<2+>}     {Fe(III)}
FeOHpp_a##      = Fe + O + H     + 2Pls ; {@FeOH^<2+>}       {Fe(III)}
FeOH4m_a##      = Fe + 4O + 4H   + Min  ; {@Fe(OH)_4^-}      {Fe(III)}
FeOHHO2p_a##    = Fe + 3O + 2H   + Pls  ; {@Fe(OH)(HO_2)^+}  {Fe(III)}
FeClpp_a##      = Fe + Cl        + 2Pls ; {@FeCl^<2+>}       {Fe(III)}
FeCl2p_a##      = Fe + 2Cl       + Pls  ; {@FeCl_2^+}        {Fe(III)}
FeBrpp_a##      = Fe + Br        + 2Pls ; {@FeBr^<2+>}       {Fe(III)}
FeBr2p_a##      = Fe + 2Br       + Pls  ; {@FeBr_2^+}        {Fe(III)}
FeFpp_a##       = Fe + F         + 2Pls ; {@FeF^<2+>}        {Fe(III)}
FeF2p_a##       = Fe + 2F        + 2Pls ; {@FeF_2^+}         {Fe(III)}
FeSO3p_a##      = Fe + 3O + S    + Pls  ; {@FeSO_3^+}        {Fe(III)}
FeSO4p_a##      = Fe + 4O + S    + Pls  ; {@FeSO_4^+}        {Fe(III)}
FeSO42m_a##     = Fe + 8O + 2S   + Min  ; {@Fe(SO_4)_2^-}    {Fe(III)}
FeOH2Fepppp_a## = 2 Fe + O + H   + 4Pls ; {@Fe(OH)_2Fe^<4+>} {Fe(III)}

{-----------------------------------------------------------------------------}
{------------------------------------ dummies --------------------------------}
{-----------------------------------------------------------------------------}

D1O_a##        = Ignore              ; {@D_1O}         {}
Nap_a##        = Ignore              ; {@Na^+}         {dummy cation}
