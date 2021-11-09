{Time-stamp: <2019-01-09 16:19:59 sander>}

{----------------------------------------------------------------------------}

{ SYNTAX AND NAMING CONVENTIONS FOR KPP SPECIES                              }
{ - Species are sorted by elements in the following order:                   }
{   O,H,N,C,F,Cl,Br,I,S,Hg                                                   }
{ - Organics are sorted by increasing number of C, H, O, N                   }
{ - All peroxides are called ROOH, all peroxy radicals are called RO2        }
{ - All species are defined here with #DEFVAR as VARIABLES. Some species     }
{   will be turned into FIXED species with #SETFIX in messy_mecca_kpp.kpp    }
{ - Lumped species start with the letter "L".                                }
{ - The maximum length for the species name is 13 (15 may also be ok?).      }
{ - The species name must not contain the underscore character "_".          }
{ - The elemental composition is needed for graphviz (spc_extract.awk) and   }
{   to check the mass balance (check_conservation.pl). There must be spaces  }
{   around the "+" sign but no spaces between a number and the element       }
{   symbol.                                                                  }
{ - The name of the species in LaTeX sytax follows after the "@" sign.       }

{----------------------------------------------------------------------------}

#DEFVAR

{----------------------------------------------------------------------------}
{--------------------------------- gas phase --------------------------------}
{----------------------------------------------------------------------------}

{------------------------------------- O ------------------------------------}

O1D             =  O                   ; {@O(^1D)}            {O singlet D}
O3P             =  O                   ; {@O(^3P)}            {O triplet P}
O2              = 2O                   ; {@O_2}               {oxygen}
O3              = 3O                   ; {@O_3}               {ozone}

{------------------------------------- H ------------------------------------}

H               =  H                   ; {@H}                 {hydrogen atom}
H2              = 2H                   ; {@H_2}               {hydrogen}
OH              =  H +  O              ; {@OH}                {hydroxyl radical}
HO2             =  H + 2O              ; {@HO_2}              {hydroperoxy radical}
H2O             = 2H +  O              ; {@H_2O}              {water}
H2O2            = 2H + 2O              ; {@H_2O_2}            {hydrogen peroxide}
H2OH2O          = 4H + 2O              ; {@(H_2O)_2}          {water dimer}

{------------------------------------- N ------------------------------------}

N               =            N         ; {@N}                 {nitrogen atom}
N2D             =            N         ; {@N(^2D)}            {N doublet D}
N2              =           2N         ; {@N_2}               {nitrogen}
NH3             = 3H      +  N         ; {@NH_3}              {ammonia}
N2O             =       O + 2N         ; {@N_2O}              {nitrous oxide}
NO              =       O +  N         ; {@NO}                {nitric oxide}
NO2             =      2O +  N         ; {@NO_2}              {nitrogen dioxide}
NO3             =      3O +  N         ; {@NO_3}              {nitrogen trioxide}
N2O5            =      5O + 2N         ; {@N_2O_5}            {dinitrogen pentoxide}
HONO            =  H + 2O +  N         ; {@HONO}              {nitrous acid}
HNO3            =  H + 3O +  N         ; {@HNO_3}             {nitric acid}
HNO4            =  H + 4O +  N         ; {@HNO_4}             {peroxynitric acid}
NH2             = 2H      +  N         ; {@NH_2}              {}
HNO             =  H +  O +  N         ; {@HNO}               {}
NHOH            = 2H +  O +  N         ; {@NHOH}              {}
NH2O            = 2H +  O +  N         ; {@NH_2O}             {}
NH2OH           = 3H +  O +  N         ; {@NH_2OH}            {}
LNITROGEN       =            N         ; {@LNITROGEN}         {lumped N species}

{------------------------------------- C ------------------------------------}

{1C (CHO)}
CH2OO           =   C +  2H +  2O      ; {@CH_2OO}            {MCM: carbonyl oxide - stabilized Criegee Intermediate}
CH2OOA          =   C +  2H +  2O      ; {@CH_2OO^*}          {MCM: carbonyl oxide - excited Criegee Intermediate}
CH3             =   C +  3H            ; {@CH_3}              {methyl radical}
CH3O            =   C +  3H +   O      ; {@CH_3O}             {MCM: methoxy radical}
CH3O2           =   C +  3H +  2O      ; {@CH_3O_2}           {MCM: methylperoxy radical}
CH3OH           =   C +  4H +   O      ; {@CH_3OH}            {MCM: methanol}
CH3OOH          =   C +  4H +  2O      ; {@CH_3OOH}           {MCM: methyl peroxide}
CH4             =   C +  4H            ; {@CH_4}              {MCM: methane}
CO              =   C       +   O      ; {@CO}                {carbon monoxide}
CO2             =   C       +  2O      ; {@CO_2}              {carbon dioxide}
HCHO            =   C +  2H +   O      ; {@HCHO}              {MCM: methanal = formaldehyde}
HCOOH           =   C +  2H +  2O      ; {@HCOOH}             {MCM: formic acid}
HOCH2O2         =   C +  3H +  3O      ; {@HOCH_2O_2}         {hydroxy methyl peroxy radical}
HOCH2OH         =   C +  4H +  2O      ; {@HOCH_2OH}          {dyhydroxy methane}
HOCH2OOH        =   C +  4H +  3O      ; {@HOCH_2OOH}         {hydroxy methyl hydroperoxide}
{1C (CHON)}
CH3NO3          =   C +  3H +  3O +  N ; {@CH_3ONO_2}         {MCM: methylnitrate}
CH3O2NO2        =   C +  3H +  4O +  N ; {@CH_3O_2NO_2}       {MCM: peroxy methylnitrate}
CH3ONO          =   C +  3H +  2O +  N ; {@CH_3ONO}           {methylnitrite}
CN              =   C             +  N ; {@CN}                {}
HCN             =   C +   H       +  N ; {@HCN}               {}
HOCH2O2NO2      =   C +  3H +  5O +  N ; {@HOCH_2O_2NO_2}     {hydroxy methyl peroxy nitrate}
NCO             =   C       +   O +  N ; {@NCO}               {}
{1C (lumped)}
LCARBON         =   C                  ; {@LCARBON}           {lumped carbon}
{2C (CHO)}
C2H2            =  2C +  2H            ; {@C_2H_2}            {MCM: ethyne}
C2H4            =  2C +  4H            ; {@C_2H_4}            {MCM: ethene}
C2H5O2          =  2C +  5H +  2O      ; {@C_2H_5O_2}         {MCM: ethylperoxy radical}
C2H5OH          =  2C +  6H +   O      ; {@C_2H_5OH}          {MCM: ethanol}
C2H5OOH         =  2C +  6H +  2O      ; {@C_2H_5OOH}         {MCM: ethyl hydro peroxide}
C2H6            =  2C +  6H            ; {@C_2H_6}            {MCM: ethane}
CH2CHOH         =  2C +  4H +   O      ; {@CH_2CHOH}          {vinyl alcohol}
CH2CO           =  2C +  2H +   O      ; {@CH2CO}             {ketene}
CH3CHO          =  2C +  4H +   O      ; {@CH_3CHO}           {MCM: acetaldehyde}
CH3CHOHO2       =  2C +  5H +  3O      ; {@CH3CHOHO2}         {}
CH3CHOHOOH      =  2C +  6H +  3O      ; {@CH3CHOHOOH}        {}
CH3CO           =  2C +  3H +  2O      ; {@CH_3C(O)}          {acetyl radical}
CH3CO2H         =  2C +  4H +  2O      ; {@CH_3COOH}          {MCM: acetic acid}
CH3CO3          =  2C +  3H +  3O      ; {@CH_3C(O)OO}        {MCM: peroxy acetyl radical}
CH3CO3H         =  2C +  4H +  3O      ; {@CH_3C(O)OOH}       {MCM: peroxy acetic acid}
ETHGLY          =  2C +  6H +  2O      ; {@ETHGLY}            {MCM: HOCH2CH2OH}
GLYOX           =  2C +  2H +  2O      ; {@GLYOX}             {MCM: CHOCHO = glyoxal}
HCOCH2O2        =  2C +  3H +  3O      ; {@HCOCH_2O_2}        {MCM}
HCOCO           =  2C +   H +  2O      ; {@HCOCO}             {MOM}
HCOCO2H         =  2C +  2H +  3O      ; {@HCOCO_2H}          {MCM: oxoethanoic acid}
HCOCO3          =  2C +   H +  4O      ; {@HCOCO_3}           {MCM}
HCOCO3H         =  2C +  2H +  4O      ; {@HCOCO_3H}          {MCM}
HOCH2CH2O       =  2C +  5H +  2O      ; {@HOCH_2CH_2O}       {MCM: (2-hydroxyethyl)oxidanyl}
HOCH2CH2O2      =  2C +  5H +  3O      ; {@HOCH_2CH_2O_2}     {MCM: (2-hydroxyethyl)dioxidanyl}
HOCH2CHO        =  2C +  4H +  2O      ; {@HOCH_2CHO}         {MCM: glycolaldehyde}
HOCH2CO         =  2C +  3H +  2O      ; {@HOCH2CO}           {}
HOCH2CO2H       =  2C +  4H +  3O      ; {@HOCH_2CO_2H}       {MCM: hydroxyethanoic acid}
HOCH2CO3        =  2C +  3H +  4O      ; {@HOCH_2CO_3}        {MCM}
HOCH2CO3H       =  2C +  4H +  4O      ; {@HOCH_2CO_3H}       {MCM}
HOCHCHO         =  2C +  3H +  2O      ; {@HOCHCHO}           {}
HOOCH2CHO       =  2C +  4H +  3O      ; {@HOOCH2CHO}         {}
HOOCH2CO2H      =  2C +  4H +  4O      ; {@HOOCH2CO2H}        {}
HOOCH2CO3       =  2C +  3H +  5O      ; {@HOOCH_2CO_3}       {}
HOOCH2CO3H      =  2C +  4H +  5O      ; {@HOOCH2CO3H}        {}
HYETHO2H        =  2C +  6H +  3O      ; {@HYETHO2H}          {MCM: HOCH2CH2OOH}
{2C (CHON)}
C2H5NO3         =  2C +  5H +  3O +  N ; {@C_2H_5ONO_2}       {MCM: ethyl nitrate}
C2H5O2NO2       =  2C +  5H +  4O +  N ; {@C_2H_5O_2NO_2}     {ethyl peroxy nitrate}
CH3CN           =  2C +  3H       +  N ; {@CH_3CN}            {acetonitrile}
ETHOHNO3        =  2C +  5H +  4O +  N ; {@ETHOHNO3}          {MCM: HOCH2CH2ONO2}
NCCH2O2         =  2C +  2H +  2O +  N ; {@NCCH_2O_2}         {}
NO3CH2CHO       =  2C +  3H +  4O +  N ; {@NO_3CH2CHO}        {MCM}
NO3CH2CO3       =  2C +  2H +  6O +  N ; {@NO_3CH2CO_3}       {MCM}
NO3CH2PAN       =  2C +  2H +  8O + 2N ; {@NO_3CH2CHO}        {MCM}
PAN             =  2C +  3H +  5O +  N ; {@PAN}               {MCM: CH3C(O)OONO2 = peroxyacetylnitrate}
PHAN            =  2C +  3H +  6O +  N ; {@PHAN}              {MCM: HOCH2C(O)OONO2}
{3C (CHO)}
ACETOL          =  3C +  6H +  2O      ; {@CH_3COCH_2OH}      {MCM: HO-CH2-CO-CH3 = hydroxy acetone}
ALCOCH2OOH      =  3C +  4H +  4O      ; {@HCOCOCH_2OOH}      {MCM}
C2H5CHO         =  3C +  6H +   O      ; {@C_2H_5CHO}         {MCM: propanal}
PROPACID        =  3C +  6H +  2O      ; {@C_2H_5CO_2H}       {MCM}
C2H5CO3         =  3C +  5H +  3O      ; {@C_2H_5CO_3}        {MCM}
PERPROACID      =  3C +  6H +  3O      ; {@C_2H_5CO_3H}       {MCM}
C33CO           =  3C +  2H +  3O      ; {@HCOCOCHO}          {MCM}
C3H6            =  3C +  6H            ; {@C_3H_6}            {MCM: propene}
C3H8            =  3C +  8H            ; {@C_3H_8}            {MCM: propane}
CH3CHCO         =  3C +  4H +   O      ; {@CH3CHCO}           {CH3CHCO}
CH3COCH2O2      =  3C +  5H +  3O      ; {@CH_3COCH_2O_2}     {MCM: peroxyradical from acetone}
CH3COCH3        =  3C +  6H +   O      ; {@CH_3COCH_3}        {MCM: acetone}
CH3COCO2H       =  3C +  4H +  3O      ; {@CH_3COCO_2H}       {MCM: pyruvic acid}
CH3COCO3        =  3C +  3H +  4O      ; {@CH_3COCO_3}        {MCM}
CH3COCO3H       =  3C +  4H +  4O      ; {@CH_3COCO_3H}       {MCM}
CHOCOCH2O2      =  3C +  3H +  4O      ; {@HCOCOCH_2O_2}      {MCM}
HCOCH2CHO       =  3C +  4H +  3O      ; {@HCOCH2CHO}         {MCM}
HCOCH2CO2H      =  3C +  4H +  4O      ; {@HCOCH2CO2H}        {MCM}
HCOCH2CO3       =  3C +  3H +  5O      ; {@HCOCH2CO3}         {MCM}
HCOCH2CO3H      =  3C +  4H +  5O      ; {@HCOCH2CO3H}        {MCM}
HCOCOCH2OOH     =  3C +  4H +  4O      ; {@HCOCOCH_2OOH}      {}
HOC2H4CO2H      =  3C +  6H +  3O      ; {@HOC2H4CO2H}        {MCM: 3-hydroxypropanoic acid}
HOC2H4CO3       =  3C +  5H +  4O      ; {@HOC_2H_4CO_3}      {MCM}
HOC2H4CO3H      =  3C +  6H +  4O      ; {@HOC2H4CO3H}        {MCM}
HOCH2COCH2O2    =  3C +  5H +  4O      ; {@HOCH2COCH2O2}      {}
HOCH2COCH2OOH   =  3C +  6H +  4O      ; {@HOCH2COCH2OOH}     {}
HOCH2COCHO      =  3C +  4H +  3O      ; {@HOCH2COCHO}        {MCM: hydroxypyruvaldehyde}
HYPERACET       =  3C +  6H +  3O      ; {@CH_3COCH_2O_2H}    {MCM: hydroperoxide from CH3COCH2O2}
HYPROPO2        =  3C +  7H +  3O      ; {@HYPROPO2}          {MCM: CH3CH(O2)CH2OH}
HYPROPO2H       =  3C +  8H +  3O      ; {@HYPROPO2H}         {MCM: CH3CH(OOH)CH2OH}
IC3H7O2         =  3C +  7H +  2O      ; {@iC_3H_7O_2}        {MCM: isopropylperoxy radical}
IC3H7OOH        =  3C +  8H +  2O      ; {@iC_3H_7OOH}        {MCM: isopropyl hydro peroxide}
IPROPOL         =  3C +  8H +   O      ; {@IPROPOL}           {MCM: isopropylic alcohol}
MGLYOX          =  3C +  4H +  2O      ; {@MGLYOX}            {MCM: CH3COCHO = methylglyoxal}
NC3H7O2         =  3C +  7H +  2O      ; {@C_3H_7O_2}         {MCM: propylperoxy radical}
NC3H7OOH        =  3C +  8H +  2O      ; {@C_3H_7OOH}         {MCM: propyl hydro peroxide}
NPROPOL         =  3C +  8H +   O      ; {@NPROPOL}           {MCM: n-propylic alcohol}
PROPENOL        =  3C +  6H +   O      ; {@CH_2CHCH_2OH}      {}
{3C (CHO) aromatics}
C32OH13CO       =  3C +  4H +  3O      ; {@C32OH13CO}         {MCM: hydroxymalonaldehyde}
C3DIALO2        =  3C +  3H +  4O      ; {@C3DIALO2}          {MCM}
C3DIALOOH       =  3C +  4H +  4O      ; {@C3DIALOOH}         {MCM}
HCOCOHCO3       =  3C +  3H +  5O      ; {@HCOCOHCO3}         {MCM}
HCOCOHCO3H      =  3C +  4H +  5O      ; {@HCOCOHCO3H}        {MCM}
METACETHO       =  3C +  4H +  3O      ; {@METACETHO}         {MCM: acetic formic anhydride}
{3C (CHON)}
C3PAN1          =  3C +  5H +  6O +  N ; {@C_3PAN1}           {MCM}
C3PAN2          =  3C +  3H +  6O +  N ; {@C_3PAN2}           {MCM}
CH3COCH2O2NO2   =  3C +  5H +  5O +  N ; {@CH_3COCH_2OONO_2}  {CH3-C(O)-CH2-OONO2}
IC3H7NO3        =  3C +  7H +  3O +  N ; {@iC_3H_7ONO_2}      {MCM: isopropyl nitrate}
NC3H7NO3        =  3C +  7H +  3O +  N ; {@C_3H_7ONO_2}       {MCM: propyl nitrate}
NOA             =  3C +  5H +  4O +  N ; {@NOA}               {MCM: CH3-CO-CH2ONO2 = nitro-oxy-acetone}
PPN             =  3C +  5H +  5O +  N ; {@PPN}               {MCM: CH3CH2C(O)OONO2}
PR2O2HNO3       =  3C +  7H +  5O +  N ; {@PR2O2HNO3}         {MCM: CH3-CH(OOH)-CH2ONO2}
PRONO3BO2       =  3C +  6H +  5O +  N ; {@PRONO3BO2}         {MCM: CH3-CH(O2)-CH2ONO2}
PROPOLNO3       =  3C +  7H +  4O +  N ; {@PROPOLNO3}         {MCM: HOCH2-CH(CH3)ONO2)}
{3C (CHON) aromatics}
HCOCOHPAN       =  3C +  3H +  7O +  N ; {@HCOCOHPAN}         {MCM}
{4C (CHO)}
BIACET          =  4C +  6H +  2O      ; {@BIACET}            {MCM: CH3-CO-CO-CH3}
BIACETO2        =  4C +  5H +  4O      ; {@CH_3COCOCH_2O_2}   {MCM}
BIACETOH        =  4C +  6H +  3O      ; {@BIACETOH}          {MCM: CH3-CO-CO-CH2OH}
BIACETOOH       =  4C +  6H +  4O      ; {@CH_3COCOCH_2OOH}   {MCM}
BUT1ENE         =  4C +  8H            ; {@BUT1ENE}           {MCM}
BUT2OLO         =  4C +  8H +  3O      ; {@BUT2OLO}           {MCM}
BUT2OLO2        =  4C +  9H +  2O      ; {@BUT2OLO2}          {MCM}
BUT2OLOOH       =  4C + 10H +  3O      ; {@BUT2OLOOH}         {MCM}
BUTENOL         =  4C +  8H +   O      ; {@BUTENOL}           {CH3CH2CHCHOH}
C312COCO3       =  4C +  3H +  5O      ; {@C312COCO3}         {MCM}
C312COCO3H      =  4C +  4H +  5O      ; {@C312COCO3H}        {MCM}
C3H7CHO         =  4C +  8H +   O      ; {@C_3H_7CHO}         {MCM: n-butanal}
C413COOOH       =  4C +  6H +  4O      ; {@C413COOOH}         {MCM}
C44O2           =  4C +  5H +  5O      ; {@C44O2}             {MCM}
C44OOH          =  4C +  6H +  5O      ; {@C44OOH}            {MCM}
C4CODIAL        =  4C +  4H +  3O      ; {@C4CODIAL}          {MCM}
CBUT2ENE        =  4C +  8H            ; {@CBUT2ENE}          {MCM}
CH3COCHCO       =  4C +  4H +  2O      ; {@CH_3COCHCO}        {}
CH3COCHO2CHO    =  4C +  5H +  4O      ; {@CH_3COCHO_2CHO}    {}
CH3COCOCO2H     =  4C +  6H +  4O      ; {@CH3COCOCO2H}       {}
CH3COOHCHCHO    =  4C +  6H +  3O      ; {@CH_3COOHCHCHO}     {}
CHOC3COO2       =  4C +  5H +  4O      ; {@CHOC3COO2}         {MCM}
CO23C3CHO       =  4C +  4H +  3O      ; {@CH_3COCOCHO}       {MCM}
CO2C3CHO        =  4C +  6H +  2O      ; {@CO2C3CHO}          {MCM: CH3COCH2CHO}
CO2H3CHO        =  4C +  5H +  3O      ; {@CO2H3CHO}          {MCM: CH3-CO-CH(OH)-CHO}
CO2H3CO2H       =  4C +  6H +  5O      ; {@CO2H3CO2H}         {}
CO2H3CO3        =  4C +  5H +  5O      ; {@CO2H3CO3}          {MCM: CH3-CO-CH(OH)-C(O)O2}
CO2H3CO3H       =  4C +  6H +  5O      ; {@CO2H3CO3H}         {MCM: CH3-CO-CH(OH)-C(O)OOH}
EZCH3CO2CHCHO   =  4C +  5H +  3O      ; {@EZCH3CO2CHCHO}     {}
EZCHOCCH3CHO2   =  4C +  5H +  3O      ; {@EZCHOCCH3CHO2}     {}
HCOCCH3CHOOH    =  4C +  6H +  3O      ; {@HCOCCH_3CHOOH}     {}
HCOCCH3CO       =  4C +  4H +  2O      ; {@HCOCCH_3CO}        {}
HCOCO2CH3CHO    =  4C +  5H +  4O      ; {@HCOCO_2CH_3CHO}    {}
HMAC            =  4C +  6H +  2O      ; {@HMAC}              {MCM: HCOC(CH3)CHOH}
HO12CO3C4       =  4C +  8H +  3O      ; {@HO12CO3C4}         {MCM: CH3-CO-CH(OH)-CH2OH}
HVMK            =  4C +  6H +  2O      ; {@HVMK}              {MCM: CH3COCHCHOH = hydroxy vinyl methyl ketone}
IBUTALOH        =  4C +  8H +  2O      ; {@IBUTALOH}          {MCM}
IBUTDIAL        =  4C +  6H +  2O      ; {@IBUTDIAL}          {MCM: HCOC(CH3)CHO}
IBUTOLBO2       =  4C +  9H +  2O      ; {@IBUTOLBO2}         {MCM}
IBUTOLBOOH      =  4C + 10H +  3O      ; {@IBUTOLBOOH}        {}
IC4H10          =  4C + 10H            ; {@iC_4H_<10>}        {MCM: (CH3)3-CH = i-butane}
IC4H9O2         =  4C +  9H +  2O      ; {@IC_4H_9O_2}        {MCM: (CH3)2-CHCH2O2 IC4H9O2}
IC4H9OOH        =  4C + 10H +  2O      ; {@IC_4H_9OOH}        {MCM: (CH3)2-CHCH2OOH MCM: IC4H9OOH}
IPRCHO          =  4C +  8H +   O      ; {@IPRCHO}            {MCM: (CH3)2CHCHO MCM : methylpropanal}
IPRCO3          =  4C +  7H +  3O      ; {@IPRCO3}            {MCM: (CH3)2CHCO3}
IPRHOCO2H       =  4C +  8H +  3O      ; {@IPRHOCO2H}         {MCM}
IPRHOCO3        =  4C +  7H +  4O      ; {@IPRHOCO3}          {MCM}
IPRHOCO3H       =  4C +  8H +  4O      ; {@IPRHOCO3H}         {MCM}
MACO2           =  4C +  5H +  2O      ; {@MACO2}             {}
MACO2H          =  4C +  6H +  2O      ; {@MACO2H}            {MCM: CH2=C(CH3)COOH = methacrylic acid}
MACO3           =  4C +  5H +  3O      ; {@MACO3}             {MCM: CH2=C(CH3)C(O)O2}
MACO3H          =  4C +  6H +  3O      ; {@MACO3H}            {MCM: CH2=C(CH3)C(O)OOH}
MACR            =  4C +  6H +   O      ; {@MACR}              {MCM: CH2=C(CH3)CHO = methacrolein}
MACRO           =  4C +  7H +  3O      ; {@MACRO}             {MCM}
MACRO2          =  4C +  7H +  4O      ; {@MACRO2}            {MCM: HOCH2C(OO)(CH3)CHO}
MACROH          =  4C +  8H +  3O      ; {@MACROH}            {MCM: HOCH2C(OH)(CH3)CHO}
MACROOH         =  4C +  8H +  4O      ; {@MACROOH}           {MCM: HOCH2C(OOH)(CH3)CHO}
MBOOO           =  4C +  8H +  3O      ; {@MBOOO}             {MCM}
MEK             =  4C +  8H +   O      ; {@MEK}               {MCM: CH3-CO-CH2-CH3 = methyl ethyl ketone}
MEPROPENE       =  4C +  8H            ; {@MEPROPENE}         {MCM}
MPROPENOL       =  4C +  8H +   O      ; {@MPROPENOL}         {(CH3)2CCHOH methylpropenol}
MVK             =  4C +  6H +   O      ; {@MVK}               {MCM: CH3-CO-CH=CH2 = methyl vinyl ketone}
NC4H10          =  4C + 10H            ; {@C_4H_<10>}         {MCM: CH3-CH2-CH2-CH3 = n-butane}
PERIBUACID      =  4C +  8H +  3O      ; {@PERIBUACID}        {MCM: (CH3)2CHCO3H}
TBUT2ENE        =  4C +  8H            ; {@TBUT2ENE}          {MCM}
TC4H9O2         =  4C +  9H +  2O      ; {@TC_4H_9O_2}        {MCM: (CH3)3-CO2}
TC4H9OOH        =  4C + 10H +  2O      ; {@TC_4H_9OOH}        {MCM: (CH3)3-COOH}
{4C (CHO) aromatics}
BZFUCO          =  4C +  4H +  4O      ; {@BZFUCO}            {MCM}
BZFUO2          =  4C +  5H +  3O      ; {@BZFUO2}            {MCM}
BZFUONE         =  4C +  4H +  2O      ; {@BZFUONE}           {MCM: 2(5H)-furanone}
BZFUOOH         =  4C +  6H +  5O      ; {@BZFUOOH}           {MCM}
CO14O3CHO       =  4C +  4H +  4O      ; {@CO14O3CHO}         {MCM}
CO14O3CO2H      =  4C +  4H +  5O      ; {@CO14O3CO2H}        {MCM}
CO2C4DIAL       =  4C +  2H +  4O      ; {@CO2C4DIAL}         {MCM: 2,3-dioxosuccinaldehyde}
EPXC4DIAL       =  4C +  4H +  3O      ; {@EPXC4DIAL}         {MCM}
EPXDLCO2H       =  4C +  4H +  4O      ; {@EPXDLCO2H}         {MCM}
EPXDLCO3        =  4C +  3H +  5O      ; {@EPXDLCO3}          {MCM}
EPXDLCO3H       =  4C +  4H +  5O      ; {@EPXDLCO3H}         {MCM}
HOCOC4DIAL      =  4C +  4H +  4O      ; {@HOCOC4DIAL}        {MCM: 2-hydroxy-3-oxosuccinaldehyde}
MALANHY         =  4C +  2H +  3O      ; {@MALANHY}           {MCM: maleic anhydride}
MALANHYO2       =  4C +  3H +  6O      ; {@MALANHYO2}         {MCM}
MALANHYOOH      =  4C +  4H +  6O      ; {@MALANHYOOH}        {MCM}
MALDALCO2H      =  4C +  4H +  3O      ; {@MALDALCO2H}        {MCM: 4-oxo-2-butenoic acid}
MALDALCO3H      =  4C +  4H +  4O      ; {@MALDALCO3H}        {MCM}
MALDIAL         =  4C +  4H +  2O      ; {@MALDIAL}           {MCM: 2-butenedial}
MALDIALCO3      =  4C +  3H +  4O      ; {@MALDIALCO3}        {MCM}
MALDIALO2       =  4C +  5H +  5O      ; {@MALDIALO2}         {MCM}
MALDIALOOH      =  4C +  6H +  5O      ; {@MALDIALOOH}        {MCM}
MALNHYOHCO      =  4C +  2H +  5O      ; {@MALNHYOHCO}        {MCM}
MECOACEOOH      =  4C +  6H +  5O      ; {@MECOACEOOH}        {MCM}
MECOACETO2      =  4C +  5H +  5O      ; {@MECOACETO2}        {MCM}
{4C (CHON)}
BUT2OLNO3       =  4C +  9H +  5O +  N ; {@BUT2OLNO3}         {MCM}
C312COPAN       =  4C +  3H +  7O +  N ; {@C312COPAN}         {MCM}
C4PAN5          =  4C +  7H +  6O +  N ; {@C4PAN5}            {MCM}
IBUTOLBNO3      =  4C +  9H +  4O +  N ; {@IBUTOLBNO3}        {MCM}
IC4H9NO3        =  4C +  9H +  3O +  N ; {@IC4H9NO3}          {MCM}
MACRNO3         =  4C +  7H +  5O +  N ; {@MACRNO3}           {MCM}
MPAN            =  4C +  5H +  5O +  N ; {@MPAN}              {MCM: CH2=C(CH3)C(O)OONO2 = peroxymethacryloyl nitrate = peroxymethacrylic nitric anhydride}
MVKNO3          =  4C +  7H +  5O +  N ; {@MVKNO3}            {MCM}
PIPN            =  4C +  7H +  5O +  N ; {@PIPN}              {MCM: (CH3)2CHCO3} 
TC4H9NO3        =  4C +  9H +  3O +  N ; {@TC4H9NO3}          {MCM}
{4C (CHON) aromatics}
EPXDLPAN        =  4C +  3H +  7O +  N ; {@EPXDLPAN}          {MCM}
MALDIALPAN      =  4C +  3H +  6O +  N ; {@MALDIALPAN}        {MCM}
NBZFUO2         =  4C +  4H +  7O +  N ; {@NBZFUO2}           {MCM}
NBZFUONE        =  4C +  3H +  6O +  N ; {@NBZFUONE}          {MCM}
NBZFUOOH        =  4C +  5H +  7O +  N ; {@NBZFUOOH}          {MCM}
NC4DCO2H        =  4C +  3H +  5O +  N ; {@NC4DCO2H}          {MCM}
{4C (CHO) (lumped)}
LBUT1ENO2       =  4C +  9H +  2O      ; {@LBUT1ENO2}         {HO3C4O2 + NBUTOLAO2}
LBUT1ENOOH      =  4C + 10H +  3O      ; {@LBUT1ENOOH}        {HO3C4OOH + NBUTOLAOOH}
LC4H9O2         =  4C +  9H +  2O      ; {@LC_4H_9O_2}        {CH3-CH2-CH(O2)-CH3 + CH3-CH2-CH2-CH2O2 = NC4H9O2 + SC4H9O2}
LC4H9OOH        =  4C + 10H +  2O      ; {@LC_4H_9OOH}        {CH3-CH2-CH(OOH)-CH3 + CH3-CH2-CH2-CH2OOH = NC4H9OOH + SC4H9OOH}
LHMVKABO2       =  4C +  7H +  4O      ; {@LHMVKABO2}         {HOCH2-CH(O2)-CO-CH3 + CH2(O2)-CH(OH)-CO-CH3}
LHMVKABOOH      =  4C +  8H +  4O      ; {@LHMVKABOOH}        {HOCH2-CH(OOH)-CO-CH3 + CH2(OOH)-CH(OH)-CO-CH3}
LMEKO2          =  4C +  7H +  3O      ; {@LMEKO2}            {CH3-CO-CH2-CH2-OO + CH3-CO-CH(O2)-CH3}
LMEKOOH         =  4C +  8H +  3O      ; {@LMEKOOH}           {CH3-CO-CH2-CH2-OOH + CH3-CO-CH(OOH)-CH3}
{4C (CHON) (lumped)}
LBUT1ENNO3      =  4C +  9H +  5O +  N ; {@LBUT1ENNO3}        {HO3C4NO3 + NBUTOLANO3}
LC4H9NO3        =  4C +  9H +  3O +  N ; {@LC4H9NO3}          {NC4H9NO3 + SC4H9NO3}
LMEKNO3         =  4C +  7H +  5O +  N ; {@LMEKNO3}           {CH3-CO-CH2-CH2-ONO2 + CH3-CO-CH(ONO2)-CH3}
{5C (CHO)}
C1ODC2O2C4OD    =  5C +  7H +  4O      ; {@C1ODC2O2C4OD}      {}
C1ODC2O2C4OOH   =  5C +  9H +  5O      ; {@C1ODC2O2C4OOH}     {}
C1ODC2OOHC4OD   =  5C +  8H +  4O      ; {@C1ODC2OOHC4OD}     {}
C1ODC3O2C4OOH   =  5C +  9H +  5O      ; {@C1ODC3O2C4OOH}     {}
C1OOHC2O2C4OD   =  5C +  9H +  5O      ; {@C1OOHC2O2C4OD}     {}
C1OOHC2OOHC4OD  =  5C + 10H +  5O      ; {@C1OOHC2OOHC4OD}    {}
C1OOHC3O2C4OD   =  5C +  9H +  5O      ; {@C1OOHC3O2C4OD}     {}
C4MDIAL         =  5C +  6H +  2O      ; {@C4MDIAL}           {MCM: 2-methyl-butenedial}
C511O2          =  5C +  7H +  4O      ; {@C511O2}            {MCM}
C511OOH         =  5C +  8H +  4O      ; {@C511OOH}           {MCM}
C512O2          =  5C +  7H +  4O      ; {@C512O2}            {MCM}
C512OOH         =  5C +  8H +  4O      ; {@C512OOH}           {MCM}
C513CO          =  5C +  6H +  4O      ; {@C513CO}            {MCM}
C513O2          =  5C +  7H +  5O      ; {@C513O2}            {MCM}
C513OOH         =  5C +  8H +  5O      ; {@C513OOH}           {MCM}
C514O2          =  5C +  7H +  4O      ; {@C514O2}            {MCM}
C514OOH         =  5C +  8H +  4O      ; {@C514OOH}           {MCM}
C59O2           =  5C +  9H +  5O      ; {@C59O2}             {MCM: HOCH2-CO-C(CH3)(O2)-CH2OH}
C59OOH          =  5C + 10H +  5O      ; {@C59OOH}            {MCM: HOCH2-CO-C(CH3)(OOH)-CH2OH}
C5H8            =  5C +  8H            ; {@C_5H_8}            {MCM: CH2=C(CH3)CH=CH2 = isoprene}
CHOC3COCO3      =  5C +  5H +  5O      ; {@CHOC3COCO3}        {MCM}
CHOC3COOOH      =  5C +  6H +  4O      ; {@CHOC3COOOH}        {MCM}
CO13C4CHO       =  5C +  6H +  3O      ; {@CO13C4CHO}         {MCM}
CO23C4CHO       =  5C +  6H +  3O      ; {@CO23C4CHO}         {MCM}
CO23C4CO3       =  5C +  5H +  5O      ; {@CO23C4CO3}         {MCM}
CO23C4CO3H      =  5C +  6H +  5O      ; {@CO23C4CO3H}        {MCM}
DB1O            =  5C +  9H +  3O      ; {@DB1O2}             {Alkoxy radical which undergoes the double H-shift predicted by T. Dibble and confirmed by F. Paulot}
DB1O2           =  5C +  9H +  4O      ; {@DB1O2}             {Peroxy radical with a vinyl alcohol part}
DB1OOH          =  5C + 10H +  4O      ; {@DB1OOH}            {}
DB2O2           =  5C +  9H +  5O      ; {@DB1O2}             {}
DB2OOH          =  5C + 10H +  5O      ; {@DB2OOH}            {}
HCOC5           =  5C +  8H +  2O      ; {@HCOC5}             {MCM: HOCH2-CO-C(CH3)=CH2}
ISOPAOH         =  5C + 10H +  2O      ; {@ISOPAOH}           {MCM: HOCH2-C(CH3)=CH-CH2OH}
ISOPBO2         =  5C +  9H +  3O      ; {@ISOPBO2}           {MCM: HOCH2-C(CH3)(O2)-CH=CH2}
ISOPBOH         =  5C + 10H +  2O      ; {@ISOPBOH}           {MCM: HOCH2-C(CH3)(OH)-CH=CH2}
ISOPBOOH        =  5C + 10H +  3O      ; {@ISOPBOOH}          {MCM: HOCH2-C(CH3)(OOH)-CH=CH2}
ISOPDO2         =  5C +  9H +  3O      ; {@ISOPDO2}           {MCM: CH2=C(CH3)CH(O2)-CH2OH}
ISOPDOH         =  5C + 10H +  2O      ; {@ISOPDOH}           {MCM: CH2=C(CH3)CH(OH)-CH2OH}
ISOPDOOH        =  5C + 10H +  3O      ; {@ISOPDOOH}          {MCM: CH2=C(CH3)CH(OOH)-CH2OH}
MBO             =  5C + 10H +   O      ; {@MBO}               {MCM: 2-methyl-3-buten-2-ol}
MBOACO          =  5C + 10H +  3O      ; {@MBOACO}            {MCM}
MBOCOCO         =  5C +  8H +  3O      ; {@MBOCOCO}           {MCM}
ME3FURAN        =  5C +  6H +   O      ; {@3METHYLFURAN}      {3-methyl-furan}
{5C aromatics (CHO)}
ACCOMECHO       =  5C +  6H +  4O      ; {@ACCOMECHO}         {MCM}
ACCOMECO3       =  5C +  5H +  6O      ; {@ACCOMECO3}         {MCM}
ACCOMECO3H      =  5C +  6H +  6O      ; {@ACCOMECO3H}        {MCM}
C24O3CCO2H      =  5C +  6H +  5O      ; {@C24O3CCO2H}        {MCM}
C4CO2DBCO3      =  5C +  3H +  5O      ; {@C4CO2DBCO3}        {MCM}
C4CO2DCO3H      =  5C +  4H +  5O      ; {@C4CO2DCO3H}        {MCM}
C5134CO2OH      =  5C +  6H +  4O      ; {@C5134CO2OH}        {MCM: 2-hydroxy-3,4-dioxopentanal}
C54CO           =  5C +  4H +  4O      ; {@C54CO}             {MCM: 2,3,4-trioxopentanal}
C5CO14O2        =  5C +  5H +  4O      ; {@C5CO14O2}          {MCM}
C5CO14OH        =  5C +  6H +  3O      ; {@C5CO14OH}          {MCM: 4-oxo-2-pentenoic acid}
C5CO14OOH       =  5C +  6H +  4O      ; {@C5CO14OOH}         {MCM}
C5DIALCO        =  5C +  4H +  3O      ; {@C5DIALCO}          {MCM}
C5DIALO2        =  5C +  5H +  4O      ; {@C5DIALO2}          {MCM}
C5DIALOOH       =  5C +  6H +  4O      ; {@C5DIALOOH}         {MCM}
C5DICARB        =  5C +  6H +  2O      ; {@C5DICARB}          {MCM: 4-oxo-2-pentenal}
C5DICARBO2      =  5C +  7H +  5O      ; {@C5DICARBO2}        {MCM: carboxy(hydroxy)acetate}
C5DICAROOH      =  5C +  8H +  5O      ; {@C5DICAROOH}        {MCM}
MC3ODBCO2H      =  5C +  6H +  3O      ; {@MC3ODBCO2H}        {MCM}
MMALANHY        =  5C +  4H +  3O      ; {@MMALANHY}          {MCM: 3-methyl-2,5-furandione}
MMALANHYO2      =  5C +  5H +  6O      ; {@MMALANHYO2}        {MCM}
MMALNHYOOH      =  5C +  6H +  6O      ; {@MMALNHYOOH}        {MCM}
TLFUO2          =  5C +  7H +  5O      ; {@TLFUO2}            {MCM}
TLFUONE         =  5C +  6H +  2O      ; {@TLFUONE}           {MCM: 5-methyl-2(5H)-furanone}
TLFUOOH         =  5C +  8H +  5O      ; {@TLFUOOH}           {MCM}
{5C (CHON)}
C4MCONO3OH      =  5C +  9H +  5O +  N ; {@C4MCONO3OH}        {MCM}
C514NO3         =  5C +  7H +  5O +  N ; {@C514NO3}           {MCM}
C5PAN9          =  5C +  5H +  7O +  N ; {@C5PAN9}            {MCM}
CHOC3COPAN      =  5C +  5H +  5O +  N ; {@CHOC3COPAN}        {MCM}
DB1NO3          =  5C +  9H +  6O +  N ; {@DB1NO3}            {}
ISOPBDNO3O2     =  5C + 10H +  7O +  N ; {@ISOPBDNO3O2}       {}
ISOPBNO3        =  5C +  9H +  4O +  N ; {@ISOPBNO3}          {MCM: HOCH2-C(CH3)(ONO2)-CH=CH2}
ISOPDNO3        =  5C +  9H +  4O +  N ; {@ISOPDNO3}          {MCM: CH2=C(CH3)CH(ONO2)-CH2OH}
NC4CHO          =  5C +  7H +  4O +  N ; {@NC4CHO}            {MCM: O2NOCH2-C(CH3)=CH-CHO}
NC4OHCO3        =  5C +  8H +  6O +  N ; {@NC4OHCO3}          {MCM}
NC4OHCO3H       =  5C +  9H +  6O +  N ; {@NC4OHCO3H}         {MCM}
NC4OHCPAN       =  5C +  8H +  8O + 2N ; {@NC4OHCPAN}         {MCM}
NISOPO2         =  5C +  8H +  5O +  N ; {@NISOPO2}           {MCM: O2NOCH2-C(CH3)=CH-CH2O2}
NISOPOOH        =  5C +  9H +  5O +  N ; {@NISOPOOH}          {MCM: O2NOCH2-C(CH3)=CH-CH2OOH}
NMBOBCO         =  5C +  9H +  5O +  N ; {@NMBOBCO}           {MCM}
{5C aromatics (CHON)}
ACCOMEPAN       =  5C +  5H +  6O +  N ; {@ACCOMEPAN}         {MCM}
C4CO2DBPAN      =  5C +  3H +  7O +  N ; {@C4CO2DBPAN}        {MCM}
C5COO2NO2       =  5C +  5H +  6O +  N ; {@C5COO2NO2}         {MCM}
NC4MDCO2H       =  5C +  5H +  5O +  N ; {@NC4MDCO2H N}       {MCM}
NTLFUO2         =  5C +  6H +  7O +  N ; {@NTLFUO2}           {MCM}
NTLFUOOH        =  5C +  7H +  6O +  N ; {@NTLFUOOH}          {MCM}
{5C (CHO) (lumped)}
LC578O2         =  5C +  9H +  5O      ; {@LC578O2}           {HOCH2-CH(OH)C(CH3)(O2)-CHO + HOCH2-C(CH3)(O2)-CH(OH)-CHO}
LC578OOH        =  5C + 10H +  5O      ; {@LC578OOH}          {HOCH2-CH(OH)C(CH3)(OOH)-CHO + HOCH2-C(CH3)(OOH)-CH(OH)-CHO}
LDISOPACO       =  5C +  9H +  2O      ; {@LISOPACO}          {}
LDISOPACO2      =  5C +  9H +  3O      ; {@LDISOPACO2}        {}
LHC4ACCHO       =  5C +  8H +  2O      ; {@LHC4ACCHO}         {HOCH2-C(CH3)=CH-CHO + HOCH2-CH=C(CH3)-CHO}
LHC4ACCO2H      =  5C +  8H +  3O      ; {@LHC4ACCO2H}        {HOCH2-C(CH3)=CH-C(O)OH + HOCH2-CH=C(CH3)-C(O)OH}
LHC4ACCO3       =  5C +  7H +  4O      ; {@LHC4ACCO3}         {HOCH2-C(CH3)=CH-C(O)O2 + HOCH2-CH=C(CH3)-C(O)O2}
LHC4ACCO3H      =  5C +  8H +  4O      ; {@LHC4ACCO3H}        {HOCH2-C(CH3)=CH-C(O)OOH + HOCH2-CH=C(CH3)-C(O)OOH}
LIEPOX          =  5C + 10H +  3O      ; {@LIEPOX}            {epoxydiol}
LISOPAB         =  5C +  9H +   O      ; {@LISOPAB}           {}
LISOPACO        =  5C +  9H +  2O      ; {@LISOPACO}          {HOCH2-C(CH3)=CH-CH2O + HOCH2-CH=C(CH3)-CH2O}
LISOPACO2       =  5C +  9H +  3O      ; {@LISOPACO2}         {HOCH2-C(CH3)=CH-CH2O2 + HOCH2-CH=C(CH3)-CH2O2}
LISOPACOOH      =  5C + 10H +  3O      ; {@LISOPACOOH}        {HOCH2-C(CH3)=CH-CH2OOH + HOCH2-CH=C(CH3)-CH2OOH}
LISOPCD         =  5C +  9H +   O      ; {@LISOPCD}           {}
LISOPEFO        =  5C +  9H +  2O      ; {@LISOPEFO}          {}
LISOPEFO2       =  5C +  9H +  3O      ; {@LISOPEFO2}         {}
LMBOABO2        =  5C + 11H +  4O      ; {@LMBOABO2}          {}
LMBOABOOH       =  5C + 12H +  4O      ; {@LMBOABOOH}         {}
LME3FURANO2     =  5C +  7H +  4O      ; {@L3METHYLFURANO2}   {hydroxy-3-methyl-furan peroxy radical}
LZCO3C23DBCOD   =  5C +  5H +  4O      ; {@LZCO3C23DBCOD}     {}
LZCO3HC23DBCOD  =  5C +  6H +  4O      ; {@LZCO3HC23DBCOD}    {C5PACALD1 + C5PACALD2}
LZCODC23DBCOOH  =  5C +  8H +  3O      ; {@LZCODC23DBCOOH}    {C5HPALD1 + C5HPALD2}
{5C (CHON) (lumped)}
LC5PAN1719      =  5C +  7H +  6O +  N ; {@LC5PAN1719}        {HOCH2-C(CH3)=CH-C(O)OONO2 + HOCH2-CH=C(CH3)C(O)OONO2}
LISOPACNO3      =  5C +  9H +  4O +  N ; {@LISOPACNO3}        {HOCH2-C(CH3)=CH-CH2ONO2 + HOCH2-CH=C(CH3)-CH2ONO2}
LISOPACNO3O2    =  5C + 10H +  7O +  N ; {@LISOPACNO3O2}      {RO2 resulting from OH-addition to ISOPANO3 and ISOPCNO3}
LMBOABNO3       =  5C + 11H +  5O +  N ; {@LMBOABNO3}         {}
LNISO3          =  5C             +  N ; {@LNISO3}            {C510O2+NC4CO3 = CHO-CH(OH)-C(CH3)(O2)-CH2ONO2 + O2NOCH2-C(CH3)=CH-C(O)O2}
LNISOOH         =  5C             +  N ; {@LNISOOH}           {CHO-CH(OH)-C(CH3)(OOH)-CH2ONO2 + O2NOCH2-C(CH3)=CH-C(O)OOH}
LNMBOABO2       =  5C +  9H +  6O +  N ; {@LNMBOABO2}         {}
LNMBOABOOH      =  5C + 10H +  6O +  N ; {@LNMBOABOOH}        {}
LZCPANC23DBCOD  =  5C +  5H +  6O +  N ; {@LZCPANC23DBCOD}    {}
{6C (CHO)}
C614CO          =  6C +  8H +  4O      ; {@C614CO}            {MCM}
C614O2          =  6C +  9H +  5O      ; {@C614O2}            {MCM}
C614OOH         =  6C + 10H +  5O      ; {@C614OOH}           {MCM}
CO235C5CHO      =  6C +  6H +  4O      ; {@CO235C5CHO}        {MCM}
CO235C6O2       =  6C +  7H +  5O      ; {@CO235C6O2}         {MCM}
CO235C6OOH      =  6C +  8H +  5O      ; {@CO235C6OOH}        {MCM}
{C6 (CHO) aromatics}
BENZENE         =  6C +  6H            ; {@BENZENE}           {MCM}
BZBIPERO2       =  6C +  7H +  5O      ; {@BZBIPERO2}         {MCM}
BZBIPEROOH      =  6C +  8H +  5O      ; {@BZBIPEROOH}        {MCM}
BZEMUCCO        =  6C +  6H +  5O      ; {@BZEMUCCO}          {MCM}
BZEMUCCO2H      =  6C +  6H +  4O      ; {@BZEMUCCO2H}        {MCM}
BZEMUCCO3       =  6C +  5H +  5O      ; {@BZEMUCCO3}         {MCM}
BZEMUCCO3H      =  6C +  6H +  5O      ; {@BZEMUCCO3H}        {MCM}
BZEMUCO2        =  6C +  7H +  6O      ; {@BZEMUCO2}          {MCM}
BZEMUCOOH       =  6C +  8H +  6O      ; {@BZEMUCOOH}         {MCM}
BZEPOXMUC       =  6C +  6H +  3O      ; {@BZEPOXMUC}         {MCM}
BZOBIPEROH      =  6C +  6H +  4O      ; {@BZOBIPEROH}        {MCM}
C5CO2DBCO3      =  6C +  5H +  5O      ; {@C5CO2DBCO3}        {MCM}
C5CO2DCO3H      =  6C +  6H +  5O      ; {@C5CO2DCO3H}        {MCM}
C5CO2OHCO3      =  6C +  5H +  6O      ; {@C5CO2OHCO3}        {MCM}
C5COOHCO3H      =  6C +  6H +  6O      ; {@C5COOHCO3H}        {MCM}
C6125CO         =  6C +  6H +  3O      ; {@C6125CO}           {MCM: 2,5-dioxo-3-hexenal}
C615CO2O2       =  6C +  7H +  4O      ; {@C615CO2O2}         {MCM}
C615CO2OOH      =  6C +  8H +  4O      ; {@C615CO2OOH}        {MCM}
C6CO4DB         =  6C +  4H +  4O      ; {@C6CO4DB}           {MCM}
C6H5O           =  6C +  5H +   O      ; {@C6H5O}             {MCM: phenyloxidanyl}
C6H5O2          =  6C +  5H +  2O      ; {@C6H5O2}            {MCM}
C6H5OOH         =  6C +  6H +  2O      ; {@C6H5OOH}           {MCM: phenyl hydroperoxide}
CATEC1O         =  6C +  5H +  2O      ; {@CATEC1O}           {MCM: 2-λ1-oxidanylphenol}
CATEC1O2        =  6C +  5H +  3O      ; {@CATEC1O2}          {MCM}
CATEC1OOH       =  6C +  6H +  3O      ; {@CATEC1OOH}         {MCM}
CATECHOL        =  6C +  4H +  2O      ; {@CATECHOL}          {MCM: catechol}
CPDKETENE       =  6C +  4H +   O      ; {@CPDKETENE}         {hv nitrophenol: cyclopentadiene ketene (Luc Vereecken's prediction)}
PBZQCO          =  6C +  4H +  4O      ; {@PBZQCO}            {MCM}
PBZQO2          =  6C +  5H +  5O      ; {@PBZQO2}            {MCM}
PBZQONE         =  6C +  4H +  2O      ; {@PBZQONE}           {MCM: 1,4-benzoquinone}
PBZQOOH         =  6C +  6H +  5O      ; {@PBZQOOH}           {MCM}
PHENO2          =  6C +  7H +  6O      ; {@PHENO2}            {MCM}
PHENOL          =  6C +  6H +   O      ; {@PHENOL}            {MCM}
PHENOOH         =  6C +  8H +  6O      ; {@PHENOOH}           {MCM}
{6C (CHON)}
C614NO3         =  6C +  9H +  6O +  N ; {@C614NO3}           {MCM}
{C6 (CHON) aromatics}
BZBIPERNO3      =  6C +  7H +  6O +  N ; {@BZBIPERNO3}        {MCM}
BZEMUCNO3       =  6C +  7H +  7O +  N ; {@BZEMUCNO3}         {MCM}
BZEMUCPAN       =  6C +  5H +  7O +  N ; {@BZEMUCPAN}         {MCM}
C5CO2DBPAN      =  6C +  5H +  7O +  N ; {@C5CO2DBPAN}        {MCM}
C5CO2OHPAN      =  6C +  5H +  8O +  N ; {@C5CO2OHPAN}        {MCM}
DNPHEN          =  6C +  4H +  5O + 2N ; {@DNPHEN}            {MCM: 2,4-dinitrophenol}
DNPHENO2        =  6C +  5H + 10O + 2N ; {@DNPHENO2}          {MCM}
DNPHENOOH       =  6C +  6H + 10O + 2N ; {@DNPHENOOH}         {MCM}
HOC6H4NO2       =  6C +  5H +  3O +  N ; {@HOC6H4NO2}         {MCM: 2-nitrophenol}
NBZQO2          =  6C +  4H +  7O +  N ; {@NBZQO2}            {MCM}
NBZQOOH         =  6C +  5H +  7O +  N ; {@NBZQOOH}           {MCM}
NCATECHOL       =  6C +  5H +  4O +  N ; {@NCATECHOL}         {MCM}
NCATECO2        =  6C +  6H +  9O +  N ; {@NCATECO2}          {MCM}
NCATECOOH       =  6C +  7H +  9O +  N ; {@NCATECOOH}         {MCM}
NCPDKETENE      =  6C +  3H +  3O +  N ; {@NCPDKETENE}        {hv nitrophenol: cyclopentadiene ketene (Luc Vereecken's prediction)}
NDNPHENO2       =  6C +  4H + 12O + 3N ; {@NDNPHENO2}         {MCM}
NDNPHENOOH      =  6C +  5H + 12O + 3N ; {@NDNPHENOOH}        {MCM}
NNCATECO2       =  6C +  5H + 11O + 2N ; {@NNCATECO2}         {MCM}
NNCATECOOH      =  6C +  6H + 11O + 2N ; {@NNCATECOOH}        {MCM}
NPHEN1O         =  6C +  4H +  3O +  N ; {@NPHEN1O}           {MCM}
NPHEN1O2        =  6C +  4H +  4O +  N ; {@NPHEN1O2}          {MCM}
NPHEN1OOH       =  6C +  5H +  4O +  N ; {@NPHEN1OOH}         {MCM}
NPHENO2         =  6C +  6H +  8O +  N ; {@NPHENO2}           {MCM}
NPHENOOH        =  6C +  7H +  8O +  N ; {@NPHENOOH}          {MCM}
{7C (CHO)}
C235C6CO3H      =  7C +  8H +  6O      ; {@C235C6CO3H}        {MCM}
C716O2          =  7C +  9H +  5O      ; {@C716O2}            {MCM}
C716OOH         =  7C + 10H +  5O      ; {@C716OOH}           {MCM}
C721O2          =  7C + 11H +  4O      ; {@C721O2}            {MCM}
C721OOH         =  7C + 12H +  4O      ; {@C721OOH}           {MCM}
C722O2          =  7C + 11H +  5O      ; {@C722O2}            {MCM}
C722OOH         =  7C + 12H +  5O      ; {@C722OOH}           {MCM}
CO235C6CHO      =  7C +  8H +  4O      ; {@CO235C6CHO}        {MCM}
CO235C6CO3      =  7C +  7H +  6O      ; {@CO235C6CO3}        {MCM}
MCPDKETENE      =  7C +  6H +  2O      ; {@MCPDKETENE}        {hv nitrophenol: cyclopentadiene ketene (Luc Vereecken's prediction)}
ROO6R3O         =  7C + 11H +  4O      ; {@ROO6R3O}           {from ref3019}
ROO6R3O2        =  7C + 11H +  5O      ; {@ROO6R3O2}          {ROO6R3OO from ref3019}
ROO6R5O2        =  7C + 11H +  7O      ; {@ROO6R5O2}          {ROO6R5OO from ref3019}
{C7 (CHO) aromatics}
BENZAL          =  7C +  6H +   O      ; {@BENZAL}            {MCM}
C6CO2OHCO3      =  7C +  7H +  6O      ; {@C6CO2OHCO3}        {MCM}
C6COOHCO3H      =  7C +  8H +  6O      ; {@C6COOHCO3H}        {MCM}
C6H5CH2O2       =  7C +  7H +  2O      ; {@C6H5CH2O2}         {MCM: benzyldioxidanyl}
C6H5CH2OOH      =  7C +  8H +  2O      ; {@C6H5CH2OOH}        {MCM: benzyl hydroperoxide}
C6H5CO3         =  7C +  5H +  3O      ; {@C6H5CO3}           {MCM}
C6H5CO3H        =  7C +  6H +  3O      ; {@C6H5CO3H}          {MCM: perbenzoic acid}
C7CO4DB         =  7C +  6H +  4O      ; {@C7CO4DB}           {MCM}
CRESO2          =  7C +  9H +  6O      ; {@CRESO2}            {MCM}
CRESOL          =  7C +  8H +   O      ; {@CRESOL}            {MCM: 2-methylphenol}
CRESOOH         =  7C + 10H +  6O      ; {@CRESOOH}           {MCM}
MCATEC1O        =  7C +  7H +  2O      ; {@MCATEC1O}          {MCM}
MCATEC1O2       =  7C +  7H +  3O      ; {@MCATEC1O2}         {MCM}
MCATEC1OOH      =  7C +  8H +  3O      ; {@MCATEC1OOH}        {MCM}
MCATECHOL       =  7C +  8H +  2O      ; {@MCATECHOL}         {MCM: 3-methylcatechol}
OXYL1O2         =  7C +  7H +  2O      ; {@OXYL1O2}           {MCM: 1-methyl-2-(oxo-λ3-oxidanyl)benzene}
OXYL1OOH        =  7C +  8H +  2O      ; {@OXYL1OOH}          {MCM}
PHCOOH          =  7C +  6H +  2O      ; {@PHCOOH}            {MCM: benzoic acid}
PTLQCO          =  7C +  6H +  4O      ; {@PTLQCO}            {MCM}
PTLQO2          =  7C +  7H +  5O      ; {@PTLQO2}            {MCM}
PTLQONE         =  7C +  6H +  2O      ; {@PTLQONE}           {MCM: 2-methyl-1,4-benzoquinone}
PTLQOOH         =  7C +  8H +  5O      ; {@PTLQOOH}           {MCM}
TLBIPERO2       =  7C +  9H +  5O      ; {@TLBIPERO2}         {MCM}
TLBIPEROOH      =  7C + 10H +  5O      ; {@TLBIPEROOH}        {MCM}
TLEMUCCO        =  7C +  8H +  5O      ; {@TLEMUCCO}          {MCM}
TLEMUCCO2H      =  7C +  8H +  4O      ; {@TLEMUCCO2H}        {MCM}
TLEMUCCO3       =  7C +  7H +  5O      ; {@TLEMUCCO3}         {MCM}
TLEMUCCO3H      =  7C +  8H +  5O      ; {@TLEMUCCO3H}        {MCM}
TLEMUCO2        =  7C +  9H +  6O      ; {@TLEMUCO2}          {MCM}
TLEMUCOOH       =  7C + 10H +  6O      ; {@TLEMUCOOH}         {MCM}
TLEPOXMUC       =  7C +  8H +  3O      ; {@TLEPOXMUC}         {MCM}
TLOBIPEROH      =  7C +  8H +  4O      ; {@TLOBIPEROH}        {MCM}
TOL1O           =  7C +  7H +   O      ; {@TOL1O}             {MCM: (2-methylphenyl)oxidanyl}
TOLUENE         =  7C +  8H            ; {@TOLUENE}           {MCM}
{7C (CHON)}
C7PAN3          =  7C +  7H +  8O +  N ; {@C7PAN3}            {MCM}
{C7 (CHON) aromatics}
C6CO2OHPAN      =  7C +  7H +  8O +  N ; {@C6CO2OHPAN}        {MCM}
C6H5CH2NO3      =  7C +  7H +  3O +  N ; {@C6H5CH2NO3}        {MCM: benzyl nitrate}
DNCRES          =  7C +  6H +  5O + 2N ; {@DNCRES}            {MCM: 2-methyl-4,6-dinitrophenol}
DNCRESO2        =  7C +  7H + 10O + 2N ; {@DNCRESO2}          {MCM}
DNCRESOOH       =  7C +  8H + 10O + 2N ; {@DNCRESOOH}         {MCM}
MNCATECH        =  7C +  7H +  4O +  N ; {@MNCATECH}          {MCM: 3-methyl-6-nitro-1,2-benzenediol}
MNCATECO2       =  7C +  8H +  9O +  N ; {@MNCATECO2}         {MCM}
MNCATECOOH      =  7C +  9H +  9O +  N ; {@MNCATECOOH}        {MCM}
MNCPDKETENE     =  7C +  5H +  3O +  N ; {@MNCPDKETENE}       {hv nitrophenol: cyclopentadiene ketene (Luc Vereecken's prediction)}
MNNCATCOOH      =  7C +  8H + 11O + 2N ; {@MNNCATCOOH}        {MCM}
MNNCATECO2      =  7C +  7H + 11O + 2N ; {@MNNCATECO2}        {MCM}
NCRES1O         =  7C +  6H +  3O +  N ; {@NCRES1O}           {MCM}
NCRES1O2        =  7C +  6H +  4O +  N ; {@NCRES1O2}          {MCM}
NCRES1OOH       =  7C +  7H +  4O +  N ; {@NCRES1OOH}         {MCM}
NCRESO2         =  7C +  8H +  8O +  N ; {@NCRESO2}           {MCM}
NCRESOOH        =  7C +  9H +  8O +  N ; {@NCRESOOH}          {MCM}
NDNCRESO2       =  7C +  6H +  2O + 3N ; {@NDNCRESO2}         {MCM}
NDNCRESOOH      =  7C +  7H + 12O + 3N ; {@NDNCRESOOH}        {MCM}
NPTLQO2         =  7C +  6H +  7O +  N ; {@NPTLQO2}           {MCM}
NPTLQOOH        =  7C +  7H +  7O +  N ; {@NPTLQOOH}          {MCM}
PBZN            =  7C +  5H +  5O +  N ; {@PBZN}              {MCM: benzoyl nitro peroxide}
TLBIPERNO3      =  7C +  9H +  6O +  N ; {@TLBIPERNO3}        {MCM}
TLEMUCNO3       =  7C +  9H +  7O +  N ; {@TLEMUCNO3}         {MCM}
TLEMUCPAN       =  7C +  7H +  7O +  N ; {@TLEMUCPAN}         {MCM}
TOL1OHNO2       =  7C +  7H +  3O +  N ; {@TOL1OHNO2}         {MCM: 2-methyl-6-nitrophenol}
{8C (CHO)}
C721CHO         =  8C + 12H +  3O      ; {@C721CHO}           {MCM}
C721CO3         =  8C + 11H +  5O      ; {@C721CO3}           {MCM}
C721CO3H        =  8C + 12H +  5O      ; {@C721CO3H}          {MCM}
C810O2          =  8C + 13H +  4O      ; {@C810O2}            {MCM}
C810OOH         =  8C + 14H +  4O      ; {@C810OOH}           {MCM}
C811O2          =  8C + 13H +  4O      ; {@C811O2}            {MCM}
C812O2          =  8C + 13H +  5O      ; {@C812O2}            {MCM}
C812OOH         =  8C + 14H +  5O      ; {@C812OOH}           {MCM}
C813O2          =  8C + 13H +  6O      ; {@C813O2}            {MCM}
C813OOH         =  8C + 14H +  5O      ; {@C813OOH}           {MCM}
C85O2           =  8C + 13H +  3O      ; {@C85O2}             {MCM}
C85OOH          =  8C + 14H +  3O      ; {@C85OOH}            {MCM}
C86O2           =  8C + 13H +  4O      ; {@C86O2}             {MCM}
C86OOH          =  8C + 14H +  4O      ; {@C86OOH}            {MCM}
C89O2           =  8C + 13H +  3O      ; {@C89O2}             {MCM}
C89OOH          =  8C + 14H +  3O      ; {@C89OOH}            {MCM}
C8BC            =  8C + 14H            ; {@C8BC}              {MCM}
C8BCCO          =  8C + 12H +  O       ; {@C8BCCO}            {MCM}
C8BCO2          =  8C + 11H +  2O      ; {@C8BCO2}            {MCM}
C8BCOOH         =  8C + 12H +  2O      ; {@C8BCOOH}           {MCM}
NORPINIC        =  8C + 12H +  4O      ; {@NORPINIC}          {MCM}
{C8 (CHO) aromatics}
EBENZ           =  8C + 10H            ; {@EBENZ}             {MCM: ethylbenzene}
STYRENE         =  8C +  8H            ; {@STYRENE}           {MCM}
STYRENO2        =  8C +  9H +  3O      ; {@STYRENO2}          {MCM}
STYRENOOH       =  8C + 10H +  3O      ; {@STYRENOOH}         {MCM}
{8C (CHON)}
C721PAN         =  8C + 11H +  7O +  N ; {@C721PAN}           {MCM}
C810NO3         =  8C + 14H +  5O +  N ; {@C810NO3}           {MCM}
C89NO3          =  8C + 13H +  4O +  N ; {@C89NO3}            {MCM}
C8BCNO3         =  8C + 11H +  3O +  N ; {@C8BCNO3}           {MCM}
{C8 (CHON) aromatics}
NSTYRENO2       =  8C +  8H +  5O +  N ; {@NSTYRENO2}         {MCM}
NSTYRENOOH      =  8C +  9H +  5O +  N ; {@NSTYRENOOH}        {MCM}
{C8 aromatics (lumped)}
LXYL            =  8C + 10H            ; {@LXYL}              {xylenes}
{9C (CHO)}
C811CO3         =  9C + 13H +  5O      ; {@C811CO3}           {MCM}
C811CO3H        =  9C + 14H +  5O      ; {@C811CO3H}          {MCM}
C85CO3          =  9C + 11H +  4O      ; {@C85CO3}            {MCM}
C85CO3H         =  9C + 12H +  4O      ; {@C85CO3H}           {MCM}
C89CO2H         =  9C + 14H +  3O      ; {@C89CO2H}           {MCM}
C89CO3          =  9C + 13H +  4O      ; {@C89CO3}            {MCM}
C89CO3H         =  9C + 14H +  4O      ; {@C89CO3H}           {MCM}
C96O2           =  9C + 15H +  3O      ; {@C96O2}             {MCM}
C96OOH          =  9C + 16H +  3O      ; {@C96OOH}            {MCM}
C97O2           =  9C + 15H +  4O      ; {@C97O2}             {MCM}
C97OOH          =  9C + 16H +  4O      ; {@C97OOH}            {MCM}
C98O2           =  9C + 15H +  5O      ; {@C98O2}             {MCM}
C98OOH          =  9C + 16H +  5O      ; {@C98OOH}            {MCM}
NOPINDCO        =  9C + 12H +  2O      ; {@NOPINDCO}          {MCM}
NOPINDO2        =  9C + 13H +  3O      ; {@NOPINDO2}          {MCM}
NOPINDOOH       =  9C + 14H +  3O      ; {@NOPINDOOH}         {MCM}
NOPINONE        =  9C + 14H +   O      ; {@NOPINONE}          {MCM}
NOPINOO         =  9C + 14H +  2O      ; {@NOPINOO}           {MCM}
NORPINAL        =  9C + 14H +  2O      ; {@NORPINAL}          {MCM: norpinaldehyde}
NORPINENOL      =  9C + 14H +  2O      ; {@NORPINENOL}        {}
PINIC           =  9C + 14H +  4O      ; {@PINIC}             {MCM: pinic acid}
{9C (CHON)}
C811PAN         =  9C + 13H +  7O +  N ; {@C811PAN}           {MCM}
C89PAN          =  9C + 13H +  5O +  N ; {@C89PAN}            {MCM}
C96NO3          =  9C + 15H +  4O +  N ; {@C96NO3}            {MCM}
C9PAN2          =  9C + 13H +  6O +  N ; {@C9PAN2}            {MCM}
{C9 aromatics (lumped)}
LTMB            =  9C + 12H            ; {@LTMB}              {trimethylbenzenes}
{10C (CHO)}
APINAOO         = 10C + 16H +  3O      ; {@APINAOO}           {stabilized APINOOA}
APINBOO         = 10C + 16H +  3O      ; {@APINBOO}           {MCM}
APINENE         = 10C + 16H            ; {@APINENE}           {MCM: alpha pinene}
BPINAO2         = 10C + 17H +  3O      ; {@BPINAO2}           {MCM}
BPINAOOH        = 10C + 18H +  3O      ; {@BPINAOOH}          {MCM}
BPINENE         = 10C + 16H            ; {@BPINENE}           {MCM: beta pinene}
C106O2          = 10C + 15H +  5O      ; {@C106O2}            {MCM}
C106OOH         = 10C + 16H +  5O      ; {@C106OOH}           {MCM}
C109CO          = 10C + 10H +  3O      ; {@C109CO}            {MCM}
C109O2          = 10C + 15H +  4O      ; {@C109O2}            {MCM}
C109OOH         = 10C + 16H +  4O      ; {@C109OOH}           {MCM}
C96CO3          = 10C + 15H +  4O      ; {@C96CO3}            {MCM}
CAMPHENE        = 10C + 16H            ; {@CAMPHENE}          {}
CARENE          = 10C + 16H            ; {@CARENE}            {3-carene}
MENTHEN6ONE     = 10C + 16H +  3O      ; {@MENTHEN6ONE}       {8-OOH-menthen-6-one, Taraborrelli, pers. comm.}
OH2MENTHEN6ONE  = 10C + 17H +  4O      ; {@2OHMENTHEN6ONE}    {2-OH-8-OOH-menthen-6-one, Taraborrelli, pers. comm.}
OHMENTHEN6ONEO2 = 10C + 17H +  5O      ; {@OHMENTHEN6ONEO2}   {2-OH-8-OOH_menthen-6-peroxy radical, Taraborrelli, pers. comm.}
PERPINONIC      = 10C + 16H +  4O      ; {@PERPINONIC}        {MCM}
PINAL           = 10C + 16H +  2O      ; {@PINAL}             {MCM: pinonaldehyde}
PINALO2         = 10C + 13H +  4O      ; {@PINALO2}           {MCM}
PINALOOH        = 10C + 14H +  4O      ; {@PINALOOH}          {MCM}
PINENOL         = 10C + 16H +  2O      ; {@PINEOL}            {}
PINONIC         = 10C + 16H +  3O      ; {@PINONIC}           {MCM: pinonic acid}
RO6R1O2         = 10C + 17H +  4O      ; {@RO6R1O2}           {cyclo-oxy peroxy radical from BPINENE, ref3019}
RO6R3O2         = 10C + 17H +  5O      ; {@RO6R3O2}           {cyclo-oxy peroxy radical from BPINENE, ref3019}
ROO6R1O2        = 10C + 17H +  5O      ; {@ROO6R1O2}          {cyclo-peroxy peroxy radical from BPINENE based on ROO6R1 from ref3019}
SABINENE        = 10C + 16H            ; {@SABINENE}          {}
{10C (CHON)}
BPINANO3        = 10C + 17H +  4O +  N ; {@BPINANO3}          {MCM}
C106NO3         = 10C + 15H +  6O +  N ; {@C106NO3}           {MCM}
C10PAN2         = 10C + 15H +  6O +  N ; {@C10PAN2}           {MCM}
PINALNO3        = 10C + 13H +  5O +  N ; {@PINALNO3}          {MCM}
RO6R1NO3        = 10C + 17H +  5O +  N ; {@RO6R1NO3}          {nitrate from cyclo-oxy peroxy radical from BPINENE, ref3019}
ROO6R1NO3       = 10C + 17H +  6O +  N ; {@ROO6R1NO3}         {nitrate from cyclo-peroxy peroxy radical from BPINENE, ref3019}
{10C (lumped)}
LAPINABNO3      = 10C + 17H +  4O +  N ; {@LAPINABNO3}        {APINANO3 + APINBNO3 lumped (ratio 1:2)}
LAPINABO2       = 10C + 17H +  3O      ; {@LAPINABO2}         {APINAO2 + APINBO2 lumped (ratio 1:2)}
LAPINABOOH      = 10C + 18H +  3O      ; {@LAPINABOOH}        {APINAOOH + APINBOOH lumped (ratio 1:2)}
LNAPINABO2      = 10C + 16H +  5O +  N ; {@LNAPINABO2}        {.65 NAPINAO2 + .35 NAPINBO2}
LNAPINABOOH     = 10C + 17H +  5O +  N ; {@LNAPINABOOH}       {.65 NAPINAOOH + .35 NAPINBOOH}
LNBPINABO2      = 10C + 16H +  5O +  N ; {@LNBPINABO2}        {.8 NBPINAO2 + .2 NBPINBO2}
LNBPINABOOH     = 10C + 17H +  5O +  N ; {@LNBPINABOOH}       {.8 NBPINAO2 + .2 NBPINBO2}
{C10 aromatics (lumped)}
LHAROM          = 11C + 14H            ; {@LHAROM}            {higher aromatics: model compound DIET35TOL(from MCM)}
{------------------------------------- F ------------------------------------}

LFLUORINE       =            F         ; {@LFLUORINE}         {lumped F species}
CHF3            =  C +  H + 3F         ; {@CHF_3}             {trifluoromethane, fluoroform = HFC-23}
CHF2CF3         = 2C +  H + 5F         ; {@CHF_2CF_3}         {pentafluoroethane = HFC-125}
CH3CF3          = 2C + 3H + 3F         ; {@CH_3CF_3}          {1,1,1-trifluoroethane = HFC-143a}
CH2F2           =  C + 2H + 2F         ; {@CH_2F_2}           {difluoromethane = HFC-32}
CH3CHF2         = 2C + 4H + 2F         ; {@CH_3CHF_2}         {1,1-difluoroethane = HFC-152a}
{------------------------------------- Cl -----------------------------------}

CCl4            =  C                    + 4Cl ; {@CCl_4}             {tetrachloro methane}
CF2Cl2          =  C               + 2F + 2Cl ; {@CF_2Cl_2}          {dichlorodifluoromethane = F12}
CF2ClCF2Cl      = 2C               + 4F + 2Cl ; {@CF_2ClCF_2Cl}      {1,1,2,2-tetrafluoro-1,2-dichloroethane = CFC-114}
CF2ClCFCl2      = 2C               + 3F + 3Cl ; {@CF_2ClCFCl_2}      {1,1,2-trifluoro-1,2,2-trichloroethane = CFC-113}
CF3CF2Cl        = 2C               + 5F +  Cl ; {@CF_3CF_2Cl}        {pentafluorochloroethane = CFC-115}
CFCl3           =  C               +  F + 3Cl ; {@CFCl_3}            {trichlorofluoromethane = F11}
CH2Cl2          =  C + 2H               + 2Cl ; {@CH_2Cl_2}          {dichloromethane}
CH2FCF3         = 2C + 2H          + 4F       ; {@CH_2FCF_3}         {1,1,1,2-tetrafluoroethane = HFC-134a}
CH3CCl3         = 2C + 3H               + 3Cl ; {@CH_3CCl_3}         {1,1,1-trichloroethane = methyl chloroform = MCF}
CH3CFCl2        = 2C + 3H          +  F + 2Cl ; {@CH_3CFCl_2}        {1,1,1-fluorodichloroethane = HCFC-141b}
CH3Cl           =  C + 3H               +  Cl ; {@CH_3Cl}            {chloromethane}
CHCl3           =  C +  H               + 3Cl ; {@CHCl_3}            {trichloromethane = chloroform}
CHF2Cl          =  C +  H          + 2F +  Cl ; {@CHF_2Cl}           {difluorochloromethane = HCFC-22}
Cl              =                          Cl ; {@Cl}                {chlorine atom}
Cl2             =                         2Cl ; {@Cl_2}              {chlorine}
Cl2O2           =           2O          + 2Cl ; {@Cl_2O_2}           {dichlorine dioxide}
ClNO2           =           2O + N      +  Cl ; {@ClNO_2}            {nitryl chloride}
ClNO3           =           3O + N      +  Cl ; {@ClNO_3}            {chlorine nitrate}
ClO             =            O          +  Cl ; {@ClO}               {chlorine oxide}
HCl             =       H               +  Cl ; {@HCl}               {hydrochloric acid}
HOCl            =       H +  O          +  Cl ; {@HOCl}              {hypochlorous acid}
OClO            =           2O          +  Cl ; {@OClO}              {chlorine dioxide}
LCHLORINE       =                          Cl ; {@LCHLORINE}         {lumped Cl species}

{------------------------------------- Br -----------------------------------}

Br              =                               Br ; {@Br}                {bromine atom}
Br2             =                              2Br ; {@Br_2}              {bromine}
BrCl            =                         Cl +  Br ; {@BrCl}              {bromine chloride}
BrNO2           =          2O + N            +  Br ; {@BrNO_2}            {nitryl bromide}
BrNO3           =          3O + N            +  Br ; {@BrNO_3}            {bromine nitrate}
BrO             =           O                +  Br ; {@BrO}               {bromine oxide}
CF2ClBr         = C               + 2F +  Cl +  Br ; {@CF_2ClBr}          {Halon 1211}
CF3Br           = C               + 3F       +  Br ; {@CF_3Br}            {Halon 1301}
CH2Br2          = C + 2H                     + 2Br ; {@CH_2Br_2}          {}
CH2ClBr         = C + 2H               +  Cl +  Br ; {@CH_2ClBr}          {}
CH3Br           = C + 3H                     +  Br ; {@CH_3Br}            {bromomethane}
CHBr3           = C +  H                     + 3Br ; {@CHBr_3}            {}
CHCl2Br         = C +  H               + 2Cl +  Br ; {@CHCl_2Br}          {}
CHClBr2         = C +  H               +  Cl + 2Br ; {@CHClBr_2}          {}
HBr             =      H                     +  Br ; {@HBr}               {hydrobromic acid}
HOBr            =      H +  O                +  Br ; {@HOBr}              {hypobromous acid}
LBROMINE        =                               Br ; {@LBROMINE}          {lumped Br species}

{------------------------------------- I ------------------------------------}

C3H7I           = 3C + 7H                    +  I ; {@CH_3CHICH_3}       {2-iodopropane}
CH2ClI          =  C + 2H          + Cl      +  I ; {@CH_2ClI}           {chloroiodomethane}
CH2I2           =  C + 2H                    + 2I ; {@CH_2I_2}           {diiodomethane}
CH3I            =  C + 3H                    +  I ; {@CH_3I}             {iodomethane}
HI              =       H                    +  I ; {@HI}                {hydrogen iodide}
HIO3            =       H + 3O               +  I ; {@HIO_3}             {}
HOI             =       H +  O               +  I ; {@HOI}               {hypoiodous acid}
I               =                               I ; {@I}                 {iodine atomic ground state}
I2              =                              2I ; {@I_2}               {molecular iodine}
I2O2            =           2O               + 2I ; {@I_2O_2}            {}
IBr             =                         Br +  I ; {@IBr}               {iodine bromide}
ICl             =                    Cl +       I ; {@ICl}               {iodine chloride}
INO2            =           2O + N           +  I ; {@INO_2}             {iodine nitrite}
INO3            =           3O + N           +  I ; {@INO_3}             {iodine nitrate}
IO              =            O               +  I ; {@IO}                {iodine monoxide radical}
IPART           =                              2I ; {@I(part)}           {iodine particles}
OIO             =           2O               +  I ; {@OIO}               {}

{------------------------------------- S ------------------------------------}

CH3SO2          =  C + 3H + 2O    + S ; {@CH_3SO_2}          {}
CH3SO3          =  C + 3H + 3O    + S ; {@CH_3SO_3}          {}
CH3SO3H         =  C + 4H + 3O    + S ; {@CH_3SO_3H}         {MSA: methane sulfonic acid}
DMS             = 2C + 6H         + S ; {@DMS}               {dimethyl sulfide}
DMSO            = 2C + 6H +  O    + S ; {@DMSO}              {dimethyl sulfoxide: CH3SOCH3}
H2SO4           =      2H + 4O    + S ; {@H_2SO_4}           {sulfuric acid}
OCS             =  C      +  O    + S ; {@OCS}               {}
S               =                   S ; {@S}                 {sulfur atomic ground state}
SF6             =              6F + S ; {@SF_6}              {sulfur hexaflouride}
SH              =       H         + S ; {@SH}                {}
SO              =            O    + S ; {@SO}                {sulfur monoxide}
SO2             =           2O    + S ; {@SO_2}              {sulfur dioxide}
SO3             =           3O    + S ; {@SO_3}              {sulfur trioxide}
LSULFUR         =                   S ; {@LSULFUR}           {lumped S species}

{--------------------------------- Hg ---------------------------------------}

Hg              = Hg                   ; {@Hg}                {}
HgO             = Hg + O               ; {@HgO}               {}
HgCl            = Hg + Cl              ; {@HgCl}              {}
HgCl2           = Hg + 2Cl             ; {@HgCl_2}            {}
HgBr            = Hg + Br              ; {@HgBr}              {}
HgBr2           = Hg + 2Br             ; {@HgBr_2}            {}
ClHgBr          = Hg + Cl + Br         ; {@ClHgBr}            {}
BrHgOBr         = Hg + O + 2Br         ; {@BrHgOBr}           {}
ClHgOBr         = Hg + O + Cl + Br     ; {@ClHgOBr}           {}

{--- mz_pj_20070209+}
{------------------------- Pseudo Aerosol -----------------------------------}
NO3m_cs         = N + 3O               ; {@NO_3^-(cs)}        {}
Hp_cs           = H                    ; {@H^+(cs)}           {}
RGM_cs          = Hg                   ; {@Hg(cs)}            {from reactive gaseous Hg}
{--- mz_pj_20070209-}

{------------------------------- Dummies ------------------------------------}

Dummy           = IGNORE               ; {@Dummy}
PRODUCTS        = IGNORE               ; {@PRODUCTS}
M               = IGNORE               ; {@M}

{ mz_pj_20070621+}
{------------------------- O3 Budget Tracers (via eval2.3.rpl) --------------}
O3s             = 3O                   ; {@O_3(s)}            {strat. ozone}
LO3s            = IGNORE               ; {@LO_3(s)}           {lost strat. ozone}
{ mz_pj_20070621-}

{ mz_rs_20100227+}
{only for MIM1, not used in MIM2:}
ISO2            = 5C + 9H + 3O         ; {@ISO2}              {isoprene (hydroxy) peroxy radicals}
ISON            = 5C +           N     ; {@ISON}              {organic nitrates from ISO2 and C5H8+NO3}
ISOOH           = 5C + 10H + 3O        ; {@ISOOH}             {isoprene (hydro) peroxides}
LHOC3H6O2       = 3C + 7H + 3O         ; {@CH_3CH(O_2)CH_2OH} {hydroxyperoxyradical from propene+OH}
LHOC3H6OOH      = 3C + 8H + 3O         ; {@CH_3CH(OOH)CH_2OH} {C3H6OHOOH = hydroxyhydroperoxides from C3H6}
MVKO2           = 4C + 7H + 4O         ; {@MVKO2}             {MVK/MACR peroxy radicals}
MVKOOH          = 4C + 8H + 4O         ; {@MVKOOH}            {MVK hydroperoxides}
NACA            = 2C + 3H + 4O + N     ; {@NACA}              {nitro-oxy acetaldehyde}
{ mz_rs_20100227-}

{ mz_ab_20100908+}
{---------------------------------- ions ------------------------------------}
Op              =  O           + Pls   ; {@O^+}               {O+}
O2p             =  2O          + Pls   ; {@O_2^+}             {O2+}
Np              =  N           + Pls   ; {@N^+}               {N+}
N2p             =  2N          + Pls   ; {@N_2^+}             {N2+}
NOp             =  O + N       + Pls   ; {@NO^+}              {NO+}
em              =                Min   ; {@e^-}               {electron}
kJmol           =  IGNORE              ; {@kJ/mol}            {released energy}
{ mz_ab_20100908-}

{ op_pj_20130723+}
{------------------------------ additional diagnostic tracers ----------------}
CFCl3_c         = C + F + 3Cl          ; {@(CFCl_3)_c}        {trichlorofluoromethane = F11}
CF2Cl2_c        = C + 2F + 2Cl         ; {@(CF_2Cl_2)_c}      {dichlorodifluoromethane = F12}
N2O_c           = O + 2N               ; {@(N_2O)_c}          {nitrous oxide}
CH3CCl3_c       = 2C + 3H + 3Cl        ; {@(CH_3CCl_3)_c}     {1,1,1-trichloroethane = methyl chloroform = MCF}
CF2ClBr_c       = Br + 2F + Cl + C     ; {@(CF_2ClBr)_c}      {Halon 1211}
CF3Br_c         = Br + 3F + C          ; {@(CF_3Br)_c}        {Halon 1301}
{ op_pj_20130723-}

{ mz_at_20131015+ needed for ORACLE.rpl}
{-----------------------Organic Condesable Gases and VOCs--------------------}
LTERP           =  IGNORE              ; {@LTERP}             {terpenes}
LALK4           =  IGNORE              ; {@LALK4}             {alkanes}
LALK5           =  IGNORE              ; {@LALK5}             {alkanes}
LARO1           =  IGNORE              ; {@LARO1}             {aromatic VOC}
LARO2           =  IGNORE              ; {@LARO2}             {aromatic VOC}
LOLE1           =  IGNORE              ; {@LOLE1}             {olefins}
LOLE2           =  IGNORE              ; {@LOLE2}             {olefins}
LfPOG02         =  IGNORE              ; {@LfPOG02}           {FF  condensable gas 2}
LfPOG03         =  IGNORE              ; {@LfPOG03}           {FF  condensable gas 3}
LfPOG04         =  IGNORE              ; {@LfPOG04}           {FF  condensable gas 4}
LfPOG05         =  IGNORE              ; {@LfPOG05}           {FF  condensable gas 5}
LbbPOG02        =  IGNORE              ; {@LbbPOG02}          {BB  condensable gas 2}
LbbPOG03        =  IGNORE              ; {@LbbPOG03}          {BB  condensable gas 3}
LbbPOG04        =  IGNORE              ; {@LbbPOG04}          {BB  condensable gas 4}
LfSOGsv01       =  IGNORE              ; {@LfSOGsv01}         {sFF condensable gas 1}
LfSOGsv02       =  IGNORE              ; {@LfSOGsv02}         {sFF condensable gas 2}
LbbSOGsv01      =  IGNORE              ; {@LbbSOGsv01}        {sBB condensable gas 1}
LbbSOGsv02      =  IGNORE              ; {@LbbSOGsv02}        {sBB condensable gas 2}
LfSOGiv01       =  IGNORE              ; {@LfSOGiv01}         {iFF condensable gas 1}
LfSOGiv02       =  IGNORE              ; {@LfSOGiv02}         {iFF condensable gas 2}
LfSOGiv03       =  IGNORE              ; {@LfSOGiv03}         {iFF condensable gas 3}
LfSOGiv04       =  IGNORE              ; {@LfSOGiv04}         {iFF condensable gas 4}
LbbSOGiv01      =  IGNORE              ; {@LbbSOGiv01}        {iBB condensable gas 1}
LbbSOGiv02      =  IGNORE              ; {@LbbSOGiv02}        {iBB condensable gas 2}
LbbSOGiv03      =  IGNORE              ; {@LbbSOGiv03}        {iBB condensable gas 3}
LbSOGv01        =  IGNORE              ; {@LbSOGv01}          {Bio condensable gas 1}
LbSOGv02        =  IGNORE              ; {@LbSOGv02}          {Bio condensable gas 2}
LbSOGv03        =  IGNORE              ; {@LbSOGv03}          {Bio condensable gas 3}
LbSOGv04        =  IGNORE              ; {@LbSOGv04}          {Bio condensable gas 4}
LbOSOGv01       =  IGNORE              ; {@LbOSOGv01}         {Bio condensable gas 1}
LbOSOGv02       =  IGNORE              ; {@LbOSOGv02}         {Bio condensable gas 2}
LbOSOGv03       =  IGNORE              ; {@LbOSOGv03}         {Bio condensable gas 3}
LaSOGv01        =  IGNORE              ; {@LaSOGv01}          {Ant condensable gas 1}
LaSOGv02        =  IGNORE              ; {@LaSOGv02}          {Ant condensable gas 2}
LaSOGv03        =  IGNORE              ; {@LaSOGv03}          {Ant condensable gas 3}
LaSOGv04        =  IGNORE              ; {@LaSOGv04}          {Ant condensable gas 4}
LaOSOGv01       =  IGNORE              ; {@LaOSOGv01}         {Ant condensable gas 1}
LaOSOGv02       =  IGNORE              ; {@LaOSOGv02}         {Ant condensable gas 2}
LaOSOGv03       =  IGNORE              ; {@LaOSOGv03}         {Ant condensable gas 3}
{ mz_at_20131015- needed for ORACLE.rpl}

{ mz_rs_20170601+ jam}
ACBZO2          =                     5H + 7C + 3O ; {@C_7H_5O_3}             {acyl peroxy radical from benzaldehyde}
ALKNO3          =               11H + 5C + 3O +  N ; {@C_5H_<11>NO_3}         {nitrate from BIGALKANE}
ALKO2           =                    11H + 5C + 2O ; {@C_5H_<11>O_2}          {peroxy radical from large alkanes}
ALKOH           =                    12H + 5C +  O ; {@C_5H_<12>O}            {alcohol from BIGALKANE}
ALKOOH          =                    12H + 5C + 2O ; {@C_5H_<12>O_2}          {peroxide from large alkanes}
BCARY           =                        24H + 15C ; {@C_<15>H_<24>}          {(1R,4E,9S)-4,11,11-trimethyl-8-methylidenebicyclo[7.2.0]undec-4-ene}
BENZO2          =                     7H + 6C + 5O ; {@C_6H_7O_5}             {peroxy radical from benzene}
BENZOOH         =                     8H + 6C + 5O ; {@C_6H_8O_5}             {peroxide from BENZO2}
BEPOMUC         =                     6H + 6C + 3O ; {@C_6H_6O_3}             {benzene eopoxy diol}
BIGALD1         =                     4H + 4C + 2O ; {@C_4H_4O_2}             {but-2-enedial}
BIGALD2         =                     6H + 5C + 2O ; {@C_5H_6O_2}             {4-oxopent-2-enal}
BIGALD3         =                     6H + 5C + 2O ; {@C_5H_6O_2}             {2-methylbut-2-enedial}
BIGALD4         =                     8H + 6C + 2O ; {@C_6H_8O_2}             {aldehyde from xylene oxidation}
BIGALKANE       =                         12H + 5C ; {@C_5H_<12>}             {large alkanes}
BIGENE          =                          8H + 4C ; {@C_4H_8}                {large alkenes}
BrONO           = IGNORE                           ; {@BrONO}
BZALD           =                     6H + 7C +  O ; {@C_7H_6O}               {benzaldehyde}
BZOO            =                     7H + 7C + 2O ; {@C_7H_7O_2}             {peroxy radical from toluene}
BZOOH           =                     8H + 7C + 2O ; {@C_7H_8O_2}             {peroxide from BZOO}
C3H7O2          =                     7H + 3C + 2O ; {@C_3H_7O_2}             {lumped peroxy radical from propane}
C3H7OOH         =                     8H + 3C + 2O ; {@C_3H_8O_2}             {lumped propyl hydro peroxide}
CFC113          =                    2C + 3F + 3Cl ; {@C_2F_3Cl_3}            {1,1,2-trichloro-1,2,2-trifluoroethane}
CFC114          =                    2C + 4F + 2Cl ; {@C_2F_4Cl_2}            {1,2-dichloro-1,1,2,2-tetrafluoro-ethane}
CFC115          =                    2C + 5F +  Cl ; {@C_2F_5Cl}              {1-chloro-1,1,2,2,2-pentafluoro-ethane}
COF2            =                      C +  O + 2F ; {@CF_2O}                 {carbonyl difluoride}
COFCL           =                C +  F +  O +  Cl ; {@CFClO}                 {carbonyl chloride fluoride}
DICARBO2        =                     5H + 5C + 4O ; {@C_5H_5O_4}             {dicarbonyl from photolysis of BIGALD2}
ELVOC           = IGNORE                           ; {@ELVOC}
ENEO2           =                     9H + 4C + 3O ; {@C_4H_9O_3}             {peroxy radical from BIGENE/OLTP}
EOOH            =                     6H + 2C + 3O ; {@C_2H_6O_3}             {2-hydroperoxyethanol}
F               =                                F ; {@F}                     {fluoride}
H1202           =                     C + 2Br + 2F ; {@CF_2Br_2}              {dibromo(difluoro)methane}
H2402           =                    2C + 2Br + 4F ; {@C_2F_4Br_2}            {1,2-dibromo-1,1,2,2-tetrafluoroethane}
HCFC141B        =               3H + 2C +  F + 2Cl ; {@C_2H_3FCl_2}           {1,1-dichloro-1-fluoroethane}
HCFC142B        =               3H + 2C + 2F +  Cl ; {@C_2H_3F_2Cl}           {1-chloro-1,1-difluoroethane}
HCFC22          =                H +  C + 2F +  Cl ; {@CHF_2Cl}               {chloro(difluoro)methane}
HF              =                           H +  F ; {@HF}                    {fluorane}
HOCH2OO         =                     3H +  C + 3O ; {@CH_3O_3}               {(hydroxymethyl)dioxidanyl}
HPALD           = IGNORE                           ; {@HPALD}
IEC1O2          =                     9H + 5C + 5O ; {@C_5H_9O_5}             {peroxy radical from LIEPOX+OH}
LIECHO          =                     8H + 5C + 3O ; {@C_5H_8O_3}             {aldehyde from LIEPOX}
LIECO3          =                     7H + 5C + 5O ; {@C_5H_7O_5}             {peroxy radical from LIECHO}
LIECO3H         =                     8H + 5C + 5O ; {@C_5H_8O_5}             {peroxide from LIECO3}
LIMON           =                        16H + 10C ; {@C_<10>H_<16>}          {1-methyl-4-prop-1-en-2-ylcyclohexene}
LISOPNO3NO3     = IGNORE                           ; {@LISOPNO3NO3}
LISOPNO3O2      = IGNORE                           ; {@LISOPNO3O2}
LISOPNO3OOH     = IGNORE                           ; {@LISOPNO3OOH}
LISOPOOHO2      = IGNORE                           ; {@LISOPOOHO2}
LISOPOOHOOH     = IGNORE                           ; {@LISOPOOHOOH}
MALO2           =                     3H + 4C + 4O ; {@C_4H_3O_4}             {peroxy radical from photolysis of BIGALD1}
MBONO3O2        =               10H + 5C + 6O +  N ; {@C_5H_<10>NO_6}         {peroxy nitrate radical from MBO+NO3}
MBOO2           =                    11H + 5C + 4O ; {@C_5H_<11>O_4}          {peroxy radical from MBO}
MBOOOH          =                    12H + 5C + 4O ; {@C_5H_<12>O_4}          {peroxide from MBO}
MDIALO2         =                     5H + 5C + 4O ; {@C_5H_5O_4}             {peroxy radical from photolysis of BIGALD3}
MEKNO3          = IGNORE                           ; {@MEKNO3}
MVKN            = IGNORE                           ; {@MVKN}
MYRC            =                        16H + 10C ; {@C_<10>H_<16>}          {2-methyl-6-methylideneocta-1,7-diene}
NTERPNO3        = IGNORE                           ; {@NTERPNO3}
NTERPO2         =              16H + 10C + 5O +  N ; {@C_<10>H_<16>NO_5}      {nitro peroxy radical from terpenes}
PACALD          = IGNORE                           ; {@PACALD}
PBZNIT          =                5H + 7C + 5O +  N ; {@C_7H_5NO_5}            {nitrate from benzaldehyde}
TEPOMUC         =                     8H + 7C + 3O ; {@C_7H_8O_3}             {epoxide from toluene}
TERP2O2         =                   15H + 10C + 4O ; {@C_<10>H_<15>O_4}       {peroxy radical from TERPROD1}
TERP2OOH        =                   16H + 10C + 4O ; {@C_<10>H_<16>O_4}       {peroxide from TERP2O2}
TERPNO3         =              17H + 10C + 4O +  N ; {@C_<10>H_<17>NO_4}      {nitrate from terpenes}
TERPO2          =                   17H + 10C + 3O ; {@C_<10>H_<17>O_3}       {peroxy radical from terpenes}
TERPOOH         =                   18H + 10C + 3O ; {@C_<10>H_<18>O_3}       {peroxide from terpenes}
TERPROD1        =                   16H + 10C + 2O ; {@C_<10>H_<16>O_2}       {terpene oxidation product C10}
TERPROD2        =                    10H + 7C + 2O ; {@C_7H_<10>O_2}          {terpene oxidation product C9}
TOLO2           =                     9H + 7C + 5O ; {@C_7H_9O_5}             {peroxy radical from toluene}
TOLOOH          =                    10H + 7C + 5O ; {@C_7H_<10>O_5}          {peroxide from toluene}
XYLENO2         =                    11H + 8C + 5O ; {@C_8H_<11>O_5}          {peroxy radical from xylene}
XYLENOOH        =                    12H + 8C + 5O ; {@C_8H_<12>O_5}          {peroxide from XYLENO2}
XYLOL           =                    10H + 8C +  O ; {@C_8H_<10>O}            {2,3-dimethylphenol}
XYLOLO2         =                    11H + 8C + 6O ; {@C_8H_<11>O_6}          {peroxy radical from xylol}
XYLOLOOH        =                    12H + 8C + 6O ; {@C_8H_<12>O_6}          {peroxide from xylol}
{ mz_rs_20170601-}

{ mz_rs_20171213+ MOZART}
O2_1D           = 2O                 ; {@O_2}               {excited molecular oxygen (singlett D state)}
O2_1S           = 2O                 ; {@O_2}               {excited molecular oxygen (singlett S state)}
ONIT            =  3C +  5H + 4O + N ; {@C_3H_5NO_4}        {organic nitrate}
C4H8            =  4C +  8H          ; {@C4H8}              {large alkenes}
C4H9O3          =  4C +  9H + 3O     ; {@C_4H_9O_3}         {peroxy radical from C4H8}
C5H12           =  5C + 12H          ; {@C5H12}             {large alkanes}
C5H11O2         =  5C + 11H + 2O     ; {@C5H11O2}           {peroxy radical from large alkanes}
C5H6O2          =  5C +  6H + 2O     ; {@C5H6O2}            {aldehyde from toluene oxidation}
HYDRALD         =  5C +  8H + 2O     ; {@C_5H_8O_2}         {lumped unsaturated hydroxycarbonyl}
ISOPO2          =  5C +  9H + 3O     ; {@C_5H_9O_3}         {lumped peroxy radical from isoprene}
C5H9O3          =  5C +  9H + 4O     ; {@C_5H_9O_4}         {peroxy radical from OH+HYDRALD}
ISOPOOH         =  5C + 10H + 3O     ; {@C_5H_10O_3}        {peroxide from isoprene}
C5H12O2         =  5C + 12H + 2O     ; {@C5H12O2}           {peroxide from large alkanes}
ONITR           =  5C +  9H + 4O + N ; {@C_5H_9NO_4}        {alkyl nitrate from ISOPO2+NO3}
C5H10O4         =  5C + 10H + 4O     ; {@C_5H_10O_4}        {peroxide from C5H9O3}
ROO6R5P         =  7C + 10H + 6O     ; {@ROO6R5P}           {from ref3019}
NH4             =        4H      + N ; {@NH_4}              {aq. ammonium ion}
SO4             = S + 4O             ; {@SO_4}              {aq. sulfate}
{ mz_rs_20171213-}

{ mz_rs_20171213+ CB05BASCOE}
HCO             =  C +   H +  O      ; {@HCO}               {CHO formyl radical}
ISPD            =  4C +  6H +  O     ; {@ISPD}              {lumped MACR MVK}
ClOO            = Cl + 2O            ; {@CLOO}              {asymmetrical chlorine dioxide radical}
Rn              = Rn                 ; {@Rn}                {radon}
Pb              = Pb                 ; {@Pb}                {lead}
XO2             = IGNORE             ; {@XO2}               {NO_to_NO2_operator}
XO2N            = IGNORE             ; {@XO2N}              {NO_to_alkyl_nitrate_operator}
ROOH            = IGNORE             ; {@ROOH}              {peroxides}
OLE             = IGNORE             ; {@OLE}               {olefins}
ROR             = IGNORE             ; {@ROR}               {organic_ethers}
ORGNTR          = IGNORE             ; {@ORGNTR}            {organic nitrates called ONIT in mocage}
ACO2            = IGNORE             ; {@ACO2}              {acetone oxidation product}
PAR             = IGNORE             ; {@PAR}               {parafins}
RXPAR           = IGNORE             ; {@RXPAR}             {olefins}
{ mz_rs_20171213-}
