{----------------------------------------------------------------------------}

{----------------------------------------------------------------------------}
{--------------------------- gas phase cont. --------------------------------}
{----------------------------------------------------------------------------}


{-------------------------- SOA(g) species ----------------------------------}

BSOV          = IGNORE                ; {@BSOV}              {SVOC,  secondary oxidized biogenic}
BLOV          = IGNORE                ; {@BLOV}              {LVOC,  secondary oxidized biogenic}
BELV          = IGNORE                ; {@BELV}              {ELVOC, secondary oxidized biogenic}
ASOV          = IGNORE                ; {@ASOV}              {SVOC,  secondary oxidized aromatic}
ALOV          = IGNORE                ; {@ALOV}              {LVOC,  secondary oxidized aromatic}
AELV          = IGNORE                ; {@AELV}              {ELVOC, secondary oxidized aromatic}
PIOV          = IGNORE                ; {@PIOV}              {IVOC,  primary emitted n-alkane}
PSOV          = IGNORE                ; {@PSOV}              {SVOC,  primary emitted n-alkane}
PELV          = IGNORE                ; {@PELV}              {ELVOC, primary emitted n-alkane}

{------------------------------------- N -------------------------------------}

N2O3          =      3O + 2N          ; {@N_2O_3}            {dinitrogen trioxide}
N2O4          =      4O + 2N          ; {@N_2O_4}            {dinitrogen tetraoxide}

{------------------------------------- C - N --------------------------------}

{1C (amines)}

H2NCHO        = C + 3H + O  + N       ; {@H2NCHO}            {formamide}
MMA           = C + 5H      + N       ; {@MMA}               {methylamine}
CH2NH         = C + 3H      + N       ; {@CH2NH}             {methanimine} 
MMAO2         = C + 4H + 2O + N       ; {@MMAO2}             {MMA-peroxyradical}
CH3NH         = C + 4H      + N       ; {@CH3NH}             {N-radical of MMA}
MMNNO2        = C + 4H + 2O + 2N      ; {@MMNNO2}            {N-nitro methylamine}
CH3NO         = C + 3H + O  + N       ; {@CH3NO}             {nitroso methane}
HNCO          = C +  H + O  + N       ; {@HNCO}              {isocyanic acid}

{2C (amines)}

MEA           = 2C + 7H + O + N       ; {@MEA}               {monoethanolamin}
MEABO2        = 2C + 6H + 3O + N      ; {@MEABO2}            {C2-amine peroxy radical}
MEABO         = 2C + 6H + 2O + N      ; {@MEABO}             {C2-amine alkoxy radical}
MEAN          = 2C + 6H + O + N       ; {@MEAN}              {N-amine radical} 
H2NCOCHO      = 2C + 3H + 2O + N      ; {@H2NCOCHO}          {2-oxo acetamide}
H2NCH2CHO     = 2C + 5H + O + N       ; {@H2NCH2CHO}         {amino acetaldehyde}
H2NCOCH2OH    = 2C + 5H + 2O + N      ; {@H2NCOCH2OH}        {2-hydroxy acetamide}
HNCHCH2OH     = 2C + 5H + O + N       ; {@HNCHCH2OH}         {ethanol imine} 
H2NCHO2CHO    = 2C + 4H + 3O + N      ; {@H2NCHO2CHO}        {amino peroxy acetaldehyde}
H2NCH2CO3     = 2C + 4H + 3O + N      ; {@H2NCH2CO3}         {C2-amino peroxy acetyl radical}
H2NCOCO3      = 2C + 2H + 4O + N      ; {@H2NCOCO3}          {amido peroxy acetyl radical}
MEANNO2       = 2C + 6H + 3O + 2N     ; {@MEANNO2}           {N-nitroamino ethanol}
MEANHA        = 2C + 6H + 4O + 2N     ; {@MEANHA}            {N-nitro hydroxyacetamide}
MEANNO        = 2C + 6H + 2O + 2N     ; {@MEANNO}            {N-nitrosoamino ethanol}
DMA           = 2C + 7H      + N      ; {@DMA}               {dimethylamine} 
CH3NCH3       = 2C + 6H      + N      ; {@CH3NCH3}           {N-radical of DMA}  
CH2NCH3       = 2C + 5H      + N      ; {@CH2NCH3}           {N-methyl methanimine}
DMAO2         = 2C + 6H + 2O + N      ; {@DMAO2}             {DMA-peroxyradical}
NDMA          = 2C + 6H + O  + 2N     ; {@NDMA}              {N-nitroso dimethylamine}
DMNNO2        = 2C + 6H + 2O + 2N     ; {@DMNNO2}            {N-nitro dimethylamine}
CH3NHCHO      = 2C + 5H + O  + N      ; {@CH3NHCHO}          {N-methyl formamide}
HOCH2CH2NO    = 2C + 5H + 2O + N      ; {@HOCH2CH2NO}        {nitroso ethanol}
H2NCOCH3      = 2C + 5H + O  + N      ; {@H2NCOCH3}          {acetamide}

{3C (amines)}

TMA           = 3C + 9H      + N      ; {@TMA}               {trimethylamine} 
TMAO2         = 3C + 8H + 2O + N      ; {@TMAO2}             {TMA-peroxyradical} 
TMAO          = 3C + 8H + O  + N      ; {@TMAO}              {alkoxy-radical of TMA}
DMNCHO        = 3C + 7H + O  + N      ; {@DMNCHO}            {N,N-dimethyl formamide}
DMNCHOO2      = 3C + 6H + 3O + N      ; {@DMNCHOO2}          {peroxyradical of N,N-dimethyl formamide}
TMADF         = 3C + 5H + 2O + N      ; {@TMADF}             {N-methyl diformamide}
HOETNHCHO     = 3C + 7H + 2O + N      ; {@HOETNHCHO}         {ethanol amide}
HOCH2CONHCHO  = 3C + 5H + 3O + N      ; {@HOCH2CONHCHO}      {hydroxyaceto formamide}
DMCNH2        = 3C + 8H + N           ; {@DMCNH2}            {amino propyl radical}
DMCOONH2      = 3C + 8H + 2O + N      ; {@DMCOONH2}          {amino propyl peroxyradical}
CH2CNH2CH3    = 3C + 7H + N           ; {@CH2CNH2CH3}        {2-amino propene}
DMCNH         = 3C + 7H + N           ; {@DMCNH}             {2-propane imine}
CH3CNH2MOH    = 3C + 8H + O  + N      ; {@CH3CNH2MOH}        {amino propanol radical}
H2NCCHOHCH3   = 3C + 7H + O  + N      ; {@H2NCCHOHCH3}       {2-aminoprop-2-en-1-ol}
HNCCH3MOH     = 3C + 7H + O  + N      ; {@HNCCH3MOH}         {2-iminopropan-1-ol}
H2NCCH2MOH    = 3C + 7H + O  + N      ; {@H2NCCH2MOH}        {2-aminoprop-1-en-1-ol}

{3C (CHON)}

IPN           = 3C + 7H + 2O + N      ; {@IPN}               {isopropyl nitrite}
CH3CHOCH3     = 3C + 7H + O           ; {@CH3CHOCH3}         {isopropyloxy radical}
MGLYOAC       = 3C + 4H + 3O          ; {@MGLYOAC}           {CH3COCOOH = methylglyoxylic acid}

{4C (amines)}

DEA           = 4C + 11H + 2O + N     ; {@DEA}               {diethanolamine}
HOETNETOH     = 4C + 10H + 2O + N     ; {@HOETNETOH}         {N-radical of DEA}  
DEAO2         = 4C + 10H + 4O + N     ; {@DEAO2}             {DEA-peroxyradical}
HOETNHCH2CHO  = 4C + 10H + 2O + N     ; {@HOETNHCH2CHO}      {ethanolamine acetaldehyde}
NDELA         = 4C + 10H + 3O + 2N    ; {@NDELA}             {N-nitroso diethanolamine}
HOCH2CHNETOH  = 4C + 9H  + 2O + N     ; {@HOCH2CHNETOH}      {DEA imine}
DEANNO2       = 4C + 10H + 4O + 2N    ; {@DEANNO2}           {N-nitro diethanolamine} 
HOCH2CONETOH  = 4C + 8H  + 3O + N     ; {@HOCH2CONETOH}      {ethanol hydroxyacetamide}
AMP           = 4C + 11H + O + N      ; {@AMP}               {2-amino-2-methyl-1-propanol}
AMPN          = 4C + 10H + O + N      ; {@AMPN}              {N-radical of AMP}
NAMP          = 4C + 10H + 2O + 2N    ; {@NAMP}              {N-nitroso AMP}
AMPNNO2       = 4C + 10H + 3O + 2N    ; {@AMPNNO2}           {N-nitro AMP}
AMPOX         = 4C + 10H + 2O + N     ; {@AMPOX}             {AMP N-oxide}
DMCNH2CHO     = 4C + 9H + O + N       ; {@DMCNH2CHO}         {2-amino-2-methyl-1-propanal}
AMPNA         = 4C + 9H + 3O + 2N     ; {@AMPNA}             {N-Nitro-2-amino-2-methyl-1-propanal}
DMCNH2CO3     = 4C + 8H + 3O + N      ; {@DMCNH2CO3}         {AMP peroxy acetyl radical}
AMPAN         = 4C + 8H + 5O + 2N     ; {@AMPAN}             {AMP PAN-type compound}
AMPO          = 4C + 10H + 2O + N     ; {@AMPO}              {AMP alkoxy radical}
DMOCNH2MOH    = 4C + 9H  + 2O + N     ; {@DMOCNH2MOH}        {2-amino-3-hydroxy-2-methylpropanal}

{5C (amines)}

DEANCHO       = 5C + 11H + 3O + N     ; {@DEANCHO}           {N-diethanol formamide}
DEANCH2O2     = 5C + 12H + 4O + N     ; {@DEANCH2O2}         {N-diethanol formamide peroxyradical}

{6C (CHO)}

TME           = 6C + 12H              ; {@TME}               {Tetramethyl ethylene}
TMEO2         = 6C + 13H + 3O         ; {@TMEO2}             {Tetramethyl ethylene peroxide}
CHEX          = 6C + 12H              ; {@CHEX}              {Cyclohexane}
CHEXO2        = 6C + 11H + 2O         ; {@CHEXO2}            {Cyclohexane peroxyradical}
CHEXO         = 6C + 11H + O          ; {@CHEXO}             {Cyclohexane alkoxyradical}
CHEXOL        = 6C + 12H + O          ; {@CHEXOL}            {Cyclohexanol}
CHEXONE       = 6C + 11H + O          ; {@CHEXONE}           {Cyclohexone}
CHEXOOH       = 6C + 12H + 2O         ; {@CHEXOOH}           {Cyclohexane hydroperoxide}

{6C (amines)}

TEA           = 6C + 15H + 3O + N     ; {@TEA}               {triethanolamin}
TEAO2         = 6C + 14H + 5O + N     ; {@TEAO2}             {TEA-peroxyradical}
TEAO          = 6C + 14H + 4O + N     ; {@TEAO}              {TEA-alkoxyradical} 
DEANCOCH2OH   = 6C + 13H + 4O + N     ; {@DEANCOCH2OH}       {N,N-diethanol hydroxyacetamide}
DEANCH2CHO    = 6C + 13H + 3O + N     ; {@DEANCH2CHO}        {N,N-diethanol acetamide}
DEANCH2COO2   = 6C + 12H + 5O + N     ; {@DEANCH2COO2}       {N,N-diethanol acetamide peroxyradical}



{------------------------------------- S ------------------------------------}

HSO3          = S + H + 3O            ; {@HSO_3}             {sulfonic acid}
CH3SOH        = C + 4H + S + O        ; {@CH_3SOH}           {MSEA}
CH3SOOH       = C + 4H + S + 2O       ; {@CH_3SOOH}          {MSIA: methane sulfinic acid}
CH3SOO2H      = C + 4H + S + 3O       ; {@CH_3SOO_2H}        {}
CH3SO4H       = C + 4H + S + 4O       ; {@CH_3SO_4H}         {}
CH3SCH2       = 2C + 5H + S           ; {@CH_3SCH_2}         {dimethyl sulfide radical}
DMSOO         = 2C + 5H + S + 2O      ; {@CH_3SCH_2OO}       {dimethyl sulfide peroxyradical}
DMSOOH        = 2C + 6H + S + 2O      ; {@CH_3SCH_2OOH}      {dimethyl sulfide hydroperoxide}
DMSOH         = 2C + 7H + S + O       ; {@DMSOH}             {dimethyl sulfhydroxide: CH3SOHCH3}
DMSOHO        = 2C + 7H + S + 2O      ; {@DMSOHO}            {}
DMSOHOO       = 2C + 7H + S + 3O      ; {@DMSOHOO}           {}
CH3SOCH2      = 2C + 5H + S + O       ; {@CH_3SOCH_2}        {dimethyl sulfoxide radical}
DMSOOO        = 2C + 5H + S + 3O      ; {@CH_3SOCH_2O_2}     {dimethyl sulfoxide peroxyradical}
DMSO2         = 2C + 6H + S + 2O      ; {@DMSO_2}            {dimethyl sulfone: CH3SO2CH3}
DMSO2O        = 2C + 6H + S + 3O      ; {@DMSO_2O}           {dimethyl sulfone oxyradical}
DMSO2OO       = 2C + 6H + S + 4O      ; {@DMSO_2OO}          {dimethyl sulfone peroxyradical}
DMSO2OOH      = 2C + 6H + S + 4O      ; {@DMSO_2OOH}         {dimethyl sulfone hydroperoxide}
CH3S          = C + 3H + S            ; {@CH_3S}             {}
CH3SO         = C + 3H + S + O        ; {@CH_3SO}            {}
CH3SOO        = C + 3H + S + 2O       ; {@CH_3SOO}           {}
CH3SOO2       = C + 3H + S + 3O       ; {@CH_3SOO_2}         {}
CH3SO4        = C + 3H + S + 4O       ; {@CH_3SO_2O_2}       {}
MSON          = C + 3H + S + 4O + N   ; {@CH_3SOONO_2}       {}
MSOON         = C + 3H + S + 5O + N   ; {@CH_3SOO_2NO_2}     {}
MSPN          = C + 3H + S + 6O + N   ; {@CH_3SO_2O_2NO_2}   {methyl sulfonyl peroxynitrate}
MSAH2O        = C + 7H + S + 4O       ; {@MSA(H_2O)}         {[MSA*H2O]: methane sulfonic acid - water cluster}
MSADMAH2O     = 3C+ 14H+ S + 4O + N   ; {@MSA(DMA)(H_2O)}    {[MSA*DMA*H2O]: methane sulfonic acid - DMA - water cluster}
MSADMA        = 3C+ 11H+ S + 3O + N   ; {@MSA(DMA)}          {[MSA*DMA]: methane sulfonic acid - DMA cluster} 
MSATMAH2O     = 4C+ 16H+ S + 4O + N   ; {@MSA(TMA)(H_2O)}    {[MSA*TMA*H2O]: methane sulfonic acid - TMA - water cluster}
MSATMA        = 4C+ 13H+ S + 3O + N   ; {@MSA(TMA)}          {[MSA*TMA]: methane sulfonic acid - TMA cluster} 

