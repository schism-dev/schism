
## State variables and sub-modules in ICM
<pre>
Core Module
     1  PB1   :  Diatom                                     g/m^3
     2  PB2   :  Green Algae                                g/m^3
     3  PB3   :  Cyanobacteria                              g/m^3
     4  RPOC  :  Refractory Particulate Organic Carbon      g/m^3
     5  LPOC  :  Labile Particulate Organic Carbon          g/m^3
     6  DOC   :  Dissolved Orgnaic Carbon                   g/m^3
     7  RPON  :  Refractory Particulate Organic Nitrogen    g/m^3
     8  LPON  :  Labile Particulate Organic Nitrogen        g/m^3
     9  DON   :  Dissolved Orgnaic Nitrogen                 g/m^3
     10 NH4   :  Ammonium Nitrogen                          g/m^3
     11 NO3   :  Nitrate Nitrogen                           g/m^3
     12 RPOP  :  Refractory Particulate Organic Phosphorus  g/m^3
     13 LPOP  :  Labile Particulate Organic Phosphorus      g/m^3
     14 DOP   :  Dissolved Orgnaic Phosphorus               g/m^3
     15 PO4   :  Total Phosphate                            g/m^3
     16 COD   :  Chemical Oxygen Demand                     g/m^3
     17 DOX   :  Dissolved Oxygen                           g/m^3
Silica Module
     1  SU    :  Particulate Biogenic Silica                g/m^3
     2  SA    :  Available Silica                           g/m^3
Zooplankton Module
     1  ZB1   :  1st zooplankton                            g/m^3
     2  ZB2   :  2nd zooplankton                            g/m^3
pH Module
     1  TIC   :  Total Inorganic Carbon                     g/m^3
     2  ALK   :  Alkalinity                                 g[CaCO3]/m^3
     3  CA    :  Dissolved Calcium                          g[CaCO3]/m^3
     4  CACO3 :  Calcium Carbonate                          g[CaCO3]/m^3
CBP Module
     1  SRPOC :  Slow Refractory Particulate Organic Carbon g/m^3
     2  SRPON :  Slow Refractory Particulate Organic Nitro. g/m^3
     3  SRPOP :  Slow Refractory Particulate Organic Phosp. g/m^3
     4  PIP   :  Particulate Inorganic Phosphate            g/m^3
SAV Module (no transport variables)
VEG Module (no transport variables)
SFM Module (no transport variables)
</pre>

## 1. Core Module

### 1.1 Mass Balance Equations of State Variables in ICM

**Note: ** <br>
>> a). settling of variables are addressed in separate section, and not included in the equations below <br>
>>  b). When iCBP=0, equations of ($\text{SRPOC},\text{SRPON},\text{SRPOP},\text{PIP}$) are omiited, 
and $(FCP_4^m,FCM_4^m,FNP_5^m,FNM_5^m,FPP_5^m,FPM_5^m)=0, (KC_S,KN_S,KP_S,KP_P)=0$ 

**Phytoplankton (PB1, PB2, PB3)**:

\begin{flalign}
  & d\text{PB}^i = \text{GP}^i-\text{MT}^i-\text{PR}^i \text{ ,   i=1,3}  \\
\end{flalign}


**Carbon (RPOC, LPOC, DOC)**:

\begin{flalign}
  & d\text{RPOC} = -KC_1 \cdot \text{RPOC}+ \sum_{m=1,3} \left( FCP_1^m \cdot \text{PR}^m + FCM_1^m \cdot \text{MT}^m \right)\\
  & d\text{LPOC} = -KC_2 \cdot \text{LPOC}+ \sum_{m=1,3} \left( FCP_2^m \cdot \text{PR}^m + FCM_2^m \cdot \text{MT}^m \right) \\
  & \begin{split}
    & d\text{DOC}=KC_1 \cdot \text{RPOC}+ KC_2 \cdot \text{LPOC} + KC_S \cdot \text{SRPOC} -K_{HR} \cdot \text{DOC} -Denit \cdot \text{DOC} \\
    & + \sum_{m=1,3} FCP_3^m \cdot \text{PR}^m + \sum_{m=1,3} \left[FCM_3^m+ \left(1-\sum_{i=1,4} FCM_i^m \right) \cdot \frac{KhDO^m}{DO+KhDO^m} \right] \cdot \text{MT}^m \\
  & \end{split} \\
  & d\text{SRPOC} = -KC_S \cdot \text{SRPOC}+ \sum_{m=1,3} \left( FCP_4^m \cdot \text{PR}^m + FCM_4^m \cdot \text{MT}^m \right)\\
\end{flalign}


**Nitrogen (RPON, LPON, DON, NH4, NO3)**:

\begin{flalign}
  & d\text{RPON} = -KN_1 \cdot \text{RPON}+ \sum_{m=1,3} n2c^m \cdot \left( FNP_1^m \cdot \text{PR}^m + FNM_1^m \cdot \text{MT}^m \right) \\
  & d\text{LPON} = -KN_2 \cdot \text{LPON}+ \sum_{m=1,3} n2c^m \cdot \left( FNP_2^m \cdot \text{PR}^m + FNM_2^m \cdot \text{MT}^m \right) \\
  & \begin{split} 
     & d\text{DON} = KN_1 \cdot \text{RPON} + KN_2 \cdot \text{LPON} + KN_S \cdot \text{SRPON} -KN_3 \cdot \text{DON} \\
     & + \sum_{m=1,3} n2c^m \cdot \left( FNP_3^m \cdot \text{PR}^m + FNM_3^m \cdot \text{MT}^m \right) 
  & \end{split} \\
  & d\text{NH4} = KN_3 \cdot \text{DON}-Nit \cdot \text{NH4}+ \sum_{m=1,3} n2c^m \cdot \left( FNP_4^m \cdot \text{PR}^m + FNM_4^m \cdot \text{MT}^m -fPN^m \cdot \text{GP}^m \right) \\
  & d\text{NO3} = Nit \cdot \text{NH4}-dn2c \cdot Denit \cdot \text{DOC}-\sum_{m=1,3} n2c^m \cdot \left( 1-fPN^m \right) \cdot \text{GP}^m  \\
  & d\text{SRPON} = -KN_S \cdot \text{SRPON}+ \sum_{m=1,3} n2c^m \cdot \left( FNP_5^m \cdot \text{PR}^m + FNM_5^m \cdot \text{MT}^m \right)\\
\end{flalign}


**Phosphorus (RPOP, LPOP, DOP, PO4)**:

\begin{flalign}
  & d\text{RPOP}= -KP_1 \cdot \text{RPOP}+ \sum_{m=1,3} p2c^m \cdot \left( FPP_1^m \cdot \text{PR}^m + FPM_1^m \cdot \text{MT}^m \right) \\
  & d\text{LPOP}= -KP_2 \cdot \text{LPOP}+ \sum_{m=1,3} p2c^m \cdot \left( FPP_2^m \cdot \text{PR}^m + FPM_2^m \cdot \text{MT}^m \right) \\
  & \begin{split} 
     & d\text{DOP} = KP_1 \cdot \text{RPOP} + KP_2 \cdot \text{LPOP} + KP_S \cdot \text{SRPOP} -KP_3 \cdot \text{DOP}  \\
     & + \sum_{m=1,3} p2c^m \cdot \left( FPP_3^m \cdot \text{PR}^m + FPM_3^m \cdot \text{MT}^m \right) 
  & \end{split} \\
  & d\text{PO4} = KP_3 \cdot \text{DOP} + KP_P \cdot \text{PIP} + \sum_{m=1,3} p2c^m \cdot \left( FPP_4^m \cdot \text{PR}^m + FPM_4^m \cdot \text{MT}^m - \text{GP}^m \right) \\
  & d\text{SRPOP} = -KP_S \cdot \text{SRPOP}+ \sum_{m=1,3} p2c^m \cdot \left( FPP_5^m \cdot \text{PR}^m + FPM_5^m \cdot \text{MT}^m \right)\\
  & d\text{PIP} = -KP_P \cdot \text{PIP} \\
\end{flalign}

**Oxygen (COD, DO)**:

\begin{flalign}
 & d\text{COD} = -K_{COD} \cdot \text{COD}  \\
 & \begin{split} \\
   & d\text{DO} = -o2n \cdot Nit \cdot \text{NH4} -o2c \cdot K_{HR} \cdot \text{DOC} -K_{COD} \cdot \text{COD} \\
   & + \sum_{m=1,3} o2c \cdot \left[ \left(1.3-0.3 \cdot fPN^m \right) \cdot \text{GP}^m -\left(1- \sum_{i=1,4} FCM_i^m \right) \cdot \frac{DO}{DO+KhDO^m} \cdot \text{MT}^m \right] \\  
 & \end{split} \\
\end{flalign}


### 1.2. Pre-Calculation

#### 1.2.1. Growth, metabolism, predation
\begin{flalign}
  & \text{GP}^i = \text{GPM}^i \cdot f(T) \cdot f(Sal) \cdot f(I) \cdot \text{min} \left[ f(N),f(P),f(S) \right] \cdot \text{PB}^i \\
  & \text{MT}^i = \text{MTR}^i \cdot \text{GP} +\text{MTB}^i \cdot \text{exp} \left[ KT_{MT}^i \cdot \left( T-T_{MT}^i \right) \right] \cdot \text{PB}^i \\
  & \text{PR}^i =
\begin{cases} 
  & \text{PRR}^i \cdot \text{exp} \left[ KT_{MT}^i \cdot \left( T-T_{MT}^i \right) \right] \cdot \text{PB}^i \text{,        iPR=0} \\ 
  & \text{PRR}^i \cdot \text{exp} \left[ KT_{MT}^i \cdot \left( T-T_{MT}^i \right) \right] \cdot \left( \text{PB}^i \right)^2 \text{,   iPR=1} \\ 
\end{cases} \\ 
\end{flalign}

\begin{flalign}
  & f(I)= \frac{\text{mLight}}{\sqrt{\text{mLight}^2+IK^2}}   \\
  & f(T)= 
\begin{cases} 
   \text{exp}\left[-KTGP_1^i \cdot \left( T-TGP^i \right)^2 \right] \text{, if }T < TGP^i \\
   \text{exp}\left[-KTGP_2^i \cdot \left( T-TGP^i \right)^2 \right] \text{, if }T \geq TGP^i \\
\end{cases} \\ 
  & f(N)= \frac{\text{DIN}}{\text{DIN}+KhN^i} \\
  & f(P)= \frac{\text{PO4d}}{\text{PO4d}+KhP^i} \\
  & f(Sal)= \frac{KhSal_i^2}{KhSal_i^2+Sal^2} \\
\end{flalign}

#### 1.2.2. Decay rates of orgnaic matter
\begin{flalign}
  & KC_i = \left( KC_i^0+KC_i^{alg} \cdot APB \right) \cdot KT_M \\
  & KN_i = \left( KN_i^0+KN_i^{alg} \cdot APB \cdot \frac{mKhN}{mKhN+DIN} \right) \cdot KT_M \\
  & KP_i = \left( KP_i^0+KP_i^{alg} \cdot APB \cdot \frac{mKhP}{mKhP+PO4_d} \right) \cdot KT_M \\
\end{flalign}

\begin{flalign}
  KT_M=\text{exp} \left[ KT_{RM}^i \cdot \left( T-T_{RM}^i \right) \right]
\end{flalign}

#### 1.2.3. Respiration, denitrification, decay of COD, nitrification
\begin{flalign}
  & K_{HR} = KC_3 \cdot \frac{DO}{KhDO_{OX}+DO} \\ 
  & K_{COD}= \frac{DO}{KhCOD+DO} \cdot KCD \cdot \text{exp}[KT_{RCOD} \cdot (T-T_{RCOD})] \\
  & Denit = an2c \cdot KC_3 \cdot \frac{KhDO_{OX}}{KhDO_{OX}+DO} \cdot \frac{NO3}{KhNO3_{dn}+NO3} \\
  & Nit = Nit^{max} \cdot \frac{DO}{KhDO_n+DO} \cdot \frac{KhNH4_n}{KhNH4_n+NH4} \cdot 
\begin{cases}
  \text{exp}[-KT_{Nit}^1 \cdot (T_{Nit}-T)] \text{, if } T<T_{Nit} \\
  \text{exp}[-KT_{Nit}^2 \cdot (T-T_{Nit})] \text{, if } T \geq T_{Nit}\\
\end{cases}
\end{flalign}

### 1.3 Light

### 1.4 Surface/bottom fluxes

#### DO reareation

### 1.4 2D spatially varying parameters

### 1.6 TSS



# [Old User Guide of ICM](http://ccrm.vims.edu/schismweb/ICM_UserGuide_v2.pdf).

