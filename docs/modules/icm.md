
\usepackage{setspace}

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
SAV Module (no transport variables)
VEG Module (no transport variables)
SFM Module (no transport variables)
</pre>

## 1. Core Module

### 1.1. Pre-Calculation

#### 1.1.1. Growth, metabolism, predation
\begin{flalign}
  & GP^i = GPM^i \cdot f(T) \cdot f(Sal) \cdot f(I) \cdot \text{min} \left[ f(N),f(P),f(S) \right] \cdot PB^i \\
  & MT^i = MTB^i \cdot \text{exp} \left[ KT_{MT}^i \cdot \left( T-T_{MT}^i \right) \right] \cdot PB^i \\
  & PR^i = PRR^i \cdot \text{exp} \left[ KT_{MT}^i \cdot \left( T-T_{MT}^i \right) \right] \cdot PB^i \\ 
\end{flalign}

\begin{flalign}
  & f(T)
\end{flalign}

#### 1.1.2. Decay rates of orgnaic matter
\begin{flalign}
  & KC^i = \left( KC_0^i+KC_{alg}^i \cdot APB \right) \cdot KT_M \\
  & KN^i = \left( KN_0^i+KN_{alg}^i \cdot APB \cdot \frac{mKhN}{mKhN+DIN} \right) \cdot KT_M \\
  & KP^i = \left( KP_0^i+KP_{alg}^i \cdot APB \cdot \frac{mKhP}{mKhP+PO4_d} \right) \cdot KT_M \\
\end{flalign}

\begin{align}
  KT_M=\text{exp} \left[ KT_{RM}^i \cdot \left( T-T_{RM}^i \right) \right]
\end{align}

#### 1.1.3. Respiration, denitrification, decay of COD, nitrification
\begin{flalign}
  & K_{HR} = KC^3 \cdot \frac{DO}{KhDO_{OX}+DO} \\ 
  & K_{COD}= \frac{DO}{KhCOD+DO} \cdot KCD \cdot \text{exp}[KT_{RCOD} \cdot (T-T_{RCOD})] \\
  & Denit = an2c \cdot KC^3 \cdot \frac{KhDO_{OX}}{KhDO_{OX}+DO} \cdot \frac{NO3}{KhNO3_{dn}+NO3} \\
  & Nit = Nit^{max} \cdot \frac{DO}{KhDO_n+DO} \cdot \frac{KhNH4_n}{KhNH4_n+NH4} \cdot 
\begin{cases}
  \text{exp}[-KT_{Nit}^1 \cdot (T_{Nit}-T)] \text{, if } T<T_{Nit} \\
  \text{exp}[-KT_{Nit}^2 \cdot (T-T_{Nit})] \text{, if } T>=T_{Nit}\\
\end{cases}
\end{flalign}

\begin{flalign}
   
\end{flalign}


### **Light**

### **Surface/bottom fluxes**

#### DO reareation

### **2D spatially varying parameters**

### **TSS**
  




# [Old User Guide of ICM](http://ccrm.vims.edu/schismweb/ICM_UserGuide_v2.pdf).
