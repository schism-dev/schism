************************************************************************
**                                                                    **
**                  Chesapeake Bay Sediment Model                     **
**                                                                    **
** Third version (speeded up) received from Fitzpatrick May 8, 1990   **
**              First modified by C. Cerco May 9, 1990                **
**           Final modifications D.M. Di Toro, Jan 27, 1992           **
** org flux to sediment ==> diagenesis ==> inorg flux to water column **
**                                                                    **
**     3-G model - G1=labile, G2=refractory, G3=slow refractory       **
**                                                                    **
************************************************************************
**                                                                    **
** Inputs                                                             **
**                                                                    **
**            Required inputs for sediment sub-model                  **
**                                                                    **
** A.  Passed to sediment subroutine from water quality subroutine    **
**                                                                    **
**   1.  Overlying water column segment volume (V1)                   **
**   2.  Overlying water column segment depth (BL(ISEG,3))            **
**   3.  Overlying water column depth                                 **
**   4.  Overlying water column segment temperature and salinity      **
**   5.  Overlying water column ammonia, nitrate, phosphate, silica   **
**       and dissolved oxygen concentrations                          **
**                                                                    **
** B. Inputs supplied via direct input to the sediment subroutine     **
**                                                                    **
**  Variable names        Description                         Units   **
**                                                                    **
**                                                                    **
**    HSEDALL       Depth of sediment layer  (h2)               cm    **
**      DIFFT       Water column-sediment layer diffusion     cm2/sec **
**                  coefficient for temperature                       **
**     SALTSW       Salinity concentration for determining      ppt   **
**                  whether methane or sulfide SOD formulation        **
**                  is to be used.  Also determines PO4 sorption.     **
**     SALTND       Determines whether fresh or saltwater             **
**                  nitrification/denitrification rates are used      **
**                                                                    **
**   Diagenesis stoichiometry                                         **
**                                                                    **
**                 Fractions of G1, G2, and G3 contained in ...       **
**                                                                    **
**   FRPPH1(3)     Algal group no.1 phosphorus                        **
**   FRPPH2(3)     Algal group no.2 phosphorus                        **
**   FRPPH3(3)     Algal group no.3 phosphorus                        **
**   FRPOP(NSED,3) Non-algal particulate organic phosphorus           **
**   FRNPH1(3)     Algal group no.1 nitrogen                          **
**   FRNPH2(3)     Algal group no.2 nitrogen                          **
**   FRNPH3(3)     Algal group no.3 nitrogen                          **
**   FRPON(NSED,3) Non-algal particulate organic nitrogen             **
**   FRCPH1(3)     Algal group no.1 carbon                            **
**   FRCPH2(3)     Algal group no.2 carbon                            **
**   FRCPH3(3)     Algal group no.3 carbon                            **
**   FRPOC(NSED,3) Non-algal particulate organic carbon               **
**                                                                    **
**   Diagenesis kinetics                                              **
**                                                                    **
**   KPDIAG(3)     Reaction rates for POP G1, G2, and G3        /day  **
**   DPTHTA(3)     Temperature thetas for POP G1, G2, and G3          **
**   KNDIAG(3)     Reaction rates for PON G1, G2, and G3        /day  **
**   DNTHTA(3)     Temperature thetas for PON G1, G2, and G3          **
**   KCDIAG(3)     Reaction rates for POC G1, G2, and G3        /day  **
**   DCTHTA(3)     Temperature thetas for POC G1, G2, and G3          **
**   KSI           Reaction rate for Part. Biogenic Si (PSi)    /day  **
**   THTASI        Temperature theta for PSi                          **
**                                                                    **
**   Solids and transport                                             **
**                                                                    **
**   VPMIX(NSED)   Particulate diffusion coefficient (Dp)   m**2/day  **
**   THTADP        Temperature theta for Dp                           **
**   VDMIX(NSED)   Porewater diffusion coefficient (Dd)     m**2/day  **
**   THTADD        Temperature theta for Dd                           **
**      M1         Concentration of solids in layer 1       kg/l      **
**      M2         Concentration of solids in layer 2       kg/l      **
**                                                                    **
**   Reaction kinetics                                                **
**                                                                    **
**    KAPPNH4F     Nitrification reaction velocity                    **
**                 for freshwater in layer 1                m/day     **
**    KAPPNH4S     Nitrification reaction velocity                    **
**                 for saltwater  in layer 1                m/day     **
**    PIENH4       Ammonia partition coefficient            L/kg      **
**    THTANH4      Theta for nitrification reaction velicities        **
**    KMNH4        Nitrification half saturation constant             **
**                 for ammonia                              mg N/m3   **
**    KMNH4O2      Nitrification half saturation constant             **
**                 for oxygen                               mg O2/L   **
**    KAPPNO3F     Denitrification reaction velocity                  **
**                 for freshwater in layer 1                m/day     **
**    KAPPNO3S     Denitrification reaction velocity                  **
**                 for saltwater  in layer 1                m/day     **
**    K2NO3        Denitrification reaction velocity                  **
**                 in layer 2                               m/day     **
**    THTANO3      Theta for denitrification                          **
**    KAPPD1       Reaction velocity for dissolved sulfide            **
**                 oxidation in layer 1                     m/day     **
**    KAPPP1       Reaction velocity for particulate sulfide          **
**                 oxidation in layer 1                     m/day     **
**    PIE1S        Partition coefficient for sulfide                  **
**                 in layer 1                               L/kg      **
**    PIE2S        Partition coefficient for sulfide                  **
**                 in layer 2                               L/kg      **
**    THTAPD1      Theta for both dissolved and particulate           **
**                 sulfide oxidation                                  **
**    KMHSO2       Sulfide oxidation normalization constant           **
**                 for oxygen                               mg O2/L   **
**    CSISAT       Saturation concentration for pore water            **
**                 silica                                   mg Si/m3  **
**    DPIE1SI      Incremental partition coefficient for              **
**                 silica in layer 1                        L/kg      **
**    PIE2SI       Partition coefficient for silica in                **
**                 layer 2                                  L/kg      **
**    KMPSI        Particulate biogenic silica half saturation        **
**                 constant for dissolution                 mg Si/m3  **
**    O2CRITSI     Critical dissolved oxygen concentration            **
**                 for layer 1 incremental silica sorption  mg O2/L   **
**    JSIDETR      Detrital biogenic silica source to                 **
**                 sediment                                 mg Si/m2-d**
**    DPIE1PO4F    Incremental partition coefficient for              **
**                 phosphate in layer 1 (freshwater)        L/kg      **
**    DPIE1PO4S    Incremental partition coefficient for              **
**                 phosphate in layer 1 (saltwater)         L/kg      **
**    PIE2PO4      Partition coefficient for phosphate                **
**                 in layer 2                               L/kg      **
**    O2CRIT       Critical dissolved oxygen concentration for        **
**                 layer 1 incremental phosphate sorption   mg O2/L   **
**    KMO2DP       Particle mixing half saturation constant           **
**                 for oxygen                               mg O2/L   **
**    TEMPBEN      Temperature at which benthic stress                **
**                 accumulation is reset to zero            deg C     **
**    KBENSTR      Decay constant for benthic stress        /day      **
**    KLBNTH       Ratio of bio-irrigation to bioturbation            **
**    DPMIN        Minimum particle diffusion coefficient   m2/day    **
**    KAPPCH4      methane oxidation reaction velocity      m/day     **
**    THTACH4      theta for methane oxidation                        **
**                                                                    **
** Output                                                             **
**                                                                    **
**    The subroutine returns fluxes for                               **
**                                                                    **
**     JSOD, JAQSOD, JCH4AQ and JCH4G  [gm o2*/m2-day]                **
**     JNH4, JPO4, JNO3 and JSI  [mg/m2-day]                          **
**                                                                    **
**    via array BFLUX in COMMON /BENTHC/                              **
**                                                                    **
************************************************************************


      SUBROUTINE SED_READ
      SAVE

      INCLUDE  'wqm_com.inc'

***** Variable declarations

      CHARACTER FRNAME(14)*24

***** Data declarations

      DATA FRNAME
     .     /'Cyanobacteria phosphorus', 'Diatom phosphorus       ',
     .      'Green algal phosphorus  ', 'Detrital org phosphorus ',
     .      'Cyanobacteria nitrogen  ', 'Diatom nitrogen         ',
     .      'Green algal nitrogen    ', 'Detrital org nitrogen   ',
     .      'Cyanobacteria carbon    ', 'Diatom carbon           ',
     .      'Green algal carbon      ', 'Detrital org carbon     ',
     .      'Diatom silica           ', 'Particulate biog silica '/

********************************************************************************
**                                  Inputs                                    **
********************************************************************************

      READ(BFI,1000,ERR=10100)  HSEDALL, INTSEDC
      READ(BFI,1010,ERR=10100)  DIFFT
      READ(BFI,1010,ERR=10100)  SALTSW, SALTND
      READ(BFI,1010,ERR=10100)  FRPPH1
      READ(BFI,1010,ERR=10100)  FRPPH2
      READ(BFI,1010,ERR=10100)  FRPPH3
      READ(BFI,1010,ERR=10100)  FRNPH1
      READ(BFI,1010,ERR=10100)  FRNPH2
      READ(BFI,1010,ERR=10100)  FRNPH3
      READ(BFI,1010,ERR=10100)  FRCPH1
      READ(BFI,1010,ERR=10100)  FRCPH2
      READ(BFI,1010,ERR=10100)  FRCPH3
      READ(BFI,1010,ERR=10100) (KPDIAG(JG),DPTHTA(JG),JG=1,3)
      READ(BFI,1010,ERR=10100) (KNDIAG(JG),DNTHTA(JG),JG=1,3)
      READ(BFI,1010,ERR=10100) (KCDIAG(JG),DCTHTA(JG),JG=1,3)
      READ(BFI,1010,ERR=10100)  KSI,THTASI
      READ(BFI,1010,ERR=10100)  M1,M2,THTADP,THTADD
      READ(BFI,1010,ERR=10100)  KAPPNH4F,KAPPNH4S,PIENH4,THTANH4,KMNH4,
     .                          KMNH4O2
      READ(BFI,1010,ERR=10100)  KAPPNO3F,KAPPNO3S,K2NO3,THTANO3
      READ(BFI,1010,ERR=10100)  KAPPD1,KAPPP1,PIE1S,PIE2S,THTAPD1,
     .                          KMHSO2
      READ(BFI,1010,ERR=10100)  CSISAT,DPIE1SI,PIE2SI,KMPSI
      READ(BFI,1010,ERR=10100)  O2CRITSI,JSIDETR
      READ(BFI,1010,ERR=10100)  DPIE1PO4F,DPIE1PO4S,PIE2PO4,O2CRIT,
     .                          KMO2DP
      READ(BFI,1010,ERR=10100)  TEMPBEN,KBENSTR,KLBNTH,DPMIN
      READ(BFI,1010,ERR=10100)  KAPPCH4,THTACH4
      READ(BFI,1022,ERR=10100)  BFI_TYPE
      IF (.NOT.((BFI_TYPE.EQ.' UNIFORM').OR.
     .            (BFI_TYPE.EQ.'  VARIED'))) THEN
          ERROR_CHECK(32) = .TRUE.
          OPEN (ERR,FILE=ERRFN,STATUS='UNKNOWN')
          IF (ERROR_CHECK(32)) WRITE (ERR,*) ERROR_MESS(32)
          IF (ERROR_CHECK(32)) WRITE (*,*)   ERROR_MESS(32)
          STOP
      END IF
      IF (BFI_TYPE.EQ.' UNIFORM') THEN
        READ(BFI,1020,ERR=10100)  WSSNET(1),WSLNET(1),WSRNET(1),
     .                            WSCNET(1),WSDNET(1),WSGNET(1),
     .                            VSED(1),VPMIX(1),VDMIX(1)
        READ(BFI,1040,ERR=10100)  FRPOP(1,2),FRPOP(1,3),FRPON(1,2),
     .                            FRPON(1,3),FRPOC(1,2),FRPOC(1,3)

        DO BB=2,NBB
          WSSNET(BB)  = WSSNET(1)
          WSLNET(BB)  = WSLNET(1)
          WSRNET(BB)  = WSRNET(1)
          WSCNET(BB)  = WSCNET(1)
          WSDNET(BB)  = WSDNET(1)
          WSGNET(BB)  = WSGNET(1)
          VSED(BB)    = VSED(1)
          VPMIX(BB)   = VPMIX(1)
          VDMIX(BB)   = VDMIX(1)
          FRPOP(BB,2) = FRPOP(1,2)
          FRPOP(BB,3) = FRPOP(1,3)
          FRPON(BB,2) = FRPON(1,2)
          FRPON(BB,3) = FRPON(1,3)
          FRPOC(BB,2) = FRPOC(1,2)
          FRPOC(BB,3) = FRPOC(1,3)
        END DO
      ELSE IF (BFI_TYPE.EQ.'  VARIED') THEN
        READ(BFI,1020,ERR=10100) (WSSNET(BB),WSLNET(BB),WSRNET(BB),
     .                            WSCNET(BB),WSDNET(BB),WSGNET(BB),
     .                            VSED(BB),VPMIX(BB),VDMIX(BB),BB=1,NBB)
        READ(BFI,1040,ERR=10100) (FRPOP(BB,2),FRPOP(BB,3),FRPON(BB,2),
     .                            FRPON(BB,3),FRPOC(BB,2),FRPOC(BB,3),
     .                            BB=1,NBB)
      END IF
 
***** Define logical variables

      STEADY_STATE_SED = INTSEDC.EQ.1

********************************************************************************
**                                 Outputs                                    **
********************************************************************************

      IF (BENTHIC_OUTPUT) THEN
        OPEN (BFO,FILE=BFOFN,STATUS='NEW')
        WRITE(BFO,2000)
        WRITE(BFO,2020) HSEDALL
        IF (STEADY_STATE_SED) THEN
          WRITE(BFO,2022)
        ELSE
          WRITE(BFO,2025)
        END IF
        WRITE(BFO,2030)
        WRITE(BFO,2040)  SSNAME(1)
        WRITE(BFO,2050) (CTEMP(BB),BB=1,NBB)
        WRITE(BFO,2060)  SSNAME(2)
        WRITE(BFO,2070) ((CPOP(BB,JG),JG=1,3),BB=1,NBB)
        WRITE(BFO,2060)  SSNAME(3)
        WRITE(BFO,2070) ((CPON(BB,JG),JG=1,3),BB=1,NBB)
        WRITE(BFO,2060)  SSNAME(4)
        WRITE(BFO,2070) ((CPOC(BB,JG),JG=1,3),BB=1,NBB)
        WRITE(BFO,2060)  SSNAME(5)
        WRITE(BFO,2070) (CPOS(BB),BB=1,NBB)
        WRITE(BFO,2040)  SSNAME(6)
        WRITE(BFO,2050) (PO4T2TM1S(BB),BB=1,NBB)
        WRITE(BFO,2040)  SSNAME(7)
        WRITE(BFO,2050) (NH4T2TM1S(BB),BB=1,NBB)
        WRITE(BFO,2040)  SSNAME(8)
        WRITE(BFO,2050) (NO3T2TM1S(BB),BB=1,NBB)
        WRITE(BFO,2040)  SSNAME(9)
        WRITE(BFO,2050) (HST2TM1S(BB),BB=1,NBB)
        WRITE(BFO,2040)  SSNAME(10)
        WRITE(BFO,2050) (SIT2TM1S(BB),BB=1,NBB)
        WRITE(BFO,2040)  SSNAME(11)
        WRITE(BFO,2050) (BENSTR1S(BB),BB=1,NBB)
        WRITE(BFO,2080)  0.0001*DIFFT
        WRITE(BFO,2090)  SALTSW,SALTND
        WRITE(BFO,2100)
        WRITE(BFO,2110)  FRNAME(1),FRPPH1
        WRITE(BFO,2110)  FRNAME(2),FRPPH2
        WRITE(BFO,2110)  FRNAME(3),FRPPH3
        WRITE(BFO,2110)  FRNAME(5),FRNPH1
        WRITE(BFO,2110)  FRNAME(6),FRNPH2
        WRITE(BFO,2110)  FRNAME(7),FRNPH3
        WRITE(BFO,2110)  FRNAME(9),FRCPH1
        WRITE(BFO,2110)  FRNAME(10),FRCPH2
        WRITE(BFO,2110)  FRNAME(11),FRCPH3
        WRITE(BFO,2250) (BBN(BB),FRPOP(BB,2),FRPOP(BB,3),FRPON(BB,2),
     .                   FRPON(BB,3),FRPOC(BB,2),FRPOC(BB,3),BB=1,NBB)
        WRITE(BFO,2120) (KPDIAG(JG),DPTHTA(JG),JG=1,3),(KNDIAG(JG),
     .                   DNTHTA(JG),JG=1,3),(KCDIAG(JG),DCTHTA(JG),
     .                   JG=1,3),KSI,THTASI
        WRITE(BFO,2170)  M1,M2,THTADP,THTADD
        WRITE(BFO,2180)  KAPPNH4F,KAPPNH4S,PIENH4,THTANH4,KMNH4,KMNH4O2
        WRITE(BFO,2190)  KAPPNO3F,KAPPNO3S,K2NO3,THTANO3
        WRITE(BFO,2200)  KAPPD1,KAPPP1,PIE1S,PIE2S,THTAPD1,KMHSO2
        WRITE(BFO,2210)  CSISAT,DPIE1SI,PIE2SI,KMPSI,O2CRITSI,JSIDETR
        WRITE(BFO,2220)  DPIE1PO4F,DPIE1PO4S,PIE2PO4,O2CRIT,KMO2DP
        WRITE(BFO,2230)  TEMPBEN,KBENSTR,KLBNTH,DPMIN
        WRITE(BFO,2240)  KAPPCH4,THTACH4
        WRITE(BFO,2130) (BBN(BB),WSSNET(BB),WSLNET(BB),WSRNET(BB),
     .                   WSCNET(BB),WSDNET(BB),WSGNET(BB),BB=1,NBB)
        WRITE(BFO,2140) (BBN(BB),VSED(BB),BB=1,NBB)
        WRITE(BFO,2150) (BBN(BB),VPMIX(BB),BB=1,NBB)
        WRITE(BFO,2160) (BBN(BB),VDMIX(BB),BB=1,NBB)
        CLOSE(BFO)
      END IF

********************************************************************************
**                             Initializations                                **
********************************************************************************

***** Convert cell heights and burial velocities to sediment units

      DIFFT = 0.0001*DIFFT
      DO 10000 BB=1,NBB
        HSED(BB) = HSEDALL*0.01
        VSED(BB) = VSED(BB)*2.73791E-5
10000 CONTINUE

***** Set sediment concentrations to initial concentrations

      DO 10010 BB=1,NBB
        POP1TM1S(BB) = CPOP(BB,1)
        POP2TM1S(BB) = CPOP(BB,2)
        POP3TM1S(BB) = CPOP(BB,3)
        PON1TM1S(BB) = CPON(BB,1)
        PON2TM1S(BB) = CPON(BB,2)
        PON3TM1S(BB) = CPON(BB,3)
        POC1TM1S(BB) = CPOC(BB,1)
        POC2TM1S(BB) = CPOC(BB,2)
        POC3TM1S(BB) = CPOC(BB,3)
        PSITM1S(BB)  = CPOS(BB)
10010 CONTINUE

***** Initialize mass balance variables

      DO 10015 BB=1,NBB
        CPIP(BB) = PO4T2TM1S(BB)
        CNO3(BB) = NO3T2TM1S(BB)
        CNH4(BB) = NH4T2TM1S(BB)
        ISEDMN   = ISEDMN+(CPON(BB,1)+CPON(BB,2)+CPON(BB,3)+CNH4(BB)
     .             +CNO3(BB))*A(VFN(1,BB))*HSED(BB)/1.E6
        ISEDMP   = ISEDMP+(CPOP(BB,1)+CPOP(BB,2)+CPOP(BB,3)+CPIP(BB))
     .             *A(VFN(1,BB))*HSED(BB)/1.E6
        ISEDMC   = ISEDMC+(CPOC(BB,1)+CPOC(BB,2)+CPOC(BB,3))
     .             *A(VFN(1,BB))*HSED(BB)/1.E6
10015 CONTINUE

***** Set up reaction rates in table look-up form

      DO 10020 JT=1,350
        TEMP         = FLOAT(JT-1)/10.+0.05
        TEMP20       = TEMP-20.
        TEMP202      = TEMP20/2.
        ZHTANH4F(JT) = KAPPNH4F*THTANH4**TEMP202
        ZHTANH4S(JT) = KAPPNH4S*THTANH4**TEMP202
        ZHTAD1(JT)   = KAPPD1*THTAPD1**TEMP202
        ZHTAP1(JT)   = KAPPP1*THTAPD1**TEMP202
        ZHTANO3F(JT) = KAPPNO3F*THTANO3**TEMP202
        ZHTANO3S(JT) = KAPPNO3S*THTANO3**TEMP202
        ZHTA2NO3(JT) = K2NO3*THTANO3**TEMP20
        ZL12NOM(JT)  = THTADD**TEMP20
        ZW12NOM(JT)  = THTADP**TEMP20
        ZHTAPON1(JT) = KPON1*THTAPON1**TEMP20
        ZHTAPON2(JT) = KPON2*THTAPON2**TEMP20
        ZHTAPON3(JT) = KPON3*THTAPON3**TEMP20
        ZHTAPOC1(JT) = KPOC1*THTAPOC1**TEMP20
        ZHTAPOC2(JT) = KPOC2*THTAPOC2**TEMP20
        ZHTAPOC3(JT) = KPOC3*THTAPOC3**TEMP20
        ZHTAPOP1(JT) = KPOP1*THTAPOP1**TEMP20
        ZHTAPOP2(JT) = KPOP2*THTAPOP2**TEMP20
        ZHTAPOP3(JT) = KPOP3*THTAPOP3**TEMP20
        ZHTASI(JT)   = KSI*THTASI**TEMP20
        ZHTACH4(JT)  = KAPPCH4*THTACH4**TEMP202
10020 CONTINUE

***** Turn off settling

      IF (.NOT.SETTLING) THEN
        DO 10030 BB=1,NBB
          WSSNET(BB) = 0.
          WSLNET(BB) = 0.
          WSRNET(BB) = 0.
          WSCNET(BB) = 0.
          WSDNET(BB) = 0.
          WSGNET(BB) = 0.
10030   CONTINUE
      END IF

***** Initialize accumulators for steady-state computations

      IF (STEADY_STATE_SED) THEN
        TINTIM = 0.
        DO 10035 BB=1,NBB
          AG3CFL(BB) = 0.
          AG3NFL(BB) = 0.
          AG3PFL(BB) = 0.
          ASDTMP(BB) = 0.
10035   CONTINUE
      END IF

***** Input FORMAT'S

 1000 FORMAT(:////8X,F8.0,I8)
 1022 FORMAT(///(8X,A8)/)
 1010 FORMAT(8F10.0)
C1020 FORMAT(:////(8X,9F8.1))
 1020 FORMAT(:/(8X,9F8.1))
 1030 FORMAT(I10)
 1040 FORMAT(:////(8X,6F8.1))
 
***** Output FORMAT'S

 2000 FORMAT(///34X,'Sediment-water column linkages and sediment ',
     .       'depths and volumes'/)
 2020 FORMAT(/' ACTIVE LAYER DEPTH ',F8.3,' CM')
 2022 FORMAT(/' STEADY-STATE VALUES OF G3 COMPUTED'/)
 2025 FORMAT(/' NO STEADY-STATE VALUES OF G3 COMPUTED'/)
 2030 FORMAT(////33X,'S E D I M E N T   I N I T I A L   C O N D I T ',
     .       'I O N S'/)
 2040 FORMAT(//25X,'Sediment initial conditions for ',A20/)
 2050 FORMAT(13X,3(7X,1PE11.4))
 2060 FORMAT(//25X,'Sediment initial conditions for ',A20/
     .       37X,'G1',22X,'G2',22X,'G3'/)
 2070 FORMAT(18X,3(2X,1PE11.4))
 2080 FORMAT(//30X,'Temperature diffusion coefficient ',E10.3,
     .       ' cm**2/sec')
 2090 FORMAT(//31X,'If salinity < ',F10.3,' ppt, methane formed',/
     .      30X,'If salinity < ',F10.3,' ppt, high nit/denit used')
 2100 FORMAT(//30X,'Particulate organic matter G-model splits'/
     .       10X,'fraction of....',5X,'recycled to',5X,'G1',5X,'G2',
     .       5X,'G3')
 2110 FORMAT(6X,A24,11X,3F7.2)
 2120 FORMAT(//30X,'Diagenesis rates (/day)  Temp corr factor'/
     .       30X,'Phosphorus'/
     .       39X,'G1',E11.3,5X,F7.3/
     .       39X,'G2',E11.3,5X,F7.3/
     .       39X,'G3',E11.3,5X,F7.3/
     .       30X,'Nitrogen'/
     .       39X,'G1',E11.3,5X,F7.3/
     .       39X,'G2',E11.3,5X,F7.3/
     .       39X,'G3',E11.3,5X,F7.3/
     .       30X,'Carbon'/
     .       39X,'G1',E11.3,5X,F7.3/
     .       39X,'G2',E11.3,5X,F7.3/
     .       39X,'G3',E11.3,5X,F7.3/
     .       30X,'Silica'/
     .       41X,E11.3,5X,F7.3)
 2130 FORMAT(//3X,'BBN',6X,'WSSNET',4X,'WSLNET',4X,'WSRNET',4X,
     .       'WSCNET',4X,'WSDNET',4X,'WSGNET'/
     .       (I7,6F10.3))
 2140 FORMAT(//31X,'Sedimentation rates (cm/yr)'/
     .       10X,8(I5,F6.2))
 2150 FORMAT(//30X,'Sediment solid-phase mixing rates (m**2/day)'/
     .       10X,8(I5,F6.2))
 2160 FORMAT(//30X,'Sediment dissolved-phase mixing rates (m**2/day)'/
     .       10X,8(I5,F6.2))
 2170 FORMAT(//35X,'Additional constants'/
     .       30X,'M1........',F8.2,' kg/l'/
     .       30X,'M2........',F8.2,' kg/l'/
     .       30X,'THTADP....',F8.3,/
     .       30X,'THTADD....',F8.3)
 2180 FORMAT(30X,'KAPPNH4F..',F8.3,' m/day'/
     .       30X,'KAPPNH4S..',F8.3,' m/day'/
     .       30X,'PIENH4....',F8.3,' l/kg'/
     .       30X,'THTANH4...',F8.3,/
     .       30X,'KMNH4.....',F8.3,' mg n/m**3'/
     .       30X,'KMNH4O2...',F8.3,' mg o2/l')
 2190 FORMAT(30X,'KAPPNO3F..',F8.3,' m/day'/
     .       30X,'KAPPNO3S..',F8.3,' m/day'/
     .       30X,'K2NO3.....',F8.3,' /day'/
     .       30X,'THTANO3...',F8.3)
 2200 FORMAT(30X,'KAPPD1....',F8.3,' m/day'/
     .       30X,'KAPPP1....',F8.3,' m/day'/
     .       30X,'PIE1S.....',F8.3,' l/kg'/
     .       30X,'PIE2S.....',F8.3,' l/kg'/
     .       30X,'THTAPD1...',F8.3,/
     .       30X,'KMHSO2....',F8.3,' mg o2/l')
 2210 FORMAT(30X,'CSISAT....',F8.1,' mg si/m**3'/
     .       30X,'DPIE1SI....',F8.3,' l/kg'/
     .       30X,'PIE2SI....',F8.3,' l/kg'/
     .       30X,'KMPSI.....',E8.2,' mg si/m**3'/
     .       30X,'O2CRITSI..',E8.2,' mg O2/L'/
     .       30X,'JSIDETR...',E8.2,' mg Si/m2-d')
 2220 FORMAT(30X,'DPIE1PO4F..',F8.3,' l/kg'/
     .       30X,'DPIE1PO4S..',F8.3,' l/kg'/
     .       30X,'PIE2PO4..',F8.3,' l/kg'/
     .       30X,'O2CRIT....',F8.3,' mg o2/l'/
     .       30X,'KMO2DP....',F8.3,' mg o2/l')
 2230 FORMAT(30X,'TEMPBEN...',F8.3,' deg c'/
     .       30X,'KBENSTR...',F8.3,' /day'/
     .       30X,'KLBNTH....',F8.3,'---'/
     .       30X,'DPMIN.....',F8.3,' m2/d')
 2240 FORMAT(30X,'KAPPCH4...',F8.3,' m/day'/
     .       30X,'THTACH4...',F8.3)
 2250 FORMAT(//31X,'G2 - G3 splits for Refractory Particulates'/
     .       '   NBB   FRG2P   FRG3P   FRG2N   FRG3N   FRG2C   FRG3C'/
     .       (I5,6F8.3))

***** Error output FORMAT's

 3000 FORMAT(///5X,'Zbrent failure, ierr=',I2,' in segment 'I5/
     .       5X,'Program termination follows diagnostic dumps')
 3010 FORMAT(/' Read error in sediment input deck')
      RETURN

********************************************************************************
**                           Sediment Calculations                            **
********************************************************************************

      ENTRY SED_CALC

***** Pass WQM time-step (in days) to sediment subr

      DLTS = DLT/86400.
      IF (STEADY_STATE_SED) TINTIM = TINTIM+DLTS

***** Initialize sediment nutrient masses

      SEDMN = 0.
      SEDMP = 0.
      SEDMC = 0.

***** Calculate fluxes

      DO 10060 BB=1,NBB
        IWC = BBN(BB)

***** Flux rate

        FLX1WC    = 1000.*WSSNET(BB)
        FLX3WC    = 1000.*WSRNET(BB)
        FLX4WC    = 1000.*WSCNET(BB)
        FLX5WC    = 1000.*WSDNET(BB)
        FLX6WC    = 1000.*WSGNET(BB)

***** Fluxes

        FLXPOP(BB,1) = FLX4WC*APC(IWC)*FRPPH1(1)*BC(IWC)
     .                +FLX5WC*APC(IWC)*FRPPH2(1)*BD(IWC)
     .                +FLX6WC*APC(IWC)*FRPPH3(1)*BG(IWC)
     .                +1000.*WSLNET(BB)*LPOP(IWC)
        FLXPOP(BB,2) = FLX4WC*APC(IWC)*FRPPH1(2)*BC(IWC)
     .                +FLX5WC*APC(IWC)*FRPPH2(2)*BD(IWC)
     .                +FLX6WC*APC(IWC)*FRPPH3(2)*BG(IWC)
     .                +FLX3WC*RPOP(IWC)*FRPOP(BB,2)/
     .                 (FRPOP(BB,2)+FRPOP(BB,3))
        FLXPOP(BB,3) = FLX4WC*APC(IWC)*FRPPH1(3)*BC(IWC)
     .                +FLX5WC*APC(IWC)*FRPPH2(3)*BD(IWC)
     .                +FLX6WC*APC(IWC)*FRPPH3(3)*BG(IWC)
     .                +FLX3WC*RPOP(IWC)*FRPOP(BB,3)/
     .                 (FRPOP(BB,2)+FRPOP(BB,3))
        FLXPON(BB,1) = FLX4WC*ANCC*FRNPH1(1)*BC(IWC)
     .                +FLX5WC*ANCD*FRNPH2(1)*BD(IWC)
     .                +FLX6WC*ANCG*FRNPH3(1)*BG(IWC)
     .                +1000.*WSLNET(BB)*LPON(IWC)
        FLXPON(BB,2) = FLX4WC*ANCC*FRNPH1(2)*BC(IWC)
     .                +FLX5WC*ANCD*FRNPH2(2)*BD(IWC)
     .                +FLX6WC*ANCG*FRNPH3(2)*BG(IWC)
     .                +FLX3WC*RPON(IWC)*FRPON(BB,2)/
     .                 (FRPON(BB,2)+FRPON(BB,3))
        FLXPON(BB,3) = FLX4WC*ANCC*FRNPH1(3)*BC(IWC)
     .                +FLX5WC*ANCD*FRNPH2(3)*BD(IWC)
     .                +FLX6WC*ANCG*FRNPH3(3)*BG(IWC)
     .                +FLX3WC*RPON(IWC)*FRPON(BB,3)/
     .                 (FRPON(BB,2)+FRPON(BB,3))
        FLXPOC(BB,1) = FLX4WC*FRCPH1(1)*BC(IWC)
     .                +FLX5WC*FRCPH2(1)*BD(IWC)
     .                +FLX6WC*FRCPH3(1)*BG(IWC)
     .                +1000.*WSLNET(BB)*LPOC(IWC)
        FLXPOC(BB,2) = FLX4WC*FRCPH1(2)*BC(IWC)
     .                +FLX5WC*FRCPH2(2)*BD(IWC)
     .                +FLX6WC*FRCPH3(2)*BG(IWC)
     .                +FLX3WC*RPOC(IWC)*FRPOC(BB,2)/
     .                 (FRPOC(BB,2)+FRPOC(BB,3))
        FLXPOC(BB,3) = FLX4WC*FRCPH1(3)*BC(IWC)+
     .                 FLX5WC*FRCPH2(3)*BD(IWC)
     .                +FLX6WC*FRCPH3(3)*BG(IWC)
     .                +FLX3WC*RPOC(IWC)*FRPOC(BB,3)/
     .                 (FRPOC(BB,2)+FRPOC(BB,3))

        FLXPOS(BB) = FLX5WC*ASCD*BD(IWC)+FLX1WC*SU(IWC)

C       *** Add adsorbed phosphate and silica to inorganic pool

        PF           = KADPO4*SSI(IWC)/(1.+KADPO4*SSI(IWC))
        PIP          = PF*(PO4(IWC)-APC(IWC)*(BC(IWC)+BD(IWC)+BG(IWC)))
        PO4T2TM1S(BB)= PO4T2TM1S(BB)+FLX1WC*PIP*DLTS/HSED(BB)
        PF           = KADSA*SSI(IWC)/(1.+KADSA*SSI(IWC))
        SIT2TM1S(BB) = SIT2TM1S(BB)+FLX1WC*PF*SA(IWC)*DLTS/HSED(BB)

C       *** Sum particulate fluxes to sediments, negative into sediments

        PPFWS(BB)=-0.001*(FLXPOP(BB,1)+FLXPOP(BB,2)+FLXPOP(BB,3))
        PNFWS(BB)=-0.001*(FLXPON(BB,1)+FLXPON(BB,2)+FLXPON(BB,3))
        PCFWS(BB)=-0.001*(FLXPOC(BB,1)+FLXPOC(BB,2)+FLXPOC(BB,3))
        PSFWS(BB)=-0.001*FLXPOS(BB)
        SSFWS(BB)=-WSSNET(BB)*SSI(IWC)

10060 CONTINUE

C     *** Accumulate fluxes for steady-state computation

      IF (STEADY_STATE_SED) THEN
        DO 10050 BB=1,NBB
          AG3CFL(BB) = AG3CFL(BB)+FLXPOC(BB,3)*DLTS
          AG3NFL(BB) = AG3NFL(BB)+FLXPON(BB,3)*DLTS
          AG3PFL(BB) = AG3PFL(BB)+FLXPOP(BB,3)*DLTS
10050   CONTINUE
      ENDIF

***** Assign previous timestep concentrations to particulate organics

      DO 10070 BB=1,NBB
        CPOP(BB,1) = POP1TM1S(BB)
        CPOP(BB,2) = POP2TM1S(BB)
        CPOP(BB,3) = POP3TM1S(BB)
        CPON(BB,1) = PON1TM1S(BB)
        CPON(BB,2) = PON2TM1S(BB)
        CPON(BB,3) = PON3TM1S(BB)
        CPOC(BB,1) = POC1TM1S(BB)
        CPOC(BB,2) = POC2TM1S(BB)
        CPOC(BB,3) = POC3TM1S(BB)
        CPOS(BB)   = PSITM1S(BB)
10070 CONTINUE

***** Update sediment concentrations

      DO 10080 BB=1,NBB

******* Assign previous timestep concentrations

        NH41TM1  = NH41TM1S(BB)
        NO31TM1  = NO31TM1S(BB)
        HS1TM1   = HS1TM1S(BB)
        SI1TM1   = SI1TM1S(BB)
        PO41TM1  = PO41TM1S(BB)
        BENSTR1  = BENSTR1S(BB)
        NH4T2TM1 = NH4T2TM1S(BB)
        NO3T2TM1 = NO3T2TM1S(BB)
        HST2TM1  = HST2TM1S(BB)
        SIT2TM1  = SIT2TM1S(BB)
        PO4T2TM1 = PO4T2TM1S(BB)
        PON1TM1  = PON1TM1S(BB)
        PON2TM1  = PON2TM1S(BB)
        PON3TM1  = PON3TM1S(BB)
        POC1TM1  = POC1TM1S(BB)
        POC1     = POC1TM1
        POC2TM1  = POC2TM1S(BB)
        POC3TM1  = POC3TM1S(BB)
        POP1TM1  = POP1TM1S(BB)
        POP2TM1  = POP2TM1S(BB)
        POP3TM1  = POP3TM1S(BB)
        PSITM1   = PSITM1S(BB)
        BFORMAX  = BFORMAXS(BB)
        ISWBEN   = ISWBENS(BB)
        H2       = HSED(BB)

******* Sedimentation, mixing rates, and sediment temperature

        W2    = VSED(BB)
        DP    = VPMIX(BB)
        DD    = VDMIX(BB)
        TEMPD = CTEMP(BB)
        STP20 = TEMPD-20.

******* Convert overlying water column concentrations into mg/m**3

        IWC  = BBN(BB)
        DF   = 1./(1.+KADPO4*SSI(IWC))
        PO4AVL = DF*(PO4(IWC)-APC(IWC)*(BC(IWC)+BD(IWC)+BG(IWC)))
        PO40 = PO4AVL*1000.
        NH40 = NH4(IWC)*1000.
        NO30 = NO3(IWC)*1000.
        DF   = 1./(1.+KADSA*SSI(IWC))
        SI0  = DF*SA(IWC)*1000.
        O20  = AMAX1(DO(IWC),0.001)
        SAL  = SALT(IWC)

******* Methane saturation

        CH4SAT = 0.099*(1.+(ZD(IWC)+BL(IWC,3)+HSED(BB))/10.)
     .           *0.9759**STP20

******* Evaluate the temperature dependent coefficients

        ITEMP    = 10.*TEMPD+1

******* Salinity dependence of nitrification and denitrification

        IF (SAL.LE.SALTND) THEN
          XAPPNH4  = ZHTANH4F(ITEMP)
          XAPP1NO3 = ZHTANO3F(ITEMP)
        ELSE
          XAPPNH4  = ZHTANH4S(ITEMP)
          XAPP1NO3 = ZHTANO3S(ITEMP)
        END IF
        XAPPD1   = ZHTAD1(ITEMP)
        XAPPP1   = ZHTAP1(ITEMP)
        XK2NO3   = ZHTA2NO3(ITEMP)*H2
        XKSI     = ZHTASI(ITEMP)*H2
        XAPPCH4  = ZHTACH4(ITEMP)
        KL12NOM  = DD/H2*ZL12NOM(ITEMP)
        W12NOM   = DP/H2*ZW12NOM(ITEMP)*POC1/1.0E5
        IF (ISWBEN.EQ.0) THEN
          IF (TEMPD.GE.TEMPBEN) THEN
            ISWBEN  = 1
            BFORMAX = 0.
          ENDIF
          BFOR = KMO2DP/(KMO2DP+O20)
        ELSE
          IF (TEMPD.LT.TEMPBEN) THEN
            ISWBEN = 0
          ENDIF
          BFORMAX = AMAX1(KMO2DP/(KMO2DP+O20),BFORMAX)
          BFOR    = BFORMAX
        ENDIF
        BENSTR = (BENSTR1+DLTS*BFOR)/(1.+KBENSTR*DLTS)
c## -- add minimum mixing term and bio-irrigation formulation
c##
c##     W12    = W12NOM*(1.-KBENSTR*BENSTR)
c##     KL12   = KL12NOM
c## -- w12min= Dpmin/h2 is minimum particle mixing
        W12MIN = DPMIN/H2
        W12    = W12NOM*(1.-KBENSTR*BENSTR)+W12MIN
c## -- klbnth is ratio of bio-irrigation to bio-particle mixing
        KL12   = KL12NOM + KLBNTH*W12

******* Lookup reaction rates

        ITEMP  = 10.*TEMPD+1
        XKPOP1 = ZHTAPOP1(ITEMP)*H2
        XKPOP2 = ZHTAPOP2(ITEMP)*H2
        XKPOP3 = ZHTAPOP3(ITEMP)*H2
        XKPON1 = ZHTAPON1(ITEMP)*H2
        XKPON2 = ZHTAPON2(ITEMP)*H2
        XKPON3 = ZHTAPON3(ITEMP)*H2
        XKPOC1 = ZHTAPOC1(ITEMP)*H2
        XKPOC2 = ZHTAPOC2(ITEMP)*H2
        XKPOC3 = ZHTAPOC3(ITEMP)*H2
        XKSI   = ZHTASI(ITEMP)*H2

******* Calculate sediment concentrations

        DOH2 = DLTS/H2
        FD2  = 1./(1.+M2*PIE2SI)
        K3   = XKSI*(CSISAT-FD2*SIT2TM1)/(PSITM1+KMPSI)
        PON1 = (FLXPON(BB,1)*DOH2+PON1TM1)/(1.+(XKPON1+W2)*DOH2)
        PON2 = (FLXPON(BB,2)*DOH2+PON2TM1)/(1.+(XKPON2+W2)*DOH2)
        PON3 = (FLXPON(BB,3)*DOH2+PON3TM1)/(1.+(XKPON3+W2)*DOH2)
        POC1 = (FLXPOC(BB,1)*DOH2+POC1TM1)/(1.+(XKPOC1+W2)*DOH2)
        POC2 = (FLXPOC(BB,2)*DOH2+POC2TM1)/(1.+(XKPOC2+W2)*DOH2)
        POC3 = (FLXPOC(BB,3)*DOH2+POC3TM1)/(1.+(XKPOC3+W2)*DOH2)
        POP1 = (FLXPOP(BB,1)*DOH2+POP1TM1)/(1.+(XKPOP1+W2)*DOH2)
        POP2 = (FLXPOP(BB,2)*DOH2+POP2TM1)/(1.+(XKPOP2+W2)*DOH2)
        POP3 = (FLXPOP(BB,3)*DOH2+POP3TM1)/(1.+(XKPOP3+W2)*DOH2)
c## -- modification for detrital Si input to sediment
c##     PSI  = (FLXPOS(BB)*DLTS/H2+PSITM1)/(1.+(K3+W2)*DLTS/H2)
        PSI  = ((FLXPOS(BB)+JSIDETR)*DOH2+PSITM1)/(1.+(K3+W2)*DOH2)

******* Assign diagenesis values for sediment model

        XJN = XKPON1*PON1+XKPON2*PON2+XKPON3*PON3
        XJC = XKPOC1*POC1+XKPOC2*POC2+XKPOC3*POC3
        XJP = XKPOP1*POP1+XKPOP2*POP2+XKPOP3*POP3

C TEMPORARY BYPASS OF FLUX ALGORITHMS
C        GO TO 66666

******* Evaluate the NH4, NO3, and SOD equations

        SOD = ZBRENT(IERR)
        IF (IERR.NE.0.AND.BENTHIC_OUTPUT) WRITE(BFO,3000) IERR,BB

******* Accumulate remaining sums for steady-state computation

        IF (STEADY_STATE_SED) THEN
          ASDTMP(BB) = ASDTMP(BB)+TEMPD*DLTS
        END IF

******* Evaluate the PO4 and Si equations

        K0H1D = 0.
        K0H1P = 0.
        KMC1  = 0.
        K1H1D = S
        K1H1P = 0.
        K2H2D = 0.
        K2H2P = 0.
        J1    = S*SI0

******* Oxygen dependency of pie1

        IF (O20.LT.O2CRITSI) THEN
          PIE1 = PIE2SI*DPIE1SI**(O20/O2CRITSI)
        ELSE
          PIE1 = PIE2SI*DPIE1SI
        ENDIF
        PIE2 = PIE2SI

******* Silica dissolution kinetics

        FD2 = 1./(1.+M2*PIE2)
        K3  = XKSI*PSI/(PSI+KMPSI)*FD2
        J2  = XKSI*PSI/(PSI+KMPSI)*CSISAT
        CALL SEDTSFNL (SI1,SI2,SIT1,SIT2,SI1TM1,SIT2TM1)
        JSI = S*(SI1-SI0)

******* Phosphate

        K0H1D = 0.
        K0H1P = 0.
        KMC1  = 0.
        K1H1D = S
        K1H1P = 0.
        K2H2D = 0.
        K2H2P = 0.
        J1    = S*PO40
        K3    = 0.
        J2    = XJP

******* Salinity dependence of PIE1

        IF (SAL.LE.SALTSW) THEN
          DPIE1PO4=DPIE1PO4F
        ELSE
          DPIE1PO4=DPIE1PO4S
        END IF

******* Oxygen dependency of pie1

        IF (O20.LT.O2CRIT) THEN
          PIE1 = PIE2PO4*DPIE1PO4**(O20/O2CRIT)
        ELSE
          PIE1 = PIE2PO4*DPIE1PO4
        ENDIF
        PIE2 = PIE2PO4
        CALL SEDTSFNL (PO41,PO42,PO4T1,PO4T2,PO41TM1,PO4T2TM1)
        JPO4 = S*(PO41-PO40)

******* Replace the t minus 1 concentrations

        NH41TM1S(BB)  = NH41
        NO31TM1S(BB)  = NO31
        HS1TM1S(BB)   = HS1
        SI1TM1S(BB)   = SI1
        PO41TM1S(BB)  = PO41
        BENSTR1S(BB)  = BENSTR
        NH4T2TM1S(BB) = NH4T2
        NO3T2TM1S(BB) = NO3T2
        HST2TM1S(BB)  = HST2
        SIT2TM1S(BB)  = SIT2
        PO4T2TM1S(BB) = PO4T2
        PON1TM1S(BB)  = PON1
        PON2TM1S(BB)  = PON2
        PON3TM1S(BB)  = PON3
        POC1TM1S(BB)  = POC1
        POC2TM1S(BB)  = POC2
        POC3TM1S(BB)  = POC3
        POP1TM1S(BB)  = POP1
        POP2TM1S(BB)  = POP2
        POP3TM1S(BB)  = POP3
        PSITM1S(BB)   = PSI
        BFORMAXS(BB)  = BFORMAX
        ISWBENS(BB)   = ISWBEN

******* Assign flux-flux results to wqm arrays

        ITEMP      = 10*TEMPD+1
        IF (SAL.LE.SALTND) THEN
          XAPP1NO3 = ZHTANO3F(ITEMP)
        ELSE
          XAPP1NO3 = ZHTANO3S(ITEMP)
        END IF
        XK2NO3     = ZHTA2NO3(ITEMP)*H2
        BENDO(BB)  = -SOD
        MTVEL(BB)  = SOD/O20
        BENNH4(BB) = JNH4/1000.
        BENNO3(BB) = JNO3/1000.
        BENPO4(BB) = JPO4/1000.
        BENDOC(BB) = 0.
C       BENCOD(BB) = JHS+JCH4AQ+JCH4G
        BENCOD(BB) = JHS+JCH4AQ 
        BENSA(BB)  = JSI/1000.
        BENDEN(BB) = (XAPP1NO3*XAPP1NO3*NO31/S+XK2NO3*NO32)/1000.

******* Fluxes due to burial of particulates

        BURIALN(BB) = (PON1+PON2+PON3+NO3T2+NH4T2)*W2
        BURIALP(BB) = (POP1+POP2+POP3+PO4T2)*W2
        BURIALC(BB) = (POC1+POC2+POC3)*W2

******* Diagenesis of carbon forms

        DIAGENC(BB) = XJC/1000.

C END OF BYPASS OF FLUX ALGORITHMS
66666   CONTINUE

******* Total sediment nutrient mass

        SEDMN = SEDMN+(PON1+PON2+PON3+NH4T2+NO3T2)*A(VFN(1,BB))*H2/1.E6
        SEDMP = SEDMP+(POP1+POP2+POP3+PO4T2)*A(VFN(1,BB))*H2/1.E6
        SEDMC = SEDMC+(POC1+POC2+POC3)*A(VFN(1,BB))*H2/1.E6
10080 CONTINUE

***** Assign concentrations to plot variables

      DO 10085 BB=1,NBB
        CPON(BB,1) = PON1TM1S(BB)
        CPON(BB,2) = PON2TM1S(BB)
        CPON(BB,3) = PON3TM1S(BB)
        CNH4(BB)   = NH4T2TM1S(BB)
        CNO3(BB)   = NO3T2TM1S(BB)
        CPOP(BB,1) = POP1TM1S(BB)
        CPOP(BB,2) = POP2TM1S(BB)
        CPOP(BB,3) = POP3TM1S(BB)
        CPIP(BB)   = PO4T2TM1S(BB)
        CPOC(BB,1) = POC1TM1S(BB)
        CPOC(BB,2) = POC2TM1S(BB)
        CPOC(BB,3) = POC3TM1S(BB)
        CPOS(BB)   = PSITM1S(BB)
10085 CONTINUE

***** Take temperature integration step

      DO 10090 BB=1,NBB
        CTEMP(BB) = CTEMP(BB)+DLT*DIFFT/HSED(BB)/HSED(BB)
     .              *(T(BBN(BB))-CTEMP(BB))
10090 CONTINUE
      RETURN

***** Compute and print out steady-state sediment concentrations

      ENTRY SED_INT

***** Compute time-average values

      DO 20000 BB=1,NBB
        AG3CFL(BB) = AG3CFL(BB)/TINTIM
        AG3NFL(BB) = AG3NFL(BB)/TINTIM
        AG3PFL(BB) = AG3PFL(BB)/TINTIM
        ASDTMP(BB) = ASDTMP(BB)/TINTIM
20000 CONTINUE

***** Compute G3 organic concentrations

      DO 20010 BB=1,NBB
        CPOC(BB,3) = AG3CFL(BB)/(KCDIAG(3)*DCTHTA(3)**(ASDTMP(BB)-20.)
     .               *HSED(BB)+VSED(BB))
        CPON(BB,3) = AG3NFL(BB)/(KNDIAG(3)*DNTHTA(3)**(ASDTMP(BB)-20.)
     .               *HSED(BB)+VSED(BB))
        CPOP(BB,3) = AG3PFL(BB)/(KPDIAG(3)*DPTHTA(3)**(ASDTMP(BB)-20.)
     .               *HSED(BB)+VSED(BB))
20010 CONTINUE


      RETURN

***** Error traps

10100 CONTINUE
      IF (BENTHIC_OUTPUT) WRITE(BFO,3010)
      STOP
      END

********************************************************************************
**                          F U N C T I O N   S E D F                         **
********************************************************************************

      FUNCTION SEDF(SOD1)
      SAVE
      INCLUDE 'wqm_com.inc'


***** Compute the NH4, NO3, and SOD fluxes

      S = SOD1/O20

***** Ammonia flux

      K0H1P = 0.
      K1H1P = 0.
      K2H2D = 0.
      K2H2P = 0.
      IF (KMNH4.NE.0.) THEN
        K0H1D = XAPPNH4**2/S*KMNH4*(O20/(KMNH4O2+O20))
        K1H1D = S
      ELSE
        K1H1D = XAPPNH4**2/S*(O20/(KMNH4O2+O20))+S
        K0H1D = 0.
      ENDIF
      J1   = S*NH40
      K3   = 0.
      J2   = XJN
      PIE1 = PIENH4
      PIE2 = PIENH4
      KMC1 = KMNH4
      CALL SEDTSFNL (NH41,NH42,NH4T1,NH4T2,NH41TM1,NH4T2TM1)
      JNH4 = S*(NH41-NH40)

***** Oxygen consumed by nitrification

      A1 = 0.0045714
      IF (KMNH4.NE.0.) THEN
        JO2NH4 = A1*K0H1D*NH41/(KMNH4+NH41TM1)
      ELSE
        JO2NH4 = A1*(K1H1D-S)*NH41
      ENDIF

***** Denitrification

      K0H1D = 0.
      K0H1P = 0.
      KMC1  = 0.
      K1H1D = XAPP1NO3**2/S+S
      K1H1P = 0.
      K2H2D = XK2NO3
      K2H2P = 0.
      IF (KMNH4.NE.0.) THEN
        J1 = S*NO30+XAPPNH4**2/S*KMNH4*(O20/(KMNH4O2+O20))*NH41
     .       /(KMNH4+NH41TM1)
      ELSE
        J1 = S*NO30+XAPPNH4**2/S*(O20/(KMNH4O2+O20))*NH41
      ENDIF
      K3   = 0.
      J2   = 0.
      PIE1 = 0.
      PIE2 = 0.
      CALL SEDTSFNL(NO31,NO32,NO3T1,NO3T2,NO31TM1,NO3T2TM1)
      JNO3 = S*(NO31-NO30)

***** Sulfide/methane oxidation

      A2      = 0.00285714
      XJCNO31 = A2*XAPP1NO3**2/S*NO31
      XJCNO3  = A2*XK2NO3*NO32

***** Add the aerobic and first anaerobic layer to keep mass balance

      XJCNO3 = XJCNO31+XJCNO3

***** Convert carbon diagenesis flux to O2 units

      XJC1 = AMAX1(2.666666666E-3*XJC-XJCNO3,0.)
         
***** Sulfide or methane in O2 equivalents

      IF (SAL.GT.SALTSW) THEN

******* Sulfide

        K0H1D = 0.
        K0H1P = 0.
        KMC1  = 0.
        K1H1D = XAPPD1**2/S*(O20/KMHSO2)+S
        K1H1P = XAPPP1**2/S*(O20/KMHSO2)
        K2H2D = 0.
        K2H2P = 0.
        J1    = 0.
        K3    = 0.
        J2    = XJC1
        PIE1  = PIE1S
        PIE2  = PIE2S
        CALL SEDTSFNL(HS1,HS2,HST1,HST2,HS1TM1,HST2TM1)
        JHS    = S*HS1
        CSOD   = (XAPPD1**2/S*FD1+XAPPP1**2/S*FP1)*(O20/KMHSO2)*HST1
        JCH4AQ = 0.
        JCH4G  = 0.
      ELSE

******* Methane

        CSODMX = SQRT(2.*KL12*CH4SAT*XJC1)
        IF (CSODMX.GT.XJC1) CSODMX = XJC1
        IF ((XAPPCH4/S).LT.80) THEN
          SECHXC = 2./(EXP(XAPPCH4/S)+EXP(-XAPPCH4/S))
        ELSE
          SECHXC = 0.
        ENDIF
        CSOD   = CSODMX*(1.-SECHXC)
        JCH4AQ = CSODMX*SECHXC
        JCH4G  = XJC1-JCH4AQ-CSOD
        JHS    = 0.
      ENDIF

***** SOD function

      SOD  = CSOD+JO2NH4
      SEDF = SOD-SOD1
      RETURN
      END

********************************************************************************
**                        F U N C T I O N   Z B R E N T                       **
********************************************************************************

      FUNCTION ZBRENT(IERR)
      EXTERNAL SEDF

***** Parameter declarations

      PARAMETER (IMAX=100,EPS=3.E-8,TOL=1.E-5,SODMIN=1.E-4,SODMAX=100.)

***** Initialize upper and lower limits for solution

      IERR = 0
      A    = SODMIN
      B    = SODMAX
      FA   = SEDF(A)
      FB   = SEDF(B)

***** Root must bracket ZBRENT

      IF (FB*FA.GT.0.) THEN
        IERR = 1
        RETURN
      ENDIF
      FC = FB
      DO 10000 I=1,IMAX
        IF (FB*FC.GT.0.) THEN
          C  = A
          FC = FA
          D  = B-A
          E  = D
        ENDIF
        IF (ABS(FC).LT.ABS(FB)) THEN
          A  = B
          B  = C
          C  = A
          FA = FB
          FB = FC
          FC = FA
        ENDIF
        TOL1 = 2.*EPS*ABS(B)+0.5*TOL
        XM   = 0.5*(C-B)
        IF (ABS(XM).LE.TOL1.OR.FB.EQ.0.) THEN
          ZBRENT = B
          RETURN
        ENDIF
        IF (ABS(E).GE.TOL1.AND.ABS(FA).GT.ABS(FB)) THEN
          S = FB/FA
          IF (A.EQ.C) THEN
            P = 2.*XM*S
            Q = 1.-S
          ELSE
            Q = FA/FC
            R = FB/FC
            P = S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q = (Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF (P.GT.0.) Q = -Q
          P = ABS(P)
          IF (2.*P.LT.MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E = D
            D = P/Q
          ELSE
            D = XM
            E = D
          ENDIF
        ELSE
          D = XM
          E = D
        ENDIF
        A  = B
        FA = FB
        IF (ABS(D).GT.TOL1) THEN
          B = B+D
        ELSE
          B = B+SIGN(TOL1,XM)
        ENDIF
        FB = SEDF(B)
10000 CONTINUE
      IERR   = 2
      ZBRENT = B
      RETURN
      END

********************************************************************************
**                    S U B R O U T I N E   S E D T S F N L                   **
********************************************************************************

      SUBROUTINE SEDTSFNL(C1S,C2S,CT1S,CT2S,C1TM1S,CT2TM1S)
      SAVE
      INCLUDE 'wqm_com.inc'

***** Initialize constants

      FD1 = 1./(1.+M1*PIE1)
      FP1 = M1*PIE1/(1.+M1*PIE1)
      FD2 = 1./(1.+M2*PIE2)
      FP2 = M2*PIE2/(1.+M2*PIE2)
      F12 = W12*FP1+KL12*FD1
      F21 = W12*FP2+KL12*FD2

***** Evaluate the MM term at time level t-1

      IF (KMC1.NE.0.) THEN
        XK0 = (K0H1D*FD1+K0H1P*FP1)/(KMC1+C1TM1S)
      ELSE
        XK0 = 0.
      ENDIF
      XK1 = XK0+K1H1D*FD1+K1H1P*FP1
      XK2 = K2H2D*FD2+K2H2P*FP2
      A11 = -F12-XK1-W2
      A21 = F12+W2
      A12 = F21
      B1  = -J1
      A22 = -F21-XK2-W2-K3-H2/DLTS
      B2  = -J2-H2/DLTS*CT2TM1S

***** Solve the 2x2 set of linear equations

      DELTA = A11*A22-A12*A21
      IF (DELTA.EQ.0.) THEN
        PRINT *,'Twod is singular: A11,A12,A21,A22'
        PRINT *,A11,A12,A21,A22
        STOP
      ENDIF

***** Assign results

      CT1S = (B1*A22-B2*A12)/DELTA
      CT2S = (B2*A11-B1*A21)/DELTA
      C1S  = FD1*CT1S
      C2S  = FD2*CT2S
      RETURN
      END
