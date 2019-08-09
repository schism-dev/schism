********************************************************************************
**                                                                            **
**                        COMMON file for CE-QUAL-IC                          **
**                                                                            **
**                                Version 1.1 (really 2.0)                    **
**                               July 28, 1992                                **
**                                                                            **
**                        Water Quality Modeling Group                        **
**                        U.S. Army Corps of Engineers                        **
**                        Waterways Experiment Station                        **
**                        Vicksburg, Mississippi 39180                        **
**                                                                            **
********************************************************************************

********************************************************************************
**                                                                            **
**                           PARAMETER definitions                            **
**                                                                            **
**    NCP   - Number of constituent state variables                           **
**
**    NBP   - Number of boxes                                                 **
**    NQFP  - Number of horizontal and vertical flow faces                    **
**    NHQP  - Number of horizontal flow faces                                 **
**    NSBP  - Number of boxes in the surface layer                            **
**    NLP   - Maximum number of layers                                        **
**    NS1P  - Number of source one (fall line load) inputs                    **
**    NS2P  - Number of source two (below fall line load) inputs              **
**    NBCP  - Number of boundary concentration inputs                         **
**    NMP   - Number of modifications to initial box concentrations           **
**    NDP   - Maximum number of days for any output                           **
**    NSAVP -
**    NFLP  - Number of files for each type of time-varying data              **
**    NOIP  -
**                                                                            **
********************************************************************************


********************************************************************************
**                        Water Quality Model Setup                           **
********************************************************************************

***** Parameter declarations

      PARAMETER (NCP=22)

      PARAMETER (NBP=10739,NQFP=26661,NHQP=18638,NSBP=2716,NLP=17,     !uches_9
     .           NS1P=20,NS2P=0,NBCP=200,NMP=0,NDP=500,                 !uches_9
     .           NSAVP=50,NFLP=20,NOIP=10)                             !uches_9

C      PARAMETER (NBP=50,NQFP=60,NHQP=60,NSBP=20,NLP=4,NS1P=20,     !30 box
C     .           NS2P=20,NBCP=30,NMP=20,NDP=500,NSAVP=50,NFLP=20,  !30 box
C     .           NOIP=10)                                          !30 box

C      PARAMETER (NBP=4073,NQFP=9874,NHQP=6530,NSBP=729,NLP=15,NS1P=600,  !CHESAPEAKE
C     .           NS2P=50,NBCP=120,NMP=30,NDP=500,NSAVP=729,NFLP=100,     !CHESAPEAKE
C     .           NOIP=10)                                                !CHESAPEAKE

***** Data type declarations

      INTEGER B,      SBN,    BB,     BU,     BBN,    VFN,    S1LB,
     .        S1LN,   S2LB,   S2LN,   SAVB,   AC
      INTEGER DIA,    HYD,    CBC,    S1,     S2,     BFI,    BFO,
     .        ATM,    STL,    AGR,    SAV,    ERR
      INTEGER HYDPTR, METPTR, CBCPTR, S1PTR,  S2PTR,  BFIPTR, KEIPTR,
     .        ATMPTR, SAVPTR
      REAL    KHRC,   KHRD,   KHRG,   KHNC,   KHND,   KHNG,   KHPC,
     .        KHPD,   KHPG,   KHSD,   KHONT,  KHOCOD, KHNNT,  KHODOC,
     .        KHTIS,  S2L,    KHNDN
      REAL    KDC,    KDCALG, KLC,    KLCALG, KRC,    KRCALG, KSUA,   !MRN 2/93
     .        KCOD,   KDN,    KDNALG, KDP,    KDPALG, KLN,    KLNALG,
     .        KLP,    KLPALG, KRN,    KRNALG, KRP,    KRPALG, NTM, 
     .        KRDO
      REAL    KHNAVG, KHPAVG
      REAL    KSO,    KSNH4,  KSNO3,  KSPO4,  KSSA,   KSDOC,  KHSO,
     .        MTCNO3
      REAL    KTNT1,  KTGC1,  KTGD1,  KTGG1,  KTNT2,  KTGC2,  KTGD2,
     .        KTGG2,  KTBC,   KTBD,   KTBG,   KTCOD,  KTMNL,  KTHDR,
     .        KTSUA
      REAL    LPOC,   LPON,   LPOP,   MNLDOC, MTVEL
      REAL    NT,     KT,     NXCBC,  JDAY,   JDAYMBL
      REAL    NXMET,  NXS1,   NXS2,   NXBFI,  NXKEI,  NXATM,  NXSAV
      REAL    KE,     KECHL,  KESS
      REAL    I0,     I0NX,   IATTOP, IATBOT, IHC,    IHD,    IHG
      REAL    MBGM,   NPP,    NH4,    NO3,    NLC,
     .        NLD,    NLG
      REAL    KADPO4, KADSA
      REAL    MNLDON, MNLDOP, NITRIF

***** Change for application on Cray or other machines

      DOUBLE PRECISION V1,       V2,      ELTMS
      DOUBLE PRECISION MAXDLTDP, MXDLTDP, NXHYDDP, ELTMS1
C     REAL             MAXDLTDP, MXDLTDP, NXHYDDP

      LOGICAL BOUNDARY_CONC,  SOURCE_ONE,       SOURCE_TWO,
     .        BENTHIC_FLUXES, SEDIMENT_CALC,    ALGAL_GROWTH,
     .        ATMOS_LOADS,    TRANSPORT_FLUXES, SUB_AQ_VEG_CALC,
     .        AVERAGE_PLOTS,  QUALITY_DIAG,     KINETIC_FLUXES,
     .        SEDIMENT_DIAG,  DIAGNOSTICS,      CONSERVE_MASS,
     .        SETTLING,       STEP_BOUNDARY,    MASS_BALANCE,
     .        NO_SEDIMENT,    ERROR_OPEN,       ERROR_CHECK
      LOGICAL POSITIVE_FLOW,  FLOW,             XY_DIFFUSION,
     .        Z_DIFFUSION,    ZGRID_HYDRO,      ASCII_HYDRO,
     .        SIGMA_HYDRO
      CHARACTER*20 SSNAME
      CHARACTER*63 ERROR_MESS
      CHARACTER*72 BFOFN
      CHARACTER*72 METFN,  S1FN,   S2FN,  HYDFN,  CBCFN,  KEIFN,
     .             ATMFN,  SAVFN,  BFIFN, ERRFN
      CHARACTER*8  KEI_TYPE,   BFI_TYPE

***** Dimension declarations

      DIMENSION ERROR_CHECK(36)
      DIMENSION POSITIVE_FLOW(NQFP)
      DIMENSION NVF(NSBP),        VFN(0:NLP,NSBP)
      DIMENSION BFLUX(NSBP,9),    BFLUXB(NSBP,7)
      DIMENSION C1(0:NBP,NCP),    C2(0:NBP,NCP),    CSTAR(0:NBP,NCP),
     .          DTC(0:NBP,NCP),   AC1(0:NBP,NCP),   CMASS(NCP)
      DIMENSION S1L(NS1P,NCP),    S2L(NS2P,NCP),    S1LB(NS1P,NCP),
     .          S2LB(NS2P,NCP),   CB(NBCP,NCP),     CBNX(NBCP,NCP),
     .          SAVB(NSAVP),      SAVAREA(NSAVP),   CBOLD(NBCP,NCP)
      DIMENSION CPOP(NSBP,3),     FLXPOP(NSBP,3),   CPON(NSBP,3),
     .          FLXPON(NSBP,3),   CPOC(NSBP,3),     FLXPOC(NSBP,3),
     .          CPOS(NSBP),       CPIP(NSBP),       CNO3(NSBP),
     .          CNH4(NSBP),       FLXPOS(NSBP),     CTEMP(NSBP),
     .          CDTEMP(NSBP),     BSVOL(NSBP),      HSED(NSBP)
      DIMENSION ATMFLXNB(NSBP),   ATMFLXPB(NSBP),   ATMFLXCB(NSBP),
     .          BENFLXPNB(NSBP),  BENFLXDNB(NSBP),  BENFLXPPB(NSBP),
     .          BENFLXDPB(NSBP),  BENFLXPCB(NSBP),  DLSEDKNB(NSBP),
     .          DLSEDKCB(NSBP),
     .          BURIALFLXNB(NSBP),BURIALFLXPB(NSBP),BURIALFLXCB(NSBP)
      DIMENSION FLUXT(0:NQFP,NCP),FLUXS(NBP,9),     AFLUX(NQFP,11)
      DIMENSION FLXTTEM(NQFP),    FLXTSAL(NQFP),    FLXTSSI(NQFP),
     .          FLXTC(NQFP),      FLXTD(NQFP),      FLXTG(NQFP),
     .          FLXTDOC(NQFP),    FLXTLPOC(NQFP),   FLXTRPOC(NQFP),
     .          FLXTNH4(NQFP),    FLXTNO3(NQFP),    FLXTDON(NQFP),
     .          FLXTLPON(NQFP),   FLXTRPON(NQFP),   FLXTPO4(NQFP),
     .          FLXTDOP(NQFP),    FLXTLPOP(NQFP),   FLXTRPOP(NQFP),
     .          FLXTCOD(NQFP),    FLXTDO(NQFP),     FLXTSU(NQFP),
     .          FLXTSA(NQFP)
      DIMENSION FLXSSSI(NBP),     FLXSC(NBP),       FLXSD(NBP),
     .          FLXSG(NBP),       FLXSPOC(NBP),     FLXSPON(NBP),
     .          FLXSPO4(NBP),     FLXSPOP(NBP),     FLXSSI(NBP)
      DIMENSION AFLXPOC(NQFP),    AFLXDOC(NQFP),    AFLXTC(NQFP),
     .          AFLXDON(NQFP),    AFLXDIN(NQFP),    AFLXPON(NQFP),
     .          AFLXTN(NQFP),     AFLXDOP(NQFP),    AFLXDIP(NQFP),
     .          AFLXPOP(NQFP),    AFLXTP(NQFP)
      DIMENSION ABENDOC(NSBP),    ABENNH4(NSBP),    ABENNO3(NSBP),
     .          ABENPO4(NSBP),    ABENCOD(NSBP),    ABENDO(NSBP),
     .          ABENSA(NSBP),     ASSFWS(NSBP),     APCFWS(NSBP),
     .          APNFWS(NSBP),     APPFWS(NSBP),     APSFWS(NSBP),
     .          ACPOC(NSBP,3),    ACPON(NSBP,3),    ACPOP(NSBP,3),
     .          ACPIP(NSBP),      ACPOS(NSBP)
      DIMENSION AAPC(0:NBP)
      DIMENSION S1LN(NCP),    S2LN(NCP),    NCB(NCP),     AC(NCP)
      DIMENSION Q(0:NQFP),    ILB(0:NQFP),  IB(0:NQFP),   JB(0:NQFP),
     .          JRB(0:NQFP),  A(0:NQFP),    DIFF(0:NQFP)
      DIMENSION BL(0:NBP,3),  V1S(0:NBP),
     .          HMV(0:NBP),   HMBV(NSBP),   ZD(0:NBP),    BU(0:NBP)
      DIMENSION V1(0:NBP),    V2(0:NBP)
      DIMENSION PNC(NBP),     PND(NBP),     PNG(NBP),     BMC(NBP),
     .          BMD(NBP),     BMG(NBP),     MNLDOC(NBP),  FTCOD(NBP),
     .          NPP(NBP),     PC(NBP),      PD(NBP),      PG(NBP),
     .          NT(NBP),      PRC(NBP),     PRD(NBP),     PRG(NBP),
     .          DENIT(NBP),   RATOX(NBP)
      DIMENSION AFIC(NBP),    ANLC(NBP),    APLC(NBP),    AFID(NBP),
     .          ANLD(NBP),    APLD(NBP),    ASLD(NBP),    AFIG(NBP),
     .          ANLG(NBP),    APLG(NBP),    ANPP(NBP),    ARESP(NBP),
     .          AKE(NBP)
      DIMENSION T(0:NBP),     SSI(0:NBP),   BC(0:NBP),    BD(0:NBP),
     .          BG(0:NBP),    DOC(0:NBP),   LPOC(0:NBP),  RPOC(0:NBP),
     .          NH4(0:NBP),   NO3(0:NBP),   DON(0:NBP),   LPON(0:NBP),
     .          RPON(0:NBP),  PO4(0:NBP),   DOP(0:NBP),   LPOP(0:NBP),
     .          RPOP(0:NBP),  COD(0:NBP),   DO(0:NBP),    SU(0:NBP),
     .          SA(0:NBP),    SALT(0:NBP)
      DIMENSION WSS(0:NBP),   WSL(0:NBP),   WSR(0:NBP),   WSC(0:NBP),
     .          WSD(0:NBP),   WSG(0:NBP),   APC(0:NBP)
      DIMENSION PMC(NBP),     PMD(NBP),     PMG(NBP),     BMRC(NBP),
     .          BMRD(NBP),    BMRG(NBP),    BPRC(NBP),    BPRD(NBP),
     .          BPRG(NBP)
      DIMENSION DTT(NBP),     DTSSI(NBP),   DTBC(NBP),    DTBD(NBP),
     .          DTBG(NBP),    DTDOC(NBP),   DTLPOC(NBP),  DTRPOC(NBP),
     .          DTNH4(NBP),   DTNO3(NBP),   DTDON(NBP),   DTLPON(NBP),
     .          DTRPON(NBP),  DTPO4(NBP),   DTDOP(NBP),   DTLPOP(NBP),
     .          DTRPOP(NBP),  DTCOD(NBP),   DTDO(NBP),    DTSU(NBP),
     .          DTSA(NBP)
      DIMENSION KE(NBP),      FIC(NBP),     FID(NBP),     FIG(NBP),
     .          NLC(NBP),     NLD(NBP),     NLG(NBP),     PLC(NBP),
     .          PLD(NBP),     PLG(NBP),     SLD(NBP),     RESP(NBP),
     .          KESS(NBP),    IATTOP(NBP),  IATBOT(NBP)
      DIMENSION FTMNL(NBP),   FTHDR(NBP)
      DIMENSION RESPC(NBP),   DLALGC(NBP)
      DIMENSION IWCMNB(NBP),  IWCMPB(NBP),  IWCMCB(NBP),  IWCMSB(NBP),
     .          WCMNB(NBP),   WCMPB(NBP),   WCMCB(NBP),   WCMSB(NBP),
     .          DLWCMNB(NBP), DLWCMPB(NBP), DLWCMCB(NCP), DLWCKMNB(NBP),
     .          DLWCKMCB(NBP),S1FLXNB(NBP), S1FLXPB(NBP), S1FLXCB(NBP),
     .          S2FLXNB(NBP), S2FLXPB(NBP), S2FLXCB(NBP)
      DIMENSION BENDOC(NSBP), BENNH4(NSBP), BENNO3(NSBP), BENDON(NSBP),
     .          BENDOP(NSBP), BENPO4(NSBP), BENCOD(NSBP), BENDO(NSBP),
     .          BENSA(NSBP),  BENDEN(NSBP)
      DIMENSION BENDOCB(NSBP),BENNH4B(NSBP),BENNO3B(NSBP),BENPO4B(NSBP),
     .          BENCODB(NSBP),BENDOB(NSBP), BENSAB(NSBP)
      DIMENSION SBN(NSBP),    BBN(NSBP),    HMSBV(NSBP),  SFA(NSBP)
      DIMENSION WSSNET(NSBP), WSLNET(NSBP), WSRNET(NSBP), WSCNET(NSBP),
     .          WSDNET(NSBP), WSGNET(NSBP), VSED(NSBP),   VPMIX(NSBP),
     .          VDMIX(NSBP),  JDIAGP(NSBP), JDIAGN(NSBP), JDIAGC(NSBP),
     .          JDIAGS(NSBP), MTVEL(NSBP)
      DIMENSION PPFWS(NSBP),  PNFWS(NSBP),  PCFWS(NSBP),  PSFWS(NSBP),
     .          SSFWS(NSBP)
      DIMENSION METFN(NFLP),  S1FN(NFLP),   S2FN(NFLP),   HYDFN(NFLP),
     .          CBCFN(NFLP),  KEIFN(NFLP),  ATMFN(NFLP),  SAVFN(NFLP),
     .          BFIFN(NFLP)
      DIMENSION FRPPH1(3),    FRPPH2(3),    FRPPH3(3),  FRPOP(NSBP,3),
     .          FRNPH1(3),    FRNPH2(3),    FRNPH3(3),  FRPON(NSBP,3),
     .          FRCPH1(3),    FRCPH2(3),    FRCPH3(3),  FRPOC(NSBP,3)
      DIMENSION SSNAME(11)
      DIMENSION ERROR_MESS(36)
      DIMENSION HQCFA(NHQP/NLP)
      DIMENSION KDC(NBP),     KDCALG(NBP),   KLC(NBP),    KLCALG(NBP),
     .          KRC(NBP),     KRCALG(NBP),   KSUA(NBP),   KCOD(NBP),
     .          KDN(NBP),     KDNALG(NBP),   KDP(NBP),    KDPALG(NBP),
     .          KLN(NBP),     KLNALG(NBP),   KLP(NBP),    KLPALG(NBP),
     .          KRN(NBP),     KRNALG(NBP),   KRP(NBP),    KRPALG(NBP), 
     .          NTM(NBP),     KRDO(NBP)
      DIMENSION ALGDOC(NBP),  ALGPOC(NBP),   HDRLPOC(NBP), HDRRPOC(NBP) 
      DIMENSION ALGNH4(NBP),  ALGNO3(NBP),   ALGDON(NBP),  ALGPON(NBP),
     .          DENNO3(NBP),  MNLDON(NBP),   HDRLPON(NBP), HDRRPON(NBP)
      DIMENSION ALGPO4(NBP),  ALGDOP(NBP),   ALGPOP(NBP),  MNLDOP(NBP),
     .          HDRLPOP(NBP), HDRRPOP(NBP)
      DIMENSION PSD(NBP),     SAP(NBP),      ALGUP(NBP),   ALGRES(NBP)
      DIMENSION DORALG(NBP),  DOPR(NBP),
     .          DCOD(NBP),    DDOC(NBP),     NITRIF(NBP)
      DIMENSION A_T(NBP),     A_PC(NBP),     ABMC(NBP),    APRC(NBP),
     .          APD(NBP),     ABMD(NBP),     APRD(NBP),    APG(NBP),
     .          ABMG(NBP),    APRG(NBP)
      DIMENSION AALGDOC(NBP), AALGPOC(NBP),  ADENIT(NBP),  AMNLDOC(NBP),
     .          AHDRLPOC(NBP),AHDRRPOC(NBP)
      DIMENSION AALGNH4(NBP), AALGNO3(NBP),  AALGDON(NBP), AALGPON(NBP),
     .          ANT(NBP),     ADENNO3(NBP),  AMNLDON(NBP),AHDRLPON(NBP),
     .          AHDRRPON(NBP)
      DIMENSION AALGPO4(NBP), AALGDOP(NBP),  AALGPOP(NBP), AMNLDOP(NBP),
     .          AHDRLPOP(NBP),AHDRRPOP(NBP)
      DIMENSION APSD(NBP),    ASAP(NBP),     AALGUP(NBP),  AALGRES(NBP)
      DIMENSION ADO(NBP),     ADORALG(NBP),  ADOPR(NBP),   ADCOD(NBP),
     .          ADDOC(NBP),   ANITRIF(NBP)
      DIMENSION HMCV(NSBP),   DLZF(NLP)
      
***** Equivalence declarations

      EQUIVALENCE (NSB,        NBB)
      EQUIVALENCE (C1(0,1),    CSTAR(0,1))
      EQUIVALENCE (T(0),       C2(0,1)),     (SALT(0),     C2(0,2)),
     .            (SSI(0),     C2(0,3)),     (BC(0),       C2(0,4)),
     .            (BD(0),      C2(0,5)),     (BG(0),       C2(0,6)),
     .            (DOC(0),     C2(0,7)),     (LPOC(0),     C2(0,8)),
     .            (RPOC(0),    C2(0,9)),     (NH4(0),      C2(0,10)),
     .            (NO3(0),     C2(0,11)),    (DON(0),      C2(0,12)),
     .            (LPON(0),    C2(0,13)),    (RPON(0),     C2(0,14)),
     .            (PO4(0),     C2(0,15)),    (DOP(0),      C2(0,16)),
     .            (LPOP(0),    C2(0,17)),    (RPOP(0),     C2(0,18)),
     .            (COD(0),     C2(0,19)),    (DO(0),       C2(0,20)),
     .            (SU(0),      C2(0,21)),    (SA(0),       C2(0,22))
      EQUIVALENCE (DTT(1),     DTC(1,1)),    (DTSSI(1),    DTC(1,3)),
     .            (DTBC(1),    DTC(1,4)),    (DTBD(1),     DTC(1,5)),
     .            (DTBG(1),    DTC(1,6)),    (DTDOC(1),    DTC(1,7)),
     .            (DTLPOC(1),  DTC(1,8)),    (DTRPOC(1),   DTC(1,9)),
     .            (DTNH4(1),   DTC(1,10)),   (DTNO3(1),    DTC(1,11)),
     .            (DTDON(1),   DTC(1,12)),   (DTLPON(1),   DTC(1,13)),
     .            (DTRPON(1),  DTC(1,14)),   (DTPO4(1),    DTC(1,15)),
     .            (DTDOP(1),   DTC(1,16)),   (DTLPOP(1),   DTC(1,17)),
     .            (DTRPOP(1),  DTC(1,18)),   (DTCOD(1),    DTC(1,19)),
     .            (DTDO(1),    DTC(1,20)),   (DTSU(1),     DTC(1,21)),
     .            (DTSA(1),    DTC(1,22))
      EQUIVALENCE (FLXTTEM(1), FLUXT(1,1)),  (FLXTSAL(1),  FLUXT(1,2)),
     .            (FLXTSSI(1), FLUXT(1,3)),  (FLXTC(1),    FLUXT(1,4)),
     .            (FLXTD(1),   FLUXT(1,5)),  (FLXTG(1),    FLUXT(1,6)),
     .            (FLXTDOC(1), FLUXT(1,7)),  (FLXTLPOC(1), FLUXT(1,8)),
     .            (FLXTRPOC(1),FLUXT(1,9)),  (FLXTNH4(1),  FLUXT(1,10)),
     .            (FLXTNO3(1), FLUXT(1,11)), (FLXTDON(1),  FLUXT(1,12)),
     .            (FLXTLPON(1),FLUXT(1,13)), (FLXTRPON(1), FLUXT(1,14)),
     .            (FLXTPO4(1), FLUXT(1,15)), (FLXTDOP(1),  FLUXT(1,16)),
     .            (FLXTLPOP(1),FLUXT(1,17)), (FLXTRPOP(1), FLUXT(1,18)),
     .            (FLXTCOD(1), FLUXT(1,19)), (FLXTDO(1),   FLUXT(1,20)),
     .            (FLXTSU(1),  FLUXT(1,21)), (FLXTSA(1),   FLUXT(1,22))
      EQUIVALENCE (FLXSSSI(1), FLUXS(1,1)),  (FLXSC(1),    FLUXS(1,2)),
     .            (FLXSD(1),   FLUXS(1,3)),  (FLXSG(1),    FLUXS(1,4)),
     .            (FLXSPOC(1), FLUXS(1,5)),  (FLXSPON(1),  FLUXS(1,6)),
     .            (FLXSPO4(1), FLUXS(1,7)),  (FLXSPOP(1),  FLUXS(1,8)),
     .            (FLXSSI(1),  FLUXS(1,9))
      EQUIVALENCE (AFLXPOC(1), AFLUX(1,1)),  (AFLXDOC(1),  AFLUX(1,2)),
     .            (AFLXTC(1),  AFLUX(1,3)),  (AFLXDON(1),  AFLUX(1,4)),
     .            (AFLXDIN(1), AFLUX(1,5)),  (AFLXPON(1),  AFLUX(1,6)),
     .            (AFLXTN(1),  AFLUX(1,7)),  (AFLXDOP(1),  AFLUX(1,8)),
     .            (AFLXDIP(1), AFLUX(1,9)),  (AFLXPOP(1),  AFLUX(1,10)),
     .            (AFLXTP(1),  AFLUX(1,11))
      EQUIVALENCE (BENDOC(1),  BFLUX(1,1)),  (BENNH4(1),   BFLUX(1,2)),
     .            (BENNO3(1),  BFLUX(1,3)),  (BENDON(1),   BFLUX(1,4)),
     .            (BENPO4(1),  BFLUX(1,5)),  (BENDOP(1),   BFLUX(1,6)),
     .            (BENCOD(1),  BFLUX(1,7)),  (BENDO(1),    BFLUX(1,8)),
     .            (BENSA(1),   BFLUX(1,9))
      EQUIVALENCE (BENDOCB(1), BFLUXB(1,1)), (BENNH4B(1),  BFLUXB(1,2)),
     .            (BENNO3B(1), BFLUXB(1,3)), (BENPO4B(1),  BFLUXB(1,4)),
     .            (BENCODB(1), BFLUXB(1,5)), (BENDOB(1),   BFLUXB(1,6)),
     .            (BENSAB(1),  BFLUXB(1,7))

***** Common declarations

      COMMON /TVDLGC/ BOUNDARY_CONC,  SOURCE_ONE,     SOURCE_TWO,
     .                BENTHIC_FLUXES, SEDIMENT_CALC,  ALGAL_GROWTH,
     .                ATMOS_LOADS,    SUB_AQ_VEG_CALC,STEP_BOUNDARY,
     .                NO_SEDIMENT
      COMMON /HYDLGC/ POSITIVE_FLOW,  FLOW,           XY_DIFFUSION,
     .                Z_DIFFUSION,    DIAGNOSTICS,    ZGRID_HYDRO,
     .                CONSERVE_MASS,  MASS_BALANCE,   ASCII_HYDRO,
     .                SIGMA_HYDRO
      COMMON /AVGLGC/ AVERAGE_PLOTS,  QUALITY_DIAG,
     .                SEDIMENT_DIAG
      COMMON /SEDLGC/ SETTLING
      COMMON /HYDROC/ Q,      A,      DIFF,   HMBV,   HMSBV,  HMV,
     .                HQCFA,  SFA,    ILB,    IB,     JB,     JRB,
     .                VFN,    NVF,    DLT,    AHMDLT, NHMDLT, NWQMR,
     .                NHMR,   ZDFMUL, HMEND,  JDAY,   NIT
      COMMON /TVDSC/  S1L,    S1LB,   S1LN,   S2L,    S2LB,   S2LN,
     .                CB,     NCB,    NXCBC,  CBNX,   TMEND,
     .                CBOLD,  NXMET,  NXS1,   NXS2,   NXBFI,  NXKEI,
     .                NXATM,  NXSAV
      COMMON /DBLPR/  V1,     V2,     ELTMS
      COMMON /AVGC1/  AC,     NAC,    AC1,    AAPC
      COMMON /AVGC2/  ABENDOC,ABENNH4,ABENNO3,ABENPO4,ABENCOD,ABENDO,
     .                ABENSA, ASSFWS, APCFWS, APNFWS, APPFWS, APSFWS,
     .                ACPOC,  ACPON,  ACPOP,  ACPIP,  ACPOS 
      COMMON /AVGC3/  AFIC,   ANLC,   APLC,   AFID,   ANLD,   APLD,
     .                ASLD,   AFIG,   ANLG,   APLG,   ANPP,   ARESP,
     .                AKE
      COMMON /AVGC4/  A_T,    A_PC,   ABMC,   APRC,   APD,    ABMD,
     .                APRD,   APG,    ABMG,   APRG,   AALGDOC,AALGPOC,
     .                ADENIT, AMNLDOC,AHDRLPOC,       AHDRRPOC,
     .                AALGNH4,AALGNO3,AALGDON,AALGPON,ANT,    ADENNO3,
     .                AMNLDON,AHDRLPON,       AHDRRPON,       AALGPO4,
     .                AALGDOP,AALGPOP,AMNLDOP,AHDRLPOP,       AHDRRPOP, 
     .                APSD,   ASAP,   AALGUP, AALGRES,ADO,    ADORALG, 
     .                ADOPR,  ADCOD,  ADDOC,  ANITRIF
      COMMON /UNITNC/ HYD,    DIA,    CBC,    S1,     S2,     BFI,
     .                MET,    BFO,    KEI,    ATM,    STL,    AGR,
     .                SAV,    KFL,    ERR
      COMMON /GEOMC/  BL,     BU,     SBN,    BBN,
     .                ZD,     NB,     NQF,    NL,     NSB,    NHQF,
     .                NSQF
      COMMON /TEMPC/  KT,     TE
      COMMON /SOLIDC/ KADPO4, KADSA
      COMMON /OXYGNC/ AOCR,   AONT,   FTCOD,  FDOP,   KRDO
      COMMON /CODMDC/ KCOD,   TRCOD,  KTCOD
      COMMON /SILICC/ ASCD,   KSUA,   FSAP,   KTSUA,  TRSUA
      COMMON /NITROC/ TMNT,   NTM,    KTNT1,  KTNT2,  KHONT,  KHNNT,
     .                ANCC,   ANCD,   ANCG,   FNIC,   KDN,    KDNALG,
     .                KLN,    KLNALG, KRN,    KRNALG, KHNAVG, FNID,
     .                FNIG,   FNDC,   FNDD,   FNDG,   FNLC,   FNLD,
     .                FNLG,   FNRC,   FNRD,   FNRG,   FNIP,   FNDP,
     .                FNLP,   FNRP,   NT,     KHNDN,  DENIT,  ANDC
      COMMON /PHOSPC/ APC,    APCMIN, APCMAX, APCRAT, FPIC,   FPID,
     .                FPIG,   FPDC,   FPDD,   FPDG,   FPLC,   FPLD,
     .                FPLG,   FPRC,   FPRD,   FPRG,   FPIP,   FPDP,
     .                FPLP,   FPRP,   KDP,    KDPALG, KLP,    KLPALG,
     .                KRP,    KRPALG, KHPAVG, PO4DMAX
      COMMON /CARBC/  FCDC,   FCDD,   FCDG,   FCDP,   FCLP,   FCRP,
     .                KTMNL,  KTHDR,  TRMNL,  TRHDR,  FTMNL,  FTHDR,
     .                KDC,    KDCALG, KLC,    KLCALG, KRC,    KRCALG,
     .                MNLDOC, RATOX,  AANOX
      COMMON /ALGAEC/ TMC,    TMD,    TMG,    KTGC1,  KTGC2,  KTGD1,
     .                KTGD2,  KTGG1,  KTGG2,  KHPC,   KHPD,   KHPG,
     .                KHSD,   BMRC,   BMRD,   BMRG,   TRC,    TRD,
     .                TRG,    KTBC,   KTBD,   KTBG,   PMC,    PMD,
     .                PMG,    MBGM,   UCM,    UDM,    UGM,    FR,
     .                PRC,    PRD,    PRG,    BPRC,   BPRD,   BPRG,
     .                SCTOX
      COMMON /SED1C/  CPOP,   FLXPOP, CPON,   FLXPON, CPOC,   FLXPOC,
     .                CPOS,   CPIP,   CNO3,   CNH4,   FLXPOS, CTEMP,
     .                CDTEMP, BSVOL,  WSSNET, WSCNET, WSDNET, WSGNET,
     .                WSLNET, WSRNET, HSED,   VSED,   VPMIX,
     .                VDMIX,  JDIAGP, JDIAGN, JDIAGC, JDIAGS, RESTART_IN
      COMMON /SED2C/  FRPPH1, FRPPH2, FRPPH3, FRPOP,  FRNPH1, FRNPH2,
     .                FRNPH3, FRPON,  FRCPH1, FRCPH2, FRCPH3, FRPOC
      COMMON /SED4C/  SSNAME
      COMMON /CONSTC/ C1,     C2,     DTC
      COMMON /FLUXC/  FLUXT,  FLUXS,  AFLUX,  TRANSPORT_FLUXES
      COMMON /PRODC/  PC,     PD,     PG
      COMMON /RESPRC/ BMC,    BMD,    BMG
      COMMON /NPREFC/ PNC,    PND,    PNG
      COMMON /SETTLC/ WSS,    WSC,    WSD,    WSG,    WSL,    WSR
      COMMON /HALFSC/ KHNC,   KHND,   KHNG,   KHRC,   KHRD,   KHRG,
     .                KHOCOD, KHODOC
      COMMON /LIGHTC/ FD,     KE,     I0,     I0NX,   IATTOP, IATBOT,
     .                KHTIS,  CEXTC,  CEXTD,  CEXTG,  KECHL,
     .                IHC,    IHD,    IHG,    KESS
      COMMON /ATMC/   PRECIP, ATMNH4, ATMNO3, ATMDON, ATMPO4, ATMDOP
      COMMON /SAVC/   NSAV,   SAVB,   SAVAREA,SAVLPOC,SAVRPOC,SAVLPON,
     .                SAVRPON,SAVLPOP,SAVRPOP,SAVDO
      COMMON /BENTHC/ KSO,    KSNH4,  KSNO3,  KSPO4,  KSSA,   KSDOC,
     .                TRSO,   TRSNH4, TRSNO3, TRSPO4, TRSSA,  TRSDOC,
     .                KHSO,   SEDNO3, MTCNO3, BFLUX,  BFLUXB, BENDEN, 
     .                MTVEL
      COMMON /DIAPLC/ FIC,    FID,    FIG,    NLC,    NLD,    NLG,
     .                PLC,    PLD,    PLG,    SLD,    RESP,   NPP
      COMMON /BENPLC/ PPFWS,  PNFWS,  PCFWS,  PSFWS,  SSFWS
      COMMON /FILENC/ METFN,  S1FN,   S2FN,   HYDFN,  CBCFN,  KEIFN,
     .                ATMFN,  SAVFN,  BFIFN,  BFOFN,  ERRFN
      COMMON /POINTC/ HYDPTR, METPTR, CBCPTR, S1PTR,  S2PTR, BFIPTR,
     .                KEIPTR, ATMPTR, SAVPTR
      COMMON /REDUCC/ REDS1C, REDS1N, REDS1P, REDS2C, REDS2N, REDS2P,
     .                REDCBC, REDCBN, REDCBP
      COMMON /MASSBC/ DLALGC, RESPC
      COMMON OLDNXCBC
      COMMON /KINFLX/ KINETIC_FLUXES
      COMMON /KINCAR/ ALGDOC, ALGPOC, HDRLPOC, HDRRPOC 
      COMMON /KINNIT/ ALGNH4, ALGNO3, ALGDON,  ALGPON,
     .                DENNO3, MNLDON, HDRLPON, HDRRPON
      COMMON /KINPHO/ ALGPO4, ALGDOP, ALGPOP,  MNLDOP,
     .                HDRLPOP, HDRRPOP
      COMMON /KINSIL/ PSD,    SAP,    ALGUP,   ALGRES
      COMMON /KINOXY/ DORALG, DOPR,
     .                DCOD,   DDOC,   NITRIF
      COMMON /SIGHYD/ HMCV,   DLZF
      COMMON /TYPES/  KEI_TYPE,    BFI_TYPE
      COMMON /ERRS1/  ERROR_OPEN,  ERROR_CHECK
      COMMON /ERRS2/  NERRS
      COMMON /ERRS3/  ERROR_MESS


********************************************************************************
**                           Sediment Model Setup                             **
********************************************************************************

***** Variable declarations

      REAL    M1,       M2,        KAPPNH4F, KMNH4,    KMNH4O2,
     .        KAPPNO3F, K2NO3,     KAPPD1,   KAPPP1,   KMHSO2,
     .        KSI,      KMPSI,     KMO2DP,   KBENSTR,
     .        KAPPCH4,  KPDIAG,    KNDIAG,   KCDIAG,   KAPPNH4S,
     .        KAPPNO3S, O2CRITSI,  JSIDETR,  KLBNTH,   DPMIN
      REAL    K0H1D,    K0H1P,     K1H1D,    K1H1P,    K2H2D,
     .        K2H2P,    K3,        J1,       J2,       KMC1,
     .        KL12,     KL12NOM,   ISWBENS
      REAL    NH41TM1S, NH4T2TM1S, NO31TM1S, NO3T2TM1S
      REAL    NH40,     NH41,      NH42,     NH4T1,     NH4T2,
     .        NH41TM1,  NH4T2TM1,  JNH4,     NO30,      NO31,
     .        NO32,     NO3T1,     NO3T2,    NO31TM1,   NO3T2TM1,
     .        JNO3,     JHS,       JSI,      JPO4,      JCH4AQ,
     .        JCH4G,    JO2NH4
      REAL    KPOP1,    KPOP2,     KPOP3,    KPON1,     KPON2,
     .        KPON3,    KPOC1,     KPOC2,    KPOC3
      REAL    ISEDMN,   ISEDMP,    ISEDMC,   S2FLXN,    S2FLXP,
     .        S2FLXC,   S2FLXNB,   S2FLXPB,  S2FLXCB,   IWCMN,
     .        IWCMP,    IWCMC
      LOGICAL BENTHIC_OUTPUT, STEADY_STATE_SED

***** Dimension declarations

      DIMENSION KPDIAG(3),       DPTHTA(3),       KNDIAG(3),
     .          DNTHTA(3),       KCDIAG(3),       DCTHTA(3) 
      DIMENSION PON1TM1S(NSBP),  PON2TM1S(NSBP),  PON3TM1S(NSBP),
     .          POC1TM1S(NSBP),  POC2TM1S(NSBP),  POC3TM1S(NSBP),
     .          POP1TM1S(NSBP),  POP2TM1S(NSBP),  POP3TM1S(NSBP),
     .          PSITM1S(NSBP),   BENSTR1S(NSBP),  BFORMAXS(NSBP),
     .          ISWBENS(NSBP)
      DIMENSION NH41TM1S(NSBP),  NH4T2TM1S(NSBP), NO31TM1S(NSBP),
     .          NO3T2TM1S(NSBP), HS1TM1S(NSBP),   HST2TM1S(NSBP),
     .          SI1TM1S(NSBP),   SIT2TM1S(NSBP),  PO41TM1S(NSBP),
     .          PO4T2TM1S(NSBP)
      DIMENSION BURIALN(NSBP),   BURIALP(NSBP),   BURIALC(NSBP),
     .          DIAGENC(NSBP)
      DIMENSION ZHTADP(350),     ZHTADD(350),     ZHTANH4F(350),
     .          ZHTANO3F(350),   ZHTA2NO3(350),   ZHTAD1(350),
     .          ZHTAP1(350),     ZHTASI(350),     ZL12NOM(350),
     .          ZW12NOM(350),    ZHTAPON1(350),   ZHTAPON2(350),
     .          ZHTAPON3(350),   ZHTAPOC1(350),   ZHTAPOC2(350),
     .          ZHTAPOC3(350),   ZHTAPOP1(350),   ZHTAPOP2(350),
     .          ZHTAPOP3(350),   ZHTACH4(350),    ZHTANH4S(350),
     .          ZHTANO3S(350)
      DIMENSION AG3CFL(NSBP),    AG3NFL(NSBP),    AG3PFL(NSBP),
     .          ASDTMP(NSBP)

***** Equivalence declarations

      EQUIVALENCE (KPDIAG(1),KPOP1),    (KPDIAG(2),KPOP2),
     .            (KPDIAG(3),KPOP3),    (KNDIAG(1),KPON1),
     .            (KNDIAG(2),KPON2),    (KNDIAG(3),KPON3),
     .            (KCDIAG(1),KPOC1),    (KCDIAG(2),KPOC2),
     .            (KCDIAG(3),KPOC3),    (DPTHTA(1),THTAPOP1),
     .            (DPTHTA(2),THTAPOP2), (DPTHTA(3),THTAPOP3),
     .            (DNTHTA(1),THTAPON1), (DNTHTA(2),THTAPON2),
     .            (DNTHTA(3),THTAPON3), (DCTHTA(1),THTAPOC1),
     .            (DCTHTA(2),THTAPOC2), (DCTHTA(3),THTAPOC3)

***** Common declarations

      COMMON /INPUTC/ M1,       M2,       DP,       W2,       DD,
     .                THTADP,   THTADD,   KAPPNH4,  PIENH4,   THTANH4,
     .                KMNH4,    KMNH4O2,  KAPP1NO3, K2NO3,    THTANO3,
     .                KAPPD1,   KAPPP1,   PIE1S,    PIE2S,    THTAPD1,
     .                KMHSO2,   KSI,      CSISAT,   DPIE1SI,  PIE2SI,
     .                H2,       THTASI,   KMPSI,    DPIE1PO4, PIE2PO4,
     .                O2CRIT,   KMO2DP,   TEMPBEN,
     .                KBENSTR,  KAPPCH4,  THTACH4,  KPDIAG,   DPTHTA,
     .                KNDIAG,   DNTHTA,   KCDIAG,   DCTHTA,   DPMIN,
     .                DPIE1PO4F,DPIE1PO4S,O2CRITSI, JSIDETR,  KLBNTH
      COMMON /LOGSC1/ BENTHIC_OUTPUT,     STEADY_STATE_SED
      COMMON /SEDPOM/ PON1TM1S, PON2TM1S, PON3TM1S, POC1TM1S, POC2TM1S,
     .                POC3TM1S, POP1TM1S, POP2TM1S, POP3TM1S, PSITM1S,
     .                BENSTR1S, BFORMAXS, ISWBENS
      COMMON /CONCC1/ NH41TM1S, NH4T2TM1S,NO31TM1S, NO3T2TM1S,HS1TM1S,
     .                HST2TM1S, SI1TM1S,  SIT2TM1S, PO41TM1S, PO4T2TM1S
      COMMON /CONCC2/ NH40,     NH41,     NH42,     NH4T1,    NH4T2,
     .                NH41TM1,  NH4T2TM1, JNH4,     NO30,     NO31,
     .                NO32,     NO3T1,    NO3T2,    NO31TM1,  NO3T2TM1,
     .                JNO3,     HS0,      HS1,      HS2,      HST1,
     .                HST2,     HS1TM1,   HST2TM1,  JHS,      SI0,
     .                SI1,      SI2,      SIT1,     SIT2,     SI1TM1,
     .                SIT2TM1,  JSI,      PO40,     PO41,     PO42,
     .                PO4T1,    PO4T2,    PO41TM1,  PO4T2TM1, JPO4,
     .                JCH4AQ,   JCH4G
      COMMON /DIAGC/  PON1,     PON1TM1,  PON2,     PON2TM1,  PON3,
     .                PON3TM1,  POC1,     POC1TM1,  POC2,     POC2TM1,
     .                POC3,     POC3TM1,  POP1,     POP1TM1,  POP2,
     .                POP2TM1,  POP3,     POP3TM1,  XJN,      XJC,
     .                XJP,      PSI,      PSITM1,   XJCNO3,   XJCO2,
     .                XJC1,     JO2NH4
      COMMON /NLPARS/ K0H1D,    K0H1P,    K1H1D,    K1H1P,    K2H2D,
     .                K2H2P,    K3,       PIE1,     PIE2,     J1,
     .                J2,       KMC1,     W12,      KL12,     TEMPD,
     .                O20,      CH4SAT,   SAL
      COMMON /TSC/    SALTSW,   DLTS,     SALTND
      COMMON /STOREC/ XAPPNH4,  XAPPD1,   XAPPP1,   XAPP1NO3, XK2NO3,
     .                XKSI,     XAPPCH4,  TEMP20,   TEMP202,  FD1,
     .                FP1,      FD2,      FP2,      SOD,      CSOD,
     .                S,        W12NOM,   BENSTR,   BENSTRS,  BENSTR1,
     .                ISWBEN,   BFORMAX,  KL12NOM
      COMMON /THETAS/ ZHTADP,   ZHTADD,   ZHTANH4,  ZHTANO3,  ZHTA2NO3,
     .                ZHTAD1,   ZHTAP1,   ZHTASI,   ZL12NOM,  ZW12NOM,
     .                ZHTAPON1, ZHTAPON2, ZHTAPON3, ZHTAPOC1, ZHTAPOC2,
     .                ZHTAPOC3, ZHTAPOP1, ZHTAPOP2, ZHTAPOP3, ZHTACH4
      COMMON /MASSGC/ SEDMN,    SEDMP,    SEDMC,    ISEDMN,   ISEDMP,
     .                ISEDMC,   BURIALN,  BURIALP,  BURIALC,  DIAGENC,
     .                S1FLXN,   S1FLXP,   S1FLXC,   S2FLXN,   S2FLXP,
     .                S2FLXC,   ATMFLXN,  ATMFLXP,  BENFLXPN, BENFLXPP,
     .                BENFLXPC, BENFLXDN, BENFLXDP, DLWCKMN,  DLWCKMC,
     .                BNDFLXN,  BNDFLXP,  BNDFLXC,  DLSEDKN,  DLSEDKC,
     .                CMASS,    IWCMN,    IWCMP,    IWCMC,    IWCMS,
     .                BURIALFLXN,  BURIALFLXP,  BURIALFLXC
      COMMON /MASSBC/ S1FLXNB,  S1FLXPB,  S1FLXCB,  S2FLXNB,  S2FLXPB,
     .                S2FLXCB,  ATMFLXNB, ATMFLXPB, ATMFLXCB, BENFLXPNB,
     .                BENFLXDNB,BENFLXPPB,BENFLXDPB,BENFLXPCB,DLWCKMNB,
     .                DLWCKMCB, DLSEDKNB, DLSEDKCB, JDAYMBL,
     .                BURIALFLXNB, BURIALFLXPB, BURIALFLXCB
      COMMON /SEDINT/ TINTIM,   AG3CFL,   AG3NFL,   AG3PFL,   ASDTMP
