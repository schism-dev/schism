      MODULE STAT_POOL 
      IMPLICIT NONE

      CHARACTER                        :: DUMP*500, DUMP2*15 
      CHARACTER                        :: DUMP3*1
      CHARACTER(LEN=30), ALLOCATABLE   :: STATIONNAMES (:)
      CHARACTER(LEN=30), ALLOCATABLE   :: STATIONNAMES_O (:)
      CHARACTER(LEN=15), ALLOCATABLE   :: DATUM_SP     (:)
      CHARACTER(LEN=15), ALLOCATABLE   :: DATUM_O      (:,:)
      CHARACTER(LEN=15), ALLOCATABLE   :: FMT_O        (:)  
      CHARACTER(LEN=15), ALLOCATABLE   :: FMT_O_F      (:)
      CHARACTER(LEN=3)                 :: CHTMP 
      CHARACTER(LEN=4)                 :: YYYY         !YYYY MM DD hh
      CHARACTER(LEN=2)                 :: MM, DD, hh, MINUTE, SEC ! MIN(x,y) is a fortran function
      CHARACTER(LEN=30)                :: FMT_ERG_HEADER 
      CHARACTER(LEN=30)                :: FMT_NAM_HEADER
      CHARACTER(LEN=30)                :: FMT_ERG_RESULT 
      CHARACTER(LEN=30)                :: FMT_ERG_ALL_STAT
      CHARACTER(LEN=30),  ALLOCATABLE  :: BNAMES(:)
      CHARACTER(LEN=2)                 :: TMPCHAR

      REAL*8, PARAMETER                :: SMALL = 0.0000000001d0 
      
      LOGICAL                          :: FOUND
      
      INTEGER  , ALLOCATABLE           :: N_OBSRV_ERR(:), N_OUT_TIME(:), N_STAT(:)
      
      INTEGER                          :: MSC_O_MAX = 47
      
      REAL*8                           :: PI, DUMP_R
      INTEGER                          :: N_DT_O_MAX, DUMP_I
      INTEGER, ALLOCATABLE             :: N_DT_O (:)
      INTEGER, ALLOCATABLE             :: MSC_O(:), N_DT_ERR(:,:)

      INTEGER                          :: DATAEND, DUMP_INT, COUNTER
      INTEGER                          :: BUOYS, BUOYS_O, N_DT_SP, ERG, MSC, MDC, ISEMAX
      
      REAL*8                             :: M0, M1, M2, INCR_DT, INTER_DT, INTER, PI2, DDIR
      REAL*8                             :: TMP1, TMP2, TMP3, TMP4, CUT_OFF, DUMMY, DX, DIFF_DX
      REAL*8                             :: DT, DIFF_DT, YINTER, MINMAXFREQ_O, FRINTF, SGHIG, SGLOW
          
      REAL*8, ALLOCATABLE                :: FREQ_SP  (:)  , DFREQ_SP(:)
      REAL*8, ALLOCATABLE                :: SPSIG(:), SPDIR(:)
      REAL*8, ALLOCATABLE                :: E_SP     (:)  , NDIR  (:), DSPR  (:), COSTH(:), SINTH(:)
      REAL*8, ALLOCATABLE                :: EMAX_O(:,:), EMAX_SP(:,:)
      REAL*8, ALLOCATABLE                :: AUX_M0_SP(:)  , AUX_M1_SP(:), AUX_M2_SP(:), DS_INCR_SP(:)
      
      REAL*8, ALLOCATABLE                :: HS_SP (:,:), TM01_SP(:,:), TM02_SP(:,:), TP_SP(:,:), TP_O(:,:)
      REAL*8, ALLOCATABLE                :: HS_SP_CLEAN (:,:), TM01_SP_CLEAN(:,:), TM02_SP_CLEAN(:,:), TP_SP_CLEAN(:,:)

      REAL*8, ALLOCATABLE                :: DIFF_HS_SP(:,:), DIFF2_HS_SP(:,:), ABS_DIFF_HS_SP(:,:)
      REAL*8, ALLOCATABLE                :: DIFF_TM01_SP(:,:), DIFF2_TM01_SP(:,:), ABS_DIFF_TM01_SP(:,:)
      REAL*8, ALLOCATABLE                :: DIFF_TM02_SP(:,:), DIFF2_TM02_SP(:,:), ABS_DIFF_TM02_SP(:,:)
      REAL*8, ALLOCATABLE                :: DIFF_TP_SP(:,:), DIFF2_TP_SP(:,:), ABS_DIFF_TP_SP(:,:)


      REAL*8     , ALLOCATABLE           :: DIFF3_HS_SP(:,:), DIFF3_TM01_SP(:,:), DIFF3_TM02_SP(:,:), DIFF3_TP_SP(:,:)
      REAL*8     , ALLOCATABLE           :: DIFF4_HS_SP(:,:), DIFF4_TM01_SP(:,:), DIFF4_TM02_SP(:,:), DIFF4_TP_SP(:,:)
      
      REAL*8     , ALLOCATABLE           :: DIFF1_ESP(:,:), DIFF2_ESP(:,:), MEAN_OSP(:,:), MEAN_BSP(:,:), BIAS_SP(:,:)
      REAL*8     , ALLOCATABLE           :: DIFF1_K(:,:), DIFF2_K(:,:), DIFF3_K(:,:), DIFF4_K(:,:), KORR_SP(:,:)
      
      REAL*8     , ALLOCATABLE           :: MEAN_HS_O(:)  , BIAS_HS_SP(:)  , MAE_HS_SP(:)
      REAL*8     , ALLOCATABLE           :: RMS_HS_SP(:)  , SCI_HS_SP(:)   , MEAN_HS_SP (:)
      REAL*8     , ALLOCATABLE           :: MEAN_TM01_O(:), BIAS_TM01_SP(:), MAE_TM01_SP(:)
      REAL*8     , ALLOCATABLE           :: RMS_TM01_SP(:), SCI_TM01_SP(:) , MEAN_TM01_SP(:)
      REAL*8     , ALLOCATABLE           :: MEAN_TM02_O(:), BIAS_TM02_SP(:), MAE_TM02_SP(:)
      REAL*8     , ALLOCATABLE           :: RMS_TM02_SP(:), SCI_TM02_SP(:) , MEAN_TM02_SP(:)  
      REAL*8     , ALLOCATABLE           :: MEAN_TP_O(:)  , BIAS_TP_SP(:)  , MAE_TP_SP(:)
      REAL*8     , ALLOCATABLE           :: RMS_TP_SP(:)  , SCI_TP_SP(:)   , MEAN_TP_SP (:)
      REAL*8     , ALLOCATABLE           :: KORR_HS(:), KORR_TM01(:), KORR_TM02(:), KORR_TP(:)

      REAL*8                             :: MEAN_O_HS_ALL  , MEAN_SP_HS_ALL
      REAL*8                             :: MEAN_O_TM01_ALL, MEAN_SP_TM01_ALL
      REAL*8                             :: MEAN_O_TM02_ALL, MEAN_SP_TM02_ALL
      REAL*8                             :: MEAN_O_TP_ALL  , MEAN_SP_TP_ALL
      REAL*8                             :: BIAS_HS_ALL  
      REAL*8                             :: BIAS_TM01_ALL
      REAL*8                             :: BIAS_TM02_ALL
      REAL*8                             :: RMS_HS_ALL  
      REAL*8                             :: RMS_TM01_ALL
      REAL*8                             :: RMS_TM02_ALL
      REAL*8                             :: MAE_HS_ALL  
      REAL*8                             :: MAE_TM01_ALL
      REAL*8                             :: MAE_TM02_ALL  
      REAL*8                             :: SCI_HS_ALL  
      REAL*8                             :: SCI_TM01_ALL
      REAL*8                             :: SCI_TM02_ALL      
      REAL*8                             :: KORR_HS_ALL  
      REAL*8                             :: KORR_TM01_ALL
      REAL*8                             :: KORR_TM02_ALL
      REAL*8                             :: BIAS_TP_ALL
      REAL*8                             :: RMS_TP_ALL
      REAL*8                             :: MAE_TP_ALL
      REAL*8                             :: SCI_TP_ALL
      REAL*8                             :: KORR_TP_ALL
      
      LOGICAL, PARAMETER               :: LWRITEALL = .TRUE.
       
      REAL*8, ALLOCATABLE                ::  E_O (:), E_OS(:,:), E_B(:,:), NDIR_B(:,:), DSPR_B(:,:)
      REAL*8, ALLOCATABLE                ::  E_OSF(:,:), E_BT(:,:), RMS_B(:,:)
      REAL*8, ALLOCATABLE                ::  AUX_M0_O(:) , AUX_M1_O(:), AUX_M2_O(:)   
      REAL*8, ALLOCATABLE                ::  HS_O (:,:), TM01_O(:,:), TM02_O(:,:)
      REAL*8, ALLOCATABLE                ::  HS_O_CLEAN (:,:), TM01_O_CLEAN(:,:), TM02_O_CLEAN(:,:)
      REAL*8, ALLOCATABLE                ::  TP_O_CLEAN (:,:)
      
      REAL*8, ALLOCATABLE                ::  FREQ_O(:,:), DFREQ_O(:,:)
      DOUBLE PRECISION, ALLOCATABLE    ::  ZEIT_SP(:)
      DOUBLE PRECISION, ALLOCATABLE    ::  ZEIT_O(:,:) 
      DOUBLE PRECISION                 ::  XMJD

      END MODULE

