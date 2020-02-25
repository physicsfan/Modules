MODULE irjj
!*******************************************************************
!                                                                  *  
!   This MODULE groups together the subroutines in the series      *
!   gg1112 - gg1234, rather than relying on calls to the           *
!   subroutines in separate files.                                 *
!                                                                  *   
!                                                                  *
!   SUBROUTINES written by  G. Gaigalas                            *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   MODULE created by A. Senchuk                    February 2020  *  
!                                                                  *
!*******************************************************************
  

  use wignerj
  USE mtjj

  IMPLICIT NONE

  !> Data needed for rangular
  INTEGER, DIMENSION(63) :: IMPTJJ, IMGTJJ, IMPNJJ, IMGNJJ
  DATA IMPNJJ/2,1,4,2*3,3*9,3*6,6*18,8*12,20*46,18*26/
  DATA IMGNJJ/2,1,5,2*3,3*11,3*8,6*25,8*17,20*63,18*45/
  DATA IMPTJJ/1,2,3,2*4,3*6,3*9,6*12,8*18,20*26,18*46/
  DATA IMGTJJ/1,2,3,2*5,3*8,3*11,6*17,8*25,20*45,18*63/

  INTEGER, DIMENSION(6) :: IMPTJJ9, IMGTJJ9, IMPNJJ9, IMGNJJ9
  DATA IMPTJJ9/301,5*302/
  DATA IMGTJJ9/301,5*306/
  DATA IMPNJJ9/302,5*301/
  DATA IMGNJJ9/306,5*301/


  INTEGER, DIMENSION(189) :: IMPTJJ11, IMGTJJ11, IMPNJJ11, IMGNJJ11
  DATA IMPTJJ11/1,6*2,8,7*9,16,8*17,25,9*26,35,10*36,  &
       46,11*47,58,12*59,71,13*72,85,14*86,100,15*101,      &
       116,16*117,133,17*134,151,18*152,170,19*171/
  DATA IMGTJJ11/1,6*7,8,7*15,16,8*24,25,9*34,35,10*45, &
       46,11*57,58,12*70,71,13*84,85,14*99,100,15*115,      &
       116,16*132,133,17*150,151,18*169,170,19*189/
  DATA IMPNJJ11/2,6*1,9,7*8,17,8*16,26,9*25,36,10*35,  &
       47,11*46,59,12*58,72,13*71,86,14*85,101,15*100,      &
       117,16*116,134,17*133,152,18*151,171,19*100/
  DATA IMGNJJ11/7,6*1,15,7*8,24,8*16,34,9*25,45,10*35, &
       57,11*46,70,12*58,84,13*71,99,14*85,115,15*100,      &
       132,16*116,150,17*133,169,18*151,189,19*100/


CONTAINS


  INTEGER FUNCTION ITJJ(IK,ID,KG,BK,BD,IBT,BT,KG1,ITP,ITG,IQ)

    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)                :: KG, IQ
    INTEGER,      INTENT(OUT)               :: KG1, ITP, ITG
    INTEGER,      INTENT(IN),  DIMENSION(7) :: IK, ID
    INTEGER,      INTENT(OUT), DIMENSION(7) :: IBT
    REAL(KIND=8), INTENT(IN),  DIMENSION(3) :: BK, BD
    REAL(KIND=8), INTENT(OUT), DIMENSION(3) :: BT
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: ITK, ITD, ITP1, ITG1, ITTK
    !-----------------------------------------------
    ITJJ=0
    IF(ID(3) > 37) RETURN
    KG1=2*KG
    IF(ITTK(ID(6),IK(6),KG1) == 0)RETURN
    ITK=IK(1)
    ITD=ID(1)
    IF(ID(3) < 9) THEN
       ITP1=IMPTJJ(ITK)
       ITP=IMPTJJ(ITD)
       IF(ITP1 /= ITP)RETURN
       ITG1=IMGTJJ(ITK)
       ITG=IMGTJJ(ITD)
    ELSEIF(ID(3) == 9) THEN
       IF(ITK > 300) THEN
          IF(ITD < 300) CALL MES(51)
          IF(ID(4) > 2) CALL MES(11)
          IF(IK(4) > 2) CALL MES(11)
          ITK=ITK-300
          ITD=ITD-300
          ITP1=IMPTJJ9(ITK)
          ITP=IMPTJJ9(ITD)
          IF(ITP1 /= ITP)RETURN
          ITG1=IMGTJJ9(ITK)
          ITG=IMGTJJ9(ITD)
       ELSE
          PRINT*, "ERROR in ITJJ"
          STOP
       END IF
    ELSE
       IF(ID(4) > 2) CALL MES(11)
       IF(IK(4) > 2) CALL MES(11)
       ITP1=IMPTJJ11(ITK)
       ITP=IMPTJJ11(ITD)
       IF(ITP1 /= ITP)RETURN
       ITG1=IMGTJJ11(ITK)
       ITG=IMGTJJ11(ITD)
    ENDIF
    IF(ITG1 /= ITG)RETURN
    ITJJ=1
    IBT(2)=ID(2)
    IBT(3)=ID(3)
    IBT(4)=ID(4)+IQ
    BT(3)=BD(3)+HALF*DBLE(IQ)
    RETURN
  END FUNCTION ITJJ



  INTEGER FUNCTION ITJJ2(IK,ID,KG1,BK,BD,IBT,BT,ITP,ITG,IQ)
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)                :: KG1, IQ
    INTEGER,      INTENT(OUT)               :: ITP, ITG
    INTEGER,      INTENT(IN),  DIMENSION(7) :: IK, ID
    INTEGER,      INTENT(OUT), DIMENSION(7) :: IBT
    REAL(kind=8), INTENT(IN),  DIMENSION(3) :: BK, BD
    REAL(kind=8), INTENT(OUT), DIMENSION(3) :: BT
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: ITK, ITD, ITP1, ITG1, ITTK
    !-----------------------------------------------
    ITJJ2=0
    IF(ID(3) > 37) RETURN
    IF(ITTK(ID(6),IK(6),KG1) == 0)RETURN
    ITK=IK(1)
    ITD=ID(1)
    IF(ID(3) < 9) THEN
       ITP1=IMPNJJ(ITK)
       ITP=IMPTJJ(ITD)
       IF(ITP1.NE.ITP)RETURN
       ITG1=IMGNJJ(ITK)
       ITG=IMGTJJ(ITD)
    ELSEIF(ID(3) == 9) THEN
       IF(ITK > 300) THEN
          IF(ITD < 300) CALL MES(52)
          IF(ID(4) > 2) CALL MES(12)
          IF(IK(4) > 2) CALL MES(12)
          ITK=ITK-300
          ITD=ITD-300
          ITP1=IMPNJJ9(ITK)
          ITP=IMPTJJ9(ITD)
          IF(ITP1 /= ITP)RETURN
          ITG1=IMGNJJ9(ITK)
          ITG=IMGTJJ9(ITD)
       ELSE
          PRINT*, "ERROR in ITJJ2"
          STOP
       END IF
    ELSE
       IF(ID(4) > 2) CALL MES(12)
       IF(IK(4) > 2) CALL MES(12)
       ITP1=IMPNJJ11(ITK)
       ITP=IMPTJJ11(ITD)
       IF(ITP1 /= ITP)RETURN
       ITG1=IMGNJJ11(ITK)
       ITG=IMGTJJ11(ITD)
    ENDIF
    IF(ITG1 /= ITG)RETURN
    ITJJ2=1
    IBT(2)=ID(2)
    IBT(3)=ID(3)
    IBT(4)=ID(4)+IQ
    BT(3)=BD(3)+HALF*DBLE(IQ)
    RETURN
  END FUNCTION ITJJ2


  INTEGER FUNCTION ITJJ3(IK,ID,KG1,BK,BD,IBT,BT,ITP,ITG,IQ)
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)                :: KG1, IQ
    INTEGER,      INTENT(OUT)               :: ITP, ITG
    INTEGER,      INTENT(IN),  DIMENSION(7) :: IK, ID
    INTEGER,      INTENT(OUT), DIMENSION(7) :: IBT
    REAL(kind=8), INTENT(IN),  DIMENSION(3) :: BK, BD
    REAL(kind=8), INTENT(OUT), DIMENSION(3) :: BT
    !      DIMENSION ID(7),IK(7),IBT(7),BT(3),BD(3),BK(3)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: ITK, ITD, ITP1, ITG1, ITTK
    !-----------------------------------------------
    ITJJ3=0
    IF(ID(3) > 37) RETURN
    IF(ITTK(ID(6),IK(6),KG1) == 0)RETURN
    ITK=IK(1)
    ITD=ID(1)
    IF(ID(3) < 9) THEN
       ITP1=IMPTJJ(ITK)
       ITP=IMPNJJ(ITD)
       IF(ITP1 /= ITP)RETURN
       ITG1=IMGTJJ(ITK)
       ITG=IMGNJJ(ITD)
    ELSEIF(ID(3) == 9) THEN
       IF(ITK > 300) THEN
          IF(ITD < 300) CALL MES(53)
          IF(ID(4) > 2) CALL MES(13)
          IF(IK(4) > 2) CALL MES(13)
          ITK=ITK-300
          ITD=ITD-300
          ITP1=IMPTJJ9(ITK)
          ITP=IMPNJJ9(ITD)
          IF(ITP1 /= ITP)RETURN
          ITG1=IMGTJJ9(ITK)
          ITG=IMGNJJ9(ITD)
       ELSE
          PRINT*, "ERROR in ITJJ3"
          STOP
       ENDIF
    ELSE
       IF(ID(4) > 2) CALL MES(13)
       IF(IK(4) > 2) CALL MES(13)
       ITP1=IMPTJJ11(ITK)
       ITP=IMPNJJ11(ITD)
       IF(ITP1 /= ITP)RETURN
       ITG1=IMGTJJ11(ITK)
       ITG=IMGNJJ11(ITD)
    ENDIF
    IF(ITG1 /= ITG)RETURN
    ITJJ3=1
    IBT(2)=ID(2)
    IBT(3)=ID(3)
    IBT(4)=ID(4)+IQ
    BT(3)=BD(3)+HALF*DBLE(IQ)
    RETURN
  END FUNCTION ITJJ3


  SUBROUTINE RMEAJJ(LL,IT,LQ,J,ITS,LQS,J1S,COEF)
    !*******************************************************************
    !   Written by  G. Gaigalas                                        *
    !   The last modification made by G. Gaigalas       October  2017  *
    !   Restructured by A. Senchuk                     September 2019  *
    !                                                                  *
    !*******************************************************************

    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)  :: LL, IT, LQ, J, ITS, LQS, J1S
    REAL(kind=8), INTENT(OUT) :: COEF
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: I1, I2, IE1,  JT, JTS
    INTEGER, DIMENSION(1,2) :: IDS3, IDV3
    INTEGER, DIMENSION(3,3) :: IDS5, IDV5
    INTEGER, DIMENSION(6,8) :: IDS7, IDV7
    !-----------------------------------------------
    DATA IDS3/12,20/
    DATA IDV3/1,1/
    DATA IDS5/-24,2*0,-30,120,-90,-54,-48,330/
    DATA IDV5/4*1,2*7,1,2*7/
    DATA IDS7/2*0,-40,3*0,-54,-33,-40,-65,30,0,88,-1,0,-1755,-40,  &
         0,198,-36,-72,4500,-234,-360,-78,156,0,108,-30,96,66,49,0,27,  &
         130,136,0,195,-104,-245,-624,1224,3*0,2720,-2448,-7752/
    DATA IDV7/6*1,7,2*1,7,2*1,7,2*1,77,11,1,7,11,1,77,2*11,35,5,1, &
         91,1,13,2*5,2*1,2*7,1,11,1,3*11,3*1,143,77,91/
    !
    COEF=ZERO
    IF(LL > 37) RETURN
    IF(LL <= 9) THEN
       IF(IT < 300) THEN
          IF(IMPTJJ(IT) /= IMPNJJ(ITS)) RETURN
          IF(IT < ITS) THEN
             JT=IT
             JTS=ITS
          ELSE
             JT=ITS
             JTS=IT
          ENDIF
       ENDIF
    ENDIF
    IF(LL == 1) THEN
       COEF=-2
    ELSEIF(LL == 3)THEN
       I1=JT-2
       I2=JTS-3
       COEF=-DSQRT(DBLE(IDS3(I1,I2))/DBLE(IDV3(I1,I2)))
    ELSEIF(LL == 5)THEN
       I1=JT-5
       I2=JTS-8
       IF(IDS5(I1,I2) >= 0) THEN
          COEF=DSQRT(DBLE(IDS5(I1,I2))/DBLE(IDV5(I1,I2)))
       ELSE
          COEF=-DSQRT(-DBLE(IDS5(I1,I2))/DBLE(IDV5(I1,I2)))
       ENDIF
    ELSEIF(LL == 7)THEN
       I1=JT-11
       I2=JTS-17
       IF(IDS7(I1,I2) >= 0) THEN
          COEF=DSQRT(DBLE(IDS7(I1,I2))/DBLE(IDV7(I1,I2)))
       ELSE
          COEF=-DSQRT(-DBLE(IDS7(I1,I2))/DBLE(IDV7(I1,I2)))
       ENDIF
    ELSEIF(LL == 9) THEN
       IF(IT > 300) THEN
          IF(ITS < 300) THEN
             WRITE(6,'(A)') ' ERROR IN  RMEAJJ '
             STOP
          ENDIF
          IF(LL == J1S) THEN
             CALL RMEAJJ11(IT,ITS,LL,COEF)
          ELSE
             CALL RMEAJJ11(ITS,IT,LL,COEF)
             IE1=LQ-LQS+J-J1S+LL-1
             IF((IE1/4)*4 /= IE1)COEF=-COEF
          ENDIF
       ELSE
          CALL RMEAJJ9(IT,LQ,J,ITS,LQS,J1S,COEF)
          IF(IT > ITS) THEN
             IF(MOD(LQ+J-LQS-J1S+LL-1,4) /= 0) COEF=-COEF
          ENDIF
          WRITE(0,'(A)') ' KLAIDA RMEAJJ SUB. '
          STOP
       END IF
    ELSE
       IF(LL == J1S) THEN
          CALL RMEAJJ11(IT,ITS,LL,COEF)
       ELSE
          CALL RMEAJJ11(ITS,IT,LL,COEF)
          IE1=LQ-LQS+J-J1S+LL-1
          IF((IE1/4)*4 /= IE1)COEF=-COEF
       END IF
    ENDIF
    IF(LL < 9) THEN
       IF(IT > ITS) THEN
          IF(MOD(LQ+J-LQS-J1S+LL-1,4) /= 0) COEF=-COEF
       ENDIF
    ENDIF
    RETURN

  CONTAINS


    SUBROUTINE RMEAJJ11(J1,J2,LL,S)

      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER,      INTENT(IN)  :: LL, J1, J2
      REAL(kind=8), INTENT(OUT) :: S
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER      :: LQ, LV, LQS, LVS, J, J1S, N, INN
      REAL(kind=8) :: A1, A4, Q, QQ, QS, QM, QMS
      !-----------------------------------------------
      S=ZERO
      CALL RUMTJJ(J1,LL,LQ,LV,J)
      CALL RUMTJJ(J2,LL,LQS,LVS,J1S)
      QQ=HALF
      N=2
      Q=HALF*DBLE(LQ)
      QS=HALF*DBLE(LQS)
      QM=-HALF*DBLE(HALF*(LL+1)-N)
      QMS=-HALF*DBLE(HALF*(LL+1)-(N-1))
      CALL C0T5S(QS,QMS,QQ,Q,QM,A4)
      IF(DABS(A4) < EPS) RETURN
      A1=DBLE(N*(LQ+1)*(J+1))
      S=DSQRT(A1)/A4
      INN=-N-1
      IF((INN/2)*2 /= INN)S=-S
      RETURN
    END SUBROUTINE RMEAJJ11


    SUBROUTINE RMEAJJ9(IT,LQ,J,ITS,LQS,J1S,COEF)
      !
      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER,      INTENT(IN)  :: IT, LQ, J, ITS, LQS, J1S
      REAL(kind=8), INTENT(OUT) :: COEF
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER :: JT, JTS
      INTEGER, DIMENSION(20) :: IDS3,IDS4,IDS5,IDS6,IDS7,IDV7,IDS8,&
           IDV8,IDS9,IDS10,IDS11,IDV11,IDS12,IDV12,IDS13,IDV13,IDS14,   &
           IDS15,IDV15,IDS16,IDV16,IDS17,IDV17,IDS18,IDV18
      !-----------------------------------------------
      DATA IDS3/2*0,660,0,1664,0,1650,130,0,816,0,1680,8*0/
      DATA IDS4/2*0,-55770,-20592,6292,-13260,0,12740,-19110,     &
           79968,-47880,-25410,-19278,7*0/
      DATA IDS5/0,-9009,-4056,-260,585,16575,0,25200,4200,-7140,  &
           -475,504,-14280,13566,-4250,5*0/
      DATA IDS6/0,-4992,2028,0,-1920,0,-12870,7350,0,8160,        &
           0,1512,0,-3648,0,-9000,4*0/
      DATA IDS7/-2184,-63,-59904,-302460,5265,848691,0,-145152,   &
           217728,1049580,287337,5184,261120,691866,-750,0,-76608,3*0/
      DATA IDV7/253,23,19481,19481,161,253253,1,36179,            &
           36179,2*36179,299,3*36179,1,3289,3*1/
      DATA IDS8/1224,2652,188598,-31824,204,-12996,0,25500,-38250,&
           -768,81396,3213,-60543,-3066144,727776,207,-41553,3*0/
      DATA IDV8/1265,115,13915,2783,23,13915,1,3*2783,13915,115,  &
           2783,236555,47311,17,21505,3*1/
      DATA IDS9/3380,-2340,5460,-1400,3150,-3570,0,9000,1500,2550,&
           -5320,-10710,1050,1140,6300,8550,-330,5750,2*0/
      DATA IDS10/0,2160,6240,0,-32,0,-21450,-9610,0,2688,         &
           0,7140,0,20160,0,6156,2*0,-10164,0/
      DATA IDS11/0,132,-52728,196,-50,-24990,0,-21160,31740,26250,&
           12920,-357,9583,-344988,5700,171,-15,-4830,-84,0/
      DATA IDV11/1,5,6655,1331,11,1331,1,4*1331,55,1331,6655,     &
           1331,11,2*121,11,1/
      DATA IDS12/2*0,-209950,77520,12920,-25688,0,-4522,2261,     &
           -48640,-285144,931,273885,-112908,-2138580,137781,6654375,  &
           -59616,284089,0/
      DATA IDV12/2*1,2*9317,231,9317,1,3993,1331,22627,9317,44,   &
           90508,429913,158389,3740,156332,39083,17765,1/
      DATA IDS13/2*0,1530,13056,29376,720,0,-1890,-315,101124,    &
           -4560,-13965,-35131,13500,-685900,-1197,-28875,-5060,759,0/
      DATA IDV13/2*1,77,3*1001,1,2*143,2431,1001,572,9724,2431,   &
           17017,2*884,221,17,1/
      DATA IDS14/4*0,22848,0,-121550,4590,0,-32832,0,45220,0,     &
           -31680,0,82764,2*0,144716,0/
      DATA IDS15/5*0,2128,0,-1938,2907,-4860,-17136,-1309,-8505,  &
           15876,420,1287,-6075,-132,-253,-650/
      DATA IDV15/5*1,143,1,3*143,2717,52,572,247,13,20,988,19,95, &
           19/
      DATA IDS16/7*0,570,95,4104,504,-1463,-39501,-60516,-1596,   &
           3933,621,-840,-805,390/
      DATA IDV16/7*1,2*13,221,65,260,884,1105,221,68,3740,2*17,   &
           11/
      DATA IDS17/9*0,-5796,17664,-5313,1771,-16632,-88,693,       &
           -94269,192500,30030,24570/
      DATA IDV17/9*1,221,1235,65,221,1615,2*17,1615,7429,323,437/
      DATA IDS18/13*0,-15000,280,-1170,-48750,-15600,3510,-33930/
      DATA IDV18/13*1,323,2*17,3553,437,19,253/
      !
      COEF=ZERO
      IF(IT < ITS) THEN
         JT=IT-25
         JTS=ITS-45
      ELSE
         JT=ITS-25
         JTS=IT-45
      ENDIF
      IF(JTS == 1) THEN
         IF(JT == 7) COEF=-DSQRT(DBLE(60))
      ELSEIF(JTS == 2) THEN
         IF(JT == 8) COEF=-DSQRT(DBLE(12))
         IF(JT == 9) COEF=-DSQRT(DBLE(8))
      ELSEIF(JTS == 3)THEN
         COEF=-DSQRT(DBLE(IDS3(JT))/DBLE(33))
      ELSEIF(JTS == 4)THEN
         IF(IDS4(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS4(JT))/DBLE(3003))
         ELSE
            COEF=-DSQRT(-DBLE(IDS4(JT))/DBLE(3003))
         ENDIF
      ELSEIF(JTS == 5)THEN
         IF(IDS5(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS5(JT))/DBLE(715))
         ELSE
            COEF=-DSQRT(-DBLE(IDS5(JT))/DBLE(715))
         ENDIF
      ELSEIF(JTS == 6)THEN
         IF(IDS6(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS6(JT))/DBLE(143))
         ELSE
            COEF=-DSQRT(-DBLE(IDS6(JT))/DBLE(143))
         ENDIF
      ELSEIF(JTS == 7)THEN
         IF(IDS7(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS7(JT))/DBLE(IDV7(JT)))
         ELSE
            COEF=-DSQRT(-DBLE(IDS7(JT))/DBLE(IDV7(JT)))
         ENDIF
      ELSEIF(JTS == 8)THEN
         IF(IDS8(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS8(JT))/DBLE(IDV8(JT)))
         ELSE
            COEF=-DSQRT(-DBLE(IDS8(JT))/DBLE(IDV8(JT)))
         ENDIF
      ELSEIF(JTS == 9)THEN
         IF(IDS9(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS9(JT))/DBLE(325))
         ELSE
            COEF=-DSQRT(-DBLE(IDS9(JT))/DBLE(325))
         ENDIF
      ELSEIF(JTS == 10)THEN
         IF(IDS10(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS10(JT))/DBLE(165))
         ELSE
            COEF=-DSQRT(-DBLE(IDS10(JT))/DBLE(165))
         ENDIF
      ELSEIF(JTS == 11)THEN
         IF(IDS11(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS11(JT))/DBLE(IDV11(JT)))
         ELSE
            COEF=-DSQRT(-DBLE(IDS11(JT))/DBLE(IDV11(JT)))
         ENDIF
      ELSEIF(JTS == 12)THEN
         IF(IDS12(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS12(JT))/DBLE(IDV12(JT)))
         ELSE
            COEF=-DSQRT(-DBLE(IDS12(JT))/DBLE(IDV12(JT)))
         ENDIF
      ELSEIF(JTS == 13)THEN
         IF(IDS13(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS13(JT))/DBLE(IDV13(JT)))
         ELSE
            COEF=-DSQRT(-DBLE(IDS13(JT))/DBLE(IDV13(JT)))
         ENDIF
      ELSEIF(JTS == 14)THEN
         IF(IDS14(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS14(JT))/DBLE(715))
         ELSE
            COEF=-DSQRT(-DBLE(IDS14(JT))/DBLE(715))
         ENDIF
      ELSEIF(JTS == 15)THEN
         IF(IDS15(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS15(JT))/DBLE(IDV15(JT)))
         ELSE
            COEF=-DSQRT(-DBLE(IDS15(JT))/DBLE(IDV15(JT)))
         ENDIF
      ELSEIF(JTS == 16)THEN
         IF(IDS16(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS16(JT))/DBLE(IDV16(JT)))
         ELSE
            COEF=-DSQRT(-DBLE(IDS16(JT))/DBLE(IDV16(JT)))
         ENDIF
      ELSEIF(JTS == 17)THEN
         IF(IDS17(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS17(JT))/DBLE(IDV17(JT)))
         ELSE
            COEF=-DSQRT(-DBLE(IDS17(JT))/DBLE(IDV17(JT)))
         ENDIF
      ELSEIF(JTS == 18)THEN
         IF(IDS18(JT) >= 0) THEN
            COEF=DSQRT(DBLE(IDS18(JT))/DBLE(IDV18(JT)))
         ELSE
            COEF=-DSQRT(-DBLE(IDS18(JT))/DBLE(IDV18(JT)))
         ENDIF
      ELSE
         WRITE(0,'(A,4I5)') ' IT ITS JT JTS= ',IT,ITS,JT,JTS
         WRITE(0,'(A)') ' ERROR IN SUB. RMEAJJ9 '
         STOP
      ENDIF
      RETURN
    END SUBROUTINE RMEAJJ9

  end subroutine RMEAJJ



  SUBROUTINE RWJJ(J,J1,J2,K1,K2,COEF)
    !*******************************************************************
    !   Written by  G. Gaigalas                                        *
    !   The last modification made by G. Gaigalas       October  2017  *
    !   Restructured by A. Senchuk                     September 2019  *
    !*******************************************************************

    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)  :: J, J1, J2, K1, K2
    REAL(KIND=8), INTENT(OUT) :: COEF
    !-----------------------------------------------
    IF(J == 1) THEN
       CALL RMEW1JJ(J1,J2,K1,K2,COEF)
    ELSEIF(J == 3) THEN
       CALL RMEW3JJ(J1,J2,K1,K2,COEF)
    ELSEIF(J == 5) THEN
       CALL RMEW5JJ(J1,J2,K1,K2,COEF)
    ELSEIF(J == 7) THEN
       CALL RMEW7JJ(J1,J2,K1,K2,COEF)
       !GG        ELSEIF(J.EQ.9) THEN
       !GG          CALL SUWJJ(K1,K2,J,J1,J2,COEF)
    ELSE
       WRITE(0,'(A,I5)') ' KLAIDA SUB. RWJJ J=',J
       STOP
    ENDIF
    RETURN

  CONTAINS


    SUBROUTINE RMEW1JJ(J1,J2,K1,K2,COEF)
      !
      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER,      INTENT(IN)  :: J1, J2, K1, K2
      REAL(KIND=8), INTENT(OUT) :: COEF
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER, DIMENSION(2) :: I10, I01
      !-----------------------------------------------
      DATA I01/6,0/
      DATA I10/0,6/
      !
      COEF=ZERO
      IF(IMPTJJ(J1) /= IMPTJJ(J2)) RETURN
      IF(J1 > 2) RETURN
      IF(K1 == 0 .AND. K2 == 0) THEN
         COEF=-DSQRT(DBLE(2))
      ELSEIF(K1 == 1 .AND. K2 == 0) THEN
         COEF=-DSQRT(DBLE(I10(J1)))
      ELSEIF(K1 == 0 .AND. K2 == 1) THEN
         COEF=-DSQRT(DBLE(I01(J1)))
      ELSE
         WRITE(0,'(A,4I5)') ' J1 J2 = ',J1,J2
         WRITE(0,'(A)') ' ERROR IN SUB. RMEW1JJ '
         STOP
      ENDIF
      RETURN
    END SUBROUTINE RMEW1JJ


    SUBROUTINE RMEW3JJ(J1,J2,K1,K2,COEF)

      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER,      INTENT(IN)  :: J1, J2, K1, K2
      REAL(KIND=8), INTENT(OUT) :: COEF
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER               :: JI1
      INTEGER, DIMENSION(3) :: I00, I10, I01
      !-----------------------------------------------
      DATA I00/16,6,10/
      DATA I10/12,12,0/
      DATA I01/12,0,12/
      !
      COEF=ZERO
      IF(IMPTJJ(J1) /= IMPTJJ(J2)) RETURN
      IF(J1 < 2 .OR. J1 > 5) RETURN
      JI1=J1-2
      IF(K1 == 0 .AND. K2 == 0) THEN
         IF(J1 /= J2) RETURN
         COEF=-DSQRT(DBLE(I00(JI1)))
      ELSEIF(K1 == 1 .AND. K2 == 0) THEN
         IF(J1 /= J2) RETURN
         COEF=-DSQRT(DBLE(I10(JI1)))
      ELSEIF(K1 == 0 .AND. K2 == 1) THEN
         IF(J1 /= J2) RETURN
         COEF=-DSQRT(DBLE(I01(JI1)))
      ELSEIF(K1 == 1 .AND. K2 == 2) THEN
         IF(J1 == 3 .AND. J2 == 3) THEN
            COEF=DSQRT(DBLE(60))
         ELSEIF(J1 == 5 .AND. J2 == 4) THEN
            COEF=DSQRT(DBLE(30))    
         ELSEIF(J1 == 4 .AND. J2 == 5) THEN
            COEF=-DSQRT(DBLE(30))
         ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 3) THEN
         IF(J1 == 3 .AND. J2 == 3) THEN
            COEF=-DSQRT(DBLE(28))
         ELSEIF(J1 == 5 .AND. J2 == 5) THEN
            COEF=DSQRT(DBLE(28))
         ENDIF
      ELSE
         WRITE(0,'(A,4I5)') ' J1 J2 = ',J1,J2
         WRITE(0,'(A)') ' ERROR IN SUB. RMEW3JJ '
         STOP
      ENDIF
      RETURN
    END SUBROUTINE RMEW3JJ

    !*******************************************************************
    !                                                                  *
    SUBROUTINE RMEW5JJ(J1,J2,K1,K2,COEF)
      !                                                                  *
      !   Written by  G. Gaigalas                                        *
      !   The last modification made by G. Gaigalas       October  2017  *
      !                                                                  *
      !*******************************************************************
      !
      implicit none
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER,      INTENT(IN)  :: J1, J2, K1, K2
      REAL(KIND=8), INTENT(OUT) :: COEF
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER                 :: JI1, JI2
      INTEGER, DIMENSION(6)   :: I00, I01
      INTEGER, DIMENSION(3,3) :: I12P,I12A,I03P, &
           I03A,I14P,I14A,I05P,I05A
      !-----------------------------------------------
      DATA I00/54,12,30,12,30,54/
      DATA I01/126,12,198,0,48,288/
      DATA I12P/420,360,270,360,2*0,-270,2*0/
      DATA I12A/0,1960,0,-1960,-1000,2430,0,2430,1980/
      DATA I03P/-882,3*0,384,-400,0,400,286/
      DATA I03A/4*0,162,-300,0,-300,22/
      DATA I14P/3780,-720,-4950,-720,2*0,4950,2*0/
      DATA I14A/2*0,3528,0,2430,1980,-3528,1980,-7722/
      DATA I05P/-1386,4*0,-440,0,440,-1430/
      DATA I05A/5*0,330,0,330,572/
      !
      COEF=ZERO
      IF(IMPTJJ(J1) /= IMPTJJ(J2)) RETURN
      IF(J1 < 6 .OR. J1 > 11) RETURN
      IF(K1 == 0 .AND. K2 == 0) THEN
         IF(J1 /= J2) RETURN
         COEF=-DSQRT(DBLE(I00(J1-5)))
      ELSEIF(K1 == 1 .AND. K2 == 0) THEN
         IF(J1 /= J2) RETURN
         IF(J1 == 6) THEN
            COEF=-DSQRT(DBLE(48))
         ELSEIF(J1 == 9) THEN
            COEF=-DSQRT(DBLE(20))
         ELSEIF(J1 == 10) THEN
            COEF=-DSQRT(DBLE(10))
         ELSEIF(J1 == 11) THEN
            COEF=-DSQRT(DBLE(18))
         ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 1) THEN
         IF(J1 /= J2) RETURN
         COEF=-DSQRT(DBLE(I01(J1-5))/DBLE(7))
      ELSEIF(K1 == 1 .AND. K2 == 2) THEN
         IF(J1 < 9) THEN
            JI1=J1-5
            JI2=J2-5
            IF(I12P(JI1,JI2) >= 0) THEN
               COEF=DSQRT(DBLE(I12P(JI1,JI2))/DBLE(7))
            ELSE
               COEF=-DSQRT(-DBLE(I12P(JI1,JI2))/DBLE(7))
            ENDIF
         ELSE
            JI1=J1-8
            JI2=J2-8
            IF(I12A(JI1,JI2) >= 0) THEN
               COEF=DSQRT(DBLE(I12A(JI1,JI2))/DBLE(49))
            ELSE
               COEF=-DSQRT(-DBLE(I12A(JI1,JI2))/DBLE(49))
            ENDIF
         ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 3) THEN
         IF(J1 < 9) THEN
            JI1=J1-5
            JI2=J2-5
            IF(I03P(JI1,JI2) >= 0) THEN
               COEF=DSQRT(DBLE(I03P(JI1,JI2))/DBLE(21))
            ELSE
               COEF=-DSQRT(-DBLE(I03P(JI1,JI2))/DBLE(21))
            ENDIF
         ELSE
            JI1=J1-8
            JI2=J2-8
            IF(I03A(JI1,JI2) >= 0) THEN
               COEF=DSQRT(DBLE(I03A(JI1,JI2))/DBLE(7))
            ELSE
               COEF=-DSQRT(-DBLE(I03A(JI1,JI2))/DBLE(7))
            ENDIF
         ENDIF
      ELSEIF(K1 == 1 .AND. K2 == 4) THEN
         IF(J1 < 9) THEN
            JI1=J1-5
            JI2=J2-5
            IF(I14P(JI1,JI2) >= 0) THEN
               COEF=DSQRT(DBLE(I14P(JI1,JI2))/DBLE(35))
            ELSE
               COEF=-DSQRT(-DBLE(I14P(JI1,JI2))/DBLE(35))
            ENDIF
         ELSE
            JI1=J1-8
            JI2=J2-8
            IF(I14A(JI1,JI2) >= 0) THEN
               COEF=DSQRT(DBLE(I14A(JI1,JI2))/DBLE(49))
            ELSE
               COEF=-DSQRT(-DBLE(I14A(JI1,JI2))/DBLE(49))
            ENDIF
         ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 5) THEN
         IF(J1 < 9) THEN
            JI1=J1-5
            JI2=J2-5
            IF(I05P(JI1,JI2) >= 0) THEN
               COEF=DSQRT(DBLE(I05P(JI1,JI2))/DBLE(21))
            ELSE
               COEF=-DSQRT(-DBLE(I05P(JI1,JI2))/DBLE(21))
            ENDIF
         ELSE
            JI1=J1-8
            JI2=J2-8
            IF(I05A(JI1,JI2) >= 0) THEN
               COEF=DSQRT(DBLE(I05A(JI1,JI2))/DBLE(7))
            ELSE
               COEF=-DSQRT(-DBLE(I05A(JI1,JI2))/DBLE(7))
            ENDIF
         ENDIF
      ELSE
         WRITE(0,'(A,4I5)') ' J1 J2 = ',J1,J2
         WRITE(0,'(A)') ' ERROR IN SUB. RMEW5JJ '
         STOP
      ENDIF
      RETURN
    END SUBROUTINE RMEW5JJ


    SUBROUTINE RMEW7JJ(J1,J2,K1,K2,COEF)
      !
      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER,      INTENT(IN)  :: J1, J2, K1, K2
      REAL(KIND=8), INTENT(OUT) :: COEF
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER, DIMENSION(14) :: I00, I10, I01
      !-----------------------------------------------
      DATA I00/32,48,128,80,96,128,20,60,20,108,36,44,156,68/
      DATA I10/6,9,120,15,18,24,2*30,0,54,2*0,78,0/
      DATA I01/10,35,168,165,286,680,0,30,10,180,60,110,546,408/
      COEF=ZERO
      IF(IMPTJJ(J1) /= IMPTJJ(J2)) RETURN
      IF(J1 < 12 .OR. J1 > 25) RETURN
      IF(K1 == 0 .AND. K2 == 0) THEN
         IF(J1 /= J2) RETURN
         COEF=-DSQRT(DBLE(I00(J1-11)))
      ELSEIF(K1 == 1 .AND. K2 == 0) THEN
         IF(J1 /= J2) RETURN
         COEF=-DSQRT(DBLE(I10(J1-11)))
      ELSEIF(K1 == 0 .AND. K2 == 1) THEN
         IF(J1 /= J2) RETURN
         COEF=-DSQRT(DBLE(I01(J1-11))/DBLE(7))
      ELSE
         CALL RMEW7BJJ(J1,J2,K1,K2,COEF)
      ENDIF
      RETURN
    END SUBROUTINE RMEW7JJ


    SUBROUTINE RMEW7BJJ(J1,J2,K1,K2,COEF)
      !
      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER,      INTENT(IN)  :: J1, J2, K1, K2
      REAL(KIND=8), INTENT(OUT) :: COEF
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER :: LL, IQ1, LV1, IL1, IQ2, LV2, IL2, JI1, JI2, IFAZ, J, L
      INTEGER, DIMENSION(6)  :: IPR1
      INTEGER, DIMENSION(8)  :: IPR2
      INTEGER, DIMENSION(21) :: I12PS,I12PV,I03PS,I03PV,I14PS, &
           I14PV,I05PS,I05PV,I16PS,I16PV,I07PS,I07PV
      INTEGER, DIMENSION(36) :: I12AS,I12AV,I03AS,I03AV,I14AS, &
           I14AV,I05AS,I05AV,I16AS,I16AV,I07AS,I07AV
      !-----------------------------------------------
      DATA IPR1/0,5,9,12,14,15/
      DATA IPR2/0,7,13,18,22,25,27,28/
      DATA I12PS/-252,1056,144,3*0,507,-88,-325,2*0,200,-520,80,  &
           0,-6125,-560,0,390,-3840,2040/
      DATA I12PV/10,70,7,3*1,10,1,14,2*1,3,21,2*1,462,11,1,539,   &
           2*49/
      DATA I12AS/0,-50,6*0,-1280,990,2640,1950,4*0,-480,4*0,-360, &
           1872,42,390,3*0,48,2*0,234,0,1040,340,0/
      DATA I12AV/8*1,4*49,4*1,49,4*1,539,49,1,11,6*1,7,1,11,7,1/
      DATA I03PS/-1188,-196,0,-234,2*0,189,0,-1911,1470,0,-56,3*0,&
           394805,5250,53760,-78,17408,12920/
      DATA I03PV/70,10,1,7,2*1,110,1,242,121,5*1,22022,121,1573,  &
           2*847,1001/
      DATA I03AS/8*0,110,0,-240,4*0,-1920,0,52,-224,2*0,32490,2*0,&
           -7644,0,12,224,2*0,-52,0,-2040,-364,0,1292/
      DATA I03AV/8*1,7,1,7,4*1,77,1,7,11,2*1,847,2*1,121,1,70,10, &
           2*1,70,1,77,121,1,77/
      DATA I14PS/0,54,-528,546,-378,0,-4335,-96,11271,3822,0,120, &
           12000,-624,-960,-30345,-210,-228480,580476,146880,-627912/
      DATA I14PV/1,2*7,2*11,1,110,11,1694,121,2*1,77,2*11,3146,   &
           121,1573,2*5929,7007/
      DATA I14AS/3*0,-90,4*0,2640,480,-360,20592,-42,390,2*0,     &
           468180,2*0,21840,0,-359424,-6750,0,917280,-10710,2*0,36,2*0,&
           -858,0,-5304,-69768,0/
      DATA I14AV/8*1,2*49,2*539,1,11,2*1,5929,2*1,847,1,5929,539, &
           1,5929,121,2*1,11,2*1,7,1,121,847,1/
      DATA I05PS/3*0,144,14,0,975,0,-50421,336,-7000,-88,3*0,     &
           -6845,-3360,14280,-1836,-103360,-28424/
      DATA I05PV/3*1,7,2*1,22,1,2002,11,143,4*1,26026,143,1859,77,&
           1001,143143/
      DATA I05AS/10*0,390,2*0,-70,3*0,32,14,2*0,-576,2*0,-210,0,  &
           63888,-154,0,-8568,-176,0,1938,1088,0,-28424/
      DATA I05AV/10*1,7,6*1,7,3*1,77,2*1,11,1,1183,13,1,169,7,1,  &
           91,11,1,1183/
      DATA I16PS/3*0,576,390,-144,0,520,-49,-12480,-408,520,      &
           -1960,-1664,3264,552250,-43520,-38760,-2652,15504,38760/
      DATA I16PV/3*1,11,77,7,1,11,2*121,11,3,33,2*11,4719,847,    &
           11011,2*121,143/
      DATA I16AS/6*0,-130,3*0,390,-48,-234,1040,340,0,3120,2*0,   &
           -1170,0,18720,-36,858,-5304,-69768,2*0,1020,4*0,-33592,     &
           31654,0/
      DATA I16AV/10*1,11,1,7,11,7,1,121,2*1,121,1,121,11,7,121,   &
           847,2*1,11,4*1,2*121,1/
      DATA I07PS/4*0,162,272,2*0,11025,4624,-1632,-120,3*0,       &
           1224510,306000,12558240,-6460,77520,-297160/
      DATA I07PV/4*1,2*7,2*1,1573,121,143,4*1,20449,11011,143143, &
           121,1573,1859/
      DATA I07AS/13*0,60,4*0,1600,0,2040,4410,2*0,18360,0,-18816, &
           -11016,0,11628,34,0,7752,9690,0,222870/
      DATA I07AV/18*1,77,1,77,121,2*1,121,1,845,455,1,1183,5,1,   &
           143,121,1,1859/
      !
      COEF = ZERO
      LL = 7
      CALL RUMTJJ(J1,LL,IQ1,LV1,IL1)
      CALL RUMTJJ(J2,LL,IQ2,LV2,IL2)
      IF(J1 > J2) THEN
         JI1=J2
         JI2=J1
         IFAZ=IL2-IL1+IQ2-IQ1
      ELSE
         JI1=J1
         JI2=J2
         IFAZ=4
      ENDIF
      IF(J1 > 17) THEN
         JI1=JI1-17
         JI2=JI2-17
         J=IPR2(JI1)+JI2
         L=2
      ELSE 
         JI1=JI1-11
         JI2=JI2-11
         L=1
         J=IPR1(JI1)+JI2
      ENDIF
      IF(K1 == 1 .AND. K2 == 2) THEN
         IF(L == 1) THEN
            IF(I12PS(J) >= 0) THEN
               COEF=DSQRT(DBLE(I12PS(J))/DBLE(I12PV(J)))
            ELSE
               COEF=-DSQRT(-DBLE(I12PS(J))/DBLE(I12PV(J)))
            ENDIF
         ELSE
            IF(I12AS(J) >= 0) THEN
               COEF=DSQRT(DBLE(I12AS(J))/DBLE(I12AV(J)))
            ELSE
               COEF=-DSQRT(-DBLE(I12AS(J))/DBLE(I12AV(J)))
            ENDIF
         ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 3) THEN
         IF(L == 1) THEN
            IF(I03PS(J) >= 0) THEN
               COEF=DSQRT(DBLE(I03PS(J))/DBLE(I03PV(J)))
            ELSE
               COEF=-DSQRT(-DBLE(I03PS(J))/DBLE(I03PV(J)))
            ENDIF
         ELSE
            IF(I03AS(J) >= 0) THEN
               COEF=DSQRT(DBLE(I03AS(J))/DBLE(I03AV(J)))
            ELSE
               COEF=-DSQRT(-DBLE(I03AS(J))/DBLE(I03AV(J)))
            ENDIF
         ENDIF
      ELSEIF(K1 == 1 .AND. K2 == 4) THEN
         IF(L == 1) THEN
            IF(I14PS(J) >= 0) THEN
               COEF=DSQRT(DBLE(I14PS(J))/DBLE(I14PV(J)))
            ELSE
               COEF=-DSQRT(-DBLE(I14PS(J))/DBLE(I14PV(J)))
            ENDIF
         ELSE
            IF(I14AS(J) >= 0) THEN
               COEF=DSQRT(DBLE(I14AS(J))/DBLE(I14AV(J)))
            ELSE
               COEF=-DSQRT(-DBLE(I14AS(J))/DBLE(I14AV(J)))
            ENDIF
         ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 5) THEN
         IF(L == 1) THEN
            IF(I05PS(J) >= 0) THEN
               COEF=DSQRT(DBLE(I05PS(J))/DBLE(I05PV(J)))
            ELSE
               COEF=-DSQRT(-DBLE(I05PS(J))/DBLE(I05PV(J)))
            ENDIF
         ELSE
            IF(I05AS(J) >= 0) THEN
               COEF=DSQRT(DBLE(I05AS(J))/DBLE(I05AV(J)))
            ELSE
               COEF=-DSQRT(-DBLE(I05AS(J))/DBLE(I05AV(J)))
            ENDIF
         ENDIF
      ELSEIF(K1 == 1 .AND. K2 == 6) THEN
         IF(L == 1) THEN
            IF(I16PS(J) >= 0) THEN
               COEF=DSQRT(DBLE(I16PS(J))/DBLE(I16PV(J)))
            ELSE
               COEF=-DSQRT(-DBLE(I16PS(J))/DBLE(I16PV(J)))
            ENDIF
         ELSE
            IF(I16AS(J) >= 0) THEN
               COEF=DSQRT(DBLE(I16AS(J))/DBLE(I16AV(J)))
            ELSE
               COEF=-DSQRT(-DBLE(I16AS(J))/DBLE(I16AV(J)))
            ENDIF
         ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 7) THEN
         IF(L == 1) THEN
            IF(I07PS(J) >= 0) THEN
               COEF=DSQRT(DBLE(I07PS(J))/DBLE(I07PV(J)))
            ELSE
               COEF=-DSQRT(-DBLE(I07PS(J))/DBLE(I07PV(J)))
            ENDIF
         ELSE
            IF(I07AS(J) >= 0) THEN
               COEF=DSQRT(DBLE(I07AS(J))/DBLE(I07AV(J)))
            ELSE
               COEF=-DSQRT(-DBLE(I07AS(J))/DBLE(I07AV(J)))
            ENDIF
         ENDIF
      ELSE
         WRITE(0,'(A,4I5)') ' J1 J2 = ',J1,J2
         WRITE(0,'(A)') ' ERROR IN SUB. RMEW7JJ '
         STOP
      ENDIF
      IF(MOD(IFAZ,4) /= 0)COEF=-COEF
      RETURN
    END SUBROUTINE RMEW7BJJ

  END SUBROUTINE RWJJ


  SUBROUTINE MES(I)
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER, INTENT(IN) :: I
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: J
    CHARACTER(LEN=10), DIMENSION(6) :: STRING5
    !
    DATA STRING5/'   I T L S',' I T L S 2',' I T L S 3', &
         'AW P 1 L S','WA P 1 L S','       W 1'/
    !-----------------------------------------------
    IF(I > 50) THEN
       J=I-50
       WRITE(6,'(A)') ' error in func./sub. ' 
       WRITE(6,'(20X,A10)') STRING5(J)
       WRITE(6,'(A)') ' susimaise f sluoksnio termu kodavimas  '
    ELSE
       WRITE(6,'(A)') ' yra daugiau nei 2 ele. sluoks. f,g,h,i,k,l,m'
       WRITE(6,'(3X,I5)') I
       IF(I == 1) THEN
          WRITE(6,'(A)') ' error in Subroutine   W 1 G '
       ELSEIF(I == 2) THEN
          WRITE(6,'(A)') ' error in Subroutine   W 1 '
       ELSEIF(I == 11) THEN
          WRITE(6,'(A)') ' error in Function     I T L S  '
       ELSEIF(I == 12) THEN
          WRITE(6,'(A)') ' error in Function     I T L S 2  '
       ELSEIF(I == 13) THEN
          WRITE(6,'(A)') ' error in Function     I T L S 3  '
       ELSEIF(I == 30) THEN
          WRITE(6,'(A)') ' error in Subroutine   A 1 A 2 A 3 A 4 L S '
       ELSEIF(I == 31) THEN
          WRITE(6,'(A)') ' error in Subroutine   A 1 A 2 L S '
       ELSEIF(I == 32) THEN
          WRITE(6,'(A)') ' error in Subroutine   A 1 A 2 W 3 L S '
       ELSEIF(I == 33) THEN
          WRITE(6,'(A)') ' error in Subroutine   A 1 A W 2 L S '
       ELSEIF(I == 34) THEN
          WRITE(6,'(A)') ' error in Subroutine   W A 1 A 2 L S '
       ELSEIF(I == 35) THEN
          WRITE(6,'(A)') ' error in Subroutine   W 1 W 2 L S '
       ELSE
          WRITE(6,'(A)') ' error in unknown Subroutine  '
       ENDIF
    ENDIF
    WRITE(6,'(A)') ' Contact to   G. Gediminas please ! ! ! '
    STOP
  END SUBROUTINE MES

END MODULE irjj
