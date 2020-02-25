MODULE redme
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
  
  USE constants
  USE wignerj
  USE mtjj
  USE irjj

  IMPLICIT NONE



CONTAINS


  SUBROUTINE A1JJ(IK,BK,ID,BD,QM1,A)
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
    REAL(kind=8), INTENT(IN), DIMENSION(3) :: BK, BD
    REAL(kind=8), INTENT(IN)               :: QM1
    REAL(kind=8), INTENT(OUT)              :: A
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER      :: IFAZ, ISUMA
    REAL(kind=8) :: AB
    !-----------------------------------------------
    A=ZERO
    IF(QM1 < EPS) THEN
       ISUMA=(ID(6)+1)*ID(4)
       AB=DBLE(ISUMA)
       A=DSQRT(AB)
       IFAZ=ID(6)+ID(3)-IK(6)+ID(4)*2
       IF((IFAZ/4)*4 /= IFAZ)A=-A
    ELSE
       ISUMA=(IK(6)+1)*IK(4)
       AB=DBLE(ISUMA)
       A=DSQRT(AB)
       IF((IK(4)/2)*2 /= IK(4))A=-A
    ENDIF
    RETURN
  END SUBROUTINE A1JJ


  SUBROUTINE AWP1(IK,BK,ID,BD,K1,BK2,QM1,QM2,QM3,AW)
    !*******************************************************************
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
    !                                                                  *
    !                    N      (j)  (k1) (k2)   N'     +-             *
    !     ELEMENT:     (j QJ ::[A  * W    ]   ::j  QJ)  -+             *
    !                                                   ++             *
    !                                                   --  B17 (2.3)  *
    !                                                                  *
    !*******************************************************************
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)               :: K1
    INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
    REAL(kind=8), INTENT(IN)               :: BK2, QM1, QM2, QM3
    REAL(kind=8), INTENT(IN), DIMENSION(3) :: BK, BD
    REAL(kind=8), INTENT(OUT)              :: AW
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: KK1, KK2, IE, IQ, IQ2, IQ3, IQM, IT, ITP, ITG, IBTT
    INTEGER,      DIMENSION(7) :: IBT
    REAL(kind=8)               :: ENQP, D1, S, SI, W
    REAL(kind=8), DIMENSION(3) :: BT
    !-----------------------------------------------
    AW=ZERO
    IF(ID(3) == 9) THEN
       IF(MAX0(IK(4),ID(4)) < 3) THEN
          IF(IK(1) < 300) CALL MES(54)
          IF(ID(1) < 300) CALL MES(54)
          CALL AWP1JJG(K1,BK2,QM1,QM2,QM3,IK,BK,ID,BD,AW)
          RETURN
       ELSE
          PRINT*, "ERROR in AWP1"
          STOP
       ENDIF
    ELSEIF(ID(3) > 9) THEN
       CALL AWP1JJG(K1,BK2,QM1,QM2,QM3,IK,BK,ID,BD,AW)
       RETURN
    ENDIF
    IF(IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0)RETURN
    ENQP=ZERO
    IQ2=QM2*TWO+QM2*TENTH
    IQ3=QM3*TWO+QM3*TENTH
    IQ=IQ2+IQ3
    KK1=K1*2
    KK2=BK2+BK2+TENTH*BK2
    IF(ITJJ2(IK,ID,KK2,BK,BD,IBT,BT,ITP,ITG,IQ) == 0)RETURN
    IQM=TWO*DABS(BT(3))+TENTH
    DO IT=ITP,ITG
       CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
       IF(IQM > IBT(7)) CYCLE
       IF(IXJTIK(IK(3),KK1,KK2,ID(6),IK(6),IBT(6)) == 0) CYCLE
       IBT(1)=IT
       BT(2)=DBLE(IBT(6))/TWO
       BT(1)=DBLE(IBT(7))/TWO
       CALL C0T5S(BT(1),BT(3),QM1,BK(1),BK(3),D1)
       IF(DABS(D1) < EPS) CYCLE
       CALL RMEAJJ(IK(3),IK(1),IK(7),IK(6),IBT(1),IBT(7),IBT(6),S)
       IF(DABS(S) < EPS) CYCLE
       CALL WJ1(IBT,BT,ID,BD,K1,QM2,QM3,W)
       D1=D1*W*S
       IF(DABS(D1) < EPS) CYCLE
       CALL SIXJ(IK(3),KK1,KK2,ID(6),IK(6),IBT(6),0,SI)
       D1=D1*SI/DSQRT(DBLE(IK(7)+1))
       ENQP=ENQP+D1
    END DO
    AW=ENQP*DSQRT(DBLE(KK2+1))
    IE=KK2+IK(6)+ID(6)+2
    IF(((IE/4)*4) /= IE)AW=-AW
    RETURN

  CONTAINS

    SUBROUTINE AWP1JJG(K1,BK2,QM1,QM2,QM3,IK,BK,ID,BD,AW)
      !
      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER,      INTENT(IN)               :: K1
      INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
      REAL(KIND=8), INTENT(IN)               :: BK2, QM1, QM2, QM3
      REAL(KIND=8), INTENT(IN), DIMENSION(3) :: BK, BD
      REAL(KIND=8), INTENT(OUT)              :: AW
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER :: KK1, KK2, IE, IQ, IQ2, IQ3, IQM, IT, ITP, ITG, IBTT
      INTEGER,      DIMENSION(7) :: IBT 
      REAL(KIND=8)               :: ENQP, D1, SI1, W
      REAL(KIND=8), DIMENSION(3) :: BT
      !-----------------------------------------------
      AW=ZERO
      IF(IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0)RETURN
      ENQP=ZERO
      KK1=K1*2
      KK2=BK2+BK2+TENTH*BK2
      IQ2=QM2*TWO+QM2*TENTH
      IQ3=QM3*TWO+QM3*TENTH
      IQ=IQ2+IQ2
      IF(ID(3) > 37) RETURN
      IF(ITJJ2(IK,ID,KK2,BK,BD,IBT,BT,ITP,ITG,IQ) == 0)RETURN
      IQM=TWO*DABS(BT(3))+TENTH
      DO IT=ITP,ITG
         CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
         IF(IQM <= IBT(7)) THEN
            IF(IXJTIK(IK(3),KK1,KK2,ID(6),IK(6),IBT(6)) /= 0) THEN
               IBT(1)=IT
               BT(2)=DBLE(IBT(6))/TWO
               BT(1)=DBLE(IBT(7))/TWO
               CALL A1JJ(IK,BK,IBT,BT,QM1,D1)
               IF(DABS(D1) > EPS) THEN
                  CALL W1JJG(K1,QM2,QM3,IBT,BT,ID,BD,W)
                  IF(DABS(W) > EPS) THEN
                     D1=D1*W
                     CALL SIXJ(IK(3),KK1,KK2,ID(6),IK(6),IBT(6),0,SI1)
                     ENQP=ENQP+D1*SI1
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      END DO
      AW=ENQP*DSQRT(DBLE((KK2+1)))
      IE=KK2+IK(6)+ID(6)
      IF(((IE/4)*4) /= IE)AW=-AW
      RETURN
    END SUBROUTINE AWP1JJG

  END SUBROUTINE AWP1

  
  INTEGER FUNCTION IZAS1(IB,QB,IK,QK)
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN) :: IB, IK
    REAL(kind=8), INTENT(IN) :: QB, QK
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: IQB, IQK
    !-----------------------------------------------
    IZAS1=0
    IQB=TWO*DABS(QB)+TENTH
    IF(IQB > IB)RETURN
    IF(MOD(IB+IQB,2) /= 0)RETURN
    IQK=TWO*DABS(QK)+TENTH
    IF(IQK > IK)RETURN
    IF(MOD(IK+IQK,2) /= 0)RETURN
    IZAS1=1
    RETURN
  END  FUNCTION IZAS1

  
  SUBROUTINE WAP1(IK,BK,ID,BD,K1,BK2,QM1,QM2,QM3,WA)
    !*******************************************************************
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
    !                                                                  *
    !                    N       (k1)   (j) (k2)   N'    +-            *
    !     ELEMENT:     (j QJ::[ W    * A   ]    ::j QJ)  -+            *
    !                                                    ++            *
    !                                                    -- B17 (2.3)  *
    !                                                                  *
    !*******************************************************************
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)               :: K1
    INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
    REAL(KIND=8), INTENT(IN)               :: BK2, QM1, QM2, QM3
    REAL(KIND=8), INTENT(IN), DIMENSION(3) :: BK, BD
    REAL(KIND=8), INTENT(OUT)              :: WA
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: KK1, KK2, IE, IQ3, IQM, IT, ITP, ITG, IBTT
    INTEGER,      DIMENSION(7) :: IBT
    REAL(KIND=8)               :: ENQP, D1, S, SI, W
    REAL(KIND=8), DIMENSION(3) :: BT
    !-----------------------------------------------
    WA=ZERO
    IF(ID(3) == 9) THEN
       IF(MAX0(IK(4),ID(4)) < 3) THEN
          IF(IK(1) < 300) CALL MES(55)
          IF(ID(1) < 300) CALL MES(55)
          CALL WAP1JJG(K1,BK2,QM1,QM2,QM3,IK,BK,ID,BD,WA)
          RETURN
       ENDIF
    ELSEIF(ID(3) > 9) THEN
       CALL WAP1JJG(K1,BK2,QM1,QM2,QM3,IK,BK,ID,BD,WA)
       RETURN
    ENDIF
    IF(IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0)RETURN
    ENQP=ZERO
    KK1=K1*2
    IQ3=QM3*TWO
    KK2=BK2+BK2+TENTH*BK2
    IF(ITJJ3(IK,ID,KK2,BK,BD,IBT,BT,ITP,ITG,IQ3) == 0)RETURN
    IQM=TWO*DABS(BT(3))+TENTH
    DO IT=ITP,ITG
       CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
       IF(IQM > IBT(7)) CYCLE
       IF(IXJTIK(KK1,IK(3),KK2,ID(6),IK(6),IBT(6)) == 0) CYCLE
       IBT(1)=IT
       BT(2)=DBLE(IBT(6))/TWO
       BT(1)=DBLE(IBT(7))/TWO
       CALL C0T5S(BD(1),BD(3),QM3,BT(1),BT(3),D1)
       IF(DABS(D1) < EPS) CYCLE
       CALL RMEAJJ(ID(3),IBT(1),IBT(7),IBT(6),ID(1),ID(7),ID(6),S)
       IF(DABS(S) < EPS) CYCLE
       D1=D1*S
       CALL WJ1(IK,BK,IBT,BT,K1,QM1,QM2,W)
       IF(DABS(W) < EPS) CYCLE
       D1=D1*W
       CALL SIXJ(KK1,IK(3),KK2,ID(6),IK(6),IBT(6),0,SI)
       D1=D1*SI/DSQRT(DBLE(IBT(7)+1))
       ENQP=ENQP+D1
    END DO
    WA=ENQP*DSQRT(DBLE(KK2+1))
    IE=KK2+IK(6)+ID(6)+2
    IF(((IE/4)*4) /= IE)WA=-WA
    RETURN
  END SUBROUTINE WAP1


  SUBROUTINE WAP1JJG(K1,BK2,QM1,QM2,QM3,IK,BK,ID,BD,WA)
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)               :: K1
    INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
    REAL(kind=8), INTENT(IN)               :: BK2, QM1, QM2, QM3
    REAL(kind=8), INTENT(IN), DIMENSION(3) :: BK, BD
    REAL(kind=8), INTENT(OUT)              :: WA
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: KK1, KK2, IE, IQ3, IQM, IT, ITP, ITG, IBTT
    INTEGER,      DIMENSION(7) :: IBT
    REAL(kind=8)               :: ENQP, D1, SI1, W
    REAL(kind=8), DIMENSION(3) :: BT
    !-----------------------------------------------
    WA=ZERO
    IF(IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0)RETURN
    ENQP=ZERO
    KK1=K1*2
    KK2=BK2+BK2+TENTH*BK2
    IQ3=QM3*TWO+QM3*TENTH
    IF(ID(3) > 37) RETURN
    IF(ITJJ3(IK,ID,KK2,BK,BD,IBT,BT,ITP,ITG,IQ3) == 0)RETURN
    IQM=TWO*DABS(BT(3))+TENTH
    DO IT=ITP,ITG
       CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
       IF(IQM > IBT(7)) CYCLE
       IF(IXJTIK(KK1,IK(3),KK2,ID(6),IK(6),IBT(6)) == 0) CYCLE
       IBT(1)=IT
       BT(2)=DBLE(IBT(6))/TWO
       BT(1)=DBLE(IBT(7))/TWO
       CALL A1JJ(IBT,BT,ID,BD,QM3,D1)
       IF(DABS(D1) < EPS) CYCLE
       CALL W1JJG(K1,QM1,QM2,IK,BK,IBT,BT,W)
       IF(DABS(W) < EPS) CYCLE
       D1=D1*W
       CALL SIXJ(KK1,IK(3),KK2,ID(6),IK(6),IBT(6),0,SI1)
       ENQP=ENQP+D1*SI1
    END DO
    WA=ENQP*DSQRT(DBLE((KK2+1)))
    IE=KK2+IK(6)+ID(6)
    IF(((IE/4)*4) /= IE)WA=-WA
    RETURN
  END SUBROUTINE WAP1JJG


  SUBROUTINE WJ1(IK,BK,ID,BD,K2,QM1,QM2,WJ)
    !*******************************************************************
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
    !                                                                  *
    !                      N       (k2)  N'     +-                     *
    !     ELEMENT:       (j  QJ:: W   ::j  QJ)  -+                     *
    !                                           ++                     *
    !                                           -- S5(1.47),(1.48),    *
    !                                                (1.49),(1.50).    *
    !                                                                  *
    !*******************************************************************
    !
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)               :: K2
    INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
    REAL(kind=8), INTENT(IN)               :: QM1, QM2
    REAL(kind=8), INTENT(IN), DIMENSION(3) :: BK, BD
    REAL(kind=8), INTENT(OUT)              :: WJ
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER      :: K1, IQ, IQM2
    REAL(kind=8) :: A, QQ, W, WK1
    !-----------------------------------------------
    WJ=ZERO
    IF(ID(3) == 9) THEN
       IF(MAX0(IK(4),ID(4)) < 3) THEN
          IF(IK(1) < 300) CALL MES(56)
          IF(ID(1) < 300) CALL MES(56)
          IQM2=QM2+QM2+QM2*EPS
          IF((ID(4)+IQM2) > 2) CALL MES(2)
          CALL W1JJG(K2,QM1,QM2,IK,BK,ID,BD,WJ)
          RETURN
       ELSE
          PRINT*, "ERROR in  WJ1"
          STOP
       ENDIF
    ELSEIF(ID(3) > 9) THEN
       IQM2=QM2+QM2+QM2*EPS
       IF((ID(4)+IQM2) > 2) CALL MES(2)
       CALL W1JJG(K2,QM1,QM2,IK,BK,ID,BD,WJ)
       RETURN
    ENDIF
    QQ=QM1+QM2
    IF(DABS(QQ) >= EPS) THEN
       IF(((K2/2)*2) /= K2)RETURN
       IQ=QQ+QQ*TENTH
       IF((IK(4)-ID(4)-2*IQ) /= 0)RETURN
       CALL C1E1SM(BD(1),BD(3),QQ,BK(1),BK(3),A)
       IF(DABS(A) < EPS)RETURN
       CALL RWJJ(IK(3),IK(1),ID(1),1,K2,W)
       WJ=A*W/DSQRT(TWO*BK(1)+ONE)
    ELSE IF(IK(4) /= ID(4))THEN
       RETURN
    ELSE IF(K2 == 0) THEN
       IF(ID(1) /= IK(1))RETURN
       IF(QM1 < EPS) THEN
          A=DBLE(ID(3)+1-ID(4))
       ELSE
          A=-DBLE(ID(4))
       END IF
       WJ=A*DSQRT(DBLE(IK(6)+1)/DBLE(IK(3)+1))
    ELSE
       K1=1
       IF(((K2/2)*2) /= K2)K1=0
       WK1=DBLE(K1)
       CALL CLE0SM(BD(1),BD(3),WK1,BK(1),BK(3),A)
       IF(DABS(A) < EPS)RETURN
       CALL RWJJ(IK(3),IK(1),ID(1),K1,K2,W)
       A=A*W
       WJ=A/DSQRT(TWO*TWO*BK(1)+TWO)
       IF(QM1 >= EPS)RETURN
       IF(((K2/2)*2) /= K2)WJ=-WJ
    END IF
    RETURN
  END SUBROUTINE WJ1


  SUBROUTINE W1JJG(K1,QM1,QM2,IK,BK,ID,BD,WW)
    !
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)               :: K1
    INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
    REAL(kind=8), INTENT(IN)               :: QM1, QM2
    REAL(kind=8), INTENT(IN), DIMENSION(3) :: BK, BD
    REAL(kind=8), INTENT(OUT)              :: WW
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: KK1,IQ,IQM,IE,IT,ITK,ITD,ITP,ITP1,ITG,ITG1,IBTT,ITTK
    INTEGER,      DIMENSION(7) :: IBT
    REAL(kind=8)               :: ENQP, D1, SI1, W
    REAL(kind=8), DIMENSION(3) :: BT
    !-----------------------------------------------
    WW=ZERO
    IF(IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0)RETURN
    ENQP=ZERO
    KK1=K1*2
    IQ=QM2*TWO+QM2*TENTH
    IF(ID(3) > 37) RETURN
    IF(ITTK(ID(6),IK(6),KK1) == 0)RETURN
    ITK=IK(1)
    ITD=ID(1)
    IF(ID(3) == 9) THEN
       IF(ID(4) > 2) CALL MES(1)
       IF(IK(4) > 2) CALL MES(1)
       ITK=ITK-300
       ITD=ITD-300
       ITP1=IMPNJJ9(ITK)
       ITP=IMPNJJ9(ITD)
       IF(ITP1 /= ITP)RETURN
       ITG1=IMGNJJ9(ITK)
       ITG=IMGNJJ9(ITD)
    ELSE
       IF(ID(4) > 2) CALL MES(1)
       IF(IK(4) > 2) CALL MES(1)
       ITP1=IMPNJJ11(ITK)
       ITP=IMPNJJ11(ITD)
       IF(ITP1 /= ITP)RETURN
       ITG1=IMGNJJ11(ITK)
       ITG=IMGNJJ11(ITD)
    ENDIF
    IF(ITG1 /= ITG)RETURN
    IBT(2)=ID(2)
    IBT(3)=ID(3)
    IBT(4)=ID(4)+IQ
    BT(3)=BD(3)+HALF*DBLE(IQ)
    IQM=TWO*DABS(BT(3))+TENTH
    DO IT=ITP,ITG
       CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
       IF(IQM <= IBT(7)) THEN
          IF(IXJTIK(IK(3),IK(3),KK1,ID(6),IK(6),IBT(6)) /= 0) THEN
             IBT(1)=IT
             BT(2)=DBLE(IBT(6))/TWO
             BT(1)=DBLE(IBT(7))/TWO
             CALL A1JJ(IK,BK,IBT,BT,QM1,D1)
             IF(DABS(D1) > EPS) THEN
                CALL A1JJ(IBT,BT,ID,BD,QM2,W)
                IF(DABS(W) > EPS) THEN
                   D1=D1*W
                   CALL SIXJ(IK(3),IK(3),KK1,ID(6),IK(6),IBT(6),0,SI1)
                   ENQP=ENQP+D1*SI1
                ENDIF
             ENDIF
          ENDIF
       ENDIF
    END DO
    WW=ENQP*DSQRT(DBLE((KK1+1)))
    IE=KK1+IK(6)+ID(6)
    IF(((IE/4)*4) /= IE)WW=-WW
    RETURN
  END SUBROUTINE W1JJG


  SUBROUTINE WW1(IK,BK,ID,BD,K2,QM1,QM2,QM3,QM4,WW)
    !*******************************************************************
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
    !                                                                  *
    !                    N      (k2)   (k2) (0)   N'     +-            *
    !     ELEMENT      (j QJ::[W   *  W    ]   ::j QJ)   -+            *
    !                                                    ++            *
    !                                                    -- B17 (2.4)  *
    !*******************************************************************
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)               :: K2
    INTEGER,      INTENT(IN), DIMENSION(7) :: IK, ID
    REAL(KIND=8), INTENT(IN)               :: QM1, QM2, QM3, QM4
    REAL(KIND=8), INTENT(IN), DIMENSION(3) :: BK, BD
    REAL(KIND=8), INTENT(OUT)              :: WW
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: KK2,KK6,IQ,IQM,IQ3,IQ4,IE,IE1,IT,ITP,ITG,IBTT
    INTEGER, DIMENSION(7)      :: IBT
    REAL(KIND=8)               :: ENQP, D1, W
    REAL(KIND=8), DIMENSION(3) :: BT
    !-----------------------------------------------
    WW=ZERO
    IF(ID(6) /= IK(6))RETURN
    IF(IZAS1(ID(7),BD(3),IK(7),BK(3)) == 0)RETURN
    ENQP=ZERO
    KK2=K2*2
    IQ3=QM3*TWO+QM3*TENTH
    IQ4=QM4*TWO+QM4*TENTH
    IQ=IQ3+IQ4
    IF(ITJJ(IK,ID,0,BK,BD,IBT,BT,KK6,ITP,ITG,IQ) == 0)RETURN
    IE1=KK2-IK(6)
    IQM=TWO*DABS(BT(3))+TENTH
    DO IT=ITP,ITG
       CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
       IF(IQM > IBT(7)) CYCLE
       IF(IXJTIK(KK2,KK2,0,ID(6),IK(6),IBT(6)) == 0) CYCLE
       IBT(1)=IT
       BT(2)=DBLE(IBT(6))/TWO
       BT(1)=DBLE(IBT(7))/TWO
       CALL WJ1(IK,BK,IBT,BT,K2,QM1,QM2,D1)
       IF(DABS(D1) < EPS) CYCLE
       CALL WJ1(IBT,BT,ID,BD,K2,QM3,QM4,W)
       IF(DABS(W) < EPS) CYCLE
       D1=D1*W
       IE=IE1+IBT(6)
       IF(((IE/4)*4) /= IE)D1=-D1
       ENQP=ENQP+D1
    END DO
    WW=ENQP/DSQRT(DBLE(KK2+1)*(IK(6)+1))
    RETURN
  END SUBROUTINE WW1


END MODULE redme
