MODULE gg
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
!   MODULE created by A. Senchuk                   September 2019  *  
!                                                                  *
!*******************************************************************
  USE redme
  
  IMPLICIT NONE

  !> 
  INTEGER, DIMENSION(7)      :: ID1,ID2,IK1,IK2,ID3,ID4,IK3,IK4
  REAL(kind=8), DIMENSION(3) :: BD1,BD2,BK1,BK2,BD3,BD4,BK3,BK4

CONTAINS

  
  SUBROUTINE GG1112(IK1,IK2,BK1,BK2,ID1,ID2,BD1,BD2,K1,  &
       QM1,QM2,QM3,QM4,WW)
    !*******************************************************************
    !                                                                  *
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
    !     ELEMENT:                                                     *
    !                                                                  *
    !       N1       (k1) (j1)(j2)  N1'       N2     (j2)   N2'     +- *
    !     (j Q J ::[W(11)*A(1)]  ::j  Q'J')*(j Q J ::A(2)::j Q'J')  -+ *
    !       1 1 1                   1  1 1    2 2 2         2 2 2   ++ *
    !                                                               -- *
    !                                                                  *
    !*******************************************************************
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)               :: K1
    INTEGER,      INTENT(IN), DIMENSION(7) :: IK1, IK2, ID1, ID2
    REAL(KIND=8), INTENT(IN)               :: QM1, QM2, QM3, QM4
    REAL(KIND=8), INTENT(IN), DIMENSION(3) :: BK1, BK2, BD1, BD2
    REAL(KIND=8), INTENT(OUT)              :: WW
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: IQMM1, IQMM2, IQMM3, IQMM4, IQMM23, KK1
    REAL(KIND=8) :: A1, S, BK, WA
    !-----------------------------------------------
    WW=ZERO
    IF(IK1(3) > 9) THEN
       IF(IK1(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   W A 1 A 2 L S '
       IF(ID1(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   W A 1 A 2 L S '
    ENDIF
    IF(IK2(3) > 9) THEN
       IF(IK2(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   W A 1 A 2 L S '
       IF(ID2(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   W A 1 A 2 L S '
    ENDIF
    IQMM1=QM1+QM1+TENTH*QM1
    IQMM2=QM2+QM2+TENTH*QM2
    IQMM3=QM3+QM3+TENTH*QM3
    IQMM23=IQMM1+IQMM2+IQMM3
    IF(IK1(4) /= (ID1(4)+IQMM23))RETURN
    IQMM4=QM4+QM4+TENTH*QM4
    IF(IK2(4) /= (ID2(4)+IQMM4))RETURN
    KK1=K1*2
    CALL C0T5S(BD2(1),BD2(3),QM4,BK2(1),BK2(3),A1)
    IF(ABS(A1) < EPS)RETURN
    CALL RMEAJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),S)
    IF(ABS(S) < EPS)RETURN
    BK=HALF*DBLE(IK2(3))
    CALL WAP1(IK1,BK1,ID1,BD1,K1,BK,QM1,QM2,QM3,WA)
    IF(ABS(WA) < EPS)RETURN
    WW=-A1*WA*S/SQRT(DBLE(IK2(7)+1))
    RETURN
  END SUBROUTINE GG1112
  
  
  
  SUBROUTINE GG1122(K1,K2,QM1,QM2,QM3,QM4,AA)
    !*******************************************************************
    !                                                                  *    
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
    !     ELEMENTS:                                                    *
    !                                                                  *
    !       N1       (k1)   N1'         N2       (k2)    N2'        +- *
    !     (j  Q J ::W(11)::j  Q'J') * (j  Q J ::W(22):: j  Q'J')    -+ *
    !       1  1 1          1  1 1      2  2 2           2  2 2     ++ *
    !                                                               -- *
    !                                                                  *
    !*******************************************************************
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)               :: K1, K2
    REAL(KIND=8), INTENT(IN)               :: QM1, QM2, QM3, QM4
    REAL(KIND=8), INTENT(OUT)              :: AA
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: IQMM1, IQMM2, IQMM3, IQMM4, IQMM12, IQMM34
    REAL(KIND=8) :: A1, W
    !-----------------------------------------------
    AA=ZERO
    IF(IK1(3) > 9) THEN
       IF(IK1(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   W 1 W 2 L S '
       IF(ID1(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   W 1 W 2 L S '
    ENDIF
    IF(IK2(3) > 9) THEN
       IF(IK2(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   W 1 W 2 L S '
       IF(ID2(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   W 1 W 2 L S '
    ENDIF
    IQMM1=QM1+QM1+TENTH*QM1
    IQMM2=QM2+QM2+TENTH*QM2
    IQMM12=IQMM1+IQMM2
    IF(IK1(4) /= (ID1(4)+IQMM12))RETURN
    IQMM3=QM3+QM3+TENTH*QM3
    IQMM4=QM4+QM4+TENTH*QM4
    IQMM34=IQMM3+IQMM4
    IF(IK2(4) /= (ID2(4)+IQMM34))RETURN
    CALL WJ1(IK1,BK1,ID1,BD1,K1,QM1,QM2,A1)
    IF(DABS(A1) < EPS)RETURN
    CALL WJ1(IK2,BK2,ID2,BD2,K2,QM3,QM4,W)
    IF(DABS(W) < EPS)RETURN
    AA=A1*W
    RETURN
  END SUBROUTINE GG1122
  


  SUBROUTINE GG12(IK1,IK2,BK1,BK2,ID1,ID2,BD1,BD2,QM1,QM2,WW)
    !*******************************************************************    
    !                                                                  *
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
    !                                                                  *
    !                   N1     (j1)  N1'      N2     (j2)  N2'      +- *
    !     ELEMENTS:   (j Q J::A(1)::j Q'J')*(j Q J::A(2)::j Q'J')*  -+ *
    !                   1 1 1        1 1 1    2 2 2        2 2 2    ++ *
    !                                                               -- *
    !                                                                  *
    !*******************************************************************
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN), DIMENSION(7) :: IK1, IK2, ID1, ID2
    REAL(KIND=8), INTENT(IN)               :: QM1, QM2
    REAL(KIND=8), INTENT(IN), DIMENSION(3) :: BK1, BK2, BD1, BD2
    REAL(KIND=8), INTENT(OUT)              :: WW
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: IQMM1, IQMM2
    REAL(KIND=8) :: A1, C, C1, S
    !-----------------------------------------------
    WW=ZERO
    IF(IK1(3) > 9) THEN
       IF(IK1(4) > 2) WRITE(6,'(A)') ' error in Subroutine   A 1 A 2 W 3 L S '
       IF(ID1(4) > 2) WRITE(6,'(A)') ' error in Subroutine   A 1 A 2 W 3 L S '
    ENDIF
    IF(IK2(3) > 9) THEN
       IF(IK2(4) > 2) WRITE(6,'(A)') ' error in Subroutine   A 1 A 2 W 3 L S '
       IF(ID2(4) > 2) WRITE(6,'(A)') ' error in Subroutine   A 1 A 2 W 3 L S '
    ENDIF
    IQMM1=QM1+QM1+TENTH*QM1
    IF(IK1(4) /= (ID1(4)+IQMM1))RETURN
    IQMM2=QM2+QM2+TENTH*QM2
    IF(IK2(4) /= (ID2(4)+IQMM2))RETURN
    CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A1)
    IF(DABS(A1) < EPS)RETURN
    CALL C0T5S(BD2(1),BD2(3),QM2,BK2(1),BK2(3),C1)
    IF(DABS(C1) < EPS)RETURN
    A1=A1*C1
    CALL RMEAJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
    IF(DABS(S) < EPS)RETURN
    CALL RMEAJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),C)
    IF(DABS(C) < EPS)RETURN
    WW=A1*S*C/DSQRT(DBLE((IK1(7)+1)*(IK2(7)+1)))
    RETURN
  END SUBROUTINE GG12



  SUBROUTINE GG1222(IK1,IK2,BK1,BK2,ID1,ID2,BD1,BD2,K1, &
       QM1,QM2,QM3,QM4,WW)
    !*******************************************************************
    !                                                                  *
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
    !                                                                  *
    !                                          N1      (j1)   N1'      *
    !     ELEMENTS:                          (j  Q J ::A(1)::j Q'J')*  *
    !                                          1  1 1         1 1 1    *
    !        N2        (j2)    (k1) (j1)  N2'                       +- *
    !     *(j  Q J ::[ A(2) * W(22) ]  ::j  Q'J')                   -+ *
    !        2  2 2                       2  2 2                    ++ *
    !                                                               -- *
    !                                                                  *
    !*******************************************************************
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)               :: K1
    INTEGER,      INTENT(IN), DIMENSION(7) :: IK1, IK2, ID1, ID2
    REAL(KIND=8), INTENT(IN)               :: QM1, QM2, QM3, QM4
    REAL(KIND=8), INTENT(IN), DIMENSION(3) :: BK1, BK2, BD1, BD2
    REAL(KIND=8), INTENT(OUT)              :: WW
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER      :: IQMM1,IQMM2,IQMM3,IQMM4,IQMM34,KK1
    REAL(KIND=8) :: A1, S, BK, AW
    !-----------------------------------------------
    WW=ZERO
    IF(IK1(3) > 9) THEN
       IF(IK1(4) > 2) WRITE(6,'(A)') ' error in Subroutine   A 1 A W 2 L S '
       IF(ID1(4) > 2) WRITE(6,'(A)') ' error in Subroutine   A 1 A W 2 L S '
    ENDIF
    IF(IK2(3) > 9) THEN
       IF(IK2(4) > 2) WRITE(6,'(A)') ' error in Subroutine   A 1 A W 2 L S '
       IF(ID2(4) > 2) WRITE(6,'(A)') ' error in Subroutine   A 1 A W 2 L S '
    ENDIF
    IQMM1=QM1+QM1+TENTH*QM1
    IF(IK1(4) /= (ID1(4)+IQMM1))RETURN
    IQMM2=QM2+QM2+TENTH*QM2
    IQMM3=QM3+QM3+TENTH*QM3
    IQMM4=QM4+QM4+TENTH*QM4
    IQMM34=IQMM2+IQMM3+IQMM4
    IF(IK2(4) /= (ID2(4)+IQMM34))RETURN
    KK1=K1*2
    CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A1)
    IF(DABS(A1) < EPS)RETURN
    CALL RMEAJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
    IF(DABS(S) < EPS)RETURN
    BK=HALF*DBLE(IK1(3))
    CALL AWP1(IK2,BK2,ID2,BD2,K1,BK,QM2,QM3,QM4,AW)
    IF(DABS(AW) < EPS)RETURN
    WW=-A1*AW*S/DSQRT(DBLE(IK1(7)+1))
    RETURN
  END SUBROUTINE GG1222



  SUBROUTINE GG1233(IK1,IK2,IK3,BK1,BK2,BK3,ID1,ID2,ID3,BD1,  &
       BD2,BD3,K1,QM1,QM2,QM3,QM4,WW)
    !*******************************************************************
    !                                                                  *
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
    !                                                                  *
    !                      N1     (j1)  N1'      N2     (j2)  N2'      *
    !     ELEMENTS:      (j Q J::A(1)::j Q'J')*(j Q J::A(2)::j Q'J')*  *
    !                      1 1 1        1 1 1    2 2 2        2 2 2    *
    !                                                                  *
    !        N3     (k1)   N3'                                      +- *
    !     *(j Q J::W(33)::j Q'J')                                   -+ *
    !        3 3 3         3 3 3                                    ++ *
    !                                                               -- *
    !                                                                  *
    !*******************************************************************
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN)               :: K1
    INTEGER,      INTENT(IN), DIMENSION(7) :: IK1,IK2,IK3,ID1,ID2,ID3
    REAL(KIND=8), INTENT(IN)               :: QM1,QM2,QM3,QM4
    REAL(KIND=8), INTENT(IN), DIMENSION(3) :: BK1,BK2,BK3,BD1,BD2,BD3
    REAL(KIND=8), INTENT(OUT)              :: WW
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER      :: IQMM1,IQMM2,IQMM3,IQMM4,IQMM34
    REAL(KIND=8) :: A1, C, C1, S, W
    !-----------------------------------------------
    WW=ZERO
    IF(IK1(3) > 9) THEN
       IF(IK1(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 W 3 L S '
       IF(ID1(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 W 3 L S '
    ENDIF
    IF(IK2(3) > 9) THEN
       IF(IK2(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 W 3 L S '
       IF(ID2(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 W 3 L S '
    ENDIF
    IF(IK3(3) > 9) THEN
       IF(IK3(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 W 3 L S '
       IF(ID3(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 W 3 L S '
    ENDIF
    IQMM1=QM1+QM1+TENTH*QM1
    IF(IK1(4) /= (ID1(4)+IQMM1))RETURN
    IQMM2=QM2+QM2+TENTH*QM2
    IF(IK2(4) /= (ID2(4)+IQMM2))RETURN
    IQMM3=QM3+QM3+TENTH*QM3
    IQMM4=QM4+QM4+TENTH*QM4
    IQMM34=IQMM3+IQMM4
    IF(IK3(4) /= (ID3(4)+IQMM34))RETURN
    CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A1)
    IF(DABS(A1) < EPS)RETURN
    CALL C0T5S(BD2(1),BD2(3),QM2,BK2(1),BK2(3),C1)
    IF(DABS(C1) < EPS)RETURN
    A1=A1*C1
    CALL RMEAJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
    IF(DABS(S) < EPS)RETURN
    CALL RMEAJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),C)
    IF(DABS(C) < EPS)RETURN
    CALL WJ1(IK3,BK3,ID3,BD3,K1,QM3,QM4,W)
    IF(DABS(W) < EPS)RETURN
    WW=A1*W*S*C/DSQRT(DBLE((IK1(7)+1)*(IK2(7)+1)))
    RETURN
  END SUBROUTINE GG1233



  SUBROUTINE GG1234(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,  &
       ID3,ID4,BD1,BD2,BD3,BD4,QM1,QM2,QM3,QM4,WW)
    !*******************************************************************
    !                                                                  *
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
    !                                                                  *
    !                      N1     (j1)  N1'      N2    (j2)   N2'      *
    !     ELEMENTS:      (j Q J::A(1)::j Q'J')*(j Q J::A(2)::j Q'J')*  *
    !                      1 1 1        1 1 1    2 2 2        2 2 2    *
    !                                                                  *
    !        N3     (j3)  N3'      N4     (j4)  N4'                 +- *
    !     *(j Q J::A(3)::j Q'J')*(j Q J::A(4)::j Q'J')              -+ *
    !        3 3 3        3 3 3    4 4 4        4 4 4               ++ *
    !                                                               -- *
    !                                                                  *
    !*******************************************************************
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER,      INTENT(IN), DIMENSION(7) :: IK1,IK2,IK3,IK4,ID1,ID2,ID3,ID4
    REAL(KIND=8), INTENT(IN)               :: QM1,QM2,QM3,QM4
    REAL(KIND=8), INTENT(IN), DIMENSION(3) :: BK1,BK2,BK3,BK4,BD1,BD2,BD3,BD4
    REAL(KIND=8), INTENT(OUT)              :: WW
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: IQMM1, IQMM2, IQMM3, IQMM4
    REAL(KIND=8) :: A1, C, C1, S, V, Z
    !-----------------------------------------------
    WW=ZERO
    IF(IK1(3) > 9) THEN
       IF(IK1(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 A 3 A 4 L S '
       IF(ID1(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 A 3 A 4 L S '
    ENDIF
    IF(IK2(3) > 9) THEN
       IF(IK2(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 A 3 A 4 L S '
       IF(ID2(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 A 3 A 4 L S '
    ENDIF
    IF(IK3(3) > 9) THEN
       IF(IK3(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 A 3 A 4 L S '
       IF(ID3(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 A 3 A 4 L S '
    ENDIF
    IF(IK4(3) > 9) THEN
       IF(IK4(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 A 3 A 4 L S '
       IF(ID4(4) > 2) WRITE(6,'(A)') ' ERROR IN SUBROUTINE   A 1 A 2 A 3 A 4 L S '
    ENDIF
    IQMM1=QM1+QM1+TENTH*QM1
    IF(IK1(4) /= (ID1(4)+IQMM1))RETURN
    IQMM2=QM2+QM2+TENTH*QM2
    IF(IK2(4) /= (ID2(4)+IQMM2))RETURN
    IQMM3=QM3+QM3+TENTH*QM3
    IF(IK3(4) /= (ID3(4)+IQMM3))RETURN
    IQMM4=QM4+QM4+TENTH*QM4
    IF(IK4(4) /= (ID4(4)+IQMM4))RETURN
    CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A1)
    IF(DABS(A1) < EPS)RETURN
    CALL C0T5S(BD2(1),BD2(3),QM2,BK2(1),BK2(3),C1)
    IF(DABS(C1) < EPS)RETURN
    A1=A1*C1
    CALL C0T5S(BD3(1),BD3(3),QM3,BK3(1),BK3(3),C1)
    IF(DABS(C1) < EPS)RETURN
    A1=A1*C1
    CALL C0T5S(BD4(1),BD4(3),QM4,BK4(1),BK4(3),C1)
    IF(DABS(C1) < EPS)RETURN
    A1=A1*C1
    CALL RMEAJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
    IF(DABS(S) < EPS)RETURN
    CALL RMEAJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),C)
    IF(DABS(C) < EPS)RETURN
    CALL RMEAJJ(IK3(3),IK3(1),IK3(7),IK3(6),ID3(1),ID3(7),ID3(6),V)
    IF(DABS(V) < EPS)RETURN
    CALL RMEAJJ(IK4(3),IK4(1),IK4(7),IK4(6),ID4(1),ID4(7),ID4(6),Z)
    IF(DABS(Z) < EPS)RETURN
    WW=A1*S*C*V*Z/DSQRT(DBLE((IK1(7)+1)*(IK2(7)+1)*(IK3(7)+1)*(IK4(7)+1)))
    RETURN
  END SUBROUTINE GG1234

END MODULE gg
