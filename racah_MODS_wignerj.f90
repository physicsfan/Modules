MODULE wignerj
!*******************************************************************
!                                                                  *  
!   This MODULE groups together the subroutines to calculate the   *
!   the Clebsch-Gordan (3J), 6J and 9J coefficients in various     *
!   circumstances, rather than relying on calls to the subroutines *
!   in separate files.                                             *
!                                                                  *   
!                                                                  *
!   SUBROUTINES written by  G. Gaigalas             December 1993  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!   MODULE created by A. Senchuk                    February 2020  *  
!                                                                  *
!*******************************************************************
  USE constants
  USE factorials

  IMPLICIT NONE


CONTAINS

  SUBROUTINE C0T5S(Q,QM,SM,C,CM,A)
    !*******************************************************************
    !                                                                  *
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
    !                                                ---          ---  *
    !                                                I  Q  1/2  C   I  *
    !     CLEBSCH - GORDAN COEFFICIENT:              I              I  *
    !                                                I  QM  SM  CM  I  *
    !                                                ---          ---  *
    !                                                                  *
    !*******************************************************************
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    REAL(KIND=8), INTENT(IN)  :: Q, QM, SM, C, CM
    REAL(KIND=8), INTENT(OUT) :: A
    !      DIMENSION GC(2)
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER                    :: IIQ, IIC, IE, ITTK
    REAL(KIND=8), DIMENSION(2) :: GC
    !-----------------------------------------------
    GC(1)=ONE
    GC(2)=-ONE
    A=ZERO
    IIQ=TWO*Q+TENTH
    IIC=TWO*C+TENTH
    IF(ITTK(IIQ,IIC,1) == 0)RETURN
    IF(DABS(QM+SM-CM) > EPS)RETURN
    IF((HALF+TENTH) < DABS(SM))RETURN
    IF((Q+TENTH) < DABS(QM))RETURN
    IF((C+TENTH) < DABS(CM))RETURN
    IE=DABS(HALF-SM)+ONE+TENTH
    IF(DABS(Q+HALF-C) < EPS) THEN
       A=DSQRT((C+GC(IE)*CM)/(TWO*C))
    ELSE
       IF(DABS(Q-HALF-C) > EPS)RETURN
       A=-GC(IE)*DSQRT((C-GC(IE)*CM+ONE)/(TWO*C+TWO))
    ENDIF
    RETURN
  END SUBROUTINE C0T5S



  SUBROUTINE CLE0SM(Q,QM,S,C,CM,A)
    !*******************************************************************
    !                                                                  *
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
    !                                                 ---         ---  *
    !                                                 I  Q   S  C   I  *
    !     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
    !                                                 I  QM  0  CM  I  *
    !                                                 ---         ---  *
    !                                                                  *
    !*******************************************************************
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    REAL(KIND=8), INTENT(IN)  :: Q, QM, S, C, CM
    REAL(KIND=8), INTENT(OUT) :: A
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: IIQ, IIC, IIS, ITTK
    !-----------------------------------------------
    A=ZERO
    IIQ=TWO*Q+TENTH
    IIC=TWO*C+TENTH
    IIS=TWO*S+TENTH
    IF(ITTK(IIQ,IIC,IIS).EQ.0)RETURN
    IF(S.LT.EPS) THEN
       IF((Q+TENTH).LT.DABS(QM))RETURN
       IF((C+TENTH).LT.DABS(CM))RETURN
       IF(DABS(Q-C).GT.EPS)RETURN
       IF(DABS(QM-CM).GT.EPS)RETURN
       A=ONE
    ELSE
       CALL C1E0SM(Q,QM,C,CM,A)
    END IF
    RETURN
  END SUBROUTINE CLE0SM



  SUBROUTINE C1E0SM(Q,QM,C,CM,A)
    !*******************************************************************
    !                                                                  *
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
    !                                                 ---         ---  *
    !                                                 I  Q   1  C   I  *
    !     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
    !                                                 I  QM  0  CM  I  *
    !                                                 ---         ---  *
    !                                                                  *
    !*******************************************************************
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    REAL(KIND=8), INTENT(IN)  :: Q, QM, C, CM
    REAL(KIND=8), INTENT(OUT) :: A
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: IIQ, IIC, IS, IG, ITTK
    !-----------------------------------------------
    A=ZERO
    IIQ=TWO*Q+TENTH
    IIC=TWO*C+TENTH
    IF(ITTK(IIQ,IIC,2) == 0)RETURN
    IF(DABS(QM-CM) > EPS) RETURN 
    IF((Q+TENTH) < DABS(QM)) RETURN
    IF((C+TENTH) < DABS(CM)) RETURN
    IF(DABS(QM) <= EPS) THEN
       IS=Q+C+ONE+TENTH
       IF((IS/2)*2 /= IS) RETURN
    END IF
    IG=Q-C+TWO+TENTH
    IF(IG <= 0) RETURN
    IF(IG > 3) RETURN
    IF (IG == 1) THEN
       A=DSQRT(((C+CM)*(C-CM))/((TWO*C-ONE)*C))
    ELSE IF (IG == 2) THEN
       A=CM/DSQRT(C*(C+ONE))
    ELSE IF (IG == 3) THEN
       A=-DSQRT(((C+CM+ONE)*(C-CM+ONE))/((C+ONE)*(TWO*C+THREE)))
    END IF
    RETURN
  END SUBROUTINE C1E0SM



  SUBROUTINE C1E1SM(Q,QM,SM,C,CM,A)
    !*******************************************************************
    !                                                                  *
    !     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
    !                                                 ---         ---  *
    !                                                 I  Q   1  C   I  *
    !     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
    !                                                 I  QM  1  CM  I  *
    !                                                 ---         ---  *
    !                                                                  *
    !*******************************************************************
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    REAL(KIND=8), INTENT(IN)  :: Q, QM, SM, C, CM
    REAL(KIND=8), INTENT(OUT) :: A
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER                    :: IE, IIQ, IIC, ITTK
    REAL(KIND=8), DIMENSION(2) :: GC
    !-----------------------------------------------
    GC(1)=ONE
    GC(2)=-ONE
    A=ZERO
    IIQ=TWO*Q+TENTH
    IIC=TWO*C+TENTH
    IF(ITTK(IIQ,IIC,2).EQ.0)RETURN
    IF(DABS(QM+SM-CM).GT.EPS)RETURN
    IF((Q+TENTH).LT.DABS(QM))RETURN
    IF((C+TENTH).LT.DABS(CM))RETURN
    IE=0
    IF(DABS(SM-ONE).LT.EPS)IE=1
    IF(DABS(SM+ONE).LT.EPS)IE=2
    IF(IE.EQ.0)RETURN
    IF(DABS(Q+ONE-C).LT.EPS) THEN
       A=DSQRT((C+GC(IE)*CM-ONE)*(C+GC(IE)*CM)/ &
            ((TWO*C-ONE)*TWO*C))
    ELSE IF(DABS(Q-C).LT.EPS) THEN
       A=-GC(IE)*DSQRT((C+GC(IE)*CM)*(C-GC(IE)*CM+ONE)/ &
            ((C+ONE)*TWO*C))
    ELSE IF(DABS(Q-ONE-C).GT.EPS) THEN
       RETURN
    ELSE
       A=DSQRT((C-GC(IE)*CM+ONE)*(C-GC(IE)*CM+TWO)/ &
            ((TWO*C+TWO)*(TWO*C+THREE)))
    END IF
    RETURN
  END SUBROUTINE C1E1SM


  SUBROUTINE DRACAH(I, J, K, L, M, N, RAC) 
    !***********************************************************************
    !   SUBROUTINE  to calculate Racah coefficients. The arguments I, J,   *
    !   K, L, M, N should be twice their actual value. Works for integer   *
    !   and  half-integer  values of  angular momenta. The routine makes   *
    !   use of the GAM  array, thus  SUBROUTINE FACTT must be called be-   *
    !   fore this routine is used.                                         *
    !                                                                      *
    !   Written by N S Scott                    Last update: 16 Oct 1992   *
    !   The last modification made by G. Gaigalas           October 2017   *
    !                                                                      *
    !***********************************************************************
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER , INTENT(IN) :: I 
    INTEGER , INTENT(IN) :: J 
    INTEGER , INTENT(IN) :: K 
    INTEGER , INTENT(IN) :: L 
    INTEGER , INTENT(IN) :: M 
    INTEGER , INTENT(IN) :: N 
    REAL(kind=8) , INTENT(OUT) :: RAC 
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: J1, J2, J3, J4, J5, J6, J7, NUMIN, NUMAX, ICOUNT, KK, KI 
    !-----------------------------------------------
    J1 = I + J + M 
    J2 = K + L + M 
    J3 = I + K + N 
    J4 = J + L + N
    RAC = 0.0d00
    IF (2*MAX(MAX(I,J),M) - J1>0 .OR. MOD(J1,2)/=0) RETURN 
    IF (2*MAX(MAX(K,L),M) - J2>0 .OR. MOD(J2,2)/=0) RETURN 
    IF (2*MAX(MAX(I,K),N) - J3>0 .OR. MOD(J3,2)/=0) RETURN 
    IF (2*MAX(MAX(J,L),N) - J4>0 .OR. MOD(J4,2)/=0) RETURN 

    J1 = J1/2 
    J2 = J2/2 
    J3 = J3/2 
    J4 = J4/2 
    J5 = (I + J + K + L)/2 
    J6 = (I + L + M + N)/2 
    J7 = (J + K + M + N)/2 
    NUMIN = MAX(MAX(MAX(J1,J2),J3),J4) + 1 
    NUMAX = MIN(MIN(J5,J6),J7) + 1 
    RAC = 1.0D00 
    ICOUNT = 0 
    !
    IF (NUMIN /= NUMAX) THEN 
       NUMIN = NUMIN + 1 
       !
       DO KK = NUMIN, NUMAX 
          KI = NUMAX - ICOUNT 
          RAC = 1.0D00 - RAC*DBLE(KI*(J5 - KI + 2)*(J6 - KI + 2)*(J7 - KI + 2&
               ))/DBLE((KI - 1 - J1)*(KI - 1 - J2)*(KI - 1 - J3)*(KI - 1 - J4)) 
          ICOUNT = ICOUNT + 1 
       END DO
       !
       NUMIN = NUMIN - 1 
    ENDIF
    RAC = RAC*(-1.0D00)**(J5 + NUMIN + 1)*EXP((GAM(NUMIN+1)-GAM(NUMIN-J1)-GAM&
         (NUMIN-J2)-GAM(NUMIN-J3)-GAM(NUMIN-J4)-GAM(J5+2-NUMIN)-GAM(J6+2-NUMIN)&
         -GAM(J7+2-NUMIN))+(GAM(J1+1-I)+GAM(J1+1-J)+GAM(J1+1-M)-GAM(J1+2)+GAM(&
         J2+1-K)+GAM(J2+1-L)+GAM(J2+1-M)-GAM(J2+2)+GAM(J3+1-I)+GAM(J3+1-K)+GAM(&
         J3+1-N)-GAM(J3+2)+GAM(J4+1-J)+GAM(J4+1-L)+GAM(J4+1-N)-GAM(J4+2))*&
         0.5D00) 
    !
    RETURN  
  END SUBROUTINE DRACAH


  INTEGER FUNCTION IXJTIK (I, J, K, L, M, N) 
    !*******************************************************************
    !     CHECKED TRIANGULAR CONDITIONS FOR 6j COEFFICIENT             *
    !                                                                  *
    !     | I/2  J/2  K/2 |            IXJTIK=1 - IF NOT SATISFY       *
    !     | L/2  M/2  N/2 |            IXJTIK=0 - IN OVER CASES        *
    !                                                                  *
    !*******************************************************************
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER  :: I, J, K, L, M, N, ITTK
    !-----------------------------------------------
    IXJTIK = 0 
    IF (ITTK(I,J,K) == 0) RETURN  
    IF (ITTK(I,M,N) == 0) RETURN  
    IF (ITTK(L,J,N) == 0) RETURN  
    IF (ITTK(L,M,K) == 0) RETURN  
    IXJTIK = 1 
    RETURN  
  END FUNCTION IXJTIK


  SUBROUTINE SIXJ(I,J,K,L,M,N,ITIK,SI)
    !*******************************************************************  
    !                                                                  *
    !     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
    !                                                                  *
    !     | I/2  J/2  K/2 |                                            *
    !     | L/2  M/2  N/2 |          (29.1A) [J.B.77]                  *
    !                                                                  *
    !   This restructured package uses the CONTAINS keyword to group   *
    !   related subroutines together into self-contained blocks,       *
    !   rather than relying on calls to subroutines in separate files. *                                 *
    !                                                                  *
    !   Written by G. Gaigalas,                                        *
    !   Vanderbilt University,  Nashville               October  1996  *
    !   The last modification made by G. Gaigalas       October  2017  *
    !   Restructed by C. Froese Fischer                 May      2019  *
    !   Reststructured by A. Senchuk                    September 2019 *
    !                                                                  *
    !*******************************************************************
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER             :: I, J, K, L, M, N
    INTEGER, INTENT(IN) :: ITIK 
    REAL(kind=8)        :: SI 
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: IFA 
    REAL(kind=8), DIMENSION(0:4,0:4,0:4,0:4,0:4,0:4) :: RACA 
    REAL(kind=8) :: UNDEF, A 
    LOGICAL :: SAVE 
    !-----------------------------------------------
    DATA RACA/ 15625*1.D-20/  
    DATA UNDEF/ 1.D-20/  
    SI = ZERO 
    IF (ITIK /= 0) THEN 
       !
       !     CHECKED TRIANGULAR CONDITIONS
       !
       IF (IXJTIK(I,J,K,L,M,N) == 0) RETURN  
    ENDIF
    SAVE = .FALSE. 
    IF (MAX0(I,J,K,L,M,N) <= 4) THEN 
       SI = RACA(I,J,K,L,M,N) 
       IF (SI == UNDEF) THEN 
          SAVE = .TRUE. 
       ELSE 
          RETURN  
       ENDIF
    ENDIF
    !
    !     CALCULATED IN CASE WHEN ONE OF PERAMETERS EQUAL 0
    !
    IF (I*J*K*L*M*N == 0) THEN 
       IF (I == 0) THEN 
          A = DBLE((M + 1)*(K + 1)) 
          IFA = L + M + K 
       ELSE IF (J == 0) THEN 
          A = DBLE((L + 1)*(K + 1)) 
          IFA = I + M + N 
       ELSE IF (K == 0) THEN 
          A = DBLE((I + 1)*(L + 1)) 
          IFA = I + M + N 
       ELSE IF (L == 0) THEN 
          A = DBLE((J + 1)*(K + 1)) 
          IFA = I + J + K 
       ELSE IF (M == 0) THEN 
          A = DBLE((I + 1)*(K + 1)) 
          IFA = I + J + K 
       ELSE 
          A = DBLE((I + 1)*(J + 1)) 
          IFA = I + J + K 
       ENDIF
       SI = ONE/DSQRT(A) 
       IF (MOD(IFA,4) /= 0) SI = -SI 
       !
       !     THE CASE 1/2
       !
    ELSE IF (MIN0(I,J,K,L,M,N) == 1) THEN 
       IF (I == 1) THEN 
          CALL SIXJ5 (M, K, L, J, N, 0, SI) 
       ELSE IF (J == 1) THEN 
          CALL SIXJ5 (I, N, M, L, K, 0, SI) 
       ELSE IF (K == 1) THEN 
          CALL SIXJ5 (I, M, N, L, J, 0, SI) 
       ELSE IF (L == 1) THEN 
          CALL SIXJ5 (J, K, I, M, N, 0, SI) 
       ELSE IF (M == 1) THEN 
          CALL SIXJ5 (I, K, J, L, N, 0, SI) 
       ELSE 
          CALL SIXJ5 (I, J, K, L, M, 0, SI) 
       ENDIF
       !
       !     THE CASE 1
       !
    ELSE IF (MIN0(I,J,K,L,M,N) == 2) THEN 
       IF (I == 2) THEN 
          CALL SIXJ1 (M, K, L, J, N, 0, SI) 
       ELSE IF (J == 2) THEN 
          CALL SIXJ1 (I, N, M, L, K, 0, SI) 
       ELSE IF (K == 2) THEN 
          CALL SIXJ1 (I, M, N, L, J, 0, SI) 
       ELSE IF (L == 2) THEN 
          CALL SIXJ1 (J, K, I, M, N, 0, SI) 
       ELSE IF (M == 2) THEN 
          CALL SIXJ1 (I, K, J, L, N, 0, SI) 
       ELSE 
          CALL SIXJ1 (I, J, K, L, M, 0, SI) 
       ENDIF
       !
       !     THE CASE 3/2
       !
    ELSE IF (MIN0(I,J,K,L,M,N) == 3) THEN 
       IF (I == 3) THEN 
          CALL SIXJ35 (M, K, L, J, N, 0, SI) 
       ELSE IF (J == 3) THEN 
          CALL SIXJ35 (I, N, M, L, K, 0, SI) 
       ELSE IF (K == 3) THEN 
          CALL SIXJ35 (I, M, N, L, J, 0, SI) 
       ELSE IF (L == 3) THEN 
          CALL SIXJ35 (J, K, I, M, N, 0, SI) 
       ELSE IF (M == 3) THEN 
          CALL SIXJ35 (I, K, J, L, N, 0, SI) 
       ELSE 
          CALL SIXJ35 (I, J, K, L, M, 0, SI) 
       ENDIF
       !
       !     THE CASE 2
       !
    ELSE IF (MIN0(I,J,K,L,M,N) == 4) THEN 
       IF (I == 4) THEN 
          CALL SIXJ2 (M, K, L, J, N, 0, SI) 
       ELSE IF (J == 4) THEN 
          CALL SIXJ2 (I, N, M, L, K, 0, SI) 
       ELSE IF (K == 4) THEN 
          CALL SIXJ2 (I, M, N, L, J, 0, SI) 
       ELSE IF (L == 4) THEN 
          CALL SIXJ2 (J, K, I, M, N, 0, SI) 
       ELSE IF (M == 4) THEN 
          CALL SIXJ2 (I, K, J, L, N, 0, SI) 
       ELSE 
          CALL SIXJ2 (I, J, K, L, M, 0, SI) 
       ENDIF
       !
       !     THE CASE 5/2
       !
    ELSE IF (MIN0(I,J,K,L,M,N) == 5) THEN 
       CALL DRACAH (I, J, M, L, K, N, SI) 
       IF (MOD(I + J + M + L,4) /= 0) SI = -SI 
       !
       !     CASES 3
       !
    ELSE IF (MIN0(I,J,K,L,M,N) == 6) THEN 
       IF (I == 6) THEN 
          CALL SIXJ3 (M, K, L, J, N, 0, SI) 
       ELSE IF (J == 6) THEN 
          CALL SIXJ3 (I, N, M, L, K, 0, SI) 
       ELSE IF (K == 6) THEN 
          CALL SIXJ3 (I, M, N, L, J, 0, SI) 
       ELSE IF (L == 6) THEN 
          CALL SIXJ3 (J, K, I, M, N, 0, SI) 
       ELSE IF (M == 6) THEN 
          CALL SIXJ3 (I, K, J, L, N, 0, SI) 
       ELSE 
          CALL SIXJ3 (I, J, K, L, M, 0, SI) 
       ENDIF
       !
       !     THE CASE 7/2
       !
    ELSE IF (MIN0(I,J,K,L,M,N) == 7) THEN 
       CALL DRACAH (I, J, M, L, K, N, SI) 
       IF (MOD(I + J + M + L,4) /= 0) SI = -SI 
       !
       !     CASES 4
       !
    ELSE IF (MIN0(I,J,K,L,M,N) == 8) THEN 
       IF (I == 8) THEN 
          CALL SIXJ4 (M, K, L, J, N, 0, SI) 
       ELSE IF (J == 8) THEN 
          CALL SIXJ4 (I, N, M, L, K, 0, SI) 
       ELSE IF (K == 8) THEN 
          CALL SIXJ4 (I, M, N, L, J, 0, SI) 
       ELSE IF (L == 8) THEN 
          CALL SIXJ4 (J, K, I, M, N, 0, SI) 
       ELSE IF (M == 8) THEN 
          CALL SIXJ4 (I, K, J, L, N, 0, SI) 
       ELSE 
          CALL SIXJ4 (I, J, K, L, M, 0, SI) 
       ENDIF
       !
       !     THE CASE 9/2
       !
    ELSE IF (MIN0(I,J,K,L,M,N) == 9) THEN 
       CALL DRACAH (I, J, M, L, K, N, SI) 
       IF (MOD(I + J + M + L,4) /= 0) SI = -SI 
       !
       !     CALCULATED OTHER CASES
       !
    ELSE 
       CALL DRACAH(I,J,M,L,K,N,SI)
       IF (MOD(I + J + M + L,4) /= 0) SI = -SI 
    ENDIF
    IF (SAVE) RACA(I,J,K,L,M,N) = SI 
    RETURN

  CONTAINS

    SUBROUTINE SIXJ5(J,K,L,M,N,ITIK,SI) 
      !*******************************************************************
      !     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
      !                                                                  *
      !     | J/2  K/2  L/2 |                                            *
      !     | M/2  N/2  1/2 |             [B.M.X. 75]                    *
      !                                                                  *
      !*******************************************************************
      !
      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER                   :: J, K, L, M, N
      INTEGER, INTENT(IN)       :: ITIK 
      REAL(KIND=8), INTENT(OUT) :: SI 
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER      :: I1 
      REAL(KIND=8) :: AS, A, B, C, AKA 
      !-----------------------------------------------
      SI = ZERO 
      IF (ITIK /= 0) THEN 
         !
         !     CHECKED TRIANGULAR CONDITIONS
         !
         IF (IXJTIK(J,K,L,M,N,1) == 0) RETURN  
      ENDIF
      I1 = (J + K + L)/2 
      AS = DBLE(I1) 
      A = DBLE(L) 
      B = DBLE(K) 
      C = DBLE(J) 
      AKA = ONE 
      IF (MOD(I1,2) /= 0) AKA = -AKA 
      IF (K < M) THEN 
         IF (J < N) THEN
            !              M > K,  J < N.
            SI = -AKA*DSQRT((AS + TWO)*(AS - A + ONE)/((B + ONE)*(B + TWO)*(C    &
                 + ONE)*(C + TWO))) 
         ELSE IF (J > N) THEN 
            !              M > K,  J > N.
            SI = AKA*DSQRT((AS - C + ONE)*(AS - B)/((B + ONE)*(B + TWO)*C*(C +   &
                 ONE))) 
         ENDIF
      ELSE IF (K > M) THEN 
         IF (J < N) THEN 
            !             M < K,  J < N.
            SI = AKA*DSQRT((AS - C)*(AS - B + ONE)/(B*(B + ONE)*(C + ONE)*(C +   &
                 TWO))) 
         ELSE IF (J > N) THEN 
            !             M < K,  J > N.
            SI = AKA*DSQRT((AS + ONE)*(AS - A)/(B*(B + ONE)*C*(C + ONE))) 
         ENDIF
      ENDIF
      RETURN  
    END SUBROUTINE SIXJ5



    SUBROUTINE SIXJ1(I,J,K,L,M,ITIK,SI) 
      !*******************************************************************
      !     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
      !                                                                  *
      !     | I/2  J/2  K/2 |                                            *
      !     | L/2  M/2   1  |               [B.M.X.  75].                *
      !                                                                  *
      !******************************************************************* 
      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER                   :: I, J, K, L, M
      INTEGER, INTENT(IN)       :: ITIK 
      REAL(KIND=8), INTENT(OUT) :: SI 
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER      :: IFA 
      REAL(KIND=8) :: AS, AKA, A, B, C 
      !-----------------------------------------------
      SI = ZERO 
      IF (ITIK /= 0) THEN 
         !
         !     CHECKED TRIANGULAR CONDITIONS
         !
         IF (IXJTIK(I,J,K,L,M,2) == 0) RETURN  
      ENDIF
      IFA = (I + J + K)/2 
      AS = DBLE(IFA) 
      AKA = ONE 
      IF (MOD(IFA,2) /= 0) AKA = -AKA 
      A = DBLE(K) 
      B = DBLE(J) 
      C = DBLE(I) 
      IF (I < M) THEN 
         IF (J < L) THEN 
            !              M > I,   L > J.
            SI = AKA*DSQRT((AS + TWO)*(AS + THREE)*(AS - A + ONE)*(AS - A + TWO    &
                 )/((B + ONE)*(B + TWO)*(B + THREE)*(C + ONE)*(C + TWO)*(C +       &
                 THREE))) 
         ELSE IF (J == L) THEN 
            !              M > I,  L = J.
            SI = (-AKA)*DSQRT(TWO*(AS + TWO)*(AS - C)*(AS - B + ONE)*(AS - A +     &
                 ONE)/(B*(B + ONE)*(B + TWO)*(C + ONE)*(C + TWO)*(C + THREE))) 
         ELSE 
            !              M > I,  L < J.
            SI = AKA*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + ONE)*(AS - B + TWO    &
                 )/((B - ONE)*B*(B + ONE)*(C + ONE)*(C + TWO)*(C + THREE))) 
         ENDIF
      ELSE IF (I == M) THEN 
         IF (J < L) THEN 
            !              M = L,  L > J.
            SI = (-AKA)*DSQRT((AS + TWO)*(AS - C + ONE)*(AS - B)*(AS - A + ONE)    &
                 *TWO/((B + ONE)*(B + TWO)*(B + THREE)*C*(C + ONE)*(C + TWO))) 
         ELSE IF (J == L) THEN 
            !              M = I,  L = J.
            SI = (-AKA)*((B*B + C*C - A*A)*HALF + B + C - A)/DSQRT(B*(B + ONE)*    &
                 (B + TWO)*C*(C + ONE)*(C + TWO)) 
         ELSE 
            !              M = I,  L < J.
            SI = AKA*DSQRT((AS + ONE)*(AS - C)*(AS - B + ONE)*(AS - A)*TWO/((B     &
                 - ONE)*B*(B + ONE)*C*(C + ONE)*(C + TWO))) 
         ENDIF
      ELSE 
         IF (J < L) THEN 
            !              M < I,   L > J.
            SI = AKA*DSQRT((AS - C + ONE)*(AS - C + TWO)*(AS - B - ONE)*(AS - B    &
                 )/((B + ONE)*(B + TWO)*(B + THREE)*(C - ONE)*C*(C + ONE))) 
         ELSE IF (J == L) THEN 
            !              M < I,   L = J.
            SI = AKA*DSQRT((AS + ONE)*(AS - C + ONE)*(AS - B)*(AS - A)*TWO/(B*(    &
                 B + ONE)*(B + TWO)*(C - ONE)*C*(C + ONE))) 
         ELSE 
            !              M < I,   L < J.
            SI = AKA*DSQRT(AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((B - ONE)*B*(    &
                 B + ONE)*(C - ONE)*C*(C + ONE))) 
         ENDIF
      ENDIF
      RETURN  
    END SUBROUTINE SIXJ1



    SUBROUTINE SIXJ35(J,K,L,M,N,ITIK,SI) 
      !*******************************************************************
      !     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
      !                                                                  *
      !     | J/2  K/2  L/2 |                                            *
      !     | M/2  N/2  3/2 |             [B.M.X. 75]                    *
      !                                                                  *
      !*******************************************************************
      !
      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER                   :: J, K, L, M, N
      INTEGER, INTENT(IN)       :: ITIK 
      REAL(KIND=8), INTENT(OUT) :: SI 
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER :: I1 
      REAL(KIND=8) :: AS, A, B, C, AKA 
      !-----------------------------------------------
      SI = ZERO 
      IF (ITIK /= 0) THEN 
         !
         !     CHECKED TRIANGULAR CONDITIONS
         !
         IF (IXJTIK(J,K,L,M,N,3) == 0) RETURN  
      ENDIF
      I1 = (J + K + L)/2 
      AS = DBLE(I1) 
      A = DBLE(L) 
      B = DBLE(J) 
      C = DBLE(K) 
      AKA = ONE 
      IF (MOD(I1,2) /= 0) AKA = -AKA 
      IF (J - N == 3) THEN 
         ! -3
         IF (K - M == 3) THEN 
            !  I                      -3/2  -3/2
            SI = AKA*DSQRT((AS - ONE)*AS*(AS + ONE)*(AS - A - TWO)*(AS - A -       &
                 ONE)*(AS - A)/((B - TWO)*(B - ONE)*B*(B + ONE)*(C - TWO)*(C -     &
                 ONE)*C*(C + ONE))) 
         ELSE IF (M - K == 3) THEN 
            !  IV  P(12)              3/2   -3/2
            SI = AKA*DSQRT((AS - C - TWO)*(AS - C - ONE)*(AS - C)*(AS - B + ONE    &
                 )*(AS - B + TWO)*(AS - B + THREE)/((C + 1)*(C + TWO)*(C + THREE)  &
                 *(C + FOUR)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
         ELSE IF (K - M == 1) THEN 
            !  II  P(12)             -1/2   -3/2
            SI = AKA*DSQRT(THREE*AS*(AS + ONE)*(AS - A - ONE)*(AS - A)*(AS - C)    &
                 *(AS - B + ONE)/((C - ONE)*C*(C + ONE)*(C + TWO)*(B - TWO)*(B -   &
                 ONE)*B*(B + ONE))) 
         ELSE IF (M - K == 1) THEN 
            !  III P(12)              1/2   -3/2
            SI = AKA*DSQRT(THREE*(AS + ONE)*(AS - A)*(AS - C - ONE)*(AS - C)*(     &
                 AS - B + ONE)*(AS - B + TWO)/(C*(C + ONE)*(C + TWO)*(C + THREE)*  &
                 (B - TWO)*(B - ONE)*B*(B + ONE))) 
         ENDIF
      ELSE IF (N - J == 3) THEN 
         !  3
         IF (K - M == 3) THEN 
            !  IV                     -3/2   3/2
            SI = AKA*DSQRT((AS - B - TWO)*(AS - B - ONE)*(AS - B)*(AS - C + ONE    &
                 )*(AS - C + TWO)*(AS - C + THREE)/((B + ONE)*(B + TWO)*(B +       &
                 THREE)*(B + FOUR)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
         ELSE IF (M - K == 3) THEN 
            !  2       pataisyta               3/2   3/2
            SI = -AKA*DSQRT((AS + TWO)*(AS + THREE)*(AS + FOUR)*(AS - A + ONE)*    &
                 (AS - A + TWO)*(AS - A + THREE)/((B + ONE)*(B + TWO)*(B + THREE)  &
                 *(B + FOUR)*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR))) 
         ELSE IF (K - M == 1) THEN 
            !  1   P(12)   pataisytas          -1/2    3/2
            SI = -AKA*DSQRT(THREE*(AS + TWO)*(AS - A + ONE)*(AS - C + ONE)*(AS     &
                 - C + TWO)*(AS - B - ONE)*(AS - B)/((C - ONE)*C*(C + ONE)*(C +    &
                 TWO)*(B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR))) 
         ELSE IF (M - K == 1) THEN 
            !  3  P(12)     taisyta           1/2    3/2
            SI = AKA*DSQRT(THREE*(AS + TWO)*(AS + THREE)*(AS - A + ONE)*(AS - A    &
                 + TWO)*(AS - B)*(AS - C + ONE)/(C*(C + ONE)*(C + TWO)*(C +        &
                 THREE)*(B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR))) 
         ENDIF
         ! -1
      ELSE IF (J - N == 1) THEN 
         IF (K - M == 3) THEN 
            !  II                   -3/2   -1/2
            SI = AKA*DSQRT((THREE*AS*(AS + ONE)*(AS - A - ONE)*(AS - A)*(AS - B    &
                 )*(AS - C + ONE))/((B - ONE)*B*(B + ONE)*(B + TWO)*(C - TWO)*(C   &
                 - ONE)*C*(C + ONE))) 
         ELSE IF (M - K == 3) THEN 
            !  1                     3/2   -1/2
            SI = -AKA*DSQRT(THREE*(AS + TWO)*(AS - A + ONE)*(AS - B + ONE)*(AS     &
                 - B + TWO)*(AS - C - ONE)*(AS - C)/((B - ONE)*B*(B + ONE)*(B +    &
                 TWO)*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR))) 
         ELSE IF (K - M == 1) THEN 
            !  V                    -1/2   -1/2
            SI = AKA*(TWO*(AS - B)*(AS - C) - (AS + TWO)*(AS - A - ONE))*DSQRT(    &
                 (AS + ONE)*(AS - A)/((B - ONE)*B*(B + ONE)*(B + TWO)*(C - ONE)*C  &
                 *(C + ONE)*(C + TWO))) 
         ELSE IF (M - K == 1) THEN 
            !  VI P(12)              1/2   -1/2
            SI = AKA*((AS - B + TWO)*(AS - C + ONE) - TWO*(AS - A + ONE)*(AS +     &
                 ONE))*DSQRT((AS - C)*(AS - B + ONE)/(C*(C + ONE)*(C + TWO)*(C +   &
                 THREE)*(B - ONE)*B*(B + ONE)*(B + TWO))) 
         ENDIF
         ! 1
      ELSE IF (N - J == 1) THEN 
         IF (K - M == 3) THEN 
            !  III                  -3/2    1/2
            SI = AKA*DSQRT(THREE*(AS + ONE)*(AS - A)*(AS - B - ONE)*(AS - B)*(     &
                 AS - C + ONE)*(AS - C + TWO)/(B*(B + ONE)*(B + TWO)*(B + THREE)*  &
                 (C - TWO)*(C - ONE)*C*(C + ONE))) 
         ELSE IF (M - K == 3) THEN 
            !  3              pataisyta       3/2    1/2
            SI = AKA*DSQRT(THREE*(AS + TWO)*(AS + THREE)*(AS - A + ONE)*(AS - A    & 
                 + TWO)*(AS - B + ONE)*(AS - C)/(B*(B + ONE)*(B + TWO)*(B +        &
                 THREE)*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR))) 
         ELSE IF (K - M == 1) THEN 
            !  VI                   -1/2    1/2
            SI = AKA*((AS - C + TWO)*(AS - B + ONE) - TWO*(AS - A + ONE)*(AS +     &
                 ONE))*DSQRT((AS - B)*(AS - C + ONE)/(B*(B + ONE)*(B + TWO)*(B +   &
                 THREE)*(C - ONE)*C*(C + ONE)*(C + TWO))) 
         ELSE IF (M - K == 1) THEN 
            !  4      pataisyta               1/2    1/2
            SI = -AKA*(TWO*(AS - B)*(AS - C) - (AS + THREE)*(AS - A))*DSQRT((AS    &
                 + TWO)*(AS - A + ONE)/(B*(B + ONE)*(B + TWO)*(B + THREE)*C*(C     &
                 + ONE)*(C + TWO)*(C + THREE))) 
         ENDIF
      ENDIF
      RETURN  
    END SUBROUTINE SIXJ35



    SUBROUTINE SIXJ2(J,K,L,M,N,ITIK,SI) 
      !*******************************************************************
      !     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
      !                                                                  *
      !     | J/2  K/2  L/2 |                                            *
      !     | M/2  N/2   2  |             [B.M.X. 75]                    *
      !                                                                  *
      !*******************************************************************

      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER                   :: J, K, L, M, N
      INTEGER, INTENT(IN)       :: ITIK 
      REAL(KIND=8), INTENT(OUT) :: SI 
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER      :: I1 
      REAL(KIND=8) :: AS, A, B, C, AKA, X1, X2, X3 
      !-----------------------------------------------
      SI = ZERO 
      IF (ITIK /= 0) THEN 
         !
         !     CHECKED TRIANGULAR CONDITIONS
         !
         IF (IXJTIK(J,K,L,M,N,4) == 0) RETURN  
      ENDIF
      I1 = (J + K + L)/2 
      AS = DBLE(I1) 
      A = DBLE(L) 
      B = DBLE(J) 
      C = DBLE(K) 
      AKA = ONE 
      IF (MOD(I1,2) /= 0) AKA = -AKA 
      IF (J - N == 4) THEN 
         ! -2
         IF (K - M == 4) THEN 
            !  I                      -2  -2
            SI = AKA*DSQRT((AS - TWO)*(AS - ONE)*AS*(AS + ONE)/((B - THREE)*(B      &
                 - TWO)*(B - ONE)*B*(B + ONE))) 
            SI = SI*DSQRT((AS - A - THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS -      &
                 A)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
         ELSE IF (M - K == 4) THEN 
            !  V   P(12)               2  -2
            SI = AKA*DSQRT((AS - C - THREE)*(AS - C - TWO)*(AS - C - ONE)*(AS       &
                 - C)/((C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR)*(C + TWO +        &
                 THREE))) 
            SI = SI*DSQRT((AS - B + ONE)*(AS - B + TWO)*(AS - B + THREE)*(AS -      &
                 B + FOUR)/((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
         ELSE IF (K - M == 2) THEN 
            !  II   P(12)              -1   -2
            SI = AKA*TWO*DSQRT((AS - ONE)*AS*(AS + ONE)/((C - TWO)*(C - ONE)*C*     &
                 (C + ONE)*(C + TWO))) 
            SI = SI*DSQRT((AS - A - TWO)*(AS - A - ONE)*(AS - A)*(AS - C)*(AS       &
                 - B + ONE)/((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
         ELSE IF (M - K == 2) THEN 
            !  IV   P(12)               1   -2
            SI = AKA*TWO*DSQRT((AS + ONE)*(AS - A)*(AS - C - TWO)*(AS - C - ONE     &
                 )*(AS - C)/(C*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR))) 
            SI = SI*DSQRT((AS - B + ONE)*(AS - B + TWO)*(AS - B + THREE)/((B -      &
                 THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
         ELSE IF (K - M == 0) THEN 
            !  III  P(12)               0   -2
            SI = AKA*DSQRT(TWO*THREE*AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((C      &
                 - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE))) 
            SI = SI*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + ONE)*(AS - B + TWO)     &
                 /((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
         ENDIF
         !  2
      ELSE IF (N - J == 4) THEN 
         IF (K - M == 4) THEN 
            !  V                      -2   2
            SI = AKA*DSQRT((AS - B - THREE)*(AS - B - TWO)*(AS - B - ONE)*(AS       &
                 - B)/((B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR)*(B + TWO +        &
                 THREE))) 
            SI = SI*DSQRT((AS - C + ONE)*(AS - C + TWO)*(AS - C + THREE)*(AS -      &
                 C + FOUR)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
         ELSE IF (M - K == 4) THEN 
            !  1                       2   2
            SI = AKA*DSQRT((AS - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(AS      &
                 - A + ONE)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*    &
                 (B + ONE))) 
            SI = SI*DSQRT((AS + THREE + TWO)*(AS + FOUR)*(AS + THREE)*(AS + TWO     &
                 )/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE))   &
                 ) 
         ELSE IF (K - M == 2) THEN 
            !  3                      -1   2
            SI = -AKA*DSQRT((AS - A + ONE)*(AS + TWO)*(AS - B - TWO)*(AS - B -      &
                 ONE)*(AS - B)/((B + TWO + THREE)*(B + FOUR)*(B + THREE)*(B + TWO   &
                 )*(B + ONE))) 
            SI = SI*TWO*DSQRT((AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((     &
                 C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO))) 
         ELSE IF (M - K == 2) THEN 
            !  2                       1   2
            SI = -AKA*DSQRT((AS - B)*(AS - C + ONE)*(AS - A + THREE)*(AS - A +      &
                 TWO)*(AS - A + ONE)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B   &
                 + TWO)*(B + ONE))) 
            SI = SI*TWO*DSQRT((AS + FOUR)*(AS + THREE)*(AS + TWO)/((C + FOUR)*(     &
                 C + THREE)*(C + TWO)*(C + ONE)*C)) 
         ELSE IF (K - M == 0) THEN 
            !  5                        0   2
            SI = AKA*DSQRT(THREE*TWO*(AS - B)*(AS - B - ONE)*(AS - C + TWO)*(AS     &
                 - C + ONE)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*    &
                 (B + ONE))) 
            SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)     &
                 /((C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE))) 
         ENDIF
      ELSE IF (J - N == 2) THEN 
         ! -1
         IF (K - M == 4) THEN 
            !  II   P(12)              -2  -1
            SI = AKA*TWO*DSQRT((AS - ONE)*AS*(AS + ONE)/((B - TWO)*(B - ONE)*B*     &
                 (B + ONE)*(B + TWO))) 
            SI = SI*DSQRT((AS - A - TWO)*(AS - A - ONE)*(AS - A)*(AS - B)*(AS       &
                 - C + ONE)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
         ELSE IF (M - K == 4) THEN 
            !  3   P(12)                2  -1
            SI = -AKA*DSQRT((AS - A + ONE)*(AS + TWO)*(AS - C - TWO)*(AS - C -      &
                 ONE)*(AS - C)/((C + TWO + THREE)*(C + FOUR)*(C + THREE)*(C + TWO   &
                 )*(C + ONE))) 
            SI = SI*TWO*DSQRT((AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((     &
                 B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO))) 
         ELSE IF (K - M == 2) THEN 
            !  VI                       -1  -1
            SI = AKA*((A + B)*(A - B + TWO) - (C - TWO)*(C - B + TWO))/DSQRT((B     &
                 - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)) 
            SI = SI*DSQRT(AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((C - TWO)*(C       &
                 - ONE)*C*(C + ONE)*(C + TWO))) 
         ELSE IF (M - K == 2) THEN 
            !  VIII P(12)              1  -1
            SI = AKA*((A + C + FOUR)*(A - C - TWO) - (B - TWO)*(B + C + FOUR))/     &
                 DSQRT(C*(C + ONE)*(C + TWO)*(C + THREE)*(C + FOUR)) 
            SI = SI*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + ONE)*(AS - B + TWO)     &
                 /((B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO))) 
         ELSE IF (K - M == 0) THEN
            !  VII  P(12)              0  -1
            SI = AKA*HALF*((A + C + TWO)*(A - C) - B*B + FOUR)/DSQRT((C - ONE)*     &
                 C*(C + ONE)*(C + TWO)*(C + THREE)) 
            SI = SI*DSQRT(THREE*TWO*(AS + ONE)*(AS - A)*(AS - C)*(AS - B + ONE)     &
                 /((B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO))) 
         ENDIF
      ELSE IF (N - J == 2) THEN 
         !  1
         IF (K - M == 4) THEN 
            !  IV                     -2   1
            SI = AKA*TWO*DSQRT((AS + ONE)*(AS - A)*(AS - B - TWO)*(AS - B - ONE     &
                 )*(AS - B)/(B*(B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR))) 
            SI = SI*DSQRT((AS - C + ONE)*(AS - C + TWO)*(AS - C + THREE)/((C -      &
                 THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
         ELSE IF (M - K == 4) THEN 
            !  2 P(12)                 2   1
            SI = -AKA*DSQRT((AS - C)*(AS - B + ONE)*(AS - A + THREE)*(AS - A +      &
                 TWO)*(AS - A + ONE)/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C   &
                 + TWO)*(C + ONE))) 
            SI = SI*TWO*DSQRT((AS + FOUR)*(AS + THREE)*(AS + TWO)/((B + FOUR)*(     &
                 B + THREE)*(B + TWO)*(B + ONE)*B)) 
         ELSE IF (K - M == 2) THEN 
            !  VIII                   -1   1
            SI = AKA*((A + B + FOUR)*(A - B - TWO) - (C - TWO)*(B + C + FOUR))/     &
                 DSQRT(B*(B + ONE)*(B + TWO)*(B + THREE)*(B + FOUR)) 
            SI = SI*DSQRT((AS - B - ONE)*(AS - B)*(AS - C + ONE)*(AS - C + TWO)     &
                 /((C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO))) 
         ELSE IF (M - K == 2) THEN 
            !  4                       1   1
            SI = AKA*(THREE*(AS - B)*(AS - C) - (AS - A)*(AS + FOUR))/DSQRT((B      &
                 + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B) 
            SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)     &
                 /((C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C)) 
         ELSE IF (K - M == 0) THEN 
            !  6 P(12)                 0   1
            SI = -AKA*((AS - B - ONE)*(AS - C) - (AS - A)*(AS + THREE))/DSQRT((     &
                 B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B) 
            SI = SI*DSQRT(THREE*TWO*(AS - B)*(AS - C + ONE)*(AS - A + ONE)*(AS      &
                 + TWO)/((C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE))) 
         ENDIF
      ELSE IF (N - J == 0) THEN 
         ! 0
         IF (K - M == 4) THEN 
            !  III                     -2   0
            SI = AKA*DSQRT(THREE*TWO*AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((B      &
                 - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE))) 
            SI = SI*DSQRT((AS - B - ONE)*(AS - B)*(AS - C + ONE)*(AS - C + TWO)     &
                 /((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
         ELSE IF (M - K == 4) THEN 
            !  5                        2   0
            SI = AKA*DSQRT(THREE*TWO*(AS - C)*(AS - C - ONE)*(AS - B + TWO)*(AS     &
                 - B + ONE)/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*    &
                 (C + ONE))) 
            SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)     &
                 /((B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE))) 
         ELSE IF (K - M == 2) THEN 
            !  VII                     -1   0
            SI = AKA*HALF*((A + B + TWO)*(A - B) - C*C + FOUR)/DSQRT((B - ONE)*&
                 B*(B + ONE)*(B + TWO)*(B + THREE)) 
            SI = SI*DSQRT(THREE*TWO*(AS + ONE)*(AS - A)*(AS - B)*(AS - C + ONE)&
                 /((C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO))) 
         ELSE IF (M - K == 2) THEN 
            !  6                        1   0
            SI = -AKA*((AS - C - ONE)*(AS - B) - (AS - A)*(AS + THREE))/DSQRT((&
                 C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C) 
            SI = SI*DSQRT(THREE*TWO*(AS - C)*(AS - B + ONE)*(AS - A + ONE)*(AS&
                 + TWO)/((B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE))) 
         ELSE IF (K - M == 0) THEN 
            !  IX                       0   0
            X1 = (AS - B)*(AS - B - ONE)*(AS - C)*(AS - C - ONE) 
            X2 = FOUR*(AS - B)*(AS - C)*(AS - A)*(AS + TWO) 
            X3 = (AS - A)*(AS - A - ONE)*(AS + THREE)*(AS + TWO) 
            SI = AKA*(X1 - X2 + X3)/DSQRT((B - ONE)*B*(B + ONE)*(B + TWO)*(B + &
                 THREE)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE)) 
         ENDIF
      ENDIF
      RETURN  
    END SUBROUTINE SIXJ2



    SUBROUTINE SIXJ3(J, K, L, M, N, ITIK, SI) 
      !*******************************************************************
      !     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
      !                                                                  *
      !     | J/2  K/2  L/2 |                                            *
      !     | M/2  N/2   2  |             [B.M.X. 75]                    *
      !                                                                  *
      !*******************************************************************

      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER  :: J, K, L, M, N
      INTEGER, INTENT(IN) :: ITIK 
      REAL(KIND=8), INTENT(OUT) :: SI 
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER      :: I1 
      REAL(KIND=8) :: AS, A, B, C, AKA 
      !-----------------------------------------------
      SI = ZERO 
      IF (ITIK /= 0) THEN 
         !
         !     CHECKED TRIANGULAR CONDITIONS
         !
         IF (IXJTIK(J,K,L,M,N,6) == 0) RETURN  
      ENDIF
      I1 = (J + K + L)/2 
      AS = DBLE(I1) 
      A = DBLE(L) 
      B = DBLE(J) 
      C = DBLE(K) 
      AKA = ONE 
      IF (MOD(I1,2) /= 0) AKA = -AKA 
      IF (J - N == 6) THEN 
         ! -3
         IF (K - M == 6) THEN 
            !                        -3  -3
            SI = AKA*DSQRT((AS + ONE)*AS*(AS - ONE)*(AS - TWO)*(AS - THREE)*(AS&
                 - FOUR)/((B + ONE)*B*(B - ONE)*(B - TWO)*(B - THREE)*(B - FOUR)&
                 *(B - THREE - TWO))) 
            SI = SI*DSQRT((AS - A - THREE - TWO)*(AS - A - FOUR)*(AS - A - &
                 THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS - A)/((C + ONE)*C*(C - &
                 ONE)*(C - TWO)*(C - THREE)*(C - FOUR)*(C - THREE - TWO))) 
         ELSE IF (M - K == 6) THEN 
            !                        3  -3
            SI = AKA*DSQRT((AS - C)*(AS - C - ONE)*(AS - C - TWO)*(AS - C - &
                 THREE)*(AS - C - FOUR)*(AS - C - THREE - TWO)/((B + ONE)*B*(B - &
                 ONE)*(B - TWO)*(B - THREE)*(B - FOUR)*(B - THREE - TWO))) 
            SI = SI*DSQRT((AS - B + THREE + THREE)*(AS - B + THREE + TWO)*(AS&
                 - B + FOUR)*(AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((C&
                 + FOUR + THREE)*(C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR&
                 )*(C + THREE)*(C + TWO)*(C + ONE))) 
         ELSE IF (K - M == 4) THEN 
            !                       -2  -3
            SI = AKA*DSQRT(TWO*THREE*(AS - C)*(AS - B + ONE)*(AS - THREE)*(AS&
                 - TWO)*(AS - ONE)*AS*(AS + ONE)/((B - THREE - TWO)*(B - FOUR)*(&
                 B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
            SI = SI*DSQRT((AS - A - FOUR)*(AS - A - THREE)*(AS - A - TWO)*(AS&
                 - A - ONE)*(AS - A)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)&
                 *C*(C + ONE)*(C + TWO))) 
         ELSE IF (M - K == 4) THEN 
            !                        2  -3
            SI = AKA*DSQRT(TWO*THREE*(AS + ONE)*(AS - A)*(AS - C - FOUR)*(AS - &
                 C - THREE)*(AS - C - TWO)*(AS - C - ONE)*(AS - C)/((C + THREE + &
                 THREE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + &
                 ONE)*C)) 
            SI = SI*DSQRT((AS - B + THREE + TWO)*(AS - B + FOUR)*(AS - B + &
                 THREE)*(AS - B + TWO)*(AS - B + ONE)/((B - THREE - TWO)*(B - &
                 FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
         ELSE IF (K - M == 2) THEN 
            !                       -1   -3
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - C)*(AS - C - ONE)*(AS - B&
                 + TWO)*(AS - B + ONE)*(AS - TWO)*(AS - ONE)*AS*(AS + ONE)/((B&
                 - THREE - TWO)*(B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B&
                 + ONE))) 
            SI = SI*DSQRT((AS - A - THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS - &
                 A)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + &
                 THREE))) 
         ELSE IF (M - K == 2) THEN 
            !                        1   -3
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS + ONE)*AS*(AS - A - ONE)*(AS&
                 - A)*(AS - C - THREE)*(AS - C - TWO)*(AS - C - ONE)*(AS - C)/((&
                 B - THREE - TWO)*(B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B&
                 + ONE))) 
            SI = SI*DSQRT((AS - B + FOUR)*(AS - B + THREE)*(AS - B + TWO)*(AS&
                 - B + ONE)/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*&
                 (C + ONE)*C*(C - ONE))) 
         ELSE IF (M - K == 0) THEN 
            !                        0   -3
            SI = AKA*TWO*DSQRT((THREE + TWO)*(AS + ONE)*AS*(AS - ONE)*(AS - A&
                 - TWO)*(AS - A - ONE)*(AS - A)/((B - THREE - TWO)*(B - FOUR)*(B&
                 - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE))) 
            SI = SI*DSQRT((AS - C - TWO)*(AS - C - ONE)*(AS - C)*(AS - B + &
                 THREE)*(AS - B + TWO)*(AS - B + ONE)/((C + FOUR)*(C + THREE)*(C&
                 + TWO)*(C + ONE)*C*(C - ONE)*(C - TWO))) 
         ENDIF
      ELSE IF (N - J == 6) THEN 
         !  3
         IF (K - M == 6) THEN 
            !                        -3  3
            SI = AKA*DSQRT((AS - B)*(AS - B - ONE)*(AS - B - TWO)*(AS - B - &
                 THREE)*(AS - B - FOUR)*(AS - B - THREE - TWO)/((C + ONE)*C*(C - &
                 ONE)*(C - TWO)*(C - THREE)*(C - FOUR)*(C - THREE - TWO))) 
            SI = SI*DSQRT((AS - C + THREE + THREE)*(AS - C + THREE + TWO)*(AS&
                 - C + FOUR)*(AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((B&
                 + FOUR + THREE)*(B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR&
                 )*(B + THREE)*(B + TWO)*(B + ONE))) 
         ELSE IF (M - K == 6) THEN 
            !                        3   3
            SI = AKA*DSQRT((AS - A + THREE + THREE)*(AS - A + THREE + TWO)*(AS&
                 - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(AS - A + ONE)/((B&
                 + SEVEN)*(B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B + &
                 THREE)*(B + TWO)*(B + ONE))) 
            SI = SI*DSQRT((AS + SEVEN)*(AS + THREE + THREE)*(AS + THREE + TWO)*&
                 (AS + FOUR)*(AS + THREE)*(AS + TWO)/((C + SEVEN)*(C + THREE + &
                 THREE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + &
                 ONE))) 
         ELSE IF (K - M == 4) THEN 
            !                       -2   3
            SI = -AKA*DSQRT(TWO*THREE*(AS - A + ONE)*(AS + TWO)*(AS - B - FOUR)&
                 *(AS - B - THREE)*(AS - B - TWO)*(AS - B - ONE)*(AS - B)/((B + &
                 SEVEN)*(B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B + &
                 THREE)*(B + TWO)*(B + ONE))) 
            SI = SI*DSQRT((AS - C + THREE + TWO)*(AS - C + FOUR)*(AS - C + &
                 THREE)*(AS - C + TWO)*(AS - C + ONE)/((C - FOUR)*(C - THREE)*(C&
                 - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO))) 
         ELSE IF (M - K == 4) THEN 
            !                        2   3
            SI = -AKA*DSQRT(TWO*THREE*(AS - B)*(AS - C + ONE)*(AS - A + THREE&
                 + TWO)*(AS - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(AS - A&
                 + ONE)/((B + SEVEN)*(B + THREE + THREE)*(B + THREE + TWO)*(B + &
                 FOUR)*(B + THREE)*(B + TWO)*(B + ONE))) 
            SI = SI*DSQRT((AS + THREE + THREE)*(AS + THREE + TWO)*(AS + FOUR)*(&
                 AS + THREE)*(AS + TWO)/((C + THREE + THREE)*(C + THREE + TWO)*(C&
                 + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C)) 
         ELSE IF (K - M == 2) THEN 
            !                       -1   3
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - A + ONE)*(AS - A + TWO)*(&
                 AS + THREE)*(AS + TWO)*(AS - B - THREE)*(AS - B - TWO)*(AS - B&
                 - ONE)*(AS - B)/((B + SEVEN)*(B + THREE + THREE)*(B + THREE + &
                 TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE))) 
            SI = SI*DSQRT((AS - C + FOUR)*(AS - C + THREE)*(AS - C + TWO)*(AS&
                 - C + ONE)/((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + &
                 TWO)*(C + THREE))) 
         ELSE IF (M - K == 2) THEN 
            !                        1   3
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - B)*(AS - B - ONE)*(AS - C&
                 + TWO)*(AS - C + ONE)*(AS - A + FOUR)*(AS - A + THREE)*(AS - A&
                 + TWO)*(AS - A + ONE)/((B + SEVEN)*(B + THREE + THREE)*(B + &
                 THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE))) 
            SI = SI*DSQRT((AS + THREE + TWO)*(AS + FOUR)*(AS + THREE)*(AS + TWO&
                 )/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*&
                 C*(C - ONE))) 
         ELSE IF (M - K == 0) THEN 
            !                        0   3
            SI = -AKA*TWO*DSQRT((THREE + TWO)*(AS - B)*(AS - B - ONE)*(AS - B&
                 - TWO)*(AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((B + &
                 SEVEN)*(B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B + &
                 THREE)*(B + TWO)*(B + ONE))) 
            SI = SI*DSQRT((AS - A + THREE)*(AS - A + TWO)*(AS - A + ONE)*(AS + &
                 FOUR)*(AS + THREE)*(AS + TWO)/((C + FOUR)*(C + THREE)*(C + TWO)*&
                 (C + ONE)*C*(C - ONE)*(C - TWO))) 
         ENDIF
      ELSE IF (J - N == 4) THEN 
         ! -2
         IF (K - M == 6) THEN 
            !                       -3  -2
            SI = AKA*DSQRT(TWO*THREE*(AS - B)*(AS - C + ONE)*(AS - THREE)*(AS&
                 - TWO)*(AS - ONE)*AS*(AS + ONE)/((C - THREE - TWO)*(C - FOUR)*(&
                 C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
            SI = SI*DSQRT((AS - A - FOUR)*(AS - A - THREE)*(AS - A - TWO)*(AS&
                 - A - ONE)*(AS - A)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)&
                 *B*(B + ONE)*(B + TWO))) 
         ELSE IF (M - K == 6) THEN 
            !                       3   -2
            SI = -AKA*DSQRT(TWO*THREE*(AS - A + ONE)*(AS + TWO)*(AS - C - FOUR)&
                 *(AS - C - THREE)*(AS - C - TWO)*(AS - C - ONE)*(AS - C)/((C + &
                 SEVEN)*(C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C + &
                 THREE)*(C + TWO)*(C + ONE))) 
            SI = SI*DSQRT((AS - B + THREE + TWO)*(AS - B + FOUR)*(AS - B + &
                 THREE)*(AS - B + TWO)*(AS - B + ONE)/((B - FOUR)*(B - THREE)*(B&
                 - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO))) 
         ELSE IF (K - M == 4) THEN 
            !                      -2  -2
            SI = AKA*((TWO + THREE)*(AS - C)*(AS - B) - (AS + TWO)*(AS - A - &
                 FOUR))*DSQRT((AS - TWO)*(AS - ONE)*AS*(AS + ONE)/((B - FOUR)*(B&
                 - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO))) 
            SI = SI*DSQRT((AS - A - THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS - &
                 A)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + &
                 TWO))) 
         ELSE IF (M - K == 4) THEN 
            !                       2  -2
            SI = -AKA*((TWO + THREE)*(AS + ONE)*(AS - A + ONE) - (AS - C + ONE)&
                 *(AS - B + THREE + TWO))*DSQRT((AS - C - THREE)*(AS - C - TWO)*(&
                 AS - C - ONE)*(AS - C)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - &
                 ONE)*B*(B + ONE)*(B + TWO))) 
            SI = SI*DSQRT((AS - B + FOUR)*(AS - B + THREE)*(AS - B + TWO)*(AS&
                 - B + ONE)/((C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C&
                 + THREE)*(C + TWO)*(C + ONE)*C)) 
         ELSE IF (K - M == 2) THEN 
            !                       -1  -2
            SI = AKA*(TWO*(AS - C - ONE)*(AS - B) - (AS + TWO)*(AS - A - THREE)&
                 )*DSQRT(TWO*(THREE + TWO)*(AS - C)*(AS - B + ONE)*(AS - ONE)*AS*&
                 (AS + ONE)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B + &
                 ONE)*(B + TWO))) 
            SI = SI*DSQRT((AS - A - TWO)*(AS - A - ONE)*(AS - A)/((C - THREE)*(&
                 C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE))) 
         ELSE IF (M - K == 2) THEN 
            !                        1  -2
            SI = -AKA*(TWO*AS*(AS - A + ONE) - (AS - C + ONE)*(AS - B + FOUR))*&
                 DSQRT(TWO*(THREE + TWO)*(AS + ONE)*(AS - A)*(AS - C - TWO)*(AS&
                 - C - ONE)*(AS - C)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)&
                 *B*(B + ONE)*(B + TWO))) 
            SI = SI*DSQRT((AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((C + &
                 THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - &
                 ONE))) 
         ELSE IF (M - K == 0) THEN 
            !                        0  -2
            SI = -AKA*((AS - ONE)*(AS - A + ONE) - (AS - C + ONE)*(AS - B + &
                 THREE))*DSQRT(TWO*THREE*(THREE + TWO)*(AS + ONE)*AS*(AS - A - &
                 ONE)*(AS - A)/((B - FOUR)*(B - THREE)*(B - TWO)*(B - ONE)*B*(B&
                 + ONE)*(B + TWO))) 
            SI = SI*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + TWO)*(AS - B + ONE)&
                 /((C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE)*(C - &
                 TWO))) 
         ENDIF
      ELSE IF (N - J == 4) THEN 
         !  2
         IF (K - M == 6) THEN 
            !                        -3  2
            SI = AKA*DSQRT(TWO*THREE*(AS + ONE)*(AS - A)*(AS - B - FOUR)*(AS - &
                 B - THREE)*(AS - B - TWO)*(AS - B - ONE)*(AS - B)/((B + THREE + &
                 THREE)*(B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + &
                 ONE)*B)) 
            SI = SI*DSQRT((AS - C + THREE + TWO)*(AS - C + FOUR)*(AS - C + &
                 THREE)*(AS - C + TWO)*(AS - C + ONE)/((C - THREE - TWO)*(C - &
                 FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
         ELSE IF (M - K == 6) THEN 
            !                        3  2
            SI = -AKA*DSQRT(TWO*THREE*(AS - C)*(AS - B + ONE)*(AS - A + THREE&
                 + TWO)*(AS - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(AS - A&
                 + ONE)/((C + SEVEN)*(C + THREE + THREE)*(C + THREE + TWO)*(C + &
                 FOUR)*(C + THREE)*(C + TWO)*(C + ONE))) 
            SI = SI*DSQRT((AS + THREE + THREE)*(AS + THREE + TWO)*(AS + FOUR)*(&
                 AS + THREE)*(AS + TWO)/((B + THREE + THREE)*(B + THREE + TWO)*(B&
                 + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B)) 
         ELSE IF (K - M == 4) THEN 
            !                       -2  2
            SI = -AKA*((TWO + THREE)*(AS + ONE)*(AS - A + ONE) - (AS - B + ONE)&
                 *(AS - C + THREE + TWO))*DSQRT((AS - B - THREE)*(AS - B - TWO)*(&
                 AS - B - ONE)*(AS - B)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - &
                 ONE)*C*(C + ONE)*(C + TWO))) 
            SI = SI*DSQRT((AS - C + FOUR)*(AS - C + THREE)*(AS - C + TWO)*(AS&
                 - C + ONE)/((B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B&
                 + THREE)*(B + TWO)*(B + ONE)*B)) 
         ELSE IF (M - K == 4) THEN 
            !                        2  2
            SI = AKA*((TWO + THREE)*(AS - B)*(AS - C) - (AS - A)*(AS + THREE + &
                 THREE))*DSQRT((AS - A + FOUR)*(AS - A + THREE)*(AS - A + TWO)*(&
                 AS - A + ONE)/((B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*&
                 (B + THREE)*(B + TWO)*(B + ONE)*B)) 
            SI = SI*DSQRT((AS + THREE + TWO)*(AS + FOUR)*(AS + THREE)*(AS + TWO&
                 )/((C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*&
                 (C + TWO)*(C + ONE)*C)) 
         ELSE IF (K - M == 2) THEN 
            !                       -1  2
            SI = AKA*(TWO*(AS - A + TWO)*(AS + ONE) - (AS - B + ONE)*(AS - C +   &
                 FOUR))*DSQRT(TWO*(THREE + TWO)*(AS - A + ONE)*(AS + TWO)*(AS - B&
                 - TWO)*(AS - B - ONE)*(AS - B)/((B + THREE + THREE)*(B + THREE  &
                 + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B)) 
            SI = SI*DSQRT((AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((C -   &
                 THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE))) 
         ELSE IF (M - K == 2) THEN 
            !                        1  2
            SI = -AKA*(TWO*(AS - B - ONE)*(AS - C) - (AS - A)*(AS + THREE + TWO  &
                 ))*DSQRT(TWO*(THREE + TWO)*(AS - B)*(AS - C + ONE)*(AS - A +    &
                 THREE)*(AS - A + TWO)*(AS - A + ONE)/((B + THREE + THREE)*(B +  &
                 THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B)) 
            SI = SI*DSQRT((AS + FOUR)*(AS + THREE)*(AS + TWO)/((C + THREE + TWO  &
                 )*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE))) 
         ELSE IF (M - K == 0) THEN 
            !                        0  2
            SI = AKA*((AS - B - TWO)*(AS - C) - (AS - A)*(AS + FOUR))*DSQRT(TWO  &
                 *THREE*(THREE + TWO)*(AS - B)*(AS - B - ONE)*(AS - C + TWO)*(AS &
                 - C + ONE)/((B + THREE + THREE)*(B + THREE + TWO)*(B + FOUR)*(B &
                 + THREE)*(B + TWO)*(B + ONE)*B)) 
            SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)  &
                 /((C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE)*(C -   &
                 TWO))) 
         ENDIF
      ELSE IF (J - N == 2) THEN 
         ! - 1
         IF (K - M == 6) THEN 
            !                       -3   -1
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - B)*(AS - B - ONE)*(AS - C   &
                 + TWO)*(AS - C + ONE)*(AS - TWO)*(AS - ONE)*AS*(AS + ONE)/((C   &
                 - THREE - TWO)*(C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C  &
                 + ONE))) 
            SI = SI*DSQRT((AS - A - THREE)*(AS - A - TWO)*(AS - A - ONE)*(AS -   &
                 A)/((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B +  &
                 THREE))) 
         ELSE IF (M - K == 6) THEN 
            !                       3    -1
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - A + ONE)*(AS - A + TWO)*(   &
                 AS + THREE)*(AS + TWO)*(AS - C - THREE)*(AS - C - TWO)*(AS - C  &
                 - ONE)*(AS - C)/((C + SEVEN)*(C + THREE + THREE)*(C + THREE +   &
                 TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE))) 
            SI = SI*DSQRT((AS - B + FOUR)*(AS - B + THREE)*(AS - B + TWO)*(AS    &
                 - B + ONE)/((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B +    &
                 TWO)*(B + THREE))) 
         ELSE IF (K - M == 4) THEN 
            !                       -2   -1
            SI = AKA*(TWO*(AS - B - ONE)*(AS - C) - (AS + TWO)*(AS - A - THREE)  &
                 )*DSQRT(TWO*(THREE + TWO)*(AS - B)*(AS - C + ONE)*(AS - ONE)*AS*&
                 (AS + ONE)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C +   &
                 ONE)*(C + TWO))) 
            SI = SI*DSQRT((AS - A - TWO)*(AS - A - ONE)*(AS - A)/((B - THREE)*(  &
                 B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE))) 
         ELSE IF (M - K == 4) THEN 
            !                       2    -1
            SI = AKA*(TWO*(AS - A + TWO)*(AS + ONE) - (AS - C + ONE)*(AS - B +   &
                 FOUR))*DSQRT(TWO*(THREE + TWO)*(AS - A + ONE)*(AS + TWO)*(AS - C&
                 - TWO)*(AS - C - ONE)*(AS - C)/((C + THREE + THREE)*(C + THREE  &
                 + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C)) 
            SI = SI*DSQRT((AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((B -   &
                 THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE))) 
         ELSE IF (K - M == 2) THEN 
            !                        -1   -1
            SI = AKA*(TWO*THREE*(AS - C)*(AS - C - ONE)*(AS - B - ONE)*(AS - B)  &
                 - TWO*FOUR*(AS - C)*(AS - B)*(AS + TWO)*(AS - A - TWO) + (AS +  &
                 THREE)*(AS + TWO)*(AS - A - THREE)*(AS - A - TWO)) 
            SI = SI*DSQRT(AS*(AS + ONE)*(AS - A - ONE)*(AS - A)/((B - THREE)*(B  &
                 - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE)*(C - THREE)* &
                 (C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE))) 
         ELSE IF (M - K == 2) THEN 
            !                         1   -1
            SI = AKA*(TWO*THREE*(AS + ONE)*AS*(AS - A + ONE)*(AS - A + TWO) -    &
                 TWO*FOUR*(AS + ONE)*(AS - A + ONE)*(AS - C + ONE)*(AS - B +     &
                 THREE) + (AS - C + ONE)*(AS - C + TWO)*(AS - B + FOUR)*(AS - B  &
                 + THREE)) 
            SI = SI*DSQRT((AS - C - ONE)*(AS - C)*(AS - B + TWO)*(AS - B + ONE)  &
                 /((B - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B +    &
                 THREE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C +  &
                 ONE)*C*(C - ONE))) 
         ELSE IF (M - K == 0) THEN 
            !                         0   -1
            SI = AKA*TWO*(AS*(AS - ONE)*(AS - A + ONE)*(AS - A + TWO) - THREE*   &
                 AS*(AS - A + ONE)*(AS - C + ONE)*(AS - B + TWO) + (AS - C + ONE)&
                 *(AS - C + TWO)*(AS - B + THREE)*(AS - B + TWO)) 
            SI = SI*DSQRT(THREE*(AS + ONE)*(AS - A)*(AS - C)*(AS - B + ONE)/((B  &
                 - THREE)*(B - TWO)*(B - ONE)*B*(B + ONE)*(B + TWO)*(B + THREE)* &
                 (C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C - ONE)*(C - TWO)&
                 )) 
         ENDIF
      ELSE IF (N - J == 2) THEN 
         ! - 1
         IF (K - M == 6) THEN 
            !                        -3   1
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS + ONE)*AS*(AS - A - ONE)*(AS  &
                 - A)*(AS - B - THREE)*(AS - B - TWO)*(AS - B - ONE)*(AS - B)/(( &
                 C - THREE - TWO)*(C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C&
                 + ONE))) 
            SI = SI*DSQRT((AS - C + FOUR)*(AS - C + THREE)*(AS - C + TWO)*(AS    &
                 - C + ONE)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)* &
                 (B + ONE)*B*(B - ONE))) 
         ELSE IF (M - K == 6) THEN 
            !                        3    1
            SI = AKA*DSQRT(THREE*(THREE + TWO)*(AS - C)*(AS - C - ONE)*(AS - B   &
                 + TWO)*(AS - B + ONE)*(AS - A + FOUR)*(AS - A + THREE)*(AS - A  &
                 + TWO)*(AS - A + ONE)/((C + SEVEN)*(C + THREE + THREE)*(C +     &
                 THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE))) 
            SI = SI*DSQRT((AS + THREE + TWO)*(AS + FOUR)*(AS + THREE)*(AS + TWO  &
                 )/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*&
                 B*(B - ONE))) 
         ELSE IF (K - M == 4) THEN 
            !                        -2   1
            SI = -AKA*(TWO*AS*(AS - A + ONE) - (AS - B + ONE)*(AS - C + FOUR))*  &
                 DSQRT(TWO*(THREE + TWO)*(AS + ONE)*(AS - A)*(AS - B - TWO)*(AS  &
                 - B - ONE)*(AS - B)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE) &
                 *C*(C + ONE)*(C + TWO))) 
            SI = SI*DSQRT((AS - C + THREE)*(AS - C + TWO)*(AS - C + ONE)/((B +   &
                 THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B -  &
                 ONE))) 
         ELSE IF (M - K == 4) THEN 
            !                         2   1
            SI = -AKA*(TWO*(AS - C - ONE)*(AS - B) - (AS - A)*(AS + THREE + TWO  &
                 ))*DSQRT(TWO*(THREE + TWO)*(AS - C)*(AS - B + ONE)*(AS - A +    &
                 THREE)*(AS - A + TWO)*(AS - A + ONE)/((C + THREE + THREE)*(C +  &
                 THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C)) 
            SI = SI*DSQRT((AS + FOUR)*(AS + THREE)*(AS + TWO)/((B + THREE + TWO  &
                 )*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE))) 
         ELSE IF (M - K == 2) THEN 
            !                         1   1
            SI = AKA*(TWO*THREE*(AS - B)*(AS - B - ONE)*(AS - C)*(AS - C - ONE)  &
                 - TWO*FOUR*(AS - B)*(AS - C)*(AS - A)*(AS + FOUR) + (AS - A)*(  &
                 AS - A - ONE)*(AS + THREE + TWO)*(AS + FOUR)) 
            SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)  &
                 /((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B&
                 *(B - ONE)*(C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C&
                 + ONE)*C*(C - ONE))) 
         ELSE IF (K - M == 2) THEN 
            !                        -1   1
            SI = AKA*(TWO*THREE*(AS + ONE)*AS*(AS - A + ONE)*(AS - A + TWO) -    &
                 TWO*FOUR*(AS + ONE)*(AS - A + ONE)*(AS - B + ONE)*(AS - C +     &
                 THREE) + (AS - B + ONE)*(AS - B + TWO)*(AS - C + FOUR)*(AS - C  &
                 + THREE)) 
            SI = SI*DSQRT((AS - B - ONE)*(AS - B)*(AS - C + TWO)*(AS - C + ONE)  &
                 /((C - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C +    &
                 THREE)*(B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B +  &
                 ONE)*B*(B - ONE))) 
         ELSE IF (M - K == 0) THEN 
            !                         0   1
            SI = -AKA*TWO*((AS - B - ONE)*(AS - B - TWO)*(AS - C)*(AS - C - ONE  &
                 ) - THREE*(AS - B - ONE)*(AS - C)*(AS - A)*(AS + THREE) + (AS - &
                 A)*(AS - A - ONE)*(AS + FOUR)*(AS + THREE)) 
            SI = SI*DSQRT(THREE*(AS - B)*(AS - C + ONE)*(AS - A + ONE)*(AS +     &
                 TWO)/((B + THREE + TWO)*(B + FOUR)*(B + THREE)*(B + TWO)*(B +   &
                 ONE)*B*(B - ONE)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C&
                 - ONE)*(C - TWO))) 
         ENDIF
      ELSE IF (J - N == 0) THEN 
         !  0
         IF (K - M == 6) THEN 
            !                       -3    0
            SI = AKA*TWO*DSQRT((THREE + TWO)*(AS + ONE)*AS*(AS - ONE)*(AS - A    &
                 - TWO)*(AS - A - ONE)*(AS - A)/((C - THREE - TWO)*(C - FOUR)*(C &
                 - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE))) 
            SI = SI*DSQRT((AS - B - TWO)*(AS - B - ONE)*(AS - B)*(AS - C +       &
                 THREE)*(AS - C + TWO)*(AS - C + ONE)/((B + FOUR)*(B + THREE)*(B &
                 + TWO)*(B + ONE)*B*(B - ONE)*(B - TWO))) 
         ELSE IF (M - K == 6) THEN 
            !                        3    0
            SI = -AKA*TWO*DSQRT((THREE + TWO)*(AS - C)*(AS - C - ONE)*(AS - C    &
                 - TWO)*(AS - B + THREE)*(AS - B + TWO)*(AS - B + ONE)/((C +     &
                 SEVEN)*(C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C +    &
                 THREE)*(C + TWO)*(C + ONE))) 
            SI = SI*DSQRT((AS - A + THREE)*(AS - A + TWO)*(AS - A + ONE)*(AS +   &
                 FOUR)*(AS + THREE)*(AS + TWO)/((B + FOUR)*(B + THREE)*(B + TWO)*&
                 (B + ONE)*B*(B - ONE)*(B - TWO))) 
         ELSE IF (K - M == 4) THEN 
            !                       -2    0
            SI = -AKA*((AS - ONE)*(AS - A + ONE) - (AS - B + ONE)*(AS - C +      &
                 THREE))*DSQRT(TWO*THREE*(THREE + TWO)*(AS + ONE)*AS*(AS - A -   &
                 ONE)*(AS - A)/((C - FOUR)*(C - THREE)*(C - TWO)*(C - ONE)*C*(C  &
                 + ONE)*(C + TWO))) 
            SI = SI*DSQRT((AS - B - ONE)*(AS - B)*(AS - C + TWO)*(AS - C + ONE)  &
                 /((B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE)*(B -   &
                 TWO))) 
         ELSE IF (M - K == 4) THEN 
            !                        2    0
            SI = AKA*((AS - C - TWO)*(AS - B) - (AS - A)*(AS + FOUR))*DSQRT(TWO  &
                 *THREE*(THREE + TWO)*(AS - C)*(AS - C - ONE)*(AS - B + TWO)*(AS &
                 - B + ONE)/((C + THREE + THREE)*(C + THREE + TWO)*(C + FOUR)*(C &
                 + THREE)*(C + TWO)*(C + ONE)*C)) 
            SI = SI*DSQRT((AS - A + TWO)*(AS - A + ONE)*(AS + THREE)*(AS + TWO)  &
                 /((B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE)*(B -   &
                 TWO))) 
         ELSE IF (K - M == 2) THEN 
            !                       -1    0
            SI = AKA*TWO*(AS*(AS - ONE)*(AS - A + ONE)*(AS - A + TWO) - THREE*   &
                 AS*(AS - A + ONE)*(AS - B + ONE)*(AS - C + TWO) + (AS - B + ONE)&
                 *(AS - B + TWO)*(AS - C + THREE)*(AS - C + TWO)) 
            SI = SI*DSQRT(THREE*(AS + ONE)*(AS - A)*(AS - B)*(AS - C + ONE)/((C  &
                 - THREE)*(C - TWO)*(C - ONE)*C*(C + ONE)*(C + TWO)*(C + THREE)* &
                 (B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE)*(B - TWO)&
                 )) 
         ELSE IF (M - K == 2) THEN 
            !                        1    0
            SI = -AKA*TWO*((AS - C - ONE)*(AS - C - TWO)*(AS - B)*(AS - B - ONE  &
                 ) - THREE*(AS - C - ONE)*(AS - B)*(AS - A)*(AS + THREE) + (AS - &
                 A)*(AS - A - ONE)*(AS + FOUR)*(AS + THREE)) 
            SI = SI*DSQRT(THREE*(AS - C)*(AS - B + ONE)*(AS - A + ONE)*(AS +     &
                 TWO)/((C + THREE + TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C +   &
                 ONE)*C*(C - ONE)*(B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B&
                 - ONE)*(B - TWO))) 
         ELSE IF (K - M == 0) THEN 
            !                        0    0
            SI = AKA*((AS - B)*(AS - B - ONE)*(AS - B - TWO)*(AS - C)*(AS - C    &
                 - ONE)*(AS - C - TWO) - THREE*THREE*(AS - B)*(AS - B - ONE)*(AS &
                 - C)*(AS - C - ONE)*(AS - A)*(AS + TWO) + THREE*THREE*(AS - B)* &
                 (AS - C)*(AS - A)*(AS - A - ONE)*(AS + THREE)*(AS + TWO) - (AS  &
                 - A)*(AS - A - ONE)*(AS - A - TWO)*(AS + FOUR)*(AS + THREE)*(AS &
                 + TWO)) 
            SI = SI/DSQRT((B + FOUR)*(B + THREE)*(B + TWO)*(B + ONE)*B*(B - ONE  &
                 )*(B - TWO)*(C + FOUR)*(C + THREE)*(C + TWO)*(C + ONE)*C*(C -   &
                 ONE)*(C - TWO)) 
         ENDIF
      ENDIF
      RETURN  
    END SUBROUTINE SIXJ3



    SUBROUTINE SIXJ4(JC,JE,JD,JB,JF,ITIK,SI) 
      !*******************************************************************
      !     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
      !                                                                  *
      !     | JC/2  JE/2  JD/2 |                                         *
      !     | JB/2  JF/2    4  |                                         *
      !                                                                  *
      !                                                                  *
      !*******************************************************************

      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER             :: JC, JE, JD, JB, JF
      INTEGER, INTENT(IN) :: ITIK 
      REAL(KIND=8)        :: SI 
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      REAL(KIND=8) :: A, C, E, D, B, F, X1, X2, X3, S2, S3 
      !-----------------------------------------------
      SI = ZERO 
      IF (ITIK /= 0) THEN 
         !
         !     CHECKED TRIANGULAR CONDITIONS
         !
         IF (IXJTIK(JC,JE,JD,JB,JF,8) == 0) RETURN  
      ENDIF
      IF (IXJTIK(JC,JE,JD,JB,JF,6) == 0) THEN 
         CALL DRACAH (JC, JE, JF, JB, JD, 8, SI) 
         IF (MOD(JC + JE + JF + JB,4) /= 0) SI = -SI 
      ELSE 
         A = THREE 
         C = DBLE(JC)*HALF 
         E = DBLE(JE)*HALF 
         D = DBLE(JD)*HALF 
         B = DBLE(JB)*HALF 
         F = DBLE(JF)*HALF 
         X1 = A*DSQRT((A+B+E+TWO)*(A-B+E+ONE)*(A+B-E+ONE)*((-           &
              A)+B+E)*(A+C+F+TWO)*(A-C+F+ONE)*(A+C-F+ONE)*(             &
              (-A)+C+F)) 
         X2 = (A+ONE)*DSQRT((A+B+E+ONE)*(A-B+E)*(A+B-E)*((-A)           &
              + B+E+ONE)*(A+C+F+ONE)*(A-C+F)*(A+C-F)*((-A)+C            &
              + F+ONE)) 
         X3 = (TWO*A+ONE)*(TWO*(A*(A+ONE)*D*(D+ONE)-B*(B+ONE)*C*(C+     &
              ONE)-E*(E+ONE)*F*(F+ONE))+(A*(A+ONE)-B*(B+ONE)-E*(E       &
              +ONE))*(A*(A+ONE)-C*(C+ONE)-F*(F+ONE))) 
         IF (DABS(X2) < EPS) THEN 
            S2 = ZERO 
         ELSE 
            CALL SIXJ2 (JC, JE, JD, JB, JF, 0, S2) 
         ENDIF
         IF (DABS(X3) < EPS) THEN 
            S3 = ZERO 
         ELSE 
            CALL SIXJ3 (JC, JE, JD, JB, JF, 0, S3) 
         ENDIF
         SI = (X3*S3 - X2*S2)/X1 
      ENDIF
      RETURN  
    END SUBROUTINE SIXJ4

  END SUBROUTINE SIXJ


  SUBROUTINE NINEJ(J1,J2,J3,L1,L2,L3,K1,K2,K3,I,INN,AA)
    !*******************************************************************
    !                                                                  *
    !     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
    !                                                                  *
    !     |  J1/2  J2/2  J3/2 |                                        *
    !     |  L1/2  L2/2  L3/2 |                                        *
    !     |  K1/2  K2/2  K3/2 |                                        *
    !                                                                  *
    !*******************************************************************
    !
    IMPLICIT NONE
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    INTEGER  :: J1, J2, J3, L1, L2, L3, K1, K2, K3
    INTEGER, INTENT(IN)  :: I 
    INTEGER, INTENT(OUT) :: INN
    REAL(kind=8)  :: AA 
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    INTEGER :: N1, N2, N3, N4, N5, N6, MAX_, MIN_, IX, ITTK 
    REAL(kind=8) :: S1, S2, S3, X 
    !-----------------------------------------------
    IF (I == 1) THEN 
       INN = 0 
       IF (ITTK(J1,J2,J3) == 0) RETURN  
       IF (ITTK(L1,L2,L3) == 0) RETURN  
       IF (ITTK(K1,K2,K3) == 0) RETURN  
       IF (ITTK(J1,L1,K1) == 0) RETURN  
       IF (ITTK(J2,L2,K2) == 0) RETURN  
       IF (ITTK(J3,L3,K3) == 0) RETURN  
       INN = 1 
       RETURN  
    ENDIF
    IF (J1*J2*J3*L1*L2*L3*K1*K2*K3 == 0) THEN 
       INN = 1 
       CALL NINEJ0 (J1, J2, J3, L1, L2, L3, K1, K2, K3, AA) 
    ELSE 
       N1 = IABS(J1 - K3) 
       N2 = IABS(L3 - J2) 
       N3 = IABS(L1 - K2) 
       N4 = IABS(J2 - L3) 
       N5 = IABS(K2 - L1) 
       N6 = IABS(J1 - K3) 
       MAX_ = MAX0(N1,N2,N3,N4,N5,N6) 
       N1 = J1 + K3 
       N2 = L3 + J2 
       N3 = J2 + L3 
       N4 = K2 + L1 
       N5 = J1 + K3 
       N6 = L1 + K2 
       MIN_ = MIN0(N1,N2,N3,N4,N5,N6) 
       INN = 1 
       AA = ZERO 
       DO IX = MAX_, MIN_, 2 
          CALL SIXJ (J1, J2, J3, L3, K3, IX, 0, S1) 
          CALL SIXJ (L1, L2, L3, J2, IX, K2, 0, S2) 
          CALL SIXJ (K1, K2, K3, IX, J1, L1, 0, S3) 
          X = S1*S2*S3*DBLE(IX + 1) 
          IF (MOD(IX,2) /= 0) X = -X 
          AA = X + AA 
       END DO
    ENDIF
    RETURN  

  CONTAINS

    SUBROUTINE NINEJ0(J1,J2,J3,L1,L2,L3,K1,K2,K3,AA) 
      !*******************************************************************
      !    THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT          *
      !                                                                  *
      !     |  J1/2  J2/2  J3/2 |                                        *
      !     |  L1/2  L2/2  L3/2 |                                        *
      !     |  K1/2  K2/2    0  |                                        *
      !                                                                  *
      !*******************************************************************
      !
      IMPLICIT NONE
      !-----------------------------------------------
      !   D u m m y   A r g u m e n t s
      !-----------------------------------------------
      INTEGER  :: J1, J2, J3, L1, L2, L3, K1, K2, K3
      REAL(KIND=8), INTENT(OUT) :: AA 
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      INTEGER      :: IFA 
      REAL(KIND=8) :: A, B 
      !-----------------------------------------------
      IF (J1 == 0) THEN 
         CALL SIXJ (L2, K2, J2, K3, L3, L1, 0, A) 
         B = DBLE((J2 + 1)*(L1 + 1)) 
         IFA = K2 + J2 + L3 + L1 
      ELSE IF (J2 == 0) THEN 
         CALL SIXJ (L3, K3, J3, K1, L1, L2, 0, A) 
         B = DBLE((J3 + 1)*(L2 + 1)) 
         IFA = K3 + J3 + L1 + L2 
      ELSE IF (J3 == 0) THEN 
         CALL SIXJ (L1, K1, J1, K2, L2, L3, 0, A) 
         B = DBLE((J1 + 1)*(L3 + 1)) 
         IFA = K1 + J1 + L2 + L3 
      ELSE IF (L1 == 0) THEN 
         CALL SIXJ (K2, J2, L2, J3, K3, K1, 0, A) 
         B = DBLE((L2 + 1)*(K1 + 1)) 
         IFA = J2 + L2 + K3 + K1 
      ELSE IF (L2 == 0) THEN 
         CALL SIXJ (K3, J3, L3, J1, K1, K2, 0, A) 
         B = DBLE((L3 + 1)*(K2 + 1)) 
         IFA = J3 + L3 + K1 + K2 
      ELSE IF (L3 == 0) THEN 
         CALL SIXJ (K1, J1, L1, J2, K2, K3, 0, A) 
         B = DBLE((L1 + 1)*(K3 + 1)) 
         IFA = J1 + L1 + K2 + K3 
      ELSE IF (K1 == 0) THEN 
         CALL SIXJ (J2, J3, J1, L3, L2, K2, 0, A) 
         B = DBLE((J1 + 1)*(K2 + 1)) 
         IFA = J3 + J1 + L2 + K2 
      ELSE IF (K2 == 0) THEN 
         CALL SIXJ (J3, J1, J2, L1, L3, K3, 0, A) 
         B = DBLE((J2 + 1)*(K3 + 1)) 
         IFA = J1 + J2 + L3 + K3 
      ELSE IF (K3 == 0) THEN 
         CALL SIXJ (J1, J2, J3, L2, L1, K1, 0, A) 
         B = DBLE((J3 + 1)*(K1 + 1)) 
         IFA = J2 + J3 + L1 + K1 
      ELSE 
         A = ZERO 
         B = ONE 
      ENDIF
      AA = A/DSQRT(B) 
      IF (MOD(IFA,4) /= 0) AA = -AA 
      RETURN  
    END SUBROUTINE NINEJ0

  END SUBROUTINE NINEJ
  
END MODULE wignerj
